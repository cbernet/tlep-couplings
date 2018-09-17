import math
import ROOT
from ROOT import RooFit, RooRealVar, RooGaussian, RooUniform, RooGenericPdf, RooArgList, RooArgSet, RooDataSet, RooProdPdf
    

class UniformConstraint(object):
    
    def __init__(self, name, expr, deps):
        self.name = name
        self.expr = expr
        self.deps = deps
        dependentlist = RooArgList()
        for dep in deps:
            dependentlist.add(dep)                
        self.pdf_yield =  RooGenericPdf(name, name, expr, dependentlist)
        pname = name + "Constraint"
        self.pdf_constraint= RooUniform(
            pname, pname, RooArgSet(self.pdf_yield)
        )
        self.pulls = ROOT.TH1F('pulls_'+name, name,
                                   1000, -10, 10)
                
    def info(self):
        print self.name
        self.pdf_yield.Print()
        self.pdf_constraint.Print()
        
    def fill_pull(self):
        pass
    
class GaussianConstraint(object):
    
    def __init__(self, name, expr, deps, mean, sigma):
        self.name = name
        self.expr = expr
        self.deps = deps
        self.mean = mean
        self.sigma = sigma
        # pdf of yield vs coupling modificators
        dependentlist = RooArgList()
        for dep in deps:
            dependentlist.add(dep)        
        self.pdf_yield = RooGenericPdf(name, name, expr, dependentlist)
        # observable (measurement)
        obsname = name + 'Obs'
        self.var_obs = RooRealVar(obsname, obsname, mean)
        # width of Gaussian pdf
        sname = name + 'Sigma'
        self.var_sigma = RooRealVar(sname, sname, sigma)
        # Gaussian pdf
        gname = name + "Constraint"
        self.pdf_constraint= RooGaussian(
            gname, gname, self.var_obs,self.pdf_yield,self.var_sigma
        )
        self.pulls = ROOT.TH1F('pulls_'+name, name,
                               1000, -10, 10)

    def fill_pull(self):      
        pull = (self.pdf_yield.getVal() - 1)/ self.var_sigma.getVal()
        self.pulls.Fill(pull)

    def info(self):
        print self.name
        self.pdf_yield.Print()
        self.var_obs.Print()
        self.var_sigma.Print()
        self.pdf_constraint.Print()
    
    
class CouplingsFitter2(object):
    
    def __init__(self):
        self.poiLabels = []
        self.BR = dict(
            b = 0.577, 
            tau = 0.063, 
            mu = 2.2e-4, 
            c = 2.91e-2, 
            g = 8.57e-2, 
            gamma = 3.82e-3, 
            W = 0.215, 
            Z = 0.0264, 
            t = 0.0
        )
        
        self.poi = dict()
        self.poilabels = dict()
        self.constraint = dict()
        self.canvases = dict()
        self._keep = []
        
        
    def addPOI(self,poi,label='',minimum =-0.3 ,maximum = 0.3):
        '''Add a parameter of interest.

        Example:
        
        addPOI('Z','Z',-0.05,0.05)
        ->
        Z[0,-0.05,0.05]
        # adds variable Z with value 0,
        # allow it to scale between -0.05 and 0.05
        '''
        self.poi[poi] = ROOT.RooRealVar(poi, poi, 0, minimum, maximum)
        if label == '':
            label = poi
        self.poilabels[poi] = label

    def createWidthDeviation(self):
        '''Compute the width deviation (denoted \kappa_H^2 by M.Peskin in arxiv 1312.4974).
        
        Note that we fit an additive modification of the coupling: (1 + dx)
        is therefore equal to kappa_x
        '''
        expr='0'
        sumBR = sum(self.BR.values())
        pwidths = []
        for dcoupling,br in self.BR.iteritems():
            pwidth = None
            if dcoupling in self.poi:
                pwidth = str(br/sumBR) + "*(1+"+dcoupling+")*(1+"+dcoupling+")"
            else:
                # using sm partial width
                pwidth = str(br/sumBR)
            pwidths.append(pwidth)
        expr = '+'.join(pwidths)
        if 'inv' in self.poi:
            expr='('+expr+')/(1.0-inv)'
        else:
            # setting invisible width to 0.
            expr='('+expr+')'            
        dependentlist = RooArgList()    
        for dep in self.poi.values():
            dependentlist.add(dep)            
        self.width = RooGenericPdf('width', 'width', expr, dependentlist)
                
           
    def addConstraint(self, name, expr, deplist, mean, sigma):
        '''Add a constraint on one of the observables
        
        For example, for ZH inclusive: 
        
        f.addConstraint('Zh','(1+Z)*(1+Z)','Z',1,0.004)  

        Z is an additive modification of the gZ coupling w/r to the standard model,
        so 1+Z = \kappa_Z
        
        Zh is the pdf of the ratio of the yield w/r to the one expected in the standard model.
        This pdf depends on Z, as (1+Z)*(1+Z).
        
        ZhObs is the measured value, here 1 so we assume that the SM yield is observed.
        
        The fit varies the parameter of interest Z, thus modifying the pdf,
        while ZhObs is fixed at 1. The likelihood of each value of Z is evaluated at ZhObs on the pdf.
        '''
        print 'constraint:', name
        print expr
        print deplist
        print mean, '+/-', sigma
        deps = self._getdeps(deplist)
        self.constraint[name] = GaussianConstraint(name, expr, deps, mean, sigma)

    def _getdeps(self, deplist):
        depnames = deplist
        try:
            depnames = deplist.split(',')
        except:
            pass
        deps = []
        for dep in depnames:
            if dep == 'width':
                deps.append(self.width)
            else:
                deps.append(self.poi[dep])
        return deps

    def addChannel(self, name, mean, sigma, prod, decay=None):
        expr_prod = '(1+{prod})*(1+{prod})'.format(prod=prod)
        expr = expr_prod
        variables = [prod]
        if decay:
            expr_decay = '(1+{decay})*(1+{decay})'.format(decay=decay)    
            expr = '{expr_prod}*{expr_decay}/width'.format(
                expr_prod=expr_prod,
                expr_decay=expr_decay,
            )
            variables.extend([decay, 'width'])
        deplist = ','.join(variables)
        self.addConstraint(name, expr, deplist, mean, sigma)
    
    def addUniformConstraint(self,name, expr):
        '''Adds a uniform constraint with pdf name, and expression expr.
        
        For example:
        
        f.addPOI('inv','inv', 0, 0.01)
        f.addUniformConstraint('Zhinv','inv') ####->Means free floating

        inv (the invisible BR) is free to float between 0 and 1%
        '''
        deps = self._getdeps([expr])        
        self.constraint[name] = UniformConstraint(name, expr, deps)

        
    def info(self):
        for name, constraint in sorted(self.constraint.iteritems()):
            constraint.info()
    

    def fit(self):
        pdfs = RooArgList()
        obsvars = RooArgSet('set')
        for constraint in self.constraint.values():
            pdfs.add(constraint.pdf_constraint)
            if hasattr(constraint, 'var_obs'):
                obsvars.add(constraint.var_obs)
        self.model = RooProdPdf('model', 'model', pdfs)
        self.data = RooDataSet('data', 'data', obsvars)
        self.data.add(obsvars)
        self.data.Print()
        self.fit_result = self.model.fitTo(self.data,
                                           ROOT.RooFit.PrintLevel(3),
                                           ROOT.RooFit.Optimize(1),
                                           ROOT.RooFit.Hesse(1),
                                           ROOT.RooFit.Minos(1),
                                           ROOT.RooFit.Strategy(2),
                                           ROOT.RooFit.Save(1))        

    def canvas(self, name, *args):
        canvas = self.canvases.setdefault(name, ROOT.TCanvas(name, name, *args))
        canvas.cd()
        return canvas

    def keep(self, obj):
        self._keep.append(obj)
        return obj

    def createSummary(self):
        #sample the covariance matrix for the width
        # ROOT.gStyle.SetOptTitle(0)
        ROOT.gStyle.SetStatW(0.4);                
        ROOT.gStyle.SetStatH(0.4);                
        
        self.graph_couplings = ROOT.TGraphAsymmErrors(len(self.poi)+2)
        
        order_BR = ['Z', 'W', 'b', 'c', 'g', 'tau', 'mu', 'gamma', 'inv']
        
        for br in order_BR:
            if not self.poi.get(br, None):
                order_BR.remove(br)

        for i, poiname in enumerate(order_BR):
            poi = self.poi.get(poiname)
            self.graph_couplings.SetPoint(i, i+0.5, poi.getVal())
            self.graph_couplings.SetPointError(i, 0.0, 0.0, -poi.getAsymErrorLo(), poi.getAsymErrorHi())
            
        print 'Sampling the covariance matrix to propagate error on width'

        self.h_width = ROOT.TH1F('h_width','width',1000,0.5,1.5)
        ntoys = 10000
        for i in range(ntoys):
            randomizedPars = self.fit_result.randomizePars()
            for j in range(0,randomizedPars.getSize()):
                self.poi[randomizedPars.at(j).GetName()].setVal(randomizedPars.at(j).getVal())
            self.h_width.Fill(self.width.getVal())
            for cstr in self.constraint.values():
                cstr.fill_pull()        
        self.graph_couplings.SetMarkerStyle(20)
        self.graph_couplings.SetLineWidth(3)
        can_couplings = self.canvas('couplings')
        self.graph_couplings.Draw("AP")
        self.graph_couplings.GetYaxis().SetTitle("68% CL on d(A) ")
        self.graph_couplings.GetXaxis().SetNdivisions(0)
        l= self.keep(ROOT.TLine())
        l.SetLineColor(ROOT.kRed)
        l.SetLineWidth(3)
        l.DrawLine(0.0,0.0,len(self.poi)+1.5,0)

        self.graph_couplings.SetPoint(len(self.poi),len(self.poi)+0.5,0.0)
        self.graph_couplings.SetPointError(len(self.poi),0.0,0.0,self.h_width.GetRMS()/self.h_width.GetMean(),
                            self.h_width.GetRMS()/self.h_width.GetMean())

        for i, poiname in enumerate(order_BR+['#Gamma_{T}']):    
            label = self.poilabels.get(poiname, poiname)
            tex = self.keep(ROOT.TLatex(i+0.5,0.95*self.graph_couplings.GetYaxis().GetXmin(),label))
            tex.Draw()

        print """
###############################################################
###############################################################
###############################################################
                         RESULTS
###############################################################
###############################################################
############################################################### 
              """

        print 'RESULTS FOR THE CONFIDENCE INTERVALS------>'
        for name in order_BR:
            poi = self.poi[name]
            poiLabel = self.poilabels.get(name, name)
            mind = poi.getAsymErrorLo() * 100
            maxd =  poi.getAsymErrorHi() * 100
            avd = abs(maxd - mind) / 2.
            # print poiLabel+':   ('+str(poi.getAsymErrorLo())+','+str(poi.getAsymErrorHi())+'), ' + str(avd)
            print '{label:10}:\t{mind:5.3f}%\t{maxd:5.3f}%\t{avd:5.3f}%'.format(label=poiLabel,
                                                                                mind=mind, 
                                                                                maxd=maxd, 
                                                                                avd=avd)

        can_gamma = self.canvas('gamma')
        self.h_width.GetXaxis().SetTitle("#Gamma_{T}")
        self.h_width.GetYaxis().SetTitle("N toys")
        self.h_width.Draw()
        print 'Relative error on the total width ',self.h_width.GetRMS()/self.h_width.GetMean()
        print 'Please check the histogram to see that the dist is Gaussian. If not the fit is biased'
        print 'The fit can be biased when floating the width sometimes.'

        can_pulls = self.canvas('pulls', 1000, 1000)
        npulls = len(self.constraint)
        nxy = int( math.ceil(math.sqrt(npulls)) )
        can_pulls.Divide(nxy, nxy)
        for i, c in enumerate(self.constraint.values()):
            can_pulls.cd(i+1)
            c.pulls.Draw()

##
##    def contour(self,var1,var2,n1=1,n2=2,n3=0,n4=0,n5=0):
##        self.w.factory('PROD::model('+','.join(self.pdfs)+')')
##        #create dataset
##        argset=ROOT.RooArgSet('set')
##        for obs in self.obs:
##            argset.add(self.w.var(obs)
##        data = ROOT.RooDataSet('data','data',argset)
##        data.add(argset)
##        data.Print()
##
##        nll = self.w.pdf('model').createNLL(data,ROOT.RooFit.Optimize(1),ROOT.RooFit.Offset(0))
##        minuit =ROOT.RooMinuit(nll)
##        return minuit.contour(self.w.var(var1),self.w.var(var2),n1,n2,n3,n4,n5)
##    #        self.w.pdf('model').fitTo(data)
##            
##    
##    def w(self):
##        return self.w
##        




        
        
