import ROOT
from ROOT import RooFit, RooRealVar, RooGaussian, RooGenericPdf, RooArgList, RooArgSet, RooDataSet, RooProdPdf

class Constraint(object):
    
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
    
    
##class POI(object):
##    def __init__(self, name, minimum, maximum, label=None):
##        self.name = name
##        self.pdf = ROOT.RooRealVar(name, name, 0, minimum, maximum)
##        if not label:
##            label = name
##        self.label = label
        

class CouplingsFitter2(object):
    
    def __init__(self):
##        self.poi = []
        self.poiLabels = []
##        self.pdfs = []
##        self.obs = []
##        self.cyields = []
##        self.cwidths = []
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
        
    def addPOI(self,poi,label='',minimum =-0.3 ,maximum = 0.3):
        '''Add a parameter of interest,

        Example:
        
        addPOI('Z','Z',-0.05,0.05)
        ->
        Z[0,-0.05,0.05]
        # adds variable gZ with value 0,
        # allow it to scale between 0.95 and 1.05      
        '''
        self.poi[poi] = ROOT.RooRealVar(poi, poi, 0, minimum, maximum)
        if label == '':
            label = poi
        self.poilabels[poi] = label

    def createWidthDeviation(self):
        '''Compute the width according to the specified coupling scaling factors. 
        '''
        expr='0'
        sumBR = sum(self.BR.values())
        pwidths = []
        for coupling,br in self.BR.iteritems():
            pwidth = None
            if coupling in self.poi:
                pwidth = str(br/sumBR) + "*(1+"+coupling+")*(1+"+coupling+")"
            else:
                # using sm partial width
                pwidth = str(br/sumBR)
            pwidths.append(pwidth)
        expr = '+'.join(pwidths)
        if 'inv' in self.poi:
            # ? 
            expr='('+expr+')/(1.0-inv*inv)'
        else:
            # setting invisible width to 0.
            expr='('+expr+')'
            
        dependentlist = RooArgList()    
        for dep in self.poi.values():
            dependentlist.add(dep)            
        self.width = RooGenericPdf('width', 'width', expr, dependentlist)
        
        
    def addConstraint(self, name, expr, deplist, mean, sigma):
        '''Add a constraint on one of the observables
        
        For example, for WH inclusive: 
        
        expr::Wh('(1+W)*(1+W)',W)
        #W is kW. the WH yield depends quadratically on kW
        #Wh is a pdf on W
        WhObs[1]
        #create a variable: observed Wh scaling factor w/r standard model, with value 1.
        RooGaussian::WhConstraint(WhObs,Wh,0.004)
        #create a Gaussian on variable WhObs, mean Wh, sigma 0.004
        
        the fit varies W, thus moving the mean of the Gaussian.
        the value of the gaussian is then evaluated at WhObs.
        '''
##        fcmd0 = 'expr::'+name+"('"+equation+"',"+dependents+')'
##        self.w.factory(fcmd0)
##        fcmd1 = name+'Obs['+str(number)+']'
##        self.w.factory(fcmd1)
##        self.obs.append(name+'Obs')
##        fcmd2 = 'RooGaussian::'+name+'Constraint('+name+'Obs,'+name+','+str(error)+')'
##        self.w.factory(fcmd2)
##        print 'Add constraint', '*' * 50
##        print fcmd0
##        print fcmd1
##        print fcmd2
##        self.pdfs.append(name+'Constraint')
##        self.cyields.append(name)
##        self.cwidths.append(error)
        depnames = deplist.split(',')
        deps = []
        for dep in depnames:
            if dep == 'width':
                deps.append(self.width)
            else:
                deps.append(self.poi[dep])
        self.constraint[name] = Constraint(name, expr, deps, mean, sigma)
        
    def info(self):
        for name, constraint in sorted(self.constraint.iteritems()):
            constraint.info()
        

##    def addUniformConstraint(self,name,equation,dependents):
##        self.w.factory('expr::'+name+"('"+equation+"',"+dependents+')')
##        self.w.factory('RooUniform::'+name+'Constraint('+name+')')
##        self.pdfs.append(name+'Constraint')
##
    
    
    def fit(self):
##        self.w.factory('PROD::model('+','.join(self.pdfs)+')')
##        #create dataset
##        argset=ROOT.RooArgSet('set')
##        for obs in self.obs:
##            argset.add(self.w.var(obs))
##        self.data = ROOT.RooDataSet('data','data',argset)
##        self.data.add(argset)
##        self.data.Print()
##        self.fitResult=self.w.pdf('model').fitTo(self.data,
##                                                 ROOT.RooFit.PrintLevel(3),
##                                                 ROOT.RooFit.Optimize(1),
##                                                 ROOT.RooFit.Hesse(1),
##                                                 ROOT.RooFit.Minos(1),
##                                                 ROOT.RooFit.Strategy(2),
##                                                 ROOT.RooFit.Save(1))
##        #        self.w.pdf('model').fitTo(data)
##
        pdfs = RooArgList()
        obsvars = RooArgSet('set')
        for constraint in self.constraint.values():
            pdfs.add(constraint.pdf_constraint)
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

    def createSummary(self):
        #sample the covariance matrix for the width
        ROOT.gStyle.SetOptTitle(0)
        graph = ROOT.TGraphAsymmErrors(len(self.poi)+2)
        
        order_BR = ['Z', 'W', 'b', 'c', 'g', 'tau', 'mu', 'gamma', 'inv']

        for i, poiname in enumerate(order_BR):
            poi = self.poi[poiname]
            graph.SetPoint(i, i+0.5, poi.getVal())
            graph.SetPointError(i, 0.0, 0.0, -poi.getAsymErrorLo(), poi.getAsymErrorHi())
            
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
        histo = self.h_width
        histo2 = self.constraint['Zh'].pulls
        
        graph.SetMarkerStyle(20)
        graph.SetLineWidth(3)
        c = ROOT.TCanvas('canvas','')
        c.cd()
        obj=[graph]
        graph.Draw("AP")
        graph.GetYaxis().SetTitle("68% CL on d(A) ")
        graph.GetXaxis().SetNdivisions(0)
        l=ROOT.TLine()
        l.SetLineColor(ROOT.kRed)
        l.SetLineWidth(3)
        l.DrawLine(0.0,0.0,len(self.poi)+1.5,0)
        obj.append(l)

        graph.SetPoint(len(self.poi),len(self.poi)+0.5,0.0)
        graph.SetPointError(len(self.poi),0.0,0.0,histo.GetRMS()/histo.GetMean(),
                            histo.GetRMS()/histo.GetMean())

        for i, poiname in enumerate(order_BR+['#Gamma_{T}']):
            label = self.poilabels.get(poiname, poiname)
            obj.append(ROOT.TLatex(i+0.5,0.95*graph.GetYaxis().GetXmin(),label))
            obj[-1].Draw()

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
        for name, poi in self.poi.iteritems():
            poiLabel = self.poilabels.get(name, name)
            print poiLabel+':   ('+str(poi.getAsymErrorLo())+','+str(poi.getAsymErrorHi())+')'

        cprime = ROOT.TCanvas('canvas2','')
        cprime.cd()
        histo.GetXaxis().SetTitle("#Gamma_{T}")
        histo.GetYaxis().SetTitle("N toys")
        histo.Draw()
        obj.append(histo)
        print 'Relative error on the total width ',histo.GetRMS()/histo.GetMean()
        print 'Please check the histogram to see that the dist is Gaussian. If not the fit is biased'
        print 'The fit can be biased when floating the width sometimes.'
        
        c2 =  ROOT.TCanvas('canvas3','')
        c2.cd()
        histo2.Draw()
        obj.append(histo2)

        return c, c2, cprime, obj

        

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




        
        
