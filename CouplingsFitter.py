from GenericFitter import *
from math import sqrt
class CouplingsFitter(GenericFitter):
    def __init__(self):
        self.BR=dict()
        self.BR['b'] = 0.577
        self.BR['tau'] = 0.063
        self.BR['mu'] = 2.2e-4
        self.BR['c'] = 2.91e-2
        self.BR['g'] = 8.57e-2
        self.BR['gamma'] = 3.82e-3
        self.BR['W'] = 0.215
        self.BR['Z'] = 0.0264
        self.BR['t'] = 0.0

        super(CouplingsFitter,self).__init__()



    def createWidthDeviation(self):
        expr='0'
        sumBR=0.0
        for coupling,br in self.BR.iteritems():
            sumBR=sumBR+br
        for coupling,br in self.BR.iteritems():
            if coupling in self.poi:
                expr=expr+'+'+str(br/sumBR)+"*(1+"+coupling+")*(1+"+coupling+")"
            else:    
                expr=expr+'+'+str(br/sumBR)
        if 'inv' in self.poi:
            expr='('+expr+')/((1.0-inv*inv))'
        else:    
            expr='('+expr+')'

        dependents =   ",".join(self.poi)
        print 'Width expression defined as ',expr,' with dependends'+dependents
 
        self.w.factory('expr::width'+"('"+expr+"',"+dependents+')')



    def createPOIs(self,inv = True):
        self.addPOI('Z','Z',-0.1,0.000001)
        self.addPOI('W','W',-0.1,0.000001)
        self.addPOI('b')
        self.addPOI('c')
        self.addPOI('g')
        self.addPOI('tau','#tau')
        self.addPOI('mu','#mu')
        self.addPOI('t')
        self.addPOI('gamma','#gamma')
        if inv:
            self.addPOI('inv','#Gamma_{inv}',0,0.1)
        self.createWidthDeviation()    



    def createSummary(self):
        #sample the covariance matrix for the width
        ROOT.gStyle.SetOptTitle(0)
        graph = ROOT.TGraphAsymmErrors(len(self.poi)+2)

        for i,poi in enumerate(self.poi):
            graph.SetPoint(i,i+0.5,self.w.var(poi).getVal())
            graph.SetPointError(i,0.0,0.0,-self.w.var(poi).getAsymErrorLo(),self.w.var(poi).getAsymErrorHi())

        print 'Sampling the covariance matrix to propagate error on width'
        histo=ROOT.TH1F('h','h',1000,0.5,1.5)
        for i in range(0,10000):
            # this samples the covariance matrix
            randomizedPars = self.fitResult.randomizePars()
            for j in range(0,randomizedPars.getSize()):
                self.w.var(randomizedPars.at(j).GetName()).setVal(randomizedPars.at(j).getVal())
            histo.Fill(self.w.function('width').getVal())    
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
        graph.SetPointError(len(self.poi),0.0,0.0,histo.GetRMS()/histo.GetMean(),histo.GetRMS()/histo.GetMean())


        for i,label in enumerate(self.poiLabels+['#Gamma_{T}']):
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
        for poi,poiLabel in zip(self.poi,self.poiLabels):
            print poiLabel+':   ('+str(self.w.var(poi).getAsymErrorLo())+','+str(self.w.var(poi).getAsymErrorHi())+')'

        cprime = ROOT.TCanvas('canvas2','')
        cprime.cd()
        histo.GetXaxis().SetTitle("#Gamma_{T}")
        histo.GetYaxis().SetTitle("N toys")
        histo.Draw()
        obj.append(histo)
        print 'Relative error on the total width ',histo.GetRMS()/histo.GetMean()
        print 'Please check the histogram to see that the dist is Gaussian. If not the fit is biased'
        print 'The fit can be biased when floating the width sometimes.'
        

        return c,cprime,obj

        
        


