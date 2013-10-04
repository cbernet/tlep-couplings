from GenericFitter import *
from math import sqrt
class CouplingsFitter(GenericFitter):
    def __init__(self):
        self.BR=dict()
#        self.BR['b'] = 0.561
#        self.BR['tau'] = 0.061
#        self.BR['mu'] = 2.14e-4
#        self.BR['c'] = 2.83e-2
#        self.BR['g'] = 8.48e-2
#        self.BR['gamma'] = 3.9e-3
#        self.BR['W'] = 0.239
#        self.BR['Z'] = 0.0289
#        self.BR['t'] = 0.0

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



    def getAnalyticalWidthError(self):
        sumBR=0.0
        widthError=0.0
        
        for coupling,br in self.BR.iteritems():
            sumBR=sumBR+br
        for coupling,br in self.BR.iteritems():
            if coupling in self.poi:
                widthError+=4*((br/sumBR)*(1+self.w.var(coupling).getVal())*self.w.var(coupling).getError())*((br/sumBR)*(1+self.w.var(coupling).getVal())*self.w.var(coupling).getError())
        if 'inv' in self.poi:
            widthError=widthError/((1+self.w.var('inv').getVal()*self.w.var('inv').getVal())*(1+self.w.var('inv').getVal()*self.w.var('inv').getVal()))

            
            #loop to get the second term
            secondterm=0.0
            for coupling,br in self.BR.iteritems():
                if coupling in self.poi:
                    secondterm+=(br/sumBR)*(1+self.w.var(coupling).getVal())*(1+self.w.var(coupling).getVal())
            secondterm=secondterm*2*self.w.var('inv').getVal()*self.w.var('inv').getError()/((1-self.w.var('inv').getVal()*self.w.var('inv').getVal())*(1-self.w.var('inv').getVal()*self.w.var('inv').getVal()))
            widthError = widthError+secondterm*secondterm

        widthError=sqrt(widthError)    
        print 'Total width (analytical)=+-',widthError    

    def getSemiAnalyticalWidthError(self):
        sumBR=0.0
        widthErrorUp=0.0
        widthErrorDown=0.0

        
        for poi in self.poi:
            preVal = self.w.var(poi).getVal()
            preWidth=self.w.function('width').getVal()
            self.w.var(poi).setVal(self.w.var(poi).getAsymErrorHi())
            widthUp=self.w.function('width').getVal()
            self.w.var(poi).setVal(self.w.var(poi).getAsymErrorLo())
            widthDown=self.w.function('width').getVal()



            if widthUp>widthDown:
                widthErrorUp+=(widthUp-preWidth)*(widthUp-preWidth)
                widthErrorDown+=(preWidth-widthDown)*(preWidth-widthDown)
                if poi=='inv':
                    print 'Gamma_{inv}:(-'+str((preWidth-widthDown))+','+str((widthUp-preWidth))+')'

            else:
                widthErrorDown+=(preWidth-widthUp)*(preWidth-widthUp)
                widthErrorUp+=(widthDown-preWidth)*(widthDown-preWidth)
                if poi=='inv':
                    print 'Gamma_{inv}:(-'+str((preWidth-widthUp))+','+str((widthDown-preWidth))+')'

            self.w.var(poi).setVal(preVal)


        print 'Total width : (-'+str(sqrt(widthErrorDown))+','+str(sqrt(widthErrorUp))+')'    
            


        
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
#            self.w.var('inv').setVal(0.0)

        self.createWidthDeviation()    



    def createSummary(self):
        #sample the covariance matrix for the width
        ROOT.gStyle.SetOptTitle(0)
        graph = ROOT.TGraphAsymmErrors(len(self.poi)+1)

        for i,poi in enumerate(self.poi):
            graph.SetPoint(i,i+0.5,self.w.var(poi).getVal())
            if poi !='inv':
                graph.SetPointError(i,0.0,0.0,-self.w.var(poi).getAsymErrorLo(),self.w.var(poi).getAsymErrorHi())


        print 'Sampling the covariance matrix to propagate error on width'
        histo=ROOT.TH1F('h','h',1000,0,2)
        for i in range(0,10000):
            randomizedPars = self.fitResult.randomizePars()
            for j in range(0,randomizedPars.getSize()):
                self.w.var(randomizedPars.at(j).GetName()).setVal(randomizedPars.at(j).getVal())
            histo.Fill(self.w.function('width').getVal())    

        self.histo=histo    
#        f=ROOT.TF1('g','gaus',0,2)
#        widthError = histo.GetRMS()/histo.GetMean()    
#        graph.SetPoint(len(self.poi),len(self.poi)+0.5,0.0)
#        graph.SetPointError(len(self.poi),0.0,0.0,widthError,widthError)


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

        for i,label in enumerate(self.poiLabels+['#Gamma_{T}']):
            obj.append(ROOT.TLatex(i+0.5,0.95*graph.GetYaxis().GetXmin(),label))
            obj[-1].Draw()


        print 'POSTFIT SUMMARY------------------------------------'    
        for poi in self.poi:
            print 'PostFit:'+poi + ' ='+str(self.w.var(poi).getVal())


        print 'MINOS SUMMARY------------------------------------'    
        for poi,poiLabel in zip(self.poi,self.poiLabels):
            print poiLabel+':   ('+str(self.w.var(poi).getAsymErrorLo())+','+str(self.w.var(poi).getAsymErrorHi())+')'


#        print 'PARABOLIC SUMMARY------------------------------------'    
#        for poi,poiLabel in zip(self.poi,self.poiLabels):
#            print poiLabel+':   '+str(self.w.var(poi).getError())



#        print 'Total width (toys)=+-',widthError    
#        self.getAnalyticalWidthError()
        self.getSemiAnalyticalWidthError()
        return c,obj

        
        


