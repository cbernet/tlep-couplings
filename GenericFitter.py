import ROOT
ROOT.gROOT.ProcessLine(".x tdrstyle.C")

class GenericFitter(object):

    def __init__(self):
        self.w = ROOT.RooWorkspace('w','w')
        self.poi = []
        self.poiLabels = []
        self.pdfs=[]
        self.obs=[]

    def addPOI(self,poi,label='',minimum =-0.3 ,maximum = 0.3):
        self.w.factory(poi+'[0,'+str(minimum)+','+str(maximum)+']')
        self.poi.append(poi)
        if label == '':
            self.poiLabels.append(poi)
        else:    
            self.poiLabels.append(label)

    def addConstraint(self,name,equation,dependents,number,error):
        self.w.factory('expr::'+name+"('"+equation+"',"+dependents+')')
        self.w.factory(name+'Obs['+str(number)+']')
        self.obs.append(name+'Obs')
        self.w.factory('RooGaussian::'+name+'Constraint('+name+'Obs,'+name+','+str(error)+')')
        self.pdfs.append(name+'Constraint')

    def addUniformConstraint(self,name,equation,dependents):
        self.w.factory('expr::'+name+"('"+equation+"',"+dependents+')')
        self.w.factory('RooUniform::'+name+'Constraint('+name+')')
        self.pdfs.append(name+'Constraint')

    def fit(self):
        self.w.factory('PROD::model('+','.join(self.pdfs)+')')
        #create dataset
        argset=ROOT.RooArgSet('set')
        for obs in self.obs:
            argset.add(self.w.var(obs))
        data = ROOT.RooDataSet('data','data',argset)
        data.add(argset)
        data.Print()
        self.fitResult=self.w.pdf('model').fitTo(data,ROOT.RooFit.PrintLevel(3),ROOT.RooFit.Optimize(1),ROOT.RooFit.Hesse(1),ROOT.RooFit.Minos(1),ROOT.RooFit.Strategy(2),ROOT.RooFit.Save(1))
#        self.w.pdf('model').fitTo(data)


    def contour(self,var1,var2,n1=1,n2=2,n3=0,n4=0,n5=0):
        self.w.factory('PROD::model('+','.join(self.pdfs)+')')
        #create dataset
        argset=ROOT.RooArgSet('set')
        for obs in self.obs:
            argset.add(self.w.var(obs))
        data = ROOT.RooDataSet('data','data',argset)
        data.add(argset)
        data.Print()

        nll = self.w.pdf('model').createNLL(data,ROOT.RooFit.Optimize(1),ROOT.RooFit.Offset(0))
        minuit =ROOT.RooMinuit(nll)
        return minuit.contour(self.w.var(var1),self.w.var(var2),n1,n2,n3,n4,n5)
#        self.w.pdf('model').fitTo(data)



        

    
    def w(self):
        return self.w
        




        
        
