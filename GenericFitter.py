import ROOT

class GenericFitter(object):
    
    def __init__(self):
        self.w = ROOT.RooWorkspace('w','w')
        self.poi = []
        self.poiLabels = []
        self.pdfs = []
        self.obs = []
        self.cyields = []
        self.cwidths = []

    def addPOI(self,poi,label='',minimum =-0.3 ,maximum = 0.3):
        '''Add a parameter of interest,

        Example:
        
        addPOI('Z','Z',-0.05,0.05)
        ->
        Z[0,-0.05,0.05]
        # adds variable gZ with value 0,
        # allow it to scale between 0.95 and 1.05      
        '''
        
        fcmd = poi+'[0,'+str(minimum)+','+str(maximum)+']'
        self.w.factory(fcmd)
        print 'Add POI', '*' * 50
        print fcmd
        self.poi.append(poi)
        if label == '':
            self.poiLabels.append(poi)
        else:    
            self.poiLabels.append(label)

    def addConstraint(self,name,equation,dependents,number,error):
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
        fcmd0 = 'expr::'+name+"('"+equation+"',"+dependents+')'
        self.w.factory(fcmd0)
        fcmd1 = name+'Obs['+str(number)+']'
        self.w.factory(fcmd1)
        self.obs.append(name+'Obs')
        fcmd2 = 'RooGaussian::'+name+'Constraint('+name+'Obs,'+name+','+str(error)+')'
        self.w.factory(fcmd2)
        print 'Add constraint', '*' * 50
        print fcmd0
        print fcmd1
        print fcmd2
        self.pdfs.append(name+'Constraint')
        self.cyields.append(name)
        self.cwidths.append(error)

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
        self.data = ROOT.RooDataSet('data','data',argset)
        self.data.add(argset)
        self.data.Print()
        self.fitResult=self.w.pdf('model').fitTo(self.data,
                                                 ROOT.RooFit.PrintLevel(3),
                                                 ROOT.RooFit.Optimize(1),
                                                 ROOT.RooFit.Hesse(1),
                                                 ROOT.RooFit.Minos(1),
                                                 ROOT.RooFit.Strategy(2),
                                                 ROOT.RooFit.Save(1))
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
        




        
        
