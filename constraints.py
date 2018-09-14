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
