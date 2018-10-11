from CouplingsFitter2 import CouplingsFitter2

import math

f= CouplingsFitter2()

lumi = 5.
base_lumi = 5.
lumi_factor = 1 / math.sqrt(lumi/base_lumi)

###Here add the Parameters of interest with a reasonable range
##############################################################
f.addPOI('Z','Z',-0.5,0.5)
f.addPOI('W','W',-0.5,0.5)
f.addPOI('b','b',-0.5,0.5)
f.addPOI('c','c',-0.5,0.5)
f.addPOI('g','g',-0.5,0.5)
f.addPOI('tau','#tau',-0.5,0.5)
# f.addPOI('t','t',-1,1)
f.addPOI('mu','#mu',-1., 1.)
f.addPOI('gamma','#gamma',-1., 1.)
f.addPOI('inv','inv', -1, 1.)
f.createWidthDeviation()    

###Here add the constraints 'name','formula','dependents',mean value ,error 
################################################
# mike added a factor of 0.95 which represents
# the improvement on the 240 measurements based on ZH
f.addConstraint('Zh','(1+Z)*(1+Z)','Z',1,0.005*lumi_factor)  
# f.addConstraint('Wh','(1+W)*(1+W)','W',1,0.004)
# f.addConstraint('Wh250','(1+W)*(1+W)','W',1,0.02)
f.addConstraint('Whbb240','(1+W)*(1+W)*(1+b)*(1+b)/width','W,b,width',1,0.031*lumi_factor)
# f.addConstraint('Whbb350','(1+W)*(1+W)*(1+b)*(1+b)/width','W,b,width',1,0.006*lumi_factor)
f.addConstraint('Zhbb','(1+Z)*(1+Z)*(1+b)*(1+b)/width','Z,b,width',1,0.003*lumi_factor)
f.addConstraint('Zhcc','(1+Z)*(1+Z)*(1+c)*(1+c)/width','Z,c,width',1,0.022*lumi_factor)  
f.addConstraint('Zhgg','(1+Z)*(1+Z)*(1+g)*(1+g)/width','Z,g,width',1,0.019*lumi_factor)
f.addConstraint('ZhWW','(1+Z)*(1+Z)*(1+W)*(1+W)/width','Z,W,width',1,0.012*lumi_factor)
f.addConstraint('Zhtautau','(1+Z)*(1+Z)*(1+tau)*(1+tau)/width','Z,tau,width',1,0.009*lumi_factor)
f.addConstraint('ZhZZ','(1+Z)*(1+Z)*(1+Z)*(1+Z)/width','Z,width',1,0.044*lumi_factor)
f.addConstraint('Zhgammagamma','(1+Z)*(1+Z)*(1+gamma)*(1+gamma)/width','Z,gamma,width',1,0.09*lumi_factor)
f.addConstraint('Zhmumu','(1+Z)*(1+Z)*(1+mu)*(1+mu)/width','Z,mu,width',1,0.19*lumi_factor)
f.addUniformConstraint('Zhinv','inv') ####->Means free floating

f.info()

print

################################################
f.fit()
f.createSummary()


