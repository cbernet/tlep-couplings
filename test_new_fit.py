from CouplingsFitter2 import CouplingsFitter2

f= CouplingsFitter2()

###Here add the Parameters of interest with a reasonable range
##############################################################
f.addPOI('Z','Z',-0.05,0.05)
f.addPOI('W','W',-0.05,0.05)
f.addPOI('b','b',-0.1,0.1)
f.addPOI('c','c',-0.5,0.5)
f.addPOI('g','g',-0.5,0.5)
f.addPOI('tau','#tau',-0.2,0.2)
# f.addPOI('t','t',-1,1)
f.addPOI('mu','#mu',-0.5,0.5)
f.addPOI('gamma','#gamma',-1,1)
f.addPOI('inv','inv', 0, 0.01)
f.createWidthDeviation()    
##
##f.addPOI('Z','Z',-2, 2)
##f.addPOI('W','W',-2, 2)
##f.addPOI('b','b',-2, 2)
##f.addPOI('c','c',-2, 2)
##f.addPOI('g','g',-2, 2)
##f.addPOI('tau','#tau',-2, 2)
##f.addPOI('mu','#mu',-2, 2)s
##f.addPOI('gamma','#gamma',-2, 2)
### f.addPOI('inv','inv',0,0.1)
##f.createWidthDeviation()    

###Here add the constraints 'name','formula','dependents',mean value ,error 
################################################
# mike added a factor of 0.95 which represents
# the improvement on the 240 measurements based on ZH
f350 = 0.95
f.addChannel('Zh', 1., 0.004*f350, prod='Z')
f.addChannel('Whbb240', 1., 0.022, prod='W', decay=[('b', 1.)])
f.addChannel('Whbb350', 1., 0.006, prod='W', decay='b')
f.addChannel('Zhbb', 1., 0.002*f350, prod='Z', decay='b')
f.addChannel('Zhcc', 1., 0.012*f350, prod='Z', decay='c')
f.addChannel('Zhgg', 1., 0.014*f350, prod='Z', decay='g')
f.addChannel('ZhWW', 1., 0.009*f350, prod='Z', decay='W')
f.addChannel('Zhtautau', 1., 0.007*f350, prod='Z', decay='tau')
f.addChannel('ZhZZ', 1., 0.031*f350, prod='Z', decay='Z')
f.addChannel('Zhgammagamma', 1., 0.03*f350, prod='Z', decay='gamma')
f.addChannel('Zhmumu', 1., 0.13*f350, prod='Z', decay='mu')
f.addUniformConstraint('Zhinv','inv') ####->Means free floating

f.info()

print

################################################
f.fit()
f.createSummary()


