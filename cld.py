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
# f.addConstraint('Zh','(1+Z)*(1+Z)','Z',1,0.0085*f350)
f.addChannel('Zh', 1., 0.0085*f350, prod='Z')
f.addConstraint('Whbb240','(1+W)*(1+W)*(1+b)*(1+b)/width','W,b,width',1,0.022)
f.addConstraint('Whbb350','(1+W)*(1+W)*(1+b)*(1+b)/width','W,b,width',1,0.006)
# f.addConstraint('Zhbb','(1+Z)*(1+Z)*(1+b)*(1+b)/width','Z,b,width',1,0.0079*f350)
f.addChannel('Zhbb', 1., 0.0079*f350, prod='Z', decay='b')
# f.addConstraint('Zhcc','(1+Z)*(1+Z)*(1+c)*(1+c)/width','Z,c,width',1,0.012*f350)  
f.addChannel('Zhcc', 1., 0.012*f350, prod='Z', decay='c')
f.addChannel('ZhllWWhad', 1., 0.015*f350, prod='Z',
             decay=[('W', 0.65),
                    ('g', 0.25),
                    ('Z', 0.1)]
             )
f.addChannel('ZhllWW1lep', 1., 0.0182*f350, prod='Z',
             decay=[('W', 0.9),
                    ('tau', 0.04),
                    ('Z', 0.06)]
             )
f.addChannel('ZhqqWW2lep', 1., 0.0243*f350, prod='Z',
             decay=[('W', 0.9),
                    ('tau', 0.07),
                    ('b', 0.03)]
             )
f.addChannel('ZhnunuWW', 1., 0.02*f350, prod='Z',
             decay=[('W', 0.6),
                    ('Z', 0.07),
                    ('g', 0.33)]
             )
f.addConstraint('Zhlltautau','(1+Z)*(1+Z)*(1+tau)*(1+tau)/width','Z,tau,width',1,0.0271*f350)
f.addConstraint('Zhqqtautau','(1+Z)*(1+Z)*(1+tau)*(1+tau)/width','Z,tau,width',1,0.0126*f350)
f.addConstraint('ZhZZ','(1+Z)*(1+Z)*(1+Z)*(1+Z)/width','Z,width',1,0.031*f350)
f.addConstraint('Zhgammagamma','(1+Z)*(1+Z)*(1+gamma)*(1+gamma)/width','Z,gamma,width',1,0.03*f350)
f.addConstraint('Zhmumu','(1+Z)*(1+Z)*(1+mu)*(1+mu)/width','Z,mu,width',1,0.13*f350)
f.addUniformConstraint('Zhinv','inv') ####->Means free floating

f.info()

print

################################################
f.fit()
f.createSummary()


