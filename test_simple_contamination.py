from CouplingsFitterTest import CouplingsFitterTest

f= CouplingsFitterTest()

###Here add the Parameters of interest with a reasonable range
##############################################################
f.addPOI('b','b',-0.1,0.1)

###Here add the constraints 'name','formula','dependents',mean value ,error 
################################################
# mike added a factor of 0.95 which represents
# the improvement on the 240 measurements based on ZH
y = 2
f.addConstraint('yb','{}*(1+b)'.format(y),
                'b', y, y*.01)  

f.info()

print

################################################
f.fit()
f.createSummary()



