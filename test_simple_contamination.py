from CouplingsFitterTest import CouplingsFitterTest

f= CouplingsFitterTest()

###Here add the Parameters of interest with a reasonable range
##############################################################
f.addPOI('b','b',-0.1,0.1)
f.addPOI('a','a',-0.1,0.1)

###Here add the constraints 'name','formula','dependents',mean value ,error 
################################################
# mike added a factor of 0.95 which represents
# the improvement on the 240 measurements based on ZH
y2 = 2
y1 = 4
b1 = 0.6 * y1
a1 = 0.4 * y1

f.addConstraint('y2','{}*(1+b)'.format(y2),
                'b', y2, y2*.01)
# f.addConstraint('y1','({b1}*(1+b)/({b1}*(1+b)+{a1}*(1+a))*(1+b) + {a1}*(1+a)/({b1}*(1+b)+{a1}*(1+a))*(1+a))*{y1}'.format(b1=b1, a1=a1, y1=y1),
#                 'b,a', y1, y1*.01)

f.addConstraint('y1','{b1}*(1+b) + {a1}*(1+a)'.format(b1=b1, a1=a1),
                 'b,a', y1, y1*.01)

#TODO
# build more complex simple model
# think about interface. formula will get complicated with width, etc. 

f.info()

print

################################################
f.fit()
f.createSummary()



