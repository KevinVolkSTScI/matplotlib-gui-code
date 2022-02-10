#! /usr/bin/env python
#
import sys
import numpy
import mpfit
import mpfitexpr

def myfunction(p, fjac=None, x=None, y=None, err=None):
    model = p[2]/numpy.exp((x-p[0])*(x-p[0])/(2.*p[1]*p[1]))
    model = model + p[5]/numpy.exp((x-p[3])*(x-p[3])/(2.*p[4]*p[4]))
    status = 0
    deviations = (y-model)/err
    return [status, deviations]

a0 = 30.0
a1 = 50.0
s1 = 5.0
s2 = 6.0
p1 = 2.0
p2 = 3.0
xvalues = numpy.arange(100001)/1000.
yvalues = p1/numpy.exp((xvalues-a0)*(xvalues-a0)/(2.*s1*s1))
yvalues = yvalues + p2/numpy.exp((xvalues-a1)*(xvalues-a1)/(2.*s2*s2))
noise = 0.03*(numpy.random.random(100001)-0.5)
yvalues = yvalues+noise
errors = yvalues*0.+0.01
startvalues = numpy.asarray([27., 3.0, 2.5, 55., 3.0, 3.5])
fa = {'x': xvalues, 'y': yvalues, 'err': errors}
results = mpfit.mpfit(myfunction, startvalues,functkw=fa)
print(results.params)
ymodel = results.params[2]/numpy.exp((xvalues-results.params[0])*(xvalues-results.params[0])/(2.*results.params[1]*results.params[1]))
ymodel = ymodel + results.params[3]/numpy.exp((xvalues-results.params[4])*(xvalues-results.params[4])/(2.*results.params[5]*results.params[5]))
outfile = open('test_fit_values.txt', 'w')
for loop in range(len(xvalues)):
    print('%7.3f %13.6e %13.6e %13.6e' % (xvalues[loop], yvalues[loop], errors[loop], ymodel[loop]), file=outfile)
outfile.close()
