'''
Created on Jan 4, 2012

@author: tel
'''

import os, sys, re
import numpy
from pylab import *
from scipy import *
from scipy import optimize
import matplotlib.pyplot as plt


class Parameter:
    def __init__(self, value):
        self.value = value

    def set(self, value):
        self.value = value

    def __call__(self):
        return self.value

def fit(function, parameters, y, x = None):
    def f(params):
        i = 0
        for p in parameters:
            p.set(params[i])
            i += 1
        return y - function(x)

    if x is None: x = arange(y.shape[0])
    p = [param() for param in parameters]
    optimize.leastsq(f, p, maxfev = 10000)


def fit_onepk(xdata, ydata):
    n0 = Parameter(1)
    pk0 = Parameter(7)
    isacid = 1
    for i, entry in enumerate(ydata):
    	if entry < 0:
    		ydata[i] = -entry
    		isacid = -1
    def fitfunc(x): return 10**(-isacid*n0()*(pk0()-x))/(1 + 10**(-isacid*n0()*(pk0()-x)))
    
    fit(fitfunc, [n0, pk0], ydata, xdata)
    return (pk0.value, n0.value*isacid)
            
def Chi2onepk(onepk, ph_array, mcce_array):
    su = 0
    pk0 = onepk[0]
    n0 = onepk[1]
    for ph, mcce_occ in zip(ph_array, mcce_array):
        fit_occ = 10**(n0*(ph-pk0))/(1 + 10**(n0*(ph-pk0)))
        diff = fit_occ - mcce_occ
        su += (diff**2)
    return su

def string_to_list(string):
	list = []
	temp = string.split()
	for entry in temp:
		if entry != '':
			list.append(float(entry))
	return list
	
def x_maker(list):
	return range(len(list))

def do_onepk_fit(x, y):
	onepk = fit_onepk(x, y)
	onepk_chi2 = Chi2onepk(onepk, x, y)
	prefixes = ['pK', 'n', '1000*chi**2']
	output = [onepk[0], onepk[1], 1000*onepk_chi2]
	for i, prefix in enumerate(prefixes):	
		print '%s: %s' % (prefix, output[i])
	return onepk

# def fit_2pk(xdata, ydata, pk1g = 4.75, n1g = 2, pk2g = 4.75, n2g =2):
#     n1 = Parameter(n1g)
#     pk1 = Parameter(pk1g)
#     n2 = Parameter(n2g)
#     pk2 = Parameter(pk2g)
#     a = Parameter(1)
#     def fitfunc(x): return a()*(10**(n1()*(x-pk1()))/(1 + 10**(n1()*(x-pk1())))) + (1 - a())*(10**(n2()*(x-pk2()))/(1 + 10**(n2()*(x-pk2()))))
#     
#     fit(fitfunc, [n1, pk1, n2, pk2, a], ydata, xdata)
#     return (pk1.value, n1.value, pk2.value, n2.value, a.value)
# 
# def Chi2twopk(twopk, ph_array, mcce_array):
#     su = 0
#     pk1 = twopk[0]
#     n1 = twopk[1]
#     pk2 = twopk[2]
#     n2 = twopk[3]
#     a= twopk[4]
#     for ph, mcce_occ in zip(ph_array, mcce_array):
#         fit_occ = a*(10**(n1*(ph-pk1))/(1 + 10**(n1*(ph-pk1)))) + (1 - a)*(10**(n2*(ph-pk2))/(1 + 10**(n2*(ph-pk2))))
#         diff = fit_occ - mcce_occ
#         su += (diff**2)
#     return su

def sum_crg_plot(x, y, pk1):
	pk0 = pk1[0]
	n0 = pk1[1]
	x_smooth = numpy.linspace(x[0], x[-1], 100)
	y_smooth = []
	for ph in x_smooth:
		y_smooth.append(10**(n0*(ph-pk0))/(1 + 10**(n0*(ph-pk0))))
	plt.scatter(x, y)
	plt.plot(x_smooth, y_smooth)
	plt.xlabel('pH')
	plt.ylabel('occupancy of ionized states')
	plt.show()
	plt.savefig("sum_crg_plot.png")

if __name__ == '__main__':
	if len(sys.argv) == 2:
		y = numpy.array(string_to_list(sys.argv[1]))
		x = numpy.array(x_maker(string_to_list(sys.argv[1])))
		do_onepk_fit(x, y)
	elif len(sys.argv) == 3:
		if sys.argv[2] == 'True' or sys.argv[2] == 'true':
			y = numpy.array(string_to_list(sys.argv[1]))
			x = numpy.array(x_maker(string_to_list(sys.argv[1])))
			onepk = do_onepk_fit(x, y)
			sum_crg_plot(x, y, onepk)
		elif sys.argv[2] == 'False' or sys.argv[2] == 'false':
			y = numpy.array(string_to_list(sys.argv[1]))
			x = numpy.array(x_maker(string_to_list(sys.argv[1])))
			onepk = do_onepk_fit(x, y)
		else:
			y = numpy.array(string_to_list(sys.argv[1]))
			x = numpy.array(string_to_list(sys.argv[2]))
			onepk = do_onepk_fit(x, y)
	elif len(sys.argv) == 4:
		if sys.argv[3] == 'True' or sys.argv[3] == 'true':
			y = numpy.array(string_to_list(sys.argv[1]))
			x = numpy.array(string_to_list(sys.argv[2]))
			onepk = do_onepk_fit(x, y)
			sum_crg_plot(x, y, onepk)
		elif sys.argv[3] == 'False' or sys.argv[3] == 'false':
			y = numpy.array(string_to_list(sys.argv[1]))
			x = numpy.array(string_to_list(sys.argv[2]))
			onepk = do_onepk_fit(x, y)		
	else:
		print "error: check arguments for format"
		print '''Example syntax:
python fitpk.py 'list of occupancies' 'list of pH values' true

as in:			
user@sibyl:$ python fitpk.py '-0.00 -0.00 -0.00 -0.01 -0.08 -0.43 -0.81 -0.96 -0.99 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00' '0     1     2     3     4     5     6     7     8     9    10    11    12    13    14'

optionally, you can also omit the list of pH values and the program will guess them:
user@sibyl:$ python fitpk.py '-0.00 -0.00 -0.00 -0.01 -0.08 -0.43 -0.81 -0.96 -0.99 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00'

as well, if you add true as the last argument, the program will plot your data and save it as a .png:
user@sibyl:$ python fitpk.py '-0.00 -0.00 -0.00 -0.01 -0.08 -0.43 -0.81 -0.96 -0.99 -1.00 -1.00 -1.00 -1.00 -1.00 -1.00' '0     1     2     3     4     5     6     7     8     9    10    11    12    13    14' true'''