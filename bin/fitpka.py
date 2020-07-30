#!/usr/bin/env python

import argparse
import math
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

PH2KCAL = 1.364


def sigmoid(x, x0, k):
    # This function changes from 0 to -1, k is the slope and x0 is the midpoint
    e = np.exp(k*(x-x0))
    y = e/(1.0+e)
    return y


def plot_res(resname, charge, t_type, t_points):

    # make parameters to convert y to sigmoid function range (0, -1)
    # charge[x] = a(exp(k(x-x0))/(1+exp(k(x-x0)) + b
    # for charge change range 1, a = 1 or -1, k is slope, x0 is midpoint, b is offset for oxidation like from +2 to +3
    deltaY = charge[-1] - charge[0]  # This determines scale up value a (always positive), 1, 2, 3 ...
    a = math.ceil(abs(deltaY))
    lowest = min(charge)  # This determines the lower end of the y value
    b = math.floor(lowest)
    # assume y=exp(k(x-x0))/(1+exp(k(x-x0)
    x = t_points
    y = [(c - b) / a for c in charge]
    xdata = np.array(x)
    ydata = np.array(y)

    # fit the function to get x0 and k
    msg = ""
    try:
        (popt, pcov) = curve_fit(sigmoid, xdata, ydata, bounds=([x[0], -4], [x[-1], 4]))
        # (popt, pcov) = curve_fit(sigmoid, xdata, ydata)
    except RuntimeError:
        msg = "Titration out of range"
    except ValueError:
        msg = "Input value not valid"

    if popt[0] < x[0] + 0.001 or popt[0] > x[-1] - 0.001:
        msg = "Titration out of range"

    # convert back to the titration curve, midpoint, curve, and error
    if msg:
        print(msg)
    else:
        chi_squared = np.sum([(sigmoid(xdata, *popt) - ydata) ** 2])
        midpoint = popt[0]
        nslope = 0
        if t_type.upper() == "PH":
            nslope = 0.4342 * popt[1]
        elif t_type.upper() == "EM":
            nslope = 0.4342 * popt[1] * 58.0
        elif t_type.upper() == "CH" or t_type.upper() == "EXTRA":
            nslope = 0.4342 * popt[1] * PH2KCAL
        else:
            print("Why am I here?")

        # plot the titration graph
        plt.plot(xdata, charge, 'kx')
        Npoints = 100
        xfit = [xdata[0] + i * (xdata[-1] - xdata[0]) / Npoints for i in range(Npoints + 1)]
        yfit = [a * sigmoid(x, popt[0], popt[1]) + b for x in xfit]

        plt.plot(xfit, yfit,
                 label='%s: pKa=%.2f,  n=%.2f,  chi^2*1000=%.2f' % (resname, midpoint, abs(nslope), chi_squared * 1000.0))

    return

def fitpka(resnames):
    # get titration type and titration points, this is x axis
    lines = open("sum_crg.out").readlines()
    headline = lines.pop(0)
    fields = headline.split()
    t_type = fields[0].strip()
    t_points = [float(x) for x in fields[1:]]

    # get residue charge, this is y axis
    plt.figure(figsize=(10, 5))
    for resname in resnames:
        found = False
        charge = []
        for line in lines:
            fields = line.split()
            if fields[0] == resname:
                found = True
                charge = [float(x) for x in fields[1:]]
                break

        if found:
            plot_res(resname, charge, t_type, t_points)
            plt.xlabel("%s" % t_type)
            plt.ylabel('Charge')
            plt.legend(loc='upper right')

    plt.show()


    return

if __name__ == "__main__":

    # Get the command arguments
    helpmsg = "Fit a titration of charged residues"
    parser = argparse.ArgumentParser(description=helpmsg)
    parser.add_argument("residue", metavar="RES", nargs="+", help="Charged residue names to plot, as in sum_crg.out")
    args = parser.parse_args()

    resnames = args.residue
    fitpka(resnames)