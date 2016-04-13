"""
A library of methods for use in analyzing data generated from the quickcyst
method of estimating the germination frequency of algal cysts

Written for python 2.7.6

Owen Coyle
ocoyle@uw.edu
12/17/14
"""

from scipy.stats import beta
from numpy.random import gamma
from scipy.misc import comb
from numpy import exp
import numpy as np
from decimal import Decimal, getcontext



def gCDF(g, k, l, n):
    """Calculate the value of the posterior CDF(g|k, l, n) at g.
        g is the germination frequency (can be an array)
        k is the number of wells with no germinated cells
        l is the expected number of cysts per well (can be an array)
        n is the number of wells run
        
        If g is an array, evaluates the CDF at multiple points g
        If l is an array, treats the array as the distribution of l and convolves
    """
    a = k
    b = n - k + 1
    #allocate space for storage
    num = np.zeros(len(g))
    denom = np.zeros(len(g))
    #freeze a beta distribution for fast computation
    rv = beta(a, b)
    for thisL in l:
        num += (1.0 - rv.cdf(exp(-thisL*g)))/thisL
        denom += (1.0 - rv.cdf(exp(-thisL)))/thisL
    return num/denom

def fZero(g, l, n, dg):
    """A helper function for the case when k = 0.
        Computes the integral of (1-e^(-l*g))^n from 0 to g
        g must be a single value.
        l can be an array, in which case it will be treated as the distribution of l
    """
    gs = np.arange(0, g, dg)
    gs = gs.reshape([len(gs), 1])
    l = l.reshape([1, len(l)])
    prod = gs.dot(l)
    p = dg*(1-exp(-prod))**n
    return sum(sum(p))



def gCDFZero(g, l, n, dg):
    """Calculate the value of the posterior CDF(g|k, l, n) at g
        for the special case where k = 1.
        l can be an array, but g must be a single value (in an array of len 1)
    """
    
    return fZero(g, l, n, dg)/fZero(1.0, l, n, dg)

def gCDFfast(g, l, rv):
    """A helper method for the efficient calculation of the invCDF of g.
        Does the same thing as gCDF, but you hand it rv, a frozen beta(a, b)
        where a = k, b = n - k + 1
        g must be a single value (in an array of len 1)
    """
    num = sum((1.0 - rv.cdf(exp(-l*g)))/l)
    denom = sum((1.0 - rv.cdf(exp(-l)))/l)
    return num/denom

def gQDF(p, k, l, n, err, dg):
    """Calculates the inverse of CDF(g|k, l, n)
        p is the value of the CDF for which we would like to find g (can be an array)
        k is the number of wells with no germinated cells
        l is the expected number of cysts per well (can be an array)
        n is the number of wells run
        err is the acceptable error in the calculation of g
        
        If p is an array, solves for the inverse CDF at multiple p
        If l is an array, treats the array as the distirbution of l and convolves
    """
    #allocate storage for the results
    g = np.zeros(len(p))
    #use different methods depending on whether k is non-zero
    if k > 0:
        a = k
        b = n - k + 1
        rv = beta(a, b)
        for i in range(len(p)):
            gLeft = 0.0
            gRight = 1.0
            while (gRight - gLeft) > err:
                gMid = 0.5*(gLeft + gRight)
                pMid = gCDFfast(gMid, l, rv)
                if pMid < p[i]:
                    gLeft = gMid
                else:
                    gRight = gMid
            g[i] = 0.5*(gLeft + gRight)
    elif k == 0:
        for i in range(len(p)):
            gLeft = 0.0
            gRight = 1.0
            while (gRight - gLeft) > err:
                gMid = 0.5*(gLeft + gRight)
                pMid = gCDFZero(gMid, l, n, dg)
                if pMid < p[i]:
                    gLeft = gMid
                else:
                    gRight = gMid
            g[i] = 0.5*(gLeft + gRight)
    return g


def distLambda(m, Nmud, Vmud, Vwell):
    """Choose m values from the distribution of lambda given Nmud, Vmud, Vwell.
        m is the number of values to pull
        Nmud is the number of cysts counted in estimating the concentration of cysts
        Vmud is the volume of mud counted over in estimating the concentrtion of cysts
        Vwell is the volume of mud placed in each well in the wellplate
        """
    return Vwell/Vmud*gamma(Nmud, 1.0, m)

def inputYN(str):
    """Ask the user for a yes or no answer by printing str.
        Keep asking until you get a y or n. Return a boolean True for y, False for n
    """
    yes = ["y", "yes", "Y", "Yes", "YES", ""]
    no = ["n", "no", "N", "No", "NO"]
    userInput = raw_input(str)
    if userInput in yes:
        keepAsking = False
        toReturn = True
    elif userInput in no:
        keepAsking = False
        toReturn = False
    else:
        keepAsking = True
    while keepAsking:
        print "Please enter y or n"
        userInput = raw_input(str)
        if userInput in yes:
            keepAsking = False
            toReturn = True
        elif userInput in no:
            keepAsking = False
            toReturn = False
        else:
            keepAsking = True
    return toReturn

def quickcystmain():
    """The main method, which queries the user and creates a detailed data report
    """
    print "----------------------------------------------------------------------------"
    print "Welcome to quickcyst.py, a program for analyzing algal cyst germination data"
    print "\nPlease start by entering the number of algal cysts you counted in estimating\nthe concentration of cysts in your mud sample"
    Nmud = input("No. of cysts counted: ")
    print "\nIn what volume of mud did you count this many cysts?"
    Vmud = input("Total vol. counted (ml) : ")
    print "\nEstimated cyst concentration = %.02f cysts/ml of mud\n" % (Nmud*1.0/Vmud)
    print "Now, how much mud did you pipette into each well?"
    Vwell = input("Vol. mud per well (ml): ")
    print "\nAvg. No. cysts/well = %.02f\n" % (1.0*Nmud*Vwell/Vmud)
    print "Ok, next how many total wells (No. wells/plate x No. plates) did you run?"
    n = input("No. wells: ")
    print "\nAlmost done. Lastly, in how many wells did you see germinated cells?"
    k = n - input("No. wells w/ germinated cells: ")
    print "\nOk, I have everything I need, would you like to use the default settings?"
    useDefault = inputYN("Use default? (y/n): ")
    if useDefault:
        m = 10000
        pBest = 0.5
        ciSize = 95
        pLower = 0.025
        pUpper = 0.975
        dg = 0.005
        print "Ok, going ahead using the default settings.\n"
        print "I will calculate my best guess for g, the germination probability of a"
        print "single algal cyst, as well the 95% confidence interval of g."
        print "I will use %d draws to estimate the convolution of the distributions" % m
        print "P(g|C) and P(C). There is a certain amount of randomness involved, so"
        print "I recommend you rerun and see if the answers I give are about the same."
        print "If they are not, I recommend you rerun and not use my default settings."
        print "If all of your wells germinated but I give you a weird answer, don't use"
        print "the default settings and try decreasing my num. integr. step size."
        print "Also don't use the default if you want to see a C.I. other than 95%"
    else:
        print "How many draws should I use to estimate the convolution?"
        m = input("No. draws for convolution (recommend > 10,000): ")
        pBest = 0.5
        print "If I need to numerically integrate over g (k = 0) what step size should I use?"
        dg = input("Num. integr. step size (recommend < 0.005): ")
        print "What confidence interval would you like to see?"
        ciSize = input("C.I.(%): ")
        pLower = (100-ciSize)/200.0
        pUpper = 1.0-(100-ciSize)/200.0
        print "Got it, this time I'll use %d draws, and calculate the %.3f%% C.I." % (m, ciSize)
    print "\nCalculating, this may take a bit...especially if all of your wells germinated\n"
    l = distLambda(m, Nmud, Vmud, Vwell)
    err = 0.00001 #no real reason to change this
    p = np.array([pLower, pBest, pUpper])
    g = gQDF(p, k, l, n, err, dg)
    gBest = g[1]
    gLower = g[0]
    gUpper = g[2]
    print ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
    print "Best guess for g: %f" % gBest
    print "%.3f%% C.I.: %f-%f" % (ciSize, gLower, gUpper)
    print ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
    #make some recommendations about how to reduce the uncertainty
    CIwUn = gUpper - gLower
    lBest = Nmud*Vwell/Vmud
    pCI = np.array([pLower, pUpper])
    gnoUn = gQDF(pCI, k, np.array([lBest]), n, err, dg)
    CInoUn = gnoUn[1]-gnoUn[0]
    percDueN = CInoUn/CIwUn
    percDueC = (1.0-percDueN)*100
    print "\nIt looks like %d%% of the uncertainty on g is due to uncertainty on C" % percDueC
    n2 = n*2
    Nmud2 = Nmud*2
    Vmud2 = Vmud*2
    l2 = distLambda(m, Nmud2, Vmud2, Vwell)
    gL2 = gQDF(pCI, k, l2, n, err, dg)
    gN2 = gQDF(pCI, 2*k, l, n2, err, dg)
    gL2N2 = gQDF(pCI, 2*k, l2, n2, err, dg)
    print "\nHere are some recommendations for how you can improve your estimate of g:"
    print "\nIf you doubled the number of cysts you counted for C from %d to %d" % (Nmud, Nmud2)
    print "your C.I. would probably be closer to %f-%f" % (gL2[0], gL2[1])
    print "\nIf you doubled the number of wells you ran from %d to %d" % (n, n2)
    print "your C.I. would probably be closer to %f-%f" % (gN2[0], gN2[1])
    print "\nIf you doubled both your C.I. would probably be closer to %f-%f" % (gL2N2[0], gL2N2[1])
    print "\nJust remember, these are only estimates, you must still run the experiments!"
    print "\nThank you for using quickCyst.py, I hope I was helpful"
    print "----------------------------------------------------------------------------"

if __name__ == "__main__":
    quickcystmain()
