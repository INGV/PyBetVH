#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" 

This file is part of pyBetVH.

"""

import numpy as np

def makeAlpha16(n, p, l, pd):
    """
    Calculate alpha values for nodes 1 to 6
  
    Variables:
  
    n: n. branches
    a: alpha
    p: prior probability
    l: equivalent number of data
    pd: past data

    """
    a = [0]*n
    a0 = l + n - 1
    for i in range(n):
        if p[i] > 0:
            a[i] = (p[i] * a0) + pd[i]
        else:
            #a[i] = 0
            a[i] = 1e-3

    return a
  


def makeAlphaBeta78(p, l):
    """
    Calculte alpha and beta values for nodes 7/8
  
    """
    eps = 0.00001

    ns = len(p)

    #a0 = 0.
    #b0 = 0.
    a0 = 1e-3
    b0 = 1e-3
  
    a = [0]*ns
    b = [0]*ns
  
    a[0] = a0+p[0]*l
    b[0] = b0+l*(1-p[0])
  
  
    #l8_betvh = a0 + b0 - 1 + l;

    for it in range(ns-1):
        if (p[it] > 0):
            #theta8_betvh = a0/(a0+b0+l) + p[it+1]/p[it] * (l / (a0+b0+l))
            ##theta8_betvh = (a0/(a0+b0)) * ((a0+b0)/(a0+b0+l)) + p[it+1]/p[it] * (l / (a0+b0+l))
            #if (theta8_betvh > 0):
                #a[it+1]=theta8_betvh*l8_betvh + theta8_betvh
                #b[it+1]=l8_betvh - a[it+1] + 1
            if (p[it+1] > 0):
                a[it+1] = a[it+1] + l*p[it+1]/p[it]
                b[it+1] = b[it+1] + l*(1-p[it+1]/p[it])
                if b[it+1] <= 0:
                    b[it+1] = 1e-3

            else:
                a[it+1] = 1e-3
                b[it+1] = 1

        else:
            #a[it+1]=0
            a[it+1] = 1e-3
            b[it+1] = 1
  
    return a, b


def theoreticalAverage(a):
    """
    """
    if (np.size(a) <= 2):
        ave = a[0]/np.sum(a)
    else:
        ave = a/np.sum(a)

    return ave

