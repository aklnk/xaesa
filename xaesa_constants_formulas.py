# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 14:49:33 2018

@author: AKalinko
"""

from numpy import asarray, exp, zeros

#constants
me = 9.10938215* 10**-31 # kg
h = 6.626070040* 10**-34 #Js = m^2 kg / s
hev =  4.135667516 * 10**-15 #eV s
c = 299792458 #m/s
hbar = 1.054571800*10**-34

def victoreen(x, a, b):
    return (a / x**4 / (hev * c)**3 + b / x**3 / (hev * c)**3 )

def windowGauss10(k, kmin, kmax):
    window = zeros(len(k))
    ka = (kmin+kmax)/2.0
    kw = (kmax-kmin)*(kmax-kmin)/9.210340372
    knp = asarray(k, float)
    wp = -(knp-ka)*(knp-ka)/kw
    exp(wp, window)
        
    return window