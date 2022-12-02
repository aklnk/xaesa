# -*- coding: utf-8 -*-
"""
Created on Fri Nov  4 07:20:37 2022

@author: A.Kuzmin
"""

import numpy as np
import matplotlib.pyplot as plt

from scipy import interpolate
from scipy.interpolate import make_lsq_spline
from scipy.interpolate import InterpolatedUnivariateSpline


#, BSpline

#from scipy.interpolate import make_interp_spline


##########################################################
#
#  Zero line for atomic absorption by LSQ spline
#  
##########################################################

#
#  Absorption coefficient mu(E) from E0 to Emax
# 
def LSQ_ZeroLine(x, y, k):

    xk = np.zeros(shape=(x.size))

    ysk = np.zeros(shape=(xk.size))

# B-spline degree. Default is cubic, k=3.
#    k = 3

# Calculate number of intervals
    nk = round(x[-1]/500.0)+1
      
    if nk <= 3: 
        nk = 3
    print(nk)
        
# Calculate k axis    
    xk = (2.0/7.62*x)**0.5
    print(xk)

# Calculate knots 
    tk = np.empty(shape=(nk-1))
    tk.fill(0)
    
    for i in range(0, nk-1):
        tk[i] = np.sqrt((2.0/7.62*(i+1)*x[-1]/(nk))) #**0.5
    print(tk)
    
    # p = np.linspace(1, nk, nk)
    # tk = (2.0/7.62*p*x[-1]/(nk))**0.5
    # print(tk)
    
    tt = np.r_[(xk[0],)*(k+1),  tk,  (xk[-1],)*(k+1)]
    # print(tt)


# Compute the (coefficients of) an LSQ B-spline.
    splk3 = make_lsq_spline(xk, y*xk*xk*xk, tt, k)

# Compute zero line
    ysk = splk3(xk)/xk/xk/xk
    
    i = 0
    while xk[i] < tk[1]:
        i = i + 1

    n = i
    
    
    xk_ = np.empty(shape=(xk.size-n))
    xk_.fill(0) 
    
    ysk_ = np.empty(shape=(xk.size-n))
    ysk_.fill(0)       
    
    for j in range(n,xk.size):
        ysk_[j-n] = ysk[j]
        xk_[j-n] = xk[j]

    return ysk_

def extend_zeroline(x, y, x_new):    

# Spline order: 1 linear, 2 quadratic, 3 cubic ... 
    order = 3

# Do extrapolation from E0 to Emin
    s = InterpolatedUnivariateSpline(x, y, k=order)

    return s(x_new)

##########################################################


# if __name__ == '__main__':


#     # Main code
     
#     nf = 999
    
#     x = np.empty(shape=(nf))
#     x.fill(0)    
    
#     xk = np.empty(shape=(nf))
#     xk.fill(0) 
    
    
    
#     y = np.empty(shape=(nf))
#     y.fill(0)
    
#     #ysk2 = np.empty(shape=(nf))
#     #ysk2.fill(0)
    
#     #yy = np.empty(shape=(nf))
#     #yy.fill(0)
    
#     #
#     # Create "Experimental data"
#     #   
#     for i in range (0,nf):
#         xk[i] = 30.0*(i+1)/(nf+1)    # k space
#         x[i] = xk[i]*xk[i]*7.62/2.0  # E space
#         y[i] = 6.0*np.sin(2.0*xk[i]*2.0)*np.exp(-2.0*0.01*0.01*xk[i]*xk[i])/(xk[i]*2.0*2.0) 
#     #    y[i] = y[i] + 2*np.sin(0.15*xk[i])       
#     #    y[i] = y[i] + 2*np.sin(0.1*xk[i])         
#         y[i] = y[i] + 2*np.sin(0.01*xk[i])   

#     #
#     # Calculate zero line for absorption coefficient mu(E)=y(x) from E0 to Emax
#     #
#     yy = LSQ_ZeroLine(x, y, 3)

#     plt.plot(x, y, 'r-', lw=1, ms=5)
    
#     #plt.plot(x, ysk2,'b-', lw=3, label='LSQ spline k^2')
    
#     plt.plot(x, yy,'k-', lw=3, label='LSQ spline k^2')

#     plt.legend(loc='best')
    
#     plt.show()

