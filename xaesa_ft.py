#import numpy as np
from numpy import sqrt, zeros, pi, arctan, asarray, arange, delete, sum, multiply, sin, cos

def FT(k, exafs, rmin, rmax, dr):
    
    con = sqrt(2 / pi)
    
    nn = len(k)
    
    rx = zeros(int((rmax - rmin) / dr), float)
    exafs_re = zeros(nn, float)
#    exafs_im = zeros(nn, float)
#    transform_re = zeros(nn, float)
#    transform_im = zeros(nn, float)
#    fourier_re = zeros(int((rmax - rmin) / dr), float)
#    fourier_im = zeros(int((rmax - rmin) / dr), float)
#    sn = zeros(nn, float)
#    cs = zeros(nn, float)
    
    exafs_re = asarray(exafs, float)

    rx = arange(rmin, rmax, dr)
    r = rx*2
    
#    dx1 = delete(k,0)
#    dx2 = delete(k, len(k)-1)
#    dx = dx1 - dx2
    dx = k[1:]-k[:-1]
    
    v1 = multiply.outer(r, asarray(k,float))
    
#    eix = exp(1j*v1)
#    cos1, sin1 = eix.real, eix.imag
    
    sin1 = sin(v1)    
    cos1 = cos(v1) 
    sin1 = sin1*(-1)
    
#    transform_re1 = exafs_re * cos1 - exafs_im * sin1
#    transform_im1 = exafs_re * sin1 + exafs_im * cos1
    transform_re1 = exafs_re * cos1
    transform_im1 = exafs_re * sin1
    

    tre11 = delete(transform_re1, 0, 1)
    tre22 = delete(transform_re1, len(transform_re1[0])-1, 1)
    
    tim11 = delete(transform_im1, 0, 1)
    tim22 = delete(transform_im1, len(transform_im1[0])-1, 1)
    
    tre111 = (tre11+tre22)/2
    tim111 = (tim11+tim22)/2

    r11 = sum(dx*tre111, axis=1)
    r22 = sum(dx*tim111, axis=1)
    
    fourier_re = r11*con
    fourier_im = r22*con
    
    return rx, fourier_re, fourier_im
    
    
def BFT(r, fre, fim, kmin, kmax, dk):
    
    print(kmin, kmax, dk)
    
    con = sqrt(2 / pi)
    
    n  = int((kmax-kmin)/dk + 1.0)
    dk = (kmax-kmin)/(n-1.0)
    bftk = []

    nn = len(r)

    bftr = zeros(n, float)
    bfti = zeros(n, float)
    #transform_re = zeros(nn, float)
    #transform_im = zeros(nn, float)
    #sn = zeros(nn, float)
    #cs = zeros(nn, float)

    bftk = arange(kmin, kmax+dk, dk)
    
    dx1 = delete(r,0)
    dx2 = delete(r, len(r)-1)
    dx = dx1 - dx2
    r5 = bftk * 2
    
    v1 = multiply.outer(r5, asarray(r,float))
    sin1 = sin(v1)
    cos1 = cos(v1) 
    transform_re1 = fre * cos1 - fim * sin1
    transform_im1 = fre * sin1 + fim * cos1
    

    tre11 = delete(transform_re1, 0, 1)
    tre22 = delete(transform_re1, len(transform_re1[0])-1, 1)
    
    tim11 = delete(transform_im1, 0, 1)
    tim22 = delete(transform_im1, len(transform_im1[0])-1, 1)
    
    tre111 = (tre11+tre22)/2
    tim111 = (tim11+tim22)/2

    r11 = sum(dx*tre111, axis=1)
    r22 = sum(dx*tim111, axis=1)
    
    bftr = r11*con
    bfti = r22*con
  
    return bftk, bftr, bfti
    
def BFTWindow(r, rmin, rmax, a):
    
    xmin = rmin
    xmax = rmax
    a1 = a
    #wind = np.zeros(len(r))
    wind = [0] * len(r)
    for i in range(0, len(r)):
        if r[i] < rmin:
            wind[i] = 0.0
        if r[i] > rmax:
            wind[i] = 0.0
        if r[i] >= rmin and r[i] <= rmax:
            if r[i]>=(xmin+a1) and r[i]<=(xmax-a1):
                wind[i] = 1.0
            if r[i] < (xmin+a1):
                wind[i] = 0.5*(1-cos(pi*((r[i]-xmin)/a1)))
            if r[i] > (xmax-a1):
                wind[i] = 0.5*(1+cos(pi*((r[i]-xmax+a1)/a1)))
    
    return wind
    
    
def GETPHASE(yr, yi):

    n1 = 0
    n2  = 0
    i = 0
    j = 0
    p = 0
#{ Phase calculation }
    if(yr[0] < 0.0):
        n1=-1
    if(yr[0] > 0.0):
        n1=1
    if( (yr[0] == 0.0) and (yr[1]<0.0) ):
        n1=1
    if( (yr[0] == 0.0) and (yr[1]>0.0) ):
        n1=-1
    n2 = n1
    p = 0
    i = 0
    
    fi = []
    for j in range(0, len(yr)):
        fi.append(0)

    for j in range(0, len(yr)):
        if (yr[j] < 0.0):
            n1=-1
        if (yr[j] > 0.0): 
            n1=+1;
        if (yr[j] == 0.0):
            n1=-n1
        if (n1 != n2):
            p=p+1
        n2 = n1
        if (yr[j] > yr[i]):
            i=j
        if (yr[j] != 0.0):
            fi[j] = arctan(yi[j]/yr[j])+pi/2 + p*pi
        else:
            fi[j] = p * pi
    
    if (sin(fi[i]) <= 0.0):
       for j in range(0, len(yr)):
           fi[j]=fi[j] + pi;
    
    return fi
