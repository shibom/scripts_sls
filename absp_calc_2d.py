#!/usr/local/bin/env python

import sys, os
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import interpolate
from mpl_toolkits.mplot3d import Axes3D

#code to calculate optimum Wavelength and crystal size for anomalous signal
#@Author by S.Basu

a = 0.22*10e-4 #um^-1A^-3 (value from Arndt 1984 paper)

def steps(start, stop, step):
    i = start
    while i <= stop:
          yield i
          i += step


def transmit_air(wave, air_path): # analytical function by me

    air = 0.0033; # cm^-1A^-3
    trans_air = [];

    for li in xrange(0,len(wave)):
        mu_air = air*(wave[li]**3)
        trans = math.exp(-(mu_air*air_path))
        trans_air.append(trans)
    return trans_air

def absorp_silicon(wave, si_thick): # analytical function by me
    si = 39.36*10e-4 #um^-1 A^-3
    abs_si = [];

    for ii in xrange(0, len(wave)):
        mu_si = si*(wave[ii]**3)
        tmp = 1 - math.exp(-(mu_si*si_thick))
        abs_si.append(tmp)
    return abs_si


def silicon_abs(wave, si_thickness): # function from Oliver

    data = np.loadtxt('si_transmission_wavelength.dat', skiprows=2)
    absp = np.empty((data.shape[0],2))
    for i in range(data.shape[0]):
        absp[i,0] = 10*data[i,0] # converting nm to Ang
        if data[i,1] > 0:
           tmp = math.log(data[i,1])
           # convert absorption into transmission
           absp[i,1] = 1.0 - math.exp((si_thickness/320)*tmp)

    f = interpolate.interp1d(absp[:,0], absp[:,1])
    return f(wave)

def fprime(wave): # function from Oliver

    data = np.loadtxt('sulphur_fprime_and_double_prime.dat')
    data = data[::-1] #interp1d function needs array ordered in ascending order..
    fpp = np.empty((data.shape[0],2))
    for i in range(np.shape(fpp)[0]):
        fpp[i,0] = 12.39852066*1e3 / data[i,0]

        fpp[i,1] = data[i,2]

    fpp_ = interpolate.interp1d(fpp[:,0], fpp[:,1])

    fpp_y = fpp_(wave)
    return fpp_y


xray = []; d_xtal = [];

for ii in steps(0.5,6,0.1):
    xray.append(ii)

for ii in steps(10,300,2):
    d_xtal.append(ii)

d_water = 50;

trans_air = transmit_air(xray, 0)

f2 = fprime(xray)

xray = np.asarray(xray); d_xtal = np.asarray(d_xtal);
X, Y = np.meshgrid(xray, d_xtal);
diff = np.empty((X.shape))

#calculating integrated intensity..

for ii in range(d_xtal.shape[0]):
    for jj in range(xray.shape[0]):
        mu = a*(xray[jj]**3)
        den = 1 - math.exp(-mu*d_xtal[ii])
        diff[ii][jj] = (d_xtal[ii]**1)*(xray[jj]**3)*math.exp(-mu*(d_xtal[ii]+d_water))*f2[jj]*trans_air[jj]/den

#fancy 2d plotting...

fig = plt.figure()
ax = fig.add_subplot(111)
ab = ax.contourf(X, Y, diff, 20, extend='both')
ax.set_xlabel('Wavelength (Ang)', fontsize=20, fontweight='bold')
ax.set_ylabel('crystal size (um)', fontsize=20, fontweight='bold')
ax.set_title('I*f'',t=1,d_water=50um,air=0cm', fontsize=20, fontweight='bold')
plt.colorbar(ab,extend='both')
ax.grid(which='both', alpha=0.5, linewidth=2)
ax.grid(True)
plt.show()
