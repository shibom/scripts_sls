#!/usr/bin/python

import sys, os
import numpy as np
import matplotlib.pyplot as plt
import math

a = 0.22*10e-4 #um^-1A^-3
t = [10, 20, 50, 75, 100, 150, 200, 250, 300, 400, 500]
lamda = [];

for ti in t:
    tmp = (2.0/(3.0*a*ti))**(1./3) 
    lamda.append(tmp)

plt.plot(t, lamda, linewidth=2)
plt.xlabel('crystal_size', fontsize=24, fontweight='bold')
plt.ylabel('wavelength', fontsize=24, fontweight='bold')
plt.show()

#xray = [0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6]
xray = []
I_scat = [];
If2 = [];
Isqf2 = []; spot = [];

def steps(start, stop, step):
    i = start
    while i <= stop:
          yield i
          i += step

def normalize(wavelength):
    mu_xtal = a*wavelength
    d = 100
    Inorm = (wavelength**3)*(d**3)/(math.exp((mu_xtal*d) - 1))
    return Inorm

#f2 = [0.0933, 0.1212, 0.1998, 0.2307, 0.4648, 0.5428, 0.6988, 0.8548, 0.9719, 1.206, 1.401, 1.518, 1.791, 1.947, 2.064, 2.259, 2.493, 2.922, 3.117, 3.390, 3.507, 3.702, 3.975]

f2 = [0.1998, 0.2307, 0.4648, 0.5428, 0.6988, 0.8548, 0.9719, 1.206, 1.401, 1.518, 1.791, 1.947, 2.064, 2.259, 2.493, 2.922, 3.117, 3.390, 3.507, 3.702, 3.975]

for li in steps(1,5,0.2):
    xray.append(li)

for li in xrange(0,len(xray)):
    t_xtal = 100
    t_air = 30e-4
    t_si = 320
    mu_air = 0.38e-6*(xray[li]**3)
    mu_si = 0.46e-2*(xray[li]**3)
    mu = a*(xray[li]**3)

#   scatter = (li**3)*(ti**3)/(math.exp((mu*t_xtal)+(mu_air*t_air)+(mu_si*t_si)) - 1)
#    Ie_f2 = ((xray[li]**3)*(t_xtal**3)/(math.exp(mu*t_xtal) - 1))*f2[li]
#    Ie_sqf2 = ((xray[li]**3)*(t_xtal**3)/(math.exp(mu*t_xtal) - 1))*(f2[li]**2)
 
    Ie = (xray[li]**3)*(t_xtal**3)/(math.exp(mu*t_xtal) - 1)
    Ie_f2 = Ie*f2[li]*4.22; #Ie_sqf2 = Ie*(f2[li]**2)
    I_scat.append(Ie)
    If2.append(Ie_f2)
    #Isqf2.append(Ie_sqf2)

    '''
    Ispot = (xray[li]**2)*(t_xtal**3)*math.exp(mu*t_xtal)
    spot.append(Ispot)
    '''
plt.plot(xray, I_scat,linewidth=2)
plt.plot(xray, If2, linewidth=2)
#plt.plot(xray, Isqf2, linewidth=2)
#plt.legend(['I_scat','If2', 'Isqf2'])
'''
plt.plot(xray, spot, linewidth=2)
'''
plt.show()

