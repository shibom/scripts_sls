#!/usr/local/bin/env python

import sys, os
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import interpolate
from mpl_toolkits.mplot3d import Axes3D

vacuum = 0;

a = 0.22*10e-4 #um^-1A^-3

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


def air_transmission(wave, air_path): # function from Oliver

    data = np.loadtxt('air_transmission_wavelength.dat', skiprows=2)
    absp = np.empty((data.shape[0],2))

    for i in range(absp.shape[0]):
        absp[i,0] = 10*data[i,0] # converting nm to Ang
        if data[i,1] > 0:
           absp[i,1] = math.exp((air_path/30) * math.log(data[i,1]))

    f = interpolate.interp1d(absp[:,0], absp[:,1])

    return f(wave)

def water_absp(wave): # function from Oliver

    data = np.loadtxt('water_attenuationlength_wavelength.dat', skiprows=2)
    absp = np.empty((data.shape[0],2))

    for i in range(absp.shape[0]):
        absp[i,0] = 10*data[i,0] # converting nm to Ang
        absp[i,1] = data[i,1]

    f = interpolate.interp1d(absp[:,0],absp[:,1])
    return f(wave)

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

def efficiency(wave, t_xtal):
    Eff = []; mu_2A = a*(2**3)
    #norm = (2**3)*((t_xtal)**3)/(math.exp(mu_2A*t_xtal) - 1)

    for li in xrange(0,len(wave)):

        mu_xtal = a*(wave[li]**3)
        Ie = (wave[li]**3)*(t_xtal**1)/(math.exp(mu_xtal*t_xtal) - 1)
        #Ie = Ie/norm
        Eff.append(Ie)
    return Eff


def fprime(wave): # function from Oliver

    data = np.loadtxt('sulphur_fprime_and_double_prime.dat')
    data = data[::-1]
    fpp = np.empty((data.shape[0],2))
    for i in range(np.shape(fpp)[0]):
        fpp[i,0] = 12.39852066*1e3 / data[i,0]

        fpp[i,1] = data[i,2]

    fpp_ = interpolate.interp1d(fpp[:,0], fpp[:,1])

    fpp_y = fpp_(wave)
    return fpp_y


def form_fac(wave): # analytical function by me

    f2 = []; sulfur = 52.62 #cm^-1 A^-3
    re = 2.82*10e-13 # cm
    density = 2.05 #gm cm^-3
    NA = 6.023*10e23; MW = 32.066
    const = (sulfur*MW) / (2*re*density*NA) #cm A^-3
    const = const*10e8 #A^-2
    print const

    for ii in xrange(0, len(wave)):
        tmp = const*(wave[ii]**2)
        f2.append(tmp)
    return f2

xray = []; xtal = [];
for li in steps(0.5,6,0.1):
    xray.append(li)

for jj in steps(1,200,2):
    xtal.append(jj)

trans_air = transmit_air(xray, 10)
abs_si = absorp_silicon(xray, 320)

#abs_si = silicon_abs(xray, 320)
#trans_air = air_transmission(xray, 20)

if vacuum == 1:
   trans_air = np.ones((len(xray)))

diff = efficiency(xray, 100)
'''
diff = np.asarray(diff); xray = np.asarray(xray); xtal = np.asarray(xtal)
X, Y = np.meshgrid(xray, xtal)
diff = diff.reshape(X.shape)
'''
#f2 = form_fac(xray)

f2 = fprime(xray)

diff_f2 = [diff[ii]*f2[ii] for ii in range(len(diff))];
diff_f2_si = [diff[ii]*f2[ii]*abs_si[ii] for ii in range(len(diff))];

diff_sqf2 = [diff[ii]*(f2[ii]**2) for ii in range(len(diff))];
diff_sqf2_air = [diff_sqf2[ii]*trans_air[ii] for ii in range(len(diff))];
diff_sqf2_si = [diff_sqf2[ii]*abs_si[ii] for ii in range(len(diff))];
diff_sqf2_air_si = [diff_sqf2[ii]*trans_air[ii]*abs_si[ii] for ii in range(len(diff))];
diff_air = [diff[ii]*trans_air[ii] for ii in range(len(diff))];
diff_air_si = [diff[ii]*trans_air[ii]*abs_si[ii] for ii in range(len(diff))];

'''
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, diff)
plt.show()
'''
plt.plot(xray, diff, linewidth=2)
plt.plot(xray, diff_f2, linewidth=2)
plt.plot(xray, diff_f2_si, linewidth=2)

#plt.plot(xray, diff_sqf2, linewidth=2)

#plt.plot(xray, diff_sqf2_air, linewidth=2)
#plt.plot(xray, diff_sqf2_si, linewidth=2)
#plt.plot(xray, diff_sqf2_air_si, linewidth=2)
#plt.plot(xray, diff_air_si, linewidth=2)

plt.legend(['Intensity', 'Intensity_f2', 'Intensity_f2_si'], loc= 'upper left')
plt.xlabel('Wavelength (Ang)')
plt.ylabel('I')

#plt.legend(['Intensity', 'Intensity*f2', 'Intensity*f2^2', 'Intensity*f2^2*air', 'Intensity*f2^2*si', 'Intensity*f2^2*air*si', 'Intensity*air*si'], loc= 'upper left')
plt.show()

#wagner paper..
diff_50 = efficiency(xray, 50)
diff_20 = efficiency(xray, 20)
diff_300 = efficiency(xray, 300)
diff_200 = efficiency(xray, 200)

#plotting only intensities at various crystal sizes in vacuum..
plt.plot(xray, diff_20)
plt.plot(xray, diff_50)
plt.plot(xray, diff)
plt.plot(xray, diff_200)
plt.plot(xray, diff_300)
plt.legend(['20um', '50um', '100um', '200um','300um'], loc= 'upper left')
plt.xlabel('Incident Wavelength(A)', fontweight = 'bold', fontsize=20)
plt.ylabel('Intensity', fontweight='bold', fontsize=20)
plt.title('I, t=1, d_water=0, air=0', fontweight='bold', fontsize = 20)
plt.grid(True)
plt.show()

#calculating anomalous intensities now..
diff_20_f2 = [diff_20[ii]*f2[ii]*trans_air[ii]*abs_si[ii]*4.22 for ii in range(len(diff))];
diff_50_f2 = [diff_50[ii]*f2[ii]*trans_air[ii]*abs_si[ii]*4.22 for ii in range(len(diff))];
diff_100_f2 = [diff[ii]*f2[ii]*trans_air[ii]*abs_si[ii]*4.22 for ii in range(len(diff))];
diff_200_f2 = [diff_200[ii]*f2[ii]*trans_air[ii]*abs_si[ii]*4.22 for ii in range(len(diff))];
diff_300_f2 = [diff_300[ii]*f2[ii]*trans_air[ii]*abs_si[ii]*4.22 for ii in range(len(diff))];

plt.plot(xray, diff_20_f2)
plt.plot(xray, diff_50_f2)
plt.plot(xray, diff_100_f2)
plt.plot(xray, diff_200_f2)
plt.plot(xray, diff_300_f2)
plt.legend(['20um', '50um', '100um', '200um','300um'], loc= 'upper left')
plt.xlabel('Incident Wavelength(A)', fontweight = 'bold', fontsize=20)
plt.ylabel(r'$\delta$ I', fontweight='bold', fontsize=20)
plt.title('I*F, t=1, d_water=0um, air=10cm', fontweight='bold', fontsize = 20)
plt.grid(True)
plt.show()
