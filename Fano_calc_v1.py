import sys, os
import math
import numpy as np
import matplotlib.pyplot as plt
import subprocess as sub
import argparse

'''
@Author S.Basu..
This program calculates an essential metric, namely <Delta_F/F>, i.e., signal from anomalous amplitude..
if bugs are encountered, email shibom.basu@psi.ch
'''

parser = argparse.ArgumentParser(prog='Fano_calc.py', usage='python %(prog)s --mtzfile mtz-filename --cell cell-parameters --high_res high-resolution')
parser.add_argument("--mtzfile", type=str,
                        help='provide me a merged and scaled mtz file with amplitudes')
parser.add_argument("--cell", type=str, nargs='+',
                        help='provide me the unit cell info..a b c alpha beta gamma')
parser.add_argument("--low_res", type=str,
                        help='provide me the low-resolution limit in Angstrom. Default is 20 A')
parser.add_argument("--high_res", type=str,
                        help='provide me the high-resolution limit in Angstrom')
parser.add_argument("--res_bins", type=str,
                        help='provide resolution bin number. Default is 20')

args = parser.parse_args()


def convert_mtz(mtzfilename, rmin, rmax):
    fh = open('mtzinp', 'w')

    fh.write('#!/bin/bash \n')
    fh.write('set e \n')
    fh.write('mtz2various HKLIN '+ mtzfilename + ' HKLOUT cns.hkl ' + '<< eof \n')
    fh.write('LABIN F(+)=F(+) SIGF(+)=SIGF(+) F(-)=F(-) SIGF(-)=SIGF(-) \n')
    fh.write('OUTPUT CNS\n')
    fh.write('RESOLUTION '+ str(rmin) + ' ' + str(rmax) +'\n')
    fh.write('END\n')
    fh.write('eof\n')
    fh.close()

    sub.call(["chmod +x mtzinp"], shell=True)
    sub.call(["./mtzinp"], shell=True)

    os.remove('mtzinp')

    return

def createDict(keys, values):

    dictionary = {}
    i = 0
    for key in keys:
        dictionary[key] = values[i]
        i += 1
        if i > len(values):
           dictionary[key] = 'None'
    return dictionary

def readcns(filename, rmin, rmax):

    convert_mtz(filename, rmin, rmax)

    fh = open('cns.hkl')
    all_lines = fh.readlines()
    select = []; start = 0

    for lines in all_lines:
        line = lines.split()
        if 'FOBS=' in line:
            select.append(line)

    os.remove('cns.hkl')

    return select

def readxds(filename):

    fh = open(filename)
    all_lines = fh.readlines()
    select = [];

    for lines in all_lines:
        if "!" in lines:
            pass
        else:
            line = lines.split()
            select.append(line)

    h = []; k = []; l = []; F = [];
    for vals in select:
        h.append(int(vals[0]))
        k.append(int(vals[1]))
        l.append(int(vals[2]))
        intensity = abs(float(vals[3]))

        F.append(math.sqrt(intensity))
    indices = zip(h,k,l)

    return indices, F

def extract_F(data_from_cns):
    h = []; k = []; l = [];
    F = []; sigF = [];

    for vals in data_from_cns:
        h.append(int(vals[1]))
        k.append(int(vals[2]))
        l.append(int(vals[3]))

        F.append(float(vals[5]))
        sigF.append(float(vals[8]))
    indices = zip(h,k,l)
    fridel = []; Fpm = [];

    for j in range(len(indices)-1):
        mate = (-indices[j][0],-indices[j][1],-indices[j][2])
        if indices[j+1] == mate:
           fridel.append((indices[j],indices[j+1]))
           Fpm.append((F[j],F[j+1]))


    return fridel, Fpm

def set_bins(nbins, reslim1, reslim2):

    stolmax3 = (1.0/(2*reslim2))**3
    stolmin3 = (1.0/(2*reslim1))**3
    stol2min = math.exp(math.log(stolmin3)*2.0/3)

    stolinc = (stolmax3 - stolmin3)/nbins
    stol2max = np.zeros((nbins)); resol_range = np.zeros((nbins));
    for i in range(nbins):
        stol3max = stolinc*i + stolmin3
        stol2max[i] = math.exp(math.log(stol3max)*2.0/3)
        resol_range[i] = math.sqrt(1/stol2max[i])/2.0

    return stol2max, resol_range

def get_bins(dstar, q_sq_array, nbins, ibin=0):

    for i in range(0,nbins):
        if dstar <= q_sq_array[i]:
           break
        else: ibin = 0
    ibin = i

    return ibin

def calc_resolution(ind_p, F, rmin, rmax, nbins, cell):

    a = float(cell[0]); b = float(cell[1]); c = float(cell[2])
    al = float(cell[3]); be = float(cell[4]); ga = float(cell[5])

    d_star = []
    cosa = abs(math.cos(al)); cosb = abs(math.cos(be)); cosg = abs(math.cos(ga));
    sina = abs(math.sin(al)); sinb = abs(math.sin(be)); sing = abs(math.sin(ga));
    cosas = (cosb*cosg-cosa)/(sinb*sing)
    cosbs = (cosa*cosg-cosb)/(sina*sing)
    cosgs=(cosa*cosb-cosg)/(sina*sinb)

    vol=a*b*c*math.sqrt(1.0-cosa**2-cosb**2-cosg**2+2*cosa*cosb*cosg)
    ast=b*c*sina/vol
    bst=a*c*sinb/vol
    cst=a*b*sing/vol

    #calculate resolution bins and put hkl reflections into correct resolution bins..

    q_sq_array, resolution = set_bins(nbins, rmin, rmax)

    putbin = []; anoFs = [];

    for ii in range(len(ind_p)):

        indx = ind_p[ii][0]
        h = indx[0]; k = indx[1]; l = indx[2]
        s2 = ((h*ast)**2 + (k*bst)**2 + (l*cst)**2 + 2.0*h*k*ast*bst*cosgs + 2.0*h*l*ast*cst*cosbs + 2.0*k*l*bst*cst*cosas)/4

        bin_num = get_bins(s2,q_sq_array, nbins) #find out which hkl reflection belongs to which resolution bin
        putbin.append(bin_num)
        ds = math.sqrt(s2)
        resol = 1/(2*ds)
        d_star.append(resol) #storing resolution associated with each hkl just in case if needed in future..

        '''
        calculate Delta_of_F for each hkl within the same loop and then store them in a list.
        At the end, create a dictionary, which will have Delta_of_Fs for each hkl as key and values will be the bin-number of those hkls.
        '''
        F_hkl = F[ii]
        delF = abs(F_hkl[0] - F_hkl[1])
        sumF = math.sqrt(abs((F_hkl[0]**2) + (F_hkl[1])**2))

        anoFs.append((delF,sumF))

    data_dict = createDict(anoFs, putbin)

    return data_dict

def calc_anomalous(data_dict, nbins):
    '''
    Calculate <deltaF_by_F> and put them into correct resolution bins. data_dict has keys as deltaF and values as bin-number of each hkls.
    Here, we need to count how many times each bin-number occurred and pull out the deltaF (i.e., keys) for each bin-number and average them off.
    '''
    factors = [];

    for i in range(nbins):
        delta = []; sums = [];
        count = 0;
        for k, v in data_dict.iteritems():
            if i == v:
               count += 1
               delta.append(k[0])
               sums.append(k[1])
        if len(sums) > 0:
            factors.append(sum(delta)/sum(sums))

    return factors

def main_calc_plot(fname, rmin, rmax, nbins, cell):

    select = readcns(fname, rmin, rmax)


    ind_p, F = extract_F(select)


    data_dict = calc_resolution(ind_p, F, rmin, rmax, nbins, cell)

    q, res = set_bins(nbins, rmin, rmax)
    anomalous_signal = calc_anomalous(data_dict, nbins)
    anom_ar = np.array(anomalous_signal)

   # print len(res)
   # print len(anomalous_signal)

    if len(anomalous_signal) < len(res):
       diff = len(res) - len(anomalous_signal)
       fit = np.polyfit(res[diff:len(res)], anom_ar, 1)
       #plt.plot(res[diff:len(res)], anomalous_signal, 'o')
       plt.plot(res[diff:len(res)], fit)
    else:
        fit = np.polyfit(res, anom_ar, 1)
        print fit.shape
        #plt.plot(res, anomalous_signal, '-o')
        plt.plot(fit)

    plt.gca().invert_xaxis()
    plt.xlabel("Resolution ($\AA$)", fontsize=14, fontweight='bold')
    plt.ylabel("<$\Delta$F/F>", fontsize=14, fontweight='bold')
    plt.show()




def main():

    if args.mtzfile is None:
       sys.exit('Need a mtz file.. --help\n')
    else:
       mtz_file = args.mtzfile

    if not len(args.cell) == 6:
           sys.exit('provide a,b,c and angles. --help\n')
    else:
        cell = args.cell

    if args.low_res is None:
       args.low_res = 20
       rmin = float(args.low_res)
    else:
       rmin=float(args.low_res);

    if args.high_res is None:
       sys.exit('provide a high-resolution limit. --help\n')
    else:
       rmax=float(args.high_res);

    if args.res_bins is None:
       args.res_bins=20

    nbins = int(args.res_bins)
    main_calc_plot(mtz_file, rmin, rmax, nbins, cell)

if __name__ == '__main__':
    main()
