#!/usr/local/bin env python

import sys, os
import numpy as np
import matplotlib.pyplot as plt
import argparse
import glob

parser = argparse.ArgumentParser()
parser.add_argument("--graphics", type=str,
                          help='Which shelx plots you want?, e.g. contrast, chi_sq, CCall_CCweak, or site_occupancy')
parser.add_argument("--input", type=str, nargs='+',
                          help='provide correct file names with path if not in the same folder.')
parser.add_argument("--row", type=str, 
                          help='Number of rows in cluster CCs.')
parser.add_argument("--col", type=str, 
                          help='Number of cols in cluster CCs.')

args = parser.parse_args()


def get_string(filename, string):
    file = open(filename, 'r')
    all_lines = file.readlines()
    found = False
    select = []
    for lines in all_lines:
        line = lines.split()
        if string in line:
            found = True
            select.append(line)
    return select

def list_file(dirname):

    file_list = [];
    os.chdir(dirname)
    for file in glob.glob("tr-*shelxd.log"):
        file_list.append(file)
    file_list.sort()
    
    return file_list

            
def plot_shelxe(data1, data2, stringX, stringY):
    
    plt.plot(data1, '-ro', linewidth=2)
    plt.plot(data2, '-o', linewidth=2)
    plt.xlabel(stringX, fontsize=24, fontweight='bold')
    plt.ylabel(stringY, fontsize=24, fontweight='bold')
    plt.legend(['Original', 'Inverted'], loc='upper left')
    plt.show()


def plot_shelxd(data1, data2, stringX, stringY):
    if args.graphics == 'site_occupancy':
       plt.plot(data1, data2, '-o', linewidth=2)
     
    if args.graphics == 'CCall_CCweak':
       plt.plot(data1, data2, 'o', linewidth=2)
    plt.xlabel(stringX, fontsize=24, fontweight='bold')
    plt.ylabel(stringY, fontsize=24, fontweight='bold')
    plt.show()

def plot_shelxc(data1, data2, stringX, stringY):

    plt.plot(data1, data2, '-o', linewidth=2)
    plt.gca().invert_xaxis()
    plt.xlabel(stringX, fontsize=24, fontweight='bold')
    plt.ylabel(stringY, fontsize=24, fontweight='bold')
    plt.show()

def get_contrast(filename):

    contrast = get_string(filename, 'Contrast')
    cont_lst = []
    for letter in contrast:
        tmp = letter[5].strip(',')
        cont_lst.append(tmp)
    return cont_lst

def get_CCs(filename):
   
    CCs = get_string(filename, 'All/Weak')
    CCall = []; CCweak = [];

    # a nasty way of sanity checks on CCs list if it got some funny string pattern 

    for line in CCs:
        sub = 'CFOM-'
        for s in line:
            if sub in s:
               CCs.remove(line)
    # hope now it should be clean..
    for line in CCs:       
        if len(line) == 15:
           
           CCall.append(line[6].strip(','))
           CCweak.append(line[8].strip(','))
        else:
           CCall.append(line[5].strip(','))
           CCweak.append(line[7].strip(','))

    return CCall, CCweak

def get_sites(filename):

    sites = get_string(filename, '1')
    cfom = get_string(filename, 'CFOM')
    for vals in cfom:
        CCall = vals[5].strip(',') + '%'
        CCweak = vals[7].strip(',')+ '%'
        fom = vals[9].strip(',')+ '%'
    
    occu = []; peaks = [];
    count = 0;
    for line in sites:
        
        tmp = line[5].strip(',')
        occu.append(tmp)
        count += 1
        peaks.append(count)

    return occu, peaks, CCall, CCweak, fom

def get_shelxc_stats(filename):
    metrics = ['Resl.','Chi-sq', '<I/sig>', '%Complete', 'Multipl.', 'R(pim)%', '<d"/sig>', 'CC(1/2)']
    stats = []; tmp = []; reso = [];
    for ele in metrics:
        line = get_string(filename, ele)
        stats.append(line)
    for vals in stats:
        for val in vals:     
            tmp.append(val)

    tmp[0].pop(0); tmp[0].pop(0);
    tmp[1].pop(0); tmp[2].pop(0);
    tmp[3].pop(0); tmp[4].pop(0);
    tmp[5].pop(0); tmp[6].pop(0);
    tmp[8].pop(0); tmp.pop(7)
    
    stat_array = np.array(tmp, dtype=float)
    
    return stat_array

def plot_cc(folder, rows, cols):
    lists = list_file(folder)
    nums = len(lists)
    fig, axs = plt.subplots(rows,cols, sharex=True, sharey=True, squeeze=True, facecolor='w', edgecolor='k')
    fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1, hspace = 0.5, wspace=0.1)

    fig.text(0.5, 0.04, 'CCweak', ha='center')
    fig.text(0.04, 0.5, 'CCall', va='center', rotation='vertical')

    axs = axs.flatten()
    print axs.shape

    for i in range(nums):
        ccall, ccweak = get_CCs(lists[i])
        axs[i].plot(ccall, ccweak, 'o')
        titles = lists[i].strip('-shelxd.log')
        axs[i].set_title(titles)

    plt.show()

def main():
   
   if args.graphics is None:
      sys.exit('tell me what you wanna plot. --help\n')
   if args.input is None:
      sys.exit('dude! provide data to plot --help\n')
   
    
   if args.graphics == 'contrast':
      if not len(args.input) > 1:
         sys.exit('Need two file names\n')
      else:
         orig = get_contrast(args.input[0])
         invert = get_contrast(args.input[1])

         plot_shelxe(orig, invert, 'Cycle', 'Contrast')

   elif args.graphics == 'CCall_CCweak':
        if not len(args.input) == 1:
           sys.exit('Need one file only \n')
        else:
           CCall, CCweak = get_CCs(args.input[0])
   #        print CCall, CCweak
           plot_shelxd(CCweak, CCall, 'CCweak', 'CCall')
   
   elif args.graphics == 'cluster':
        if not len(args.input) == 1:
           sys.exit('Need a folder name \n')
        
        elif (len(args.row) or len(args.col)) == None:
             sys.exit('provide me the number of rows and columns for cluster plotting \n')
        else:
             plot_cc(args.input[0], int(args.row), int(args.col))

   elif args.graphics == 'site_occupancy':
        if not len(args.input) == 1:
           sys.exit('Need one file only.. *.res file \n')
        else:
          
           occu, peaks, CCall, CCweak, fom  = get_sites(args.input[0])
           print('CCall: %s, CCweak: %s, CFOM: %s\n' % (CCall, CCweak, fom))
           plot_shelxd(peaks, occu, 'Peak Number', 'Site Occupancy')

   elif args.graphics == 'chi_sq':
        if not len(args.input) == 1:
           sys.exit('Need one file only.. shelxc.log file \n')
        else:
           stat_array = get_shelxc_stats(args.input[0])
           plot_shelxc(stat_array[0,:], stat_array[1,:], 'Resolution', 'Chi^2')
   
   elif args.graphics == 'I/sig':
        if not len(args.input) == 1:
           sys.exit('Need one file only.. shelxc.log file \n')
        else:
           stat_array = get_shelxc_stats(args.input[0])
           plot_shelxc(stat_array[0,:], stat_array[2,:], 'Resolution', 'I/sig')

   elif args.graphics == 'completeness':
        if not len(args.input) == 1:
           sys.exit('Need one file only.. shelxc.log file \n')
        else:
           stat_array = get_shelxc_stats(args.input[0])
           plot_shelxc(stat_array[0,:], stat_array[3,:], 'Resolution', 'Completeness(%)')

   elif args.graphics == 'Multiplicity':
        if not len(args.input) == 1:
           sys.exit('Need one file only.. shelxc.log file \n')
        else:
           stat_array = get_shelxc_stats(args.input[0])
           plot_shelxc(stat_array[0,:], stat_array[4,:], 'Resolution', 'Multiplicity')

   elif args.graphics == 'Rpim':
        if not len(args.input) == 1:
           sys.exit('Need one file only.. shelxc.log file \n')
        else:
           stat_array = get_shelxc_stats(args.input[0])
           plot_shelxc(stat_array[0,:], stat_array[5,:], 'Resolution', 'Rpim')

   elif args.graphics == 'Anomalous_signal':
        if not len(args.input) == 1:
           sys.exit('Need one file only.. shelxc.log file \n')
        else:
           stat_array = get_shelxc_stats(args.input[0])
           plot_shelxc(stat_array[0,:], stat_array[6,:], 'Resolution', 'd/sig')

   elif args.graphics == 'CC1/2':
        if not len(args.input) == 1:
           sys.exit('Need one file only.. shelxc.log file \n')
        else:
           stat_array = get_shelxc_stats(args.input[0])
           plot_shelxc(stat_array[0,:], stat_array[7,:], 'Resolution', 'CC(1/2)')


if __name__ == '__main__':
   main()
       


