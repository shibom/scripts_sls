#!/usr/local/bin/env python

import sys, os
import matplotlib.pyplot as plt

 
def extract_xscale(filename):
    file = open(filename, 'r')
    all_lines = file.readlines()
    file.close()
    #look for indexes of the specific lines
    
    start = all_lines.index(' SUBSET OF INTENSITY DATA WITH SIGNAL/NOISE >= -3.0 AS FUNCTION OF RESOLUTION\n')
    end = all_lines.index(' ========== STATISTICS OF INPUT DATA SET ==========\n')
    # read only the selected chunk ..
    select = all_lines[start:end]
    # prune a bit of junks and then store everything else in a list
    useful = [];
    for lines in select:
       if 'SUBSET' in lines:
          pass
       elif 'RESOLUTION' in lines or 'LIMIT' in lines:
             pass
       elif 'total' in lines:
            pass
       else:
           line = lines.split()
           if len(line) > 0:
              useful.append(line)   
    return useful

def get_stats(data_list, min_range, max_range):
    reso = []; Isig = []; Ano = []; Rfactor = [];
    norm = min_range + 2
    for ii in xrange(min_range, max_range):
        reso.append(float(data_list[ii][0]))
        Isig.append(float(data_list[ii][8]))
        Ano.append(float(data_list[ii][12]))
        Rfactor.append(float(data_list[ii][5].strip('%')))

    return reso, Isig, Ano, Rfactor

def wrap(file_list):
    in_list = open(file_list)
    for line in in_list:
        line = line.strip('\n')
        # extract data structure from each xscale file in the list
        data = extract_xscale(line)
        # plot them sequentially..
        reso, Isig, Ano, Rf = get_stats(data, 0, 19)

        plt.subplot(3,1,1)
        plt.plot(reso, Isig, '-o', linewidth=2)
        plt.xlabel('Resolution')
        plt.ylabel('I/sig(I)')
        plt.gca().invert_xaxis()

        plt.subplot(3,1,2)
        plt.plot(reso, Ano, '-o', linewidth=2)
        plt.xlabel('Resolution')
        plt.ylabel('Anomalous_sig')
        plt.gca().invert_xaxis()

        plt.subplot(3,1,3)
        plt.plot(reso, Rf, '-o', linewidth=2)
        plt.xlabel('Resolution')
        plt.ylabel('Rfactor')
        plt.gca().invert_xaxis()

    
    plt.show()

 
wrap('list')

