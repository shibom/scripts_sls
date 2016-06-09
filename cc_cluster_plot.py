#!/usr/local/bin/env python

import sys, os
import numpy as np
import matplotlib.pyplot as plt
import plot_shelx as psx
import glob

def list_file(dirname):

    file_list = [];
    os.chdir(dirname)
    for file in glob.glob("tr-*shelxd.log"):
        file_list.append(file)
    file_list.sort()
    print file_list
    return file_list

def plot_ano():
    lists = list_file("shelx_v1")
    nums = len(lists)
    fig, axs = plt.subplots(3,2, sharex=True, figsize=(15, 6), facecolor='w', edgecolor='k')
    fig.subplots_adjust(hspace = .5, wspace=.001)
    axs[0,0].invert_xaxis()

    axs = axs.ravel()
    for i in range(len(lists)):
        data = psx.get_shelxc_stats(lists[i])
        axs[i].plot(data[0,:], data[6,:], '-o', linewidth=2)
        plt.gca().invert_xaxis()
        axs[i].set_xlabel('Resolution', fontsize=24, fontweight='bold')
        axs[i].set_ylabel('d/sig', fontsize=24, fontweight='bold')
        axs[i].set_title(lists[i])

    plt.show()     

def plot_cc(folder):
    lists = list_file(folder)
    nums = len(lists)
    fig, axs = plt.subplots(7,4, sharex=True, sharey=True, squeeze=True, facecolor='w', edgecolor='k')
    fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1, hspace = 0.5, wspace=0.1)
    
    fig.text(0.5, 0.04, 'CCweak', ha='center')
    fig.text(0.04, 0.5, 'CCall', va='center', rotation='vertical')

    axs = axs.flatten()
    print axs.shape

    for i in range(nums):
        ccall, ccweak = psx.get_CCs(lists[i])
        axs[i].plot(ccall, ccweak, 'o')
       # axs[i].set_xlabel('CCweak', fontsize=24, fontweight='bold')
       # axs[i].set_ylabel('CCall', fontsize=24, fontweight='bold')
        titles = lists[i].strip('-shelxd.log')
        axs[i].set_title(titles)
    
   # fig.tight_layout()
    plt.show()

folder = "/sls/X06DA/data/e11206/Data10/shibom/20160416_20S/merge/XSCALE_21_09/shelx_v1"
plot_cc(folder)

