#!/usr/local/bin/env python

import os, sys
import subprocess as sub
import shutil
import multiprocessing as mp


def shelxc_infile(ref_file, outfile, cell, symm, sites):
    fh = open(outfile, 'w')
    fh.write('SAD ' + ref_file + '\n')
    fh.write('CELL ' + cell + '\n')
    fh.write('SPAG ' + symm + '\n')
    fh.write('FIND ' + str(sites) + '\n')
    fh.write('SFAC S\n')
    fh.write('MAXM 1000\n')
    fh.close()

def run_shelxc(project, inp):
    sub.call(["shelxc "+ project + "< "+inp + "| tee " + project + "-shelxc.log"], shell=True)


def create_copies(ifh, ofh):
    ifh = ifh + '_fa.hkl'
    ofh = ofh + '_fa.hkl'
    shutil.copyfile(ifh, ofh)
#    shutil.move(ifh, ofh)


def shelxd_input(filename, resolution, iters):
    '''function that modifies .ins file for running shelxd.
    '''
    
    os.rename(filename, filename+'~')
    new = open(filename, 'w')
    
    source = open(filename+'~', 'r')
    all_lines = source.readlines()
    source.close()
    all_lines.pop(7); all_lines.pop(9); #remove SHEL and NTRY lines..

    for ii in xrange(len(all_lines)):
        line = all_lines[ii].split()
        new.write(all_lines[ii])
        if 'FIND' in line:
            
            new.write('SHEL 999 ' + str(resolution)+'\n')
        
        if 'PATS' in line:
            
            new.write('NTRY '+ str(iters) + '\n')
    new.close()
   
    sub.call(["rm "+filename+"~"], shell=True)


def create_ins(base, cell, symm, site, reso, iters):
    name = base + '.ins'
    new = open(name, 'w')
    new.write('TITL '+ name +' SAD in '+ symm + '\n')
    new.write('CELL 0.98000 ' + cell +'\n')
    new.write('LATT -1' + '\n')
    new.write('SYMM -X, 1/2+Y, -Z' + '\n')
    new.write('SFAC S' + '\n')
    new.write('UNIT 64' + '\n')
    new.write('FIND ' + str(site) + '\n')
    new.write('SHEL 999 ' + str(reso) + '\n')
    new.write('MIND -1.5 -0.1' + '\n') 
    new.write('PATS' + '\n')
    new.write('NTRY ' + str(iters) + '\n')
    new.write('SEED 1' + '\n')
    new.write('HKLF 3' + '\n')
    new.write('END' + '\n')
    new.close()


def run_shelxd(filename, outname):
    sub.call(["shelxd -t2 " + filename + " | tee " + outname], shell=True)

def get_20percent(total_sites):
    val = total_sites*0.2
    return round(val)

def reso_range(start, stop, step):
    i = start
    while i <= stop:
          yield i
          i += step
    

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


def par_shelxc(ref_file, cell, symm, tot_sites):
    proc = [];
    #project = 'tr-'+str(sites)
    for sites in tot_sites:
        project = 'tr-'+str(sites)
        inp = str(sites)+'.inp'
        shelxc_infile(ref_file, inp, cell, symm, sites)
        proc.append(mp.Process(target=run_shelxc, args=(project,inp)))
    for p in proc:
        p.start()
    for p in proc:
        p.join()

def par_shelxd(tot_sites):

    proc = []; #pool = mp.Pool(processes=4)
    for sites in tot_sites:
        for res in reso_range(reso, reso+2, 0.5):
            for trial in ntry:
                base_c = 'tr-'+str(sites)
                project = base_c +'-'+str(res)+'-'+str(trial)
                create_copies(base_c, project)
                create_ins((project+'_fa'), cell, symm, sites, res, trial)
                name = project + '_fa'
                outname = project +'-shelxd.log'
                #pool = mp.Pool(processes=4)
                proc.append(mp.Process(target=run_shelxd, args=(name, outname)))
              #  pool.apply(run_shelxd, args=(name, outname))
    #pool.close()

    for p in proc:
        p.start()

    for p in proc:
        p.join()

ref_file = sys.argv[1]
#total_sites = int(sys.argv[2])
tot_sites = [50, 75, 100, 125, 140]
reso = 3.8
ntry = [5000, 10000]

cell = "134.58   301.96   144.35  90.000 113.372  90.000"
symm = "P21"

par_shelxc(ref_file, cell, symm, tot_sites)
par_shelxd(tot_sites)



   
