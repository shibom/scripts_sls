
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


def shelxd_input(outname, inname, resolution, iters, emins):
    '''function that modifies .ins file for running shelxd.
    '''

    #os.rename(filename, filename+'~')
    new = open(outname, 'w')

    source = open(inname, 'r')
    all_lines = source.readlines()
    source.close()
    count_symm = 0
    for lines in all_lines:
        each_line = lines.split()
        if 'SYMM' in each_line:
            count_symm += 1;
    if count_symm == 1:
        all_lines.pop(7); all_lines.pop(9);
    elif count_symm == 2:
          all_lines.pop(8); all_lines.pop(10)
    elif count_symm == 3:
          all_lines.pop(9); all_lines.pop(11);
    elif  count_symm == 4:
          all_lines.pop(10); all_lines.pop(12);
    elif  count_symm == 4:
          all_lines.pop(10); all_lines.pop(12);
    elif  count_symm == 11:
          all_lines.pop(17); all_lines.pop(19);
   # all_lines.pop(11); all_lines.pop(13); #remove SHEL and NTRY lines..buggy. case when multiple SYM lines appear.
    #all_lines.pop(9); all_lines.pop(11); #remove SHEL and NTRY lines..

    for ii in xrange(len(all_lines)):
        line = all_lines[ii].split()
        new.write(all_lines[ii])
        if 'FIND' in line:

            new.write('SHEL 999 ' + str(resolution)+'\n')
            new.write('ESEL ' + str(emins)+ '\n')
        if 'PATS' in line:

            new.write('NTRY '+ str(iters) + '\n')
    new.close()

    #sub.call(["rm "+filename+"~"], shell=True)


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
    sub.call(["shelxd " + filename + " | tee " + outname], shell=True)

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
        run_shelxc(project,inp)


def par_shelxd(tot_sites):


    for sites in tot_sites:
        for res in reso:
            for emin in eval:
                for trial in ntry:
                    base_c = 'tr-'+str(sites)
                    project = base_c +'-'+str(res)+'-'+str(emin)+'-'+str(trial)
                    create_copies(base_c, project)
                    #create_ins((project+'_fa'), cell, symm, sites, res, trial)
                    inname = base_c + '_fa' + '.ins'
                    name = project + '_fa'
                    outname = name + '.ins'
                    logname = project +'-shelxd.log'
                    shelxd_input(outname, inname, res, trial, emin)
                    run_shelxd(name, logname)


ref_file = sys.argv[1]
tot_sites = [int(sys.argv[2])]
reso = [float(sys.argv[3])]; eval = [float(sys.argv[4])]
ntry = [int(sys.argv[5])]

#cell = "104.83   158.82   180.22  90.000  90.000  90.000"
cell = "103.77   155.79   179.72  90.000  90.000  90.000"
symm = "P212121"

par_shelxc(ref_file, cell, symm, tot_sites)
par_shelxd(tot_sites)
