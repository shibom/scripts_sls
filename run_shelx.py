import os, sys
import subprocess as sub
import shutil
import fnmatch

class prepfiles(object):
    def __init__(self, ref_file, symm):
        self.ref_file = ref_file;
        if not self.ref_file.endswith('.sca'):
           self.cell = sub.check_output(["grep CELL_CONSTANT "+self.ref_file + " | awk '{print $2,$3,$4,$5,$6,$7}'" ], shell=True)
        else:
           #self.cell = sub.check_output(["cat "+self.ref_file + " | awk 'NR>2 && NR<4 {print $1, $2, $3, $4, $5, $6}'"], shell=True) 
           #print self.cell
           self.cell = "105.6469  160.4173  179.6551   90.0000   90.0000   90.0000" 
        self.symm = symm;

    def shelxc_infile(self, outfile):
        fh = open(outfile, 'w')
        fh.write('SAD ' + self.ref_file + '\n')
        fh.write('CELL ' + self.cell + '\n')
        fh.write('SPAG ' + self.symm + '\n')
        fh.write('MAXM 1000\n')
        fh.close()

    def get_20percent(self, total_sites):
        val = total_sites*0.2
        return int(val)

    def custom_range(self, start, stop, step):
        i = start
        while i <= stop:
              yield i
              i += step

    def create_paramspace(self, resolution, total_sites):
        self.sites = []; self.reso = []; self.ntries = [10000,30000];
        #self.sites = [32,36,40,44,46,50,55,57,60,65,70,74,78,80]; 
        self.super_sulfur = True;
         
        init_site = self.get_20percent(total_sites)
        for ii in self.custom_range(init_site, total_sites, 5):
            self.sites.append(ii)
        
        for ii in self.custom_range(resolution, resolution+1, 0.1):
            self.reso.append(ii)

        self.emins = [1.2,1.3,1.4,1.5]

    def create_copies(self, ifh, ofh):
        ifh = ifh + '_fa.hkl';
        ofh = ofh + '_fa.hkl'
        shutil.copyfile(ifh, ofh)

    def shelxd_input(self, outname, inname, site, resolution, iters, emins):
        '''function that modifies .ins file for running shelxd.
        '''

        new = open(outname, 'w')
        select = [];   keys = ['FIND','SHEL', 'NTRY', 'MIND', 'NTPR']
        source = open(inname, 'r')
        all_lines = source.readlines()
        source.close()
       # all_lines.pop(11); all_lines.pop(13); #remove SHEL and NTRY lines..buggy. case when multiple SYM lines appear.
    #    all_lines.pop(7); all_lines.pop(9); #remove SHEL and NTRY lines..

        for line in all_lines:

            if not any(string in line for string in keys):
               select.append(line)

        all_lines = select

        for ii in xrange(len(all_lines)):
            line = all_lines[ii]
            new.write(line)
            line = line.split()

            if 'UNIT' in line:
                new.write('FIND ' + str(site)+'\n')
                new.write('SHEL 999 ' + str(resolution)+'\n')
                new.write('ESEL ' + str(emins)+ '\n')
                new.write('MIND -3.5 3.0 \n')
                new.write('NTPR 500\n')
                if self.super_sulfur:
                   new.write('DSUL ' + str(site)+'\n')
            if 'PATS' in line:

                new.write('NTRY '+ str(iters) + '\n')
        new.close()

class running(object):
    def __init__(self,filename,symm,res,tot_sites):
        self.filename = filename;
        self.res = float(res); self.tot_sites = int(tot_sites);
        self.pp = prepfiles(self.filename, symm)
        self.pp.create_paramspace(self.res,self.tot_sites)

    def run_shelxc(self):
        project = 'tr'; infile = 'tr.inp';
        self.pp.shelxc_infile(infile);
        sub.call(["shelxc "+ project + "< "+infile + "| tee " + project + "-shelxc.log"], shell=True)

    def write_script(self, project, logout):
        fh = open("run_cluster", 'w')
        fh.write("#!/bin/bash\n")
        fh.write("\n\n")
        fh.write("shelxd %s | tee %s" %(project, logout))
        fh.close()

        sub.call(["chmod a+x run_cluster"], shell=True)

    def prep_shelxd(self):
        base_c = 'tr'; inname_d_4m_c = base_c + '_fa.ins'
        for site in self.pp.sites:
            for r in self.pp.reso:
                for emin in self.pp.emins:
                    for ntry in self.pp.ntries:
                        project = base_c+'-'+str(site)+'-'+str(r)+'-'+str(emin)+'-'+str(ntry)
                        name = project + '_fa';
                        outname = name + '.ins'
                        logname = project +'-shelxd.log'
                        self.pp.create_copies(base_c, project)
                        self.pp.shelxd_input(outname, inname_d_4m_c,site,r,ntry,emin)
                        self.write_script(name, logname)
                        sub.call(["sbatch -p day -J shelx run_cluster"], shell=True)

def main():
    job = running(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])
    job.run_shelxc()
    job.prep_shelxd()

if __name__ == '__main__':
    main()
