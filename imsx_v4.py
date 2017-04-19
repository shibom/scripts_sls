'''
This code belongs to SLS PX beamlines. It automatically looks for IMISX or serial crystallography data collected 
at Eiger or Pilatus detector, processes them in parallelized fashion, and runs XSCALE using ISa based selection.
In future, It will provide other types of data visualization and quality checks.
This project was initiated by R. Warshamanage, then developed and maintained by S. Basu
Report bugs at shibom.basu@psi.ch
'''
import os, sys, errno
import glob
import h5py
import numpy as np
import multiprocessing as mp
import argparse
'''
============First provide helps/options for users=====================
'''
parser = argparse.ArgumentParser()
parser.add_argument("--well_path", type=str, nargs='+',
                       help="provide path for each well, containing minisets folder or provide path for parent folder, containing all wells, e.g. your/path/to/parent/<well-id> or your/path/to/parent")
parser.add_argument("--output_dir", type=str,
                        help="provide path where processing stuffs will be dumped using identical directory tree")
parser.add_argument("--Beam", type=str, help="Beamline ID needs to be specified, eg. PXI or PXII\n")

parser.add_argument("--total_degree", type=str, help="provide angular range to process")

parser.add_argument("--SG_num", type=str, default="0",
                        help="optionally, Space-group number can be specified, default is 0")
parser.add_argument("--cell", type=str, default="70 70 30 90 90 90",
                        help="optionally, unit-cell can be specified as 70 70 30 90 90 90; otherwise it will try to determine by itself")

parser.add_argument("--highres", type=str, default="2.5", 
                         help="optionally high-resolution limit can be given, default: 2.5")
parser.add_argument("--friedel", type=str, default="FALSE", help="optionally, it can be changed to true..")

parser.add_argument("--refs", type=str, help='optionally, reference data set for indexing can be provided..')

parser.add_argument("--strong_pixel", type=str)
parser.add_argument("--min_pix_spot", type=str)
parser.add_argument("--merge_paths", type=str, nargs='+')
parser.add_argument("--ISa_cutoff", type=str)

args = parser.parse_args()
'''
============================================
'''
class Xtal(object):
    '''class Xtal defines each crystal as object to run xds job, also it prepares XDS.INP files for each detector type.
       kwargs option is used for optional parameters, e.g., space group, reference file, resolution cut off so on.
    '''
    def __init__(self,xtalImgPath,xtalProcessPath,xtalNum, BL, tot_angle=None, **kwargs):
        
        
        self.xtalimgpath = xtalImgPath
        self.xtalprocesspath = xtalProcessPath
        self.xtalnum = xtalNum
        self.xtalname = None
        self.wavelength = None
        self.det_dis = None
        self.beamx = None
        self.beamy = None
        self.osci = None
        self.osci_range = tot_angle
        self.datrange_str = None
        self.bgrange_str = None
        self.sprange_str = None
        self.beamline = BL
        self.SG = kwargs.get('SG', '0')
        self.cell = kwargs.get('cell', "70 70 30 90 90 90")
        self.res_cut = kwargs.get('res_cut', "2.5")
        self.idx_res = kwargs.get('idx_res', "5.0")
        self.friedel = kwargs.get('friedel', "FALSE")
        self.refdata = kwargs.get('refdata', " ")
        self.strong_pixel = kwargs.get('strong_pixel', '6.0')
        self.min_pix_spot = kwargs.get('min_pix_spot', '3')
        self.content = {}
	
    def read_masterh5(self, master_file):
    
        #read master file headers for xds.inp preparation
        header = h5py.File(master_file, 'r')
	#beam center x and y
        beamx = header['/entry/instrument/detector/beam_center_x']
        beamx = np.array(beamx)
        self.beamx = beamx
        beamy = header['/entry/instrument/detector/beam_center_y']
        beamy = np.array(beamy)
        self.beamy = beamy 
	
	#wavelength and detector distance
        wave = header['/entry/instrument/beam/incident_wavelength']
        wave = np.array(wave)
        self.wavelength = wave
        detZ = header['/entry/instrument/detector/detector_distance']
        detZ = np.array(detZ)
        self.det_dis = round(detZ*1e3) # convert distance into millimeter
	
        #omega and oscillation
        omega = header['/entry/sample/goniometer/omega_increment']
        omega = np.array(omega)
        self.osci = omega

	#number of images
        if self.osci_range is None:
           nimages = header['/entry/instrument/detector/detectorSpecific/nimages']
           nimages = np.array(nimages)
	else:
           nimages = round(float(self.osci_range)/self.osci)
           nimages = int(nimages)
        data_start = self.xtalnum*nimages + 1; data_end = (self.xtalnum+1)*nimages
        self.datrange_str = (str(data_start), str(data_end))
        self.sprange_str = self.datrange_str
	self.bgrange_str = (str(data_start), str(data_start+10))
	
	
	# xtal filename template with path   
        name = os.path.basename(master_file)
        name = name.split( 'master' )
	#self.xtalname = self.xtalimgpath + name[0]+"??????.h5"
        img = name[0]+"??????.h5"
        self.xtalname = os.path.join(self.xtalimgpath, img)

    def read_cbf(self, headerfile):
        #read cbf file header and store the info in dictionary and later prepare xds.inp

        cmd = "head -35 "+headerfile+" > filehead.txt"
        os.system(cmd)
        self.xtalname = headerfile[:-9]+"?????.cbf"
        keys_head = ["Beam_xy","Wavelength", "Detector_distance","Angle_increment"]
        fh = open('filehead.txt', 'r')
        all_lines = fh.readlines()
        fh.close()

        for lines in all_lines:
            if any(key in lines for key in keys_head):
               line = lines.split()
               if line[1] == 'Beam_xy':
                  self.content['ORGX'] = str(line[2].strip('(').strip(','))
                  self.content['ORGY'] = str(line[3].strip(')'))
               elif line[1] == 'Wavelength':
                    self.content['X-RAY_WAVELENGTH'] = line[2]
               elif line[1] == 'Detector_distance':
                    self.content["DETECTOR_DISTANCE"] = str(float(line[2])*1e3)
               else:
                    self.content["OSCILLATION_RANGE"] = line[2]

                    self.datrange_end = int(int(self.osci_range)/float(line[2]))
                    self.datrange_str = (1, str(self.datrange_end))
                    self.sprange_str = self.datrange_str

        self.content["NAME_TEMPLATE_OF_DATA_FRAMES"] = self.xtalname
        self.content["DATA_RANGE"] = str(self.datrange_str[0])+"  "+str(self.datrange_str[1])
        self.content["SPOT_RANGE"] = str(self.sprange_str[0])+"  "+str(self.sprange_str[1])

        self.content["SPACE_GROUP_NUMBER"] = self.SG
        self.content["UNIT_CELL_CONSTANTS"] = self.cell
        self.content["INCLUDE_RESOLUTION_RANGE"] = str('50  '+ str(self.res_cut))
        self.content["FRIEDEL'S_LAW"] = self.friedel
        self.content["MINIMUM_NUMBER_OF_PIXELS_IN_A_SPOT"] = self.min_pix_spot
        self.content["STRONG_PIXEL"] = self.strong_pixel

        if self.refdata != " ":
           os.chdir(self.xtalprocesspath)
           ref_link = "reference.HKL"
           os.symlink(self.refdata, ref_link)
           self.content["REFERENCE_DATA_SET"] = ref_link
        else:
           self.content["REFERENCE_DATA_SET"] = self.refdata

        
 
    def locatextalpath(self):
        #Locate and sanity check if the data exists and then read the headers/master files

        if not os.path.exists(self.xtalimgpath):
	   print 'Error: path does not exist\n' 
	   sys.exit()
        if self.beamline == "PXI":
           master_file = glob.glob(os.path.join(self.xtalimgpath,"*_master.h5"))[0]
           try:
             self.read_masterh5(master_file)
           except OSError:
               print "master file may not exist\n"

        elif self.beamline == "PXII":
           cbf_header = glob.glob(os.path.join(self.xtalimgpath,"*_00001.cbf"))[0]
           try:
             self.read_cbf(cbf_header)
           except OSError:
               print "cbf may not be collected yet\n"

    def create_eigerInp(self):
        #create xds.inp for Eiger data..

        try:
          os.chdir(self.xtalprocesspath)
        except OSError:
            print "xtal process folder may not be created yet\n"

        if not os.path.isfile("XDS.INP"):
           eiger = open("XDS.INP", 'w')
           eiger.write("! Settings for EIGER @ SLS generated by IMISX automation \n")
           eiger.write("!====================== JOB CONTROL PARAMETERS ===============================\n")
           eiger.write("JOB=XYCORR INIT COLSPOT IDXREF !DEFPIX INTEGRATE CORRECT\n")
           eiger.write("!JOB=DEFPIX INTEGRATE CORRECT\n")
           eiger.write("MAXIMUM_NUMBER_OF_JOBS=1\n")
           eiger.write("MAXIMUM_NUMBER_OF_PROCESSORS=12\n")
           eiger.write("\n\n")
           eiger.write("! for this experiment:\n")

           for k, v in sorted(self.content.items()):
	       eiger.write("%s= %s\n" %(k,v))
           
           eiger.write("!REIDX=   0  0 -1  0  0 -1  0  0 -1  0  0  0\n")
           eiger.write("REFINE(IDXREF)=BEAM AXIS ORIENTATION CELL !POSITION\n")
           eiger.write("REFINE(INTEGRATE)= BEAM ORIENTATION CELL !AXIS POSITION\n")
           eiger.write("REFINE(CORRECT)= BEAM ORIENTATION CELL AXIS !POSITION\n")
           eiger.write("STRICT_ABSORPTION_CORRECTION=FALSE\n")
           eiger.write("TEST_RESOLUTION_RANGE=10.0 5.0\n")
           eiger.write("TRUSTED_REGION=0.00 1.2\n")
           eiger.write("VALUE_RANGE_FOR_TRUSTED_DETECTOR_PIXELS=4000. 30000.\n")
           eiger.write("CORRECTIONS=DECAY MODULATION ABSORP\n")
           
           eiger.write("BACKGROUND_PIXEL=6.0\n")
           
           eiger.write("MINIMUM_FRACTION_OF_INDEXED_SPOTS=0.10\n")
           eiger.write("SEPMIN=4\n")
           eiger.write("CLUSTER_RADIUS=2\n")
           eiger.write("NUMBER_OF_PROFILE_GRID_POINTS_ALONG_ALPHA/BETA=21\n\n")
           eiger.write("! parameters specifically for this detector and beamline:\n\n")
           eiger.write("DETECTOR=EIGER NX=4150 NY=4371 QX=0.075  QY=0.075\n")
           eiger.write("MINIMUM_VALID_PIXEL_VALUE=0  OVERLOAD=3000000\n")
           eiger.write("DIRECTION_OF_DETECTOR_X-AXIS=1 0 0\n")
           eiger.write("DIRECTION_OF_DETECTOR_Y-AXIS=0 1 0\n")
           eiger.write("INCIDENT_BEAM_DIRECTION=0 0 1\n")
           eiger.write("ROTATION_AXIS=1 0 0\n")
           eiger.write("FRACTION_OF_POLARIZATION=0.99\n")
           eiger.write("POLARIZATION_PLANE_NORMAL=0 1 0\n")
           eiger.write("SENSOR_THICKNESS=0.45\n")
           eiger.write("!EXCLUSION OF HORIZONTAL DEAD AREAS OF THE EIGER 16M DETECTOR + ONE PIXEL ON EACH SIDE\n")
           eiger.write("UNTRUSTED_RECTANGLE=    0 4151    513  553\n")
           eiger.write("UNTRUSTED_RECTANGLE=    0 4151   1064 1104\n")
           eiger.write("UNTRUSTED_RECTANGLE=    0 4151   1615 1655\n")
           eiger.write("UNTRUSTED_RECTANGLE=    0 4151   2166 2206\n")
           eiger.write("UNTRUSTED_RECTANGLE=    0 4151   2717 2757\n")
           eiger.write("UNTRUSTED_RECTANGLE=    0 4151   3268 3308\n")
           eiger.write("UNTRUSTED_RECTANGLE=    0 4151   3819 3859\n")
           eiger.write("!EXCLUSION OF VERTICAL DEAD AREAS OF THE EIGER 16M DETECTOR + ONE PIXEL ON EACH SIDE\n")
           eiger.write("UNTRUSTED_RECTANGLE= 1029 1042      0 4372\n")
           eiger.write("UNTRUSTED_RECTANGLE= 2069 2082      0 4372\n")
           eiger.write("UNTRUSTED_RECTANGLE= 3109 3122      0 4372\n\n")
        
           eiger.close()
        else:
            pass

    def mod_eigerInp(self):
        os.system("cat 'XDS.INP' | sed 's/!JOB/JOB/g' > mod.inp")
        os.system("cat mod.inp | sed 's/INCLUDE_RESOLUTION_RANGE= 50  5.0/INCLUDE_RESOLUTION_RANGE=50  "+self.res_cut+"/g' > XDS.INP")
        #nasty patch to re-run xds with integrate and correct steps with full resolution. Need to be organized better
       
    def create_pilatusInp(self):
        #create xds.inp for Pilatus 6M data..
        try:
          os.chdir(self.xtalprocesspath)
        except OSError:
          print "xtal process folder may not be created yet\n"

        if not os.path.isfile("XDS.INP"):
           pilatus = open("XDS.INP", 'w')
           pilatus.write("! Settings for PILATUS @ SLS generated by IMISX automation \n")
           pilatus.write("!====================== JOB CONTROL PARAMETERS ===============================\n")
           pilatus.write("JOB=XYCORR INIT COLSPOT IDXREF DEFPIX INTEGRATE CORRECT\n")
           pilatus.write("MAXIMUM_NUMBER_OF_JOBS=4\n")
           pilatus.write("MAXIMUM_NUMBER_OF_PROCESSORS=1\n")
           pilatus.write("\n\n")
           pilatus.write("! for this experiment:\n")

           for k, v in sorted(self.content.items()):
	       pilatus.write("%s= %s\n" %(k,v))
           
           pilatus.write("!REIDX=   0  0 -1  0  0 -1  0  0 -1  0  0  0\n")
           pilatus.write("REFINE(IDXREF)=BEAM AXIS ORIENTATION CELL !POSITION\n")
           pilatus.write("REFINE(INTEGRATE)= BEAM ORIENTATION CELL !AXIS POSITION\n")
           pilatus.write("REFINE(CORRECT)= BEAM ORIENTATION CELL AXIS !POSITION\n")
           pilatus.write("STRICT_ABSORPTION_CORRECTION=FALSE\n")
           pilatus.write("TRUSTED_REGION=0.00 1.15\n")
           pilatus.write("VALUE_RANGE_FOR_TRUSTED_DETECTOR_PIXELS=8000. 30000.\n")
           pilatus.write("CORRECTIONS=DECAY MODULATION ABSORP\n")
           
           #pilatus.write("BACKGROUND_RANGE=2 10\n")
           
           pilatus.write("MINIMUM_FRACTION_OF_INDEXED_SPOTS=0.20\n")
           pilatus.write("NUMBER_OF_PROFILE_GRID_POINTS_ALONG_ALPHA/BETA=13\n\n")
           pilatus.write("! parameters specifically for this detector and beamline:\n\n")
           pilatus.write("DETECTOR=PILATUS NX=2463 NY=2527 QX=0.172  QY=0.172 !PILATUS 6M\n")
           pilatus.write("MINIMUM_VALID_PIXEL_VALUE=0  OVERLOAD=1048500\n")
           pilatus.write("DIRECTION_OF_DETECTOR_X-AXIS=1 0 0\n")
           pilatus.write("DIRECTION_OF_DETECTOR_Y-AXIS=0 1 0\n")
           pilatus.write("INCIDENT_BEAM_DIRECTION=0 0 1\n")
           pilatus.write("ROTATION_AXIS=1 0 0\n")
           pilatus.write("FRACTION_OF_POLARIZATION=0.99\n")
           pilatus.write("POLARIZATION_PLANE_NORMAL=0 1 0\n")
           pilatus.write("SENSOR_THICKNESS=0.32\n")
           pilatus.write("UNTRUSTED_RECTANGLE= 487  495     1 2527\n")
           pilatus.write("UNTRUSTED_RECTANGLE= 981  989     1 2527\n")
           pilatus.write("UNTRUSTED_RECTANGLE=1475 1483     1 2527\n")
           pilatus.write("UNTRUSTED_RECTANGLE=1969 1977     1 2527\n")
           pilatus.write("UNTRUSTED_RECTANGLE=   1 2463   195  213\n")
           pilatus.write("UNTRUSTED_RECTANGLE=   1 2463   407  425\n")
           pilatus.write("UNTRUSTED_RECTANGLE=   1 2463   619  637\n")
           pilatus.write("UNTRUSTED_RECTANGLE=   1 2463   831  849\n")
           pilatus.write("UNTRUSTED_RECTANGLE=   1 2463  1043 1061\n")
           pilatus.write("UNTRUSTED_RECTANGLE=   1 2463  1255 1273\n")
           pilatus.write("UNTRUSTED_RECTANGLE=   1 2463  1467 1485\n")
           pilatus.write("UNTRUSTED_RECTANGLE=   1 2463  1679 1697\n")
           pilatus.write("UNTRUSTED_RECTANGLE=   1 2463  1891 1909\n")
           pilatus.write("UNTRUSTED_RECTANGLE=   1 2463  2103 2121\n")
           pilatus.write("UNTRUSTED_RECTANGLE=   1 2463  2315 2333\n")
           
           pilatus.close()
        else:
           pass

    def update_eiger(self):
    
    #this method is to update/create xds.inp file based on each xtal dataset info..
    
     
	self.locatextalpath()
	
	self.content["NAME_TEMPLATE_OF_DATA_FRAMES"] = self.xtalname
	self.content["OSCILLATION_RANGE"] = str(self.osci)
	self.content["DATA_RANGE"] = str(self.datrange_str[0])+"  "+str(self.datrange_str[1])
	self.content["SPOT_RANGE"] = str(self.sprange_str[0])+"  "+str(self.sprange_str[1])
        #self.content["BACKGROUND_RANGE"] = str(self.bgrange_str[0])+"  "+str(self.bgrange_str[1])
	self.content["X-RAY_WAVELENGTH"] = str(self.wavelength)
	self.content["ORGX"] = str(self.beamx)
	self.content["ORGY"] = str(self.beamy)
	self.content["DETECTOR_DISTANCE"] = str(self.det_dis)
        self.content["SPACE_GROUP_NUMBER"] = self.SG
        self.content["UNIT_CELL_CONSTANTS"] = self.cell
        self.content["INCLUDE_RESOLUTION_RANGE"] = str('50  '+ str(self.res_cut))
        self.content["FRIEDEL'S_LAW"] = self.friedel
        self.content["MINIMUM_NUMBER_OF_PIXELS_IN_A_SPOT"] = self.min_pix_spot
        self.content["STRONG_PIXEL"] = self.strong_pixel

        if not self.refdata is " ":
           os.chdir(self.xtalprocesspath)
           ref_link = "reference.HKL"
           os.symlink(self.refdata, ref_link)
           self.content["REFERENCE_DATA_SET"] = ref_link
        else:
           self.content["REFERENCE_DATA_SET"] = " ";

        #self.create_eigerInp()
        self.mod_eigerInp()
        
    def idx_eiger(self):
    
    #this method is to run xds up to IDXREF step with low-res cut-off of 5A.
    
     
	self.locatextalpath()
	
	self.content["NAME_TEMPLATE_OF_DATA_FRAMES"] = self.xtalname
	self.content["OSCILLATION_RANGE"] = str(self.osci)
	self.content["DATA_RANGE"] = str(self.datrange_str[0])+"  "+str(self.datrange_str[1])
	self.content["SPOT_RANGE"] = str(self.sprange_str[0])+"  "+str(self.sprange_str[1])
        #self.content["BACKGROUND_RANGE"] = str(self.bgrange_str[0])+"  "+str(self.bgrange_str[1])
	self.content["X-RAY_WAVELENGTH"] = str(self.wavelength)
	self.content["ORGX"] = str(self.beamx)
	self.content["ORGY"] = str(self.beamy)
	self.content["DETECTOR_DISTANCE"] = str(self.det_dis)
        self.content["SPACE_GROUP_NUMBER"] = self.SG
        self.content["UNIT_CELL_CONSTANTS"] = self.cell
        self.content["INCLUDE_RESOLUTION_RANGE"] = str('50  '+ str(self.idx_res))
        self.content["FRIEDEL'S_LAW"] = self.friedel
        self.content["MINIMUM_NUMBER_OF_PIXELS_IN_A_SPOT"] = self.min_pix_spot
        self.content["STRONG_PIXEL"] = self.strong_pixel

        if not self.refdata is " ":
           os.chdir(self.xtalprocesspath)
           ref_link = "reference.HKL"
           os.symlink(self.refdata, ref_link)
           self.content["REFERENCE_DATA_SET"] = ref_link
        else:
           self.content["REFERENCE_DATA_SET"] = " ";

        self.create_eigerInp()

    def update_pilatus(self):

        try:
          self.locatextalpath()
        except OSError:
            pass
 
        if not self.content:
            print "cbf header has not been read\n"
        else:
            self.create_pilatusInp()

        
    def runxds(self):
        #method to call create necessary files, locate data and finally run xds for each xtal.
        if self.beamline == "PXI":
           self.idx_eiger()
           os.system("xds_par >/dev/null")
           self.update_eiger()
 
        elif self.beamline == "PXII":
             os.chdir(self.xtalprocesspath)
        os.system("xds_par >/dev/null")
        if not os.path.isfile("XDS_ASCII.HKL"):
           print "xtal: %d is failed from folder %s\n" %(self.xtalnum, self.xtalprocesspath),
        else:
           print "xtal: %d is processed\n" %(self.xtalnum),

class process(object):
     '''
       Process Class handles folder structures, data finders, and launching multiple xds jobs for xtals. 
       It also provide methods for running XSCALE and ISa based selection.
       
     '''
     def __init__(self, data_dir, output, tot_angle=None):
         self.data_folder = data_dir
	 self.process_folder = None
	 self.output = os.path.join(output, "proc")
         if not os.path.exists(self.output):
            os.makedirs(self.output, 0755)
         else:
            pass 
	 self.nxtals = None
	 self.mininame = None
         self.dataName = None
         self.ISa_th = 4.0
         self.total_range = tot_angle
         self.xscale_file_list = [];
         self.xtals_lists = []
     
     def find_HKLs(self, **kwargs):
         #method to look up for XDS_ASCII.HKL files produced by xds and create a list for XSCALE.

         hklpath = kwargs.get('hklpath', None)

         try:
           os.chdir(self.output)
         except OSError:
             print "check if the output folder exists\n"
         if not kwargs:
           for parent, dirs, files in os.walk(self.output):
             for fh in files:
                 if fh == "XDS_ASCII.HKL":
                    HKLpath = os.path.join(parent,fh)
                    self.xscale_file_list.append(HKLpath)
                 else:
                    pass
         else:
            for parent, dirs, files in os.walk(hklpath):
             for fh in files:
                 if fh == "XDS_ASCII.HKL":
                    HKLpath = os.path.join(parent,fh)
                    self.xscale_file_list.append(HKLpath)
                 else:
                    pass
 
     def read_xscale(self):
        #method to read xscale.LP file and print the statistics for user..

        if os.path.isfile("XSCALE.LP"):
           fh = open("XSCALE.LP", 'r')
           all_lines = fh.readlines()
           fh.close()
           try:
             start = all_lines.index(' SUBSET OF INTENSITY DATA WITH SIGNAL/NOISE >= -3.0 AS FUNCTION OF RESOLUTION\n')
             end = all_lines.index(' ========== STATISTICS OF INPUT DATA SET ==========\n')
             useful = all_lines[start+1:end-1]
           except ValueError, IndexError:
               print "there was no XDS_ASCII files, so nothing to XSCALE. Check input/output directory naming\n"
           for lines in useful:
               print "%s" %lines,

     def Isa_select(self, ISa_th):
        #method to perform ISa based selection of 'good' HKL files for next round of XSCALEing
        if os.path.isfile("XSCALE.LP"):
         fh = open("XSCALE.LP", 'r')
         all_lines = fh.readlines()
         fh.close()
         #ISa_th = kwargs.get('ISa_th', "4.0")
         self.selected_files = [];
         try:
           start = all_lines.index('     a        b          ISa    ISa0   INPUT DATA SET\n')
         # select the table with ISa values from the file..
           Isa_list = all_lines[start+1:start+len(self.xscale_file_list)]
         except ValueError, IndexError:
              print "check if XSCALE ran properly or XSCALE.INP screwed up \n"

         if len(Isa_list) > 0:
          for lines in Isa_list:
              line = lines.split()
              try:
                if float(line[2]) > float(ISa_th):
                   self.selected_files.append(line[4])
              except IndexError:
                   pass
         

     def run_xscale(self, **kwargs):
        #method to run xscale, it calls Isa_select and find_HKLs methods internally and run two rounds of XSCALing
        
        merge_paths = kwargs.get('merge_paths', [])
        isa_cut = kwargs.get('isa_cut', '4.0')
        
        if len(merge_paths) > 0:
           for path in merge_paths:
               self.find_HKLs(path=kwargs['hklpath'])
        else:
           self.find_HKLs()
        
        if len(self.xscale_file_list) > 0:
           try:
             os.makedirs("multi", 0755)
           except OSError:
               #print "either there's no HKL files or no permission to create multi folder\n"
               pass

           os.chdir("multi")

           fh = open("XSCALE.INP",'w')
           fh.write("OUTPUT_FILE=XSCALE.HKL\n")
           fh.write("SAVE_CORRECTION_IMAGES=FALSE\n")
           fh.write("FRIEDEL'S_LAW=TRUE\n")
           fh.write("\n\n")
           for ii in range (len(self.xscale_file_list)):
             try:
               linkname = "xtal_"+str(ii)+".HKL"
               os.symlink(self.xscale_file_list[ii], linkname)
             except OSError, e:
                  if e.errno == errno.EEXIST:
                     os.remove(linkname)
                     os.symlink(self.xscale_file_list[ii], linkname)
                  else:
                     raise e

             fh.write("INPUT_FILE="+linkname+"\n")
             fh.write("MINIMUM_I/SIGMA=0\n")

           fh.close()
        print "Running 1st round of xscale\n"
        os.system("xscale_par >/dev/null")

        self.read_xscale()
        
        self.Isa_select(isa_cut)
        os.system("cp XSCALE.LP nofilter.LP")
        os.system("cp XSCALE.HKL nofilter.HKL")

        if len(self.selected_files) > 1:
           fh = open("XSCALE.INP",'w')
           fh.write("OUTPUT_FILE=XSCALE.HKL\n")
           fh.write("SAVE_CORRECTION_IMAGES=FALSE\n")
           fh.write("FRIEDEL'S_LAW=TRUE\n")
           fh.write("\n\n")
           for files in self.selected_files:
               fh.write("INPUT_FILE="+files+"\n")
               fh.write("MINIMUM_I/SIGMA=0\n")
           fh.close()

        print "running xscale after ISa selection\n\n"

        os.system("xscale_par > /dev/null")

        print "%d out of %d crystals selected\n" %(len(self.selected_files), len(self.xscale_file_list))

        self.read_xscale()

class Eigerh5(process):
    #child class of Process for Eiger data; folder structure and parallization protocol are bit different from Pilatus.

    def __init__(self,data_dir, output, tot_angle=None):
        process.__init__(self,data_dir, output, tot_angle=None)
        self.miniset_folders = [];
        self.beamline = "PXI"
    
    def Lookup(self):
        for ii in range(len(self.data_folder)):
            for dirs in sorted(glob.glob(os.path.join(self.data_folder[ii], "minisets*"))):
                self.miniset_folders.append(dirs)

    def get_xtals_eiger(self, **kwargs):

        for ii in range(len(self.data_folder)):
            if self.data_folder[ii].endswith('*'):
               parent_dir = self.data_folder[ii][:-1]
               if not os.path.exists(parent_dir):
                  print "Error: data directory does not exist!\n"
	          sys.exit()
            else:
               if not os.path.exists(self.data_folder[ii]):
	          print "Error: data directory does not exist!\n"
	          sys.exit()
        
        try:
           os.chdir(self.output)
        except OSError:
            print "output folder may not exist. check!\n"
            sys.exit()

        self.Lookup()
        if len(self.miniset_folders) > 0:
          for i in range(len(self.miniset_folders)):
            self.xtal_each_miniset = sorted(glob.glob(os.path.join(self.miniset_folders[i], "*data*.h5")))
            try:
              self.xtal_each_miniset.pop()
            except IndexError:
               print "all data is not saved yet, will process in next round\n"
               pass
            self.nxtals = len(self.xtal_each_miniset)
	    print "%d xtals in %s miniset \n" %(self.nxtals, self.miniset_folders[i]),
            self.mininame = os.path.basename(self.miniset_folders[i])
            dir_von_miniset = os.path.dirname(self.miniset_folders[i])
            self.dataName = os.path.basename(dir_von_miniset)
	    self.process_folder = os.path.join(self.output, self.dataName,self.mininame)

            if not os.path.exists(self.process_folder):
                print "creating processing directory %s\n" %(self.process_folder),
	        os.makedirs(self.process_folder, 0755)
            else:
                pass

            os.chdir(self.process_folder)
            for k in range(self.nxtals):
                xtal_process_path = self.process_folder +'/xtal_'+str(k)
                if not os.path.exists(xtal_process_path):
                   os.mkdir(xtal_process_path, 0755)
                   os.chdir(xtal_process_path)
                   if not kwargs:
                      xtalobj = Xtal(self.miniset_folders[i],xtal_process_path,k, self.beamline, self.total_range)
                   else:
                      xtalobj = Xtal(self.miniset_folders[i],xtal_process_path,k, self.beamline, self.total_range, SG=kwargs['SG'], cell=kwargs['cell'], res_cut=kwargs['res_cut'], friedel=kwargs['friedel'], refdata=kwargs['refdata'], strong_pixel=kwargs['strong_pixel'], min_pix_spot=kwargs['min_pix_spot'])
                   #xtalobj.update_eiger()
                   self.xtals_lists.append(xtalobj)
                    
                else:
                   print "folder may exist, skipping it %s\n" %xtal_process_path,  
                

            os.chdir(self.output)

          print "\n I gathered %d xtals from data directory\n\n" %(len(self.xtals_lists))


    def run_eiger(self):
          
       job_cnt = 0
       if len(self.xtals_lists) > 0:
        for j in range(len(self.xtals_lists)): 
       
          proc = [];
          for i in range (0,4):
            try:
              jobid = mp.Process(target=self.xtals_lists[(j*4)+i].runxds)
              proc.append(jobid)
                    
            except IndexError:
                pass
          for p in proc:
              p.start()
            
          for p in proc:
              p.join()
             
          job_cnt += 4
          print "%d crystals has been processed\n" %job_cnt

class PilatusCBF(process):
    #child class of Process for Pilatus data; folder structure and parallization protocol are bit different from Eiger.

    def __init__(self, data_dir, output, tot_angle="10"):
        process.__init__(self, data_dir, output, tot_angle="10")
        self.beamline = "PXII"
        self.minixtal_folders = [];
        self.total_range = tot_angle;
    
    def Lookup(self):
        for ii in range(len(self.data_folder)):
            for xtal_dirs in sorted(glob.glob(os.path.join(self.data_folder[ii], "minisets/xtal*"))):
                self.minixtal_folders.append(xtal_dirs)
        print "%d xtals found so far\n" %(len(self.minixtal_folders))
    
    
    def get_xtals_cbf(self, **kwargs):

        for ii in range(len(self.data_folder)):
            if self.data_folder[ii].endswith('*'):
               parent_dir = self.data_folder[ii][:-1]
               if not os.path.exists(parent_dir):
                  print "Error: data directory does not exist!\n"
	          sys.exit()
            else:
               if not os.path.exists(self.data_folder[ii]):
	          print "Error: data directory does not exist!\n"
	          sys.exit()

        try:
           os.chdir(self.output)
        except OSError:
            print "output folder may not exist. check!\n"
            sys.exit()

        self.Lookup()
        if len(self.minixtal_folders) > 0:
           for i in range(len(self.minixtal_folders)):
               self.mininame = os.path.basename(self.minixtal_folders[i])
               dir_von_miniset = os.path.dirname(self.minixtal_folders[i])
               dir_well = os.path.dirname(dir_von_miniset)
               self.dataName = os.path.basename(dir_well)
               dir_puck = os.path.dirname(dir_well)
               puckname = os.path.basename(dir_puck)
               self.process_folder = os.path.join(self.output, puckname, self.dataName, self.mininame)
 
               try:
                  os.makedirs(self.process_folder, 0755)
                  print "creating process directory %s\n" %self.process_folder,
                  os.chdir(self.process_folder)
               
               
                  if not kwargs:
                     xtalobj = Xtal(self.minixtal_folders[i], self.process_folder, i, self.beamline, self.total_range)
                  else:
                     xtalobj = Xtal(self.minixtal_folders[i],self.process_folder, i, self.beamline, self.total_range, SG=kwargs['SG'], cell=kwargs['cell'], res_cut=kwargs['res_cut'], friedel=kwargs['friedel'], refdata=kwargs['refdata'], strong_pixel=kwargs['strong_pixel'], min_pix_spot=kwargs['min_pix_spot'])

                  xtalobj.update_pilatus() #run update method first before parallization to avoid racing symtoms 
                  self.xtals_lists.append(xtalobj)
               except OSError as e:
                  print "folder may exist %s" %self.process_folder
                  if e.errno == errno.EEXIST:
                     pass
                  else:
                     raise e

           os.chdir(self.output)
           print "\nI gathered %d xtals in total\n\n" %(len(self.xtals_lists))

    def run_pilatus(self):
       job_cnt = 0
       if len(self.xtals_lists) > 0:
         for j in range(len(self.xtals_lists)): 
         
           proc = [];
           for i in range (0,10):
            
            try:
              jobid = mp.Process(target=self.xtals_lists[(j*10)+i].runxds)
              proc.append(jobid)
                    
            except IndexError:
                pass
           for p in proc:
              p.start()
            
           for p in proc:
              p.join()
             
           job_cnt += 10
           print "%d crystals has been processed\n" %job_cnt
       else:
          pass   

'''
=============================Main code to control and call classes==================
'''                           
def main():

   if args.well_path is None:
      sys.exit("tell me where the data is located \n")
   elif len(args.well_path) > 0:
      data_path = args.well_path  # It takes a list type of data folders, you can even use bash wild-card
       
   if args.output_dir is None:
      print "output directory path is not provided \n"
      cwd = os.getcwd()
      
      print "Dumping output at the current working directory %s \n" %cwd
      args.output_dir = cwd

   if args.total_degree is None:
      print "better if total angle range is specified\n"
   else:
      osci_range = args.total_degree


   if args.Beam == "PXI":
      proc_well = Eigerh5(data_path, args.output_dir, osci_range)
      if any (c is None for c in (args.SG_num, args.cell, args.highres, args.friedel, args.refs, args.strong_pixel, args.min_pix_spot)):
           
          proc_well.get_xtals_eiger()
      else:
          proc_well.get_xtals_eiger(SG=args.SG_num, cell=args.cell, res_cut=args.highres, friedel=args.friedel, refdata=args.refs, strong_pixel=args.strong_pixel, min_pix_spot=args.min_pix_spot)
      
      proc_well.run_eiger()
   elif args.Beam == "PXII":
      proc_well = PilatusCBF(data_path, args.output_dir, osci_range)

      if any (c is None for c in (args.SG_num, args.cell, args.highres, args.friedel, args.refs, args.strong_pixel, args.min_pix_spot)):
           
          proc_well.get_xtals_cbf()
      else:
         proc_well.get_xtals_cbf(SG=args.SG_num, cell=args.cell, res_cut=args.highres, friedel=args.friedel, refdata=args.refs, strong_pixel=args.strong_pixel, min_pix_spot=args.min_pix_spot)

      proc_well.run_pilatus()

   if args.merge_paths is None:
      proc_well.run_xscale(isa_cut=args.ISa_cutoff)
   else:
      proc_well.run_xscale(merge_paths=args.merge_paths)
 
if __name__ == '__main__':
   main()




