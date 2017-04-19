# -*- coding: utf-8 -*-
"""
Created on Tue Aug 5 14:52:14 2016

@author: shibom
"""

import sys, os, time
import PyQt4.QtGui as Qt
import PyQt4.QtCore as QC
import subprocess as sub
import plot_utils as psx


def writer(filename, symm, sites, resolution, emins, ntries):
    ofh = open('run_cluster','w')
    ofh.write("#!/bin/bash \n")
    ofh.write('\n\n')

    ofh.write('''\
python shelx_batch.py --hklfile %s --symm %s\
 --sites %s --resolution %s --emins %s --ntries %s\
''' %(filename,symm,sites,resolution,emins,ntries))

    ofh.close()
    sub.call(['chmod +x run_cluster'],shell=True)

def mail_write(email_id,output):
    fh = open("tmp_mail.sh", 'w')
    fh.write("#!/bin/bash\n")
    fh.write('\n\n')
    fh.write('sleep 2\n')
    fh.write('python /mnt/das-gpfs/work/p11206/shibom/20S/scripts/plot_cluster.py %s\n' %output)
    fh.write('sleep 2\n')
    fh.write('cd %s\n' %output)
    fh.write('echo "All jobs finished, check output" | mail -s "shelx jobs status" -a CCplot.pdf %s\n' %email_id)
    fh.write('sleep 2\n')
    fh.close()

    sub.call(['chmod +x tmp_mail.sh'], shell=True)


class Tabs(Qt.QTabWidget):
    def __init__(self):
        Qt.QTabWidget.__init__(self)
        self.tab1 = Scale_Merge()
        self.tab2 = Substructure()
        self.tab3 = Qt.QMainWindow()

        self.addTab(self.tab1, "scale_merge")
        self.addTab(self.tab2, "Substructure")
        self.addTab(self.tab3, "Model_build")
        self.setWindowTitle("Native-SAD ui")

class Scale_Merge(Qt.QMainWindow):
    def __init__(self):
        Qt.QMainWindow.__init__(self)
        self.layout = Qt.QVBoxLayout()
        self.layout.addWidget(MainLayout(adding_files()))
        self.layout.addWidget(MainLayout(scale_merge()))
        self.layout.addWidget(MainLayout(quick_anom()))
        self.centralWidget = Qt.QWidget()
        self.centralWidget.setLayout(self.layout)

        self.setCentralWidget(self.centralWidget)

class Substructure(Qt.QMainWindow):
    def __init__(self):
        Qt.QMainWindow.__init__(self)
        self.layout = Qt.QVBoxLayout()
        self.layout.addWidget(MainLayout(shelx_ui()))
        self.layout.addWidget(MainLayout(plotting_ui()))
        self.centralWidget = Qt.QWidget()
        self.centralWidget.setLayout(self.layout)

        self.setCentralWidget(self.centralWidget)

class MainLayout(Qt.QMainWindow):
    def __init__(self, classname):
        Qt.QMainWindow.__init__(self)


        self.scrollLayout = Qt.QFormLayout()
        self.scrollwidget = Qt.QWidget()
        self.scrollwidget.setLayout(self.scrollLayout)

        self.scrollArea = Qt.QScrollArea()
        self.scrollArea.setWidgetResizable(True)
        self.scrollArea.setWidget(self.scrollwidget)

        # main layout
        self.mainLayout = Qt.QVBoxLayout()

        # add all main to the main vLayout

        self.mainLayout.addWidget(self.scrollArea)
        self.add_widget(classname)
        #self.scrollLayout.addRow(plotting_ui())

        # central widget
        self.centralWidget = Qt.QWidget()

        self.centralWidget.setLayout(self.mainLayout)

        # set central widget
        self.setCentralWidget(self.centralWidget)

    def add_widget(self, classname):
        self.splitter = Qt.QSplitter()  #Nice bordering provides..
        self.splitter.addWidget(classname)
        self.scrollLayout.addRow(self.splitter)


    @classmethod
    def closing(MainLayout):
        MainLayout.close()
        print "App is closing"

class scale_merge(Qt.QWidget):
    HKLs = [];
    def __init__(self):
        Qt.QWidget.__init__(self)
        self.layout = Qt.QHBoxLayout()
        
        self.morefile_btn = widgets.createButton("Add more", self.add_filewidget)
        self.morefile_btn.setToolTip("Click to add more HKL files for merging")
        self.xscale_btn = widgets.createButton("xscale", self.merging)
        self.layout.addWidget(self.morefile_btn)
        self.layout.addWidget(self.xscale_btn)
        self.setLayout(self.layout)

    def add_filewidget(self):
        adding_files.File_layer.addWidget(files_widget("Unmerged HKL"))
        
        

    def merging(self):
        print "will merge with xscale"

class adding_files(Qt.QWidget):
    File_layer = Qt.QVBoxLayout()
    
    def __init__(self):
        Qt.QWidget.__init__(self)
        
        self.setWindowTitle("Add your HKL files")
        self.unmerged_data = files_widget("Unmerged HKL")
        
        self.__class__.File_layer.addWidget(self.unmerged_data)
        self.setLayout(self.__class__.File_layer)

class quick_anom(Qt.QWidget):
    def __init__(self):
        Qt.QWidget.__init__(self)
        layout = Qt.QVBoxLayout()
        self.get_file = files_widget("Reflection File")
        self.get_file.setToolTip("Reflection file must be SHELX supported format")
        self.sg = widgets("Space-Group")
        self.sg.setToolTip("Please use name, not space-group number.")
        self.xc_xcp = two_buttons("SHELXC", "SHELXC_plot")
        self.xc_xcp.Button1.clicked.connect(lambda: self.shelxc_run())
        self.xc_xcp.Button2.clicked.connect(lambda: self.shelxc_plot())
        layout.addWidget(self.get_file)
        layout.addWidget(self.sg)
        layout.addWidget(self.xc_xcp)
        self.setLayout(layout)

    def shelxc_run(self):
        inp = open("quicky.inp", 'w')
        if self.get_file.fname != None:
           inp.write("SAD "+self.get_file.fname+ '\n')
           cell = sub.check_output(["grep CELL "+self.get_file.fname + " | awk '{print $2,$3,$4,$5,$6,$7}'" ], shell=True)
           print cell
           inp.write("CELL "+ cell + '\n')
        else: 
           print "provide filename\n"
        if self.sg.value != None:
           inp.write("SPAG "+ self.sg.value + '\n')
           inp.write("MAXM 1000 \n")
        else:
           print "provide space group name, not number\n"
        inp.close()

        sub.call(["shelxc quick < quicky.inp" +' | tee quick.log'], shell=True)

    def shelxc_plot(self):
        try:
           self.ap = Qt.QDialog()
           self.ap.ui = psx.Main('quick.log', 'shelxc')
           self.ap.ui.show()
        except NameError:
             pass

class shelx_ui(Qt.QWidget):

    def __init__(self):
        Qt.QWidget.__init__(self)
        self.setWindowTitle("SHELX-GUI")
        # shelx_ui class should have layouts in vertical mode..
        layout = Qt.QVBoxLayout()
        #create instance for file loading..
        self.get_file = files_widget("Reflection-File")
        self.get_file.setToolTip("Reflection file must be SHELX supported format")
        #connect to enter key to make loaded file effective..
        self.get_file.file_box.returnPressed.connect(lambda: self.show_cell())

        #create Instance for output directory..
        self.outdir = dir_widget("Output-dir")
        self.outdir.setToolTip("A 'shelx_grid' directory with all files will be created under this output directory")

        #add file and output directory widgets vertically to the layout..
        layout.addWidget(self.get_file)
        layout.addWidget(self.outdir)
        #create Instances of widgets class for each different parameter..
        self.cell = widgets("Unit-cell")
        self.cell.setToolTip("cell value will pop up automatically, don't need to do anything")
        self.sg = widgets("Space-Group")
        self.sg.setToolTip("Please use name, not space-group number.")
        self.resol = widgets("Resolution")
        self.resol.setToolTip("provide values as comma separated, same for emins, sites")

        self.emins = widgets("Emins-values")
        self.sites = widgets("# Sites")
        self.tries = widgets("# trials")
        self.tries.setToolTip("please provide at least two numbers as comma separated")
        self.emails = widgets("email_id")
        self.emails.setToolTip("provide email-id within quote-sign '', job status will be reported once finished. ")

        #create submit-cancel buttons using Instance of two_buttons class..

        self.sub_abort = two_buttons("Submit", "Cancel")
        #add job submission functionality to "submit" and "cancel" buttons..
        self.sub_abort.Button1.clicked.connect(lambda: self.job_run())
        self.sub_abort.Button2.clicked.connect(lambda: self.closing())

        #add all other widgets to the vertical layout..

        layout.addWidget(self.cell)
        layout.addWidget(self.sg)
        layout.addWidget(self.resol)
        layout.addWidget(self.emins)
        layout.addWidget(self.sites)
        layout.addWidget(self.tries)
        layout.addWidget(self.emails)
        layout.addWidget(self.sub_abort)
        self.setLayout(layout)
        

    def show_cell(self):
        if self.get_file.fname != None:
            cell_value = sub.check_output(["grep CELL_CONSTANT "+self.get_file.fname + " | awk '{print $2,$3,$4,$5,$6,$7}'" ], shell=True)
            self.cell.textbox.setText(cell_value)
        else:
            msgbox = Qt.QMessageBox()
            msgbox.setText("Error: Load a reflection file first")
            msgbox.exec_()

    def job_run(self):
         #create an output directory for shelx jobs..

        if self.outdir.dir is None:
           self.path = "./shelx_grid"
        else:
           self.path = self.outdir.dir + "/shelx_grid";
           print self.path
        try:
          os.mkdir(self.path, 0755)
        except OSError:  #catching error related to overriding existing folder and impose overriding..
           pass

        print("we will dump all results in %s directory path" %self.path)

        os.chdir( self.path )
        sub.call(['cp /mnt/das-gpfs/work/p11206/shibom/20S/scripts/shelx_batch.py .'],shell=True)

        job_cnt = 0
        try:
           for site in self.sites.value:
               for res in self.resol.value:
                   for emin in self.emins.value:
                       for val in self.tries.value:
                           writer(self.get_file.fname,self.sg.value,site,res,emin,val)
                           sub.call(['sbatch -p day -J shelx run_cluster'], shell=True)
                           job_cnt += 1;
        except TypeError:
            pass
        print "%d Jobs have been submitted" %job_cnt
        self.keep_session()

        if self.emails.value is not None:
            self.tracker()

    def tracker(self):
        mail_write(self.emails.value,self.path)
        sub.call(['sbatch -d singleton -J shelx -t 300:00 tmp_mail.sh'], shell=True)

        '''
        job_status = 2;
        user = sub.check_output(['who am i'], shell=True)
        while job_status > 1:
            job_no = sub.check_output(['squeue -u '+user +' | wc -l'], shell=True)
            job_status = int(job_no)
            time.sleep(1800)

        if job_status == 1:
            send_email(self.emails.value)
        '''
    def keep_session(self):
        fh = open("session.param", 'w')
        
        fh.write("Reflection filename= %s\n" %self.get_file.fname)
        fh.write("Output directory= %s\n" %self.outdir.dir)
        fh.write("Space group= %s \n" %self.sg.value)
        fh.write("Sites= %s \n" %self.sites.textbox.text())
        fh.write("Resolution= %s\n" %self.resol.textbox.text())
        fh.write("Emins= %s\n" %self.emins.textbox.text())
        fh.write("#trials= %s\n" %self.tries.textbox.text())
        fh.write("Email-id= %s\n" %self.emails.value)
        fh.write('\n\n')
        fh.close()

    def closing(self):
        username = sub.check_output(["whoami"], shell=True)
        sub.call(["scancel -u "+username], shell=True)

        print "Job submission is cancelled.."

class Ui_Dialog(Qt.QDialog):
    def __init__(self, dbConnection):
        Qt.QDialog.__init__(self)
        global c
        c = dbConnection

class plotting_ui(Qt.QWidget):
    def __init__(self):
       Qt.QWidget.__init__(self)
       self.setWindowTitle("Plot-me")

       layout = Qt.QVBoxLayout()
       self.datafile = files_widget("Input file")
       self.datafile.setToolTip("Browse filename .log or .res based on graphics")
       layout.addWidget(self.datafile)
       self.datafolder = dir_widget("folder")
       self.datafolder.setToolTip("Browse folder name for CC cluster plot")
       layout.addWidget(self.datafolder)

       #self.graphics = widgets("graphics")
       self.graphics = dropmenu("graphics")
       self.graphics.setToolTip("mention the plot type and hit <enter>")
       self.rows = widgets("Rows")
       self.cols = widgets("Columns")
       self.plot_clear_btn = two_buttons("Plot", "Clear")
       self.plot_clear_btn.Button1.clicked.connect(self.call_plotter)
       self.plot_clear_btn.Button2.clicked.connect(self.cleaner)

       layout.addWidget(self.graphics)
       layout.addWidget(self.rows)
       layout.addWidget(self.cols)
       layout.addWidget(self.plot_clear_btn)
       self.setLayout(layout)

    def call_plotter(self):
        if self.graphics.value == 'cluster':
           try:
               self.ap = Qt.QDialog()
               self.ap.ui = psx.Main(self.datafolder.dir, self.graphics.value, self.rows.value, self.cols.value)
               self.ap.ui.show()
           except NameError:
               pass

        if self.graphics.value == 'CCall_CCweak':
           try:
               self.ap = Qt.QDialog()
               self.ap.ui = psx.Main(self.datafile.fname, self.graphics.value)
               self.ap.ui.show()
           except NameError:
               pass
        if self.graphics.value == 'shelxc':
           try:
              self.ap = Qt.QDialog()
              self.ap.ui = psx.Main(self.datafile.fname, self.graphics.value)
              self.ap.ui.show()
           except NameError:
               pass
        if self.graphics.value == 'site_occupancy':
           try:
              self.ap = Qt.QDialog()
              self.ap.ui = psx.Main(self.datafile.fname, self.graphics.value)
              self.ap.ui.show()
           except NameError:
               pass

    '''
    def call_plotter(self):
        if self.graphics.value == 'CCall_CCweak':
            try:
                CCall, CCweak = psx.get_CCs(self.datafile.fname)
                psx.plot_shelxd(CCweak, CCall, 'CCweak', 'CCall')
            except NameError:
                pass
        elif self.graphics.value == 'site_occupancy':
             occu, peaks, CCall, CCweak, fom  = psx.get_sites(self.datafile.fname)
             print('CCall: %s, CCweak: %s, CFOM: %s\n' % (CCall, CCweak, fom))
             psx.plot_shelxd(peaks, occu, 'Peak Number', 'Site Occupancy')

        elif self.graphics.value == 'cluster':
             if self.datafile.fname != None:
                print('Need a folder name \n')

             elif (self.rows.value or self.cols.value) == None:
                 print('provide me the number of rows and columns for cluster plotting \n')
             else:
                 psx.plot_cc(self.datafile.fname, self.rows.value, self.cols.value)

        elif self.graphics.value == 'chi_sq':
            stat_array = psx.get_shelxc_stats(self.datafile.fname)
            psx.plot_shelxc(stat_array[0,:], stat_array[1,:], 'Resolution', 'Chi^2')

        elif self.graphics.value == 'I/sig':
            stat_array = psx.get_shelxc_stats(self.datafile.fname)
            psx.plot_shelxc(stat_array[0,:], stat_array[2,:], 'Resolution', 'I/sig')

        elif self.graphics.value == 'completeness':
            stat_array = psx.get_shelxc_stats(self.datafile.fname)
            psx.plot_shelxc(stat_array[0,:], stat_array[3,:], 'Resolution', 'Completeness(%)')

        elif self.graphics.value == 'multiplicity':
            stat_array = psx.get_shelxc_stats(self.datafile.fname)
            psx.plot_shelxc(stat_array[0,:], stat_array[4,:], 'Resolution', 'Multiplicity')
        elif self.graphics.value == 'Rpim':
            stat_array = psx.get_shelxc_stats(self.datafile.fname)
            psx.plot_shelxc(stat_array[0,:], stat_array[5,:], 'Resolution', 'Rpim')

        elif self.graphics.value == 'Anomalous_signal':
            stat_array = psx.get_shelxc_stats(self.datafile.fname)
            psx.plot_shelxc(stat_array[0,:], stat_array[6,:], 'Resolution', 'd/sig')

        elif self.graphics.value == 'CC1/2':
            stat_array = psx.get_shelxc_stats(self.datafile.fname)
            psx.plot_shelxc(stat_array[0,:], stat_array[7,:], 'Resolution', 'CC(1/2)')
    '''
    def cleaner(self):
        self.datafile.file_box.setText(" ")
        self.graphics.textbox.setText(" ")
        self.rows.textbox.setText(" ")
        self.cols.textbox.setText(" ")



class widgets(Qt.QWidget):

    def __init__(self, strings):
        Qt.QWidget.__init__(self)
        self.label = strings
        self.value = None

        layout = Qt.QHBoxLayout()
        layout.addWidget(widgets.createLabel(self.label))
        self.textbox = widgets.createBox()
        self.textbox.returnPressed.connect(lambda: self.getter())
        layout.addWidget(self.textbox)

        self.setLayout(layout)

    @classmethod
    def createLabel(widgets, string1):
        lab = Qt.QLabel(string1)
        font = lab.font()
        font.setBold(True)
        font.setPointSize(14)

        lab.setFont(font)
        return lab

    @classmethod
    def createBox(widgets):
        box = Qt.QLineEdit()
        box.resize(100,20)

        return box

    @classmethod
    def createButton(widgets, btn_name, func, *args):
        button = Qt.QPushButton(btn_name)
        button.connect(button, QC.SIGNAL("clicked()"), func)
        return button


    def getter(self):

        try:
            self.value = eval(str(self.textbox.text()))

        except NameError:
            self.value = str(self.textbox.text())
        except SyntaxError:
              pass

        print self.value

class dropmenu(Qt.QWidget):
    def __init__(self, string1):
        Qt.QWidget.__init__(self)
        self.label = string1;
        self.value = None

        layout = Qt.QHBoxLayout()
        layout.addWidget(widgets.createLabel(self.label))
        self.comb = Qt.QComboBox(self)
        self.comb.addItems([' ','CCall_CCweak', 'cluster', 'shelxc', 'site_occupancy'])
        self.comb.activated[str].connect(self.get_text_combo)
        layout.addWidget(self.comb)
        self.setLayout(layout)

    def get_text_combo(self):

        try:
            self.value = unicode(self.comb.currentText())  
        except:
            pass
        print self.value

class files_widget(Qt.QWidget):
    def __init__(self, string1):
        Qt.QWidget.__init__(self)
        self.fname = None;

        self.layout = Qt.QHBoxLayout()
        self.layout.addWidget(widgets.createLabel(string1))
        self.file_box = widgets.createBox()

        self.brw_btn = widgets.createButton("Browse", self.Loader)
        self.rm_btn = widgets.createButton("Remove", self.cleaner)
        self.layout.addWidget(self.file_box)
        self.layout.addWidget(self.brw_btn)
        self.layout.addWidget(self.rm_btn)
        self.setLayout(self.layout)

    def Loader(self):
        self.file_box.setFocus()
        filename = Qt.QFileDialog.getOpenFileName(self, 'Open File', '~/')
        self.file_box.setText(filename)
        try:
            self.fname = unicode(self.file_box.text())
            print "filename received %s" %self.fname
            if len(self.fname) > 72:
               msgbox = Qt.QMessageBox()
               msgbox.setText("Error: shelx hates filename with >72 character!")
               msgbox.exec_()
              
        except:
            pass

    def cleaner(self):

        self.file_box.setText(None)
        if self.fname != None:
            self.fname = None

class dir_widget(Qt.QWidget):
    def __init__(self,dirname):
        Qt.QWidget.__init__(self)
        self.dir = None;

        layout = Qt.QHBoxLayout()
        layout.addWidget(widgets.createLabel(dirname))
        self.dir_box = widgets.createBox()
        self.brw_btn = widgets.createButton("Browse", self.Load_dir)
        self.rm_btn = widgets.createButton("Remove", self.cleaner)
        layout.addWidget(self.dir_box)
        layout.addWidget(self.brw_btn)
        layout.addWidget(self.rm_btn)
        self.setLayout(layout)

    def Load_dir(self):
        self.dir_box.setFocus()
        directory = Qt.QFileDialog.getExistingDirectory(self)
        self.dir_box.setText(directory)
        try:
            self.dir = unicode(self.dir_box.text())
            print "output directory name: %s" %self.dir
        except:
            pass

    def cleaner(self):

        self.dir_box.setText(None)
        if self.dir != None:
            self.dir = None

class two_buttons(Qt.QWidget):

    def __init__(self, butn1, butn2):
        Qt.QWidget.__init__(self)
        layout = Qt.QHBoxLayout()
        self.Button1 = Qt.QPushButton(butn1)
        self.Button2 = Qt.QPushButton(butn2)

        layout.addWidget(self.Button1)
        layout.addWidget(self.Button2)
        self.setLayout(layout)


def main():
    app = Qt.QApplication(sys.argv)
    #w = MainLayout()
    w = Tabs()
    w.show()
    return app.exec_()

if __name__ == '__main__':
   main()
