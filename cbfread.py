import sys, os
import numpy as np
from iotbx.xds import xds_cbf
import matplotlib.pyplot as plt
import glob
import h5py
import h5read

def list_file(dirname):

    file_list = [];
    os.chdir(dirname)
    for file in glob.glob("*.cbf"):
        file_list.append(file)
    file_list.sort()
   
    return file_list

#def read_bg(fname):
    '''
    fh = xds_cbf.reader()
    fh.read_file("BKGPIX.cbf")
    bg = fh.get_data()
    
    
    f = h5py.File(fname, 'r')
    bg = f["/data/data"]
    bg = np.copy(bg)
    '''
#    return bg

def powder():
    
    lists = list_file(sys.argv[1])
    summed = np.empty((1679,1475))
    #summed = np.zeros((2476525))
    avg_img = np.empty((len(lists),1475))
    bg = h5read.read_bg()

    for j in range(len(lists)):
        '''
        handle = ImageFactory(name)
        handle.read()
        image = handle.get_raw_data()
        '''
        handle = xds_cbf.reader()
        handle.read_file(lists[j])
        image = handle.get_data()
        image = image - bg
        for i in range(image.shape[1]):
            avg_img[j][i] = np.sum(image[:,i])        
        
        
    return avg_img

summed = powder()
summed[summed > 20000] = 0
summed[summed < 0 ] = 0
    
f = h5py.File(sys.argv[2], "w")
print "writing h5 file \n"
data = f.create_group("data")
data.create_dataset("data", data=summed)
f.close()
cols = summed.shape[1]
rad = [];
for i in range(cols):
    tmp = np.mean(summed[:,i])
    rad.append(tmp)

plt.plot(rad)
plt.xlabel("#pixels", fontweight='bold', fontsize = 14)
plt.ylabel("average_intensity", fontweight='bold', fontsize = 14)
plt.show()

