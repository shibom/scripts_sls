#include "stdio.h"
#include "string.h"
#include "ctype.h"
#include <stdint.h>
#include <iostream>
#include <hdf5.h>

const unsigned int MaxStrLen = 255;
const unsigned int MaxHeader = 10000;
const float Threshold = 0;
const float PanelGapValue = -1e10;

#define dtype float

//extern "C" bool ReadCBFfile(char* fname, dtype** outArray, int* dims, float* pixel,
//                 float* expo, float* waveLen, float* dist, float* beamxy,
//                 float* flux, unsigned short** badmask)
extern "C" dtype* ReadCBFfile(char* fname)
{
  //int* dims; float* pixel; float* expo; float* waveLen; float* dist;
  //float* beamxy; float* flux; unsigned short** badmask;
  float pixel[2] = {0,0};
  float expo = 0;
  float waveLen = 0;
  float dist = 0;
  float beamxy[2] = {0,0};
  float flux = 0;
  int dims[3] = {0,0,0};
  dtype* outArray;

  FILE* inpCBF = fopen(fname,"rb");
  if (inpCBF == NULL) exit(-1);
  char* headerAr = new char[MaxHeader];
  signed char aByte;
  size_t cInd = 0;
  // Searching for the data and reading the header
  bool notCBF = false;
  while (!feof(inpCBF))
  {
    fread(&aByte,1,1,inpCBF);
    if (aByte==0xC)
    { unsigned char ch3[3];
      fread(&ch3,3,1,inpCBF);
      if (ch3[0]==0x1A && ch3[1]==0x04 && ch3[2]==0xD5) break;
    }
    //here header just copied to a big array of chars here?
    if (cInd<MaxHeader)
    { headerAr[cInd] = aByte;
      cInd++;
    } else
    { printf("Big Header! Is it realy CBF?\n");
      notCBF = true;
      break;
    }

  }
  if (feof(inpCBF) || notCBF)
  { printf("Data not found in %s\n",fname);
    exit(-1);
  }
  // Here analyse the header
  // (can be combined with the previous one)
  char str1[MaxStrLen];
  int headerLen = cInd;
  cInd = 0;
  while (cInd<headerLen)
  { int sPo = 0;
    while (headerAr[cInd]!='\n' && cInd<headerLen && sPo<MaxStrLen-1)
    { str1[sPo] = headerAr[cInd];
      cInd++;
      sPo++;
    }
    str1[sPo]=0;
    cInd++;
    sscanf(str1,"# Pixel_size %f m x %f m", pixel,pixel+1);
    sscanf(str1,"# Exposure_time %f", expo);
    sscanf(str1,"# Wavelength %f", waveLen);
    sscanf(str1,"# Detector_distance %f", dist);
    sscanf(str1,"# Beam_xy (%f, %f) pixels", beamxy, beamxy+1);
    sscanf(str1,"# Flux %f", flux);
    sscanf(str1,"X-Binary-Size-Fastest-Dimension: %d", dims);
    sscanf(str1,"X-Binary-Size-Second-Dimension: %d", dims+1);
  }
  size_t totLen = dims[0]*dims[1];
  if (totLen<1)
  { printf("Some Dimentions are 0!\n");
    exit(-1);
  }
  delete [] headerAr;

  outArray = new dtype[totLen];
  //int* outAr = *outArray;
  //std::vector<int> outArray(totLen);
//mask  *badmask = new unsigned short[totLen];
//mask  unsigned short* badMa = *badmask;

  dtype cVal = 0;                       // what about float? Make more universal?
  // The main reading of the CBS file
  signed short aWord;
  int anInt = sizeof(char);
  cInd = 0;
  while (!feof(inpCBF) && cInd<totLen)
  { fread(&aByte,1,1,inpCBF);
    if (aByte == -128)
    { fread(&aWord,2,1,inpCBF);
      if (aWord == -32768)
      { fread(&anInt,4,1,inpCBF);
        cVal += anInt;
      } else cVal += aWord;
    } else cVal += aByte;
    outArray[cInd] = (dtype)cVal;
//doesn't work :(    if (cVal<0) outAr[cInd] = PanelGapValue;
//mask    if (cVal<Threshold-1e-10)
//mask      badMa[cInd] = 0;  //mask negative intensities
//mask    else
//mask      badMa[cInd] = 1;
    cInd++;
  }
  fclose(inpCBF);

  return outArray;
}

extern "C" bool write_h5(char* fname, dtype* data_array){
  #ifdef _HDF5_H
        hid_t out_type_id = H5T_NATIVE_FLOAT;
        hid_t file_id = H5Fcreate(fname,  H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        hid_t dataspace_id, dataset_id;
        int ndims = 2; hsize_t dimsH[2];
        dimsH[0] = 1679;
        dimsH[1] = 1475;
        hid_t gid2 = H5Gcreate1(file_id,"data",0);
        dataspace_id = H5Screate_simple(ndims, dimsH, NULL);
        bool compress = false;
        hid_t dcpl = H5Pcreate (H5P_DATASET_CREATE);
        hsize_t chunk[2] = {64, 64};
        if (compress)
        { H5Pset_deflate (dcpl, 9);
          H5Pset_chunk (dcpl, 2, chunk);
        }



        dataset_id = H5Dcreate(file_id, "/data/data", out_type_id, dataspace_id,
                     H5P_DEFAULT, dcpl, H5P_DEFAULT);
        H5Dwrite(dataset_id, out_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_array);
        H5Dclose(dataset_id);
        H5Sclose(dataspace_id);
        H5Gclose(gid2);
        H5Fclose(file_id);
        return true;
  #endif
        return false;

}
extern "C" int main(int argc, char* argv[])
{
  dtype* sum = NULL;
  sum = new dtype[2476525];
  FILE* StList; size_t numFi=0; char str1[MaxStrLen];
  char fName[MaxStrLen];
  strcpy(fName, argv[1]);

  StList = fopen(fName, "rt");
  char* filesList= new char[MaxStrLen*numFi];

  while (!feof(StList))
  { fgets(str1,MaxStrLen,StList);
      if (feof(StList)) break;
      sscanf(str1,"%s",str1);
      sprintf(&filesList[numFi*MaxStrLen],"%s\0",str1);
      numFi++;
  }
  printf("Found %ld files in the list %s\n",numFi,fName);
  fclose(StList);
  for (size_t ifi=0; ifi<numFi; ifi++)
  {
  char justName[MaxStrLen];
  char fName1[MaxStrLen];
  char extention[MaxStrLen];
  char comment[MaxStrLen];
  strncpy(fName1,&filesList[ifi*MaxStrLen],MaxStrLen);
  strcpy(justName,(strrchr(fName1,'/')!=NULL?strrchr(fName1,'/')+1:(strrchr(fName1,'\\')!=NULL?strrchr(fName1,'\\')+1:fName1)));
  if (strrchr(justName,'.')!=NULL) strrchr(justName,'.')[0] = 0;

  strcpy(extention,(strrchr(fName1,'.')!=NULL?strrchr(fName1,'.')+1:""));
  for (int i=0; i<strlen(extention); i++)
  { extention[i] = tolower(extention[i]);
    if (extention[i]==' ' || extention[i] == '\\' || extention[i] == '/' || extention[i] == '\n' || extention[i] == '\t')
    { extention[i] = 0;
      break;
    }
  }
  unsigned short *badMask;
  //char* fName1;
  //fName1 = "testshot_1_00001.cbf";

  dtype* data = NULL;
  data = ReadCBFfile(fName1);//, dims, pixel, &expo, &waveLen, &dist, beamxy, &flux, &badMask);
  for (int i=0; i<2476525; i++){
    sum[i] += data[i];
  }
  #ifdef _HDF5_H
        write_h5("out.h5", sum);
  #endif
}

}
