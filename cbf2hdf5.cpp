
// cbf2hdf5.cpp
//
// Convert a Pilatus .CBF file into an HDF5 file
// Also a simple hit-finder
//
// Copyright � 2014 Deutsches Elektronen-Synchrotron DESY,
//                  a research centre of the Helmholtz Association.
//
// Author: Oleksandr Yefanov (oleksandr.yefanov@desy.de)

#include "stdio.h"
#include "string.h"
#include "ctype.h"
#include <hdf5.h>    //can be commented to use .cbf -> .raw/.tif converter
#include "stdlib.h"
#include <stdint.h>

#ifdef _OPENMP
#include <omp.h>
omp_lock_t  lock;
#endif

const unsigned int MaxStrLen = 255;
const unsigned int MaxHeader = 10000;
const float Threshold = 0;
const float PanelGapValue = -1e10;

#define dtype float

bool ReadCBFfile(char* fname, dtype** outArray, int* dims, float* pixel,
                 float* expo, float* waveLen, float* dist, float* beamxy,
                 float* flux, unsigned short** badmask)
{
  FILE* inpCBF = fopen(fname,"rb");
  if (inpCBF == NULL) return false;
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
    return false;
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
    return false;
  }
  delete headerAr;

  *outArray = new dtype[totLen];
  dtype* outAr = *outArray;

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
    outAr[cInd] = (dtype)cVal;
//doesn't work :(    if (cVal<0) outAr[cInd] = PanelGapValue;
//mask    if (cVal<Threshold-1e-10)
//mask      badMa[cInd] = 0;  //mask negative intensities
//mask    else
//mask      badMa[cInd] = 1;
    cInd++;
  }
  fclose(inpCBF);
  
  return true;
}

bool TIFFloatWriter(char* fnam, float* data, int dx, int dy, char* comment)
{ char aChar;
  short aWord;
  int anInt;
  int commentLen = strlen(comment)+1;
  FILE* tifF = fopen(fnam,"wb");
  if (tifF == NULL) return false;
  //Header
  fputc(0x49,tifF);                // I can check endians - still from CASS
  fputc(0x49,tifF);
  aWord = 42;
  fwrite(&aWord,sizeof(aWord),1,tifF);
  anInt = 8; //offset
  fwrite(&anInt,sizeof(anInt),1,tifF);

  int numFields = 9;
  aWord = numFields; // Num of fields in IFD
  fwrite(&aWord,sizeof(aWord),1,tifF);

  aWord = 254; //NewSubfileType ???
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 4; // type int
  fwrite(&aWord,sizeof(aWord),1,tifF);
  anInt = 1; // num elements
  fwrite(&anInt,sizeof(anInt),1,tifF);
  anInt = 0;
  fwrite(&anInt,sizeof(anInt),1,tifF);

  aWord = 256; // Width
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 4; // type int
  fwrite(&aWord,sizeof(aWord),1,tifF);
  anInt = 1; // num elements
  fwrite(&anInt,sizeof(anInt),1,tifF);
  anInt = dx; // the width
  fwrite(&anInt,sizeof(anInt),1,tifF);

  aWord = 257; // Height
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 4; // type int
  fwrite(&aWord,sizeof(aWord),1,tifF);
  anInt = 1; // num elements
  fwrite(&anInt,sizeof(anInt),1,tifF);
  anInt = dy; // the height
  fwrite(&anInt,sizeof(anInt),1,tifF);

  aWord = 258; // bits per sample
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 3; // type short
  fwrite(&aWord,sizeof(aWord),1,tifF);
  anInt = 1; // num elements
  fwrite(&anInt,sizeof(anInt),1,tifF);
  aWord = 32; //
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 0; // 12 byte
  fwrite(&aWord,sizeof(aWord),1,tifF);

  aWord = 262; // Photometric Interpretation ???
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 3; // type short
  fwrite(&aWord,sizeof(aWord),1,tifF);
  anInt = 1; // num elements
  fwrite(&anInt,sizeof(anInt),1,tifF);
  aWord = 1; // ?
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 0; // 12 byte
  fwrite(&aWord,sizeof(aWord),1,tifF);

  aWord = 273; // offset to the data
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 4; // type int
  fwrite(&aWord,sizeof(aWord),1,tifF);
  anInt = 1; // num elements
  fwrite(&anInt,sizeof(anInt),1,tifF);
  anInt = 14+numFields*12+commentLen; // calculate!
  fwrite(&anInt,sizeof(anInt),1,tifF);

  aWord = 279; // data size
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 4; // type int
  fwrite(&aWord,sizeof(aWord),1,tifF);
  anInt = 1; // num elements
  fwrite(&anInt,sizeof(anInt),1,tifF);
  anInt = dx*dy*sizeof(float); // calculate!
  fwrite(&anInt,sizeof(anInt),1,tifF);

  aWord = 339; // Sample Format
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 3; // type short
  fwrite(&aWord,sizeof(aWord),1,tifF);
  anInt = 1; // num elements
  fwrite(&anInt,sizeof(anInt),1,tifF);
  aWord = 3; // IEEE floar
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 0; // 12 byte
  fwrite(&aWord,sizeof(aWord),1,tifF);

  aWord = 305; // Software
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 2; // type string
  fwrite(&aWord,sizeof(aWord),1,tifF);
  anInt = commentLen; // the length
  fwrite(&anInt,sizeof(anInt),1,tifF);
  anInt = 14+numFields*12; // offset - easy to calc
  fwrite(&anInt,sizeof(anInt),1,tifF);

  anInt = 0; // next IFD - not present
  fwrite(&anInt,sizeof(anInt),1,tifF);

  //write comment
  for (int i=0; i<commentLen; i++)
    fputc(comment[i],tifF);
//  fprintf(tifF,"%s",comment);

  // Just saving the data
  size_t numel = (size_t)dx * (size_t)dy;
  for (size_t i=0; i<numel; i++)
    fwrite(&data[i],sizeof(float),1,tifF);

  fclose(tifF);
  return true;
}

bool TIFFloatReader(char* fnam, float** data, int* dx, int* dy, char* comment)
{ char aChar;
  short aWord;
  int anInt;
  FILE* tifF = fopen(fnam,"rb");
  if (tifF==NULL) return false;
  fread(&anInt,4,1,tifF); //Header - could be analysed for little/big endians
  fread(&anInt,4,1,tifF); // offset of IFD
  fseek(tifF,anInt,SEEK_SET);

  fread(&aWord,2,1,tifF); // number of fields
  int numFields = aWord;
  if (numFields<3) return false;

  struct tFields
  { short tag;
    short type;
    int numEl;
    short valS1;
    short valS2;
    int valI;
    char valC[4];
    float valF;
  };

  tFields*  aField = new tFields[numFields];
  for (int i=0; i<numFields; i++)  //Reading fields in IFD
  {
    fread(&aField[i].tag,2,1,tifF);
    fread(&aField[i].type,2,1,tifF);
    fread(&aField[i].numEl,4,1,tifF);
    if (aField[i].type == 3)
    { fread(&aField[i].valS1,2,1,tifF);
      fread(&aField[i].valS2,2,1,tifF);
    } else if (aField[i].type == 4 || aField[i].type == 2)
      fread(&aField[i].valI,4,1,tifF);
      else if (aField[i].type == 1)
      fread(aField[i].valC,1,4,tifF);
      else if (aField[i].type == 11)
      fread(&aField[i].valF,4,1,tifF);
  }

  int imOffset = 0;
  int dataSize = 0;
  float dataType = 0;
  float rightType = true;
  for (int i=0; i<numFields; i++)  //Analysing fields in IFD
  { if (aField[i].tag == 256) *dx = aField[i].valI;
    else if (aField[i].tag == 257) *dy = aField[i].valI;
    else if (aField[i].tag == 273) imOffset = aField[i].valI;
    else if (aField[i].tag == 258)
      rightType = (aField[i].valS1==32);
    else if (aField[i].tag == 339)
      dataType = aField[i].valS1;
    else if (aField[i].tag == 279)
      dataSize = aField[i].valI;
    else if (aField[i].tag == 305) // This is string reader
    { fseek(tifF,aField[i].valI,SEEK_SET);
      if (aField[i].numEl>MaxStrLen) aField[i].numEl=MaxStrLen;
      for (size_t k=0; k<aField[i].numEl; k++)
        fread(&comment[k],1,1,tifF);
//      aField[i].valI;
    }
  }
  if (!rightType || (dataType!=2 && dataType!=3))
  { printf("Supported data types are: float or int (both 32bit)!\n");
    return false;
  }

  // creating array for the image
  if (imOffset<8 || *dx<1 || *dy<1) return false; //???
  size_t numEl = *dx * *dy;
  if (numEl<1) return false;
#ifdef linux
  *data = (float*)valloc(numEl*sizeof(float));
#else   //090615
  *data = new float[numEl];
#endif
  float* _Ar = *data;

  // Reading the actual image
  fseek(tifF,imOffset,SEEK_SET);
  if (dataType==3)
    fread(_Ar,sizeof(float),numEl,tifF);
  else if (dataType==2)
  { int32_t aval;
    for (int i=0; i<numEl; i++)
    { fread(&aval,4,1,tifF);
      _Ar[i] = (float)aval;
    }
  }

  delete aField;
  fclose(tifF);
  return true;
}

bool RAWFloatReader(char* fnam, float** data, int* numEl)
{
  FILE* rawF = fopen(fnam,"rb");
  if (rawF==NULL) return false;
#ifdef linux
  *data = (float*)valloc(*numEl*sizeof(float));
#else   //090615
  *data = new float[*numEl];
#endif
  float* _Ar = *data;

  int i = 0;
  while (!feof(rawF))
  { fread(&_Ar[i],sizeof(float),1,rawF);
    i++;
    if (i>=*numEl) break;
  }
  *numEl = i;
  fclose(rawF);
  return true;
}

bool WriteHDF5file(char* fname, dtype* outArray, int* dims, float* pixel,
                float* expo, float* waveLen, float* dist, float* beamxy,
                float* flux, bool compress, unsigned short* badmask)
{
#ifdef _HDF5_H
      hid_t out_type_id = H5T_NATIVE_FLOAT;
      hid_t file_id = H5Fcreate(fname,  H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      hid_t gid = H5Gcreate1(file_id,"LCLS",0);
      hid_t dataspace_id, dataset_id;
      int ndims = 1;
      hsize_t dimsH[2];
      dimsH[0] = 1;
      dimsH[1] = 0;
      dataspace_id = H5Screate_simple(1, dimsH, NULL);

      dataset_id = H5Dcreate1(gid, "detectorPosition", out_type_id, dataspace_id, H5P_DEFAULT);
      if(H5Dwrite(dataset_id, out_type_id, H5S_ALL, H5S_ALL,H5P_DEFAULT, dist)< 0)
        printf("Error writing 1D data to file\n");
      H5Dclose(dataset_id);

      dataset_id = H5Dcreate1(gid, "photon_wavelength_A", out_type_id, dataspace_id, H5P_DEFAULT);
      if(H5Dwrite(dataset_id, out_type_id, H5S_ALL, H5S_ALL,H5P_DEFAULT, waveLen)< 0)
        printf("Error writing 1D data to file\n");
      H5Dclose(dataset_id);

      dataset_id = H5Dcreate1(gid, "exposure_s", out_type_id, dataspace_id, H5P_DEFAULT);
      if(H5Dwrite(dataset_id, out_type_id, H5S_ALL, H5S_ALL,H5P_DEFAULT, expo)< 0)
        printf("Error writing 1D data to file\n");
      H5Dclose(dataset_id);

      dataset_id = H5Dcreate1(gid, "pixelX_m", out_type_id, dataspace_id, H5P_DEFAULT);
      if(H5Dwrite(dataset_id, out_type_id, H5S_ALL, H5S_ALL,H5P_DEFAULT, &pixel[0])< 0)
        printf("Error writing 1D data to file\n");
      H5Dclose(dataset_id);

      dataset_id = H5Dcreate1(gid, "pixelY_m", out_type_id, dataspace_id, H5P_DEFAULT);
      if(H5Dwrite(dataset_id, out_type_id, H5S_ALL, H5S_ALL,H5P_DEFAULT, &pixel[1])< 0)
        printf("Error writing 1D data to file\n");
      H5Dclose(dataset_id);

      dataset_id = H5Dcreate1(gid, "flux_ph_s", out_type_id, dataspace_id, H5P_DEFAULT);
      if(H5Dwrite(dataset_id, out_type_id, H5S_ALL, H5S_ALL,H5P_DEFAULT, flux)< 0)
        printf("Error writing 1D data to file\n");
      H5Dclose(dataset_id);

      dataset_id = H5Dcreate1(gid, "beamX_px", out_type_id, dataspace_id, H5P_DEFAULT);
      if(H5Dwrite(dataset_id, out_type_id, H5S_ALL, H5S_ALL,H5P_DEFAULT, &beamxy[0])< 0)
        printf("Error writing 1D data to file\n");
      H5Dclose(dataset_id);

      dataset_id = H5Dcreate1(gid, "beamY_px", out_type_id, dataspace_id, H5P_DEFAULT);
      if(H5Dwrite(dataset_id, out_type_id, H5S_ALL, H5S_ALL,H5P_DEFAULT, &beamxy[1])< 0)
        printf("Error writing 1D data to file\n");
      H5Dclose(dataset_id);
      H5Sclose(dataspace_id);

      ndims = 2;
      dimsH[0] = dims[1];
      dimsH[1] = dims[0];
      hid_t gid2 = H5Gcreate1(file_id,"data",0);
      dataspace_id = H5Screate_simple(ndims, dimsH, NULL);
      //COMPRESSION
      hid_t dcpl = H5Pcreate (H5P_DATASET_CREATE);
//      hsize_t chunk[2] = {dimsH[0], dimsH[1]};
//#ifdef ZLIB_H
//CASS    hsize_t chunk[2] = {40,2};
      hsize_t chunk[2] = {64, 64};
      if (compress)
      { H5Pset_deflate (dcpl, 9);
        H5Pset_chunk (dcpl, 2, chunk);
      }
//#endif
 //     H5Pset_szip (dcpl, H5_SZIP_NN_OPTION_MASK, 8);
 //     H5Pset_chunk (dcpl, 2, chunk);

      dataset_id = H5Dcreate(file_id, "/data/data", out_type_id, dataspace_id,
                   H5P_DEFAULT, dcpl, H5P_DEFAULT);
//      dataset_id = H5Dcreate1(file_id, "/data/data", out_type_id, dataspace_id, H5P_DEFAULT);
      if(H5Dwrite(dataset_id, out_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, outArray)< 0)
        printf("Error writing 2D data to file\n");
      H5Dclose(dataset_id);
      H5Sclose(dataspace_id);
      H5Gclose(gid2);

      //Here badmask

      H5Gclose(gid);
      H5Fclose(file_id);
      return true;
#endif
      return false;
}


bool ReadHDF5file(char* filename, char* fieldname, dtype** outArray, int* dims)
{
#ifdef _HDF5_H
	// Open the file
	hid_t file_id;
	file_id = H5Fopen(filename,H5F_ACC_RDONLY,H5P_DEFAULT);
//?	file_id = H5Fopen(filename,H5F_ACC_RDONLY,faplist_id);
	if(file_id < 0){
		printf("ERROR: Could not open file %s\n",filename);
		return false;
	}

	// Open the dataset
	hid_t dataset_id;
	hid_t dataspace_id;
	dataset_id = H5Dopen1(file_id, fieldname);
	if(dataset_id < 0){
		printf("ERROR: Could not open the data field %s\n",fieldname);
		return false;
	}

	dataspace_id = H5Dget_space(dataset_id);

	// Test if 2D data
	int ndims;
	ndims = H5Sget_simple_extent_ndims(dataspace_id);

	// Get dimensions of data set (nx, ny, nn)
	hsize_t* dimsl = new hsize_t[ndims];
	H5Sget_simple_extent_dims(dataspace_id,dimsl,NULL);
	for(int i = 0;(i<ndims&&i<3);i++)
	  dims[i] = dimsl[ndims-1-i];  //!!!!!!!! NOT SURE
	size_t nn = 1;
	for(int i = 0;i<ndims;i++)
		nn *= dimsl[i];

	// Create space for the new data
    dtype* data = *outArray;
    if (data!=NULL) delete data;//free(data);
    *outArray = new dtype[nn];
    data = *outArray;

	hid_t		datatype_id;
	H5T_class_t dataclass;
	size_t size;
	datatype_id =  H5Dget_type(dataset_id);
	dataclass = H5Tget_class(datatype_id);
	size = H5Tget_size(datatype_id);
    int rrr = sizeof(int);
	if(dataclass == H5T_FLOAT){
		if (size == sizeof(float)) {
			float* buffer = (float *) calloc(nn, sizeof(float));
			H5Dread(dataset_id, datatype_id, H5S_ALL,H5S_ALL, H5P_DEFAULT, buffer);
			for(long i=0; i<nn; i++)
				data[i] = buffer[i];
			free(buffer);
		}
		else if (size == sizeof(double)) {
			double* buffer = (double *) calloc(nn, sizeof(double));
			H5Dread(dataset_id, datatype_id, H5S_ALL,H5S_ALL, H5P_DEFAULT, buffer);
			for(long i=0; i<nn; i++)
				data[i] = buffer[i];
			free(buffer);
		}
		else {
			printf("2dData::readHDF5: unknown floating point type, size=%i\n",(int) size);
			return false;
		}
	}
	else if(dataclass == H5T_INTEGER){
		if (size == sizeof(short)) {
			short* buffer = (short*) calloc(nn, sizeof(short));
			H5Dread(dataset_id, datatype_id, H5S_ALL,H5S_ALL, H5P_DEFAULT, buffer);
			for(long i=0; i<nn; i++)
				data[i] = buffer[i];
			free(buffer);
		}
		else if (size == sizeof(int)) {
			int* buffer = (int *) calloc(nn, sizeof(int));
			H5Dread(dataset_id, datatype_id, H5S_ALL,H5S_ALL, H5P_DEFAULT, buffer);
			for(long i=0; i<nn; i++)
				data[i] = buffer[i];
			free(buffer);
		}
		else if (size == sizeof(long)) {
			long* buffer = (long *) calloc(nn, sizeof(long));
			H5Dread(dataset_id, datatype_id, H5S_ALL,H5S_ALL, H5P_DEFAULT, buffer);
			for(long i=0; i<nn; i++)
				data[i] = buffer[i];
			free(buffer);
		}
		else {
			printf("2dData::readHDF5: unknown integer type, size=%i\n",(int) size);
			return false;
		}
	}
	else {
		printf("2dData::readHDF5: unknown HDF5 data type\n");
		return false;
	}

	// Close and cleanup
	H5Dclose(dataset_id);


	// Cleanup stale IDs
	hid_t ids[256];
	int n_ids = H5Fget_obj_ids(file_id, H5F_OBJ_ALL, 256, ids);
	for (long i=0; i<n_ids; i++ ) {

		hid_t id;
		H5I_type_t type;

		id = ids[i];
		type = H5Iget_type(id);

		if ( type == H5I_GROUP )
			H5Gclose(id);
		if ( type == H5I_DATASET )
			H5Dclose(id);
		if ( type == H5I_DATASPACE )
			H5Sclose(id);
		//if ( type == H5I_DATATYPE )
		//	H5Dclose(id);
	}

	H5Fclose(file_id);
    return true;
#endif
    return false;
}


#define elem_type float
#define ELEM_SWAP(a,b) { register elem_type t=(a);(a)=(b);(b)=t; }
elem_type quick_select(elem_type arr[], int n)
{
    int low, high ;
    int median;
    int middle, ll, hh;

    low = 0 ; high = n-1 ; median = (low + high) / 2;
    for (;;) {
        if (high <= low) /* One element only */
            return arr[median] ;

        if (high == low + 1) {  /* Two elements only */
            if (arr[low] > arr[high])
                ELEM_SWAP(arr[low], arr[high]) ;
            return arr[median] ;
        }

    /* Find median of low, middle and high items; swap into position low */
    middle = (low + high) / 2;
    if (arr[middle] > arr[high])    ELEM_SWAP(arr[middle], arr[high]) ;
    if (arr[low] > arr[high])       ELEM_SWAP(arr[low], arr[high]) ;
    if (arr[middle] > arr[low])     ELEM_SWAP(arr[middle], arr[low]) ;

    /* Swap low item (now in position middle) into position (low+1) */
    ELEM_SWAP(arr[middle], arr[low+1]) ;

    /* Nibble from each end towards middle, swapping items when stuck */
    ll = low + 1;
    hh = high;
    for (;;) {
        do ll++; while (arr[low] > arr[ll]) ;
        do hh--; while (arr[hh]  > arr[low]) ;

        if (hh < ll)
        break;

        ELEM_SWAP(arr[ll], arr[hh]) ;
    }

    /* Swap middle item (in position low) back into correct position */
    ELEM_SWAP(arr[low], arr[hh]) ;

    /* Re-set active partition */
    if (hh <= median)
        low = ll;
        if (hh >= median)
        high = hh - 1;
    }
}
#undef ELEM_SWAP
#undef elem_type

int LocBackSub(dtype* theArray, int* dims, int* box)
{
  float* arr1 = new float[(2*box[0]+1)*(2*box[1]+1)];
  for (int xi=0; xi<dims[0]; xi++)
    for (int yi=0; yi<dims[1]; yi++)
    { int nuel = 0;
      if (theArray[xi+dims[0]*yi]<0) continue;
      for (int bxi=-box[0]; bxi<=box[0]; bxi++)
        for (int byi=-box[1]; byi<=box[1]; byi++)
          if (xi+bxi>=0 && xi+bxi<dims[0] && yi+byi>=0 && yi+byi<dims[1])
            if (theArray[xi+bxi+dims[0]*(yi+byi)]>=0)
            { arr1[nuel] = theArray[xi+bxi+dims[0]*(yi+byi)];
              nuel++;
            }
      if (nuel>0) theArray[xi+dims[0]*yi] -= quick_select(arr1, nuel);
    }
  delete arr1;
  return 1;
}

int HitFinder(dtype* theArray, int* dims, float threshold, int numconnected, float* averint)
{ int numpx = 1;
  for (int i=0; i<2; i++) if (dims[i]>1) numpx *= dims[i];
  if (numpx<2) return 0;

//  if (locBg>0)
//  { int box[2] = {2*locBg+1, 2*locBg+1};
//  { int box[2] = {locBg, locBg};
//    LocBackSub(theArray, dims, box);
//  }

  int numfound = 0;
  if (numconnected<2)
  { *averint = 0;
    for (int i=0; i<numpx; i++)
      if (theArray[i]>threshold)
      { *averint += theArray[i];
        numfound++;
      }
    *averint /= (float)numfound;
  } else  // connected peaks
  { *averint = 0;
    int* considered = new int[numpx];
    struct TPoint
    { int x;
      int y;
      float I;
    };
    int maxNumPo = numpx/10;
    TPoint* setP = new TPoint[maxNumPo];

    for (int i=0; i<numpx; i++) considered[i] = 0;

    for (int xi=0; xi<dims[0]; xi++)
      for (int yi=0; yi<dims[1]; yi++)
        if (considered[xi+dims[0]*yi]==0)
          if (theArray[xi+dims[0]*yi]>threshold)
          { int numpoints = 0;
            setP[numpoints].x=xi;
            setP[numpoints].y=yi;
            setP[numpoints].I=theArray[xi+dims[0]*yi];
            numpoints++;
//            while (numcon<cpoint) // here loop over setP ???
            for (int cpoint=0; cpoint<numpoints; cpoint++) // here loop over setP ???
            { considered[setP[cpoint].x+dims[0]*setP[cpoint].y]=1;
              for (int bxi=-1; bxi<=1; bxi++)
                for (int byi=-1; byi<=1; byi++)
                  if (setP[cpoint].x+bxi>=0 && setP[cpoint].x+bxi<dims[0] &&
                      setP[cpoint].y+byi>=0 && setP[cpoint].y+byi<dims[1])
                    if (considered[setP[cpoint].x+bxi+dims[0]*(setP[cpoint].y+byi)]==0)
                      if (theArray[setP[cpoint].x+bxi+dims[0]*(setP[cpoint].y+byi)]>threshold)
                      { int icp;
                        for (icp=0; icp<numpoints; icp++)
                          if (setP[icp].x==setP[cpoint].x+bxi && setP[icp].y==setP[cpoint].y+byi) break;
                        if (icp<numpoints) continue;
                        setP[numpoints].x=setP[cpoint].x+bxi;
                        setP[numpoints].y=setP[cpoint].y+byi;
                        setP[numpoints].I=theArray[setP[numpoints].x+dims[0]*setP[numpoints].y];
//                        *averint += setP[cpoint].I;
                        numpoints++;
                        if (numpoints>=maxNumPo) break;
                      }
            }
            if (numpoints>=numconnected) numfound++;
          } else considered[xi+dims[0]*yi]=1;
    delete considered;
    delete setP;
  }

  return numfound;
}

int main(int argc, char* argv[])
{
  //list threshold=10.1 numpeaks=20 numconnected=5 locbgrad=2 tif save
  printf("Converter .cbf/.raw/.tif/.h5 -> .h5/.raw/.tif (float 32bit)\n");
  printf("Also simple hitfinder. Achtung! Detector position in hdf5 in [mm]\n");
//  printf("Achtung! Now (from 19.03.13) detector position in hdf5 in [mm]!\n");
//#ifdef ZLIB_H
//  printf("Using gzip compression\n");
//#endif

  if (argc<2)
  { printf("   cbf2hdf5 filelist [raw] [tif] [datafield=XXX] [locbgrad=K]\n");
    printf("           [threshold=N] [numpeaks=M] [numconnected=L] [save] [dims=x,y,z]\n");
    printf("Instead of filelist can be just a single file - like fileinp.cbf\n");
    printf("The output file(s) will be stored in curent folder\n");
    printf("Option \"raw\" - the output is binary, \"tif\" - 32 bit TIFF, otherwise .h5\n");
    printf("By default input/output datafield=/data/data, but can be changed for input\n");
    printf("For local background correction set its radius: locbgrad=K\n");
    printf("If input files are RAW (32bit float), set its dimensions dims=x,y,z\n");
//140404    printf("If option \"comp\" is set the output .h5 will be gzip compressed\n");
//140404    printf("  but it may not work in multithreaded mode\n");
//140404    printf("For OMP use (in zsh), for example, \"export OMP_NUM_THREADS=10\"\n");
//140404    printf("  it can speed up convertion 3 times.\n");
    printf("HITFINDING. A hit when \"numpeaks\" peaks with \"numconnected\" connected\n");
    printf("  pixels with values > \"threshold\". To save the patterns add \"save\".\n");
    printf("\nEXAMPLE: cbf2hdf5 files.lst threshold=10.1 numpeaks=20 numconnected=5 locbgrad=2\n");

#ifndef _HDF5_H
    printf("\nAchtung! Compiled without HDF5 library -\n");
    printf("         only conversion to .raw/.tif is possible.\n");
#endif

    return 0;
  }
  char fName[MaxStrLen];
  char datafield[MaxStrLen] = "/data/data";
//  strcpy(datafield,"/entry/instrument/detector/data");
  bool compress = false;
  bool convToRAW = false;
  bool convToTIF = false;
  bool saveHDF = false;
  int doHitFind = 0;
  float threshold = 0;
  int numpeaks = 0;
  int numconnected = 0;
  int locBg = 0;
  int inputDataType = 0;// 0 - CBF, 1 - HDF, 2(?) - TIF, 3 - RAW
  int dimsI[3] = {0,0,0};
  strcpy(fName,argv[1]);
//printf("len %d\n",strlen(strrchr(fName,'.')));
//printf("argv1 \"%s\"\n",argv[1]);
//+  if (argc>2) if (argv[2][0]=='r' && argv[2][1]=='a' && argv[2][2]=='w') convToRAW = true;
  for (int i=2; i<argc; i++)
  { if (strncmp(argv[i],"raw",3)==0) convToRAW = true;
    if (strncmp(argv[i],"tif",3)==0) convToTIF = true;
    if (strncmp(argv[i],"comp",4)==0) compress = true;
    if (strncmp(argv[i],"save",4)==0) saveHDF = true;
    if (strncmp(argv[i],"threshold",9)==0) threshold = atof(strchr(argv[i],'=')+1);
    if (strncmp(argv[i],"numpeaks",8)==0) numpeaks = atoi(strchr(argv[i],'=')+1);
    if (strncmp(argv[i],"numconnected",11)==0) numconnected = atoi(strchr(argv[i],'=')+1);
    if (strncmp(argv[i],"locbgrad",8)==0) locBg = atoi(strchr(argv[i],'=')+1);
    if (strncmp(argv[i],"datafield",9)==0) strcpy(datafield,(strchr(argv[i],'=')+1));
    if (strncmp(argv[i],"dims",4)==0)
    { char* initpo = strchr(argv[i],'=')+1;
      int ci = 0;
      while (strlen(initpo)>0)
      { char* finpo = strchr(initpo,',');
        if (finpo!=NULL) finpo[0] = 0;
        dimsI[ci] = atof(initpo);
        if (finpo!=NULL) initpo=finpo+1;
        else initpo[0] = 0;
        ci++;
      }
      printf("dims=%d,%d,%d\n",dimsI[0],dimsI[1],dimsI[2]);
    }
  }
  if (convToRAW) printf("Converting to RAW files\n");
  else if (convToTIF) printf("Converting to TIF (32 bit) files\n");
  else if (compress) printf("Saving compressed H5 files\n");

  if (threshold>0 && numpeaks>0)
  { doHitFind = 1;
    printf("Doing hit-finding with threshold of %0.2f and minimul %d peaks\n", threshold, numpeaks);
  } else
    saveHDF = true;

#ifndef _HDF5_H
  if ((!convToRAW && !convToTIF) && doHitFind!=1)
  { printf("HDF5 library not found - only conversion to .raw is possible.\n");
    printf("  Run the program in the form: cbf2hdf5 filelist raw\n");
    return 0;
  }
#endif
  char str1[MaxStrLen];

  bool convSingleF = false;

  if (strrchr(fName,'.')!=NULL)
  { strcpy(str1,strrchr(fName,'.')+1);
    if (tolower(str1[0])=='c' && tolower(str1[1])=='b' && tolower(str1[2])=='f')
      convSingleF = true;
    if (tolower(str1[0])=='h' && tolower(str1[1])=='5')
      convSingleF = true;
  }

  strcpy(str1,(strrchr(fName,'/')!=NULL?strrchr(fName,'/')+1:(strrchr(fName,'\\')!=NULL?strrchr(fName,'\\')+1:fName)));
  if (convSingleF) printf("converting a file %s to %s\n",str1,(convToRAW?".raw":(convToTIF?".tif":".h5")));
  else printf("converting a list of files %s to %s\n",str1,(convToRAW?".raw":(convToTIF?".tif":".h5")));
  FILE* stList;
  if (!convSingleF)
  { stList = fopen(fName,"rt");
    if (stList==NULL)
    { printf("List of files %s not found!\n",fName);
      return 0;
    }
  }

//printf("1\n");

  size_t numFi = 0;
  if (!convSingleF)
  { while (!feof(stList))
    { fgets(str1,MaxStrLen,stList);
      numFi++;
    }
    rewind(stList);
  }
  else
    numFi=1;

  char* filesList= new char[MaxStrLen*numFi];
  if (!convSingleF)
  { numFi=0;
    while (!feof(stList))
    { fgets(str1,MaxStrLen,stList);
      if (feof(stList)) break;
      sscanf(str1,"%s",str1);
      sprintf(&filesList[numFi*MaxStrLen],"%s\0",str1);
      numFi++;
    }
    printf("Found %ld files in the list %s\n",numFi,fName);
    fclose(stList);
  } else
    strcpy(filesList,fName);

//  int numHitPatterns = 0;

  int* numHit = new int[numFi];
  float* averIntHit = new float[numFi];
  for (size_t ifi=0; ifi<numFi; ifi++)
  { numHit[ifi] = 0;
    averIntHit[ifi] = 0;
  }
  int numProc0th = 0;



//140404 - OMP disabled for now
#undef _OPENMP

#ifdef _OPENMP
  printf ("Found %d threads\n", omp_get_max_threads());
//?  omp_set_dynamic(0);
//    omp_set_num_threads(NUM_THREADS);
//?  omp_init_lock(&lock);
#endif

//  while (contRead)
#ifdef _OPENMP
//130715
#pragma omp parallel for
#endif

  for (size_t ifi=0; ifi<numFi; ifi++)
  {
    dtype* outArray1 = NULL;
    unsigned short *badMask;
    char fileoutname[MaxStrLen];
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
    if (strncmp(extention,"cbf",3)==0) inputDataType = 0;
    else if (strncmp(extention,"h5",2)==0 || strncmp(extention,"hdf",3)==0) inputDataType = 1;
    else if (strncmp(extention,"tif",3)==0) inputDataType = 2;
    else if (strncmp(extention,"raw",3)==0 || strncmp(extention,"bin",3)==0 ||
             strncmp(extention,"dat",3)==0 || strlen(extention)==0)
      inputDataType = 3;
    else
    { printf("File extension of %s is not recognized, skipping\n", fName1);
      continue;
    }

    float pixel[2] = {0,0};
    float expo = 0;
    float waveLen = 0;
    float dist = 0;
    float beamxy[2] = {0,0};
    float flux = 0;
    int dims[3] = {0,0,0};

    bool readOK = true;
    if (inputDataType==0) readOK = ReadCBFfile(fName1, &outArray1, dims, pixel, &expo, &waveLen, &dist, beamxy, &flux, &badMask);
    else if (inputDataType==1) readOK = ReadHDF5file(fName1, datafield, &outArray1, dims);
    else if (inputDataType==2) readOK = TIFFloatReader(fName1, &outArray1, &dims[0], &dims[1], comment);
    else if (inputDataType==3)
    {
      if (!(dimsI[0]>0 && dimsI[1]>0))
      { printf("For raw files you must set dimentions dims=x,y,z!\n");
        continue;
      }

      int fulldims = 1;
      for (int i=0; i<3; i++)
        if (dimsI[i]>0) fulldims *= dimsI[i];
      int numEl = fulldims;

      for (int i=0; i<3; i++) dims[i]=dimsI[i];

      readOK = RAWFloatReader(fName1, &outArray1, &numEl);
      if (numEl != fulldims) printf("Strange, read %d elements instead of expected %d\n",numEl,fulldims);
    }
    else readOK = false;
    if (!readOK || dims[0]<2 || dims[1]<1)
    { printf("The file %s couldn't be read!\n", fName1);
      continue;
    }

    if (locBg>0)
//  { int box[2] = {2*locBg+1, 2*locBg+1};
    { int box[2] = {locBg, locBg};
      LocBackSub(outArray1, dims, box);
    }

    if (doHitFind>0)
    { //int numFoundPeaks = 0;
      //float averint = 0;
//      if (ifi % 100 == 0) printf("%d / %ld patterns with hits\n",numHitPatterns,ifi);
      numHit[ifi] = HitFinder(outArray1, dims, threshold, numconnected, &averIntHit[ifi]);
//      if (numFoundPeaks >= numpeaks) isHit[ifi] = true;
//131213      if (numFoundPeaks < numpeaks) continue;
//131213      fprintf(listOfHits,"%s\n",fName1);
//131213      fprintf(listOfHitsInfo,"%s\t%d\t%0.2f\n",fName1,numFoundPeaks,averint);
//131213      numHitPatterns++;
    }

    if ((saveHDF && doHitFind==0) || (saveHDF && doHitFind>0 && numHit[ifi]>=numpeaks))
    { dist *= 1000; // in mm - needed for crystfel
      if (convToRAW)
      { //Save the RAW - make a choise through command line
        sprintf(fileoutname, "%s_%dx%d.raw",justName,dims[0],dims[1]);
//      sprintf(fileoutname, "%s_%dx%d.raw",ExtractFileNameNoExtC(ExtractFileNameC(fName1)),dims[0],dims[1]);
        FILE* outRAW = fopen(fileoutname,"wb");
        for (size_t i=0; i<dims[0]*dims[1]; i++)
          fwrite(&outArray1[i],sizeof(dtype),1,outRAW);
        fclose(outRAW);
      } else if (convToTIF)
      { char str1[MaxStrLen];
        sprintf(str1, "Pixel=%0.2fx%0.2fum, expo=%0.4fs, dist=%0.2fmm, Lambda=%0.4fA, beamXY=%0.1fx%0.1fpx",pixel[0]*1e6,pixel[1]*1e6,expo,dist*1000,waveLen,beamxy[0],beamxy[1]);
        sprintf(fileoutname, "%s.tif",justName);
        TIFFloatWriter(fileoutname, outArray1, dims[0], dims[1], str1);
      } else
      { //Here save an HDF5
        sprintf(fileoutname, "%s.h5", justName);
//      sprintf(fileoutname, "%s.h5",ExtractFileNameNoExtC(ExtractFileNameC(fName1)));
//#ifdef _OPENMP
//      omp_set_lock(&lock);
//#endif
#ifdef _HDF5_H
        WriteHDF5file(fileoutname, outArray1, dims, pixel, &expo, &waveLen, &dist, beamxy, &flux, compress, badMask);
#endif
//#ifdef _OPENMP
//      omp_unset_lock(&lock);
//#endif
      }
    }
    delete outArray1;

#ifdef _OPENMP
    if (omp_get_thread_num()==0)
    { numProc0th++;
      if (numProc0th % 10 == 0)
        printf("  %d patterns found in 1 thread\n",numProc0th);//,(strncmp(SNum2, "Cell parameters",15)==0?" indexed":""));
    }
//+    printf("(debug)Thread %d executes file %ld\n", omp_get_thread_num(), ifi);
#endif
#ifndef _OPENMP
    if (ifi % 100 == 0)
      printf("  %d patterns found\n",ifi);//,(strncmp(SNum2, "Cell parameters",15)==0?" indexed":""));
#endif

  }
  if (doHitFind>0)
  { int totNumHits = 0;
    char fName1[MaxStrLen];
    FILE* listOfHits = fopen("listofhits.txt","wt");
    FILE* listOfHitsInfo = fopen("listofpeaksinfo.txt","wt");
    for (size_t ifi=0; ifi<numFi; ifi++)
      if (numHit[ifi]>=numpeaks)
      {
        strncpy(fName1,&filesList[ifi*MaxStrLen],MaxStrLen);
        fprintf(listOfHits,"%s\n",fName1);
        fprintf(listOfHitsInfo,"%s\t%d\t%0.2f\n",fName1,numHit[ifi],averIntHit[ifi]);
        totNumHits++;
      }
    printf("Total number of patterns with hits: %d\n",totNumHits);
    fclose(listOfHits);
    fclose(listOfHitsInfo);
  } else
    printf("Converted %ld files.\n",numFi);
  return 1;
}
//---------------------------------------------------------------------------
