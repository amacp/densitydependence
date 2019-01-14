#include "WolframLibrary.h"
#include <iostream>

EXTERN_C DLLEXPORT mint WolframLibrary_getVersion(){return WolframLibraryVersion;}
EXTERN_C DLLEXPORT int WolframLibrary_initialize( WolframLibraryData libData) {return 0;}
EXTERN_C DLLEXPORT void WolframLibrary_uninitialize( WolframLibraryData libData) {}

EXTERN_C DLLEXPORT int myFunction(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res){

 int err; // error code

 MTensor m1; // input tensor
 MTensor m2; // output tensor

 mint const* dims; // dimensions of the tensor

 mreal *data1; // actual data of the input tensor
 mreal *data2; // data for the output tensor

 mint i; // bean counters
 mint j;

 m1 = MArgument_getMTensor(Args[0]);
 dims = libData->MTensor_getDimensions(m1);
 err = libData->MTensor_new(MType_Real, 2, dims,&m2);
 data1 = libData->MTensor_getRealData(m1);
 data2 = libData->MTensor_getRealData(m2);

// C code
for(i = 0; i < dims[0]; i++) {
  for(j = 0; j < dims[1]; j++) {
   data2[i*dims[1]+j] = 1;
  }
 }

 MArgument_setMTensor(Res, m2);
 return LIBRARY_NO_ERROR;

}