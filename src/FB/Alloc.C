// #################################################################
//
// FFV : Frontflow / violet
//
// Copyright (c) All right reserved. 2012
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
//
// #################################################################

/** 
 * @file Alloc.C
 * @brief FlowBase Alloc class
 * @author kero
 */

#include "Alloc.h"


// データ領域をアロケートする（Scalar4:REAL_TYPE）
float* Alloc::alloc_Float_S4D(const int* sz, const int gc, const int dnum, unsigned long &mc)
{
  if ( !sz ) return NULL;
  
  size_t dims[3], nx;
  
  dims[0] = (size_t)(sz[0] + 2*gc);
  dims[1] = (size_t)(sz[1] + 2*gc); 
  dims[2] = (size_t)(sz[2] + 2*gc); 
  
  nx = dim[0] * dim[1] * dim[2] * (size_t)dnum;
  
  if ( nx > LONG_MAX ) {
    printf("Error : Allocation size overflow\n");
    return NULL;
  }
  
  float* var = new float[nx];
  
  memset(var, 0, sizeof(float)*nx);
  
  mc += (unsigned long)( nx*sizeof(float) );
  
  return var;
}


// データ領域をアロケートする（Scalar:float）
float* Alloc::alloc_Float_S3D(const int* sz, const int gc, unsigned long &mc)
{
  if ( !sz ) return NULL;
  
  size_t dims[3], nx;
  
  dims[0] = (size_t)(sz[0] + 2*gc); 
  dims[1] = (size_t)(sz[1] + 2*gc); 
  dims[2] = (size_t)(sz[2] + 2*gc); 
  
  nx = dim[0] * dim[1] * dim[2];
  
  if ( nx > LONG_MAX ) {
    printf("Error : Allocation size overflow\n");
    return NULL;
  }
  
  float* var = new float[nx];
  
  memset(var, 0, sizeof(float)*nx);
  
  mc += (unsigned long)( nx*sizeof(float) );
  
  return var;
}


// データ領域をアロケートする（Scalar:REAL_TYPE）
REAL_TYPE* Alloc::alloc_Real_S3D(const int* sz, const int gc, unsigned long &mc)
{
  if ( !sz ) return NULL;
  
  size_t dims[3], nx;
  
  dims[0] = (size_t)(sz[0] + 2*gc); 
  dims[1] = (size_t)(sz[1] + 2*gc); 
  dims[2] = (size_t)(sz[2] + 2*gc); 
  
  nx = dim[0] * dim[1] * dim[2];
  
  if ( nx > LONG_MAX ) {
    printf("Error : Allocation size overflow\n");
    return NULL;
  }
  
  REAL_TYPE* var = new REAL_TYPE[nx];
  
  memset(var, 0, sizeof(REAL_TYPE)*nx);
  
  mc += (unsigned long)( nx*sizeof(REAL_TYPE) );
  
  return var;
}


// データ領域をアロケートする（Scalar:int）
int* Alloc::alloc_Int_S3D(const int* sz, const int gc, unsigned long &mc)
{
  if ( !sz ) return NULL;
  
  size_t dims[3], nx;
  
  dims[0] = (size_t)(sz[0] + 2*gc); 
  dims[1] = (size_t)(sz[1] + 2*gc); 
  dims[2] = (size_t)(sz[2] + 2*gc); 
  
  nx = dim[0] * dim[1] * dim[2];
  
  if ( nx > LONG_MAX ) {
    printf("Error : Allocation size overflow\n");
    return NULL;
  }
  
  int* var = new int[nx];
  
  memset(var, 0, sizeof(int)*nx);

  mc += (unsigned long)( nx*sizeof(int) );
  
  return var;
}


// データ領域をアロケートする（Scalar:unsigned）
unsigned* Alloc::alloc_Uint_S3D(const int* sz, const int gc, unsigned long &mc)
{
  if ( !sz ) return NULL;
  
  size_t dims[3], nx;
  
  dims[0] = (size_t)(sz[0] + 2*gc); 
  dims[1] = (size_t)(sz[1] + 2*gc); 
  dims[2] = (size_t)(sz[2] + 2*gc); 
  
  nx = dim[0] * dim[1] * dim[2];
  
  if ( nx > LONG_MAX ) {
    printf("Error : Allocation size overflow\n");
    return NULL;
  }
  
  unsigned* var = new unsigned[nx];
  
  memset(var, 0, sizeof(unsigned)*nx);
  
  mc += (unsigned long)( nx*sizeof(unsigned) );
  
  return var;
}


// データ領域をアロケートする（Vector:REAL_TYPE）
REAL_TYPE* Alloc::alloc_Real_V3D(const int* sz, const int gc, unsigned long &mc)
{
  if ( !sz ) return NULL;
  
  size_t dims[3], nx;
  
  dims[0] = (size_t)(sz[0] + 2*gc); 
  dims[1] = (size_t)(sz[1] + 2*gc); 
  dims[2] = (size_t)(sz[2] + 2*gc); 
  
  nx = dim[0] * dim[1] * dim[2] * 3;
  
  if ( nx > LONG_MAX ) {
    printf("Error : Allocation size overflow\n");
    return NULL;
  }
  
  REAL_TYPE* var = new REAL_TYPE[nx];
  
  memset(var, 0, sizeof(REAL_TYPE)*nx);
  
  mc += (unsigned long)( nx*sizeof(REAL_TYPE) );
  
  return true;
}
