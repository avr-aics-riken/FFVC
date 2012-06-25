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



// データ領域をアロケートする（Scalar:float）
float* Alloc::Float_S3D(const int* sz, const int gc, float &mc)
{
  if ( !sz ) return NULL;
  
  size_t dims[3], nx;
  
  dims[0] = (size_t)(sz[0] + 2*gc); 
  dims[1] = (size_t)(sz[1] + 2*gc); 
  dims[2] = (size_t)(sz[2] + 2*gc); 
  
  nx = dims[0] * dims[1] * dims[2];
  
  if ( nx > LONG_MAX ) {
    printf("Error : Allocation size overflow\n");
    return NULL;
  }
  
  float* var = new float[nx];
  
  memset(var, 0, sizeof(float)*nx);
  
  mc += (float)( nx*sizeof(float) );
  
  return var;
}



// データ領域をアロケートする（Scalar4:REAL_TYPE）
float* Alloc::Float_S4D(const int* sz, const int gc, const int dnum, float &mc)
{
  if ( !sz ) return NULL;
  
  size_t dims[3], nx;
  
  dims[0] = (size_t)(sz[0] + 2*gc);
  dims[1] = (size_t)(sz[1] + 2*gc); 
  dims[2] = (size_t)(sz[2] + 2*gc); 
  
  nx = dims[0] * dims[1] * dims[2] * (size_t)dnum;
  
  if ( nx > LONG_MAX ) {
    printf("Error : Allocation size overflow\n");
    return NULL;
  }
  
  float* var = new float[nx];
  
  memset(var, 0, sizeof(float)*nx);
  
  mc += (float)( nx*sizeof(float) );
  
  return var;
}



// データ領域をアロケートする（Scalar:int）
int* Alloc::Int_S3D(const int* sz, const int gc, float &mc)
{
  if ( !sz ) return NULL;
  
  size_t dims[3], nx;
  
  dims[0] = (size_t)(sz[0] + 2*gc); 
  dims[1] = (size_t)(sz[1] + 2*gc); 
  dims[2] = (size_t)(sz[2] + 2*gc); 
  
  nx = dims[0] * dims[1] * dims[2];
  
  if ( nx > LONG_MAX ) {
    printf("Error : Allocation size overflow\n");
    return NULL;
  }
  
  int* var = new int[nx];
  
  memset(var, 0, sizeof(int)*nx);

  mc += (float)( nx*sizeof(int) );
  
  return var;
}


// データ領域をアロケートする（Scalar:REAL_TYPE）
REAL_TYPE* Alloc::Real_S3D(const int* sz, const int gc, float &mc)
{
  if ( !sz ) return NULL;
  
  size_t dims[3], nx;
  
  dims[0] = (size_t)(sz[0] + 2*gc); 
  dims[1] = (size_t)(sz[1] + 2*gc); 
  dims[2] = (size_t)(sz[2] + 2*gc); 
  
  nx = dims[0] * dims[1] * dims[2];
  
  if ( nx > LONG_MAX ) {
    printf("Error : Allocation size overflow\n");
    return NULL;
  }
  
  REAL_TYPE* var = new REAL_TYPE[nx];
  
  memset(var, 0, sizeof(REAL_TYPE)*nx);
  
  mc += (float)( nx*sizeof(REAL_TYPE) );
  
  return var;
}


// データ領域をアロケートする（Vector:REAL_TYPE）
REAL_TYPE* Alloc::Real_V3D(const int* sz, const int gc, float &mc)
{
  if ( !sz ) return NULL;
  
  size_t dims[3], nx;
  
  dims[0] = (size_t)(sz[0] + 2*gc); 
  dims[1] = (size_t)(sz[1] + 2*gc); 
  dims[2] = (size_t)(sz[2] + 2*gc); 
  
  nx = dims[0] * dims[1] * dims[2] * 3;
  
  if ( nx > LONG_MAX ) {
    printf("Error : Allocation size overflow\n");
    return NULL;
  }
  
  REAL_TYPE* var = new REAL_TYPE[nx];
  
  memset(var, 0, sizeof(REAL_TYPE)*nx);
  
  mc += (float)( nx*sizeof(REAL_TYPE) );
  
  return var;
}


// データ領域をアロケートする（Scalar:unsigned）
unsigned* Alloc::Uint_S3D(const int* sz, const int gc, float &mc)
{
  if ( !sz ) return NULL;
  
  size_t dims[3], nx;
  
  dims[0] = (size_t)(sz[0] + 2*gc); 
  dims[1] = (size_t)(sz[1] + 2*gc); 
  dims[2] = (size_t)(sz[2] + 2*gc); 
  
  nx = dims[0] * dims[1] * dims[2];
  
  if ( nx > LONG_MAX ) {
    printf("Error : Allocation size overflow\n");
    return NULL;
  }
  
  unsigned* var = new unsigned[nx];
  
  memset(var, 0, sizeof(unsigned)*nx);
  
  mc += (float)( nx*sizeof(unsigned) );
  
  return var;
}
