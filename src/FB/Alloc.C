//##################################################################################
//
// Flow Base class
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2015 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2016 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
// Copyright (c) 2016-2020 Research Institute for Information Research, Kyushu University.
// All rights reserved.
//
//##################################################################################


/** 
 * @file   Alloc.C
 * @brief  FlowBase Alloc class
 * @author aics
 */

#include "Alloc.h"
#include <stdio.h>


// #################################################################
// データ領域をアロケートする（Scalar:double）
double* Alloc::Double_S3D(const int* sz, const int gc)
{
  if ( !sz ) return NULL;
  
  size_t dims[3], nx;
  
  dims[0] = (size_t)(sz[0] + 2*gc); 
  dims[1] = (size_t)(sz[1] + 2*gc); 
  dims[2] = (size_t)(sz[2] + 2*gc); 
  
  nx = dims[0] * dims[1] * dims[2];
  
  double* var = new double[nx];
  
  memset(var, 0, sizeof(double)*nx);
  
  return var;
}


// #################################################################
// データ領域をアロケートする（Scalar:float）
float* Alloc::Float_S3D(const int* sz, const int gc)
{
  if ( !sz ) return NULL;
  
  size_t dims[3], nx;
  
  dims[0] = (size_t)(sz[0] + 2*gc);
  dims[1] = (size_t)(sz[1] + 2*gc);
  dims[2] = (size_t)(sz[2] + 2*gc);
  
  nx = dims[0] * dims[1] * dims[2];
  
  float* var = new float[nx];
  
  memset(var, 0, sizeof(float)*nx);
  
  return var;
}


// #################################################################
// データ領域をアロケートする（Scalar4:float）
float* Alloc::Float_S4D(const int* sz, const int gc, const int dnum)
{
  if ( !sz ) return NULL;
  
  size_t dims[3], nx;
  
  dims[0] = (size_t)(sz[0] + 2*gc);
  dims[1] = (size_t)(sz[1] + 2*gc); 
  dims[2] = (size_t)(sz[2] + 2*gc); 
  
  nx = dims[0] * dims[1] * dims[2] * (size_t)dnum;
  
  float* var = new float[nx];
  
  memset(var, 0, sizeof(float)*nx);
  
  return var;
}


// #################################################################
// データ領域をアロケートする（Scalar:int）
int* Alloc::Int_S3D(const int* sz, const int gc)
{
  if ( !sz ) return NULL;
  
  size_t dims[3], nx;
  
  dims[0] = (size_t)(sz[0] + 2*gc); 
  dims[1] = (size_t)(sz[1] + 2*gc); 
  dims[2] = (size_t)(sz[2] + 2*gc); 
  
  nx = dims[0] * dims[1] * dims[2];
  
  int* var = new int[nx];
  
  memset(var, 0, sizeof(int)*nx);
  
  return var;
}


// #################################################################
// データ領域をアロケートする（Scalar:REAL_TYPE）
REAL_TYPE* Alloc::Real_S3D(const int* sz, const int gc)
{
  if ( !sz ) return NULL;
  
  size_t dims[3], nx;
  
  dims[0] = (size_t)(sz[0] + 2*gc); 
  dims[1] = (size_t)(sz[1] + 2*gc); 
  dims[2] = (size_t)(sz[2] + 2*gc); 
  
  nx = dims[0] * dims[1] * dims[2];
  
  REAL_TYPE* var = new REAL_TYPE[nx];
  
  memset(var, 0, sizeof(REAL_TYPE)*nx);
  
  return var;
}

// #################################################################
// データ領域をアロケートする（Scalar4:REAL_TYPE）
REAL_TYPE* Alloc::Real_S4D(const int* sz, const int gc, const int dnum)
{
  if ( !sz ) return NULL;
  
  size_t dims[3], nx;
  
  dims[0] = (size_t)(sz[0] + 2*gc);
  dims[1] = (size_t)(sz[1] + 2*gc);
  dims[2] = (size_t)(sz[2] + 2*gc);
  
  nx = dims[0] * dims[1] * dims[2] * (size_t)dnum;
  
  REAL_TYPE* var = new REAL_TYPE[nx];
  
  memset(var, 0, sizeof(REAL_TYPE)*nx);
  
  return var;
}

// #################################################################
// データ領域をアロケートする（Vector:REAL_TYPE）
REAL_TYPE* Alloc::Real_V3D(const int* sz, const int gc)
{
  if ( !sz ) return NULL;
  
  size_t dims[3], nx;
  
  dims[0] = (size_t)(sz[0] + 2*gc); 
  dims[1] = (size_t)(sz[1] + 2*gc); 
  dims[2] = (size_t)(sz[2] + 2*gc); 
  
  nx = dims[0] * dims[1] * dims[2] * 3;
  
  REAL_TYPE* var = new REAL_TYPE[nx];
  
  memset(var, 0, sizeof(REAL_TYPE)*nx);
  
  return var;
}


// #################################################################
// データ領域をアロケートする（Scalar:unsigned）
unsigned* Alloc::Uint_S3D(const int* sz, const int gc)
{
  if ( !sz ) return NULL;
  
  size_t dims[3], nx;
  
  dims[0] = (size_t)(sz[0] + 2*gc); 
  dims[1] = (size_t)(sz[1] + 2*gc); 
  dims[2] = (size_t)(sz[2] + 2*gc); 
  
  nx = dims[0] * dims[1] * dims[2];
  
  unsigned* var = new unsigned[nx];
  
  memset(var, 0, sizeof(unsigned)*nx);
  
  return var;
}
