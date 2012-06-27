/*
 * CPM - Cartesian Partition Manager
 *
 * Copyright (C) 2012 University of Tokyo.
 *
 */

/**
 * @file   cpm_ParaManager_Alloc.cpp
 * パラレルマネージャクラスのソースファイル
 * @author University of Tokyo
 * @date   2012/05/31
 */
#include <stdlib.h>
#include "cpm_ParaManager.h"

////////////////////////////////////////////////////////////////////////////////
// 配列確保
REAL_TYPE*
cpm_ParaManager::AllocRealS4D( int nmax, int vc, int procGrpNo )
{
  const int *sz = GetLocalVoxelSize( procGrpNo );
  if( !sz ) return NULL;

  size_t nw = size_t(sz[0]+2*vc) * size_t(sz[1]+2*vc) * size_t(sz[2]+2*vc) * size_t(nmax);
  return new REAL_TYPE[nw];
}

////////////////////////////////////////////////////////////////////////////////
// 配列確保
double*
cpm_ParaManager::AllocDoubleS4D( int nmax, int vc, int procGrpNo )
{
  const int *sz = GetLocalVoxelSize( procGrpNo );
  if( !sz ) return NULL;

  size_t nw = size_t(sz[0]+2*vc) * size_t(sz[1]+2*vc) * size_t(sz[2]+2*vc) * size_t(nmax);
  return new double[nw];
}

////////////////////////////////////////////////////////////////////////////////
// 配列確保
float*
cpm_ParaManager::AllocFloatS4D( int nmax, int vc, int procGrpNo )
{
  const int *sz = GetLocalVoxelSize( procGrpNo );
  if( !sz ) return NULL;

  size_t nw = size_t(sz[0]+2*vc) * size_t(sz[1]+2*vc) * size_t(sz[2]+2*vc) * size_t(nmax);
  return new float[nw];
}

////////////////////////////////////////////////////////////////////////////////
// 配列確保
int*
cpm_ParaManager::AllocIntS4D( int nmax, int vc, int procGrpNo )
{
  const int *sz = GetLocalVoxelSize( procGrpNo );
  if( !sz ) return NULL;

  size_t nw = size_t(sz[0]+2*vc) * size_t(sz[1]+2*vc) * size_t(sz[2]+2*vc) * size_t(nmax);
  return new int[nw];
}

////////////////////////////////////////////////////////////////////////////////
// 配列確保
REAL_TYPE*
cpm_ParaManager::AllocRealS3D( int vc, int procGrpNo )
{
  return AllocRealS4D(1, vc, procGrpNo);
}

////////////////////////////////////////////////////////////////////////////////
// 配列確保
double*
cpm_ParaManager::AllocDoubleS3D( int vc, int procGrpNo )
{
  return AllocDoubleS4D(1, vc, procGrpNo);
}

////////////////////////////////////////////////////////////////////////////////
// 配列確保
float*
cpm_ParaManager::AllocFloatS3D( int vc, int procGrpNo )
{
  return AllocFloatS4D(1, vc, procGrpNo);
}

////////////////////////////////////////////////////////////////////////////////
// 配列確保
int*
cpm_ParaManager::AllocIntS3D( int vc, int procGrpNo )
{
  return AllocIntS4D(1, vc, procGrpNo);
}

////////////////////////////////////////////////////////////////////////////////
// 配列確保
REAL_TYPE*
cpm_ParaManager::AllocRealV3D( int vc, int procGrpNo )
{
  return AllocRealS4D(3, vc, procGrpNo);
}

////////////////////////////////////////////////////////////////////////////////
// 配列確保
double*
cpm_ParaManager::AllocDoubleV3D( int vc, int procGrpNo )
{
  return AllocDoubleS4D(3, vc, procGrpNo);
}

////////////////////////////////////////////////////////////////////////////////
// 配列確保
float*
cpm_ParaManager::AllocFloatV3D( int vc, int procGrpNo )
{
  return AllocFloatS4D(3, vc, procGrpNo);
}

////////////////////////////////////////////////////////////////////////////////
// 配列確保
int*
cpm_ParaManager::AllocIntV3D( int vc, int procGrpNo )
{
  return AllocIntS4D(3, vc, procGrpNo);
}

////////////////////////////////////////////////////////////////////////////////
// 配列確保
REAL_TYPE*
cpm_ParaManager::AllocRealV3DEx( int vc, int procGrpNo )
{
  return AllocRealS4D(3, vc, procGrpNo);
}

////////////////////////////////////////////////////////////////////////////////
// 配列確保
double*
cpm_ParaManager::AllocDoubleV3DEx( int vc, int procGrpNo )
{
  return AllocDoubleS4D(3, vc, procGrpNo);
}

////////////////////////////////////////////////////////////////////////////////
// 配列確保
float*
cpm_ParaManager::AllocFloatV3DEx( int vc, int procGrpNo )
{
  return AllocFloatS4D(3, vc, procGrpNo);
}

////////////////////////////////////////////////////////////////////////////////
// 配列確保
int*
cpm_ParaManager::AllocIntV3DEx( int vc, int procGrpNo )
{
  return AllocIntS4D(3, vc, procGrpNo);
}

////////////////////////////////////////////////////////////////////////////////
// 配列確保
REAL_TYPE*
cpm_ParaManager::AllocRealS4DEx( int nmax, int vc, int procGrpNo )
{
  return AllocRealS4D(nmax, vc, procGrpNo);
}

////////////////////////////////////////////////////////////////////////////////
// 配列確保
double*
cpm_ParaManager::AllocDoubleS4DEx( int nmax, int vc, int procGrpNo )
{
  return AllocDoubleS4D(nmax, vc, procGrpNo);
}

////////////////////////////////////////////////////////////////////////////////
// 配列確保
float*
cpm_ParaManager::AllocFloatS4DEx( int nmax, int vc, int procGrpNo )
{
  return AllocFloatS4D(nmax, vc, procGrpNo);
}

////////////////////////////////////////////////////////////////////////////////
// 配列確保
int*
cpm_ParaManager::AllocIntS4DEx( int nmax, int vc, int procGrpNo )
{
  return AllocIntS4D(nmax, vc, procGrpNo);
}

