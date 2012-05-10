#ifndef _SKL_FB_ALLOC_H_
#define _SKL_FB_ALLOC_H_

/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file Alloc.h
//@brief FlowBase Files class Header
//@author keno, FSI Team, VCAD, RIKEN

#include "FBDefine.h"
#include "Skl.h"
#include "SklSolverBase.h"
#include "parallel/SklParaComponent.h"
#include "FB_Ffunc.h"
#include "limits.h"


class Alloc {
  
public:
  Alloc() {}
  ~Alloc() {}
  
public:
  bool alloc_Int_S1D   (SklSolverBase* obj, SklScalar<int>*           &dc_var, const char* label, unsigned  sz, unsigned gc, 
                        int init,      unsigned long &mc, int para_key=-1, int m_procGrp=0);
  bool alloc_Int_S3D   (SklSolverBase* obj, SklScalar3D<int>*         &dc_var, const char* label, unsigned* sz, unsigned gc, 
                        int init,      unsigned long &mc, int para_key=-1, int m_procGrp=0);
	bool alloc_Real_S3D  (SklSolverBase* obj, SklScalar3D<REAL_TYPE>*   &dc_var, const char* label, unsigned* sz, unsigned gc, 
                        REAL_TYPE init, unsigned long &mc, int para_key=-1, int m_procGrp=0);
  bool alloc_Real_V3DEx(SklSolverBase* obj, SklVector3DEx<REAL_TYPE>* &dc_var, const char* label, unsigned* sz, unsigned gc, 
                        REAL_TYPE init, unsigned long &mc, int para_key=-1, int m_procGrp=0);
  bool alloc_Uint_S1D  (SklSolverBase* obj, SklScalar<unsigned>*      &dc_var, const char* label, unsigned  sz, unsigned gc, 
                        unsigned init, unsigned long &mc, int para_key=-1, int m_procGrp=0);
  bool alloc_Uint_S3D  (SklSolverBase* obj, SklScalar3D<unsigned>*    &dc_var, const char* label, unsigned* sz, unsigned gc, 
                        unsigned init, unsigned long &mc, int para_key=-1, int m_procGrp=0);
  bool alloc_Float_S3D (SklSolverBase* obj, SklScalar3D<float>* &dc_var, const char* label, unsigned* sz, unsigned gc, 
                        float init, unsigned long &mc, int para_key=-1, int m_procGrp=0);
  bool alloc_Float_S4DEx(SklSolverBase* obj, SklScalar4DEx<float>* &dc_var, const char* label, unsigned* sz, unsigned gc, unsigned dNum,
                         float init, unsigned long &mc, int para_key=-1, int m_procGrp=0);
};

#endif // _SKL_FB_ALLOC_H_
