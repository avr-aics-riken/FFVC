#ifndef _SKL_FB_FILE_IO_H_
#define _SKL_FB_FILE_IO_H_

/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2011
 *
 */

//@file FileIO.h
//@brief FlowBase FileIO class Header
//@author keno, FSI Team, VCAD, RIKEN

#include "Skl.h"
#include "SklSolverBase.h"
#include "FBDefine.h"
#include "mydebug.h"
#include "SklUtil.h"
#include <fstream>
#include "Parallel_node.h"
#include "fileio/SklSbxDataSet.h"
#include "parallel/SklParaComponent.h"
#include <math.h>
#include <time.h>
#include <fcntl.h>
#include "FBUtility.h"
#include "Core_Util.h"

#ifndef _WIN32
#include <unistd.h>
#include <strings.h>
#else
#include "sph_win32_util.h"
#endif
#include <sys/types.h>

#if defined(IA32_LINUX) || defined(IA64_LINUX) || defined(SGI_ALTIX)
#include <sys/stat.h>
#endif

#ifdef MacOSX
#include <sys/uio.h>
#endif

class FileIO : public Parallel_Node {
  
public:
  FileIO() {}
  ~FileIO() {}
  
protected:
  /// sphファイルのスカラ/ベクトルの種別
  enum sv_type {
    kind_scalar=1,
    kind_vector
  };
  
  void loadSBXfile (SklSolverBase* obj, FILE* fp, const char* name, unsigned* size, unsigned guide, SklScalar3D<unsigned char>* data) const;
  void loadSBXfile (SklSolverBase* obj, FILE* fp, const char* name, unsigned* size, unsigned guide, SklScalar3D<int>* mid_data, 
                    SklScalar3D<SKL_REAL>* vol_data = NULL) const;
  
public:
  void cnv_Div           (SklScalar3D<SKL_REAL>* dst, const SklScalar3D<SKL_REAL>* src, const SKL_REAL coef, SKL_REAL& flop) const;
  void cnv_P_D2ND        (SklScalar3D<SKL_REAL>* dst, const SKL_REAL Base_prs, const SKL_REAL Ref_rho, const SKL_REAL Ref_v, const unsigned mode) const;
  void cnv_P_ND2D        (SklScalar3D<SKL_REAL>* dst, const SklScalar3D<SKL_REAL>* src, 
                          const SKL_REAL Base_prs, const SKL_REAL Ref_rho, const SKL_REAL Ref_v, const unsigned mode, SKL_REAL& flop) const ;
  void cnv_T_D2ND        (SklScalar3D<SKL_REAL>* dst, const SKL_REAL Base_tmp, const SKL_REAL Diff_tmp, const unsigned Unit) const;
  void cnv_T_ND2D        (SklScalar3D<SKL_REAL>* dst, const SklScalar3D<SKL_REAL>* src, const SKL_REAL Base_tmp, const SKL_REAL Diff_tmp, 
                          const unsigned Unit, SKL_REAL& flop) const;
  void cnv_TP_ND2D       (SklScalar3D<SKL_REAL>* dst, const SklScalar3D<SKL_REAL>* src, const SKL_REAL Ref_rho, const SKL_REAL Ref_v, SKL_REAL& flop) const;
  void cnv_V_D2ND        (SklVector3D<SKL_REAL>* dst, const SKL_REAL Ref_v) const;
  void cnv_V_D2ND        (SklVector3DEx<SKL_REAL>* dst, const SKL_REAL Ref_v) const;
  void cnv_V_ND2D        (SklVector3DEx<SKL_REAL>* dst, const SklVector3D<SKL_REAL>* src, const SKL_REAL v00[3], const SKL_REAL Ref_v, 
                          SKL_REAL& flop, unsigned stepAvr=1) const;
  void cnv_V_ND2D        (SklVector3DEx<SKL_REAL>* dst, const SklVector3DEx<SKL_REAL>* src, const SKL_REAL v00[3], const SKL_REAL Ref_v, 
                          SKL_REAL& flop, unsigned stepAvr=1) const;
  void loadSphScalar3D   (SklSolverBase* obj, FILE* fp, const char* fname, const unsigned* size, const unsigned guide, 
                          SklScalar3D<SKL_REAL>* dc_s, int& step, SKL_REAL& time, unsigned Dmode) const;
  void loadSphScalar4D   (SklSolverBase* obj, FILE* fp, const char* fname, const unsigned* size, const unsigned guide, 
                          SklScalar4D<SKL_REAL>* dc_s, int& step, SKL_REAL& time, unsigned Dmode) const;
  void loadSphVector3D   (SklSolverBase* obj, FILE* fp, const char* fname, const unsigned* size, const unsigned guide, 
                          SklVector3D<SKL_REAL>* dc_v, int& step, SKL_REAL& time, SKL_REAL *v00, unsigned Dmode) const;
  void loadSphVector3D   (SklSolverBase* obj, FILE* fp, const char* fname, const unsigned* size, const unsigned guide, 
                          SklVector3DEx<SKL_REAL>* dc_v, int& step, SKL_REAL& time, SKL_REAL *v00, unsigned Dmode) const;
  void loadSphScalar3DAvr(SklSolverBase* obj, FILE* fp, const char* fname, const unsigned* size, const unsigned guide, 
                          SklScalar3D<SKL_REAL>* dc_s, int& step, SKL_REAL& time, unsigned Dmode) const;
  void loadSphVector3DAvr(SklSolverBase* obj, FILE* fp, const char* fname, const unsigned* size, const unsigned guide, 
                          SklVector3D<SKL_REAL>* dc_v, int& step, SKL_REAL& time, SKL_REAL *v00, unsigned Dmode) const;
  void loadSphVector3DAvr(SklSolverBase* obj, FILE* fp, const char* fname, const unsigned* size, const unsigned guide, 
                          SklVector3DEx<SKL_REAL>* dc_v, int& step, SKL_REAL& time, SKL_REAL *v00, unsigned Dmode) const;
  void readSBX           (SklSolverBase* obj, FILE* fp, const char* file_attr, unsigned* size, unsigned guide, 
                          SklScalar3D<unsigned char>* dc_mid) const;
  void readSBX           (SklSolverBase* obj, FILE* fp, const char* mid_str, unsigned* size, unsigned guide, 
                          SklScalar3D<int>* dc_mid, SklScalar3D<SKL_REAL>* dc_vol = NULL) const;
  void readSVX           (SklSolverBase* obj, FILE* fp, const char* fname, unsigned* size, unsigned guide, 
                          SklScalar3D<int>* dc_mid, bool vf_mode=false, SklScalar3D<SKL_REAL>* dc_ws=NULL) const;
  void writeRawSPH       (const SKL_REAL *vf, const unsigned* size, const unsigned gc, const SKL_REAL* org, const SKL_REAL* ddx, 
                          const unsigned m_ModePrecision) const;
};
#endif // _SKL_FB_FILE_IO_H_
