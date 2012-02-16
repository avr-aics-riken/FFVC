/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

/**
 @file SklSolverCBC.h
 @brief SklSolverCBC class Header
 @author keno, FSI Team, VCAD, RIKEN
 
 @note 
 - 有効セル数の変数をC++からFortranに渡すとき，unsignedから直接SKL_REAL*にキャスト(SKL_REAL*)&varとすると値がおかしくなる．
 (int*)&varでintで渡して，fortran内で型変換して利用する．
 */

#ifndef _SKL_SOLVER_CBC_CLASS_H_
#define _SKL_SOLVER_CBC_CLASS_H_

#include "Skl.h"
#include "SklSolverBase.h"
#include "SklSolverCBCDefine.h"
#include "FortranFuncCBC.h"
#include "SklTiming.h"
#include "FB_Ffunc.h"
#include "FBDefine.h"
#include "Control.h"
#include "ParseBC.h"
#include "ParseMat.h"
#include "FBUtility.h"
#include "VoxInfo.h"
#include "SklUtil.h"
#include "IP_Users.h"
#include "IterationCtl.h"
#include "History.h"
#include "CBC_SetBC.h"
#include "Parallel_node.h"
#include "Alloc.h"
#include "FileIO.h"
#include "Monitor.h"
#include "PerfMonitor.h"
#include "DTcntl.h"

#include "IP_Duct.h"
#include "IP_PPLT2D.h"
#include "IP_SHC1D.h"
#include "IP_PMT.h"
#include "IP_Rect.h"
#include "IP_Step.h"
#include "IP_Cylinder.h"
#include "IP_Polygon.h"

// K用のプロファイラ
#ifdef __K_FPCOLL
#include "fjcoll.h"
#endif

// Polylib
#include "Polylib.h"
#include "MPIPolylib.h"
#include "file_io/PolylibConfig.h"
#include "file_io/TriMeshIO.h"
using namespace PolylibNS;

// Cutlib
#include "Cutlib.h"
using namespace cutlib;

class SklSolverCBC : public SklSolverBase {

protected:
  int* GC_bv; // グローバルなコンポーネントBV
  
public:
  // Timing
  int cm_mode;
  
  // Global variables
  int para_key;
  unsigned G_Fcell, G_Acell, G_Wcell;
  unsigned G_size[3];
  SKL_REAL G_Lbx[3], G_org[3];
  
  Control         C;
  SetBC3D         BC;
  ItrCtl          IC[ItrCtl::ic_END];
  Alloc           A;
  FileIO          F;
  //Core_Utility    CU;
  MonitorList     MO;
  PerfMonitor     PM;
  DTcntl          DT;
  ReferenceFrame  RF;
  
  History*				H;
  Intrinsic*      Ex; // pointer to a base class
  CompoList*      cmp;
  MaterialList*   mat;

  // Polylib
  MPIPolylib*     PL;
  POLYLIB_STAT	  poly_stat;

  // Cutlib
  CutPos32Array* cutPos;
  CutBid8Array*  cutBid;
  
  // for parallel
  Parallel_Info pn;
  
  char *m_condition;  // for file name of text output
  char *m_log;
  
  // These pointers are to be used for arguments from C/C++ to fortran
  int *ix, ixc;
  int *jx, jxc;
  int *kx, kxc;
  int *gc;
  SKL_REAL *dh;
  SKL_REAL *dh0;
  SKL_REAL *x0;
  SKL_REAL *y0;
  SKL_REAL *z0;
  SKL_REAL v00[4];
  
  int sz[3];
  unsigned size[3];
  unsigned guide;
  
  // 履歴情報のファイルポインタ
  FILE *mp;    // 標準出力
  FILE *fp_b;  // 基本情報
  FILE *fp_w;  // 壁面情報
  FILE *fp_c;  // コンポーネント情報
  FILE *fp_d;  // 流量収支情報
  FILE *fp_i;  // 反復履歴情報
  
  SKL_REAL checkTime;
  SKL_REAL convergence_prev, convergence_rate;
  SKL_REAL range_Ut[2], range_Yp[2];
  
  // (3, ix+guide*2, jx+guide*2, kx+guide*2)
  SklVector3DEx<SKL_REAL>   *dc_v;
  SklVector3DEx<SKL_REAL>   *dc_vc;
  SklVector3DEx<SKL_REAL>   *dc_v0;
  SklVector3DEx<SKL_REAL>   *dc_wv;
  SklVector3DEx<SKL_REAL>   *dc_abf;
  SklVector3DEx<SKL_REAL>   *dc_vf0;
  
  // (3, ix, jx, kx)
  SklVector3DEx<SKL_REAL>   *dc_av;
  
  // (3, ix+guide*2, jx+guide*2, kx+guide*2)
  SklVector3DEx<SKL_REAL> *dc_wvex;
  SklVector3DEx<SKL_REAL> *dc_qbc;
  
  // (ix+guide*2, jx+guide*2, kx+guide*2)
  SklScalar3D<int>        *dc_mid;
  SklScalar3D<unsigned>   *dc_bcd;
  SklScalar3D<unsigned>   *dc_bcp;
  SklScalar3D<unsigned>   *dc_bcv;
  SklScalar3D<unsigned>   *dc_bh1;
  SklScalar3D<unsigned>   *dc_bh2;
  SklScalar3D<SKL_REAL>   *dc_ws;
  SklScalar3D<SKL_REAL>   *dc_p;
  SklScalar3D<SKL_REAL>   *dc_wk2;
  SklScalar3D<SKL_REAL>   *dc_dp;
  SklScalar3D<SKL_REAL>   *dc_p0;
  SklScalar3D<SKL_REAL>   *dc_t;
  SklScalar3D<SKL_REAL>   *dc_t0;
  SklScalar3D<SKL_REAL>   *dc_wkj;
  SklScalar3D<SKL_REAL>   *dc_vt;
  SklScalar3D<SKL_REAL>   *dc_vof;
  
  // (ix, jx, kx)
  SklScalar3D<SKL_REAL>   *dc_ap;
  SklScalar3D<SKL_REAL>   *dc_at;
  
  // out file object
  SklVoxDataSet *m_outPrs;
  SklVoxDataSet *m_outUVW;
  SklVoxDataSet *m_outTmp;
  SklVoxDataSet *m_outAvrPrs;
  SklVoxDataSet *m_outAvrUVW;
  SklVoxDataSet *m_outAvrTmp;
  SklVoxDataSet *m_outVrt;
  SklVoxDataSet *m_outTP;
  SklVoxDataSet *m_outVOF;
  SklVoxDataSet *m_outI2VGT;
  SklVoxDataSet *m_outHlcty;
  SklVoxDataSet *m_outDiv;

  // (6, ix+guide*2, jx+guide*2, kx+guide*2)
  float* cut; // Cutlibで確保する配列のポインタを受け取る

  // Fluid cell
  SklScalar<int> *dc_index3;      // index test; 
  SklScalar<unsigned> *dc_index;  // index test; 
  
  // for tuning
  int cf_sz[3];  //buffer size
  SKL_REAL *cf_x; //comm buffer of each dir
  SKL_REAL *cf_y;
  SKL_REAL *cf_z;
  
  // timing management
  unsigned ModeTiming;
  
  // プロファイラ用のラベル
  char tm_label_ptr[tm_END][TM_LABEL_MAX];
  
protected:
  SklSolverCBC();
  
  //@fn int getGlCompoBV_st_x(unsigned odr)
  //@brief コンポーネントのBV情報st_xを返す
  inline int getGlCompoBV_st_x(unsigned odr) {
    return ( GC_bv[6*odr+0] );
  }
  
  //@fn int getGlCompoBV_st_y(unsigned odr)
  //@brief コンポーネントのBV情報st_yを返す
  inline int getGlCompoBV_st_y(unsigned odr) {
    return ( GC_bv[6*odr+1] );
  }
  
  //@fn int getGlCompoBV_st_z(unsigned odr)
  //@brief コンポーネントのBV情報st_zを返す
  inline int getGlCompoBV_st_z(unsigned odr) {
    return ( GC_bv[6*odr+2] );
  }
  
  //@fn int getGlCompoBV_ed_x(unsigned odr)
  //@brief コンポーネントのBV情報st_xを返す
  inline int getGlCompoBV_ed_x(unsigned odr) {
    return ( GC_bv[6*odr+3] );
  }
  
  //@fn int getGlCompoBV_ed_y(unsigned odr)
  //@brief コンポーネントのBV情報ed_yを返す
  inline int getGlCompoBV_ed_y(unsigned odr) {
    return ( GC_bv[6*odr+4] );
  }
  
  //@fn int getGlCompoBV_ed_z(unsigned odr)
  //@brief コンポーネントのBV情報ed_zを返す
  inline int getGlCompoBV_ed_z(unsigned odr) {
    return ( GC_bv[6*odr+5] );
  }
  
public:
  SklSolverCBC                      (int sType);
  virtual ~SklSolverCBC             (void);
  virtual int SklSolverInitialize   (void);
  virtual int SklSolverLoop         (const unsigned int step);
  virtual bool SklSolverPost        (void);
  virtual bool SklSolverUsage       (const char* cmd);
  
  virtual void connectExample       (Control* Cref);
  virtual void getXMLExample        (Control* Cref);
  virtual void set_timing_label     (void);
  virtual void VoxelInitialize      (void);
  
  bool hasLinearSolver (unsigned L);
  
  float min_distance     (float* cut, FILE* fp);
  SKL_REAL Norm_Poisson     (ItrCtl* IC);
  
  void allocArray_AB2       (unsigned long &total);
  void allocArray_average   (unsigned long &total, FILE* fp);
  void allocArray_Collocate (unsigned long &total);
  void allocArray_forcing   (unsigned long &total);
  void allocArray_heat      (unsigned long &total);
  void allocArray_index     (unsigned long &total);
  void allocArray_index3    (unsigned long &total);
  void allocArray_interface (unsigned long &total);
  void allocArray_Jacobi    (unsigned long &total);
  void allocArray_LES       (unsigned long &total);
  void allocArray_main      (unsigned long &total);
  void allocArray_prep      (unsigned long &total, unsigned long &prep);
  void allocArray_RK        (unsigned long &total);
  
  void AverageOutput        (unsigned mode, SKL_REAL& flop);
  void Averaging_Time       (SKL_REAL& flop);
  void Averaging_Space      (SKL_REAL* avr, SKL_REAL& flop);
  void VoxEncode            (VoxInfo* Vinfo, ParseMat* M, int* mid, CutPos32Array* cutPos);
  void VoxScan              (VoxInfo* Vinfo, ParseBC* B, int* mid, FILE* fp);
  //void CN_Itr               (ItrCtl* IC);
  void DomainMonitor        (BoundaryOuter* ptr, Control* R, SKL_REAL& flop);
  void FileOutput           (unsigned mode, SKL_REAL& flop);
  void gather_DomainInfo    (void);
  void getEnlargedIndex     (int& m_st, int& m_ed, unsigned st_i, unsigned len, unsigned m_x, unsigned dir, int m_id);
  void getGlobalCmpIdx      (VoxInfo* Vinfo);
  void getLocalCmpIdx       (void);
  void IF_TRP_VOF           (void);
  void load_Restart_avr_file(FILE* fp);
  void load_Restart_file    (FILE* fp);
  void LS_Binary            (ItrCtl* IC, SKL_REAL b2);
  void LS_Planar            (ItrCtl* IC, SKL_REAL b2);
  void prepOutput           (void);
  void setBCinfo            (ParseBC* B);
  void setEnsComponent      (void);
  void setIDtables          (ParseBC* B, FILE* fp, FILE* mp);
  void setMaterialList      (ParseBC* B, ParseMat* M, FILE* mp, FILE* fp);
  void set_label            (unsigned key, char* label, PerfMonitor::Type type, bool exclusive=true);
  void set_Parallel_Info    (void);
  void setup_CutInfo4IP     (unsigned long& m_prep, unsigned long& m_total, FILE* fp);
  void setup_Polygon2CutInfo(unsigned long& m_prep, unsigned long& m_total, FILE* fp);
  void setVOF               (SKL_REAL* vof, unsigned* bx);
  void swap_ptr_SKL_REAL    (SKL_REAL* a, SKL_REAL* b);
  
  void NS_FS_E_CBC          (void);
  void NS_FS_E_CDS          (void);
  
  void PS_E_CBC             (void);
  void PS_EE_EI_CBC         (void);
  
  SKL_REAL PSOR(SKL_REAL* p, SKL_REAL* src0, SKL_REAL* src1, unsigned* bp, ItrCtl* IC, SKL_REAL& flop);
  SKL_REAL PSOR2sma_core(SKL_REAL* p, int ip, int color, SKL_REAL* src0, SKL_REAL* src1, unsigned* bp, ItrCtl* IC, SKL_REAL& flop);
  SKL_REAL count_comm_size (unsigned sz[3], unsigned guide) const;
  
  // CBC_Heat.C
	SKL_REAL ps_Diff_SM_EE    (SKL_REAL* t, SKL_REAL dt, SKL_REAL* qbc, unsigned* bh2, SKL_REAL* ws, SKL_REAL& flop);
  SKL_REAL ps_Diff_SM_PSOR  (SKL_REAL* t, SKL_REAL& b2, SKL_REAL dt, SKL_REAL* qbc, unsigned* bh2, SKL_REAL* ws, ItrCtl* IC, SKL_REAL& flop);

  void Buoyancy         (SKL_REAL* v, SKL_REAL dgr, SKL_REAL* t, unsigned* bd, SKL_REAL& flop);
  void ps_ConvectionEE  (SKL_REAL* tc, SKL_REAL dt, unsigned* bd, SKL_REAL* t0, SKL_REAL& flop);
  void ps_LS            (ItrCtl* IC);
  
  //@fn 時刻をRFクラスからv00[4]にコピーする
  //@param time 設定する時刻
  void copyV00fromRF(double m_time) {
    double g[4];
    RF.setV00(m_time);
    RF.copyV00(g);
    for (int i=0; i<4; i++) v00[i]=(SKL_REAL)g[i];
  }
  
  //@fn プロファイラのラベル取り出し
  //@param 格納番号
  inline const char* get_tm_label(unsigned key) {
    return (const char*)tm_label_ptr[key];
  }
  
  //@fn タイミング測定開始
  //@param 格納番号
  inline void TIMING_start(unsigned key) {
    TIMING__ PM.start(key);
    start_collection( get_tm_label(key) );
  }
  
  //@fn タイミング測定終了
  //@param 格納番号
  //@param[in] flopPerTask 「タスク」あたりの計算量/通信量(バイト) (ディフォルト0)
  //@param[in] iterationCount  実行「タスク」数 (ディフォルト1)
  inline void TIMING_stop(unsigned key, SKL_REAL flopPerTask=0.0, unsigned iterationCount=1) {
    stop_collection( get_tm_label(key) );
    TIMING__ PM.stop(key, flopPerTask, iterationCount);
  }

#ifndef __K_FPCOLL
  //@fn K用プロファイラのスタブ -D__K_FPCOLLがオプション指定されないと、こちらが有効
  void start_collection(const char*) {}
  void stop_collection (const char*) {}
#endif
  
};

#endif // _SKL_SOLVER_CBC_CLASS_H_
