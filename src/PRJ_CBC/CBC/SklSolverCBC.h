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
 - 有効セル数の変数をC++からFortranに渡すとき，unsignedから直接REAL_TYPE*にキャスト(REAL_TYPE*)&varとすると値がおかしくなる．
 (int*)&varでintで渡して，fortran内で型変換して利用する．
 */

#ifndef _SKL_SOLVER_CBC_CLASS_H_
#define _SKL_SOLVER_CBC_CLASS_H_

#include "FBDefine.h"
#include "Skl.h"
#include "SklSolverBase.h"
#include "SklSolverCBCDefine.h"
#include "FortranFuncCBC.h"
#include "SklTiming.h"
#include "FB_Ffunc.h"
#include "Control.h"
#include "ParseBC.h"
#include "ParseMat.h"
#include "FBUtility.h"
#include "VoxInfo.h"
#include "SklUtil.h"
#include "History.h"
#include "CBC_SetBC.h"
#include "Parallel_node.h"
#include "Alloc.h"
#include "FileIO.h"
#include "Monitor.h"
#include "PerfMonitor.h"

#include "IP_Duct.h"
#include "IP_PPLT2D.h"
#include "IP_SHC1D.h"
#include "IP_PMT.h"
#include "IP_Rect.h"
#include "IP_Step.h"
#include "IP_Cylinder.h"
#include "IP_Polygon.h"
#include "IP_Sphere.h"

// K/FX10 profiler
#if defined __K_FPCOLL
#include "fjcoll.h"
#elif defined __FX_FAPP
#include "/fj_tool/fapp.h"
#endif

// Polylib
#include "Polylib.h"
#include "MPIPolylib.h"
#include "file_io/PolylibConfig.h"
#include "file_io/TriMeshIO.h"

// Cutlib
#include "Cutlib.h"

#include "omp.h"

using namespace PolylibNS;
using namespace cutlib;

class SklSolverCBC : public SklSolverBase {

protected:
  int* compo_global_bbox; // グローバルなコンポーネントBbox 表示に利用
  
public:
  // Timing
  int cm_mode;
  
  // Global variables
  int para_key;
  unsigned G_Fcell, G_Acell, G_Wcell;
  unsigned G_size[3];
  REAL_TYPE G_Lbx[3], G_org[3];
  
  Control         C;
  SetBC3D         BC;
  ItrCtl          IC[ItrCtl::ic_END];
  Alloc           A;
  FileIO          F;
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
  CutBid5Array*  cutBid;
  
  // for parallel
  Parallel_Info pn;
  
  char *m_condition;  // for file name of text output
  char *m_log;
  
  // These pointers are to be used for arguments from C/C++ to fortran
  int *ix, ixc;
  int *jx, jxc;
  int *kx, kxc;
  int *gc;
  REAL_TYPE *dh;
  REAL_TYPE *dh0;
  REAL_TYPE *x0;
  REAL_TYPE *y0;
  REAL_TYPE *z0;
  REAL_TYPE v00[4];
  
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
  
  REAL_TYPE checkTime;
  REAL_TYPE convergence_prev, convergence_rate;
  REAL_TYPE range_Ut[2], range_Yp[2];
  
  // (3, ix+guide*2, jx+guide*2, kx+guide*2)
  SklVector3DEx<REAL_TYPE> *dc_v;
  SklVector3DEx<REAL_TYPE> *dc_vc;
  SklVector3DEx<REAL_TYPE> *dc_v0;
  SklVector3DEx<REAL_TYPE> *dc_wv;
  SklVector3DEx<REAL_TYPE> *dc_abf;
  SklVector3DEx<REAL_TYPE> *dc_vf0;
  SklVector3DEx<REAL_TYPE> *dc_av;
  SklVector3DEx<REAL_TYPE> *dc_wvex;
  SklVector3DEx<REAL_TYPE> *dc_qbc;
  
  // (ix+guide*2, jx+guide*2, kx+guide*2)
  SklScalar3D<int>        *dc_mid;
  SklScalar3D<unsigned>   *dc_bcd;
  SklScalar3D<unsigned>   *dc_bcp;
  SklScalar3D<unsigned>   *dc_bcv;
  SklScalar3D<unsigned>   *dc_bh1;
  SklScalar3D<unsigned>   *dc_bh2;
  SklScalar3D<REAL_TYPE>  *dc_ws;
  SklScalar3D<REAL_TYPE>  *dc_p;
  SklScalar3D<REAL_TYPE>  *dc_wk2;
  SklScalar3D<REAL_TYPE>  *dc_dp;
  SklScalar3D<REAL_TYPE>  *dc_p0;
  SklScalar3D<REAL_TYPE>  *dc_t;
  SklScalar3D<REAL_TYPE>  *dc_t0;
  SklScalar3D<REAL_TYPE>  *dc_vt;
  SklScalar3D<REAL_TYPE>  *dc_vof;
  SklScalar3D<REAL_TYPE>  *dc_ap;
  SklScalar3D<REAL_TYPE>  *dc_at;
  SklScalar3D<float>      *dc_cvf;
  
  // コンポーネントワーク配列のアドレス管理
  REAL_TYPE** component_array;
  
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
  int* cut_id;   // カット点のID
  
  // for tuning
  int cf_sz[3];  //buffer size
  REAL_TYPE *cf_x; //comm buffer of each dir
  REAL_TYPE *cf_y;
  REAL_TYPE *cf_z;
  
  // timing management
  unsigned ModeTiming;
  
  // プロファイラ用のラベル
  char tm_label_ptr[tm_END][TM_LABEL_MAX];
  
protected:
  SklSolverCBC();
  
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
  
  bool chkMediumConsistency(void);
  bool hasLinearSolver (unsigned L);
  
  float min_distance     (float* cut, FILE* fp);
  REAL_TYPE Norm_Poisson     (ItrCtl* IC);
  
  void allocArray_AB2       (unsigned long &total);
  void allocArray_average   (unsigned long &total, FILE* fp);
  void allocArray_Collocate (unsigned long &total);
  void allocArray_compoVF   (unsigned long &prep, unsigned long &total);
  void allocArray_forcing   (unsigned long &total);
  void allocArray_heat      (unsigned long &total);
  void allocArray_interface (unsigned long &total);
  void allocArray_LES       (unsigned long &total);
  void allocArray_main      (unsigned long &total);
  void allocArray_prep      (unsigned long &prep, unsigned long &total);
  void allocArray_RK        (unsigned long &total);
  void allocComponentArray  (unsigned long& m_prep, unsigned long& m_total, FILE* fp);
  
  void AverageOutput        (unsigned mode, REAL_TYPE& flop);
  void Averaging_Time       (REAL_TYPE& flop);
  void Variation_Space      (REAL_TYPE* avr, REAL_TYPE& flop);
  void VoxEncode            (VoxInfo* Vinfo, ParseMat* M, int* mid, float* vf);
  void VoxScan              (VoxInfo* Vinfo, ParseBC* B, int* mid, FILE* fp);
  //void CN_Itr               (ItrCtl* IC);
  void DomainMonitor        (BoundaryOuter* ptr, Control* R, REAL_TYPE& flop);
  void EnlargeIndex         (int& m_st, int& m_ed, const unsigned st_i, const unsigned len, const unsigned m_x, const unsigned dir, const int m_id=0);
  void FileOutput           (unsigned mode, REAL_TYPE& flop);
  void fixed_parameters     (void);
  void gather_DomainInfo    (void);
  void getXML_Monitor       (SklSolverConfig* CF, MonitorList* M);
  void getXML_Mon_Line      (MonitorList* M, const CfgElem *elmL2, REAL_TYPE from[3], REAL_TYPE to[3], int& nDivision);
  void getXML_Mon_Pointset  (MonitorList* M, const CfgElem *elmL2, vector<MonitorCompo::MonitorPoint>& pointSet);
  void IF_TRP_VOF           (void);
  void load_Restart_avr_file(FILE* fp, REAL_TYPE& flop);
  void load_Restart_file    (FILE* fp, REAL_TYPE& flop);
  void LS_Binary            (ItrCtl* IC, REAL_TYPE b2);
  void LS_Planar            (ItrCtl* IC, REAL_TYPE b2);
  void prepOutput           (void);
  void Pressure_Shift       (REAL_TYPE* p);
  void resizeBVcell         (const int* st, const int* ed, const unsigned n, const unsigned* bx);
  void resizeBVface         (const int* st, const int* ed, const unsigned n, const unsigned* bx);
  void resizeCompoBV        (const unsigned* bd, const unsigned* bv, const unsigned* bh1, const unsigned* bh2, const unsigned kos, const bool isHeat);
  void setBbox_of_VIBC_Cut  (VoxInfo* Vinfo, const unsigned* bv);
  void setBCinfo            (ParseBC* B);
  void setComponentVF       (float* cvf);
  void setEnsComponent      (void);
  void setGlobalCmpIdx      (void);
  void setIDtables          (ParseBC* B, FILE* fp, FILE* mp);
  void setLocalCmpIdx_Binary(void);
  void setMaterialList      (ParseBC* B, ParseMat* M, FILE* mp, FILE* fp);
  void set_label            (unsigned key, char* label, PerfMonitor::Type type, bool exclusive=true);
  void set_Parallel_Info    (void);
  void setParallelism       (void);
  void setShapeMonitor      (int* mid);
  void setup_CutInfo4IP     (unsigned long& m_prep, unsigned long& m_total, FILE* fp);
  void setup_Polygon2CutInfo(unsigned long& m_prep, unsigned long& m_total, FILE* fp);
  void setVIBC_from_Cut     (VoxInfo* Vinfo, const unsigned* bv);
  void setVOF               (REAL_TYPE* vof, unsigned* bx);
  void write_distance       (float* cut);
  
  void NS_FS_E_CBC          (void);
  void NS_FS_E_CDS          (void);
  
  void PS_E_CBC             (void);
  void PS_EE_EI_CBC         (void);
  
  REAL_TYPE PSOR(REAL_TYPE* p, REAL_TYPE* src0, REAL_TYPE* src1, unsigned* bp, ItrCtl* IC, REAL_TYPE& flop);
  REAL_TYPE PSOR2sma_core(REAL_TYPE* p, int ip, int color, REAL_TYPE* src0, REAL_TYPE* src1, unsigned* bp, ItrCtl* IC, REAL_TYPE& flop);
  REAL_TYPE count_comm_size (unsigned sz[3], unsigned guide) const;
  
  // CBC_Heat.C
	REAL_TYPE ps_Diff_SM_EE    (REAL_TYPE* t, REAL_TYPE dt, REAL_TYPE* qbc, unsigned* bh2, REAL_TYPE* ws, REAL_TYPE& flop);
  REAL_TYPE ps_Diff_SM_PSOR  (REAL_TYPE* t, REAL_TYPE& b2, REAL_TYPE dt, REAL_TYPE* qbc, unsigned* bh2, REAL_TYPE* ws, ItrCtl* IC, REAL_TYPE& flop);

  void Buoyancy         (REAL_TYPE* v, REAL_TYPE dgr, REAL_TYPE* t, unsigned* bd, REAL_TYPE& flop);
  void ps_ConvectionEE  (REAL_TYPE* tc, REAL_TYPE dt, unsigned* bd, REAL_TYPE* t0, REAL_TYPE& flop);
  void ps_LS            (ItrCtl* IC);
  
  //@fn 時刻をRFクラスからv00[4]にコピーする
  //@param time 設定する時刻
  void copyV00fromRF(double m_time) {
    double g[4];
    RF.setV00(m_time);
    RF.copyV00(g);
    for (int i=0; i<4; i++) v00[i]=(REAL_TYPE)g[i];
  }
  
  //@fn プロファイラのラベル取り出し
  //@param 格納番号
  inline const char* get_tm_label(unsigned key) {
    return (const char*)tm_label_ptr[key];
  }
  
  //@fn タイミング測定開始
  //@param 格納番号
  inline void TIMING_start(unsigned key) {
    // Intrinsic profiler
    TIMING__ PM.start(key);
    
    // Venus FX profiler
#if defined __K_FPCOLL
    start_collection( get_tm_label(key) );
#elif defined __FX_FAPP
    fapp_start( get_tm_label(key), 0, 0);
#endif
  }
  
  //@fn タイミング測定終了
  //@param 格納番号
  //@param[in] flopPerTask 「タスク」あたりの計算量/通信量(バイト) (ディフォルト0)
  //@param[in] iterationCount  実行「タスク」数 (ディフォルト1)
  inline void TIMING_stop(unsigned key, REAL_TYPE flopPerTask=0.0, unsigned iterationCount=1) {
    // Venus FX profiler
#if defined __K_FPCOLL
    stop_collection( get_tm_label(key) );
#elif defined __FX_FAPP
    fapp_stop( get_tm_label(key), 0, 0);
#endif
    
    // Intrinsic profiler
    TIMING__ PM.stop(key, flopPerTask, iterationCount);
  }
  
  //@fn void normalizeCord(REAL_TYPE x[3])
  //@brief 座標値を無次元化する
  void normalizeCord(REAL_TYPE x[3]) {
    x[0] /= C.RefLength;
    x[1] /= C.RefLength;
    x[2] /= C.RefLength;
  }
  
};

#endif // _SKL_SOLVER_CBC_CLASS_H_
