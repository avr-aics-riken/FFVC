#ifndef _SKL_FB_DT_CNTL_H_
#define _SKL_FB_DT_CNTL_H_

/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file DTcntl.h
//@brief FlowBase DTcntl class Header
//@author keno, FSI Team, VCAD, RIKEN

#include <string>

#include "FBDefine.h"
#include "stdio.h"

class DTcntl {//: public Parallel_Node {
public:
  enum dt_Type {
    dt_direct=1,      // 入力値がΔt
    dt_cfl_ref_v,     // dt < c dx/U0
    dt_cfl_max_v,     // dt < c dx/Umax
    dt_dfn,           // 拡散数制限
    dt_cfl_dfn_ref_v, // dt = min( c dx/U0, diffusion number )
    dt_cfl_dfn_max_v, // dt = min( c dx/Umax, diffusion number )
    dt_cfl_max_v_cp   // 圧縮性　dt < (cfl+soundSpeed) dx/Umax
  };
  
private:
  unsigned scheme;  // Δtのスキーム種別
  unsigned KOS;     // Kind of Solver
  unsigned mode;    // 入力パラメータの次元モード（無次元/有次元）
  double   CFL;     // Δtの決定に使うCFLなど
  double   deltaT;  // Δt（無次元）
  double   dh;      // 格子幅（無次元）
  double   Reynolds;// レイノルズ数
  double   Peclet;  // ペクレ数
  
public:
  DTcntl() {
    scheme = 0;
    CFL = 0.0;
    KOS = 0;
    mode = 0;
    deltaT = 0.0;
    dh = Reynolds = Peclet = 0.0;
  }
  ~DTcntl() {}
   
public:
  //@fn unsigned get_Scheme(void)
  unsigned get_Scheme(void) const { return scheme; };
  
  //@fn double get_dt(void)
  double get_DT(void) const { return deltaT; };
  
  //@fn Sdouble get_CFL(void)
  double get_CFL(void) const { return CFL; };

  double dtCFL  (const double Uref) const;
  double dtDFN  (const double coef) const;
  
  bool chkDtSelect(void) const;
  bool set_Scheme (const char* str, const double val);
  
  unsigned set_DT (const double vRef);
  void set_Vars   (const unsigned m_kos, const unsigned m_mode, const double m_dh, const double re, const double pe);
};

#endif // _SKL_FB_DT_CNTL_H_
