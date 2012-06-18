#ifndef _FB_SETBC_H_
#define _FB_SETBC_H_

// #################################################################
//
// CAERU Library
//
// Copyright (c) All right reserved. 2012
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
//
// #################################################################

//@file SetBC.h
//@brief FlowBase SetBC class Header
//@author keno, FSI Team, VCAD, RIKEN

#include "Skl.h"
#include "SklSolverBase.h"
#include "FB_Define.h"
#include "FBUtility.h"
#include "BndOuter.h"
#include "Control.h"
#include "Component.h"
#include "Medium.h"
#include "SklUtil.h"
#include "Intrinsic.h"
#include "Parallel_node.h"

class SetBC : public Parallel_Node {
protected:
  int *ix, *jx, *kx, *gc;
  int dim_sz[3], ixc, jxc, kxc;
  
  REAL_TYPE dh, accel, Dp1, Dp2, mach, BasePrs;
  REAL_TYPE RefV, RefL, DiffTemp, BaseTemp, Peclet, Reynolds, rei, pei;
  REAL_TYPE Lbx[3], Rayleigh, Grashof, Prandtl;
  REAL_TYPE rho, nyu, cp, lambda, beta;
  
  unsigned imax, jmax, kmax, guide, size[3];
  unsigned Example, NoBC, Unit_Temp, Unit_Prs;
  bool     inout_flag, isCDS;
  
  BoundaryOuter   obc[NOFACE];
  CompoList*      cmp;
  MediumList*     mat;
  Intrinsic*      Ex;
  
  void getOuterLoopIdx (int face, int* st, int* ed);
  
  /// 外部境界の種類
  enum obc_kind {
    id_specvel = 0,
    id_outflow,
    id_wall,
    id_symmetric,
    id_periodic
  };
  
public:
  SetBC() {
    imax = jmax = kmax = guide = 0;
    ix = jx = kx = gc = NULL;
    ixc = jxc = kxc = 0;
    dh = rei = accel = Dp1 = Dp2 = mach = RefV = RefL = DiffTemp = BaseTemp = pei = 0.0;
    rho = nyu = cp = lambda = beta = BasePrs = 0.0;
    Peclet = Reynolds = Rayleigh = Grashof = Prandtl = 0.0;
    Example = NoBC = Unit_Temp = Unit_Prs = 0;
    inout_flag = isCDS = false;
    
    for (int i=0; i<3; i++) {
      dim_sz[i] = 0;
      size[i]   = 0.0;
      Lbx[i]    = 0.0;
    }
    
    cmp = NULL;
    mat = NULL;
    Ex  = NULL;
  }
  virtual ~SetBC() {}

public:
	REAL_TYPE getVrefOut (REAL_TYPE tm);
  
  void setControlVars (Control* Cref, MediumList* mat, CompoList* cmp, ReferenceFrame* RF, Intrinsic* ExRef=NULL);
  void setWorkList    (CompoList* m_CMP, MediumList* m_MAT);
  
  //@fn bool has_InOut(void)
  //@brief 流入出境界条件が設定されている場合にtrueを返す
  bool has_InOut(void) { return inout_flag; }
  
  //@fn BoundaryOuter* get_OBC_Ptr(void)
  //@brief 外部境界リストのポインタを返す
  BoundaryOuter* get_OBC_Ptr(void) { 
    return obc;
  }
  
  //@fn BoundaryOuter* get_OBC_Ptr(int face)
  //@brief 引数の外部境界面の外部境界リストのポインタを返す
  //@param face 面番号
  BoundaryOuter* get_OBC_Ptr(int face) { 
    return &obc[face];
  }
  
};

#endif // _FB_SETBC_H_
