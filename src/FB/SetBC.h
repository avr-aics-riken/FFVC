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

/** 
 * @file   SetBC.h
 * @brief  FlowBase SetBC class Header
 * @author kero
 */

#include "cpm_Define.h"
#include "cpm_ParaManager.h"

#include "DomainInfo.h"
#include "FB_Define.h"
#include "BndOuter.h"
#include "Control.h"
#include "Component.h"
#include "Medium.h"
#include "Intrinsic.h"


class SetBC : public DomainInfo {
protected:
  int *ix, *jx, *kx, *gc; ///< Fortranの引数用のポインタ
  
  REAL_TYPE accel, Dp1, Dp2, mach, BasePrs;
  REAL_TYPE RefV, RefL, DiffTemp, BaseTemp, Peclet, Reynolds, rei, pei;
  REAL_TYPE Lbx[3], Rayleigh, Grashof, Prandtl;
  REAL_TYPE rho, nyu, cp, lambda, beta;
  
  int Example, NoBC, Unit_Temp, Unit_Prs;
  bool     inout_flag, isCDS;
  
  BoundaryOuter   obc[NOFACE];
  CompoList       *cmp;
  MediumList      *mat;
  Intrinsic       *Ex;
  
  /** 外部境界の種類 */
  enum obc_kind 
  {
    id_specvel = 0,
    id_outflow,
    id_wall,
    id_symmetric,
    id_periodic
  };

  
public:
  /** コンストラクタ */
  SetBC() {
    ix = jx = kx = gc = NULL;
    rei = accel = Dp1 = Dp2 = mach = RefV = RefL = DiffTemp = BaseTemp = pei = 0.0;
    rho = nyu = cp = lambda = beta = BasePrs = 0.0;
    Peclet = Reynolds = Rayleigh = Grashof = Prandtl = 0.0;
    Example = NoBC = Unit_Temp = Unit_Prs = 0;
    inout_flag = isCDS = false;
    
    cmp = NULL;
    mat = NULL;
    Ex  = NULL;
  }
  
  /**　デストラクタ */
  virtual ~SetBC() {}
  
protected:
  /**
   * @brief 外部境界処理用のループインデクスを取得する
   * @param [in] face 外部境界面番号
   * @param [out] st  開始インデクス
   * @param [out] ed  終了インデクス
   */
  void getOuterLoopIdx(const int face, int* st, int* ed);
  

public:
  /**
   @brief 静止座標系のときの流出速度制御の値を計算する
   @param [in] tm 時刻
   @retval 流出境界速度
   @todo experimental
   */
	REAL_TYPE getVrefOut (const REAL_TYPE tm);
  
  
  /** 
   * @brief クラスに必要な変数のコピー
   * @param [in] Cref       Controlクラス
   * @param [in] mat        MediumListクラス
   * @param [in] cmp        Componentクラス
   * @param [in] RF         ReferenceFrameクラス
   * @param [in] ExRef      Intrinsicクラス
   */
  void setControlVars(Control* Cref, MediumList* mat, CompoList* cmp, ReferenceFrame* RF, Intrinsic* ExRef=NULL);
  
  
  /** 
   * @brief クラスのポインタコピー
   * @param [in] m_CMP        CompoListクラス
   * @param [in] m_MAT        MediumListクラス
   */
  void importCMP_MAT(CompoList* m_CMP, MediumList* m_MAT);
  
  
  //@brief 流入出境界条件が設定されている場合にtrueを返す
  bool has_InOut(void) 
  { 
    return inout_flag; 
  }
  
  
  /** 外部境界リストのポインタを返す */
  BoundaryOuter* export_OBC()
  { 
    return obc;
  }
  
  
  /** 引数の外部境界面の外部境界リストのポインタを返す
   * @param [in] face 面番号
   */
  BoundaryOuter* export_OBC(const int face)
  { 
    return &obc[face];
  }
  
};

#endif // _FB_SETBC_H_
