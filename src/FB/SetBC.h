#ifndef _FB_SETBC_H_
#define _FB_SETBC_H_

//##################################################################################
//
// Flow Base class
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2013 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2013 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

/** 
 * @file   SetBC.h
 * @brief  FlowBase SetBC class Header
 * @author kero
 */

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
  
  REAL_TYPE accel, Dp1, Dp2, mach, BasePrs;
  REAL_TYPE RefV, RefL, DiffTemp, BaseTemp, Peclet, Reynolds, rei, pei;
  REAL_TYPE Lbx[3], Rayleigh, Grashof, Prandtl;
  REAL_TYPE rho, nyu, cp, lambda, beta;
  
  int NoCompo;     ///< コンポーネント数
  int NoMedium;    ///< 媒質数
  int Example;     ///<
  int Unit_Temp;   ///<
  int Unit_Prs;    ///<
  bool inout_flag; ///<
  
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
    rei = accel = Dp1 = Dp2 = mach = RefV = RefL = DiffTemp = BaseTemp = pei = 0.0;
    rho = nyu = cp = lambda = beta = BasePrs = 0.0;
    Peclet = Reynolds = Rayleigh = Grashof = Prandtl = 0.0;
    Example = Unit_Temp = Unit_Prs = NoCompo = NoMedium = 0;
    inout_flag = false;
    
    cmp = NULL;
    mat = NULL;
    Ex  = NULL;
  }
  
  /**　デストラクタ */
  virtual ~SetBC() {}
  

public:
  
  /**
   * @brief 静止座標系のときの流出速度制御の値を計算する
   * @param [in] tm 時刻
   * @retval 流出境界速度
   * @todo experimental
   */
	REAL_TYPE getVrefOut (const REAL_TYPE tm);
  
  
  /**
   * @brief クラスに必要な変数のコピー
   * @param [in,out] Cref       Controlクラス
   * @param [in]     mat        MediumListクラス
   * @param [in]     RF         ReferenceFrameクラス
   * @param [in,out] ExRef      Intrinsicクラス
   */
  void setControlVars(Control* Cref, const MediumList* mat, const ReferenceFrame* RF, Intrinsic* ExRef);
  
  
  /**
   * @brief クラスのポインタコピー
   * @param [in] m_CMP        CompoListクラス
   * @param [in] m_MAT        MediumListクラス
   */
  void importCMP_MAT(CompoList* m_CMP, MediumList* m_MAT);
  
  
  // @brief 流入出境界条件が設定されている場合にtrueを返す
  bool has_InOut() 
  { 
    return inout_flag; 
  }
  
  
  /** 
   * @brief 外部境界リストのポインタを返す 
   */
  BoundaryOuter* exportOBC()
  { 
    return obc;
  }
  
  
  /** 
   * @brief 引数の外部境界面の外部境界リストのポインタを返す
   * @param [in] face 面番号
   */
  BoundaryOuter* exportOBC(const int face)
  { 
    return &obc[face];
  }
  
};

#endif // _FB_SETBC_H_
