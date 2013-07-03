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
 * @file   SetBC.C
 * @brief  FlowBase SetBC class
 * @author kero
 */

#include <math.h>
#include "SetBC.h"

// #################################################################
/**
 * @brief クラスのポインタコピー
 * @param [in] m_CMP        CompoListクラス
 * @param [in] m_MAT        MediumListクラス
 */
void SetBC::importCMP_MAT(CompoList* m_CMP, MediumList* m_MAT)
{
  if ( !m_CMP ) Exit(0);
  cmp = m_CMP;
  
  if ( !m_MAT ) Exit(0);
  mat = m_MAT;
}


// #################################################################
/**
 * @brief クラスに必要な変数のコピー
 * @param [in] Cref       Controlクラス
 * @param [in] mat        MediumListクラス
 * @param [in] cmp        Componentクラス
 * @param [in] RF         ReferenceFrameクラス
 * @param [in] ExRef      Intrinsicクラス
 */
void SetBC::setControlVars(Control* Cref, MediumList* mat, CompoList* cmp, ReferenceFrame* RF, Intrinsic* ExRef)
{
  Reynolds  = Cref->Reynolds;
  rei       = Cref->getRcpReynolds();
  Peclet    = Cref->Peclet;
  Example   = Cref->Mode.Example;
  accel     = (REAL_TYPE)RF->getAccel();
  mach      = Cref->Mach;
  RefV      = Cref->RefVelocity;
  RefL      = Cref->RefLength;
  DiffTemp  = Cref->DiffTemp;
	BaseTemp  = Cref->BaseTemp;
	Unit_Temp = Cref->Unit.Temp;
  Unit_Prs  = Cref->Unit.Prs;
	BasePrs   = Cref->BasePrs;
  NoBC      = Cref->NoBC;
  Rayleigh  = Cref->Rayleigh;
  Grashof   = Cref->Grashof;
  Prandtl   = Cref->Prandtl;
  
  isCDS = Cref->isCDS();
  
  if ( Cref->isHeatProblem() ) pei = Cref->getRcpPeclet();

  Dp1       = Cref->Domain_p1;
  Dp2       = Cref->Domain_p2;
  
  Ex = ExRef;
  
  // get reference values >> 媒質はIDがユニークに定まる

  for (int n=Cref->NoBC+1; n<=Cref->NoCompo; n++) {
    if ( cmp[n].getMatOdr() == Cref->RefMat ) {
      
      int m = cmp[n].getMatOdr();
      
      if ( mat[m].getState() == FLUID ) {
        rho    = mat[m].P[p_density];
        nyu    = mat[m].P[p_kinematic_viscosity];
        cp     = mat[m].P[p_specific_heat];
        lambda = mat[m].P[p_thermal_conductivity];
        beta   = mat[m].P[p_vol_expansion]; // can be replaced by 1/K in the case of gas
      }
      else {
        rho    = mat[m].P[p_density];
        cp     = mat[m].P[p_specific_heat];
        lambda = mat[m].P[p_thermal_conductivity];
      }
    }
  }

}


// #################################################################
/**
 * @brief 静止座標系のときの流出速度制御の値を計算する
 * @param [in] tm 時刻
 * @retval 流出境界速度
 * @todo experimental
 */
REAL_TYPE SetBC::getVrefOut(const REAL_TYPE tm)
{
	REAL_TYPE u0;
	
	if ( accel == 0.0 ) {
		u0=0.0;
	}
	else {
    const REAL_TYPE c_pai = (REAL_TYPE)(2.0*asin(1.0));
		u0 = 0.5*(1.0-cos(2.0*c_pai*tm/accel));
		if ( tm > accel ) u0 = 0.0;
	}
	
  return u0;
}


// #################################################################
/**
 * @brief 外部境界処理用のループインデクスを取得する
 * @param [in] face 外部境界面番号
 * @param [out] st  開始インデクス
 * @param [out] ed  終了インデクス
 */
void SetBC::getOuterLoopIdx(const int face, int* st, int* ed)
{
  switch (face) 
  {
    case X_MINUS:
      st[0] = 1;         ed[0] = 1;
      st[1] = 1;         ed[1] = size[1];
      st[2] = 1;         ed[2] = size[2];
      break;
      
    case X_PLUS:
      st[0] = size[0]; ed[0] = size[0];
      st[1] = 1;         ed[1] = size[1];
      st[2] = 1;         ed[2] = size[2];
      break;
      
    case Y_MINUS:
      st[0] = 1;         ed[0] = size[0];
      st[1] = 1;         ed[1] = 1;
      st[2] = 1;         ed[2] = size[2];
      break;
      
    case Y_PLUS:
      st[0] = 1;         ed[0] = size[0];
      st[1] = size[1];   ed[1] = size[1];
      st[2] = 1;         ed[2] = size[2];
      break;
      
    case Z_MINUS:
      st[0] = 1;         ed[0] = size[0];
      st[1] = 1;         ed[1] = size[1];
      st[2] = 1;         ed[2] = 1;
      break;
      
    case Z_PLUS:
      st[0] = 1;         ed[0] = size[0];
      st[1] = 1;         ed[1] = size[1];
      st[2] = size[2];   ed[2] = size[2];
      break;
  }
}
