/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2010
 *
 */

//@file SetBC.C
//@brief FlowBase SetBC class Header
//@author keno, FSI Team, VCAD, RIKEN

#include <math.h>
#include "SetBC.h"

/**
 @fn void SetBC::setWorklist(CompoList* m_CMP, MaterialList* m_MAT)
 @brief CompoListをポイント
 @param m_CMP
 @param m_MAT 
 */
void SetBC::setWorkList(CompoList* m_CMP, MaterialList* m_MAT)
{
  if ( !m_CMP ) assert(0);
  cmp = m_CMP;
  
  if ( !m_MAT ) assert(0);
  mat = m_MAT;
}

/**
 @fn void SetBC::setControlVars(Control* Cref, MaterialList* mat, CompoList* cmp, Intrinsic* ExRef)
 @brief 必要な値のコピー
 */
void SetBC::setControlVars(Control* Cref, MaterialList* mat, CompoList* cmp, ReferenceFrame* RF, Intrinsic* ExRef)
{
  guide          = Cref->guide;
  imax = size[0] = Cref->imax;
  jmax = size[1] = Cref->jmax;
  kmax = size[2] = Cref->kmax;
  gc = (int*)&guide;
  ix = (int*)&(Cref->imax);
  jx = (int*)&(Cref->jmax);
  kx = (int*)&(Cref->kmax);
  ixc = (int)imax;
  jxc = (int)jmax;
  kxc = (int)kmax;
  dh        = Cref->dh;
  Reynolds  = Cref->Reynolds;
  rei       = Cref->getRcpReynolds();
  Peclet    = Cref->Peclet;
  Example   = Cref->Mode.Example;
  accel     = (SKL_REAL)RF->getAccel();
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
  for (int i=0; i<3; i++) Lbx[i] = Cref->Lbx[i];
  
  if ( Cref->isHeatProblem() ) pei = Cref->getRcpPeclet();

  Dp1       = Cref->Domain_p1;
  Dp2       = Cref->Domain_p2;
  
  dim_sz[0] = (int)imax;
  dim_sz[1] = (int)jmax;
  dim_sz[2] = (int)kmax;
  
  Ex = ExRef;
  
  // get reference values >> 媒質はIDがユニークに定まる
  unsigned m;
  for (unsigned n=Cref->NoBC+1; n<=Cref->NoCompo; n++) {
    if ( cmp[n].getID() == Cref->RefID ) {
      m = cmp[n].getMatOdr();
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

/**
 @fn SKL_REAL SetBC::getVrefOut(SKL_REAL tm)
 @brief 静止座標系のときの流出速度制御の値を計算する
 @param tm 時刻
 @retval 流出境界速度
 @todo experimental
 */
SKL_REAL SetBC::getVrefOut(SKL_REAL tm)
{
	SKL_REAL u0;
	
	if ( accel == 0.0 ) {
		u0=0.0;
	}
	else {
    const SKL_REAL c_pai = (SKL_REAL)(2.0*asin(1.0));
		u0 = 0.5*(1.0-cos(2.0*c_pai*tm/accel));
		if ( tm > accel ) u0 = 0.0;
	}
	
  return u0;
}

/**
 @fn void SetBC::getOuterLoopIdx(int face, int* st, int* ed)
 @brief 外部境界処理用のループインデクスを取得する
 @param face 外部境界面番号
 @param st 開始インデクス
 @param ed 終了インデクス
 */
void SetBC::getOuterLoopIdx(int face, int* st, int* ed)
{
  switch (face) {
    case X_MINUS:
      st[0] = 1;         ed[0] = 1;
      st[1] = 1;         ed[1] = dim_sz[1];
      st[2] = 1;         ed[2] = dim_sz[2];
      break;
      
    case X_PLUS:
      st[0] = dim_sz[0]; ed[0] = dim_sz[0];
      st[1] = 1;         ed[1] = dim_sz[1];
      st[2] = 1;         ed[2] = dim_sz[2];
      break;
      
    case Y_MINUS:
      st[0] = 1;         ed[0] = dim_sz[0];
      st[1] = 1;         ed[1] = 1;
      st[2] = 1;         ed[2] = dim_sz[2];
      break;
      
    case Y_PLUS:
      st[0] = 1;         ed[0] = dim_sz[0];
      st[1] = dim_sz[1]; ed[1] = dim_sz[1];
      st[2] = 1;         ed[2] = dim_sz[2];
      break;
      
    case Z_MINUS:
      st[0] = 1;         ed[0] = dim_sz[0];
      st[1] = 1;         ed[1] = dim_sz[1];
      st[2] = 1;         ed[2] = 1;
      break;
      
    case Z_PLUS:
      st[0] = 1;         ed[0] = dim_sz[0];
      st[1] = 1;         ed[1] = dim_sz[1];
      st[2] = dim_sz[2]; ed[2] = dim_sz[2];
      break;
  }
}

/**
 @fn void SetBC3D::flipDir_OBC(unsigned* bv, Control* C)
 @brief 外部境界に流入出境界条件が設定されている場合にフラグを立てておく
 */
void SetBC::set_InOut_flag(void)
{
  for (int face=0; face<NOFACE; face++) {
    
    // 計算領域の最外郭領域でないときに，境界処理をスキップ，次のface面を評価
    if( pn.nID[face] >= 0 ) continue;
    
    if ( obc[face].get_BCtype() == OBC_IN_OUT ) inout_flag = true;
  }
}
