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
 * @file SetBC.C
 * @brief FlowBase SetBC class
 * @author kero
 */

#include <math.h>
#include "SetBC.h"


// 作業用ポインタのコピー
void SetBC::setWorkList(CompoList* m_CMP, MediumList* m_MAT)
{
  if ( !m_CMP ) Exit(0);
  cmp = m_CMP;
  
  if ( !m_MAT ) Exit(0);
  mat = m_MAT;
}


// CPMのインポート
void SetBC::importCPM(cpm_ParaManager* m_paraMngr)
{
  if ( !m_paraMngr ) Exit(0);
  paraMngr = m_paraMngr;
}


// 必要な値のコピー
void SetBC::setControlVars(Control* Cref, MediumList* mat, CompoList* cmp, ReferenceFrame* RF, Intrinsic* ExRef)
{
  guide          = Cref->guide;
  imax = size[0] = Cref->imax;
  jmax = size[1] = Cref->jmax;
  kmax = size[2] = Cref->kmax;
  gc = &guide;
  ix = &(Cref->imax);
  jx = &(Cref->jmax);
  kx = &(Cref->kmax);
  ixc = imax;
  jxc = jmax;
  kxc = kmax;
  dh        = Cref->dh;
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
  
  for (int i=0; i<3; i++) Lbx[i] = Cref->Lbx[i];
  
  if ( Cref->isHeatProblem() ) pei = Cref->getRcpPeclet();

  Dp1       = Cref->Domain_p1;
  Dp2       = Cref->Domain_p2;
  
  dim_sz[0] = imax;
  dim_sz[1] = jmax;
  dim_sz[2] = kmax;
  
  Ex = ExRef;
  
  // get reference values >> 媒質はIDがユニークに定まる
  int m;
  for (int n=Cref->NoBC+1; n<=Cref->NoCompo; n++) {
    if ( cmp[n].getMatOdr() == Cref->RefMat ) {
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


// 静止座標系のときの流出速度制御の値を計算する
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


// 外部境界処理用のループインデクスを取得する
void SetBC::getOuterLoopIdx(const int face, int* st, int* ed)
{
  switch (face) 
  {
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
