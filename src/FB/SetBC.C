//##################################################################################
//
// Flow Base class
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2014 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

/** 
 * @file   SetBC.C
 * @brief  FlowBase SetBC class
 * @author aics
 */

#include "SetBC.h"
#include <math.h>

// #################################################################
// 無次元の媒質情報をコピー
void SetBC::copyNDmatTable(const int m_compo, const double* m_mtbl, const REAL_TYPE* m_vtbl)
{
  // Fortran用のデータ保持配列 >> m_mtbl[C.NoCompo+1][3]のイメージ
  if ( !(mtbl = new double[3*(m_compo+1)]) ) Exit(0);
  
  for (int i=0; i<3*(m_compo+1); i++)
  {
    mtbl[i] = m_mtbl[i];
  }
  
  if ( !(vtbl = new REAL_TYPE[7*(m_compo+1)]) ) Exit(0);
  
  for (int i=0; i<7*(m_compo+1); i++)
  {
    vtbl[i] = m_vtbl[i];
  }
}


// #################################################################
// 作業用ポインタのコピー
void SetBC::importCMP_MAT(CompoList* m_CMP, MediumList* m_MAT)
{
  if ( !m_CMP ) Exit(0);
  cmp = m_CMP;
  
  if ( !m_MAT ) Exit(0);
  mat = m_MAT;
}


// #################################################################
// クラスに必要な変数のコピー
void SetBC::setControlVars(Control* Cref, const MediumList* mat, const ReferenceFrame* RF, Intrinsic* ExRef)
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
  Unit_Prs  = Cref->Unit.Prs;
	BasePrs   = Cref->BasePrs;
  Rayleigh  = Cref->Rayleigh;
  Grashof   = Cref->Grashof;
  Prandtl   = Cref->Prandtl;
  NoCompo   = Cref->NoCompo;
  NoMedium  = Cref->NoMedium;
  //rho       = Cref->RefDensity;
  
  if ( Cref->isHeatProblem() ) pei = Cref->getRcpPeclet();

  Dp1       = Cref->Domain_p1;
  Dp2       = Cref->Domain_p2;
  
  Ex = ExRef;
  
  // get reference values >> 媒質はIDがユニークに定まる

  for (int n=1; n<=NoMedium; n++)
  {
    if ( n == Cref->RefMat )
    {
      if ( mat[n].getState() == FLUID )
      {
        rho_0    = mat[n].P[p_density];
        cp_0     = mat[n].P[p_specific_heat];
        lambda_0 = mat[n].P[p_thermal_conductivity];
      }
      else
      {
        rho_0    = mat[n].P[p_density];
        cp_0     = mat[n].P[p_specific_heat];
        lambda_0 = mat[n].P[p_thermal_conductivity];
      }
    }
  }

}


// #################################################################
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
