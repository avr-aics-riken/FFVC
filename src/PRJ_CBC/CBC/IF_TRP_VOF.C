/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2010
 *
 */

//@file IF_TRP_VOF.C
//@brief SklSolverCBC class
//@author keno, FSI Team, VCAD, RIKEN

#include "SklSolverCBC.h"

/**
 @fn void SklSolverCBC::IF_TRP_VOF(void)
 @brief Interfaceの移流方程式を解く
 */
void SklSolverCBC::IF_TRP_VOF(void)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  
  // local variables
  SKL_REAL tm = SklGetTotalTime();  /// 計算開始からの積算時刻
  SKL_REAL dt = SklGetDeltaT();     /// 時間積分幅
  SKL_REAL flop_count=0.0;          /// 浮動小数演算数
  SKL_REAL nrm=0.0;
  
  // Dimensions
  SKL_REAL *v;    
  SKL_REAL *vof;  
  SKL_REAL *ws;   
  unsigned *bcv;  

  v = vof = NULL;
  ws = NULL;
  bcv = NULL;
  dt = NULL;

  // point Data
  if( !(ws  = dc_ws->GetData()) )   assert(0);
  if( !(bcv = dc_bcv->GetData()) )  assert(0);
  if( !(vof = dc_vof->GetData()) )  assert(0);
  if( !(v   = dc_vf0->GetData()) )  assert(0);
  
  // convection
  TIMING_start(tm_vof_cnv);
  vof_uwd_(vof, sz, gc, v00, &dt, dh, v, ws, (int*)bcv, &flop_count);
  TIMING_stop(tm_vof_cnv, flop_count);
  
  TIMING_start(tm_vof_cnv_comm);
  dc_vof->CommBndCell(guide);
  TIMING_stop(tm_vof_cnv_comm, 0.0);
  
  // variable range cutoff
  //TIMING_start(tm_vof_range);
  //CU.CutOffRange (t, cmp, bcp, &C);
  //TIMING_stop(tm_vof_range, flop_task);
}
