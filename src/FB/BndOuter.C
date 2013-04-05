// #################################################################
//
// CAERU Library
//
// Copyright (c) 2012-2013  All right reserved.
//
// Institute of Industrial Science, University of Tokyo, Japan. 
//
// #################################################################

//@file   BndOuter.C
//@brief  FlowBase BoundaryOuter class
//@author kero

#include "BndOuter.h"


// #################################################################
/**
 * @brief 計算領域の流入出量と有効セル数を保持する
 * @param [in] vv   保持する値
 */
void BoundaryOuter::setDomainV(const REAL_TYPE* vv)
{
  dm[0] = vv[0];
  dm[1] = vv[1];
}

// #################################################################
/**
 * @brief 計算領域の流入出量を保持する
 * @param [in] vv   保持する値
 */
void BoundaryOuter::setDomainMF(const REAL_TYPE vv)
{
  dm[0] = vv;
}


// #################################################################
/**
 * @brief メンバー変数のコピー
 * @param [in] src   BoundaryOuterクラス
 */
void BoundaryOuter::dataCopy(BoundaryOuter* src)
{
  BCclass   = src->BCclass;
  wallType  = src->wallType;
  drv_dir   = src->drv_dir;
  drv_lid   = src->drv_lid;
  gc_medium = src->gc_medium;
  v_profile = src->v_profile;
  Face_mode = src->Face_mode;
  hType     = src->hType;
  HTref     = src->HTref;
  HTmode    = src->HTmode;
  Prdc_mode = src->Prdc_mode;
  pType     = src->pType;
  valid_cell= src->valid_cell;
  var1      = src->var1;
  var2      = src->var2;
  label     = src->label;
  alias     = src->alias;
  
  p         = src->p;

  for (int i=0; i<3; i++) {
    nv[i] = src->nv[i];
  }
  
  for (int i=0; i<5; i++) {
    ca[i] = src->ca[i];
    cb[i] = src->cb[i];
  }
}
