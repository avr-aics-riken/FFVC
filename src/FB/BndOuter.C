//##################################################################################
//
// Flow Base class
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2015 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2016 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
// Copyright (c) 2016-2020 Research Institute for Information Research, Kyushu University.
// All rights reserved.
//
//##################################################################################

//@file   BndOuter.C
//@brief  FlowBase BoundaryOuter class
//@author aics

#include "BndOuter.h"


// #################################################################
// メンバー変数のコピー
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
  alias     = src->alias;
  p         = src->p;
  ptr_cmp   = src->ptr_cmp;

  for (int i=0; i<3; i++)
  {
    nv[i] = src->nv[i];
  }
  
  for (int i=0; i<5; i++)
  {
    ca[i] = src->ca[i];
    cb[i] = src->cb[i];
  }
}


// #################################################################
// 計算領域の流入出量を保持する
void BoundaryOuter::setDomainMF(const REAL_TYPE vv)
{
  dm[0] = vv;
}


// #################################################################
// 計算領域の流入出量と有効セル数を保持する
void BoundaryOuter::setDomainV(const REAL_TYPE* vv)
{
  dm[0] = vv[0];
  dm[1] = vv[1];
}
