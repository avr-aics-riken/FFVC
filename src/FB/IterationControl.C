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
 * @file   IterationControl.C
 * @brief  FlowBase IterationCtl class
 * @author kero
 */

#include "IterationControl.h"


// #################################################################
// 基本変数のコピー
void IterationCtl::copy(IterationCtl* src)
{
  if ( !src ) return;
  
  valid        = src->valid;
  alias        = src->alias;
  LinearSolver = src->LinearSolver;
  MaxIteration = src->MaxIteration;
  NormType     = src->NormType;
  eps          = src->eps;
  Sync         = src->Sync;
  omg          = src->omg;
}


// #################################################################
// ノルムのラベルを返す
std::string IterationCtl::getNormString()
{
  std::string nrm;
	
	if (NormType == v_div_dbg)
  {
    nrm = "Max. Norm : Divergence of velocity with Monitoring  ### Forced to be selected since Iteration Log is specified ###";
  }
  else if (NormType == v_div_max)
  {
    nrm = "Max. Norm : Divergence of velocity";
  }
  else if (NormType == dx_b)
  {
    nrm = "dx_b : Increment of vector x divided by RHS vector b";
  }
  else if (NormType == r_b)
  {
    nrm = "r_b  : Residual vector divided by RHS vector b";
  }
	else if (NormType == r_r0)
  {
    nrm = "r_r0 : Residual vector divided by initial residual vector";
  }
	
  return nrm;
}
