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
//##################################################################################

/**
 * @file   IterationControl.C
 * @brief  FlowBase IterationCtl class
 * @author aics
 */

#include "IterationControl.h"


// #################################################################
// 基本変数のコピー
void IterationCtl::copy(IterationCtl* src)
{
  if ( !src ) return;
  
  residual     = src->residual;
  error        = src->error;
  eps_res      = src->eps_res;
  eps_err      = src->eps_err;
  omg          = src->omg;
  ErrNorm      = src->ErrNorm;
  ResNorm      = src->ResNorm;
  MaxIteration = src->MaxIteration;
  LinearSolver = src->LinearSolver;
  LoopCount    = src->LoopCount;
  Sync         = src->Sync;
  alias        = src->alias;
  precondition = src->precondition;
  InnerItr     = src->InnerItr;
  smoother     = src->smoother;
}


// #################################################################
// 固有パラメータを取得
bool IterationCtl::getInherentPara(TextParser* tpCntl, const string base)
{
  switch (LinearSolver)
  {
      case JACOBI:
      getParaJacobi(tpCntl, base);
      break;
      
      case SOR:
      getParaJacobi(tpCntl, base);
      break;
      
      case SOR2SMA:
      getParaSOR2(tpCntl, base);
      break;
      
      //case GMRES:
      //getParaGmres(tpCntl, base);
      //break;
      
      //case PCG:
      //getParaPCG(tpCntl, base);
      //break;
      
      case BiCGSTAB:
      getParaBiCGSTAB(tpCntl, base);
      break;
      
      default:
      return false;
  }
  
  return true;
}


// #################################################################
/**
 * @brief Gmres反復固有のパラメータを指定する
 * @param [in] tpCntl TextParser pointer
 * @param [in] base   ラベル
 */
void IterationCtl::getParaGmres(TextParser* tpCntl, const string base)
{
  ;
}


// #################################################################
/**
 * @brief Jacobi反復固有のパラメータを指定する
 * @param [in] tpCntl TextParser pointer
 * @param [in] base   ラベル
 */
void IterationCtl::getParaJacobi(TextParser* tpCntl, const string base)
{
  string str, label;
  double tmp=0.0; // 加速係数
  
  label = base + "/Omega";
  
  if ( !(tpCntl->getInspectedValue(label, tmp )) )
  {
    Exit(0);
  }
  omg = tmp;
  
}


// #################################################################
/**
 * @brief BiCGSTAB反復固有のパラメータを指定する
 * @param [in] tpCntl TextParser pointer
 * @param [in] base   ラベル
 */
void IterationCtl::getParaBiCGSTAB(TextParser* tpCntl, const string base)
{
  string str, label;
  
  label = base + "/Preconditioner";
  
  if ( !tpCntl->chkLabel(label) )
  {
    Exit(0);
  }
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
    Exit(0);
  }
  else
  {
    if      ( !strcasecmp(str.c_str(), "none") )
    {
      precondition = OFF;
    }
    else if ( !strcasecmp(str.c_str(), "sor") )
    {
      precondition = ON;
      smoother = SOR;
    }
    else if ( !strcasecmp(str.c_str(), "sor2sma") )
    {
      precondition = ON;
      smoother = SOR2SMA;
    }
    else
    {
      Exit(0);
    }
  }
  
  if ( precondition == OFF ) return;

  
  int ct = 0;
  label = base + "/InnerIteration";
  if ( !(tpCntl->getInspectedValue(label, ct )) )
  {
    Exit(0);
  }
  InnerItr = ct;
  
  getParaSOR2(tpCntl, base);

}



// #################################################################
/**
 * @brief RB-SOR反復固有のパラメータを指定する
 * @param [in] tpCntl TextParser pointer
 * @param [in] base   ラベル
 */
void IterationCtl::getParaSOR2(TextParser* tpCntl, const string base)
{
  string str, label;
  
  getParaJacobi(tpCntl, base);
  
  label = base + "/commMode";
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
    Exit(0);
  }
  if ( !strcasecmp(str.c_str(), "sync") )
  {
    setSyncMode(comm_sync);
  }
  else if ( !strcasecmp(str.c_str(), "async") )
  {
    setSyncMode(comm_async);
  }
  else
  {
    Exit(0);
  }
  
}


// #################################################################
// 残差ノルムのラベルを返す
string IterationCtl::getResNormString()
{
  string nrm;
	
  if (ResNorm == nrm_r_b)
  {
    nrm = "r_b  : Residual vector divided by constant vector b";
  }
  else if (ResNorm == nrm_r_x)
  {
    nrm = "r_x  : Residual vector divided by solution vector x";
  }
	else if (ResNorm == nrm_r_r0)
  {
    nrm = "r_r0 : Residual vector divided by initial residual vector";
  }
	
  return nrm;
}


// #################################################################
// 誤差ノルムのラベルを返す
string IterationCtl::getErrNormString()
{
  string nrm;
  
  if (ErrNorm == nrm_dx)
  {
    nrm = "dx   : Increment vector x(k+1)-x(k)";
  }
  else if (ErrNorm == nrm_dx_x)
  {
    nrm = "dx_x : Increment vector x(k+1)-x(k) divided by solution vector x(k)";
  }
  else if (ErrNorm == nrm_div_max)
  {
    nrm = "div_max : Max(Divergence)";
  }
  else if (ErrNorm == nrm_div_l2)
  {
    nrm = "div_L2 : |Divergence|_L2";
  }
  else
  {
    Exit(0);
  }
  
  return nrm;
}


// #################################################################
// 誤差ノルムのタイプを保持
bool IterationCtl::setErrType(TextParser* tpCntl, const string label)
{
  string str;
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
    return false;
  }
  
  if ( !strcasecmp(str.c_str(), "DeltaXbyX") )
  {
    ErrNorm = nrm_dx_x;
  }
  else if ( !strcasecmp(str.c_str(), "DeltaX") )
  {
    ErrNorm = nrm_dx;
  }
  else
  {
    return false;
  }
  
  return true;
}


// #################################################################
// 線形ソルバの種類を設定する
bool IterationCtl::setLS(const string str)
{
	if     ( !strcasecmp(str.c_str(), "SOR") )          LinearSolver = SOR;
  else if( !strcasecmp(str.c_str(), "SOR2SMA") )      LinearSolver = SOR2SMA;
  else if( !strcasecmp(str.c_str(), "JACOBI") )       LinearSolver = JACOBI;
  else if( !strcasecmp(str.c_str(), "GMRES") )        LinearSolver = GMRES;
  else if( !strcasecmp(str.c_str(), "PCG") )          LinearSolver = PCG;
  else if( !strcasecmp(str.c_str(), "BiCGstab") )     LinearSolver = BiCGSTAB;
  else
  {
    return false;
  }
	
  return true;
}
