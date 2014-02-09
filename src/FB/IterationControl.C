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
  
  valid        = src->valid;
  alias        = src->alias;
  LinearSolver = src->LinearSolver;
  MaxIteration = src->MaxIteration;
  NormType     = src->NormType;
  eps          = src->eps;
  Sync         = src->Sync;
  omg          = src->omg;
  Naive        = src->Naive;
  Bit3option   = src->Bit3option;
}


// #################################################################
// 固有パラメータを取得
bool IterationCtl::getInherentPara(TextParser* tpCntl, const string base, int& m_naive)
{
  switch (LinearSolver)
  {
    case JACOBI:
      getParaJacobi(tpCntl, base);
      break;
      
    case SOR:
      getParaSOR(tpCntl, base);
      break;
      
    case SOR2SMA:
      getParaSOR2(tpCntl, base);
      m_naive = Naive;
      break;
      
    case RBGS:
      getParaRBGS(tpCntl, base);
      break;
      
    case GMRES:
      getParaGmres(tpCntl, base);
      break;
      
    case PCG:
      getParaPCG(tpCntl, base);
      break;
      
    case PBiCGSTAB:
      getParaPBiCGSTAB(tpCntl, base);
      break;
      
    case VP_ITERATION:
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
  setOmega(tmp);
  
}


// #################################################################
/**
 * @brief PCG反復固有のパラメータを指定する
 * @param [in] tpCntl TextParser pointer
 * @param [in] base   ラベル
 */
void IterationCtl::getParaPCG(TextParser* tpCntl, const string base)
{
  ;
}


// #################################################################
/**
 * @brief PBiCGSTAB反復固有のパラメータを指定する
 * @param [in] tpCntl TextParser pointer
 * @param [in] base   ラベル
 */
void IterationCtl::getParaPBiCGSTAB(TextParser* tpCntl, const string base)
{
  ;
}


// #################################################################
/**
 * @brief RBGS反復固有のパラメータを指定する
 * @param [in] tpCntl TextParser pointer
 * @param [in] base   ラベル
 */
void IterationCtl::getParaRBGS(TextParser* tpCntl, const string base)
{
  ;
}


// #################################################################
/**
 * @brief SOR反復固有のパラメータを指定する
 * @param [in] tpCntl TextParser pointer
 * @param [in] base   ラベル
 */
void IterationCtl::getParaSOR(TextParser* tpCntl, const string base)
{
  getParaJacobi(tpCntl, base);
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
  
  // not mandatory
  label = base + "/NaiveImplementation";
  
  if ( tpCntl->chkLabel(label) )
  {
    if ( !(tpCntl->getInspectedValue(label, str )) )
    {
      Exit(0);
    }
    else
    {
      if ( !strcasecmp(str.c_str(), "on") ) Naive = ON;
    }
  }
  
  // Bit3Test (NOT mandatory)
  label = "/Bit3option";
  
  if ( tpCntl->chkLabel(label) )
  {
    if ( !(tpCntl->getInspectedValue(label, str )) )
    {
      Exit(0);
    }
    else
    {
      if ( !strcasecmp(str.c_str(), "on") ) Bit3option = ON;
    }
  }
  
}



// #################################################################
// ノルムのラベルを返す
string IterationCtl::getNormString()
{
  string nrm;
	
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


// #################################################################
// 線形ソルバの種類を設定する
bool IterationCtl::setLS(const string str)
{
	if     ( !strcasecmp(str.c_str(), "SOR") )          LinearSolver = SOR;
  else if( !strcasecmp(str.c_str(), "SOR2SMA") )      LinearSolver = SOR2SMA;
  else if( !strcasecmp(str.c_str(), "JACOBI") )       LinearSolver = JACOBI;
  else if( !strcasecmp(str.c_str(), "GMRES") )        LinearSolver = GMRES;
  else if( !strcasecmp(str.c_str(), "RBGS") )         LinearSolver = RBGS;
  else if( !strcasecmp(str.c_str(), "PCG") )          LinearSolver = PCG;
  else if( !strcasecmp(str.c_str(), "PBiCGSTAB") )    LinearSolver = PBiCGSTAB;
  else if( !strcasecmp(str.c_str(), "VPiteration") )  LinearSolver = VP_ITERATION;
  else
  {
    return false;
  }
	
  return true;
}
