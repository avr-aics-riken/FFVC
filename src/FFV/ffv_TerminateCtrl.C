//##################################################################################
//
// FFV-C : Frontflow / violet Cartesian
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2015 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2015 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

/** 
 * @file   ffvTerminateCtrl.C
 * @brief  FFV Terminate Control Class
 * @author aics
 */

#ifndef _NO_CPMLIB_
#include "cpm_ParaManager.h"
#endif

#include "ffv_TerminateCtrl.h"


// スタティック変数の実体
bool FFV_TerminateCtrl::m_ffvTerminateFlag;

// #################################################################
// シグナルハンドラの初期化
void FFV_TerminateCtrl::initialize()
{
  // 中断フラグにfalseをセット
  FFV_TerminateCtrl::m_ffvTerminateFlag = false;

  // シグナルハンドラをセット
  //   SIGINT  : キーボードからの割り込み (Interrupt)
  //   SIGABRT : abort(3) からの中断 (Abort) シグナル
  signal(SIGINT, FFV_TerminateCtrl::terminateSignalHandler);
  signal(SIGABRT, FFV_TerminateCtrl::terminateSignalHandler);
}

// #################################################################
// シグナル発生時のコールバック関数
// 中断フラグにtrueをセット
void FFV_TerminateCtrl::terminateSignalHandler(int signum)
{
#ifndef _NO_CPMLIB_
  // 並列実行時はAbortをコール
  cpm_ParaManager* paraMngr = cpm_ParaManager::get_instance();
  if( paraMngr )
  {
    if( paraMngr->IsParallel() && paraMngr->GetNumRank() > 1 ) paraMngr->Abort(signum);
  }
#endif

  // 中断フラグにtrueをセット
  FFV_TerminateCtrl::m_ffvTerminateFlag = true;
}
 
// #################################################################
// 中断フラグの取得
bool FFV_TerminateCtrl::getTerminateFlag()
{
  // 中断フラグ
  return FFV_TerminateCtrl::m_ffvTerminateFlag;
}

