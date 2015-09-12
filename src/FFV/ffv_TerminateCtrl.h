#ifndef _FFV_TERMINATE_CTRL_
#define _FFV_TERMINATE_CTRL_

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
 * @file   ffvTerminateCtrl.h
 * @brief  FFV Terminate Control Class Header
 * @author aics
 */

#include <signal.h>

class FFV_TerminateCtrl
{
protected:
  /** コンストラクタ */
  FFV_TerminateCtrl(){
  };

  /** デストラクタ */
  ~FFV_TerminateCtrl(){};

  /** 中断フラグ(true:中断、false:継続) */
  static bool m_ffvTerminateFlag;

  /**
   * @brief シグナル発生時のコールバック関数
   * @param [in] signum 発生したシグナル
   */
  static void terminateSignalHandler(int signum);

public:

  /** シグナルハンドラの初期化 */
  static void initialize();

  /** 中断フラグの取得 */
  static bool getTerminateFlag();

};

#endif /* _FFV_TERMINATE_CTRL_ */

