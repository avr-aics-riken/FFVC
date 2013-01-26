#ifndef _FFV_TERMINATE_CTRL_
#define _FFV_TERMINATE_CTRL_

// #################################################################
//
// FFV : Frontflow / violet Cartesian
//
// Copyright (c) 2012-2013 All right reserved.
//
// Institute of Industrial Science, University of Tokyo, Japan. 
//
// #################################################################

/** 
 * @file   ffvTerminateCtrl.h
 * @brief  FFV Terminate Control Class Header
 * @author kero
 */

#include <signal.h>

class FFV_TerminateCtrl
{
protected:
  /** コンストラクタ */
  FFV_TerminateCtrl(){};

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

