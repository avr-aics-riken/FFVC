#ifndef _FB_INTVL_MNGR_H_
#define _FB_INTVL_MNGR_H_

// #################################################################
//
// CAERU Library
//
// Copyright (c) 2012-2013  All right reserved.
//
// Institute of Industrial Science, University of Tokyo, Japan. 
//
// #################################################################

/** 
 * @file   Interval_Mngr.h
 * @brief  FlowBase Interval_Manager class Header 出力ファイルおよび計算時間のタイミングを無次元で管理する
 * @author kero
 */

#include "math.h"
#include "stdio.h"
#include "FB_Define.h"

class Interval_Manager {
protected:
  int id;              ///< 管理対象を表すID
  int mode;            ///< 出力指定モード
  int start_step;      ///< 開始ステップ数
  int intvl_step;      ///< ステップ数指定のインターバル
  int next_step;       ///< ステップ数指定の場合の次の出力ステップ
  int m_count;         ///< セッション内のインターバル数のカウント
  double   intvl_tm;   ///< 時刻指定のインターバル（無次元）
  double   start_tm;   ///< 開始時刻（無次元）
  double   next_tm;    ///< 時刻指定の場合の次の出力時刻（無次元）
  double   delta_t;    ///< 時間積分幅（無次元）
  bool     step_flag;  ///< 1ステップの間に複数回のコールを許すためのフラグ 

public:
  /** 管理対象のリスト */
  enum flush_trigger {
    tg_compute,  ///< セッションの計算時間
    tg_console,  ///< コンソール出力
    tg_history,  ///< ファイル出力
    tg_instant,  ///< 瞬時値出力
    tg_average,  ///< 平均値出力
    tg_sampled,  ///< サンプリング出力 
    tg_accelra,  ///< 加速時間
    tg_plot3d,   ///< PLOT3D瞬時値出力
    tg_END
  };
  
  /** ファイル出力タイミングのモード指定 */
  enum type_IO_spec {
    By_step=1, ///< ステップ数指定
    By_time    ///< 時刻指定
  };
  
  Interval_Manager() {
    mode       = 0;
    intvl_step = 0;
    next_step  = 0;
    m_count    = 1;
    start_step = 0;
    next_tm    = 0.0;
    intvl_tm   = 0.0;
    delta_t    = 0.0;
    start_tm   = 0.0;
    step_flag  = false;
  }
  ~Interval_Manager() {}
  
public:
  
  // トリガーを初期化する
  bool initTrigger(const int stp, const double tm, const double m_dt, const int m_id);
  
  
  // 指定の時刻になったかどうかを判断する
  bool isTriggered(const int stp, const double tm, bool d_flag=false);
 
  
  // インターバル値をセットする
  void setInterval(const double m_interval);
  
  
  // 開始時刻とインターバル時刻を無次元化する
  void normalizeTime(const double scale);
  
  
  // セッションの開始時刻をセットする
  void setStart(const double m_start);
  
  
  // 指定モードをステップにする
  void setMode_Step() 
  {
    mode = By_step;
  }
  

  // 指定モードをステップにする
  void setMode_Time() 
  {
    mode = By_time;
  }
  
  
  // 各タイムステップの最初にstep_flagをリセットする 
  void resetTrigger() 
  {
    step_flag = false;
  }
  

  // インターバル（ステップ）を返す
  int getIntervalStep() const 
  {
    return intvl_step;
  }
  

  // インターバル（時刻）を返す
  double getIntervalTime() const 
  {
    return intvl_tm;
  }
  
  
  // 開始ステップを返す
  int getStartStep() const
  {
    return start_step;
  }
  
  
  // 開始時刻を返す
  double getStartTime() const
  {
    return start_tm;
  }
  

  // インターバル指定がステップの場合trueを返す
  bool isStep() const 
  {
    return (mode==By_step) ? true : false;
  }
  
  // 開始時刻を過ぎているかどうかを判断
  bool isStarted(const unsigned m_step, const double m_time);

};

#endif // _FB_INTVL_MNGR_H_
