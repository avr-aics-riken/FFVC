#ifndef _FB_INTVL_MNGR_H_
#define _FB_INTVL_MNGR_H_

/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file Interval_Mngr.h
//@brief FlowBase Interval_Manager class Header
//@author keno, FSI Team, VCAD, RIKEN
//@note 出力ファイルおよび計算時間のタイミングを管理する

#include "math.h"
#include "stdio.h"
#include "FB_Define.h"

class Interval_Manager {
protected:
  unsigned id;         /// 管理対象を表すID
  unsigned mode;       /// 出力指定モード
  unsigned intvl_step; /// ステップ数指定のインターバル
  unsigned next_step;  /// ステップ数指定の場合の次の出力ステップ
  unsigned m_count;    /// セッション内のインターバル数のカウント
  double   intvl_tm;   /// 時刻指定のインターバル（無次元） tg_avstartの場合には，スタート開始時刻として扱う
  double   next_tm;    /// 時刻指定の場合の次の出力時刻（無次元）
  double   delta_t;    /// 時間積分幅（無次元）
  double   init_tm;    /// セッションの初期時刻（無次元）
  bool     step_flag;  /// 1ステップの間に複数回のコールを許すためのフラグ 
  bool     init_state; /// 初期化確認フラグ

public:
  // 管理対象のリスト
  enum flush_trigger {
    tg_compute,  // セッションの計算時間
    tg_console,  // コンソール出力
    tg_history,  // ファイル出力
    tg_instant,  // 瞬時値出力
    tg_average,  // 平均値出力
    tg_sampled,  // サンプリング出力 
    tg_accelra,  // 加速時間
    tg_avstart,  // 平均開始時間
    tg_END
  };
  
  /// ファイル出力タイミングのモード指定
  enum type_IO_spec {
    By_step=1, // ステップ数指定
    By_time    // 時刻指定
  };
  
  Interval_Manager() {
    id         = 0;
    mode       = 0;
    intvl_step = 0;
    next_step  = 0;
    m_count    = 1;
    next_tm    = 0.0;
    intvl_tm   = 0.0;
    delta_t    = 0.0;
    init_tm    = 0.0;
    step_flag  = false;
    init_state = false;
  }
  ~Interval_Manager() {}
  
public:
  bool initTrigger(const unsigned stp, const double tm, const double m_dt, const unsigned m_id, const double tscale=0.0);
  bool isTriggered(const unsigned stp, const double tm, bool d_flag=false);
 
  void setInterval      (const double m_interval);
  void normalizeInterval(const double scale);
  void setTime_init     (const double m_tm);
  
  //@fn void setMode_Step(void)
  //@brief 指定モードをステップにする
  void setMode_Step(void) {
    mode = By_step;
  }
  
  //@fn void setMode_Time(void)
  //@brief 指定モードをステップにする
  void setMode_Time(void) {
    mode = By_time;
  }
  
  //@fn void resetTrigger(void)
  //@brief 各タイムステップの最初にstep_flagをリセットする
  void resetTrigger(void) {
    step_flag = false;
  }
  
  //@fn unsigned getMIntervalStep(void) const
  //@brief インターバル（ステップ）を返す
  unsigned getIntervalStep(void) const {
    return intvl_step;
  }
  
  //@fn double getIntervalTime(void) const
  //@brief インターバル（時刻）を返す
  double getIntervalTime(void) const {
    return intvl_tm;
  }
  
  //@fn bool isStep(void) const 
  //@brief インターバル指定がステップの場合trueを返す
  bool isStep(void) const {
    return (mode==By_step) ? true : false;
  }
};
#endif // _FB_INTVL_MNGR_H_
