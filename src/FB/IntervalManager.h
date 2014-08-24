#ifndef _FB_INTERVAL_MNGR_H_
#define _FB_INTERVAL_MNGR_H_

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

// 使用例
// 1. 初期化処理
//   (1) インスタンス
//     IntervalManager intvl;
//   (2) モードをセット
//     intvl.setMode( IntervalManager::By_step or By_time );
//   (3) インターバルステップ(時刻)のセット
//     intvl.setInterval( (double)stepInterval );
//   (4) 開始ステップ(時刻)のセット(未指定のとき0)
//     intvl.setStart( (double)startStep );
//   (5) 最終ステップ(時刻)のセット(指定の場合、計算の最終ステップは必ず出力する)
//     intvl.setLast( (double)lastStep );
//   (6) トリガーの初期化
//     intvl.initTrigger( step, time, dt );
// 2.ループ処理
//   (1) 出力ステップ(時刻)かどうか判定
//     bool out = intvl.isTriggerd( step, time );

#include <stdlib.h>

class IntervalManager
{
  
// enum
public:
  /// モード
  enum type_IO_spec
  {
    noset  =-1,
    By_step= 1,
    By_time= 2,
  };

  unsigned restartStep; ///< By_time指定時のリスタートステップの指定。浮動小数点では誤差がでるための対応策。
  
  
// メンバ変数
protected:
  type_IO_spec m_mode; ///< 時制モード(noset:未定義、By_step:ステップ間隔指定、By_time:時刻間隔指定)

  unsigned    m_base_step;  ///< インターバルの基点となるステップ(By_stepのとき有効)
  unsigned    m_intvl_step; ///< ステップ間隔(By_stepのとき有効)
  unsigned    m_start_step; ///< 開始ステップ(By_stepのとき有効)
  unsigned    m_last_step;  ///< 最終ステップ(By_stepのとき有効)

  double m_base_time;       ///< インターバルの基点となる時刻(By_timeのとき有効)
  double m_intvl_time;      ///< 時刻間隔(By_timeのとき有効)
  double m_start_time;      ///< 開始時刻(By_timeのとき有効)
  double m_last_time;       ///< 最終時刻(By_timeのとき有効)
  double m_dt;              ///< 計算の時間間隔Δt

  
// メンバ関数
public:
  /// コンストラクタ
  IntervalManager()
  {
    m_mode = noset;
    m_base_step  = 0;
    m_intvl_step = 0;
    m_start_step = 0;
    m_last_step  = 0;
    m_base_time  = 0.0;
    m_intvl_time = 0.0;
    m_start_time = 0.0;
    m_last_time  = -1.0;
    m_dt         = 0.0;
    restartStep  = 0;
  }

  /// デストラクタ
  ~IntervalManager()
  {
  }

  
  /*
   * @brief インターバル(ステップ)の取得
   * @retval インターバルステップ
   */
  unsigned getIntervalStep()
  {
    return m_intvl_step;
  }

  
  /*
   * @brief インターバル(時間)の取得
   * @retval インターバル時間
   */
  double getIntervalTime()
  {
    return m_intvl_time;
  }

  
  /*
   * @brief 時制モードを取得
   * @retval 時制モード(noset:未定義、By_step:ステップ間隔指定、By_time:時刻間隔指定)
   */
  type_IO_spec getMode()
  {
    return m_mode;
  }
  
  
  /*
   * @brief 終了ステップを取得
   * @retval 終了ステップ
   */
  unsigned getLastStep()
  {
    return m_last_step;
  }
  
  
  /*
   * @brief 終了時刻を取得
   * @retval 終了時刻
   */
  double getLastTime()
  {
    return m_last_time;
  }
  
  
  /*
   * @brief 開始ステップを取得
   * @retval 開始ステップ
   */
  unsigned getStartStep()
  {
    return m_start_step;
  }
  
  
  /*
   * @brief 開始時刻を取得
   * @retval 開始時刻
   */
  double getStartTime()
  {
    return m_start_time;
  }
  
  
  /*
   * @brief トリガーの初期化
   * @retval true-成功 / false-時制モードが指定されていない
   * @param [in] step   基準ステップ
   * @param [in] time   基準時刻
   * @param [in] dt     時間積分幅Δt
   * @param [in] d_flag デバッグフラグ
   */
  bool initTrigger( const unsigned step, const double time, const double dt, bool d_flag=false )
  {
    m_dt = dt;
    
    if( m_mode == By_step )
    {
      m_base_step = step;
      
      if( d_flag )
      {
        unsigned next_step = calcNextStep(step);
        if( next_step == m_base_step ) next_step += m_intvl_step;
        
        printf("IntervalManager : mode=By_step, base_step=%d, interval=%d, start=%d, first_step=%d"
               , m_base_step, m_intvl_step, m_start_step, next_step);
        
        if( m_last_step > 0 )
        {
          printf(", last_step=%d\n", m_last_step);
        }
        else
        {
          printf("\n");
        }
      }
    }
    else if( m_mode == By_time )
    {
      m_base_time = time;
      
      if( d_flag )
      {
        double next_time = calcNextTime(time);
        
        if( next_time == m_base_time ) next_time += m_intvl_time;
        printf("IntervalManager : mode=By_time, base_time=%e, delta_t=%e, interval=%e, start=%e, first_time=%e"
               , m_base_time, m_dt, m_intvl_time, m_start_time, next_time);
        if( m_last_time > 0.0 )
        {
          printf(", last_time=%e\n", m_last_time);
        }
        else
        {
          printf("\n");
        }
      }
    }
    else
    {
      return false;
    }
    return true;
  }
  
  
  /*
   * @brief 最終ステップ/時刻かどうか
   * @retval true-最終ステップ/時刻 / false-最終ステップ/時刻でない
   * @param [in] step 現在のステップ
   * @param [in] time 現在時刻
   */
  bool isLast( const unsigned step, const double time )
  {
    if( m_mode == By_step )
    {
      if( m_last_step <= 0 ) return false;
      if( m_last_step == step ) return true;
    }
    else if( m_mode == By_time )
    {
      if( m_last_time <= 0 ) return false;
      if( time <= m_last_time && m_last_time < time + m_dt ) return true;
    }
    return false;
  }
  
  
  /*
   * @brief 開始しているかをチェック
   * @retval true-開始している / false-開始していない
   * @param [in] step 現在ステップ
   * @param [in] time 現在時刻
   */
  bool isStarted( const unsigned step, const double time )
  {
    if( m_mode == By_step )
    {
      if( m_start_step <= step ) return true;
    }
    else if( m_mode == By_time )
    {
      if( m_start_time <= time ) return true;
    }
    return false;
  }
  
  
  /*
   * @brief タイミングの評価
   * @retval true-出力タイミングである / false-出力タイミングではない
   * @param [in] step       現在のステップ
   * @param [in] time       現在の時刻
   * @param [in] forced_out 強制出力モード
   * @param [in] d_flag     デバッグフラグ
   * @note 出力ステップ(時刻)であっても、同じステップ(時刻)に既にisTriggeredがコールされている場合はfalseを返す(重複出力対応)
   */
  bool isTriggered( const unsigned step, const double time, bool forced_out=false, bool d_flag=false )
  {
    // 強制出力フラグ
    if( forced_out )
    {
      return true;
    }
    
    // nosetのときは必ずtrue(上位プログラムで判定しているものとする)
    if( m_mode == noset )
    {
      return true;
    }
    
    // 開始しているか
    if( !isStarted(step, time) )
    {
      return false;
    }
    
    // ステップ指定のとき
    if( m_mode == By_step )
    {
      // intervalが正常か
      if( m_intvl_step <= 0 )
      {
        return false;
      }
      // 次の出力ステップかどうか
      if( step == calcNextStep( step ) )
      {
        if( d_flag )
        {
          printf("IntervalManager::isTriggerd : step=%d\n", step);
        }
        return true;
      }
      
      // 最終ステップかどうか
      if( isLast( step, time ) )
      {
        if( d_flag )
        {
          printf("IntervalManager::isTriggerd : last step=%d\n", step);
        }
        return true;
      }
    }
    // 時刻指定のとき
    else if( m_mode == By_time )
    {
      // intervalが正常か
      if( m_intvl_time <= 0.0 )
      {
        return false;
      }
      
      // 次の出力時刻
      double next_time = calcNextTime( time );
      
      // 次の出力時刻かどうか
      if( time <= next_time && next_time < time + m_dt )
      {
        if( d_flag )
        {
          printf("IntervalManager::isTriggerd : time=%e\n", time);
        }
        return true;
      }
      
      // 最終時刻かどうか
      if( isLast( step, time ) )
      {
        if( d_flag )
        {
          printf("IntervalManager::isTriggerd : last time=%e\n", time);
        }
        return true;
      }
    }
    return false;
  }
  
  
  /*
   * @brief 時刻モードがステップのとき，インターバルの計算に使われる全ての時間を scale で無次元化する
   * @param [in] scale 無次元化の代表時間
   */
  void normalizeTime( const double scale )
  {
    if( m_mode != By_time )
    {
      return;
    }
    normalizeBaseTime(scale);
    normalizeIntervalTime(scale);
    normalizeStartTime(scale);
    normalizeLastTime(scale);
    normalizeDelteT(scale);
    return;
  }
  
  
  /*
   * @brief インターバル値のセット
   * @retval true-成功 / false-時制モードが指定されていない
   * @param [in] interval 制御間隔
   */
  bool setInterval( const double interval )
  {
    if( m_mode == By_step )
    {
      m_intvl_step = (unsigned)interval;
    }
    else if( m_mode == By_time )
    {
      m_intvl_time = interval;
    }
    else
    {
      return false;
    }
    return true;
  }
  
  
  /*
   * @brief 最終ステップ/時刻のセット
   * @retval true-成功 / false-時制モードが指定されていない
   * @param [in] last 終了時刻
   */
  bool setLast( const double last )
  {
    if( m_mode == By_step )
    {
      m_last_step = (unsigned)last;
    }
    else if( m_mode == By_time )
    {
      m_last_time = last;
    }
    else
    {
      return false;
    }
    return true;
  }

  
  /*
   * @brief 時制モードをセット
   * @param [in] mode 時制モード(noset:未定義、By_step:ステップ間隔指定、By_time:時刻間隔指定)
   */
  void setMode( const type_IO_spec mode )
  {
    m_mode = mode;
  }
  
  
  /*
   * @brief 開始ステップ/時刻のセット
   * @retval true-成功 / false-時制モードが指定されていない
   * @param [in] 開始時刻
   */
  bool setStart( const double start )
  {
    if( m_mode == By_step )
    {
      m_start_step = (unsigned)start;
    }
    else if( m_mode == By_time )
    {
      m_start_time = start;
    }
    else
    {
      return false;
    }
    return true;
  }
  
  
  
protected:
  
  /*
   * @brief 次の出力ステップを計算
   * @retval 次の出力ステップ数
   * @param [in] step 現在のステップ
   */
  unsigned calcNextStep( const unsigned step )
  {
    unsigned s_step = (m_start_step > step) ? m_start_step : step;
    unsigned inc    = ((s_step-m_base_step)%m_intvl_step==0) ? 0 : 1;
    unsigned next_step = unsigned((s_step-m_base_step)/m_intvl_step+inc) * m_intvl_step + m_base_step;
    return next_step;
  }

  /*
   * @brief 次の出力時刻を計算
   * @retval 次の出力時刻
   * @param [in] time 現在時刻
   */
  double calcNextTime( const double time )
  {
    //unsigned s_time = (m_start_time > time) ? m_start_time : time;
    double s_time = (m_start_time > time) ? m_start_time : time;
    unsigned inc    = (dmod(s_time-m_base_time, m_intvl_time)==0.0) ? 0 : 1;
    double next_time = unsigned((s_time-m_base_time)/m_intvl_time+inc) * m_intvl_time + m_base_time;
    return next_time;
  }

  /*
   * @brief 実数の余り(fortranのmodと同じ)
   * @retval 剰余 
   * @param [in] a 被除数
   * @param [in] b 除数
   */
  static double dmod( double a, double b )
  {
    return a - int(a/b) * b;
  }

 
  /*
   * @brief 基準時刻を scale で無次元化する
   * @param [in] scale 無次元化の代表時間
   */
  void normalizeBaseTime(const double scale)
  {
    m_base_time /= scale;
  }
  
  /*
   * @brief インターバルを scale で無次元化する
   * @param [in] scale 無次元化の代表時間
   */
  void normalizeIntervalTime(const double scale)
  {
    m_intvl_time /= scale;
  }
  
  /*
   * @brief 開始時刻を scale で無次元化する
   * @param [in] scale 無次元化の代表時間
   */
  void normalizeStartTime(const double scale)
  {
    m_start_time /= scale;
  }
  
  /*
   * @brief 終了時刻を scale で無次元化する
   * @param [in] scale 無次元化の代表時間
   */
  void normalizeLastTime(const double scale)
  {
    m_last_time /= scale;
  }
  
  /*
   * @brief 計算の時間間隔Δtを scale で無次元化する
   * @param [in] scale 無次元化の代表時間
   */
  void normalizeDelteT(const double scale)
  {
    m_dt /= scale;
  }
  
};

#endif  /* _FB_INTERVAL_MNGR_H_ */

