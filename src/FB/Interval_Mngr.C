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
 * @file   Interval_Mngr.C
 * @brief  FlowBase Interval_Manager class
 * @author kero
 */

#include "Interval_Mngr.h"


// #################################################################
/**
 * @brief 指定の時刻になったかどうかを判断する
 * @retval 出力タイミングの場合，trueを返す
 * @param [in] stp    現在ステップ
 * @param [in] tm     現時刻
 * @param [in] d_flag 表示用フラグ（デバッグ）
 * @note 指定時刻を過ぎて，かつ1時刻前が指定時刻に満たない場合が出力タイミングとなる
 */
bool Interval_Manager::isTriggered(const int stp, const double tm, bool d_flag)
{
  if (d_flag) printf("tm=%f ntm=%f nst=%d itv=%d\n",tm, next_tm, next_step, intvl_step);
  
    
  if ( step_flag )
  {
    return true;
  }
  else
  {
    if (mode == By_step)
    {
      if ( (stp >= next_step) && ((stp-1) < next_step) )
      {
        m_count++;
        next_step = start_step + m_count * intvl_step;
        step_flag = true;
        return true;
      }
    }
    else
    {
      if ( (tm >= next_tm) && ((tm-delta_t) < next_tm) )
      {
        m_count++;
        next_tm = start_tm + (double)(m_count) * intvl_tm;
        step_flag = true;
        return true;
      }
    }
  }
  
  return false;
}
  

// #################################################################
/**
 * @brief トリガーを初期化する
 * @param [in] stp    現在ステップ
 * @param [in] tm     現時刻（無次元）
 * @param [in] m_dt   時間積分幅（無次元）
 * @param [in] m_id   管理対象を示すID
 * @param [in] tscale タイムスケール
 */
bool Interval_Manager::initTrigger(const int stp, const double tm, const double m_dt, const int m_id, const double tscale)
{
  delta_t = m_dt;
  
  
  if ( mode == By_step )
  {
    if (intvl_step == 0) return false;

    if ( stp < (start_step + intvl_step) )
    {
      next_step = start_step + intvl_step;
    }
    else
    {
      next_step = start_step + (stp/intvl_step + 1)*intvl_step;
    }
  }
  else // by_time
  {
    if ( intvl_tm == 0.0 ) return false;

    if ( tm < (start_tm + intvl_tm) )
    {
      next_tm = start_tm + intvl_tm;
    }
    else
    {
      next_tm = start_tm + floor(tm/intvl_tm + 1.0)*intvl_tm; 
    }
  }
  
  // 全計算時間は，step/timeの両方で制御
  if ( m_id == tg_compute )
  {
    if ( mode == By_step )
    {
      intvl_tm = (double)(intvl_step+1) * delta_t;
    }
    else
    {
      intvl_step = (int)(intvl_tm/delta_t + 1);
    }
  }
  
  return true;
}
 

// #################################################################
/**
 * @brief インターバル値をセットする
 * @param [in] m_interval インターバル値
 */
void Interval_Manager::setInterval(const double m_interval)
{
  if (mode == By_step)
  {
    intvl_step = (int)m_interval;
  }
  else
  {
    intvl_tm = (double)m_interval;
  }

}


// #################################################################
/**
 * @brief インターバル値を無次元化する
 * @param [in] scale  時間スケール
 * @note BY_stepの場合には変化なし
 */
void Interval_Manager::normalizeTime(const double scale)
{
  if (mode == By_time)
  {
    intvl_tm /= scale;
    start_tm /= scale;
  }
}


// #################################################################
/**
 * @brief セッションの開始時刻をセットする
 * @param [in] m_start 無次元時刻 or ステップ
 */
void Interval_Manager::setStart(const double m_start)
{
  if (mode == By_step)
  {
    start_step = (int)m_start;
  }
  else
  {
    start_tm = (double)m_start;
  }
}

// #################################################################
// 開始時刻が過ぎている場合、true
bool Interval_Manager::isStarted(const double m_time, const unsigned m_step)
{
  if (mode == By_step)
  {
    if ( start_step <= (int)m_step ) return true;
  }
  else
  {
    if ( start_tm <= (double)m_time ) return true;
  }
  return false;
}
