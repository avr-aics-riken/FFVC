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
 * @file   Interval_Mngr.C
 * @brief  FlowBase Interval_Manager class
 * @author kero
 */

#include "Interval_Mngr.h"
  

// #################################################################
// トリガーを初期化する
bool Interval_Manager::initTrigger(const int stp, const double tm, const double m_dt, const int m_id)
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
      next_step = start_step + (stp/intvl_step)*intvl_step;
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
      next_tm = start_tm + floor(tm/intvl_tm)*intvl_tm;
    }
  }
  
  // 全計算時間は，step/timeの両方で制御 指定時刻を越えるように+1
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
// セッションの開始時刻が過ぎているかを判断する
bool Interval_Manager::isStarted(const unsigned m_step, const double m_time)
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


// #################################################################
// 指定の時刻になったかどうかを判断する
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
      if ( ((stp+1) > next_step) && (stp <= next_step) )
      {
        m_count++;
        next_step = start_step + m_count * intvl_step;
        step_flag = true;
        return true;
      }
    }
    else
    {
      if ( ((tm+delta_t) > next_tm) && (tm <= next_tm) )
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
// 開始時刻とインターバル時刻を無次元化する
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
 * @brief セッションの開始時刻をセットする
 * @param [in] m_step 無次元評価ステップ
 * @param [in] m_time 無次元評価時刻
 */
void Interval_Manager::setStart(const unsigned m_step, const double m_time)
{
  if (mode == By_step)
  {
    start_step = (int)m_step;
  }
  else
  {
    start_tm = (double)m_time;
  }
}
