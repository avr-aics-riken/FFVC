// #################################################################
//
// CAERU Library
//
// Copyright (c) All right reserved. 2012
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
//
// #################################################################

/** 
 * @file Interval_Mngr.C
 * @brief FlowBase Interval_Manager class
 * @author kero
 */

#include "Interval_Mngr.h"


// 指定の時刻になったかどうかを判断する
bool Interval_Manager::isTriggered(const unsigned stp, const double tm, bool d_flag) 
{
  if (d_flag) printf("tm=%f ntm=%f nst=%d itv=%d\n",tm, next_tm, next_step, intvl_step);
  
  if ( step_flag ) {
    return true;
  }
  else {
    if (mode == By_step) {
      if ( (stp >= next_step) && ((stp-1) < next_step) ) {
        next_step += intvl_step;
        step_flag = true;
        return true;
      }
    }
    else {
      if ( (tm >= next_tm) && ((tm-delta_t) < next_tm) ) {
        //next_tm += intvl_tm;
        m_count++;
        next_tm = init_tm + (double)(m_count) * intvl_tm;
        step_flag = true;
        return true;
      }
    }
  }
  
  return false;
}
  

// トリガーを初期化する
bool Interval_Manager::initTrigger(const unsigned stp, const double tm, const double m_dt, const unsigned m_id, const double tscale)
{
  delta_t = m_dt;
  id = m_id;
  
  // 平均値トリガーの場合には，初期化は一度だけを保証し，初期化済みの場合はなにもしない
  if ( id == tg_average ) {
    if ( init_state == true ) return true;
  }
  
  if ( mode == By_step ) {
    if (intvl_step == 0) return false;

    if ( stp < intvl_step ) {
      next_step = intvl_step;
    }
    else {
      next_step = (stp/intvl_step + 1)*intvl_step;
    }
  }
  else {
    if ( intvl_tm == 0.0 ) return false;

    if ( tm < intvl_tm ) {
      next_tm = intvl_tm; // 最小インターバル未満は次インターバル
    }
    else {
      next_tm = floor(tm/intvl_tm + 1.0)*intvl_tm; 
    }
    
    // フレームワークに対するダミー用の値をいれておく >> Control::tell_Interval_2_Sphere()
    if ( intvl_step == 0 ) {
      if ( id == tg_compute ) {
        if ( tscale == 0.0 ) Exit(0);
        intvl_step = intvl_tm/(delta_t*tscale);
      }
      else {
        intvl_step = 1; // 必ず1を入れておくこと．フレームワーク側のタイミング管理を常にtrueにするため．
      }
    }
  }
  
  init_state = true;
  return true;
}
 

// インターバル値をセットする
void Interval_Manager::setInterval(const double m_interval) 
{
  if (mode == By_step) {
    intvl_step = (unsigned)m_interval;
  }
  else {
    intvl_tm = (double)m_interval;
  }
}


// インターバル値を無次元化する
void Interval_Manager::normalizeInterval(const double scale) 
{
  intvl_tm /= scale;
}


// セッションの開始時刻をセットする
void Interval_Manager::setTime_init(const double m_tm) 
{
  init_tm = m_tm;
}
