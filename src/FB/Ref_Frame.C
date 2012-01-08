/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2010
 *
 */

//@file Ref_Frame.C
//@brief FlowBase ReferenceFrame class
//@author keno, FSI Team, VCAD, RIKEN

#include "Ref_Frame.h"

/**
 @fn void ReferenceFrame::setFrame(const unsigned m_frame)
 @brief 変数をセットする
 @param m_frame 参照フレームの種類
 */
void ReferenceFrame::setFrame(const unsigned m_frame)
{
  Frame = m_frame;
}

/**
 @fn void ReferenceFrame::setAccel(const double m_timeAccel)
 @brief 変数をセットする
 @param m_timeAccel 加速時間（無次元）
 */
void ReferenceFrame::setAccel(const double m_timeAccel)
{
  TimeAccel = m_timeAccel;
}

/**
 @fn void ReferenceFrame::setGridVel(const double* m_Gvel)
 @brief 変数をセットする
 @param m_Gvel 格子速度成分の単位方向ベクトル
 */
void ReferenceFrame::setGridVel(const double* m_Gvel)
{
  GridVel[0] = m_Gvel[0];
  GridVel[1] = m_Gvel[1];
  GridVel[2] = m_Gvel[2];
}

/**
 @fn void setV00(const double time, const bool init)
 @brief 参照速度を計算する
 @param time 時刻（無次元）
 @param init フラグ
 @todo 回転
 */
void ReferenceFrame::setV00(const double time, const bool init) 
{
  if (init == true) {
    printf("\tReference velocity is set to 1.0 for checking.\n\n");
    v00[0]=1.0;
  }
  else {
    if ( TimeAccel == 0.0 ) {
      v00[0] = 1.0;
    }
    else {
      const double c_pai = (double)(2.0*asin(1.0));
      v00[0] = 0.5*(1.0-cos(c_pai*time/(TimeAccel)));
      if ( time > (TimeAccel) ) v00[0] = 1.0;
    }
  }
  
  double u0 = v00[0];
  
  switch (Frame) {
    case frm_static:
      v00[1] = 0.0;
      v00[2] = 0.0;
      v00[3] = 0.0;
      break;
      
    case frm_translation:
      v00[1] = u0*GridVel[0];  // v0x
      v00[2] = u0*GridVel[1];  // v0y
      v00[3] = u0*GridVel[2];  // v0z
      break;
      
    case frm_rotation:
      break;
  }
  
}
