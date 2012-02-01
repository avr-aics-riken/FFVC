#ifndef _SKL_FB_REF_FRAME_H_
#define _SKL_FB_REF_FRAME_H_

/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2010
 *
 */

//@file Ref_Frame.h
//@brief FlowBase ReferenceFrame class Header
//@author keno, FSI Team, VCAD, RIKEN

#include <math.h>
#include <stdio.h>

class ReferenceFrame {
  
protected:
  unsigned Frame;    /// 参照座標系
  double TimeAccel;  /// 加速時間（無次元）
  double v00[4];     /// 参照速度（無次元
  double GridVel[3]; /// 座標系の移動速度（無次元）

public:
  /// 参照系の定義
  enum frame_type {
    frm_static,
    frm_translation,
    frm_rotation
  };
  
  ReferenceFrame(){
    Frame     = 0;
    TimeAccel = 0.0;
    v00[0] = v00[1] = v00[2] = v00[3] = 0.0;
    GridVel[0] = GridVel[1] = GridVel[2] = 0.0;
  }
  ~ReferenceFrame() {}

  void setAccel  (const double m_timeAccel);
  void setFrame  (const unsigned m_frame);
  void setGridVel(const double* m_Gvel);
  void setV00    (const double time, const bool init=false);
  
  //@fn unsigned getFrame(void) const
  //@brief Frameを返す
  unsigned getFrame(void) const {
    return Frame;
  }
  
  //@fn double getAccel(void) const
  //@brief Frameを返す
  double getAccel(void) const {
    return TimeAccel;
  }
  
  //@fn void copyV00(double* m_v0) const
  //@brief v00をコピーする
  void copyV00(double* m_v0) const {
    for (int i=0; i<4; i++) m_v0[i] = v00[i];
  }
  
  //@fn void copyGridVel(double* m_gv) const
  //@brief GridVelをコピーする
  void copyGridVel(double* m_gv) const {
    for (int i=0; i<3; i++) m_gv[i] = GridVel[i];
  }
};

#endif // _SKL_FB_REF_FRAME_H_
