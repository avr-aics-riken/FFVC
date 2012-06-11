#ifndef _IP_SHERE_H_
#define _IP_SHERE_H_

/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file IP_Sphere.h
//@brief IP_Sphere class Header
//@author keno, FSI Team, VCAD, RIKEN

#include "Intrinsic.h"
#include "IP_Define.h"

class IP_Sphere : public Intrinsic {
protected:
  REAL_TYPE radius;     ///< 球の半径
  REAL_TYPE drv_length; ///< ドライバの長さ
  int drv_mode;         ///< ドライバのON/OFF
  FB::Vec3f pch;        ///< セル幅
  FB::Vec3f org;        ///< 計算領域の基点
  FB::Vec3f wth;        ///< 計算領域の大きさ
  FB::Vec3f box_min;    ///< Bounding boxの最小値
  FB::Vec3f box_max;    ///< Bounding boxの最大値
  FB::Vec3i box_st;     ///< Bounding boxの始点インデクス
  FB::Vec3i box_ed;     ///< Bounding boxの終点インデクス
  
  std::string m_fluid;
  std::string m_solid;
  std::string m_driver;
  std::string m_driver_face;
  
public:
  IP_Sphere(){
    radius = 0.0;
    drv_length = 0.0;
    drv_mode = OFF;
  }
  ~IP_Sphere() {}

public:
  virtual bool getTP(Control* R, TPControl* tpCntl);
  
  virtual void setDomain(Control* R, unsigned sz[3], REAL_TYPE org[3], REAL_TYPE wth[3], REAL_TYPE pch[3]);
  virtual void printPara(FILE* fp, Control* R);
  virtual const char* getExampleName(void) {
    return ("Sphere");
  }
  
  virtual void setup(int* mid, Control* R, REAL_TYPE* G_org);
  virtual void setup_cut(int* mid, Control* R, REAL_TYPE* G_org, float* cut);
  
  FB::Vec3i find_index(const FB::Vec3f p, const FB::Vec3f ol);
  
  float cut_line(const FB::Vec3f b, const int dir, const float r, const float dh);
};
#endif // _IP_SHERE_H_
