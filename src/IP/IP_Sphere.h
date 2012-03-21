#ifndef _SKL_IP_SHERE_H_
#define _SKL_IP_SHERE_H_

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
#include "vec3.h"

class IP_Sphere : public Intrinsic {
protected:
  REAL_TYPE offset;     ///< オフセット距離
  REAL_TYPE radius;     ///< 球の半径
  REAL_TYPE drv_length; ///< ドライバの長さ
  FB::Vec3f pch;        ///< セル幅
  FB::Vec3f org;        ///< 計算領域の基点
  FB::Vec3f wth;        ///< 計算領域の大きさ
  FB::Vec3f ctr;        ///< 計算領域のセンター
  FB::Vec3f os;         ///< 球の中心座標
  FB::Vec3i sz;         ///< 領域分割数
  FB::Vec3f box_min;    ///< Bounding boxの最小値
  FB::Vec3f box_max;    ///< Bounding boxの最大値
  
public:
  IP_Shere(){
    offset = 0.0;
    radius = 0.0;
    drv_length = 0.0;
  }
  ~IP_Shere() {}

public:
  virtual bool getXML(SklSolverConfig* CF, Control* R);
  
  virtual void setDomain(Control* R, unsigned sz[3], REAL_TYPE org[3], REAL_TYPE wth[3], REAL_TYPE pch[3]);
  virtual void setup(int* mid, Control* R, REAL_TYPE* G_org);
  virtual void printPara(FILE* fp, Control* R);
  
  virtual const char* getExampleName(void) {
    return ("Sphere");
  }
};
#endif // _SKL_IP_SHERE_H_
