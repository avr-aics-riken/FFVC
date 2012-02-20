#ifndef _SKL_FB_FRACTION_H_
#define _SKL_FB_FRACTION_H_

/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file CompoFraction.h
//@brief Component Fraction class Header
//@author keno, AICS, RIKEN

//#include "SklUtil.h"
//#include "vec3.h"
#include "basic_func.h"
#include "vec3.h"
#include "FBDefine.h"
#include "math.h"

class CompoFraction {
public:
  
  
public:
  /// ディフォルトコンストラクタ
  CompoFraction() {}
  
  
  /// デストラクタ
  virtual ~CompoFraciton() {}
  
  
protected:
  Vec3f get_alpha_beta(const Vec3f w);
  Vec3f transform(const Vec3f angle, const Vec3f u);
  
  //@fn get_gamma()
  //@brief アフィン変換の角度γを求める
  void get_gamma();
  
  //@fn trans_rotate()
  //@brief 回転変換
  Vec3f trans_rotate();
  
  //@fn Vec3f trans_shift(const Ve3f a, const Vec3f b)
  //@brief bが原点となる平行移動変換
  Vec3f trans_shift(const Ve3f a, const Vec3f b) {
    return (a-b);
  }
  
};

#endif // _SKL_FB_FRACTION_H_
