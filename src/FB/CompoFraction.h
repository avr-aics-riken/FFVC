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

// テンポラリに使用
#define SKL_REAL float

#include "vec3.h"
#include "FBDefine.h"

class CompoFraction {
public:
  /// サンプリング形状
  enum Shapes {
    RECT,      ///< 矩形領域
    CYLINDER,  ///< 円筒領域
  };
  
protected:
  // ローカル計算領域情報
  unsigned size[3];  ///< セル(ローカル)サイズ
  unsigned guide;    ///< ガイドセル数
  Vec3f pch;         ///< セル幅
  Vec3f org;         ///< 基点
   
  int division;      ///< 細分化の分割数
  
  /// 形状パラメータ 各コンポーネント毎に異なる
  Shapes smode; ///< 形状モード
  Vec3f dir;    ///< 矩形の方向規定の参照ベクトル
  Vec3f nv;     ///< 法線方向ベクトル（流出方向）
  Vec3f ctr;    ///< 前面の中心座標
  float depth;  ///< 厚さ
  float radius; ///< 円筒の半径
  float width;  ///< 矩形の幅（dir方向）
  float height; ///< 矩形の高さ
  
  Vec3f angle;  ///< アフィン変換の回転角度
  
public:
  /// デフォルトコンストラクタ
  CompoFraction() {}
  
  /// コンストラクタ
  ///   @param[in] size,guide ローカルセル数，ガイドセル数
  ///   @param[in] crd  モニタ点座標
  ///   @param[in] org,pch  ローカル領域基点座標，セル幅
  ///   @param[in] div サブディビジョンの分割数
  CompoFraction(unsigned size[], unsigned guide, Vec3f pch, Vec3f org, int div) {
    this->size[0]  = size[0];
    this->size[1]  = size[1];
    this->size[2]  = size[2];
    this->guide    = guide;
    this->pch      = pch;
    this->org      = org;
    this->division = div;
  }
  
  /// デストラクタ
  ~CompoFraction() {}
  
protected:
  void get_angle(void);
  void subdivision     (int st[], int ed[], float* vf);
  void vertex8         (int st[], int ed[], float* vf);
  
  float judge_cylinder (Vec3f p);
  float judge_rect     (Vec3f p);
  
  Vec3f transform      (const Vec3f angle, const Vec3f u);
  
  /// セルインデックスを(1,0,0)シフト
  Vec3f shift_f1(Vec3f index, float h) { return Vec3f(index.x+h, index.y  , index.z  ); }
  
  /// セルインデックスを(0,1,0)シフト
  Vec3f shift_f2(Vec3f index, float h) { return Vec3f(index.x  , index.y+h, index.z  ); }
  
  /// セルインデックスを(1,1,0)シフト
  Vec3f shift_f3(Vec3f index, float h) { return Vec3f(index.x+h, index.y+h, index.z  ); }
  
  /// セルインデックスを(0,0,1)シフト
  Vec3f shift_f4(Vec3f index, float h) { return Vec3f(index.x  , index.y  , index.z+h); }
  
  /// セルインデックスを(1,0,1)シフト
  Vec3f shift_f5(Vec3f index, float h) { return Vec3f(index.x+h, index.y  , index.z+h); }
  
  /// セルインデックスを(0,1,1)シフト
  Vec3f shift_f6(Vec3f index, float h) { return Vec3f(index.x  , index.y+h, index.z+h); }
  
  /// セルインデックスを(1,1,1)シフト
  Vec3f shift_f7(Vec3f index, float h) { return Vec3f(index.x+h, index.y+h, index.z+h); }
  
public:
  void get_fraction(int st[], int ed[], float* vf);
  
};

#endif // _SKL_FB_FRACTION_H_
