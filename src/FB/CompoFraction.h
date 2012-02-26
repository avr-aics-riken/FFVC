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


#include "FBDefine.h"
#include "FBUtility.h"
#include "vec3.h"

class CompoFraction {
public:
  /// サンプリング形状
  enum Shapes {
    RECT,      ///< 矩形領域
    CYLINDER,  ///< 円筒領域
  };
  
protected:
  // ローカル計算領域情報
  int size[3];     ///< セル(ローカル)サイズ
  int guide;       ///< ガイドセル数
  int division;    ///< 細分化の分割数
  FB::Vec3f pch;   ///< セル幅
  FB::Vec3f org;   ///< 計算領域の基点
  FB::Vec3f angle; ///< アフィン変換の回転角度
  
  // 形状パラメータ
  Shapes smode;    ///< 形状モード
  float depth;     ///< 厚さ
  float width;     ///< 矩形の幅（dir方向)
  float height;    ///< 矩形の高さ
  float radius;    ///< 半径
  FB::Vec3f nv;    ///< 法線方向ベクトル（流出方向）
  FB::Vec3f center;///< 形状の中心座標（前面の中心位置）
  FB::Vec3f dir;   ///< 矩形の方向規定の参照ベクトル
  
public:
  /// デフォルトコンストラクタ
  CompoFraction() {}
  
  /// コンストラクタ
  ///   @param[in] size,guide ローカルセル数，ガイドセル数
  ///   @param[in] crd  モニタ点座標
  ///   @param[in] org,pch  ローカル領域基点座標，セル幅
  ///   @param[in] div サブディビジョンの分割数
  CompoFraction(unsigned size[], unsigned guide, FB::Vec3f pch, FB::Vec3f org, int div) {
    this->size[0]  = (int)size[0];
    this->size[1]  = (int)size[1];
    this->size[2]  = (int)size[2];
    this->guide    = (int)guide;
    this->pch      = pch;
    this->org      = org;
    this->division = div;
  }
  
  /// デストラクタ
  ~CompoFraction() {}
  
protected:
  void get_angle       (void);
  void subdivision     (int st[], int ed[], float* vf);
  void vertex8         (int st[], int ed[], float* vf);
  
  FB::Vec3f transform  (const FB::Vec3f angle, const FB::Vec3f u);
  
  /// セルインデックスを(1,0,0)シフト
  FB::Vec3f shift_f1(const FB::Vec3f index, const float h) { return FB::Vec3f(index.x+h, index.y  , index.z  ); }
  
  /// セルインデックスを(0,1,0)シフト
  FB::Vec3f shift_f2(const FB::Vec3f index, const float h) { return FB::Vec3f(index.x  , index.y+h, index.z  ); }
  
  /// セルインデックスを(1,1,0)シフト
  FB::Vec3f shift_f3(const FB::Vec3f index, const float h) { return FB::Vec3f(index.x+h, index.y+h, index.z  ); }
  
  /// セルインデックスを(0,0,1)シフト
  FB::Vec3f shift_f4(const FB::Vec3f index, const float h) { return FB::Vec3f(index.x  , index.y  , index.z+h); }
  
  /// セルインデックスを(1,0,1)シフト
  FB::Vec3f shift_f5(const FB::Vec3f index, const float h) { return FB::Vec3f(index.x+h, index.y  , index.z+h); }
  
  /// セルインデックスを(0,1,1)シフト
  FB::Vec3f shift_f6(const FB::Vec3f index, const float h) { return FB::Vec3f(index.x  , index.y+h, index.z+h); }
  
  /// セルインデックスを(1,1,1)シフト
  FB::Vec3f shift_f7(const FB::Vec3f index, const float h) { return FB::Vec3f(index.x+h, index.y+h, index.z+h); }
  
  //@fn inline float judge_cylinder(const FB::Vec3f p)
  //@brief 円筒形状の内外判定
  //@param p テスト点座標
  //@ret 内部のときに1.0を返す
  inline float judge_cylinder(const FB::Vec3f p) {
    FB::Vec3f q = transform(angle, p);
    
    if (q.z > depth) return 0.0;
    float r = sqrtf(q.x*q.x + q.y*q.y);
    
    return (r<=radius) ? 1.0 : 0.0;
  }
  
  //@fn inline float judge_rect(const FB::Vec3f p)
  //@brief 矩形形状の内外判定
  //@param p テスト点座標
  //@ret 内部のときに1.0を返す
  inline float judge_rect(const FB::Vec3f p) {
    FB::Vec3f q = transform(angle, p);
    
    if (q.z > depth) return 0.0;
    if (q.x > fabs(0.5*width)) return 0.0;
    if (q.y > fabs(0.5*height)) return 0.0;
    
    return 1.0;
  }
  
public:
  void get_fraction  (int st[], int ed[], float* vf);
  void setShapeParam (FB::Vec3f m_nv, FB::Vec3f m_ctr, FB::Vec3f m_dir, float m_depth, float m_width, float m_height);
  void setShapeParam (FB::Vec3f m_nv, FB::Vec3f m_ctr, float m_depth, float m_radius);
  
};

#endif // _SKL_FB_FRACTION_H_
