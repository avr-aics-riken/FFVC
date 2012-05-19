#ifndef _FB_FRACTION_H_
#define _FB_FRACTION_H_

/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file CompoFraction.h
//@brief Component Fraction class Header
//@author keno, AICS, RIKEN


#include "FB_Define.h"
#include "FBUtility.h"
#include "vec3.h"
#include <stdio.h>

class CompoFraction {
  
protected:
  // ローカル計算領域情報
  int size[3];     ///< セル(ローカル)サイズ
  int guide;       ///< ガイドセル数
  int division;    ///< 細分化の分割数
  FB::Vec3f pch;   ///< セル幅
  FB::Vec3f org;   ///< 計算領域の基点
  FB::Vec3f angle; ///< 変換の回転角度
  
  // 形状パラメータ
  int smode;         ///< 形状モード
  float depth;       ///< 厚さ
  float width;       ///< 矩形の幅（dir方向)
  float height;      ///< 矩形の高さ
  float r_fan;       ///< ファン半径
  float r_boss;      ///< ボス半径
  FB::Vec3f nv;      ///< 法線方向ベクトル（流出方向）
  FB::Vec3f center;  ///< 形状の中心座標（前面の中心位置）
  FB::Vec3f dir;     ///< 矩形の方向規定の参照ベクトル
  FB::Vec3f box_min; ///< Bounding boxの最小値
  FB::Vec3f box_max; ///< Bounding boxの最大値
  
public:
  /// デフォルトコンストラクタ
  CompoFraction() {}
  
  /// コンストラクタ
  ///   @param[in] size,guide ローカルセル数，ガイドセル数
  ///   @param[in] crd  モニタ点座標
  ///   @param[in] org,pch  ローカル領域基点座標，セル幅
  ///   @param[in] div サブディビジョンの分割数
  CompoFraction(unsigned size[], unsigned guide, float pch[], float org[], int div) {
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
  float bbox_rect_cylinder(FB::Vec3f& mn, FB::Vec3f& mx);
  float bbox_circ_cylinder(FB::Vec3f& mn, FB::Vec3f& mx);
  
  void find_index(int* w, const FB::Vec3f p);
  
  inline FB::Vec3f shift_f1 (const FB::Vec3f index, const float h) { return FB::Vec3f(index.x+h, index.y  , index.z  ); } /// インデックスを(1,0,0)シフト
  inline FB::Vec3f shift_f2 (const FB::Vec3f index, const float h) { return FB::Vec3f(index.x  , index.y+h, index.z  ); } /// インデックスを(0,1,0)シフト
  inline FB::Vec3f shift_f3 (const FB::Vec3f index, const float h) { return FB::Vec3f(index.x+h, index.y+h, index.z  ); } /// インデックスを(1,1,0)シフト
  inline FB::Vec3f shift_f4 (const FB::Vec3f index, const float h) { return FB::Vec3f(index.x  , index.y  , index.z+h); } /// インデックスを(0,0,1)シフト
  inline FB::Vec3f shift_f5 (const FB::Vec3f index, const float h) { return FB::Vec3f(index.x+h, index.y  , index.z+h); } /// インデックスを(1,0,1)シフト
  inline FB::Vec3f shift_f6 (const FB::Vec3f index, const float h) { return FB::Vec3f(index.x  , index.y+h, index.z+h); } /// インデックスを(0,1,1)シフト
  inline FB::Vec3f shift_f7 (const FB::Vec3f index, const float h) { return FB::Vec3f(index.x+h, index.y+h, index.z+h); } /// インデックスを(1,1,1)シフト
  
  //@fn inline void get_min()
  inline void get_min(FB::Vec3f& mn, const FB::Vec3f p) {
    mn.x = (mn.x < p.x) ? mn.x : p.x;
    mn.y = (mn.y < p.y) ? mn.y : p.y;
    mn.z = (mn.z < p.z) ? mn.z : p.z;
  }
  
  //@fn inline void get_max()
  inline void get_max(FB::Vec3f& mx, const FB::Vec3f p) {
    mx.x = (mx.x > p.x) ? mx.x : p.x;
    mx.y = (mx.y > p.y) ? mx.y : p.y;
    mx.z = (mx.z > p.z) ? mx.z : p.z;
  }
  
  //@fn inline float judge_cylinder(const FB::Vec3f p)
  //@brief 円筒形状の内外判定
  //@param p テスト点座標
  //@ret 内部のときに1.0を返す
  //@note 186 flop
  inline float judge_cylinder(const FB::Vec3f p) {
    FB::Vec3f q = rotate(angle, p-center); // 181 flop
    
    if ( (q.z < 0.0) || (q.z > depth)  ) return 0.0;
    float r = sqrtf(q.x*q.x + q.y*q.y);
    
    return ( (r<=r_fan) && (r>=r_boss) ) ? 1.0 : 0.0;
  }
  
  //@fn inline float judge_rect(const FB::Vec3f p)
  //@brief 矩形形状の内外判定
  //@param p テスト点座標
  //@ret 内部のときに1.0を返す
  //@note 183 flop
  inline float judge_rect(const FB::Vec3f p) {
    FB::Vec3f q = rotate(angle, p-center); // 181 flop
    
    if ( (q.z < 0.0) || (q.z > depth)  ) return 0.0;
    if ( fabs(q.x) > 0.5*width )  return 0.0;
    if ( fabs(q.y) > 0.5*height ) return 0.0;
    
    return 1.0;
  }
  
  //@fn inline FB::Vec3f CompoFraction::rotate(const Vec3f p, const Vec3f u)
  //@brief 回転ベクトルp(alpha, beta, gamma)でベクトルuを回転する
  //@param p 回転角度
  //@param u 方向ベクトル
  //@ret 角度
  //@note sin, cos = 5flop, dot=5flop, total 181flop
  inline FB::Vec3f rotate(const FB::Vec3f p, const FB::Vec3f u)
  {
    FB::Vec3f a, b, c;
    
    // line vector expression
    a.x =  cos(p.y)*cos(p.z);
    a.y =  sin(p.x)*sin(p.y)*cos(p.z) - cos(p.x)*sin(p.z);
    a.z =  cos(p.x)*sin(p.y)*cos(p.z) + sin(p.x)*sin(p.z);
    
    b.x =  cos(p.y)*sin(p.z);
    b.y =  sin(p.x)*sin(p.y)*sin(p.z) + cos(p.x)*cos(p.z);
    b.z =  cos(p.x)*sin(p.y)*sin(p.z) - sin(p.x)*cos(p.z);
    
    c.x = -sin(p.y);
    c.y =  sin(p.x)*cos(p.y);
    c.z =  cos(p.x)*cos(p.y);
    
    return FB::Vec3f( dot(a, u), dot(b, u), dot(c, u) );
  }
  
  //@fn inline FB::Vec3f CompoFraction::rotate_inv(const Vec3f p, const Vec3f u)
  //@brief 回転ベクトルp(alpha, beta, gamma)に対して，-pでベクトルuを回転する
  //@param p 回転角度
  //@param u 方向ベクトル
  //@ret 角度
  //@note sin, cos = 5flop, dot=5flop, total 181flop
  inline FB::Vec3f rotate_inv(const FB::Vec3f p, const FB::Vec3f u)
  {
    FB::Vec3f a, b, c;
    
    // line vector expression
    a.x =  cos(p.y)*cos(p.z);
    a.y =  cos(p.y)*sin(p.z);
    a.z = -sin(p.y);
    
    b.x =  sin(p.x)*sin(p.y)*cos(p.z) - cos(p.x)*sin(p.z);
    b.y =  sin(p.x)*sin(p.y)*sin(p.z) + cos(p.x)*cos(p.z);
    b.z =  sin(p.x)*cos(p.y);
    
    c.x =  cos(p.x)*sin(p.y)*cos(p.z) + sin(p.x)*sin(p.z);
    c.y =  cos(p.x)*sin(p.y)*sin(p.z) - sin(p.x)*cos(p.z);
    c.z =  cos(p.x)*cos(p.y);
    
    return FB::Vec3f( dot(a, u), dot(b, u), dot(c, u) );
  }
  
public:
  float get_BboxArea (void);
  
  void bbox_index(int* st, int* ed);
  void get_angle     (void);
  void setShapeParam (const float m_nv[3], const float m_ctr[3], const float m_dir[3], const float m_depth, const float m_width, const float m_height);
  void setShapeParam (const float m_nv[3], const float m_ctr[3], const float m_depth, const float m_r_fan, const float m_r_boss=0.0f);
  void subdivision   (const int st[], const int ed[], float* vf, double& flop);
  void vertex8       (const int st[], const int ed[], float* vf, double& flop);
};


class ShapeMonitor : public CompoFraction {

public:
  ShapeMonitor() {}
  
  /// コンストラクタ
  ///   @param[in] size,guide ローカルセル数，ガイドセル数
  ///   @param[in] crd  モニタ点座標
  ///   @param[in] org,pch  ローカル領域基点座標，セル幅
  ShapeMonitor(unsigned size[], unsigned guide, float pch[], float org[]) {
    this->size[0]  = (int)size[0];
    this->size[1]  = (int)size[1];
    this->size[2]  = (int)size[2];
    this->guide    = (int)guide;
    this->pch      = pch;
    this->org      = org;
  }
  
  ~ShapeMonitor() {}
  
public:
  void setID(const int st[], const int ed[], int* mid, int m_id);
};


#endif // _FB_FRACTION_H_
