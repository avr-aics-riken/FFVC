#ifndef _FB_FRACTION_H_
#define _FB_FRACTION_H_

//##################################################################################
//
// Flow Base class
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2014 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

/** 
 * @file   CompoFraction.h
 * @brief  Component Fraction class Header
 * @author aics
 */

#include "FB_Define.h"
#include "Vec3.h"
#include <stdio.h>

using namespace Vec3class;

class CompoFraction {
  
protected:
  // ローカル計算領域情報
  int size[3];     ///< セル(ローカル)サイズ
  int guide;       ///< ガイドセル数
  int division;    ///< 細分化の分割数
  Vec3r pch;       ///< セル幅
  Vec3r org;       ///< 計算領域の基点
  Vec3r angle;     ///< 変換の回転角度
  
  // 形状パラメータ
  int smode;         ///< 形状モード
  REAL_TYPE depth;   ///< 厚さ
  REAL_TYPE width;   ///< 矩形の幅（dir方向)
  REAL_TYPE height;  ///< 矩形の高さ
  REAL_TYPE r_fan;   ///< ファン半径
  REAL_TYPE r_boss;  ///< ボス半径
  Vec3r nv;          ///< 法線方向ベクトル（流出方向）
  Vec3r center;      ///< 形状の中心座標（前面の中心位置）
  Vec3r dir;         ///< 矩形の方向規定の参照ベクトル
  Vec3r box_min;     ///< Bounding boxの最小値
  Vec3r box_max;     ///< Bounding boxの最大値
  
public:
  /** デフォルトコンストラクタ */
  CompoFraction() {}
  
  /** コンストラクタ
   * @param [in] size   ローカルセル数
   * @param [in] guide  ガイドセル数
   * @param [in] pch    ローカル領域セル幅
   * @param [in] org    ローカル領域基点座標
   * @param [in] div    サブディビジョンの分割数
   */
  CompoFraction(const int* size, const int guide, const REAL_TYPE* pch, const REAL_TYPE* org, const int div)
  {
    this->size[0]  = size[0];
    this->size[1]  = size[1];
    this->size[2]  = size[2];
    this->guide    = guide;
    this->pch.x    = pch[0];
    this->pch.y    = pch[1];
    this->pch.z    = pch[2];
    this->org.x    = org[0];
    this->org.y    = org[1];
    this->org.z    = org[2];
    this->division = div;
  }
  
  /** デストラクタ */
  virtual ~CompoFraction() {}
  
protected:
  /**
   * @brief 矩形領域のbboxを計算
   * @param [in] mn bboxの最小位置
   * @param [in] mx bboxの最大位置
   * @erturn 投影面積
   */
  REAL_TYPE bbox_rect_cylinder(Vec3r& mn, Vec3r& mx);
  
  
  /**
   * @brief 円筒領域のbboxを計算
   * @param [in] mn bboxの最小位置
   * @param [in] mx bboxの最大位置
   * @erturn 投影面積
   */
  REAL_TYPE bbox_circ_cylinder(Vec3r& mn, Vec3r& mx);
  
  
  /**
   * @brief 点pの属するセルインデクスを求める
   * @param [out] w インデクス
   * @param [in]  p テストする座標値
   * @note Fortran index
   */
  void find_index(int* w, const Vec3r p);
  
  
  /** インデックスを(1,0,0)シフト
   * @param [in] index 元のインデクス
   * @param [in] h     シフト幅
   */
  inline Vec3r shift_f1 (const Vec3r index, const REAL_TYPE h) 
  { 
    return Vec3r(index.x+h, index.y  , index.z  ); 
  }
  
  /** インデックスを(0,1,0)シフト
   * @param [in] index 元のインデクス
   * @param [in] h     シフト幅
   */
  inline Vec3r shift_f2 (const Vec3r index, const REAL_TYPE h) 
  { 
    return Vec3r(index.x  , index.y+h, index.z  ); 
  }
  
  /** インデックスを(1,1,0)シフト
   * @param [in] index 元のインデクス
   * @param [in] h     シフト幅
   */
  inline Vec3r shift_f3 (const Vec3r index, const REAL_TYPE h) 
  { 
    return Vec3r(index.x+h, index.y+h, index.z  ); 
  }
  
  /** インデックスを(0,0,1)シフト
   * @param [in] index 元のインデクス
   * @param [in] h     シフト幅
   */
  inline Vec3r shift_f4 (const Vec3r index, const REAL_TYPE h) 
  {
    return Vec3r(index.x  , index.y  , index.z+h); 
  }
  
  /** インデックスを(1,0,1)シフト
   * @param [in] index 元のインデクス
   * @param [in] h     シフト幅
   */
  inline Vec3r shift_f5 (const Vec3r index, const REAL_TYPE h) 
  { 
    return Vec3r(index.x+h, index.y  , index.z+h); 
  }
  
  /** インデックスを(0,1,1)シフト
   * @param [in] index 元のインデクス
   * @param [in] h     シフト幅
   */
  inline Vec3r shift_f6 (const Vec3r index, const REAL_TYPE h)
  { 
    return Vec3r(index.x  , index.y+h, index.z+h); 
  }
  
  /** インデックスを(1,1,1)シフト
   * @param [in] index 元のインデクス
   * @param [in] h     シフト幅
   */
  inline Vec3r shift_f7 (const Vec3r index, const REAL_TYPE h) 
  {
    return Vec3r(index.x+h, index.y+h, index.z+h); 
  }

  
  /**
   * @brief 円筒形状の内外判定
   * @param [in] p テスト点座標
   * @return 内部のときに1.0を返す
   * @note 186 flop
   */
  inline REAL_TYPE judge_cylinder(const Vec3r p) 
  {
    Vec3r q = rotate(angle, p-center); // 181 flop
    
    if ( (q.z < 0.0) || (q.z > depth)  ) return 0.0;
    REAL_TYPE r = sqrtf(q.x*q.x + q.y*q.y);
    
    return ( (r<=r_fan) && (r>=r_boss) ) ? 1.0 : 0.0;
  }
  
  
  /**
   * @brief 矩形形状の内外判定
   * @param [in] p テスト点座標
   * @return 内部のときに1.0を返す
   * @note 183 flop
   */
  inline REAL_TYPE judge_rect(const Vec3r p) 
  {
    Vec3r q = rotate(angle, p-center); // 181 flop
    
    if ( (q.z < 0.0) || (q.z > depth)  ) return 0.0;
    if ( fabs(q.x) > 0.5*width )  return 0.0;
    if ( fabs(q.y) > 0.5*height ) return 0.0;
    
    return 1.0;
  }
  
  
  /**
   * @brief 回転ベクトルp(alpha, beta, gamma)でベクトルuを回転する
   * @param [in] p 回転角度
   * @param [in] u 方向ベクトル
   * @return 角度
   * @note sin, cos = 5flop, dot=5flop, total 181flop
   */
  inline Vec3r rotate(const Vec3r p, const Vec3r u)
  {
    Vec3r a, b, c;
    
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
    
    return Vec3r( dot(a, u), dot(b, u), dot(c, u) );
  }
  
  
  /**
   * @brief 回転ベクトルp(alpha, beta, gamma)に対して，-pでベクトルuを回転する
   * @param [in] p 回転角度
   * @param [in] u 方向ベクトル
   * @return 角度
   * @note sin, cos = 5flop, dot=5flop, total 181flop
   */
  inline Vec3r rotate_inv(const Vec3r p, const Vec3r u)
  {
    Vec3r a, b, c;
    
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
    
    return Vec3r( dot(a, u), dot(b, u), dot(c, u) );
  }
  
  
public:
  
  /** 
   * @brief 形状のbboxと投影面積を求める
   * @return 投影面積
   */
  REAL_TYPE get_BboxArea ();
  
  
  /**
   * @brief コンポーネントの属するセルインデクスを求める
   * @param [out] st 開始インデクス
   * @param [out] ed 終了インデクス
   */
  void bbox_index(int* st, int* ed);
  
  /**
   * @brief 指定法線nvがz軸の方向ベクトルに向かう回転角を計算する
   */
  void get_angle();
  
  
  /**
   * @brief 矩形の形状パラメータをセットする
   * @param [in] m_nv     法線ベクトル
   * @param [in] m_ctr    中心座標
   * @param [in] m_dir    方向ベクトル
   * @param [in] m_depth  厚み
   * @param [in] m_width  幅
   * @param [in] m_height 高さ
   */
  void setShapeParam(const REAL_TYPE m_nv[3],
                     const REAL_TYPE m_ctr[3],
                     const REAL_TYPE m_dir[3],
                     const REAL_TYPE m_depth,
                     const REAL_TYPE m_width,
                     const REAL_TYPE m_height);
  
  /**
   * @brief 円筒の形状パラメータをセットする
   * @param [in] m_nv     法線ベクトル
   * @param [in] m_ctr    中心座標
   * @param [in] m_depth  厚み
   * @param [in] m_r_fan  外径
   * @param [in] m_r_boss 内径
   */
  void setShapeParam(const REAL_TYPE m_nv[3],
                     const REAL_TYPE m_ctr[3],
                     const REAL_TYPE m_depth,
                     const REAL_TYPE m_r_fan,
                     const REAL_TYPE m_r_boss=0.0f);
  
  
  /**
   * @brief 体積率が(0,1)の間のセルに対してサブディビジョンを実施
   * @param [in]     st    開始インデクス
   * @param [in]     ed    終了インデクス
   * @param [in,out] vf    フラクション
   * @param [in,out] flop  浮動小数点演算数
   */
  void subdivision(const int st[], const int ed[], float* vf, double& flop);
  
  
  /**
   * @brief セルの8頂点の内外判定を行い，0, 1, otherに分類
   * @param [in]     st    開始インデクス
   * @param [in]     en    終了インデクス
   * @param [in,out] vf    フラクション
   * @param [in,out] flop  浮動小数点演算数
   */
  void vertex8(const int st[], const int ed[], float* vf, double& flop);
  
  
  /**
   * @brief ベクトルの最小成分
   * @param [in,out] mn 比較して小さい成分
   * @param [in]     p  参照ベクトル
   */
  static inline void get_min(Vec3r& mn, const Vec3r p)
  {
    mn.x = (mn.x < p.x) ? mn.x : p.x;
    mn.y = (mn.y < p.y) ? mn.y : p.y;
    mn.z = (mn.z < p.z) ? mn.z : p.z;
  }
  
  
  /**
   * @brief ベクトルの最大値成分
   * @param [in,out] mx 比較して大きい成分
   * @param [in]     p  参照ベクトル
   */
  static inline void get_max(Vec3r& mx, const Vec3r p)
  {
    mx.x = (mx.x > p.x) ? mx.x : p.x;
    mx.y = (mx.y > p.y) ? mx.y : p.y;
    mx.z = (mx.z > p.z) ? mx.z : p.z;
  }
  
};



class ShapeMonitor : public CompoFraction {

public:
  /** デフォルトコンストラクタ */
  ShapeMonitor() {}
  
  
  /** コンストラクタ
   * @param [in] size   ローカルセル数
   * @param [in] guide  ガイドセル数
   * @param [in] pch    ローカル領域セル幅
   * @param [in] org    ローカル領域基点座標
   */
  ShapeMonitor(const int* size, const int guide, const REAL_TYPE* pch, const REAL_TYPE* org) 
  {
    this->size[0]  = size[0];
    this->size[1]  = size[1];
    this->size[2]  = size[2];
    this->guide    = guide;
    this->pch.x    = pch[0];
    this->pch.y    = pch[1];
    this->pch.z    = pch[2];
    this->org.x    = org[0];
    this->org.y    = org[1];
    this->org.z    = org[2];
  }
  
  
   /** デストラクタ */
  virtual ~ShapeMonitor() {}
  
public:
  /**
   * @brief セルの8頂点の内外判定より50%以上のセルにIDを設定する
   * @param [in]     st   開始インデクス
   * @param [in]     en   終了インデクス
   * @param [in,out] bcd  BCindex B
   * @param [in]     id   指定ID
   */
  void setID(const int st[], const int ed[], int* bcd, const int m_id);
};


#endif // _FB_FRACTION_H_
