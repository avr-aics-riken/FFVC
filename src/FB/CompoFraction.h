#ifndef _FB_FRACTION_H_
#define _FB_FRACTION_H_

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
 * @file   CompoFraction.h
 * @brief  Component Fraction class Header
 * @author kero
 */

#include <stdio.h>
#include "FB_Define.h"
#include "vec3.h"


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
  /** デフォルトコンストラクタ */
  CompoFraction() {}
  
  /** コンストラクタ
   * @param [in] size   ローカルセル数
   * @param [in] guide  ガイドセル数
   * @param [in] crd    モニタ点座標
   * @param [in] org    ローカル領域基点座標
   * @param [in] pch    ローカル領域セル幅
   * @param [in] div    サブディビジョンの分割数
   */
  CompoFraction(const int* size, const int guide, const float* pch, const float* org, const int div)
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
  float bbox_rect_cylinder(FB::Vec3f& mn, FB::Vec3f& mx);
  
  
  /**
   * @brief 円筒領域のbboxを計算
   * @param [in] mn bboxの最小位置
   * @param [in] mx bboxの最大位置
   * @erturn 投影面積
   */
  float bbox_circ_cylinder(FB::Vec3f& mn, FB::Vec3f& mx);
  
  
  /**
   * @brief 点pの属するセルインデクスを求める
   * @param [out] w インデクス
   * @param [in]  p テストする座標値
   * @note Fortran index
   */
  void find_index(int* w, const FB::Vec3f p);
  
  
  /** インデックスを(1,0,0)シフト
   * @param [in] index 元のインデクス
   * @param [in] h     シフト幅
   */
  inline FB::Vec3f shift_f1 (const FB::Vec3f index, const float h) 
  { 
    return FB::Vec3f(index.x+h, index.y  , index.z  ); 
  }
  
  /** インデックスを(0,1,0)シフト
   * @param [in] index 元のインデクス
   * @param [in] h     シフト幅
   */
  inline FB::Vec3f shift_f2 (const FB::Vec3f index, const float h) 
  { 
    return FB::Vec3f(index.x  , index.y+h, index.z  ); 
  }
  
  /** インデックスを(1,1,0)シフト
   * @param [in] index 元のインデクス
   * @param [in] h     シフト幅
   */
  inline FB::Vec3f shift_f3 (const FB::Vec3f index, const float h) 
  { 
    return FB::Vec3f(index.x+h, index.y+h, index.z  ); 
  }
  
  /** インデックスを(0,0,1)シフト
   * @param [in] index 元のインデクス
   * @param [in] h     シフト幅
   */
  inline FB::Vec3f shift_f4 (const FB::Vec3f index, const float h) 
  {
    return FB::Vec3f(index.x  , index.y  , index.z+h); 
  }
  
  /** インデックスを(1,0,1)シフト
   * @param [in] index 元のインデクス
   * @param [in] h     シフト幅
   */
  inline FB::Vec3f shift_f5 (const FB::Vec3f index, const float h) 
  { 
    return FB::Vec3f(index.x+h, index.y  , index.z+h); 
  }
  
  /** インデックスを(0,1,1)シフト
   * @param [in] index 元のインデクス
   * @param [in] h     シフト幅
   */
  inline FB::Vec3f shift_f6 (const FB::Vec3f index, const float h)
  { 
    return FB::Vec3f(index.x  , index.y+h, index.z+h); 
  }
  
  /** インデックスを(1,1,1)シフト
   * @param [in] index 元のインデクス
   * @param [in] h     シフト幅
   */
  inline FB::Vec3f shift_f7 (const FB::Vec3f index, const float h) 
  {
    return FB::Vec3f(index.x+h, index.y+h, index.z+h); 
  }
  
  
  /**
   * @brief ベクトルの最小成分
   * @param [in/out] mn 比較して小さい成分
   * @param [in]     p  参照ベクトル
   */
  inline void get_min(FB::Vec3f& mn, const FB::Vec3f p) 
  {
    mn.x = (mn.x < p.x) ? mn.x : p.x;
    mn.y = (mn.y < p.y) ? mn.y : p.y;
    mn.z = (mn.z < p.z) ? mn.z : p.z;
  }
  
  
  /**
   * @brief ベクトルの最大値成分
   * @param [in/out] mx 比較して大きい成分
   * @param [in]     p  参照ベクトル
   */
  inline void get_max(FB::Vec3f& mx, const FB::Vec3f p) 
  {
    mx.x = (mx.x > p.x) ? mx.x : p.x;
    mx.y = (mx.y > p.y) ? mx.y : p.y;
    mx.z = (mx.z > p.z) ? mx.z : p.z;
  }
  
  /**
   * @brief 円筒形状の内外判定
   * @param [in] p テスト点座標
   * @return 内部のときに1.0を返す
   * @note 186 flop
   */
  inline float judge_cylinder(const FB::Vec3f p) 
  {
    FB::Vec3f q = rotate(angle, p-center); // 181 flop
    
    if ( (q.z < 0.0) || (q.z > depth)  ) return 0.0;
    float r = sqrtf(q.x*q.x + q.y*q.y);
    
    return ( (r<=r_fan) && (r>=r_boss) ) ? 1.0 : 0.0;
  }
  
  
  /**
   * @brief 矩形形状の内外判定
   * @param [in] p テスト点座標
   * @return 内部のときに1.0を返す
   * @note 183 flop
   */
  inline float judge_rect(const FB::Vec3f p) 
  {
    FB::Vec3f q = rotate(angle, p-center); // 181 flop
    
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
  
  
  /**
   * @brief 回転ベクトルp(alpha, beta, gamma)に対して，-pでベクトルuを回転する
   * @param [in] p 回転角度
   * @param [in] u 方向ベクトル
   * @return 角度
   * @note sin, cos = 5flop, dot=5flop, total 181flop
   */
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
  
  /** 
   * @brief 形状のbboxと投影面積を求める
   * @return 投影面積
   */
  float get_BboxArea ();
  
  
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
  void setShapeParam(const float m_nv[], const float m_ctr[], const float m_dir[], const float m_depth, const float m_width, const float m_height);
  
  /**
   * @brief 円筒の形状パラメータをセットする
   * @param [in] m_nv     法線ベクトル
   * @param [in] m_ctr    中心座標
   * @param [in] m_depth  厚み
   * @param [in] m_r_fan  外径
   * @param [in] m_r_boss 内径
   */
  void setShapeParam(const float m_nv[], const float m_ctr[], const float m_depth, const float m_r_fan, const float m_r_boss=0.0f);
  
  
  /**
   * @brief 体積率が(0,1)の間のセルに対してサブディビジョンを実施
   * @param [in]     st    開始インデクス
   * @param [in]     ed    終了インデクス
   * @param [in/out] vf    フラクション
   * @param [in/out] flop  浮動小数点演算数
   */
  void subdivision(const int st[], const int ed[], float* vf, double& flop);
  
  
  /**
   * @brief セルの8頂点の内外判定を行い，0, 1, otherに分類
   * @param [in]     st    開始インデクス
   * @param [in]     en    終了インデクス
   * @param [in/out] vf    フラクション
   * @param [in/out] flop  浮動小数点演算数
   */
  void vertex8(const int st[], const int ed[], float* vf, double& flop);
};



class ShapeMonitor : public CompoFraction {

public:
  /** デフォルトコンストラクタ */
  ShapeMonitor() {}
  
  
  /** コンストラクタ
   * @param [in] size   ローカルセル数
   * @param [in] guide  ガイドセル数
   * @param [in] crd    モニタ点座標
   * @param [in] org    ローカル領域基点座標
   * @param [in] pch    ローカル領域セル幅
   */
  ShapeMonitor(const int* size, const int guide, const float* pch, const float* org) 
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
   * @param [in/out] mid  IDの配列
   * @param [in]     id   指定ID
   */
  void setID(const int st[], const int ed[], int* mid, const int m_id);
};


#endif // _FB_FRACTION_H_
