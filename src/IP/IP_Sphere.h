#ifndef _IP_SHERE_H_
#define _IP_SHERE_H_

// #################################################################
//
// FFV : Frontflow / violet Cartesian
//
// Copyright (c) 2012-2013  All right reserved.
//
// Institute of Industrial Science, University of Tokyo, Japan. 
//
// #################################################################

/** 
 * @file   IP_Sphere.h
 * @brief  IP_Sphere class Header
 * @author kero
 */

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
  
  std::string m_fluid;       ///< 流体のラベル
  std::string m_solid;       ///< 固体のラベル
  std::string m_driver;      ///< ドライバ部分のラベル
  std::string m_driver_face; ///< ドライバ指定面のラベル
  
public:
  /** コンストラクタ */
  IP_Sphere() {
    radius = 0.0;
    drv_length = 0.0;
    drv_mode = OFF;
  }
  
  /**　デストラクタ */
  ~IP_Sphere() {}

public:
  // パラメータを取得する
  virtual bool getTP(Control* R, TPControl* tpCntl);
  
  
  // 領域パラメータを設定する
  virtual void setDomainParameter(Control* R, const int* sz, REAL_TYPE* m_org, REAL_TYPE* m_reg, REAL_TYPE* m_pch);
  
  
  // パラメータの表示
  virtual void printPara(FILE* fp, const Control* R);
  
  
  /** 
   * @brief 例題の名称を返す
   */
  virtual const char* getExampleName() 
  {
    return ("Sphere");
  }
  
  
  // 計算領域のセルIDを設定する（バイナリー）
  virtual void setup(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, MediumList* mat);
  
  
  
  // 計算領域のセルIDとカット情報を設定する
  void setup_cut(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, MediumList* mat, float* cut);
  
  
  // 点pの属するセルインデクスを求める
  FB::Vec3i find_index(const FB::Vec3f p, const FB::Vec3f ol);
  
  
  // 交点計算
  float cut_line(const FB::Vec3f b, const int dir, const float r, const float dh);
};
#endif // _IP_SHERE_H_
