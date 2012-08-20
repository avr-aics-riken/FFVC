#ifndef _IP_SHERE_H_
#define _IP_SHERE_H_

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
 * @file IP_Sphere.h
 * @brief IP_Sphere class Header
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
  /** パラメータをロード
   * @param [in] R      Controlクラス
   * @param [in] tpCntl テキストパーサクラス
   * @return true-成功, false-エラー
   */
  virtual bool getTP(Control* R, TPControl* tpCntl);
  
  
  /** 領域を設定する
   * @param [in]     R     Controlクラスのポインタ
   * @param [in]     sz    分割数
   * @param [in,out] m_org 計算領域の基点
   * @param [in,out] m_reg 計算領域のbounding boxサイズ
   * @param [in,out] m_pch セル幅
   */
  virtual void setDomain(Control* R, const int* sz, REAL_TYPE* m_org, REAL_TYPE* m_reg, REAL_TYPE* m_pch);
  
  
  /**
   @brief パラメータの表示
   @param[in] fp ファイルポインタ
   @param[in] R  コントロールクラスのポインタ
   */
  virtual void printPara(FILE* fp, const Control* R);
  
  
  /** 
   @brief 例題の名称を返す
   */
  virtual const char* getExampleName(void) 
  {
    return ("Sphere");
  }
  
  
  /** 矩形の計算領域のセルIDを設定する
   * @param[in,out] mid   媒質情報の配列
   * @param[in]     R     Controlクラスのポインタ
   * @param[in]     G_org グローバルな原点（無次元）
   * @param[in]     Nmax  Controlクラスのポインタ
   * @param[in]     mat   MediumListクラスのポインタ
   */
  virtual void setup(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, MediumList* mat);
  
  
  /**
   @brief 計算領域のセルIDとカット情報を設定する
   @param[in,out] mid    IDの配列
   @param[in]     R      Controlクラスのポインタ
   @param[in]     G_org  グローバルな原点（無次元）
   @param[in]     Nmax   
   @param[in]     mat    MediumListクラスのポインタ
   @param[out]    cut    カット情報
   */
  virtual void setup_cut(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, MediumList* mat, float* cut);
  
  
  /** 点pの属するセルインデクスを求める
   * @param[in] p   探索座標
   * @param[in] ol  基点座標
   * @return 
   */
  FB::Vec3i find_index(const FB::Vec3f p, const FB::Vec3f ol);
  
  
  /**
   @brief 交点の無次元距離を計算する
   @param[in] p 基点座標
   @param[in] dir テスト方向
   @param[in] r radius
   @param[in] dh 格子幅
   @return 交点距離
   */
  float cut_line(const FB::Vec3f b, const int dir, const float r, const float dh);
};
#endif // _IP_SHERE_H_
