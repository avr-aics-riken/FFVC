#ifndef _IP_SHERE_H_
#define _IP_SHERE_H_

//##################################################################################
//
// FFV-C : Frontflow / violet Cartesian
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
 * @file   IP_Sphere.h
 * @brief  IP_Sphere class Header
 * @author aics
 */

#include "../FB/Intrinsic.h"
#include "IP_Define.h"
#include "Vec3.h"
#include "FBUtility.h"


using namespace Vec3class;

class IP_Sphere : public Intrinsic {
protected:
  REAL_TYPE radius;     ///< 球の半径
  REAL_TYPE drv_length; ///< ドライバの長さ
  int drv_mode;         ///< ドライバのON/OFF
  
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

  
  
protected:
  
  /**
   * @brief 交点の無次元距離を計算する
   * @param [in] p 基点座標
   * @param [in] dir テスト方向
   * @param [in] r radius
   * @param [in] dh 格子幅
   * @return 交点距離
   */
  float cut_line(const Vec3f b, const int dir, const float r, const float dh);
  
  
  /**
   * @brief 点pの属するセルインデクスを求める
   * @param [in] p   探索座標
   * @param [in] ol  基点座標
   * @param [in] pch 格子幅
   * @return cell index
   */
  Vec3i find_index(const Vec3f p, const Vec3f ol, const Vec3f pch);
  
  
  /**
   * @brief パラメータをロード
   * @param [in] R      Controlクラス
   * @param [in] tpCntl テキストパーサクラス
   * @return true-成功, false-エラー
   */
  virtual bool getTP(Control* R, TextParser* tpCntl);
  
  
  /**
   * @brief パラメータの表示
   * @param [in] fp ファイルポインタ
   * @param [in] R  コントロールクラスのポインタ
   */
  virtual void printPara(FILE* fp, const Control* R);
  
  
  /**
   * @brief 計算領域のセルIDとカット情報を設定する
   * @param [in,out] bcd      BCindex B
   * @param [in]     R        Controlクラスのポインタ
   * @param [in]     NoMedium 媒質数
   * @param [in]     mat      MediumListクラスのポインタ
   * @param [in]     cut      交点情報
   * @param [in]     bid      境界ID
   */
  virtual void setup(int* bcd, Control* R, const int NoMedium, const MediumList* mat, float* cut, int* bid);
  
};
#endif // _IP_SHERE_H_
