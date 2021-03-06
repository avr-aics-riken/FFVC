#ifndef _IP_SHERE_H_
#define _IP_SHERE_H_

//##################################################################################
//
// FFV-C : Frontflow / violet Cartesian
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2015 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2016 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
// Copyright (c) 2016-2020 Research Institute for Information Research, Kyushu University.
// All rights reserved.
//
//##################################################################################

/*
 * @file   IP_Sphere.h
 * @brief  IP_Sphere class Header
 * @author aics
 */

#include "Intrinsic.h"
#include "IP_Define.h"
#include "common/Vec3.h" // defined in Polylib
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
   * @return 交点距離
   */
  REAL_TYPE cut_line(const Vec3r b, const int dir, const REAL_TYPE r);


  /**
   * @brief 点pの属するセルインデクスを求める
   * @param [in] p   探索座標
   * @param [in] ol  基点座標
   * @param [in] pch 格子幅
   * @return cell index
   */
  Vec3i find_index(const Vec3r p, const Vec3r ol, const Vec3r pch);


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
   * @brief 計算領域のセルIDを設定する
   * @param [in,out] bcd   　　BCindex B
   * @param [in]     R     　　Controlクラスのポインタ
   * @param [in]     NoMedium 媒質数
   * @param [in]     mat   　　MediumListクラスのポインタ
   * @param [in]     NoCompo  コンポーネント数
   * @param [in]     cmp      CompoListクラスのポインタ
   * @param [out]    cut      カット情報
   * @param [out]    bid      境界ID
   */
  virtual void setup(int* bcd,
                     Control* R,
                     const int NoMedium,
                     const MediumList* mat,
                     const int NoCompo,
                     const CompoList* cmp,
                     int* cutL,
                     int* cutU,
                     int* bid);

};
#endif // _IP_SHERE_H_
