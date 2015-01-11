#ifndef _IP_CYL_H_
#define _IP_CYL_H_

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
// Copyright (c) 2012-2015 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

/** 
 * @file   IP_Cylinder.h
 * @brief  IP_Cylinder class Header
 * @author aics
 */

#include "Intrinsic.h"
#include "IP_Define.h"
#include "Vec3.h"
#include "FBUtility.h"

class CYL {
public:
  int ens;                 ///< cylinder1の有無
  int shape;               ///< 1-circular, 2-rect
  REAL_TYPE length_z;      ///< Z方向の長さ
  REAL_TYPE length_y;      ///< Y方向の長さ
  REAL_TYPE length_x;      ///< Y方向の長さ
  REAL_TYPE radius;        ///< 半径
  REAL_TYPE x;             ///< cylinderの中心位置座標
  REAL_TYPE y;             ///<
  std::string solid;       ///< 固体のラベル
  
  CYL(){
    ens      = 0;
    shape    = 0;
    length_x = 0.0;
    length_y = 0.0;
    length_z = 0.0;
    radius   = 0.0;
    x        = 0.0;
    y        = 0.0;
  }
  ~CYL() {}
};


class IP_Cylinder : public Intrinsic {
protected:
  
  REAL_TYPE drv_length;      ///< ドライバの長さ
  int drv_mode;              ///< ドライバのON/OFF
  std::string m_driver;      ///< ドライバ部分のラベル
  std::string m_driver_face; ///< ドライバ指定面のラベル
  int num_cyls;              ///< Cylinderの個数
  CYL cyl1;
  CYL cyl2;

  
public:
  /** コンストラクタ */
  IP_Cylinder(){
    num_cyls = 0;
    drv_mode = 0;
    drv_length = 0.0;
  }
  
  /**　デストラクタ */
  ~IP_Cylinder() {}

protected:
  
  
  /**
   * @brief 交点の無次元距離を計算する
   * @param [in] p   基点座標
   * @param [in] dir テスト方向
   * @param [in] r   radius
   * @param [in] px  格子幅
   * @return 交点距離
   */
  REAL_TYPE cut_line_2d(const Vec3r p, const int dir, const REAL_TYPE r, const REAL_TYPE px);
  
  
  /**
   * @brief 点pの属するセルインデクスを求める
   * @param [in] p   探索座標
   * @param [in] ol  基点座標
   * @param [in] pch 格子幅
   * @return cell index
   */
  Vec3i find_index(const Vec3r p, const Vec3r ol, const Vec3r pch);
  
  
  /**
   * @brief 円柱を設定
   * @param [in]     R     　　 Controlクラスのポインタ
   * @param [in]     pos_x     中心座標(x,y)
   * @param [in]     pos_y     中心座標(x,y)
   * @param [in]     radius    半径
   * @param [in]     len_z     円柱のZ方向の長さ
   * @param [in]     mid_solid 固体媒質ID
   * @param [out]    cut       カット情報
   * @param [out]    bid       境界ID
   */
  void setCircle(Control* R,
                 const REAL_TYPE pos_x,
                 const REAL_TYPE pos_y,
                 const REAL_TYPE radius,
                 const REAL_TYPE len_z,
                 const int mid_solid,
                 long long* cut,
                 int* bid);
  
  
  /**
   * @brief Rectangularを設定
   * @param [in]     R     　　 Controlクラスのポインタ
   * @param [in]     pos_x     角柱の中心座標(x,y)
   * @param [in]     pos_y     角柱の中心座標(x,y)
   * @param [in]     len_x     角柱のX方向の長さ
   * @param [in]     len_y     角柱のY方向の長さ
   * @param [in]     len_z     角柱のZ方向の長さ
   * @param [in]     mid_solid 固体媒質ID
   * @param [out]    cut       カット情報
   * @param [out]    bid       境界ID
   */
  void setRect(Control* R,
               const REAL_TYPE pos_x,
               const REAL_TYPE pos_y,
               const REAL_TYPE len_x,
               const REAL_TYPE len_y,
               const REAL_TYPE len_z,
               const int mid_solid,
               long long* cut,
               int* bid);
  
  
public:
  
  /**
   * @brief Cylinderの個数を返す
   */
  int get_num_cyls() const
  {
    return num_cyls;
  }
  
  
  /**
   * @brief パラメータをロード
   * @apram [in] R      Controlクラス
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
   * @brief Cylinderの計算領域のセルIDを設定する
   * @param [in,out] bcd   　　BCindex B
   * @param [in]     R     　　Controlクラスのポインタ
   * @param [in]     NoMedium 媒質数
   * @param [in]     mat   　　MediumListクラスのポインタ
   * @param [out]    cut      カット情報
   * @param [out]    bid      境界ID
   */
  virtual void setup(int* bcd, Control* R, const int NoMedium, const MediumList* mat, long long* cut, int* bid);
  
};
#endif // _IP_CYL_H_
