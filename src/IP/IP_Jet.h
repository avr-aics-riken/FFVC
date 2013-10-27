#ifndef _IP_JET_H_
#define _IP_JET_H_

//##################################################################################
//
// FFV-C : Frontflow / violet Cartesian
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2013 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2013 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

/**
 * @file   IP_Jet.h
 * @brief  IP_Jet class Header
 * @author kero
 */

#include "../FB/Intrinsic.h"
#include "IP_Define.h"

class IP_Jet : public Intrinsic {

protected:
  /*
  typedef struct {
    unsigned shape;     ///< 形状
    int direction;      ///< 方向
    REAL_TYPE diameter; ///< 直径
    REAL_TYPE length;   ///< 長さ
  } Driver_property;
  
  Driver_property driver;    ///< ドライバの特性
   */
  
  // パラメータは無次元で保持
  int pat_1;                 ///< ジェットのパターン (Ring1 on/off)
  int pat_2;                 ///< ジェットのパターン (Ring2 on/off)
  REAL_TYPE r1i, r1o;        ///< Ring1の内外径 [m]
  REAL_TYPE r2i, r2o;        ///< Ring2の内外径 [m]
  REAL_TYPE omg1;            ///< Ring1の角速度（符号が正のとき、x軸に向かって右ねじ）
  REAL_TYPE omg2;            ///< Ring2の角速度 [rad/s]
  REAL_TYPE q1, q2;          ///< Ring1, 2の流入流量 [m^3/s]
  REAL_TYPE n1, n2;          ///< Ring1, 2の回転数 [1/s]
  REAL_TYPE a1, a2;          ///< Ring1, 2の面積 [m^2]
  std::string m_fluid;       ///< 流体のラベル
  std::string m_solid;       ///< 固体のラベル
  
  /*
   
   0           r1i   r1o  r2i   r2o
   |------------+-----+----+-----+-------> Radius
   |    Solid   |  F  |  S |  F  |  Solid
   
   */
  
public:
  /** コンストラクタ */
  IP_Jet(){
    pat_1 = -1;
    pat_2 = -1;
    r1i  = 0.0;
    r1o  = 0.0;
    r2i  = 0.0;
    r2o  = 0.0;
    omg1 = 0.0;
    omg2 = 0.0;
    q1   = 0.0;
    q2   = 0.0;
    n1   = 0.0;
    n2   = 0.0;
    a1   = 0.0;
    a2   = 0.0;
/*
    driver.shape = 0;
    driver.direction = -1;
    driver.diameter = 0.0;
    driver.length = 0.0;
 */
  }
  
  /**　デストラクタ */
  ~IP_Jet() {}
  
  
  
public:

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
   * @brief 外部境界の設定
   * @param [in]     face     面番号
   * @param [in,out] bcd      BCindexx B
   * @param [in]     R        Controlクラスのポインタ
   * @param [in]     G_org    グローバルな原点（無次元）
   * @param [in]     NoMedium 媒質数
   * @param [in]     mat      MediumListクラスのポインタ
   * @param [out]    cut      カット情報
   * @param [out]    bid      境界ID
   */
  virtual void setOBC(const int face,
                      int* bcd,
                      Control* R,
                      REAL_TYPE* G_org,
                      const int NoMedium,
                      const MediumList* mat,
                      float* cut,
                      int* bid);
  
  
  /**
   * @brief Jetの流入境界条件による発散値の修正
   * @param [in,out] div   発散値
   * @param [in]     face  面番号
   * @param [in,out] vf    セルフェイス速度
   * @param [out]    sum   sum[0] 無次元流入量, sum[1] 無次元平均速度のもと
   * @param [in,out] flop  flop count
   */
  void divJetInflow(REAL_TYPE* div, const int face, REAL_TYPE* vf, REAL_TYPE* sum, double& flop);
  
  
  /**
   * @brief Jetの流入境界条件をガイドセルに代入
   * @param [in,out] v     セルセンター速度
   * @param [in]     face  面番号
   */
  void vobcJetInflowGC(REAL_TYPE* v, const int face);
  
  
  /**
   * @brief Jetの流入境界条件　Xマイナス方向のみ
   * @param [in,out] wv    疑似速度
   * @param [in]     face  面番号
   * @param [in]     rei   レイノルズ数の逆数
   * @param [in]     v0    速度ベクトル（n-step）
   * @param [in,out] flop  flop count
   */
  void vobc_pv_JetInflow(REAL_TYPE* wv,
                         const int face, 
                         const REAL_TYPE rei,
                         const REAL_TYPE* v0,
                         double& flop);
  
};

#endif // _IP_JET_H_
