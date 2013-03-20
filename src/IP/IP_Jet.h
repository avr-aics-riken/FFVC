#ifndef _IP_JET_H_
#define _IP_JET_H_

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
 * @file   IP_Jet.h
 * @brief  IP_Jet class Header
 * @author kero
 */

#include "Intrinsic.h"
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
  int mode;                  ///< モード（2D or 3D）
  int pat_1;                 ///< ジェットのパターン (Ring1 on/off)
  int pat_2;                 ///< ジェットのパターン (Ring2 on/off)
  REAL_TYPE r1i, r1o;        ///< Ring1の内外径 [m]
  REAL_TYPE r2i, r2o;        ///< Ring2の内外径 [m]
  REAL_TYPE omg1;            ///< Ring1の角速度（符号が正のとき、x軸に向かって右ねじ）
  REAL_TYPE omg2;            ///< Ring2の角速度 [rad/s]
  REAL_TYPE q1, q2;          ///< Ring1, 2の流入流量 [m^3/s]
  REAL_TYPE n1, n2;          ///< Ring1, 2の回転数 [1/s]
  REAL_TYPE a1, a2;          ///< Ring1, 2の面積 [m^2]
  REAL_TYPE RefV;            ///< 代表速度 [m/s]
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
    mode = 0;
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
    RefV = 0.0;
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

  // パラメータを取得する
  virtual bool getTP(Control* R, TPControl* tpCntl);
  
  
  // Rectの領域情報を設定する
  virtual void setDomain(Control* R, const int* sz, REAL_TYPE* org, REAL_TYPE* reg, REAL_TYPE* pch);
  
  
  // 計算領域のセルIDを設定する
  virtual void setup(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, MediumList* mat);
  
  
  // パラメータ表示
  virtual void printPara(FILE* fp, const Control* R);
  
  
  // 例題の名称を返す
  virtual const char* getExampleName() 
  {
    return ("Jet");
  }
  
  
  
  // Jetの流入境界条件による発散値の修正
  REAL_TYPE divJetInflow(REAL_TYPE* div, const int* bv, REAL_TYPE* vf, double& flop);
  
  
  // Jetの流入境界条件をガイドセルに代入
  void vobcJetInflowGC(REAL_TYPE* v);
  
  
  // 流束型流入境界条件
  void vobc_pv_JetInflow(REAL_TYPE* wv,
                         const REAL_TYPE rei,
                         const REAL_TYPE* v0,
                         const int* bv,
                         double& flop);
  
};

#endif // _IP_JET_H_
