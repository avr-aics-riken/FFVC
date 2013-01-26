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
 @file   IP_Jet.h
 @brief  IP_Jet class Header
 @author kero
 */

#include "Intrinsic.h"
#include "IP_Define.h"

class IP_Jet : public Intrinsic {

public:
  int mode;                  ///< モード（2D or 3D）
  int pat_0;                 ///< ジェットのパターン (コア　 on/off)
  int pat_1;                 ///< ジェットのパターン (Ring1 on/off)
  int pat_2;                 ///< ジェットのパターン (Ring2 on/off)
  REAL_TYPE r0;              ///< コア半径
  REAL_TYPE r1i, r1o;        ///< Ring1の内外径
  REAL_TYPE r2i, r2o;        ///< Ring2の内外径
  REAL_TYPE omg0;            ///< コアの角速度
  REAL_TYPE omg1;            ///< Ring1の角速度
  REAL_TYPE omg2;            ///< Ring2の角速度
  REAL_TYPE RefV;            ///< 代表速度
  std::string m_fluid;       ///< 流体のラベル
  std::string m_solid;       ///< 固体のラベル
  
  /*
   
   0       r0    r1i   r1o  r2i   r2o
   |--------+-----+-----+----+-----+-------> Radius
   | Fluid  |  S  |  F  |  S |  F  |  Solid
   
   */
  
public:
  /** コンストラクタ */
  IP_Jet(){
    mode = 0;
    r0   = 0.0;
    r1i  = 0.0;
    r1o  = 0.0;
    r2i  = 0.0;
    r2o  = 0.0;
    omg0 = 0.0;
    omg1 = 0.0;
    omg2 = 0.0;
    RefV = 0.0;
    pat_0 = -1;
    pat_1 = -1;
    pat_2 = -1;
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
  virtual bool printPara(FILE* fp, const Control* R);
  
  
  // 例題の名称を返す
  virtual const char* getExampleName(void) 
  {
    return ("Jet");
  }
  
};

#endif // _IP_JET_H_
