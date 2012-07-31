#ifndef _IP_POLYGON_H_
#define _IP_POLYGON_H_

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
 @file IP_Polygon.h
 @brief IP_Polygon class Header
 @author kero
 */

#include "Intrinsic.h"
#include "IP_Define.h"

class IP_Polygon : public Intrinsic {
protected:
  
public:
  /** コンストラクタ */
  IP_Polygon() {}
  
  /**　デストラクタ */
  ~IP_Polygon() {}

public:
  /** 領域を設定する
   * @param [in]     R   Controlクラスのポインタ
   * @param [in]     sz  分割数
   * @param [in,out] org 計算領域の基点
   * @param [in,out] reg 計算領域のbounding boxサイズ
   * @param [in,out] pch セル幅
   */
  virtual void setDomain(Control* R, const int* sz, REAL_TYPE* org, REAL_TYPE* reg, REAL_TYPE* pch);
  
  
  /** 矩形の計算領域のセルIDを設定する
   * @param [in,out] mid   媒質情報の配列
   * @param [in]     R     Controlクラスのポインタ
   * @param [in]     G_org グローバルな原点（無次元）
   * @param [in]     Nmax  Controlクラスのポインタ
   * @param [in]     mat   MediumListクラスのポインタ
   */
  virtual void setup(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, MediumList* mat);
  
  
  /** 
   @brief 例題の名称を返す
   */
  virtual const char* getExampleName(void)
  {
    return ("Polygon");
  }
};
#endif // _IP_POLYGON_H_
