#ifndef _IP_RSP_H_
#define _IP_RSP_H_

// #################################################################
//
// CAERU Library
//
// Copyright (c) 2012-2013  All right reserved.
//
// Institute of Industrial Science, University of Tokyo, Japan. 
//
// #################################################################

//@file IP_RSP.h
//@brief IP_RSP class Header
//@author keno, FSI Team, VCAD, RIKEN

#include "Intrinsic.h"
#include "IP_Define.h"

class IP_RSP : public Intrinsic {
public:
  /** コンストラクタ */
  IP_RSP(){}
  
  /**　デストラクタ */
  ~IP_RSP() {}
  
protected:

public:
  /** 領域を設定する
   * @param [in]     R   Controlクラスのポインタ
   * @param [in]     sz  分割数
   * @param [in,out] org 計算領域の基点
   * @param [in,out] reg 計算領域のbounding boxサイズ
   * @param [in,out] pch セル幅
   */
  void setDomain(Control* R, const int* sz, REAL_TYPE* org, REAL_TYPE* reg, REAL_TYPE* pch);
  
  
  /** 矩形の計算領域のセルIDを設定する
   * @param [in,out] mid   媒質情報の配列
   * @param [in]     R     Controlクラスのポインタ
   * @param [in]     G_org グローバルな原点（無次元）
   * @param [in]     Nmax  Controlクラスのポインタ
   * @param [in]     mat   MediumListクラスのポインタ
   */
  void setup(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, MediumList* mat);
  
  
  /** 
   @brief 例題の名称を返す
   */
  const char* getExampleName(void) 
  {
    return ("Rayleigh's Problem");
  }
};
#endif // _IP_RSP_H_
