#ifndef _IP_RSP_H_
#define _IP_RSP_H_

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
  
  // 領域パラメータを設定する
  virtual void setDomainParameter(Control* R, const int* sz, REAL_TYPE* org, REAL_TYPE* reg, REAL_TYPE* pch);
  
  
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
