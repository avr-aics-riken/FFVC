#ifndef _IP_RSP_H_
#define _IP_RSP_H_

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
   * @param [in,out] bcd      BCindex
   * @param [in]     R        Controlクラスのポインタ
   * @param [in]     NoMedium 媒質数
   * @param [in]     mat      MediumListクラスのポインタ
   * @param [in]     cut      交点情報
   * @param [in]     bid      境界ID
   */
  void setup(int* bcd, Control* R, const int NoMedium, const MediumList* mat);

};
#endif // _IP_RSP_H_
