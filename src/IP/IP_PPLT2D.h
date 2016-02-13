#ifndef _IP_PPLT2D_H_
#define _IP_PPLT2D_H_

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
//##################################################################################

/** 
 * @file   IP_PPLT2D.h
 * @brief  IP_PPLT2D class Header
 * @author aics
 */

#include "Intrinsic.h"
#include "IP_Define.h"

class IP_PPLT2D : public Intrinsic {
  
public:
  /** コンストラクタ */
  IP_PPLT2D(){}
  
  /**　デストラクタ */
  ~IP_PPLT2D() {}
  
protected:

public:
  
  /*
   * @brief パラメータをロード
   * @param [in] R      Controlクラス
   * @param [in] tpCntl テキストパーサクラス
   * @return true-成功, false-エラー
   */
  virtual bool getTP(Control* R, TextParser* tpCntl);
  
  
  /* @brief 領域パラメータを設定する
   * @param [in]     R   Controlクラスのポインタ
   * @param [in]     sz  分割数
   * @param [in,out] org 計算領域の基点
   * @param [in,out] reg 計算領域のbounding boxサイズ
   * @param [in,out] pch セル幅
   */
  virtual void setDomainParameter(Control* R, const int* sz, REAL_TYPE* org, REAL_TYPE* reg, REAL_TYPE* pch);
  
};
#endif // _IP_PPLT2D_H_
