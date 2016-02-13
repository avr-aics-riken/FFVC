#ifndef _IP_RECT_H_
#define _IP_RECT_H_

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
 * @file   IP_Rect.h
 * @brief  IP_Rect class Header
 * @author aics
 */

#include "Intrinsic.h"
#include "IP_Define.h"

class IP_Rect : public Intrinsic {
 
public:
  /** コンストラクタ */
  IP_Rect() {}
  
  /**　デストラクタ */
  ~IP_Rect() {}

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
  
};
#endif // _IP_RECT_H_
