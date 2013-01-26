#ifndef _IP_RECT_H_
#define _IP_RECT_H_

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
 * @file   IP_Rect.h
 * @brief  IP_Rect class Header
 * @author kero
 */

#include "Intrinsic.h"
#include "IP_Define.h"

class IP_Rect : public Intrinsic {
protected:
  unsigned even;       ///< 偶数分割のチェック
  
public:
  std::string m_fluid; ///< 流体のラベル
  std::string m_solid; ///< 固体のラベル
  
public:
  /** コンストラクタ */
  IP_Rect() {
    even = OFF;
  }
  
  /**　デストラクタ */
  ~IP_Rect() {}

public:
  
  // パラメータを取得する
  virtual bool getTP(Control* R, TPControl* tpCntl);
  
  
  // Rectの領域情報を設定する
  virtual void setDomain(Control* R, const int* sz, REAL_TYPE* org, REAL_TYPE* reg, REAL_TYPE* pch);
  
  
  // 計算領域のセルIDを設定する
  virtual void setup(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, MediumList* mat);
  
  
  // 例題の名称を返す
  virtual const char* getExampleName(void) 
  {
    return ("Rectangular");
  }
};
#endif // _IP_RECT_H_
