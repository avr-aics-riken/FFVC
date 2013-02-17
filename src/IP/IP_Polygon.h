#ifndef _IP_POLYGON_H_
#define _IP_POLYGON_H_

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
 @file   IP_Polygon.h
 @brief  IP_Polygon class Header
 @author kero
 */

#include "Intrinsic.h"
#include "IP_Define.h"

#ifdef __ARCH_BG
#include <string.h> // for memset()
#endif

class IP_Polygon : public Intrinsic {
protected:
  
public:
  /** コンストラクタ */
  IP_Polygon() {}
  
  /**　デストラクタ */
  ~IP_Polygon() {}

public:

  // 領域情報を設定する
  virtual void setDomain(Control* R, const int* sz, REAL_TYPE* org, REAL_TYPE* reg, REAL_TYPE* pch);
  
  
  // 計算領域のセルIDを設定する
  virtual void setup(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, MediumList* mat);
  
  
  /** 
   * @brief 例題の名称を返す
   */
  virtual const char* getExampleName()
  {
    return ("Polygon");
  }
};
#endif // _IP_POLYGON_H_
