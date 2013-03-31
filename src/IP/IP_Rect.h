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
 
public:
  /** コンストラクタ */
  IP_Rect() {}
  
  /**　デストラクタ */
  ~IP_Rect() {}

public:
  
  virtual bool getTP(Control* R, TPControl* tpCntl);
  
  virtual void printPara(FILE* fp, const Control* R);
  
  virtual void setup(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, MediumList* mat);
  
};
#endif // _IP_RECT_H_
