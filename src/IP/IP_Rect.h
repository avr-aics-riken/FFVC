#ifndef _IP_RECT_H_
#define _IP_RECT_H_

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
  
  virtual void setup(int* mid, Control* R, REAL_TYPE* G_org, const int NoMedium, const MediumList* mat);
  
};
#endif // _IP_RECT_H_
