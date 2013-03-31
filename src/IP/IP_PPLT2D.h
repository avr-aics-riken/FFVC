#ifndef _IP_PPLT2D_H_
#define _IP_PPLT2D_H_

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
 * @file   IP_PPLT2D.h
 * @brief  IP_PPLT2D class Header
 * @author kero
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
  
  virtual bool getTP(Control* R, TPControl* tpCntl);
  
  virtual void setDomainParameter(Control* R, const int* sz, REAL_TYPE* org, REAL_TYPE* reg, REAL_TYPE* pch);
  
  virtual void setup(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, MediumList* mat);
  
};
#endif // _IP_PPLT2D_H_
