#ifndef _IP_PMT_H_
#define _IP_PMT_H_

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
 * @file   IP_PMT.h
 * @brief  IP_PMT class Header
 * @author kero
 */

#include "Intrinsic.h"
#include "IP_Define.h"

class IP_PMT : public Intrinsic {
public:
  /** コンストラクタ */
  IP_PMT(){}
  
  /**　デストラクタ */
  ~IP_PMT() {}

public:
  std::string m_fluid; ///< 流体のラベル
  std::string m_solid; ///< 固体のラベル
  
protected:

public:
  
  virtual bool getTP(Control* R, TPControl* tpCntl);
  
  virtual void setDomainParameter(Control* R, const int* sz, REAL_TYPE* org, REAL_TYPE* reg, REAL_TYPE* pch);
  
  virtual void setup(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, MediumList* mat);
  
};
#endif // _IP_PMT_H_
