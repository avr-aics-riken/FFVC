#ifndef _IP_SHC1D_H_
#define _IP_SHC1D_H_

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
 @file   IP_SHC1D.h
 @brief  IP_SHC1D class Header
 @author kero
 */

#include "Intrinsic.h"
#include "IP_Define.h"

class IP_SHC1D : public Intrinsic {
  
public:   
  std::string m_inactive;   ///< 固体で不活性セルのラベル
  
public:
  /** コンストラクタ */
  IP_SHC1D(){}
  
  /**　デストラクタ */
  ~IP_SHC1D() {}

public:

  virtual bool getTP(Control* R, TPControl* tpCntl);
  
  virtual void setDomainParameter(Control* R, const int* sz, REAL_TYPE* org, REAL_TYPE* reg, REAL_TYPE* pch);
  
  virtual void setup(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, MediumList* mat);
  
  
  
  // 計算領域のセルIDとカット情報を設定する
  void setup_bc(int* bid);
  
};

#endif // _IP_SHC1D_H_
