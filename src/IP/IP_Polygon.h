#ifndef _IP_POLYGON_H_
#define _IP_POLYGON_H_

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
  
  virtual void setup(int* mid, Control* R, REAL_TYPE* G_org, const int NoMedium, const MediumList* mat);
  
};
#endif // _IP_POLYGON_H_
