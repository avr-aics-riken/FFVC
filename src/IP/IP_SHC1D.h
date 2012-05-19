#ifndef _IP_SHC1D_H_
#define _IP_SHC1D_H_

/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

/**
 @file IP_SHC1D.h
 @brief IP_SHC1D class Header
 @author keno, FSI Team, RIKEN
 */

#include "Intrinsic.h"
#include "IP_Define.h"

class IP_SHC1D : public Intrinsic {
public:
  IP_SHC1D(){}
  ~IP_SHC1D() {}
  
protected:

public:
  virtual void setDomain(Control* R, unsigned sz[3], REAL_TYPE org[3], REAL_TYPE wth[3], REAL_TYPE pch[3]);
  virtual void setup(int* mid, Control* R, REAL_TYPE* G_org);
  
  virtual const char* getExampleName(void) {
    return ("Steady 1D Heat Conduction");
  }
};

#endif // _IP_SHC1D_H_
