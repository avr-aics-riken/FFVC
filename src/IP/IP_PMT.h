#ifndef _SKL_IP_PMT_H_
#define _SKL_IP_PMT_H_

/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2011
 *
 */

//@file IP_PMT.h
//@brief IP_PMT class Header
//@author keno, FSI Team, VCAD, RIKEN

#include "Intrinsic.h"
#include "IP_Define.h"

class IP_PMT : public Intrinsic {
public:
  IP_PMT(){}
  ~IP_PMT() {}
  
protected:

public:
  virtual void setDomain(Control* R, unsigned sz[3], SKL_REAL org[3], SKL_REAL wth[3], SKL_REAL pch[3]);
  virtual void setup(int* mid, Control* R, SKL_REAL* G_org);
  
  virtual const char* getExampleName(void) {
    return ("Performance Test");
  }
};
#endif // _SKL_IP_PMT_H_
