#ifndef _IP_RSP_H_
#define _IP_RSP_H_

/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file IP_RSP.h
//@brief IP_RSP class Header
//@author keno, FSI Team, VCAD, RIKEN

#include "../FB/Intrinsic.h"
#include "IP_Define.h"

class IP_RSP : public Intrinsic {
public:
  IP_RSP(){}
  ~IP_RSP() {}
  
protected:

public:
  void setDomain(Control* R, unsigned sz[3], REAL_TYPE org[3], REAL_TYPE wth[3], REAL_TYPE pch[3]);
  void setup(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, MediumList* mtl);
  
  const char* getExampleName(void) {
    return ("Rayleigh's Problem");
  }
};
#endif // _IP_RSP_H_
