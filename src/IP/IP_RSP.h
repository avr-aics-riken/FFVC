#ifndef _SKL_IP_RSP_H_
#define _SKL_IP_RSP_H_

/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2010
 *
 */

//@file IP_RSP.h
//@brief IP_RSP class Header
//@author keno, FSI Team, VCAD, RIKEN

#include "Intrinsic.h"
#include "IP_Define.h"

class IP_RSP : public Intrinsic {
public:
  IP_RSP(){}
  ~IP_RSP() {}
  
protected:

public:
  void setDomain(Control* R, unsigned sz[3], SKL_REAL org[3], SKL_REAL wth[3], SKL_REAL pch[3]);
  void setup(int* mid, Control* R, SKL_REAL* G_org);
  
  const char* getExampleName(void) {
    return ("Rayleigh's Problem");
  }
};
#endif // _SKL_IP_RSP_H_
