#ifndef _SKL_IP_POLYGON_H_
#define _SKL_IP_POLYGON_H_

/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2011
 *
 */

//@file IP_Polygon.h
//@brief IP_Polygon class Header
//@author keno, FSI Team, VCAD, RIKEN

#include "Intrinsic.h"
#include "IP_Define.h"

class IP_Polygon : public Intrinsic {
protected:
  
public:
  IP_Polygon(){
  }
  ~IP_Polygon() {}

public:
  virtual void setDomain(Control* R, unsigned sz[3], SKL_REAL org[3], SKL_REAL wth[3], SKL_REAL pch[3]);
  virtual void setup(int* mid, Control* R, SKL_REAL* G_org);
  
  virtual const char* getExampleName(void) {
    return ("Polygon");
  }
};
#endif // _SKL_IP_POLYGON_H_
