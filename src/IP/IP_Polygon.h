#ifndef _IP_POLYGON_H_
#define _IP_POLYGON_H_

/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file IP_Polygon.h
//@brief IP_Polygon class Header
//@author keno, FSI Team, VCAD, RIKEN

#include "../FB/Intrinsic.h"
#include "IP_Define.h"

class IP_Polygon : public Intrinsic {
protected:
  
public:
  IP_Polygon(){
  }
  ~IP_Polygon() {}

public:
  virtual void setDomain(Control* R, unsigned sz[3], REAL_TYPE org[3], REAL_TYPE wth[3], REAL_TYPE pch[3]);
  virtual void setup(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, MediumList* mat);
  
  virtual const char* getExampleName(void) {
    return ("Polygon");
  }
};
#endif // _IP_POLYGON_H_
