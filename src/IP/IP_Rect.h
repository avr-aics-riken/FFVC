#ifndef _SKL_IP_RECT_H_
#define _SKL_IP_RECT_H_

/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2011
 *
 */

//@file IP_Rect.h
//@brief IP_Rect class Header
//@author keno, FSI Team, VCAD, RIKEN

#include "Intrinsic.h"
#include "IP_Define.h"

class IP_Rect : public Intrinsic {
protected:
  unsigned even;
  
public:
  IP_Rect(){
    even = OFF;
  }
  ~IP_Rect() {}

public:
  virtual bool getXML(SklSolverConfig* CF, Control* R);
  
  virtual void setDomain(Control* R, unsigned sz[3], SKL_REAL org[3], SKL_REAL wth[3], SKL_REAL pch[3]);
  virtual void setup(int* mid, Control* R, SKL_REAL* G_org);
  
  virtual const char* getExampleName(void) {
    return ("Rectangular");
  }
};
#endif // _SKL_IP_RECT_H_
