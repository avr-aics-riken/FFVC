#ifndef _SKL_IP_CYL_H_
#define _SKL_IP_CYL_H_

/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2011
 *
 */

//@file IP_Cylinder.h
//@brief IP_Cylinder class Header
//@author keno, FSI Team, VCAD, RIKEN

#include "Intrinsic.h"
#include "IP_Define.h"

class IP_Cylinder : public Intrinsic {
protected:
  unsigned mode;
  SKL_REAL width;
  SKL_REAL height;
  SKL_REAL drv_length;
  
public:
  IP_Cylinder(){
    mode   = 0;
    width  = 0.0;
    height = 0.0;
    drv_length = 0.0;
  }
  ~IP_Cylinder() {}

public:
  virtual bool getXML(SklSolverConfig* CF, Control* R);
  
  virtual void setDomain(Control* R, unsigned sz[3], SKL_REAL org[3], SKL_REAL wth[3], SKL_REAL pch[3]);
  virtual void setup(int* mid, Control* R, SKL_REAL* G_org);
  virtual void printPara(FILE* fp, Control* R);
  
  virtual const char* getExampleName(void) {
    return ("Back Step");
  }
};
#endif // _SKL_IP_CYL_H_
