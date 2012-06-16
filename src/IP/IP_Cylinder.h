#ifndef _IP_CYL_H_
#define _IP_CYL_H_

/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file IP_Cylinder.h
//@brief IP_Cylinder class Header
//@author keno, FSI Team, VCAD, RIKEN

#include "../FB/Intrinsic.h"
#include "IP_Define.h"

class IP_Cylinder : public Intrinsic {
protected:
  unsigned mode;
  REAL_TYPE width;
  REAL_TYPE height;
  REAL_TYPE drv_length;
  std::string m_fluid;
  std::string m_solid;
  std::string m_driver;
  std::string m_driver_face;
  
public:
  IP_Cylinder(){
    mode   = 0;
    width  = 0.0;
    height = 0.0;
    drv_length = 0.0;
  }
  ~IP_Cylinder() {}

public:
  virtual bool getTP(Control* R, TPControl* tpCntl);
  
  virtual void setDomain(Control* R, unsigned sz[3], REAL_TYPE org[3], REAL_TYPE wth[3], REAL_TYPE pch[3]);
  virtual void setup(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, MediumList* mtl);
  virtual void printPara(FILE* fp, Control* R);
  
  virtual const char* getExampleName(void) {
    return ("Back Step");
  }
};
#endif // _IP_CYL_H_
