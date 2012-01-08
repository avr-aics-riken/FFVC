#ifndef _SKL_IP_DUCT_H_
#define _SKL_IP_DUCT_H_

/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2011
 *
 */

//@file IP_Duct.h
//@brief IP_Duct class Header
//@author keno, FSI Team, VCAD, RIKEN

#include "Intrinsic.h"
#include "IP_Define.h"

class IP_Duct : public Intrinsic {
protected:
  typedef struct {
    unsigned shape;
    int direction;
    SKL_REAL diameter;
    SKL_REAL length;
  } Driver_property;
  
  Driver_property driver;
  
  enum shape_type {
    id_circular = 1,
    id_rectangular
  };
  
public:
  IP_Duct(){
    driver.shape = 0;
    driver.direction = -1;
    driver.diameter = 0.0;
    driver.length = 0.0;
  }
  ~IP_Duct() {}

public:
  virtual bool getXML(SklSolverConfig* CF, Control* R);
  
  virtual void setDomain(Control* R, unsigned sz[3], SKL_REAL org[3], SKL_REAL wth[3], SKL_REAL pch[3]);
  virtual void setup(int* mid, Control* R, SKL_REAL* G_org);
  virtual void printPara(FILE* fp, Control* R);
  
  virtual const char* getExampleName(void) {
    return ("Duct");
  }
};
#endif // _SKL_IP_DUCT_H_
