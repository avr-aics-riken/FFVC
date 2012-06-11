#ifndef _IP_DUCT_H_
#define _IP_DUCT_H_

/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
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
    REAL_TYPE diameter;
    REAL_TYPE length;
  } Driver_property;
  
  Driver_property driver;
  
  std::string m_fluid;
  std::string m_solid;
  std::string m_driver;
  std::string m_driver_face;
  
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
  virtual bool getTP(Control* R, TPControl* tpCntl);
  
  virtual void setDomain(Control* R, unsigned sz[3], REAL_TYPE org[3], REAL_TYPE wth[3], REAL_TYPE pch[3]);
  virtual void setup(int* mid, Control* R, REAL_TYPE* G_org);
  virtual void printPara(FILE* fp, Control* R);
  
  virtual const char* getExampleName(void) {
    return ("Duct");
  }
};
#endif // _IP_DUCT_H_
