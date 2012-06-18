#ifndef _IP_PMT_H_
#define _IP_PMT_H_

/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file IP_PMT.h
//@brief IP_PMT class Header
//@author keno, FSI Team, VCAD, RIKEN

#include "../FB/Intrinsic.h"
#include "IP_Define.h"

class IP_PMT : public Intrinsic {
public:
  IP_PMT(){}
  ~IP_PMT() {}

public:
  std::string m_fluid;
  std::string m_solid;
  
protected:

public:
  virtual bool getTP(Control* R, TPControl* tpCntl);
  
  virtual void setDomain(Control* R, unsigned sz[3], REAL_TYPE org[3], REAL_TYPE wth[3], REAL_TYPE pch[3]);
  virtual void setup(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, MediumList* mat);
  
  virtual const char* getExampleName(void) {
    return ("Performance Test");
  }
};
#endif // _IP_PMT_H_
