#ifndef _IP_SHC1D_H_
#define _IP_SHC1D_H_

/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

/**
 @file IP_SHC1D.h
 @brief IP_SHC1D class Header
 @author keno, FSI Team, RIKEN
 */

#include "../FB/Intrinsic.h"
#include "IP_Define.h"

class IP_SHC1D : public Intrinsic {
  
public:
  std::string m_inactive;
  std::string m_fluid;
  std::string m_fin;
  std::string m_isothermal;
  std::string m_adiabatic;
  
public:
  IP_SHC1D(){}
  ~IP_SHC1D() {}
  
protected:

public:
  virtual bool getTP(Control* R, TPControl* tpCntl);
  
  virtual void setDomain(Control* R, unsigned sz[3], REAL_TYPE org[3], REAL_TYPE wth[3], REAL_TYPE pch[3]);
  virtual void setup(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, MediumList* mtl);
  
  virtual const char* getExampleName(void) {
    return ("Steady 1D Heat Conduction");
  }
};

#endif // _IP_SHC1D_H_
