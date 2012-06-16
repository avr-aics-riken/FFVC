#ifndef _IP_FREEJET_H_
#define _IP_FREEJET_H_

/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

/**
 @file IP_FreeJet.h
 @brief IP_FreeJet class Header
 @author keno, FSI Team, RIKEN
 */

#include "../FB/Intrinsic.h"
#include "IP_Define.h"

class FreeJet : public Intrinsic {
public:
  FreeJet(){}
  ~FreeJet() {}
  
protected:
  virtual bool printPara(FILE* fp, Control* R);
  
  virtual void setup(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, IDtable* itbl);
  
public:
  virtual bool getParaXML(SklSolverConfig* CF, Control* R);
  virtual bool initVars(Control* R);
  virtual bool setDomain(Control* R, unsigned* G_size, REAL_TYPE* G_org, REAL_TYPE* G_Lbx);
  
  virtual const char* getExampleName(void) {
    return ("FreeJet");
  }
  
  virtual void PostInit(REAL_TYPE &checkTime, Control* R);
};

#endif // _IP_FREEJET_H_
