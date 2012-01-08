#ifndef _SKL_IP_FREEJET_H_
#define _SKL_IP_FREEJET_H_

/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2010
 *
 */

/**
 @file IP_FreeJet.h
 @brief IP_FreeJet class Header
 @author keno, FSI Team, RIKEN
 */

#include "Intrinsic.h"
#include "IP_Define.h"

class FreeJet : public Intrinsic {
public:
  FreeJet(){}
  ~FreeJet() {}
  
protected:
  virtual bool getXML(SklSolverConfig* CF, Control* R);
  virtual bool printPara(FILE* fp, Control* R);
  
  virtual void setup(int* mid, Control* R, SKL_REAL* G_org);
  
public:
  virtual bool getParaXML(SklSolverConfig* CF, Control* R);
  virtual bool initVars(Control* R);
  virtual bool setDomain(Control* R, unsigned* G_size, SKL_REAL* G_org, SKL_REAL* G_Lbx);
  
  virtual const char* getExampleName(void) {
    return ("FreeJet");
  }
  
  virtual void PostInit(SKL_REAL &checkTime, Control* R);
};

#endif // _SKL_IP_FREEJET_H_
