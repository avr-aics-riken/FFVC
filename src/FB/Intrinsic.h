#ifndef _SKL_FB_INTRNSC_H_
#define _SKL_FB_INTRNSC_H_

/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2010
 *
 */

//@file Intrinsic.h
//@brief FlowBase Intrinsic class Header
//@author keno, FSI Team, VCAD, RIKEN

#include "Skl.h"
#include "SklSolverBase.h"
#include "FBDefine.h"
#include "BndOuter.h"
#include "mydebug.h"
#include "Control.h"
#include "Component.h"
#include "FBUtility.h"
#include <math.h>
#include <fstream>
#include "Parallel_node.h"
#include "SklUtil.h"
extern SklParaComponent* ParaCmpo;

enum Intrinsic_class {
  id_Users = 0,
  id_Duct,
  id_PPLT2D,
  id_SHC1D,
  id_PMT,
  id_Rect,
  id_Cylinder,
  id_Step,
  id_Polygon
};

class Intrinsic : public Parallel_Node {

public:
  unsigned size[3];
  unsigned imax, jmax, kmax, guide;
  SKL_REAL RefL;
  
  enum dim_mode {
    dim_2d = 1,
    dim_3d
  };
  
  Intrinsic() {
    for (unsigned i=0; i<3; i++) size[i]=0.0;
    imax = jmax = kmax = guide = 0;
    RefL = 0.0;
  }
  virtual ~Intrinsic() {}
    
public:
  virtual bool getXML(SklSolverConfig* CF, Control* R) { return true; };
  
  virtual const char* getExampleName(void) { return NULL; };
  
  virtual void genVFfromBcx(SKL_REAL* VF, unsigned* bx);
  virtual void initCond(SKL_REAL* v, SKL_REAL* p) {};
  virtual void PostInit(SKL_REAL &checkTime, Control* R) {};
  virtual void printExample(FILE* fp, const char* str);
  virtual void printParaInfo(FILE* mp, FILE* fp, Control* R);
  virtual void printPara(FILE* fp, Control* R);
  virtual void setDomain(Control* R, unsigned sz[3], SKL_REAL org[3], SKL_REAL wth[3], SKL_REAL pch[3]) {};
  virtual void setup(int* mid, Control* R, SKL_REAL* G_org) {};
  virtual void writeSVX(SKL_REAL *vf, int *id, Control* R);
  virtual void writeSVX(int *id, Control* R);
  
  void setControlVars(Control* R);
};
#endif // _SKL_FB_INTRNSC_H_