#ifndef _FB_INTRNSC_H_
#define _FB_INTRNSC_H_

/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file Intrinsic.h
//@brief FlowBase Intrinsic class Header
//@author keno, FSI Team, VCAD, RIKEN

#include "Skl.h"
#include "SklSolverBase.h"
#include "FB_Define.h"
#include "BndOuter.h"
#include "Control.h"
#include "Component.h"
#include "FBUtility.h"
#include <math.h>
#include <fstream>
#include "Parallel_node.h"
#include "SklUtil.h"
#include "vec3.h"

// TextParser 
#include "TPControl.h"

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
  id_Polygon,
  id_Sphere
};

class Intrinsic : public Parallel_Node {

public:
  unsigned size[3];
  unsigned imax, jmax, kmax, guide;
  REAL_TYPE RefL;
  
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
  
  virtual bool getTP(Control* R, TPControl* tpCntl) { return true; };
  
  virtual void initCond(REAL_TYPE* v, REAL_TYPE* p) {};
  virtual void PostInit(REAL_TYPE &checkTime, Control* R) {};
  virtual void printExample(FILE* fp, const char* str);
  virtual void printParaInfo(FILE* mp, FILE* fp, Control* R);
  virtual void printPara(FILE* fp, Control* R);
  virtual void setDomain(Control* R, unsigned sz[3], REAL_TYPE org[3], REAL_TYPE wth[3], REAL_TYPE pch[3]) {};
  virtual void setup(int* mid, Control* R, REAL_TYPE* G_org) {};
  virtual void setup_cut(int* mid, Control* R, REAL_TYPE* G_org, float* cut) {};
  
  void setControlVars(Control* R);
  void writeSVX(REAL_TYPE *vf, int *id, Control* R);
  void writeSVX(int *id, Control* R);
  
};


class IP_Users : public Intrinsic {
public:
  IP_Users(){}
  ~IP_Users() {}
  
protected:
  
public:
  // @fn const char* getExampleName(void)
  // @brief ユーザー例題の名称を返す
  const char* getExampleName(void) {
    return ("User's problem");
  }
};

#endif // _FB_INTRNSC_H_
