/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2010
 *
 */

#ifndef _SKL_SOLVER_CBC_FACTORY_H_
#define _SKL_SOLVER_CBC_FACTORY_H_

#include "SklSolvFactoryBase.h"

class SklFactoryCBC : public SklSolvFactoryBase {
public:
  SklFactoryCBC(){;};
  virtual ~SklFactoryCBC(){;};
  virtual SklSolverBase* CreateSolvObj(int sType);
  virtual void DestorySolvObj(int sType, SklSolverBase* obj);
};

#endif // _SKL_SOLVER_CBC_FACTORY_H_

