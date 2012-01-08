/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2010
 *
 */

#include <stdlib.h>
#include "SklGloval.h"
#include "SklSolverBase.h"
#include "SklSolverType.h"
#include "SklFactoryCBC.h"

void
SklDeleteSolverObject(int sType, SklSolverBase* obj) {
  if( !obj ) return;
  switch (sType) {
  case 0:
    break;
  case SKL_CBC: {
      SklFactoryCBC factory;
      factory.DestorySolvObj(sType, obj);
    } break;
  default:
    break;
  }
}

