/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2010
 *
 */

#include "SklFactoryCBC.h"
#include "SklSolverCBC.h"

SklSolverBase* SklFactoryCBC::CreateSolvObj(int sType)
{
  return (SklSolverBase*)new SklSolverCBC(sType);
}

void SklFactoryCBC::DestorySolvObj(int sType, SklSolverBase* obj)
{
  SklSolverCBC* solv = 
           dynamic_cast<SklSolverCBC*>(obj);
  if( solv ) delete solv;
  else       delete obj;
}

