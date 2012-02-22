/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file SklSolverCBCUsage.C
//@brief SklSolverCBC class
//@author keno, FSI Team, VCAD, RIKEN

#include "SklSolverCBC.h"

//@fn bool SklSolverCBC::SklSolverUsage(const char* cmd)
//@brief ソルバークラスの説明
//@param cmd ソルバークラス名
bool
SklSolverCBC::SklSolverUsage(const char* cmd) {
  if( !cmd ) return false;
  
  printf("\nSolver Name\n");
  printf("\t%s::CBC      version %5.2f with FlowBase(%5.2f)\n", cmd, (REAL_TYPE)VERS_CBC*0.01, (REAL_TYPE)FB_VERS*0.01);
  printf("Features\n");
  printf("\tThree-dimensional Cartesian flow solver for unsteady flow on Collocated Mesh\n");
  printf("\t2010-     keno@FSI, VCAD System Research Program, RIKEN\n");
  printf("\n");

  return true;
}
