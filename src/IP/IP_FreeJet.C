/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

/**
 @file IP_FreeJet.C
 @brief IP_FreeJet class
 @author keno, FSI Team, RIKEN
 */

#include "FreeJet.h"

bool FreeJet::getXML(SklSolverConfig* CF, Control* R)
{
  const char *keyword=NULL;
  ParseSteer Tree(CF);
  
  if ( !(keyword=Tree.getParam("Example")) ) return false;
  
  if ( !strcasecmp(keyword, "FreeJet") )  R->Mode.Buoyancy = FREEJET;
  else {
    stamped_printf("\tInvalid keyword is described for Example definition\n");
    return false;
  }
  return true;
}

bool FreeJet::setDomain(Control* R, unsigned* G_size, REAL_TYPE* G_org, REAL_TYPE* G_Lbx)
{
  REAL_TYPE px, py, pz;
  px = py = pz = 0.0;
  
  // check Dimension
  if (R->NoDimension != 3 ) {
    printf("\tThe size of dimension must be 3 for Free Jet example\n");
    return false;
  }
  
  // Setting depends on Example,  INTRINSIC
  px         = 1.0 / (REAL_TYPE)imax;
  py         = 1.0 / (REAL_TYPE)jmax;
  pz         = 1.0 / (REAL_TYPE)kmax;
  R->dh      = px;      // Non-dimensional values
  R->dx[0]   = px;
  R->dx[1]   = py;
  R->dx[2]   = pz;
  R->org[0]  = -0.5;
  R->org[1]  = -0.5;
  R->org[2]  = -0.5;
  break;
  
  return true;  
}

bool FreeJet::initVars(Control* R)
{
  R->RefLength   = 1.0;  // Non-dimensional length
  R->RefVelocity = 1.0;
  
  return true;
}

bool FreeJet::getParaXML(SklSolverConfig* CF, Control* R)
{
  ParsePara Tree(CF);
  
  if ( !Tree.IsSetElem("Intrinsic_Examples") )  return false;
  if ( !Tree.getEParam("Plot_Interval", R->PlotIntvl)) return false;
  return true;
}

void FreeJet::PostInit(REAL_TYPE &checkTime, Control* R)
{
  checkTime = R->PlotIntvl;  // first check time
}

bool FreeJet::printPara(FILE* fp, Control* R)
{
  if ( !fp ) {
    stamped_printf("\tFail to write into file\n");
    return false;
  }
  fprintf(fp,"   ---> Intrinsic\n");
  fprintf(fp,"\tInterval for Plot    [-]          : %12.5e\n", R->PlotIntvl);
  fprintf(fp,"\n");
  
  fflush(fp);
  
  return true;
}

void FreeJet::setup(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, MediumList* mtl)
{
  unsigned i0, j0, j1, jd;
  unsigned m;
  REAL_TYPE  fj0, fj1;
  int i,j,k;
  
  // Inner
  for (k=0; k<=(int)(kmax+1); k++) {
    for (j=0; j<=(int)(jmax+1); j++) {
      for (i=0; i<=(int)(imax+1); i++) {
        m = FBUtility::getFindexS3D(size, guide, i, j, k);
        mid[m] = 1;
      }
    }
  }
  
  fj0 = -0.5*Domain_p1/RefLength - org[1];  //  -D/2
  fj1 =  0.5*Domain_p1/RefLength - org[1];  //  +D/2
  j0 = (int)( floor( fj0/dx[1]+0.5 ) );
  j1 = (int)( floor( fj1/dx[1]+0.5 ) );
  printf("\tJet index\n");
  printf("\tj0=%d j1=%d [Fortran index]\n", j0, j1);
  
  // X_MINUS
  for (k=0; k<(int)(kmax+2*guide); k++) {
    for (j=0; j<(int)(jmax+2*guide); j++) {
      for (i=0; i<=(int)(guide-1); i++) {
        m = FBUtility::getFindexS3D(size, guide, i, j, k);
        mid[m] = 600;
      }
    }
  }
  
  // X_MINUS - Nozzle
  for (k=0; k<(int)(kmax+2*guide); k++) {
    for (j=(int)(j0+guide-1); j<=(int)(j1+guide-1); j++) {
      for (i=0; i<=(int)(guide-1); i++) {
        m = FBUtility::getFindexS3D(size, guide, i, j, k);
        mid[m] = 0;
      }
    }
  }

}
