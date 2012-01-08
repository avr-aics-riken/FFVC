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

extern int SklSolvTypeNum;

const char* SklSolvIdentifier[SKL_SOLV_TYPE_NUM] = {
    "Unknown"
    ,"CBC"
};

bool
CreateSolverNameList() {
  SklSolvTypeNum = SKL_SOLV_TYPE_NUM;
  if( !(SklSolverName = new char*[SklSolvTypeNum]) )
    return false;
  int i;
  for(i=0; i<SklSolvTypeNum; i++){
    if( !(SklSolverName[i] = new char[SKL_MAX_STR_LENGTH]) )
      return false;
    strncpy(SklSolverName[i],
            SklSolvIdentifier[i],
            strlen(SklSolvIdentifier[i])+1 );
  }
  return true;
}

SklSolverBase*
SklCreateSolverObject(int sType) {
  SklSolverBase* obj = NULL;
  switch (sType){
  case SKL_UNKNOWN:
    break;
  case SKL_CBC: {
      SklFactoryCBC factory;
      obj = factory.CreateSolvObj(sType);
    } break;
  default:
    break;
  }

  return obj;
}

