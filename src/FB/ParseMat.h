#ifndef _FB_PARSE_M_H_
#define _FB_PARSE_M_H_

/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file ParseMat.h
//@brief FlowBase ParseMat class Header
//@author keno, FSI Team, VCAD, RIKEN

#include "Skl.h"
#include "SklSolverBase.h"
#include "FB_Define.h"
#include "Control.h"
#include "Material.h"
#include "IDtable.h"
#include "Component.h"
#include "FBUtility.h"
#include "config/SklSolverConfig.h"

using namespace SklCfg;  // to use SklSolverConfig* cfg

class  ParseMat {
private:
  bool ChkList[property_END];  // # of parameters in MaterialList must be less than # of property_END
  
  unsigned NoCompo, NoBaseMat, NoBC, NoID, Unit_Temp;
  unsigned NoMaterial, NoFluid, NoSolid, KOS;
  unsigned imax, jmax, kmax, guide, size[3];
  REAL_TYPE BaseTemp, DiffTemp;
  
  SklSolverConfig*  CF;  // for XML parsing
  MaterialList*     mat;
  MaterialList*     BaseMat;      // for internal use of this class
  IDtable*          iTable;
  
public:
  ParseMat() {
    NoCompo      = 0;
    NoBC         = 0;
    NoID         = 0;
    NoBaseMat    = 0;
    NoMaterial   = 0;
    NoFluid      = 0;
    NoSolid      = 0;
    imax         = 0;
    jmax         = 0;
    kmax         = 0;
    Unit_Temp    = 0;
    KOS          = 0;
    BaseTemp     = 0.0;
    DiffTemp     = 0.0;
    CF           = NULL;
    mat          = NULL;
    BaseMat      = NULL;
    iTable       = NULL;
    for (int i=0; i<LABEL; i++) ChkList[i]=false;
    for (unsigned i=0; i<3; i++) size[i]=0.0;
  }
  ~ParseMat() {
    if (BaseMat) delete [] BaseMat;
  }
  
protected:
  bool chkList4Solver      (int m);
  bool copyMaterials       (unsigned& odr, unsigned id);
  bool isIDinList          (unsigned RefID);
  
  unsigned getMatIDinTable (unsigned id);
  unsigned getOdrInMatList (unsigned id, unsigned Msize, MaterialList* mlist);
  unsigned missingMessage  (int m, unsigned key);
  
  void chkList             (FILE* fp, CompoList* compo, unsigned basicEq);
  void getPvalue           (const CfgParam* p, REAL_TYPE &value);
  void printMatList        (FILE* fp, unsigned Max, MaterialList* mlist);
  void printRelation       (FILE* fp, CompoList* compo);
  void readMedium          (const CfgElem *elmL2, int m);
  
  //@fn bool isHeatProblem(void) const
  bool isHeatProblem(void) const {
    return ( ( KOS != FLOW_ONLY ) ? true : false );
  }
    
public:
  bool chkStateList          (CompoList* compo);
  bool receiveCfgPtr         (SklSolverConfig* cfg);
  
  void chkList               (FILE* mp, FILE* fp, CompoList* compo, unsigned basicEq);
  void chkState_Mat_Cmp      (CompoList* compo);
  void getXMLmaterial        (void);
  void makeLinkCmpMat        (CompoList* compo);
  void makeMaterialList      (void);
  void printMaterialList     (FILE* mp, FILE* fp);
  void setControlVars        (Control* Cref, IDtable* itbl, MaterialList* m_mat, SklSolverConfig* cfg);
  
  // ----------> debug function
  void dbg_printBaseMaterialList (FILE* mp, FILE* fp);
  void dbg_printRelation         (FILE* mp, FILE* fp, CompoList* compo);
  
};

#endif // _FB_PARSE_M_H_
