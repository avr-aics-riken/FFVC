#ifndef _SKL_FB_PARA_BC_H_
#define _SKL_FB_PARA_BC_H_

/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2010
 *
 */

//@file ParseBC.h
//@brief FlowBase ParseBC class Header
//@author keno, FSI Team, VCAD, RIKEN

#include "Skl.h"
#include "SklSolverBase.h"
#include "FBDefine.h"
#include "FBUtility.h"
#include "BndOuter.h"
#include "mydebug.h"
#include "string.h"
#include "Control.h"
#include "Component.h"
#include "IDtable.h"
#include "SklUtil.h"
#include "config/SklSolverConfig.h"
#include "Parallel_node.h"
#include "parallel/SklParaComponent.h"
#include "Material.h"
using namespace SklCfg;  // to use SklSolverConfig* cfg

class ParseBC : public Parallel_Node {
private:
  SKL_REAL RefVelocity, BaseTemp, DiffTemp, RefDensity, RefSpecificHeat;
  SKL_REAL RefLength, BasePrs;
  SKL_REAL rho, nyu, cp, lambda, beta; // 無次元化の参照値
  
  unsigned ix, jx, kx, guide;
  unsigned KindOfSolver;
  unsigned NoCompo, NoBC, NoID, NoBaseBC;
  unsigned Unit_Param, monitor, Unit_Temp, Unit_Prs, Mode_Gradp;
  
  SklSolverConfig*  CF;  // for XML parsing
  CompoList*        compo;
  BoundaryOuter*    bc;
  BoundaryOuter*    BaseBc;
	
	bool HeatProblem;
  char OBCname[NOFACE][LABEL];
	
public:
	IDtable*          iTable;

public:
  ParseBC(){
    ix = jx = kx = guide = 0;
    KindOfSolver = 0;
    BaseTemp = DiffTemp = BasePrs = 0.0;
    RefVelocity = RefDensity = RefSpecificHeat = RefLength = 0.0;
    NoCompo = NoBC = NoID = NoBaseBC = monitor = 0;
    Unit_Param = 0;
    Unit_Temp = 0;
    Unit_Prs = 0;
    Mode_Gradp = 0;
    HeatProblem = false;
    CF = NULL;
    bc = NULL;
    compo = NULL;
    BaseBc = NULL;
    iTable = NULL;
    for (int i=0; i<NOFACE; i++) {
      strcpy(&OBCname[i][0], "\0");
    }
  }
  ~ParseBC() {
    if (iTable) delete [] iTable;
    if (BaseBc) delete [] BaseBc;
  }
  
protected:
  bool chkID              (void);
  bool isAlterOBC         (void);
  bool isComponent        (unsigned label);
  bool isCompoTransfer    (unsigned label);
  bool isIDinTable        (int candidate);
  
  int getNBC              (void);
  int getStateinTable     (unsigned id);
  
  unsigned getNoMaterial      (void);
  unsigned getXML_Vel_profile (const CfgElem *elmL, const char* err_str);
  unsigned oppositDir         (unsigned dir);
  unsigned scanXMLmodel       (void);
  
  void chkKeywordIBC        (const char *keyword, unsigned m);
  void chkKeywordOBC        (const char *keyword, unsigned m);
  void getFan               (const CfgElem *elmL, unsigned n);
  void getDarcy             (const CfgElem *elmL, unsigned n);
  void get_Vec              (const CfgElem *elmL, unsigned n, const char* str, SKL_REAL* v);
  void getUnitVec           (SKL_REAL* v);
  void getXML_Cell_Monitor  (const CfgElem *elmL, unsigned n, Control* C);
  void getXML_IBC_Adiabatic (const CfgElem *elmL, unsigned n);
  void getXML_IBC_CnstTemp  (const CfgElem *elmL, unsigned n);
  void getXML_IBC_IBM_DF    (const CfgElem *elmL, unsigned n);
  void getXML_IBC_HeatFlux  (const CfgElem *elmL, unsigned n);
  void getXML_IBC_HeatSrc   (const CfgElem *elmL, unsigned n);
  void getXML_IBC_HT_B      (const CfgElem *elmL, unsigned n);
  void getXML_IBC_HT_N      (const CfgElem *elmL, unsigned n);
  void getXML_IBC_HT_S      (const CfgElem *elmL, unsigned n);
  void getXML_IBC_HT_SF     (const CfgElem *elmL, unsigned n);
  void getXML_IBC_HT_SN     (const CfgElem *elmL, unsigned n);
  void getXML_IBC_IsoTherm  (const CfgElem *elmL, unsigned n);
  void getXML_IBC_Outflow   (const CfgElem *elmL, unsigned n);
  void getXML_IBC_Periodic  (const CfgElem *elmL, unsigned n);
  void getXML_IBC_PrsLoss   (const CfgElem *elmL, unsigned n);
  void getXML_IBC_Radiant   (const CfgElem *elmL, unsigned n);
  void getXML_IBC_SpecVel   (const CfgElem *elmL, unsigned n);
  void getXML_OBC_HT        (const CfgElem *elmL, unsigned n, const char* kind);
  void getXML_OBC_InOut     (const CfgElem *elmL, unsigned n);
  void getXML_OBC_Outflow   (const CfgElem *elmL, unsigned n);
  void getXML_OBC_Periodic  (const CfgElem *elmL, unsigned n);
  void getXML_OBC_SpecVH    (const CfgElem *elmL, unsigned n);
  void getXML_OBC_Trcfree   (const CfgElem *elmL, unsigned n);
  void getXML_OBC_Wall      (const CfgElem *elmL, unsigned n);
  void getXML_Vel_Params    (const CfgElem *elmL, unsigned type, SKL_REAL* ca, SKL_REAL vel, const char* err_str);
  void printBaseOBC         (FILE* fp);
  void printCompo           (FILE* fp, SKL_REAL* nv, int* ci, MaterialList* mat);
  void printFaceOBC         (FILE* fp, SKL_REAL* G_Lbx);
  void printOBC             (FILE* fp, BoundaryOuter* ref, SKL_REAL* G_Lbx, unsigned face);
  void set_Deface           (const CfgElem *elmL, unsigned n, const char* str);
  
  const CfgElem* selectFace (int face, const CfgElem* elmTop);
  
  string getOBCstr          (unsigned id);
  
  //@fn int getGlCompoBV_st_x(unsigned odr, int* gci)
  //@brief コンポーネントのBV情報st_xを返す
  int getGlCompoBV_st_x(unsigned odr, int* gci) {
    return ( gci[6*odr+0] );
  }
  
  //@fn int getGlCompoBV_st_y(unsigned odr, int* gci)
  //@brief コンポーネントのBV情報st_yを返す
  int getGlCompoBV_st_y(unsigned odr, int* gci) {
    return ( gci[6*odr+1] );
  }
  
  //@fn int getGlCompoBV_st_z(unsigned odr, int* gci)
  //@brief コンポーネントのBV情報st_zを返す
  int getGlCompoBV_st_z(unsigned odr, int* gci) {
    return ( gci[6*odr+2] );
  }
  
  //@fn int getGlCompoBV_ed_x(unsigned odr, int* gci)
  //@brief コンポーネントのBV情報st_xを返す
  int getGlCompoBV_ed_x(unsigned odr, int* gci) {
    return ( gci[6*odr+3] );
  }
  
  //@fn int getGlCompoBV_ed_y(unsigned odr, int* gci)
  //@brief コンポーネントのBV情報ed_yを返す
  int getGlCompoBV_ed_y(unsigned odr, int* gci) {
    return ( gci[6*odr+4] );
  }
  
  //@fn int getGlCompoBV_ed_z(unsigned odr, int* gci)
  //@brief コンポーネントのBV情報ed_zを返す
  int getGlCompoBV_ed_z(unsigned odr, int* gci) {
    return ( gci[6*odr+5] );
  }
  
public:
  bool isIDinCompo        (unsigned candidate_id, unsigned now);
  bool isIDinCompo        (unsigned candidate_id, int def, unsigned now);
  
  int get_BCval_int       (const CfgElem *elmL, const char* key);
  
  SKL_REAL get_BCval_real (const CfgElem *elmL, const char* key);
  
  unsigned count_Outer_Cell_ID (unsigned* medium);
  
  void chkBCconsistency   (unsigned kos);
  void getXML_Model       (void);
  void getXML_Phase       (void);
  void loadOuterBC        (void);
  void printCompoInfo     (FILE* mp, FILE* fp, SKL_REAL* nv, int* ci, MaterialList* mat);
  void printOBCinfo       (FILE* mp, FILE* fp, SKL_REAL* G_Lbx);
  void printTable         (FILE* fp);
  void receiveCompoPtr    (CompoList* CMP);
  void receiveCfgPtr      (SklSolverConfig* cfg);
  void setMedium          (Control* Cref);
  void setCompoList       (Control* C);
  void setControlVars     (Control* Cref);
  void getXML_Medium_InitTemp(void);
  void setObcPtr          (BoundaryOuter* ptr);
  void setRefMedium       (MaterialList* mat, Control* Cref);
  void setRefValue        (MaterialList* mat, CompoList* cmp, Control* C);
  
  IDtable* get_IDtable_Ptr(void);
  
  //@fn BoundaryOuter* get_BaseBC_Ptr(void)
  //@brief 基本境界リストのポインタを返す
  BoundaryOuter* get_BaseBC_Ptr(void) { 
    return BaseBc;
  }
  
  //@fn unsigned get_NoBaseBC(void)
  //@brief 基本境界条件数を返す
  unsigned get_NoBaseBC(void) {
    return NoBaseBC;
  }
};

#endif // _SKL_FB_PARA_BC_H_