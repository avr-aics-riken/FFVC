#ifndef _FB_PARA_BC_H_
#define _FB_PARA_BC_H_

/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file ParseBC.h
//@brief FlowBase ParseBC class Header
//@author keno, FSI Team, VCAD, RIKEN

#include "Skl.h"
#include "SklSolverBase.h"
#include "FB_Define.h"
#include "FBUtility.h"
#include "BndOuter.h"
#include "string.h"
#include "Control.h"
#include "Component.h"
#include "IDtable.h"
#include "config/SklSolverConfig.h"
#include "parallel/SklParaComponent.h"
#include "Material.h"
#include "vec3.h"
#include "Parallel_node.h"
#include "Intrinsic.h"

// TextParser 
#include "TPControl.h"

using namespace SklCfg;  // to use SklSolverConfig* cfg

class ParseBC : public Parallel_Node {
private:
  REAL_TYPE RefVelocity, BaseTemp, DiffTemp, RefDensity, RefSpecificHeat;
  REAL_TYPE RefLength, BasePrs;
  REAL_TYPE rho, nyu, cp, lambda, beta; // 無次元化の参照値
  
  unsigned ix, jx, kx, guide;
  unsigned KindOfSolver;
  unsigned NoCompo, NoBC, NoID, NoBaseBC;
  unsigned Unit_Param, monitor, Unit_Temp, Unit_Prs, Mode_Gradp;
  bool isCDS;
  
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
    HeatProblem = isCDS = false;
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
  void dbg_printBaseOBC     (FILE* fp);
  void getDarcy             (const CfgElem *elmL, unsigned n);
  void get_NV               (const CfgElem *elmL, unsigned n, const char* str, REAL_TYPE* v);
  void get_Dir              (const CfgElem *elmL, unsigned n, const char* str, REAL_TYPE* v);
  void get_Center           (const CfgElem *elmL, unsigned n, const char* str, REAL_TYPE* v);
  void getUnitVec           (REAL_TYPE* v);
  void getXML_IBC_Adiabatic (const CfgElem *elmL, unsigned n);
  void getXML_IBC_CnstTemp  (const CfgElem *elmL, unsigned n);
  void getXML_IBC_Fan       (const CfgElem *elmL, unsigned n);
  void getXML_IBC_IBM_DF    (const CfgElem *elmL, unsigned n);
  void getXML_IBC_HeatFlux  (const CfgElem *elmL, unsigned n);
  void getXML_IBC_HeatSrc   (const CfgElem *elmL, unsigned n);
  void getXML_IBC_HT_B      (const CfgElem *elmL, unsigned n);
  void getXML_IBC_HT_N      (const CfgElem *elmL, unsigned n);
  void getXML_IBC_HT_S      (const CfgElem *elmL, unsigned n);
  void getXML_IBC_HT_SF     (const CfgElem *elmL, unsigned n);
  void getXML_IBC_HT_SN     (const CfgElem *elmL, unsigned n);
  void getXML_IBC_IsoTherm  (const CfgElem *elmL, unsigned n);
  void getXML_IBC_Monitor   (const CfgElem *elmL, unsigned n, Control* C);
  void getXML_IBC_Outflow   (const CfgElem *elmL, unsigned n);
  void getXML_IBC_Periodic  (const CfgElem *elmL, unsigned n);
  void getXML_IBC_PrsLoss   (const CfgElem *elmL, unsigned n);
  void getXML_IBC_Radiant   (const CfgElem *elmL, unsigned n);
  void getXML_IBC_SpecVel   (const CfgElem *elmL, unsigned n);
  void getXML_OBC_FarField  (const CfgElem *elmL, unsigned n);
  void getXML_OBC_HT        (const CfgElem *elmL, unsigned n, const char* kind);
  void getXML_OBC_InOut     (const CfgElem *elmL, unsigned n);
  void getXML_OBC_Outflow   (const CfgElem *elmL, unsigned n);
  void getXML_OBC_Periodic  (const CfgElem *elmL, unsigned n);
  void getXML_OBC_SpecVH    (const CfgElem *elmL, unsigned n);
  void getXML_OBC_Trcfree   (const CfgElem *elmL, unsigned n);
  void getXML_OBC_Wall      (const CfgElem *elmL, unsigned n);
  void getXML_Vel_Params    (const CfgElem *elmL, unsigned type, REAL_TYPE* ca, REAL_TYPE vel, const char* err_str);
  void printCompo           (FILE* fp, REAL_TYPE* nv, int* ci, MaterialList* mat);
  void printFaceOBC         (FILE* fp, REAL_TYPE* G_Lbx);
  void printOBC             (FILE* fp, BoundaryOuter* ref, REAL_TYPE* G_Lbx, unsigned face);
  void set_Deface           (const CfgElem *elmL, unsigned n, const char* str);
  
  const CfgElem* selectFace (int face, const CfgElem* elmTop);
  
  string getOBCstr          (unsigned id);
  
  //@fn int getCmpGbbox_st_x(unsigned odr, int* gci)
  //@brief コンポーネントのBbox情報st_xを返す
  int getCmpGbbox_st_x(unsigned odr, int* gci) {
    return ( gci[6*odr+0] );
  }
  
  //@fn int getCmpGbbox_st_y(unsigned odr, int* gci)
  //@brief コンポーネントのBbox情報st_yを返す
  int getCmpGbbox_st_y(unsigned odr, int* gci) {
    return ( gci[6*odr+1] );
  }
  
  //@fn int getCmpGbbox_st_z(unsigned odr, int* gci)
  //@brief コンポーネントのBbox情報st_zを返す
  int getCmpGbbox_st_z(unsigned odr, int* gci) {
    return ( gci[6*odr+2] );
  }
  
  //@fn int getCmpGbbox_ed_x(unsigned odr, int* gci)
  //@brief コンポーネントのBbox情報st_xを返す
  int getCmpGbbox_ed_x(unsigned odr, int* gci) {
    return ( gci[6*odr+3] );
  }
  
  //@fn int getCmpGbbox_ed_y(unsigned odr, int* gci)
  //@brief コンポーネントのBbox情報ed_yを返す
  int getCmpGbbox_ed_y(unsigned odr, int* gci) {
    return ( gci[6*odr+4] );
  }
  
  //@fn int getCmpGbbox_ed_z(unsigned odr, int* gci)
  //@brief コンポーネントのBbox情報ed_zを返す
  int getCmpGbbox_ed_z(unsigned odr, int* gci) {
    return ( gci[6*odr+5] );
  }
  
  //@fn void copyVec(REAL_TYPE* dst, REAL_TYPE* src)
  //@brief ベクトルのコピー
  void copyVec(REAL_TYPE* dst, REAL_TYPE* src) {
    dst[0] = src[0];
    dst[1] = src[1];
    dst[2] = src[2];
  }
  
public:
  bool isIDinCompo        (unsigned candidate_id, unsigned now);
  bool isIDinCompo        (unsigned candidate_id, int def, unsigned now);
  
  int get_BCval_int       (const CfgElem *elmL, const char* key);
  
  REAL_TYPE get_BCval_real(const CfgElem *elmL, const char* key);
  
  void count_Outer_Cell_ID(int* medium);
  void chkBCconsistency   (unsigned kos);
  void getXML_Model       (void);
  void getXML_Phase       (void);
  void loadOuterBC        (void);
  void printCompoInfo     (FILE* mp, FILE* fp, REAL_TYPE* nv, int* ci, MaterialList* mat);
  void printOBCinfo       (FILE* mp, FILE* fp, REAL_TYPE* G_Lbx);
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
  
  
  //@fn IDtable* ParseBC::get_IDtable_Ptr(void)
  IDtable* get_IDtable_Ptr(void) {
    return iTable;
  }
  
  //@fn BoundaryOuter* get_BaseBC_Ptr(void)
  BoundaryOuter* get_BaseBC_Ptr(void) { 
    return BaseBc;
  }
  
  //@fn unsigned get_NoBaseBC(void)
  unsigned get_NoBaseBC(void) {
    return NoBaseBC;
  }
  
  
  //for text parser
  
protected:
  
  // TPControl
  TPControl* tpCntl;
  
  
  int getTP_NBC              (void);
  unsigned scanTPmodel       (void);
  unsigned getTP_Vel_profile (const string label_base);
  void getTP_IBC_Adiabatic (const string label_base, unsigned n);
  void getTP_IBC_CnstTemp  (const string label_base, unsigned n);
  void getTP_IBC_Fan       (const string label_base, unsigned n);
  void getTP_IBC_IBM_DF    (const string label_base, unsigned n);
  void getTP_IBC_HeatFlux  (const string label_base, unsigned n);
  void getTP_IBC_HeatSrc   (const string label_base, unsigned n);
  void getTP_IBC_HT_B      (const string label_base, unsigned n);
  void getTP_IBC_HT_N      (const string label_base, unsigned n);
  void getTP_IBC_HT_S      (const string label_base, unsigned n);
  void getTP_IBC_HT_SF     (const string label_base, unsigned n);
  void getTP_IBC_HT_SN     (const string label_base, unsigned n);
  void getTP_IBC_IsoTherm  (const string label_base, unsigned n);
  void getTP_IBC_Monitor   (const string label_base, unsigned n, Control* C);
  void getTP_IBC_Outflow   (const string label_base, unsigned n);
  void getTP_IBC_Periodic  (const string label_base, unsigned n);
  void getTP_IBC_PrsLoss   (const string label_base, unsigned n);
  void getTP_IBC_Radiant   (const string label_base, unsigned n);
  void getTP_IBC_SpecVel   (const string label_base, unsigned n);
  void getTP_OBC_FarField  (const string label_base, unsigned n);
  void getTP_OBC_HT        (const string label_base, unsigned n, string kind);
  //void getTP_OBC_InOut     (const string label_base, unsigned n);
  void getTP_OBC_Outflow   (const string label_base, unsigned n);
  void getTP_OBC_Periodic  (const string label_base, unsigned n);
  void getTP_OBC_SpecVH    (const string label_base, unsigned n);
  void getTP_OBC_Trcfree   (const string label_base, unsigned n);
  void getTP_OBC_Wall      (const string label_base, unsigned n);
  void getTP_Vel_Params    (const string label_base, unsigned type, REAL_TYPE* ca, REAL_TYPE vel);
  void getTP_Darcy         (const string label_base, unsigned n);
  void setTP_Deface        (const string label_base, unsigned n);
  void getTP_NV            (const string label_base, unsigned n, REAL_TYPE* v);
  void getTP_Dir              (const string label_base, unsigned n, REAL_TYPE* v);
  void getTP_Center           (const string label_base, unsigned n, REAL_TYPE* v);
  
public:
  
  
  int getTP_BCval_int(const string label);
  REAL_TYPE getTP_BCval_real(const string label);
  
  // Medium Table
  int nMedium_TableTP;
  int nMedium_TableDB;
  MediumTableInfo *MTITP;//Medium Table <--- textparser
  MediumTableInfo *MTIDB;//Medium Table <--- database
  
  void TPloadOuterBC        (void);
  void TPsetCompoList       (Control* C);
  void TPsetControlVars     (Control* Cref);
  void TPsetObcPtr        (BoundaryOuter* ptr);
  
  void getTP_Model       (Control* C);
  void getTP_Phase       (void);
  void getTP_Medium_InitTemp();
  
  void setMediumPoint(
                      int m_nMedium_TableTP,
                      int m_nMedium_TableDB,
                      MediumTableInfo *m_MTITP,
                      MediumTableInfo *m_MTIDB);
  
  bool receive_TP_Ptr(TPControl* tp);
  
  unsigned TPcount_Outer_Cell_ID (unsigned* medium);
  
};

#endif // _FB_PARA_BC_H_
