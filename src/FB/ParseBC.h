#ifndef _FB_PARA_BC_H_
#define _FB_PARA_BC_H_

// #################################################################
//
// CAERU Library
//
// Copyright (c) All right reserved. 2012
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
//
// #################################################################

//@file ParseBC.h
//@brief FlowBase ParseBC class Header
//@author keno

#include "Skl.h"
#include "SklSolverBase.h"
#include "FB_Define.h"
#include "FBUtility.h"
#include "BndOuter.h"
#include "string.h"
#include "Control.h"
#include "Component.h"
#include "IDtable.h"
#include "parallel/SklParaComponent.h"
#include "Medium.h"
#include "vec3.h"
#include "Parallel_node.h"
#include "Intrinsic.h"
#include "TPControl.h"

class ParseBC : public Parallel_Node {
private:
  REAL_TYPE RefVelocity, BaseTemp, DiffTemp, RefDensity, RefSpecificHeat;
  REAL_TYPE RefLength, BasePrs;
  REAL_TYPE rho, nyu, cp, lambda, beta; // 無次元化の参照値
  
  unsigned ix, jx, kx, guide;
  unsigned KindOfSolver;
  unsigned NoCompo, NoBC, NoMedium;
  unsigned Unit_Param, monitor, Unit_Temp, Unit_Prs, Mode_Gradp;
  bool isCDS;
  
  CompoList*     compo;
  BoundaryOuter* bc;
	
	bool HeatProblem;
  char OBCname[NOFACE][LABEL];
	
public:
	IDtable* iTable;
  
  // Medium Table
  MediumTableInfo *MTITP; //Medium Table <--- textparser

public:
  ParseBC(){
    ix = jx = kx = guide = 0;
    KindOfSolver = 0;
    BaseTemp = DiffTemp = BasePrs = 0.0;
    RefVelocity = RefDensity = RefSpecificHeat = RefLength = 0.0;
    NoCompo = NoBC = NoMedium = monitor = 0;
    Unit_Param = 0;
    Unit_Temp = 0;
    Unit_Prs = 0;
    Mode_Gradp = 0;
    HeatProblem = isCDS = false;
    bc = NULL;
    compo = NULL;
    iTable = NULL;
    for (int i=0; i<NOFACE; i++) {
      strcpy(&OBCname[i][0], "\0");
    }
  }
  ~ParseBC() {
    if (iTable) delete [] iTable;
  }
  
protected:
  bool chkID              (void);
  bool isComponent        (unsigned label);
  bool isCompoTransfer    (unsigned label);
  bool isIDinTable        (int candidate);
  
  int getStateinTable     (unsigned id);
  
  unsigned oppositDir         (unsigned dir);
  
  void chkKeywordIBC        (const char *keyword, unsigned m);
  void chkKeywordOBC        (const char *keyword, unsigned m);
  void dbg_printBaseOBC     (FILE* fp);
  void get_NV               (const std::string label_base, unsigned n, REAL_TYPE* v);
  void get_Dir              (const std::string label_base, unsigned n, REAL_TYPE* v);
  void getUnitVec           (REAL_TYPE* v);
  void printCompo           (FILE* fp, REAL_TYPE* nv, int* ci, MediumList* mat);
  void printFaceOBC         (FILE* fp, REAL_TYPE* G_Lbx);
  void printOBC             (FILE* fp, BoundaryOuter* ref, REAL_TYPE* G_Lbx, unsigned face);
  void set_Deface           (const CfgElem *elmL, unsigned n, const char* str);
  
  const CfgElem* selectFace (int face, const CfgElem* elmTop);
  
  std::string getOBCstr     (unsigned id);
  
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
  bool receive_TP_Ptr     (TPControl* tp);
  
  int getNoLocalBC        (void);
  
  void construct_iTable   (Control* C);
  void chkBCconsistency   (unsigned kos);
  void setMediumPoint     (MediumTableInfo *m_MTITP);
  void printCompoInfo     (FILE* mp, FILE* fp, REAL_TYPE* nv, int* ci, MediumList* mat);
  void printOBCinfo       (FILE* mp, FILE* fp, REAL_TYPE* G_Lbx);
  void printTable         (FILE* fp);
  void receiveCompoPtr    (CompoList* CMP);
  void setMedium          (Control* Cref);
  void setControlVars     (Control* Cref);
  void setRefMedium       (MediumList* mat, Control* Cref);
  void setRefValue        (MediumList* mat, CompoList* cmp, Control* C);
  
  IDtable* get_IDtable_Ptr(void) {
    return iTable;
  }
  
  
  //for text parser
  
protected:
  
  // TPControl
  TPControl* tpCntl;
  
  unsigned getTP_Vel_profile (const std::string label_base);
  void getTP_IBC_Adiabatic (const std::string label_base, unsigned n);
  void getTP_IBC_CnstTemp  (const std::string label_base, unsigned n);
  void getTP_IBC_Fan       (const std::string label_base, unsigned n);
  void getTP_IBC_IBM_DF    (const std::string label_base, unsigned n);
  void getTP_IBC_HeatFlux  (const std::string label_base, unsigned n);
  void getTP_IBC_HeatSrc   (const std::string label_base, unsigned n);
  void getTP_IBC_HT_B      (const std::string label_base, unsigned n);
  void getTP_IBC_HT_N      (const std::string label_base, unsigned n);
  void getTP_IBC_HT_S      (const std::string label_base, unsigned n);
  void getTP_IBC_HT_SF     (const std::string label_base, unsigned n);
  void getTP_IBC_HT_SN     (const std::string label_base, unsigned n);
  void getTP_IBC_IsoTherm  (const std::string label_base, unsigned n);
  void getTP_IBC_Monitor   (const std::string label_base, unsigned n, Control* C);
  void getTP_IBC_Outflow   (const std::string label_base, unsigned n);
  void getTP_IBC_Periodic  (const std::string label_base, unsigned n);
  void getTP_IBC_PrsLoss   (const std::string label_base, unsigned n);
  void getTP_IBC_Radiant   (const std::string label_base, unsigned n);
  void getTP_IBC_SpecVel   (const std::string label_base, unsigned n);
  void getTP_OBC_FarField  (const std::string label_base, unsigned n);
  void getTP_OBC_HT        (const std::string label_base, unsigned n, string kind);
  void getTP_OBC_Outflow   (const std::string label_base, unsigned n);
  void getTP_OBC_Periodic  (const std::string label_base, unsigned n);
  void getTP_OBC_SpecVH    (const std::string label_base, unsigned n);
  void getTP_OBC_Trcfree   (const std::string label_base, unsigned n);
  void getTP_OBC_Wall      (const std::string label_base, unsigned n);
  void getTP_Vel_Params    (const std::string label_base, unsigned type, REAL_TYPE* ca, REAL_TYPE vel);
  void getTP_Darcy         (const std::string label_base, unsigned n);
  void setTP_Deface        (const std::string label_base, unsigned n);
  void getTP_Center        (const std::string label_base, unsigned n, REAL_TYPE* v);
  
public:
  
  int getTP_BCval_int(const std::string label);
  REAL_TYPE getTP_BCval_real(const std::string label);
  
  void TPloadOuterBC        (void);
  void TPsetCompoList       (Control* C);
  void TPsetObcPtr        (BoundaryOuter* ptr);
  void getTP_Phase       (void);
  void getTP_Medium_InitTemp(void);
  
  void TPcount_Outer_Cell_ID (int* medium);
  
};

#endif // _FB_PARA_BC_H_
