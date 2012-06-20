#ifndef _FB_PARSE_M_H_
#define _FB_PARSE_M_H_

// #################################################################
//
// CAERU Library
//
// Copyright (c) All right reserved. 2012
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
//
// #################################################################

//@file ParseMat.h
//@brief FlowBase ParseMat class Header
//@author kero

#include "FB_Define.h"
#include "Medium.h"
#include "Component.h"
#include "TPControl.h"
#include "FBUtility.h"

class  ParseMat {
private:
  bool ChkList[property_END];  // # of parameters in MediumList must be less than # of property_END
  
  int NoCompo;
  int NoBC;
  int Unit_Temp;
  int KOS;
  
  MediumTableInfo *MTITP;
  
public:
  ParseMat() {
    NoCompo    = 0;
    NoBC       = 0;
    Unit_Temp  = 0;
    KOS        = 0;
    MTITP      = NULL;
    for (int i=0; i<property_END; i++) ChkList[i]=false;
  }
  ~ParseMat() {
    if ( MTITP ) delete [] MTITP;
  }
  
protected:
  TPControl* tpCntl;
  
  bool chkDuplicateLabel (MediumList* mat, const int n, const std::string m_label);
  bool chkList4Solver    (MediumList* mat, const int m);
  
  int missingMessage     (MediumList* mat, const int m, const int key);
  
  void printRelation     (FILE* fp, CompoList* compo, MediumList* mat);
  void copyProperty      (MediumList* mat, const int n);
    
public:
  bool makeMediumList    (MediumList* mat, const int NoMedium);
  bool receive_TP_Ptr    (TPControl* tp);
  
  int get_MediumTable    (void);
  
  void chkList           (FILE* fp, CompoList* compo, const int basicEq);
  void chkState_Mat_Cmp  (CompoList* compo, MediumList* mat, FILE* fp);
  //void makeLinkCmpMat       (CompoList* compo);
  void printMatList      (FILE* fp, MediumList* mat, const int NoMedium);
  void setControlVars    (const int m_NoCompo,
                          const int m_NoBC,
                          const int m_Unit_Temp,
                          const int m_KOS);
  
  MediumTableInfo* export_MTI(void) {
    return MTITP;
  }
  
  // ----------> debug function
  void dbg_printRelation (FILE* mp, FILE* fp, CompoList* compo, MediumList* mat); 

};

#endif // _FB_PARSE_M_H_
