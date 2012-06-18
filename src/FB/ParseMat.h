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
//@author keno

#include "FB_Define.h"
#include "Control.h"
#include "Medium.h"
#include "Component.h"
#include "Parallel_node.h"
#include "TPControl.h"

class  ParseMat : public Parallel_Node {
private:
  bool ChkList[property_END];  // # of parameters in MediumList must be less than # of property_END
  
  int NoMedium;    // 媒質数
  
  unsigned NoCompo;
  unsigned NoBC;
  unsigned Unit_Temp;
  unsigned KOS;
  
  MediumTableInfo *MTITP;
  
public:
  ParseMat() {
    NoCompo    = 0;
    NoBC       = 0;
    NoMedium   = 0;
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
  
  void chkList           (FILE* fp, CompoList* compo, unsigned basicEq);
  void printMatList      (FILE* fp, MediumList* mat);
  void printRelation     (FILE* fp, CompoList* compo, MediumList* mat);
  void copyProperty      (MediumList* mat, const int n);
    
public:
  bool chkStateList      (CompoList* compo, MediumList* mat);
  bool receive_TP_Ptr    (TPControl* tp);
  
  int get_MediumTable    (void);
  
  void chkList           (FILE* mp, FILE* fp, CompoList* compo, unsigned basicEq);
  void chkState_Mat_Cmp  (CompoList* compo, MediumList* mat, FILE* fp);
  
  //void makeLinkCmpMat       (CompoList* compo);
  void makeMediumList    (MediumList* mat);
  void printMediumList   (FILE* mp, FILE* fp, MediumList* mat);
  void setControlVars    (Control* Cref);
  
  
  // ----------> debug function
  void dbg_printRelation         (FILE* mp, FILE* fp, CompoList* compo, MediumList* mat); 

};

#endif // _FB_PARSE_M_H_
