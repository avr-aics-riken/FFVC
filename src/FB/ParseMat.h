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
  
  unsigned NoCompo, NoBC, Unit_Temp;
  
  MediumList* mat;

  int NoMedium;
  MediumTableInfo *MTITP;
  
public:
  ParseMat() {
    NoCompo    = 0;
    NoBC       = 0;
    NoMedium   = 0;
    Unit_Temp  = 0;
    mat        = NULL;
    MTITP      = NULL;
    for (int i=0; i<property_END; i++) ChkList[i]=false;
  }
  ~ParseMat() {}
  
protected:
  TPControl* tpCntl;
  
  bool chkDuplicateLabel   (const int n, std::string m_label);
  bool chkList4Solver      (const int m);
  
  int missingMessage       (const int m, const int key);
  
  void chkList             (FILE* fp, CompoList* compo, unsigned basicEq);
  void printMatList        (FILE* fp, int Max, MediumList* mlist);
  void printRelation       (FILE* fp, CompoList* compo);
  void storeProperty       (const int n, const int matid);
    
public:
  bool chkStateList      (CompoList* compo);
  bool receive_TP_Ptr    (TPControl* tp);
  
  int get_MediumTable    (void);
  
  void chkList           (FILE* mp, FILE* fp, CompoList* compo, unsigned basicEq);
  void chkState_Mat_Cmp  (CompoList* compo, FILE* fp);
  
  //void makeLinkCmpMat       (CompoList* compo);
  void makeMediumList    (MediumList* m_mat);
  void printMediumList   (FILE* mp, FILE* fp);
  void setControlVars    (Control* Cref);
  void set_matMedium     (const unsigned n, const int matid);
  
  
  // ----------> debug function
  void dbg_printRelation         (FILE* mp, FILE* fp, CompoList* compo); 

};

#endif // _FB_PARSE_M_H_
