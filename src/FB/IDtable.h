#ifndef _FB_ID_TABLE_H_
#define _FB_ID_TABLE_H_

// #################################################################
//
// CAERU Library
//
// Copyright (c) All right reserved. 2012
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
//
// #################################################################

//@file IDtable.h
//@brief FlowBase IDtable class Header
//@author keno

#include "FB_Define.h"

class IDtable {
private:
  int   state;
  int   mat_id;
  char  label[LABEL];
  
public:
  IDtable() {
    state = 0;
    mat_id = 0;
    for (unsigned n=0; n<LABEL; n++) label[n]='\0';
  }
  ~IDtable() {}
  
public:
    
  char* getLabel(void) { return label; }
  
  int getMatID(void) { return mat_id; }
  
  int getState(void) { return state; }
  
  void setMatID(int i) {
    mat_id = i;
  }
  
  void setLabel(const char* pnt) {
    strcpy(label, pnt);
  }
  
  void setState(int s) {
    state = s;
  }
};

#endif // _FB_ID_TABLE_H_
