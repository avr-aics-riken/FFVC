#ifndef _SKL_FB_IDTABLE_H_
#define _SKL_FB_IDTABLE_H_

/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file IDtable.h
//@brief FlowBase IDtable class Header
//@author keno, FSI Team, VCAD, RIKEN

#include "FBDefine.h"

class IDtable {
private:
  int       state;
  unsigned  id, mat_id;
  char      label[LABEL];
  
public:
  IDtable() {
    state = 0;
    id = mat_id = 0;
    for (unsigned n=0; n<LABEL; n++) label[n]='\0';
  }
  ~IDtable() {}
  
public:
    
  /**
   @fn char* getLabel(void)
   @brief labelを返す
   @retval labelのアドレス
   */
  char* getLabel(void) { return label; }
  
  /**
   @fn unsigned getID(void)
   @brief idを返す
   @retval id
   */
  unsigned getID(void) { return id; }
  
  /**
   @fn unsigned getMatID(void)
   @brief mat_idを返す
   @retval mat_id
   */
  unsigned getMatID(void) { return mat_id; }
  
  /**
   @fn int getState(void)
   @brief stateを返す
   @retval state
   */
  int getState(void) { return state; }
  
  /**
   @fn void setID(unsigned i)
   @brief iをidstateにセットする
   */
  void setID(unsigned i) {
    id = i;
  }
  
  /**
   @fn void setMatID(unsigned i)
   @brief iをidstateにセットする
   */
  void setMatID(unsigned i) {
    mat_id = i;
  }
  
  /**
   @fn void setLabel(const char* pnt)
   @brief ラベル名をセットする
   */
  void setLabel(const char* pnt) {
    strcpy(label, pnt);
  }
  
  /**
   @fn void IDtable::setState(int s)
   @brief sをstateにセットする
   */
  void setState(int s) {
    state = s;
  }
};

#endif // _SKL_FB_IDTABLE_H_
