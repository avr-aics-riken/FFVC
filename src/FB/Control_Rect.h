#ifndef _SKL_FB_CONTROL_RECT_H_
#define _SKL_FB_CONTROL_RECT_H_

/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file Control_Rect.h
//@brief FlowBase ControlRect class Header
//@author keno, FSI Team, VCAD, RIKEN

#include "Control.h"

class ControlRect : public Control {
public:
  unsigned jdim, kdim, ldim;
  int GhostPoints, InterfacialPoints;
  REAL_TYPE xdist, ydist, zdist;
  
public:
  ControlRect(){
    jdim = kdim = ldim = 0;
    GhostPoints = InterfacialPoints = 0;
    xdist = ydist = zdist = 0.0;
  }
  ~ControlRect() {}
  
protected:
  virtual void printSteerConditions(FILE* fp, ItrCtl* IC, DTcntl* DT, ReferenceFrame* RF);
  
public:
  void getXMLEquallySpace(void);
  void getXMLGhostPoint(void);
  void getXMLInterfacialpoint(void);
  void getXMLTemporarySize(void);
  
  virtual void setInitialConditions(void);
};

#endif // _SKL_FB_CONTROL_RECT_H_
