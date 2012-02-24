/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file Control_Rect.C
//@brief FlowBase ControlRect class
//@author keno, FSI Team, VCAD, RIKEN

#include <math.h>
#include "Control_Rect.h"

/**
 @fn void ControlRect::getXMLTemporarySize(void)
 @brief テンポラリのサイズを与える
 */
void ControlRect::getXMLTemporarySize(void)
{  
  ParseSteer Tree(CF);
  if ( !Tree.IsSetElem("Temporary_Grid_Size") ) Exit(0);
  if ( !Tree.getEParam("jdim", jdim))   Exit(0);
  if ( !Tree.getEParam("kdim", kdim))   Exit(0);
  if ( !Tree.getEParam("ldim", ldim))   Exit(0);
}

/**
 @fn void ControlRect::getXMLInterfacialpoint(void)
 @brief インターフェイシャルポイントの数を取得する
 */
void ControlRect::getXMLInterfacialpoint(void)
{
  ParseSteer Tree(CF);
  if ( !Tree.getParam("No_Interfacial_Point", InterfacialPoints)) Exit(0);
}

/**
 @fn void ControlRect::getXMLGhostPoint(void)
 @brief ゴーストポイントの数を取得する
 */
void ControlRect::getXMLGhostPoint(void)
{
  ParseSteer Tree(CF);
  if ( !Tree.getParam("No_Ghost_Point", GhostPoints)) Exit(0);
}

/**
 @fn void ControlRect::getXMLGhostPoint(void)
 @brief 等分割領域の大きさを取得
 */
void ControlRect::getXMLEquallySpace(void)
{
  ParseSteer Tree(CF);
  if ( !Tree.IsSetElem("Equally_Spaced_Region") ) Exit(0);
  if ( !Tree.getEParam("x_dist", xdist))   Exit(0);
  if ( !Tree.getEParam("y_dist", ydist))   Exit(0);
  if ( !Tree.getEParam("z_dist", zdist))   Exit(0);
}

/**
 @fn void ControlRect::printSteerConditions(FILE* fp, ItrCtl* IC, DTcntl* DT, ReferenceFrame* RF)
 @brief 制御パラメータSTEERの表示
 */
void ControlRect::printSteerConditions(FILE* fp, ItrCtl* IC, DTcntl* DT, ReferenceFrame* RF)
{
  Control::printSteerConditions(fp, IC, DT, RF);
}

/**
 @fn void ControlRect::setInitialConditions(void)
 @brief 初期条件を設定する
 */
void ControlRect::setInitialConditions(void)
{
  REAL_TYPE a = (SpecificHeatRatio-1.0) + 0.5* (iv.VecU*iv.VecU + iv.VecV*iv.VecV + iv.VecW*iv.VecW);
  
  if ( Unit.Param == DIMENSIONAL ) {
    iv.Density  /= RefDensity;
    iv.Pressure /= SpecificHeatRatio;
    iv.Energy   /= a;
    iv.VecU     /= RefVelocity;
    iv.VecV     /= RefVelocity;
    iv.VecW     /= RefVelocity;
  }
}
