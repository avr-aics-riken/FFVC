#ifndef _FB_POLYGRP_H_
#define _FB_POLYGRP_H_

//##################################################################################
//
// Flow Base class
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2014 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

/**
 * @file   PolyProperty.h
 * @brief  FlowBase PolygonProperty class Header
 * @author aics
 */

#include <string>
#include "FB_Define.h"
#include "Vec3.h"

using namespace std;
using namespace Vec3class;


// #################################################################
class PolygonProperty {
  
private:
  int l_ntria;       ///< ローカルなポリゴン数
  int g_ntria;       ///< グローバルなポリゴン数
  int m_id;          ///< ID
  REAL_TYPE l_area;  ///< ローカルな面積
  REAL_TYPE g_area;  ///< グローバルな面積
  string group;      ///< ポリゴングループ名
  string material;   ///< Mediumtable[@]のalias
  string bc;         ///< BCのラベル
  Vec3r bx_min;      ///< Bboxの最小値
  Vec3r bx_max;      ///< Bboxの最大値
  
public:
  PolygonProperty() {
    l_area = 0.0;
    l_ntria= 0;
    g_area = 0.0;
    g_ntria= 0;
  }
  
  ~PolygonProperty() {}
  
  string getGroup() const
  {
    return group;
  }
  
  int getID() const
  {
    return m_id;
  }
  
  string getMaterial() const
  {
    return material;
  }
  
  string getBClabel() const
  {
    return bc;
  }
  
  int getLntria() const
  {
    return l_ntria;
  }
  
  REAL_TYPE getLarea() const
  {
    return l_area;
  }
  
  int getGntria() const
  {
    return g_ntria;
  }
  
  REAL_TYPE getGarea() const
  {
    return g_area;
  }
  
  void setGroup(string key)
  {
    group = key;
  }
  
  void setID(int key)
  {
    m_id = key;
  }
  
  void setMaterial(string key)
  {
    material = key;
  }
  void setBClabel(string key)
  {
    bc = key;
  }
  
  void setLntria(int val)
  {
    l_ntria = val;
  }
  
  void setLarea(REAL_TYPE val)
  {
    l_area = val;
  }
  
  void setGntria(int val)
  {
    g_ntria = val;
  }
  
  void setGarea(REAL_TYPE val)
  {
    g_area = val;
  }
  
  Vec3r getBboxMax() const
  {
    return bx_max;
  }
  
  Vec3r getBboxMin() const
  {
    return bx_min;
  }
  
  void setBboxMax(const Vec3r bmax)
  {
    bx_max = bmax;
  }
  
  void setBboxMin(const Vec3r bmin)
  {
    bx_min = bmin;
  }
  
};

#endif // _FB_POLYGRP_H_
