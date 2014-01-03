//##################################################################################
//
// FFV-C : Frontflow / violet Cartesian
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
 * @file   IP_PMT.C
 * @brief  IP_PMT class
 * @author kero
 */

#include "IP_PMT.h"


// #################################################################
//f パラメータをロード
bool IP_PMT::getTP(Control* R, TextParser* tpCntl)
{
  std::string str;
  std::string label;

  
  // 媒質指定
  label="/IntrinsicExample/FluidMedium";
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  m_fluid = str;
  
  label="/IntrinsicExample/SolidMedium";
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  m_solid = str;
  
  return true;
}


// #################################################################
// 領域パラメータを設定する
void IP_PMT::setDomainParameter(Control* R, const int* sz, REAL_TYPE* org, REAL_TYPE* reg, REAL_TYPE* pch)
{  
  // forced
  if (R->Unit.Param != NONDIMENSIONAL)
  {
    Hostonly_ printf("\tError : PMT class is designed for only non-dimensional parameter\n");
    Exit(0);
  }
  
  // 性能測定モードをOnにする
  R->Hide.PM_Test = ON;

  
  // 偶数チェック
  even = ON;
  
  org[0] = -0.5*reg[0];
  org[1] = -0.5*reg[1];
  org[2] = -0.5*reg[2];
}
