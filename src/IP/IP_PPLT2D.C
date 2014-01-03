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
 * @file   IP_PPLT2D.C
 * @brief  IP_PPLT2D class
 * @author kero
 */

#include "IP_PPLT2D.h"


// #################################################################
//パラメータをロード
bool IP_PPLT2D::getTP(Control* R, TextParser* tpCntl)
{
  std::string str;
  std::string label;
  
  // Dimension
  mode = dim_2d;
  
  // 媒質指定
  label = "/IntrinsicExample/FluidMedium";
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  m_fluid = str;
  
  label = "/IntrinsicExample/SolidMedium";
  
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
void IP_PPLT2D::setDomainParameter(Control* R, const int* sz, REAL_TYPE* org, REAL_TYPE* reg, REAL_TYPE* pch)
{
  // forced
  if (R->Unit.Param != NONDIMENSIONAL)
  {
    Hostonly_ printf("\tError : PPLT2D class is designed for only non-dimensional parameter\n");
    Exit(0);
  }
  
  // 2D
  mode = dim_2d;
  
  // 領域アスペクト比のチェック
  if ( sz[0] != sz[1]*2 )
  {
    Hostonly_ printf("\tError : The number of division must be 2:1 (=X:Y)\n");
    Exit(0);
  }

  // 二次元の領域設定
  reg[0] =  2.0;
  reg[1] =  1.0;
  org[0] = -1.0;
  org[1] =  0.0;
  
  pch[0] = reg[0] / (REAL_TYPE)sz[0];
  pch[1] = reg[1] / (REAL_TYPE)sz[1];
  
  // Z方向は3層に合わせて調整
  pch[2] =  pch[0];
  reg[2] =  pch[0]*3.0;
  org[2] = -pch[2]*1.5;
}
