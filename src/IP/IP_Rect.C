//##################################################################################
//
// FFV-C : Frontflow / violet Cartesian
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2015 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2016 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
// Copyright (c) 2016-2020 Research Institute for Information Research, Kyushu University.
// All rights reserved.
//
//##################################################################################

/*
 * @file   IP_Rect.C
 * @brief  IP_Rect class
 * @author aics
 */

#include "IP_Rect.h"


// #################################################################
// パラメータをロード
bool IP_Rect::getTP(Control* R, TextParser* tpCntl)
{
  std::string str;
  std::string label;
  
  // 2D or 3D mode
  label="/IntrinsicExample/Dimension";
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  
  if     ( !strcasecmp(str.c_str(), "2d") )
  {
    mode = dim_2d;
  }
  else if( !strcasecmp(str.c_str(), "3d") )
  {
    mode = dim_3d;
  }
  else
  {
    Hostonly_ stamped_printf("\tParsing error : Invalid '%s'\n", label.c_str());
    return false;
  }
  
  
  // 分割数の偶数チェックオプション
  label="/IntrinsicExample/CheckEven";

  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }

  if     ( !strcasecmp(str.c_str(), "yes") )
  {
    even = ON;
  }
  else if( !strcasecmp(str.c_str(), "no") )
  {
    even = OFF;
  }
  else
  {
    Hostonly_ stamped_printf("\tParsing error : Invalid '%s'\n", label.c_str());
    return false;
  }
  
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
// パラメータの表示
void IP_Rect::printPara(FILE* fp, const Control* R)
{
  if ( !fp ) {
    stamped_printf("\tFail to write into file\n");
    Exit(0);
  }
  
  fprintf(fp,"\n----------\n\n");
  fprintf(fp,"\n\t>> Intrinsic Rectangular Class Parameters\n\n");
  
  fprintf(fp,"\tDimension Mode                     :  %s\n", (mode==dim_2d)?"2 Dimensional":"3 Dimensional");
  fprintf(fp,"\tCheck even number                  :  %s\n", (even == ON)?"Yes":"No");

}
