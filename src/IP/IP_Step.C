//##################################################################################
//
// FFV-C : Frontflow / violet Cartesian
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2013 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2013 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

/** 
 * @file   IP_Step.C
 * @brief  IP_Step class
 * @author kero
 */

#include "IP_Step.h"


// #################################################################
// パラメータをロード
bool IP_Step::getTP(Control* R, TextParser* tpCntl)
{
  std::string str;
  std::string label;
  REAL_TYPE ct;
  
  // 2D or 3D mode
  label = "/IntrinsicExample/Dimension";
  
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
  
  // x-dir step
  label = "/IntrinsicExample/StepLength";
  if ( !(tpCntl->getInspectedValue(label, ct )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  else
  {
	  width = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL;
  }
  
  // z-dir step
  label = "/IntrinsicExample/StepHeight";
  if ( !(tpCntl->getInspectedValue(label, ct )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  else
  {
	  height = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL;
  }
  
  // ドライバの設定 値が正の値のとき，有効．ゼロの場合はドライバなし
  label = "/IntrinsicExample/DriverLength";
  if ( tpCntl->getInspectedValue(label, ct ) )
  {
    drv_length = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL;
  }
  else
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  
  if ( drv_length < 0.0 )
  {
    Hostonly_ stamped_printf("\tError : Value of '%s' must be positive.\n", label.c_str());
    return false;
  }
  
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
  
  // Only driver is specified
  if ( drv_length > 0.0 )
  {
    label = "/IntrinsicExample/DriverMedium";
    if ( !(tpCntl->getInspectedValue(label, str )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      return false;
    }
    m_driver = str;
    
    label = "/IntrinsicExample/DriverFaceMedium";
    if ( !(tpCntl->getInspectedValue(label, str )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      return false;
    }
    m_driver_face = str;
  }
  
  return true;
}


// #################################################################
// パラメータの表示
void IP_Step::printPara(FILE* fp, const Control* R)
{
  if ( !fp )
  {
    stamped_printf("\tFail to write into file\n");
    Exit(0);
  }
  
  fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
  fprintf(fp,"\n\t>> Intrinsic Backstep Class Parameters\n\n");
  
  fprintf(fp,"\tDimension Mode                     :  %s\n", (mode==dim_2d)?"2 Dimensional":"3 Dimensional");
  fprintf(fp,"\tStep Width (x-dir.)    [m] / [-]   : %12.5e / %12.5e\n", width, width/RefL);
  fprintf(fp,"\tStep Height(z-dir.)    [m] / [-]   : %12.5e / %12.5e\n", height, height/RefL);
  
  if ( drv_length > 0.0 )
  {
    fprintf(fp,"\tDriver Length        [m] / [-]   : %12.5e / %12.5e\n", drv_length, drv_length/RefL);
  }
}


// #################################################################
// 計算領域のセルIDを設定する
void IP_Step::setup(int* bcd, Control* R, REAL_TYPE* G_org, const int NoMedium, const MediumList* mat, float* cut, int* bid)
{
  int mid_fluid, mid_solid, mid_driver, mid_driver_face;
  
  // 流体
  if ( (mid_fluid = R->findIDfromLabel(mat, NoMedium, m_fluid)) == 0 )
  {
    Hostonly_ printf("\tLabel '%s' is not listed in MediumList\n", m_fluid.c_str());
    Exit(0);
  }
  
  // 固体
  if ( (mid_solid = R->findIDfromLabel(mat, NoMedium, m_solid)) == 0 )
  {
    Hostonly_ printf("\tLabel '%s' is not listed in MediumList\n", m_solid.c_str());
    Exit(0);
  }
  
  if ( drv_length > 0.0 )
  {
    // ドライバ部
    if ( (mid_driver = R->findIDfromLabel(mat, NoMedium, m_driver)) == 0 )
    {
      Hostonly_ printf("\tLabel '%s' is not listed in MediumList\n", m_driver.c_str());
      Exit(0);
    }
    
    // ドライバ流出面
    if ( (mid_driver_face = R->findIDfromLabel(mat, NoMedium, m_driver_face)) == 0 )
    {
      Hostonly_ printf("\tLabel '%s' is not listed in MediumList\n", m_driver_face.c_str());
      Exit(0);
    }
  }


  REAL_TYPE ox, oy, oz, dh;
  
  // ノードローカルの無次元値
  ox = origin[0];
  oy = origin[1];
  oz = origin[2];
  dh = deltaX;
  
  // ローカルにコピー
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  // length, widthなどは有次元値
  REAL_TYPE len = G_origin[0] + (drv_length+width)/R->RefLength; // グローバルな無次元位置
  REAL_TYPE ht  = G_origin[1] + height/R->RefLength;

  
  // ドライバ部分　X-面からドライバ長さより小さい領域
  if ( drv_length > 0.0 )
  {
    
#pragma omp parallel for firstprivate(ix, jx, kx, gd, mid_driver, ox, dh, len) schedule(static)
    for (int k=1; k<=kx; k++) {
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
          REAL_TYPE x = ox + 0.5*dh + dh*(i-1);
          if ( x < len ) bcd[m] |= mid_driver;
        }
      }
    }
    
  }
  
  // ドライバの下流面にIDを設定
  if ( drv_length > 0.0 )
  {
    
#pragma omp parallel for firstprivate(ix, jx, kx, gd, mid_driver, mid_fluid, mid_driver_face) schedule(static)
    for (int k=1; k<=kx; k++) {
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i,   j, k, ix, jx, kx, gd);
          size_t m1= _F_IDX_S3D(i+1, j, k, ix, jx, kx, gd);
          if ( (DECODE_CMP(bcd[m])  == mid_driver) && (DECODE_CMP(bcd[m1]) == mid_fluid) ) {
            bcd[m] |= mid_driver_face;
          }
        }
      }
    }
    
  }

  // ステップ部分を上書き
#pragma omp parallel for firstprivate(ix, jx, kx, gd, mid_solid, ox, oy, dh, len, ht) schedule(static)
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        REAL_TYPE x = ox + 0.5*dh + dh*(i-1);
        REAL_TYPE y = oy + 0.5*dh + dh*(j-1);
        if ( (x < len) && (y < ht) )
        {
          bcd[m] |= mid_solid;
        }
      }
    }
  }
  
}
