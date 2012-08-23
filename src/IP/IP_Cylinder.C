// #################################################################
//
// FFV : Frontflow / violet Cartesian
//
// Copyright (c) 2012 All right reserved.
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
//
// #################################################################

/** 
 * @file   IP_Cylinder.C
 * @brief  IP_Cylinder class
 * @author kero
 */

#include "IP_Cylinder.h"


// #################################################################
// パラメータを取得する
bool IP_Cylinder::getTP(Control* R, TPControl* tpCntl)
{
  std::string str;
  std::string label;
  REAL_TYPE ct;
  
  // 2D or 3D mode
  label = "/Parameter/Intrinsic_Example/mode";
  
  if ( !(tpCntl->GetValue(label, &str )) )
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
  label = "/Parameter/Intrinsic_Example/width";
  
  if ( !(tpCntl->GetValue(label, &ct )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  else
  {
	  width = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL;
  }
  
  // z-dir step
  label="/Parameter/Intrinsic_Example/height";
  
  if ( !(tpCntl->GetValue(label, &ct )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  else
  {
	  height = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL;
  }
  
  // ドライバの設定 値が正の値のとき，有効．ゼロの場合はドライバなし
  label = "/Parameter/Intrinsic_Example/Driver";
  
  if ( tpCntl->GetValue(label, &ct ) )
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
    Hostonly_ stamped_printf("\tError : Value of 'Driver' in 'Intrinsic_Example' must be positive.\n");
    return false;
  }
  
  if ( drv_length > 0.0 )
  {
    drv_mode = ON;
  }
  else{
    drv_mode = OFF;
  }
  
  
  // 媒質指定
  label = "/Parameter/Intrinsic_Example/Fluid_medium";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  m_fluid = str;
  
  
  label = "/Parameter/Intrinsic_Example/Solid_medium";
  
  if ( !(tpCntl->GetValue(label, &str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  m_solid = str;
  
  if ( drv_mode == ON )
  {
    label = "/Parameter/Intrinsic_Example/driver_medium";
    
    if ( !(tpCntl->GetValue(label, &str )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      return false;
    }
    m_driver = str;
    
    
    label = "/Parameter/Intrinsic_Example/driver_face_medium";
    
    if ( !(tpCntl->GetValue(label, &str )) )
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
void IP_Cylinder::printPara(FILE* fp, const Control* R)
{
  if ( !fp )
  {
    stamped_printf("\tFail to write into file\n");
    Exit(0);
  }
  
  fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
  fprintf(fp,"\n\t>> Intrinsic Duct Parameters\n\n");
  
  fprintf(fp,"\tDimension Mode                     :  %s\n", (mode==dim_2d)?"2 Dimensional":"3 Dimensional");
  fprintf(fp,"\tStep Width (x-dir.)    [m] / [-]   : %12.5e / %12.5e\n", width, width/RefL);
  fprintf(fp,"\tStep Height(z-dir.)    [m] / [-]   : %12.5e / %12.5e\n", height, height/RefL);
  
  if ( drv_length > 0.0 )
  {
    fprintf(fp,"\tDriver Length        [m] / [-]   : %12.5e / %12.5e\n", drv_length, drv_length/RefL);
  }
}


// #################################################################
// Cylinderの領域情報を設定する
void IP_Cylinder::setDomain(Control* R, const int* sz, REAL_TYPE* org, REAL_TYPE* reg, REAL_TYPE* pch)
{
  RefL = R->RefLength;
  
  reg[0] = pch[0]*(REAL_TYPE)sz[0];
  reg[1] = pch[1]*(REAL_TYPE)sz[1];
  reg[2] = pch[2]*(REAL_TYPE)sz[2];
  
  // チェック
  if ( (pch[0] != pch[1]) || (pch[1] != pch[2]) )
  {
    Hostonly_ printf("Error : 'VoxelPitch' in each direction must be same.\n");
    Exit(0);
  }
  
  if ( ((int)(reg[0]/pch[0]) != sz[0]) ||
       ((int)(reg[1]/pch[1]) != sz[1]) ||
       ((int)(reg[2]/pch[2]) != sz[2]) )
  {
    Hostonly_ printf("Error : Invalid parameters among 'VoxelSize', 'VoxelPitch', and 'VoxelWidth' in DomainInfo section.\n");
    Exit(0);
  }

  // 次元とサイズ
  if (mode == dim_2d)
  {
    if (size[2] != 3)
    {
      Hostonly_ printf("Error : VoxelSize kmax must be 3 if 2-dimensional.\n");
    }
  }
}


// #################################################################
// Cylinderの計算領域のセルIDを設定する
void IP_Cylinder::setup(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, MediumList* mat)
{
  int mid_fluid;        ///< 流体
  int mid_solid;        ///< 固体
  int mid_driver;       ///< ドライバ部
  int mid_driver_face;  ///< ドライバ流出面
  REAL_TYPE x, y, z, dh, len, ht;
  REAL_TYPE ox, oy, oz, Lx, Ly, Lz;
  
  
  // ローカルにコピー
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  // 隣接ランクのIDを取得
  // nID[6]
  const int* nID = paraMngr->GetNeighborRankID();
  
  ox = origin[0];
  oy = origin[1];
  oz = origin[2];
  Lx = region[0];
  Ly = region[1];
  Lz = region[2];
  dh = deltaX;

  len= drv_length/R->RefLength;
  ht = height/R->RefLength;
  
  if ( (mid_fluid = R->find_ID_from_Label(mat, Nmax, m_fluid)) == 0 )
  {
    Hostonly_ printf("\tLabel '%s' is not listed in MediumList\n", m_fluid.c_str());
    Exit(0);
  }
  
  if ( (mid_solid = R->find_ID_from_Label(mat, Nmax, m_solid)) == 0 )
  {
    Hostonly_ printf("\tLabel '%s' is not listed in MediumList\n", m_solid.c_str());
    Exit(0);
  }
  
  
  // Initialize  内部領域をfluidにしておく
#pragma omp parallel for firstprivate(ix, jx, kx, gd, mid_fluid) schedule(static)
  for (int k=1; k<=kx; k++) { 
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        mid[m] = mid_fluid;
      }
    }
  }
  
  
  // driver設定 iff ドライバ長が正の場合
  if ( drv_mode == OFF ) return;
  
  
  if ( (mid_driver = R->find_ID_from_Label(mat, Nmax, m_driver)) == 0 )
  {
    Hostonly_ printf("\tLabel '%s' is not listed in MediumList\n", m_driver.c_str());
    Exit(0);
  }
  
  if ( (mid_driver_face = R->find_ID_from_Label(mat, Nmax, m_driver_face)) == 0 )
  {
    Hostonly_ printf("\tLabel '%s' is not listed in MediumList\n", m_driver_face.c_str());
    Exit(0);
  }
  
  // ドライバ部分
  if ( drv_length > 0.0 )
  {
    if ( nID[X_MINUS] < 0 )
    {
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          for (int i=1; i<=ix; i++) {
            size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
            x = ox + 0.5*dh + dh*(i-1);
            if ( x < ox+len ) mid[m] = mid_driver;
          }
        }
      }
    }     
  }
  
  
  // ドライバの下流面にIDを設定
  if ( drv_length > 0.0 )
  {
    for (int k=1; k<=kx; k++) {
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i,   j, k, ix, jx, kx, gd);
          size_t m1= _F_IDX_S3D(i+1, j, k, ix, jx, kx, gd);
          if ( (mid[m] == mid_driver) && (mid[m1] == mid_fluid) )
          {
            mid[m] = mid_driver_face;
          }
        }
      }
    }    
  }

  // ステップ部分を上書き
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        x = ox + 0.5*dh + dh*(i-1);
        y = oy + 0.5*dh + dh*(j-1);
        if ( (x < ox+len) && (y < oy+ht) )
        {
          mid[m] = mid_solid;
        }
      }
    }
  }
}
