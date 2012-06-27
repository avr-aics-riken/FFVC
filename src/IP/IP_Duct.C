// #################################################################
//
// CAERU Library
//
// Copyright (c) All right reserved. 2012
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
//
// #################################################################

/** 
 * @file IP_Duct.C
 * @brief IP_Duct class
 * @author kero
 */

#include "IP_Duct.h"


//@brief パラメータを取得する
bool IP_Duct::getTP(Control* R, TPControl* tpCntl)
{
  std::string str;
  std::string label;
  REAL_TYPE ct;
  
  // Shape
  label="/Parameter/Intrinsic_Example/Shape";
  if ( !(tpCntl->GetValue(label, &str )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Shape' in 'Intrinsic_Example'\n");
    return false;
  }
  
  if ( !strcasecmp(str.c_str(), "circular") ) {
    driver.shape = id_circular;
  }
  else if ( !strcasecmp(str.c_str(), "rectangualr") ) {
    driver.shape = id_rectangular;
  }
  else {
    Hostonly_ stamped_printf("\tParsing error : Invalid shape in 'Intrinsic_Example'\n");
    return false;
  }
  
  // Diameter
  label="/Parameter/Intrinsic_Example/Diameter";
  if ( !(tpCntl->GetValue(label, &ct )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Diameter' in 'Intrinsic_Example'\n");
    return false;
  }
  else{
	  driver.diameter = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL;
  }
  
  // periodic
  label="/Parameter/Intrinsic_Example/Direction";
  if ( !(tpCntl->GetValue(label, &str )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Direction' in 'Intrinsic_Example'\n");
    return false;
  }
  if ( !strcasecmp(str.c_str(), "X_minus")) {
    driver.direction = X_MINUS;
  }
  else if ( !strcasecmp(str.c_str(), "X_plus")) {
    driver.direction = X_PLUS;
  }
  else if ( !strcasecmp(str.c_str(), "Y_minus")) {
    driver.direction = Y_MINUS;
  }
  else if ( !strcasecmp(str.c_str(), "Y_plus")) {
    driver.direction = Y_PLUS;
  }
  else if ( !strcasecmp(str.c_str(), "Z_minus")) {
    driver.direction = Z_MINUS;
  }
  else if ( !strcasecmp(str.c_str(), "Z_plus")) {
    driver.direction = Z_PLUS;
  }
  else {
    Hostonly_ stamped_printf("\tParsing error : Invalid value of 'Direction' in 'Intrinsic_Example'\n");
    return false;
  }     
  
  // ドライバの設定 値が正の値のとき，有効．ゼロの場合はドライバなし
  label="/Parameter/Intrinsic_Example/Driver";
  if ( tpCntl->GetValue(label, &ct ) ) {
    driver.length = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL;
  }
  else {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Driver' in 'Intrinsic_Example'\n");
    return false;
  }
  
  if ( driver.length < 0.0 ) {
    Hostonly_ stamped_printf("\tError : Value of 'Driver' in 'Intrinsic_Example' must be positive.\n");
    return false;
  }
  
  // 媒質指定
  label="/Parameter/Intrinsic_Example/Fluid_medium";
  if ( !(tpCntl->GetValue(label, &str )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Fluid_medium' in 'Intrinsic_Example'\n");
    return false;
  }
  m_fluid = str;
  
  label="/Parameter/Intrinsic_Example/Solid_medium";
  if ( !(tpCntl->GetValue(label, &str )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Solid_medium' in 'Intrinsic_Example'\n");
    return false;
  }
  m_solid = str;
  
  label="/Parameter/Intrinsic_Example/driver_medium";
  if ( !(tpCntl->GetValue(label, &str )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'driver_medium' in 'Intrinsic_Example'\n");
    return false;
  }
  m_driver = str;
  
  label="/Parameter/Intrinsic_Example/driver_face_medium";
  if ( !(tpCntl->GetValue(label, &str )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'driver_face_medium' in 'Intrinsic_Example'\n");
    return false;
  }
  m_driver_face = str;
  
  return true;
}


// パラメータの表示
void IP_Duct::printPara(FILE* fp, Control* R)
{
  if ( !fp ) {
    stamped_printf("\tFail to write into file\n");
    Exit(0);
  }
  
  std::string dir;
  switch (driver.direction) {
    case X_MINUS:
    case X_PLUS:
      dir = "X dir.";
      break;
      
    case Y_MINUS:
    case Y_PLUS:
      dir = "Y dir.";
      break;
      
    case Z_MINUS:
    case Z_PLUS:
      dir = "Z dir.";
      break;
  }
  
  fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
  fprintf(fp,"\n\t>> Intrinsic Duct Parameters\n\n");
  
  fprintf(fp,"\tShape of cross-section                :  %s\n", (driver.shape==id_circular)?"Circular":"Rectangular");
  fprintf(fp,"\tDiameter                  [m] / [-]   : %12.5e / %12.5e\n", driver.diameter, driver.diameter/RefL);
  fprintf(fp,"\tDirection of periodic BC              :  %s\n", dir.c_str());
  if ( driver.length > 0.0 ) {
    fprintf(fp,"\twith Driver                           :  %s\n", FBUtility::getDirection(driver.direction).c_str());
    fprintf(fp,"\t     Driver Length        [m] / [-]   : %12.5e / %12.5e\n", driver.length, driver.length/RefL);
  }
  
}


// Ductの領域情報を設定する
void IP_Duct::setDomain(Control* R, const int* sz, const REAL_TYPE* org, const REAL_TYPE* reg, const REAL_TYPE* pch)
{
  reg[0] = pch[0]*(REAL_TYPE)sz[0];
  reg[1] = pch[1]*(REAL_TYPE)sz[1];
  reg[2] = pch[2]*(REAL_TYPE)sz[2];
  
  // チェック
  if ( (pch[0] != pch[1]) || (pch[1] != pch[2]) ) {
    Hostonly_ printf("Error : 'VoxelPitch' in each direction must be same.\n");
    Exit(0);
  }
  if ( ((int)(reg[0]/pch[0]) != sz[0]) ||
       ((int)(reg[1]/pch[1]) != sz[1]) ||
       ((int)(reg[2]/pch[2]) != sz[2]) ) {
    Hostonly_ printf("Error : Invalid parameters among 'VoxelSize', 'VoxelPitch', and 'VoxelWidth' in DomainInfo section.\n");
    Exit(0);
  }

}


// Ductの計算領域のセルIDを設定する
void IP_Duct::setup(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, MediumList* mat)
{
  int mid_fluid=1;        /// 流体
  int mid_solid=2;        /// 固体
  int mid_driver=3;       /// ドライバ部
  int mid_driver_face=4;  /// ドライバ流出面
  REAL_TYPE x, y, z, dh, r, len;
  REAL_TYPE ox, oy, oz, Lx, Ly, Lz;
  
  size_t m;
  
  // ローカルにコピー
  int imax = size[0];
  int jmax = size[1];
  int kmax = size[2];
  int gd = guide;
  
  // 隣接ランクのIDを取得 nID[6]
  const int* nID = paraMngr->GetNeighborRankID();
  
  ox = R->org[0];
  oy = R->org[1];
  oz = R->org[2];
  Lx = R->Lbx[0];
  Ly = R->Lbx[1];
  Lz = R->Lbx[2];
  dh = R->dh;
  r  = driver.diameter/R->RefLength * 0.5;
  len= driver.length/R->RefLength;
  
  // Initialize  内部領域をsolidにしておく
  for (int k=1; k<=kmax; k++) {
    for (int j=1; j<=jmax; j++) {
      for (int i=1; i<=imax; i++) {
        m = FBUtility::getFindexS3D(size, guide, i, j, k);
        mid[m] = mid_solid;
      }
    }
  }
  
  // Inner
  if (driver.shape == id_rectangular ) { // 矩形管の場合は内部は全て流体
    for (int k=1; k<=kmax; k++) {
      for (int j=1; j<=jmax; j++) {
        for (int i=1; i<=imax; i++) {
          m = FBUtility::getFindexS3D(size, guide, i, j, k);
          mid[m] = mid_fluid;
        }
      }
    }
  }
  else { // 円管の場合，半径以下のセルを流体にする（ノードにかかわらず）
    switch (driver.direction) {
      case X_MINUS:
      case X_PLUS:
        for (int k=1; k<=kmax; k++) {
          for (int j=1; j<=jmax; j++) {
            for (int i=1; i<=imax; i++) {
              m = FBUtility::getFindexS3D(size, guide, i, j, k);
              y = oy + 0.5*dh + dh*(j-1);
              z = oz + 0.5*dh + dh*(k-1);
              if ( (y-r)*(y-r)+(z-r)*(z-r) <= r*r ) mid[m] = mid_fluid;
            }
          }
        }
        break;
        
      case Y_MINUS:
      case Y_PLUS:
        for (int k=1; k<=kmax; k++) {
          for (int j=1; j<=jmax; j++) {
            for (int i=1; i<=imax; i++) {
              m = FBUtility::getFindexS3D(size, guide, i, j, k);
              x = ox + 0.5*dh + dh*(i-1);
              z = oz + 0.5*dh + dh*(k-1);
              if ( (x-r)*(x-r)+(z-r)*(z-r) <= r*r ) mid[m] = mid_fluid;
            }
          }
        }
        break;
        
      case Z_MINUS:
      case Z_PLUS:
        for (int k=1; k<=kmax; k++) {
          for (int j=1; j<=jmax; j++) {
            for (int i=1; i<=imax; i++) {
              m = FBUtility::getFindexS3D(size, guide, i, j, k);
              x = ox + 0.5*dh + dh*(i-1);
              y = oy + 0.5*dh + dh*(j-1);
              if ( (x-r)*(x-r)+(y-r)*(y-r) <= r*r ) mid[m] = mid_fluid;
            }
          }
        }
        break;
    }   
  }
  
  // ドライバ部分
  if ( driver.length > 0.0 ) {
    
    switch (driver.direction) {
      case X_MINUS:
        if ( nID[driver.direction] < 0 ) {
          for (int k=1; k<=kmax; k++) {
            for (int j=1; j<=jmax; j++) {
              for (int i=1; i<=imax; i++) {
                m = FBUtility::getFindexS3D(size, guide, i, j, k);
                x = ox + 0.5*dh + dh*(i-1);
                y = oy + 0.5*dh + dh*(j-1);
                z = oz + 0.5*dh + dh*(k-1);
                if ( x < ox+len ) {
                  if ( driver.shape == id_circular ) {
                    if ( (y-r)*(y-r)+(z-r)*(z-r) <= r*r ) mid[m] = mid_driver; // 半径以内をドライバIDにする
                  }
                  else {
                    mid[m] = mid_driver;
                  }
                }
              }
            }
          }
        }        
        break;
        
      case X_PLUS:
        if ( nID[driver.direction] < 0 ) {
          for (int k=1; k<=kmax; k++) {
            for (int j=1; j<=jmax; j++) {
              for (int i=1; i<=imax; i++) {
                m = FBUtility::getFindexS3D(size, guide, i, j, k);
                x = ox + 0.5*dh + dh*(i-1);
                y = oy + 0.5*dh + dh*(j-1);
                z = oz + 0.5*dh + dh*(k-1);
                if ( x > ox+Lx-len ) {
                  if ( driver.shape == id_circular ) {
                    if ( (y-r)*(y-r)+(z-r)*(z-r) <= r*r ) mid[m] = mid_driver;
                  }
                  else {
                    mid[m] = mid_driver;
                  }
                }                
              }
            }
          }
        }
        break;
        
      case Y_MINUS:
        if ( nID[driver.direction] < 0 ) {
          for (int k=1; k<=kmax; k++) {
            for (int j=1; j<=jmax; j++) {
              for (int i=1; i<=imax; i++) {
                m = FBUtility::getFindexS3D(size, guide, i, j, k);
                x = ox + 0.5*dh + dh*(i-1);
                y = oy + 0.5*dh + dh*(j-1);
                z = oz + 0.5*dh + dh*(k-1);
                if ( y < oy+len ) {
                  if ( driver.shape == id_circular ) {
                    if ( (x-r)*(x-r)+(z-r)*(z-r) <= r*r ) mid[m] = mid_driver;
                  }
                  else {
                    mid[m] = mid_driver;
                  }
                }
              }
            }
          }
        }        
        break;
        
      case Y_PLUS:
        if ( nID[driver.direction] < 0 ) {
          for (int k=1; k<=kmax; k++) {
            for (int j=1; j<=jmax; j++) {
              for (int i=1; i<=imax; i++) {
                m = FBUtility::getFindexS3D(size, guide, i, j, k);
                x = ox + 0.5*dh + dh*(i-1);
                y = oy + 0.5*dh + dh*(j-1);
                z = oz + 0.5*dh + dh*(k-1);
                if ( y > oy+Ly-len ) {
                  if ( driver.shape == id_circular ) {
                    if ( (x-r)*(x-r)+(z-r)*(z-r) <= r*r ) mid[m] = mid_driver;
                  }
                  else {
                    mid[m] = mid_driver;
                  }
                }
              }
            }
          }
        }
        break;
        
      case Z_MINUS:
        if ( nID[driver.direction] < 0 ) {
          for (int k=1; k<=kmax; k++) {
            for (int j=1; j<=jmax; j++) {
              for (int i=1; i<=imax; i++) {
                m = FBUtility::getFindexS3D(size, guide, i, j, k);
                x = ox + 0.5*dh + dh*(i-1);
                y = oy + 0.5*dh + dh*(j-1);
                z = oz + 0.5*dh + dh*(k-1);
                if ( z < oz+len ) {
                  if ( driver.shape == id_circular ) {
                    if ( (x-r)*(x-r)+(y-r)*(y-r) <= r*r ) mid[m] = mid_driver;
                  }
                  else {
                    mid[m] = mid_driver;
                  }
                }
              }
            }
          }
        }
        break;
        
      case Z_PLUS:
        if ( nID[driver.direction] < 0 ) {
          for (int k=1; k<=kmax; k++) {
            for (int j=1; j<=jmax; j++) {
              for (int i=1; i<=imax; i++) {
                m = FBUtility::getFindexS3D(size, guide, i, j, k);
                x = ox + 0.5*dh + dh*(i-1);
                y = oy + 0.5*dh + dh*(j-1);
                z = oz + 0.5*dh + dh*(k-1);
                if ( z > oz+Lz-len ) {
                  if ( driver.shape == id_circular ) {
                    if ( (x-r)*(x-r)+(y-r)*(y-r) <= r*r ) mid[m] = mid_driver;
                  }
                  else {
                    mid[m] = mid_driver;
                  }
                }
              }
            }
          }
        }
        break;
    }    
  }
  
  // ドライバの下流面にIDを設定
  if ( driver.length > 0.0 ) {
    
    size_t m1;
    
    switch (driver.direction) 
    {
      case X_MINUS:
        for (int k=1; k<=kmax; k++) {
          for (int j=1; j<=jmax; j++) {
            for (int i=1; i<=imax; i++) {
              m = FBUtility::getFindexS3D(size, guide, i,   j, k);
              m1= FBUtility::getFindexS3D(size, guide, i+1, j, k);
              if ( (mid[m] == mid_driver) && (mid[m1] == mid_fluid) ) 
              {
                mid[m] = mid_driver_face;
              }
            }
          }
        }        
        break;
        
      case X_PLUS:
        for (int k=1; k<=kmax; k++) {
          for (int j=1; j<=jmax; j++) {
            for (int i=1; i<=imax; i++) {
              m = FBUtility::getFindexS3D(size, guide, i,   j, k);
              m1= FBUtility::getFindexS3D(size, guide, i-1, j, k);
              if ( (mid[m] == mid_driver) && (mid[m1] == mid_fluid) ) 
              {
                mid[m] = mid_driver_face;
              }
            }
          }
        }
        break;
        
      case Y_MINUS:
        for (int k=1; k<=kmax; k++) {
          for (int j=1; j<=jmax; j++) {
            for (int i=1; i<=imax; i++) {
              m = FBUtility::getFindexS3D(size, guide, i, j,   k);
              m1= FBUtility::getFindexS3D(size, guide, i, j+1, k);
              if ( (mid[m] == mid_driver) && (mid[m1] == mid_fluid) ) 
              {
                mid[m] = mid_driver_face;
              }
            }
          }
        }        
        break;
        
      case Y_PLUS:
        for (int k=1; k<=kmax; k++) {
          for (int j=1; j<=jmax; j++) {
            for (int i=1; i<=imax; i++) {
              m = FBUtility::getFindexS3D(size, guide, i, j,   k);
              m1= FBUtility::getFindexS3D(size, guide, i, j-1, k);
              if ( (mid[m] == mid_driver) && (mid[m1] == mid_fluid) ) 
              {
                mid[m] = mid_driver_face;
              }
            }
          }
        }
        break;
        
      case Z_MINUS:
        for (int k=1; k<=kmax; k++) {
          for (int j=1; j<=jmax; j++) {
            for (int i=1; i<=imax; i++) {
              m = FBUtility::getFindexS3D(size, guide, i, j, k);
              m1= FBUtility::getFindexS3D(size, guide, i, j, k+1);
              if ( (mid[m] == mid_driver) && (mid[m1] == mid_fluid) ) 
              {
                mid[m] = mid_driver_face;
              }
            }
          }
        }
        break;
        
      case Z_PLUS:
        for (int k=1; k<=kmax; k++) {
          for (int j=1; j<=jmax; j++) {
            for (int i=1; i<=imax; i++) {
              m = FBUtility::getFindexS3D(size, guide, i, j, k);
              m1= FBUtility::getFindexS3D(size, guide, i, j, k-1);
              if ( (mid[m] == mid_driver) && (mid[m1] == mid_fluid) ) 
              {
                mid[m] = mid_driver_face;
              }
            }
          }
        }
        break;
    }    
  }
  
}
