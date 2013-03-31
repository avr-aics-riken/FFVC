// #################################################################
//
// FFV : Frontflow / violet Cartesian
//
// Copyright (c) 2012-2013  All right reserved.
//
// Institute of Industrial Science, University of Tokyo, Japan. 
//
// #################################################################

/** 
 * @file IP_Duct.C
 * @brief IP_Duct class
 * @author kero
 */

#include "IP_Duct.h"



// #################################################################
/*
 * @brief パラメータをロード
 * @param [in] R      Controlクラス
 * @param [in] tpCntl テキストパーサクラス
 * @return true-成功, false-エラー
 */
bool IP_Duct::getTP(Control* R, TPControl* tpCntl)
{
  std::string str;
  std::string label;
  REAL_TYPE ct;
  
  // Shape
  label="/Parameter/IntrinsicExample/Shape";
  if ( !(tpCntl->GetValue(label, &str )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  
  if ( !strcasecmp(str.c_str(), "circular") ) {
    driver.shape = id_circular;
  }
  else if ( !strcasecmp(str.c_str(), "rectangualr") ) {
    driver.shape = id_rectangular;
  }
  else {
    Hostonly_ stamped_printf("\tParsing error : Invalid shape in '%s'\n", label.c_str());
    return false;
  }
  
  // Diameter
  label="/Parameter/IntrinsicExample/Diameter";
  if ( !(tpCntl->GetValue(label, &ct )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  else{
	  driver.diameter = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL;
  }
  
  // periodic
  label="/Parameter/IntrinsicExample/Direction";
  if ( !(tpCntl->GetValue(label, &str )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  if ( !strcasecmp(str.c_str(), "Xminus")) {
    driver.direction = X_MINUS;
  }
  else if ( !strcasecmp(str.c_str(), "Xplus")) {
    driver.direction = X_PLUS;
  }
  else if ( !strcasecmp(str.c_str(), "Yminus")) {
    driver.direction = Y_MINUS;
  }
  else if ( !strcasecmp(str.c_str(), "Yplus")) {
    driver.direction = Y_PLUS;
  }
  else if ( !strcasecmp(str.c_str(), "Zminus")) {
    driver.direction = Z_MINUS;
  }
  else if ( !strcasecmp(str.c_str(), "Zplus")) {
    driver.direction = Z_PLUS;
  }
  else {
    Hostonly_ stamped_printf("\tParsing error : Invalid value of '%s'\n", label.c_str());
    return false;
  }     
  
  // ドライバの設定 値が正の値のとき，有効．ゼロの場合はドライバなし
  label="/Parameter/IntrinsicExample/Driver";
  if ( tpCntl->GetValue(label, &ct ) ) {
    driver.length = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL;
  }
  else {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  
  if ( driver.length < 0.0 ) {
    Hostonly_ stamped_printf("\tError : Value of 'Driver' must be positive.\n");
    return false;
  }
  
  // 媒質指定
  label="/Parameter/IntrinsicExample/FluidMedium";
  if ( !(tpCntl->GetValue(label, &str )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  m_fluid = str;
  
  label="/Parameter/IntrinsicExample/SolidMedium";
  if ( !(tpCntl->GetValue(label, &str )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  m_solid = str;
  
  label="/Parameter/IntrinsicExample/DriverMedium";
  if ( !(tpCntl->GetValue(label, &str )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  m_driver = str;
  
  label="/Parameter/IntrinsicExample/DriverFaceMedium";
  if ( !(tpCntl->GetValue(label, &str )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  m_driver_face = str;
  
  return true;
}


// #################################################################
/**
 * @brief パラメータの表示
 * @param [in] fp ファイルポインタ
 * @param [in] R  コントロールクラスのポインタ
 */
void IP_Duct::printPara(FILE* fp, const Control* R)
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
  fprintf(fp,"\n\t>> Intrinsic Duct Class Parameters\n\n");
  
  fprintf(fp,"\tShape of cross-section                :  %s\n", (driver.shape==id_circular)?"Circular":"Rectangular");
  fprintf(fp,"\tDiameter                  [m] / [-]   : %12.5e / %12.5e\n", driver.diameter, driver.diameter/RefL);
  fprintf(fp,"\tDirection of periodic BC              :  %s\n", dir.c_str());
  if ( driver.length > 0.0 ) {
    fprintf(fp,"\twith Driver                           :  %s\n", FBUtility::getDirection(driver.direction).c_str());
    fprintf(fp,"\t     Driver Length        [m] / [-]   : %12.5e / %12.5e\n", driver.length, driver.length/RefL);
  }
  
}


// #################################################################
/* @brief 領域パラメータを設定する
 * @param [in]     R   Controlクラスのポインタ
 * @param [in]     sz  分割数
 * @param [in,out] org 計算領域の基点
 * @param [in,out] reg 計算領域のbounding boxサイズ
 * @param [in,out] pch セル幅
 */
void IP_Duct::setDomainParameter(Control* R, const int* sz, REAL_TYPE* org, REAL_TYPE* reg, REAL_TYPE* pch)
{
  RefL = R->RefLength;
  
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


// #################################################################
/*
 * @brief Ductの計算領域のセルIDを設定する
 * @param [in,out] mid   媒質情報の配列
 * @param [in]     R     Controlクラスのポインタ
 * @param [in]     G_org グローバルな原点（無次元）
 * @param [in]     Nmax  Controlクラスのポインタ
 * @param [in]     mat   MediumListクラスのポインタ
 */
void IP_Duct::setup(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, MediumList* mat)
{
  int mid_fluid=1;        /// 流体
  int mid_solid=2;        /// 固体
  int mid_driver=3;       /// ドライバ部
  int mid_driver_face=4;  /// ドライバ流出面
  REAL_TYPE x, y, z, dh, r, len;
  REAL_TYPE ox, oy, oz, Lx, Ly, Lz;
  
  // ローカルにコピー
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  // 隣接ランクのIDを取得 nID[6]
  const int* nID = paraMngr->GetNeighborRankID();
  
  ox = origin[0];
  oy = origin[1];
  oz = origin[2];
  Lx = region[0];
  Ly = region[1];
  Lz = region[2];
  dh = deltaX;
  r  = driver.diameter/R->RefLength * 0.5;
  len= driver.length/R->RefLength;
  
  // Initialize  内部領域をsolidにしておく
#pragma omp parallel for firstprivate(ix, jx, kx, gd, mid_solid) schedule(static)
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        mid[m] = mid_solid;
      }
    }
  }
  
  // Inner
  if (driver.shape == id_rectangular ) { // 矩形管の場合は内部は全て流体
#pragma omp parallel for firstprivate(ix, jx, kx, gd, mid_fluid) \
schedule(static)
    for (int k=1; k<=kx; k++) {
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
          mid[m] = mid_fluid;
        }
      }
    }
  }
  else { // 円管の場合，半径以下のセルを流体にする（ノードにかかわらず）
    switch (driver.direction) {
      case X_MINUS:
      case X_PLUS:
        for (int k=1; k<=kx; k++) {
          for (int j=1; j<=jx; j++) {
            for (int i=1; i<=ix; i++) {
              size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
              y = oy + 0.5*dh + dh*(j-1);
              z = oz + 0.5*dh + dh*(k-1);
              if ( (y-r)*(y-r)+(z-r)*(z-r) <= r*r ) mid[m] = mid_fluid;
            }
          }
        }
        break;
        
      case Y_MINUS:
      case Y_PLUS:
        for (int k=1; k<=kx; k++) {
          for (int j=1; j<=jx; j++) {
            for (int i=1; i<=ix; i++) {
              size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
              x = ox + 0.5*dh + dh*(i-1);
              z = oz + 0.5*dh + dh*(k-1);
              if ( (x-r)*(x-r)+(z-r)*(z-r) <= r*r ) mid[m] = mid_fluid;
            }
          }
        }
        break;
        
      case Z_MINUS:
      case Z_PLUS:
        for (int k=1; k<=kx; k++) {
          for (int j=1; j<=jx; j++) {
            for (int i=1; i<=ix; i++) {
              size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
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
          for (int k=1; k<=kx; k++) {
            for (int j=1; j<=jx; j++) {
              for (int i=1; i<=ix; i++) {
                size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
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
          for (int k=1; k<=kx; k++) {
            for (int j=1; j<=jx; j++) {
              for (int i=1; i<=ix; i++) {
                size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
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
          for (int k=1; k<=kx; k++) {
            for (int j=1; j<=jx; j++) {
              for (int i=1; i<=ix; i++) {
                size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
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
          for (int k=1; k<=kx; k++) {
            for (int j=1; j<=jx; j++) {
              for (int i=1; i<=ix; i++) {
                size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
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
          for (int k=1; k<=kx; k++) {
            for (int j=1; j<=jx; j++) {
              for (int i=1; i<=ix; i++) {
                size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
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
          for (int k=1; k<=kx; k++) {
            for (int j=1; j<=jx; j++) {
              for (int i=1; i<=ix; i++) {
                size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
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
    
    switch (driver.direction) 
    {
      case X_MINUS:
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
        break;
        
      case X_PLUS:
        for (int k=1; k<=kx; k++) {
          for (int j=1; j<=jx; j++) {
            for (int i=1; i<=ix; i++) {
              size_t m = _F_IDX_S3D(i,   j, k, ix, jx, kx, gd);
              size_t m1= _F_IDX_S3D(i-1, j, k, ix, jx, kx, gd);
              if ( (mid[m] == mid_driver) && (mid[m1] == mid_fluid) ) 
              {
                mid[m] = mid_driver_face;
              }
            }
          }
        }
        break;
        
      case Y_MINUS:
        for (int k=1; k<=kx; k++) {
          for (int j=1; j<=jx; j++) {
            for (int i=1; i<=ix; i++) {
              size_t m = _F_IDX_S3D(i, j,   k, ix, jx, kx, gd);
              size_t m1= _F_IDX_S3D(i, j+1, k, ix, jx, kx, gd);
              if ( (mid[m] == mid_driver) && (mid[m1] == mid_fluid) ) 
              {
                mid[m] = mid_driver_face;
              }
            }
          }
        }        
        break;
        
      case Y_PLUS:
        for (int k=1; k<=kx; k++) {
          for (int j=1; j<=jx; j++) {
            for (int i=1; i<=ix; i++) {
              size_t m = _F_IDX_S3D(i, j,   k, ix, jx, kx, gd);
              size_t m1= _F_IDX_S3D(i, j-1, k, ix, jx, kx, gd);
              if ( (mid[m] == mid_driver) && (mid[m1] == mid_fluid) ) 
              {
                mid[m] = mid_driver_face;
              }
            }
          }
        }
        break;
        
      case Z_MINUS:
        for (int k=1; k<=kx; k++) {
          for (int j=1; j<=jx; j++) {
            for (int i=1; i<=ix; i++) {
              size_t m = _F_IDX_S3D(i, j, k,   ix, jx, kx, gd);
              size_t m1= _F_IDX_S3D(i, j, k+1, ix, jx, kx, gd);
              if ( (mid[m] == mid_driver) && (mid[m1] == mid_fluid) ) 
              {
                mid[m] = mid_driver_face;
              }
            }
          }
        }
        break;
        
      case Z_PLUS:
        for (int k=1; k<=kx; k++) {
          for (int j=1; j<=jx; j++) {
            for (int i=1; i<=ix; i++) {
              size_t m = _F_IDX_S3D(i, j, k,   ix, jx, kx, gd);
              size_t m1= _F_IDX_S3D(i, j, k-1, ix, jx, kx, gd);
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
