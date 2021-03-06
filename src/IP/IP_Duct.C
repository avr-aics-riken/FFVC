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
 * @file IP_Duct.C
 * @brief IP_Duct class
 * @author aics
 */

#include "IP_Duct.h"



// #################################################################
// パラメータをロード
bool IP_Duct::getTP(Control* R, TextParser* tpCntl)
{
  std::string str;
  std::string label;
  REAL_TYPE ct;
  
  
  // Shape
  label="/IntrinsicExample/Shape";
  if ( !(tpCntl->getInspectedValue(label, str )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  
  if ( !strcasecmp(str.c_str(), "circular") ) {
    driver.shape = id_circular;
  }
  else if ( !strcasecmp(str.c_str(), "rectangular") ) {
    driver.shape = id_rectangular;
  }
  else {
    Hostonly_ stamped_printf("\tParsing error : Invalid shape in '%s'\n", label.c_str());
    return false;
  }
  
  
  // Diameter
  label="/IntrinsicExample/Diameter";
  if ( !(tpCntl->getInspectedValue(label, ct )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  else{
	  driver.diameter = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL;
  }
  
  
  // 媒質指定
  label="/IntrinsicExample/FluidMedium";
  if ( !(tpCntl->getInspectedValue(label, str )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  m_fluid = str;
  
  label="/IntrinsicExample/SolidMedium";
  if ( !(tpCntl->getInspectedValue(label, str )) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  m_solid = str;

  
  
  // ドライバの設定 値が正の値のとき，有効．ゼロの場合はドライバなし
  label="/IntrinsicExample/DriverLength";
  if ( tpCntl->getInspectedValue(label, ct ) ) {
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
  
  // Only driver is specified
  if ( driver.length > 0.0 )
  {
    label="/IntrinsicExample/DriverMedium";
    if ( !(tpCntl->getInspectedValue(label, str )) ) {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      return false;
    }
    m_driver = str;
    
    label="/IntrinsicExample/DriverFaceMedium";
    if ( !(tpCntl->getInspectedValue(label, str )) ) {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      return false;
    }
    m_driver_face = str;
    
    label="/IntrinsicExample/DriverDirection";
    if ( !(tpCntl->getInspectedValue(label, str )) ) {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      return false;
    }
    if ( !strcasecmp(str.c_str(), "Xminus")) {
      driver.direction = X_minus;
    }
    else if ( !strcasecmp(str.c_str(), "Xplus")) {
      driver.direction = X_plus;
    }
    else if ( !strcasecmp(str.c_str(), "Yminus")) {
      driver.direction = Y_minus;
    }
    else if ( !strcasecmp(str.c_str(), "Yplus")) {
      driver.direction = Y_plus;
    }
    else if ( !strcasecmp(str.c_str(), "Zminus")) {
      driver.direction = Z_minus;
    }
    else if ( !strcasecmp(str.c_str(), "Zplus")) {
      driver.direction = Z_plus;
    }
    else {
      Hostonly_ stamped_printf("\tParsing error : Invalid value of '%s'\n", label.c_str());
      return false;
    }
  }
  
  
  return true;
}


// #################################################################
//  パラメータの表示
void IP_Duct::printPara(FILE* fp, const Control* R)
{
  if ( !fp )
  {
    stamped_printf("\tFail to write into file\n");
    Exit(0);
  }
  
  std::string dir;
  switch (driver.direction)
  {
    case X_minus:
    case X_plus:
      dir = "X dir.";
      break;
      
    case Y_minus:
    case Y_plus:
      dir = "Y dir.";
      break;
      
    case Z_minus:
    case Z_plus:
      dir = "Z dir.";
      break;
  }
  
  fprintf(fp,"\n----------\n\n");
  fprintf(fp,"\n\t>> Intrinsic Duct Class Parameters\n\n");
  
  fprintf(fp,"\tShape of cross-section                :  %s\n", (driver.shape==id_circular)?"Circular":"Rectangular");
  fprintf(fp,"\tDiameter                  [m] / [-]   : %12.5e / %12.5e\n", driver.diameter, driver.diameter/RefL);
  fprintf(fp,"\tDirection of periodic BC              :  %s\n", dir.c_str());
  
  if ( driver.length > 0.0 )
  {
    fprintf(fp,"\twith Driver                           :  %s\n", FBUtility::getDirection(driver.direction).c_str());
    fprintf(fp,"\t     Driver Length        [m] / [-]   : %12.5e / %12.5e\n", driver.length, driver.length/RefL);
  }
  
}



// #################################################################
// Ductの計算領域のセルIDを設定する
void IP_Duct::setup(int* bcd,
                    Control* R,
                    const int NoMedium,
                    const MediumList* mat,
                    const int NoCompo,
                    const CompoList* cmp,
                    int* cutL,
                    int* cutU,
                    int* bid)
{
  int mid_fluid;        /// 流体
  int mid_solid;        /// 固体
  int mid_driver;       /// ドライバ部
  int mid_driver_face;  /// ドライバ流出面
  REAL_TYPE x, y, z, r, len;
  REAL_TYPE ox, oy, oz, Lx, Ly, Lz;
  
  // ローカルにコピー
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  // 隣接ランクのIDを取得 nID[6]
  const int* nID = paraMngr->GetNeighborRankID(procGrp);
  
  ox = origin[0];
  oy = origin[1];
  oz = origin[2];
  Lx = region[0];
  Ly = region[1];
  Lz = region[2];

  r  = driver.diameter/R->RefLength * 0.5;
  len= driver.length/R->RefLength;
  
  REAL_TYPE dx = pitch[0];
  REAL_TYPE dy = pitch[1];
  REAL_TYPE dz = pitch[2];
  
  
  // 媒質設定
  if ( (mid_fluid = FBUtility::findIDfromLabel(mat, NoMedium, m_fluid)) == 0 )
  {
    Hostonly_ printf("\tLabel '%s' is not listed in MediumList\n", m_fluid.c_str());
    Exit(0);
  }
  
  if ( (mid_solid = FBUtility::findIDfromLabel(mat, NoMedium, m_solid)) == 0 )
  {
    Hostonly_ printf("\tLabel '%s' is not listed in MediumList\n", m_solid.c_str());
    Exit(0);
  }
  
  // 矩形管の場合は内部は全て流体なので，処理の必要なし
  
  
  // 円管の場合にはsphereを参考にカット情報を入れる
  // ドライバー情報はフィルの後にしなければいけないか？　@todo
  // いろんな組み込み問題をまとめて処理する方がよい
  

  // 円管の場合，半径以下のセルを流体にする（ノードにかかわらず）
  if (driver.shape == id_circular) {
    switch (driver.direction) {
      case X_minus:
      case X_plus:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, dy, dz, mid_fluid) \
schedule(static)
        for (int k=1; k<=kx; k++) {
          for (int j=1; j<=jx; j++) {
            for (int i=1; i<=ix; i++) {
              size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
              y = oy + 0.5*dy + dy*(j-1);
              z = oz + 0.5*dz + dz*(k-1);
              if ( (y-r)*(y-r)+(z-r)*(z-r) <= r*r ) bcd[m] |= mid_fluid;
            }
          }
        }
        break;
        
      case Y_minus:
      case Y_plus:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, dx, dz, mid_fluid) \
schedule(static)
        for (int k=1; k<=kx; k++) {
          for (int j=1; j<=jx; j++) {
            for (int i=1; i<=ix; i++) {
              size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
              x = ox + 0.5*dx + dx*(i-1);
              z = oz + 0.5*dz + dz*(k-1);
              if ( (x-r)*(x-r)+(z-r)*(z-r) <= r*r ) bcd[m] |= mid_fluid;
            }
          }
        }
        break;
        
      case Z_minus:
      case Z_plus:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, dx, dy, mid_fluid) \
schedule(static)
        for (int k=1; k<=kx; k++) {
          for (int j=1; j<=jx; j++) {
            for (int i=1; i<=ix; i++) {
              size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
              x = ox + 0.5*dx + dx*(i-1);
              y = oy + 0.5*dy + dy*(j-1);
              if ( (x-r)*(x-r)+(y-r)*(y-r) <= r*r ) bcd[m] |= mid_fluid;
            }
          }
        }
        break;
    }   
  }
  
  // ドライバ部分
  if ( driver.length > 0.0 ) {
    
    switch (driver.direction) {
      case X_minus:
        if ( nID[driver.direction] < 0 ) {
          for (int k=1; k<=kx; k++) {
            for (int j=1; j<=jx; j++) {
              for (int i=1; i<=ix; i++) {
                size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
                x = ox + 0.5*dx + dx*(i-1);
                y = oy + 0.5*dy + dy*(j-1);
                z = oz + 0.5*dz + dz*(k-1);
                if ( x < ox+len ) {
                  if ( driver.shape == id_circular ) {
                    if ( (y-r)*(y-r)+(z-r)*(z-r) <= r*r ) bcd[m] |= mid_driver; // 半径以内をドライバIDにする
                  }
                  else {
                    bcd[m] |= mid_driver;
                  }
                }
              }
            }
          }
        }        
        break;
        
      case X_plus:
        if ( nID[driver.direction] < 0 ) {
          for (int k=1; k<=kx; k++) {
            for (int j=1; j<=jx; j++) {
              for (int i=1; i<=ix; i++) {
                size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
                x = ox + 0.5*dx + dx*(i-1);
                y = oy + 0.5*dy + dy*(j-1);
                z = oz + 0.5*dz + dz*(k-1);
                if ( x > ox+Lx-len ) {
                  if ( driver.shape == id_circular ) {
                    if ( (y-r)*(y-r)+(z-r)*(z-r) <= r*r ) bcd[m] |= mid_driver;
                  }
                  else {
                    bcd[m] |= mid_driver;
                  }
                }                
              }
            }
          }
        }
        break;
        
      case Y_minus:
        if ( nID[driver.direction] < 0 ) {
          for (int k=1; k<=kx; k++) {
            for (int j=1; j<=jx; j++) {
              for (int i=1; i<=ix; i++) {
                size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
                x = ox + 0.5*dx + dx*(i-1);
                y = oy + 0.5*dy + dy*(j-1);
                z = oz + 0.5*dz + dz*(k-1);
                if ( y < oy+len ) {
                  if ( driver.shape == id_circular ) {
                    if ( (x-r)*(x-r)+(z-r)*(z-r) <= r*r ) bcd[m] |= mid_driver;
                  }
                  else {
                    bcd[m] |= mid_driver;
                  }
                }
              }
            }
          }
        }        
        break;
        
      case Y_plus:
        if ( nID[driver.direction] < 0 ) {
          for (int k=1; k<=kx; k++) {
            for (int j=1; j<=jx; j++) {
              for (int i=1; i<=ix; i++) {
                size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
                x = ox + 0.5*dx + dx*(i-1);
                y = oy + 0.5*dy + dy*(j-1);
                z = oz + 0.5*dz + dz*(k-1);
                if ( y > oy+Ly-len ) {
                  if ( driver.shape == id_circular ) {
                    if ( (x-r)*(x-r)+(z-r)*(z-r) <= r*r ) bcd[m] |= mid_driver;
                  }
                  else {
                    bcd[m] |= mid_driver;
                  }
                }
              }
            }
          }
        }
        break;
        
      case Z_minus:
        if ( nID[driver.direction] < 0 ) {
          for (int k=1; k<=kx; k++) {
            for (int j=1; j<=jx; j++) {
              for (int i=1; i<=ix; i++) {
                size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
                x = ox + 0.5*dx + dx*(i-1);
                y = oy + 0.5*dy + dy*(j-1);
                z = oz + 0.5*dz + dz*(k-1);
                if ( z < oz+len ) {
                  if ( driver.shape == id_circular ) {
                    if ( (x-r)*(x-r)+(y-r)*(y-r) <= r*r ) bcd[m] |= mid_driver;
                  }
                  else {
                    bcd[m] |= mid_driver;
                  }
                }
              }
            }
          }
        }
        break;
        
      case Z_plus:
        if ( nID[driver.direction] < 0 ) {
          for (int k=1; k<=kx; k++) {
            for (int j=1; j<=jx; j++) {
              for (int i=1; i<=ix; i++) {
                size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
                x = ox + 0.5*dx + dx*(i-1);
                y = oy + 0.5*dy + dy*(j-1);
                z = oz + 0.5*dz + dz*(k-1);
                if ( z > oz+Lz-len ) {
                  if ( driver.shape == id_circular ) {
                    if ( (x-r)*(x-r)+(y-r)*(y-r) <= r*r ) bcd[m] |= mid_driver;
                  }
                  else {
                    bcd[m] |= mid_driver;
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
      case X_minus:
        for (int k=1; k<=kx; k++) {
          for (int j=1; j<=jx; j++) {
            for (int i=1; i<=ix; i++) {
              size_t m = _F_IDX_S3D(i,   j, k, ix, jx, kx, gd);
              size_t m1= _F_IDX_S3D(i+1, j, k, ix, jx, kx, gd);
              if ( (DECODE_CMP( bcd[m] )  == mid_driver) && (DECODE_CMP( bcd[m1] ) == mid_fluid) )
              {
                bcd[m] |= mid_driver_face;
              }
            }
          }
        }        
        break;
        
      case X_plus:
        for (int k=1; k<=kx; k++) {
          for (int j=1; j<=jx; j++) {
            for (int i=1; i<=ix; i++) {
              size_t m = _F_IDX_S3D(i,   j, k, ix, jx, kx, gd);
              size_t m1= _F_IDX_S3D(i-1, j, k, ix, jx, kx, gd);
              if ( (DECODE_CMP( bcd[m] ) == mid_driver) && (DECODE_CMP( bcd[m1] )  == mid_fluid) ) 
              {
                bcd[m] |= mid_driver_face;
              }
            }
          }
        }
        break;
        
      case Y_minus:
        for (int k=1; k<=kx; k++) {
          for (int j=1; j<=jx; j++) {
            for (int i=1; i<=ix; i++) {
              size_t m = _F_IDX_S3D(i, j,   k, ix, jx, kx, gd);
              size_t m1= _F_IDX_S3D(i, j+1, k, ix, jx, kx, gd);
              if ( (DECODE_CMP( bcd[m] ) == mid_driver) && (DECODE_CMP( bcd[m1] )  == mid_fluid) ) 
              {
                bcd[m] |= mid_driver_face;
              }
            }
          }
        }        
        break;
        
      case Y_plus:
        for (int k=1; k<=kx; k++) {
          for (int j=1; j<=jx; j++) {
            for (int i=1; i<=ix; i++) {
              size_t m = _F_IDX_S3D(i, j,   k, ix, jx, kx, gd);
              size_t m1= _F_IDX_S3D(i, j-1, k, ix, jx, kx, gd);
              if ( (DECODE_CMP( bcd[m] ) == mid_driver) && (DECODE_CMP( bcd[m1] )  == mid_fluid) ) 
              {
                bcd[m] |= mid_driver_face;
              }
            }
          }
        }
        break;
        
      case Z_minus:
        for (int k=1; k<=kx; k++) {
          for (int j=1; j<=jx; j++) {
            for (int i=1; i<=ix; i++) {
              size_t m = _F_IDX_S3D(i, j, k,   ix, jx, kx, gd);
              size_t m1= _F_IDX_S3D(i, j, k+1, ix, jx, kx, gd);
              if ( (DECODE_CMP( bcd[m] ) == mid_driver) && (DECODE_CMP( bcd[m1] )  == mid_fluid) ) 
              {
                bcd[m] |= mid_driver_face;
              }
            }
          }
        }
        break;
        
      case Z_plus:
        for (int k=1; k<=kx; k++) {
          for (int j=1; j<=jx; j++) {
            for (int i=1; i<=ix; i++) {
              size_t m = _F_IDX_S3D(i, j, k,   ix, jx, kx, gd);
              size_t m1= _F_IDX_S3D(i, j, k-1, ix, jx, kx, gd);
              if ( (DECODE_CMP( bcd[m] ) == mid_driver) && (DECODE_CMP( bcd[m1] )  == mid_fluid) ) 
              {
                bcd[m] |= mid_driver_face;
              }
            }
          }
        }
        break;
    }    
  }
  
}
