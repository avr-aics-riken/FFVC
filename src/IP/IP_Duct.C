/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2011
 *
 */

//@file IP_Duct.C
//@brief IP_Duct class
//@author keno, FSI Team, VCAD, RIKEN

#include "IP_Duct.h"

/**
 @fn bool IP_Duct::getXML(SklSolverConfig* CF, Control* R)
 @brief Ductのパラメータを取得する
 @param CF コンフィギュレーションツリー
 @param R Controlクラスのポインタ
 */
bool IP_Duct::getXML(SklSolverConfig* CF, Control* R)
{
  const CfgElem *elemTop=NULL, *elmL1=NULL;
  const char *str=NULL;
  SKL_REAL ct=0.0;
  
  if ( !(elemTop = CF->GetTop(PARAMETER)) ) return false;
  
  if( !(elmL1 = elemTop->GetElemFirst("Intrinsic_Example")) ) {
    Hostonly_ stamped_printf("\tParsing error : Missing the section of 'Intrinsic_Example'\n");
    return false;
  }
  
  // 断面形状
  if ( !elmL1->GetValue("Shape", &str) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Shape' in 'Intrinsic_Example'\n");
    return false;
  }
  if ( !strcasecmp(str, "circular") ) {
    driver.shape = id_circular;
  }
  else if ( !strcasecmp(str, "rectangualr") ) {
    driver.shape = id_rectangular;
  }
  else {
    Hostonly_ stamped_printf("\tParsing error : Invalid shape in 'Intrinsic_Example'\n");
    return false;
  }
  
  // 直径
  if ( elmL1->GetValue(CfgIdt("Diameter"), &ct) ) {
    driver.diameter = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL;
  }
  else {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Diameter' in 'intrinsic_Example'\n");
    return false;
  }
  
  // 周期境界の方向
  if ( elmL1->GetValue("Direction", &str) ) {
    if ( !strcasecmp(str, "X_minus")) {
      driver.direction = X_MINUS;
    }
    else if ( !strcasecmp(str, "X_plus")) {
      driver.direction = X_PLUS;
    }
    else if ( !strcasecmp(str, "Y_minus")) {
      driver.direction = Y_MINUS;
    }
    else if ( !strcasecmp(str, "Y_plus")) {
      driver.direction = Y_PLUS;
    }
    else if ( !strcasecmp(str, "Z_minus")) {
      driver.direction = Z_MINUS;
    }
    else if ( !strcasecmp(str, "Z_plus")) {
      driver.direction = Z_PLUS;
    }
    else {
      Hostonly_ stamped_printf("\tParsing error : Invalid value of 'Direction' in 'Intrinsic_Example'\n");
      return false;
    }      
  }
  else {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Direction' in 'Intrinsic_Example'\n");
    return false;
  }
  
  // ドライバの設定 値が正の値のとき，有効．ゼロの場合はドライバなし
  if ( elmL1->GetValue(CfgIdt("Driver"), &ct) ) {
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
  
  return true;
}

/**
 @fn void IP_Duct::printPara(FILE* fp, Control* R)
 @brief パラメータの表示
 @param fp ファイルポインタ
 @param R コントロールクラスのポインタ
 */
void IP_Duct::printPara(FILE* fp, Control* R)
{
  if ( !fp ) {
    stamped_printf("\tFail to write into file\n");
    assert(0);
  }
  
  string dir;
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

/**
 @fn bool IP_Duct::setDomain(Control* R, unsigned sz[3], SKL_REAL org[3], SKL_REAL wth[3], SKL_REAL pch[3])
 @brief Ductの領域情報を設定する
 @param R Controlクラスのポインタ
 @param sz グローバル計算領域のセルサイズ
 @param org グローバル計算領域の基点
 @param wth グローバル計算領域のbounding boxサイズ
 @param pch セルピッチ
 */
void IP_Duct::setDomain(Control* R, unsigned sz[3], SKL_REAL org[3], SKL_REAL wth[3], SKL_REAL pch[3])
{
  wth[0] = pch[0]*(SKL_REAL)sz[0];
  wth[1] = pch[1]*(SKL_REAL)sz[1];
  wth[2] = pch[2]*(SKL_REAL)sz[2];
  
  // チェック
  if ( (pch[0] != pch[1]) || (pch[1] != pch[2]) ) {
    Hostonly_ printf("Error : 'VoxelPitch' in each direction must be same.\n");
    assert(0);
  }
  if ( ((unsigned)(wth[0]/pch[0]) != sz[0]) ||
       ((unsigned)(wth[1]/pch[1]) != sz[1]) ||
       ((unsigned)(wth[2]/pch[2]) != sz[2]) ) {
    Hostonly_ printf("Error : Invalid parameters among 'VoxelSize', 'VoxelPitch', and 'VoxelWidth' in DomainInfo section.\n");
    assert(0);
  }

}

/**
 @fn void IP_Duct::setup(int* mid, Control* R, SKL_REAL* G_org)
 @brief Ductの計算領域のセルIDを設定する
 @param mid IDの配列
 @param R Controlクラスのポインタ
 @param G_org グローバルな原点（無次元）
 */
void IP_Duct::setup(int* mid, Control* R, SKL_REAL* G_org)
{
  int i,j,k, gd;
  int mid_fluid=1;        /// 流体
  int mid_solid=2;        /// 固体
  int mid_driver=3;       /// ドライバ部
  int mid_driver_face=4;  /// ドライバ流出面
  unsigned m;
  SKL_REAL x, y, z, dh, r, len;
  SKL_REAL ox, oy, oz, Lx, Ly, Lz;
  
  ox = R->org[0];
  oy = R->org[1];
  oz = R->org[2];
  Lx = R->Lbx[0];
  Ly = R->Lbx[1];
  Lz = R->Lbx[2];
  dh = R->dh;
  gd = (int)guide;
  r  = driver.diameter/R->RefLength * 0.5;
  len= driver.length/R->RefLength;
  
  // Initialize  内部領域をsolidにしておく
  for (k=1; k<=(int)kmax; k++) { 
    for (j=1; j<=(int)jmax; j++) {
      for (i=1; i<=(int)imax; i++) {
        m = SklUtil::getFindexS3D(size, guide, i, j, k);
        mid[m] = mid_solid;
      }
    }
  }
  
  // Inner
  if (driver.shape == id_rectangular ) { // 矩形管の場合は内部は全て流体
    for (k=1; k<=(int)kmax; k++) {
      for (j=1; j<=(int)jmax; j++) {
        for (i=1; i<=(int)imax; i++) {
          m = SklUtil::getFindexS3D(size, guide, i, j, k);
          mid[m] = mid_fluid;
        }
      }
    }
  }
  else { // 円管の場合，半径以下のセルを流体にする（ノードにかかわらず）
    switch (driver.direction) {
      case X_MINUS:
      case X_PLUS:
        for (k=1; k<=(int)kmax; k++) {
          for (j=1; j<=(int)jmax; j++) {
            for (i=1; i<=(int)imax; i++) {
              m = SklUtil::getFindexS3D(size, guide, i, j, k);
              y = oy + 0.5*dh + dh*(j-1);
              z = oz + 0.5*dh + dh*(k-1);
              if ( (y-r)*(y-r)+(z-r)*(z-r) <= r*r ) mid[m] = mid_fluid;
            }
          }
        }
        break;
        
      case Y_MINUS:
      case Y_PLUS:
        for (k=1; k<=(int)kmax; k++) {
          for (j=1; j<=(int)jmax; j++) {
            for (i=1; i<=(int)imax; i++) {
              m = SklUtil::getFindexS3D(size, guide, i, j, k);
              x = ox + 0.5*dh + dh*(i-1);
              z = oz + 0.5*dh + dh*(k-1);
              if ( (x-r)*(x-r)+(z-r)*(z-r) <= r*r ) mid[m] = mid_fluid;
            }
          }
        }
        break;
        
      case Z_MINUS:
      case Z_PLUS:
        for (k=1; k<=(int)kmax; k++) {
          for (j=1; j<=(int)jmax; j++) {
            for (i=1; i<=(int)imax; i++) {
              m = SklUtil::getFindexS3D(size, guide, i, j, k);
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
        if ( pn.nID[driver.direction] < 0 ) {
          for (k=1; k<=(int)kmax; k++) {
            for (j=1; j<=(int)jmax; j++) {
              for (i=1; i<=(int)imax; i++) {
                m = SklUtil::getFindexS3D(size, guide, i, j, k);
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
        if ( pn.nID[driver.direction] < 0 ) {
          for (k=1; k<=(int)kmax; k++) {
            for (j=1; j<=(int)jmax; j++) {
              for (i=1; i<=(int)imax; i++) {
                m = SklUtil::getFindexS3D(size, guide, i, j, k);
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
        if ( pn.nID[driver.direction] < 0 ) {
          for (k=1; k<=(int)kmax; k++) {
            for (j=1; j<=(int)jmax; j++) {
              for (i=1; i<=(int)imax; i++) {
                m = SklUtil::getFindexS3D(size, guide, i, j, k);
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
        if ( pn.nID[driver.direction] < 0 ) {
          for (k=1; k<=(int)kmax; k++) {
            for (j=1; j<=(int)jmax; j++) {
              for (i=1; i<=(int)imax; i++) {
                m = SklUtil::getFindexS3D(size, guide, i, j, k);
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
        if ( pn.nID[driver.direction] < 0 ) {
          for (k=1; k<=(int)kmax; k++) {
            for (j=1; j<=(int)jmax; j++) {
              for (i=1; i<=(int)imax; i++) {
                m = SklUtil::getFindexS3D(size, guide, i, j, k);
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
        if ( pn.nID[driver.direction] < 0 ) {
          for (k=1; k<=(int)kmax; k++) {
            for (j=1; j<=(int)jmax; j++) {
              for (i=1; i<=(int)imax; i++) {
                m = SklUtil::getFindexS3D(size, guide, i, j, k);
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
    
    unsigned m1;
    switch (driver.direction) {
      case X_MINUS:
        for (k=1; k<=(int)kmax; k++) {
          for (j=1; j<=(int)jmax; j++) {
            for (i=1; i<=(int)imax; i++) {
              m = SklUtil::getFindexS3D(size, guide, i,   j, k);
              m1= SklUtil::getFindexS3D(size, guide, i+1, j, k);
              if ( (mid[m] == mid_driver) && (mid[m1] == mid_fluid) ) {
                mid[m] = mid_driver_face;
              }
            }
          }
        }        
        break;
        
      case X_PLUS:
        for (k=1; k<=(int)kmax; k++) {
          for (j=1; j<=(int)jmax; j++) {
            for (i=1; i<=(int)imax; i++) {
              m = SklUtil::getFindexS3D(size, guide, i,   j, k);
              m1= SklUtil::getFindexS3D(size, guide, i-1, j, k);
              if ( (mid[m] == mid_driver) && (mid[m1] == mid_fluid) ) {
                mid[m] = mid_driver_face;
              }
            }
          }
        }
        break;
        
      case Y_MINUS:
        for (k=1; k<=(int)kmax; k++) {
          for (j=1; j<=(int)jmax; j++) {
            for (i=1; i<=(int)imax; i++) {
              m = SklUtil::getFindexS3D(size, guide, i, j,   k);
              m1= SklUtil::getFindexS3D(size, guide, i, j+1, k);
              if ( (mid[m] == mid_driver) && (mid[m1] == mid_fluid) ) {
                mid[m] = mid_driver_face;
              }
            }
          }
        }        
        break;
        
      case Y_PLUS:
        for (k=1; k<=(int)kmax; k++) {
          for (j=1; j<=(int)jmax; j++) {
            for (i=1; i<=(int)imax; i++) {
              m = SklUtil::getFindexS3D(size, guide, i, j,   k);
              m1= SklUtil::getFindexS3D(size, guide, i, j-1, k);
              if ( (mid[m] == mid_driver) && (mid[m1] == mid_fluid) ) {
                mid[m] = mid_driver_face;
              }
            }
          }
        }
        break;
        
      case Z_MINUS:
        for (k=1; k<=(int)kmax; k++) {
          for (j=1; j<=(int)jmax; j++) {
            for (i=1; i<=(int)imax; i++) {
              m = SklUtil::getFindexS3D(size, guide, i, j, k);
              m1= SklUtil::getFindexS3D(size, guide, i, j, k+1);
              if ( (mid[m] == mid_driver) && (mid[m1] == mid_fluid) ) {
                mid[m] = mid_driver_face;
              }
            }
          }
        }
        break;
        
      case Z_PLUS:
        for (k=1; k<=(int)kmax; k++) {
          for (j=1; j<=(int)jmax; j++) {
            for (i=1; i<=(int)imax; i++) {
              m = SklUtil::getFindexS3D(size, guide, i, j, k);
              m1= SklUtil::getFindexS3D(size, guide, i, j, k-1);
              if ( (mid[m] == mid_driver) && (mid[m1] == mid_fluid) ) {
                mid[m] = mid_driver_face;
              }
            }
          }
        }
        break;
    }    
  }
  
}