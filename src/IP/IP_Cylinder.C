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
// Copyright (c) 2012-2015 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

/** 
 * @file   IP_Cylinder.C
 * @brief  IP_Cylinder class
 * @author aics
 */

#include "IP_Cylinder.h"


using namespace Vec3class;


// #################################################################
// XY面内の交点の無次元距離を計算する
REAL_TYPE IP_Cylinder::cut_line_2d(const Vec3r p, const int dir, const REAL_TYPE r, const REAL_TYPE px)
{
  REAL_TYPE x, y, s, c;
  
  x = p.x;
  y = p.y;
  
  s = 0.0;
  
  // 基点座標の符号で交点座標を判断
  switch (dir)
  {
    case 1: // X-
    case 2: // X+
      c = sqrtf(r*r - y*y);
      if ( x < 0.0 ) c *= -1.0;
      s = fabs(c-x);
      break;
      
    case 3: // Y-
    case 4: // Y+
      c = sqrtf(r*r - x*x);
      if ( y < 0.0 ) c *= -1.0;
      s = fabs(c-y);
      break;
      
    default:
      Exit(0);
      break;
  }
  
  return s/px;
}


// #################################################################
/**
 * @brief 点pの属するセルインデクスを求める
 * @param [in] p   探索座標
 * @param [in] ol  基点座標
 * @param [in] pch 格子幅
 * @return cell index
 * @note Fortran index
 */
Vec3i IP_Cylinder::find_index(const Vec3r p, const Vec3r ol, const Vec3r pch)
{
  Vec3r q = (p-ol)/pch;
  Vec3i idx( ceil(q.x), ceil(q.y), ceil(q.z) );
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  
  if ( idx.x < 1 ) idx.x = 1;
  if ( idx.y < 1 ) idx.y = 1;
  if ( idx.z < 1 ) idx.z = 1;
  if ( idx.x > ix ) idx.x = ix;
  if ( idx.y > jx ) idx.y = jx;
  if ( idx.z > kx ) idx.z = kx;
  
  return idx;
}


// #################################################################
// パラメータをロード
bool IP_Cylinder::getTP(Control* R, TextParser* tpCntl)
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
  
  
  // Cylinder 1
  label = "/IntrinsicExample/Cylinder1/UseCylinder";
  
  if ( !(tpCntl->getInspectedValue(label, str)) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  if     ( !strcasecmp(str.c_str(), "yes") )
  {
    cyl1.ens = ON;
    num_cyls++;
  }
  else if( !strcasecmp(str.c_str(), "no") )
  {
    cyl1.ens = OFF;
  }
  else {
    Hostonly_ stamped_printf("\tParsing error : Invalid '%s'\n", label.c_str());
    return false;
  }
  
  
  if ( cyl1.ens == ON )
  {
    // shape
    label = "/IntrinsicExample/Cylinder1/Shape";
    
    if ( !(tpCntl->getInspectedValue(label, str )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      return false;
    }
    
    if     ( !strcasecmp(str.c_str(), "circular") )
    {
      cyl1.shape = id_circular;
    }
    else if( !strcasecmp(str.c_str(), "rectangular") )
    {
      cyl1.shape = id_rectangular;
    }
    else
    {
      Hostonly_ stamped_printf("\tParsing error : Invalid '%s'\n", label.c_str());
      return false;
    }
    
    if ( cyl1.shape == id_rectangular)
    {
      label = "/IntrinsicExample/Cylinder1/LengthX";
      if ( !(tpCntl->getInspectedValue(label, ct )) )
      {
        Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
        return false;
      }
      else
      {
        cyl1.length_x = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL;
      }
      
      
      label="/IntrinsicExample/Cylinder1/LengthY";
      if ( !(tpCntl->getInspectedValue(label, ct )) )
      {
        Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
        return false;
      }
      else
      {
        cyl1.length_y = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL;
      }
    }
    else
    {
      // radius
      label = "/IntrinsicExample/Cylinder1/Radius";
      
      if ( !(tpCntl->getInspectedValue(label, ct )) )
      {
        Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
        return false;
      }
      else
      {
        cyl1.radius = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL;
      }
    }
    
    if ( mode == dim_3d )
    {
      label="/IntrinsicExample/Cylinder1/LengthZ";
      if ( !(tpCntl->getInspectedValue(label, ct )) )
      {
        Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
        return false;
      }
      else
      {
        cyl1.length_z = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL;
      }
    }
    
    // position
    label = "/IntrinsicExample/Cylinder1/PositionX";
    
    if ( !(tpCntl->getInspectedValue(label, ct )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      return false;
    }
    else
    {
      cyl1.x = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL;
    }
    
    label = "/IntrinsicExample/Cylinder1/PositionY";
    
    if ( !(tpCntl->getInspectedValue(label, ct )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      return false;
    }
    else
    {
      cyl1.y = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL;
    }
    
    // 固体媒質
    label = "/IntrinsicExample/Cylinder1/SolidMedium";
    
    if ( !(tpCntl->getInspectedValue(label, str )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      return false;
    }
    cyl1.solid = str;
  }
  

  // Cylinder 2
  label = "/IntrinsicExample/Cylinder2/UseCylinder";
  
  if ( !(tpCntl->getInspectedValue(label, str)) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  if     ( !strcasecmp(str.c_str(), "yes") )
  {
    cyl2.ens = ON;
    num_cyls++;
  }
  else if( !strcasecmp(str.c_str(), "no") )
  {
    cyl2.ens = OFF;
  }
  else {
    Hostonly_ stamped_printf("\tParsing error : Invalid '%s'\n", label.c_str());
    return false;
  }
  
  
  if ( cyl2.ens == ON )
  {
    // shape
    label = "/IntrinsicExample/Cylinder2/Shape";
    
    if ( !(tpCntl->getInspectedValue(label, str )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      return false;
    }
    
    if     ( !strcasecmp(str.c_str(), "circular") )
    {
      cyl2.shape = id_circular;
    }
    else if( !strcasecmp(str.c_str(), "rectangular") )
    {
      cyl2.shape = id_rectangular;
    }
    else
    {
      Hostonly_ stamped_printf("\tParsing error : Invalid '%s'\n", label.c_str());
      return false;
    }
    
    if ( cyl2.shape == id_rectangular )
    {
      label = "/IntrinsicExample/Cylinder2/LengthX";
      if ( !(tpCntl->getInspectedValue(label, ct )) )
      {
        Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
        return false;
      }
      else
      {
        cyl2.length_x = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL;
      }
      
      label="/IntrinsicExample/Cylinder2/LengthY";
      if ( !(tpCntl->getInspectedValue(label, ct )) )
      {
        Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
        return false;
      }
      else
      {
        cyl2.length_y = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL;
      }
    }
    else
    {
      // radius
      label = "/IntrinsicExample/Cylinder2/Radius";
      
      if ( !(tpCntl->getInspectedValue(label, ct )) )
      {
        Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
        return false;
      }
      else
      {
        cyl2.radius = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL;
      }
    }
    
    if ( mode == dim_3d )
    {
      label="/IntrinsicExample/Cylinder2/LengthZ";
      if ( !(tpCntl->getInspectedValue(label, ct )) )
      {
        Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
        return false;
      }
      else
      {
        cyl2.length_z = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL;
      }
    }
    
    // position
    label = "/IntrinsicExample/Cylinder2/PositionX";
    
    if ( !(tpCntl->getInspectedValue(label, ct )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      return false;
    }
    else
    {
      cyl2.x = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL;
    }
    
    label = "/IntrinsicExample/Cylinder2/PositionY";
    
    if ( !(tpCntl->getInspectedValue(label, ct )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      return false;
    }
    else
    {
      cyl2.y = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL;
    }
    
    // 固体媒質
    label = "/IntrinsicExample/Cylinder2/SolidMedium";
    
    if ( !(tpCntl->getInspectedValue(label, str )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      return false;
    }
    cyl2.solid = str;
  }
  

  
  
  // ドライバの設定 値が正の値のとき，有効．ゼロの場合はドライバなし
  label = "/IntrinsicExample/Driver";
  
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
    Hostonly_ stamped_printf("\tError : Value of 'Driver' in 'IntrinsicExample' must be positive.\n");
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
  label = "/IntrinsicExample/FluidMedium";
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  m_fluid = str;
  
  
  
  
  if ( drv_mode == ON )
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
void IP_Cylinder::printPara(FILE* fp, const Control* R)
{
  if ( !fp )
  {
    stamped_printf("\tFail to write into file\n");
    Exit(0);
  }
  
  fprintf(fp,"\n----------\n\n");
  fprintf(fp,"\n\t>> Intrinsic Cylinder Class Parameters\n\n");
  
  fprintf(fp,"\tDimension Mode                     :  %s\n\n", (mode==dim_2d)?"2 Dimensional":"3 Dimensional");
  
  
  REAL_TYPE gap_x;           ///< cylinder間の中心距離 x方向
  REAL_TYPE gap_y;           ///< cylinder間の中心距離 y方向
  
  if ( (cyl1.ens == ON) && (cyl2.ens == ON) )
  {
    // gap X
    REAL_TYPE gx = cyl2.x - cyl1.x;
    REAL_TYPE gy = cyl2.y - cyl1.y;
    
    if ( R->Unit.Param == DIMENSIONAL )
    {
      gap_x = gx;
      gap_y = gy;
    }
    else
    {
      gap_x = gx * RefL;
      gap_y = gy * RefL;
    }
  }
  
  // Cylinder 1
  if ( cyl1.ens == ON )
  {
    if ( cyl1.shape == id_rectangular )
    {
      fprintf(fp,"\tCylinder 1  Shape                           : Rectangular\n");
      fprintf(fp,"\t            Length X            [m] / [-]   : %12.5e / %12.5e\n", cyl1.length_x, cyl1.length_x/RefL);
      fprintf(fp,"\t            Length Y            [m] / [-]   : %12.5e / %12.5e\n", cyl1.length_y, cyl1.length_y/RefL);
    }
    else
    {
      fprintf(fp,"\tCylinder 1  Shape                           : Circular\n");
      fprintf(fp,"\t            Radius              [m] / [-]   : %12.5e / %12.5e\n", cyl1.radius, cyl1.radius/RefL);
    }
    
    if ( mode == dim_3d )
    {
      fprintf(fp,"\t            Length Z            [m] / [-]   : %12.5e / %12.5e\n", cyl1.length_z, cyl1.length_z/RefL);
    }
    
    fprintf(fp,"\t            Position X          [m] / [-]   : %12.5e / %12.5e\n", cyl1.x, cyl1.x/RefL);
    fprintf(fp,"\t            Position Y          [m] / [-]   : %12.5e / %12.5e\n", cyl1.y, cyl1.y/RefL);
    fprintf(fp,"\t            Solid Medium                    : %s\n\n", cyl1.solid.c_str());
  }
  
  
  // Cylinder 2
  if ( cyl2.ens == ON )
  {
    if ( cyl2.shape == id_rectangular )
    {
      fprintf(fp,"\tCylinder 2  Shape                           : Rectangular\n");
      fprintf(fp,"\t            Length X            [m] / [-]   : %12.5e / %12.5e\n", cyl2.length_x, cyl2.length_x/RefL);
      fprintf(fp,"\t            Length Y            [m] / [-]   : %12.5e / %12.5e\n", cyl2.length_y, cyl2.length_y/RefL);
    }
    else
    {
      fprintf(fp,"\tCylinder 2  Shape                           : Circular\n");
      fprintf(fp,"\t            Radius              [m] / [-]   : %12.5e / %12.5e\n", cyl2.radius, cyl2.radius/RefL);
    }
    
    if ( mode == dim_3d )
    {
      fprintf(fp,"\t            Length Z            [m] / [-]   : %12.5e / %12.5e\n", cyl2.length_z, cyl2.length_z/RefL);
    }
    
    fprintf(fp,"\t            Position X          [m] / [-]   : %12.5e / %12.5e\n", cyl2.x, cyl2.x/RefL);
    fprintf(fp,"\t            Position Y          [m] / [-]   : %12.5e / %12.5e\n", cyl2.y, cyl2.y/RefL);
    fprintf(fp,"\t            Solid Medium                    : %s\n\n", cyl2.solid.c_str());
  }
  
  if ( cyl2.ens == ON )
  {
    fprintf(fp,"\tGap between cylinders in X-dir. [m] / [-]   : %12.5e / %12.5e\n", gap_x, gap_x/RefL);
    fprintf(fp,"\tGap between cylinders in Y-dir. [m] / [-]   : %12.5e / %12.5e\n", gap_y, gap_y/RefL);
  }


  
  
  
  if ( drv_length > 0.0 )
  {
    fprintf(fp,"\tDriver Length        [m] / [-]   : %12.5e / %12.5e\n", drv_length, drv_length/RefL);
  }
  
}


// #################################################################
/**
 * @brief 円柱を設定
 * @param [in]     R     　　 Controlクラスのポインタ
 * @param [in]     pos_x     中心座標(x,y) （無次元）
 * @param [in]     pos_y     中心座標(x,y) （無次元）
 * @param [in]     radius    半径 （無次元）
 * @param [in]     len_z     円柱のZ方向の長さ （無次元）
 * @param [in]     mid_solid 固体媒質ID
 * @param [out]    cut       カット情報
 * @param [out]    bid       境界ID
 */
void IP_Cylinder::setCircle(Control* R,
                            const REAL_TYPE pos_x,
                            const REAL_TYPE pos_y,
                            const REAL_TYPE radius,
                            const REAL_TYPE len_z,
                            const int mid_solid,
                            long long* cut,
                            int* bid)
{
  Vec3r pch(pitch);      ///< 無次元格子幅
  Vec3r org(origin);     ///< ノードローカル基点座標の無次元値

  
  REAL_TYPE ph = deltaX; ///< 格子幅（無次元）
  REAL_TYPE cx = pos_x;  ///< 円柱の中心座標（無次元）
  REAL_TYPE cy = pos_y;  ///< 円柱の中心座標（無次元）
  REAL_TYPE rs = radius; ///< 円柱の半径（無次元）
  
  int mid_s = mid_solid;
  
  // 球のbbox
  Vec3r box_min;  ///< Bounding boxの最小値
  Vec3r box_max;  ///< Bounding boxの最大値
  Vec3i box_st;   ///< Bounding boxの始点インデクス
  Vec3i box_ed;   ///< Bounding boxの終点インデクス
  box_min.assign( -(REAL_TYPE)rs+cx, -(REAL_TYPE)rs+cy, 0.0 ); // Z方向はダミー
  box_max.assign(  (REAL_TYPE)rs+cx,  (REAL_TYPE)rs+cy, 0.0 );
  box_st = find_index(box_min, org, pch);
  box_ed = find_index(box_max, org, pch);
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  
  // ローカルな無次元基点座標
  REAL_TYPE ox = origin[0];
  REAL_TYPE oy = origin[1];
  REAL_TYPE oz = origin[2];
  
  // グローバルな無次元座標
  REAL_TYPE ze;
  
  if ( mode == dim_2d )
  {
    ze = oz + region[2] + ph; // Z方向には全セルを対象, dhは安全係数
  }
  else
  {
    ze = oz + len_z; // Z-面からの距離
  }
  
  
  // カット情報
  Vec3r p[5];
  Vec3r base;  // セルのシフト量
  Vec3r b;     // セルセンタ座標
  REAL_TYPE lb[5]; // 内外判定フラグ
  
  for (int k=1; k<=kx; k++) {
    for (int j=box_st.y-1; j<=box_ed.y+1; j++) {
      for (int i=box_st.x-1; i<=box_ed.x+1; i++) {
        
        REAL_TYPE z = org.z + 0.5*ph + ph*(REAL_TYPE)(k-1);
        
        if ( z <= ze )
        {
          base.assign((REAL_TYPE)i-0.5, (REAL_TYPE)j-0.5, (REAL_TYPE)k-0.5);
          b = org + base*ph;
          
          p[0].assign(b.x   , b.y   , b.z   ); // p
          p[1].assign(b.x-ph, b.y   , b.z   ); // w
          p[2].assign(b.x+ph, b.y   , b.z   ); // e
          p[3].assign(b.x   , b.y-ph, b.z   ); // s
          p[4].assign(b.x   , b.y+ph, b.z   ); // n
          
          // (cx, cy, *)が球の中心
          for (int l=0; l<5; l++) {
            REAL_TYPE x = p[l].x - cx;
            REAL_TYPE y = p[l].y - cy;
            REAL_TYPE rr = sqrt(x*x + y*y);
            lb[l] = ( rr <= rs ) ? -1.0 : 1.0; // 内側がマイナス
          }
          
          // cut test  注意！　インデクスが1-4
          for (int l=1; l<=4; l++)
          {
            if ( lb[0]*lb[l] < 0.0 )
            {
              REAL_TYPE s = cut_line_2d(p[0], l, rs, ph);
              size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
              
              int r = (int)quantize9(s);
              setBit10(cut[m], r, l-1);
              setBit5(bid[m], mid_s, l-1);
              
              int rr = (int)quantize9(1.0-s);
              size_t m1;
              
              switch (l-1) {
                case X_minus:
                  m1 = _F_IDX_S3D(i-1, j, k, ix, jx, kx, gd);
                  setBit5(bid[m1], mid_s, X_plus);
                  setBit10(cut[m1], rr, X_plus);
                  break;
                  
                case X_plus:
                  m1 = _F_IDX_S3D(i+1, j, k, ix, jx, kx, gd);
                  setBit5(bid[m1], mid_s, X_minus);
                  setBit10(cut[m1], rr, X_minus);
                  break;
                  
                case Y_minus:
                  m1 = _F_IDX_S3D(i, j-1, k, ix, jx, kx, gd);
                  setBit5(bid[m1], mid_s, Y_plus);
                  setBit10(cut[m1], rr, Y_plus);
                  break;
                  
                case Y_plus:
                  m1 = _F_IDX_S3D(i, j+1, k, ix, jx, kx, gd);
                  setBit5(bid[m1], mid_s, Y_minus);
                  setBit10(cut[m1], rr, Y_minus);
                  break;
              }

            }
          }
        } // if z branch

      }
    }
  }
  
  
  // Z+方向
#pragma omp parallel for firstprivate(ix, jx, kx, gd, mid_s, ox, oy, oz, ph, ze, rs, cx, cy) schedule(static)
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        REAL_TYPE x = ox + 0.5*ph + ph*(i-1); // position of cell center
        REAL_TYPE y = oy + 0.5*ph + ph*(j-1);
        REAL_TYPE z = oz + 0.5*ph + ph*(k-1);
        REAL_TYPE s = (ze - z)/ph;
        
        REAL_TYPE xx = x - cx;
        REAL_TYPE yy = y - cy;
        
        if ( (xx*xx + yy*yy) <= rs*rs )
        {
          
          if ( (z <= ze) && (ze < z+ph) )
          {
            setBit5(bid[m], mid_s, Z_plus);
            int r = (int)quantize9(s);
            setBit10(cut[m], r, Z_plus);

            size_t m1 = _F_IDX_S3D(i, j, k+1, ix, jx, kx, gd);
            setBit5(bid[m1], mid_s, Z_minus);
            int rr = (int)quantize9(1.0-s);
            setBit10(cut[m1], rr, Z_minus);
          }
          else if ( (z-ph < ze) && (ze < z) )
          {
            setBit5(bid[m], mid_s, Z_minus);
            int r = (int)quantize9(-s);
            setBit10(cut[m], r, Z_minus);

            size_t m1 = _F_IDX_S3D(i, j, k-1, ix, jx, kx, gd);
            setBit5(bid[m1], mid_s, Z_plus);
            int rr = (int)quantize9(1.0+s);
            setBit10(cut[m1], rr, Z_plus);
          }
          
        }
        
      }
    }
  }
  
  
}


// #################################################################
/**
 * @brief Rectangularを設定
 * @param [in]     R     　　 Controlクラスのポインタ
 * @param [in]     pos_x     角柱の中心座標(x,y) （無次元）
 * @param [in]     pos_y     角柱の中心座標(x,y)　（無次元）
 * @param [in]     len_x     角柱のX方向の長さ　（無次元）
 * @param [in]     len_y     角柱のY方向の長さ　（無次元）
 * @param [in]     len_z     角柱のZ方向の長さ　（無次元）
 * @param [in]     mid_solid 固体媒質ID
 * @param [out]    cut       カット情報
 * @param [out]    bid       境界ID
 */
void IP_Cylinder::setRect(Control* R,
                          const REAL_TYPE pos_x,
                          const REAL_TYPE pos_y,
                          const REAL_TYPE len_x,
                          const REAL_TYPE len_y,
                          const REAL_TYPE len_z,
                          const int mid_solid,
                          long long* cut,
                          int* bid)
{
  // ローカル
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  REAL_TYPE dh = deltaX;
  
  int mid_s = mid_solid;
  
  // ローカルな無次元基点座標
  REAL_TYPE ox = origin[0];
  REAL_TYPE oy = origin[1];
  REAL_TYPE oz = origin[2];
  
  // グローバルな無次元座標　角柱の頂点
  REAL_TYPE xs = pos_x - 0.5*len_x;
  REAL_TYPE xe = pos_x + 0.5*len_x;
  REAL_TYPE ys = pos_y - 0.5*len_y;
  REAL_TYPE ye = pos_y + 0.5*len_y;
  REAL_TYPE ze;
  
  if ( mode == dim_2d )
  {
    ze = origin[2] + region[2] + dh; // Z方向には全セルを対象, dhは安全係数
  }
  else
  {
    ze = origin[2] + len_z; // Z-面からの距離
  }
  

  /*
    <-dh->
      1     2    ...
   |     |     |     |  *  |
   +-----+-----+-----+-----+
   ^                    ^
   ox                  x_i
   
   x_i = ox + 0.5*dh + dh*(i-1)
   
   
   
   <- s ->
   |     |         |
   +---------------+
  x_i    c        x_{i+1}
   
   s = (c-x_i)/dh > 0
   
   d_i^+     = s
   d_{i+1}^- = 1-s
   b_i^+     = solid
   b_{i+1}^- = solid
   
   
   
           <-  s  ->
   |      |        |
   +---------------+
 x_{i-1}  c       x_i
   
   s = (c-x_i)/dh < 0
   
   d_{i-1}^+ = 1+s
   d_i^-     = -s
   b_{i-1}^+ = solid
   b_i^-     = solid
   
   */
  
  
  // Y+ face of Rect
#pragma omp parallel for firstprivate(ix, jx, kx, gd, mid_s, ox, oy, oz, dh, xs, xe, ye, ze) schedule(static)
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        REAL_TYPE x = ox + 0.5*dh + dh*(i-1); // position of cell center
        REAL_TYPE y = oy + 0.5*dh + dh*(j-1);
        REAL_TYPE z = oz + 0.5*dh + dh*(k-1);
        REAL_TYPE s = (ye - y)/dh;
        
        if ( z <= ze )
        {
          if ( (xs <= x) && (x <= xe) )
          {
            if ( (y <= ye) && (ye < y+dh) )
            {
              setBit5(bid[m], mid_s, Y_plus);
              int r = (int)quantize9(s);
              setBit10(cut[m], r, Y_plus);
              
              size_t m1 = _F_IDX_S3D(i, j+1, k, ix, jx, kx, gd);
              setBit5(bid[m1], mid_s, Y_minus);
              int rr = (int)quantize9(1.0-s);
              setBit10(cut[m1], rr, Y_minus);
            }
            else if ( (y-dh < ye) && (ye < y) )
            {
              setBit5(bid[m], mid_s, Y_minus);
              int r = (int)quantize9(-s);
              setBit10(cut[m], r, Y_minus);
              
              size_t m1 = _F_IDX_S3D(i, j-1, k, ix, jx, kx, gd);
              setBit5(bid[m1], mid_s, Y_plus);
              int rr = (int)quantize9(1.0+s);
              setBit10(cut[m1], rr, Y_plus);
            }
          }
        } // Z
        
      }
    }
  }

  // Y- face of Rect
#pragma omp parallel for firstprivate(ix, jx, kx, gd, mid_s, ox, oy, oz, dh, xs, xe, ys, ze) schedule(static)
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        REAL_TYPE x = ox + 0.5*dh + dh*(i-1); // position of cell center
        REAL_TYPE y = oy + 0.5*dh + dh*(j-1);
        REAL_TYPE z = oz + 0.5*dh + dh*(k-1);
        REAL_TYPE s = (ys - y)/dh;
        
        if ( z <= ze )
        {
          if ( (xs <= x) && (x <= xe) )
          {
            if ( (y <= ys) && (ys < y+dh) )
            {
              setBit5(bid[m], mid_s, Y_plus);
              int r = (int)quantize9(s);
              setBit10(cut[m], r, Y_plus);

              size_t m1 = _F_IDX_S3D(i, j+1, k, ix, jx, kx, gd);
              setBit5(bid[m1], mid_s, Y_minus);
              int rr = (int)quantize9(1.0-s);
              setBit10(cut[m1], rr, Y_minus);
            }
            else if ( (y-dh < ys) && (ys < y) )
            {
              setBit5(bid[m], mid_s, Y_minus);
              int r = (int)quantize9(-s);
              setBit10(cut[m], r, Y_minus);
              
              size_t m1 = _F_IDX_S3D(i, j-1, k, ix, jx, kx, gd);
              setBit5(bid[m1], mid_s, Y_plus);
              int rr = (int)quantize9(1.0+s);
              setBit10(cut[m1], rr, Y_plus);
            }
          }
        } // Z
        
      }
    }
  }

  // X+ face of Rect
#pragma omp parallel for firstprivate(ix, jx, kx, gd, mid_s, ox, oy, oz, dh, xe, ys, ye, ze) schedule(static)
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        REAL_TYPE x = ox + 0.5*dh + dh*(i-1); // position of cell center
        REAL_TYPE y = oy + 0.5*dh + dh*(j-1);
        REAL_TYPE z = oz + 0.5*dh + dh*(k-1);
        REAL_TYPE s = (xe - x)/dh;
        
        if ( z <= ze )
        {
          if ( (ys <= y) && (y <= ye) )
          {
            if ( (x <= xe) && (xe < x+dh) )
            {
              setBit5(bid[m], mid_s, X_plus);
              int r = (int)quantize9(s);
              setBit10(cut[m], r, X_plus);

              size_t m1 = _F_IDX_S3D(i+1, j, k, ix, jx, kx, gd);
              setBit5(bid[m1], mid_s, X_minus);
              int rr = (int)quantize9(1.0-s);
              setBit10(cut[m1], rr, X_minus);
            }
            else if ( (x-dh < xe) && (xe < x) )
            {
              setBit5(bid[m], mid_s, X_minus);
              int r = (int)quantize9(-s);
              setBit10(cut[m], r, X_minus);

              size_t m1 = _F_IDX_S3D(i-1, j, k, ix, jx, kx, gd);
              setBit5(bid[m1], mid_s, X_plus);
              int rr = (int)quantize9(1.0+s);
              setBit10(cut[m1], r, X_plus);
            }
          }
        } // Z
        
      }
    }
  }

  // X- face of Rect
#pragma omp parallel for firstprivate(ix, jx, kx, gd, mid_s, ox, oy, oz, dh, xs, ys, ye, ze) schedule(static)
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        REAL_TYPE x = ox + 0.5*dh + dh*(i-1); // position of cell center
        REAL_TYPE y = oy + 0.5*dh + dh*(j-1);
        REAL_TYPE z = oz + 0.5*dh + dh*(k-1);
        REAL_TYPE s = (xs - x)/dh;
        
        if ( z <= ze )
        {
          if ( (ys <= y) && (y <= ye) )
          {
            if ( (x <= xs) && (xs < x+dh) )
            {
              setBit5(bid[m], mid_s, X_plus);
              int r = (int)quantize9(s);
              setBit10(cut[m], r, X_plus);

              size_t m1 = _F_IDX_S3D(i+1, j, k, ix, jx, kx, gd);
              setBit5(bid[m1], mid_s, X_minus);
              int rr = (int)quantize9(1.0-s);
              setBit10(cut[m1], rr, X_minus);
            }
            else if ( (x-dh < xs) && (xs < x) )
            {
              setBit5(bid[m], mid_s, X_minus);
              int r = (int)quantize9(-s);
              setBit10(cut[m], r, X_minus);
              
              size_t m1 = _F_IDX_S3D(i-1, j, k, ix, jx, kx, gd);
              setBit5(bid[m1], mid_s, X_plus);
              int rr = (int)quantize9(1.0+s);
              setBit10(cut[m1], r, X_plus);
            }
          }
        } // Z
        
      }
    }
  }
  
}


// #################################################################
// Cylinderの計算領域のカット情報を設定する
void IP_Cylinder::setup(int* bcd, Control* R, const int NoMedium, const MediumList* mat, long long* cut, int* bid)
{
  int mid_fluid;        ///< 流体
  int mid_solid1;       ///< 固体1
  int mid_solid2;       ///< 固体2
  
  // 流体
  if ( (mid_fluid = FBUtility::findIDfromLabel(mat, NoMedium, m_fluid)) == 0 )
  {
    Hostonly_ printf("\tLabel '%s' is not listed in MediumList\n", m_fluid.c_str());
    Exit(0);
  }
  
  // 固体
  if (cyl1.ens==ON)
  {
    if ( (mid_solid1 = FBUtility::findIDfromLabel(mat, NoMedium, cyl1.solid)) == 0 )
    {
      Hostonly_ printf("\tLabel '%s' is not listed in MediumList\n", cyl1.solid.c_str());
      Exit(0);
    }
  }
  
  if (cyl2.ens==ON)
  {
    if ( (mid_solid2 = FBUtility::findIDfromLabel(mat, NoMedium, cyl2.solid)) == 0 )
    {
      Hostonly_ printf("\tLabel '%s' is not listed in MediumList\n", cyl2.solid.c_str());
      Exit(0);
    }
  }
  
  
  // 媒質数チェック
  if ( (cyl1.ens==ON) && (cyl2.ens==ON) )
  {
    if ( NoMedium < 3 )
    {
      Hostonly_ printf("\t The numbe of medium must be more than 3 for two cylinders\n");
      Exit(0);
    }
  }
  
  
  if ( cyl1.ens == ON )
  {
    REAL_TYPE pos_x = cyl1.x/R->RefLength;
    REAL_TYPE pos_y = cyl1.y/R->RefLength;
    REAL_TYPE len_x = cyl1.length_x/R->RefLength;
    REAL_TYPE len_y = cyl1.length_y/R->RefLength;
    REAL_TYPE len_z = cyl1.length_z/R->RefLength;
    REAL_TYPE rad   = cyl1.radius/R->RefLength;
    
    if ( cyl1.shape == id_rectangular )
    {
      setRect(R, pos_x, pos_y, len_x, len_y, len_z, mid_solid1, cut, bid);
    }
    else
    {
      setCircle(R, pos_x, pos_y, rad, len_z, mid_solid1, cut, bid);
    }
  }
  
  
  if ( cyl2.ens == ON )
  {
    REAL_TYPE pos_x = cyl2.x/R->RefLength;
    REAL_TYPE pos_y = cyl2.y/R->RefLength;
    REAL_TYPE len_x = cyl2.length_x/R->RefLength;
    REAL_TYPE len_y = cyl2.length_y/R->RefLength;
    REAL_TYPE len_z = cyl2.length_z/R->RefLength;
    REAL_TYPE rad   = cyl2.radius/R->RefLength;
    
    if ( cyl2.shape == id_rectangular )
    {
      setRect(R, pos_x, pos_y, len_x, len_y, len_z, mid_solid2, cut, bid);
    }
    else
    {
      setCircle(R, pos_x, pos_y, rad, len_z, mid_solid2, cut, bid);
    }
  }


  
  
  /* driver設定 iff ドライバ長が正の場合
  if ( drv_mode == OFF ) return;
  
  
  if ( (mid_driver = FBUtility::findIDfromLabel(mat, NoMedium, m_driver)) == 0 )
  {
    Hostonly_ printf("\tLabel '%s' is not listed in MediumList\n", m_driver.c_str());
    Exit(0);
  }
  
  if ( (mid_driver_face = FBUtility::findIDfromLabel(mat, NoMedium, m_driver_face)) == 0 )
  {
    Hostonly_ printf("\tLabel '%s' is not listed in MediumList\n", m_driver_face.c_str());
    Exit(0);
  }
  
  // ドライバ部分
  if ( drv_length > 0.0 )
  {
    if ( nID[X_minus] < 0 )
    {
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          for (int i=1; i<=ix; i++) {
            size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
            x = ox + 0.5*dh + dh*(i-1);
            if ( x < ox+len ) bcd[m] |= mid_driver;
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
          if ( ( DECODE_CMP( bcd[m] ) == mid_driver) && (DECODE_CMP( bcd[m1] ) == mid_fluid) )
          {
            bcd[m] |= mid_driver_face;
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
          bcd[m] |= mid_solid;
        }
      }
    }
  }
   */
}
