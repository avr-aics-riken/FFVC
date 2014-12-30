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
 * @file   IP_Sphere.C
 * @brief  IP_Sphere class
 * @author aics
 */

#include "IP_Sphere.h"

using namespace Vec3class;


// #################################################################
// 交点の無次元距離を計算する
REAL_TYPE IP_Sphere::cut_line(const Vec3r p, const int dir, const REAL_TYPE r, const REAL_TYPE dh)
{
  REAL_TYPE x, y, z, s;
  REAL_TYPE c;
  
  x = p.x;
  y = p.y;
  z = p.z;
  
  s = 0.0;
  
  // 基点座標の符号で交点座標を判断
  switch (dir)
  {
    case 1: // X-
      c = sqrtf(r*r - y*y - z*z);
      if ( x < 0.0 ) c *= -1.0;
      s = fabs(c-x);
      break;
      
    case 2: // X+
      c = sqrtf(r*r - y*y - z*z);
      if ( x < 0.0 ) c *= -1.0;
      s = fabs(c-x);
      break;
      
    case 3: // Y-
      c = sqrtf(r*r - x*x - z*z);
      if ( y < 0.0 ) c *= -1.0;
      s = fabs(c-y);
      break;
      
    case 4: // Y+
      c = sqrtf(r*r - x*x - z*z);
      if ( y < 0.0 ) c *= -1.0;
      s = fabs(c-y);
      break;
      
    case 5: // Z-
      c = sqrtf(r*r - x*x - y*y);
      if ( z < 0.0 ) c *= -1.0;
      s = fabs(c-z);
      break;
      
    case 6: // Z+
      c = sqrtf(r*r - x*x - y*y);
      if ( z < 0.0 ) c *= -1.0;
      s = fabs(c-z);
      break;
      
    default:
      Exit(0);
      break;
  }
  
  return s/dh;
}



// #################################################################
//  点pの属するセルインデクスを求める
// Fortran index
Vec3i IP_Sphere::find_index(const Vec3r p, const Vec3r ol, const Vec3r pch)
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
bool IP_Sphere::getTP(Control* R, TextParser* tpCntl)
{
  std::string str;
  std::string label;
  REAL_TYPE ct;
  
  // radius
  label = "/IntrinsicExample/Radius";
  
  if ( !(tpCntl->getInspectedValue(label, ct)) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  else
  {
	  radius = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL;
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
  if ( !(tpCntl->getInspectedValue(label, str)) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  m_fluid = str;
  
  
  label = "/IntrinsicExample/SolidMedium";
  if ( !(tpCntl->getInspectedValue(label, str)) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  m_solid = str;
  
  
  if (drv_length > 0.0 )
  {
    label = "/IntrinsicExample/DriverMedium";
    if ( !(tpCntl->getInspectedValue(label, str)) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      return false;
    }
    m_driver = str;
    
    
    label = "/IntrinsicExample/DriverFaceMedium";
    if ( !(tpCntl->getInspectedValue(label, str)) )
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
void IP_Sphere::printPara(FILE* fp, const Control* R)
{
  if ( !fp )
  {
    stamped_printf("\tFail to write into file\n");
    Exit(0);
  }
  
  fprintf(fp,"\n----------\n\n");
  fprintf(fp,"\n\t>> Intrinsic Sphere Class Parameters\n\n");
  
  fprintf(fp,"\tRadius of Sphere       [m] / [-]   : %12.5e / %12.5e\n", radius, radius/RefL);
  
  if ( drv_mode == ON )
  {
    fprintf(fp,"\tDriver Length        [m] / [-]   : %12.5e / %12.5e\n", drv_length, drv_length/RefL);
  }
}


// #################################################################
// 計算領域のセルIDとカット情報を設定する
void IP_Sphere::setup(int* bcd, Control* R, const int NoMedium, const MediumList* mat, long long* cut, int* bid)
{
  int mid_fluid;        /// 流体
  int mid_solid;        /// 固体
  int mid_driver;       /// ドライバ部
  int mid_driver_face;  /// ドライバ流出面
  
  REAL_TYPE ph = deltaX;
  REAL_TYPE rs = radius/R->RefLength;
  
  // ノードローカルの無次元値
  REAL_TYPE Lx = region[0];
  REAL_TYPE Ly = region[1];
  REAL_TYPE Lz = region[2];
  REAL_TYPE dh = deltaX;
  
  Vec3r pch(pitch);        ///< セル幅
  Vec3r org(origin);       ///< 計算領域の基点
  
  // 球のbbox
  Vec3r box_min;    ///< Bounding boxの最小値
  Vec3r box_max;    ///< Bounding boxの最大値
  Vec3i box_st;     ///< Bounding boxの始点インデクス
  Vec3i box_ed;     ///< Bounding boxの終点インデクス
  
  box_min = - rs;
  box_max = + rs;
  box_st = find_index(box_min, org, pch);
  box_ed = find_index(box_max, org, pch);
  
  
  // ローカルにコピー
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
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
  
  
  
  // カット情報
  Vec3r p[7];
  Vec3r base, b;
  REAL_TYPE r_min=10.0, r_max=0.0;
  REAL_TYPE lb[7];
  
  for (int k=box_st.z-2; k<=box_ed.z+2; k++) {
    for (int j=box_st.y-2; j<=box_ed.y+2; j++) {
      for (int i=box_st.x-2; i<=box_ed.x+2; i++) {

        base.assign((REAL_TYPE)i-0.5, (REAL_TYPE)j-0.5, (REAL_TYPE)k-0.5);
        b = org + base*ph;
        
        p[0].assign(b.x   , b.y   , b.z   ); // p
        p[1].assign(b.x-ph, b.y   , b.z   ); // w 
        p[2].assign(b.x+ph, b.y   , b.z   ); // e
        p[3].assign(b.x   , b.y-ph, b.z   ); // s
        p[4].assign(b.x   , b.y+ph, b.z   ); // n
        p[5].assign(b.x   , b.y   , b.z-ph); // b
        p[6].assign(b.x   , b.y   , b.z+ph); // t
        
        // (0.0, 0.0, 0.0)が球の中心
        for (int l=0; l<7; l++) {
          lb[l] = ( p[l].length() <= rs ) ? -1.0 : 1.0; // 内側がマイナス
        }
        
        
        // cut test // 注意！　インデクスが1-6
        for (int l=1; l<=6; l++) {
          if ( lb[0]*lb[l] < 0.0 )
          {
            REAL_TYPE s = cut_line(p[0], l, rs, ph);
            size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
            int r = quantize9(s);
            setBit10(cut[m], r, l-1);
            setBit5(bid[m], mid_solid, l-1);
            
            int rr = quantize9(1.0-s);
            size_t m1;
            
            switch (l-1) {
              case X_minus:
                m1 = _F_IDX_S3D(i-1, j, k, ix, jx, kx, gd);
                setBit5(bid[m1], mid_solid, X_plus);
                setBit10(cut[m1], rr, X_plus);
                break;
                
              case X_plus:
                m1 = _F_IDX_S3D(i+1, j, k, ix, jx, kx, gd);
                setBit5(bid[m1], mid_solid, X_minus);
                setBit10(cut[m1], rr, X_minus);
                break;
                
              case Y_minus:
                m1 = _F_IDX_S3D(i, j-1, k, ix, jx, kx, gd);
                setBit5(bid[m1], mid_solid, Y_plus);
                setBit10(cut[m1], rr, Y_plus);
                break;
                
              case Y_plus:
                m1 = _F_IDX_S3D(i, j+1, k, ix, jx, kx, gd);
                setBit5(bid[m1], mid_solid, Y_minus);
                setBit10(cut[m1], rr, Y_minus);
                break;
                
              case Z_minus:
                m1 = _F_IDX_S3D(i, j, k-1, ix, jx, kx, gd);
                setBit5(bid[m1], mid_solid, Z_plus);
                setBit10(cut[m1], rr, Z_plus);;
                break;
                
              case Z_plus:
                m1 = _F_IDX_S3D(i, j, k+1, ix, jx, kx, gd);
                setBit5(bid[m1], mid_solid, Z_minus);
                setBit10(cut[m1], rr, Z_minus);
                break;
            }
            
            //printf("(%2d %2d %2d) %2d %d %f\n", i,j,k,mid_solid, l-1, s);
            r_min = (std::min)(r_min, s);
            r_max = (std::max)(r_max, s);
          }
        }
      }
    }
  }
  
  Hostonly_
  {
    printf("\n\tCut info. for Sphere\n");
    printf("\tmin. cut = %f\n", r_min);
    printf("\tmax. cut = %f\n", r_max);
  }

  
  // driver設定 iff ドライバ長が正の場合
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
  
  REAL_TYPE x, y, z;
  REAL_TYPE ox, oy, oz;
  
  // lengthは有次元値
  REAL_TYPE len = G_origin[0] + (drv_length)/R->RefLength; // グローバルな無次元位置
  
  // ドライバ部分　X-面からドライバ長さより小さい領域
  if ( drv_length > 0.0 )
  {
    for (int k=1; k<=kx; k++) {
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
          x = ox + 0.5*dh + dh*(i-1);
          if ( x < len ) bcd[m] |= mid_driver;
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
          if ( (DECODE_CMP(bcd[m]) == mid_driver) && (DECODE_CMP(bcd[m1]) == mid_fluid) )
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
        if ( x < len )
        {
          bcd[m] |= mid_solid;
        }
      }
    }
  }
}
