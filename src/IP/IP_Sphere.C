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
 * @file   IP_Sphere.C
 * @brief  IP_Sphere class
 * @author kero
 */

#include "IP_Sphere.h"


// #################################################################
//@brief パラメータを取得する
bool IP_Sphere::getTP(Control* R, TPControl* tpCntl)
{
  std::string str;
  std::string label;
  REAL_TYPE ct;
  
  // radius
  label = "/Parameter/Intrinsic_Example/radius";
  
  if ( !tpCntl->GetValue(label, &ct) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  else
  {
	  radius = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL;
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
  
  // 媒質指定
  label = "/Parameter/Intrinsic_Example/Fluid_medium";
  if ( !tpCntl->GetValue(label, &str) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  m_fluid = str;
  
  
  label = "/Parameter/Intrinsic_Example/Solid_medium";
  if ( !tpCntl->GetValue(label, &str) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  m_solid = str;
  
  
  if (drv_length < 0.0 )
  {
    label = "/Parameter/Intrinsic_Example/driver_medium";
    if ( !tpCntl->GetValue(label, &str) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      return false;
    }
    m_driver = str;
    
    
    label = "/Parameter/Intrinsic_Example/driver_face_medium";
    if ( !tpCntl->GetValue(label, &str) )
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
  
  fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
  fprintf(fp,"\n\t>> Intrinsic Sphere Parameters\n\n");
  
  fprintf(fp,"\tRadius of Sphere       [m] / [-]   : %12.5e / %12.5e\n", radius, radius/RefL);
  
  if ( drv_mode == ON )
  {
    fprintf(fp,"\tDriver Length        [m] / [-]   : %12.5e / %12.5e\n", drv_length, drv_length/RefL);
  }
}


// #################################################################
// 領域情報を設定する
void IP_Sphere::setDomain(Control* R, const int* sz, REAL_TYPE* m_org, REAL_TYPE* m_reg, REAL_TYPE* m_pch)
{
  RefL = R->RefLength;
  
  pch.x = (float)m_pch[0];
  pch.y = (float)m_pch[1];
  pch.z = (float)m_pch[2];
  
  org.x = (float)m_org[0];
  org.y = (float)m_org[1];
  org.z = (float)m_org[2];
  
  // チェック
  if ( (pch.x != pch.y) || (pch.y != pch.z) )
  {
    Hostonly_ printf("Error : 'VoxelPitch' in each direction must be same.\n");
    Exit(0);
  }
  
  // 領域サイズ
  wth.x = pch.x * (float)size[0];
  wth.y = pch.y * (float)size[1];
  wth.z = pch.z * (float)size[2];
  
  m_reg[0] = wth.x;
  m_reg[1] = wth.y;
  m_reg[2] = wth.z;

}


// #################################################################
// 点pの属するセルインデクスを求める
// Fortran index
FB::Vec3i IP_Sphere::find_index(const FB::Vec3f p, const FB::Vec3f ol)
{
  FB::Vec3f q = (p-ol)/pch;
  FB::Vec3i idx( ceil(q.x), ceil(q.y), ceil(q.z) );
  
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
// 計算領域のセルIDを設定する（バイナリー）
void IP_Sphere::setup(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, MediumList* mat)
{
  int mid_fluid;        /// 流体
  int mid_solid;        /// 固体
  int mid_driver;       /// ドライバ部
  int mid_driver_face;  /// ドライバ流出面

  REAL_TYPE x, y, z, dh, len;
  REAL_TYPE ox, oy, oz, Lx, Ly, Lz;
  REAL_TYPE ox_g, oy_g, oz_g;
  
  FB::Vec3f base, b, t, org_l;
  REAL_TYPE ph = pch.x;
  REAL_TYPE r;
  REAL_TYPE rs = radius/R->RefLength;
  
  // ノードローカルの無次元値
  org_l.x = (float)origin[0];
  org_l.y = (float)origin[1];
  org_l.z = (float)origin[2];
  Lx = region[0];
  Ly = region[1];
  Lz = region[2];
  dh = deltaX;
  ox_g = G_origin[0];
  oy_g = G_origin[1];
  oz_g = G_origin[2];
  
  // ローカルにコピー
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  //printf("%d : ox = %e, oy = %e, oz = %e\n", pn.myrank, org_l.x, org_l.y, org_l.z);
  //printf("%d : Lx = %e, Ly = %e, Lz = %e\n", pn.myrank, Lx, Ly, Lz);
  //printf("%d : oxG= %e, oyG= %e, ozG= %e\n", pn.myrank, ox_g, oy_g, oz_g);
  
  // 球のbbox 球の中心座標はゼロ
  box_min = - rs;
  box_max = + rs;
  box_st = find_index(box_min, org_l);
  box_ed = find_index(box_max, org_l);
  
  //printf("%d : st_x = %d, st_y = %d, st_z = %d\n", pn.myrank, box_st.x, box_st.y, box_st.z);
  //printf("%d : ed_x = %d, ed_y = %d, ed_z = %d\n", pn.myrank, box_ed.x, box_ed.y, box_ed.z);
  
  // 媒質設定
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
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd, mid_fluid) schedule(static)
  for (int k=1-gd; k<=kx+gd; k++) {
    for (int j=1-gd; j<=jx+gd; j++) {
      for (int i=1-gd; i<=ix+gd; i++) {
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        mid[m] = mid_fluid;
      }
    }
  }
  
  // 球内部
  for (int k=box_st.z; k<=box_ed.z; k++) { 
    for (int j=box_st.y; j<=box_ed.y; j++) {
      for (int i=box_st.x; i<=box_ed.x; i++) {
        
        base.assign((float)i-0.5, (float)j-0.5, (float)k-0.5);
        b = org_l + base*ph;
        r = b.length();
        
        if ( r <= rs )
        {
          size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
          mid[m] = mid_solid;
        }
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
  
  // @todo 以下は確認のこと
  
  // lengthは有次元値
  len = ox_g + (drv_length)/R->RefLength; // グローバルな無次元位置
  
  // ドライバ部分　X-面からドライバ長さより小さい領域
  if ( drv_length > 0.0 )
  {
    for (int k=1; k<=kx; k++) {
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
          x = ox + 0.5*dh + dh*(i-1);
          if ( x < len ) mid[m] = mid_driver;
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
        if ( x < len )
        {
          mid[m] = mid_solid;
        }
      }
    }
  }
}


// #################################################################
// 計算領域のセルIDとカット情報を設定する
void IP_Sphere::setup_cut(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, MediumList* mat, float* cut)
{
  int mid_fluid;        /// 流体
  int mid_solid;        /// 固体
  int mid_driver;       /// ドライバ部
  int mid_driver_face;  /// ドライバ流出面

  REAL_TYPE x, y, z, dh, len;
  REAL_TYPE ox, oy, oz, Lx, Ly, Lz;
  REAL_TYPE ox_g, oy_g, oz_g;
  
  FB::Vec3f base, b, t, org_l;
  REAL_TYPE ph = pch.x;
  REAL_TYPE r;
  REAL_TYPE rs = radius/R->RefLength;
  
  // ノードローカルの無次元値
  org_l.x = (float)origin[0];
  org_l.y = (float)origin[1];
  org_l.z = (float)origin[2];
  Lx = region[0];
  Ly = region[1];
  Lz = region[2];
  dh = deltaX;

  ox_g = G_origin[0];
  oy_g = G_origin[1];
  oz_g = G_origin[2];
  
  // 球のbbox
  box_min = - rs;
  box_max = + rs;
  box_st = find_index(box_min, org_l);
  box_ed = find_index(box_max, org_l);
  
  // ローカルにコピー
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  // 媒質設定
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
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd, mid_fluid) schedule(static)
  for (int k=1-gd; k<=kx+gd; k++) {
    for (int j=1-gd; j<=jx+gd; j++) {
      for (int i=1-gd; i<=ix+gd; i++) {
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        mid[m] = mid_fluid;
      }
    }
  }
  
  // 球内部
  
  for (int k=box_st.z; k<=box_ed.z; k++) { 
    for (int j=box_st.y; j<=box_ed.y; j++) {
      for (int i=box_st.x; i<=box_ed.x; i++) {
       
        base.assign((float)i-0.5, (float)j-0.5, (float)k-0.5);
        b = org_l + base*ph;
        r = b.length();
        
        if ( r <= rs )
        {
          size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
          mid[m] = mid_solid;
        }
      }
    }
  }
  
  // カット情報
  FB::Vec3f p[7];
  float lb[7], s, r_min=10.0, r_max=0.0;
  
  for (int k=box_st.z; k<=box_ed.z; k++) { 
    for (int j=box_st.y; j<=box_ed.y; j++) {
      for (int i=box_st.x; i<=box_ed.x; i++) {
        
        base.assign((float)i-0.5, (float)j-0.5, (float)k-0.5);
        b = org + base*ph;
        
        p[0].assign(b.x   , b.y   , b.z   ); // p
        p[1].assign(b.x-ph, b.y   , b.z   ); // w 
        p[2].assign(b.x+ph, b.y   , b.z   ); // e
        p[3].assign(b.x   , b.y-ph, b.z   ); // s
        p[4].assign(b.x   , b.y+ph, b.z   ); // n
        p[5].assign(b.x   , b.y   , b.z-ph); // b
        p[6].assign(b.x   , b.y   , b.z+ph); // t
        
        for (int l=0; l<7; l++) {
          lb[l] = ( p[l].length() <= rs ) ? -1.0 : 1.0; // 内側がマイナス
        }
        
        // cut test
        for (int l=1; l<=6; l++) {
          if ( lb[0]*lb[l] < 0.0 )
          {
            s = cut_line(p[0], l, rs, ph);
            size_t m = _F_IDX_S4DEX(l-1, i, j, k, 6, ix, jx, kx, gd); // 注意！　インデクスが1-6
            cut[m] = s;
            r_min = std::min(r_min, s);
            r_max = std::max(r_max, s);
          }
        }
      }
    }
  }
  
  Hostonly_
  {
    printf("\nCut info. for Sphere\n");
    printf("\tmin. cut = %f\n", r_min);
    printf("\tmax. cut = %f\n", r_max);
  }

  
  // driver設定 iff ドライバ長が正の場合
  if ( drv_mode == OFF ) return;
  
  // lengthは有次元値
  len = ox_g + (drv_length)/R->RefLength; // グローバルな無次元位置
  
  // ドライバ部分　X-面からドライバ長さより小さい領域
  if ( drv_length > 0.0 )
  {
    for (int k=1; k<=kx; k++) {
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
          x = ox + 0.5*dh + dh*(i-1);
          if ( x < len ) mid[m] = mid_driver;
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
        if ( x < len )
        {
          mid[m] = mid_solid;
        }
      }
    }
  }
}


// #################################################################
// 交点計算
float IP_Sphere::cut_line(const FB::Vec3f p, const int dir, const float r, const float dh)
{
  float x, y, z, s;
  float xc, yc, zc;

  x = p.x;
  y = p.y;
  z = p.z;
  
  s = 0.0;
  
  // 基点座標の符号で交点座標を判断
  switch (dir)
  {
    case 1: // X-
      xc = sqrtf(r*r - y*y - z*z);
      if ( x < 0.0 ) xc *= -1.0;
      s = fabs(xc-x);
      break;
      
    case 2: // X+
      xc = sqrtf(r*r - y*y - z*z);
      if ( x < 0.0 ) xc *= -1.0;
      s = fabs(xc-x);
      break;
      
    case 3: // Y-
      yc = sqrtf(r*r - x*x - z*z);
      if ( y < 0.0 ) yc *= -1.0;
      s = fabs(yc-y);
      break;
      
    case 4: // Y+
      yc = sqrtf(r*r - x*x - z*z);
      if ( y < 0.0 ) yc *= -1.0;
      s = fabs(yc-y);
      break;
      
    case 5: // Z-
      zc = sqrtf(r*r - x*x - y*y);
      if ( z < 0.0 ) zc *= -1.0;
      s = fabs(zc-z);
      break;
      
    case 6: // Z+
      zc = sqrtf(r*r - x*x - y*y);
      if ( z < 0.0 ) zc *= -1.0;
      s = fabs(zc-z);
      break;
      
    default:
      Exit(0);
      break;
  }
  
  return s/dh;
}
