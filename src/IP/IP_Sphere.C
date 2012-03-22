/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file IP_Shere.C
//@brief IP_Shere class
//@author keno, FSI Team, VCAD, RIKEN

#include "IP_Sphere.h"

/**
 @fn bool IP_Sphere::getXML(SklSolverConfig* CF, Control* R)
 @brief パラメータを取得する
 @param CF コンフィギュレーションツリー
 @param R Controlクラスのポインタ
 */
bool IP_Sphere::getXML(SklSolverConfig* CF, Control* R)
{
  const CfgElem *elemTop=NULL, *elmL1=NULL;
  REAL_TYPE ct=0.0;
  
  if ( !(elemTop = CF->GetTop(PARAMETER)) ) return false;
  
  if( !(elmL1 = elemTop->GetElemFirst("Intrinsic_Example")) ) {
    Hostonly_ stamped_printf("\tParsing error : Missing the section of 'Intrinsic_Example'\n");
    return false;
  }
  
  // offset
  if ( elmL1->GetValue(CfgIdt("offset"), &ct) ) {
    offset = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL;
  }
  else {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Offset' in 'Intrinsic_Example'\n");
    return false;
  }

  // radius
  if ( elmL1->GetValue(CfgIdt("radius"), &ct) ) {
    radius = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL;
  }
  else {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Radius' in 'Intrinsic_Example'\n");
    return false;
  }
  
  // ドライバの設定 値が正の値のとき，有効．ゼロの場合はドライバなし
  if ( elmL1->GetValue(CfgIdt("Driver"), &ct) ) {
    drv_length = ( R->Unit.Param == DIMENSIONAL ) ? ct : ct * RefL;
  }
  else {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Driver' in 'Intrinsic_Example'\n");
    return false;
  }
  if ( drv_length < 0.0 ) {
    Hostonly_ stamped_printf("\tError : Value of 'Driver' in 'Intrinsic_Example' must be positive.\n");
    return false;
  }
  
  return true;
}

/**
 @fn void IP_Sphere::printPara(FILE* fp, Control* R)
 @brief パラメータの表示
 @param fp ファイルポインタ
 @param R コントロールクラスのポインタ
 */
void IP_Sphere::printPara(FILE* fp, Control* R)
{
  if ( !fp ) {
    stamped_printf("\tFail to write into file\n");
    Exit(0);
  }
  
  fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
  fprintf(fp,"\n\t>> Intrinsic Sphere Parameters\n\n");
  
  fprintf(fp,"\tOffset                 [m] / [-]   : %12.5e / %12.5e\n", offset, offset/RefL);
  fprintf(fp,"\tRadius of Sphere       [m] / [-]   : %12.5e / %12.5e\n", radius, radius/RefL);
  if ( drv_length > 0.0 ) {
    fprintf(fp,"\tDriver Length        [m] / [-]   : %12.5e / %12.5e\n", drv_length, drv_length/RefL);
  }
}

/**
 @fn bool IP_Sphere::setDomain(Control* R, unsigned m_sz[3], REAL_TYPE m_org[3], REAL_TYPE m_wth[3], REAL_TYPE m_pch[3])
 @brief 領域情報を設定する
 @param R Controlクラスのポインタ
 @param[in]  m_sz グローバル計算領域のセルサイズ
 @param[out] m_org グローバル計算領域の基点
 @param[out] m_wth グローバル計算領域のbounding boxサイズ
 @param[in]  m_pch セルピッチ
 */
void IP_Sphere::setDomain(Control* R, unsigned m_sz[3], REAL_TYPE m_org[3], REAL_TYPE m_wth[3], REAL_TYPE m_pch[3])
{
  pch = m_pch;
  
  // チェック
  if ( (pch.x != pch.y) || (pch.y != pch.z) ) {
    Hostonly_ printf("Error : 'VoxelPitch' in each direction must be same.\n");
    Exit(0);
  }
  
  // 領域サイズ
  wth.x = pch.x * (float)m_sz[0];
  wth.y = pch.y * (float)m_sz[1];
  wth.z = pch.z * (float)m_sz[2];
  
  m_wth[0] = wth.x;
  m_wth[1] = wth.y;
  m_wth[2] = wth.z;
  
  // 領域の基点
  org = - 0.5*wth;
  m_org[0] = org.x;
  m_org[1] = org.y;
  m_org[2] = org.z;

}

//@fn FB::Vec3i IP_Sphere::find_index(const FB::Vec3f p)
//@brief 点pの属するセルインデクスを求める
//@note Fortran index
FB::Vec3i IP_Sphere::find_index(const FB::Vec3f p)
{
  FB::Vec3f q = (p-org)/pch;

  return FB::Vec3i( ceil(q.x), ceil(q.y), ceil(q.z) );
}

/**
 @fn void IP_Sphere::setup(int* mid, Control* R, REAL_TYPE* G_org)
 @brief 計算領域のセルIDを設定する
 @param mid IDの配列
 @param R Controlクラスのポインタ
 @param G_org グローバルな原点（無次元）
 */
void IP_Sphere::setup(int* mid, Control* R, REAL_TYPE* G_org)
{
  int i,j,k, gd;
  int mid_fluid=1;        /// 流体
  int mid_solid=2;        /// 固体
  int mid_driver=3;       /// ドライバ部
  int mid_driver_face=4;  /// ドライバ流出面
  unsigned m;
  REAL_TYPE x, y, z, dh, len;
  REAL_TYPE ox, oy, oz, Lx, Ly, Lz;
  REAL_TYPE ox_g, oy_g, oz_g;
  
  FB::Vec3f base, b, t;
  REAL_TYPE ph = pch.x;
  REAL_TYPE r;
  REAL_TYPE rs = radius/R->RefLength;
  
  // ノードローカルの無次元値
  ox = R->org[0];
  oy = R->org[1];
  oz = R->org[2];
  Lx = R->Lbx[0];
  Ly = R->Lbx[1];
  Lz = R->Lbx[2];
  dh = R->dh;
  gd = (int)guide;
  ox_g = G_org[0];
  oy_g = G_org[1];
  oz_g = G_org[2];
  
  // 領域の基点
  org.x -= offset;
  G_org[0] = org.x;
  
  
  // 領域中心座標
  ctr.x = -offset;
  ctr.y = 0.0;
  ctr.z = 0.0;
  
  // 球の中心座標
  os.x = 0.0;
  os.y = 0.0;
  os.z = 0.0;
  
  // 球のbbox
  box_min = os - rs;
  box_max = os + rs;
  box_st = find_index(box_min);
  box_ed = find_index(box_max);
  
  /*
   printf("pitch : %f %f %f\n", pch.x, pch.y, pch.z);
   printf("width : %f %f %f\n", wth.x, wth.y, wth.z);
   printf("center: %f %f %f\n", ctr.x, ctr.y, ctr.z);
   printf("radius: %f\n", radius);
   printf("offset: %f\n", offset);
   printf("origin: %f %f %f\n", org.x, org.y, org.z);
   printf("b_min : %f %f %f\n", box_min.x, box_min.y, box_min.z);
   printf("b_max : %f %f %f\n", box_max.x, box_max.y, box_max.z);
   printf("index : %d %d %d - %d %d %d\n", box_st.x, box_st.y, box_st.z, box_ed.x, box_ed.y, box_ed.z);
   */
  
  // 媒質設定
  size_t m_nx = (imax+2*guide) * (jmax+2*guide) * (kmax+2*guide);
  
  for (size_t n=0; n<m_nx; n++) { 
    mid[n] = mid_fluid;
  }
  
  // 球内部
  for (k=box_st.z; k<=box_ed.z; k++) { 
    for (j=box_st.y; j<=box_ed.y; j++) {
      for (i=box_st.x; i<=box_ed.x; i++) {
        
        base.assign((float)i-0.5, (float)j-0.5, (float)k-0.5);
        b = org + base*ph - os;
        r = b.length();
        
        if ( r <= rs ) {
          m = FBUtility::getFindexS3D(size, guide, i, j, k);
          mid[m] = mid_solid;
        }
      }
    }
  }

  
  // driver設定 iff ドライバ長が正の場合
  if ( drv_length < 0.0 ) return;
  
  // lengthは有次元値
  len = ox_g + (drv_length)/R->RefLength; // グローバルな無次元位置
  
  // ドライバ部分　X-面からドライバ長さより小さい領域
  if ( drv_length > 0.0 ) {
    for (k=1; k<=(int)kmax; k++) {
      for (j=1; j<=(int)jmax; j++) {
        for (i=1; i<=(int)imax; i++) {
          m = FBUtility::getFindexS3D(size, guide, i, j, k);
          x = ox + 0.5*dh + dh*(i-1);
          if ( x < len ) mid[m] = mid_driver;
        }
      }
    }  
  }
  
  // ドライバの下流面にIDを設定
  if ( drv_length > 0.0 ) {
    unsigned m1;
    for (k=1; k<=(int)kmax; k++) {
      for (j=1; j<=(int)jmax; j++) {
        for (i=1; i<=(int)imax; i++) {
          m = FBUtility::getFindexS3D(size, guide, i,   j, k);
          m1= FBUtility::getFindexS3D(size, guide, i+1, j, k);
          if ( (mid[m] == mid_driver) && (mid[m1] == mid_fluid) ) {
            mid[m] = mid_driver_face;
          }
        }
      }
    }    
  }
  
  // ステップ部分を上書き
  for (k=1; k<=(int)kmax; k++) {
    for (j=1; j<=(int)jmax; j++) {
      for (i=1; i<=(int)imax; i++) {
        m = FBUtility::getFindexS3D(size, guide, i, j, k);
        x = ox + 0.5*dh + dh*(i-1);
        y = oy + 0.5*dh + dh*(j-1);
        if ( x < len ) {
          mid[m] = mid_solid;
        }
      }
    }
  }
}

/**
 @fn void IP_Sphere::setup_cut(int* mid, Control* R, REAL_TYPE* G_org, float* cut)
 @brief 計算領域のセルIDとカット情報を設定する
 @param mid IDの配列
 @param R Controlクラスのポインタ
 @param G_org グローバルな原点（無次元）
 @param cut 
 */
void IP_Sphere::setup_cut(int* mid, Control* R, REAL_TYPE* G_org, float* cut)
{
  int i,j,k, gd;
  int mid_fluid=1;        /// 流体
  int mid_solid=2;        /// 固体
  int mid_driver=3;       /// ドライバ部
  int mid_driver_face=4;  /// ドライバ流出面
  unsigned m;
  REAL_TYPE x, y, z, dh, len;
  REAL_TYPE ox, oy, oz, Lx, Ly, Lz;
  REAL_TYPE ox_g, oy_g, oz_g;
  
  FB::Vec3f base, b, t;
  REAL_TYPE ph = pch.x;
  REAL_TYPE r;
  REAL_TYPE rs = radius/R->RefLength;
  
  // ノードローカルの無次元値
  ox = R->org[0];
  oy = R->org[1];
  oz = R->org[2];
  Lx = R->Lbx[0];
  Ly = R->Lbx[1];
  Lz = R->Lbx[2];
  dh = R->dh;
  gd = (int)guide;
  ox_g = G_org[0];
  oy_g = G_org[1];
  oz_g = G_org[2];

  // 領域の基点
  org.x -= offset;
  G_org[0] = org.x;
  
  
  // 領域中心座標
  ctr.x = -offset;
  ctr.y = 0.0;
  ctr.z = 0.0;
  
  // 球の中心座標
  os.x = 0.0;
  os.y = 0.0;
  os.z = 0.0;
  
  // 球のbbox
  box_min = os - rs;
  box_max = os + rs;
  box_st = find_index(box_min);
  box_ed = find_index(box_max);

  /*
  printf("pitch : %f %f %f\n", pch.x, pch.y, pch.z);
  printf("width : %f %f %f\n", wth.x, wth.y, wth.z);
  printf("center: %f %f %f\n", ctr.x, ctr.y, ctr.z);
  printf("radius: %f\n", radius);
  printf("offset: %f\n", offset);
  printf("origin: %f %f %f\n", org.x, org.y, org.z);
  printf("b_min : %f %f %f\n", box_min.x, box_min.y, box_min.z);
  printf("b_max : %f %f %f\n", box_max.x, box_max.y, box_max.z);
  printf("index : %d %d %d - %d %d %d\n", box_st.x, box_st.y, box_st.z, box_ed.x, box_ed.y, box_ed.z);
  */
  
  // 媒質設定
  size_t m_nx = (imax+2*guide) * (jmax+2*guide) * (kmax+2*guide);
  
  for (size_t n=0; n<m_nx; n++) { 
    mid[n] = mid_fluid;
  }
  
  // 球内部
  for (k=box_st.z; k<=box_ed.z; k++) { 
    for (j=box_st.y; j<=box_ed.y; j++) {
      for (i=box_st.x; i<=box_ed.x; i++) {
       
        base.assign((float)i-0.5, (float)j-0.5, (float)k-0.5);
        b = org + base*ph - os;
        r = b.length();
        
        if ( r <= rs ) {
          m = FBUtility::getFindexS3D(size, guide, i, j, k);
          mid[m] = mid_solid;
        }
      }
    }
  }
  
  // カット情報
  FB::Vec3f p[7];
  float lb[7], s, r_min=10.0, r_max=0.0;
  
  for (k=box_st.z; k<=box_ed.z; k++) { 
    for (j=box_st.y; j<=box_ed.y; j++) {
      for (i=box_st.x; i<=box_ed.x; i++) {
        
        base.assign((float)i-0.5, (float)j-0.5, (float)k-0.5);
        b = org + base*ph - os;
        
        p[0].assign(b.x   , b.y   , b.z   ); // p
        p[1].assign(b.x-ph, b.y   , b.z   ); // w 
        p[2].assign(b.x+ph, b.y   , b.z   ); // e
        p[3].assign(b.x   , b.y-ph, b.z   ); // s
        p[4].assign(b.x   , b.y+ph, b.z   ); // n
        p[5].assign(b.x   , b.y   , b.z-ph); // b
        p[6].assign(b.x   , b.y   , b.z+ph); // t
        
        for (int l=0; l<7; l++) {
          lb[l] = ( p[l].length() <= rs ) ? -1.0 : 1.0;
        }
        
        // cut test
        for (int l=1; l<=6; l++) {
          if ( lb[0]*lb[l] < 0.0 ) {
            s = cut_line(p[0], l, rs, ph);
            //printf("%8.5f %8.5f %8.5f > %3.1f : %8.5f %8.5f %8.5f > %3.1f : %d : %f\n", b.x, b.y, b.z, lb[0], p[l].x, p[l].y, p[l].z, lb[l], l, s);
            if ( r_min > s ) r_min = s;
            if ( r_max < s ) r_max = s;
          }
        }
      }
    }
  }
  
  printf("\nCut info. for Sphere\n");
  printf("\tmin. cut = %f\n", r_min);
  printf("\tmax. cut = %f\n", r_max);


  
  // driver設定 iff ドライバ長が正の場合
  if ( drv_length < 0.0 ) return;
  
  // lengthは有次元値
  len = ox_g + (drv_length)/R->RefLength; // グローバルな無次元位置
  
  // ドライバ部分　X-面からドライバ長さより小さい領域
  if ( drv_length > 0.0 ) {
    for (k=1; k<=(int)kmax; k++) {
      for (j=1; j<=(int)jmax; j++) {
        for (i=1; i<=(int)imax; i++) {
          m = FBUtility::getFindexS3D(size, guide, i, j, k);
          x = ox + 0.5*dh + dh*(i-1);
          if ( x < len ) mid[m] = mid_driver;
        }
      }
    }  
  }
  
  // ドライバの下流面にIDを設定
  if ( drv_length > 0.0 ) {
    unsigned m1;
    for (k=1; k<=(int)kmax; k++) {
      for (j=1; j<=(int)jmax; j++) {
        for (i=1; i<=(int)imax; i++) {
          m = FBUtility::getFindexS3D(size, guide, i,   j, k);
          m1= FBUtility::getFindexS3D(size, guide, i+1, j, k);
          if ( (mid[m] == mid_driver) && (mid[m1] == mid_fluid) ) {
            mid[m] = mid_driver_face;
          }
        }
      }
    }    
  }

  // ステップ部分を上書き
  for (k=1; k<=(int)kmax; k++) {
    for (j=1; j<=(int)jmax; j++) {
      for (i=1; i<=(int)imax; i++) {
        m = FBUtility::getFindexS3D(size, guide, i, j, k);
        x = ox + 0.5*dh + dh*(i-1);
        y = oy + 0.5*dh + dh*(j-1);
        if ( x < len ) {
          mid[m] = mid_solid;
        }
      }
    }
  }
}

/**
 @fn float IP_Sphere::cut_line(const FB::Vec3f p, const int dir, const float r, const float dh)
 @brief 計算領域のセルIDを設定する
 @param p 基点座標
 @param dir テスト方向
 @param r radius
 @param dh 格子幅
 */
float IP_Sphere::cut_line(const FB::Vec3f p, const int dir, const float r, const float dh)
{
  float x, y, z, s1, s2, e1, e2, e3, a, b, c;
  
  // unit vector
  switch (dir) {
    case 1:
      e1 =-1.0;
      e2 = 0.0;
      e3 = 0.0;
      break;
      
    case 2:
      e1 = 1.0;
      e2 = 0.0;
      e3 = 0.0;
      break;
      
    case 3:
      e1 = 0.0;
      e2 =-1.0;
      e3 = 0.0;
      break;
      
    case 4:
      e1 = 0.0;
      e2 = 1.0;
      e3 = 0.0;
      break;
      
    case 5:
      e1 = 0.0;
      e2 = 0.0;
      e3 =-1.0;
      break;
      
    case 6:
      e1 = 0.0;
      e2 = 0.0;
      e3 = 1.0;
      break;
      
    default:
      Exit(0);
      break;
  }
  
  x = p.x;
  y = p.y;
  z = p.z;
  
  a = 1.0;
  b = 2.0 * (x*e1 + y*e2 + z*e3);
  c = x*x + y*y + z*z - r*r;
  
  s1 = 0.5 * (-b + sqrtf(b*b-4.0*a*c) );
  s2 = 0.5 * (-b - sqrtf(b*b-4.0*a*c) );
  
  // [0, 1]の実根を拾う
  int f=0;
  int d=0;
  if ( (s1>=0.0) && (s1<dh) ) { f++; d=1; }
  if ( (s2>=0.0) && (s2<dh) ) { f++; d=(d==0)?2:3; }
  
  if ( f == 0 ) {
    printf("Solution error %f %f : (%f, %f, %f) %d\n" ,s1, s2, x, y, z, dir);
    float g1 = (x+s1*e1)*(x+s1*e1) + (y+s1*e2)*(y+s1*e2) + (z+s1*e3)*(z+s1*e3);
    float g2 = (x+s2*e1)*(x+s2*e1) + (y+s2*e2)*(y+s2*e2) + (z+s2*e3)*(z+s2*e3);
    printf("g1= %f  g2= %f\n", g1, g2);
    Exit(0);
  }
  
  if ( d == 1 ) {
    return s1/dh;
  }
  else if (d == 2 ) {
    return s2/dh;
  }
  else {
    return (s1<s2) ? s1/dh : s2/dh;
  }
}
