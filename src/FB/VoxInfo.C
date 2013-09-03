//##################################################################################
//
// Flow Base class
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
 * @file   VoxInfo.C
 * @brief  FlowBase VoxInfo class
 * @author kero
 */

#include <set>
#include <algorithm>
#include <map>
#include "VoxInfo.h"


// #################################################################
// 外部境界に接するガイドセルのmid[]にIDをエンコードする（内部周期境界の場合）
void VoxInfo::adjMediumPrdc_Inner(int* mid, CompoList* cmp)
{
  int st[3], ed[3];
  
  for (int n=1; n<=NoCompo; n++)
  {
    if ( !cmp[n].isKindMedium() )
    {
      cmp[n].getBbox(st, ed);
      int dir = cmp[n].getPeriodicDir();
      
      if ( cmp[n].getType() == PERIODIC )
      {
        if ( numProc > 1 )
        {
          Hostonly_ printf("Error : Inner Periodic condition is limited to use for serial execution on a temporary\n.");
          Exit(0);
        }
        copyID_Prdc_Inner(mid, st, ed, n, dir);
      }
    }
  }
  
}


// #################################################################
// コンポーネントのbid情報から近似的な面積と法線を計算する
void VoxInfo::calCompoArea(const float dhd, const int* bid, CompoList* cmp, FILE* fp)
{
  Hostonly_
  {
    printf("\n---------------------------------------------------------------------------\n\n");
    printf("\t>> Comparison of Component normal and area\n\n");
    printf("\t  No.  Normal                 Original  /                      Approximated  : Area [m*m] Org./  Approx.     Err.[%%]\n");
    
    fprintf(fp, "\n---------------------------------------------------------------------------\n\n");
    fprintf(fp, "\t>> Comparison of Component normal and area\n\n");
  }
  
  int st[3], ed[3]; // コンポーネントのBbox
  int ss[3];        // 各軸方向の面積
  float nv[3];      // 近似計算した法線

  for (int n=1; n<=NoCompo; n++)
  {
    // 媒質のときにはスキップ
    if ( cmp[n].isKindMedium() ) continue;
    
    float vec[3] = { (float)cmp[n].nv[0], (float)cmp[n].nv[1], (float)cmp[n].nv[2] };
    
    // コンポーネントのBbox情報
    cmp[n].getBbox(st, ed);
    
    
    // 各方向の張るスクリーンに要素を投影し，面積を求める
    projectCompo(st, ed, n, bid, ss);
    
    
    // 無次元面積
    float aa = sqrt( ss[0]*ss[0] + ss[1]*ss[1] + ss[2]*ss[2] );
    
    if ( aa != 0.0 )
    {
      nv[0] = (float)ss[0] / aa;
      nv[1] = (float)ss[1] / aa;
      nv[2] = (float)ss[2] / aa;
    }
    else
    {
      nv[0] = nv[1] = nv[2] = 0.0;
    }
    
    float ap = aa * dhd * dhd; // dhd は有次元値，近似的な面積
    float ao = (float)cmp[n].area;
    
    Hostonly_
    {
      printf("\t %3d (%10.3e %10.3e %10.3e) / (%10.3e %10.3e %10.3e) :  %12.5e  / %12.5e  %6.2f\n",
             n, vec[0], vec[1], vec[2], nv[0], nv[1], nv[2], ao, ap, fabs(ap-ao)/ao*100.0);
      fprintf(fp,"\t %3d (%10.3e %10.3e %10.3e) / (%10.3e %10.3e %10.3e) :  %12.5e  / %12.5e  %6.2f\n",
              n, vec[0], vec[1], vec[2], nv[0], nv[1], nv[2], ao, ap, fabs(ap-ao)/ao*100.0);
    }
    
    cmp[n].area = (REAL_TYPE)ap;
    
    // 法線はオリジナルを使う
    //cmp[n].nv[0]= (REAL_TYPE)nv[0];
    //cmp[n].nv[1]= (REAL_TYPE)nv[1];
    //cmp[n].nv[2]= (REAL_TYPE)nv[2];

  }
  
}



// #################################################################
// dst[]にsrc[]のstate, activeビットの情報をコピーする
void VoxInfo::copyBCIbase(int* dst, const int* src)
{
  size_t nx = (size_t)(size[0]+2*guide) * (size_t)(size[1]+2*guide) * (size_t)(size[2]+2*guide);
  
  // 30, 31ビットのみコピー
  for (size_t m=0; m<nx; m++)
  {
    dst[m] = src[m] & 0xc0000000;
  }  
}


// #################################################################
/**
 * @brief 外部境界に接するガイドセルのmid[]にIDを内部周期境界からコピーする
 * @param [in,out] mid   ID配列のデータクラス
 * @param [in]     m_st  コンポーネントのbbox始点
 * @param [in]     m_ed  コンポーネントのbbox終点
 * @param [in]     m_id  対象のID
 * @param [in]     m_dir ドライバの方向
 */
void VoxInfo::copyID_Prdc_Inner(int* mid, const int* m_st, const int* m_ed, const int m_id, const int m_dir)
{
  size_t m0, m1, m2;
  
  int st[3] = {m_st[0], m_st[1], m_st[2]};
  int ed[3] = {m_ed[0], m_ed[1], m_ed[2]};
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int dir = m_dir;
  int id  = m_id;
  int gd  = guide;
  
  switch (dir) {
    case X_MINUS:
      if ( nID[dir] < 0 )
      {
        int i = st[0];
        for (int k=st[2]; k<=ed[2]; k++) {
          for (int j=st[1]; j<=ed[1]; j++) {
            m1 = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
            if ( mid[m1] == id )
            {
              for (int ii=1-gd; ii<=0; ii++) {
                m0 = _F_IDX_S3D(ii,   j, k, ix, jx, kx, gd);
                m2 = _F_IDX_S3D(ii+i, j, k, ix, jx, kx, gd);
                mid[m0] = mid[m2];
              }
            }
          }
        }
      }
      break;
      
    case X_PLUS:
      if ( nID[dir] < 0 )
      {
        int i = st[0];
        for (int k=st[2]; k<=ed[2]; k++) {
          for (int j=st[1]; j<=ed[1]; j++) {
            m1 = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
            if ( mid[m1] == id )
            {
              for (int ii=ix+1; ii<=ix+gd; ii++) {
                m0 = _F_IDX_S3D(ii,        j, k, ix, jx, kx, gd);
                m2 = _F_IDX_S3D(ii+i-ix-1, j, k, ix, jx, kx, gd);
                mid[m0] = mid[m2];
              }
            }
          }
        }
      }
      break;
      
    case Y_MINUS:
      if ( nID[dir] < 0 )
      {
        int j = st[1];
        for (int k=st[2]; k<=ed[2]; k++) {
          for (int i=st[0]; i<=ed[0]; i++) {
            m1 = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
            if ( mid[m1] == id )
            {
              for (int jj=1-gd; jj<=0; jj++) {
                m0 = _F_IDX_S3D(i, jj,   k, ix, jx, kx, gd);
                m2 = _F_IDX_S3D(i, jj+j, k, ix, jx, kx, gd);
                mid[m0] = mid[m2];
              }
            }
          }
        }
      }
      break;
      
    case Y_PLUS:
      if ( nID[dir] < 0 )
      {
        int j = st[1];
        for (int k=st[2]; k<=ed[2]; k++) {
          for (int i=st[0]; i<=ed[0]; i++) {
            m1 = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
            if ( mid[m1] == id )
            {
              for (int jj=jx+1; jj<=jx+gd; jj++) {
                m0 = _F_IDX_S3D(i, jj,        k, ix, jx, kx, gd);
                m2 = _F_IDX_S3D(i, jj+j-jx-1, k, ix, jx, kx, gd);
                mid[m0] = mid[m2];
              }
            }
          }
        }
      }
      break;
      
    case Z_MINUS:
      if ( nID[dir] < 0 )
      {
        int k = st[2];
        for (int j=st[1]; j<=ed[1]; j++) {
          for (int i=st[0]; i<=ed[0]; i++) {
            m1 = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
            if ( mid[m1] == id )
            {
              for (int kk=1-gd; kk<=0; kk++) {
                m0 = _F_IDX_S3D(i, j, kk,   ix, jx, kx, gd);
                m2 = _F_IDX_S3D(i, j, kk+k, ix, jx, kx, gd);
                mid[m0] = mid[m2];
              }
            }
          }
        }
      }
      break;
      
    case Z_PLUS:
      if ( nID[dir] < 0 )
      {
        int k = st[2];
        for (int j=st[1]; j<=ed[1]; j++) {
          for (int i=st[0]; i<=ed[0]; i++) {
            m1 = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
            if ( mid[m1] == id )
            {
              for (int kk=kx+1; kk<=kx+gd; kk++) {
                m0 = _F_IDX_S3D(i, j, kk,        ix, jx, kx, gd);
                m2 = _F_IDX_S3D(i, j, kk+k-kx-1, ix, jx, kx, gd);
                mid[m0] = mid[m2];
              }
            }
          }
        }
      }
      break;
  }  
  
}


// #################################################################
// mid[]内にあるm_idのセルを数える
// painted : ID=0以外でペイント済みを求める(true), m_idのセルをカウント(false)
unsigned long VoxInfo::countCell(const int* mid, bool painted, int m_id)
{
  unsigned long c=0;
  int id = m_id;
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  bool sw = painted;
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd, id, sw) schedule(static) reduction(+:c)
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        
        if ( sw )
        {
          if ( mid[m] != id ) c++;
        }
        else
        {
          if ( mid[m] == id ) c++;
        }
      }
    }
  }
  
  if ( numProc > 1 )
  {
    unsigned long c_tmp = c;
    if ( paraMngr->Allreduce(&c_tmp, &c, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
  }
  
  return c;
}


// #################################################################
// セルの状態をカウントして，その個数をLcell, Gcellに保持する
void VoxInfo::countCellState(unsigned long& Lcell, unsigned long& Gcell, const int* bx, const int state)
{
  unsigned long cell=0;    // local
  unsigned long g_cell=0;  // global 
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  int st = state;
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd, st) schedule(static) reduction(+:cell)
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);

        if ( st == SOLID)
        {
          if (!IS_FLUID(bx[m]) ) cell++;  //  IS_Fluid() > 0=SOLID, 1=FLUID
        }
        else
        {
          if ( IS_FLUID(bx[m]) ) cell++;
        }
      }
    }
  }
  
  g_cell = cell;
  Lcell  = cell;
  
  if ( numProc > 1 )
  {
    unsigned long c_tmp = g_cell;
    if ( paraMngr->Allreduce(&c_tmp, &g_cell, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
  }
  
  Gcell = g_cell;

}


// #################################################################
// 計算領域の外部境界でガイドセル1層と内側の両方が流体セル数の場合にカウントする
void VoxInfo::countOpenAreaOfDomain(const int* bx, REAL_TYPE* OpenArea)
{
  unsigned g;
  unsigned m_area[NOFACE];
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  for (int i=0; i<NOFACE; i++)
  {
    m_area[i] = 0.0;
  }
  
  // described in Fortran index
  
  // X_MINUS
  if( nID[X_MINUS] < 0 )
  {
    g=0;
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static) reduction(+:g)
    for (int k=1; k<=kx; k++) {
      for (int j=1; j<=jx; j++) {
        size_t m0 = _F_IDX_S3D(0, j, k, ix, jx, kx, gd);
        size_t m1 = _F_IDX_S3D(1, j, k, ix, jx, kx, gd);
        if ( (  FLUID == IS_FLUID(bx[m0]))
            && (FLUID == IS_FLUID(bx[m1])) ) g++;
      }
    }
    m_area[X_MINUS] = g;
  }
  
  // X_PLUS
  if( nID[X_PLUS] < 0 )
  {
    g=0;
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static) reduction(+:g)
    for (int k=1; k<=kx; k++) {
      for (int j=1; j<=jx; j++) {
        size_t m0 = _F_IDX_S3D(ix+1, j, k, ix, jx, kx, gd);
        size_t m1 = _F_IDX_S3D(ix,   j, k, ix, jx, kx, gd);
        if ( (  FLUID == IS_FLUID(bx[m0]))
            && (FLUID == IS_FLUID(bx[m1])) ) g++;
      }
    }
    m_area[X_PLUS] = g;
  }
  
  // Y_MINUS
  if( nID[Y_MINUS] < 0 )
  {
    g=0;
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static) reduction(+:g)
    for (int k=1; k<=kx; k++) {
      for (int i=1; i<=ix; i++) {
        size_t m0 = _F_IDX_S3D(i, 0, k, ix, jx, kx, gd);
        size_t m1 = _F_IDX_S3D(i, 1, k, ix, jx, kx, gd);
        if ( (  FLUID == IS_FLUID(bx[m0]))
            && (FLUID == IS_FLUID(bx[m1])) ) g++;
      }
    }
    m_area[Y_MINUS] = g;
  }
  
  // Y_PLUS
  if( nID[Y_PLUS] < 0 )
  {
    g=0;
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static) reduction(+:g)
    for (int k=1; k<=kx; k++) {
      for (int i=1; i<=ix; i++) {
        size_t m0 = _F_IDX_S3D(i, jx+1, k, ix, jx, kx, gd);
        size_t m1 = _F_IDX_S3D(i, jx,   k, ix, jx, kx, gd);
        if ( (  FLUID == IS_FLUID(bx[m0]))
            && (FLUID == IS_FLUID(bx[m1])) ) g++;
      }
    }
    m_area[Y_PLUS] = g;
  }
  
  // Z_MINUS
  if( nID[Z_MINUS] < 0 )
  {
    g=0;
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static) reduction(+:g)
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        size_t m0 = _F_IDX_S3D(i, j, 0, ix, jx, kx, gd);
        size_t m1 = _F_IDX_S3D(i, j, 1, ix, jx, kx, gd);
        if ( (  FLUID == IS_FLUID(bx[m0]))
            && (FLUID == IS_FLUID(bx[m1])) ) g++;
      }
    }
    m_area[Z_MINUS] = g;
  }
  
  // Z_PLUS
  if( nID[Z_PLUS] < 0 )
  {
    g=0;
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static) reduction(+:g)
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        size_t m0 = _F_IDX_S3D(i, j, kx,   ix, jx, kx, gd);
        size_t m1 = _F_IDX_S3D(i, j, kx+1, ix, jx, kx, gd);
        if ( (  FLUID == IS_FLUID(bx[m0]))
            && (FLUID == IS_FLUID(bx[m1])) ) g++;
      }
    }
    m_area[Z_PLUS] = g;
  }
  
  // 外部境界面の面素がunsignedの値域を超えることはないと仮定
  if ( numProc > 1 )
  {
    unsigned tmp[NOFACE];
    for (int i=0; i<NOFACE; i++) tmp[i] = m_area[i];
    if ( paraMngr->Allreduce(tmp, m_area, NOFACE, MPI_SUM) != CPM_SUCCESS ) Exit(0);
  }
  
  for (int i=0; i<NOFACE; i++) OpenArea[i] = (REAL_TYPE)m_area[i];
}


// #################################################################
/**
 * @brief 外部境界面の有効セル数をカウントする
 * @param [in] face 外部境界面番号
 * @param [in] bv   BCindex V
 * @param [in] typ  BC class
 * @note 外部境界面の両側のセルがFluidのときのみカウント
 */
unsigned long VoxInfo::countValidCellOBC(const int face, const int* bv, const int typ)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int sz[3] = {ix, jx, kx};
  int gd = guide;
  
  unsigned long g=0;
  
  if( nID[face] < 0 ) // 外部境界をもつノードのみ
  {
    switch (face)
    {
      case X_MINUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static) reduction(+:g)
        for (int k=1; k<=kx; k++) {
          for (int j=1; j<=jx; j++) {
            size_t m1 = _F_IDX_S3D(1, j, k, ix, jx, kx, gd);
            size_t m2 = _F_IDX_S3D(0, j, k, ix, jx, kx, gd);
            int s1 = bv[m1];
            int s2 = bv[m2];
            if ( IS_FLUID(s1) && IS_FLUID(s2) ) g++;
          }
        }
        break;
        
      case X_PLUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static) reduction(+:g)
        for (int k=1; k<=kx; k++) {
          for (int j=1; j<=jx; j++) {
            size_t m1 = _F_IDX_S3D(ix,   j, k, ix, jx, kx, gd);
            size_t m2 = _F_IDX_S3D(ix+1, j, k, ix, jx, kx, gd);
            int s1 = bv[m1];
            int s2 = bv[m2];
            if ( IS_FLUID(s1) && IS_FLUID(s2) ) g++;
          }
        }
        break;
        
      case Y_MINUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static) reduction(+:g)
        for (int k=1; k<=kx; k++) {
          for (int i=1; i<=ix; i++) {
            size_t m1 = _F_IDX_S3D(i, 1, k, ix, jx, kx, gd);
            size_t m2 = _F_IDX_S3D(i, 0, k, ix, jx, kx, gd);
            int s1 = bv[m1];
            int s2 = bv[m2];
            if ( IS_FLUID(s1) && IS_FLUID(s2) ) g++;
          }
        }
        break;
        
      case Y_PLUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static) reduction(+:g)
        for (int k=1; k<=kx; k++) {
          for (int i=1; i<=ix; i++) {
            size_t m1 = _F_IDX_S3D(i, jx,   k, ix, jx, kx, gd);
            size_t m2 = _F_IDX_S3D(i, jx+1, k, ix, jx, kx, gd);
            int s1 = bv[m1];
            int s2 = bv[m2];
            if ( IS_FLUID(s1) && IS_FLUID(s2) ) g++;
          }
        }
        break;
        
      case Z_MINUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static) reduction(+:g)
        for (int j=1; j<=jx; j++) {
          for (int i=1; i<=ix; i++) {
            size_t m1 = _F_IDX_S3D(i, j, 1, ix, jx, kx, gd);
            size_t m2 = _F_IDX_S3D(i, j, 0, ix, jx, kx, gd);
            int s1 = bv[m1];
            int s2 = bv[m2];
            if ( IS_FLUID(s1) && IS_FLUID(s2) ) g++;
          }
        }
        break;
        
      case Z_PLUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static) reduction(+:g)
        for (int j=1; j<=jx; j++) {
          for (int i=1; i<=ix; i++) {
            size_t m1 = _F_IDX_S3D(i, j, kx,   ix, jx, kx, gd);
            size_t m2 = _F_IDX_S3D(i, j, kx+1, ix, jx, kx, gd);
            int s1 = bv[m1];
            int s2 = bv[m2];
            if ( IS_FLUID(s1) && IS_FLUID(s2) ) g++;
          }
        }
        break;
    } // end of switch
  }

  
  if ( numProc > 1 )
  {
    unsigned long tmp = g;
    if ( paraMngr->Allreduce(&tmp, &g, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
  }
  
  return g;
}



// #################################################################
// BCindexIDを表示する（デバッグ用）
void VoxInfo::dbg_chkBCIndexD(const int* bcd, const char* fname)
{
  int s;
  size_t m;
  FILE *fp=NULL;
  
  if ( !(fp=fopen(fname,"w")) )
  {
    Hostonly_ printf("\tSorry, can't open '%s', write failed.\n", fname);
    Exit(0);
  }
  Hostonly_ printf("\n\tCheck BC Index '%s' for check\n", fname);
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  // ガイドセルを含む全領域
  for (int k=1-gd; k<=kx+gd; k++) {
    for (int j=1-gd; j<=jx+gd; j++) {
      for (int i=1-gd; i<=ix+gd; i++) {
        m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        s = bcd[m];
        Hostonly_ fprintf(fp, "[%4d %4d %4d], state=%1d: cmp=%3d vf=%3d force=%3d\n",
                          i, j, k, IS_FLUID(s), 
                          DECODE_CMP(s), 
                          DECODE_VF(s), 
                          (s>>FORCING_BIT)&0x1);
      }
    }
  }
  fflush(fp);
  fclose(fp);
}


// #################################################################
//BCindexH1を表示する（デバッグ用）
void VoxInfo::dbg_chkBCIndexH1(const int* bh, const char* fname)
{
  FILE *fp=NULL;
  
  if ( !(fp=fopen(fname,"w")) )
  {
    Hostonly_ printf("\tSorry, can't open '%s', write failed.\n", fname);
    Exit(0);
  }
  Hostonly_ printf("\n\tCheck BC Index '%s' for check\n", fname);
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  // ガイドセルを含む全領域
  for (int k=1-gd; k<=kx+gd; k++) {
    for (int j=1-gd; j<=jx+gd; j++) {
      for (int i=1-gd; i<=ix+gd; i++) {
        
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        int s = bh[m];
        
        Hostonly_ fprintf(fp, "[%4d %4d %4d], state=%1d: HBC(ewnstb) [%2d %2d %2d %2d %2d %2d] : idx=%d\n",
                          i, j, k, IS_FLUID(s),
                          GET_FACE_BC(s, BC_FACE_E), // 5-bit
                          GET_FACE_BC(s, BC_FACE_W),
                          GET_FACE_BC(s, BC_FACE_N),
                          GET_FACE_BC(s, BC_FACE_S),
                          GET_FACE_BC(s, BC_FACE_T),
                          GET_FACE_BC(s, BC_FACE_B),
                          s);
        
      }
    }
  }
  fflush(fp);
  fclose(fp);
}


// #################################################################
//BCindexH2を表示する（デバッグ用）
void VoxInfo::dbg_chkBCIndexH2(const int* bh, const char* fname)
{
  FILE *fp=NULL;
  
  if ( !(fp=fopen(fname,"w")) )
  {
    Hostonly_ printf("\tSorry, can't open '%s', write failed.\n", fname);
    Exit(0);
  }
  Hostonly_ printf("\n\tCheck BC Index '%s' for check\n", fname);
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  // ガイドセルを含む全領域
  for (int k=1-gd; k<=kx+gd; k++) {
    for (int j=1-gd; j<=jx+gd; j++) {
      for (int i=1-gd; i<=ix+gd; i++) {
        
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        int s = bh[m];
        
        Hostonly_ fprintf(fp, "[%4d %4d %4d], state=%1d: A(ewnstb) [%2d %2d %2d %2d %2d %2d] : G[%2d %2d %2d %2d %2d %2d] : DIAG = %3d : ORDER = %3d\n",
                          i, j, k, IS_FLUID(s),
                          (s>>ADIABATIC_E)&0x1, // 1-bit
                          (s>>ADIABATIC_W)&0x1,
                          (s>>ADIABATIC_N)&0x1,
                          (s>>ADIABATIC_S)&0x1,
                          (s>>ADIABATIC_T)&0x1,
                          (s>>ADIABATIC_B)&0x1,
                          (s>>GMA_E)&0x1,
                          (s>>GMA_W)&0x1,
                          (s>>GMA_N)&0x1,
                          (s>>GMA_S)&0x1,
                          (s>>GMA_T)&0x1,
                          (s>>GMA_B)&0x1,
                          (s>>H_DIAG)&0x7, // 3-bit
                          DECODE_CMP(s)
                          );
      }
    }
  }
  fflush(fp);
  fclose(fp);
}



// #################################################################
// BCindexPを表示する（デバッグ用）
void VoxInfo::dbg_chkBCIndexP(const int* bcd, const int* bcp, const char* fname)
{
  size_t m;
  int q, s, d;
  FILE *fp=NULL;
  
  if ( !(fp=fopen(fname,"w")) ) {
    Hostonly_ printf("\tSorry, can't open '%s', write failed.\n", fname);
    Exit(0);
  }
  
  Hostonly_ printf("\n\tCheck BC Index '%s' for check\n", fname);
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  // ガイドセルを含む全領域
  for (int k=1-gd; k<=kx+gd; k++) {
    for (int j=1-gd; j<=jx+gd; j++) {
      for (int i=1-gd; i<=ix+gd; i++) {
        m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        d = bcd[m];
        s = bcp[m];
        q = DECODE_CMP(d);
        Hostonly_ fprintf(fp, "[%4d %4d %4d], cmp[%3d], state=%1d: D(ewnstb) [%d %d %d %d %d %d] N [%d %d %d %d %d %d] NDAG [%d %d %d %d %d %d] DIAG=%1d :idx=%d\n", 
                          i, j, k, q, IS_FLUID(d), 
                          (s>>BC_D_E)&0x1,
                          (s>>BC_D_W)&0x1,
                          (s>>BC_D_N)&0x1,
                          (s>>BC_D_S)&0x1,
                          (s>>BC_D_T)&0x1,
                          (s>>BC_D_B)&0x1,
                          (s>>BC_N_E)&0x1,
                          (s>>BC_N_W)&0x1,
                          (s>>BC_N_N)&0x1,
                          (s>>BC_N_S)&0x1,
                          (s>>BC_N_T)&0x1,
                          (s>>BC_N_B)&0x1,
                          (s>>BC_NDAG_E)&0x1,
                          (s>>BC_NDAG_W)&0x1,
                          (s>>BC_NDAG_N)&0x1,
                          (s>>BC_NDAG_S)&0x1,
                          (s>>BC_NDAG_T)&0x1,
                          (s>>BC_NDAG_B)&0x1,
                          (s>>BC_DIAG)&0x7,
                          s);
      }
    }
  }
  fflush(fp);
  fclose(fp);
}



// #################################################################
//BCindexを表示する（デバッグ用）
void VoxInfo::dbg_chkBCIndexV(const int* bcv, const char* fname)
{
  size_t m;
  int s;
  FILE *fp=NULL;
  
  if ( !(fp=fopen(fname,"w")) )
  {
    Hostonly_ printf("\tSorry, can't open '%s', write failed.\n", fname);
    Exit(0);
  }
  Hostonly_ printf("\n\tCheck BC Index '%s' for check\n", fname);
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  // ガイドセルを含む全領域
  for (int k=1-gd; k<=kx+gd; k++) {
    for (int j=1-gd; j<=jx+gd; j++) {
      for (int i=1-gd; i<=ix+gd; i++) {
        m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        s = bcv[m];
        Hostonly_ fprintf(fp, "[%4d %4d %4d], state=%1d: VBC(ewnstb) [%2d %2d %2d %2d %2d %2d] : idx=%d\n", 
                          i, j, k, IS_FLUID(s), 
                          GET_FACE_BC(s, BC_FACE_E), // 5-bit
                          GET_FACE_BC(s, BC_FACE_W),
                          GET_FACE_BC(s, BC_FACE_N),
                          GET_FACE_BC(s, BC_FACE_S),
                          GET_FACE_BC(s, BC_FACE_T),
                          GET_FACE_BC(s, BC_FACE_B),
                          s);
      }
    }
  }
  fflush(fp);
  fclose(fp);
}



// #################################################################
/**
 * @brief BCindexにそのセルが計算に有効(active)かどうかをエンコードする
 * @param Lcell[out] ノードローカルの有効セル数
 * @param Gcell[out] グローバルの有効セル数
 * @param bx BCindex ID
 * @param KOS 解くべき方程式の種類 KIND_OF_SOLVER
 * @note
   - IS_FLUID returns true if FLUID
   - 設定は内部領域のみ
 */
void VoxInfo::encActive(unsigned long& Lcell, unsigned long& Gcell, int* bx, const int KOS)
{
  unsigned long c=0;
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  switch ( KOS )
  {
    case FLOW_ONLY:
    case THERMAL_FLOW:
    case THERMAL_FLOW_NATURAL:
      
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static) reduction(+:c)
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          for (int i=1; i<=ix; i++) {
            size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
            int s = bx[m];
            
            if ( IS_FLUID( s ) )
            {
              s = onBit( s, ACTIVE_BIT );
              c++;
            }
            else
            {
              s = offBit( s, ACTIVE_BIT );
            }
            bx[m] = s;
          }
        }
      }
      break;
      
    case SOLID_CONDUCTION:
      
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static) reduction(+:c)
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          for (int i=1; i<=ix; i++) {
            size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
            int s = bx[m];
            
            if ( !IS_FLUID( s ) )
            {
              s = onBit( s, ACTIVE_BIT );
              c++;
            }
            else
            {
              s = offBit( s, ACTIVE_BIT );
            }
            bx[m] = s;
          }
        }
      }
      break;
      
    case CONJUGATE_HEAT_TRANSFER:
      
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static) reduction(+:c)
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          for (int i=1; i<=ix; i++) {
            size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
            bx[m] = onBit( bx[m], ACTIVE_BIT );
            c++;
          }
        }
      }
      break;
  }
  
  Lcell = c;
  
  
  if ( numProc > 1 )
  {
    unsigned long c_tmp = c;
    if ( paraMngr->Allreduce(&c_tmp, &c, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
  }
  
  Gcell = c;
}


// #################################################################
/**
 * @brief 外部境界に接するセルに対称境界面の断熱マスクをセットする
 * @param [in]     face 外部境界面番号
 * @param [in,out] bh2  BCindex H2
 */
void VoxInfo::encAmaskSymtrcBC(const int face, int* bh2)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  switch (face)
  {
    case X_MINUS:
      if( nID[face] < 0 )
      {
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
        for (int k=1; k<=kx; k++) {
          for (int j=1; j<=jx; j++) {
            size_t m = _F_IDX_S3D(1, j, k, ix, jx, kx, gd); // 最外層のID
            bh2[m] = offBit(bh2[m], ADIABATIC_W);           // 断熱ビットを0にする
          }
        }
      }
      break;
      
    case X_PLUS:
      if( nID[face] < 0 )
      {
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
        for (int k=1; k<=kx; k++) {
          for (int j=1; j<=jx; j++) {
            size_t m = _F_IDX_S3D(ix, j, k, ix, jx, kx, gd);
            bh2[m] = offBit(bh2[m], ADIABATIC_E);
          }
        }
      }
      break;
      
    case Y_MINUS:
      if( nID[face] < 0 )
      {
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
        for (int k=1; k<=kx; k++) {
          for (int i=1; i<=ix; i++) {
            size_t m = _F_IDX_S3D(i, 1, k, ix, jx, kx, gd);
            bh2[m] = offBit(bh2[m], ADIABATIC_S);
          }
        }
      }
      break;
      
    case Y_PLUS:
      if( nID[face] < 0 )
      {
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
        for (int k=1; k<=kx; k++) {
          for (int i=1; i<=ix; i++) {
            size_t m = _F_IDX_S3D(i, jx, k, ix, jx, kx, gd);
            bh2[m] = offBit(bh2[m], ADIABATIC_N);
          }
        }
      }
      break;
      
    case Z_MINUS:
      if( nID[face] < 0 )
      {
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
        for (int j=1; j<=jx; j++) {
          for (int i=1; i<=ix; i++) {
            size_t m = _F_IDX_S3D(i, j, 1, ix, jx, kx, gd);
            bh2[m] = offBit(bh2[m], ADIABATIC_B);
          }
        }
      }
      break;
      
    case Z_PLUS:
      if( nID[face] < 0 )
      {
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
        for (int j=1; j<=jx; j++) {
          for (int i=1; i<=ix; i++) {
            size_t m = _F_IDX_S3D(i, j, kx, ix, jx, kx, gd);
            bh2[m] = offBit(bh2[m], ADIABATIC_T);
          }
        }
      }
      break;
  }
}


// #################################################################
/**
 * @brief セルの各面を調べ，境界条件が設定されていれば，ビットをON
 * @param [in,out] bh1 BCindex H1
 * @param [in,out] bh2 BCindex H2
 * @note 対角要素の係数をエンコードする
 */
void VoxInfo::encHbit(int* bh1, int* bh2)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  // 初期化
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
  for (int k=0; k<=kx+1; k++) {
    for (int j=0; j<=jx+1; j++) {
      for (int i=0; i<=ix+1; i++) {
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        bh2[m] |= ( 0x3f << GMA_W ); // 6bitまとめて(1)で初期化
      }
    }
  }
  
  // 境界条件ビットの設定
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        int s1 = bh1[m];
        int s2 = bh2[m];
        
        if ( GET_FACE_BC(s1, BC_FACE_W) != 0 ) s2 = offBit( s2, GMA_W );
        if ( GET_FACE_BC(s1, BC_FACE_E) != 0 ) s2 = offBit( s2, GMA_E );
        if ( GET_FACE_BC(s1, BC_FACE_S) != 0 ) s2 = offBit( s2, GMA_S );
        if ( GET_FACE_BC(s1, BC_FACE_N) != 0 ) s2 = offBit( s2, GMA_N );
        if ( GET_FACE_BC(s1, BC_FACE_B) != 0 ) s2 = offBit( s2, GMA_B );
        if ( GET_FACE_BC(s1, BC_FACE_T) != 0 ) s2 = offBit( s2, GMA_T );
        
        bh2[m] = s2;
      }
    }
  }
  
  // 対角要素の係数のチェックとエンコード
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        int s2= bh2[m];
        
        int s_w = BIT_SHIFT(s2, GMA_W); // 対角要素の係数：0 or 1
        int s_e = BIT_SHIFT(s2, GMA_E);
        int s_s = BIT_SHIFT(s2, GMA_S);
        int s_n = BIT_SHIFT(s2, GMA_N);
        int s_b = BIT_SHIFT(s2, GMA_B);
        int s_t = BIT_SHIFT(s2, GMA_T);
        
        int a_w = BIT_SHIFT(s2, ADIABATIC_W);
        int a_e = BIT_SHIFT(s2, ADIABATIC_E);
        int a_s = BIT_SHIFT(s2, ADIABATIC_S);
        int a_n = BIT_SHIFT(s2, ADIABATIC_N);
        int a_b = BIT_SHIFT(s2, ADIABATIC_B);
        int a_t = BIT_SHIFT(s2, ADIABATIC_T);
        
        int ss = s_w * a_w
               + s_e * a_e
               + s_s * a_s
               + s_n * a_n
               + s_b * a_b
               + s_t * a_t;
        
        bh2[m] = s2 | (ss << H_DIAG);
        
        if ( (ss == 0) && BIT_SHIFT(s2, ACTIVE_BIT) )
        {
          Hostonly_ printf("\tDiagonal element is zero at (%d,%d,%d) : (Gamma:wesnbt)[%1d %1d %1d %1d %1d %1d] (A:wesnbt)[%1d %1d %1d %1d %1d %1d] : Act= %d  State= %d\n",
                           i,j,k,
                           s_w, s_e, s_s, s_n, s_b, s_t, 
                           a_w, a_e, a_s, a_n, a_b, a_t,
                           BIT_SHIFT(s2, ACTIVE_BIT),
                           BIT_SHIFT(s2, STATE_BIT) );
        }
        
      }
    }
  }
  
  // ゼロ割防止のためのダミー係数 >> 全領域
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
  for (int k=1-gd; k<=kx+gd; k++) {
    for (int j=1-gd; j<=jx+gd; j++) {
      for (int i=1-gd; i<=ix+gd; i++) {
        
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        int s2 = bh2[m];
        
        if ( ((s2>>H_DIAG) & 0x7) == 0 )  // 0x7 = 3 bit
        {
          bh2[m] = s2 | (0x1<<H_DIAG);
        }
      }
    }
  }
  
}



// #################################################################
/**
 * @brief mat[]/cmp[]の格納順をbx[]へエンコードする
 * @param [in]     order エンコードする格納順
 * @param [in]     mid   セルID配列
 * @param [in,out] bx    BCindex ID/H2
 * @param [in,out] cmp   CompoList
 * @retval エンコードした個数
 * @note mid[]が指定されたidならば，CompoListのエントリをbx[]エンコードする
 *       対象範囲をサブドメイン内であること．拡大すると，並列時に余計にカウントしてしまう
 */
unsigned long VoxInfo::encodeOrder(const int order, const int* mid, int* bx, CompoList* cmp)
{
  unsigned long g=0;
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  int odr = order;
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd, odr) schedule(static) reduction(+:g)
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        
        if ( mid[m] == odr )
        {
          bx[m] |= odr; // bx[m]の下位6bitにエントリをエンコード  >> ParseBC:sertControlVars()でビット幅をチェック
          g++;
        }
      }
    }
  }
  
  if ( g > 0 )
  {
    cmp->setEnsLocal(ON);
  }
  
  if ( numProc > 1 )
  {
    unsigned long tmp = g;
    if ( paraMngr->Allreduce(&tmp, &g, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
  }

  return g;
}




// #################################################################
/**
 * @brief ディリクレ条件とノイマン条件の排他性をチェックし，反復行列の非対角要素/対角要素の係数をエンコードする
 * @param bx BCindex P
 * @note
   - ディリクレ条件とノイマン条件の排他性のチェック
   - 非対角要素と対角要素の係数をエンコードする
 */
void VoxInfo::encPbit(int* bx)
{
  size_t m, flag;
  int coef;
  int s_e, s_w, s_n, s_s, s_t, s_b, ss;
  int d_e, d_w, d_n, d_s, d_t, d_b;
  int s;
  bool exclusive;
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  // ディリクレ条件とノイマン条件の排他性のチェック
  exclusive = true;
  
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        s = bx[m];
        flag = 0;
        
        // ノイマン条件：値がゼロのとき，BCがセットされている
        s_e = BIT_SHIFT(s, BC_N_E); 
        s_w = BIT_SHIFT(s, BC_N_W);
        s_n = BIT_SHIFT(s, BC_N_N);
        s_s = BIT_SHIFT(s, BC_N_S);
        s_t = BIT_SHIFT(s, BC_N_T);
        s_b = BIT_SHIFT(s, BC_N_B);
        
        // ディリクレ条件：値がゼロのとき，BCがセットされている
        d_e = BIT_SHIFT(s, BC_D_E); 
        d_w = BIT_SHIFT(s, BC_D_W);
        d_n = BIT_SHIFT(s, BC_D_N);
        d_s = BIT_SHIFT(s, BC_D_S);
        d_t = BIT_SHIFT(s, BC_D_T);
        d_b = BIT_SHIFT(s, BC_D_B);
        
        // 非対角要素の係数をエンコード
        if ( (s_e * d_e) == 0 ) s = offBit( s, BC_NDAG_E );
        if ( (s_w * d_w) == 0 ) s = offBit( s, BC_NDAG_W );
        if ( (s_n * d_n) == 0 ) s = offBit( s, BC_NDAG_N );
        if ( (s_s * d_s) == 0 ) s = offBit( s, BC_NDAG_S );
        if ( (s_t * d_t) == 0 ) s = offBit( s, BC_NDAG_T );
        if ( (s_b * d_b) == 0 ) s = offBit( s, BC_NDAG_B );
        
        bx[m] = s;
        
        if ( (s_e==0) && (d_e==0) ) flag++;
        if ( (s_w==0) && (d_w==0) ) flag++;
        if ( (s_n==0) && (d_n==0) ) flag++;
        if ( (s_s==0) && (d_s==0) ) flag++;
        if ( (s_t==0) && (d_t==0) ) flag++;
        if ( (s_b==0) && (d_b==0) ) flag++;
        
        if ( flag != 0)
        {
          Hostonly_ printf("\tDirichlet and Neumann BC are specified on the same face in cell (%d,%d,%d)\n", i,j,k);
          exclusive = false;
        }
      }
    }
  }
  if ( !exclusive ) Exit(0);
  
  // 対角要素の係数のチェックとエンコード
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        s = bx[m];
        
        s_e = BIT_SHIFT(s, BC_N_E); // 0-Neumann / 1-normal
        s_w = BIT_SHIFT(s, BC_N_W);
        s_n = BIT_SHIFT(s, BC_N_N);
        s_s = BIT_SHIFT(s, BC_N_S);
        s_t = BIT_SHIFT(s, BC_N_T);
        s_b = BIT_SHIFT(s, BC_N_B);
        
        ss = s_e + s_w + s_n + s_s + s_t + s_b;
        bx[m] = s | (ss<<BC_DIAG);
        
        if ( (ss == 0) && (TEST_BIT(s,ACTIVE_BIT)) )
        {
          Hostonly_ printf("\tError : Coefficient of diagonal element is zero at (%d,%d,%d) : (wesnbt)[%1d %1d %1d %1d %1d %1d]\n", i,j,k,
                           IS_FLUID(bx[_F_IDX_S3D(i-1, j,   k,   ix, jx, kx, gd)]),
                           IS_FLUID(bx[_F_IDX_S3D(i+1, j,   k,   ix, jx, kx, gd)]),
                           IS_FLUID(bx[_F_IDX_S3D(i,   j-1, k,   ix, jx, kx, gd)]),
                           IS_FLUID(bx[_F_IDX_S3D(i,   j+1, k,   ix, jx, kx, gd)]),
                           IS_FLUID(bx[_F_IDX_S3D(i,   j,   k-1, ix, jx, kx, gd)]),
                           IS_FLUID(bx[_F_IDX_S3D(i,   j,   k+1, ix, jx, kx, gd)]) );
        }
        
      }
    }
  }
  
  // ゼロ割防止のためのダミー係数 >> 全領域
  size_t nx = (size[0]+2*guide) * (size[1]+2*guide) * (size[2]+2*guide);
  
  for (m=0; m<nx; m++)
  {
    s = bx[m];
    if ( ((s>>BC_DIAG) & 0x7) == 0 ) // 0x7 = 3 bit
    {
      bx[m] = s | (0x6<<BC_DIAG);
    }
  }
}



// #################################################################
/**
 * @brief 圧力のディリクレ境界ビットをエンコードする
 * @retval エンコードしたセル数
 * @param [in]     order  cmp[]のエントリ番号
 * @param [in]     id     CellID
 * @param [in]     mid    ボクセル配列
 * @param [in,out] bcd    BCindex ID
 * @param [in,out] bcp    BCindex P
 * @param [in]     deface 面を指定するid
 * @note
   - 対象セルが流体セルの場合，隣接する面にDirichletフラグをエンコードする
   - 同種のBCは1セルに一つだけ
 */
unsigned long VoxInfo::encPbit_D_IBC(const int order, 
                                     const int id, 
                                     const int* mid, 
                                     int* bcd, 
                                     int* bcp, 
                                     const int deface)
{
  unsigned long g=0;
  
  int s_e, s_w, s_n, s_s, s_t, s_b;
  int c_p, c_e, c_w, c_n, c_s, c_t, c_b;
  int s, d, m;
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  int idd = id;
  
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        d = bcd[m];
        s = bcp[m];
        c_p = mid[m];
        
        if ( IS_FLUID( s ) && (c_p == deface) ) // 対象セルが流体セル，かつdefaceである
        {
#include "FindexS3D.h"
          
          c_e = mid[m_e];
          c_w = mid[m_w];
          c_n = mid[m_n];
          c_s = mid[m_s];
          c_t = mid[m_t];
          c_b = mid[m_b];
          
          s_e = bcp[m_e];
          s_w = bcp[m_w];
          s_n = bcp[m_n];
          s_s = bcp[m_s];
          s_t = bcp[m_t];
          s_b = bcp[m_b];
          
          // X_MINUS
          if ( c_w == idd ) // Wセルが指定IDかどうかを先にチェック
          {
            if ( !((nID[X_MINUS] < 0) && (i == 1)) ) // 続いて領域チェック 外部境界面は除外
            {
              d |= order; // dにorderをエンコード
              s = offBit( s, BC_D_W );
              g++;
            }
          }
          
          // X_PLUS
          if ( c_e == idd )
          {
            if ( !((nID[X_PLUS] < 0) && (i == ix)) )
            {
              d |= order;
              s = offBit( s, BC_D_E );
              g++;
            }
          }
          
          // Y_MINUS
          if ( c_s == idd )
          {
            if ( !((nID[Y_MINUS] < 0) && (j == 1)) )
            {
              d |= order;
              s = offBit( s, BC_D_S );
              g++;
            }
          }
          
          // Y_PLUS
          if ( c_n == idd )
          {
            if ( !((nID[Y_PLUS] < 0) && (j == jx)) )
            {
              d |= order;
              s = offBit( s, BC_D_N );
              g++;
            }
          }
          
          // Z_MINUS
          if ( c_b == idd )
          {
            if ( !((nID[Z_MINUS] < 0) && (k == 1)) )
            {
              d |= order;
              s = offBit( s, BC_D_B );
              g++;
            }
          }
          
          // Z_PLUS
          if ( c_t == idd )
          {
            if ( !((nID[Z_PLUS] < 0) && (k == kx)) )
            {
              d |= order;
              s = offBit( s, BC_D_T );
              g++;
            }
          }
          
          bcd[m] = d;
          bcp[m] = s;
        }
      }
    }
  }
  
  if ( numProc > 1 )
  {
    unsigned long tmp = g;
    if ( paraMngr->Allreduce(&tmp, &g, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
  }
  
  return g;
}



// #################################################################
/**
 * @brief 圧力のノイマン境界ビットをエンコードする（バイナリボクセル）
 * @param [in,out] bx BCindex P
 * @retval 固体表面セル数
 * @note
 *  - 流体セルのうち，固体セルに隣接する面のノイマンフラグをゼロにする．ただし，内部領域のみ．
 *  - 固体セルに隣接する流体セルに方向フラグを付与する．全内部領域．
 *  - 収束判定の有効フラグをエンコードする
 */
unsigned long VoxInfo::encPbit_N_Binary(int* bx)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  // ノイマンフラグ
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        size_t m_p = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        int s = bx[m_p];
        
        if ( IS_FLUID( s ) )
        {
#include "FindexS3D.h"
          
          // X_MINUS
          if ( !IS_FLUID(bx[m_w]) )
          {
            s = offBit( s, BC_N_W );
          }
          
          // X_PLUS
          if ( !IS_FLUID(bx[m_e]) )
          {
            s = offBit( s, BC_N_E );
          }
          
          // Y_MINUS
          if ( !IS_FLUID(bx[m_s]) )
          {
            s = offBit( s, BC_N_S );
          }
          
          // Y_PLUS
          if ( !IS_FLUID(bx[m_n]) )
          {
            s = offBit( s, BC_N_N );
          }
          
          // Z_MINUS
          if ( !IS_FLUID(bx[m_b]) )
          {
            s = offBit( s, BC_N_B );
          }
          
          // Z_PLUS
          if ( !IS_FLUID(bx[m_t]) )
          {
            s = offBit( s, BC_N_T );
          }
          
          bx[m_p] = s;
        }
      }
    }
  }

  // wall locationフラグ
  unsigned long c = 0;
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static) reduction(+:c)
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        size_t m_p = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        int s = bx[m_p];
        
        if ( IS_FLUID( s ) )
        {
#include "FindexS3D.h"
          
          // X_MINUS
          if ( !IS_FLUID(bx[m_w]) ) { s = onBit( s, FACING_W ); c++; }
          
          // X_PLUS
          if ( !IS_FLUID(bx[m_e]) ) { s = onBit( s, FACING_E ); c++; }
          
          // Y_MINUS
          if ( !IS_FLUID(bx[m_s]) ) { s = onBit( s, FACING_S ); c++; }
          
          // Y_PLUS
          if ( !IS_FLUID(bx[m_n]) ) { s = onBit( s, FACING_N ); c++; }
          
          // Z_MINUS
          if ( !IS_FLUID(bx[m_b]) ) { s = onBit( s, FACING_B ); c++; }
          
          // Z_PLUS
          if ( !IS_FLUID(bx[m_t]) ) { s = onBit( s, FACING_T ); c++; }
          
          // 収束判定の有効フラグ，計算内部領域の流体セルのみ
          s = onBit(s, VLD_CNVG);
          
          bx[m_p] = s;
        }
      }
    }
  }

  if ( numProc > 1 )
  {
    unsigned long tmp = c;
    if ( paraMngr->Allreduce(&tmp, &c, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
  }
  
  return c;
}


// #################################################################
/**
 * @brief 圧力のノイマン境界ビットをエンコードする（カット）
 * @param [in,out] bx          BCindex P
 * @param [in]     bid         カット点のID情報
 * @param [in]     cut         距離情報
 * @param [in]     convergence カットのあるセルは収束判定をしないオプション（trueの時）
 * @retval 固体表面セル数
 * @note
 *   - 流体セルのうち，固体セルに隣接する面のノイマンフラグをゼロにする．ただし，内部領域のみ．
 *   - 固体セルに隣接する流体セルに方向フラグを付与する．全内部領域．
 *   - 収束判定の有効フラグをカット情報からエンコードする
 */
unsigned long VoxInfo::encPbit_N_Cut(int* bx, const int* bid, const float* cut, const bool convergence)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  // ノイマンフラグ
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        size_t m_p = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        int s = bx[m_p];
        
        if ( IS_FLUID( s ) )
        {
          int qq = bid[m_p];
          
          // 隣接セルの方向に対するカットの有無>> 0ならばカット無し
          int qw = getFaceBID(0, qq);
          int qe = getFaceBID(1, qq);
          int qs = getFaceBID(2, qq);
          int qn = getFaceBID(3, qq);
          int qb = getFaceBID(4, qq);
          int qt = getFaceBID(5, qq);
          
          // X_MINUS
          if (qw != 0)  // 交点があるなら壁面なのでノイマン条件をセット
          {
            s = offBit( s, BC_N_W );
          }
          
          // X_PLUS
          if (qe != 0)
          {
            s = offBit( s, BC_N_E );
          }
          
          // Y_MINUS
          if (qs != 0)
          {
            s = offBit( s, BC_N_S );
          }
          
          // Y_PLUS
          if (qn != 0)
          {
            s = offBit( s, BC_N_N );
          }
          
          // Z_MINUS
          if (qb != 0)
          {
            s = offBit( s, BC_N_B );
          }
          
          // Z_PLUS
          if (qt != 0)
          {
            s = offBit( s, BC_N_T );
          }
          
          bx[m_p] = s;
        }
      }
    }
  }
  
  // wall locationフラグ
  unsigned long c = 0;
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static) reduction(+:c)
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        size_t m_p = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        int s = bx[m_p];
        
        if ( IS_FLUID( s ) )
        {
          int qq = bid[m_p];

          int qw = getFaceBID(0, qq);
          int qe = getFaceBID(1, qq);
          int qs = getFaceBID(2, qq);
          int qn = getFaceBID(3, qq);
          int qb = getFaceBID(4, qq);
          int qt = getFaceBID(5, qq);

          
          if ( qw != 0 )
          {
            s = onBit(s, FACING_W); c++;
          }
          if ( qe != 0 )
          {
            s = onBit(s, FACING_E); c++;
          }
          if ( qs != 0 )
          {
            s = onBit(s, FACING_S); c++;
          }
          if ( qn != 0 )
          {
            s = onBit(s, FACING_N); c++;
          }
          if ( qb != 0 )
          {
            s = onBit(s, FACING_B); c++;
          }
          if ( qt != 0 )
          {
            s = onBit(s, FACING_T); c++;
          }
          
          bx[m_p] = s;
        }
      }
    }
  }
  
  if ( numProc > 1 )
  {
    unsigned long tmp = c;
    if ( paraMngr->Allreduce(&tmp, &c, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
  }
  
  
  // 収束判定の有効フラグ
  unsigned long g=0;
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static) reduction(+:g)
  
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        size_t m_p = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        int s = bx[m_p];
        
        if ( IS_FLUID( s ) )
        {
          size_t m = _F_IDX_S4DEX(0, i, j, k, 6, ix, jx, kx, gd);
          const float* pos = &cut[m];
          
          float q0 = floor(pos[0]);
          float q1 = floor(pos[1]);
          float q2 = floor(pos[2]);
          float q3 = floor(pos[3]);
          float q4 = floor(pos[4]);
          float q5 = floor(pos[5]);
          
          // 全周カットがあるセルは孤立セルとして固体セルへ変更
          if ( (q0+q1+q2+q3+q4+q5) == 0.0 )
          {
            s = offBit(s, VLD_CNVG);    // Out of scope
            s = offBit(s, STATE_BIT );  // Solid
            s = offBit(s, ACTIVE_BIT ); // Inactive
            g++;
          }
          else
          {
            s = onBit(s, VLD_CNVG);
          }
          
          bx[m_p] = s;
        }
      }
    }
  }
  
  
  if ( numProc > 1 )
  {
    unsigned long tmp = g;
    if ( paraMngr->Allreduce(&tmp, &g, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
  }
  
  Hostonly_ printf("\tThe number of cells which are changed to INACTIVE and SOLID because of all faces are cut = %ld\n\n", g);
  
  
  // カットのあるセルの収束判定をしないオプション
  if ( convergence )
  {
    g = 0;

#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static) reduction(+:g)
    
    for (int k=1; k<=kx; k++) {
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          size_t m_p = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
          int s = bx[m_p];
          
          if ( IS_FLUID( s ) )
          {
            size_t m = _F_IDX_S4DEX(0, i, j, k, 6, ix, jx, kx, gd);
            const float* pos = &cut[m];
            
            float q0 = floor(pos[0]);
            float q1 = floor(pos[1]);
            float q2 = floor(pos[2]);
            float q3 = floor(pos[3]);
            float q4 = floor(pos[4]);
            float q5 = floor(pos[5]);
            
            // いずれかのセルがひとつでもカットがある場合
            if ( (q0+q1+q2+q3+q4+q5) != 6.0 )
            {
              s = offBit(s, VLD_CNVG);    // Out of scope  @todo check
              g++;
            }
            
            bx[m_p] = s;
          }
        }
      }
    }
    
    if ( numProc > 1 )
    {
      unsigned long tmp = g;
      if ( paraMngr->Allreduce(&tmp, &g, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
    }
    
    Hostonly_ printf("\tThe number of cells which are excluded to convergence judgement by cut = %ld\n\n", g);
    
  }
  
  return c;
}



// #################################################################
//計算領域内部のコンポーネントのNeumannフラグをbcp[]にエンコードする
unsigned long VoxInfo::encPbit_N_IBC(const int order, 
                                     const int* mid, 
                                     int* bcd, 
                                     int* bcp,
                                     const float* vec,
                                     const int bc_dir)
{
  unsigned long g=0;
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  int odr = order;
  
  FB::Vec3f nv(vec);
  
  // 反対方向のとき符号反転 > 内積が負のときに対象位置となる
  if ( bc_dir == CompoList::opposite_direction )
  {
    nv *= -1.0;
  }
  
  FB::Vec3f e_w(-1.0,  0.0,  0.0);
  FB::Vec3f e_e(+1.0,  0.0,  0.0);
  FB::Vec3f e_s( 0.0, -1.0,  0.0);
  FB::Vec3f e_n( 0.0, +1.0,  0.0);
  FB::Vec3f e_b( 0.0,  0.0, -1.0);
  FB::Vec3f e_t( 0.0,  0.0, +1.0);
  
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd, odr) \
            firstprivate(e_w, e_e, e_s, e_n, e_b, e_t, nv) \
            schedule(static) reduction(+:g)
  
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        int d = bcd[m];
        int s = bcp[m];
        
        if ( IS_FLUID( s ) ) // 対象セルが流体セル
        {

#include "FindexS3D.h"

          // X_MINUS
          if ( (mid[m_w] == odr) && (dot(e_w, nv) < 0.0) )
          {
            d |= odr; // dにエントリをエンコード
            s = offBit( s, BC_N_W );
            g++;
          }
          
          // X_PLUS
          if ( (mid[m_e] == odr) && (dot(e_e, nv) < 0.0) )
          {
            d |= odr;
            s = offBit( s, BC_N_E );
            g++;
          }
          
          // Y_MINUS
          if ( (mid[m_s] == odr) && (dot(e_s, nv) < 0.0) )
          {
            d |= odr;
            s = offBit( s, BC_N_S );
            g++;
          }
          
          // Y_PLUS
          if ( (mid[m_n] == odr) && (dot(e_n, nv) < 0.0) )
          {
            d |= odr;
            s = offBit( s, BC_N_N );
            g++;
          }
          
          // Z_MINUS
          if ( (mid[m_b] == odr) && (dot(e_b, nv) < 0.0) )
          {
            d |= odr;
            s = offBit( s, BC_N_B );
            g++;
          }
          
          // Z_PLUS
          if ( (mid[m_t] == odr) && (dot(e_t, nv) < 0.0) )
          {
            d |= odr;
            s = offBit( s, BC_N_T );
            g++;
          }
          
          bcd[m] = d;
          bcp[m] = s;
        }
      }
    }
  }
  
  if ( numProc > 1 )
  {
    unsigned long tmp = g;
    if ( paraMngr->Allreduce(&tmp, &g, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
  }

  return g;
}


// #################################################################
/**
 * @brief 外部境界に接するセルにおいて，bx[]に圧力境界条件keyに対応するビットフラグを設定する
 * @param [in]     face 外部境界面番号
 * @param [in,out] bx   BCindex P
 * @param [in]     key  Dirichlet or Neumann
 * @param [in]     dir  壁面の場合(true)，方向フラグをON
 * @note
   - 流体セルに対してのみ，1-Normal, 0-BC
   - 固体セルに隣接する面のノイマンフラグをゼロにし，方向フラグを立てる
 */
void VoxInfo::encPbitOBC(const int face, int* bx, const string key, const bool dir)
{
  if( nID[face] >= 0 ) return;
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  

  switch (face)
  {
    case X_MINUS:
      
      if ("Neumann"==key)
      {
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
        for (int k=1; k<=kx; k++) {
          for (int j=1; j<=jx; j++) {
            size_t m = _F_IDX_S3D(1, j, k, ix, jx, kx, gd);
            int s = bx[m];
            if ( IS_FLUID(s) )
            {
              if (dir) s = onBit( s, FACING_W );
              bx[m] = offBit( s, BC_N_W );
            }
          }
        }
      }
      else
      {
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
        for (int k=1; k<=kx; k++) {
          for (int j=1; j<=jx; j++) {
            size_t m = _F_IDX_S3D(1, j, k, ix, jx, kx, gd);
            int s = bx[m];
            if ( IS_FLUID(s) )
            {
              if (dir) s = onBit( s, FACING_W );
              bx[m] = offBit( s, BC_D_W );
            }
          }
        }
      }
      break;
      
    case X_PLUS:
      
      if ("Neumann"==key)
      {
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
        for (int k=1; k<=kx; k++) {
          for (int j=1; j<=jx; j++) {
            size_t m = _F_IDX_S3D(ix, j, k, ix, jx, kx, gd);
            int s = bx[m];
            if ( IS_FLUID(s) )
            {
              if (dir) s = onBit( s, FACING_E );
              bx[m] = offBit( s, BC_N_E );
            }
          }
        }
      }
      else
      {
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
        for (int k=1; k<=kx; k++) {
          for (int j=1; j<=jx; j++) {
            size_t m = _F_IDX_S3D(ix, j, k, ix, jx, kx, gd);
            int s = bx[m];
            if ( IS_FLUID(s) )
            {
              if (dir) s = onBit( s, FACING_E );
              bx[m] = offBit( s, BC_D_E );
            }
          }
        }
      }
      break;
      
    case Y_MINUS:
      
      if ("Neumann"==key)
      {
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
        for (int k=1; k<=kx; k++) {
          for (int i=1; i<=ix; i++) {
            size_t m = _F_IDX_S3D(i, 1, k, ix, jx, kx, gd);
            int s = bx[m];
            if ( IS_FLUID(s) )
            {
              if (dir) s = onBit( s, FACING_S );
              bx[m] = offBit( s, BC_N_S );
            }
          }
        }
      }
      else
      {
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
        for (int k=1; k<=kx; k++) {
          for (int i=1; i<=ix; i++) {
            size_t m = _F_IDX_S3D(i, 1, k, ix, jx, kx, gd);
            int s = bx[m];
            if ( IS_FLUID(s) )
            {
              if (dir) s = onBit( s, FACING_S );
              bx[m] = offBit( s, BC_D_S );
            }
          }
        }
      }
      break;
      
    case Y_PLUS:
      
      if ("Neumann"==key)
      {
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
        for (int k=1; k<=kx; k++) {
          for (int i=1; i<=ix; i++) {
            size_t m = _F_IDX_S3D(i, jx, k, ix, jx, kx, gd);
            int s = bx[m];
            if ( IS_FLUID(s) )
            {
              if (dir) s = onBit( s, FACING_N );
              bx[m] = offBit( s, BC_N_N );
            }
          }
        }
      }
      else
      {
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
        for (int k=1; k<=kx; k++) {
          for (int i=1; i<=ix; i++) {
            size_t m = _F_IDX_S3D(i, jx, k, ix, jx, kx, gd);
            int s = bx[m];
            if ( IS_FLUID(s) )
            {
              if (dir) s = onBit( s, FACING_N );
              bx[m] = offBit( s, BC_D_N );
            }
          }
        }
      }
      break;
      
    case Z_MINUS:
      
      if ("Neumann"==key)
      {
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
        for (int j=1; j<=jx; j++) {
          for (int i=1; i<=ix; i++) {
            size_t m = _F_IDX_S3D(i, j, 1, ix, jx, kx, gd);
            int s = bx[m];
            if ( IS_FLUID(s) )
            {
              if (dir) s = onBit( s, FACING_B );
              bx[m] = offBit( s, BC_N_B );
            }
          }
        }
      }
      else
      {
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
        for (int j=1; j<=jx; j++) {
          for (int i=1; i<=ix; i++) {
            size_t m = _F_IDX_S3D(i, j, 1, ix, jx, kx, gd);
            int s = bx[m];
            if ( IS_FLUID(s) )
            {
              if (dir) s = onBit( s, FACING_B );
              bx[m] = offBit( s, BC_D_B );
            }
          }
        }
      }
      break;
      
    case Z_PLUS:
      
      if ("Neumann"==key)
      {
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
        for (int j=1; j<=jx; j++) {
          for (int i=1; i<=ix; i++) {
            size_t m = _F_IDX_S3D(i, j, kx, ix, jx, kx, gd);
            int s = bx[m];
            if ( IS_FLUID(s) )
            {
              if (dir) s = onBit( s, FACING_T );
              bx[m] = offBit( s, BC_N_T );
            }
          }
        }
      }
      else
      {
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
        for (int j=1; j<=jx; j++) {
          for (int i=1; i<=ix; i++) {
            size_t m = _F_IDX_S3D(i, j, kx, ix, jx, kx, gd);
            int s = bx[m];
            if ( IS_FLUID(s) )
            {
              if (dir) s = onBit( s, FACING_T );
              bx[m] = offBit( s, BC_D_T );
            }
          }
        }
      }
      break;
  } // end of switch
}


// #################################################################
/**
 * @brief 熱境界条件のBCエントリをエンコードする
 * @retval エンコードしたセル数
 * @param [in]     order  CompoListのエントリ
 * @param [in]     bid    カットID
 * @param [in,out] bh1    BCindex H1
 * @param [in,out] bh2    BCindex H2
 * @param [in]     flag   true-非断熱(1), false-断熱(0)
 * @param [in]     target 判定対象セルの状態 (0-SOLID / 1-FLUID, others)
 * @param [in,out] cmp    CompoList
 * @note 計算領域内のF-S/S-F界面を想定．S-S界面の場合，両方のセルにBCが設定される
 */
unsigned long VoxInfo::encQface(const int order,
                                const int* bid,
                                int* bh1,
                                int* bh2,
                                const bool flag,
                                const int target,
                                CompoList* cmp)
{
  unsigned long g=0;
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  int odr= order;
  bool es= flag;
  int state = target;
  
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd, odr, es, state) schedule(static) reduction(+:g)
  
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        
        size_t m  = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        
        int d  = bid[m];
        int s1 = bh1[m];
        int s2 = bh2[m];
        
        bool mode = true; // default stete==ANY_STATEの場合，true
        
        // ターゲットの状態の判定
        if ( state == FLUID )
        {
          mode = IS_FLUID(s2); // ターゲットが流体セルの場合にtrue
        }
        else // SOLID
        {
          mode = !IS_FLUID(s2); // ターゲットが流体セルでない場合にtrue
        }
        
        // 6面のいずれかにカットがあり，ターゲットの状態であるセルが候補
        if ( TEST_BC(d) && mode )
        {
          int b0 = getFaceBID(0, d);
          int b1 = getFaceBID(1, d);
          int b2 = getFaceBID(2, d);
          int b3 = getFaceBID(3, d);
          int b4 = getFaceBID(4, d);
          int b5 = getFaceBID(5, d);
          
          
          // X-
          if ( b0 == odr )
          {
            s1 |= (odr << BC_FACE_W);  // エントリをエンコード
            s2 = ( es == true ) ? onBit(s2, ADIABATIC_W) : offBit( s2, ADIABATIC_W );
            g++;
          }
          
          // X+
          if ( b1 == odr )
          {
            s1 |= (odr << BC_FACE_E);
            s2 = ( es == true ) ? onBit(s2, ADIABATIC_E) : offBit( s2, ADIABATIC_E );
            g++;
          }
          
          // Y-
          if ( b2 == odr )
          {
            s1 |= (odr<< BC_FACE_S);
            s2 = ( es == true ) ? onBit(s2, ADIABATIC_S) : offBit( s2, ADIABATIC_S );
            g++;
          }
          
          // Y+
          if ( b3 == odr )
          {
            s1 |= (odr << BC_FACE_N);
            s2 = ( es == true ) ? onBit(s2, ADIABATIC_N) : offBit( s2, ADIABATIC_N );
            g++;
          }
          
          // Z-
          if ( b4 == odr )
          {
            s1 |= (odr << BC_FACE_B);
            s2 = ( es == true ) ? onBit(s2, ADIABATIC_B) : offBit( s2, ADIABATIC_B );
            g++;
          }
          
          // Z+
          if ( b5 == odr )
          {
            s1 |= (odr << BC_FACE_T);
            s2 = ( es == true ) ? onBit(s2, ADIABATIC_T) : offBit( s2, ADIABATIC_T );
            g++;
          }
          
          bh1[m] = s1;
          bh2[m] = s2;
        }
        
      }
    }
  }
  
  // Local
  if ( g > 0 )
  {
    cmp->setEnsLocal(ON);
  }
  
  
  if ( numProc > 1 )
  {
    unsigned long tmp = g;
    if ( paraMngr->Allreduce(&tmp, &g, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
  }
  
  return g;
}




// #################################################################
/**
 @brief 速度指定，流出型の境界条件の計算に必要な情報をエンコードする
 @param order cmp[]のエントリ番号
 @param id  セルID
 @param mid ボクセル配列
 @param bcd BCindex ID
 @param bh1 BCindex H1
 @param bh2 BCindex H2
 @param deface 面を指定するid
 @note
 - encVbit_IBC()でセル数はカウント済み
 - 流体かつ指定IDであり，defaceが設定されたセルと挟まれる面に対して適用
 - 対象面の断熱ビットを非断熱にする
 - encVbit_IBC()では，最外層を除外しているが，熱では除外しない
 */
void VoxInfo::encQfaceSVO(int order, int id, int* mid, int* bcd, int* bh1, int* bh2, int deface)
{
  int idd;

  int c_p, c_e, c_w, c_n, c_s, c_t, c_b;
  int s1, s2, d;
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  idd = id;
  
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        size_t m_p = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        c_p = mid[m_p];
        
        d  = bcd[m_p];
        s1 = bh1[m_p];
        s2 = bh2[m_p];
        
        if ( (c_p == deface) && IS_FLUID(s1) ) // defaceセルかつ流体の時に，以下を評価
        {
#include "FindexS3D.h"
          
          c_e = mid[m_e];
          c_w = mid[m_w];
          c_n = mid[m_n];
          c_s = mid[m_s];
          c_t = mid[m_t];
          c_b = mid[m_b];
          
          // X-
          if ( c_w == idd ) // m_wセルの属性チェックはしていない
          {
            d |= order; // エントリをエンコード
            s1 |= (order << BC_FACE_W);  // エントリをエンコード
            s2 = onBit(s2, ADIABATIC_W); // 断熱ビットを非断熱(1)
          }
          
          // X+
          if ( c_e == idd )
          {
            d |= order;
            s1 |= (order << BC_FACE_E);
            s2 = onBit(s2, ADIABATIC_E);
          }
          
          // Y-
          if ( c_s == idd )
          {
            d |= order;
            s1 |= (order << BC_FACE_S);
            s2 = onBit(s2, ADIABATIC_S);
          }
          
          // Y+
          if ( c_n == idd )
          {
            d |= order;
            s1 |= (order << BC_FACE_N);
            s2 = onBit(s2, ADIABATIC_N);
          }
          
          // Z-
          if ( c_b == idd )
          {
            d |= order;
            s1 |= (order << BC_FACE_B);
            s2 = onBit(s2, ADIABATIC_B);
          }
          
          // Z+
          if ( c_t == idd )
          {
            d |= order;
            s1 |= (order << BC_FACE_T);
            s2 = onBit(s2, ADIABATIC_T);
          }
          
          bcd[m_p] = d;
          bh1[m_p] = s1;
          bh2[m_p] = s2;
        }
      } // i-loop
    }
  }
  
  // iddセルのSTATE_BITをFLUIDに変更
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        size_t m_p = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        c_p = mid[m_p];
        
        s1  = bh1[m_p];
        
        if ( (c_p == deface) && IS_FLUID(s1) ) // 指定セルかつ流体の時に，以下を評価
        {
#include "FindexS3D.h"
          
          c_e = mid[m_e];
          c_w = mid[m_w];
          c_n = mid[m_n];
          c_s = mid[m_s];
          c_t = mid[m_t];
          c_b = mid[m_b];
          
          // X-
          if ( c_w == idd )
          {
            bh1[m_w] = onBit( bh1[m_w], STATE_BIT );
          }
          
          // X+
          if ( c_e == idd )
          {
            bh1[m_e] = onBit( bh1[m_e], STATE_BIT );
          }
          
          // Y-
          if ( c_s == idd )
          {
            bh1[m_s] = onBit( bh1[m_s], STATE_BIT );
          }
          
          // Y+
          if ( c_n == idd )
          {
            bh1[m_n] = onBit( bh1[m_n], STATE_BIT );
          }
          
          // Z-
          if ( c_b == idd )
          {
            bh1[m_b] = onBit( bh1[m_b], STATE_BIT );
          }
          
          // Z+
          if ( c_t == idd )
          {
            bh1[m_t] = onBit( bh1[m_t], STATE_BIT );
          }
        }
      } // i-loop
    }
  }
}




// #################################################################
/**
 * @brief bv[]にVBCの境界条件に必要な情報をエンコードし，流入流出の場合にbp[]の隣接壁の方向フラグを除く．
 *        境界条件指定キーセルのSTATEを流体に変更する
 * @retval エンコードしたセル数
 * @param [in]     order  cmp[]のエントリ番号
 * @param [in,out] bv     BCindex V
 * @param [in,out] bp     BCindex P
 * @param [in]     bid    カット点ID
 * @param [in]     vec    法線ベクトル
 * @param [in]     bc_dir 境界条件の方向 same_direction=1, opposite_direction=2
 * @param [in,out] cmp    CompoList
 * @note 指定法線とセルのカット方向ベクトルの内積で判断，vspecとoutflowなのでbp[]のVBC_UWDにマスクビットを立てる
 */
unsigned long VoxInfo::encVbitIBC(const int order,
                                  int* bv,
                                  int* bp,
                                  const int* bid,
                                  const float* vec,
                                  const int bc_dir,
                                  CompoList* cmp)
{
  unsigned long g=0;
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  int odr = order;
  
  FB::Vec3f nv(vec);
  
  // 反対方向のとき符号反転 > 内積が負のときに対象位置となる
  if ( bc_dir == CompoList::opposite_direction )
  {
    nv *= -1.0;
  }
  
  FB::Vec3f e_w(-1.0,  0.0,  0.0);
  FB::Vec3f e_e(+1.0,  0.0,  0.0);
  FB::Vec3f e_s( 0.0, -1.0,  0.0);
  FB::Vec3f e_n( 0.0, +1.0,  0.0);
  FB::Vec3f e_b( 0.0,  0.0, -1.0);
  FB::Vec3f e_t( 0.0,  0.0, +1.0);
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd, odr) \
        firstprivate(e_w, e_e, e_s, e_n, e_b, e_t, nv) \
        schedule(static) reduction(+:g)

  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        
        size_t mp = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        int bd = bid[mp];
        
        if ( TEST_BC(bd) ) // 6方向のうちいずれかにカットがある
        {
          int s = bv[mp];
          int q = bp[mp];
          
          if ( IS_FLUID(s) ) // 流体セルがテスト対象
          {
            // 各方向のID
            int id_w = getFaceBID(0, bd);
            int id_e = getFaceBID(1, bd);
            int id_s = getFaceBID(2, bd);
            int id_n = getFaceBID(3, bd);
            int id_b = getFaceBID(4, bd);
            int id_t = getFaceBID(5, bd);
            
            // X-
            if ( (id_w == odr) && (dot(e_w, nv) < 0.0) )
            {
              s |= (odr << BC_FACE_W);
              q = offBit(q, FACING_W);
              q = onBit(q, VBC_UWD);
              g++;
            }
            
            // X+
            if ( (id_e == odr) && (dot(e_e, nv) < 0.0) )
            {
              s |= (odr << BC_FACE_E);
              q = offBit(q, FACING_E);
              q = onBit(q, VBC_UWD);
              g++;
            }
            
            // Y-
            if ( (id_s == odr) && (dot(e_s, nv) < 0.0) )
            {
              s |= (odr << BC_FACE_S);
              q = offBit(q, FACING_S);
              q = onBit(q, VBC_UWD);
              g++;
            }
            
            // Y+
            if ( (id_n == odr) && (dot(e_n, nv) < 0.0) )
            {
              s |= (odr << BC_FACE_N);
              q = offBit(q, FACING_N);
              q = onBit(q, VBC_UWD);
              g++;
            }
            
            // Z-
            if ( (id_b == odr) && (dot(e_b, nv) < 0.0) )
            {
              s |= (odr << BC_FACE_B);
              q = offBit(q, FACING_B);
              q = onBit(q, VBC_UWD);
              g++;
            }
            
            // Z+
            if ( (id_t == odr) && (dot(e_t, nv) < 0.0) )
            {
              s |= (odr << BC_FACE_T);
              q = offBit(q, FACING_T);
              q = onBit(q, VBC_UWD);
              g++;
            }
          } // if fluid
          
          bv[mp] = s;
          bp[mp] = q;
          
        } // if TEST_BC()
        
      } // i-loop
    }
  }
  
  
  // Local
  if ( g > 0 )
  {
    cmp->setEnsLocal(ON);
  }

// ########## check for debug
# if 0
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        
        size_t mp = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        int s = bv[mp];
        
        int m_flag = 0;
        
        if ( GET_FACE_BC(s, BC_FACE_W) == odr ) m_flag = 1;
        if ( GET_FACE_BC(s, BC_FACE_E) == odr ) m_flag = 2;
        if ( GET_FACE_BC(s, BC_FACE_S) == odr ) m_flag = 3;
        if ( GET_FACE_BC(s, BC_FACE_N) == odr ) m_flag = 4;
        if ( GET_FACE_BC(s, BC_FACE_B) == odr ) m_flag = 5;
        if ( GET_FACE_BC(s, BC_FACE_T) == odr ) m_flag = 6;
        
        if ( m_flag != 0 ) printf("%3d %3d %3d >> %d\n", i, j, k, m_flag);
      }
    }
  }
#endif
// ##########
  
  if ( numProc > 1 )
  {
    unsigned long tmp = g;
    if ( paraMngr->Allreduce(&tmp, &g, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
  }
  
  return g;
}



// #################################################################
/**
 * @brief 外部境界に接するセルにおいて，各種速度境界条件に対応する媒質をチェックし，bv[]にビットフラグをセットする
 * @param [in]     face    外部境界面番号
 * @param [in,out] bv      BCindex V
 * @param [in]     key     fluid or solid　指定するBCが要求するガイドセルの状態 >> エラーチェックに使う
 * @param [in]     enc_sw  trueのとき，エンコードする．falseの場合にはガイドセルの状態チェックのみ
 * @param [in]     chk     ガイドセルの状態をチェックするかどうかを指定
 * @param [in]     bp      BCindex P
 * @param [in]     enc_uwd trueのとき，1次精度のスイッチオン
 * @note
  - 外部境界条件の実装には，流束型とディリクレ型の2種類がある．
  - setMediumOnGC()でガイドセル上のIDを指定済み．指定BCとの適合性をチェックする
 */
void VoxInfo::encVbitOBC(const int face, int* bv, const string key, const bool enc_sw, const string chk, int* bp, bool enc_uwd)
{
  if ( nID[face] >= 0 ) return;
  
  int sw, cw;
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  bool enc = enc_sw;
  bool uwd = enc_uwd;
  
  ( !strcasecmp("fluid", key.c_str()) ) ? sw=1 : sw=0;
  ( !strcasecmp("check", chk.c_str()) ) ? cw=1 : cw=0;
  
  switch (face)
  {
    case X_MINUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, sw, cw, enc, uwd) schedule(static)
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          size_t m = _F_IDX_S3D(1, j, k, ix, jx, kx, gd);
          size_t mt= _F_IDX_S3D(0, j, k, ix, jx, kx, gd);
          
          int s = bv[m];
          int z = bv[mt];
          int q = bp[mt];
          
          if ( IS_FLUID(s) )
          {
            // エンコード処理
            if ( enc )
            {
              bv[m]  = s | (OBC_MASK << BC_FACE_W); // OBC_MASK==31 外部境界条件のフラグ
            }
            
            // 外部境界で安定化のため，スキームを1次精度にする
            if ( uwd )
            {
              bp[mt] = onBit(q, VBC_UWD);
            }
            
            // チェック
            if ( cw == 1 )
            {
              if ( sw==1 ) // ガイドセルが流体であることを要求
              {
                if ( !IS_FLUID(z) )
                {
                  Hostonly_ printf("Error : Fluid cell at Outer boundary X- is required\n");
                  Exit(0);
                }
              }
              else // ガイドセルが固体であることを要求
              {
                if ( IS_FLUID(z) )
                {
                  Hostonly_ printf("Error : Solid cell at Outer boundary X- is required\n");
                  Exit(0);
                }
              }
            }
            
          }// IS_FLUID
        }
      }
      break;
      
      
    case X_PLUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, sw, cw, enc, uwd) schedule(static)
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          size_t m = _F_IDX_S3D(ix,   j, k, ix, jx, kx, gd);
          size_t mt= _F_IDX_S3D(ix+1, j, k, ix, jx, kx, gd);
          
          int s = bv[m];
          int z = bv[mt];
          int q = bp[mt];
          
          if ( IS_FLUID(s) )
          {
            // エンコード処理
            if ( enc )
            {
              bv[m]  = s | (OBC_MASK << BC_FACE_E);
            }
            
            // 外部境界で安定化のため，スキームを1次精度にする
            if ( uwd )
            {
              bp[mt] = onBit(q, VBC_UWD);
            }
            
            // チェック
            if ( cw == 1 )
            {
              if ( sw==1 )
              {
                if ( !IS_FLUID(z) )
                {
                  Hostonly_ printf("Error : Fluid cell at Outer boundary X+ is required\n");
                  Exit(0);
                }
              }
              else
              {
                if ( IS_FLUID(z) )
                {
                  Hostonly_ printf("Error : Solid cell at Outer boundary X+ is required\n");
                  Exit(0);
                }
              }
            }
            
          }// IS_FLUID
        }
      }
      break;
      
      
    case Y_MINUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, sw, cw, enc, uwd) schedule(static)
      for (int k=1; k<=kx; k++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, 1, k, ix, jx, kx, gd);
          size_t mt= _F_IDX_S3D(i, 0, k, ix, jx, kx, gd);
          
          int s = bv[m];
          int z = bv[mt];
          int q = bp[mt];
          
          if ( IS_FLUID(s) )
          {
            // エンコード処理
            if ( enc )
            {
              bv[m]  = s | (OBC_MASK << BC_FACE_S);
            }
            
            // 外部境界で安定化のため，スキームを1次精度にする
            if ( uwd )
            {
              bp[mt] = onBit(q, VBC_UWD);
            }
            
            // チェック
            if ( cw == 1 )
            {
              if ( sw==1 )
              {
                if ( !IS_FLUID(z) )
                {
                  Hostonly_ printf("Error : Fluid cell at Outer boundary Y- is required %d %d %d\n", i, 0, k);
                  Exit(0);
                }
              }
              else
              {
                if ( IS_FLUID(z) )
                {
                  Hostonly_ printf("Error : Solid cell at Outer boundary Y- is required\n");
                  Exit(0);
                }
              }
            }
            
          }// IS_FLUID
        }
      }
      break;
      
      
    case Y_PLUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, sw, cw, enc, uwd) schedule(static)
      for (int k=1; k<=kx; k++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, jx,   k, ix, jx, kx, gd);
          size_t mt= _F_IDX_S3D(i, jx+1, k, ix, jx, kx, gd);
          
          int s = bv[m];
          int z = bv[mt];
          int q = bp[mt];
          
          if ( IS_FLUID(s) )
          {
            // エンコード処理
            if ( enc )
            {
              bv[m]  = s | (OBC_MASK << BC_FACE_N);
            }
            
            // 外部境界で安定化のため，スキームを1次精度にする
            if ( uwd )
            {
              bp[mt] = onBit(q, VBC_UWD);
            }
            
            // チェック
            if ( cw == 1 )
            {
              if ( sw==1 )
              {
                if ( !IS_FLUID(z) )
                {
                  Hostonly_ printf("Error : Fluid cell at Outer boundary Y+ is required\n");
                  Exit(0);
                }
              }
              else
              {
                if ( IS_FLUID(z) )
                {
                  Hostonly_ printf("Error : Solid cell at Outer boundary Y+ is required\n");
                  Exit(0);
                }
              }
            }
            
          }// IS_FLUID
        }
      }
      break;
      
      
    case Z_MINUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, sw, cw, enc, uwd) schedule(static)
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, j, 1, ix, jx, kx, gd);
          size_t mt= _F_IDX_S3D(i, j, 0, ix, jx, kx, gd);
          
          int s = bv[m];
          int z = bv[mt];
          int q = bp[mt];
          
          if ( IS_FLUID(s) )
          {
            // エンコード処理
            if ( enc )
            {
              bv[m]  = s | (OBC_MASK << BC_FACE_B);
            }
            
            // 外部境界で安定化のため，スキームを1次精度にする
            if ( uwd )
            {
              bp[mt] = onBit(q, VBC_UWD);
            }
            
            // チェック
            if ( cw == 1 )
            {
              if ( sw==1 )
              {
                if ( !IS_FLUID(z) )
                {
                  Hostonly_ printf("Error : Fluid cell at Outer boundary Z- is required\n");
                  Exit(0);
                }
              }
              else
              {
                if ( IS_FLUID(z) )
                {
                  Hostonly_ printf("Error : Solid cell at Outer boundary Z- is required\n");
                  Exit(0);
                }
              }
            }
            
          }// IS_FLUID
        }
      }
      break;
      
      
    case Z_PLUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, sw, cw, enc, uwd) schedule(static)
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, j, kx,   ix, jx, kx, gd);
          size_t mt= _F_IDX_S3D(i, j, kx+1, ix, jx, kx, gd);
          
          int s = bv[m];
          int z = bv[mt];
          int q = bp[mt];
          
          if ( IS_FLUID(s) )
          {
            // エンコード処理
            if ( enc )
            {
              bv[m]  = s | (OBC_MASK << BC_FACE_T);
            }
            
            // 外部境界で安定化のため，スキームを1次精度にする
            if ( uwd )
            {
              bp[mt] = onBit(q, VBC_UWD);
            }
            
            // チェック
            if ( cw == 1 )
            {
              if ( sw==1 )
              {
                if ( !IS_FLUID(z) )
                {
                  Hostonly_ printf("Error : Fluid cell at Outer boundary Z+ is required\n");
                  Exit(0);
                }
              }
              else
              {
                if ( IS_FLUID(z) )
                {
                  Hostonly_ printf("Error : Solid cell at Outer boundary Z+ is required\n");
                  Exit(0);
                }
              }
            }
            
          }// IS_FLUID
        }
      }
      break;
      
  } // end of switch

}



// #################################################################
// カットID情報に基づく流体媒質のフィルを実行
// Symmetric fillにより反復回数を減少
unsigned long VoxInfo::fillByBid(int* bid, int* mid, float* cut, const int tgt_id, unsigned long& substituted, const int* m_list)
{
  int tg = tgt_id;       ///< FLUID ID
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  unsigned long filled   = 0; ///< 流体IDでペイントされた数
  unsigned long replaced = 0; ///< 固体IDで置換された数

  int list[64];
  for (int i=0; i<64; i++) list[i] = m_list[i];
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd, tg, list) \
schedule(static) reduction(+:filled) reduction(+:replaced)
  
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        
#include "fill_bid.h"
        
      }
    }
  }
  
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd, tg, list) \
schedule(static) reduction(+:filled) reduction(+:replaced)
  
  for (int k=kx; k>=1; k--) {
    for (int j=jx; j>=1; j--) {
      for (int i=ix; i>=1; i--) {
        
#include "fill_bid.h"

      }
    }
  }
  
  if ( numProc > 1 )
  {
    unsigned long tmp = filled;
    if ( paraMngr->Allreduce(&tmp, &filled, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
    
    tmp = replaced;
    if ( paraMngr->Allreduce(&tmp, &replaced, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
  }
  
  substituted = replaced;

  return filled;
}

// #################################################################
// 媒質ID情報に基づく流体媒質のフィルを実行
// Symmetric fillにより反復回数を減少
unsigned long VoxInfo::fillByMid(int* bid, int* mid, float* cut, const int tgt_id, const int* m_list)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  int tg = tgt_id;       ///< FLUID ID
  unsigned long replaced = 0; ///< 固体IDで置換された数
  
  int list[64];
  for (int i=0; i<64; i++) list[i] = m_list[i];
  
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd, tg, list) schedule(static) reduction(+:replaced)
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
       
#include "fill_mid.h"
        
      }
    }
  }
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd, tg, list) schedule(static) reduction(+:replaced)
  
  for (int k=kx; k>=1; k--) {
    for (int j=jx; j>=1; j--) {
      for (int i=ix; i>=1; i--) {
        
#include "fill_mid.h"
        
      }
    }
  }
  
  if ( numProc > 1 )
  {
    unsigned long tmp = replaced;
    if ( paraMngr->Allreduce(&tmp, &replaced, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
  }
  
  return replaced;
}


// #################################################################
// targetセルの周囲の固体最頻値でフィルを実行
unsigned long VoxInfo::fillByModalSolid(int* mid, const int target, const int fluid_id)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  int tg = target;
  int fid = fluid_id;
  unsigned long c = 0; /// painted count
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd, tg, fid) schedule(static) reduction(+:c)
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        
        size_t m_p = _F_IDX_S3D(i  , j  , k  , ix, jx, kx, gd);
        size_t m_e = _F_IDX_S3D(i+1, j,   k,   ix, jx, kx, gd);
        size_t m_w = _F_IDX_S3D(i-1, j,   k,   ix, jx, kx, gd);
        size_t m_n = _F_IDX_S3D(i,   j+1, k,   ix, jx, kx, gd);
        size_t m_s = _F_IDX_S3D(i,   j-1, k,   ix, jx, kx, gd);
        size_t m_t = _F_IDX_S3D(i,   j,   k+1, ix, jx, kx, gd);
        size_t m_b = _F_IDX_S3D(i,   j,   k-1, ix, jx, kx, gd);
        
        int qq = mid[m_p];
        int qw = mid[m_w];
        int qe = mid[m_e];
        int qs = mid[m_s];
        int qn = mid[m_n];
        int qb = mid[m_b];
        int qt = mid[m_t];
        
        // 対象セルがtargetの場合
        if ( qq == tg )
        {
          mid[m_p] = find_mode_id(fid, qw, qe, qs, qn, qb, qt);
          c++;
        }
      }
    }
  }
  
  if ( numProc > 1 )
  {
    unsigned long c_tmp = c;
    if ( paraMngr->Allreduce(&c_tmp, &c, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
  }
  
  return c;
}


// #################################################################
// 内部フィルを実行
unsigned long VoxInfo::fillReplace(int* mid, const int target, const int fill_id)
{
  int sd = fill_id;
  int tg = target;
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  unsigned long c = 0; /// painted count
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd, sd, tg) schedule(static) reduction(+:c)
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        
        // 対象セルがtargetの場合
        if ( mid[m] == tg )
        {
          mid[m] = sd;
          c++;
        }
      }
    }
  }
  
  if ( numProc > 1 )
  {
    unsigned long c_tmp = c;
    if ( paraMngr->Allreduce(&c_tmp, &c, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
  }
  
  return c;
}



// #################################################################
// シード点をペイントする
// ヒントとして与えられた外部境界面に接するセルにおいて，確実に流体セルであるセルをフィルする
// もし，ヒント面に固体候補があれば、ぬれ面はフィルしない
// @attention 外部境界面にカットを設定している場合，フィルされない．>> 外部境界面にカットを設定するのはあとのフェイズ
unsigned long VoxInfo::fillSeed(int* mid, const int face, const int target, const float* cut)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  int tg = target;     ///< FLUID ID
  unsigned long c = 0;
  
  
  switch (face)
  {
    case X_MINUS:
      if ( nID[face] < 0 )
      {
#pragma omp parallel for firstprivate(ix, jx, kx, gd, tg) schedule(static) reduction(+:c)
        for (int k=1; k<=kx; k++) {
          for (int j=1; j<=jx; j++) {
            size_t m = _F_IDX_S3D(1, j, k, ix, jx, kx, gd);
            size_t mp = _F_IDX_S4DEX(0, 1, j, k, 6, ix, jx, kx, gd);
            const float* pos = &cut[mp];
            int flag = 0;
            
            // 対象セルの周囲6方向にセル内のカットがあるかを調べる
            //if ( (pos[0] - ROUND_EPS) <= 0.5) flag++; // x-
            if ( (pos[1] - ROUND_EPS) <= 0.5) flag++; // x+
            if ( (pos[2] - ROUND_EPS) <= 0.5) flag++; // y-
            if ( (pos[3] - ROUND_EPS) <= 0.5) flag++; // y+
            if ( (pos[4] - ROUND_EPS) <= 0.5) flag++; // z-
            if ( (pos[5] - ROUND_EPS) <= 0.5) flag++; // z+
            
            // カットがなく，未ペイントのセルの場合に流体セルとしてペイントする
            if ( (flag == 0) && (mid[m] == 0) )
            {
              mid[m] = tg;
              c++;
            }
            
          }
        }
      }
      break;
      
    case X_PLUS:
      if ( nID[face] < 0 )
      {
#pragma omp parallel for firstprivate(ix, jx, kx, gd, tg) schedule(static) reduction(+:c)
        for (int k=1; k<=kx; k++) {
          for (int j=1; j<=jx; j++) {
            size_t m = _F_IDX_S3D(ix, j, k, ix, jx, kx, gd);
            size_t mp = _F_IDX_S4DEX(0, ix, j, k, 6, ix, jx, kx, gd);
            const float* pos = &cut[mp];
            int flag = 0;
            
            if ( (pos[0] - ROUND_EPS) <= 0.5) flag++; // x-
            //if ( (pos[1] - ROUND_EPS) <= 0.5) flag++; // x+
            if ( (pos[2] - ROUND_EPS) <= 0.5) flag++; // y-
            if ( (pos[3] - ROUND_EPS) <= 0.5) flag++; // y+
            if ( (pos[4] - ROUND_EPS) <= 0.5) flag++; // z-
            if ( (pos[5] - ROUND_EPS) <= 0.5) flag++; // z+
            
            if ( (flag == 0) && (mid[m] == 0) )
            {
              mid[m] = tg;
              c++;
            }
          }
        }
      }
      break;
      
    case Y_MINUS:
      if ( nID[face] < 0 )
      {
#pragma omp parallel for firstprivate(ix, jx, kx, gd, tg) schedule(static) reduction(+:c)
        for (int k=1; k<=kx; k++) {
          for (int i=1; i<=ix; i++) {
            size_t m = _F_IDX_S3D(i, 1, k, ix, jx, kx, gd);
            size_t mp = _F_IDX_S4DEX(0, i, 1, k, 6, ix, jx, kx, gd);
            const float* pos = &cut[mp];
            int flag = 0;
            
            if ( (pos[0] - ROUND_EPS) <= 0.5) flag++; // x-
            if ( (pos[1] - ROUND_EPS) <= 0.5) flag++; // x+
            //if ( (pos[2] - ROUND_EPS) <= 0.5) flag++; // y-
            if ( (pos[3] - ROUND_EPS) <= 0.5) flag++; // y+
            if ( (pos[4] - ROUND_EPS) <= 0.5) flag++; // z-
            if ( (pos[5] - ROUND_EPS) <= 0.5) flag++; // z+
            
            if ( (flag == 0) && (mid[m] == 0) )
            {
              mid[m] = tg;
              c++;
            }
          }
        }
      }
      break;
      
    case Y_PLUS:
      if ( nID[face] < 0 )
      {
#pragma omp parallel for firstprivate(ix, jx, kx, gd, tg) schedule(static) reduction(+:c)
        for (int k=1; k<=kx; k++) {
          for (int i=1; i<=ix; i++) {
            size_t m = _F_IDX_S3D(i, jx, k, ix, jx, kx, gd);
            size_t mp = _F_IDX_S4DEX(0, i, jx, k, 6, ix, jx, kx, gd);
            const float* pos = &cut[mp];
            int flag = 0;
            
            if ( (pos[0] - ROUND_EPS) <= 0.5) flag++; // x-
            if ( (pos[1] - ROUND_EPS) <= 0.5) flag++; // x+
            if ( (pos[2] - ROUND_EPS) <= 0.5) flag++; // y-
            //if ( (pos[3] - ROUND_EPS) <= 0.5) flag++; // y+
            if ( (pos[4] - ROUND_EPS) <= 0.5) flag++; // z-
            if ( (pos[5] - ROUND_EPS) <= 0.5) flag++; // z+
            
            if ( (flag == 0) && (mid[m] == 0) )
            {
              mid[m] = tg;
              c++;
            }
            
          }
        }
      }
      break;
      
    case Z_MINUS:
      if ( nID[face] < 0 )
      {
#pragma omp parallel for firstprivate(ix, jx, kx, gd, tg) schedule(static) reduction(+:c)
        for (int j=1; j<=jx; j++) {
          for (int i=1; i<=ix; i++) {
            size_t m = _F_IDX_S3D(i, j, 1, ix, jx, kx, gd);
            size_t mp = _F_IDX_S4DEX(0, i, j, 1, 6, ix, jx, kx, gd);
            const float* pos = &cut[mp];
            int flag = 0;
            
            if ( (pos[0] - ROUND_EPS) <= 0.5) flag++; // x-
            if ( (pos[1] - ROUND_EPS) <= 0.5) flag++; // x+
            if ( (pos[2] - ROUND_EPS) <= 0.5) flag++; // y-
            if ( (pos[3] - ROUND_EPS) <= 0.5) flag++; // y+
            //if ( (pos[4] - ROUND_EPS) <= 0.5) flag++; // z-
            if ( (pos[5] - ROUND_EPS) <= 0.5) flag++; // z+
            
            if ( (flag == 0) && (mid[m] == 0) )
            {
              mid[m] = tg;
              c++;
            }
            
          }
        }
      }
      break;
      
    case Z_PLUS:
      if ( nID[face] < 0 )
      {
#pragma omp parallel for firstprivate(ix, jx, kx, gd, tg) schedule(static) reduction(+:c)
        for (int j=1; j<=jx; j++) {
          for (int i=1; i<=ix; i++) {
            size_t m = _F_IDX_S3D(i, j, kx, ix, jx, kx, gd);
            size_t mp = _F_IDX_S4DEX(0, i, j, kx, 6, ix, jx, kx, gd);
            const float* pos = &cut[mp];
            int flag = 0;
            
            if ( (pos[0] - ROUND_EPS) <= 0.5) flag++; // x-
            if ( (pos[1] - ROUND_EPS) <= 0.5) flag++; // x+
            if ( (pos[2] - ROUND_EPS) <= 0.5) flag++; // y-
            if ( (pos[3] - ROUND_EPS) <= 0.5) flag++; // y+
            if ( (pos[4] - ROUND_EPS) <= 0.5) flag++; // z-
            //if ( (pos[5] - ROUND_EPS) <= 0.5) flag++; // z+
            
            if ( (flag == 0) && (mid[m] == 0) )
            {
              mid[m] = tg;
              c++;
            }
            
          }
        }
      }
      break;
      
  } // end of switch
  
  
  if ( numProc > 1 )
  {
    unsigned long c_tmp = c;
    if ( paraMngr->Allreduce(&c_tmp, &c, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
  }
  
  return c;
}


// #################################################################
// 孤立したゼロIDのセルを隣接セルIDで埋める
unsigned long VoxInfo::fillIsolatedEmptyCell(int* mid, const int fluid_id)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  int fid = fluid_id;
  unsigned long replaced=0;
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd, fid) schedule(static) reduction(+:replaced)
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        size_t m_p = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        int s = mid[m_p];
        
        if ( s == 0 )
        {
#include "FindexS3D.h"
          
          int qe = mid[m_e];
          int qw = mid[m_w];
          int qn = mid[m_n];
          int qs = mid[m_s];
          int qt = mid[m_t];
          int qb = mid[m_b];
          
          // 隣接セルがすべて固体の場合でなく，かつempty(=0)でもない
          if ( (qw != fid) && (qe != fid) &&
               (qs != fid) && (qn != fid) &&
               (qb != fid) && (qt != fid) &&
              ( qw * qe * qs * qn * qb * qt != 0 )
              )
          {
            mid[m_p] = find_mode_id(fid, qw, qe, qs, qn, qb, qt);
            replaced++;
          }
        }
      }
    }
  }
  
  if ( numProc > 1 )
  {
    unsigned long c_tmp = replaced;
    if ( paraMngr->Allreduce(&c_tmp, &replaced, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
  }
  
  return replaced;
}


// #################################################################
// 孤立した流体セルを探し，周囲の固体媒質で置換，BCindexを修正する
unsigned long VoxInfo::findIsolatedFcell(int* bx, const int fluid_id)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  int fid = fluid_id;
  unsigned long c = 0;
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd, fid) schedule(static) reduction(+:c)
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        size_t m_p = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        int s = bx[m_p];
        
        if ( IS_FLUID(s) )
        {

#include "FindexS3D.h"
          
          int bb[6];
          bb[0] = bx[m_e];
          bb[1] = bx[m_w];
          bb[2] = bx[m_n];
          bb[3] = bx[m_s];
          bb[4] = bx[m_t];
          bb[5] = bx[m_b];
          
          // 隣接セルがすべて固体の場合
          if ( !IS_FLUID(bb[0]) && !IS_FLUID(bb[1]) &&
               !IS_FLUID(bb[2]) && !IS_FLUID(bb[3]) &&
               !IS_FLUID(bb[4]) && !IS_FLUID(bb[5]) )
          {
            
            // 最頻値を求める
            int qe = DECODE_CMP( bb[0] );
            int qw = DECODE_CMP( bb[1] );
            int qn = DECODE_CMP( bb[2] );
            int qs = DECODE_CMP( bb[3] );
            int qt = DECODE_CMP( bb[4] );
            int qb = DECODE_CMP( bb[5] );
          
            // 媒質オーダーの変更
            s |= find_mode_id(fid, qw, qe, qs, qn, qb, qt);
            
            // 固体セルへ状態を変更　
            bx[m_p] = offBit( s, STATE_BIT );
          }
        }
      }
    }
  }
  
  if ( numProc > 1 )
  {
    unsigned long tmp = c;
    if ( paraMngr->Allreduce(&tmp, &c, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
  }
  
  return c;
}



// #################################################################
// IBCのbboxを取得する
bool VoxInfo::findLBCbbox(const int tgt, const int* bid, const float* cut, int* st, int* ed, const int policy)
{
  int qw, qe, qs, qn, qb, qt;
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  int tg = tgt;
  int plcy = policy;
  
  st[0] = ix;
  st[1] = jx;
  st[2] = kx;
  ed[0] = 0;
  ed[1] = 0;
  ed[2] = 0;

#pragma omp parallel firstprivate(ix, jx, kx, gd, tg, plcy) shared(st, ed)
  {
    int s[3] = {ix, jx, kx};
    int e[3] = {0, 0, 0};
    
#pragma omp for schedule(static)
    for (int k=1; k<=kx; k++) {
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          
          size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
          int qq = bid[m];
          
          if ( TEST_BC(qq) ) // 6方向のいずれかにIDがある
          {
            size_t mp = _F_IDX_S4DEX(0, i, j, k, 6, ix, jx, kx, gd);
            const float* pos = &cut[mp];
            
            int paint=0;
            int bx;
            
            for (int l=0; l<6; l++)
            {
              bx = (qq >> l*5) & MASK_5; // カット面のID，カットがなければゼロ
              
              if ( plcy == NEIGHBOR_CELL ) // CELL_MONITORで隣接セルが指定された場合
              {
                if ( bx == tg ) paint++;
              }
              else // CELL_MONITORとそれ以外で，カットを含むセルのみ
              {
                if ( ((pos[l] - ROUND_EPS) <= 0.5) && ( bx == tg ) ) paint++; // セル内部にカットが存在する
              }

            }
            
            // min, max
            if ( paint > 0 )
            {
              int tmp[3] = {i, j, k};
              FB::vec3_min(s, s, tmp);
              FB::vec3_max(e, e, tmp);
            }
            
          }

        }
      }
    }
    
#pragma omp critical
    {
      FB::vec3_min(st, st, s);
      FB::vec3_max(ed, ed, e);
    }
    
  }
  
  int len[3];
  len[0] = ed[0] - st[0] + 1;
  len[1] = ed[1] - st[1] + 1;
  len[2] = ed[2] - st[2] + 1;
  
  // 各方向長さが全て1以上の場合に，コンポーネントが存在する
  if ( (len[0] > 0) && (len[1] > 0) && (len[2] > 0) )
  {
    return true;
  }
  else
  {
    st[0] = 0;
    st[1] = 0;
    st[2] = 0;
    ed[0] = 0;
    ed[1] = 0;
    ed[2] = 0;
    return false;
  }
}


// #################################################################
// BCindexにActive/Inactiveをエンコード
unsigned long VoxInfo::flipInactive(unsigned long& L,
                                    unsigned long& G,
                                    const int id,
                                    const int* mid,
                                    int* bx,
                                    CompoList* cmp)
{
  int c_p, s;
  size_t m;
  unsigned long c=0;
  
  int idd = id;
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        c_p = mid[m];
        s = bx[m];
        
        if ( c_p == idd )
        {
          if ( TEST_BIT(s, ACTIVE_BIT) ) // 活性化してある場合に，不活性にし，カウント
          {
            s = offBit( s, ACTIVE_BIT );
            bx[m] = s;
            c++;
          }
        }
        
      }
    }
  }
  
  if ( c > 0 )
  {
    cmp->setEnsLocal(ON);
  }
  
  L = c;
  
  if ( numProc > 1 )
  {
    unsigned long tmp = c;
    MPI_Allreduce(&tmp, &c, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  }
  
  G = c;
  
  return c;
}


// #################################################################
/* @brief 各方向の張るスクリーンにコンポーネント要素を投影し，面積を求める
 * @param [in]  st     コンポーネントのbbox始点
 * @param [in]  ed     コンポーネントのbbox終点
 * @param [in]  target サーチ対象ID
 * @param [in]  bid    カットID配列
 * @param [out] ss     各軸方向の投影面積
 */
void VoxInfo::projectCompo(const int* st, const int* ed, const int target, const int* bid, int* ss)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  int tg = target;
  
  int ist = st[0];
  int jst = st[1];
  int kst = st[2];
  int ied = ed[0];
  int jed = ed[1];
  int ked = ed[2];
  printf("%d : %d %d %d %d %d %d\n", myRank, ist,ied, jst,jed, kst,ked);
  // 各方向のスクリーンの大きさ
  int lx = ied - ist + 1;
  int ly = jed - jst + 1;
  int lz = ked - kst + 1;

  // screen
  int* sx = new int [ly*lz]; memset(sx, 0, sizeof(int)*ly*lz);
  int* sy = new int [lx*lz]; memset(sy, 0, sizeof(int)*lx*lz);
  int* sz = new int [lx*ly]; memset(sz, 0, sizeof(int)*lx*ly);
  
  // projection
  // 各方向のスクリーンに対して，列にひとつでもコンポがあればフラグをたてる（重複はカウントしない）
#pragma omp parallel for firstprivate(ix, jx, kx, gd, tg) \
firstprivate(ist, jst, kst, ied, jed ,ked, lx, ly, lz) \
schedule(static)
  for (int k=kst; k<=ked; k++) {
    for (int j=jst; j<=jed; j++) {
      for (int i=ist; i<=ied; i++) {
        
        size_t m  = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        int d = bid[m];
        
        // 6面のいずれかにカット
        if ( TEST_BC(d) )
        {
          int b0 = getFaceBID(0, d);
          int b1 = getFaceBID(1, d);
          int b2 = getFaceBID(2, d);
          int b3 = getFaceBID(3, d);
          int b4 = getFaceBID(4, d);
          int b5 = getFaceBID(5, d);
          
          // スクリーン座標
          int ii = i - ist;
          int jj = j - jst;
          int kk = k - kst;

          // X-dir
          if ( (b0 == tg) || (b1 == tg) ) sx[ kk * ly + jj ] = 1;
          
          // Y-dir
          if ( (b2 == tg) || (b3 == tg) ) sy[ kk * lx + ii ] = 1;
          
          // Z-dir
          if ( (b4 == tg) || (b5 == tg) ) sz[ jj * lx + ii ] = 1;
        }
      }
    }
  }
  
  // count 小サイズのループでリダクションだけなのでスレッド化しない
  int ax = 0;
  int ay = 0;
  int az = 0;
  
  for (int k=0; k<lz; k++) {
    for (int j=0; j<ly; j++) {
      ax += sx[ k * ly + j ];
    }
  }
  
  for (int k=0; k<lz; k++) {
    for (int i=0; i<lx; i++) {
      ay += sy[ k * lx + i ];
    }
  }
  
  for (int j=0; j<ly; j++) {
    for (int i=0; i<lx; i++) {
      az += sz[ j * lx + i ];
    }
  }
  
  ss[0] = ax;
  ss[1] = ay;
  ss[2] = az;
  
  if ( numProc > 1 )
  {
    int tmp[3] = {ax, ay, az};
    if ( paraMngr->Allreduce(tmp, ss, 3, MPI_SUM) != CPM_SUCCESS ) Exit(0);
  }
  
  if ( sx ) { delete [] sx; sx=NULL; }
  if ( sy ) { delete [] sy; sy=NULL; }
  if ( sz ) { delete [] sz; sz=NULL; }
  
}


// #################################################################
// cellで保持されるボクセルid配列をスキャンし，coloList[]に登録する
void VoxInfo::scanCell(int* mid, const int* colorList, const int ID_replace, FILE* fp)
{
  int target;
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int sz[3] = {ix, jx, kx};
  int gd = guide;
  
  
  // ID[0]を置換するオプションが指定されている場合（ID_replaceに値がセット）
  if ( ID_replace != 0 )
  {
    target = ID_replace;
    Hostonly_
    {
      printf (   "\n\tID[0] is replaced by ID[%d]\n", target);
      fprintf(fp,"\n\tID[0] is replaced by ID[%d]\n", target);
    }
    
#pragma omp parallel for firstprivate(ix, jx, kx, gd, target) schedule(static)
    for (int k=1; k<=kx; k++) {
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
          if ( mid[m] == 0 ) mid[m] = target;
        }
      }
    }
  }
  
  int r = myRank;
  
  // 内部領域に対して，マイナスとゼロをチェック スレッド化しない
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        int tg = mid[m];
        
        if ( tg<=0 )
        {
          stamped_printf (   "\tVoxel data includes non-positive ID [%d] at (%d, %d, %d) in Rank %d\n", tg, i, j, k, r);
          stamped_fprintf(fp,"\tVoxel data includes non-positive ID [%d] at (%d, %d, %d) in Rank %d\n", tg, i, j, k, r);
          Exit(0);
        }
      }
    }
  }
  

  // 内部領域の媒質IDに対して、カウント
  
  int* colorSet = new int[NoCompo+1];
  
  for (int i=1; i<=NoCompo; i++) colorSet[i]=0;
  
  // colorSet[] ローカルなIDのカウント
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        int tg = mid[m];
        
        if ( tg<1 || tg>NoCompo )
        {
          Hostonly_
          {
            stamped_printf (   "\tError : Scanned ID is out of range.\n");
            stamped_fprintf(fp,"\tError : Scanned ID is out of range.\n");
          }
          Exit(0);
        }
        
        colorSet[tg] = 1;
      }
    }
  }
  
  
  // 外部領域の媒質IDをcolorSetに追加する
  for (int i=1; i<=NoCompo; i++)
  {
    int tg = colorList[i];
    
    if ( tg != 0 ) colorSet[i] = 1;
  }

  

// ##########
#if 0
  cout << "color list" << endl;
  for (int i=1; i<=NoCompo; i++)
  {
    printf("\t%3d : %d\n", i, colorSet[i]);
  }
  
#endif
// ##########

	// colorSet[] の集約
  if ( numProc > 1 )
  {
    int* clist = new int[NoCompo+1];
    for (int i=0; i<NoCompo+1; i++) clist[i] = colorSet[i];
    if ( paraMngr->Allreduce(clist, colorSet, NoCompo+1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
    delete [] clist;
  }
  // この時点で、存在するIDの数は最大でもnumProc個 >> int の範囲内

  
  
  // ボクセルをスキャンしたIDの数と境界条件の数，含まれるIDのリストを表示する
  Hostonly_
  {
    printf("\n---------------------------------------------------------------------------\n\n");
    printf("\t>> Information of Scanned Voxel ID\n\n");
    for (int i=1; i<=NoCompo; i++)
    {
      if (colorSet[i] != 0) printf("\tID = [%d]\n", i);
    }
    fflush(stdout);
    
    fprintf(fp, "\n---------------------------------------------------------------------------\n\n");
    fprintf(fp,"\t>> Information of Scanned Voxel ID\n\n");
    for (int i=1; i<=NoCompo; i++)
    {
      if (colorSet[i] != 0) fprintf(fp,"\tID = [%d]\n", i);
    }
    fflush(fp);
  }
  
  if ( colorSet ) delete [] colorSet;
  
}



// #################################################################
/**
 * @brief 不活性セルに対する断熱マスクをエンコード
 * @param [in,out] bh  BCindex H2
 */
void VoxInfo::setAmaskInactive(int* bh)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        
        int s = bh[m];
        
        // Inactive指定の場合には，セルの全周を断熱にする
        if ( !TEST_BIT(s, ACTIVE_BIT) )
        {
          s = offBit( s, ADIABATIC_W );
          s = offBit( s, ADIABATIC_E );
          s = offBit( s, ADIABATIC_S );
          s = offBit( s, ADIABATIC_N );
          s = offBit( s, ADIABATIC_B );
          s = offBit( s, ADIABATIC_T );
        }
        bh[m] = s;
      }
    }
  }
}



// #################################################################
/**
 * @brief KOSがSOLID_CONDUCTIONの場合の断熱マスクの処理
 * @param [in,out] bh BCindex H2
 * @note S-F面のS側を断熱にする
 */
void VoxInfo::setAmaskSolid(int* bh)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
  
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        size_t m_p = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        int s   = bh[m_p];
        
        if ( !IS_FLUID(s) ) // pセルが固体セルの場合が対象
        {
#include "FindexS3D.h"
          
          int s_e = bh[m_e];
          int s_w = bh[m_w];
          int s_n = bh[m_n];
          int s_s = bh[m_s];
          int s_t = bh[m_t];
          int s_b = bh[m_b];
          
          // X-
          if ( IS_FLUID(s_w) ) // S-F界面であること
          {
            s = offBit( s, ADIABATIC_W );
          }
          
          // X+
          if ( IS_FLUID(s_e) )
          {
            s = offBit( s, ADIABATIC_E );
          }
          
          // Y-
          if ( IS_FLUID(s_s) )
          {
            s = offBit( s, ADIABATIC_S );
          }
          
          // Y+
          if ( IS_FLUID(s_n) )
          {
            s = offBit( s, ADIABATIC_N );
          }
          
          // Z-
          if ( IS_FLUID(s_b) )
          {
            s = offBit( s, ADIABATIC_B );
          }
          
          // Z+
          if ( IS_FLUID(s_t) )
          {
            s = offBit( s, ADIABATIC_T );
          }
        }
        bh[m_p] = s;
      }
    }
  }
}


// #################################################################
/**
 * @brief KOSがTHERMAL_FLOW, THERMAL_FLOW_NATURALの場合の断熱マスクの処理
 * @param [in,out] bh BCindex H2
 * @note S-F面のF側を断熱にする
 */
void VoxInfo::setAmaskThermal(int* bh)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int sz[3] = {ix, jx, kx};
  int gd = guide;
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
  
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        size_t m_p = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        int s   = bh[m_p];
        
        if ( IS_FLUID(s) ) // pセルが流体セルの場合が対象
        {
#include "FindexS3D.h"
          
          int s_e = bh[m_e];
          int s_w = bh[m_w];
          int s_n = bh[m_n];
          int s_s = bh[m_s];
          int s_t = bh[m_t];
          int s_b = bh[m_b];
          
          // X-
          if ( !IS_FLUID(s_w) ) // S-F界面であること
          {
            s = offBit( s, ADIABATIC_W );
          }
          
          // X+
          if ( !IS_FLUID(s_e) )
          {
            s = offBit( s, ADIABATIC_E );
          }
          
          // Y-
          if ( !IS_FLUID(s_s) )
          {
            s = offBit( s, ADIABATIC_S );
          }
          
          // Y+
          if ( !IS_FLUID(s_n) )
          {
            s = offBit( s, ADIABATIC_N );
          }
          
          // Z-
          if ( !IS_FLUID(s_b) )
          {
            s = offBit( s, ADIABATIC_B );
          }
          
          // Z+
          if ( !IS_FLUID(s_t) )
          {
            s = offBit( s, ADIABATIC_T );
          }
        }
        bh[m_p] = s;
      }
    }
  }
}


// #################################################################
/**
 * @brief bx[]に各境界条件の共通のビット情報をエンコードする
 * @param [in,out] bx    BCindex ID
 * @param [in,out] mid   ID配列
 * @param [in]     cvf   コンポーネントの体積率
 * @param [in]     mat   MediumList
 * @param [in,out] cmp   CompoList
 * @param [in,out] Lcell ノードローカルの有効セル数
 * @param [in,out] Gcell グローバルの有効セル数
 * @param [in]     KOS   解くべき方程式の種類 KIND_OF_SOLVER
 */
void VoxInfo::setBCIndexBase(int* bx, int* mid, const float* cvf, const MediumList* mat, CompoList* cmp,
                              unsigned long& Lcell, unsigned long& Gcell, const int KOS)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  size_t nx = (ix+2*gd) * (jx+2*gd) * (kx+2*gd); // ガイドセルを含む全領域を対象にする
  
  // セルの状態を流体で初期化
  for (size_t m=0; m<nx; m++)
  {
    bx[m] = onBit( bx[m], STATE_BIT );
  }
  
  
  // 状態のエンコード
  for (size_t m=0; m<nx; m++)
  {
    int odr = mid[m];
    int s = bx[m];
    
    if ( mat[odr].getState() == FLUID )
    {
      s = onBit( s, STATE_BIT );
    }
    else  // SOLID
    {
      s = offBit( s, STATE_BIT );
    }
    bx[m] = s;
  }


  // 後の境界条件処理でエンコードしない種類の格納順をエンコード
  for (int n=1; n<=NoCompo; n++)
  {
    switch ( cmp[n].getType() )
    {
      case 0: // Medium
      case OBSTACLE:
      case CELL_MONITOR:
        cmp[n].setElement( encodeOrder(n, mid, bx, &cmp[n]) );
    }
  }
  
  
  // BCIndexにそのセルが計算に有効(active)かどうかをエンコードする．KindOfSolverによって異なる
  encActive(Lcell, Gcell, bx, KOS);
  
  
  // Inactive指定のセルを不活性にする
  unsigned long m_L, m_G;
  
  for (int n=1; n<=NoCompo; n++)
  {
    if ( cmp[n].getType() == INACTIVE )
    {
      m_L = m_G = 0;
      cmp[n].setElement( flipInactive(m_L, m_G, n, mid, bx, &cmp[n]) );
      Lcell -= m_L;
      Gcell -= m_G;
    }
  }
  
}



// #################################################################
// 温度境界条件のビット情報をエンコードする
void VoxInfo::setBCIndexH(int* bh1, int* bh2, int* mid, SetBC* BC, const int kos, CompoList* cmp, float* cut, int* bid)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  
  // 断熱マスクを非断熱(1)に初期化する
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
  for (int k=0; k<=kx+1; k++) {
    for (int j=0; j<=jx+1; j++) {
      for (int i=0; i<=ix+1; i++) {
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        bh2[m] |= ( 0x3f << ADIABATIC_W ); // 6bitまとめて初期化
      }
    }
  }
  
  
  // THERMAL_FLOW, THERMAL_FLOW_NATURAL, SOLID_CONDUCTIONのときに，デフォルトとしてSolid-Fluid面を断熱にする
  switch (kos)
  {
    case THERMAL_FLOW:
    case THERMAL_FLOW_NATURAL:
      setAmaskThermal(bh2);
      break;
      
    case SOLID_CONDUCTION:
      setAmaskSolid(bh2);
      break;
  }
  
  
  // 対称境界面に断熱マスクをセット
  for (int face=0; face<NOFACE; face++)
  {
    if ( BC->exportOBC(face)->getClass() == OBC_SYMMETRIC )
    {
      encAmaskSymtrcBC(face, bh2);
    }
  }
  
  
  // 不活性セルの場合の断熱マスク処理
  setAmaskInactive(bh2);
  
  
  
  // bh2の下位5ビットにはコンポーネントのエントリをエンコード
  for (int n=1; n<=NoCompo; n++)
  {
    
    switch ( cmp[n].getType() )
    {
      case ADIABATIC:
        if ( (kos == THERMAL_FLOW) || (kos == THERMAL_FLOW_NATURAL) )
        {
          cmp[n].setElement( encQface(n, bid, bh1, bh2, false, FLUID, &cmp[n]) );
        }
        else if ( (kos == SOLID_CONDUCTION) )
        {
          cmp[n].setElement( encQface(n, bid, bh1, bh2, false, SOLID, &cmp[n]) );
        }
        else if ( (kos == CONJUGATE_HEAT_TRANSFER) )
        {
          cmp[n].setElement( encQface(n, bid, bh1, bh2, false, ANY_STATE, &cmp[n]) );
        }
        else
        {
          Exit(0);
        }
        break;
        
        
      case HEATFLUX:
      case ISOTHERMAL:
        if ( (kos == THERMAL_FLOW) || (kos == THERMAL_FLOW_NATURAL) )
        {
          cmp[n].setElement( encQface(n, bid, bh1, bh2, true, FLUID, &cmp[n]) );
        }
        else if ( (kos == SOLID_CONDUCTION) )
        {
          cmp[n].setElement( encQface(n, bid, bh1, bh2, true, SOLID, &cmp[n]) );
        }
        else if ( (kos == CONJUGATE_HEAT_TRANSFER) )
        {
          cmp[n].setElement( encQface(n, bid, bh1, bh2, true, ANY_STATE, &cmp[n]) );
        }
        else
        {
          Exit(0);
        }
        break;
        
        
      case SPEC_VEL_WH: // 要素数については，setBCIndexV()でカウントしているので，不要
      case OUTFLOW:
        //encQfaceSVO(n, n, mid, bcd, bh1, bh2, deface);
        cmp[n].setEnsLocal(ON);
        break;
        
        
      case TRANSFER:
        switch ( cmp[n].getHtype() )
      {
          case HT_N:
            //
            break;
            
          case HT_S:
          case HT_SN:
          case HT_SF:
            if ( (kos == THERMAL_FLOW) || (kos == THERMAL_FLOW_NATURAL) )
            {
              cmp[n].setElement( encQface(n, bid, bh1, bh2, true, FLUID, &cmp[n]) );
            }
            else
            {
              Hostonly_ printf("\tHeat Transfer(S, SF, SN) can be specified only in the case of 'ThermalFlow' or 'ThermalFlowNatural'\n");
              Exit(0);
            }
            break;
          
          case HT_B:
            if ( (kos == THERMAL_FLOW) || (kos == THERMAL_FLOW_NATURAL) )
            {
              Hostonly_ printf("\tHeat Transfer(B) can be specified only in the case of 'ConjugateHeatTransfer' or 'SolidConduction'\n");
              Exit(0);
            }
            else if ( (kos == SOLID_CONDUCTION) )
            {
              cmp[n].setElement( encQface(n, bid, bh1, bh2, true, SOLID, &cmp[n]) );
            }
            else if ( (kos == CONJUGATE_HEAT_TRANSFER) )
            {
              cmp[n].setElement( encQface(n, bid, bh1, bh2, true, ANY_STATE, &cmp[n]) );
            }
            else
            {
              Exit(0);
            }
            break;
        }
        break;
        
        
      case RADIANT:
        break;
        
        
        // Q BC at Volume; idのガイドセルチェックなし
      case HEAT_SRC:
      case CNST_TEMP:
        cmp[n].setElement( encodeOrder(n, mid, bh2, &cmp[n]) );
        break;
    }
    
  }// end loop
  
  // set gamma coef. for Heat BC
  encHbit(bh1, bh2);
  
}


// #################################################################
// 圧力境界条件のビット情報をエンコードする
unsigned long VoxInfo::setBCIndexP(int* bcd, int* bcp, int* mid, SetBC* BC, CompoList* cmp, int icls, const float* cut, const int* bid, const bool isBinary)
{
  unsigned long surface = 0;
  

  // 初期化 @note ビットを1に初期化する．初期化範囲はガイドセルを含む全領域．セルフェイスの射影処理で必要．
  size_t mx = (size[0]+2*guide) * (size[1]+2*guide) * (size[2]+2*guide);
  
#pragma omp parallel for firstprivate(mx) schedule(static)
  for (size_t m=0; m<mx; m++)
  {
    bcp[m] |= ( 0x3ffff << BC_NDAG_W ); // BC_NDAG_W〜BC_D_Tまで18bitまとめて1に初期化
  }

  
  // 計算領域内の壁面のNeumannBCのマスク処理と固体に隣接するFセルに方向フラグをエンコードし，表面セル数を返す
  if ( isBinary )
  {
    surface = encPbit_N_Binary(bcp);
    //surface = encPbit_N_Cut(bcp, bid, cut, false);
  }
  else
  {
    surface = encPbit_N_Cut(bcp, bid, cut, false);
  }
  

  
  // 外部境界のビットフラグをエンコード
  BoundaryOuter* m_obc=NULL;
  int F;
  
  for (int face=0; face<NOFACE; face++)
  {
    m_obc = BC->exportOBC(face);
    F = m_obc->getClass();
    
    
    switch ( F )
    {
      case OBC_WALL:
      case OBC_SYMMETRIC:
        encPbitOBC(face, bcp, "Neumann", true);
        break;
        
      case OBC_SPEC_VEL:
        encPbitOBC(face, bcp, "Neumann", false);
        break;
        
      case OBC_TRC_FREE:
        encPbitOBC(face, bcp, "Dirichlet", false);
        break;
        
      case OBC_OUTFLOW:
      case OBC_FAR_FIELD:
        if ( m_obc->get_pType() == P_DIRICHLET )
        {
          encPbitOBC(face, bcp, "Dirichlet", false);
        }
        else
        {
          encPbitOBC(face, bcp, "Neumann", false);
        }
        break;
        
      case OBC_INTRINSIC:
        if ( (icls == id_Jet) && (face==0) )
        {
          encPbitOBC(face, bcp, "Neumann", false);
        }
        break;
        
      // 外部境界値を指定
      case OBC_PERIODIC:
        break;
        
      default:
        Exit(0);
    }
  }

  // 内部境界のコンポーネントのエンコード
  
  for (int n=1; n<=NoCompo; n++)
  {
    int m_dir = cmp[n].getBClocation();
    float vec[3] = { (float)cmp[n].nv[0], (float)cmp[n].nv[1], (float)cmp[n].nv[2] };
    
    switch ( cmp[n].getType() )
    {
      case SPEC_VEL:
      case SPEC_VEL_WH:
      case OUTFLOW:
        encPbit_N_IBC(n, mid, bcd, bcp, vec, m_dir);
        break;
    }
  }


  // 全周Neumannフラグのセルと排他性をチェックし，反復行列の非対角要素/対角要素をエンコードする
  encPbit(bcp);



// ########## debug
#if 0
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  size_t m;
  float w, q;
  int i=1;
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
      w = GET_SHIFT_F(bcp[m], BC_N_W);
      q = GET_SHIFT_F(bcp[m], BC_NDAG_W);
      printf("(1, %3d, %3d) N=%f Coef_W=%f\n", j,k,w,q);
    }
  }
#endif
// ##########
  
  return surface;
}



// #################################################################
// bv[]に境界条件のビット情報をエンコードする
void VoxInfo::setBCIndexV(int* bv, const int* mid, int* bp, SetBC* BC, CompoList* cmp, int icls, float* cut, int* bid)
{
  // ガイドセルの媒質情報をチェックし，流束形式のBCの場合にビットフラグをセット
  BoundaryOuter* m_obc=NULL;
  int F;
  
  // 外部境界
  for (int face=0; face<NOFACE; face++)
  {
    m_obc = BC->exportOBC(face);
    F = m_obc->getClass();
    
    switch ( F )
    {
      case OBC_WALL:
        encVbitOBC(face, bv, "solid", true, "check", bp); // 流束形式
        break;
        
      case OBC_SPEC_VEL:
        encVbitOBC(face, bv, "fluid", true, "check", bp); // 流束形式
        break;
      
      case OBC_TRC_FREE:
      case OBC_OUTFLOW:
      case OBC_FAR_FIELD:
      case OBC_SYMMETRIC:
        encVbitOBC(face, bv, "fluid", false, "check", bp); // 境界値指定
        break;
        
      case OBC_INTRINSIC:
        if ( (icls == id_Jet) && (face==0) )
        {
          encVbitOBC(face, bv, "fluid", true, "nocheck", bp); // 流束形式
        }
        break;
        
      case OBC_PERIODIC:
        encVbitOBC(face, bv, "fluid", false, "nocheck", bp); // 境界値指定，内部セルの状態をコピーするので，ガイドセル状態のチェックなし
        break;
        
      default:
        Exit(0);
    }
    
    // 有効セル数をカウントし，集約
    m_obc->setValidCell( countValidCellOBC(face, bv, F) );
  }
  
  // 内部境界のコンポーネントのエンコード
  int m_dir;
  float vec[3];
  
  for (int n=1; n<=NoCompo; n++)
  {
    int m_dir = cmp[n].getBClocation(); // same_direction=1, opposite_direction=2
    float vec[3] = { (float)cmp[n].nv[0], (float)cmp[n].nv[1], (float)cmp[n].nv[2] };
    
    switch ( cmp[n].getType() )
    {
      case SPEC_VEL:
      case SPEC_VEL_WH:
      case OUTFLOW:
        cmp[n].setElement( encVbitIBC(n, bv, bp, bid, vec, m_dir, &cmp[n]) );
        break;
    }
  }
  
}



// #################################################################
// bx[]のコンポーネントエントリを参照して体積率を計算し，圧力損失コンポーネントの場合にはビットを立てる
void VoxInfo::setCmpFraction(CompoList* cmp, int* bx, const float* vf)
{
	size_t m;
  int st[3], ed[3];
  int f, s;
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  for (int n=1; n<=NoCompo; n++)
  {
    if ( cmp[n].isVFraction() ) // 対象のコンポーネント
    {
      cmp[n].getBbox(st, ed);
      
      for (int k=st[2]; k<=ed[2]; k++) {
        for (int j=st[1]; j<=ed[1]; j++) {
          for (int i=st[0]; i<=ed[0]; i++) {
            m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd); //FBUtility::getFindexS3D(sz, gd, i, j, k);
            s = bx[m];
            if ( ( s & MASK_6) == n )
            {
              f = floorf(vf[m]*255.0 + 0.5); // 0.0<vf<1.0 の四捨五入 > 8ビットで量子化
              if ( f > 255 ) Exit(0);
              s |= (f << TOP_VF); // 8ビットで量子化
              if ( cmp[n].isFORCING() ) s = onBit(s, FORCING_BIT);
              bx[m] = s;
            }
          }
        }
      }
      
    }
  }
}


// #################################################################
// コンポーネントの操作に必要な定数の設定
void VoxInfo::setControlVars(const int m_NoCompo, Intrinsic* ExRef)
{
  NoCompo = m_NoCompo;
  Ex      = ExRef;
}


// #################################################################
/**
 @brief コンポーネントに接するdef_faceのIDをもつセルを不活性セルにする
 @param id 対象BCのセルID
 @param def 面指定ID
 @param mid ボクセル配列
 @param bh1 BCindex H1
 @param bh2 BCindex H2
 */
void VoxInfo::setInactive_Compo(int id, int def, int* mid, int* bh1, int* bh2)
{
  int c_p, c_e, c_w, c_n, c_s, c_t, c_b;
  int idd = id;
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        size_t m_p = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        c_p = mid[m_p];
        
        if ( c_p == idd ) // BCの対象セル
        {
#include "FindexS3D.h"
          
          c_e = mid[m_e];
          c_w = mid[m_w];
          c_n = mid[m_n];
          c_s = mid[m_s];
          c_t = mid[m_t];
          c_b = mid[m_b];
          
          // X-
          if ( c_w == def )
          {
            bh1[m_w] = offBit( bh1[m_w], ACTIVE_BIT ); // defセルを不活性化
            bh2[m_w] = offBit( bh2[m_w], ACTIVE_BIT );
          }
          
          // X+
          if ( c_e == def )
          {
            bh1[m_e] = offBit( bh1[m_e], ACTIVE_BIT );
            bh2[m_e] = offBit( bh2[m_e], ACTIVE_BIT );
          }
          
          // Y-
          if ( c_s == def )
          {
            bh1[m_s] = offBit( bh1[m_s], ACTIVE_BIT );
            bh2[m_s] = offBit( bh2[m_s], ACTIVE_BIT );
          }
          
          // Y+
          if ( c_n == def )
          {
            bh1[m_n] = offBit( bh1[m_n], ACTIVE_BIT );
            bh2[m_n] = offBit( bh2[m_n], ACTIVE_BIT );
          }
          
          // Z-
          if ( c_b == def )
          {
            bh1[m_b] = offBit( bh1[m_b], ACTIVE_BIT );
            bh2[m_b] = offBit( bh2[m_b], ACTIVE_BIT );
          }
          
          // Z+
          if ( c_t == def )
          {
            bh1[m_t] = offBit( bh1[m_t], ACTIVE_BIT );
            bh2[m_t] = offBit( bh2[m_t], ACTIVE_BIT );
          }
        }
        
      }
    }
  }
}


// #################################################################
// 計算領域外部のガイドセルに媒質IDをエンコードする（周期境界以外の場合）
void VoxInfo::setMediumOnGC (const int face, int* mid, const int c_id)
{
  if ( nID[face] >= 0 ) return;
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  int tgt= c_id;
  
  switch (face)
  {
    case X_MINUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, tgt) schedule(static)
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          size_t m = _F_IDX_S3D(0, j, k, ix, jx, kx, gd);
          mid[m] = tgt;
        }
      }
      break;
      
    case X_PLUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, tgt) schedule(static)
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          size_t m = _F_IDX_S3D(ix+1, j, k, ix, jx, kx, gd);
          mid[m] = tgt;
        }
      }
      break;
      
    case Y_MINUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, tgt) schedule(static)
      for (int k=1; k<=kx; k++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, 0, k, ix, jx, kx, gd);
          mid[m] = tgt;
        }
      }
      break;
      
    case Y_PLUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, tgt) schedule(static)
      for (int k=1; k<=kx; k++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, jx+1, k, ix, jx, kx, gd);
          mid[m] = tgt;
        }
      }
      break;
      
    case Z_MINUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, tgt) schedule(static)
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, j, 0, ix, jx, kx, gd);
          mid[m] = tgt;
        }
      }
      break;
      
    case Z_PLUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, tgt) schedule(static)
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, j, kx+1, ix, jx, kx, gd);
          mid[m] = tgt;
        }
      }
      break;
  } // end of switch

  
}


// #################################################################
// 計算領域外部のガイドセルに媒質IDをエンコードする（周期境界の場合）
void VoxInfo::setMediumOnGCperiodic (const int face, int* mid, const int prdc_mode)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  // 内部周期境界の場合には，別メソッド
  if ( prdc_mode == BoundaryOuter::prdc_Driver ) return;
  
    
  if ( numProc > 1 )
  {
    switch (face)
    {
      case X_MINUS:
        if ( paraMngr->PeriodicCommS3D(mid, ix, jx, kx, gd, 1, X_DIR, PLUS2MINUS) != CPM_SUCCESS ) Exit(0);
        break;
        
      case X_PLUS:
        if ( paraMngr->PeriodicCommS3D(mid, ix, jx, kx, gd, 1, X_DIR, MINUS2PLUS) != CPM_SUCCESS ) Exit(0);
        break;
        
      case Y_MINUS:
        if ( paraMngr->PeriodicCommS3D(mid, ix, jx, kx, gd, 1, Y_DIR, PLUS2MINUS) != CPM_SUCCESS ) Exit(0);
        break;
        
      case Y_PLUS:
        if ( paraMngr->PeriodicCommS3D(mid, ix, jx, kx, gd, 1, Y_DIR, MINUS2PLUS) != CPM_SUCCESS ) Exit(0);
        break;
        
      case Z_MINUS:
        if ( paraMngr->PeriodicCommS3D(mid, ix, jx, kx, gd, 1, Z_DIR, PLUS2MINUS) != CPM_SUCCESS ) Exit(0);
        break;
        
      case Z_PLUS:
        if ( paraMngr->PeriodicCommS3D(mid, ix, jx, kx, gd, 1, Z_DIR, MINUS2PLUS) != CPM_SUCCESS ) Exit(0);
        break;
    }
  }
  else // 逐次
  {
    if ( nID[face] >= 0 ) return;
    
    switch (face)
    {
      case X_MINUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
        for (int k=1; k<=kx; k++) {
          for (int j=1; j<=jx; j++) {
            for (int i=1-gd; i<=0; i++) {
              size_t m0 = _F_IDX_S3D(i,    j, k, ix, jx, kx, gd);
              size_t m1 = _F_IDX_S3D(i+ix, j, k, ix, jx, kx, gd);
              mid[m0] = mid[m1];
            }
          }
        }
        break;
        
      case X_PLUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
        for (int k=1; k<=kx; k++) {
          for (int j=1; j<=jx; j++) {
            for (int i=ix+1; i<=ix+gd; i++) {
              size_t m0 = _F_IDX_S3D(i,    j, k, ix, jx, kx, gd);
              size_t m1 = _F_IDX_S3D(i-ix, j, k, ix, jx, kx, gd);
              mid[m0] = mid[m1];
            }
          }
        }
        break;
        
      case Y_MINUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
        for (int k=1; k<=kx; k++) {
          for (int j=1-gd; j<=0; j++) {
            for (int i=1; i<=ix; i++) {
              size_t m0 = _F_IDX_S3D(i, j,    k, ix, jx, kx, gd);
              size_t m1 = _F_IDX_S3D(i, j+jx, k, ix, jx, kx, gd);
              mid[m0] = mid[m1];
            }
          }
        }
        break;
        
      case Y_PLUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
        for (int k=1; k<=kx; k++) {
          for (int j=jx+1; j<=jx+gd; j++) {
            for (int i=1; i<=ix; i++) {
              size_t m0 = _F_IDX_S3D(i, j,    k, ix, jx, kx, gd);
              size_t m1 = _F_IDX_S3D(i, j-jx, k, ix, jx, kx, gd);
              mid[m0] = mid[m1];
            }
          }
        }
        break;
        
      case Z_MINUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
        for (int k=1-gd; k<=0; k++) {
          for (int j=1; j<=jx; j++) {
            for (int i=1; i<=ix; i++) {
              size_t m0 = _F_IDX_S3D(i, j, k,    ix, jx, kx, gd);
              size_t m1 = _F_IDX_S3D(i, j, k+kx, ix, jx, kx, gd);
              mid[m0] = mid[m1];
            }
          }
        }
        break;
        
      case Z_PLUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
        for (int k=kx+1; k<=kx+gd; k++) {
          for (int j=1; j<=jx; j++) {
            for (int i=1; i<=ix; i++) {
              size_t m0 = _F_IDX_S3D(i, j, k,    ix, jx, kx, gd);
              size_t m1 = _F_IDX_S3D(i, j, k-kx, ix, jx, kx, gd);
              mid[m0] = mid[m1];
            }
          }
        }
        break;
    }
  }
  
}



// #################################################################
// 対象セルの場合，セルモニターのIDをmid[]にセットする
// @note チェック半径は1なのでモニタ領域は切断面の両側になる
unsigned long VoxInfo::setMonitorCellID(int* mid, const int* bid, const float* cut, const int target, const int fluid, const int policy)
{
  unsigned long c=0;
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  int tg = target;
  
  
  if ( policy == SINGLE_CELL )
  {
#pragma omp parallel for firstprivate(ix, jx, kx, gd, tg) \
schedule(static) reduction(+:c)
    for (int k=1; k<=kx; k++) {
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          
          size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
          int bd = bid[m];
          
          if ( TEST_BC(bd) && (mid[m] == fluid) ) // 6面のいずれかにIDがあり，流体セルの場合
          {
            size_t mp = _F_IDX_S4DEX(0, i, j, k, 6, ix, jx, kx, gd);
            const float* pos = &cut[mp];
            
            int flag = false;
            for (int i=0; i<6; i++)
            {
              if ( ((pos[i] - ROUND_EPS) <= 0.5) && ( ((bd >> i*5) & MASK_5) == tg) ) flag = true; // セル内部にカットが存在する
            }
            
            if ( flag )
            {
              mid[m] = tg;
              c++;
            }
          }
        }
      }
    }
  }
  else if (policy == NEIGHBOR_CELL)
  {
#pragma omp parallel for firstprivate(ix, jx, kx, gd, tg) \
schedule(static) reduction(+:c)
    for (int k=1; k<=kx; k++) {
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          
          size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
          int bd = bid[m];
          
          if ( TEST_BC(bd) && (mid[m] == fluid) ) // 6面のいずれかにIDがあり，流体セルの場合
          {
            int flag = false;
            for (int i=0; i<6; i++)
            {
              if ( ((bd >> i*5) & MASK_5) == tg ) flag = true;
            }
            
            if ( flag )
            {
              mid[m] = tg;
              c++;
            }
          }
        }
      }
    }
  }
  else
  {
    Exit(0);
  }
  
  
  
  unsigned long cl = c;
  
  if ( numProc > 1 )
  {
    unsigned long tmp = cl;
    if ( paraMngr->Allreduce(&tmp, &cl, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
  }
  
  return cl;
}


// #################################################################
// Cell_Monitorで指定するIDでモニタ部分を指定するための準備 (SHAPE_BOX, SHAPE_CYLINDER)
void VoxInfo::setMonitorShape(int* mid, const int n, ShapeMonitor* SM, CompoList* cmp, const REAL_TYPE RefL)
{
  int f_st[3], f_ed[3];
  
  float ctr[3];
  float nv[3];
  float dr[3];
  float m_depth;
  float m_radius;
  float m_width;
  float m_height;
  
  
  // 形状パラメータのセット
  nv[0]   = cmp[n].nv[0];
  nv[1]   = cmp[n].nv[1];
  nv[2]   = cmp[n].nv[2];
  
  dr[0]   = cmp[n].dr[0];
  dr[1]   = cmp[n].dr[1];
  dr[2]   = cmp[n].dr[2];
  
  ctr[0]  = cmp[n].oc[0]/(float)RefL;
  ctr[1]  = cmp[n].oc[1]/(float)RefL;
  ctr[2]  = cmp[n].oc[2]/(float)RefL;
  
  m_depth = cmp[n].depth/(float)RefL;
  m_radius= cmp[n].shp_p1/(float)RefL;
  m_width = cmp[n].shp_p1/(float)RefL;
  m_height= cmp[n].shp_p2/(float)RefL;
  
  
  switch ( cmp[n].get_Shape() )
  {
    case SHAPE_BOX:
      SM->setShapeParam(nv, ctr, dr, m_depth, m_width, m_height);
      break;
      
    case SHAPE_CYLINDER:
      SM->setShapeParam(nv, ctr, m_depth, m_radius);
      break;
      
    default:
      Exit(0);
      break;
  }
  
  // 回転角度の計算
  SM->get_angle();
  
  // bboxと投影面積の計算
  cmp[n].area = SM->get_BboxArea();
  
  // インデクスの計算 > あとでresize
  SM->bbox_index(f_st, f_ed);
  
  // インデクスのサイズ登録
  cmp[n].setBbox(f_st, f_ed);
  
  SM->setID(f_st, f_ed, mid, n);
  
}



// #################################################################
// 外部境界のガイドセルが固体の場合に距離情報をセット
void VoxInfo::setOBCcut(SetBC* BC, float* cut, int* bid)
{
  float pos=0.5f;
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  for (int face=0; face<NOFACE; face++)
  {
    if( nID[face] >= 0 ) continue;
    
    BoundaryOuter* m_obc = BC->exportOBC(face);
    int F  = m_obc->getClass();
    int id = m_obc->getGuideMedium();
    
    
    if ( F == OBC_WALL )
    {
      switch (face)
      {
        case X_MINUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, pos, id) schedule(static)
          for (int k=1; k<=kx; k++) {
            for (int j=1; j<=jx; j++) {
              size_t m = _F_IDX_S4DEX(X_MINUS, 1, j, k, 6, ix, jx, kx, gd);
              cut[m] = pos;
              size_t l = _F_IDX_S3D(1  , j  , k  , ix, jx, kx, gd);
              int q = bid[l];
              setFaceBID(q, X_MINUS, id);
              bid[l] = q;
            }
          }
          break;
          
        case X_PLUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, pos, id) schedule(static)
          for (int k=1; k<=kx; k++) {
            for (int j=1; j<=jx; j++) {
              size_t m = _F_IDX_S4DEX(X_PLUS, ix, j, k, 6, ix, jx, kx, gd);
              cut[m] = pos;
              size_t l = _F_IDX_S3D(ix  , j  , k  , ix, jx, kx, gd);
              int q = bid[l];
              setFaceBID(q, X_PLUS, id);
              bid[l] = q;
            }
          }
          break;
          
        case Y_MINUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, pos, id) schedule(static)
          for (int k=1; k<=kx; k++) {
            for (int i=1; i<=ix; i++) {
              size_t m = _F_IDX_S4DEX(Y_MINUS, i, 1, k, 6, ix, jx, kx, gd);
              cut[m] = pos;
              size_t l = _F_IDX_S3D(i  , 1  , k  , ix, jx, kx, gd);
              int q = bid[l];
              setFaceBID(q, Y_MINUS, id);
              bid[l] = q;
            }
          }
          break;
          
        case Y_PLUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, pos, id) schedule(static)
          for (int k=1; k<=kx; k++) {
            for (int i=1; i<=ix; i++) {
              size_t m = _F_IDX_S4DEX(Y_PLUS, i, jx, k, 6, ix, jx, kx, gd);
              cut[m] = pos;
              size_t l = _F_IDX_S3D(i  , jx  , k  , ix, jx, kx, gd);
              int q = bid[l];
              setFaceBID(q, Y_PLUS, id);
              bid[l] = q;
            }
          }
          break;
          
        case Z_MINUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, pos, id) schedule(static)
          for (int j=1; j<=jx; j++) {
            for (int i=1; i<=ix; i++) {
              size_t m = _F_IDX_S4DEX(Z_MINUS, i, j, 1, 6, ix, jx, kx, gd);
              cut[m] = pos;
              size_t l = _F_IDX_S3D(i  , j  , 1  , ix, jx, kx, gd);
              int q = bid[l];
              setFaceBID(q, Z_MINUS, id);
              bid[l] = q;
            }
          }
          break;
          
        case Z_PLUS:
#pragma omp parallel for firstprivate(ix, jx, kx, gd, pos, id) schedule(static)
          for (int j=1; j<=jx; j++) {
            for (int i=1; i<=ix; i++) {
              size_t m = _F_IDX_S4DEX(Z_PLUS, i, j, kx, 6, ix, jx, kx, gd);
              cut[m] = pos;
              size_t l = _F_IDX_S3D(i  , j  , kx  , ix, jx, kx, gd);
              int q = bid[l];
              setFaceBID(q, Z_PLUS, id);
              bid[l] = q;
            }
          }
          break;
      }
      
    }
  }
  
}




// #################################################################
// ボクセルモデルにカット情報から得られた固体情報を転写する
unsigned long VoxInfo::SolidFromCut(int* mid, const int* bid, const float* cut, CompoList* cmp)
{
  unsigned long c=0;
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  // チェック用のリスト
  int* list = new int[NoCompo+1];
  for (int n=1; n<=NoCompo; n++) list[n] = 0;
  
  //  対象とするIDとして，非媒質(BC+OBSTACLE)かつ固体属性をリストアップ
  for (int n=1; n<=NoCompo; n++)
  {
    if ( !cmp[n].isKindMedium() && !cmp[n].isFluid() ) list[n] = n;
  }
  
  
  /*　本来は i+1/2 上にあるはずのカットが，浮動小数点の丸め誤差の影響で不確定になる問題
   
      i           i+1/2          i+1
   ---+-------------+-------------+---
                  ^   ^
                  <-+->
                  2*eps
   
   安全側に評価する場合，d_i^+ の評価のとき epsの幅だけ拡張して0.5+epsをセル内と考える
   同様に，d_(i+1)^- の評価のときにも 0.5+epsで抑える
   
   d - 0.5の計算は桁落ちする可能性があるが，ゼロに近くなる?
   */
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd, list) \
            schedule(static) reduction(+:c)
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        int bd = bid[m];
        int bx = 0;
        
        if ( TEST_BC(bd) ) // 6面のいずれかにIDがある
        {
          size_t mp = _F_IDX_S4DEX(0, i, j, k, 6, ix, jx, kx, gd);
          const float* pos = &cut[mp];
          int tgt = -1;
          
          for (int l=0; l<6; l++)
          {
            if ( (pos[l] - 0.5) <= ROUND_EPS ) // セル内部にカットが存在する
            {
              bx = (bd >> l*5)  & MASK_5; // カット面のID，カットがなければゼロ
              tgt = max(tgt, bx); // 6方向のカット面のIDのうち最大のものを使う
            }
          }
          
          if ( tgt > 0 )
          {
            // 対象IDの確認
            int flag = 0;
            for (int n=1; n<=NoCompo; n++)
            {
              if ( list[n] == tgt ) flag++;
            }
            
            if ( flag > 0 )
            {
              mid[m] = tgt;
              c++;
            }
            
          }
          
        } // TEST_BC
        
      }
    }
  }
  
  unsigned long cl = c;
  
  if ( numProc > 1 )
  {
    unsigned long tmp = cl;
    if ( paraMngr->Allreduce(&tmp, &cl, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
  }
  
  if ( list ) delete [] list;
  
  return cl;
}
