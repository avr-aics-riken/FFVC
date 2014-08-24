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
 * @file   ffv_Geometry.C
 * @brief  FFV geometry related functions
 * @author aics
 */

#include "ffv.h"



// #################################################################
/* @brief d_mid[]の対象IDに対して、d_pvf[]に指定値を代入する
 * @param [in] target 対象ID
 * @param [in] value  指定値
 * @retval 指定したセルの数
 */
unsigned long FFV::assignVF(const int target, const REAL_TYPE value)
{
  int ix, jx, kx, gd;

  ix = size[0];
  jx = size[1];
  kx = size[2];
  gd = guide;

  
  int tid = target;
  REAL_TYPE val = value;
  unsigned long c = 0; /// painted count
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd, tid, val) schedule(static) reduction(+:c)
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        
        if ( d_mid[m] == tid )
        {
          d_pvf[m] = val;
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
// ポリゴングループの座標値からboxを計算する
void FFV::calcBboxFromPolygonGroup()
{
  Vec3f m_min, m_max;
  Vec3f t1(poly_org), t2(poly_dx), t3;
  
  t3.assign((float)size[0]*t2.x, (float)size[1]*t2.y, (float)size[2]*t2.z);
  
  // サブドメインの1層外側までをサーチ対象とする
  m_min = t1 - t2;
  m_max = t1 + t3 + t2;
  //printf("Search area Bbox min : %f %f %f\n", m_min.x, m_min.y, m_min.z);
  //printf("Search area Bbox max : %f %f %f\n", m_max.x, m_max.y, m_max.z);
  //printf("\tBounding box for polygon group\n");
  
  // ポリゴン情報へのアクセス
  vector<PolygonGroup*>* pg_roots = PL->get_root_groups();
  vector<PolygonGroup*>::iterator it;
  
  // ポリゴングループのループ
  int m=0;
  for (it = pg_roots->begin(); it != pg_roots->end(); it++)
  {
    std::string m_pg = (*it)->get_name();     // グループラベル
    int ntria = (*it)->get_group_num_tria();  // ローカルのポリゴン数
    
    // 対象ポリゴンがある場合のみ
    if ( ntria > 0 )
    {
      // false; ポリゴンが一部でもかかればピックアップ
      vector<Triangle*>* trias = PL->search_polygons(m_pg, m_min, m_max, false);
      
      Vec3f *p;
      Vec3f bbox_min( 1.0e6,  1.0e6,  1.0e6);
      Vec3f bbox_max(-1.0e6, -1.0e6, -1.0e6);
      unsigned c=0;
      vector<Triangle*>::iterator it2;
      
      for (it2 = trias->begin(); it2 != trias->end(); it2++)
      {
        Vertex** org = (*it2)->get_vertex();
        Vec3f p[3];
        p[0] = *(org[0]);
        p[1] = *(org[1]);
        p[2] = *(org[2]);
        
        CompoFraction::get_min(bbox_min, p[0]);
        CompoFraction::get_min(bbox_min, p[1]);
        CompoFraction::get_min(bbox_min, p[2]);
        
        CompoFraction::get_max(bbox_max, p[0]);
        CompoFraction::get_max(bbox_max, p[1]);
        CompoFraction::get_max(bbox_max, p[2]);
        
#if 0
        printf("%d : p0=(%10.3e %10.3e %10.3e)  p1=(%10.3e %10.3e %10.3e) p2=(%10.3e %10.3e %10.3e) n=(%10.3e %10.3e %10.3e)\n", c++,
               p[0].x, p[0].y, p[0].z,
               p[1].x, p[1].y, p[1].z,
               p[2].x, p[2].y, p[2].z,
               n.x, n.y, n.z);
#endif
      }
      
      PG[m].setBboxMin(bbox_min);
      PG[m].setBboxMax(bbox_max);
      
#if 0
      printf("[%d] %20s : (%10.3e %10.3e %10.3e) - (%10.3e %10.3e %10.3e)\n",
             myRank, m_pg.c_str(), bbox_min.x, bbox_min.y, bbox_min.z,
             bbox_max.x, bbox_max.y, bbox_max.z);
#endif
      
      delete trias; // 後始末
    }
    else // ntria == 0
    {
      Vec3f dummy(0.0, 0.0, 0.0);
      PG[m].setBboxMin(dummy);
      PG[m].setBboxMax(dummy);
    }
    
    m++;
  }
  
  //printf("R[%d] : (%10.3e %10.3e %10.3e) - (%10.3e %10.3e %10.3e)\n",
  //       myRank, origin[0], origin[1], origin[2],
  //       origin[0] + region[0], origin[1] + region[1], origin[2] + region[2]);
  
  // 領域内に収まっているかどうかをチェック >> ポリゴンは少しでも触れれば対象となり、領域外にはみ出すことがある
  for (int i=0; i<C.num_of_polygrp; i++)
  {
    Vec3f b_min = PG[i].getBboxMin();
    Vec3f b_max = PG[i].getBboxMax();
    
    float f_min[3];
    float f_max[3];
    f_min[0] = b_min.x;
    f_min[1] = b_min.y;
    f_min[2] = b_min.z;
    f_max[0] = b_max.x;
    f_max[1] = b_max.y;
    f_max[2] = b_max.z;
    
#if 0
    printf("B[%d] %20s : (%10.3e %10.3e %10.3e) - (%10.3e %10.3e %10.3e)\n",
           myRank, PG[i].getGroup().c_str(), b_min.x, b_min.y, b_min.z,
           b_max.x, b_max.y, b_max.z);
#endif
    
    for (int q=0; q<3; q++)
    {
      if ( f_min[q] < origin[q] ) f_min[q] = origin[q];
      float tmp = origin[q] + region[q];
      if ( f_max[q] > tmp ) f_max[q] = tmp;
    }
    
    Vec3f rmin(f_min);
    Vec3f rmax(f_max);
    
    PG[i].setBboxMin(rmin);
    PG[i].setBboxMax(rmax);
    
#if 0
    printf("A[%d] %20s : (%10.3e %10.3e %10.3e) - (%10.3e %10.3e %10.3e)\n",
           myRank, PG[i].getGroup().c_str(), rmin.x, rmin.y, rmin.z,
           rmax.x, rmax.y, rmax.z);
#endif
  }
  
}



// #################################################################
/* @brief ポリゴンの場合のフィル操作
 * @param [in] fp    ファイルポインタ
 */
void FFV::fill(FILE* fp)
{
  
  // フィル媒質のチェック
  if ( cmp[C.FillID].getState() != FLUID )
  {
    Hostonly_ printf("\tSpecified Medium of filling fluid is not FLUID\n");
    Exit(0);
  }
  
  
  
  unsigned long target_count; ///< フィルの対象となるセル数
  unsigned long replaced;     ///< 置換された数
  unsigned long filled;       ///< FLUIDでフィルされた数
  unsigned long sum_replaced; ///< 置換された数の合計
  unsigned long sum_filled;   ///< FLUIDでフィルされた数の合計
  
  
  
  // 最初にフィル対象のセル数を求める >> 全計算内部セル数
  unsigned long total_cell = (unsigned long)size[0] * (unsigned long)size[1] * (unsigned long)size[2];
  
  
  if ( numProc > 1 )
  {
    unsigned long tmp_fc = total_cell;
    if ( paraMngr->Allreduce(&tmp_fc, &total_cell, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
  }
  
  
  target_count = total_cell;
  
  Hostonly_
  {
    printf(    "\tFill initialize -----\n\n");
    fprintf(fp,"\tFill initialize -----\n\n");
    
    printf    ("\t\tInitial target count   = %16ld\n", target_count);
    fprintf(fp,"\t\tInitial target count   = %16ld\n", target_count);
  }
  
  
  // 定義点上に交点がある場合の処理 >> カットするポリゴンのエントリ番号でフィルする
  unsigned long fill_cut = V.fillCutOnCellCenter(d_bcd, d_bid, d_cut);
  target_count -= fill_cut;
  
  Hostonly_
  {
    if ( fill_cut > 0 )
    {
      printf(    "\t\tFill center cut        = %16ld\n", fill_cut);
      fprintf(fp,"\t\tFill center cut        = %16ld\n", fill_cut);
      
      
      printf    ("\t\tTarget count           = %16ld\n\n", target_count);
      fprintf(fp,"\t\tTarget count           = %16ld\n\n", target_count);
    }
  }
  
  if ( numProc > 1 )
  {
    if ( paraMngr->BndCommS3D(d_bid, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
  }
  
  
  Hostonly_
  {
    printf(    "\n\tFill -----\n\n");
    printf(    "\t\tFilling Fluid Medium   : %s\n", mat[C.FillID].getAlias().c_str());
    printf(    "\t\tHint of Seeding Dir.   : %s\n", FBUtility::getDirection(C.FillSeedDir).c_str());
    printf(    "\t\tFill Seed Medium       : %s\n", mat[C.SeedID].getAlias().c_str());
    printf(    "\t\tFill Control (X, Y, Z) : (%s, %s, %s)\n\n",
           ( !C.FillSuppress[0] ) ? "Suppress" : "Fill",
           ( !C.FillSuppress[1] ) ? "Suppress" : "Fill",
           ( !C.FillSuppress[2] ) ? "Suppress" : "Fill");
    
    fprintf(fp,"\n\tFill -----\n\n");
    fprintf(fp,"\t\tFilling Fluid Medium   : %s\n", mat[C.FillID].getAlias().c_str());
    fprintf(fp,"\t\tHint of Seeding Dir.   : %s\n", FBUtility::getDirection(C.FillSeedDir).c_str());
    fprintf(fp,"\t\tFill Seed Medium       : %s\n", mat[C.SeedID].getAlias().c_str());
    fprintf(fp,"\t\tFill Control (X, Y, Z) : (%s, %s, %s)\n\n",
            ( !C.FillSuppress[0] ) ? "Suppress" : "Fill",
            ( !C.FillSuppress[1] ) ? "Suppress" : "Fill",
            ( !C.FillSuppress[2] ) ? "Suppress" : "Fill");
  }
  
  
  // ヒントが与えられている場合
  filled = V.fillSeedBcd(d_bcd, C.FillSeedDir, C.SeedID, d_bid);
  
  if ( numProc > 1 )
  {
    if ( paraMngr->BndCommS3D(d_bid, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
  }
  
  if ( filled == 0 )
  {
    Hostonly_
    {
      printf(    "No cells to paint\n");
      fprintf(fp,"No cells to paint\n");
    }
    Exit(0);
  }
  
  Hostonly_
  {
    printf(    "\t\tPainted cells          = %16ld\n", filled);
    fprintf(fp,"\t\tPainted cells          = %16ld\n", filled);
  }
  
  // ペイントされたシードセル数をターゲットから差し引く
  target_count -= filled;
  
  
  
  Hostonly_
  {
    printf(    "\t\tRemaining target cells = %16ld\n\n", target_count);
    fprintf(fp,"\t\tRemaining target cells = %16ld\n\n", target_count);
  }
  
  //Ex->writeSVX(d_bcd, &C); Exit(0);
  
  
  // 隣接する流体セルと接続しており，かつ固体セルに挟まれていないセルのみペイントする
  
  int c=-1; // iteration
  sum_replaced = 0;
  sum_filled = 0;
  
  while (target_count > 0) {
    
    // SeedIDで指定された媒質でフィルする．FLUID/SOLIDの両方のケースがある
    unsigned long fs;
    filled = V.fillByBid(d_bid, d_bcd, d_cut, C.SeedID, C.FillSuppress, fs);
    replaced = fs;
    
    if ( numProc > 1 )
    {
      if ( paraMngr->BndCommS4DEx(d_cut, 6, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
      if ( paraMngr->BndCommS3D(d_bcd, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
      if ( paraMngr->BndCommS3D(d_bid, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
    }
    
    target_count -= filled;
    target_count -= replaced;
    sum_filled   += filled;
    sum_replaced += replaced;
    
    if ( filled <= 0 ) break; // フィル対象がなくなったら終了
    c++;
  }
  
  
  Hostonly_
  {
    printf(    "\t\tBID Iteration          = %5d\n", c+1);
    fprintf(fp,"\t\tBID Iteration          = %5d\n", c+1);
    printf(    "\t\t    Filled by [%02d]     = %16ld\n", C.SeedID, sum_filled);
    fprintf(fp,"\t\t    Filled by [%02d]     = %16ld\n", C.SeedID, sum_filled);
    printf(    "\t\t    SOLID replaced     = %16ld\n", sum_replaced);
    fprintf(fp,"\t\t    SOLID replaced     = %16ld\n", sum_replaced);
    printf(    "\t\t    Remaining cell     = %16ld\n\n", target_count);
    fprintf(fp,"\t\t    Remaining cell     = %16ld\n\n", target_count);
  }
  
#if 0
  for (int k=1; k<=size[2]; k++) {
    for (int j=1; j<=size[1]; j++) {
      for (int i=1; i<=size[0]; i++) {
        size_t m = _F_IDX_S3D(i  , j  , k  , size[0], size[1], size[2], guide);
        if ( DECODE_CMP(d_bcd[m]) == 0 )
        {
          printf("(%d %d %d) = %d\n", i,j,k, DECODE_CMP(d_bcd[m]) );
        }
      }
    }
  }
  Ex->writeSVX(d_bcd, &C);
#endif
  
  if ( target_count == 0 ) return;
  
  
  
  
  // 未ペイントのセルを最頻値IDでペイントする -------------------
  
  // 未ペイント（ID=0）のセルを検出
  unsigned long upc = V.countCellB(d_bcd, 0);
  
  if ( upc == 0 )
  {
    Hostonly_
    {
      printf(    "\t\tUnpainted cell         = %16ld\n\n", upc);
      fprintf(fp,"\t\tUnpainted cell         = %16ld\n\n", upc);
    }
  }
  
  
  
  // SeedMediumと反対の媒質に変更
  c = -1;
  sum_replaced = 0;
  int fill_mode = -1;
  if ( cmp[C.SeedID].getState() == FLUID )
  {
    fill_mode = SOLID;
  }
  else
  {
    fill_mode = FLUID;
  }
  
  // 未ペイントのセルに対して、指定媒質でフィルする
  while ( target_count > 0 ) {
    
    if ( fill_mode == SOLID )
    {
      replaced = V.fillByModalSolid(d_bcd, C.FillID, d_bid);
    }
    else
    {
      replaced = V.fillByFluid(d_bcd, C.FillID, d_bid);
    }
    
    if ( numProc > 1 )
    {
      if ( paraMngr->BndCommS3D(d_bcd, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
    }
    
    target_count -= replaced;
    sum_replaced += replaced;
    
    if ( replaced <= 0 ) break;
    c++;
  }
  
  
  Hostonly_
  {
    printf(    "\t\tFinal Filling Iteration= %5d\n", c+1);
    fprintf(fp,"\t\tFinal Filling Iteration= %5d\n", c+1);
    printf(    "\t\t   Filled by %s     = %16ld\n\n", (fill_mode==FLUID)?"FLUID":"SOLID", sum_replaced);
    fprintf(fp,"\t\t   Filled by %s     = %16ld\n\n", (fill_mode==FLUID)?"FLUID":"SOLID", sum_replaced);
  }
  
  
  
#if 0
  for (int k=1; k<=size[2]; k++) {
    for (int j=1; j<=size[1]; j++) {
      for (int i=1; i<=size[0]; i++) {
        size_t m = _F_IDX_S3D(i  , j  , k  , size[0], size[1], size[2], guide);
        if ( DECODE_CMP(d_bcd[m]) == 0 )
        {
          printf("(%d %d %d) = %d\n", i,j,k, DECODE_CMP(d_bcd[m]) );
        }
      }
    }
  }
  Ex->writeSVX(d_bcd, &C);
#endif
  

  
  // ID=0をカウント
  upc = V.countCellB(d_bcd, 0);
  
  if ( upc != 0 )
  {
    Hostonly_
    {
      printf(    "\tFill operation is done, but still remains %ld unpainted cells.\n", upc);
      fprintf(fp,"\tFill operation is done, but still remains %ld unpainted cells.\n", upc);
    }
    Exit(0);
  }
  
}



// #################################################################
/**
 * @brief list[]内の最頻値IDを求める
 * @param [in] m_sz  配列のサイズ
 * @param [in] list  ID配列
 * @param [in] m_noc NoCompo
 * @note 候補がない場合には、0が戻り値
 */
int FFV::find_mode(const int m_sz, const int* list, const int m_noc)
{
  int key[CMP_BIT_W]; ///< ID毎の頻度 @note ffv_Initialize() >> fill()でif ( C.NoCompo+1 > CMP_BIT_W )をチェック
  memset(key, 0, sizeof(int)*CMP_BIT_W);
  
  
  for (int l=0; l<m_sz; l++) key[ list[l] ]++;
  
  
  int mode = 0; // サーチの初期値，IDの大きい方から
  int z = 0;    // 最頻値のID
  
  for (int l=m_noc; l>=1; l--)
  {
    if ( key[l] > mode )
    {
      mode = key[l];
      z = l;
    }
  }
  
  return z;
}


// #################################################################
// @brief セルに含まれるポリゴンを探索し、d_midに記録
unsigned long FFV::findPolygonInCell()
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  // d_midを-1で初期化
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static) collapse(3)
  for (int k=1-gd; k<=kx+gd; k++) {
    for (int j=1-gd; j<=jx+gd; j++) {
      for (int i=1-gd; i<=ix+gd; i++) {
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        d_mid[m] = -1;
      }
    }
  }
  
  
  // ポリゴングループ毎にアクセス
  vector<PolygonGroup*>* pg_roots = PL->get_root_groups();
  vector<PolygonGroup*>::iterator it;
  
  int odr = 0;
  for (it = pg_roots->begin(); it != pg_roots->end(); it++)
  {
    std::string m_pg = (*it)->get_name();     // グループラベル
    std::string m_bc = (*it)->get_type();     // 境界条件ラベル
    int ntria = (*it)->get_group_num_tria();  // ローカルのポリゴン数
    
#if 0
    Hostonly_ printf("\n%s : %s\n", m_pg.c_str(), m_bc.c_str() );
#endif
    
    // 対象ポリゴンがある場合のみ
    if ( ntria > 0 )
    {
      // Monitor属性のポリゴンはスキップ
      if ( strcasecmp(m_bc.c_str(), "monitor"))
      {
        // ポリゴングループのインデクス
        int wmin[3];
        int wmax[3];
        findIndex( PG[odr].getBboxMin(), wmin );
        findIndex( PG[odr].getBboxMax(), wmax );
        
        printf("[%d] (%3d %3d %3d) - (%3d %3d %3d) %s \n", myRank, wmin[0], wmin[1], wmin[2], wmax[0], wmax[1], wmax[2], m_pg.c_str());
        
        //#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static) reduction(+:c) collapse(3), thread safe?
        for (int k=wmin[2]; k<=wmax[2]; k++) {
          for (int j=wmin[1]; j<=wmax[1]; j++) {
            for (int i=wmin[0]; i<=wmax[0]; i++) {
              
              Vec3f bx_min(poly_org[0]+poly_dx[0]*(float)(i-1),
                           poly_org[1]+poly_dx[1]*(float)(j-1),
                           poly_org[2]+poly_dx[2]*(float)(k-1)); // セルBboxの対角座標
              Vec3f bx_max(poly_org[0]+poly_dx[0]*(float)i,
                           poly_org[1]+poly_dx[1]*(float)j,
                           poly_org[2]+poly_dx[2]*(float)k);     // セルBboxの対角座標
              
              vector<Triangle*>* trias = PL->search_polygons(m_pg, bx_min, bx_max, false); // false; ポリゴンが一部でもかかる場合
              int polys = trias->size();
              
              if (polys>0)
              {
                // IDの返却用の配列
                int* ary = new int[polys];
                unsigned c=0;
                vector<Triangle*>::iterator it2;
                
                for (it2 = trias->begin(); it2 != trias->end(); it2++)
                {
                  ary[c] = (*it2)->get_exid();
                  c++;
                }
                
                int z = find_mode(polys, ary, C.NoCompo);
                if ( z == 0 ) Exit(0);
                //printf("(%3d %3d %3d) = %2d [ ", i,j,k, z);
                //for (int l=0; l<polys; l++) printf("%d ", ary[l]);
                //printf("]\n");
                size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
                d_mid[m] = z;
                
                if ( ary ) delete [] ary;
              }
              
              //後始末
              delete trias;
            }
          }
        }
        
      } // skip monitor
    } // ntria
    
    odr++;
  }
  
  
  delete pg_roots;
  
  
  // ID=0 をチェック
  unsigned long c = V.countCellM(d_mid, 0, true);
  
  if ( c != 0 )
  {
    Hostonly_ stamped_printf("\tID=0 was found in water-tight fill process. : %d\n", c);
    Exit(0);
  }
  
  
  // count the number of replaced cells >> except ID=-1
  c = V.countCellM(d_mid, -1, false);
  
  return c;
}



// #################################################################
/* @brief サブセルID配列をrefIDでペイント
 * @param [in,out] svf       サブセルの体積率
 * @param [in,out] smd       サブセルID配列
 * @param [in]     sdv       サブセル分割数
 * @param [in]     dir       ペイント開始方向
 * @param [in]     refID     ペイントID
 * @param [in]     refVf     ペイントする体積率
 * @note single threadで実行
 */
int FFV::SubCellFill(float* svf,
                     int* smd,
                     const int sdv,
                     const int dir,
                     const int refID,
                     const float refVf
                     )
{
  int filled = 0;
  
  if ( dir == X_minus )
  {
    for (int i=1; i<=sdv; i++) {
      for (int k=1; k<=sdv; k++) {
        for (int j=1; j<=sdv; j++) {
          size_t m_p = _F_IDX_S3D(i  , j  , k  , sdv, sdv, sdv, 1);
          size_t m_e = _F_IDX_S3D(i+1, j,   k,   sdv, sdv, sdv, 1);
          size_t m_w = _F_IDX_S3D(i-1, j,   k,   sdv, sdv, sdv, 1);
          size_t m_n = _F_IDX_S3D(i,   j+1, k,   sdv, sdv, sdv, 1);
          size_t m_s = _F_IDX_S3D(i,   j-1, k,   sdv, sdv, sdv, 1);
          size_t m_t = _F_IDX_S3D(i,   j,   k+1, sdv, sdv, sdv, 1);
          size_t m_b = _F_IDX_S3D(i,   j,   k-1, sdv, sdv, sdv, 1);
          
          if ( smd[m_p] == -1 )
          {
            int ff = 0;
            if ( smd[m_w] == refID ) ff++;
            if ( smd[m_e] == refID ) ff++;
            if ( smd[m_s] == refID ) ff++;
            if ( smd[m_n] == refID ) ff++;
            if ( smd[m_b] == refID ) ff++;
            if ( smd[m_t] == refID ) ff++;
            
            if ( ff>0 )
            {
              smd[m_p] = refID;
              svf[m_p] = refVf;
              filled++;
            }
          }
        }
      }
    }
  }
  else if ( dir == X_plus)
  {
    for (int i=sdv; i>=1; i--) {
      for (int k=1; k<=sdv; k++) {
        for (int j=1; j<=sdv; j++) {
          size_t m_p = _F_IDX_S3D(i  , j  , k  , sdv, sdv, sdv, 1);
          size_t m_e = _F_IDX_S3D(i+1, j,   k,   sdv, sdv, sdv, 1);
          size_t m_w = _F_IDX_S3D(i-1, j,   k,   sdv, sdv, sdv, 1);
          size_t m_n = _F_IDX_S3D(i,   j+1, k,   sdv, sdv, sdv, 1);
          size_t m_s = _F_IDX_S3D(i,   j-1, k,   sdv, sdv, sdv, 1);
          size_t m_t = _F_IDX_S3D(i,   j,   k+1, sdv, sdv, sdv, 1);
          size_t m_b = _F_IDX_S3D(i,   j,   k-1, sdv, sdv, sdv, 1);
          
          if ( smd[m_p] == -1 )
          {
            int ff = 0;
            if ( smd[m_w] == refID ) ff++;
            if ( smd[m_e] == refID ) ff++;
            if ( smd[m_s] == refID ) ff++;
            if ( smd[m_n] == refID ) ff++;
            if ( smd[m_b] == refID ) ff++;
            if ( smd[m_t] == refID ) ff++;
            
            if ( ff>0 )
            {
              smd[m_p] = refID;
              svf[m_p] = refVf;
              filled++;
            }
          }
        }
      }
    }
  }
  else if ( dir == Y_minus)
  {
    for (int j=1; j<=sdv; j++) {
      for (int k=1; k<=sdv; k++) {
        for (int i=1; i<=sdv; i++) {
          size_t m_p = _F_IDX_S3D(i  , j  , k  , sdv, sdv, sdv, 1);
          size_t m_e = _F_IDX_S3D(i+1, j,   k,   sdv, sdv, sdv, 1);
          size_t m_w = _F_IDX_S3D(i-1, j,   k,   sdv, sdv, sdv, 1);
          size_t m_n = _F_IDX_S3D(i,   j+1, k,   sdv, sdv, sdv, 1);
          size_t m_s = _F_IDX_S3D(i,   j-1, k,   sdv, sdv, sdv, 1);
          size_t m_t = _F_IDX_S3D(i,   j,   k+1, sdv, sdv, sdv, 1);
          size_t m_b = _F_IDX_S3D(i,   j,   k-1, sdv, sdv, sdv, 1);
          
          if ( smd[m_p] == -1 )
          {
            int ff = 0;
            if ( smd[m_w] == refID ) ff++;
            if ( smd[m_e] == refID ) ff++;
            if ( smd[m_s] == refID ) ff++;
            if ( smd[m_n] == refID ) ff++;
            if ( smd[m_b] == refID ) ff++;
            if ( smd[m_t] == refID ) ff++;
            
            if ( ff>0 )
            {
              smd[m_p] = refID;
              svf[m_p] = refVf;
              filled++;
            }
          }
        }
      }
    }
  }
  else if ( dir == Y_plus)
  {
    for (int j=sdv; j>=1; j--) {
      for (int k=1; k<=sdv; k++) {
        for (int i=1; i<=sdv; i++) {
          size_t m_p = _F_IDX_S3D(i  , j  , k  , sdv, sdv, sdv, 1);
          size_t m_e = _F_IDX_S3D(i+1, j,   k,   sdv, sdv, sdv, 1);
          size_t m_w = _F_IDX_S3D(i-1, j,   k,   sdv, sdv, sdv, 1);
          size_t m_n = _F_IDX_S3D(i,   j+1, k,   sdv, sdv, sdv, 1);
          size_t m_s = _F_IDX_S3D(i,   j-1, k,   sdv, sdv, sdv, 1);
          size_t m_t = _F_IDX_S3D(i,   j,   k+1, sdv, sdv, sdv, 1);
          size_t m_b = _F_IDX_S3D(i,   j,   k-1, sdv, sdv, sdv, 1);
          
          if ( smd[m_p] == -1 )
          {
            int ff = 0;
            if ( smd[m_w] == refID ) ff++;
            if ( smd[m_e] == refID ) ff++;
            if ( smd[m_s] == refID ) ff++;
            if ( smd[m_n] == refID ) ff++;
            if ( smd[m_b] == refID ) ff++;
            if ( smd[m_t] == refID ) ff++;
            
            if ( ff>0 )
            {
              smd[m_p] = refID;
              svf[m_p] = refVf;
              filled++;
            }
          }
        }
      }
    }
  }
  else if ( dir == Z_minus)
  {
    for (int k=1; k<=sdv; k++) {
      for (int j=1; j<=sdv; j++) {
        for (int i=1; i<=sdv; i++) {
          size_t m_p = _F_IDX_S3D(i  , j  , k  , sdv, sdv, sdv, 1);
          size_t m_e = _F_IDX_S3D(i+1, j,   k,   sdv, sdv, sdv, 1);
          size_t m_w = _F_IDX_S3D(i-1, j,   k,   sdv, sdv, sdv, 1);
          size_t m_n = _F_IDX_S3D(i,   j+1, k,   sdv, sdv, sdv, 1);
          size_t m_s = _F_IDX_S3D(i,   j-1, k,   sdv, sdv, sdv, 1);
          size_t m_t = _F_IDX_S3D(i,   j,   k+1, sdv, sdv, sdv, 1);
          size_t m_b = _F_IDX_S3D(i,   j,   k-1, sdv, sdv, sdv, 1);
          
          if ( smd[m_p] == -1 )
          {
            int ff = 0;
            if ( smd[m_w] == refID ) ff++;
            if ( smd[m_e] == refID ) ff++;
            if ( smd[m_s] == refID ) ff++;
            if ( smd[m_n] == refID ) ff++;
            if ( smd[m_b] == refID ) ff++;
            if ( smd[m_t] == refID ) ff++;
            
            if ( ff>0 )
            {
              smd[m_p] = refID;
              svf[m_p] = refVf;
              filled++;
            }
          }
        }
      }
    }
  }
  else if ( dir == Z_plus)
  {
    for (int k=sdv; k>=1; k--) {
      for (int j=1; j<=sdv; j++) {
        for (int i=1; i<=sdv; i++) {
          size_t m_p = _F_IDX_S3D(i  , j  , k  , sdv, sdv, sdv, 1);
          size_t m_e = _F_IDX_S3D(i+1, j,   k,   sdv, sdv, sdv, 1);
          size_t m_w = _F_IDX_S3D(i-1, j,   k,   sdv, sdv, sdv, 1);
          size_t m_n = _F_IDX_S3D(i,   j+1, k,   sdv, sdv, sdv, 1);
          size_t m_s = _F_IDX_S3D(i,   j-1, k,   sdv, sdv, sdv, 1);
          size_t m_t = _F_IDX_S3D(i,   j,   k+1, sdv, sdv, sdv, 1);
          size_t m_b = _F_IDX_S3D(i,   j,   k-1, sdv, sdv, sdv, 1);
          
          if ( smd[m_p] == -1 )
          {
            int ff = 0;
            if ( smd[m_w] == refID ) ff++;
            if ( smd[m_e] == refID ) ff++;
            if ( smd[m_s] == refID ) ff++;
            if ( smd[m_n] == refID ) ff++;
            if ( smd[m_b] == refID ) ff++;
            if ( smd[m_t] == refID ) ff++;
            
            if ( ff>0 )
            {
              smd[m_p] = refID;
              svf[m_p] = refVf;
              filled++;
            }
          }
        }
      }
    }
  }
  
  return filled;
}


// #################################################################
/* @brief サブセルのポリゴン含有テスト
 * @param [in,out] svf       サブセルの体積率
 * @param [in,out] smd       サブセルID配列
 * @param [in]     sdv       サブセル分割数
 * @param [in]     ip,jp,kp  プライマリセルインデクス
 * @param [in]     pch       サブセルの格子幅
 * @param [in]     m_pg      ポリゴングループ名
 * @retval ポリゴンを含むセル数
 * @note single threadで実行
 */
int FFV::SubCellIncTest(float* svf,
                        int* smd,
                        const int sdv,
                        const int ip,
                        const int jp,
                        const int kp,
                        const Vec3f pch,
                        const string m_pg
                        )
{
  // プライマリセルの基点（有次元）
  Vec3f o(poly_org[0]+poly_dx[0]*(float)(ip-1),
          poly_org[1]+poly_dx[1]*(float)(jp-1),
          poly_org[2]+poly_dx[2]*(float)(kp-1));
  
  int pic = 0;
  
  // サブセルの体積率を評価（ポリゴンをもつサブセルのみ0.5）
  for (int k=1; k<=sdv; k++) {
    for (int j=1; j<=sdv; j++) {
      for (int i=1; i<=sdv; i++) {
        
        // サブセルのbbox
        Vec3f b1(o.x + pch.x * (float)(i-1),
                 o.y + pch.y * (float)(j-1),
                 o.z + pch.z * (float)(k-1));
        Vec3f b2(o.x + pch.x * (float)i,
                 o.y + pch.y * (float)j,
                 o.z + pch.z * (float)k);
        
        vector<Triangle*>* trias = PL->search_polygons(m_pg, b1, b2, false); // false; ポリゴンが一部でもかかる場合
        int polys = trias->size();
        
        if (polys>0)
        {
          // IDの返却用の配列
          int* ary = new int[polys];
          unsigned c=0;
          vector<Triangle*>::iterator it2;
          
          for (it2 = trias->begin(); it2 != trias->end(); it2++)
          {
            ary[c] = (*it2)->get_exid();
            c++;
          }
          
          int z = find_mode(polys, ary, C.NoCompo);
          if ( z == 0 ) Exit(0);
          size_t m = _F_IDX_S3D(i, j, k, sdv, sdv, sdv, 1);
          svf[m] = 0.5;
          smd[m] = z;
          pic++;
          //printf("(%3d %3d %3d) = %2d [ ", i,j,k, z);
          //for (int l=0; l<polys; l++) printf("%d ", ary[l]);
          //printf("]\n");
          
          if ( ary ) delete [] ary;
        }
        
        delete trias;
      }
    }
  }
  return pic;
}

  
// #################################################################
/* @brief サブセル分割
 * @param [in,out] svf       サブセルの体積率
 * @param [in,out] smd       サブセルID配列
 * @param [in]     sdv       サブセル分割数
 * @param [in]     ip,jp,kp  プライマリセルインデクス
 * @note single threadで実行
 */
void FFV::SubDivision(float* svf,
                      int* smd,
                      const int sdv,
                      const int ip,
                      const int jp,
                      const int kp
                      )
{
  // 外縁部にポリゴンがない面を探す
  int face_flag = 0;
  int c;
  int sum = 0;

  
  // X-
  c = 0;
  for (int k=1; k<=sdv; k++) {
    for (int j=1; j<=sdv; j++) {
      size_t m = _F_IDX_S3D(1, j, k, sdv, sdv, sdv, 1);
      if ( svf[m] > 0.0 ) c++; // ポリゴンを含むセルがある
    }
  }
  printf("X- : %d\n", c);
  if (c==0) face_flag |= (0x1 << X_minus);
  sum += c;
  
  // X+
  c = 0;
  for (int k=1; k<=sdv; k++) {
    for (int j=1; j<=sdv; j++) {
      size_t m = _F_IDX_S3D(sdv, j, k, sdv, sdv, sdv, 1);
      if ( svf[m] > 0.0 ) c++;
    }
  }
  printf("X+ : %d\n", c);
  if (c==0) face_flag |= (0x1 << X_plus);
  sum += c;
  
  // Y-
  c = 0;
  for (int k=1; k<=sdv; k++) {
    for (int i=1; i<=sdv; i++) {
      size_t m = _F_IDX_S3D(i, 1, k, sdv, sdv, sdv, 1);
      if ( svf[m] > 0.0 ) c++;
    }
  }
  printf("Y- : %d\n", c);
  if (c==0) face_flag |= (0x1 << Y_minus);
  sum += c;
  
  // Y+
  c = 0;
  for (int k=1; k<=sdv; k++) {
    for (int i=1; i<=sdv; i++) {
      size_t m = _F_IDX_S3D(i, sdv, k, sdv, sdv, sdv, 1);
      if ( svf[m] > 0.0 ) c++;
    }
  }
  printf("Y+ : %d\n", c);
  if (c==0) face_flag |= (0x1 << Y_plus);
  sum += c;

  // Z-
  c = 0;
  for (int j=1; j<=sdv; j++) {
    for (int i=1; i<=sdv; i++) {
      size_t m = _F_IDX_S3D(i, j, 1, sdv, sdv, sdv, 1);
      if ( svf[m] > 0.0 ) c++;
    }
  }
  printf("Z- : %d\n", c);
  if (c==0) face_flag |= (0x1 << Z_minus);
  sum += c;
  
  // Z+
  c = 0;
  for (int j=1; j<=sdv; j++) {
    for (int i=1; i<=sdv; i++) {
      size_t m = _F_IDX_S3D(i, j, sdv, sdv, sdv, sdv, 1);
      if ( svf[m] > 0.0 ) c++;
    }
  }
  printf("Z+ : %d\n", c);
  if (c==0) face_flag |= (0x1 << Z_plus);
  sum += c;

  printf("face_flag=%d\n", face_flag);
  
  // プライマリセルの隣接ID（確定済み）を参照
  // 確定済みのセルは、W-TプロセスでFillID, SeedIDでペイントしたセル
  int refID = -1;
  int fillDir = -1;
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  // code block
  {
    size_t m;
    int r;
    printf("%d %d\n", C.SeedID, C.FillID);
    if ( TEST_BIT(face_flag, X_minus) )
    {
      
      m = _F_IDX_S3D(ip-1, jp, kp, ix, jx, kx, gd);
      r = d_mid[m];
      if ( (r == C.SeedID) || (r == C.FillID) )
      {
        refID = r;
        fillDir = 0;
      }
    }
    else if ( TEST_BIT(face_flag, X_plus) )
    {
      m = _F_IDX_S3D(ip+1, jp, kp, ix, jx, kx, gd);
      r = d_mid[m];
      if ( (r == C.SeedID) || (r == C.FillID) )
      {
        refID = r;
        fillDir = 1;
      }
    }
    else if ( TEST_BIT(face_flag, Y_minus) )
    {
      m = _F_IDX_S3D(ip, jp-1, kp, ix, jx, kx, gd);
      r = d_mid[m];
      if ( (r == C.SeedID) || (r == C.FillID) )
      {
        refID = r;
        fillDir = 2;
      }
    }
    else if ( TEST_BIT(face_flag, Y_plus) )
    {
      m = _F_IDX_S3D(ip, jp+1, kp, ix, jx, kx, gd);
      r = d_mid[m];
      if ( (r == C.SeedID) || (r == C.FillID) )
      {
        refID = r;
        fillDir = 3;
      }
    }
    else if ( TEST_BIT(face_flag, Z_minus) )
    {
      m = _F_IDX_S3D(ip, jp, kp-1, ix, jx, kx, gd);
      r = d_mid[m];
      if ( (r == C.SeedID) || (r == C.FillID) )
      {
        refID = r;
        fillDir = 4;
      }
    }
    else if ( TEST_BIT(face_flag, Z_plus) )
    {
      m = _F_IDX_S3D(ip, jp, kp+1, ix, jx, kx, gd);
      r = d_mid[m];
      if ( (r == C.SeedID) || (r == C.FillID) )
      {
        refID = r;
        fillDir = 5;
      }
    }
    
    Hostonly_
    {
      printf("\tP-cell (%3d %3d %3d) dir= %d  refID= %d : %d\n", ip, jp, kp, fillDir, refID, face_flag);
    }
    
    // 6面ともポリゴンがある
    if (fillDir == -1) Exit(0);
    
    // check
    if ( refID < 0 ) Exit(0);
  }

  
  // fill
  int target_count = 0; ///< フィルの対象となるセル数

  target_count = sdv*sdv*sdv;
  target_count -= sum;
  
  
  // 開始面をペイント
  c = 0;
  if ( TEST_BIT(face_flag, X_minus) )
  {
    for (int k=1; k<=kx; k++) {
      for (int j=1; j<=jx; j++) {
        size_t m = _F_IDX_S3D(1, j, k, sdv, sdv, sdv, 1);
        
        if ( smd[m] == -1 )
        {
          smd[m] = refID;
          c++;
        }
      }
    }
  }
  else if ( TEST_BIT(face_flag, X_plus) )
  {
    for (int k=1; k<=kx; k++) {
      for (int j=1; j<=jx; j++) {
        size_t m = _F_IDX_S3D(sdv, j, k, sdv, sdv, sdv, 1);
        
        if ( smd[m] == -1 )
        {
          smd[m] = refID;
          c++;
        }
      }
    }
  }
  else if ( TEST_BIT(face_flag, Y_minus) )
  {
    for (int k=1; k<=kx; k++) {
      for (int i=1; i<=ix; i++) {
        size_t m = _F_IDX_S3D(i, 1, k, sdv, sdv, sdv, 1);
        
        if ( smd[m] == -1 )
        {
          smd[m] = refID;
          c++;
        }
      }
    }
  }
  else if ( TEST_BIT(face_flag, Y_plus) )
  {
    for (int k=1; k<=kx; k++) {
      for (int i=1; i<=ix; i++) {
        size_t m = _F_IDX_S3D(i, sdv, k, sdv, sdv, sdv, 1);
        
        if ( smd[m] == -1 )
        {
          smd[m] = refID;
          c++;
        }
      }
    }
  }
  else if ( TEST_BIT(face_flag, Z_minus) )
  {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        size_t m = _F_IDX_S3D(i, j, 1, sdv, sdv, sdv, 1);
        
        if ( smd[m] == -1 )
        {
          smd[m] = refID;
          c++;
        }
      }
    }
  }
  else if ( TEST_BIT(face_flag, Z_plus) )
  {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        size_t m = _F_IDX_S3D(i, j, sdv, sdv, sdv, sdv, 1);
        
        if ( smd[m] == -1 )
        {
          smd[m] = refID;
          c++;
        }
      }
    }
  }
  
  target_count -= c;
  
  int filled = 0;       ///< フィルされた数
  int sum_filled = 0;   ///< フィルされた数の合計
  float refVf = (mat[refID].getState()==FLUID) ? 1.0 : 0.0; ///< refIDの体積率
  
  c = 0;
  while (target_count > 0) {
    
    filled = SubCellFill(svf, smd, sdv, fillDir, refID, refVf);
    
    target_count -= filled;
    sum_filled   += filled;
    c++;
    printf("\titr=%3d %d\n", c+1, filled);
    
    if ( filled <= 0 ) break; // フィル対象がなくなったら終了
  }
  
  // 未ペイント(-1)をチェック
  int flag = 0;
  for (int k=1; k<=sdv; k++) {
    for (int j=1; j<=sdv; j++) {
      for (int i=1; i<=sdv; i++) {
        size_t m = _F_IDX_S3D(i, j, k, sdv, sdv, sdv, 1);
        if ( smd[m] == -1 ) flag++;
      }
    }
  }
  
  // 未ペイント部分がある場合
  if ( flag > 0 )
  {
    // 逆の属性で残りをフィル
    int paintID = (refID == C.SeedID) ? C.FillID : C.SeedID;
    float paintVf = (mat[refID].getState()==FLUID) ? 0.0 : 1.0;
    
    for (int k=1; k<=sdv; k++) {
      for (int j=1; j<=sdv; j++) {
        for (int i=1; i<=sdv; i++) {
          size_t m = _F_IDX_S3D(i, j, k, sdv, sdv, sdv, 1);
          if ( smd[m] == -1 )
          {
            smd[m] = paintID;
            svf[m] = paintVf;
          }
        }
      }
    }
  }
  else // 未ペイント部分なし >> 残りの数はPolyIDの数と等しくなければならない
  {
    if (flag != target_count)
    {
      printf("unresolved cells\n");
      Exit(0);
    }
  }
  
  
  
  //if ( flag > 0 )
  //{
  //  printf("Subcell is not completely filled.\n");
  //  Exit(0);
  //}

  
  // サブセルの体積率からプライマリセルの体積率を求める
  float sff = 0.0;
  
  for (int k=1; k<=sdv; k++) {
    for (int j=1; j<=sdv; j++) {
      for (int i=1; i<=sdv; i++) {
        size_t m = _F_IDX_S3D(i, j, k, sdv, sdv, sdv, 1);
        sff += svf[m];
      }
    }
  }

  float ff = 1.0/(float)(sdv*sdv*sdv);
  d_pvf[_F_IDX_S3D(ip, jp, kp, ix, jx, kx, gd)] = sff*ff;
  
}



// #################################################################
/* @brief サブサンプリング
 * @param [in] fp    ファイルポインタ
 */
void FFV::SubSampling(FILE* fp)
{
  unsigned long target_count = 0; ///< フィルの対象となるセル数
  unsigned long replaced = 0;     ///< 置換された数
  unsigned long filled = 0;       ///< フィルされた数
  unsigned long sum_replaced = 0; ///< 置換された数の合計
  unsigned long sum_filled = 0;   ///< FLUIDでフィルされた数の合計
  
  int ix, jx, kx, gd;
  ix = size[0];
  jx = size[1];
  kx = size[2];
  gd = guide;
  
  // -1.0で初期化
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        d_pvf[m] = -1.0;
      }
    }
  }
  
  // W-T processのInner/Outer fillで確定したセルに対応する体積率を代入初期化
  Hostonly_
  {
    printf(    "\tS-S initialize -----\n\n");
    fprintf(fp,"\tS-S initialize -----\n\n");
  }
  
  // Outer fill => SeedID
  REAL_TYPE tmp = (mat[C.SeedID].getState()==FLUID) ? 1.0 : 0.0;
  filled = assignVF(C.SeedID, tmp);
  
  Hostonly_
  {
    printf    ("\t\tOuter assigned Vf(%3.1f) = %16ld  (%s)\n", tmp, filled, mat[C.SeedID].getAlias().c_str());
    fprintf(fp,"\t\tOuter assigned Vf(%3.1f) = %16ld  (%s)\n", tmp, filled, mat[C.SeedID].getAlias().c_str());
  }
  
  // Inner fill => FillID
  tmp = (mat[C.FillID].getState()==FLUID) ? 1.0 : 0.0;
  filled = assignVF(C.FillID, tmp);
  
  Hostonly_
  {
    printf    ("\t\tInner assigned Vf(%3.1f) = %16ld  (%s)\n", tmp, filled, mat[C.FillID].getAlias().c_str());
    fprintf(fp,"\t\tInner assigned Vf(%3.1f) = %16ld  (%s)\n", tmp, filled, mat[C.FillID].getAlias().c_str());
  }
  
  // 分割数
  int sdv = 10; //C.Hide.Subdivision;
  
  // Local
  REAL_TYPE m_org[3], m_pit[3];
  int m_siz[3];
  
  for (int i=0; i<3; i++)
  {
    m_siz[i] = sdv;
  }

  
  
  
  // sub-cell work array
  size_t nx = (sdv+2)*(sdv+2)*(sdv+2);
  float* svf = new float [nx]; // guide cell is 1 for each dir.
  int* smd = new int [nx];
  
#pragma omp parallel for firstprivate(sdv) schedule(static) collapse(3)
  for (int k=0; k<=sdv+1; k++) {
    for (int j=0; j<=sdv+1; j++) {
      for (int i=0; i<=sdv+1; i++) {
        size_t m = _F_IDX_S3D(i, j, k, sdv, sdv, sdv, 1);
        svf[m] = 0.0;
        smd[m] = -1;
      }
    }
  }
  

  
  // Poly IDsに対して、サブセルテスト
  
//#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        
        for (int kk=0; kk<=sdv+1; kk++) {
          for (int jj=0; jj<=sdv+1; jj++) {
            for (int ii=0; ii<=sdv+1; ii++) {
              size_t m = _F_IDX_S3D(ii, jj, kk, sdv, sdv, sdv, 1);
              svf[m] = 0.0;
              smd[m] = -1;
            }
          }
        }
        
        // サブセルのファイル出力ヘッダ
        for (int l=0; l<3; l++) m_pit[l] = poly_dx[l]/(float)sdv;
        
        m_org[0] = poly_org[0]+poly_dx[0]*(float)(i-1);
        m_org[1] = poly_org[1]+poly_dx[1]*(float)(j-1);
        m_org[2] = poly_org[2]+poly_dx[2]*(float)(k-1);
        
        
        
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        int q = d_mid[m];
        
        // PolyIDsの場合のみ
        if ( (q != C.FillID) && (q != C.SeedID) )
        {
          vector<PolygonGroup*>* pg_roots = PL->get_root_groups();
          vector<PolygonGroup*>::iterator it;
          
          // ポリゴングループのループ
          for (it = pg_roots->begin(); it != pg_roots->end(); it++)
          {
            string m_pg = (*it)->get_name();     // グループラベル
            string m_bc = (*it)->get_type();     // 境界条件ラベル
            int ntria = (*it)->get_group_num_tria();  // ローカルのポリゴン数
            
            // 対象ポリゴンがある場合のみ
            if ( ntria > 0 )
            {
              // Monitor属性のポリゴンはスキップ
              if ( strcasecmp(m_bc.c_str(), "monitor"))
              {
                // サブセルのポリゴン含有テスト
                int cp = SubCellIncTest(svf, smd, sdv, i, j, k, m_pit, m_pg);
                printf("%3d %3d %3d : %3d\n", i,j,k,cp);
              }
            }
          } // Polygon Group
          
          F.writeSVX(smd, i, j, k, m_siz, 1, m_pit, m_org);
          
          /*
          for (it = pg_roots->begin(); it != pg_roots->end(); it++)
          {
            string m_pg = (*it)->get_name();     // グループラベル
            string m_bc = (*it)->get_type();     // 境界条件ラベル
            int ntria = (*it)->get_group_num_tria();  // ローカルのポリゴン数
            
            // 対象ポリゴンがある場合のみ
            if ( ntria > 0 )
            {
              // Monitor属性のポリゴンはスキップ
              if ( strcasecmp(m_bc.c_str(), "monitor"))
              {
                SubDivision(svf, smd, sdv, i, j, k);
              }
            }
          } // Polygon Group
          */
          SubDivision(svf, smd, sdv, i, j, k);
          
          delete pg_roots;
        } // polyIDs
        
      }
    }
  }
  
  
  
  /*
  // ポリゴングループ毎に細分化テスト
  vector<PolygonGroup*>* pg_roots = PL->get_root_groups();
  vector<PolygonGroup*>::iterator it;
  
  // サブセルのポリゴン含有テスト
  int odr = 0;
  for (it = pg_roots->begin(); it != pg_roots->end(); it++)
  {
    std::string m_pg = (*it)->get_name();     // グループラベル
    std::string m_bc = (*it)->get_type();     // 境界条件ラベル
    int ntria = (*it)->get_group_num_tria();  // ローカルのポリゴン数
    
    // 対象ポリゴンがある場合のみ
    if ( ntria > 0 )
    {
      // Monitor属性のポリゴンはスキップ
      if ( strcasecmp(m_bc.c_str(), "monitor"))
      {
        // ポリゴングループのBboxインデクス
        int wmin[3];
        int wmax[3];
        findIndex( PG[odr].getBboxMin(), wmin );
        findIndex( PG[odr].getBboxMax(), wmax );
        
//#pragma omp parallel for firstprivate(sdv) schedule(dynamic) reduction(+:c) collapse(3)
        for (int k=wmin[2]; k<=wmax[2]; k++) {
          for (int j=wmin[1]; j<=wmax[1]; j++) {
            for (int i=wmin[0]; i<=wmax[0]; i++) {
              
              int cp = SubCellIncTest(svf, smd, sdv, i, j, k, m_pg);
              //if ( cp>0 ) printf("%3d %3d %3d : %3d\n", i,j,k,cp);
            }
          }
        }
        
      } // skip monitor
    } // ntria
    
    odr++;
  }
  
  
  //
  odr = 0;
  for (it = pg_roots->begin(); it != pg_roots->end(); it++)
  {
    std::string m_pg = (*it)->get_name();     // グループラベル
    std::string m_bc = (*it)->get_type();     // 境界条件ラベル
    int ntria = (*it)->get_group_num_tria();  // ローカルのポリゴン数
    
    // 対象ポリゴンがある場合のみ
    if ( ntria > 0 )
    {
      // Monitor属性のポリゴンはスキップ
      if ( strcasecmp(m_bc.c_str(), "monitor"))
      {
        // ポリゴングループのBboxインデクス
        int wmin[3];
        int wmax[3];
        findIndex( PG[odr].getBboxMin(), wmin );
        findIndex( PG[odr].getBboxMax(), wmax );
        
//#pragma omp parallel for firstprivate(sdv) schedule(dynamic) collapse(3)
        for (int k=wmin[2]; k<=wmax[2]; k++) {
          for (int j=wmin[1]; j<=wmax[1]; j++) {
            for (int i=wmin[0]; i<=wmax[0]; i++) {
              
              SubDivision(svf, smd, sdv, i, j, k, m_pg);
              
            }
          }
        }
        
      } // skip monitor
    } // ntria
    
    odr++;
  }
  
  delete pg_roots;
  */
  
  if ( svf ) delete [] svf;
  if ( smd ) delete [] smd;
  
  F.writeRawSPH(d_pvf, size, guide, 0, m_org, m_pit, sizeof(REAL_TYPE));
  Exit(0);

}



// #################################################################
/* @brief ポリゴンの水密化
 * @param [in] fp    ファイルポインタ
 * @note ここまで、d_bcdにはsetMonitorList()でモニタIDが入っている
 */
void FFV::WaterTightening(FILE* fp)
{
  unsigned long target_count = 0; ///< フィルの対象となるセル数
  unsigned long replaced = 0;     ///< 置換された数
  unsigned long filled = 0;       ///< FLUIDでフィルされた数
  unsigned long sum_replaced = 0; ///< 置換された数の合計
  unsigned long sum_filled = 0;   ///< FLUIDでフィルされた数の合計
  
  
  // フィル媒質のチェック
  if ( cmp[C.FillID].getState() != FLUID )
  {
    Hostonly_ printf("\tSpecified Medium of filling fluid is not FLUID\n");
    Exit(0);
  }
  
  
  // 最初にフィル対象のセル数を求める >> 全計算内部セル数
  unsigned long total_cell = (unsigned long)size[0] * (unsigned long)size[1] * (unsigned long)size[2];
  
  
  if ( numProc > 1 )
  {
    unsigned long tmp_fc = total_cell;
    if ( paraMngr->Allreduce(&tmp_fc, &total_cell, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
  }
  
  
  target_count = total_cell;
  
  Hostonly_
  {
    printf(    "\tW-T initialize -----\n\n");
    fprintf(fp,"\tW-T initialize -----\n\n");
    
    printf    ("\t\tInitial target count   = %16ld\n", target_count);
    fprintf(fp,"\t\tInitial target count   = %16ld\n", target_count);
    
    printf(    "\t\tFilling Fluid Medium   : %s\n", mat[C.FillID].getAlias().c_str());
    printf(    "\t\tHint of Seeding Dir.   : %s\n", FBUtility::getDirection(C.FillSeedDir).c_str());
    printf(    "\t\tFill Seed Medium       : %s\n", mat[C.SeedID].getAlias().c_str());
    printf(    "\t\tFill Control (X, Y, Z) : (%s, %s, %s)\n\n",
           ( !C.FillSuppress[0] ) ? "Suppress" : "Fill",
           ( !C.FillSuppress[1] ) ? "Suppress" : "Fill",
           ( !C.FillSuppress[2] ) ? "Suppress" : "Fill");
    
    fprintf(fp,"\t\tFilling Fluid Medium   : %s\n", mat[C.FillID].getAlias().c_str());
    fprintf(fp,"\t\tHint of Seeding Dir.   : %s\n", FBUtility::getDirection(C.FillSeedDir).c_str());
    fprintf(fp,"\t\tFill Seed Medium       : %s\n", mat[C.SeedID].getAlias().c_str());
    fprintf(fp,"\t\tFill Control (X, Y, Z) : (%s, %s, %s)\n\n",
            ( !C.FillSuppress[0] ) ? "Suppress" : "Fill",
            ( !C.FillSuppress[1] ) ? "Suppress" : "Fill",
            ( !C.FillSuppress[2] ) ? "Suppress" : "Fill");
  }
  
  
  
  
  // セルに含まれるポリゴンを探索し、d_midに記録
  Hostonly_
  {
    printf(    "\tPaint cells that contain polygons -----\n\n");
    fprintf(fp,"\tPaint cells that contain polygons -----\n\n");
  }
  
  sum_replaced = findPolygonInCell();
  
  target_count -= sum_replaced;
  
  if ( numProc > 1 )
  {
    if ( paraMngr->BndCommS3D(d_mid, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
  }
  
  
  Hostonly_ // Up to here, d_mid={-1, Poly-IDs}
  {
    printf    ("\t\t# of cells touch Polys = %16ld\n", sum_replaced);
    fprintf(fp,"\t\t# of cells touch Polys = %16ld\n", sum_replaced);
    printf(    "\t\tRemaining target cells = %16ld\n\n", target_count);
    fprintf(fp,"\t\tRemaining target cells = %16ld\n\n", target_count);
  }
  
  
  
  
  // ヒントが与えられている場合
  Hostonly_
  {
    printf(    "\tHint of filling -----\n\n");
    fprintf(fp,"\thint of filling -----\n\n");
  }
  
  filled = V.fillSeedMid(d_mid, C.FillSeedDir, C.SeedID);
  target_count -= filled;
  
  if ( filled == 0 )
  {
    Hostonly_
    {
      printf(    "No cells painted\n");
      fprintf(fp,"No cells painted\n");
    }
    Exit(0);
  }
  
  if ( numProc > 1 )
  {
    if ( paraMngr->BndCommS3D(d_mid, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
  }
  
  
  Hostonly_ // Up to here, d_mid={-1, Poly-IDs, SeedID}
  {
    printf(    "\t\tPainted cells          = %16ld  (%s)\n", filled, mat[C.SeedID].getAlias().c_str());
    fprintf(fp,"\t\tPainted cells          = %16ld  (%s)\n", filled, mat[C.SeedID].getAlias().c_str());
    printf(    "\t\tRemaining target cells = %16ld\n\n", target_count);
    fprintf(fp,"\t\tRemaining target cells = %16ld\n\n", target_count);
  }
  
  
  
  
  // ターゲットセルが-1で隣接するセルがSeedIDであるセルをペイントする
  Hostonly_
  {
    printf(    "\tFill outside of objects -----\n\n");
    fprintf(fp,"\tFill outside of objects -----\n\n");
  }
  
  int c=0; // iteration
  sum_filled = 0;
  
  while (target_count > 0) {
    
    // SeedIDで指定された媒質でフィルする
    filled = V.fillByMid(d_mid, C.SeedID);
    
    if ( numProc > 1 )
    {
      if ( paraMngr->BndCommS3D(d_mid, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
    }
    
    target_count -= filled;
    sum_filled   += filled;
    c++;
    
    if ( filled <= 0 ) break; // フィル対象がなくなったら終了
  }
  
  
  Hostonly_
  {
    printf(    "\t\tIteration              = %5d\n", c);
    fprintf(fp,"\t\tIteration              = %5d\n", c);
    printf(    "\t\t    Filled cells       = %16ld  (%s)\n", sum_filled, mat[C.SeedID].getAlias().c_str());
    fprintf(fp,"\t\t    Filled cells       = %16ld  (%s)\n", sum_filled, mat[C.SeedID].getAlias().c_str());
    printf(    "\t\t    Remaining cells    = %16ld\n\n", target_count);
    fprintf(fp,"\t\t    Remaining cells    = %16ld\n\n", target_count);
  }
  
  
  
  
  
  // Solid部分を埋める
  Hostonly_
  {
    printf(    "\tFill inside of objects -----\n\n");
    fprintf(fp,"\tFill inside of objects -----\n\n");
  }
  
  if ( target_count != 0 ) // >> target_countがゼロの場合、中実部分がない
  {
    // 未ペイント（ID=-1）のセルを検出
    unsigned long upc = V.countCellM(d_mid, -1);
    
    if ( target_count != upc )
    {
      Hostonly_
      {
        printf(    "\t\tUnpainted cell         = %16ld\n", upc);
        fprintf(fp,"\t\tUnpainted cell         = %16ld\n", upc);
      }
      Exit(0);
    }
    
    
    // 未ペイントのセルに対して、指定媒質でフィルする
    c = 0;
    sum_replaced = 0;
    
    while ( target_count > 0 ) {
      
      replaced = V.fillByID(d_mid, C.FillID);
      
      if ( numProc > 1 )
      {
        if ( paraMngr->BndCommS3D(d_mid, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
      }
      
      target_count -= replaced;
      sum_replaced += replaced;
      c++;
      
      if ( replaced <= 0 ) break;
    }
    
    Hostonly_
    {
      printf(    "\t\tIteration              = %5d\n", c);
      fprintf(fp,"\t\tIteration              = %5d\n", c);
      printf(    "\t\tFilled cells           = %16ld  (%s)\n\n", sum_replaced, mat[C.FillID].getAlias().c_str());
      fprintf(fp,"\t\tFilled cells           = %16ld  (%s)\n\n", sum_replaced, mat[C.FillID].getAlias().c_str());
    }
    
    
    
    // ID=-1をカウントしてチェック
    upc = V.countCellM(d_mid, -1);
    
    if ( upc != 0 )
    {
      Hostonly_
      {
        printf(    "\tFill operation is done, but still remains %ld unpainted cells.\n", upc);
        fprintf(fp,"\tFill operation is done, but still remains %ld unpainted cells.\n", upc);
      }
      Exit(0);
    }
  }
  
  
  //Ex->writeSVX(d_mid, &C);
  //Exit(0);
  
}

