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
 * @file   Geometry.C
 * @brief  FFV geometry related functions
 * @author aics
 */

#include "Geometry.h"

// #################################################################
/* @brief d_mid[]がtargetであるセルに対して、d_pvf[]に指定値valueを代入する
 * @param [in]  target  キーID
 * @param [in]  value   代入値
 * @param [in]  d_mid   識別子配列
 * @param [out] d_pvf   体積率
 * @retval 値を代入したセルの数
 */
unsigned long Geometry::assignVF(const int target, const REAL_TYPE value, const int* d_mid, REAL_TYPE* d_pvf)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

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
void Geometry::calcBboxFromPolygonGroup()
{
  Vec3<REAL_TYPE> m_min;
  Vec3<REAL_TYPE> m_max;
  Vec3<REAL_TYPE> t1(m_poly_org);
  Vec3<REAL_TYPE> t2(m_poly_dx);
  Vec3<REAL_TYPE> t3;
  
  t3.assign((REAL_TYPE)size[0]*t2.x, (REAL_TYPE)size[1]*t2.y, (REAL_TYPE)size[2]*t2.z);
  
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
      
      Vec3<REAL_TYPE> *p;
      Vec3<REAL_TYPE> bbox_min( 1.0e6,  1.0e6,  1.0e6);
      Vec3<REAL_TYPE> bbox_max(-1.0e6, -1.0e6, -1.0e6);
      unsigned c=0;
      vector<Triangle*>::iterator it2;
      
      for (it2 = trias->begin(); it2 != trias->end(); it2++)
      {
        Vertex** org = (*it2)->get_vertex();
        Vec3<REAL_TYPE> p[3];
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
      Vec3<REAL_TYPE> dummy(0.0, 0.0, 0.0);
      PG[m].setBboxMin(dummy);
      PG[m].setBboxMax(dummy);
    }
    
    m++;
  }
  
  //printf("R[%d] : (%10.3e %10.3e %10.3e) - (%10.3e %10.3e %10.3e)\n",
  //       myRank, origin[0], origin[1], origin[2],
  //       origin[0] + region[0], origin[1] + region[1], origin[2] + region[2]);
  
  // 領域内に収まっているかどうかをチェック >> ポリゴンは少しでも触れれば対象となり、領域外にはみ出すことがある
  for (int i=0; i<m_NoPolyGrp; i++)
  {
    Vec3<REAL_TYPE> b_min = PG[i].getBboxMin();
    Vec3<REAL_TYPE> b_max = PG[i].getBboxMax();
    
    REAL_TYPE f_min[3];
    REAL_TYPE f_max[3];
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
      REAL_TYPE tmp = origin[q] + region[q];
      if ( f_max[q] > tmp ) f_max[q] = tmp;
    }
    
    Vec3<REAL_TYPE>  rmin(f_min);
    Vec3<REAL_TYPE>  rmax(f_max);
    
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
/**
 * @brief セル数をカウント
 * @param [in] bcd     BCindex B
 * @param [in] m_id    検査するID
 * @param [in] painted m_id以外のセル数をカウント(false), m_idのセル数をカウント(true)
 * @param [in] Dsize   サイズ
 * @note painted : m_id以外のセル数をカウント(false), m_idのセル数をカウント(true)
 */
unsigned long Geometry::countCellB(const int* bcd, const int m_id, const bool painted, const int* Dsize)
{
  int ix, jx, kx, gd;
  
  if ( !Dsize )
  {
    ix = size[0];
    jx = size[1];
    kx = size[2];
    gd = guide;
  }
  else // ASD module用
  {
    ix = Dsize[0];
    jx = Dsize[1];
    kx = Dsize[2];
    gd = 1;
  }
  
  unsigned long c=0;
  int id = m_id;
  
  if ( painted )
  {
#pragma omp parallel for firstprivate(ix, jx, kx, gd, id) schedule(static) reduction(+:c)
    for (int k=1; k<=kx; k++) {
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
          if ( DECODE_CMP(bcd[m]) == id ) c++;
        }
      }
    }
  }
  else
  {
#pragma omp parallel for firstprivate(ix, jx, kx, gd, id) schedule(static) reduction(+:c)
    for (int k=1; k<=kx; k++) {
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
          if ( DECODE_CMP(bcd[m]) != id ) c++;
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
/**
 * @brief セル数をカウント
 * @param [in] mid     work array
 * @param [in] m_id    検査するID
 * @param [in] painted m_id以外のセル数をカウント(false), m_idのセル数をカウント(true)
 * @param [in] Dsize   サイズ
 * @note painted : m_id以外のセル数をカウント(false), m_idのセル数をカウント(true)
 */
unsigned long Geometry::countCellM(const int* mid, const int m_id, const bool painted, const int* Dsize)
{
  int ix, jx, kx, gd;
  
  if ( !Dsize )
  {
    ix = size[0];
    jx = size[1];
    kx = size[2];
    gd = guide;
  }
  else // ASD module用
  {
    ix = Dsize[0];
    jx = Dsize[1];
    kx = Dsize[2];
    gd = 1;
  }
  
  unsigned long c=0;
  int id = m_id;
  
  if ( painted )
  {
#pragma omp parallel for firstprivate(ix, jx, kx, gd, id) schedule(static) reduction(+:c)
    for (int k=1; k<=kx; k++) {
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
          if ( mid[m] == id ) c++;
        }
      }
    }
  }
  else
  {
#pragma omp parallel for firstprivate(ix, jx, kx, gd, id) schedule(static) reduction(+:c)
    for (int k=1; k<=kx; k++) {
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
          if ( mid[m] != id ) c++;
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
/* @brief フィル操作
 * @param [in]      fp     ファイルポインタ
 * @param [in]      cmp    CompoList class
 * @param [in]      mat    MediumList
 * @param [in,out]  d_bcd  BCindex ID
 * @param [in]      d_cut  距離情報
 * @param [in]      d_bid  BC情報
 */
void Geometry::fill(FILE* fp, CompoList* cmp, MediumList* mat, int* d_bcd, float* d_cut, int* d_bid)
{

  // フィル媒質のチェック
  if ( cmp[FillID].getState() != FLUID )
  {
    Hostonly_ printf("\tSpecified Medium of filling fluid is not FLUID\n");
    Exit(0);
  }
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  
  unsigned long target_count; ///< フィルの対象となるセル数
  unsigned long replaced;     ///< 置換された数
  unsigned long filled;       ///< FLUIDでフィルされた数
  unsigned long sum_replaced; ///< 置換された数の合計
  unsigned long sum_filled;   ///< FLUIDでフィルされた数の合計
  

  
  // 最初にフィル対象のセル数を求める >> 全計算内部セル数
  unsigned long total_cell = (unsigned long)ix * (unsigned long)jx * (unsigned long)kx;
  
  
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
  unsigned long fill_cut = fillCutOnCellCenter(d_bcd, d_bid, d_cut);
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
    if ( paraMngr->BndCommS3D(d_bid, ix, jx, kx, gd, gd) != CPM_SUCCESS ) Exit(0);
  }
  
  
  Hostonly_
  {
    printf(    "\n\tFill -----\n\n");
    printf(    "\t\tFilling Fluid Medium   : %s\n", mat[FillID].getAlias().c_str());
    printf(    "\t\tHint of Seeding Dir.   : %s\n", FBUtility::getDirection(FillSeedDir).c_str());
    printf(    "\t\tFill Seed Medium       : %s\n", mat[SeedID].getAlias().c_str());
    printf(    "\t\tFill Control (X, Y, Z) : (%s, %s, %s)\n\n",
           ( !FillSuppress[0] ) ? "Suppress" : "Fill",
           ( !FillSuppress[1] ) ? "Suppress" : "Fill",
           ( !FillSuppress[2] ) ? "Suppress" : "Fill");
    
    fprintf(fp,"\n\tFill -----\n\n");
    fprintf(fp,"\t\tFilling Fluid Medium   : %s\n", mat[FillID].getAlias().c_str());
    fprintf(fp,"\t\tHint of Seeding Dir.   : %s\n", FBUtility::getDirection(FillSeedDir).c_str());
    fprintf(fp,"\t\tFill Seed Medium       : %s\n", mat[SeedID].getAlias().c_str());
    fprintf(fp,"\t\tFill Control (X, Y, Z) : (%s, %s, %s)\n\n",
            ( !FillSuppress[0] ) ? "Suppress" : "Fill",
            ( !FillSuppress[1] ) ? "Suppress" : "Fill",
            ( !FillSuppress[2] ) ? "Suppress" : "Fill");
  }
  
  
  // ヒントが与えられている場合
  filled = fillSeedBcd(d_bcd, FillSeedDir, SeedID, d_bid);
  
  if ( numProc > 1 )
  {
    if ( paraMngr->BndCommS3D(d_bid, ix, jx, kx, gd, gd) != CPM_SUCCESS ) Exit(0);
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
    filled = fillByBid(d_bid, d_bcd, d_cut, SeedID, FillSuppress, fs);
    replaced = fs;
    
    if ( numProc > 1 )
    {
      if ( paraMngr->BndCommS4DEx(d_cut, 6, ix, jx, kx, gd, gd) != CPM_SUCCESS ) Exit(0);
      if ( paraMngr->BndCommS3D(d_bcd, ix, jx, kx, gd, gd) != CPM_SUCCESS ) Exit(0);
      if ( paraMngr->BndCommS3D(d_bid, ix, jx, kx, gd, gd) != CPM_SUCCESS ) Exit(0);
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
    printf(    "\t\t    Filled by [%02d]     = %16ld\n", SeedID, sum_filled);
    fprintf(fp,"\t\t    Filled by [%02d]     = %16ld\n", SeedID, sum_filled);
    printf(    "\t\t    SOLID replaced     = %16ld\n", sum_replaced);
    fprintf(fp,"\t\t    SOLID replaced     = %16ld\n", sum_replaced);
    printf(    "\t\t    Remaining cell     = %16ld\n\n", target_count);
    fprintf(fp,"\t\t    Remaining cell     = %16ld\n\n", target_count);
  }
  
#if 0
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        size_t m = _F_IDX_S3D(i  , j  , k  , ix, jx, kx, gd);
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
  unsigned long upc = countCellB(d_bcd, 0);
  
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
  if ( cmp[SeedID].getState() == FLUID )
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
      replaced = fillByModalSolid(d_bcd, FillID, d_bid);
    }
    else
    {
      replaced = fillByFluid(d_bcd, FillID, d_bid);
    }
    
    if ( numProc > 1 )
    {
      if ( paraMngr->BndCommS3D(d_bcd, ix, jx, kx, gd, gd) != CPM_SUCCESS ) Exit(0);
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
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        size_t m = _F_IDX_S3D(i  , j  , k  , ix, jx, kx, gd);
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
  upc = countCellB(d_bcd, 0);
  
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
/* 
 * @brief 未ペイントセルをフィル
 * @param [in,out] bcd      BCindex B
 * @param [in]     fluid_id フィルをする流体ID
 * @param [in]     bid      境界ID
 * @param [in]     Dsize    サイズ
 * @retval 置換されたセル数
 */
unsigned long Geometry::fillByFluid(int* bcd, const int fluid_id, const int* bid, const int* Dsize)
{
  int ix, jx, kx, gd;
  
  if ( !Dsize )
  {
    ix = size[0];
    jx = size[1];
    kx = size[2];
    gd = guide;
  }
  else // ASD module用
  {
    ix = Dsize[0];
    jx = Dsize[1];
    kx = Dsize[2];
    gd = 1;
  }
  
  int fid = fluid_id;
  unsigned long c = 0; /// painted count
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd, fid) schedule(static) reduction(+:c)
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        
        size_t m_p = _F_IDX_S3D(i  , j  , k  , ix, jx, kx, gd);
        int dd = DECODE_CMP( bcd[m_p] );
        
        // 対象セルが未ペイントの場合
        if ( dd == 0 )
        {
          setBitID(bcd[m_p], fid);
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
/* @brief 未ペイントセルをフィル
 * @param [in,out] mid      work array
 * @param [in]     target   フィルをするID
 * @param [in]     Dsize    サイズ
 * @retval 置換されたセル数
 */
unsigned long Geometry::fillByID(int* mid, const int target, const int* Dsize)
{
  int ix, jx, kx, gd;
  
  if ( !Dsize )
  {
    ix = size[0];
    jx = size[1];
    kx = size[2];
    gd = guide;
  }
  else // ASD module用
  {
    ix = Dsize[0];
    jx = Dsize[1];
    kx = Dsize[2];
    gd = 1;
  }
  
  int fid = target;
  unsigned long c = 0; /// painted count
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd, fid) schedule(static) reduction(+:c)
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        
        if ( mid[m] == -1 )
        {
          mid[m] = fid;
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
/**
 * @brief 流体媒質のフィルをbid情報を元に実行
 * @param [in,out] bid         境界ID（5ビット幅x6方向）
 * @param [in,out] bcd         BCindex B
 * @param [in,out] cut         カット情報
 * @param [in]     tgt_id      フィルする流体IDのエントリ
 * @param [in]     suppress    各軸方向のフィル抑止モード（Periodic, Symmetric時の対策）
 * @param [out]    substituted 固体IDに置換された数
 * @param [in]     Dsize       サイズ
 * @note Symmetric fillにより反復回数を減少
 */
unsigned long Geometry::fillByBid (int* bid, int* bcd, float* cut, const int tgt_id, const int* suppress, unsigned long& substituted, const int* Dsize)
{
  int ix, jx, kx, gd;
  
  if ( !Dsize )
  {
    ix = size[0];
    jx = size[1];
    kx = size[2];
    gd = guide;
  }
  else // ASD module用
  {
    ix = Dsize[0];
    jx = Dsize[1];
    kx = Dsize[2];
    gd = 1;
  }
  
  int tg = tgt_id;
  unsigned long filled   = 0; ///< 流体IDでペイントされた数
  unsigned long replaced = 0; ///< 固体IDで置換された数
  
  // 隣接サブドメインのランク番号
  int sdw = nID[X_minus];
  int sde = nID[X_plus];
  int sds = nID[Y_minus];
  int sdn = nID[Y_plus];
  int sdb = nID[Z_minus];
  int sdt = nID[Z_plus];
  
  int mode_x = suppress[0]; // if 0, suppress connectivity evaluation
  int mode_y = suppress[1];
  int mode_z = suppress[2];
  
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd, tg, sdw, sde, sds, sdn, sdb, sdt, mode_x, mode_y, mode_z) \
schedule(static) reduction(+:filled) reduction(+:replaced)
  
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        
#include "fill_bid.h"
        
      }
    }
  }
  
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd, tg, sdw, sde, sds, sdn, sdb, sdt, mode_x, mode_y, mode_z) \
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
/**
 * @brief 流体媒質のフィルをbid情報を元に実行
 * @param [in,out] mid         work array
 * @param [in]     tgt_id      フィルする流体IDのエントリ
 * @param [in]     Dsize       サイズ
 * @note Symmetric fillにより反復回数を減少
 */
unsigned long Geometry::fillByMid(int* mid, const int tgt_id, const int* Dsize)
{
  int ix, jx, kx, gd;
  
  if ( !Dsize )
  {
    ix = size[0];
    jx = size[1];
    kx = size[2];
    gd = guide;
  }
  else // ASD module用
  {
    ix = Dsize[0];
    jx = Dsize[1];
    kx = Dsize[2];
    gd = 1;
  }
  
  int tg = tgt_id;
  unsigned long filled = 0; ///< tgt_idでペイントされた数
  
  
  /// 検査対象セル{-1}の隣接6方向を見て、tgt_idと同じIDがあれば対象セルをtgt_idでペイント
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd, tg) \
schedule(static) reduction(+:filled) collapse(3)
  
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
        
        if ( mid[m_p] == -1 )
        {
          int ff = 0;
          
          if ( mid[m_w] == tg ) ff++;
          if ( mid[m_e] == tg ) ff++;
          if ( mid[m_s] == tg ) ff++;
          if ( mid[m_n] == tg ) ff++;
          if ( mid[m_b] == tg ) ff++;
          if ( mid[m_t] == tg ) ff++;
          
          if ( ff>0 )
          {
            mid[m_p] = tg;
            filled++;
          }
        }
      }
    }
  }
  
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd, tg) \
schedule(static) reduction(+:filled) collapse(3)
  
  for (int k=kx; k>=1; k--) {
    for (int j=jx; j>=1; j--) {
      for (int i=ix; i>=1; i--) {
        
        size_t m_p = _F_IDX_S3D(i  , j  , k  , ix, jx, kx, gd);
        size_t m_e = _F_IDX_S3D(i+1, j,   k,   ix, jx, kx, gd);
        size_t m_w = _F_IDX_S3D(i-1, j,   k,   ix, jx, kx, gd);
        size_t m_n = _F_IDX_S3D(i,   j+1, k,   ix, jx, kx, gd);
        size_t m_s = _F_IDX_S3D(i,   j-1, k,   ix, jx, kx, gd);
        size_t m_t = _F_IDX_S3D(i,   j,   k+1, ix, jx, kx, gd);
        size_t m_b = _F_IDX_S3D(i,   j,   k-1, ix, jx, kx, gd);
        
        if ( mid[m_p] == -1 )
        {
          int ff = 0;
          
          if ( mid[m_w] == tg ) ff++;
          if ( mid[m_e] == tg ) ff++;
          if ( mid[m_s] == tg ) ff++;
          if ( mid[m_n] == tg ) ff++;
          if ( mid[m_b] == tg ) ff++;
          if ( mid[m_t] == tg ) ff++;
          
          if ( ff>0 )
          {
            mid[m_p] = tg;
            filled++;
          }
        }
        
      }
    }
  }
  
  if ( numProc > 1 )
  {
    unsigned long tmp = filled;
    if ( paraMngr->Allreduce(&tmp, &filled, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
  }
  
  return filled;
}


// #################################################################
/* @brief 未ペイントセルを周囲の媒質IDの固体最頻値でフィル
 * @param [in,out] bcd      BCindex B
 * @param [in]     fluid_id フィルをする流体ID
 * @param [in]     bid      境界ID
 * @retval 置換されたセル数
 * @note 周囲の媒質IDの固体最頻値がゼロの場合には，境界IDで代用
 */
unsigned long Geometry::fillByModalSolid(int* bcd, const int fluid_id, const int* bid)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  int fid = fluid_id;
  unsigned long c = 0; /// painted count
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd, fid) schedule(static) reduction(+:c)
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
        
        int dd = DECODE_CMP( bcd[m_p] );
        
        int qw = DECODE_CMP( bcd[m_w] );
        int qe = DECODE_CMP( bcd[m_e] );
        int qs = DECODE_CMP( bcd[m_s] );
        int qn = DECODE_CMP( bcd[m_n] );
        int qb = DECODE_CMP( bcd[m_b] );
        int qt = DECODE_CMP( bcd[m_t] );
        
        
        // 対象セルが未ペイントの場合
        if ( dd == 0 )
        {
          int sd = FBUtility::find_mode_id(fid, qw, qe, qs, qn, qb, qt, m_NoCompo);
          
          // 周囲の媒質IDの固体最頻値がゼロの場合
          if ( sd == 0 )
          {
            int qq = bid[m_p];
            qw = getBit5(qq, 0);
            qe = getBit5(qq, 1);
            qs = getBit5(qq, 2);
            qn = getBit5(qq, 3);
            qb = getBit5(qq, 4);
            qt = getBit5(qq, 5);
            sd = FBUtility::find_mode_id(fid, qw, qe, qs, qn, qb, qt, m_NoCompo);
            if ( sd == 0 ) Exit(0); // 何かあるはず
          }
          setBitID(bcd[m_p], sd);
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
/**
 * @brief 交点が定義点にある場合にそのポリゴンのエントリ番号でフィルする
 * @param [in,out] bcd    BCindex B
 * @param [in]     bid    境界ID（5ビット幅x6方向）
 * @param [in]     cut    カット情報
 * @param [in]     Dsize  サイズ
 * @retval フィルされたセル数
 */
unsigned long Geometry::fillCutOnCellCenter(int* bcd, const int* bid, const float* cut, const int* Dsize)
{
  int ix, jx, kx, gd;
  
  if ( !Dsize )
  {
    ix = size[0];
    jx = size[1];
    kx = size[2];
    gd = guide;
  }
  else // ASD module用
  {
    ix = Dsize[0];
    jx = Dsize[1];
    kx = Dsize[2];
    gd = 1;
  }
  
  unsigned long c = 0;
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static) reduction(+:c)
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
        
        int qq = bid[m_p];
        
        // 隣接セルの方向に対する境界ID
        int qw = getBit5(qq, X_minus);
        int qe = getBit5(qq, X_plus);
        int qs = getBit5(qq, Y_minus);
        int qn = getBit5(qq, Y_plus);
        int qb = getBit5(qq, Z_minus);
        int qt = getBit5(qq, Z_plus);
        
        const float* pos = &cut[ _F_IDX_S4DEX(0, i, j, k, 6, ix, jx, kx, gd) ];
        
        // いずれかの方向で交点が定義点上の場合
        if ( pos[0]*pos[1]*pos[2]*pos[3]*pos[4]*pos[5] == 0.0 )
        {
          //printf("%d %d %d : %f %f %f %f %f %f : %d\n",i,j,k,pos[0],pos[1],pos[2],pos[3],pos[4],pos[5], qw);
          setBitID(bcd[m_p], qw); // qwで代表
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
/**
 * @brief シード点をbcd[]にペイントする
 * @param [in,out] bcd    BCindex B
 * @param [in]     face   ヒント面
 * @param [in]     target ペイントするIDのエントリ
 * @param [in]     bid    境界ID
 * @param [in]     Dsize  サイズ
 * @note ヒントとして与えられた外部境界面に接するセルにおいて，確実に流体セルであるセルをフィルする
 *       もし，外部境界面以外に固体候補があれば、ぬれ面はフィルしない
 */
unsigned long Geometry::fillSeedBcd(int* bcd, const int face, const int target, const int* bid, const int* Dsize)
{
  int ix, jx, kx, gd;
  
  if ( !Dsize )
  {
    ix = size[0];
    jx = size[1];
    kx = size[2];
    gd = guide;
  }
  else // ASD module用
  {
    ix = Dsize[0];
    jx = Dsize[1];
    kx = Dsize[2];
    gd = 1;
  }
  
  
  int tg = target;     ///< order of FLUID ID
  unsigned long c = 0;
  
  // 各外部境界面以外のフェイスにカットがあるかどうかを検査するためのマスク
  unsigned mask0 = 0x3fffffff;
  unsigned mx_w  = mask0 & ( ~(0x1f << (X_minus*5)) );
  unsigned mx_e  = mask0 & ( ~(0x1f << (X_plus *5)) );
  unsigned mx_s  = mask0 & ( ~(0x1f << (Y_minus*5)) );
  unsigned mx_n  = mask0 & ( ~(0x1f << (Y_plus *5)) );
  unsigned mx_b  = mask0 & ( ~(0x1f << (Z_minus*5)) );
  unsigned mx_t  = mask0 & ( ~(0x1f << (Z_plus *5)) );
  
  
  switch (face)
  {
    case X_minus:
      if ( nID[face] < 0 )
      {
#pragma omp parallel for firstprivate(ix, jx, kx, gd, tg, mx_w) schedule(static) reduction(+:c)
        for (int k=1; k<=kx; k++) {
          for (int j=1; j<=jx; j++) {
            size_t m = _F_IDX_S3D(1, j, k, ix, jx, kx, gd);
            int s = bid[m];
            
            // 対象セルの周囲に指定方向以外にカットがなく，未ペイントのセルの場合
            if ( ((s & mx_w) == 0) && (DECODE_CMP(bcd[m]) == 0) )
            {
              setBitID(bcd[m], tg);
              c++;
            }
          }
        }
      }
      break;
      
    case X_plus:
      if ( nID[face] < 0 )
      {
#pragma omp parallel for firstprivate(ix, jx, kx, gd, tg, mx_e) schedule(static) reduction(+:c)
        for (int k=1; k<=kx; k++) {
          for (int j=1; j<=jx; j++) {
            size_t m = _F_IDX_S3D(ix, j, k, ix, jx, kx, gd);
            int s = bid[m];
            if ( ((s & mx_e) == 0) && (DECODE_CMP(bcd[m]) == 0) )
            {
              setBitID(bcd[m], tg);
              c++;
            }
          }
        }
      }
      break;
      
    case Y_minus:
      if ( nID[face] < 0 )
      {
#pragma omp parallel for firstprivate(ix, jx, kx, gd, tg, mx_s) schedule(static) reduction(+:c)
        for (int k=1; k<=kx; k++) {
          for (int i=1; i<=ix; i++) {
            size_t m = _F_IDX_S3D(i, 1, k, ix, jx, kx, gd);
            int s = bid[m];
            
            if ( ((s & mx_s) == 0) && (DECODE_CMP(bcd[m]) == 0) )
            {
              setBitID(bcd[m], tg);
              c++;
            }
          }
        }
      }
      break;
      
    case Y_plus:
      if ( nID[face] < 0 )
      {
#pragma omp parallel for firstprivate(ix, jx, kx, gd, tg, mx_n) schedule(static) reduction(+:c)
        for (int k=1; k<=kx; k++) {
          for (int i=1; i<=ix; i++) {
            size_t m = _F_IDX_S3D(i, jx, k, ix, jx, kx, gd);
            int s = bid[m];
            
            if ( ((s & mx_n) == 0) && (DECODE_CMP(bcd[m]) == 0) )
            {
              setBitID(bcd[m], tg);
              c++;
            }
          }
        }
      }
      break;
      
    case Z_minus:
      if ( nID[face] < 0 )
      {
#pragma omp parallel for firstprivate(ix, jx, kx, gd, tg, mx_b) schedule(static) reduction(+:c)
        for (int j=1; j<=jx; j++) {
          for (int i=1; i<=ix; i++) {
            size_t m = _F_IDX_S3D(i, j, 1, ix, jx, kx, gd);
            int s = bid[m];
            
            if ( ((s & mx_b) == 0) && (DECODE_CMP(bcd[m]) == 0) )
            {
              setBitID(bcd[m], tg);
              c++;
            }
          }
        }
      }
      break;
      
    case Z_plus:
      if ( nID[face] < 0 )
      {
#pragma omp parallel for firstprivate(ix, jx, kx, gd, tg, mx_t) schedule(static) reduction(+:c)
        for (int j=1; j<=jx; j++) {
          for (int i=1; i<=ix; i++) {
            size_t m = _F_IDX_S3D(i, j, kx, ix, jx, kx, gd);
            int s = bid[m];
            
            if ( ((s & mx_t) == 0) && (DECODE_CMP(bcd[m]) == 0) )
            {
              setBitID(bcd[m], tg);
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
/**
 * @brief シード点をmid[]にペイントする
 * @param [in,out] mid    work array
 * @param [in]     face   ヒント面
 * @param [in]     target ペイントするIDのエントリ
 * @param [in]     Dsize  サイズ
 * @note ヒントとして与えられた外部境界面に接するセルにおいて，-1のセルをtargetでペイントする
 */
unsigned long Geometry::fillSeedMid(int* mid, const int face, const int target, const int* Dsize)
{
  int ix, jx, kx, gd;
  
  if ( !Dsize )
  {
    ix = size[0];
    jx = size[1];
    kx = size[2];
    gd = guide;
  }
  else // ASD module用
  {
    ix = Dsize[0];
    jx = Dsize[1];
    kx = Dsize[2];
    gd = 1;
  }
  
  
  int tg = target;
  unsigned long c = 0;
  
  
  switch (face)
  {
    case X_minus:
      if ( nID[face] < 0 )
      {
#pragma omp parallel for firstprivate(ix, jx, kx, gd, tg) schedule(static) reduction(+:c)
        for (int k=1; k<=kx; k++) {
          for (int j=1; j<=jx; j++) {
            size_t m = _F_IDX_S3D(1, j, k, ix, jx, kx, gd);
            
            if ( mid[m] == -1 )
            {
              mid[m] = tg;
              c++;
            }
          }
        }
      }
      break;
      
    case X_plus:
      if ( nID[face] < 0 )
      {
#pragma omp parallel for firstprivate(ix, jx, kx, gd, tg) schedule(static) reduction(+:c)
        for (int k=1; k<=kx; k++) {
          for (int j=1; j<=jx; j++) {
            size_t m = _F_IDX_S3D(ix, j, k, ix, jx, kx, gd);
            
            if ( mid[m] == -1 )
            {
              mid[m] = tg;
              c++;
            }
          }
        }
      }
      break;
      
    case Y_minus:
      if ( nID[face] < 0 )
      {
#pragma omp parallel for firstprivate(ix, jx, kx, gd, tg) schedule(static) reduction(+:c)
        for (int k=1; k<=kx; k++) {
          for (int i=1; i<=ix; i++) {
            size_t m = _F_IDX_S3D(i, 1, k, ix, jx, kx, gd);
            
            if ( mid[m] == -1 )
            {
              mid[m] = tg;
              c++;
            }
          }
        }
      }
      break;
      
    case Y_plus:
      if ( nID[face] < 0 )
      {
#pragma omp parallel for firstprivate(ix, jx, kx, gd, tg) schedule(static) reduction(+:c)
        for (int k=1; k<=kx; k++) {
          for (int i=1; i<=ix; i++) {
            size_t m = _F_IDX_S3D(i, jx, k, ix, jx, kx, gd);
            
            if ( mid[m] == -1 )
            {
              mid[m] = tg;
              c++;
            }
          }
        }
      }
      break;
      
    case Z_minus:
      if ( nID[face] < 0 )
      {
#pragma omp parallel for firstprivate(ix, jx, kx, gd, tg) schedule(static) reduction(+:c)
        for (int j=1; j<=jx; j++) {
          for (int i=1; i<=ix; i++) {
            size_t m = _F_IDX_S3D(i, j, 1, ix, jx, kx, gd);
            
            if ( mid[m] == -1 )
            {
              mid[m] = tg;
              c++;
            }
          }
        }
      }
      break;
      
    case Z_plus:
      if ( nID[face] < 0 )
      {
#pragma omp parallel for firstprivate(ix, jx, kx, gd, tg) schedule(static) reduction(+:c)
        for (int j=1; j<=jx; j++) {
          for (int i=1; i<=ix; i++) {
            size_t m = _F_IDX_S3D(i, j, kx, ix, jx, kx, gd);
            
            if ( mid[m] == -1 )
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
/**
 * @brief list[]内の最頻値IDを求める
 * @param [in] m_sz  配列のサイズ
 * @param [in] list  ID配列
 * @note 候補がない場合には、0が戻り値
 */
int Geometry::find_mode(const int m_sz, const int* list)
{
  int key[CMP_BIT_W]; ///< ID毎の頻度 @note ffv_Initialize() >> fill()でif ( NoCompo+1 > CMP_BIT_W )をチェック
  memset(key, 0, sizeof(int)*CMP_BIT_W);
  
  
  for (int l=0; l<m_sz; l++) key[ list[l] ]++;
  
  
  int mode = 0; // サーチの初期値，IDの大きい方から
  int z = 0;    // 最頻値のID
  
  for (int l=m_NoCompo; l>=1; l--)
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
// @param [in,out] d_mid  識別子配列
unsigned long Geometry::findPolygonInCell(int* d_mid)
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
    string m_pg = (*it)->get_name();     // グループラベル
    string m_bc = (*it)->get_type();     // 境界条件ラベル
    int ntria = (*it)->get_group_num_tria();  // ローカルのポリゴン数
    
#if 0
    Hostonly_ printf("\n%s : %s\n", m_pg.c_str(), m_bc.c_str() );
#endif
    
    printf("\n");
    
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
              
              Vec3<REAL_TYPE> bx_min(m_poly_org[0]+m_poly_dx[0]*(REAL_TYPE)(i-1),
                                     m_poly_org[1]+m_poly_dx[1]*(REAL_TYPE)(j-1),
                                     m_poly_org[2]+m_poly_dx[2]*(REAL_TYPE)(k-1)); // セルBboxの対角座標
              Vec3<REAL_TYPE> bx_max(m_poly_org[0]+m_poly_dx[0]*(REAL_TYPE)i,
                                     m_poly_org[1]+m_poly_dx[1]*(REAL_TYPE)j,
                                     m_poly_org[2]+m_poly_dx[2]*(REAL_TYPE)k);     // セルBboxの対角座標
              
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
                
                int z = find_mode(polys, ary);
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
  unsigned long c = countCellM(d_mid, 0, true);
  
  if ( c != 0 )
  {
    Hostonly_ stamped_printf("\tID=0 was found in water-tight fill process. : %d\n", c);
    Exit(0);
  }
  
  
  // count the number of replaced cells >> except ID=-1
  c = countCellM(d_mid, -1, false);
  
  return c;
}


// #################################################################
// フィルパラメータを取得
void Geometry::getFillParam(TextParser* tpCntl)
{
  string str;
  string label;
  
  // フィルの媒質指定
  label = "/GeometryModel/FillMedium";
  if ( !(tpCntl->getInspectedValue(label, str)) )
  {
    Hostonly_ printf("\tParsing error in '%s'\n", label.c_str());
    Exit(0);
  }
  FillMedium = str;
  
  
  // 流体セルのフィルの開始面指定
  label = "/GeometryModel/HintOfFillSeedDirection";
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : Invalid value in '%s'\n", label.c_str());
    Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "xminus" ) ) FillSeedDir = X_minus;
  else if( !strcasecmp(str.c_str(), "xplus"  ) ) FillSeedDir = X_plus;
  else if( !strcasecmp(str.c_str(), "yminus" ) ) FillSeedDir = Y_minus;
  else if( !strcasecmp(str.c_str(), "yplus"  ) ) FillSeedDir = Y_plus;
  else if( !strcasecmp(str.c_str(), "zminus" ) ) FillSeedDir = Z_minus;
  else if( !strcasecmp(str.c_str(), "zplus"  ) ) FillSeedDir = Z_plus;
  else
  {
    FillSeedDir = X_minus;
    Hostonly_ printf("\tDefault 'X_minus' is set for Hint Of FillSeed direction\n");
  }
  
  
  // ヒントに使うフィルの媒質指定
  label = "/GeometryModel/HintOfFillSeedMedium";
  if ( !(tpCntl->getInspectedValue(label, str)) )
  {
    Hostonly_ printf("\tParsing error in '%s'\n", label.c_str());
    Exit(0);
  }
  SeedMedium = str;
  
  
  // フィル方向制御 (NOT mandatory)
  string dir[3];
  label = "/GeometryModel/FillDirectionControl";
  if ( tpCntl->chkLabel(label) )
  {
    if ( !(tpCntl->getInspectedVector(label, dir, 3)) )
    {
      Hostonly_ printf("\tParsing error in '%s'\n", label.c_str());
      Exit(0);
    }
    else
    {
      FillSuppress[0] = ( !strcasecmp(dir[0].c_str(), "fill" ) ) ? 1 : 0;
      FillSuppress[1] = ( !strcasecmp(dir[1].c_str(), "fill" ) ) ? 1 : 0;
      FillSuppress[2] = ( !strcasecmp(dir[2].c_str(), "fill" ) ) ? 1 : 0;
    }
  }
  
}


// #################################################################
void Geometry::Initialize(MPIPolylib* PL,
                          PolygonProperty* PG,
                          const REAL_TYPE* poly_org,
                          const REAL_TYPE* poly_dx,
                          const int NoCompo,
                          const int NoPolyGrp)
{
  this->PL = PL;
  this->PG = PG;
  
  for ( int i=0; i<3; i++)
  {
    m_poly_org[i] = poly_org[i];
    m_poly_dx[i]  = poly_dx[i];
  }

  m_NoCompo   = NoCompo;
  m_NoPolyGrp = NoPolyGrp;
}



// #################################################################
/* @brief サブセルID配列をrefIDでペイント
 * @param [in,out] svf       サブセルの体積率
 * @param [in,out] smd       サブセルID配列
 * @param [in]     dir       ペイント開始方向
 * @param [in]     refID     ペイントID
 * @param [in]     refVf     ペイントする体積率
 * @note single threadで実行
 */
int Geometry::SubCellFill(REAL_TYPE* svf,
                     int* smd,
                     const int dir,
                     const int refID,
                     const REAL_TYPE refVf)
{
  int filled = 0;
  int sdv = NumSuvDiv; // OpenMPのfirstprivateで使うためローカル変数にコピー
  
  if ( dir == X_minus )
  {
    
#pragma omp parallel for firstprivate(sdv, refID, refVf) reduction(+:filled) schedule(static) collapse(3)
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
    
#pragma omp parallel for firstprivate(sdv, refID, refVf) reduction(+:filled) schedule(static) collapse(3)
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
    
#pragma omp parallel for firstprivate(sdv, refID, refVf) reduction(+:filled) schedule(static) collapse(3)
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

#pragma omp parallel for firstprivate(sdv, refID, refVf) reduction(+:filled) schedule(static) collapse(3)
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
    
#pragma omp parallel for firstprivate(sdv, refID, refVf) reduction(+:filled) schedule(static) collapse(3)
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
    
#pragma omp parallel for firstprivate(sdv, refID, refVf) reduction(+:filled) schedule(static) collapse(3)
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
 * @param [in]     ip,jp,kp  プライマリセルインデクス
 * @param [in]     pch       サブセルの格子幅
 * @param [in]     m_pg      ポリゴングループ名
 * @retval ポリゴンを含むセル数
 * @note 呼び出し元がスレッド並列の場合、single threadで実行
 */
int Geometry::SubCellIncTest(REAL_TYPE* svf,
                             int* smd,
                             const int ip,
                             const int jp,
                             const int kp,
                             const Vec3<REAL_TYPE> pch,
                             const string m_pg)
{
  // プライマリセルの基点（有次元）
  Vec3<REAL_TYPE> o(m_poly_org[0]+m_poly_dx[0]*(REAL_TYPE)(ip-1),
                    m_poly_org[1]+m_poly_dx[1]*(REAL_TYPE)(jp-1),
                    m_poly_org[2]+m_poly_dx[2]*(REAL_TYPE)(kp-1));
  
  int pic = 0;
  
  int sdv = NumSuvDiv; // OpenMPのfirstprivateで使うためローカル変数にコピー
  
  // サブセルの体積率を評価（ポリゴンをもつサブセルのみ0.5）
  for (int k=1; k<=sdv; k++) {
    for (int j=1; j<=sdv; j++) {
      for (int i=1; i<=sdv; i++) {
        
        // サブセルのbbox
        Vec3<REAL_TYPE> b1(o.x + pch.x * (REAL_TYPE)(i-1),
                           o.y + pch.y * (REAL_TYPE)(j-1),
                           o.z + pch.z * (REAL_TYPE)(k-1));
        Vec3<REAL_TYPE> b2(o.x + pch.x * (REAL_TYPE)i,
                           o.y + pch.y * (REAL_TYPE)j,
                           o.z + pch.z * (REAL_TYPE)k);
        
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
          
          int z = find_mode(polys, ary);
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
 * @param [in]     ip,jp,kp  プライマリセルインデクス
 * @param [in]     d_mid     識別子配列
 * @param [in]     mat       MediumList
 * @param [out]    d_pvf     体積率
 * @note 呼び出し先でスレッド化している場合には、single threadで実行
 */
void Geometry::SubDivision(REAL_TYPE* svf,
                           int* smd,
                           const int ip,
                           const int jp,
                           const int kp,
                           const int* d_mid,
                           const MediumList* mat,
                           REAL_TYPE* d_pvf)
{
  // 外縁部にポリゴンがない面を探す
  int face_flag = 0;
  int c;
  int sum = 0;

  int sdv = NumSuvDiv; // OpenMPのfirstprivateで使うためローカル変数にコピー
  
  // X-
  c = 0;
  
#pragma omp parallel for firstprivate(sdv) reduction(+:c) schedule(static) collapse(2)
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
  
#pragma omp parallel for firstprivate(sdv) reduction(+:c) schedule(static) collapse(2)
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
  
#pragma omp parallel for firstprivate(sdv) reduction(+:c) schedule(static) collapse(2)
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
  
#pragma omp parallel for firstprivate(sdv) reduction(+:c) schedule(static) collapse(2)
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
  
#pragma omp parallel for firstprivate(sdv) reduction(+:c) schedule(static) collapse(2)
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
  
#pragma omp parallel for firstprivate(sdv) reduction(+:c) schedule(static) collapse(2)
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
    printf("%d %d\n", SeedID, FillID);
    if ( TEST_BIT(face_flag, X_minus) )
    {
      
      m = _F_IDX_S3D(ip-1, jp, kp, ix, jx, kx, gd);
      r = d_mid[m];
      if ( (r == SeedID) || (r == FillID) )
      {
        refID = r;
        fillDir = 0;
      }
    }
    else if ( TEST_BIT(face_flag, X_plus) )
    {
      m = _F_IDX_S3D(ip+1, jp, kp, ix, jx, kx, gd);
      r = d_mid[m];
      if ( (r == SeedID) || (r == FillID) )
      {
        refID = r;
        fillDir = 1;
      }
    }
    else if ( TEST_BIT(face_flag, Y_minus) )
    {
      m = _F_IDX_S3D(ip, jp-1, kp, ix, jx, kx, gd);
      r = d_mid[m];
      if ( (r == SeedID) || (r == FillID) )
      {
        refID = r;
        fillDir = 2;
      }
    }
    else if ( TEST_BIT(face_flag, Y_plus) )
    {
      m = _F_IDX_S3D(ip, jp+1, kp, ix, jx, kx, gd);
      r = d_mid[m];
      if ( (r == SeedID) || (r == FillID) )
      {
        refID = r;
        fillDir = 3;
      }
    }
    else if ( TEST_BIT(face_flag, Z_minus) )
    {
      m = _F_IDX_S3D(ip, jp, kp-1, ix, jx, kx, gd);
      r = d_mid[m];
      if ( (r == SeedID) || (r == FillID) )
      {
        refID = r;
        fillDir = 4;
      }
    }
    else if ( TEST_BIT(face_flag, Z_plus) )
    {
      m = _F_IDX_S3D(ip, jp, kp+1, ix, jx, kx, gd);
      r = d_mid[m];
      if ( (r == SeedID) || (r == FillID) )
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
    
#pragma omp parallel for firstprivate(sdv, refID, jx, kx) reduction(+:c) schedule(static) collapse(2)
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
    
#pragma omp parallel for firstprivate(sdv, refID, jx, kx) reduction(+:c) schedule(static) collapse(2)
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
    
#pragma omp parallel for firstprivate(sdv, refID, ix, kx) reduction(+:c) schedule(static) collapse(2)
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

#pragma omp parallel for firstprivate(sdv, refID, ix, kx) reduction(+:c) schedule(static) collapse(2)
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
    
#pragma omp parallel for firstprivate(sdv, refID, ix, jx) reduction(+:c) schedule(static) collapse(2)
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
    
#pragma omp parallel for firstprivate(sdv, refID, ix, jx) reduction(+:c) schedule(static) collapse(2)
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
  REAL_TYPE refVf = (mat[refID].getState()==FLUID) ? 1.0 : 0.0; ///< refIDの体積率
  
  c = 0;
  while (target_count > 0) {
    
    filled = SubCellFill(svf, smd, fillDir, refID, refVf);
    
    target_count -= filled;
    sum_filled   += filled;
    c++;
    printf("\titr=%3d %d\n", c+1, filled);
    
    if ( filled <= 0 ) break; // フィル対象がなくなったら終了
  }
  
  // 未ペイント(-1)をチェック
  int flag = 0;
  
#pragma omp parallel for firstprivate(sdv) reduction(+:flag) schedule(static) collapse(3)
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
    int paintID = (refID == SeedID) ? FillID : SeedID;
    REAL_TYPE paintVf = (mat[refID].getState()==FLUID) ? 0.0 : 1.0;
    
#pragma omp parallel for firstprivate(sdv, paintID, paintVf) schedule(static) collapse(3)
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
  REAL_TYPE sff = 0.0;
  
#pragma omp parallel for firstprivate(sdv) reduction(+:sff) schedule(static) collapse(3)
  for (int k=1; k<=sdv; k++) {
    for (int j=1; j<=sdv; j++) {
      for (int i=1; i<=sdv; i++) {
        size_t m = _F_IDX_S3D(i, j, k, sdv, sdv, sdv, 1);
        sff += svf[m];
      }
    }
  }

  REAL_TYPE ff = 1.0/(REAL_TYPE)(sdv*sdv*sdv);
  d_pvf[_F_IDX_S3D(ip, jp, kp, ix, jx, kx, gd)] = sff*ff;
  
}



// #################################################################
// サブサンプリング
void Geometry::SubSampling(FILE* fp, MediumList* mat, int* d_mid, REAL_TYPE* d_pvf)
{
  unsigned long target_count = 0; ///< フィルの対象となるセル数
  unsigned long replaced = 0;     ///< 置換された数
  unsigned long filled = 0;       ///< フィルされた数
  unsigned long sum_replaced = 0; ///< 置換された数の合計
  unsigned long sum_filled = 0;   ///< FLUIDでフィルされた数の合計
  

  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
   int sdv = NumSuvDiv; // OpenMPのfirstprivateで使うためローカル変数にコピー
  
  
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
  REAL_TYPE tmp = (mat[SeedID].getState()==FLUID) ? 1.0 : 0.0;
  filled = assignVF(SeedID, tmp, d_mid, d_pvf);
  
  Hostonly_
  {
    printf    ("\t\tOuter assigned Vf(%3.1f) = %16ld  (%s)\n", tmp, filled, mat[SeedID].getAlias().c_str());
    fprintf(fp,"\t\tOuter assigned Vf(%3.1f) = %16ld  (%s)\n", tmp, filled, mat[SeedID].getAlias().c_str());
  }
  
  // Inner fill => FillID
  tmp = (mat[FillID].getState()==FLUID) ? 1.0 : 0.0;
  filled = assignVF(FillID, tmp, d_mid, d_pvf);
  
  Hostonly_
  {
    printf    ("\t\tInner assigned Vf(%3.1f) = %16ld  (%s)\n", tmp, filled, mat[FillID].getAlias().c_str());
    fprintf(fp,"\t\tInner assigned Vf(%3.1f) = %16ld  (%s)\n", tmp, filled, mat[FillID].getAlias().c_str());
  }
  
  
  // Local
  REAL_TYPE m_org[3], m_pit[3];
  int m_siz[3];
  
  for (int i=0; i<3; i++)
  {
    m_siz[i] = sdv;
  }

  
  
  
  // sub-cell work array
  size_t nx = (sdv+2)*(sdv+2)*(sdv+2);
  REAL_TYPE* svf = new REAL_TYPE [nx]; // guide cell is 1 for each dir.
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
  // この三重ループはスレッド化しない >> スレッド化する場合には、ループ内で呼び出しているメソッドをチェックすること
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        
#pragma omp parallel for firstprivate(sdv) schedule(static)
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
        for (int l=0; l<3; l++) m_pit[l] = m_poly_dx[l]/(REAL_TYPE)sdv;
        
        m_org[0] = m_poly_org[0]+m_poly_dx[0]*(REAL_TYPE)(i-1);
        m_org[1] = m_poly_org[1]+m_poly_dx[1]*(REAL_TYPE)(j-1);
        m_org[2] = m_poly_org[2]+m_poly_dx[2]*(REAL_TYPE)(k-1);
        
        
        
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        int q = d_mid[m];
        
        // PolyIDsの場合のみ
        if ( (q != FillID) && (q != SeedID) )
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
                int cp = SubCellIncTest(svf, smd, i, j, k, m_pit, m_pg);
                printf("%3d %3d %3d : %3d\n", i,j,k,cp);
              }
            }
          } // Polygon Group
          
          //writeSVX(smd, i, j, k, m_siz, 1, m_pit, m_org);
          
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
          SubDivision(svf, smd, i, j, k, d_mid, mat, d_pvf);
          
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
  
  //F->writeRawSPH(d_pvf, size, gd, 0, m_org, m_pit, sizeof(REAL_TYPE));

}



// #################################################################
// ポリゴンの水密化
// @note ここまで、d_bcdにはsetMonitorList()でモニタIDが入っている
void Geometry::WaterTightening(FILE* fp, CompoList* cmp, MediumList* mat, int* d_mid)
{
  unsigned long target_count = 0; ///< フィルの対象となるセル数
  unsigned long replaced = 0;     ///< 置換された数
  unsigned long filled = 0;       ///< FLUIDでフィルされた数
  unsigned long sum_replaced = 0; ///< 置換された数の合計
  unsigned long sum_filled = 0;   ///< FLUIDでフィルされた数の合計
  
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  
  // フィル媒質のチェック
  if ( cmp[FillID].getState() != FLUID )
  {
    Hostonly_ printf("\tSpecified Medium of filling fluid is not FLUID\n");
    Exit(0);
  }
  
  
  // 最初にフィル対象のセル数を求める >> 全計算内部セル数
  unsigned long total_cell = (unsigned long)ix * (unsigned long)jx * (unsigned long)kx;
  
  
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
    
    printf(    "\t\tFilling Fluid Medium   : %s\n", mat[FillID].getAlias().c_str());
    printf(    "\t\tHint of Seeding Dir.   : %s\n", FBUtility::getDirection(FillSeedDir).c_str());
    printf(    "\t\tFill Seed Medium       : %s\n", mat[SeedID].getAlias().c_str());
    printf(    "\t\tFill Control (X, Y, Z) : (%s, %s, %s)\n\n",
           ( !FillSuppress[0] ) ? "Suppress" : "Fill",
           ( !FillSuppress[1] ) ? "Suppress" : "Fill",
           ( !FillSuppress[2] ) ? "Suppress" : "Fill");
    
    fprintf(fp,"\t\tFilling Fluid Medium   : %s\n", mat[FillID].getAlias().c_str());
    fprintf(fp,"\t\tHint of Seeding Dir.   : %s\n", FBUtility::getDirection(FillSeedDir).c_str());
    fprintf(fp,"\t\tFill Seed Medium       : %s\n", mat[SeedID].getAlias().c_str());
    fprintf(fp,"\t\tFill Control (X, Y, Z) : (%s, %s, %s)\n\n",
            ( !FillSuppress[0] ) ? "Suppress" : "Fill",
            ( !FillSuppress[1] ) ? "Suppress" : "Fill",
            ( !FillSuppress[2] ) ? "Suppress" : "Fill");
  }
  
  
  
  
  // セルに含まれるポリゴンを探索し、d_midに記録
  Hostonly_
  {
    printf(    "\tPaint cells that contain polygons -----\n\n");
    fprintf(fp,"\tPaint cells that contain polygons -----\n\n");
  }
  
  sum_replaced = findPolygonInCell(d_mid);
  
  target_count -= sum_replaced;
  
  if ( numProc > 1 )
  {
    if ( paraMngr->BndCommS3D(d_mid, ix, jx, kx, gd, gd) != CPM_SUCCESS ) Exit(0);
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
  
  filled = fillSeedMid(d_mid, FillSeedDir, SeedID);
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
    if ( paraMngr->BndCommS3D(d_mid, ix, jx, kx, gd, gd) != CPM_SUCCESS ) Exit(0);
  }
  
  
  Hostonly_ // Up to here, d_mid={-1, Poly-IDs, SeedID}
  {
    printf(    "\t\tPainted cells          = %16ld  (%s)\n", filled, mat[SeedID].getAlias().c_str());
    fprintf(fp,"\t\tPainted cells          = %16ld  (%s)\n", filled, mat[SeedID].getAlias().c_str());
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
    filled = fillByMid(d_mid, SeedID);
    
    if ( numProc > 1 )
    {
      if ( paraMngr->BndCommS3D(d_mid, ix, jx, kx, gd, gd) != CPM_SUCCESS ) Exit(0);
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
    printf(    "\t\t    Filled cells       = %16ld  (%s)\n", sum_filled, mat[SeedID].getAlias().c_str());
    fprintf(fp,"\t\t    Filled cells       = %16ld  (%s)\n", sum_filled, mat[SeedID].getAlias().c_str());
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
    unsigned long upc = countCellM(d_mid, -1);
    
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
      
      replaced = fillByID(d_mid, FillID);
      
      if ( numProc > 1 )
      {
        if ( paraMngr->BndCommS3D(d_mid, ix, jx, kx, gd, gd) != CPM_SUCCESS ) Exit(0);
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
      printf(    "\t\tFilled cells           = %16ld  (%s)\n\n", sum_replaced, mat[FillID].getAlias().c_str());
      fprintf(fp,"\t\tFilled cells           = %16ld  (%s)\n\n", sum_replaced, mat[FillID].getAlias().c_str());
    }
    
    
    
    // ID=-1をカウントしてチェック
    upc = countCellM(d_mid, -1);
    
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

