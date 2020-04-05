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
 * @file   Geometry.C
 * @brief  FFV geometry related functions
 * @author aics
 */

#include "Geometry.h"
#include "FBUtility.h"
#include <limits.h>
#include <algorithm>
#include <numeric>

// quantizeCut のデバッグ出力
#define PROBE_DEBUG_UPDATE 0


// #################################################################
/**
 * @brief 各ノードのラベルのユニーク性を担保する
 * @param [out]    labelTop  各ランクのラベルリストの先頭ラベル
 * @param [out]    labelsz   ラベルの数
 * @param [out]    tbl       ローカルのラベルリスト
 * @param [in,out] mid       識別子配列
 * @param [in]     Dsize     配列サイズ
 */
void Geometry::assureUniqueLabel(int* labelTop, int* labelsz, vector<int>& tbl, int* mid, const int* Dsize)
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



  // 各ランクでラベルのリストを作成 逐次実行
#pragma omp single
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {

        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        const int dd = mid[m];

        if ( dd > 0 ) addLabel2List(tbl, dd);

      }
    }
  }


  /* DEBUG
  for (vector<int>::iterator it=tbl.begin(); it != tbl.end(); ++it )
  {
    fprintf(fpc, "%3d %3d vector\n", myRank, *it);
  }
  fprintf(sfpc, "\n");
  for (int i=0; i<tbl.size(); i++) {
    fprintf(fpc, "%3d %3d v[]\n", myRank, tbl[i]);
  }
   */


  // 各ランクのラベルの数
  int LabelSize = (int)tbl.size();
  //fprintf(fpc, "Rank = %6d : size = %d\n", myRank, LabelSize);



  // ラベル数の収集
  if ( numProc > 1 )
  {
    if ( paraMngr->Allgather(&LabelSize, 1, labelsz, 1, procGrp) != CPM_SUCCESS ) Exit(0);
  }
  else // serial
  {
    labelsz[0] = LabelSize;
  }

  /* DEBUG
  Hostonly_
  {
    printf("No of processes = %d\n", numProc);
    for (int i=0; i<numProc; i++)
    {
      printf("Rank = %6d : size = %d\n", i, labelsz[i]);
    }
    printf("\n");
  }
  */



  // 各ランクに割り振るラベルの開始番号を計算 >> labelTop[]に登録

  // 開始ランクを求める >> ラベルがないランクをスキップ
  int st = 0;
  for (int i=0; i<numProc; i++)
  {
    // ラベルがあれば1
    if (labelsz[i] != 0)
    {
      st = i;
      break;
    }
  }

  // Hostonly_ printf("start rank = %d\n\n", st);

  // エラーチェック
  if ( st > numProc-1 )
  {
    Hostonly_
    {
      fprintf(fpc, "Error : Start rank of gathering label is greather than %d\n", numProc);
    }
    Exit(0);
  }


  // ラベルの積算開始は1
  labelTop[st] = 1;

  // 続いて、各プロセスの開始数をセット
  for (int i=st+1; i<numProc; i++)
  {
    labelTop[i] = labelTop[i-1] + labelsz[i-1];
  }

  /* DEBUG
  Hostonly_
  {
    for (int i=0; i<numProc; i++)
    {
      fprintf(fpc, "original >> rank=%5d : num. labels=%6d : begin=%6d\n", i, labelsz[i], labelTop[i]);
    }
    printf("\n");
  }
  */


  // ラベルの変更
  // 先頭ラベルから個数分だけインクリメント

  // ラベルが存在する場合のみ実行
  if ( LabelSize > 0 )
  {

    // 各ランクの担当部分について、ラベルのふり直し
    int c = 0;
    int id= labelTop[myRank]; // 先頭のラベル

    for (vector<int>::iterator it=tbl.begin(); it != tbl.end(); ++it )
    {
      int target = *it;
      int replace= id + c;

      //fprintf(fpc, "replaced >> rank=%5d : %d :from=%6d  to=%6d\n", myRank, c, target, replace);
      c++;

      // Paint
#pragma omp parallel for firstprivate(ix, jx, kx, gd, target, replace) schedule(static)
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          for (int i=1; i<=ix; i++) {

            size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
            if ( mid[m] == target ) mid[m] = replace;
          }
        }
      }

    } // loop it

  }

  // ガイドセルの情報交換 >> makeConnectList()
  if ( numProc > 1 )
  {
    if ( paraMngr->BndCommS3D(mid, ix, jx, kx, gd, gd, procGrp) != CPM_SUCCESS ) Exit(0);
  }
}



#ifdef DISABLE_MPI
// #################################################################
/* @brief ポリゴングループの座標値からboxを計算する
 * @param [in]     PL         Polylibのインスタンス
 * @param [in,out] PG         PolygonPropertyクラス
 * @param [in]     NoPolyGrp  ポリゴングループ数
 */
void Geometry::calcBboxFromPolygonGroup(Polylib* PL,
                                        PolygonProperty* PG,
                                        const int m_NoPolyGrp)
{
  // 有次元空間でサーチ
  Vec3r m_min;
  Vec3r m_max;
  Vec3r t1(originD);
  Vec3r t2(pitchD);
  Vec3r t3;

  t3.assign((REAL_TYPE)size[0]*t2.x, (REAL_TYPE)size[1]*t2.y, (REAL_TYPE)size[2]*t2.z);

  // サブドメインの1層外側までをサーチ対象とする
  m_min = t1 - t2;
  m_max = t1 + t3 + t2;

#if 0 // debug
  printf("poly : %f %f %f\n", t1.x, t1.y, t1.z);
  printf("dx   : %f %f %f\n", t2.x, t2.y, t2.z);
  printf("Search area Bbox min : %f %f %f\n", m_min.x, m_min.y, m_min.z);
  printf("Search area Bbox max : %f %f %f\n", m_max.x, m_max.y, m_max.z);
#endif



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

      Vec3r *p;
      Vec3r bbox_min( 1.0e6,  1.0e6,  1.0e6);
      Vec3r bbox_max(-1.0e6, -1.0e6, -1.0e6);
      unsigned c=0;
      vector<Triangle*>::iterator it2;

      for (it2 = trias->begin(); it2 != trias->end(); it2++)
      {
        Vertex** org = (*it2)->get_vertex();
        Vec3r p[3];
        p[0] = *(org[0]);
        p[1] = *(org[1]);
        p[2] = *(org[2]);

        get_min(bbox_min, p[0]);
        get_min(bbox_min, p[1]);
        get_min(bbox_min, p[2]);

        get_max(bbox_max, p[0]);
        get_max(bbox_max, p[1]);
        get_max(bbox_max, p[2]);

#if 0
        /*
        printf("%d : p0=(%10.3e %10.3e %10.3e)  p1=(%10.3e %10.3e %10.3e) p2=(%10.3e %10.3e %10.3e) n=(%10.3e %10.3e %10.3e)\n", c++,
               p[0].x, p[0].y, p[0].z,
               p[1].x, p[1].y, p[1].z,
               p[2].x, p[2].y, p[2].z,
               n.x, n.y, n.z);
        */
        printf("%d : p0=(%10.3e %10.3e %10.3e)  p1=(%10.3e %10.3e %10.3e) p2=(%10.3e %10.3e %10.3e)\n", c++,
               p[0].x, p[0].y, p[0].z,
               p[1].x, p[1].y, p[1].z,
               p[2].x, p[2].y, p[2].z);
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
      Vec3r dummy(0.0, 0.0, 0.0);
      PG[m].setBboxMin(dummy);
      PG[m].setBboxMax(dummy);
    }

    m++;
  }

  //printf("R[%d] : (%10.3e %10.3e %10.3e) - (%10.3e %10.3e %10.3e)\n",
  //       myRank, originD[0], originD[1], originD[2],
  //       originD[0] + region[0], originD[1] + region[1], originD[2] + region[2]);

  // 領域内に収まっているかどうかをチェック >> ポリゴンは少しでも触れれば対象となり、領域外にはみ出すことがある
  for (int i=0; i<m_NoPolyGrp; i++)
  {
    Vec3r b_min = PG[i].getBboxMin();
    Vec3r b_max = PG[i].getBboxMax();

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
      if ( f_min[q] < originD[q] ) f_min[q] = originD[q];
      REAL_TYPE tmp = originD[q] + regionD[q];
      if ( f_max[q] > tmp ) f_max[q] = tmp;
    }

    Vec3r rmin(f_min);
    Vec3r rmax(f_max);

    PG[i].setBboxMin(rmin);
    PG[i].setBboxMax(rmax);

#if 0
    printf("A[%d] %20s : (%10.3e %10.3e %10.3e) - (%10.3e %10.3e %10.3e)\n",
           myRank, PG[i].getGroup().c_str(), rmin.x, rmin.y, rmin.z,
           rmax.x, rmax.y, rmax.z);
#endif
  }

}


#else // DISABLE_MPI
// #################################################################
/* @brief ポリゴングループの座標値からboxを計算する
 * @param [in]     PL         MPIPolylibのインスタンス
 * @param [in,out] PG         PolygonPropertyクラス
 * @param [in]     NoPolyGrp  ポリゴングループ数
 */
void Geometry::calcBboxFromPolygonGroup(MPIPolylib* PL,
                                        PolygonProperty* PG,
                                        const int m_NoPolyGrp)
{
  // 有次元空間でサーチ
  Vec3r m_min;
  Vec3r m_max;
  Vec3r t1(originD);
  Vec3r t2(pitchD);
  Vec3r t3;

  t3.assign((REAL_TYPE)size[0]*t2.x, (REAL_TYPE)size[1]*t2.y, (REAL_TYPE)size[2]*t2.z);

  // サブドメインの1層外側までをサーチ対象とする
  m_min = t1 - t2;
  m_max = t1 + t3 + t2;

#if 0 // debug
  printf("poly : %f %f %f\n", t1.x, t1.y, t1.z);
  printf("dx   : %f %f %f\n", t2.x, t2.y, t2.z);
  printf("Search area Bbox min : %f %f %f\n", m_min.x, m_min.y, m_min.z);
  printf("Search area Bbox max : %f %f %f\n", m_max.x, m_max.y, m_max.z);
#endif



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

      Vec3r bbox_min( 1.0e6,  1.0e6,  1.0e6);
      Vec3r bbox_max(-1.0e6, -1.0e6, -1.0e6);
      vector<Triangle*>::iterator it2;

      for (it2 = trias->begin(); it2 != trias->end(); it2++)
      {
        Vertex** org = (*it2)->get_vertex();
        Vec3r p[3];
        p[0] = *(org[0]);
        p[1] = *(org[1]);
        p[2] = *(org[2]);

        get_min(bbox_min, p[0]);
        get_min(bbox_min, p[1]);
        get_min(bbox_min, p[2]);

        get_max(bbox_max, p[0]);
        get_max(bbox_max, p[1]);
        get_max(bbox_max, p[2]);

#if 0
        /*
        printf("p0=(%10.3e %10.3e %10.3e)  p1=(%10.3e %10.3e %10.3e) p2=(%10.3e %10.3e %10.3e) n=(%10.3e %10.3e %10.3e)\n",
               p[0].x, p[0].y, p[0].z,
               p[1].x, p[1].y, p[1].z,
               p[2].x, p[2].y, p[2].z,
               n.x, n.y, n.z);
        */
        printf("p0=(%10.3e %10.3e %10.3e)  p1=(%10.3e %10.3e %10.3e) p2=(%10.3e %10.3e %10.3e)\n",
               p[0].x, p[0].y, p[0].z,
               p[1].x, p[1].y, p[1].z,
               p[2].x, p[2].y, p[2].z);
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
      Vec3r dummy(0.0, 0.0, 0.0);
      PG[m].setBboxMin(dummy);
      PG[m].setBboxMax(dummy);
    }

    m++;
  }

  //printf("R[%d] : (%10.3e %10.3e %10.3e) - (%10.3e %10.3e %10.3e)\n",
  //       myRank, originD[0], originD[1], originD[2],
  //       originD[0] + region[0], originD[1] + region[1], originD[2] + region[2]);

  // 領域内に収まっているかどうかをチェック >> ポリゴンは少しでも触れれば対象となり、領域外にはみ出すことがある
  for (int i=0; i<m_NoPolyGrp; i++)
  {
    Vec3r b_min = PG[i].getBboxMin();
    Vec3r b_max = PG[i].getBboxMax();

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
      if ( f_min[q] < originD[q] ) f_min[q] = originD[q];
      REAL_TYPE tmp = originD[q] + regionD[q];
      if ( f_max[q] > tmp ) f_max[q] = tmp;
    }

    Vec3r rmin(f_min);
    Vec3r rmax(f_max);

    PG[i].setBboxMin(rmin);
    PG[i].setBboxMax(rmax);

#if 0
    printf("A[%d] %20s : (%10.3e %10.3e %10.3e) - (%10.3e %10.3e %10.3e)\n",
           myRank, PG[i].getGroup().c_str(), rmin.x, rmin.y, rmin.z,
           rmax.x, rmax.y, rmax.z);
#endif
  }

}
#endif // DISABLE_MPI





// #################################################################
/**
 * @brief セル数をカウント
 * @param [in] mid     work array
 * @param [in] m_id    検査するID
 * @apram [in] mode    "global" / "local"
 * @param [in] painted m_id以外のセル数をカウント(false), m_idのセル数をカウント(true)
 * @param [in] Dsize   サイズ
 * @note painted : m_id以外のセル数をカウント(false), m_idのセル数をカウント(true)
 */
unsigned long Geometry::countCellM(const int* mid,
                                   const int m_id,
                                   const string mode,
                                   const bool painted,
                                   const int* Dsize)
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


  if ( !strcasecmp("global", mode.c_str()) )
  {
    if ( numProc > 1 )
    {
      unsigned long c_tmp = c;
      if ( paraMngr->Allreduce(&c_tmp, &c, 1, MPI_SUM, procGrp) != CPM_SUCCESS ) Exit(0);
    }
  }

  return c;
}


// #################################################################
/*
 * @brief bid情報を元にフラッドフィル
 * @param [in,out] bcd         BCindex B
 * @param [in]     bid         交点ID（5ビット幅x6方向）
 * @param [in]     mode        フィルモード (SOLID | FLUID)
 * @param [in]     Dsize       サイズ
 * @note Symmetric fillにより反復回数を減少
 */
unsigned long Geometry::fillbyBid(int* bcd,
                                  const int* bid,
                                  const int mode,
                                  const int* Dsize)
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

  // 隣接サブドメインのランク番号
  int sdw = nID[X_minus];
  int sde = nID[X_plus];
  int sds = nID[Y_minus];
  int sdn = nID[Y_plus];
  int sdb = nID[Z_minus];
  int sdt = nID[Z_plus];

  // 対称面と周期境界の場合の対応
  int mode_x = FillSuppress[0]; // if 0, suppress connectivity evaluation
  int mode_y = FillSuppress[1];
  int mode_z = FillSuppress[2];

  unsigned long filled   = 0; ///< ペイントされた数

  int fill_mode = mode;


#pragma omp parallel for reduction(+:filled)

  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {

#include "fill_bid_naive.h"

      }
    }
  }


#pragma omp parallel for reduction(+:filled)

  for (int k=kx; k>=1; k--) {
    for (int j=jx; j>=1; j--) {
      for (int i=ix; i>=1; i--) {

#include "fill_bid_naive.h"

      }
    }
  }

  if ( numProc > 1 )
  {
    unsigned long tmp = filled;
    if ( paraMngr->Allreduce(&tmp, &filled, 1, MPI_SUM, procGrp) != CPM_SUCCESS ) Exit(0);
  }

  return filled;
}




// #################################################################
/* @brief カット情報に基づくフィル操作
 * @param [in,out]  d_bcd      BCindex ID
 * @param [in,out]  d_bid      交点ID情報
 * @param [in,out]  d_mid      work array
 * @param [in,out]  d_cut      交点距離
 */
bool Geometry::fillbyCut(int* d_bcd,
                         int* d_bid,
                         int* d_mid,
                         int* d_cutL,
                         int* d_cutU)
{

  /*
   * 先に流体領域をフィルし、残った部分を固体領域としてフィルする
   * d_bid[] <= cmp[]のオーダーインデクス
   * d_bcd[] <= mat[]のオーダーインデクス
   */

  int ix = size[0];
  int jx = size[1];
  int kx = size[2];


  // 全計算内部セル数
  unsigned long total_cell = ix * jx * kx;

  if ( numProc > 1 )
  {
    unsigned long tmp_fc = total_cell;
    if ( paraMngr->Allreduce(&tmp_fc, &total_cell, 1, MPI_SUM, procGrp) != CPM_SUCCESS ) Exit(0);
  }


  // フィル対象のセル数
  unsigned long target_count = countCellB(d_bcd, 0);


  Hostonly_
  {
    fprintf(fpc,"\tFill initialize -----\n\n");
    fprintf(fpc,"\t\tTotal cell count       = %16ld\n", total_cell);
    fprintf(fpc,"\t\tInitial target count   = %16ld\n", target_count);

    fprintf(fpc,"\n\tFill -----\n\n");
  }


  // FLUIDでフィル ---------------------
  if ( !fillConnected(d_bcd, d_bid, FLUID, target_count) ) return false;

  // 対象セルがなければ終了
  if ( target_count == 0 ) return true;




  // 連結領域を同定し、ペイントする --------------------
  if ( !identifyConnectedRegion(d_mid, d_bcd, d_bid, target_count, total_cell-target_count) ) return false;



  /* 同種のセルでカットをもつ境界の修正
  unsigned modc = mergeSameSolidCell(d_bid, d_cutL, d_cutU, d_bcd);

  Hostonly_
  {
    fprintf(fpc,"\t\tModified connectivity  = %16d\n", modc);
  }
*/

  // チェック
  unsigned long upc = countCellB(d_bcd, 0);

  if ( upc != 0 )
  {
    Hostonly_
    {
      fprintf(fpc,"\n\tFill operation is done, but still remains %ld unpainted cells.\n\n", upc);
    }
    return false;
  }

  return true;
}





// #################################################################
/* @brief 連続領域のフィル
 * @param [in,out]  d_bcd        BCindex ID
 * @param [in]      d_bid        交点ID情報
 * @param [in]      fill_mode    フィルモード (SOLID | FLUID)
 * @param [in,out]  target_count ペイント対象のセル数
 * @retval success => true
 */
bool Geometry::fillConnected(int* d_bcd,
                             const int* d_bid,
                             const int fill_mode,
                             unsigned long& target_count)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  unsigned long filled=0;       ///< フィルされた数


  // 媒質のヒントを与える
  for (int m=0; m<NoHint; m++)
  {
    // MediumListの格納番号
    int key = FBUtility::findIDfromLabel(mat, NoMedium, fill_table[m].medium);

    if ( mat[key].getState() == fill_mode )
    {
      filled = fillSeedBcdOuter(d_bcd, fill_table[m].dir, key, d_bid);

      if ( numProc > 1 )
      {
        unsigned long c_tmp = filled;
        if ( paraMngr->Allreduce(&c_tmp, &filled, 1, MPI_SUM, procGrp) != CPM_SUCCESS ) Exit(0);
      }

    }
  }


  if ( numProc > 1 )
  {
    if ( paraMngr->BndCommS3D(d_bcd, ix, jx, kx, gd, gd, procGrp) != CPM_SUCCESS ) Exit(0);
  }



  Hostonly_
  {
    fprintf(fpc,"\t\tPainted %s cells by hint       = %16ld\n", (fill_mode==FLUID)?"FLUID":"SOLID", filled);
  }


  // ペイントされたシードセル数をターゲットから差し引く
  target_count -= filled;


  Hostonly_
  {
    fprintf(fpc,"\t\tRemaining cells to paint          = %16ld\n\n", target_count);
  }



  // 隣接するセルと同じfill_mode属性で接続している場合に隣接IDでフィル

  int c=0;
  unsigned long sum_filled = 0;   ///< フィルされた数の合計

  while (target_count > 0) {

    // 隣接媒質でフィルする
    filled = fillbyBid(d_bcd, d_bid, fill_mode);

    if ( numProc > 1 )
    {
      if ( paraMngr->BndCommS3D(d_bcd, ix, jx, kx, gd, gd, procGrp) != CPM_SUCCESS ) Exit(0);
    }

    target_count -= filled;
    sum_filled   += filled;

    c++;

    if ( filled <= 0 ) break; // フィル対象がなくなったら終了
  }


  Hostonly_
  {
    fprintf(fpc,"\t\tConnected fill Iteration          = %5d\n", c);
    fprintf(fpc,"\t\t               Filled by %s    = %16ld\n", (fill_mode==FLUID)?"FLUID":"SOLID", sum_filled);
    fprintf(fpc,"\t\t               Remaining cells    = %16ld\n\n", target_count);
  }

  return true;
}



// #################################################################
/* @brief 未ペイントセルをpaintIDで連結フィルする
 * @param [in,out] mid      識別子配列
 * @param [in]     bid      交点ID
 * @param [in]     targetID ペイントするID
 * @param [in]     Dsize    サイズ
 * @note ペイントしたセル数
 * @note 逐次処理
 * @todo 効率的な並列アルゴリズムがあれば
 */
unsigned long Geometry::fillConnected4ID(int* mid,
                                         const int* bid,
                                         const int paintID,
                                         const int* Dsize)
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



  // シードセルからフィル
  int pid = paintID;
  unsigned long c = 0;

#pragma omp single
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

    // 未ペイント
    if ( mid[m_p] == 0 )
    {
      int qq = bid[m_p];

      // 各方向がpaintIDで、かつ交点がないとき => 連結
      if ( ( mid[m_w] == pid && 0 == getBit5(qq, 0) ) ||
           ( mid[m_e] == pid && 0 == getBit5(qq, 1) ) ||
           ( mid[m_s] == pid && 0 == getBit5(qq, 2) ) ||
           ( mid[m_n] == pid && 0 == getBit5(qq, 3) ) ||
           ( mid[m_b] == pid && 0 == getBit5(qq, 4) ) ||
           ( mid[m_t] == pid && 0 == getBit5(qq, 5) ) )
      {
        mid[m_p] = pid;
        c++;
      }
    }
  }}}

  return c;
}




// #################################################################
/**
 * @brief 流体セルをマイナス値でワーク配列にコピー
 * @param [in,out] mid         work array
 * @param [in]     bcd         BCindex B
 * @note bcd[]には流体セルのIDのみが記録されている，それ以外はゼロを仮定
 * @retval 流体セルの個数
 */
unsigned long Geometry::fillFluidRegion(int* mid, const int* bcd)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  // zero initialize
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
  for (int k=1-gd; k<=kx+gd; k++) {
    for (int j=1-gd; j<=jx+gd; j++) {
      for (int i=1-gd; i<=ix+gd; i++) {
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        mid[m] = 0;
      }
    }
  }


  unsigned long c=0;

  // 流体セルを-1でペイント
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static) reduction(+:c)
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {

        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);

        if (DECODE_CMP( bcd[m] ) > 0)
        {
          mid[m] = -1;
          c++;
        }

      }
    }
  }

  if ( numProc > 1 )
  {
    unsigned long c_tmp = c;
    if ( paraMngr->Allreduce(&c_tmp, &c, 1, MPI_SUM, procGrp) != CPM_SUCCESS ) Exit(0);

    if ( paraMngr->BndCommS3D(mid, ix, jx, kx, gd, gd, procGrp) != CPM_SUCCESS ) Exit(0);
  }


  return c;
}



// #################################################################
/*
 * @brief bcd[]の外層にシードIDをペイントする
 * @param [in,out] bcd    BCindex B
 * @param [in]     face   ヒント面
 * @param [in]     target ペイントするIDのエントリ
 * @param [in]     bid    境界ID
 * @param [in]     Dsize  サイズ
 * @note ヒントとして与えられた外部境界面に接するセルにおいて，確実に流体セルであるセルをフィルする
 *       もし，外部境界面以外に固体候補があれば、ぬれ面はフィルしない
 */
unsigned long Geometry::fillSeedBcdOuter(int* bcd, const int face, const int target, const int* bid, const int* Dsize)
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
#pragma omp parallel for reduction(+:c)
        for (int k=1; k<=kx; k++) {
          for (int j=1; j<=jx; j++) {
            size_t m = _F_IDX_S3D(1, j, k, ix, jx, kx, gd);
            int s = bid[m];

            // 対象セルの周囲に指定方向以外にカットがなく，未ペイントのセルの場合
            if ( ((s & mx_w) == 0) && (DECODE_CMP(bcd[m]) == 0) )
            {
              setMediumID(bcd[m], tg);
              c++;
            }
          }
        }
      }
      break;

    case X_plus:
      if ( nID[face] < 0 )
      {
#pragma omp parallel for reduction(+:c)
        for (int k=1; k<=kx; k++) {
          for (int j=1; j<=jx; j++) {
            size_t m = _F_IDX_S3D(ix, j, k, ix, jx, kx, gd);
            int s = bid[m];
            if ( ((s & mx_e) == 0) && (DECODE_CMP(bcd[m]) == 0) )
            {
              setMediumID(bcd[m], tg);
              c++;
            }
          }
        }
      }
      break;

    case Y_minus:
      if ( nID[face] < 0 )
      {
#pragma omp parallel for reduction(+:c)
        for (int k=1; k<=kx; k++) {
          for (int i=1; i<=ix; i++) {
            size_t m = _F_IDX_S3D(i, 1, k, ix, jx, kx, gd);
            int s = bid[m];

            if ( ((s & mx_s) == 0) && (DECODE_CMP(bcd[m]) == 0) )
            {
              setMediumID(bcd[m], tg);
              c++;
            }
          }
        }
      }
      break;

    case Y_plus:
      if ( nID[face] < 0 )
      {
#pragma omp parallel for reduction(+:c)
        for (int k=1; k<=kx; k++) {
          for (int i=1; i<=ix; i++) {
            size_t m = _F_IDX_S3D(i, jx, k, ix, jx, kx, gd);
            int s = bid[m];

            if ( ((s & mx_n) == 0) && (DECODE_CMP(bcd[m]) == 0) )
            {
              setMediumID(bcd[m], tg);
              c++;
            }
          }
        }
      }
      break;

    case Z_minus:
      if ( nID[face] < 0 )
      {
#pragma omp parallel for reduction(+:c)
        for (int j=1; j<=jx; j++) {
          for (int i=1; i<=ix; i++) {
            size_t m = _F_IDX_S3D(i, j, 1, ix, jx, kx, gd);
            int s = bid[m];

            if ( ((s & mx_b) == 0) && (DECODE_CMP(bcd[m]) == 0) )
            {
              setMediumID(bcd[m], tg);
              c++;
            }
          }
        }
      }
      break;

    case Z_plus:
      if ( nID[face] < 0 )
      {
#pragma omp parallel for reduction(+:c)
        for (int j=1; j<=jx; j++) {
          for (int i=1; i<=ix; i++) {
            size_t m = _F_IDX_S3D(i, j, kx, ix, jx, kx, gd);
            int s = bid[m];

            if ( ((s & mx_t) == 0) && (DECODE_CMP(bcd[m]) == 0) )
            {
              setMediumID(bcd[m], tg);
              c++;
            }
          }
        }
      }
      break;

  } // end of switch

  return c;
}



// #################################################################
/* @brief 未ペイントセルのシードセルをひとつ探す
 * @param [in,out] mid      識別子配列
 * @param [in]     clabel   key ID
 * @param [in]     Dsize    サイズ
 * @note 候補となる最初のセルを見つけるため、逐次実行
 */
bool Geometry::findKeyCell(int* mid,
                           const int label,
                           const int* Dsize)
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

  int id = label;
  int c = 0;
  int flag = 0;

#pragma omp single
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {

        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);

        if ( mid[m]==0 && flag==0 )
        {
          mid[m] = id;
          flag = 1;
          c++;
        }

      }
    }
  }

  // シードセルがない => 全てのmid[] != 0
  if ( c == 0 ) return false;

  return true;
}


// #################################################################
/**
 * @brief list[]内の最頻値IDを求める
 * @param [in] m_sz      配列のサイズ
 * @param [in] list      ID配列
 * @note 候補がない場合には、0が戻り値
 */
int Geometry::find_mode(const int m_sz, const int* list)
{
  int key[CMP_BIT_W]; ///< ID毎の頻度 @note ffv_Initialize() >> fill()でif ( NoCompo+1 > CMP_BIT_W )をチェック
  memset(key, 0, sizeof(int)*CMP_BIT_W);


  for (int l=0; l<m_sz; l++) key[ list[l] ]++;


  int mode = 0; // サーチの初期値，IDの大きい方から
  int z = 0;    // 最頻値のID

  for (int l=NoCompo; l>=1; l--)
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
/**
 * @brief 各ランクのラベル情報を集める
 * @param [out] buffer       ラベル情報バッファ
 * @param [in]  width_label  バッファの最大要素数
 * @param [in]  cnct         接続リスト
 */
bool Geometry::gatherLabels(int* buffer,
                            const int width_label,
                            const vector< vector<int> > cnct)
{
  // 作業用バッファ
  int* wk = NULL;
  if ( !(wk = new int[width_label]) )  Exit(0);
  for (int i=0; i<width_label; i++) wk[i] = 0;


  int c = 0;
  for (int i=0; i<(int)cnct.size(); i++) // ローカルランクがもつラベルの数
  {
    for (int j=0; j<(int)cnct[i].size(); j++) // 各ラベルの接続情報の数
    {
      wk[c] = cnct[i][j];
      c++;
    }
  }
  if ( c > width_label ) { fprintf(fpc, "rank= %d c= %d\n", myRank, c); return false; }


  // ラベル数の収集
  if ( numProc > 1 )
  {
    if ( paraMngr->Allgather(wk, width_label, buffer, width_label, procGrp) != CPM_SUCCESS ) Exit(0);
  }
  else // serial
  {
    for (int i=0; i< width_label; i++) buffer[i] = wk[i];
  }

  if ( wk )  { delete [] wk; wk=NULL; }

  return true;
}



// #################################################################
/**
 * @brief 各ランクの最頻値情報を集める
 * @param [in,out]  mode      最頻値リスト
 */
void Geometry::gatherModes(vector<unsigned long>& mode)
{
  // 作業用バッファ
  unsigned long* wk = NULL;
  if ( !(wk = new unsigned long[NoMedium+1]) )  Exit(0);
  for (int i=0; i<NoMedium+1; i++) wk[i] = 0;


  unsigned long* buffer = NULL;
  if ( !(buffer = new unsigned long[ numProc*(NoMedium+1) ]) )   Exit(0);
  for (int i=0; i<numProc*(NoMedium+1); i++) buffer[i] = 0;

  // copy
  for (int i=0; i<NoMedium+1; i++)
  {
    wk[i] = mode[i];
  }


  // ラベル数の収集
  if ( numProc > 1 )
  {
    if ( paraMngr->Allgather(wk, NoMedium+1, buffer, NoMedium+1, procGrp) != CPM_SUCCESS ) Exit(0);
  }
  else // serial
  {
    for (int i=0; i<NoMedium+1; i++) buffer[i] = wk[i];
  }


  // zero初期化
  for (int i=0; i<NoMedium+1; i++) mode[i] = 0;

  // 全体の最頻値を求める
  for (int j=0; j<NoMedium+1; j++)
  {
    unsigned long c= 0;
    for (int i=0; i<numProc; i++)
    {
      c += buffer[ i*(NoMedium+1) + j ];
    }
    mode[j] = c;
  }


  if ( wk )  { delete [] wk; wk=NULL; }
  if ( buffer )  { delete [] buffer; buffer=NULL; }

}



// #################################################################
/**
 * @brief 各ランクの接続ルール数の情報を集める
 * @param [out] buffer      ルール要素数バッファ
 * @param [in]  width_rule  最大ルール数
 * @param [in]  cnct        接続リスト
 */
bool Geometry::gatherRules(int* buffer,
                           const int width_rule,
                           const vector< vector<int> > cnct)
{
  // 作業用バッファ
  int* wk = NULL;
  if ( !(wk = new int[width_rule]) )  Exit(0);
  for (int i=0; i<width_rule; i++) wk[i] = 0;

  if ( (int)cnct.size() > width_rule )
  {
    fprintf(fpc, "rank= %d rule= %d < %ld\n", myRank, width_rule, cnct.size());
    return false;
  }

  for (int i=0; i<(int)cnct.size(); i++) // ローカルランクがもつルールの数
  {
    wk[i] = cnct[i].size();
  }


  if ( numProc > 1 )
  {
    if ( paraMngr->Allgather(wk, width_rule, buffer, width_rule, procGrp) != CPM_SUCCESS ) Exit(0);
  }
  else // serial
  {
    for (int i=0; i< width_rule; i++) buffer[i] = wk[i];
  }

  if ( wk )  { delete [] wk; wk=NULL; }

  return true;
}




// #################################################################
/*
 * @brief フィルパラメータを取得
 * @param [in] tpCntl       TextParser
 * @param [in] Unit         DIMENSIONAL | NON_DIMENSIONAL
 * @param [in] m_RefL       代表長さ
 * @param [in] m_Nomedium   媒質数
 * @param [in] m_mat        MediumList
 * @param [in] m_fp         condition.txt
 */
void Geometry::getFillParam(TextParser* tpCntl,
                            const int Unit,
                            const REAL_TYPE m_RefL,
                            const int m_NoMedium,
                            const MediumList* m_mat,
                            FILE* m_fp)
{
  string str;
  string label_base, label_leaf, label;

  // ファイルポインタのコピー
  fpc = m_fp;

  NoMedium = m_NoMedium;

  mat = m_mat;

  label_base = "/FillHint";
  if ( (NoHint = tpCntl->countLabels(label_base)) < 1 )
  {
    Hostonly_ printf("\tParsing error in '%s' : No labels\n", label_base.c_str());
    Exit(0);
  }

  // フィルテーブルの作成
  fill_table = new KindFill[NoHint];


  for (int m=0; m<NoHint; m++)
  {
    if ( !(tpCntl->getNodeStr(label_base, m+1, str)) )
    {
      stamped_printf("\tParsing error : No Leaf Node \n");
      Exit(0);
    }

    fill_table[m].identifier = str;


    label_leaf = label_base + "/" + str;

    label = label_leaf + "/direction";
    if ( !(tpCntl->getInspectedValue(label, str )) )
    {
      Hostonly_ stamped_printf("\tParsing error : No '%s'\n", label.c_str());
      Exit(0);
    }
    if     ( !strcasecmp(str.c_str(), "xminus" ) ) fill_table[m].dir = X_minus;
    else if( !strcasecmp(str.c_str(), "xplus"  ) ) fill_table[m].dir = X_plus;
    else if( !strcasecmp(str.c_str(), "yminus" ) ) fill_table[m].dir = Y_minus;
    else if( !strcasecmp(str.c_str(), "yplus"  ) ) fill_table[m].dir = Y_plus;
    else if( !strcasecmp(str.c_str(), "zminus" ) ) fill_table[m].dir = Z_minus;
    else if( !strcasecmp(str.c_str(), "zplus"  ) ) fill_table[m].dir = Z_plus;
    else
    {
      fill_table[m].dir = X_minus;
      Hostonly_ printf("\tDefault 'X_minus' is set for Hint of FillSeed direction\n");
    }


    label = label_leaf + "/Medium";
    if ( !(tpCntl->getInspectedValue(label, str)) )
    {
      Hostonly_ stamped_printf("\tParsing error : No '%s'\n", label.c_str());
      Exit(0);
    }
    fill_table[m].medium = str;

  }


  Hostonly_
  {
    fprintf(fpc,"\n----------\n");
    fprintf(fpc,"\n\t>> Fill Hint\n\n");
    fprintf(fpc,"\t No.        Kind       Medium   Direction / Coordinate(ND)\n");
    fprintf(fpc,"\t--------------------------------------------------------------------------\n");
    for (int m=0; m<NoHint; m++)
    {
      fprintf(fpc,"\t%3d %12s   ", m+1, fill_table[m].medium.c_str());
      fprintf(fpc,"Direction = %s\n", FBUtility::getDirection(fill_table[m].dir).c_str());
    }
    fprintf(fpc,"\n----------\n");
  }


  // FillMediumがMediumList中にあるかどうかをチェック
  for (int m=0; m<NoHint; m++)
  {

    if ( FBUtility::findIDfromLabel(mat, NoMedium, fill_table[m].medium) == 0 )
    {
      Hostonly_
      {
        printf("/FillHint/%s/Medium = \"%s\" is not listed in MediumTable.\n",
               fill_table[m].identifier.c_str(),
               fill_table[m].medium.c_str());
      }
      Exit(0);
    }
  }


}



// #################################################################
/**
 * @brief ラベル情報集約のため、バッファサイズを計算する
 * @param [in] cnct  ローカルの接続ルールテーブル
 */
int Geometry::getNumBuffer(const vector< vector<int> > cnct)
{
  int width_label = 0;

  // 各ランク内の接続ルール要素の総和
  for (int i=0; i<(int)cnct.size(); i++)
  {
    width_label += cnct[i].size();
  }


  //printf("r=%d : No rule = %d : sep = %d : width_label=%d\n", myRank, cnct.size(), num_sep, width_label);

  if ( numProc > 1 )
  {
    int tmp = width_label;
    if ( paraMngr->Allreduce(&tmp, &width_label, 1, MPI_MAX, procGrp) != CPM_SUCCESS ) Exit(0);
  }

  return width_label;
}


// #################################################################
/*
 * @brief 連結領域を同定する
 * @param [in] mid          work array
 * @param [in] bcd          BCindex B
 * @param [in] bid          交点ID
 * @param [in] paintable    固体領域でペイントすべきセル数
 * @param [in] filled_fluid 流体でペイント済みの数
 */
bool Geometry::identifyConnectedRegion(int* mid,
                                       int* bcd,
                                       const int* bid,
                                       const unsigned long paintable,
                                       const unsigned long filled_fluid)
{
  // mid[]をゼロで初期化し、流体領域を-1でペイント
  // ここまでの処理でbcd[]には流体セルのIDのみが記録されている，それ以外はゼロ
  if ( fillFluidRegion(mid, bcd) != filled_fluid )
  {
    Hostonly_
    {
      fprintf(fpc, "\tFilled cells of fluid is not consistent.\n");
      return false;
    }
  }


  // このセクションは、逐次モードで実行

  int label = 1;

  while ( true )
  {

    // 未ペイントセルのシードセルをひとつ探す、falseは all mid[] != 0 を意味する
    if ( !findKeyCell(mid, label) ) break;

    // labelでペイントできなくなる（fillConnected4IDがゼロを返す）とbreak
    while ( true )
    {
      // 未ペイントセルをlabelで連結フィルする
      if ( fillConnected4ID(mid, bid, label) == 0 ) break;
    }

    // この時点で、labelでペイントできるセルは全てペイントしたので、labelをインクリメント
    label++;


    // 全てペイントされているなら、終了
    if ( countCellM(mid, 0, "local") == 0 ) break;
  }


  // この時点で、mid[]の内点には全て非ゼロの値が入る
  // 流体部分は-1、それ以外がlabel
  // labelは各ランク内ではユニークであるが、ランク間では重複している


  // ラベルを整理する

  // 各ランクのラベルリストの先頭
  int* labelTop = NULL;
  if ( !(labelTop = new int[numProc]) )   Exit(0);
  for (int i=0; i<numProc; i++) labelTop[i] = 0;


  // 各ランクのラベルの数 => 接続ルールの数
  int* labelSize = NULL;
  if ( !(labelSize = new int[numProc]) )   Exit(0);
  for (int i=0; i<numProc; i++) labelSize[i] = 0;


  // 各ランクのラベル保持コンテナ
  vector<int> local_Label;


  // 各ランクのラベルをユニークに定める
  assureUniqueLabel(labelTop, labelSize, local_Label, mid);



  // 接続リスト
  vector< vector<int> > cnnctTable;
  cnnctTable.resize( labelSize[myRank] ); // 各ランクの接続リストの値を保持するラベルの数


  // ローカルプロセスの接続ルールを作成
  makeConnectList(cnnctTable, labelTop, local_Label, mid, bid);


  /* DEBUG
  for (int i=0; i<cnnctTable.size(); i++)
  {
    for (int j=0; j<cnnctTable[i].size(); j++)
    {
      printf("rank= %d : key = %d , value = %d\n",
             myRank, labelTop[myRank]+i, cnnctTable[i][j]);  << key & values
    }
  }
   */

  // バッファ幅を決める => 各ランクでもつルールの要素数の最大値
  int width_label = getNumBuffer(cnnctTable);


  // バッファ 要素（ラベル）を各ランク毎にパックして集める
  int* packedLabels = NULL;
  if ( !(packedLabels = new int[numProc*width_label]) )   Exit(0);
  for (int i=0; i<numProc*width_label; i++) packedLabels[i] = 0;



  // ルール数バッファの幅
  int width_rule = 0;
  for (int i=0; i<numProc; i++)
  {
    if (width_rule < labelSize[i]) width_rule = labelSize[i];
  }
  //Hostonly_ printf("max rule = %d\n", width_rule);


  // バッファ 要素数を各ランク毎にパックして集める
  int* packedRules = NULL;
  if ( !(packedRules = new int[numProc*width_rule]) )   Exit(0);
  for (int i=0; i<numProc*width_rule; i++) packedRules[i] = 0;


  // ラベル情報の集約 プロセス間で通信
  if ( !gatherLabels(packedLabels, width_label, cnnctTable) ) return false;

  /* DEBUG
  Hostonly_ {
    printf("Labels\n");
    for (int i=0; i<numProc; i++) {
      printf("Rank %d : ", i);
      for (int j=0; j<width_label; j++) {
        printf(" %d", packedLabels[i*width_label+j]);
      }
      printf("\n");
    }
  }
  paraMngr->Barrier();
   */


  // ルール数の集約 プロセス間で通信
  if ( !gatherRules(packedRules, width_rule, cnnctTable) ) return false;


  /* DEBUG
  Hostonly_ {
    printf("Num of Rules\n");
    for (int i=0; i<numProc; i++) {
      printf("Rank %d : ", i);
      for (int j=0; j<width_rule; j++) {
        printf(" %d", packedRules[i*width_rule+j]);
      }
      printf("\n");
    }
  }
  paraMngr->Barrier();
   */


  // 接続ルール数の総和
  int NumRules = 0;
  for (int i=0; i<numProc; i++) NumRules += labelSize[i];
  //printf("Num of rules = %d\n", NumRule);


  // 全ての接続ルール
  vector< vector<int> > connectRules;
  connectRules.resize(NumRules);


  // 全接続ルールを作成
  if ( !makeWholeRules(connectRules,
                       labelSize,
                       packedLabels,
                       width_label,
                       packedRules,
                       width_rule) ) return false;


  // 最終のラベルと接続ルール
  vector<int> finalLabels;
  vector< vector<int> > finalRules;


  // ラベルを縮約
  reduceLabels(finalLabels, finalRules, connectRules);


  // DEBUG
  Hostonly_
  {
    fprintf(fpc, "\n ================\n");
    fprintf(fpc, "\n Reduced Labels\n\n");
    fprintf(fpc, "  key : values\n");
    for (int i=0; i<(int)finalRules.size(); i++) // ルール数
    {
      fprintf(fpc, " %4d :", finalLabels[i]);

      for (int j=0; j<(int)finalRules[i].size(); j++) // 各ルールの持つ要素数
      {
        fprintf(fpc, " %4d", finalRules[i][j]);
      }
      fprintf(fpc, "\n");
    }
    fprintf(fpc, "\n ================\n\n");
  }


  // ユニークにラベリング
  paintConnectedLabel(mid, finalLabels, finalRules);


  // ラベルに対応するインデクス
  vector<int> matIndex(finalLabels.size(), 0);


  // ラベルとmat[]インデクスの対応を決める
  if ( !makeHistgramByModalCut(matIndex, finalLabels, mid, bid) )
  {
    Hostonly_
    {
      fprintf(fpc, "\tFailed to perform makeHistgramByModalCut()\n");
      return false;
    }
  }


  // 対応づけしたmat[]インデクスによりラベル部分をペイントする
  paintCellByUniqueLabel(bcd, mid, finalLabels, matIndex);



  // cleanup
  if ( labelTop )     { delete [] labelTop;     labelTop     =NULL; }
  if ( labelSize )    { delete [] labelSize;    labelSize    =NULL; }
  if ( packedLabels ) { delete [] packedLabels; packedLabels =NULL; }
  if ( packedRules )  { delete [] packedRules;  packedRules  =NULL; }

  return true;
}



// #################################################################
/**
 * @brief 各ラベルの接続リストを作成
 * @param [in,out] cnct     接続リスト
 * @param [in]     ltop     各ランクのラベルリストの先頭
 * @param [in]     localTbl ローカルのラベルリスト
 * @param [in]     mid      識別子配列
 * @param [in]     bid      交点ID
 * @param [in]     Dsize    配列サイズ
 */
void Geometry::makeConnectList(vector< vector<int> >& cnct,
                               const int* ltop,
                               const vector<int> localTbl,
                               const int* mid,
                               const int* bid,
                               const int* Dsize)
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

  // 担当ランクの開始ラベル
  int st = ltop[myRank];

  // 各ラベルについて、接続リストを作成
  // 各ランク内については、単連結領域毎にラベルが割り当てられているので、領域境界のみでリストを作成可能
  for (int l=0; l<(int)cnct.size(); l++)
  {
    // 対象ラベル
    int id = st + l;


    // X-
#pragma omp single
    for (int k=1; k<=kx; k++) {
      for (int j=1; j<=jx; j++) {

        size_t m_p = _F_IDX_S3D(1, j, k, ix, jx, kx, gd);
        size_t m_w = _F_IDX_S3D(0, j, k, ix, jx, kx, gd);

        const int dd  = mid[m_p];
        const int d_w = mid[m_w];
        const int qw  = getBit5(bid[m_p], 0);

        // 対象ラベル、テスト方向にカットがなく、SOLIDラベルで、自ラベル（対象ラベル）でない
        if ( dd == id && qw==0 && d_w > 0 && d_w != id ) addLabel2List(cnct[l], d_w);
      }
    }


    // X+
#pragma omp single
    for (int k=1; k<=kx; k++) {
      for (int j=1; j<=jx; j++) {

        size_t m_p = _F_IDX_S3D(ix  , j, k, ix, jx, kx, gd);
        size_t m_e = _F_IDX_S3D(ix+1, j, k, ix, jx, kx, gd);

        const int dd  = mid[m_p];
        const int d_e = mid[m_e];
        const int qe  = getBit5(bid[m_p], 1);

        if ( dd == id && qe==0 && d_e > 0 && d_e != id ) addLabel2List(cnct[l], d_e);
      }
    }


    // Y-
#pragma omp single
    for (int k=1; k<=kx; k++) {
      for (int i=1; i<=ix; i++) {

        size_t m_p = _F_IDX_S3D(i, 1, k, ix, jx, kx, gd);
        size_t m_s = _F_IDX_S3D(i, 0, k, ix, jx, kx, gd);

        const int dd  = mid[m_p];
        const int d_s = mid[m_s];
        const int qs  = getBit5(bid[m_p], 2);

        if ( dd == id && qs==0 && d_s > 0 && d_s != id ) addLabel2List(cnct[l], d_s);
      }
    }


    // Y+
#pragma omp single
    for (int k=1; k<=kx; k++) {
      for (int i=1; i<=ix; i++) {

        size_t m_p = _F_IDX_S3D(i, jx,   k, ix, jx, kx, gd);
        size_t m_n = _F_IDX_S3D(i, jx+1, k, ix, jx, kx, gd);

        const int dd  = mid[m_p];
        const int d_n = mid[m_n];
        const int qn  = getBit5(bid[m_p], 3);

        if ( dd == id && qn==0 && d_n > 0 && d_n != id ) addLabel2List(cnct[l], d_n);
      }
    }


    // Z-
#pragma omp single
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {

        size_t m_p = _F_IDX_S3D(i, j, 1, ix, jx, kx, gd);
        size_t m_b = _F_IDX_S3D(i, j, 0, ix, jx, kx, gd);

        const int dd  = mid[m_p];
        const int d_b = mid[m_b];
        const int qb  = getBit5(bid[m_p], 4);

        if ( dd == id && qb==0 && d_b > 0 && d_b != id ) addLabel2List(cnct[l], d_b);
      }
    }


    // Z+
#pragma omp single
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {

        size_t m_p = _F_IDX_S3D(i, j, kx,   ix, jx, kx, gd);
        size_t m_t = _F_IDX_S3D(i, j, kx+1, ix, jx, kx, gd);

        const int dd  = mid[m_p];
        const int d_t = mid[m_t];
        const int qt  = getBit5(bid[m_p], 5);

        if ( dd == id && qt==0 && d_t > 0 && d_t != id ) addLabel2List(cnct[l], d_t);
      }
    }

  } // loop cnct

}



// #################################################################
/* @brief ラベルとmat[]インデクスの対応を決める
 * @param [in]     matIndex  ラベルに対応するインデクス
 * @param [in]     labels    ラベルリスト
 * @param [in]     mid       ラベル配列
 * @param [in]     bid       境界ID
 * @param [in]     Dsize     配列サイズ
 * @retval true -> success
 */
bool Geometry::makeHistgramByModalCut(vector<int>& matIndex,
                                      const vector<int> labels,
                                      const int* mid,
                                      const int* bid,
                                      const int* Dsize)
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

  // 頻度
  vector<unsigned long> mode(NoMedium+1, 0);

  int m_Medium = NoMedium;



  for (int l=0; l<(int)labels.size(); l++)
  {
    int target = labels[l];

    // 頻度積算をクリア
    for (int i=0; i<NoMedium+1; i++) mode[i] = 0;

    int flag = 0;
    //Hostonly_ fprintf(fpc, "label[%d] = %d\n", l, target);

#pragma omp parallel for firstprivate(ix, jx, kx, gd, target, m_Medium) \
        schedule(static) reduction(+:flag)
    for (int k=1; k<=kx; k++) {
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {

          size_t m = _F_IDX_S3D(i  , j  , k  , ix, jx, kx, gd);

          int qq = bid[m];
          int qw = getBit5(qq, 0);
          int qe = getBit5(qq, 1);
          int qs = getBit5(qq, 2);
          int qn = getBit5(qq, 3);
          int qb = getBit5(qq, 4);
          int qt = getBit5(qq, 5);

          // 対象セルが対象ラベルの場合、かつ、交点をもつ場合
          if ( mid[m] == target  &&  qw+qe+qs+qn+qb+qt > 0 )
          {
            // 交点IDの最頻値
            int sd = FBUtility::find_mode_id(qw, qe, qs, qn, qb, qt, NoCompo);
            if ( sd == 0 ) flag = 1;

            // 交点IDの媒質番号を得る
            int key = cmp[sd].getMatodr();
            if ( key == 0 )
            {
              printf("(%3d %3d %3d) sd=%2d key=%2d : %2d %2d %2d %2d %2d %2d\n", i,j,k,sd, key, qw, qe, qs, qn, qb, qt);
              flag++;
            }

            if ( cmp[key].getState() == FLUID ) // mat[]代用
            {
              printf("Modal Cut is FLUID mat[%d] (%d, %d, %d)\n", key, i, j, k);
              flag++;
            }

            if ( key < 1  ||  key > m_Medium )
            {
              printf("Out of range for mat[%d] at (%d, %d, %d)\n", key, i, j, k);
              flag++;
            }

            mode[key]++;
          }
        }
      }
    } // OMP


    // error check
    if ( numProc > 1 )
    {
      int c_tmp = flag;
      if ( paraMngr->Allreduce(&c_tmp, &flag, 1, MPI_SUM, procGrp) != CPM_SUCCESS ) Exit(0);
    }

    if ( flag > 0 )
    {
      fprintf(fpc, "rank=%d flag=%d\n", myRank, flag);
      return false;
    }


    /* DEBUG
    Hostonly_ fprintf(fpc, "          label : mode << LOCAL\n");
    for (int i=1; i<=NoMedium; i++) fprintf(fpc,"[rank=%d]  medium=%d %5d\n", myRank, i, mode[i]);
    paraMngr->Barrier();
     */

    // 全プロセスで集約 >> modeに全プロセスでの最頻値が返る
    gatherModes(mode);

    /*
    for (int j=0; j<NoMedium+1; j++)
    {
      Hostonly_ fprintf(fpc, "%d : %ld\n", j, mode[j]);
    }
     */

    // 最大値のmat[]インデクス
    int loc = 0;
    unsigned long max_mode = 0;
    for (int j=0; j<NoMedium+1; j++)
    {
      if (max_mode < mode[j])
      {
        max_mode = mode[j];
        loc = j;
      }
    }

    if ( loc < 0  ||  loc > NoMedium )
    {
      fprintf(fpc, "Out of range of max key = %d\n", loc);
      return false;
    }

    // インデクスの保存
    matIndex[l] = loc;

    /* DEBUG
    Hostonly_
    {
      fprintf(fpc, "label : loc. of max : mode\n");
      fprintf(fpc, "%5d :        %4d : %ld\n", l, loc, mode[loc]);
    }
     */


  } // label


  return true;
}



// #################################################################
/**
 * @brief 全接続ルールを作成
 * @param [in,out] connectRules 全ランクの接続ルール
 * @param [in]     labelSize    各ランクの接続ルール数
 * @param [in]     pckdLabel    ラベル情報のリスト
 * @param [in]     width_label  ラベル情報のバッファ幅
 * @param [in]     pckdRules    接続ルール数のリスト
 * @param [in]     width_rule   接続ルール数のバッファ幅
 */
bool Geometry::makeWholeRules(vector< vector<int> >& connectRules,
                              const int* labelSize,
                              const int* pckdLabel,
                              const int width_label,
                              const int* pckdRules,
                              const int width_rule)
{
  // 接続ルールをコンテナに入れる
  int c = 0;
  for (int i=0; i<numProc; i++)
  {
    // パックされているラベルをmでとりだす
    int m = 0;

    for (int j=0; j<labelSize[i]; j++) // 各ランクのルール数
    {
      int sz = pckdRules[ i * width_rule + j ]; // 各ルールの持つ要素数

      for (int l=0; l<sz; l++)
      {
        int elm = pckdLabel[ i * width_label + m ];
        if ( elm == 0 )
        {
          fprintf(fpc, "Zero element appears in constructing the whole connection rules at rank=%d #rule=%d > %dth\n", i, j, l);
        }
        addLabel2List(connectRules[c], elm);
        m++;
      }

      c++;

    }
  }

  // チェック
  if ( c != (int)connectRules.size() )
  {
    Hostonly_ {
      fprintf(fpc, "Number of rules does not agree %d\n", c);
      return false;
    }
  }


  // 昇順にならべる
  for (int i=0; i<(int)connectRules.size(); i++)
  {
    sort(connectRules[i].begin(), connectRules[i].end());
  }

  // 重複ルールの除去
  for (int i=0; i<(int)connectRules.size(); i++)
  {
    int key = i+1;

    // copy
    vector<int> wk;

    for (int j=0; j<(int)connectRules[i].size(); j++)
    {
      int lbl = connectRules[i][j];

      // key>lblのときは重複しているので追加しない
      if (key < lbl) wk.push_back(lbl);
    }

    // 入れ替え
    connectRules[i].clear();

    for (int j=0; j<(int)wk.size(); j++)
    {
      connectRules[i].push_back( wk[j] );
    }
  }



  // DEBUG
  Hostonly_
  {
    fprintf(fpc, "\n ================\n");
    fprintf(fpc, "\n Connection Rules (ascending order)\n\n");
    fprintf(fpc, "  key  : values\n");
    for (int i=0; i<(int)connectRules.size(); i++) // ルール数
    {
      fprintf(fpc, " %4d :", i+1);

      for (int j=0; j<(int)connectRules[i].size(); j++) // 各ルールの持つ要素数
      {
        fprintf(fpc, " %4d", connectRules[i][j]);
      }
      fprintf(fpc, "\n");
    }
    fprintf(fpc, "\n ================\n\n");
  }

  return true;
}



// #################################################################
/**
 * @brief 同種のセルでカットをもつ境界の修正
 * @param [in,out]  bid      カット点のID情報
 * @param [in,out]  cut      距離情報
 * @param [in]      bcd      d_bcd
 * @param [in]      Dsize    サイズ
 */
unsigned Geometry::mergeSameSolidCell(int* bid,
                                      int* cutL,
                                      int* cutU,
                                      const int* bcd,
                                      const int* Dsize)
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

  unsigned fc = 0;

#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {

        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);

        int zp = DECODE_CMP( bcd[m] );

        // SOLIDに対してチェック
        if ( cmp[zp].getState() == SOLID ) // cmp[]はmat[]の代用
        {
          int qq = bid[m];
          //int qw = getBit5(qq, 0);
          //int qe = getBit5(qq, 1);
          //int qs = getBit5(qq, 2);
          //int qn = getBit5(qq, 3);
          //int qb = getBit5(qq, 4);
          //int qt = getBit5(qq, 5);

          int posL = cutL[m];
          int posU = cutU[m];

          size_t m_e = _F_IDX_S3D(i+1, j,   k,   ix, jx, kx, gd);
          size_t m_w = _F_IDX_S3D(i-1, j,   k,   ix, jx, kx, gd);
          size_t m_n = _F_IDX_S3D(i,   j+1, k,   ix, jx, kx, gd);
          size_t m_s = _F_IDX_S3D(i,   j-1, k,   ix, jx, kx, gd);
          size_t m_t = _F_IDX_S3D(i,   j,   k+1, ix, jx, kx, gd);
          size_t m_b = _F_IDX_S3D(i,   j,   k-1, ix, jx, kx, gd);


          int zw = DECODE_CMP( bcd[m_w] );
          int ze = DECODE_CMP( bcd[m_e] );
          int zs = DECODE_CMP( bcd[m_s] );
          int zn = DECODE_CMP( bcd[m_n] );
          int zb = DECODE_CMP( bcd[m_b] );
          int zt = DECODE_CMP( bcd[m_t] );


          // X-方向が同じ媒質でカットがある場合
          if ( ensCutL(posL, X_minus)  &&  zp == zw )
          {
            setUncutL9(posL, X_minus);
            setBit5(qq, 0, X_minus);
            fc++;
            //fprintf(fpc, "rank %d : (%3d %3d %3d) dir=%d id=%2d : %3d : %d >> %3d : %d\n",
            //       myRank,i,j,k,0,zp,
            //       getBitL9(tmp, 0), ensCut(tmp, 0),
            //       getBitL9(pos, 0), ensCut(pos, 0) );
          }

          if ( ensCutL(posL, X_plus)  &&  zp == ze )
          {
            setUncutL9(posL, X_plus);
            setBit5(qq, 0, X_plus);
            fc++;
            //fprintf(fpc, "rank %d : (%3d %3d %3d) dir=%d id=%2d : %3d : %d >> %3d : %d\n",
            //       myRank,i,j,k,1,zp,
            //       getBit9(tmp, 1), ensCut(tmp, 1),
            //       getBit9(pos, 1), ensCut(pos, 1) );
          }

          if ( ensCutL(posL, Y_minus)  &&  zp == zs )
          {
            setUncutL9(posL, Y_minus);
            setBit5(qq, 0, Y_minus);
            fc++;
            //fprintf(fpc, "rank %d : (%3d %3d %3d) dir=%d id=%2d : %3d : %d >> %3d : %d\n",
            //       myRank,i,j,k,2,zp,
            //       getBit9(tmp, 2), ensCut(tmp, 2),
            //       getBit9(pos, 2), ensCut(pos, 2) );
          }

          if ( ensCutU(posU, Y_plus)  &&  zp == zn )
          {
            setUncutU9(posU, Y_plus);
            setBit5(qq, 0, Y_plus);
            fc++;
            //fprintf(fpc, "rank %d : (%3d %3d %3d) dir=%d id=%2d : %3d : %d >> %3d : %d\n",
            //       myRank,i,j,k,3,zp,
            //       getBit9(tmp, 3), ensCut(tmp, 3),
            //       getBit9(pos, 3), ensCut(pos, 3) );
          }

          if ( ensCutU(posU, Z_minus)  &&  zp == zb )
          {
            setUncutU9(posU, Z_minus);
            setBit5(qq, 0, Z_minus);
            fc++;
            //fprintf(fpc, "rank %d : (%3d %3d %3d) dir=%d id=%2d : %3d : %d >> %3d : %d\n",
            //       myRank,i,j,k,4,zp,
            //       getBit9(tmp, 4), ensCut(tmp, 4),
            //       getBit9(pos, 4), ensCut(pos, 4) );
          }

          if ( ensCutU(posU, Z_plus)  &&  zp == zt )
          {
            setUncutU9(posU, Z_plus);
            setBit5(qq, 0, Z_plus);
            fc++;
            //fprintf(fpc, "rank %d : (%3d %3d %3d) dir=%d id=%2d : %3d : %d >> %3d : %d\n",
            //       myRank,i,j,k,5,zp,
            //       getBit9(tmp, 5), ensCut(tmp, 5),
            //       getBit9(pos, 5), ensCut(pos, 5) );
          }

          cutL[m] = posL;
          cutU[m] = posU;
          bid[m] = qq;

        } // SOLID

      }
    }
  }

  if ( numProc > 1 )
  {
    unsigned c_tmp = fc;
    if ( paraMngr->Allreduce(&c_tmp, &fc, 1, MPI_SUM, procGrp) != CPM_SUCCESS ) Exit(0);
  }

  return fc;
}



// #################################################################
/* @brief 距離の最小値を求める
 * @param [in,out] cut カットの配列
 * @param [in]     bid 境界IDの配列
 */
void Geometry::minDistance(const int* cutL, const int* cutU, const int* bid)
{
  int global_min = 1024;
  int local_min = 1024;

  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

#pragma omp parallel firstprivate(ix, jx, kx, gd)
  {
    int th_min = 1024;

#pragma omp for schedule(static)
    for (int k=1; k<=kx; k++) {
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {

          size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
          int bd = bid[m];

          if ( IS_CUT(bd) ) // カットがあるか，IDによる判定
          {
            for (int n=0; n<3; n++)
            {
              int tmp = getBitL9(cutL[m], n);
              if ( (th_min > tmp) && tmp != 0) th_min = tmp;
            }
            for (int n=3; n<6; n++)
            {
              int tmp = getBitU9(cutU[m], n);
              if ( (th_min > tmp) && tmp != 0) th_min = tmp;
            }
          }

        }
      }
    }

#pragma omp critical
    {
      local_min = (std::min)(local_min, th_min);
    }
  }

  global_min = local_min;


  if ( numProc > 1 )
  {
    int tmp = global_min;
    if ( paraMngr->Allreduce(&tmp, &global_min, 1, MPI_MIN, procGrp) != CPM_SUCCESS ) Exit(0);
  }

  Hostonly_
  {
    fprintf(fpc, "\n\tMinimum non-dimensional distance except on a center = %e\n\n", (REAL_TYPE)global_min * DELTA_QT);
  }

}



// #################################################################
/* @brief 対応づけしたmat[]インデクスによりラベル部分をペイントする
 * @param [in,out] bcd      BCindex B
 * @param [in]     mid      識別子配列
 * @param [in]     labels   ラベルリスト
 * @param [in]     matIndex mat[]インデクス
 * @param [in]     Dsize    サイズ
 */
void Geometry::paintCellByUniqueLabel(int* bcd,
                                      const int* mid,
                                      const vector<int> labels,
                                      const vector<int> matIndex,
                                      const int* Dsize)
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


  for (int l=0; l<(int)labels.size(); l++)
  {
    int ref = labels[l];
    int pid = matIndex[l];


#pragma omp parallel for schedule(static) firstprivate(ix, jx, kx, gd, ref, pid)
    for (int k=1; k<=kx; k++) {
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {

          size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);

          if ( mid[m] == ref ) setMediumID(bcd[m], pid);
        }
      }
    }

  }


  // SYNC
  if ( numProc > 1 )
  {
    if ( paraMngr->BndCommS3D(bcd, ix, jx, kx, gd, gd, procGrp) != CPM_SUCCESS ) Exit(0);
  }
}



// #################################################################
/* @brief 接続している領域を指定ラベルでペイント
 * @param [in,out] mid      識別子配列
 * @param [in]     labels   ラベルリスト
 * @param [in]     rules    ルールリスト
 * @param [in]     Dsize    サイズ
 */
void Geometry::paintConnectedLabel(int* mid,
                                   const vector<int> labels,
                                   const vector< vector<int> > rules,
                                   const int* Dsize)
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


  for (int l=0; l<(int)rules.size(); l++)
  {
    int replace = labels[l];

    for (int m=0; m<(int)rules[l].size(); m++)
    {
      int target = rules[l][m];

      //Hostonly_ printf("target = %d : to be replaced %d\n", target, replace);

      unsigned long c = 0;

#pragma omp parallel for schedule(static) \
        firstprivate(ix, jx, kx, gd, target, replace) reduction(+:c)
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          for (int i=1; i<=ix; i++) {

            size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);

            if ( mid[m] == target )
            {
              mid[m] = replace;
              c++;
            }
          }
        }
      }

      //printf("replaced[%3d] %12ld rank[%d]\n", target, c, myRank);
    }
  }


  // SYNC
  if ( numProc > 1 )
  {
    if ( paraMngr->BndCommS3D(mid, ix, jx, kx, gd, gd, procGrp) != CPM_SUCCESS ) Exit(0);
  }
}



// #################################################################
/*
 * @brief 6方向にカットのあるセルの状態をsolidへ変更
 * @param [in,out]     bid       交点ID
 * @param [in,out]     cutL      距離情報 dir=0-2
 * @param [in,out]     cutU      距離情報 dir=3-5
 */
void Geometry::paintIsolatedCell(int* bid, int*cutL, int* cutU)
{
#define PROBE_DEBUG_ISOLATED 0
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

#if 0
  printf("DEBUG ========================\n");
  mark();
  Vec3r org(originD);
  Vec3r pch(pitchD);
  vector<Vec3r> xyz;
#endif
  
  unsigned replaced=0;

#pragma omp parallel for reduction(+:replaced)
  for (int k=1; k<=kx; k++) {
  for (int j=1; j<=jx; j++) {
  for (int i=1; i<=ix; i++) {
    size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);

    int qq = bid[m];
    int qw = getBit5(qq, 0);
    int qe = getBit5(qq, 1);
    int qs = getBit5(qq, 2);
    int qn = getBit5(qq, 3);
    int qb = getBit5(qq, 4);
    int qt = getBit5(qq, 5);

    // 全隣接方向に交点がある場合
    if ( qw * qe * qs * qn * qb * qt != 0 )
    {
      bid[m]  = onBit(qq, STATE_BIT);
      cutL[m] = onBit(cutL[m], STATE_BIT);
      cutU[m] = onBit(cutU[m], STATE_BIT);
      replaced++;
#if PROBE_DEBUG_ISOLATED
      Vec3r oi( (REAL_TYPE)i-0.5, (REAL_TYPE)j-0.5, (REAL_TYPE)k-0.5 );
      Vec3r o( org + pch * oi );
      xyz.push_back(o);
#endif
    }

  }}}

  if ( numProc > 1 )
  {
    unsigned c_tmp = replaced;
    if ( paraMngr->Allreduce(&c_tmp, &replaced, 1, MPI_SUM, procGrp) != CPM_SUCCESS ) Exit(0);
  }
  
  printf("\tNumber of Isolated cells by Cut >> replaced to SOLID = %d\n\n", replaced);
  
  
#if PROBE_DEBUG_ISOLATED
  std::ofstream ofs("all-cut-pnt.scab", std::ios::out | std::ios::binary);
  if (!ofs) {
    printf("\tCan't open 'all-cut-pnt.scab' file\n");
    return;
  }

  unsigned stp = 0;
  double tm = 0.0;
  unsigned csz = (unsigned)xyz.size();
  unsigned nvar = 1;
  
  printf("all cut points %d\n", csz);
  
  ofs.write((char*)&stp, sizeof(unsigned));
  ofs.write((char*)&tm,  sizeof(float));
  ofs.write((char*)&csz, sizeof(unsigned));
  ofs.write((char*)&nvar, sizeof(unsigned));
  float dummy = 1.0;
  
  for (unsigned i=0; i<csz; i++)
  {
    ofs.write((char*)&xyz[i].x, sizeof(float));
    ofs.write((char*)&xyz[i].y, sizeof(float));
    ofs.write((char*)&xyz[i].z, sizeof(float));
    ofs.write((char*)&dummy, sizeof(float));
  }
  
  ofs.close();
#endif
}



#ifdef DISABLE_MPI
// #################################################################
/*
 * @brief 交点計算を行い、量子化する
 * @param [out]    cut      量子化した交点
 * @param [out]    bid      交点ID
 * @param [in]     PL       Polylibインスタンス
 * @param [in]     PG       PolygonPropertyクラス
 * @param [in]     Dsize    サイズ
 */
void Geometry::quantizeCut(int* cutL,
                           int* cutU,
                           int* bid,
                           Polylib* PL,
                           PolygonProperty* PG,
                           const int* Dsize)
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


  Vec3r org(originD);
  Vec3r pch(pitchD);


  // ポリゴン情報へのアクセス
  vector<PolygonGroup*>* pg_roots = PL->get_root_groups();
  vector<PolygonGroup*>::iterator it;

  // ポリゴングループのループ
  int odr = 0;
  unsigned count = 0;
  for (it = pg_roots->begin(); it != pg_roots->end(); it++)
  {
    string m_pg = (*it)->get_name();     // グループラベル
    string m_bc = (*it)->get_type();     // 境界条件ラベル

    // サブドメイン内にポリゴンが存在する場合のみ処理する
    int ntria = (*it)->get_group_num_tria();  // ローカルのポリゴン数

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

        // 1層外側まで拡大
        for (int i=0; i<3; i++)
        {
          wmin[i] -= 1;
          wmax[i] += 1;
        }

#if 0
        printf("\t\t[rank=%6d] (%4d %4d %4d) - (%4d %4d %4d) %s \n",
               myRank, wmin[0], wmin[1], wmin[2], wmax[0], wmax[1], wmax[2], m_pg.c_str());
#endif

        // ポリゴングループの存在するbbox内のセルに対して交点計算
        for (int k=wmin[2]; k<=wmax[2]; k++) {
          for (int j=wmin[1]; j<=wmax[1]; j++) {
            for (int i=wmin[0]; i<=wmax[0]; i++) {

              // position of cell center
              Vec3r base((REAL_TYPE)i-0.5, (REAL_TYPE)j-0.5, (REAL_TYPE)k-0.5);

              // セルセンターを中心に2*pitchD[]の矩形領域
              Vec3r ctr(org + base * pch);

              // 1%マージンを与える
              Vec3r bx_min(ctr - pch * 1.01);
              Vec3r bx_max(ctr + pch * 1.01);


              vector<Triangle*>* trias = PL->search_polygons(m_pg, bx_min, bx_max, false); // false; ポリゴンが一部でもかかる場合
              int polys = trias->size();


              if (polys>0)
              {
                vector<Triangle*>::iterator it2;

                for (it2 = trias->begin(); it2 != trias->end(); it2++)
                {
                  Vertex** tmp = (*it2)->get_vertex();
                  Vec3r p[3];
                  p[0] = *(tmp[0]);
                  p[1] = *(tmp[1]);
                  p[2] = *(tmp[2]);

                  // Polygon ID
                  int poly_id = (*it2)->get_exid();

#if 0
                  printf("[%3d %3d %3d] (%f %f %f), (%f %f %f), (%f %f %f)\n", i,j,k,
                         p[0].x, p[0].y, p[0].z,
                         p[1].x, p[1].y, p[1].z,
                         p[2].x, p[2].y, p[2].z);
#endif

                  size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
                  int bb = bid[m];
                  int cL = cutL[m];
                  int cU = cutU[m];

                  // 各方向の交点を評価、短い距離を記録する。新規記録の場合のみカウント
                  count += updateCut(ctr, X_minus, p[0], p[1], p[2], cL, cU, bb, poly_id);
                  count += updateCut(ctr, X_plus,  p[0], p[1], p[2], cL, cU, bb, poly_id);
                  count += updateCut(ctr, Y_minus, p[0], p[1], p[2], cL, cU, bb, poly_id);
                  count += updateCut(ctr, Y_plus,  p[0], p[1], p[2], cL, cU, bb, poly_id);
                  count += updateCut(ctr, Z_minus, p[0], p[1], p[2], cL, cU, bb, poly_id);
                  count += updateCut(ctr, Z_plus,  p[0], p[1], p[2], cL, cU, bb, poly_id);

                  bid[m] = bb;
                  cutL[m] = cL;
                  cutU[m] = cU;
                } // trias

              } // polys

              //後始末
              delete trias;
            }
          }
        } // for i,j,k

      } // skip monitor
    } // ntria

    odr ++;
  }

  Hostonly_ fprintf(fpc, "\tquantize cut = %d\n\n", count);

  delete pg_roots;


  if ( numProc > 1 )
  {
    if ( paraMngr->BndCommS3D(bid, ix, jx, kx, gd, gd, procGrp) != CPM_SUCCESS ) Exit(0);
    if ( paraMngr->BndCommS3D(cut, ix, jx, kx, gd, gd, procGrp) != CPM_SUCCESS ) Exit(0);
  }

  
  // 全方向カットのあるセルをSOLIDへ変更
  paintIsolatedCell(bid, cutL, cutU);

  // 最小値カット
  minDistance(cut, bid);

}



#else // DISABLE_MPI
// #################################################################
/*
 * @brief 交点計算を行い、量子化する
 * @param [out]    cut      量子化した交点
 * @param [out]    bid      交点ID
 * @param [in]     PL       Polylibインスタンス
 * @param [in]     PG       PolygonPropertyクラス
 * @param [in]     Dsize    サイズ
 */
void Geometry::quantizeCut(int* cutL,
                           int* cutU,
                           int* bid,
                           MPIPolylib* PL,
                           PolygonProperty* PG,
                           const int* Dsize)
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


  Vec3r org(originD);
  Vec3r pch(pitchD);


  // ポリゴン情報へのアクセス
  vector<PolygonGroup*>* pg_roots = PL->get_root_groups();
  vector<PolygonGroup*>::iterator it;

  // ポリゴングループのループ
  int odr = 0;
  unsigned count = 0;
  for (it = pg_roots->begin(); it != pg_roots->end(); it++)
  {
    string m_pg = (*it)->get_name();     // グループラベル
    string m_bc = (*it)->get_type();     // 境界条件ラベル

    // サブドメイン内にポリゴンが存在する場合のみ処理する
    int ntria = (*it)->get_group_num_tria();  // ローカルのポリゴン数

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

        // 1層外側まで拡大
        for (int i=0; i<3; i++)
        {
          wmin[i] -= 1;
          wmax[i] += 1;
        }

        // ポリゴングループの存在するbbox内のセルに対して交点計算
        for (int k=wmin[2]; k<=wmax[2]; k++) {
        for (int j=wmin[1]; j<=wmax[1]; j++) {
        for (int i=wmin[0]; i<=wmax[0]; i++) {

          // position of cell center
          Vec3r base((REAL_TYPE)i-0.5, (REAL_TYPE)j-0.5, (REAL_TYPE)k-0.5);

          // セルセンターを中心に2*pitchD[]の矩形領域
          Vec3r ctr(org + base * pch);

          // 1%マージンを与える
          Vec3r bx_min(ctr - pch * 1.01);
          Vec3r bx_max(ctr + pch * 1.01);


          vector<Triangle*>* trias = PL->search_polygons(m_pg, bx_min, bx_max, false); // false; ポリゴンが一部でもかかる場合
          int polys = trias->size();


          if (polys>0)
          {
            vector<Triangle*>::iterator it2;

            for (it2 = trias->begin(); it2 != trias->end(); it2++)
            {
              Vertex** tmp = (*it2)->get_vertex();
              Vec3r p[3];
              p[0] = *(tmp[0]);
              p[1] = *(tmp[1]);
              p[2] = *(tmp[2]);

              // Polygon ID
              int poly_id = (*it2)->get_exid();

              size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
              int bb = bid[m];
              int cL = cutL[m];
              int cU = cutU[m];

              // 各方向の交点を評価、短い距離を記録する。新規記録の場合のみカウント
              count += updateCut(ctr, X_minus, p[0], p[1], p[2], cL, cU, bb, poly_id);
              count += updateCut(ctr, X_plus,  p[0], p[1], p[2], cL, cU, bb, poly_id);
              count += updateCut(ctr, Y_minus, p[0], p[1], p[2], cL, cU, bb, poly_id);
              count += updateCut(ctr, Y_plus,  p[0], p[1], p[2], cL, cU, bb, poly_id);
              count += updateCut(ctr, Z_minus, p[0], p[1], p[2], cL, cU, bb, poly_id);
              count += updateCut(ctr, Z_plus,  p[0], p[1], p[2], cL, cU, bb, poly_id);

              bid[m] = bb;
              cutL[m] = cL;
              cutU[m] = cU;
            } // trias

          } // polys

          //後始末
          delete trias;
        }}} // for i,j,k
      } // skip monitor
    } // ntria

    odr ++;
  }

  Hostonly_ {
    fprintf(fpc, "\tquantize cut = %d\n\n", count);
    printf("\tquantize cut = %d\n\n", count);
  }

  delete pg_roots;


  if ( numProc > 1 )
  {
    if ( paraMngr->BndCommS3D(bid, ix, jx, kx, gd, gd, procGrp) != CPM_SUCCESS ) Exit(0);
    if ( paraMngr->BndCommS3D(cutL, ix, jx, kx, gd, gd, procGrp) != CPM_SUCCESS ) Exit(0);
    if ( paraMngr->BndCommS3D(cutU, ix, jx, kx, gd, gd, procGrp) != CPM_SUCCESS ) Exit(0);
  }



  // 全方向カットのあるセルをSOLIDへ変更
  paintIsolatedCell(bid, cutL, cutU);

  // 最小値カット
  minDistance(cutL, cutU, bid);
  
#if PROBE_DEBUG_UPDATE // 法線チェック updateCut()のデバッグチェックも有効にする
  
  std::ofstream ofs("cut-pnt.scab", std::ios::out | std::ios::binary);
  if (!ofs) {
    printf("\tCan't open 'sdf-normal.scab' file\n");
    return;
  }
  
  unsigned stp = 0;
  double tm = 0.0;
  unsigned csz = (unsigned)prb_dbg.size();
  unsigned nvar = 3;
  
  printf("cut-pnt.scab : cut points %d was written\n", csz);
  
  ofs.write((char*)&stp, sizeof(unsigned));
  ofs.write((char*)&tm,  sizeof(float));
  ofs.write((char*)&csz, sizeof(unsigned));
  ofs.write((char*)&nvar, sizeof(unsigned));
  
  for (unsigned i=0; i<csz; i++)
  {
    ofs.write((char*)&prb_dbg[i].cp.x, sizeof(float));
    ofs.write((char*)&prb_dbg[i].cp.y, sizeof(float));
    ofs.write((char*)&prb_dbg[i].cp.z, sizeof(float));
    ofs.write((char*)&prb_dbg[i].nv.x, sizeof(float));
    ofs.write((char*)&prb_dbg[i].nv.y, sizeof(float));
    ofs.write((char*)&prb_dbg[i].nv.z, sizeof(float));
  }
  
  ofs.close();
#endif

}
#endif // DISABLE_MPI


// #################################################################
/**
 * @brief ラベルを縮約する
 * @param [out] finalLabels  縮約したラベルリスト
 * @param [out] finalRules   縮約したルールリスト
 * @param [in]  connectRules 全ランクの接続ルール
 */
void Geometry::reduceLabels(vector<int>& finalLabels,
                            vector< vector<int> >& finalRules,
                            vector< vector<int> > connectRules)
{
  // 接続ルールの総数
  int NumRules = connectRules.size();


  // Valid flag
  //      >0 -- テスト対象
  //     ==0 -- 非対象
  vector<int> valid(NumRules, 1);

  //for (int i=0; i<NumRules; i++) valid[i] = connectRules[i].size();

  int no_of_keys = 0;

  // マスターキーのループ
  for (int i=0; i<NumRules; i++)
  {
    const int key = i+1;

    if ( valid[i] > 0 )
    {
      no_of_keys++;

      // 最終リストにラベルkeyを追加
      finalLabels.push_back(key);

      // key[i]に属するラベルを作成
      registerRule(finalRules, no_of_keys, i, connectRules, valid);

      valid[i] = 0;
    } // valid

  } // NumRules


  // チェック
  for (int i=0; i<NumRules; i++)
  {
    if ( valid[i] !=0 )
    {
      Hostonly_ fprintf(fpc, "Error : Valid flag[%d] = %d\n", i, valid[i]);
      Exit(0);
    }
  }

  // 昇順ソート
  for (int i=0; i<(int)finalRules.size(); i++)
  {
    sort( finalRules[i].begin(), finalRules[i].end() );
  }

  // zero要素の処理
  for (int i=0; i<(int)finalRules.size(); i++)
  {
    if ( finalRules[i].size() == 0 )
    {
      finalRules[i].push_back( finalLabels[i] );
    }
  }


}



// #################################################################
/* @brief ラベルを登録しながら、ルールを縮約する
 * @param [in,out] finalRules   テストするルールリスト
 * @param [in]     reg_keys     ラベルを登録するマスターキーの数
 * @param [in]     test_keyodr  テストするconnectRulesのkey番号
 * @param [in]     connectRules 全ランクの接続ルール
 * @param [in,out] valid        有効フラグ
 */
void Geometry::registerRule(vector< vector<int> >& finalRules,
                            const int reg_keys,
                            const int test_keyodr,
                            vector< vector<int> >& connectRules,
                            vector<int>& valid)
{
  if ( valid[test_keyodr] == 0 ) return;

  const int key = test_keyodr + 1;

  // 登録ルールfinalRulesの配列を確保
  if ( (int)finalRules.size() != reg_keys ) { finalRules.resize(reg_keys); }

  // test_keyのルールの各ラベル
  for (int j=0; j<(int)connectRules[test_keyodr].size(); j++)
  {
    int lbl = connectRules[test_keyodr][j];

    // lbl>0 で key<lblのときのみ、マスターキーにラベルを登録
    if (lbl > 0 && key < lbl)
    {
      addLabel2List( finalRules[reg_keys-1], lbl );
    }

    // スキップモードに変更
    if (lbl > 0)
    {
      connectRules[test_keyodr][j] = -lbl;
    }

    // 再帰処理
    registerRule(finalRules, reg_keys, lbl-1, connectRules, valid);
  }

  valid[test_keyodr] = 0;
}




// #################################################################
/*
 * @brief ポリゴンと線分の交点計算
 * @param [in]   ray_o      レイの開始点
 * @param [in]   dir        レイの方向ベクトル（単位ベクトル）
 * @param [in]   v0         ポリゴンの頂点
 * @param [in]   v1         ポリゴンの頂点
 * @param [in]   v2         ポリゴンの頂点
 * @param [out]  pRetT      交点距離
 * @param [out]  pRetU      交点情報
 * @param [out]  pRetV      交点情報
 * @retval true -> 交点あり
 * @note Tomas Möller and Ben Trumbore.
         Journal of Graphics Tools, 2(1):21--28, 1997.
         http://www.graphics.cornell.edu/pubs/1997/MT97.html
 */
bool Geometry::TriangleIntersect(const Vec3r ray_o,
                                 const Vec3r dir,
                                 const Vec3r v0,
                                 const Vec3r v1,
                                 const Vec3r v2,
                                 REAL_TYPE& t,
                                 REAL_TYPE& u,
                                 REAL_TYPE& v)
{
  const REAL_TYPE m_eps = REAL_TYPE_EPSILON*2.0; //ROUND_EPS;
  Vec3r tvec, qvec;

  Vec3r e1 = v1 - v0;
  Vec3r e2 = v2 - v0;
  Vec3r pvec = cross(dir, e2);

  REAL_TYPE det = dot(e1, pvec);

  if (det > m_eps)
  {
    tvec = ray_o - v0;
    u = dot(tvec, pvec);
    if (u < 0.0 || u > det) return false;

    qvec = cross(tvec, e1);
    v = dot(dir, qvec);
    if (v < 0.0 || u + v > det) return false;
  }
  else if (det < -m_eps)
  {
    tvec = ray_o - v0;
    u = dot(tvec, pvec);
    if (u > 0.0 || u < det) return false;

    qvec = cross(tvec, e1);
    v = dot(dir, qvec);
    if (v > 0.0 || u + v < det) return false;

  }
  else
  {
    return false;
  }

  REAL_TYPE inv_det = 1.0 / det;

  t = dot(e2, qvec);

  t *= inv_det;
  u *= inv_det;
  v *= inv_det;

  return true;
}



/*
 * @brief 交点情報をアップデート
 * @param [in]     ray_o  レイの始点
 * @param [in]     dir    レイの方向
 * @param [in]     v0     テストする三角形の頂点
 * @param [in]     v1     テストする三角形の頂点
 * @param [in]     v2     テストする三角形の頂点
 * @param [in,out] cut    量子化交点距離情報
 * @param [in,out] bid    交点ID情報
 * @param [in]     pid    polygon id
 * @retval 新規交点の数
 * @note 短い距離を記録
 */
unsigned Geometry::updateCut(const Vec3r ray_o,
                             const int dir,
                             const Vec3r v0,
                             const Vec3r v1,
                             const Vec3r v2,
                             int& cutL,
                             int& cutU,
                             int& bid,
                             const int pid)
{
  // 単位方向ベクトルと格子幅
  Vec3r d;
  getDirVec(d, dir);
  REAL_TYPE pit= pitchD[0];

  // カウンタ
  unsigned count = 0;

  // 交点計算
  REAL_TYPE t, u, v;
  if ( !TriangleIntersect(ray_o, d, v0, v1, v2, t, u, v) ) return 0;

  // 格子幅で正規化
  REAL_TYPE tn = t / pit;

  if ( tn < 0.0 || 1.0 <= tn ) return 0;


  // 9bit幅の量子化
  int r = quantize9(tn);

  bool record = false;
  


  // 交点が記録されていない場合 >> 新規記録
  int e = (dir<3) ? ensCutL(cutL, dir) : ensCutU(cutU, dir);
  if ( e == 0 )
  {
    record = true;
    count = 1;
  }
  else // 交点が既に記録されている場合 >> 短い方を記録
  {
    int b = (dir<3) ? getBitL9(cutL, dir) : getBitU9(cutU, dir);
    if (r < b) record = true;
  }


  if ( record )
  {
    setBit5(bid, pid, dir);
    (dir<3) ? setCutL9(cutL, r, dir) : setCutU9(cutU, r, dir);
    
#if PROBE_DEBUG_UPDATE
    cut_probe p;
    p.cp = ray_o + d * (tn*pit);
    prb_dbg.push_back(p);
#endif
    /*
    printf("%10.6f %6d dir=%d : %10.6f %6d\n", tn, r, dir,
           (dir<3) ? getCutL9(cutL, dir) : getCutU9(cutU, dir),
           (dir<3) ? getBitL9(cutL, dir) : getBitU9(cutU, dir) ); */
  }

  return count;
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
          setMediumID(bcd[m_p], fid);
          c++;
        }
      }
    }
  }

  if ( numProc > 1 )
  {
    unsigned long c_tmp = c;
    if ( paraMngr->Allreduce(&c_tmp, &c, 1, MPI_SUM, procGrp) != CPM_SUCCESS ) Exit(0);
  }

  return c;
}






// #################################################################
/*
 * @brief ポリゴンからsdfに距離を転写し初期情報を作成
 * @param [in,out] sdf   SDF
 * @param [in]     PL    Polylibインスタンス
 * @param [out]    nrm   法線
 * @param [in]     Dsize サイズ
 * @note nrmのSTATEビットをSDFの初期化マスクとして使う
 */
void Geometry::polygon2sdf(REAL_TYPE* sdf,
                           MPIPolylib* PL,
                           int* nrm,
                           const int* Dsize)
{
#define SDF_DEBUG_NORMAL 0
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

  // SDFの初期値 >> FLT_MAX
  size_t size_n_cell = (ix+2*gd) * (jx+2*gd) * (kx+2*gd);

#pragma omp parallel for
  for (size_t i=0; i<size_n_cell; i++) {
    sdf[i] = FLT_MAX;
  }
  
  // 法線チェック
#if SDF_DEBUG_NORMAL
  vector<Vec3r> xyz;
  vector<Vec3r> nv;
#endif
  
  // 有次元空間でサーチ
  Vec3r bx_min;
  Vec3r bx_max;
  Vec3r org(originD);
  Vec3r pch(pitchD);
  Vec3r t3;
  t3.assign((REAL_TYPE)size[0]*pch.x, (REAL_TYPE)size[1]*pch.y, (REAL_TYPE)size[2]*pch.z);
  
  // サブドメインの2層外側（ガイドセル分）までをサーチ対象とする
  bx_min = org - pch - pch;
  bx_max = org + t3 + pch + pch;
  
  // 無次元距離の基準 格子幅 >> 計算座標なので格子幅が1.0となる
  REAL_TYPE dh = pitchD[0];
  
  // SDFの幅
  REAL_TYPE bnd = (REAL_TYPE)SEARCH_WIDTH;

  // ポリゴン情報へのアクセス
  vector<PolygonGroup*>* pg_roots = PL->get_root_groups();
  vector<PolygonGroup*>::iterator it;

  // ポリゴングループのループ
  for (it = pg_roots->begin(); it != pg_roots->end(); it++)
  {
    string m_pg = (*it)->get_name();     // グループラベル
    string m_bc = (*it)->get_type();     // 境界条件ラベル
    int ntria = (*it)->get_group_num_tria();  // ローカルのポリゴン数

    // ポリゴンが存在しなければ，スキップ
    if ( ntria == 0 ) continue;

    // Monitor属性のポリゴンはスキップ
    if ( !strcasecmp(m_bc.c_str(), "monitor")) continue;
    
    vector<Triangle*>* trias = PL->search_polygons(m_pg, bx_min, bx_max, false); // false; ポリゴンが一部でもかかる場合
            
    // ポリゴンがなければ，スキップ
    if ( trias->size() == 0 ) continue;
    
    
    vector<Triangle*>::iterator it2;
                
    // ポリゴングループに含まれるパッチへのアクセス
    for (it2 = trias->begin(); it2 != trias->end(); it2++)
    {
      Vertex** tmp = (*it2)->get_vertex();
      Vec3r p[3];
      p[0] = *(tmp[0]);
      p[1] = *(tmp[1]);
      p[2] = *(tmp[2]);
        
      // 対象パッチ一つのb-box
      Vec3r bbox_min(FLT_MAX);
      Vec3r bbox_max(-FLT_MAX);
      
      for (int i=0; i<3; i++) get_min(bbox_min, p[i]);
      for (int i=0; i<3; i++) get_max(bbox_max, p[i]);
      
      // SDFの初期値を与える格子点をポリゴンb-boxから片側5点分に設定
      Vec3r s_min(bbox_min.x-pch.x*bnd,
                  bbox_min.y-pch.y*bnd,
                  bbox_min.z-pch.z*bnd);
      Vec3r s_max(bbox_max.x+pch.x*bnd,
                  bbox_max.y+pch.y*bnd,
                  bbox_max.z+pch.z*bnd);
      
      // インデクス範囲（セルセンタ位置）
      Vec3r tt;
      tt = (s_min-org)/pch;
      Vec3i ids( (int)(tt.x+0.5), (int)(tt.y+0.5), (int)(tt.z+0.5) );
      
      tt = (s_max-org)/pch;
      Vec3i ide( (int)(tt.x+0.5)+1, (int)(tt.y+0.5)+1, (int)(tt.z+0.5)+1 );
      
      // 範囲チェック - 仮想セル2層
      if (ids.x < -1) ids.x = -1;
      if (ids.y < -1) ids.y = -1;
      if (ids.z < -1) ids.z = -1;
      if (ide.x > ix+2) ide.x = ix+2;
      if (ide.y > jx+2) ide.y = jx+2;
      if (ide.z > kx+2) ide.z = kx+2;
      
      Vec3r e1 = p[1] - p[0];
      Vec3r e2 = p[2] - p[0];
      Vec3r w = cross(e1, e2);
      
      // 対象範囲をテスト
      for (int k=ids.z; k<=ide.z; k++) {
      for (int j=ids.y; j<=ide.y; j++) {
      for (int i=ids.x; i<=ide.x; i++) {
        
        // d ベクトルを求める
        Vec3r oi( (REAL_TYPE)i-0.5, (REAL_TYPE)j-0.5, (REAL_TYPE)k-0.5 );
        Vec3r o( org + pch * oi );
        Vec3r s = p[0] - o;
        // 本来はcos(\theta)の符号で判断だが，内積の符号で代用 / ( s.length() * w.length() );
        Vec3r d = (dot(s, w) < 0) ? -w : w;
        d.normalize();
        
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        
        // 交点判定 qは無次元距離
        REAL_TYPE t, u, v, q;
        
        if ( TriangleIntersect(o, d, p[0], p[1], p[2], t, u, v) )
        {
          // 三角形内部
          q = t / dh;
          if (sdf[m] > q)
          {
            sdf[m] = q;
            // ポリゴンから流体側への法線
            setNrm9(nrm[m], -d);
            
#if SDF_DEBUG_NORMAL
            nv.push_back(-d);
            //nv.push_back( getNrm9(nrm[m]) );
            xyz.push_back(o);
#endif
          }
        }
        else
        {
          // 三角形外部
          int c;
          REAL_TYPE ds[3];
          ds[0] = distance(o, p[0]);
          ds[1] = distance(o, p[1]);
          ds[2] = distance(o, p[2]);
          c = (ds[1] < ds[2]) ? 1 : 2;
          if (ds[0] < ds[c]) c=0;
          t = ds[c];
          q = t / dh;
          if (sdf[m] > q)
          {
            sdf[m] = q;
            d = o-p[c];
            d.normalize();
            setNrm9(nrm[m], d);
            
#if SDF_DEBUG_NORMAL
            nv.push_back(d);
            xyz.push_back(o);
#endif
          }
        }
        
        // SDFを記録したセルのSTATE_BITを立てる
        nrm[m] = onBit(nrm[m] ,STATE_BIT);
      }}}

    } // it2

    //後始末
    delete trias;

  } // it

  delete pg_roots;

  
  if ( numProc > 1 )
  {
    // 最小値を通信するように考える
    if ( paraMngr->BndCommS3D(sdf, ix, jx, kx, gd, gd, procGrp) != CPM_SUCCESS ) Exit(0);
    
    // こちらは普通の同期
    if ( paraMngr->BndCommS3D(nrm, ix, jx, kx, gd, gd, procGrp) != CPM_SUCCESS ) Exit(0);
  }
  
  
  
#if SDF_DEBUG_NORMAL
  /*
    FILE* fp = NULL;
    if ( !(fp=fopen("normal.pwn", "w")) )
    {
      stamped_printf("\tSorry, can't open 'normal.pwn' file. Write failed.\n");
      exit(0);
    }
  
    fprintf(fp,"%d\n", xyz.size());
    for (vector<Vec3r>::iterator it=xyz.begin(); it != xyz.end(); ++it )
    {
      fprintf(fp,"%f %f %f\n", (*it).x, (*it).y, (*it).z);
    }
    for (vector<Vec3r>::iterator it=nv.begin(); it != nv.end(); ++it )
    {
      fprintf(fp,"%f %f %f\n", (*it).x, (*it).y, (*it).z);
    }

  if ( fp ) fclose(fp);
   */
  
  std::ofstream ofs("sdf-normal.scab", std::ios::out | std::ios::binary);
  if (!ofs) {
    printf("\tCan't open 'sdf-normal.scab' file at rank %d\n", myRank);
    return;
  }
  
  unsigned stp = 0;
  double tm = 0.0;
  unsigned csz = (unsigned)xyz.size();
  unsigned nvar = 3;
  
  //printf("cut points %d was written\n", csz);
  
  ofs.write((char*)&stp, sizeof(unsigned));
  ofs.write((char*)&tm,  sizeof(float));
  ofs.write((char*)&csz, sizeof(unsigned));
  ofs.write((char*)&nvar, sizeof(unsigned));
  
  for (unsigned i=0; i<csz; i++)
  {
    ofs.write((char*)&xyz[i].x, sizeof(float));
    ofs.write((char*)&xyz[i].y, sizeof(float));
    ofs.write((char*)&xyz[i].z, sizeof(float));
    ofs.write((char*)&nv[i].x, sizeof(float));
    ofs.write((char*)&nv[i].y, sizeof(float));
    ofs.write((char*)&nv[i].z, sizeof(float));
  }
  
  ofs.close();
#endif

}


// #################################################################
/*
 * @brief 外部境界条件の距離導入
 * @param [in,out] sdf   SDF
 * @param [out]    nrm   法線
 * @param [in]     oflag 外部境界条件フラグ
 */
void Geometry::outerSDF(REAL_TYPE* sdf,
                        int* nrm,
                        const int* oflag)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  // 計算空間座標なので，格子幅が1.0


  if ( oflag[X_minus] == 1 )
  {
    Vec3r d(1.0, 0.0, 0.0);
#pragma omp parallel for
    for (int i=1; i<=SEARCH_WIDTH; i++) {
      for (int k=1; k<=kx; k++) {
      for (int j=1; j<=jx; j++) {
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        REAL_TYPE q = (REAL_TYPE)i - 0.5;
        if (sdf[m] > q) {
          sdf[m] = q;
          setNrm9(nrm[m], d);
          nrm[m] = onBit(nrm[m] ,STATE_BIT);
        }
      }}
    }
  }
  
  
  if ( oflag[X_plus] == 1 )
  {
    Vec3r d(-1.0, 0.0, 0.0);
#pragma omp parallel for
    for (int i=ix; i>=ix-SEARCH_WIDTH+1; i--) {
      for (int k=1; k<=kx; k++) {
      for (int j=1; j<=jx; j++) {
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        REAL_TYPE q = (REAL_TYPE)(ix-i) + 0.5;
        if (sdf[m] > q) {
          sdf[m] = q;
          setNrm9(nrm[m], d);
          nrm[m] = onBit(nrm[m] ,STATE_BIT);
        }
      }}
    }
  }
  
  if ( oflag[Y_minus] == 1 )
  {
    Vec3r d(0.0, 1.0, 0.0);
#pragma omp parallel for
    for (int j=1; j<=SEARCH_WIDTH; j++) {
      for (int k=1; k<=kx; k++) {
      for (int i=1; i<=ix; i++) {
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        REAL_TYPE q = (REAL_TYPE)j - 0.5;
        if (sdf[m] > q) {
          sdf[m] = q;
          setNrm9(nrm[m], d);
          nrm[m] = onBit(nrm[m] ,STATE_BIT);
        }
      }}
    }
  }
  
  if ( oflag[Y_plus] == 1 )
  {
    Vec3r d(0.0, -1.0, 0.0);
#pragma omp parallel for
    for (int j=jx; j>=jx-SEARCH_WIDTH+1; j--) {
      for (int k=1; k<=kx; k++) {
      for (int i=1; i<=ix; i++) {
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        REAL_TYPE q = (REAL_TYPE)(jx-j) + 0.5;
        if (sdf[m] > q) {
          sdf[m] = q;
          setNrm9(nrm[m], d);
          nrm[m] = onBit(nrm[m] ,STATE_BIT);
        }
      }}
    }
  }
  
  if ( oflag[Z_minus] == 1 )
  {
    Vec3r d(0.0, 0.0, 1.0);
#pragma omp parallel for
    for (int k=1; k<=SEARCH_WIDTH; k++) {
      for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        REAL_TYPE q = (REAL_TYPE)k - 0.5;
        if (sdf[m] > q) {
          sdf[m] = q;
          setNrm9(nrm[m], d);
          nrm[m] = onBit(nrm[m] ,STATE_BIT);
        }
      }}
    }
  }
  
  if ( oflag[Z_plus] == 1 )
  {
    Vec3r d(0.0, 0.0, -1.0);
#pragma omp parallel for
    for (int k=kx; k>=kx-SEARCH_WIDTH+1; k--) {
      for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        REAL_TYPE q = (REAL_TYPE)(kx-k) + 0.5;
        if (sdf[m] > q) {
          sdf[m] = q;
          setNrm9(nrm[m], d);
          nrm[m] = onBit(nrm[m] ,STATE_BIT);
        }
      }}
    }
  }
  
  
  if ( numProc > 1 )
  {
    // 最小値を通信するように考える
    if ( paraMngr->BndCommS3D(sdf, ix, jx, kx, gd, gd, procGrp) != CPM_SUCCESS ) Exit(0);
    
    // こちらは普通の同期
    if ( paraMngr->BndCommS3D(nrm, ix, jx, kx, gd, gd, procGrp) != CPM_SUCCESS ) Exit(0);
  }
  
}

// #################################################################
/*
 * @brief 指定幅でトリム
 * @param [in,out] sdf   SDF
 * @param [in,out] nrm   法線
 * @param [in]     delta 指定幅
 */
void Geometry::trimNarrowBand(REAL_TYPE* sdf, int* nrm, const REAL_TYPE delta)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
#pragma omp parallel for
  for (int k=1; k<=kx; k++) {
  for (int j=1; j<=jx; j++) {
  for (int i=1; i<=ix; i++) {
    size_t mp = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
    if (sdf[mp] > delta)
    {
      sdf[mp] = FLT_MAX;
      nrm[mp] = offBit(nrm[mp] ,STATE_BIT);
    }
  }}}
  
  if ( numProc > 1 )
  {
    // 最小値を通信するように考える
    if ( paraMngr->BndCommS3D(sdf, ix, jx, kx, gd, gd, procGrp) != CPM_SUCCESS ) Exit(0);
    
    // こちらは普通の同期
    if ( paraMngr->BndCommS3D(nrm, ix, jx, kx, gd, gd, procGrp) != CPM_SUCCESS ) Exit(0);
  }
}



// #################################################################
/*
 * @brief レイヤーmarchingによるSDFの作成
 * @param [in,out] fn   SDF
 * @param [in]     fo   ワーク
 * @param [in]     nrm  法線
 * @param [in]     wk   ワーク
 * @param [in]     fv   ワーク
 */
bool Geometry::marchingSDF(REAL_TYPE* fn, REAL_TYPE* fo, int* nrm, int* wk, REAL_TYPE* fv)
{
  int num_layer = 0;
  
  
  // 1. フロントをサーチし，フロントレイヤー番号を記録
  
  if ( -1 == (num_layer = generateLayer(wk, nrm)) ) return false;
  
  Hostonly_ printf("\tSDF : Number of layers = %d\n", num_layer);
  

  // 2. ポリゴンから直接法線を計算していない部分について，法線を計算
  
  // fo[]にスムージング前の初期値を入れる
  setLayerInit(fn, fo, nrm, wk);
  
  // 非マスク部分のレイヤー場を平滑化
  smoothingScalar(nrm, fo, "Layer smoothing");
  printf("\n");
  
  // 平滑化したレイヤー場から法線を計算（非マスク部分）
  estimateNV(nrm, fo);
  
  // 非マスク部の法線を平滑化
  smoothingNV(nrm, fv);
  printf("\n");
  
  // 3. フロントのレイヤー順に距離場を計算
  
  tracingSDF(fn, nrm, wk, num_layer);
  printf("\n");

  return true;
}



// #################################################################
/*
 * @brief レイヤーの生成
 * @param [in,out]   wk   ワーク
 * @param [in]       nrm  法線
 */
int Geometry::generateLayer(int* wk, const int* nrm)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  int lx = std::max( ix, std::max(jx, kx))*5;
  int q = 0;
  size_t size_n_cell = (ix+2*gd) * (jx+2*gd) * (kx+2*gd);
  
  // -1で初期化
#pragma omp parallel for
  for (size_t i=0; i<size_n_cell; i++) {
    wk[i] = -1;
  }
  
  
  // ポリゴンから法線を直接転写した点は0を設定
#pragma omp parallel for
  for (int k=1; k<=kx; k++) {
  for (int j=1; j<=jx; j++) {
  for (int i=1; i<=ix; i++) {
    size_t mp = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
    if ( TEST_BIT(nrm[mp], STATE_BIT) ) wk[mp] = 0;
  }}}
  
  
  
  // フロントを追跡し，レイヤーを同定
  for (int layer=1; layer<=lx; layer++)
  {
    int a = 0;
    
#pragma omp parallel for reduction(+:q) reduction(+:a)
    for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
    for (int i=1; i<=ix; i++) {
      size_t mp = _F_IDX_S3D(i  , j  , k  , ix, jx, kx, gd);
      size_t me = _F_IDX_S3D(i+1, j  , k  , ix, jx, kx, gd);
      size_t mw = _F_IDX_S3D(i-1, j  , k  , ix, jx, kx, gd);
      size_t mn = _F_IDX_S3D(i  , j+1, k  , ix, jx, kx, gd);
      size_t ms = _F_IDX_S3D(i  , j-1, k  , ix, jx, kx, gd);
      size_t mt = _F_IDX_S3D(i  , j  , k+1, ix, jx, kx, gd);
      size_t mb = _F_IDX_S3D(i  , j  , k-1, ix, jx, kx, gd);
      
      int c = 0;
      int tgt = layer-1;
      
      if ( wk[mp] == -1 )
      {
        if ( wk[me] == tgt ) {
          c++;
        }
        if ( wk[mw] == tgt ) {
          c++;
        }
        if ( wk[mn] == tgt ) {
          c++;
        }
        if ( wk[ms] == tgt ) {
          c++;
        }
        if ( wk[mt] == tgt ) {
          c++;
        }
        if ( wk[mb] == tgt ) {
          c++;
        }
        
        if (c>0) {
          wk[mp] = layer;
          q++;
        }
        else {
          a++;
        }
      }

    }}}
    
    if (a == ix*jx*kx) break;
  } // layer
  
  // 最大レイヤー数
  int num_layer = 0;
  
  for (int k=1; k<=kx; k++) {
  for (int j=1; j<=jx; j++) {
  for (int i=1; i<=ix; i++) {
    size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
    if (wk[m] > num_layer) num_layer = wk[m];
  }}}
  
  if (num_layer > lx) {
    Hostonly_ stamped_printf("Number of layer is greater than IterationMax\n");
    return -1;
  }
  
  // check
  for (int k=1; k<=kx; k++) {
  for (int j=1; j<=jx; j++) {
  for (int i=1; i<=ix; i++) {
    size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
    if (wk[m] == -1)
    {
      Hostonly_ stamped_printf("still unprocessed point (%d %d %d)\n", i,j,k);
      return -1;
    }
  }}}
  
  return num_layer;
}




// #################################################################
/*
 * @brief レイヤーの初期値
 * @param [in,out] fn   SDF
 * @param [in]     fo   ワーク
 * @param [in]     nrm  法線
 * @param [in]     wk   ワーク
 */
void Geometry::setLayerInit(REAL_TYPE* fn, REAL_TYPE* fo, const int* nrm, const int* wk)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  size_t size_n_cell = (ix+2*gd) * (jx+2*gd) * (kx+2*gd);
  
  // 初期化されていないセルをゼロにする
#pragma omp parallel for
  for (size_t i=0; i<size_n_cell; i++) {
    if ( !TEST_BIT(nrm[i], STATE_BIT) ) fn[i] = 0.0;
  }
  
  // floatにコピー
#pragma omp parallel for
  for (size_t i=0; i<size_n_cell; i++) {
    fo[i] = (REAL_TYPE)wk[i];
  }

}



// #################################################################
/*
 * @brief レイヤー場のスムージング
 * @param [in]     nrm  法線
 * @param [in,out] fo   ワーク
 */
void Geometry::smoothingScalar(int* nrm, REAL_TYPE* fo, const string str)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  int Itrmax = 500;
  const REAL_TYPE eps = 1.0e-6;
  REAL_TYPE tr = 0.0;
  
  // トレース
#pragma omp parallel for reduction(+:tr)
  for (int k=1; k<=kx; k++) {
  for (int j=1; j<=jx; j++) {
  for (int i=1; i<=ix; i++) {
    size_t mp = _F_IDX_S3D(i  , j  , k  , ix, jx, kx, gd);
    REAL_TYPE b = GET_SHIFT_F(nrm[mp], STATE_BIT);
    
    tr += fabs(fo[mp])*(1.0-b);
  }}}
  
  
  
  int l;
  REAL_TYPE res;
  
  for (l=1; l<=Itrmax; l++)
  {
    res = 0.0;

    // BC
#pragma omp parallel for
    for (int j=1; j<=jx; j++) {
    for (int i=1; i<=ix; i++) {
      size_t m1 = _F_IDX_S3D(i  , j  , 1   , ix, jx, kx, gd);
      size_t m2 = _F_IDX_S3D(i  , j  , kx  , ix, jx, kx, gd);
      size_t mt = _F_IDX_S3D(i  , j  , kx+1, ix, jx, kx, gd);
      size_t mb = _F_IDX_S3D(i  , j  , 0   , ix, jx, kx, gd);
      fo[mt] = fo[m2];
      fo[mb] = fo[m1];
    }}
         
#pragma omp parallel for
    for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      size_t m1 = _F_IDX_S3D(1   , j  , k  , ix, jx, kx, gd);
      size_t m2 = _F_IDX_S3D(ix  , j  , k  , ix, jx, kx, gd);
      size_t me = _F_IDX_S3D(ix+1, j  , k  , ix, jx, kx, gd);
      size_t mw = _F_IDX_S3D(0   , j  , k  , ix, jx, kx, gd);
      fo[me] = fo[m2];
      fo[mw] = fo[m1];
    }}
          
#pragma omp parallel for
    for (int k=1; k<=kx; k++) {
    for (int i=1; i<=ix; i++) {
      size_t m1 = _F_IDX_S3D(i  , 1   , k  , ix, jx, kx, gd);
      size_t m2 = _F_IDX_S3D(i  , jx  , k  , ix, jx, kx, gd);
      size_t mn = _F_IDX_S3D(i  , jx+1, k  , ix, jx, kx, gd);
      size_t ms = _F_IDX_S3D(i  , 0   , k  , ix, jx, kx, gd);
      fo[mn] = fo[m2];
      fo[ms] = fo[m1];
    }}
    
#pragma omp parallel for reduction(+:res)
    for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
    for (int i=1; i<=ix; i++) {
      size_t mp = _F_IDX_S3D(i  , j  , k  , ix, jx, kx, gd);
      size_t me = _F_IDX_S3D(i+1, j  , k  , ix, jx, kx, gd);
      size_t mw = _F_IDX_S3D(i-1, j  , k  , ix, jx, kx, gd);
      size_t mn = _F_IDX_S3D(i  , j+1, k  , ix, jx, kx, gd);
      size_t ms = _F_IDX_S3D(i  , j-1, k  , ix, jx, kx, gd);
      size_t mt = _F_IDX_S3D(i  , j  , k+1, ix, jx, kx, gd);
      size_t mb = _F_IDX_S3D(i  , j  , k-1, ix, jx, kx, gd);
      
      REAL_TYPE b = GET_SHIFT_F(nrm[mp], STATE_BIT);
      REAL_TYPE s = (fo[me]+fo[mw]+fo[mn]
                    +fo[ms]+fo[mt]+fo[mb]) / 6.0;
      REAL_TYPE ds = (s - fo[mp])*(1.0-b);
      fo[mp] += ds;
      res += ds*ds;
    }}}
    
    REAL_TYPE rr = sqrt(res)/tr;

    if ( rr < eps )
    {
      break;
    }
  } // iteration loop
  
  Hostonly_ printf("\t\t %s : Iteration=%4d : res= %e\n", str.c_str(), l, sqrt(res)/tr);
  
}



// #################################################################
/*
 * @brief 法線の推定
 * @param [in,out] nrm  法線
 * @param [in]     fo   レイヤー
 */
void Geometry::estimateNV(int* nrm, const REAL_TYPE* fo)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  const REAL_TYPE eps = 1.0e-5;
  
#pragma omp parallel for
  for (int k=1; k<=kx; k++) {
  for (int j=1; j<=jx; j++) {
  for (int i=1; i<=ix; i++) {
    size_t mp = _F_IDX_S3D(i  , j  , k  , ix, jx, kx, gd);
    size_t me = _F_IDX_S3D(i+1, j  , k  , ix, jx, kx, gd);
    size_t mw = _F_IDX_S3D(i-1, j  , k  , ix, jx, kx, gd);
    size_t mn = _F_IDX_S3D(i  , j+1, k  , ix, jx, kx, gd);
    size_t ms = _F_IDX_S3D(i  , j-1, k  , ix, jx, kx, gd);
    size_t mt = _F_IDX_S3D(i  , j  , k+1, ix, jx, kx, gd);
    size_t mb = _F_IDX_S3D(i  , j  , k-1, ix, jx, kx, gd);

    // 非マスク部分のみ
    if ( 0 == BIT_SHIFT(nrm[mp], STATE_BIT) )
    {
      REAL_TYPE px = fo[me] - fo[mw];
      REAL_TYPE py = fo[mn] - fo[ms];
      REAL_TYPE pz = fo[mt] - fo[mb];
      REAL_TYPE u  = px / (fabs(px) + eps);
      REAL_TYPE v  = py / (fabs(py) + eps);
      REAL_TYPE w  = pz / (fabs(pz) + eps);

      Vec3r nv(u, v, w);
      nv.normalize();
      setNrm9(nrm[mp], nv);
    }

  }}}

}




// #################################################################
/*
 * @brief 法線のスムージング
 * @param [in]     nrm  法線
 * @param [in,out] fv   ワーク
 */
void Geometry::smoothingNV(int* nrm, REAL_TYPE* fv)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  size_t sn = (ix+2*gd) * (jx+2*gd) * (kx+2*gd);
  
  // fv[]に法線をコピー
  getNVfromIdx(nrm, fv);
  
  Hostonly_ printf("\tSDF : Norml smoothing\n");
  
  // nv.x
  smoothingScalar(nrm, &fv[0], "NV.x smoothing ");
  
  // nv.y
  smoothingScalar(nrm, &fv[sn], "NV.y smoothing ");
  
  // nv.z
  smoothingScalar(nrm, &fv[sn*2], "NV.z smoothing ");
  
  
  // nrmへのストア
      
#pragma omp parallel for
  for (int k=1; k<=kx; k++) {
  for (int j=1; j<=jx; j++) {
  for (int i=1; i<=ix; i++) {
    size_t mp = _F_IDX_S3D(i, j, k,    ix, jx, kx, gd);
    size_t m0 = _F_IDX_V3D(i, j, k, 0, ix, jx, kx, gd);
    size_t m1 = _F_IDX_V3D(i, j, k, 1, ix, jx, kx, gd);
    size_t m2 = _F_IDX_V3D(i, j, k, 2, ix, jx, kx, gd);
      
    Vec3r nv(fv[m0], fv[m1], fv[m2]);
    nv.normalize();
    setNrm9(nrm[mp], nv);
  }}}
}



// #################################################################
/*
 * @brief 法線の推定
 * @param [in]     nrm  法線
 * @param [out]    f    ワーク
 */
void Geometry::getNVfromIdx(const int* nrm, REAL_TYPE* f)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  for (int l=0; l<3; l++) {
    
#pragma omp parallel for
    for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
    for (int i=1; i<=ix; i++) {
      size_t mp = _F_IDX_S3D(i, j, k,    ix, jx, kx, gd);
      size_t m0 = _F_IDX_V3D(i, j, k, 0, ix, jx, kx, gd);
      size_t m1 = _F_IDX_V3D(i, j, k, 1, ix, jx, kx, gd);
      size_t m2 = _F_IDX_V3D(i, j, k, 2, ix, jx, kx, gd);
      
      Vec3r nv = getNrm9(nrm[mp]);
      nv.normalize();
      f[m0] = nv.x;
      f[m1] = nv.y;
      f[m2] = nv.z;
    }}}
  }
  
  return;
}


// #################################################################
/*
 * @brief SDFのレイヤー方向へのマーチング
 * @param [in,out] fn     SDF
 * @param [in]     nrm    法線
 * @param [in]     wk     レイヤー配列
 * @param [in]     nLayer レイヤー数
 * @note 最大10回の積分を想定
 *       |u|~1なのでu.x~1, u.x*dh がほぼ1格子幅分の移動量
 *       数回積分する程度
 */
void Geometry::tracingSDF(REAL_TYPE* fn, int* nrm, const int* wk, const int nLayer)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  int lx = ix + jx + kx;

  // 計算空間座標で考える
  Hostonly_ printf("\tSDF : Tracing distance\n");
  
  for (int l=1; l<=nLayer; l++)
  {
    int c = 0;
    
#pragma omp parallel for reduction(+:c)
    for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
    for (int i=1; i<=ix; i++) {
      size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
            
      if (wk[m] == l && !TEST_BIT(nrm[m], STATE_BIT) )
      {
        Vec3r nv = -getNrm9( nrm[m] ).normalize(); // 法線の逆方向に
        Vec3r o( (REAL_TYPE)i, (REAL_TYPE)j, (REAL_TYPE)k );
        Vec3r p;         // 次の積分点ベクトル
        Vec3r v;         // 次の積分点での法線ベクトル
        REAL_TYPE d=0.0; // 次の積分点までの距離（積分距離は積算していく）
        REAL_TYPE s=0.0; // 次の積分点でのSDF値

        int q;
        for (q=1; q<=lx; q++)
        {
          if ( getDistance(o, nv, nrm, fn, d, s, p, v) )
          {
            fn[m] = d + s;
            nrm[m] = onBit(nrm[m], STATE_BIT);
            break;
          }
          else // 候補でなければ，積分を続ける
          {
            o = p;
            nv= v;
          }
        }
        if  ( q==lx ) c++;
      }
    }}}
    
    if ( c>0 ) printf("Error : tracking at rank=%d\n", myRank);
  } // nLayer
  
}



// #################################################################
/*
 * @brief 開始座標点oからn方向へ移動した点のSDF値とその距離の和
 * @param [in]  o   開始座標点
 * @param [in]  n   法線
 * @param [in]  nrm ビット配列
 * @param [in]  f   SDF
 * @param [out] ds  移動距離（積分距離は積算していく）
 * @param [out] s   積分点のSDF
 * @param [out] p   積分点の位置
 * @param [out] v   積分点の放線
 */
bool Geometry::getDistance(const Vec3r o,
                           const Vec3r n,
                           const int* nrm,
                           const REAL_TYPE* f,
                           REAL_TYPE& ds,
                           REAL_TYPE& s,
                           Vec3r& p,
                           Vec3r& v)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  REAL_TYPE dt = 0.5;
  
  Vec3r d(n.x*dt, n.y*dt, n.z*dt);
  p = o + d;
  ds += d.length(); // 積算する
  
  // 点pが属するセルセンタインデクスの基点
  int i = (int)(p.x);
  int j = (int)(p.y);
  int k = (int)(p.z);
  
  Vec3r g;
  g.x = o.x - (REAL_TYPE)i + d.x;
  g.y = o.y - (REAL_TYPE)j + d.y;
  g.z = o.z - (REAL_TYPE)k + d.z;
  
  size_t m0 = _F_IDX_S3D(i  , j  , k  , ix, jx, kx, gd);
  size_t m1 = _F_IDX_S3D(i+1, j  , k  , ix, jx, kx, gd);
  size_t m2 = _F_IDX_S3D(i  , j+1, k  , ix, jx, kx, gd);
  size_t m3 = _F_IDX_S3D(i+1, j+1, k  , ix, jx, kx, gd);
  size_t m4 = _F_IDX_S3D(i  , j  , k+1, ix, jx, kx, gd);
  size_t m5 = _F_IDX_S3D(i+1, j  , k+1, ix, jx, kx, gd);
  size_t m6 = _F_IDX_S3D(i  , j+1, k+1, ix, jx, kx, gd);
  size_t m7 = _F_IDX_S3D(i+1, j+1, k+1, ix, jx, kx, gd);
  
  REAL_TYPE q[8];
  q[0] = f[ m0 ];
  q[1] = f[ m1 ];
  q[2] = f[ m2 ];
  q[3] = f[ m3 ];
  q[4] = f[ m4 ];
  q[5] = f[ m5 ];
  q[6] = f[ m6 ];
  q[7] = f[ m7 ];
  
  // 全ての点のSDF値が既に確定した値から計算できる
  if (   TEST_BIT(nrm[m0], STATE_BIT)
      && TEST_BIT(nrm[m1], STATE_BIT)
      && TEST_BIT(nrm[m2], STATE_BIT)
      && TEST_BIT(nrm[m3], STATE_BIT)
      && TEST_BIT(nrm[m4], STATE_BIT)
      && TEST_BIT(nrm[m5], STATE_BIT)
      && TEST_BIT(nrm[m6], STATE_BIT)
      && TEST_BIT(nrm[m7], STATE_BIT) )
  {
    s = getInterp(g, q);
    return true;
  }

  // 内挿できない場合
  Vec3r w[8];
  w[0] = getNrm9( nrm[ m0 ] );
  w[1] = getNrm9( nrm[ m1 ] );
  w[2] = getNrm9( nrm[ m2 ] );
  w[3] = getNrm9( nrm[ m3 ] );
  w[4] = getNrm9( nrm[ m4 ] );
  w[5] = getNrm9( nrm[ m5 ] );
  w[6] = getNrm9( nrm[ m6 ] );
  w[7] = getNrm9( nrm[ m7 ] );
    
  for (int l=0; l<8; l++) w[l].normalize();
  
  // 法線の逆方向
  v = -getInterp(g, w);
    
  return false;
}


// #################################################################
/*
 * @brief 格子点を走査し，交点の法線方向の参照点を求める
 * @param [in]     cut    量子化した交点
 * @param [in]     PL     Polylibインスタンスx
 * @param [in,out] bx     BCindex P
 * @param [in]     inner  内部のみの指定(true, default=false)
 */
int Geometry::getRefPointOnCut(const int* cutL,
                               const int* cutU,
                               MPIPolylib* PL,
                               int* bx,
                               bool inner)
{
#define PROBE_DEBUG_NEUMANN 1
  int ix, jx, kx, gd;

  ix = size[0];
  jx = size[1];
  kx = size[2];
  gd = guide;

  Vec3r org(originD);
  Vec3r pch(pitchD);
  Vec3r b(0.5);
  
  REAL_TYPE dh = pitchD[0];
  
#if PROBE_DEBUG_NEUMANN
  printf("DEBUG ========================\n");
  mark();
  vector<cut_probe> nrf;
#endif

  // ポリゴン情報へのアクセス
  vector<PolygonGroup*>* pg_roots = PL->get_root_groups();
  
  int Ncount = 0;

  int is,js,ks,ie,je,ke;
  
  if (inner) {
    is = js = ks = 2;
    ie = ix-1;
    je = jx-1;
    ke = kx-1;
  }
  else {
    is = js = ks = 1;
    ie = ix;
    je=  jx;
    ke = kx;
  }
  
  // 処理開始
  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
  for (int i=is; i<=ie; i++) {
    
    size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
    int posL = cutL[m];
    int posU = cutU[m];
    
    // FLUIDのセルのみが対象
    if ( !IS_FLUID(posL) ) continue;
    
    // oは有次元の(i,j,k)のセルセンター座標値
    Vec3r oi( (REAL_TYPE)i-0.5, (REAL_TYPE)j-0.5, (REAL_TYPE)k-0.5 );
    Vec3r o( org + oi * pch );
    
    // ポリゴンの探索範囲を基点を中心に10%の範囲 + 各方向
    Vec3r bx_min(o - pch * 0.1);
    Vec3r bx_max(o + pch * 0.1);
    
    
    for (int dir=0; dir<6; dir++)
    {
      // 登録済みの量子化交点距離
      int rc = (dir<3) ? getBitL9(posL, dir) : getBitU9(posU, dir);
      
      // dir方向の単位ベクトル
      Vec3r dvec;
      getDirVec(dvec, dir);
      
      // 交点の参照点と距離 （戻り値受け取り ）
      cut_probe cpr, nonref;
      cpr.dist = FLT_MAX;
      
      // 各方向に記録済みの交点があるとき
      int e = (dir<3) ? ensCutL(posL, dir) : ensCutU(posU, dir);
      if ( 1==e )
      {
        // 交点が在る場合にtrue
        bool flag=false;
        
        // ポリゴングループのループ
        vector<PolygonGroup*>::iterator it;
        for (it = pg_roots->begin(); it != pg_roots->end(); it++)
        {
          string m_pg = (*it)->get_name();     // グループラベル
          string m_bc = (*it)->get_type();     // 境界条件ラベル

          // サブドメイン内にポリゴンが存在する場合のみ処理する
          int ntria = (*it)->get_group_num_tria();  // ローカルのポリゴン数

          // ポリゴンが存在しなければ，スキップ
          if ( ntria == 0 ) continue;
          
          // Monitor属性のポリゴンはスキップ
          if ( !strcasecmp(m_bc.c_str(), "monitor")) continue;
          
          // 各方向に拡げたbbox
          getBboxCutDir(i, j, k, dir, dh, bx_min, bx_max);
          
          // false; ポリゴンがb-boxに一部でもかかる場合
          vector<Triangle*>* trias = PL->search_polygons(m_pg, bx_min, bx_max, false);
                  
          // ポリゴンがなければ，スキップ
          if ( trias->size() == 0 ) continue;
          
          // 対象ポリゴンから記録した最短距離の交点のポリゴンを選択し，法線を得る
          Vec3r nv;
          vector<Triangle*>::iterator it2;
          
          
          bool cflag = false;
          for (it2 = trias->begin(); it2 != trias->end(); it2++)
          {
            Vertex** tmp = (*it2)->get_vertex();
            Vec3r p[3];
            p[0] = *(tmp[0]);
            p[1] = *(tmp[1]);
            p[2] = *(tmp[2]);

            // 交点計算
            REAL_TYPE t, u, v, tn;
            if ( !TriangleIntersect(o, dvec, p[0], p[1], p[2], t, u, v) ) continue;
            tn = t / dh;
            if ( tn < 0.0 || 1.0 <= tn ) continue;
            
            // ポリゴンの交点距離が登録済みの距離と同じ場合
            if ( rc == quantize9( tn ) )
            {
              // ポリゴンの法線ベクトル
              Vec3r e1 = p[1] - p[0];
              Vec3r e2 = p[2] - p[0];
              Vec3r w = cross(e1, e2).normalize();
              
              // nvはポリゴンから基準点方向のベクトル（ポリゴン表面から流体側へ）
              nv = (dot(dvec, w) < 0) ? w : -w;
              
              // 法線を量子化した値に変換
              int nml=0;
              setNrm9(nml, nv);
              nv = getNrm9(nml);
              
              cflag = true;
              break; // it2 のループを抜ける
            }
          }
          delete trias;
          
          // 交点があるのに法線がえられなければエラー
          if ( !cflag ) return -1;
          
          
          // 参照点座標を計算（計算空間座標）
          Vec3r og( (REAL_TYPE)i, (REAL_TYPE)j, (REAL_TYPE)k);
          REAL_TYPE r = (dir<3) ? getCutL9(posL, dir) : getCutU9(posU, dir);
          
          // 内挿可能な参照点が得られた場合、flag=true で、cprに値が返る
          // 内挿できない点の場合にはnonrefに値がもどる
          if ( getMinRefPoint(og, dir, dvec, nv, r, cutL, cutU, cpr, nonref) ) {
            flag = true;
            break; // 見つかったら　pg_rootのループを抜ける
          }
          
        } // pg_root
        
        // 内挿可能な参照点がある場合にはリストに追加、得られない場合ノイマン条件を課す
        if (flag) {
          prb.push_back(cpr);
        }
        else {
          // 点(i,j,k)のdir方向のノイマンフラグをセット
          bx[m] = offBit( bx[m], BC_N_W+dir );
          Ncount++;
          
#if PROBE_DEBUG_NEUMANN
          nrf.push_back(nonref);
#endif
        } // flag
        
      } // ensCut()
    } // dir
  }}}
  
  Hostonly_ printf("Number of Neumann flag = %d\n", Ncount);

  delete pg_roots;
  
#if PROBE_DEBUG_NEUMANN
  std::ofstream ofs("neumann-pnt.scab", std::ios::out | std::ios::binary);
  if (!ofs) {
    printf("\tCan't open 'neumann-pnt.scab' file at rank=%d\n", myRank);
    Exit(0);
  }
  
  unsigned stp = 0;
  double tm = 0.0;
  unsigned csz = (unsigned)nrf.size();
  unsigned nvar = 3;
  
  //printf("Neumann flag points %d\n", csz);
  
  ofs.write((char*)&stp, sizeof(unsigned));
  ofs.write((char*)&tm,  sizeof(float));
  ofs.write((char*)&csz, sizeof(unsigned));
  ofs.write((char*)&nvar, sizeof(unsigned));
  
  for (unsigned i=0; i<csz; i++)
  {
    Vec3r p( org + (nrf[i].cp - b) * pch );
    ofs.write((char*)&p.x, sizeof(float));
    ofs.write((char*)&p.y, sizeof(float));
    ofs.write((char*)&p.z, sizeof(float));
    
    //Vec3r q( nrf[i].refp - nrf[i].cp );
    Vec3r q( org + (nrf[i].refp - b) * pch );
    ofs.write((char*)&q.x, sizeof(float));
    ofs.write((char*)&q.y, sizeof(float));
    ofs.write((char*)&q.z, sizeof(float));
    //printf("%f %f %f / %f %f %f\n", p.x, p.y, p.z, q.x, q.y, q.z);
  }
  
  ofs.close();
#endif
  
  return (int)prb.size();
}



// #################################################################
/*
 * @brief
 * @param [in]  i,j,k   インデクス
 * @param [in]  dir     方向
 * @param [in]  dh      格子幅
 * @param [out] bx_min  Bboxの最小値
 * @param [out] bx_max  Bboxの最大値
 */
void Geometry::getBboxCutDir(const int i, const int j, const int k, const int dir, const REAL_TYPE dh,
                             Vec3r& bx_min, Vec3r& bx_max)
{
  switch (dir)
  {
    case X_minus:
      bx_min.x -= dh;
      break;

    case X_plus:
      bx_max.x += dh;
      break;

    case Y_minus:
      bx_min.y -= dh;
      break;

    case Y_plus:
      bx_max.y += dh;
      break;

    case Z_minus:
      bx_min.z -= dh;
      break;

    case Z_plus:
      bx_max.z += dh;
      break;
  }
}
 

// #################################################################
/*
 * @brief 既に登録された交点距離と同じポリゴンの法線を得る
 * @param [in]     trias   候補ポリゴンリストへのポインタ（有次元）
 * @param [in]     ray_o   交点計算の基点位置ベクトル（有次元）
 * @param [in]     dir     テスト方向
 * @param [in]     d       単位方向ベクトル
 * @param [in]     dh      格子幅（有次元）
 * @param [in]     rc      登録済みの量子化交点距離
 * @param [out]    nv      ポリゴンから流体側への単位法線ベクトル
 * @retval true 交点が記録された距離と一致
 */
bool Geometry::getNV4CalculatedCut(vector<Triangle*>* trias,
                                   const Vec3r ray_o,
                                   const int dir,
                                   const Vec3r d,
                                   const REAL_TYPE dh,
                                   const int rc,
                                   Vec3r& nv)
{
  vector<Triangle*>::iterator it2;

  for (it2 = trias->begin(); it2 != trias->end(); it2++)
  {
    Vertex** tmp = (*it2)->get_vertex();
    Vec3r p[3];
    p[0] = *(tmp[0]);
    p[1] = *(tmp[1]);
    p[2] = *(tmp[2]);

    // 交点計算
    REAL_TYPE t, u, v, tn;
    if ( !TriangleIntersect(ray_o, d, p[0], p[1], p[2], t, u, v) ) continue;
    tn = t / dh;
    if ( tn < 0.0 || 1.0 <= tn ) continue;
    
    //printf("%d %d\n", rc, quantize9( tn ));
    // ポリゴンの交点距離を登録済みのものと比較
    if ( rc != quantize9( tn ) ) continue;
      
    // 対象ポリゴンが登録済みのものと同じ場合
    // ポリゴンの法線ベクトル
    Vec3r e1 = p[1] - p[0];
    Vec3r e2 = p[2] - p[0];
    Vec3r w = cross(e1, e2).normalize();
    
    // nvはポリゴンから基準点方向のベクトル（ポリゴン表面から流体側へ）
    nv = (dot(d, w) < 0) ? w : -w;
    return true;
  }

  return false;
}


// #################################################################
/*
 * @brief 方向ベクトルを返す
 * @param [in]     dir     テスト方向
 */
void Geometry::getDirVec(Vec3r& d, const int dir)
{
  switch (dir)
  {
    case X_minus:
      d.assign(-1.0, 0.0, 0.0);
      break;

    case X_plus:
      d.assign(1.0, 0.0, 0.0);
      break;

    case Y_minus:
      d.assign(0.0, -1.0, 0.0);
      break;

    case Y_plus:
      d.assign(0.0, 1.0, 0.0);
      break;

    case Z_minus:
      d.assign(0.0, 0.0, -1.0);
      break;

    case Z_plus:
      d.assign(0.0, 0.0, 1.0);
      break;
  }
}


// #################################################################
/*
 * @brief 交点の法線方向の内挿可能な参照点を計算する
 * @param [in]     o      交点計算の基点位置ベクトル（計算空間）
 * @param [in]     dir    交点のあるテスト方向
 * @param [in]     dvec   テスト方向の単位ベクトル
 * @param [in]     nv     ポリゴンから流体側への単位法線ベクトル
 * @param [in]     r      テスト方向の交点距離（計算空間 0<r<1）
 * @param [in]     cut    cut配列
 * @param [out]    cpr    内挿可能な交点の参照点座標と距離
 * @param [out]    nonref 内挿できない交点の参照点座標と距離
 */
bool Geometry::getMinRefPoint(const Vec3r o,
                              const int dir,
                              const Vec3r dvec,
                              const Vec3r nv,
                              const REAL_TYPE r,
                              const int* cutL,
                              const int* cutU,
                              cut_probe& cpr,
                              cut_probe& nonref)
{
  // 方向ベクトルと法線が直交 >> 交点なし
  if ( fabs(dot(nv, dvec)) < 1.0e-6 ) return false;
  
  Vec3r aa = o + dvec * r;
  Vec3r bb = aa + nv*DIST_AB;
  
  // 交点が登録されればc>0
  int c=0;
  
  cut_probe s1;
  
  // 最小距離を選ぶための初期値
  s1.dist = FLT_MAX;
  
  for (int k=-2; k<=2; k++) {
  for (int j=-2; j<=2; j++) {
  for (int i=-2; i<=2; i++) {
    
    Vec3r ss( (REAL_TYPE)i, (REAL_TYPE)j, (REAL_TYPE)k  );
    Vec3r tt = o + ss;
    
    if ( check_index_region(tt) )
    {
      // X, Y, Z方向
      for (int l=1; l<6; l+=2) {

        // 内挿可能な交点
        if ( decideInterp(tt, aa, bb, nv, l, cutL, cutU, s1) ) {
          if ( s1.dist<cpr.dist ) {
            cpr = s1;
            c++;
          }
        }
        // 内挿できない交点
        else {
          nonref = s1;
        }
      }
    } // check_index_region
  }}}

  return (c==0) ? false : true;
}

// #################################################################
/*
 * @brief 参照点の補間可能性をチェック
 * @param [in]     o      交差判定平面の基点位置ベクトル（計算空間）
 * @param [in]     aa     点Aの位置ベクトル
 * @param [in]     bb     点Bの位置ベクトル
 * @param [in]     nv     ポリゴンから流体側への単位法線ベクトル
 * @param [in]     face   交差判定平面の方向 (X_plus, Y_plus, Z_plus)
 * @param [in]     cut    cut配列
 * @param [out]    cpr    参照点座標と距離
 * @ret 内挿可能であればtrue
 */
bool Geometry::decideInterp(const Vec3r o,
                            const Vec3r aa,
                            const Vec3r bb,
                            const Vec3r nv,
                            const int face,
                            const int* cutL,
                            const int* cutU,
                            cut_probe& cpr)
{
  // 線分aa-bbの平面による内分比
  REAL_TYPE ratio;
  
  // 平面と線分に交点があれば処理
  Vec3r d;
  getDirVec(d, face);
  
  if ( divideLineByPlane(o, d, aa, bb, ratio) ) {
    REAL_TYPE f = ratio*DIST_AB;
    Vec3r rp = aa + nv*f;
    
    cpr.dist = f;
    cpr.refp = rp;
    cpr.nv = nv;
    cpr.cp = aa;
    
    // 内挿可能ならばtrue
    if ( judgeInterp(rp, face, cutL, cutU) ) return true;
  }
  return false;
}


// #################################################################
/*
 * @brief 平面と線分ABの交差判定
 * @param [in]     o      交差判定平面の基点位置ベクトル（計算空間）
 * @param [in]     n      交差判定平面の法線ベクトル
 * @param [in]     a      点Aの位置ベクトル
 * @param [in]     b      点Bの位置ベクトル
 * @param [out]    r      ABを平面で分割する内分比
 */
bool Geometry::divideLineByPlane(const Vec3r o,
                                 const Vec3r n,
                                 const Vec3r a,
                                 const Vec3r b,
                                 REAL_TYPE& r)
{
  Vec3r oa = a - o;
  Vec3r ob = b - o;
  
  REAL_TYPE f1 = dot(oa, n);
  REAL_TYPE f2 = dot(ob, n);
  r = fabs(f1);
  
  if ( f1 * f2 < 0.0 ) {
    r /= ( r + fabs(f2) );
    return true;
  }

  return false;
}


// #################################################################
/*
 * @brief 参照点の補間可能性をチェック
 * @param [in]     p      参照点座標（計算空間）
 * @param [in]     dir    方向（補間面の指定: X, Y, Z面）
 * @param [in]     cut    cut配列
 */
bool Geometry::judgeInterp(const Vec3r p,
                           const int dir,
                           const int* cutL,
                           const int* cutU)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  Vec3r q = floorVec(p);
  Vec3r zz = p - q;
  int i = (int)q.x;
  int j = (int)q.y;
  int k = (int)q.z;

  int posL_p = cutL[_F_IDX_S3D(i  , j  , k  , ix, jx, kx, gd)];
  int posU_p = cutU[_F_IDX_S3D(i  , j  , k  , ix, jx, kx, gd)];
  int posL_n = cutL[_F_IDX_S3D(i  , j+1, k  , ix, jx, kx, gd)];
  int posU_n = cutU[_F_IDX_S3D(i  , j+1, k  , ix, jx, kx, gd)];
  int posL_t = cutL[_F_IDX_S3D(i  , j  , k+1, ix, jx, kx, gd)];
  int posU_t = cutU[_F_IDX_S3D(i  , j  , k+1, ix, jx, kx, gd)];
  int posL_e = cutL[_F_IDX_S3D(i+1, j  , k  , ix, jx, kx, gd)];
  int posU_e = cutU[_F_IDX_S3D(i+1, j  , k  , ix, jx, kx, gd)];

  int cpz = ensCutU (posU_p, Z_plus);
  int cnz = ensCutU (posU_n, Z_plus);
  int cpn = ensCutU (posU_p, Y_plus);
  int ctn = ensCutU (posU_t, Y_plus);
  int cez = ensCutU (posU_e, Z_plus);
  int cpe = ensCutL (posL_p, X_plus);
  int cte = ensCutL (posL_t, X_plus);
  int cen = ensCutU (posU_e, Y_plus);
  int cne = ensCutL (posL_n, X_plus);
  
  switch (dir) {
    case X_minus:
    case X_plus:
      zz.x = 0.0;
      if      ( zz.y==0.0 && cpz==0 ) return true;
      else if ( zz.y==1.0 && cnz==0 ) return true;
      else if ( zz.z==0.0 && cpn==0 ) return true;
      else if ( zz.z==1.0 && ctn==0 ) return true;
      else if ( cpz==0 && cnz==0 && cpn==0 && ctn==0 ) return true;
      else return false;
      break;
      
    case Y_minus:
    case Y_plus:
      zz.y = 0.0;
      if      ( zz.x==0.0 && cpz==0 ) return true;
      else if ( zz.x==1.0 && cez==0 ) return true;
      else if ( zz.z==0.0 && cpe==0 ) return true;
      else if ( zz.z==1.0 && cte==0 ) return true;
      else if ( cpz==0 && cez==0 && cpe==0 && cte==0 ) return true;
      else return false;
      break;
      
    case Z_minus:
    case Z_plus:
      zz.z = 0.0;
      if      ( zz.x==0.0 && cpn==0 ) return true;
      else if ( zz.x==1.0 && cen==0 ) return true;
      else if ( zz.y==0.0 && cpe==0 ) return true;
      else if ( zz.y==1.0 && cne==0 ) return true;
      else if ( cpn==0 && cen==0 && cpe==0 && cne==0 ) return true;
      else return false;
      break;
  }

  return false;
}


// #################################################################
/*
 * @brief プローブの格納順序をメモリアクセス順にソートし、配列にコピー
 * @param [out]    cf   参照点座標（計算空間）
 * @param [out]    ds   交点から参照点までの距離（計算空間）
 */
void Geometry::sortProbe(REAL_TYPE* cf, REAL_TYPE* ds)
{
  
  //for (int i=0; i<(int)prb.size(); i++)
  //  printf("%6d %7.5f %7.5f %7.5f %7.5f\n", i,prb[i].dist,prb[i].refp.x,prb[i].refp.y,prb[i].refp.z);
    
  //sort(prb.begin(), prb.end(), cmp_k);
  //sort(prb.begin(), prb.end(), cmp_j);
  //sort(prb.begin(), prb.end(), cmp_i);
  
  REAL_TYPE rx = 0.0;
  for (int i=0; i<(int)prb.size(); i++) {
    //printf("%6d %7.5f %7.5f %7.5f %7.5f\n", i,prb[i].dist,prb[i].refp.x,prb[i].refp.y,prb[i].refp.z);
    if (rx < prb[i].dist) rx = prb[i].dist;
  }
  printf("Max length of  probe = %f\n", rx);
  
  for (int i=0; i<(int)prb.size(); i++)
  {
    ds[i]   = prb[i].dist;
    cf[i  ] = prb[i].refp.x;
    cf[i+1] = prb[i].refp.y;
    cf[i+2] = prb[i].refp.z;
  }
  
}


// #################################################################
/*
 * @brief プローブの格納順序をメモリアクセス順にソートし、配列にコピー
 * @param [out]    cf   参照点座標（計算空間）
 * @param [out]    ds   交点から参照点までの距離（計算空間）
 */
void Geometry::writeProbe()
{
  Vec3r org(originD);
  Vec3r pch(pitchD);
  Vec3r b(0.5);
  
  unsigned stp = 0;
  double tm = 0.0;
  unsigned csz = (unsigned)prb.size();
  unsigned nvar = 3;
  
  
  // 交点と法線 >> V-Isio plotArrows
  std::ofstream ofs("cut-nml.scab", std::ios::out | std::ios::binary);
  if (!ofs) {
    printf("\tCan't open 'cut-nml.scab' file at rank=%d\n", myRank);
    return;
  }
  ofs.write((char*)&stp, sizeof(unsigned));
  ofs.write((char*)&tm,  sizeof(float));
  ofs.write((char*)&csz, sizeof(unsigned));
  ofs.write((char*)&nvar, sizeof(unsigned));
  
  for (unsigned i=0; i<csz; i++)
  {
    Vec3r p( org + (prb[i].cp - b) * pch );
    
    ofs.write((char*)&p.x, sizeof(float));
    ofs.write((char*)&p.y, sizeof(float));
    ofs.write((char*)&p.z, sizeof(float));
    ofs.write((char*)&prb[i].nv.x, sizeof(float));
    ofs.write((char*)&prb[i].nv.y, sizeof(float));
    ofs.write((char*)&prb[i].nv.z, sizeof(float));
  }
  ofs.close();
  
  
  // 交点と参照点 >> V-Isio plotArrows
  std::ofstream ofs2("cut-ref.scab", std::ios::out | std::ios::binary);
  if (!ofs2) {
    printf("\tCan't open 'cut-ref.scab' file at rank=%d\n", myRank);
    return;
  }
  ofs2.write((char*)&stp, sizeof(unsigned));
  ofs2.write((char*)&tm,  sizeof(float));
  ofs2.write((char*)&csz, sizeof(unsigned));
  ofs2.write((char*)&nvar, sizeof(unsigned));
  
  for (unsigned i=0; i<csz; i++)
  {
    Vec3r p( org + (prb[i].cp - b) * pch );
    ofs2.write((char*)&p.x, sizeof(float));
    ofs2.write((char*)&p.y, sizeof(float));
    ofs2.write((char*)&p.z, sizeof(float));
    
    Vec3r q( prb[i].refp - prb[i].cp );
    // Vec3r q( org + (prb[i].refp - b) * pch );
    ofs2.write((char*)&q.x, sizeof(float));
    ofs2.write((char*)&q.y, sizeof(float));
    ofs2.write((char*)&q.z, sizeof(float));
  }
  ofs2.close();
}


// #################################################################
/*
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
#pragma omp parallel for reduction(+:c)
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
#pragma omp parallel for reduction(+:c)
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
    if ( paraMngr->Allreduce(&c_tmp, &c, 1, MPI_SUM, procGrp) != CPM_SUCCESS ) Exit(0);
  }

  return c;
}

