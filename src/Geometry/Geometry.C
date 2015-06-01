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
 * @file   Geometry.C
 * @brief  FFV geometry related functions
 * @author aics
 */

#include "Geometry.h"
#include "FBUtility.h"
#include <limits.h>

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
/**
 * @brief 各ノードのラベルのユニーク性を担保する
 * @param [in,out] tbl   ラベルリスト
 * @param [in,out] mid   識別子配列
 * @param [in]     Dsize 配列サイズ
 */
void Geometry::assureUniqueLabel(vector<int> tbl, int* mid, const int* Dsize)
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
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        const int dd = mid[m];
        if ( dd > 0 )
        {
          addLabel2List(tbl, dd);
        }
        
      }
    }
  }
  
  // 確認
  //for (vector<int>::iterator it=tbl.begin(); it != tbl.end(); ++it )
  //{
  //  fprintf(stderr, "%3d %3d\n", myRank, *it);
  //}
  
  int LabelSize = (int)tbl.size();
  //fprintf(stderr, "Rank = %6d : size = %d\n", myRank, LabelSize);
  
  
  // ラベル数を集める
  int* count = NULL;
  if ( !(count = new int[numProc]) )   Exit(0);
  for (int i=0; i<numProc; i++) count[i] = 0;
  
  // ラベル数の収集
  if ( numProc > 1 )
  {
    if ( paraMngr->Allgather(&LabelSize, 1, count, 1) != CPM_SUCCESS ) Exit(0);
  }
  else // serial
  {
    count[0] = LabelSize;
  }
  
  /*
  Hostonly_
  {
    printf("No of processes = %d\n", numProc);
    for (int i=0; i<numProc; i++)
    {
      printf("Rank = %6d : size = %d\n", i, count[i]);
    }
    printf("\n");
  }
  */
  
  // 各ランクに割り振る開始ラベルを計算 >> accum[]に開始番号がはいる
  int* accum = NULL;
  if ( !(accum = new int[numProc]) )   Exit(0);
  for (int i=0; i<numProc; i++) accum[i] = 0;
  
  
  // 開始ランクを求める
  int st = 0;
  for (int i=0; i<numProc; i++)
  {
    // ラベルがあれば1
    if (count[i] != 0)
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
      fprintf(stderr, "Error : Start rank of gathering label is greather than %d\n", numProc);
    }
    Exit(0);
  }
  
  
  // ラベルの積算開始は1
  accum[st] = 1;

  // 続いて、各プロセスの開始数をセット
  for (int i=st+1; i<numProc; i++)
  {
    accum[i] = accum[i-1] + count[i-1];
  }

  /*
  Hostonly_
  {
    for (int i=0; i<numProc; i++)
    {
      fprintf(stderr, "original >> rank=%5d : num. labels=%6d : begin=%6d\n", i, count[i], accum[i]);
    }
    printf("\n");
  }
  */
  
  // ワーク配列
  int* wrk = NULL;
  if ( !(wrk = Alloc::Int_S3D(size, guide)) ) Exit(0);
  
  
  // 参照用コピー
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        wrk[m] = mid[m];
        
      }
    }
  }
  
  // ラベル変更テーブル
  int* Dtable=NULL;
  
  // ラベルが存在する場合のみ実行
  if ( LabelSize > 0 )
  {
    if ( !(Dtable = new int[LabelSize*2]) )   Exit(0);
    for (int i=0; i<LabelSize*2; i++) Dtable[i]=0;
    
    // searchParint()で割り当てたラベル
    int c = 0;
    for (vector<int>::iterator it=tbl.begin(); it != tbl.end(); ++it )
    {
      Dtable[c*2] = *it;
      c++;
    }
    
    // 更新するラベル 各ランク内で昇順、全ランクを通して連続
    for (int i=0; i<LabelSize; i++)
    {
      Dtable[i*2+1] = accum[myRank] + i;
    }
    
    /* 確認
    for (int i=0; i<LabelSize; i++)
    {
      fprintf(stderr, "replaced >> rank=%5d : from=%6d  to=%6d\n", myRank, Dtable[i*2], Dtable[i*2+1]);
    }
     */
    
    // ラベルのふり直し
    for (int l=0; l<LabelSize; l++)
    {
      const int target = Dtable[l*2];
      const int replace= Dtable[l*2+1];
      
#pragma omp parallel for firstprivate(ix, jx, kx, gd, target, replace) schedule(static)
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          for (int i=1; i<=ix; i++) {
            
            size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
            const int dd = wrk[m];
            
            if ( dd == target )
            {
              mid[m] = replace;
            }
            
          }
        }
      }
      
    } // loop l
    
    // release
    if ( wrk ) delete [] wrk; wrk=NULL;
    
    
    
    // ラベルのアップデート
    vector<int> v;
    
    // 置換
    for (int l=0; l<LabelSize; l++)
    {
      v.push_back( Dtable[l*2+1] );
    }
    
    tbl.clear();
    tbl.resize(v.size());
    copy(v.begin(), v.end(), tbl.begin());
    
    /* 確認
    for (vector<int>::iterator it=tbl.begin(); it != tbl.end(); ++it )
    {
      fprintf(stderr, ">>> %3d %3d\n", myRank, *it);
    }
     */
  }
  
  
  if( Dtable ) { delete [] Dtable; Dtable=NULL; }
  if( count )  { delete [] count;  count=NULL; }
  if( accum )  { delete [] accum;  accum=NULL; }
}



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


// #################################################################
/**
 * @brief ガイドセルのIDをbcdからmidに転写
 * @param [in]  face  外部境界面方向
 * @param [in]  bcd   BCindex ID
 * @param [out] mid   識別子配列
 */
void Geometry::copyIDonGuide(const int face, const int* bcd, int* mid)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  
  switch (face)
  {
    case X_minus:
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          size_t m = _F_IDX_S3D(0, j, k, ix, jx, kx, gd);
          mid[m] = DECODE_CMP(bcd[m]);
        }
      }
      break;
      
      
    case X_plus:
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
      for (int k=1; k<=kx; k++) {
        for (int j=1; j<=jx; j++) {
          size_t m = _F_IDX_S3D(ix+1, j, k, ix, jx, kx, gd);
          mid[m] = DECODE_CMP(bcd[m]);
        }
      }
      break;
      
      
    case Y_minus:
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
      for (int k=1; k<=kx; k++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, 0, k, ix, jx, kx, gd);
          mid[m] = DECODE_CMP(bcd[m]);
        }
      }
      break;
      
      
    case Y_plus:
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
      for (int k=1; k<=kx; k++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, jx+1, k, ix, jx, kx, gd);
          mid[m] = DECODE_CMP(bcd[m]);
        }
      }
      break;
      
      
    case Z_minus:
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, j, 0, ix, jx, kx, gd);
          mid[m] = DECODE_CMP(bcd[m]);
        }
      }
      break;
      
    case Z_plus:
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static)
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          size_t m = _F_IDX_S3D(i, j, kx+1, ix, jx, kx, gd);
          mid[m] = DECODE_CMP(bcd[m]);
        }
      }
      break;
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
 * @param [in]      fp         ファイルポインタ
 * @param [in,out]  d_bcd      BCindex ID
 * @param [in]      d_bid      交点ID情報
 * @param [in]      m_NoMedium 媒質数
 * @param [in]      mat        MediumList
 * @param [in]      m_NoCompo  コンポーネント数
 * @param [in]      cmp        CompoList
 * @param [in,out]  d_mid      work array
 */
bool Geometry::fill(FILE* fp,
                    int* d_bcd,
                    const int* d_bid,
                    const int m_NoMedium,
                    const MediumList* mat,
                    const int m_NoCompo,
                    const CompoList* cmp,
                    int* d_mid)
{

  /*
   * 先に流体領域をフィルし、残った部分を固体領域としてフィルする
   * d_bid[] <= cmp[]のオーダーインデクス
   * d_bcd[] <= mat[]のオーダーインデクス
   */
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  
  // 全計算内部セル数
  unsigned long total_cell = (unsigned long)ix * (unsigned long)jx * (unsigned long)kx;
  
  if ( numProc > 1 )
  {
    unsigned long tmp_fc = total_cell;
    if ( paraMngr->Allreduce(&tmp_fc, &total_cell, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
  }
  
  
  // フィル対象のセル数
  unsigned long target_count = countCellB(d_bcd, 0);

  
  Hostonly_
  {
    fprintf(fp,"\tFill initialize -----\n\n");
    fprintf(fp,"\t\tTotal cell count       = %16ld\n", total_cell);
    fprintf(fp,"\t\tInitial target count   = %16ld\n", target_count);
    
    fprintf(fp,"\n\tFill -----\n\n");
    fprintf(fp,"\t\tFill mode for each dir.: %d %d %d\n", FillSuppress[0], FillSuppress[1], FillSuppress[2]);
  }
  
  
  // FLUIDでフィル ---------------------
  if ( !fill_connected(fp, d_bcd, d_bid, mat, m_NoMedium, FLUID, target_count) ) return false;
  
  
  if ( target_count == 0 ) return true;
  
  
  // 連結領域の同定 --------------------
  if ( !identifyConnectRegion(d_mid, d_bcd, d_bid, target_count, total_cell-target_count) ) {
    mark();
    return false;
  }
  
  
  // 交点IDをもつd_bcd[]の未ペイントのセルを、最頻値交点IDの媒質でペイントする -------------------
  
  unsigned long tmp;
  
  if ( !fillByModalCutID(d_bcd, d_bid, m_NoCompo, cmp, tmp) )
  {
    Hostonly_
    {
      fprintf(fp, "\tFailed to perform fillByModalCutID()\n");
      return false;
    }
  }
  
  target_count -= tmp;
  
  if ( tmp > 0 )
  {
    Hostonly_
    {
      fprintf(fp,"\t\tPainted by cut ID's medium        = %16ld\n", tmp);
      fprintf(fp,"\t\t               Remaining cells    = %16ld\n\n", target_count);
    }
  }
  
  
  
  // 未ペイントのセルをヒントのSOLIDでペイントする -------------------
  
  // 未ペイント（ID=0）のセルを検出
  unsigned long upc = countCellB(d_bcd, 0);
  
  Hostonly_
  {
    fprintf(fp, "\tUnpainted cells after fillByModalCutID() = %ld\n", upc);
  }
  
  if ( upc == 0 ) return true;
  
  
  // SOLIDでフィル
  //if ( !fill_connected(fp, d_bcd, d_bid, mat, m_NoMedium, SOLID, target_count) ) return false;
  
  fillByModalSolid(d_bcd, d_bid, m_NoCompo, cmp);
  
  // 未ペイント(ID=0)をカウント
  upc = countCellB(d_bcd, 0);
  
  
  // 未ペイントセルがある場合は、全周カットの可能性
  tmp = replaceIsolatedCell(d_bcd, d_bid, m_NoCompo, cmp);
  
  if ( tmp > 0 )
  {
    Hostonly_
    {
      fprintf(fp,"\t\tReplaced isolated cell = %16ld\n", tmp);
    }
  }
  
  
  // チェック
  upc = countCellB(d_bcd, 0);
  
  if ( upc != 0 )
  {
    Hostonly_
    {
      fprintf(fp,"\n\tFill operation is done, but still remains %ld unpainted cells.\n\n", upc);
    }
    return false;
  }
  
  return true;
}


// #################################################################
/* @brief 連続領域のフィル
 * @param [in]      fp         ファイルポインタ
 * @param [in,out]  d_bcd      BCindex ID
 * @param [in]      d_bid      交点ID情報
 * @param [in]      mat        MediumList
 * @param [in]      m_NoMedium 媒質数
 * @param [in]      fill_mode  フィルモード (SOLID | FLUID)
 * @param [in,out]  target_count ペイント対象のセル数
 * @retval success => true
 */
bool Geometry::fill_connected(FILE* fp,
                              int* d_bcd,
                              const int* d_bid,
                              const MediumList* mat,
                              const int m_NoMedium,
                              const int fill_mode,
                              unsigned long& target_count)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  unsigned long filled;       ///< フィルされた数

  
  // 媒質のヒントを与える
  for (int m=0; m<NoHint; m++)
  {
    // MediumListの格納番号
    int key = FBUtility::findIDfromLabel(mat, m_NoMedium, fill_table[m].medium);
    
    if ( mat[key].getState() == fill_mode )
    {
      
      if ( fill_table[m].kind == kind_outerface )
      {
        filled = fillSeedBcdOuter(d_bcd, fill_table[m].dir, key, d_bid);
      }
      else
      {
        filled = fillSeedBcdInner(d_bcd, fill_table[m].point, key);
      }
      
      if ( numProc > 1 )
      {
        unsigned long c_tmp = filled;
        if ( paraMngr->Allreduce(&c_tmp, &filled, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
      }
  
      // 内部点の場合にはfilled==1のはず
      if ( fill_table[m].kind == kind_point )
      {
        if ( filled != 1 )
        {
          Hostonly_
          {
            fprintf(fp,"\tIn case of Inner Fill Hint, filled cell must be 1.\n");
          }
          Exit(0);
        }
      }
      
    }
  }
  

  if ( numProc > 1 )
  {
    if ( paraMngr->BndCommS3D(d_bcd, ix, jx, kx, gd, gd) != CPM_SUCCESS ) Exit(0);
  }

  
  if ( filled == 0 )
  {
    Hostonly_
    {
      fprintf(fp,"\tNo cells painted by %s\n", (fill_mode==FLUID)?"FLUID":"SOLID");
    }
    return false;
  }
  
  Hostonly_
  {
    fprintf(fp,"\t\tPainted %s cells by hint       = %16ld\n", (fill_mode==FLUID)?"FLUID":"SOLID", filled);
  }
  
  
  // ペイントされたシードセル数をターゲットから差し引く
  target_count -= filled;
  
  
  Hostonly_
  {
    fprintf(fp,"\t\tRemaining cells to paint          = %16ld\n\n", target_count);
  }
  
  
  
  // 隣接するセルと同じfill_mode属性で接続している場合に隣接IDでフィル
  
  int c=0;
  unsigned long sum_filled = 0;   ///< フィルされた数の合計
  
  while (target_count > 0) {
    
    // 隣接媒質でフィルする
    filled = fillByBid(d_bcd, d_bid, mat, fill_mode);
    
    if ( numProc > 1 )
    {
      if ( paraMngr->BndCommS3D(d_bcd, ix, jx, kx, gd, gd) != CPM_SUCCESS ) Exit(0);
    }
    
    target_count -= filled;
    sum_filled   += filled;
    
    c++;
    
    if ( filled <= 0 ) break; // フィル対象がなくなったら終了
  }
  
  
  Hostonly_
  {
    fprintf(fp,"\t\tConnected fill Iteration          = %5d\n", c);
    fprintf(fp,"\t\t               Filled by %s    = %16ld\n", (fill_mode==FLUID)?"FLUID":"SOLID", sum_filled);
    fprintf(fp,"\t\t               Remaining cells    = %16ld\n\n", target_count);
  }
  
  return true;
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
 * @brief bid情報を元にフラッドフィル
 * @param [in,out] bcd         BCindex B
 * @param [in]     bid         交点ID（5ビット幅x6方向）
 * @param [in]     mat         MediumList
 * @param [in]     mode        フィルモード (SOLID | FLUID)
 * @param [in]     Dsize       サイズ
 * @note Symmetric fillにより反復回数を減少
 */
unsigned long Geometry::fillByBid(int* bcd,
                                  const int* bid,
                                  const MediumList* mat,
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
  
  int mode_x = FillSuppress[0]; // if 0, suppress connectivity evaluation
  int mode_y = FillSuppress[1];
  int mode_z = FillSuppress[2];
  
  unsigned long filled   = 0; ///< ペイントされた数
  
  int fill_mode = mode;
  
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd, sdw, sde, sds, sdn, sdb, sdt) \
                         firstprivate(mode_x, mode_y, mode_z, fill_mode) \
schedule(static) reduction(+:filled)
  
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        
#include "fill_bid_naive.h"
        
      }
    }
  }
  
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd, sdw, sde, sds, sdn, sdb, sdt) \
                         firstprivate(mode_x, mode_y, mode_z, fill_mode) \
schedule(static) reduction(+:filled)
  
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
    if ( paraMngr->Allreduce(&tmp, &filled, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
  }
  
  return filled;
}


// #################################################################
/**
 * @brief 流体媒質のフィルをbid情報を元に実行
 * @param [in,out] mid         work array
 * @param [in]     tgt_id      フィルするID
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
  
  
  // findPolygonInCell()により、既に、midにはポリゴンIDがストアされている
  // 検査対象セル{-1}の隣接6方向を見て、tgt_idと同じIDがあれば対象セルをtgt_idでペイント
  
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
/* @brief 未ペイントセルを周囲の交点IDの最頻値の媒質でフィル
 * @param [in,out] bcd       BCindex B
 * @param [in]     bid       境界ID
 * @param [in]     m_NoCompo コンポーネント数
 * @param [in]     cmp       CompoList
 * @param [out]    replaced  置換されたセル数
 * @retval true -> success
 */
bool Geometry::fillByModalCutID(int* bcd,
                                const int* bid,
                                const int m_NoCompo,
                                const CompoList* cmp,
                                unsigned long& replaced)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  int flag = -1;
  
  unsigned long c = 0; /// painted count
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd) shared(flag) schedule(static) reduction(+:c)
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
        
        // 対象セルが未ペイントの場合、かつ、交点をもつ場合
        if ( DECODE_CMP(bcd[m]) == 0  &&  qw+qe+qs+qn+qb+qt > 0 )
        {
          // 交点IDの最頻値
          int sd = FBUtility::find_mode_id(qw, qe, qs, qn, qb, qt, m_NoCompo);
          if ( sd == 0 ) flag = 1;
          
          // 交点IDの媒質番号を得る
          int key = cmp[sd].getMatodr();
          if ( key == 0 )
          {
            printf("(%3d %3d %3d) sd=%2d key=%2d : %2d %2d %2d %2d %2d %2d\n", i,j,k,sd, key, qw, qe, qs, qn, qb, qt);
            flag = 1;
          }

          setMediumID(bcd[m], key);
          c++;
        }
      }
    }
  }
  
  if ( flag == 1 ) return false;
  
  
  if ( numProc > 1 )
  {
    unsigned long c_tmp = c;
    if ( paraMngr->Allreduce(&c_tmp, &c, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
  }
  
  replaced = c;
  
  return true;
}


// #################################################################
/* @brief 未ペイントセルを周囲の媒質IDの固体最頻値でフィル
 * @param [in,out] bcd       BCindex B
 * @param [in]     bid       境界ID
 * @param [in]     m_NoCompo コンポーネント数
 * @param [in]     cmp       CompoList
 * @retval 置換されたセル数
 * @note 周囲の媒質IDの固体最頻値がゼロの場合には，境界IDで代用
 */
unsigned long Geometry::fillByModalSolid(int* bcd,
                                         const int* bid,
                                         const int m_NoCompo,
                                         const CompoList* cmp)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  unsigned long c = 0; /// painted count
  
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
          int sd = FBUtility::find_mode_id(qw, qe, qs, qn, qb, qt, m_NoCompo);
          int key = sd;
          
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
            sd = FBUtility::find_mode_id(qw, qe, qs, qn, qb, qt, m_NoCompo);
            if ( sd == 0 ) Exit(0); // 何かあるはず
            
            // 交点IDの媒質番号を得る
            key = cmp[sd].getMatodr();
            if ( key == 0 )
            {
              //printf("sd=%2d key=%2d : %2d %2d %2d %2d %2d %2d\n", sd, key, qw, qe, qs, qn, qb, qt);
              Exit(0);
            }
          }

          setMediumID(bcd[m_p], key);
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
 * @brief 流体セルをマイナス値でワーク配列にコピー
 * @param [in,out] mid         work array
 * @param [in]     bcd         BCindex B
 * @note bcd[]には流体セルのIDが記録されている，それ以外はゼロ
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
    if ( paraMngr->Allreduce(&c_tmp, &c, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
  }
  
  return c;
}


// #################################################################
/* @brief サブセルの未ペイントセルを周囲の媒質IDの固体最頻値でフィル
 * @param [in,out] mid       識別子配列
 * @param [in]     fluid_id  フィルをする流体ID
 * @param [in]     m_NoCompo コンポーネント数
 * @retval 置換されたセル数
 */
unsigned long Geometry::fillSubCellByModalSolid(int* smd,
                                                const int m_NoCompo,
                                                REAL_TYPE* svf,
                                                const MediumList* mat)
{
  int sdv = NumSuvDiv;
  int fid = FillID;
  
  unsigned long c = 0; /// painted count
  
#pragma omp parallel for firstprivate(sdv, fid) schedule(static) reduction(+:c)
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
        
        int dd = smd[m_p];
        
        int qw = smd[m_w];
        int qe = smd[m_e];
        int qs = smd[m_s];
        int qn = smd[m_n];
        int qb = smd[m_b];
        int qt = smd[m_t];
        
        
        // 対象セルが未ペイントの場合
        if ( dd == -1 )
        {
          int sd = FBUtility::find_mode_id(fid, qw, qe, qs, qn, qb, qt, m_NoCompo);
          
          // 周囲の媒質IDの固体最頻値がゼロの場合
          if ( sd == 0 ) Exit(0); // 何かあるはず
          smd[m_p] = sd;
          svf[m_p] = (mat[sd].getState()==FLUID) ? 1.0 : 0.0;
          c++;
        }
      }
    }
  }
  
  return c;
}


// #################################################################
/* @brief サブセルのSolid部分の値を代入
 * @param [in,out] smd    サブセル識別子配列
 * @param [in,out] svf    サブセル体積率
 * @param [in]     m_NoCompo コンポーネント数
 * @retval 置換されたセル数
 */
unsigned long Geometry::fillSubCellSolid(int* smd, REAL_TYPE* svf)
{
  int sd  = SeedID;
  int sdv = NumSuvDiv;
  
  unsigned long c = 0; /// painted count
  
#pragma omp parallel for firstprivate(sdv, sd) schedule(static) reduction(+:c)
  for (int k=1; k<=sdv; k++) {
    for (int j=1; j<=sdv; j++) {
      for (int i=1; i<=sdv; i++) {
        
        size_t m = _F_IDX_S3D(i  , j  , k  , sdv, sdv, sdv, 1);;

        if ( smd[m] != sd )
        {
          if (smd[m] == -1) Exit(0);
          svf[m] = 0.5;
          c++;
        }
      }
    }
  }
  
  return c;
}


// #################################################################
/**
 * @brief bcd[]の内部セルにシードIDをペイントする
 * @param [in,out] bcd    BCindex B
 * @param [in]     p      有次元座標値
 * @param [in]     target ペイントするIDのエントリ
 * @param [in]     Dsize  サイズ
 */
unsigned long Geometry::fillSeedBcdInner(int* bcd, const REAL_TYPE p[3], const int target, const int* Dsize)
{
  
  unsigned long filled = 0;
  
  // 領域内でなければ、フィルしない
  if (p[0]<originD[0] || p[0]>originD[0]+regionD[0] ||
      p[1]<originD[1] || p[1]>originD[1]+regionD[1] ||
      p[2]<originD[2] || p[2]>originD[2]+regionD[2] )
  {
    return filled;
  }
  
  /*
  printf("\n\t\tfill point [rank=%5d] = (%f, %f, %f) : region (%f %f %f)-(%f %f %f)\n\n",
          myRank, p[0], p[1], p[2],
          originD[0], originD[1], originD[2],
          originD[0]+regionD[0], originD[1]+regionD[1], originD[2]+regionD[2]
          );
  */
  
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
  
  int w[3];
  findIndex(p, w);
  
  // 計算したインデクスが内部領域から外れる場合にはフィルしない
  if (w[0]<1 || w[0]>ix||
      w[1]<1 || w[1]>jx ||
      w[2]<1 || w[2]>kx )
  {
    //printf("\trank=%d : %d %d %d\n", myRank, w[0], w[1], w[2]);
    return filled;
  }
  

  size_t m = _F_IDX_S3D(w[0], w[1], w[2], ix, jx, kx, gd);
  
  // 未ペイントのセルの場合
  if ( DECODE_CMP(bcd[m]) == 0 )
  {
    setMediumID(bcd[m], target);
    filled = 1;
  }
  
  return filled;
}


// #################################################################
/**
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
#pragma omp parallel for firstprivate(ix, jx, kx, gd, tg, mx_w) schedule(static) reduction(+:c)
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
#pragma omp parallel for firstprivate(ix, jx, kx, gd, tg, mx_e) schedule(static) reduction(+:c)
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
#pragma omp parallel for firstprivate(ix, jx, kx, gd, tg, mx_s) schedule(static) reduction(+:c)
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
#pragma omp parallel for firstprivate(ix, jx, kx, gd, tg, mx_n) schedule(static) reduction(+:c)
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
#pragma omp parallel for firstprivate(ix, jx, kx, gd, tg, mx_b) schedule(static) reduction(+:c)
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
#pragma omp parallel for firstprivate(ix, jx, kx, gd, tg, mx_t) schedule(static) reduction(+:c)
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
#pragma omp parallel for firstprivate(ix, jx, kx, gd, tg) schedule(static) reduction(+:c) collapse(2)
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
#pragma omp parallel for firstprivate(ix, jx, kx, gd, tg) schedule(static) reduction(+:c) collapse(2)
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
#pragma omp parallel for firstprivate(ix, jx, kx, gd, tg) schedule(static) reduction(+:c) collapse(2)
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
#pragma omp parallel for firstprivate(ix, jx, kx, gd, tg) schedule(static) reduction(+:c) collapse(2)
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
#pragma omp parallel for firstprivate(ix, jx, kx, gd, tg) schedule(static) reduction(+:c) collapse(2)
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
#pragma omp parallel for firstprivate(ix, jx, kx, gd, tg) schedule(static) reduction(+:c) collapse(2)
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
 * @brief ペイント対象のセル数をカバーするbboxを探す
 * @param [in,out] mid      work array
 * @param [out]    bbox     bounding box
 * @param [in]     Dsize    サイズ
 * @retval 対象セル数
 * @note 各ランクで逐次実行
 */
unsigned long Geometry::findBboxforSeeding(int* mid, int* bbox, const int* Dsize)
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
  
  // ノードローカルの値
  unsigned long c = 0;
  
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        
        if ( mid[m] == 0 )
        {
          c++;
          if ( bbox[0] > i ) bbox[0] = i;
          if ( bbox[1] < i ) bbox[1] = i;
          if ( bbox[2] > j ) bbox[2] = j;
          if ( bbox[3] < j ) bbox[3] = j;
          if ( bbox[4] > k ) bbox[4] = k;
          if ( bbox[5] < k ) bbox[5] = k;
        }
        
      }
    }
  }
  
  return c;
}



// #################################################################
/**
 * @brief list[]内の最頻値IDを求める
 * @param [in] m_sz      配列のサイズ
 * @param [in] list      ID配列
 * @param [in] m_NoCompo コンポーネント数
 * @note 候補がない場合には、0が戻り値
 */
int Geometry::find_mode(const int m_sz, const int* list, const int m_NoCompo)
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
/**
 * @brief サブセル内の最頻値IDを求める
 * @param [in] smd      ID配列
 * @param [in] m_NoCompo コンポーネント数
 * @note 候補がない場合には、0が戻り値
 */
int Geometry::find_mode_smd(const int* smd, const int m_NoCompo)
{
  int key[CMP_BIT_W]; ///< ID毎の頻度 @note ffv_Initialize() >> fill()でif ( NoCompo+1 > CMP_BIT_W )をチェック
  memset(key, 0, sizeof(int)*CMP_BIT_W);
  
  int sdv = NumSuvDiv;
  
  // smd[]は-1の可能性もある
#pragma omp parallel for firstprivate(sdv) schedule(static) collapse(3)
  for (int k=1; k<=sdv; k++) {
    for (int j=1; j<=sdv; j++) {
      for (int i=1; i<=sdv; i++) {
        size_t m = _F_IDX_S3D(i, j, k, sdv, sdv, sdv, 1);
        if (smd[m] >= 0)
        {
          key[ smd[m] ]++;
        }
      }
    }
  }
  
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
/* @brief セルに含まれるポリゴンを探索し、d_midに記録
 * @param [in,out] d_mid  識別子配列
 * @param [in]     PL     MPIPolylibのインスタンス
 * @param [in]     PG     PolygonPropertyクラス
 * @param [in] m_NoCompo  コンポーネント数
 */
unsigned long Geometry::findPolygonInCell(int* d_mid, MPIPolylib* PL, PolygonProperty* PG, const int m_NoCompo)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  
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
        
#if 0
        printf("\t\t[rank=%6d] (%4d %4d %4d) - (%4d %4d %4d) %s \n",
               myRank, wmin[0], wmin[1], wmin[2], wmax[0], wmax[1], wmax[2], m_pg.c_str());
#endif
        
        for (int k=wmin[2]; k<=wmax[2]; k++) {
          for (int j=wmin[1]; j<=wmax[1]; j++) {
            for (int i=wmin[0]; i<=wmax[0]; i++) {
              
              Vec3r bx_min(originD[0]+pitchD[0]*(REAL_TYPE)(i-1),
                           originD[1]+pitchD[1]*(REAL_TYPE)(j-1),
                           originD[2]+pitchD[2]*(REAL_TYPE)(k-1)); // セルBboxの対角座標
              Vec3r bx_max(originD[0]+pitchD[0]*(REAL_TYPE)i,
                           originD[1]+pitchD[1]*(REAL_TYPE)j,
                           originD[2]+pitchD[2]*(REAL_TYPE)k);     // セルBboxの対角座標
              
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
                
                int z = find_mode(polys, ary, m_NoCompo);
                if ( z == 0 ) Exit(0);

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
  printf("\n");
  
  delete pg_roots;
  
  
  // ID=0 をチェック
  unsigned long c = countCellM(d_mid, 0, true);
  
  if ( c != 0 )
  {
    Hostonly_ stamped_printf("\tID=0 was found in water-tight fill process. : %ld\n", c);
    Exit(0);
  }
  
  
  // count the number of replaced cells >> except ID=-1
  c = countCellM(d_mid, -1, false);
  
  return c;
}


// #################################################################
/**
 * @brief ペイントするシードセルを探す
 * @param [in,out] mid      work array
 * @param [in]     bbox     探索範囲
 * @param [in]     counter  リピート数
 * @param [in]     Dsize    サイズ
 * @retval シードがあれば1
 * @note ローカルノードに対するスレッド処理
 */
int Geometry::findSeedCells(int* mid,
                            const int* bbox,
                            const int counter,
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
  
  int ist = bbox[0];
  int ied = bbox[1];
  int jst = bbox[2];
  int jed = bbox[3];
  int kst = bbox[4];
  int ked = bbox[5];
  
  // 最大スレッド数 * リピート数(counter=0~)で基点をつくる
  int mt = omp_get_max_threads() * counter;
  
  int flag = 0;
  int c = 0;
  
#pragma omp parallel for schedule(static) reduction(+:c) \
firstprivate(ix, jx, kx, gd, mt, ist, ied, jst, jed, kst, ked, flag)
  for (int k=kst; k<=ked; k++) {
    for (int j=jst; j<=jed; j++) {
      for (int i=ist; i<=ied; i++) {
        
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        
        // 未ペイントセルの場合
        if ( mid[m]==0 && flag==0 )
        {
          mid[m] = mt + omp_get_thread_num() + 1;
          flag = 1;
          c++;
        }
        
      }
    }
  }
  
  return c;
}




// #################################################################
/*
 * @brief フィルパラメータを取得
 * @param [in] tpCntl       TextParser
 * @param [in] fp           file pointer to "condition.txt"
 * @param [in] Unit         DIMENSIONAL | NON_DIMENSIONAL
 * @param [in] RefL         代表長さ
 * @param [in] m_Nomedium   媒質数
 * @param [in] mat          MediumList
 */
void Geometry::getFillParam(TextParser* tpCntl,
                            FILE* fp,
                            const int Unit,
                            const REAL_TYPE RefL,
                            const int m_NoMedium,
                            const MediumList* mat)
{
  string str;
  string label_base, label_leaf, label;
  
  
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
    
    label = label_leaf + "/kind";
    if ( !(tpCntl->getInspectedValue(label, str )) )
    {
      Hostonly_ stamped_printf("\tParsing error : No '%s'\n", label.c_str());
      Exit(0);
    }
    
    
    if ( !strcasecmp(str.c_str(), "outerface"))
    {
      fill_table[m].kind = kind_outerface;
    }
    else if ( !strcasecmp(str.c_str(), "point"))
    {
      fill_table[m].kind = kind_point;
    }
    else
    {
      Hostonly_ stamped_printf("\tError : Invalid keyword [%s]\n", str.c_str());
      Exit(0);
    }
    
    
    if ( fill_table[m].kind == kind_outerface )
    {
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
        Hostonly_ printf("\tDefault 'X_minus' is set for Hint Of FillSeed direction\n");
      }
    }
    else
    {
      label = label_leaf + "/coordinate";

      REAL_TYPE tmp[3];
      if ( !tpCntl->getInspectedVector(label, tmp, 3) )
      {
        Hostonly_ stamped_printf("\tParsing error : No '%s'\n", label.c_str());
        Exit(0);
      }
      
      // 無次元の場合に有次元化する
      if (Unit == NONDIMENSIONAL )
      {
        // 有次元化
        for (int i=0; i<3; i++) {
          tmp[i] *= RefL;
        }
      }
      
      fill_table[m].point[0] = tmp[0];
      fill_table[m].point[1] = tmp[1];
      fill_table[m].point[2] = tmp[2];
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
    fprintf(fp,"\n----------\n");
    fprintf(fp,"\n\t>> Fill Hint\n\n");
    fprintf(fp,"\t No.        Kind       Medium   Direction / Coordinate(ND)\n");
    fprintf(fp,"\t--------------------------------------------------------------------------\n");
    for (int m=0; m<NoHint; m++)
    {
      fprintf(fp,"\t%3d %12s %12s   ", m+1, (fill_table[m].kind == kind_outerface)?"OuterFace":"Point", fill_table[m].medium.c_str());
      if (fill_table[m].kind == kind_outerface)
      {
        fprintf(fp,"Direction = %s\n", FBUtility::getDirection(fill_table[m].dir).c_str());
      }
      else
      {
        fprintf(fp,"(%12.6e, %12.6e, %12.6e)\n", fill_table[m].point[0], fill_table[m].point[1], fill_table[m].point[2]);
      }
    }
    fprintf(fp,"\n----------\n");
  }
  
  
  // FillMediumがMediumList中にあるかどうかをチェック
  for (int m=0; m<NoHint; m++)
  {
    
    if ( FBUtility::findIDfromLabel(mat, m_NoMedium, fill_table[m].medium) == 0 )
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
  
  
  
  /* フィル方向制御 (NOT mandatory)
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
  */
}


// #################################################################
/**
 * @brief 連結領域を同定する
 * @param [in] mid          work array
 * @param [in] bcd          BCindex B
 * @param [in] bid          交点ID
 * @param [in] paintable    固体領域でペイントすべきセル数
 * @param [in] filled_fluid 流体でペイント済みの数
 * @param [in] Dsize        配列サイズ
 */
bool Geometry::identifyConnectRegion(int* mid,
                                     const int* bcd,
                                     const int* bid,
                                     const unsigned long paintable,
                                     const unsigned long filled_fluid,
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
  
  
  // d_midをゼロで初期化し、流体領域を-1でペイント
  if ( fillFluidRegion(mid, bcd) != filled_fluid )
  {
    Hostonly_
    {
      fprintf(stderr, "\tFilled cells of fluid is not consistent.\n");
      return false;
    }
  }

  
  // シード点Bbox
  int bbox[6];
  bbox[0] = INT_MAX; // x-
  bbox[1] = INT_MIN; // x+
  bbox[2] = INT_MAX; // y-
  bbox[3] = INT_MIN; // y+
  bbox[4] = INT_MAX; // z-
  bbox[5] = INT_MIN; // z+
  
  // 各ランク毎のempty cellの数
  unsigned long nc = findBboxforSeeding(mid, bbox);
  
  //fprintf(stderr, "empty cell in rank [%6d] = %10ld : range = (%d - %d, %d - %d, %d - %d)\n",
  //       myRank, nc, bbox[0], bbox[1], bbox[2], bbox[3], bbox[4], bbox[5]);
  
  
  int counter = 0;

  while ( nc > 0 )
  {
    // 未ペイントセルをペイントするユニークラベルを割り振る
    int unique_id = (unsigned long)findSeedCells(mid, bbox, counter);
    
    nc -= unique_id;
    
    // 未ペイントセルをユニークラベルでペイントする
    nc -= searchPaint(mid, bbox, bid);
    
    // fprintf(stderr, "\tIteration=%3d / Empty cell=%16ld / unique_id = %d\n", counter+1, nc, unique_id);
    
    counter++;
    
    if ( counter > 1000 )
    {
      Hostonly_ fprintf(stderr, "May be error: Fill iteration is greater than 1000\n");
      return false;
    }
  }
  
  // この時点で、d_mid[]は全て非ゼロの値が入る。
  // 流体部分は-1、それ以外がunique_id
  // unique_idは各ランクでユニークであるが、全プロセスではユニークになっていない
  
  
  // 各ランクのラベル保持コンテナ
  vector<int> labelTable;
  
  while ( true )
  {
    // 各ランクのラベルをユニークに定める
    assureUniqueLabel(labelTable, mid);
    
    if ( numProc > 1 )
    {
      if ( paraMngr->BndCommS3D(mid, ix, jx, kx, gd, gd) != CPM_SUCCESS ) Exit(0);
    }
    
    // 接続リスト
    vector< vector<int> > cnnctTable;
    cnnctTable.resize(labelTable.size()); // 接続リストの値をいれる大きさ
    
    
    // 接続グラフを作成
    if ( makeConnectList(cnnctTable, labelTable, mid, bid) == 0 ) break;
    
    
    // 各ラベルの接続リストのうち、自ラベルも含めて最小値のラベルでペイント
    for (int l=0; l<cnnctTable.size(); l++)
    {
      vector<int>::iterator itr_min;
      itr_min = min_element(cnnctTable[l].begin(), cnnctTable[l].end());
      const int tmp = *itr_min;                    // 接続リストの最小値
      const int myl = labelTable[l];               // 自ラベル
      const int lst_min = (myl > tmp) ? tmp : myl; // 最小値
      
      paintConnectedLabel(mid, cnnctTable[l], lst_min);
      
      if ( numProc > 1 )
      {
        if ( paraMngr->BndCommS3D(mid, ix, jx, kx, gd, gd) != CPM_SUCCESS ) Exit(0);
      }
    }
    
  }
  
  
  return true;
}


// #################################################################
/**
 * @brief 各ラベルの接続リストを作成
 * @param [in,out] cnct   接続リスト
 * @param [in,out] label  ラベルリスト
 * @param [in]     mid    識別子配列
 * @param [in]     bid    交点ID
 * @param [in]     Dsize  配列サイズ
 */
int Geometry::makeConnectList(vector< vector<int> >& cnct,
                              vector<int>& label,
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

  
  // 各ラベルについて、接続リストを作成
  for (int l=0; l<cnct.size(); l++)
  {
    const int key = label[l];
    
    // serial loop
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
          
          // 隣接セル方向に対する境界ID
          int qw = getBit5(qq, 0);
          int qe = getBit5(qq, 1);
          int qs = getBit5(qq, 2);
          int qn = getBit5(qq, 3);
          int qb = getBit5(qq, 4);
          int qt = getBit5(qq, 5);
          
          const int dd = mid[m_p];
          
          // SOLIDセルの場合
          if ( dd > 0 )
          {
            int d_w = mid[m_w];
            int d_e = mid[m_e];
            int d_s = mid[m_s];
            int d_n = mid[m_n];
            int d_b = mid[m_b];
            int d_t = mid[m_t];
            
            // テスト方向にカットがなく、自セルのラベルでない
            if ( d_w != dd && qw==0 ) addLabel2List(cnct[l], d_w);
            if ( d_e != dd && qe==0 ) addLabel2List(cnct[l], d_e);
            if ( d_s != dd && qs==0 ) addLabel2List(cnct[l], d_s);
            if ( d_n != dd && qn==0 ) addLabel2List(cnct[l], d_n);
            if ( d_b != dd && qb==0 ) addLabel2List(cnct[l], d_b);
            if ( d_t != dd && qt==0 ) addLabel2List(cnct[l], d_t);
          }
          
        }
      }
    }
  }
  
  int total=0;
  
  for (int l=0; l<cnct.size(); l++)
  {
    total += cnct[l].size();
  }
  
  return total;
}


// #################################################################
/* @brief 距離の最小値を求める
 * @param [in,out] cut カットの配列
 * @param [in]     bid 境界IDの配列
 * @param [in]     fp  file pointer
 */
void Geometry::minDistance(const long long* cut, const int* bid, FILE* fp)
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
          
          if ( TEST_BC(bd) ) // カットがあるか，IDによる判定
          {
            const long long c = cut[m];
            
            for (int n=0; n<6; n++)
            {
              int tmp = getBit9(c, n);
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
    if ( paraMngr->Allreduce(&tmp, &global_min, 1, MPI_MIN) != CPM_SUCCESS ) Exit(0);
  }
  
  Hostonly_
  {
    fprintf(fp, "\n\tMinimum non-dimensional distance except on a center = %e\n\n", (REAL_TYPE)global_min/(REAL_TYPE)QT_9); // 9bit幅
  }
  
}



// #################################################################
/**
 * @brief 交点が定義点にある場合にそのポリゴンの媒質IDでフィルし、反対側を修正する
 * @param [in,out] bcd       BCindex B
 * @param [in,out] bid       境界ID（5ビット幅x6方向）
 * @param [in,out] cut       カット情報
 * @param [in]     m_NoCompo コンポーネント数
 * @param [in]     cmp       CompoList
 * @param [out]    fillcut   定義点上の交点IDでフィルした数
 * @param [out]    modopp    対向点の修正数

 * @param [in]     Dsize     サイズ
 * @retval フィルされたセル数
 */
void Geometry::paintCutOnPoint(int* bcd,
                               int* bid,
                               long long* cut,
                               const int m_NoCompo,
                               const CompoList* cmp,
                               unsigned long& fillcut,
                               unsigned long& modopp,
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
  
  
  unsigned long fc = 0;
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static) reduction(+:fc)
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        
        size_t m = _F_IDX_S3D(i  , j  , k  , ix, jx, kx, gd);
        
        int qq = bid[m];
        long long pos = cut[m];
        
        // いずれかの方向が定義点上の場合（交点があり、かつ距離がゼロ）
        if ( chkZeroCut(pos, X_minus) ||
             chkZeroCut(pos, X_plus ) ||
             chkZeroCut(pos, Y_minus) ||
             chkZeroCut(pos, Y_plus ) ||
             chkZeroCut(pos, Z_minus) ||
             chkZeroCut(pos, Z_plus ) )
        {
          // 交点IDの最頻値
          int sd = FBUtility::find_mode_id(getBit5(qq, X_minus),
                                           getBit5(qq, X_plus),
                                           getBit5(qq, Y_minus),
                                           getBit5(qq, Y_plus),
                                           getBit5(qq, Z_minus),
                                           getBit5(qq, Z_plus),
                                           m_NoCompo);
#if 0
          printf("(%3d %3d %3d) %3d %3d %3d %3d %3d %3d : %d %d %d %d %d %d\n",
                 i,j,k,
                 getBit9(pos, 0),
                 getBit9(pos, 1),
                 getBit9(pos, 2),
                 getBit9(pos, 3),
                 getBit9(pos, 4),
                 getBit9(pos, 5),
                 ensCut(pos, X_minus),
                 ensCut(pos, X_plus),
                 ensCut(pos, Y_minus),
                 ensCut(pos, Y_plus),
                 ensCut(pos, Z_minus),
                 ensCut(pos, Z_plus)
                 );
#endif
          if ( sd == 0 ) Exit(0);
          
          // 交点IDの媒質番号を得る
          int key = cmp[sd].getMatodr();
          if ( key == 0 )
          {
            //printf("sd=%2d key=%2d : %2d %2d %2d %2d %2d %2d\n", sd, key, qw, qe, qs, qn, qb, qt);
            Exit(0);
          }
          
          // セルを固体にする
          setMediumID(bcd[m], key);
          
          
          // 6方向とも交点ゼロにする
          setCut9(pos, 0, X_minus);
          setCut9(pos, 0, X_plus);
          setCut9(pos, 0, Y_minus);
          setCut9(pos, 0, Y_plus);
          setCut9(pos, 0, Z_minus);
          setCut9(pos, 0, Z_plus);
          cut[m] = pos;
          
          setBit5(qq, sd, X_minus);
          setBit5(qq, sd, X_plus);
          setBit5(qq, sd, Y_minus);
          setBit5(qq, sd, Y_plus);
          setBit5(qq, sd, Z_minus);
          setBit5(qq, sd, Z_plus);
          bid[m] = qq;
          
          fc++;
        }
        
      }
    }
  }
  
  if ( numProc > 1 )
  {
    unsigned long c_tmp = fc;
    if ( paraMngr->Allreduce(&c_tmp, &fc, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
  }
  
  fillcut = fc;
  
  
  if ( numProc > 1 )
  {
    if ( paraMngr->BndCommS3D(bcd, ix, jx, kx, gd, gd) != CPM_SUCCESS ) Exit(0);
    if ( paraMngr->BndCommS3D(bid, ix, jx, kx, gd, gd) != CPM_SUCCESS ) Exit(0);
    if ( paraMngr->BndCommS3D(cut, ix, jx, kx, gd, gd) != CPM_SUCCESS ) Exit(0);
  }
  

  
  
  // 反対側の修正
  fc = 0;
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static) reduction(+:fc)
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        
        size_t m_p = _F_IDX_S3D(i  , j  , k  , ix, jx, kx, gd);
        size_t m_w = _F_IDX_S3D(i-1, j,   k,   ix, jx, kx, gd);
        size_t m_e = _F_IDX_S3D(i+1, j,   k,   ix, jx, kx, gd);
        size_t m_s = _F_IDX_S3D(i,   j-1, k,   ix, jx, kx, gd);
        size_t m_n = _F_IDX_S3D(i,   j+1, k,   ix, jx, kx, gd);
        size_t m_b = _F_IDX_S3D(i,   j,   k-1, ix, jx, kx, gd);
        size_t m_t = _F_IDX_S3D(i,   j,   k+1, ix, jx, kx, gd);
        
        const int qq = bid[m_p];
        const long long pos = cut[m_p];
        
        // 定義点上に交点があり、反対側のセルから見て交点がない場合
        if ( chkZeroCut(pos, X_minus) > 0 && !ensCut(cut[m_w], X_plus) )
        {
          setBit5(bid[m_w], getBit5(qq, X_minus), X_plus);
          setCut9(cut[m_w], QT_9, X_plus);
          fc++;
        }
        
        if ( chkZeroCut(pos, X_plus) > 0  && !ensCut(cut[m_e], X_minus) )
        {
          setBit5(bid[m_e], getBit5(qq, X_plus), X_minus);
          setCut9(cut[m_e], QT_9, X_minus);
          fc++;
        }
        
        if ( chkZeroCut(pos, Y_minus) > 0 && !ensCut(cut[m_s], Y_plus) )
        {
          setBit5(bid[m_s], getBit5(qq, Y_minus), Y_plus);
          setCut9(cut[m_s], QT_9, Y_plus);
          fc++;
        }
        
        if ( chkZeroCut(pos, Y_plus) > 0  && !ensCut(cut[m_n], Y_minus) )
        {
          setBit5(bid[m_n], getBit5(qq, Y_plus), Y_minus);
          setCut9(cut[m_n], QT_9, Y_minus);
          fc++;
        }
        
        if ( chkZeroCut(pos, Z_minus) > 0 && !ensCut(cut[m_b], Z_plus) )
        {
          setBit5(bid[m_b], getBit5(qq, Z_minus), Z_plus);
          setCut9(cut[m_b], QT_9, Z_plus);
          fc++;
        }
        
        if ( chkZeroCut(pos, Z_plus) > 0  && !ensCut(cut[m_t], Z_minus) )
        {
          setBit5(bid[m_t], getBit5(qq, Z_plus), Z_minus);
          setCut9(cut[m_t], QT_9, Z_minus);
          fc++;
        }
      }
    }
  }
  
  if ( numProc > 1 )
  {
    unsigned long c_tmp = fc;
    if ( paraMngr->Allreduce(&c_tmp, &fc, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
  }
  
  modopp = fc;
  
  
  if ( numProc > 1 )
  {
    if ( paraMngr->BndCommS3D(bid, ix, jx, kx, gd, gd) != CPM_SUCCESS ) Exit(0);
    if ( paraMngr->BndCommS3D(cut, ix, jx, kx, gd, gd) != CPM_SUCCESS ) Exit(0);
  }
  
}



// #################################################################
/* @brief 接続している領域を指定ラベルでペイント
 * @param [in,out] mid      識別子配列
 * @param [in]     label    ラベルリスト
 * @param [in]     replace  置換するラベル
 * @param [in]     Dsize    サイズ
 */
void Geometry::paintConnectedLabel(int* mid, const vector<int>& label, const int replace, const int* Dsize)
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

  int lsize = label.size();
  int* lbl = new int[lsize];
  for (int l=0; l<lsize; l++) lbl[l] = label[l];
  
  int rep = replace;
  
#pragma omp parallel for schedule(static) firstprivate(ix, jx, kx, gd, lsize, rep)
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        
        for (int l=0; l<lsize; l++)
        {
          if ( lbl[l] == mid[m] ) mid[m] = rep;
        }
        
      }
    }
  }
  
  if ( lbl ) delete [] lbl; lbl = NULL;
}



// #################################################################
/**
 * @brief 交点計算を行い、量子化する
 * @param [in]      fp        ファイルポインタ
 * @param [out]    cut        量子化した交点
 * @param [out]    bid        交点ID
 * @param [in,out] bcd        セルID
 * @param [in]     m_NoCompo  コンポーネント数
 * @param [in]     cmp        CompoList
 * @param [in]     PL         Polylibインスタンス
 * @param [in]     PG         PolygonPropertyクラス
 */
void Geometry::quantizeCut(FILE* fp,
                           long long* cut,
                           int* bid,
                           int* bcd,
                           const int m_NoCompo,
                           const CompoList* cmp,
                           MPIPolylib* PL,
                           PolygonProperty* PG)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  
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
                  long long cc = cut[m];
                  
                  // 各方向の交点を評価、短い距離を記録する。新規記録の場合のみカウント
                  count += updateCut(ctr, X_minus, p[0], p[1], p[2], cc, bb, poly_id);
                  count += updateCut(ctr, X_plus,  p[0], p[1], p[2], cc, bb, poly_id);
                  count += updateCut(ctr, Y_minus, p[0], p[1], p[2], cc, bb, poly_id);
                  count += updateCut(ctr, Y_plus,  p[0], p[1], p[2], cc, bb, poly_id);
                  count += updateCut(ctr, Z_minus, p[0], p[1], p[2], cc, bb, poly_id);
                  count += updateCut(ctr, Z_plus,  p[0], p[1], p[2], cc, bb, poly_id);
                  
                  bid[m] = bb;
                  cut[m] = cc;
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
  
  Hostonly_ fprintf(fp, "\tquantize cut = %d\n\n", count);
  
  delete pg_roots;
  
  
  if ( numProc > 1 )
  {
    if ( paraMngr->BndCommS3D(bid, ix, jx, kx, gd, gd) != CPM_SUCCESS ) Exit(0);
    if ( paraMngr->BndCommS3D(cut, ix, jx, kx, gd, gd) != CPM_SUCCESS ) Exit(0);
  }
  

  
  // 修正
  
  // 定義点上に交点がある場合の処理 >> カットするポリゴンの媒質番号でフィルする
  unsigned long fill_cut, cm;
  
  paintCutOnPoint(bcd, bid, cut, m_NoCompo, cmp, fill_cut, cm);

  Hostonly_
  {
    if ( fill_cut > 0 )
    {
      fprintf(fp,"\n\tPaint cells which have cut on a center  = %16ld\n", fill_cut);
      fprintf(fp,  "\tModify neighbor cells owing to painting = %16ld\n", cm);
    }
  }
  
  // 最小値カット
  minDistance(cut, bid, fp);
  
}


// #################################################################
/** 
 * @brief 6方向にカットのあるセルを交点の媒質IDでフィルする
 * @param [in,out] bcd       BCindex
 * @param [in]     bid       交点ID
 * @param [in]     m_NoCompo コンポーネント数
 * @param [in]     cmp       CompoList
 */
unsigned long Geometry::replaceIsolatedCell(int* bcd,
                                            const int* bid,
                                            const int m_NoCompo,
                                            const CompoList* cmp)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  unsigned long replaced=0;
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd) schedule(static) reduction(+:replaced)
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        
        if ( DECODE_CMP(bcd[m]) == 0 )
        {
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
            // 交点IDの最頻値
            int sd = FBUtility::find_mode_id(qw, qe, qs, qn, qb, qt, m_NoCompo);
            if ( sd==0 ) Exit(0);
            
            // MediumListの格納番号
            int key = cmp[sd].getMatodr();
            //printf("%d %d %d : key id =%d\n", i,j,k,key);
            
            if ( key == 0 ) Exit(0);
            setMediumID(bcd[m], key);
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
/* @brief 未ペイントセルをシードフィル
 * @param [in,out] mid    識別子配列
 * @param [in]     bbox   対象領域
 * @param [in]     bid    交点ID
 * @param [in]     Dsize  サイズ
 * @retval 置換されたセル数
 * @note findSeedCells()と同じスレッド処理を行うこと
 */
unsigned long Geometry::searchPaint(int* mid, const int* bbox, const int* bid, const int* Dsize)
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
  
  int ist = bbox[0];
  int ied = bbox[1];
  int jst = bbox[2];
  int jed = bbox[3];
  int kst = bbox[4];
  int ked = bbox[5];
  
  unsigned long c = 0;
  
#pragma omp parallel for schedule(static) reduction(+:c)\
            firstprivate(ix, jx, kx, gd, ist, ied, jst, jed, kst, ked)
  for (int k=kst; k<=ked; k++) {
    for (int j=jst; j<=jed; j++) {
      for (int i=ist; i<=ied; i++) {
        
        size_t m_p = _F_IDX_S3D(i  , j  , k  , ix, jx, kx, gd);
        size_t m_e = _F_IDX_S3D(i+1, j,   k,   ix, jx, kx, gd);
        size_t m_w = _F_IDX_S3D(i-1, j,   k,   ix, jx, kx, gd);
        size_t m_n = _F_IDX_S3D(i,   j+1, k,   ix, jx, kx, gd);
        size_t m_s = _F_IDX_S3D(i,   j-1, k,   ix, jx, kx, gd);
        size_t m_t = _F_IDX_S3D(i,   j,   k+1, ix, jx, kx, gd);
        size_t m_b = _F_IDX_S3D(i,   j,   k-1, ix, jx, kx, gd);
        
        int qq = bid[m_p];
        int dd = 0;
        
        // 未ペイント
        if ( mid[m_p] == 0 )
        {
          int d_w = mid[m_w];
          int d_e = mid[m_e];
          int d_s = mid[m_s];
          int d_n = mid[m_n];
          int d_b = mid[m_b];
          int d_t = mid[m_t];
          
          if ( d_w != 0 && 0 == getBit5(qq, 0) )
          {
            dd = d_w;
          }
          else if ( d_e != 0 && 0 == getBit5(qq, 1) )
          {
            dd = d_e;
          }
          else if ( d_s != 0 && 0 == getBit5(qq, 2) )
          {
            dd = d_s;
          }
          else if ( d_n != 0 && 0 == getBit5(qq, 3) )
          {
            dd = d_n;
          }
          else if ( d_b != 0 && 0 == getBit5(qq, 4) )
          {
            dd = d_b;
          }
          else if ( d_t != 0 && 0 == getBit5(qq, 5) )
          {
            dd = d_t;
          }
        }
        
        if ( dd != 0 )
        {
          mid[m_p] = dd;
          c++;
        }
        
      }
    }
  }
  
  return c;
}



// #################################################################
/* @brief サブセルID配列をrefIDでペイント
 * @param [in,out] svf       サブセルの体積率
 * @param [in,out] smd       サブセルID配列
 * @param [in]     dir       ペイント開始方向
 * @param [in]     refID     ペイントID
 * @param [in]     refVf     ペイントする体積率
 * @note 呼び出し元がスレッド化してある場合には、single threadで実行
 */
int Geometry::SubCellFill(REAL_TYPE* svf,
                     int* smd,
                     const int dir,
                     const int refID,
                     const REAL_TYPE refVf)
{
  int filled = 0;
  int sdv = NumSuvDiv; // OpenMPのfirstprivateで使うためローカル変数にコピー
  int rD = refID;
  int rV = refVf;
  
  if ( dir == X_minus )
  {
    
#pragma omp parallel for firstprivate(sdv, rD, rV) reduction(+:filled) schedule(static) collapse(3)
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
            if ( smd[m_w] == rD ) ff++;
            if ( smd[m_e] == rD ) ff++;
            if ( smd[m_s] == rD ) ff++;
            if ( smd[m_n] == rD ) ff++;
            if ( smd[m_b] == rD ) ff++;
            if ( smd[m_t] == rD ) ff++;
            
            if ( ff>0 )
            {
              smd[m_p] = rD;
              svf[m_p] = rV;
              filled++;
            }
          }
        }
      }
    }
    
  }
  else if ( dir == X_plus)
  {
    
#pragma omp parallel for firstprivate(sdv, rD, rV) reduction(+:filled) schedule(static) collapse(3)
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
            if ( smd[m_w] == rD ) ff++;
            if ( smd[m_e] == rD ) ff++;
            if ( smd[m_s] == rD ) ff++;
            if ( smd[m_n] == rD ) ff++;
            if ( smd[m_b] == rD ) ff++;
            if ( smd[m_t] == rD ) ff++;
            
            if ( ff>0 )
            {
              smd[m_p] = rD;
              svf[m_p] = rV;
              filled++;
            }
          }
        }
      }
    }
    
  }
  else if ( dir == Y_minus)
  {
    
#pragma omp parallel for firstprivate(sdv, rD, rV) reduction(+:filled) schedule(static) collapse(3)
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
            if ( smd[m_w] == rD ) ff++;
            if ( smd[m_e] == rD ) ff++;
            if ( smd[m_s] == rD ) ff++;
            if ( smd[m_n] == rD ) ff++;
            if ( smd[m_b] == rD ) ff++;
            if ( smd[m_t] == rD ) ff++;
            
            if ( ff>0 )
            {
              smd[m_p] = rD;
              svf[m_p] = rV;
              filled++;
            }
          }
        }
      }
    }
    
  }
  else if ( dir == Y_plus)
  {

#pragma omp parallel for firstprivate(sdv, rD, rV) reduction(+:filled) schedule(static) collapse(3)
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
            if ( smd[m_w] == rD ) ff++;
            if ( smd[m_e] == rD ) ff++;
            if ( smd[m_s] == rD ) ff++;
            if ( smd[m_n] == rD ) ff++;
            if ( smd[m_b] == rD ) ff++;
            if ( smd[m_t] == rD ) ff++;
            
            if ( ff>0 )
            {
              smd[m_p] = rD;
              svf[m_p] = rV;
              filled++;
            }
          }
        }
      }
    }
    
  }
  else if ( dir == Z_minus)
  {
    
#pragma omp parallel for firstprivate(sdv, rD, rV) reduction(+:filled) schedule(static) collapse(3)
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
            if ( smd[m_w] == rD ) ff++;
            if ( smd[m_e] == rD ) ff++;
            if ( smd[m_s] == rD ) ff++;
            if ( smd[m_n] == rD ) ff++;
            if ( smd[m_b] == rD ) ff++;
            if ( smd[m_t] == rD ) ff++;
            
            if ( ff>0 )
            {
              smd[m_p] = rD;
              svf[m_p] = rV;
              filled++;
            }
          }
        }
      }
    }
    
  }
  else if ( dir == Z_plus)
  {
    
#pragma omp parallel for firstprivate(sdv, rD, rV) reduction(+:filled) schedule(static) collapse(3)
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
            if ( smd[m_w] == rD ) ff++;
            if ( smd[m_e] == rD ) ff++;
            if ( smd[m_s] == rD ) ff++;
            if ( smd[m_n] == rD ) ff++;
            if ( smd[m_b] == rD ) ff++;
            if ( smd[m_t] == rD ) ff++;
            
            if ( ff>0 )
            {
              smd[m_p] = rD;
              svf[m_p] = rV;
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
 * @param [in]     PL        MPIPolylibのインスタンス
 * @param [in]     m_NoCompo コンポーネント数
 * @retval ポリゴンを含むセル数
 * @note 呼び出し元がスレッド並列の場合、single threadで実行
 */
int Geometry::SubCellIncTest(REAL_TYPE* svf,
                             int* smd,
                             const int ip,
                             const int jp,
                             const int kp,
                             const Vec3r pch,
                             const string m_pg,
                             MPIPolylib* PL,
                             const int m_NoCompo)
{
  // プライマリセルの基点（有次元）
  Vec3r o(originD[0]+pitchD[0]*(REAL_TYPE)(ip-1),
          originD[1]+pitchD[1]*(REAL_TYPE)(jp-1),
          originD[2]+pitchD[2]*(REAL_TYPE)(kp-1));
  
  int pic = 0;
  
  int sdv = NumSuvDiv; // OpenMPのfirstprivateで使うためローカル変数にコピー
  
  // サブセルの体積率を評価（ポリゴンをもつサブセルのみ0.5）
  for (int k=1; k<=sdv; k++) {
    for (int j=1; j<=sdv; j++) {
      for (int i=1; i<=sdv; i++) {
        
        // サブセルのbbox
        Vec3r b1(o.x + pch.x * (REAL_TYPE)(i-1),
                 o.y + pch.y * (REAL_TYPE)(j-1),
                 o.z + pch.z * (REAL_TYPE)(k-1));
        
        Vec3r b2(o.x + pch.x * (REAL_TYPE)i,
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
          
          int z = find_mode(polys, ary, m_NoCompo);
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
 * @param [in]     m_NoCompo コンポーネント数
 * @note 呼び出し先でスレッド化している場合には、single threadで実行
 */
void Geometry::SubDivision(REAL_TYPE* svf,
                           int* smd,
                           const int ip,
                           const int jp,
                           const int kp,
                           int* d_mid,
                           const MediumList* mat,
                           REAL_TYPE* d_pvf,
                           const int m_NoCompo)
{
  // 外縁部にポリゴンがない面を探す
  int face_flag = 0;

  int sdv = NumSuvDiv; // OpenMPのfirstprivateで使うためローカル変数にコピー
  
  // SubCellIncTest()処理で、svfには、ポリゴンを持つサブセルに0.5が代入される
  // smdには、そのサブセルに存在するポリゴンIDの最頻値が代入されている
  
  // セルの6面をについて、ポリゴンを含むセルがあるかどうかを調べる
  // 調べた方向にポリゴンを含むセルがない場合にface_flagに方向を示すビットをエンコードする
  
  // X-
  int c = 0;
  
#pragma omp parallel for firstprivate(sdv) reduction(+:c) schedule(static) collapse(2)
  for (int k=1; k<=sdv; k++) {
    for (int j=1; j<=sdv; j++) {
      size_t m = _F_IDX_S3D(1, j, k, sdv, sdv, sdv, 1);
      if ( svf[m] > 0.0 ) c++; // ポリゴンを含むセルがある
    }
  }
  
  printf("X- : %d\n", c);
  if (c==0) face_flag |= (0x1 << X_minus);
  
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

  
  
  // プライマリセルの隣接ID（確定済み）を参照
  // 確定済みのセルは、SeedFillingプロセスでSeedIDでペイントしたセル
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
    
    // face_flagには空いた方向の面にビットフラグが設定されている
    // 外部境界面に接する場合には、空いている面側の参照すべきID(=r)がSeedIDでなく、
    // ガイドセルのIDとなる場合もある
    
    if ( TEST_BIT(face_flag, X_minus) )
    {
      m = _F_IDX_S3D(ip-1, jp, kp, ix, jx, kx, gd);
      r = d_mid[m];
      if ( r == SeedID )
      {
        refID = r;
        fillDir = X_minus;
      }
    }
    
    else if ( TEST_BIT(face_flag, X_plus) )
    {
      m = _F_IDX_S3D(ip+1, jp, kp, ix, jx, kx, gd);
      r = d_mid[m];
      if ( r == SeedID )
      {
        refID = r;
        fillDir = X_plus;
      }
    }
    
    else if ( TEST_BIT(face_flag, Y_minus) )
    {
      m = _F_IDX_S3D(ip, jp-1, kp, ix, jx, kx, gd);
      r = d_mid[m];
      if ( r == SeedID )
      {
        refID = r;
        fillDir = Y_minus;
      }
    }
    
    else if ( TEST_BIT(face_flag, Y_plus) )
    {
      m = _F_IDX_S3D(ip, jp+1, kp, ix, jx, kx, gd);
      r = d_mid[m];
      if ( r == SeedID )
      {
        refID = r;
        fillDir = Y_plus;
      }
    }
    
    else if ( TEST_BIT(face_flag, Z_minus) )
    {
      m = _F_IDX_S3D(ip, jp, kp-1, ix, jx, kx, gd);
      r = d_mid[m];
      if ( r == SeedID )
      {
        refID = r;
        fillDir = Z_minus;
      }
    }
    
    else if ( TEST_BIT(face_flag, Z_plus) )
    {
      m = _F_IDX_S3D(ip, jp, kp+1, ix, jx, kx, gd);
      r = d_mid[m];
      if ( r == SeedID )
      {
        refID = r;
        fillDir = Z_plus;
      }
    }
    
    
    Hostonly_
    {
      printf("\tP-cell (%3d %3d %3d) dir= %d  refID= %d : %d\n", ip, jp, kp, fillDir, refID, face_flag);
    }
    
  }

  
  if (fillDir == -1)
  {
    int s = find_mode_smd(smd, m_NoCompo);
    if ( s == 0 )
    {
      printf("Polygon ID in array smd is zero!\n");
      Exit(0);
    }
    size_t m = _F_IDX_S3D(ip, jp, kp, ix, jx, kx, gd);
    d_mid[m] = s;
    d_pvf[m] = (mat[s].getState()==FLUID) ? 1.0 : 0.0; ///< refIDの体積率
  }
  else
  {
    int target_count = sdv*sdv*sdv; ///< フィルの対象となるセル数
    
    if (refID == -1) Exit(0);
    
    
    // 開始面をペイント
    c = 0;
    
    if ( TEST_BIT(face_flag, X_minus) )
    {
#pragma omp parallel for firstprivate(sdv, refID) reduction(+:c) schedule(static) collapse(2)
      for (int k=1; k<=sdv; k++) {
        for (int j=1; j<=sdv; j++) {
          size_t m = _F_IDX_S3D(1, j, k, sdv, sdv, sdv, 1);
          if ( smd[m] == -1 )
          {
            smd[m] = refID;
            c++;
          }
          
          // カウントしない
          size_t m0 = _F_IDX_S3D(0, j, k, sdv, sdv, sdv, 1);
          if ( smd[m0] == -1 )
          {
            smd[m0] = refID;
          }
        }
      }
    }
    
    else if ( TEST_BIT(face_flag, X_plus) )
    {
#pragma omp parallel for firstprivate(sdv, refID) reduction(+:c) schedule(static) collapse(2)
      for (int k=1; k<=sdv; k++) {
        for (int j=1; j<=sdv; j++) {
          size_t m = _F_IDX_S3D(sdv, j, k, sdv, sdv, sdv, 1);
          if ( smd[m] == -1 )
          {
            smd[m] = refID;
            c++;
          }
          
          size_t m0 = _F_IDX_S3D(sdv+1, j, k, sdv, sdv, sdv, 1);
          if ( smd[m0] == -1 )
          {
            smd[m0] = refID;
          }
        }
      }
    }
    
    else if ( TEST_BIT(face_flag, Y_minus) )
    {
#pragma omp parallel for firstprivate(sdv, refID) reduction(+:c) schedule(static) collapse(2)
      for (int k=1; k<=sdv; k++) {
        for (int i=1; i<=sdv; i++) {
          size_t m = _F_IDX_S3D(i, 1, k, sdv, sdv, sdv, 1);
          if ( smd[m] == -1 )
          {
            smd[m] = refID;
            c++;
          }
          
          size_t m0 = _F_IDX_S3D(i, 0, k, sdv, sdv, sdv, 1);
          if ( smd[m0] == -1 )
          {
            smd[m0] = refID;
          }
        }
      }
    }
    
    else if ( TEST_BIT(face_flag, Y_plus) )
    {
#pragma omp parallel for firstprivate(sdv, refID) reduction(+:c) schedule(static) collapse(2)
      for (int k=1; k<=sdv; k++) {
        for (int i=1; i<=sdv; i++) {
          size_t m = _F_IDX_S3D(i, sdv, k, sdv, sdv, sdv, 1);
          if ( smd[m] == -1 )
          {
            smd[m] = refID;
            c++;
          }
          
          size_t m0 = _F_IDX_S3D(i, sdv+1, k, sdv, sdv, sdv, 1);
          if ( smd[0] == -1 )
          {
            smd[m0] = refID;
          }
        }
      }
    }
    
    else if ( TEST_BIT(face_flag, Z_minus) )
    {
#pragma omp parallel for firstprivate(sdv, refID) reduction(+:c) schedule(static) collapse(2)
      for (int j=1; j<=sdv; j++) {
        for (int i=1; i<=sdv; i++) {
          size_t m = _F_IDX_S3D(i, j, 1, sdv, sdv, sdv, 1);
          if ( smd[m] == -1 )
          {
            smd[m] = refID;
            c++;
          }
          
          size_t m0 = _F_IDX_S3D(i, j, 0, sdv, sdv, sdv, 1);
          if ( smd[m0] == -1 )
          {
            smd[m0] = refID;
          }
        }
      }
    }
    
    else if ( TEST_BIT(face_flag, Z_plus) )
    {
#pragma omp parallel for firstprivate(sdv, refID) reduction(+:c) schedule(static) collapse(2)
      for (int j=1; j<=sdv; j++) {
        for (int i=1; i<=sdv; i++) {
          size_t m = _F_IDX_S3D(i, j, sdv, sdv, sdv, sdv, 1);
          if ( smd[m] == -1 )
          {
            smd[m] = refID;
            c++;
          }
          
          size_t m0 = _F_IDX_S3D(i, j, sdv+1, sdv, sdv, sdv, 1);
          if ( smd[m0] == -1 )
          {
            smd[m0] = refID;
          }
        }
      }
    }
    
    target_count -= c;
    
    int filled = 0;       ///< フィルした数
    int sum_filled = 0;   ///< フィルした数の合計
    REAL_TYPE refVf = (mat[refID].getState()==FLUID) ? 1.0 : 0.0; ///< refIDの体積率
    
    c = 0;
    while (target_count > 0) {
      
      // 未ペイントで隣接セルがrefIDの場合、refID, refVfを代入
      filled = SubCellFill(svf, smd, fillDir, refID, refVf);
      
      target_count -= filled;
      sum_filled   += filled;
      c++;
      printf("\titr=%3d %8d %8d %8d\n", c, filled, sum_filled, target_count);
      
      if ( filled <= 0 ) break; // フィル対象がなくなったら、残りはPolyID
    }
    
    
    c = 0;
    while (target_count > 0) {
      
      // 未ペイントで隣接セルがrefIDの場合、refID, refVfを代入
      filled = fillSubCellByModalSolid(smd, m_NoCompo, svf, mat);
      
      target_count -= filled;
      c++;
      printf("\titr=%3d %8d %8d\n", c, filled, target_count);
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
    
    // 未ペイント部分がある場合 >> Solidが残り
    if ( flag > 0 )
    {
      filled = fillSubCellSolid(smd, svf);
    }
    if ( filled != flag )
    {
      printf("\tfilled(%d) does not agree wtih #flag(%d)\n", filled, flag);
      Exit(0);
    }


    
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
  
}



// #################################################################
/* @brief sub-sampling
 * @param [in]  fp        ファイルポインタ
 * @param [in]  mat       MediumList
 * @param [in]  d_mid     識別子配列
 * @param [out] d_pvf     体積率
 * @param [in]  PL        MPIPolylibのインスタンス
 * @param [in]  m_NoCompo コンポーネント数
 */
void Geometry::SubSampling(FILE* fp,
                           MediumList* mat,
                           int* d_mid,
                           REAL_TYPE* d_pvf,
                           MPIPolylib* PL,
                           const int m_NoCompo)
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
  
  // SeedFillingで確定したセルに対応する体積率を代入初期化、それ以外は-1.0
  // Outer fill => SeedID
  REAL_TYPE tmp = (mat[SeedID].getState()==FLUID) ? 1.0 : 0.0;
  filled = assignVF(SeedID, tmp, d_mid, d_pvf);
  
  Hostonly_
  {
    printf    ("\t\tVolume fraction for Seed ID (%3.1f) = %16ld  (%s)\n", tmp, filled, mat[SeedID].alias.c_str());
    fprintf(fp,"\t\tVolume fraction for Seed ID (%3.1f) = %16ld  (%s)\n", tmp, filled, mat[SeedID].alias.c_str());
  }
  
  /* Inner fill => FillID
  tmp = (mat[FillID].getState()==FLUID) ? 1.0 : 0.0;
  filled = assignVF(FillID, tmp, d_mid, d_pvf);
  
  Hostonly_
  {
    printf    ("\t\tVolume fraction for Fill ID (%3.1f) = %16ld  (%s)\n", tmp, filled, mat[FillID].alias.c_str());
    fprintf(fp,"\t\tVolume fraction for Fill ID (%3.1f) = %16ld  (%s)\n", tmp, filled, mat[FillID].alias.c_str());
  }
  */
  

  
  int sdv = NumSuvDiv; // OpenMPのfirstprivateで使うためローカル変数にコピー
  
  // sub-cell work array
  size_t nx = (sdv+2)*(sdv+2)*(sdv+2);
  REAL_TYPE* svf = new REAL_TYPE [nx]; // guide cell is 1 for each dir.
  int* smd = new int [nx];

  
  // Poly IDsに対して、サブセルテスト
  // SeedFilling()でd_midに確定したID(SeedID)が入っている
  // この三重ループはスレッド化しない >> スレッド化する場合には、ループ内で呼び出しているメソッドをチェックすること
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        
#pragma omp parallel for firstprivate(sdv) schedule(static) collapse(3)
        for (int kk=0; kk<=sdv+1; kk++) {
          for (int jj=0; jj<=sdv+1; jj++) {
            for (int ii=0; ii<=sdv+1; ii++) {
              size_t m = _F_IDX_S3D(ii, jj, kk, sdv, sdv, sdv, 1);
              svf[m] = 0.0;
              smd[m] = -1;
            }
          }
        }
        
        REAL_TYPE m_org[3], m_pit[3];
        
        // サブセルのファイル出力ヘッダ
        for (int l=0; l<3; l++) m_pit[l] = pitchD[l]/(REAL_TYPE)sdv;
        
        m_org[0] = originD[0]+pitchD[0]*(REAL_TYPE)(i-1);
        m_org[1] = originD[1]+pitchD[1]*(REAL_TYPE)(j-1);
        m_org[2] = originD[2]+pitchD[2]*(REAL_TYPE)(k-1);
        
        
        
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        int q = d_mid[m];
        
        // SeedID以外
        if ( q != SeedID )
        {
          vector<PolygonGroup*>* pg_roots = PL->get_root_groups();
          vector<PolygonGroup*>::iterator it;
          
          // ポリゴングループのループ
          for (it = pg_roots->begin(); it != pg_roots->end(); it++)
          {
            string m_pg = (*it)->get_name();          // グループラベル
            string m_bc = (*it)->get_type();          // 境界条件ラベル
            int ntria = (*it)->get_group_num_tria();  // ローカルのポリゴン数
            
            // 対象ポリゴンがある場合のみ
            if ( ntria > 0 )
            {
              // Monitor属性のポリゴンはスキップ
              if ( strcasecmp(m_bc.c_str(), "monitor"))
              {
                // サブセルのポリゴン含有テスト
                // svfには、ポリゴンを持つサブセルに0.5が代入される
                // smdには、そのサブセルに存在するポリゴンIDの最頻値
                int cp = SubCellIncTest(svf, smd, i, j, k, m_pit, m_pg, PL, m_NoCompo);
                //printf("%3d %3d %3d : %3d\n", i,j,k,cp);
              }
            }
          } // Polygon Group
          
          SubDivision(svf, smd, i, j, k, d_mid, mat, d_pvf, m_NoCompo);
          
          delete pg_roots;
        }
        
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
/* @brief シード点によるフィル
 * @param [in]  fp        ファイルポインタ
 * @param [in]  cmp       CompoList class
 * @param [in]  mat       MediumList
 * @param [in]  d_mid     識別子配列
 * @param [in]  PL        MPIPolylibのインスタンス
 * @param [in]  PG        PolygonPropertyクラス
 * @param [in]  m_NoCompo コンポーネント数
 * @note ここまで、d_bcdにはsetMonitorList()でモニタIDが入っている
 */
void Geometry::SeedFilling(FILE* fp,
                           CompoList* cmp,
                           MediumList* mat,
                           int* d_mid,
                           MPIPolylib* PL,
                           PolygonProperty* PG,
                           const int m_NoCompo)
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
    printf    ("\t\tInitial target count   = %16ld\n\n", target_count);
    fprintf(fp,"\t\tInitial target count   = %16ld\n\n", target_count);
  }
  
  
  // セルに含まれるポリゴンを探索し、d_midに記録
  Hostonly_
  {
    printf(    "\t\tPaint cells that contain polygons -----\n");
    fprintf(fp,"\t\tPaint cells that contain polygons -----\n");
  }
  
  sum_replaced = findPolygonInCell(d_mid, PL, PG, m_NoCompo);
  
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
  
  
  
  
  // デフォルトのヒントはX-
  Hostonly_
  {
    printf(    "\tHint of filling -----\n\n");
    fprintf(fp,"\tHint of filling -----\n\n");
    printf(    "\t\tSeeding Dir.           : %s\n", FBUtility::getDirection(FillSeedDir).c_str());
    fprintf(fp,"\t\tSeeding Dir.           : %s\n", FBUtility::getDirection(FillSeedDir).c_str());
    printf(    "\t\tFill medium of SEED    : %s (%d)\n", mat[SeedID].alias.c_str(), SeedID);
    fprintf(fp,"\t\tFill medium of SEED    : %s (%d)\n", mat[SeedID].alias.c_str(), SeedID);
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
    printf(    "\t\tPainted cells          = %16ld  (%s)\n", filled, mat[SeedID].alias.c_str());
    fprintf(fp,"\t\tPainted cells          = %16ld  (%s)\n", filled, mat[SeedID].alias.c_str());
    printf(    "\t\tRemaining target cells = %16ld\n\n", target_count);
    fprintf(fp,"\t\tRemaining target cells = %16ld\n\n", target_count);
  }
  
  
  
  
  // 未ペイントのターゲットセル(d_mid==-1)を、SeedIDでペイントする
  Hostonly_
  {
    printf(    "\tFill from outside by Seed ID -----\n\n");
    fprintf(fp,"\tFill from outside by Seed ID -----\n\n");
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
    printf(    "\t\t    Filled cells       = %16ld  (%s)\n", sum_filled, mat[SeedID].alias.c_str());
    fprintf(fp,"\t\t    Filled cells       = %16ld  (%s)\n", sum_filled, mat[SeedID].alias.c_str());
    printf(    "\t\t    Remaining cells    = %16ld\n\n", target_count);
    fprintf(fp,"\t\t    Remaining cells    = %16ld\n\n", target_count);
  }
  
  
  
  if ( target_count != 0 )
  {
    Hostonly_
    {
      printf(    "\tFill cells, which are inside of objects -----\n\n");
      fprintf(fp,"\tFill cells, which are inside of objects -----\n\n");
    }
    
    // 未ペイント（ID=-1）のセルを検出
    unsigned long upc = countCellM(d_mid, -1);
    
    if ( upc > 0 )
    {
      Hostonly_
      {
        printf(    "\t\tUnpainted cell         = %16ld\n", upc);
        fprintf(fp,"\t\tUnpainted cell         = %16ld\n", upc);
      }
    }
    
    if ( target_count != upc )
    {
      Exit(0);
    }
    
    
    // 未ペイントのセルに対して、最頻値の媒質でフィルする
    c = 0;
    sum_replaced = 0;
    
    while ( target_count > 0 ) {
      
      // @todo アルゴリズム再考
      // replaced = fillByModalSolid(d_mid, FillID, m_NoCompo);
      
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
      printf(    "\t\tFilled cells           = %16ld\n\n", sum_replaced);
      fprintf(fp,"\t\tFilled cells           = %16ld\n\n", sum_replaced);
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
  
}


// #################################################################
/**
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


/**
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
                             long long& cut,
                             int& bid,
                             const int pid)
{
  // 単位方向ベクトルと格子幅
  Vec3r d;
  REAL_TYPE pit;
  
  switch (dir)
  {
    case X_minus:
      d.assign(-1.0, 0.0, 0.0);
      pit = pitchD[0];
      break;
      
    case X_plus:
      d.assign(1.0, 0.0, 0.0);
      pit = pitchD[0];
      break;
      
    case Y_minus:
      d.assign(0.0, -1.0, 0.0);
      pit = pitchD[1];
      break;
      
    case Y_plus:
      d.assign(0.0, 1.0, 0.0);
      pit = pitchD[1];
      break;
      
    case Z_minus:
      d.assign(0.0, 0.0, -1.0);
      pit = pitchD[2];
      break;
      
    case Z_plus:
      d.assign(0.0, 0.0, 1.0);
      pit = pitchD[2];
      break;
  }
  
  // カウンタ
  unsigned count = 0;
  
  // 交点計算
  REAL_TYPE t, u, v;
  if ( !TriangleIntersect(ray_o, d, v0, v1, v2, t, u, v) ) return 0;

  /* 交点 >> 必要な場合に使う
   px = d.x * t + ray_o.x;
   py = d.y * t + ray_o.y;
   pz = d.z * t + ray_o.z;
   
   // 法線
   float fDat = 1.0 - u - v;
   n1 = n1 * fDat;
   n2 = n2 * U;
   n3 = n3 * V;
   n = n1 + n2 + n3;
   */

  // 格子幅で正規化
  REAL_TYPE tn = t / pit;
  //printf("t = %f\n", t);
  
  if ( tn < 0.0 || 1.0 < tn ) return 0;

  
  // 9bit幅の量子化
  int r = quantize9(tn);
  
  bool record = false;
  
  
  // 交点が記録されていない場合 >> 新規記録
  if ( ensCut(cut, dir) == 0 )
  {
    record = true;
    count = 1;
  }
  else // 交点が既に記録されている場合 >> 短い方を記録
  {
    if ( r < getBit9(cut, dir)) record = true;
  }
  
  
  if ( record )
  {
    setBit5(bid, pid, dir);
    setCut9(cut, r, dir);
    //printf("%10.6f %6d dir=%d id=%d\n", tn, r, dir, pid);
  }
  
  return count;
}

// #################################################################
/**
 * @brief セル数をカウント（デバッグ用）
 * @param [in] bcd     BCindex B
 * @param [in] m_id    検査するID
 * @param [in] bid     交点ID
 * @param [in] Dsize   サイズ
 * @note painted : m_id以外のセル数をカウント(false), m_idのセル数をカウント(true)
 */
unsigned long Geometry::debug_countCellB(const int* bcd, const int m_id, const int* bid, const int* Dsize)
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
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd, id) schedule(static) reduction(+:c)
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        if ( DECODE_CMP(bcd[m]) == id )
        {
          size_t m_e = _F_IDX_S3D(i+1, j,   k,   ix, jx, kx, gd);
          size_t m_w = _F_IDX_S3D(i-1, j,   k,   ix, jx, kx, gd);
          size_t m_n = _F_IDX_S3D(i,   j+1, k,   ix, jx, kx, gd);
          size_t m_s = _F_IDX_S3D(i,   j-1, k,   ix, jx, kx, gd);
          size_t m_t = _F_IDX_S3D(i,   j,   k+1, ix, jx, kx, gd);
          size_t m_b = _F_IDX_S3D(i,   j,   k-1, ix, jx, kx, gd);
          int zp = DECODE_CMP( bcd[m] );
          int zw = DECODE_CMP( bcd[m_w] );
          int ze = DECODE_CMP( bcd[m_e] );
          int zs = DECODE_CMP( bcd[m_s] );
          int zn = DECODE_CMP( bcd[m_n] );
          int zb = DECODE_CMP( bcd[m_b] );
          int zt = DECODE_CMP( bcd[m_t] );
          int qq = bid[m];
          int qw = getBit5(qq, 0);
          int qe = getBit5(qq, 1);
          int qs = getBit5(qq, 2);
          int qn = getBit5(qq, 3);
          int qb = getBit5(qq, 4);
          int qt = getBit5(qq, 5);
          Hostonly_ printf("(%d %d %d) : bid %d %d %d %d %d %d : bcd %d %d %d %d %d %d %d\n",
                           i,j,k,
                           qw, qe, qs, qn, qb, qt,
                           zp, zw, ze, zs, zn, zb, zt);
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
