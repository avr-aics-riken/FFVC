
// unused
/*

   // セルに含まれるポリゴンを探索し、d_midに記録
 #ifndef DISABLE_MPI
   unsigned long findPolygonInCell(int* d_mid,
                                   MPIPolylib* PL,
                                   PolygonProperty* PG);
 #else
   unsigned long findPolygonInCell(int* d_mid,
                                   Polylib* PL,
                                   PolygonProperty* PG);
 #endif
 
 // d_mid[]がtargetであるセルに対して、d_pvf[]に指定値valueを代入する
 unsigned long assignVF(const int target,
                        const REAL_TYPE value,
                        const int* d_mid,
                        REAL_TYPE* d_pvf);
 
 
 // シード点をmid[]にペイントする
 unsigned long fillSeedMid(int* mid,
                           const int face,
                           const int target,
                           const int* Dsize=NULL);
 
 
 // ペイント対象のセルを含むbboxを探す
 unsigned long findBboxforSeeding(const int* mid, int* bbox, const int* Dsize=NULL);
 
 
 // 未ペイントセルをtargetでフィル
 unsigned long fillByID(int* mid,
                        const int target,
                        const int* Dsize=NULL);
 
 
 // ガイドセルのIDをbcdからmidに転写
 void copyIDonGuide(const int face, const int* bcd, int* mid);
 
 
     // サブサンプリング
 #ifndef DISABLE_MPI
   void SubSampling(int* d_mid,
                    REAL_TYPE* d_pvf,
                    MPIPolylib* PL);
 #else
   void SubSampling(int* d_mid,
                    REAL_TYPE* d_pvf,
                    Polylib* PL);
 #endif
 
 
 // @brief 再分割数を設定
 // @param [in] num 再分割数
 void setSubDivision(int num)
 {
   if ( num < 1 ) Exit(0);
   NumSuvDiv = num;
 }
 
 // セル数をカウント（デバッグ用）
 unsigned long debug_countCellB(const int* bcd,
                                const int m_id,
                                const int* bid,
                                const int* Dsize=NULL);

   // サブセルのポリゴン含有テスト
 #ifndef DISABLE_MPI
   int SubCellIncTest(REAL_TYPE* svf,
                      int* smd,
                      const int ip,
                      const int jp,
                      const int kp,
                      const Vec3r pch,
                      const string m_pg,
                      MPIPolylib* PL);
 #else
 int SubCellIncTest(REAL_TYPE* svf,
                    int* smd,
                    const int ip,
                    const int jp,
                    const int kp,
                    const Vec3r pch,
                    const string m_pg,
                    Polylib* PL);
 #endif
 
 
   // ポリゴンの水密化
 #ifndef DISABLE_MPI
   void SeedFilling(int* d_mid,
                    MPIPolylib* PL,
                    PolygonProperty* PG);
 #else
   void SeedFilling(int* d_mid,
                    Polylib* PL,
                    PolygonProperty* PG);
 #endif
 
 
 // サブセルのペイント
 int SubCellFill(REAL_TYPE* svf,
                 int* smd,
                 const int dir,
                 const int refID,
                 const REAL_TYPE refVf);




 // sub-division
 void SubDivision(REAL_TYPE* svf,
                  int* smd,
                  const int ip,
                  const int jp,
                  const int kp,
                  int* d_mid,
                  REAL_TYPE* d_pvf);

 
 
 
 // 流体媒質のフィルをbid情報を元に実行
 unsigned long fillByMid(int* mid,
                         const int tgt_id,
                         const int* Dsize=NULL);
 
 // 未ペイントセルを周囲のbidの固体最頻値でフィル
 bool fillByModalSolid(int* bcd,
                       const int* bid,
                       unsigned long& replaced);
 
 
 // サブセルのSolid部分の値を代入
 unsigned long fillSubCellSolid(int* smd, REAL_TYPE* svf);
 
 // サブセルの未ペイントセルを周囲の媒質IDの固体最頻値でフィル
 unsigned long fillSubCellByModalSolid(int* smd, REAL_TYPE* svf);
 
 // サブセル内の最頻値IDを求める
 int find_mode_smd(const int* smd);
 
 // 固体領域をスレッド毎にIDでペイント
 void paintSolidRegion(int* mid, const int* Dsize=NULL);
 
 
 // 6方向にカットのあるセルを交点媒質でフィルする
 unsigned long replaceIsolatedCell(int* bcd, const int* bid);
 
 
 // @brief インデックスを(1,0,0)シフト
 // @param [in] index 元のインデクス
 // @param [in] h     シフト幅
 //
 inline Vec3r shift_E(const Vec3r index, const Vec3r h)
 {
   return Vec3r(index.x+h.x, index.y, index.z);
 }


 // @brief インデックスを(-1,0,0)シフト
 // @param [in] index 元のインデクス
 //  @param [in] h     シフト幅
 inline Vec3r shift_W(const Vec3r index, const Vec3r h)
 {
   return Vec3r(index.x-h.x, index.y, index.z  );
 }


 // @brief インデックスを(0,1,0)シフト
 // @param [in] index 元のインデクス
 // @param [in] h     シフト幅
 inline Vec3r shift_N(const Vec3r index, const Vec3r h)
 {
   return Vec3r(index.x, index.y+h.y, index.z);
 }



  * @brief インデックスを(0,-1,0)シフト
  * @param [in] index 元のインデクス
  * @param [in] h     シフト幅
 inline Vec3r shift_S(const Vec3r index, const Vec3r h)
 {
   return Vec3r(index.x, index.y-h.y, index.z);
 }



  * @brief インデックスを(0,0,1)シフト
  * @param [in] index 元のインデクス
  * @param [in] h     シフト幅
 inline Vec3r shift_T(const Vec3r index, const Vec3r h)
 {
   return Vec3r(index.x, index.y, index.z+h.z);
 }


  * @brief インデックスを(0,0,-1)シフト
  * @param [in] index 元のインデクス
  * @param [in] h     シフト幅
 inline Vec3r shift_B(const Vec3r index, const Vec3r h)
 {
   return Vec3r(index.x, index.y, index.z-h.z);
 }

 // bcd[]の内部セルにシードIDをペイントする
 unsigned long fillSeedBcdInner(int* bcd,
                                const REAL_TYPE p[3],
                                const int target,
                                const int* Dsize=NULL);
 
 
 
 
 */



// #################################################################
/*
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





#ifdef DISABLE_MPI
// #################################################################
/* @brief セルに含まれるポリゴンを探索し、d_midに記録
 * @param [in,out] d_mid  識別子配列
 * @param [in]     PL     Polylibのインスタンス
 * @param [in]     PG     PolygonPropertyクラス
 */
unsigned long Geometry::findPolygonInCell(int* d_mid, Polylib* PL, PolygonProperty* PG)
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

                int z = find_mode(polys, ary);
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
  unsigned long c = countCellM(d_mid, 0, "global", true);

  if ( c != 0 )
  {
    Hostonly_ stamped_printf("\tID=0 was found in water-tight fill process. : %ld\n", c);
    Exit(0);
  }


  // count the number of replaced cells >> except ID=-1
  c = countCellM(d_mid, -1, "global", false);

  return c;
}


#else // DISABLE_MPI
// #################################################################
/* @brief セルに含まれるポリゴンを探索し、d_midに記録
 * @param [in,out] d_mid  識別子配列
 * @param [in]     PL     MPIPolylibのインスタンス
 * @param [in]     PG     PolygonPropertyクラス
 */
unsigned long Geometry::findPolygonInCell(int* d_mid, MPIPolylib* PL, PolygonProperty* PG)
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

                int z = find_mode(polys, ary);
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
  unsigned long c = countCellM(d_mid, 0, "global", true);

  if ( c != 0 )
  {
    Hostonly_ stamped_printf("\tID=0 was found in water-tight fill process. : %ld\n", c);
    Exit(0);
  }


  // count the number of replaced cells >> except ID=-1
  c = countCellM(d_mid, -1, "global", false);

  return c;
}
#endif // DISABLE_MPI




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
    if ( paraMngr->Allreduce(&c_tmp, &c, 1, MPI_SUM, procGrp) != CPM_SUCCESS ) Exit(0);
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
    if ( paraMngr->Allreduce(&c_tmp, &c, 1, MPI_SUM, procGrp) != CPM_SUCCESS ) Exit(0);
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
unsigned long Geometry::findBboxforSeeding(const int* mid, int* bbox, const int* Dsize)
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


  if ( numProc > 1 )
  {
    unsigned long c_tmp = c;
    if ( paraMngr->Allreduce(&c_tmp, &c, 1, MPI_SUM, procGrp) != CPM_SUCCESS ) Exit(0);
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
    if ( paraMngr->Allreduce(&c_tmp, &c, 1, MPI_SUM, procGrp) != CPM_SUCCESS ) Exit(0);
  }

  return c;
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
    filled = fillByBid(d_bcd, d_bid, fill_mode);

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
/*
 * @brief STATEでFluidセル数をカウント
 * @param [in] cutL    cut  cutL/Uどちらでもよい
 * @param [in] Dsize   サイズ
 */
unsigned long Geometry::countCell(const int* cut, const int* Dsize)
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

#pragma omp parallel for reduction(+:c)
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        if ( IS_FLUID(cut[m]) ) c++;
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
    if ( paraMngr->Allreduce(&c_tmp, &c, 1, MPI_SUM, procGrp) != CPM_SUCCESS ) Exit(0);
  }

  return c;
}




// #################################################################
/**
 * @brief 6方向にカットのあるセルを交点の媒質IDでフィルする
 * @param [in,out] bcd       BCindex
 * @param [in]     bid       交点ID
 */
unsigned long Geometry::replaceIsolatedCell(int* bcd, const int* bid)
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
            int sd = FBUtility::find_mode_id(qw, qe, qs, qn, qb, qt, NoCompo);
            if ( sd==0 ) Exit(0);

            // MediumListの格納番号
            int key = cmp[sd].getMatodr();
            printf("rank %d : (%d %d %d) : key id =%d\n", myRank, i,j,k,key);

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
    if ( paraMngr->Allreduce(&c_tmp, &replaced, 1, MPI_SUM, procGrp) != CPM_SUCCESS ) Exit(0);
  }

  return replaced;
}



// #################################################################
/**
 * @brief サブセル内の最頻値IDを求める
 * @param [in] smd      ID配列
 * @note 候補がない場合には、0が戻り値
 */
int Geometry::find_mode_smd(const int* smd)
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
/* @brief サブセルの未ペイントセルを周囲の媒質IDの固体最頻値でフィル
 * @param [in,out] mid       識別子配列
 * @param [in]     fluid_id  フィルをする流体ID
 * @retval 置換されたセル数
 */
unsigned long Geometry::fillSubCellByModalSolid(int* smd, REAL_TYPE* svf)
{
  int sdv = NumSuvDiv;

  unsigned long c = 0; /// painted count

#pragma omp parallel for firstprivate(sdv) schedule(static) reduction(+:c)
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
          int sd = FBUtility::find_mode_id(cmp, qw, qe, qs, qn, qb, qt, NoCompo);

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
/* @brief 未ペイントセルを周囲の媒質IDの固体最頻値でフィル
 * @param [in,out] bcd       BCindex B
 * @param [in]     bid       境界ID
 * @param [out]    replaced  置換されたセル数
 * @retval true -> success
 * @note 周囲の媒質IDの固体最頻値がゼロの場合には，境界IDで代用
 */
bool Geometry::fillByModalSolid(int* bcd,
                                const int* bid,
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
          int sd = FBUtility::find_mode_id(cmp, qw, qe, qs, qn, qb, qt, NoCompo);
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
            sd = FBUtility::find_mode_id(cmp, qw, qe, qs, qn, qb, qt, NoCompo);
            if ( sd == 0 )
            {
              flag = 1;
            }
            else
            {
              // 交点IDの媒質番号を得る
              key = cmp[sd].getMatodr();
              if ( key == 0 )
              {
                printf("rank=%d : sd=%2d key=%2d : %2d %2d %2d %2d %2d %2d\n", myRank, sd, key, qw, qe, qs, qn, qb, qt);
                flag = 1;
              }
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
    if ( paraMngr->Allreduce(&c_tmp, &c, 1, MPI_SUM, procGrp) != CPM_SUCCESS ) Exit(0);
  }

  replaced = c;

  if ( flag == 1 ) return false;

  return true;
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
    if ( paraMngr->Allreduce(&tmp, &filled, 1, MPI_SUM, procGrp) != CPM_SUCCESS ) Exit(0);
  }

  return filled;
}




// #################################################################
/* @brief サブセル分割
 * @param [in,out] svf       サブセルの体積率
 * @param [in,out] smd       サブセルID配列
 * @param [in]     ip,jp,kp  プライマリセルインデクス
 * @param [in]     d_mid     識別子配列
 * @param [out]    d_pvf     体積率
 * @note 呼び出し先でスレッド化している場合には、single threadで実行
 */
void Geometry::SubDivision(REAL_TYPE* svf,
                           int* smd,
                           const int ip,
                           const int jp,
                           const int kp,
                           int* d_mid,
                           REAL_TYPE* d_pvf)
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
    int s = find_mode_smd(smd);
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
      filled = fillSubCellByModalSolid(smd, svf);

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





#ifdef DISABLE_MPI
// #################################################################
/* @brief シード点によるフィル
 * @param [in]  d_mid     識別子配列
 * @param [in]  PL        MPIPolylibのインスタンス
 * @param [in]  PG        PolygonPropertyクラス
 * @note ここまで、d_bcdにはsetMonitorList()でモニタIDが入っている
 */
void Geometry::SeedFilling(int* d_mid,
                           Polylib* PL,
                           PolygonProperty* PG)
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
    Hostonly_ fprintf(fpc, "\tSpecified Medium of filling fluid is not FLUID\n");
    Exit(0);
  }


  // 最初にフィル対象のセル数を求める >> 全計算内部セル数
  unsigned long total_cell = (unsigned long)ix * (unsigned long)jx * (unsigned long)kx;


  if ( numProc > 1 )
  {
    unsigned long tmp_fc = total_cell;
    if ( paraMngr->Allreduce(&tmp_fc, &total_cell, 1, MPI_SUM, procGrp) != CPM_SUCCESS ) Exit(0);
  }


  target_count = total_cell;

  Hostonly_
  {
    fprintf(fpc,"\t\tInitial target count   = %16ld\n\n", target_count);
  }


  // セルに含まれるポリゴンを探索し、d_midに記録
  Hostonly_
  {
    fprintf(fpc,"\t\tPaint cells that contain polygons -----\n");
  }

  sum_replaced = findPolygonInCell(d_mid, PL, PG);

  target_count -= sum_replaced;

  if ( numProc > 1 )
  {
    if ( paraMngr->BndCommS3D(d_mid, ix, jx, kx, gd, gd, procGrp) != CPM_SUCCESS ) Exit(0);
  }


  Hostonly_ // Up to here, d_mid={-1, Poly-IDs}
  {
    fprintf(fpc,"\t\t# of cells touch Polys = %16ld\n", sum_replaced);
    fprintf(fpc,"\t\tRemaining target cells = %16ld\n\n", target_count);
  }




  // デフォルトのヒントはX-
  Hostonly_
  {
    fprintf(fpc,"\tHint of filling -----\n\n");
    fprintf(fpc,"\t\tSeeding Dir.           : %s\n", FBUtility::getDirection(FillSeedDir).c_str());
    fprintf(fpc,"\t\tFill medium of SEED    : %s (%d)\n", mat[SeedID].alias.c_str(), SeedID);
  }

  filled = fillSeedMid(d_mid, FillSeedDir, SeedID);
  target_count -= filled;

  if ( filled == 0 )
  {
    Hostonly_
    {
      fprintf(fpc,"No cells painted\n");
    }
    Exit(0);
  }

  if ( numProc > 1 )
  {
    if ( paraMngr->BndCommS3D(d_mid, ix, jx, kx, gd, gd, procGrp) != CPM_SUCCESS ) Exit(0);
  }


  Hostonly_ // Up to here, d_mid={-1, Poly-IDs, SeedID}
  {
    fprintf(fpc,"\t\tPainted cells          = %16ld  (%s)\n", filled, mat[SeedID].alias.c_str());
    fprintf(fpc,"\t\tRemaining target cells = %16ld\n\n", target_count);
  }




  // 未ペイントのターゲットセル(d_mid==-1)を、SeedIDでペイントする
  Hostonly_
  {
    fprintf(fpc,"\tFill from outside by Seed ID -----\n\n");
  }

  int c=0; // iteration
  sum_filled = 0;

  while (target_count > 0) {

    // SeedIDで指定された媒質でフィルする
    filled = fillByMid(d_mid, SeedID);

    if ( numProc > 1 )
    {
      if ( paraMngr->BndCommS3D(d_mid, ix, jx, kx, gd, gd, procGrp) != CPM_SUCCESS ) Exit(0);
    }

    target_count -= filled;
    sum_filled   += filled;
    c++;

    if ( filled <= 0 ) break; // フィル対象がなくなったら終了
  }


  Hostonly_
  {
    fprintf(fpc,"\t\tIteration              = %5d\n", c);;
    fprintf(fpc,"\t\t    Filled cells       = %16ld  (%s)\n", sum_filled, mat[SeedID].alias.c_str());
    fprintf(fpc,"\t\t    Remaining cells    = %16ld\n\n", target_count);
  }



  if ( target_count != 0 )
  {
    Hostonly_
    {
      fprintf(fpc,"\tFill cells, which are inside of objects -----\n\n");
    }

    // 未ペイント（ID=-1）のセルを検出
    unsigned long upc = countCellM(d_mid, -1, "global");

    if ( upc > 0 )
    {
      Hostonly_
      {
        fprintf(fpc,"\t\tUnpainted cell         = %16ld\n", upc);
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
      // replaced = fillByModalSolid(d_mid, FillID);

      if ( numProc > 1 )
      {
        if ( paraMngr->BndCommS3D(d_mid, ix, jx, kx, gd, gd, procGrp) != CPM_SUCCESS ) Exit(0);
      }

      target_count -= replaced;
      sum_replaced += replaced;
      c++;

      if ( replaced <= 0 ) break;
    }

    Hostonly_
    {
      fprintf(fpc,"\t\tIteration              = %5d\n", c);
      fprintf(fpc,"\t\tFilled cells           = %16ld\n\n", sum_replaced);
    }



    // ID=-1をカウントしてチェック
    upc = countCellM(d_mid, -1, "global");

    if ( upc != 0 )
    {
      Hostonly_
      {
        fprintf(fpc,"\tFill operation is done, but still remains %ld unpainted cells.\n", upc);
      }
      Exit(0);
    }
  }

}


#else // DISABLE_MPI
// #################################################################
/* @brief シード点によるフィル
 * @param [in]  d_mid     識別子配列
 * @param [in]  PL        MPIPolylibのインスタンス
 * @param [in]  PG        PolygonPropertyクラス
 * @note ここまで、d_bcdにはsetMonitorList()でモニタIDが入っている
 */
void Geometry::SeedFilling(int* d_mid,
                           MPIPolylib* PL,
                           PolygonProperty* PG)
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
    Hostonly_ fprintf(fpc, "\tSpecified Medium of filling fluid is not FLUID\n");
    Exit(0);
  }


  // 最初にフィル対象のセル数を求める >> 全計算内部セル数
  unsigned long total_cell = (unsigned long)ix * (unsigned long)jx * (unsigned long)kx;


  if ( numProc > 1 )
  {
    unsigned long tmp_fc = total_cell;
    if ( paraMngr->Allreduce(&tmp_fc, &total_cell, 1, MPI_SUM, procGrp) != CPM_SUCCESS ) Exit(0);
  }


  target_count = total_cell;

  Hostonly_
  {
    fprintf(fpc,"\t\tInitial target count   = %16ld\n\n", target_count);
  }


  // セルに含まれるポリゴンを探索し、d_midに記録
  Hostonly_
  {
    fprintf(fpc,"\t\tPaint cells that contain polygons -----\n");
  }

  sum_replaced = findPolygonInCell(d_mid, PL, PG);

  target_count -= sum_replaced;

  if ( numProc > 1 )
  {
    if ( paraMngr->BndCommS3D(d_mid, ix, jx, kx, gd, gd, procGrp) != CPM_SUCCESS ) Exit(0);
  }


  Hostonly_ // Up to here, d_mid={-1, Poly-IDs}
  {
    fprintf(fpc,"\t\t# of cells touch Polys = %16ld\n", sum_replaced);
    fprintf(fpc,"\t\tRemaining target cells = %16ld\n\n", target_count);
  }




  // デフォルトのヒントはX-
  Hostonly_
  {
    fprintf(fpc,"\tHint of filling -----\n\n");
    fprintf(fpc,"\t\tSeeding Dir.           : %s\n", FBUtility::getDirection(FillSeedDir).c_str());
    fprintf(fpc,"\t\tFill medium of SEED    : %s (%d)\n", mat[SeedID].alias.c_str(), SeedID);
  }

  filled = fillSeedMid(d_mid, FillSeedDir, SeedID);
  target_count -= filled;

  if ( filled == 0 )
  {
    Hostonly_
    {
      fprintf(fpc,"No cells painted\n");
    }
    Exit(0);
  }

  if ( numProc > 1 )
  {
    if ( paraMngr->BndCommS3D(d_mid, ix, jx, kx, gd, gd, procGrp) != CPM_SUCCESS ) Exit(0);
  }


  Hostonly_ // Up to here, d_mid={-1, Poly-IDs, SeedID}
  {
    fprintf(fpc,"\t\tPainted cells          = %16ld  (%s)\n", filled, mat[SeedID].alias.c_str());
    fprintf(fpc,"\t\tRemaining target cells = %16ld\n\n", target_count);
  }




  // 未ペイントのターゲットセル(d_mid==-1)を、SeedIDでペイントする
  Hostonly_
  {
    fprintf(fpc,"\tFill from outside by Seed ID -----\n\n");
  }

  int c=0; // iteration
  sum_filled = 0;

  while (target_count > 0) {

    // SeedIDで指定された媒質でフィルする
    filled = fillByMid(d_mid, SeedID);

    if ( numProc > 1 )
    {
      if ( paraMngr->BndCommS3D(d_mid, ix, jx, kx, gd, gd, procGrp) != CPM_SUCCESS ) Exit(0);
    }

    target_count -= filled;
    sum_filled   += filled;
    c++;

    if ( filled <= 0 ) break; // フィル対象がなくなったら終了
  }


  Hostonly_
  {
    fprintf(fpc,"\t\tIteration              = %5d\n", c);;
    fprintf(fpc,"\t\t    Filled cells       = %16ld  (%s)\n", sum_filled, mat[SeedID].alias.c_str());
    fprintf(fpc,"\t\t    Remaining cells    = %16ld\n\n", target_count);
  }



  if ( target_count != 0 )
  {
    Hostonly_
    {
      fprintf(fpc,"\tFill cells, which are inside of objects -----\n\n");
    }

    // 未ペイント（ID=-1）のセルを検出
    unsigned long upc = countCellM(d_mid, -1, "global");

    if ( upc > 0 )
    {
      Hostonly_
      {
        fprintf(fpc,"\t\tUnpainted cell         = %16ld\n", upc);
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
      // replaced = fillByModalSolid(d_mid, FillID);

      if ( numProc > 1 )
      {
        if ( paraMngr->BndCommS3D(d_mid, ix, jx, kx, gd, gd, procGrp) != CPM_SUCCESS ) Exit(0);
      }

      target_count -= replaced;
      sum_replaced += replaced;
      c++;

      if ( replaced <= 0 ) break;
    }

    Hostonly_
    {
      fprintf(fpc,"\t\tIteration              = %5d\n", c);
      fprintf(fpc,"\t\tFilled cells           = %16ld\n\n", sum_replaced);
    }



    // ID=-1をカウントしてチェック
    upc = countCellM(d_mid, -1, "global");

    if ( upc != 0 )
    {
      Hostonly_
      {
        fprintf(fpc,"\tFill operation is done, but still remains %ld unpainted cells.\n", upc);
      }
      Exit(0);
    }
  }

}
#endif // DISABLE_MPI




#ifdef DISABLE_MPI
// #################################################################
/* @brief サブセルのポリゴン含有テスト
 * @param [in,out] svf       サブセルの体積率
 * @param [in,out] smd       サブセルID配列
 * @param [in]     ip,jp,kp  プライマリセルインデクス
 * @param [in]     pch       サブセルの格子幅
 * @param [in]     m_pg      ポリゴングループ名
 * @param [in]     PL        Polylibのインスタンス
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
                             Polylib* PL)
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

#else // DISABLE_MPI
// #################################################################
/* @brief サブセルのポリゴン含有テスト
 * @param [in,out] svf       サブセルの体積率
 * @param [in,out] smd       サブセルID配列
 * @param [in]     ip,jp,kp  プライマリセルインデクス
 * @param [in]     pch       サブセルの格子幅
 * @param [in]     m_pg      ポリゴングループ名
 * @param [in]     PL        MPIPolylibのインスタンス
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
                             MPIPolylib* PL)
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
#endif // DISABLE_MPI



#ifdef DISABLE_MPI
// #################################################################
/* @brief sub-sampling
 * @param [in]  d_mid     識別子配列
 * @param [out] d_pvf     体積率
 * @param [in]  PL        Polylibのインスタンス
 */
void Geometry::SubSampling(int* d_mid,
                           REAL_TYPE* d_pvf,
                           Polylib* PL)
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
    fprintf(fpc,"\t\tVolume fraction for Seed ID (%3.1f) = %16ld  (%s)\n", tmp, filled, mat[SeedID].alias.c_str());
  }

  /* Inner fill => FillID
   tmp = (mat[FillID].getState()==FLUID) ? 1.0 : 0.0;
   filled = assignVF(FillID, tmp, d_mid, d_pvf);

   Hostonly_
   {
   fprintf(fpc,"\t\tVolume fraction for Fill ID (%3.1f) = %16ld  (%s)\n", tmp, filled, mat[FillID].alias.c_str());
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
                int cp = SubCellIncTest(svf, smd, i, j, k, m_pit, m_pg, PL);
                //printf("%3d %3d %3d : %3d\n", i,j,k,cp);
              }
            }
          } // Polygon Group

          SubDivision(svf, smd, i, j, k, d_mid, d_pvf);

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

  //F->writeRawSPH(d_pvf, size, gd, 0, m_org, m_pit, sizeof(REAL_TYPE), "field");

}

#else // DISABLE_MPI
// #################################################################
/* @brief sub-sampling
 * @param [in]  d_mid     識別子配列
 * @param [out] d_pvf     体積率
 * @param [in]  PL        MPIPolylibのインスタンス
 */
void Geometry::SubSampling(int* d_mid,
                           REAL_TYPE* d_pvf,
                           MPIPolylib* PL)
{
  unsigned long filled = 0;       ///< フィルされた数

  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;


  // -1.0で初期化
#pragma omp parallel for
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
    fprintf(fpc,"\t\tVolume fraction for Seed ID (%3.1f) = %16ld  (%s)\n", tmp, filled, mat[SeedID].alias.c_str());
  }

  /* Inner fill => FillID
   tmp = (mat[FillID].getState()==FLUID) ? 1.0 : 0.0;
   filled = assignVF(FillID, tmp, d_mid, d_pvf);

   Hostonly_
   {
   fprintf(fpc,"\t\tVolume fraction for Fill ID (%3.1f) = %16ld  (%s)\n", tmp, filled, mat[FillID].alias.c_str());
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

#pragma omp parallel for collapse(3)
        for (int kk=0; kk<=sdv+1; kk++) {
          for (int jj=0; jj<=sdv+1; jj++) {
            for (int ii=0; ii<=sdv+1; ii++) {
              size_t m = _F_IDX_S3D(ii, jj, kk, sdv, sdv, sdv, 1);
              svf[m] = 0.0;
              smd[m] = -1;
            }
          }
        }

        REAL_TYPE m_pit[3];

        // サブセルのファイル出力ヘッダ
        for (int l=0; l<3; l++) m_pit[l] = pitchD[l]/(REAL_TYPE)sdv;

        //REAL_TYPE m_org[3];
        //m_org[0] = originD[0]+pitchD[0]*(REAL_TYPE)(i-1);
        //m_org[1] = originD[1]+pitchD[1]*(REAL_TYPE)(j-1);
        //m_org[2] = originD[2]+pitchD[2]*(REAL_TYPE)(k-1);



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
                SubCellIncTest(svf, smd, i, j, k, m_pit, m_pg, PL);
              }
            }
          } // Polygon Group

          SubDivision(svf, smd, i, j, k, d_mid, d_pvf);

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

  //F->writeRawSPH(d_pvf, size, gd, 0, m_org, m_pit, sizeof(REAL_TYPE), "field");

}
#endif // DISABLE_MPI


// #################################################################
/* @brief 交点情報のグリフを生成する
 * @param [in] cut   カットの配列
 * @param [in] bid   境界IDの配列
 * @param [in] fp    file pointer
 * @param [in] m_st  範囲指定　デバッグ用
 * @param [in] m_ed  範囲指定　デバッグ用
 */
void FFV::generateGlyph(const int* cutL, const int* cutU, const int* bid, FILE* fp, int* m_st, int* m_ed)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  // 交点の総数を求める
  unsigned global_cut=0;  /// 全カット数
  unsigned local_cut=0;   /// 担当プロセスのカット数
  unsigned g=0;

  int st[3], ed[3];

  if ( m_st == NULL )
  {
    st[0] = 1;
    st[1] = 1;
    st[2] = 1;
    ed[0] = ix;
    ed[1] = jx;
    ed[2] = kx;
  }
  else
  {
    st[0] = m_st[0];
    st[1] = m_st[1];
    st[2] = m_st[2];
    ed[0] = m_ed[0];
    ed[1] = m_ed[1];
    ed[2] = m_ed[2];
  }


  for (int k=st[2]; k<=ed[2]; k++) {
    for (int j=st[1]; j<=ed[1]; j++) {
      for (int i=st[0]; i<=ed[0]; i++) {

        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        int qq = bid[m];

        if ( IS_CUT(qq) ) // カットがあるか，IDによる判定
        {
          if (getBit5(qq, 0) != 0) g++;
          if (getBit5(qq, 1) != 0) g++;
          if (getBit5(qq, 2) != 0) g++;
          if (getBit5(qq, 3) != 0) g++;
          if (getBit5(qq, 4) != 0) g++;
          if (getBit5(qq, 5) != 0) g++;
        }

      }
    }
  }

  global_cut = local_cut = g;

  if ( numProc > 1 )
  {
    unsigned tmp = global_cut;
    if ( paraMngr->Allreduce(&tmp, &global_cut, 1, MPI_SUM, procGrp) != CPM_SUCCESS ) Exit(0);
  }

  Hostonly_
  {
    printf("\tGlyph : Number of Cut points = %u\n", global_cut);
    fprintf(fp, "\tGlyph : Number of Cut points = %u\n", global_cut);
  }


  // 格子幅（有次元）
  float m_pch[3] = {
    (float)pitch[0]*C.RefLength,
    (float)pitch[1]*C.RefLength,
    (float)pitch[2]*C.RefLength
  };

  // サブドメインの起点座標（有次元）
  float m_org[3] = {
    (float)origin[0]*C.RefLength,
    (float)origin[1]*C.RefLength,
    (float)origin[2]*C.RefLength
  };


  // ポリゴンをストアする配列を確保
  Glyph glyph(m_pch, m_org, local_cut, myRank);


  // グリフの生成モード
  bool inner_only = false;
  if (C.Hide.GlyphOutput == 2) inner_only=true;

  // カット点毎にグリフのポリゴン要素を生成し，配列にストア
  Vec3i idx;

  for (int k=st[2]; k<=ed[2]; k++) {
    for (int j=st[1]; j<=ed[1]; j++) {
      for (int i=st[0]; i<=ed[0]; i++) {

        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        int qq = bid[m];

        if ( IS_CUT(qq) ) // カットがあるか，IDによる判定
        {
          const int posL = cutL[m];
          const int posU = cutU[m];

          idx.assign(i, j, k);

          if ( inner_only )
          {
            if (i != 1 )
            {
              if (getBit5(qq, 0) != 0) glyph.generateVertex(idx, posL, X_minus, qq);
            }
            if ( i != ix )
            {
              if (getBit5(qq, 1) != 0) glyph.generateVertex(idx, posL, X_plus, qq);
            }
            if ( j != 1 )
            {
              if (getBit5(qq, 2) != 0) glyph.generateVertex(idx, posL, Y_minus, qq);
            }
            if ( j != jx )
            {
              if (getBit5(qq, 3) != 0) glyph.generateVertex(idx, posU, Y_plus, qq);
            }
            if ( k != 1 )
            {
              if (getBit5(qq, 4) != 0) glyph.generateVertex(idx, posU, Z_minus, qq);
            }
            if ( k != kx )
            {
              if (getBit5(qq, 5) != 0) glyph.generateVertex(idx, posU, Z_plus, qq);
            }
          }
          else
          {
            if (getBit5(qq, 0) != 0) glyph.generateVertex(idx, posL, X_minus, qq);
            if (getBit5(qq, 1) != 0) glyph.generateVertex(idx, posL, X_plus , qq);
            if (getBit5(qq, 2) != 0) glyph.generateVertex(idx, posL, Y_minus, qq);
            if (getBit5(qq, 3) != 0) glyph.generateVertex(idx, posU, Y_plus , qq);
            if (getBit5(qq, 4) != 0) glyph.generateVertex(idx, posU, Z_minus, qq);
            if (getBit5(qq, 5) != 0) glyph.generateVertex(idx, posU, Z_plus , qq);
          }
        }

      }
    }
  }


  // ポリゴンの出力
  glyph.writeBinary("CutGlyph");

}
