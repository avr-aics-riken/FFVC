//##################################################################################
//
// FFV-C ASD module : Frontflow / violet Cartesian Active SubDomain
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
 * @file   ASDmodule.C
 * @brief  asdの関数
 * @author aics
 */

#include "ASDmodule.h"



//##################################################################################
// Active SubDomainを作成し，統計情報を表示する
void ASD::evaluateASD(int argc, char **argv)
{
  // 分割数がコマンドラインで陽に与えられている場合，不必要な表示はしない
  bool flag = false;
  
  unsigned int NumDivSubDomain = 0;
  
  if (argc == 6) {
    flag = true;
    G_division[0] = atoi(argv[3]);
    G_division[1] = atoi(argv[4]);
    G_division[2] = atoi(argv[5]);
  }
  else if (argc == 4) {
    flag = true;
    NumDivSubDomain = (unsigned int)atoi(argv[3]);
  }

  // TextPaserをインスタンス
  TextParser tpCntl;
  

  // 利用ライブラリのバージョン番号取得
  string ver_Poly= PL->getVersionInfo();
  string ver_TP  = tpCntl.getVersionInfo();
  
  if ( !flag ) {
    printf("\n\t>> Library Information\n\n");
    printf("\t     Polylib    Version %s\n", ver_Poly.c_str());
    printf("\t     TextParser Version %s\n", ver_TP.c_str());
    printf("\n");
  }
  
  
  // 入力ファイルの指定
  string input_file = argv[2];
  

  
  // パラメータのロードと保持
  {
    int ierror=0;
    
    if ( (ierror = tpCntl.read(input_file)) != TP_NO_ERROR )
    {
      stamped_printf("\tError at reading '%s' file : %d\n", input_file.c_str(), ierror);
      Exit(0);
    }
  }

  
  getDomainInfo(&tpCntl, flag);
  
  // サブドメイン分割数のみ指定された、自動分割の場合
  if (argc == 4)
  {
    unsigned Sdiv[3]={0, 0, 0};
    
    if ( divPolicy == DIV_VOX_CUBE)
    {
      if ( !DecideDivPatternCube(NumDivSubDomain, (unsigned*)size, Sdiv) )
      {
        printf("\tError at calculating division pattern\n");
        Exit(0);
      }
    }
    else // DIV_COMM_SIZE
    {
      if ( !DecideDivPatternCommSize(NumDivSubDomain, (unsigned*)size, Sdiv) )
      {
        printf("\tError at calculating division pattern\n");
        Exit(0);
      }
    }

    G_division[0] = (int)Sdiv[0];
    G_division[1] = (int)Sdiv[1];
    G_division[2] = (int)Sdiv[2];
  }
  
  
  // プロセス分割数のチェック
  if ( (G_division[0]<=0) || (G_division[1]<=0) || (G_division[2]<=0) )
  {
    printf("ERROR : Subdomain (%d, %d, %d)\n", G_division[0], G_division[1], G_division[2] );
    Exit(0);
  }

  
  // 領域情報(serial)
  for (int i=0; i<6; i++)
  {
    nID[i] = -2;
  }

  
  // domain info
  if ( !flag ) {
    printf("\n");
    
    printf("\t(ix, jx, kx)   = %13d %13d %13d\n",
           size[0],
           size[1],
           size[2]);
    
    printf("\t(dx, dy, dz)   = %13.6e %13.6e %13.6e\n",
           pitch[0],
           pitch[1],
           pitch[2]);
    
    printf("\t(ox, oy, oz)   = %13.6e %13.6e %13.6e\n",
           G_origin[0],
           G_origin[1],
           G_origin[2]);
    
    printf("\t(Lx, Ly, Lz)   = %13.6e %13.6e %13.6e\n",
           G_region[0],
           G_region[1],
           G_region[2]);
    
    printf("\t(Dx, Dy, Dz)   = %13d %13d %13d\n",
           G_division[0],
           G_division[1],
           G_division[2]);
    printf("\n");
  }
  
  
  sd_rgn[0] = G_region[0] / (REAL_TYPE)G_division[0];
  sd_rgn[1] = G_region[1] / (REAL_TYPE)G_division[1];
  sd_rgn[2] = G_region[2] / (REAL_TYPE)G_division[2];

  if ( !flag ) {
    printf("\tSubdomain size = %12.4e %12.4e %12.4e\n\n", sd_rgn[0], sd_rgn[1], sd_rgn[2]);
    
    printf("\n----------\n\n");
    printf("\t>> Polylib configuration\n\n");
  }



  
  // 計算モデルの入力ソース情報を取得
  string str, label;
  label = "/DomainInfo/Source";
  if ( !(tpCntl.getInspectedValue(label, str )) )
  {
    stamped_printf("\tParsing error : Invalid char* value in '%s'\n", label.c_str());
    Exit(0);
  } 
  
  PL = MPIPolylib::get_instance();

  
  setupPolygonASD(str, flag);
  
  
  REAL_TYPE *pos_x=NULL; 
  REAL_TYPE *pos_y=NULL; 
  REAL_TYPE *pos_z=NULL;
  pos_x = new REAL_TYPE[G_division[0]];
  pos_y = new REAL_TYPE[G_division[1]];
  pos_z = new REAL_TYPE[G_division[2]];

  createSubdomainTable(pos_x, pos_y, pos_z);
  

  // Active subdomain array
  ActiveSubdomain sd(G_division[0], G_division[1], G_division[2]);
  
  
  
  // cut & memory allocation for d_bcd, d_bid, d_cut
  CalculateCut();
  
  
  Geometry GM;
  
  // fill
  fill(flag, &GM);

  
  // サブドメインにかかるポリゴンがあれば活性にする >> フィルが必要
  int aa = active(pos_x, pos_y, pos_z, sd.get_ptr());
  
  if ( !flag ) printf("\tSubdomain touching polygons = %d\n\n", aa);
  
  
  // Active flag
  int ac=0;
  int dvx = G_division[0]+2*guide;
  int dvy = G_division[1]+2*guide;
  int dvz = G_division[2]+2*guide;
  int ix  = G_division[0];
  int jx  = G_division[1];
  int kx  = G_division[2];
  int gd  = guide;
  int mdf = md_fluid;
  unsigned char* usd = sd.get_ptr();
  
  // from sd to bcd
//#pragma omp parallel for firstprivate(ix, jx, kx, gd, mdf) schedule(static) reduction(+:ac)
#pragma omp parallel for firstprivate(ix, jx, kx, gd, mdf) schedule(dynamic) collapse(3) reduction(+:ac)
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        size_t mc= _F_IDX_S3D(i, j, k, ix, jx, kx, 0);
        if ( usd[mc] == 1 ) {
          d_bcd[m] = mdf;
          ac++;
        }
      }
    }
  }
  
  ac = 0;
  // from bcd to sd
//#pragma omp parallel for firstprivate(ix, jx, kx, gd, mdf) schedule(static) reduction(+:ac)
#pragma omp parallel for firstprivate(ix, jx, kx, gd, mdf) schedule(dynamic) collapse(3) reduction(+:ac)
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        size_t mc= _F_IDX_S3D(i, j, k, ix, jx, kx, 0);
        if ( d_bcd[m] == mdf ) {
          usd[mc] = 1;
          ac++;
        }
        else{
          usd[mc] = 0;
        }
      }
    }
  }
 
  
  if ( !flag ) {
    printf("\tActive SubDomains  = %d / %d\n", ac, G_division[0] * G_division[1] * G_division[2]);
  }
  
  
  
  /*
   label = "/DomainInfo/ActiveSubdomainFile";
   if ( !(tpCntl.getInspectedValue(label, str )) )
   {
   stamped_printf("\tParsing error : Invalid char* value in '%s'\n", label.c_str());
   Exit(0);
   }
   */
  
  
  // subdomain
  if ( !strcasecmp(out_sub.c_str(), "no" ) )
  {
    ;
  }
  else
  {
    if ( !sd.SaveActiveSubdomain(out_sub) ) {
      printf("File write error\n");
      Exit(0);
    }
    printf("\n\tsaved '%s'\n", out_sub.c_str());
  }
  
  
  
  // for V-Xgen Debug
  if ( !strcasecmp(out_svx.c_str(), "no" ) )
  {
    ;
  }
  else
  {
    if ( !sd.writeSVX(out_svx, sd_rgn, G_origin)) {
      printf("SVX file write error\n");
      Exit(0);
    }
    printf("\tsaved '%s'\n\n", out_svx.c_str());
  }
  

  
  
  // 測定モード
  REAL_TYPE x = (REAL_TYPE)size[0] / (REAL_TYPE)G_division[0];
  REAL_TYPE y = (REAL_TYPE)size[1] / (REAL_TYPE)G_division[1];
  REAL_TYPE z = (REAL_TYPE)size[2] / (REAL_TYPE)G_division[2];
  REAL_TYPE sv_ratio = 2.0*(x*y + y*z + z*x) / (x*y*z);
  REAL_TYPE load = x*y*z * (REAL_TYPE)ac;
  REAL_TYPE subload = x*y*z;
  
  int xxx = G_division[0] * G_division[1] * G_division[2];
  printf("Division %d %d %d\n", G_division[0], G_division[1], G_division[2]);
  printf("Pitch %13.6e %13.6e %13.6e\n", pitch[0], pitch[1], pitch[2]);
  printf("SubDomain_Size %d %d %d\n", (int)x, (int)y, (int)z);
  printf("Active/Total %d %d %e\n", ac, xxx, (REAL_TYPE)ac/(REAL_TYPE)xxx);
  printf("Surface/Volume %e\n", sv_ratio);
  printf("N-Active-Workload-SV %d %d %e %e %e %e\n", xxx, ac, (REAL_TYPE)ac/(REAL_TYPE)xxx, load, subload, sv_ratio);
  
  if ( pos_x )  delete [] pos_x;
  if ( pos_y )  delete [] pos_y;
  if ( pos_z )  delete [] pos_z;
  
}



// #################################################################
// active subdomain flag
int ASD::active(const REAL_TYPE* px,
                const REAL_TYPE* py,
                const REAL_TYPE* pz,
                unsigned char* sd)
{
  int dvx = G_division[0];
  int dvy = G_division[1];
  int dvz = G_division[2];

  REAL_TYPE lx = sd_rgn[0];
  REAL_TYPE ly = sd_rgn[1];
  REAL_TYPE lz = sd_rgn[2];

  vector<PolygonGroup*>* pg_roots = PL->get_root_groups();
  vector<PolygonGroup*>::iterator it;
 
  for (it = pg_roots->begin(); it != pg_roots->end(); it++)
  {
    string m_pg = (*it)->get_name();     // グループラベル
    findPolygon(px, py, pz, sd, m_pg); 
  }
 
  delete pg_roots;
  
  
  // count
  int c=0;
  
//#pragma omp parallel for firstprivate(dvx, dvy, dvz) schedule(static) reduction(+:c)
#pragma omp parallel for firstprivate(dvx, dvy, dvz) schedule(dynamic) collapse(3)  reduction(+:c)
  for (int k=0; k<dvz; k++) {
    for (int j=0; j<dvy; j++) {
      for (int i=0; i<dvx; i++) {
        
        size_t m = i + j * dvx + k * dvx * dvy;
        if ( sd[m]>0 ) c++;
      }
    }
  }
  
  return c;
}



// #################################################################
// position of min/max for each subdomain
void ASD::createSubdomainTable(REAL_TYPE* p_x, REAL_TYPE* p_y, REAL_TYPE* p_z)
{
  int div_x = G_division[0];
  int div_y = G_division[1];
  int div_z = G_division[2];
  
  REAL_TYPE ox = G_origin[0];
  REAL_TYPE oy = G_origin[1];
  REAL_TYPE oz = G_origin[2];
  
  REAL_TYPE lx = sd_rgn[0];
  REAL_TYPE ly = sd_rgn[1];
  REAL_TYPE lz = sd_rgn[2];
  
#pragma omp parallel for firstprivate(div_x, lx, ox) schedule(static)
  for (int i=1; i<=div_x; i++) {
    p_x[i-1] = (REAL_TYPE)(i-1) * lx + ox;
  }
  
#pragma omp parallel for firstprivate(div_y, ly, oy) schedule(static)
  for (int j=1; j<=div_y; j++) {
    p_y[j-1] = (REAL_TYPE)(j-1) * ly + oy;
  }
  
#pragma omp parallel for firstprivate(div_z, lz, oz) schedule(static)
  for (int k=1; k<=div_z; k++) {
    p_z[k-1] = (REAL_TYPE)(k-1) * lz + oz;
  }
  
}



// #################################################################
// @brief カット計算
void ASD::CalculateCut()
{
  size_t n_cell[3];
  
  double cut_org[3];
  double cut_dx[3];
  
  cut_dx[0]  = (double)sd_rgn[0];
  cut_dx[1]  = (double)sd_rgn[1];
  cut_dx[2]  = (double)sd_rgn[2];
  cut_org[0] = (double)G_origin[0];
  cut_org[1] = (double)G_origin[1];
  cut_org[2] = (double)G_origin[2];
  
  for ( int i=0; i<3; i++)
  {
    n_cell[i] = (size_t)(G_division[i] + 2*guide);  // サブドメイン分割数+ガイドセル両側
    cut_org[i] -= cut_dx[i]*(double)guide;    // ガイドセルを含む領域のマイナス側頂点の座標
  }
  
  
  // アロケート
  size_t size_n_cell = n_cell[0] * n_cell[1] * n_cell[2];
  d_cut = new long long[size_n_cell*6];
  memset(d_cut, 0, sizeof(float)*size_n_cell*6);
  
  d_bid = new int[size_n_cell];
  memset(d_bid, 0, sizeof(int)*size_n_cell);
  
  d_bcd = new int[size_n_cell];
  memset(d_bcd, 0, sizeof(int)*size_n_cell);
  
  
  
  // 交点計算
  //GM.quantizeCut(d_cut, d_bid, PL, PG);

  
  
  
#if 0
  int ix = G_division[0];
  int jx = G_division[1];
  int kx = G_division[2];
  int gd = guide;
  
  FILE *fp=NULL;
  
  if ( !(fp=fopen("cutinfo.txt","w")) )
  {
    Hostonly_ printf("\tSorry, can't open 'cutinfo.txt', write failed.\n");
    Exit(0);
  }
  
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        int bd = d_bid[mb];
        
        int b0 = getBit5(bd, X_minus); // (bd >> 0)  & MASK_5;
        int b1 = getBit5(bd, X_plus);  // (bd >> 5)  & MASK_5;
        int b2 = getBit5(bd, Y_minus); // (bd >> 10) & MASK_5;
        int b3 = getBit5(bd, Y_plus);  // (bd >> 15) & MASK_5;
        int b4 = getBit5(bd, Z_minus); // (bd >> 20) & MASK_5;
        int b5 = getBit5(bd, Z_plus);  // (bd >> 25) & MASK_5;
        
        //fprintf(fp, "%d %d %d %d %d %d : ", b0, b1, b2, b3, b4, b5);
        
        fprintf(fp, "%5d %5d %5d : %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %d %d %d %d %d %d\n", i,j,k,
                d_cut[m], // x- w
                d_cut[m], // x+ e
                d_cut[m], // y- s
                d_cut[m], // y+ n
                d_cut[m], // z- b
                d_cut[m], // z+ t
                b0, b1, b2, b3, b4, b5);
        
      }
    }
  }
#endif

  
}


// #################################################################
// @brief フィル
void ASD::fill(bool disp_flag, Geometry* GM)
{
  unsigned long target_count; ///< フィルの対象となるセル数
  unsigned long replaced;     ///< 置換された数
  unsigned long filled;       ///< FLUIDでフィルされた数
  unsigned long sum_replaced; ///< 置換された数の合計
  unsigned long sum_filled;   ///< FLUIDでフィルされた数の合計

  
  // 最初にフィル対象のセル数を求める >> 全計算内部セル数
  unsigned long total_cell = (unsigned long)G_division[0] * (unsigned long)G_division[1] * (unsigned long)G_division[2];
  
  
  target_count = total_cell;
  
  if ( !disp_flag )
  {
    printf(    "\tFill initialize -----\n\n");
    printf    ("\t\tInitial target count   = %16ld\n", target_count);
  }
  
  
  // 定義点上に交点がある場合の処理 >> カットするポリゴンのエントリ番号でフィルする
  unsigned long fill_cut = GM->fillCutOnCellCenter(d_bcd, d_bid, d_cut, G_division);
  target_count -= fill_cut;
  
  if ( !disp_flag )
  {
    if ( fill_cut > 0 )
    {
      printf(    "\t\tFill center cut        = %16ld\n", fill_cut);
      printf    ("\t\tTarget count           = %16ld\n\n", target_count);
    }
  }
  
  
  //int FillSuppress[3] = {1, 1, 1};
  
  
  if ( !disp_flag )
  {
    printf(    "\n\tFill -----\n\n");
    printf(    "\t\tFilling Fluid Medium   : Solid\n");
    printf(    "\t\tFill Seed Medium       : Solid\n");
    printf(    "\t\tFill Seed Direction    : %d\n", FillSeedDir);
    //printf(    "\t\tFill Control (X, Y, Z) : (%s, %s, %s)\n\n",
    //       ( !FillSuppress[0] ) ? "Suppress" : "Fill",
    //       ( !FillSuppress[1] ) ? "Suppress" : "Fill",
    //       ( !FillSuppress[2] ) ? "Suppress" : "Fill");
  }
  
  
  filled = GM->fillSeedBcd(d_bcd, FillSeedDir, md_solid, d_bid, G_division);

  
  if ( filled == 0 )
  {
    if ( !disp_flag )
    {
      printf(    "No cells to paint\n");
    }
    Exit(0);
  }
  
  if ( !disp_flag )
  {
    printf(    "\t\tPainted cells          = %16ld\n", filled);
  }
  
  // ペイントされたシードセル数をターゲットから差し引く
  target_count -= filled;
  
  
  
  if ( !disp_flag )
  {
    printf(    "\t\tRemaining target cells = %16ld\n\n", target_count);
  }
  
  
  
  // 隣接する流体セルと接続しており，かつ固体セルに挟まれていないセルのみペイントする
  
  int c=-1; // iteration
  sum_replaced = 0;
  sum_filled = 0;
  
  while (target_count > 0) {
    
    // SeedIDで指定された媒質でフィルする．FLUID/SOLIDの両方のケースがある
    unsigned long fs;
    int cs=1; // @todo 適当なので正しく動かない　修正の必要あり
    filled = GM->fillByBid(d_bid, d_bcd, d_cut, md_solid, fs, cs, G_division);
    replaced = fs;
    
    target_count -= filled;
    target_count -= replaced;
    sum_filled   += filled;
    sum_replaced += replaced;
    
    if ( filled <= 0 ) break; // フィル対象がなくなったら終了
    c++;
  }
  
  
  if ( !disp_flag )
  {
    printf(    "\t\tBID Iteration          = %5d\n", c+1);
    printf(    "\t\t    Filled by Solid    = %16ld\n", sum_filled);
    printf(    "\t\t    SOLID replaced     = %16ld\n", sum_replaced);
    printf(    "\t\t    Remaining cell     = %16ld\n\n", target_count);
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
  
  
  
  // 未ペイント（ID=0）のセルを検出
  unsigned long upc = GM->countCellB(d_bcd, 0, true, G_division);
  
  if ( upc == 0 )
  {
    if ( !disp_flag )
    {
      printf(    "\t\tUnpainted cell         = %16ld\n\n", upc);
    }
  }
  
  
  // 未ペイントのセルに対して、指定媒質でフィルする
  while ( target_count > 0 ) {
    
    replaced = GM->fillByFluid(d_bcd, md_fluid, d_bid, G_division);
    
    target_count -= replaced;
    sum_replaced += replaced;
    
    if ( replaced <= 0 ) break;
    c++;
  }
  
  
  if ( !disp_flag )
  {
    printf(    "\t\tFinal Filling Iteration= %5d\n", c+1);
    printf(    "\t\t   Filled by fluid = %16ld\n\n", sum_replaced);
  }
  
  
  
  // ID=0をカウント
  upc = GM->countCellB(d_bcd, 0, true, G_division);
  
  if ( upc != 0 )
  {
    if ( !disp_flag )
    {
      printf(    "\tFill operation is done, but still remains %ld unpainted cells.\n", upc);
    }
    Exit(0);
  }
  
}



// #################################################################
// サブドメイン内に含まれるポリゴンリストを検索し，フラグを立てる
void ASD::findPolygon(const REAL_TYPE* px,
                      const REAL_TYPE* py,
                      const REAL_TYPE* pz,
                      unsigned char* sd,
                      const string label)
{
  int dvx = G_division[0];
  int dvy = G_division[1];
  int dvz = G_division[2];
  
  REAL_TYPE lx = sd_rgn[0];
  REAL_TYPE ly = sd_rgn[1];
  REAL_TYPE lz = sd_rgn[2];
  
//#pragma omp parallel for firstprivate(dvx, dvy, dvz, lx, ly, lz) schedule(static)
#pragma omp parallel for firstprivate(dvx, dvy, dvz, lx, ly, lz) schedule(dynamic) collapse(3)
  for (int k=1; k<=dvz; k++) {
    for (int j=1; j<=dvy; j++) {
      for (int i=1; i<=dvx; i++) {
        Vec3r pos_min(px[i-1],    py[j-1],    pz[k-1]);
        Vec3r pos_max(px[i-1]+lx, py[j-1]+ly, pz[k-1]+lz);
        
        // false; ポリゴンが一部でもかかる場合
        vector<Triangle*>* trias = PL->search_polygons(label, pos_min, pos_max, false);
        int polys = trias->size();
        
        if (polys>0) {
          size_t m = i-1 + (j-1) * dvx + (k-1) * dvx * dvy;
          sd[m] = 1; //Active
        }
        
        delete trias;
      }
    }
  }
  
}



// #################################################################
// @brief グローバルな領域情報を取得
// @param [in] tpCntl テキストパーサー
// @param [in] flag   表示フラグ
void ASD::getDomainInfo(TextParser* tp, bool flag)
{
  if ( !tp ) Exit(0);
  
  string label, str;
  
  // G_origin　必須
  label = "/DomainInfo/GlobalOrigin";
  
  if ( !tp->getInspectedVector(label, G_origin, 3) )
  {
    cout << "ERROR : in parsing [" << label << "]" << endl;
    Exit(0);
  }
  
  
  // G_region 必須
  label = "/DomainInfo/GlobalRegion";
  
  if ( !tp->getInspectedVector(label, G_region, 3) )
  {
    cout << "ERROR : in parsing [" << label << "]" << endl;
    Exit(0);
  }
  
  if ( (G_region[0]>0.0) && (G_region[1]>0.0) && (G_region[2]>0.0) )
  {
    ; // skip
  }
  else
  {
    cout << "ERROR : in parsing [" << label << "]" << endl;
    Exit(0);
  }
  
  
  // G_voxel
  label = "/DomainInfo/GlobalVoxel";
  
  if ( !tp->getInspectedVector(label, size, 3) ) {
    ;
  }
  
  
  
  if ( (size[0]>0) && (size[1]>0) && (size[2]>0) )
  {
    
    pitch[0] = G_region[0] / size[0];
    pitch[1] = G_region[1] / size[1];
    pitch[2] = G_region[2] / size[2];
    
    // 等方性チェック
    REAL_TYPE p1 = pitch[0] - pitch[1];
    REAL_TYPE p2 = pitch[1] - pitch[2];
    if ( (p1 > SINGLE_EPSILON) || (p2 > SINGLE_EPSILON) )
    {
      printf("\tGlobal Pitch must be same in all direction (%14.6e, %14.6e, %14.6e)\n", pitch[0], pitch[1], pitch[2]);
      Exit(0);
    }
    pitch[2] = pitch[1] = pitch[0];
    
  }
  else
  {
    printf("ERROR : in parsing [%s] >> (%d, %d, %d)\n", label.c_str(), size[0], size[1], size[2] );
    Exit(0);
  }
  
  
  
  // G_division オプション >>  コマンドラインで指定された場合はスキップ
  if ( !flag ) {
    label = "/DomainInfo/GlobalDivision";
    
    if ( !tp->getInspectedVector(label, G_division, 3) )
    {
      cout << "\tGlobalDivision is required." << endl;
    }
  }
  
  
  // 領域分割ポリシー
  label = "/DomainInfo/DivisionPolicy";
  
  if ( !(tp->getInspectedValue(label, str)) )
  {
    Hostonly_ stamped_printf("\tParsing error : Invalid value in '%s'\n", label.c_str());
    Exit(0);
  }
  else
  {
    if     ( !strcasecmp(str.c_str(), "cube" ) )           divPolicy = DIV_VOX_CUBE;
    else if( !strcasecmp(str.c_str(), "communication" ) )  divPolicy = DIV_COMM_SIZE;
    else
    {
      Hostonly_ printf("\tInvalid string '%s'\n", str.c_str());
      Exit(0);
    }
  }
  
  
  
  // 流体セルのフィルの開始面指定
  label = "/DomainInfo/HintOfFillSeedDirection";
  
  if ( !(tp->getInspectedValue(label, str)) )
  {
    Hostonly_ stamped_printf("\tParsing error : Invalid value in '%s'\n", label.c_str());
    Exit(0);
  }
  else
  {
    if     ( !strcasecmp(str.c_str(), "xminus" ) ) FillSeedDir = 0;
    else if( !strcasecmp(str.c_str(), "xplus" ) )  FillSeedDir = 1;
    else if( !strcasecmp(str.c_str(), "yminus" ) ) FillSeedDir = 2;
    else if( !strcasecmp(str.c_str(), "yplus" ) )  FillSeedDir = 3;
    else if( !strcasecmp(str.c_str(), "zminus" ) ) FillSeedDir = 4;
    else if( !strcasecmp(str.c_str(), "zplus" ) )  FillSeedDir = 5;
    else
    {
      FillSeedDir = 0;
      Hostonly_ printf("\tDefault 'X_minus' is set for Hint Of FillSeed direction\n");
    }
  }
  
  // output
  label = "/DomainInfo/outputSVX";
  
  if ( !(tp->getInspectedValue(label, str)) )
  {
    cout << "ERROR : in parsing [" << label << "]" << endl;
    Exit(0);
  }
  out_svx = str;
  
  
  label = "/DomainInfo/outputSubdomain";
  
  if ( !(tp->getInspectedValue(label, str)) )
  {
    cout << "ERROR : in parsing [" << label << "]" << endl;
    Exit(0);
  }
  out_sub = str;
  
}



// #################################################################
// @brief 幾何形状情報のロード
// @param [in] fname Polylibファイル名
// @param [in] flag  表示フラグ
void ASD::setupPolygonASD(const string fname, bool flag)
{
  
  POLYLIB_STAT poly_stat;  ///< Polylibの戻り値
  
  unsigned poly_gc[3] = {0, 0, 0};
  REAL_TYPE poly_org[3];
  REAL_TYPE poly_dx[3];
  
  poly_dx[0]  = pitch[0];
  poly_dx[1]  = pitch[1];
  poly_dx[2]  = pitch[2];
  poly_org[0] = G_origin[0];
  poly_org[1] = G_origin[1];
  poly_org[2] = G_origin[2];
   
  poly_stat = PL->init_parallel_info(MPI_COMM_WORLD,
                                     poly_org,         // 自ランクの基点座標
                                     (unsigned*)size,  // 自ランクの分割数
                                     poly_gc,          // ガイドセル数
                                     poly_dx           // 格子幅
                                     );
   

  poly_stat = PL->load_rank0( fname );
 
  if ( poly_stat != PLSTAT_OK ) {
    printf ("\tpolylib->load_rank0() failed.");
    Exit(0);
  }
 
 
  // 階層情報表示 debug brief hierarchy
// ##########
#if 0
  PL->show_group_hierarchy();
#endif
// ##########
 
 
  // ポリゴン情報へのアクセス
  vector<PolygonGroup*>* pg_roots = PL->get_root_groups();
 
 
  // Polygon Groupの数
  int num_of_polygrp = pg_roots->size();
 
   if ( !flag ) {
     printf(     "\n\tNumber of Polygon Group = %d\n\n", num_of_polygrp);
     printf(     "\t   Polygon Group Label       Medium Alias              Local BC      Element          Area\n");
     printf(     "\t   ---------------------------------------------------------------------------------------\n");
   }

 
 
  vector<PolygonGroup*>::iterator it;
  const int mat_id = 2;
 
  // ポリゴン情報の表示
  for (it = pg_roots->begin(); it != pg_roots->end(); it++)
  {
    string m_pg = (*it)->get_name();     // グループラベル
    string m_mat = (*it)->get_label();   // 媒質ラベル
    string m_bc = (*it)->get_type();     // 境界条件ラベル
    int ntria= (*it)->get_group_num_tria();   // ローカルのポリゴン数
    REAL_TYPE area = (*it)->get_group_area(); // ローカルのポリゴン面積
    
    // PolygonにIDを割り当てる
    poly_stat = (*it)->set_all_exid_of_trias(mat_id);
    
    if ( poly_stat != PLSTAT_OK )
    {
      if ( !flag )
      {
        printf(     "\tError : Polylib::set_all_exid_of_trias()\n");
        Exit(0);
      }
    }
    
    if ( !flag ) {
      printf("\t  %20s %18s  %20s %12d  %e\n",
             m_pg.c_str(),
             m_mat.c_str(),
             m_bc.c_str(),
             ntria,
             area);
    }
 
 // ########## show corrdinates and area
 #if 0
    PL->show_group_info(m_pg); //debug
 #endif
 // ##########
 
  }
 
  delete pg_roots;
 
  if ( !flag ) printf("\n");
 
 }


// #################################################################
// @brief 通信コストの計算
// @param [in] iDiv    サブドメイン分割数（x方向）
// @param [in] jDiv    サブドメイン分割数（y方向）
// @param [in] kDiv    サブドメイン分割数（z方向）
// @param [in] voxSize 全領域のボクセル分割数
// @ret 全領域の通信面の和
inline
unsigned long long
ASD::CalcCommSize(const unsigned long long iDiv,
                  const unsigned long long jDiv,
                  const unsigned long long kDiv,
                  const unsigned long long voxSize[3])
{
  if( (iDiv==0) || (jDiv==0) || (kDiv==0) ) return 0;
  if( !voxSize ) return 0;
  
  unsigned long long Len[3];
  Len[0] = voxSize[0] / iDiv; if( Len[0] == 0 ) return 0;
  Len[1] = voxSize[1] / jDiv; if( Len[1] == 0 ) return 0;
  Len[2] = voxSize[2] / kDiv; if( Len[2] == 0 ) return 0;
  
  unsigned long long commFace[3];
  if( iDiv != 1) commFace[0] = Len[1]*Len[2]*(iDiv-1);
  else commFace[0] = 0;
  
  if( jDiv != 1) commFace[1] = Len[2]*Len[0]*(jDiv-1);
  else commFace[1] = 0;
  
  if( kDiv != 1) commFace[2] = Len[0]*Len[1]*(kDiv-1);
  else commFace[2] = 0;
  
  return (commFace[0] + commFace[1] + commFace[2]);
}


// #################################################################
/** @brief 並列プロセス数からI,J,K方向の分割数を取得する
 *  @note 通信面のトータルサイズが小さい分割パターンを採用する
 *  @param [in]  divNum  ランク数
 *  @param [in]  voxSize 空間全体のボクセル数
 *  @param [out] divPttn 領域分割数
 *  @return              終了コード(CPM_SUCCESS=正常終了)
 */
bool ASD::DecideDivPatternCommSize(const unsigned int divNum,
                                   const unsigned int voxSize[3],
                                   unsigned int divPttn[3])
{
  if ( !voxSize || !divPttn )
  {
    return false;
  }
  
  if( (voxSize[0]==0) || (voxSize[1]==0) || (voxSize[2]==0) )
  {
    return false;
  }
  
  if ( divNum <= 1 )
  {
    divPttn[0] = divPttn[1] = divPttn[2] = 1;
    return true;
  }
  
  divPttn[0] = divPttn[1] = divPttn[2] = 0;
  
  unsigned long long minCommSize = 0;
  unsigned long long divNumll = divNum;
  unsigned long long voxSizell[3] = {0, 0, 0};
  unsigned long long divPttnll[3] = {0, 0, 0};
  voxSizell[0] = voxSize[0];
  voxSizell[1] = voxSize[1];
  voxSizell[2] = voxSize[2];
  
  bool flag = false;
  unsigned long long i, j, k;
  for(i=1; i<=divNumll; i++)
  {
    if( divNumll%i != 0 ) continue;
    if( voxSizell[0] < i ) break;
    
    unsigned long long jmax = divNumll/i;
    
    for (j=1; j<=jmax; j++)
    {
      if ( jmax%j != 0 ) continue;
      if( voxSizell[1] < j ) break;
      
      k = jmax/j;
      if( voxSizell[2] < k ) continue;
      
      
      unsigned long long commSize;
      if ( (commSize=CalcCommSize(i, j, k, voxSizell)) == 0 ) break;
      
      if ( !flag )
      {
        divPttnll[0] = i; divPttnll[1] = j; divPttnll[2] = k;
        minCommSize = commSize;
        flag = true;
      }
      else if( commSize < minCommSize )
      {
        divPttnll[0] = i; divPttnll[1] = j; divPttnll[2] = k;
        minCommSize = commSize;
      }
    }
  }
  
  divPttn[0] = (unsigned)divPttnll[0];
  divPttn[1] = (unsigned)divPttnll[1];
  divPttn[2] = (unsigned)divPttnll[2];
  
  if ( (divPttn[0]==0) || (divPttn[1]==0) || (divPttn[2]==0) )
  {
    return false;
  }
  
  return true;
}



// #################################################################
/** @brief 並列プロセス数からI,J,K方向の分割数を取得する
 *  @note １つのサブドメインが立方体に一番近い分割パターンを採用する
 *  @param [in]  divNum  ランク数
 *  @param [in]  voxSize 空間全体のボクセル数
 *  @param [out] divPttn 領域分割数
 *  @return              終了コード(CPM_SUCCESS=正常終了)
 */
bool ASD::DecideDivPatternCube(const unsigned int divNum,
                               const unsigned int voxSize[3],
                               unsigned int divPttn[3])
{
  if ( !voxSize || !divPttn )
  {
    return false;
  }
  
  if ( (voxSize[0]==0) || (voxSize[1]==0) || (voxSize[2]==0) )
  {
    return false;
  }
  
  if ( divNum <= 1 )
  {
    divPttn[0] = divPttn[1] = divPttn[2] = 1;
    return true;
  }
  
  divPttn[0] = divPttn[1] = divPttn[2] = 0;
  
  unsigned long long minVoxDiff = 0;
  unsigned long long divNumll = divNum;
  unsigned long long voxSizell[3] = {0, 0, 0};
  unsigned long long divPttnll[3] = {0, 0, 0};
  voxSizell[0] = voxSize[0];
  voxSizell[1] = voxSize[1];
  voxSizell[2] = voxSize[2];
  
  bool flag = false;
  unsigned long long i, j, k;
  for(i=1; i<=divNumll; i++)
  {
    if( divNumll%i != 0 ) continue;
    if( voxSizell[0] < i ) break;
    
    unsigned long long jmax = divNumll/i;
    
    for(j=1; j<=jmax; j++)
    {
      if( jmax%j != 0 ) continue;
      if( voxSizell[1] < j ) break;
      
      k = jmax/j;
      if( voxSizell[2] < k ) continue;
      
      long long voxDiff;
      if( (voxDiff=CheckCube(i, j, k, voxSizell)) < 0 ) break;
      
      if( !flag )
      {
        divPttnll[0] = i; divPttnll[1] = j; divPttnll[2] = k;
        minVoxDiff = voxDiff;
        flag = true;
      }
      else if( voxDiff < minVoxDiff )
      {
        divPttnll[0] = i; divPttnll[1] = j; divPttnll[2] = k;
        minVoxDiff = voxDiff;
      }
    }
  }
  
  divPttn[0] = (unsigned)divPttnll[0];
  divPttn[1] = (unsigned)divPttnll[1];
  divPttn[2] = (unsigned)divPttnll[2];
  
  if ( (divPttn[0]==0) || (divPttn[1]==0) || (divPttn[2]==0) )
  {
    return false;
  }
  
  return true;
}

// #################################################################
/** I,J,K分割を行った時のI,J,Kボクセル数の最大/最小の差を取得する
 *  @param [in] iDiv    i方向領域分割数
 *  @param [in] jDiv    j方向領域分割数
 *  @param [in] kDiv    k方向領域分割数
 *  @param [in] voxSize 空間全体のボクセル数
 *  @retval 0以上        I,J,Kボクセル数の最大/最小の差
 *  @retval 負値         領域分割不可のパターン
 */
long long ASD::CheckCube(const unsigned long long iDiv,
                         const unsigned long long jDiv,
                         const unsigned long long kDiv,
                         const unsigned long long voxSize[3])
{
  if ( (iDiv==0) || (jDiv==0) || (kDiv==0) ) return -1;
  if ( !voxSize ) return -1;
  
  unsigned long long Len[3];
  Len[0] = voxSize[0] / iDiv; if( Len[0] == 0 ) return -1;
  Len[1] = voxSize[1] / jDiv; if( Len[1] == 0 ) return -1;
  Len[2] = voxSize[2] / kDiv; if( Len[2] == 0 ) return -1;
  
  unsigned long long minVox = (Len[0]<Len[1]?Len[0]:Len[1]);
  minVox = (minVox<Len[2]?minVox:Len[2]);
  
  unsigned long long maxVox = (Len[0]>Len[1]?Len[0]:Len[1]);
  maxVox = (maxVox>Len[2]?maxVox:Len[2]);
  
  return (maxVox-minVox);
}
