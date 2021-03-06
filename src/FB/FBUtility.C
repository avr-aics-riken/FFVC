//##################################################################################
//
// Flow Base class
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
 * @file   FBUtility.C
 * @brief  FlowBase FBUtility class
 * @author aics
 */


#include "FBUtility.h"


// #################################################################
/**
 * @brief ディレクトリがなければ作成、既存なら何もしない（単一ディレクトリ）
 * @param [in] path ディレクトリパス
 */
int FBUtility::c_mkdir(const char* path)
{
  // 標準ライブラリ呼び出し
  // パーミッションはumaskで指定
  umask(022);
  
  int ret = mkdir(path, 0777); // rwx
  
  if ( 0 != ret )
  {
    // 既存以外のエラー
    if ( EEXIST != errno )
    {
      printf( "\tError(errno)=[%s]\n", strerror(errno) );
      return 0;
    }
  }
  
  return 1;
}


// #################################################################
// ファイル出力時，発散値を計算する
void FBUtility::cnv_Div(REAL_TYPE* dst, REAL_TYPE* src, int* sz, int gc)
{
  copyS3D(dst, sz, gc, src, 1.0);
}


// #################################################################
// 無次元内部エネルギーから有次元/無次元温度への変換
void FBUtility::convArrayIE2Tmp(REAL_TYPE* dst, const int* size, const int guide, const REAL_TYPE* src, const int* bd, const double* mtbl, const REAL_TYPE Base_tmp, const REAL_TYPE Diff_tmp, const int mode, double& flop)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  int zz = mode; // dst[]の次元　1=dimensional, 0=non-dimensional
  
  REAL_TYPE dp = fabs(Diff_tmp);
  REAL_TYPE bt = Base_tmp;
  
  flop += (double)ix * (double)jx * (double)kx * 11.0 + 1.0;
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd, dp, bt, zz) schedule(static) collapse(2)
  for (int k=1-gd; k<=kx+gd; k++) {
    for (int j=1-gd; j<=jx+gd; j++) {
      for (int i=1-gd; i<=ix+gd; i++) {
        
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        int l = bd[m] & MASK_5;
        REAL_TYPE rho = mtbl[3*l+0];
        REAL_TYPE cp  = mtbl[3*l+1];
        REAL_TYPE tn = src[m] / (rho * cp);
        dst[m] = (zz==1) ? (tn * dp + bt) : tn;
      }
    }
  }
  
}


// #################################################################
// 圧力値を有次元から無次元へ変換
void FBUtility::convArrayPrsD2ND(REAL_TYPE* dst, const int* size, const int guide, const REAL_TYPE Base_prs, const REAL_TYPE Ref_rho, const REAL_TYPE Ref_v, double& flop)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  REAL_TYPE dp = 1.0 / (Ref_rho * Ref_v * Ref_v);
  REAL_TYPE bp = Base_prs;

  flop += (double)ix * (double)jx * (double)kx * 2.0 + 10.0;
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd, dp, bp) schedule(static) collapse(2)
  for (int k=1-gd; k<=kx+gd; k++) {
    for (int j=1-gd; j<=jx+gd; j++) {
      for (int i=1-gd; i<=ix+gd; i++) {
        
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        dst[m] = ( dst[m] - bp ) * dp;
      }
    }
  }
}


// #################################################################
// 圧力値を無次元から有次元へ変換
void FBUtility::convArrayPrsND2D(REAL_TYPE* dst, const int* size, const int guide, const REAL_TYPE* src, const REAL_TYPE Base_prs, const REAL_TYPE Ref_rho, const REAL_TYPE Ref_v, double& flop)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  REAL_TYPE dp = Ref_rho * Ref_v * Ref_v;
  REAL_TYPE bp = Base_prs;
  
  flop += (double)ix * (double)jx * (double)kx * 2.0 + 3.0;
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd, dp, bp) schedule(static) collapse(2)
  for (int k=1-gd; k<=kx+gd; k++) {
    for (int j=1-gd; j<=jx+gd; j++) {
      for (int i=1-gd; i<=ix+gd; i++) {
        
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        dst[m] = src[m] * dp + bp;
      }
    }
  }
}


// #################################################################
// 有次元/無次元温度から無次元内部エネルギーへの変換
// mtbl[]は無次元パラメータ
void FBUtility::convArrayTmp2IE(REAL_TYPE* dst, const int* size, const int guide, REAL_TYPE* src, const int* bd, const double* mtbl, const REAL_TYPE Base_tmp, const REAL_TYPE Diff_tmp, const int mode, double& flop)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  int zz = mode; // src[]の次元　1=dimensional, 0=non-dimensional
  
  REAL_TYPE dp = 1.0 / fabs(Diff_tmp);
  REAL_TYPE bt = Base_tmp;
  
  flop += (double)ix * (double)jx * (double)kx * 4.0 + 9.0;
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd, dp, bt, zz) schedule(static) collapse(2)
  for (int k=1-gd; k<=kx+gd; k++) {
    for (int j=1-gd; j<=jx+gd; j++) {
      for (int i=1-gd; i<=ix+gd; i++) {
        
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        REAL_TYPE tn = (zz==1) ? (src[m]-bt)*dp : src[m];
        
        int l = bd[m] & MASK_5;
        REAL_TYPE rho = mtbl[3*l+0];
        REAL_TYPE cp  = mtbl[3*l+1];
        dst[m] = rho * cp * tn;
      }
    }
  }
}


// #################################################################
// 全圧データについて，無次元から有次元単位に変換する
void FBUtility::convArrayTpND2D(REAL_TYPE* src, const int* size, const int guide, const REAL_TYPE Ref_rho, const REAL_TYPE Ref_v)
{
  REAL_TYPE cf = Ref_rho * Ref_v * Ref_v;
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd, cf) schedule(static) collapse(2)
  for (int k=1-gd; k<=kx+gd; k++) {
    for (int j=1-gd; j<=jx+gd; j++) {
      for (int i=1-gd; i<=ix+gd; i++) {
        
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        src[m] *= cf;
      }
    }
  }
}


// #################################################################
// S3D配列のコピー
void FBUtility::copyS3D(REAL_TYPE* dst, const int* size, const int guide, const REAL_TYPE* src, const REAL_TYPE scale)
{
  REAL_TYPE s = scale;
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd, s) schedule(static) collapse(2)
  for (int k=1-gd; k<=kx+gd; k++) {
    for (int j=1-gd; j<=jx+gd; j++) {
      for (int i=1-gd; i<=ix+gd; i++) {
        
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        dst[m] = s * src[m];
      }
    }
  }
}


// #################################################################
// V3D配列のコピー
void FBUtility::copyV3D(REAL_TYPE* dst, const int* size, const int guide, const REAL_TYPE* src, const REAL_TYPE scale)
{
  REAL_TYPE s = scale;
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
    for (int l=0; l<3; l++) {
      
#pragma omp parallel for firstprivate(ix, jx, kx, gd, s, l) schedule(static) collapse(2)
      for (int k=1-gd; k<=kx+gd; k++) {
        for (int j=1-gd; j<=jx+gd; j++) {
          for (int i=1-gd; i<=ix+gd; i++) {
          
          size_t m = _F_IDX_V3D(i, j, k, l, ix, jx, kx, gd);
          dst[m] = s * src[m];
        }
      }
    }
  }
}


// #################################################################
// MediumList中に登録されているkeyに対するIDを返す。発見できない場合はzero
int FBUtility::findIDfromLabel(const MediumList* mat, const int Nmax, const std::string key)
{
  std::string str = key;
  
  for (int i=1; i<=Nmax; i++)
  {
    if ( !strcasecmp(str.c_str(), mat[i].alias.c_str()) ) return i;
  }
  
  return 0;
}


// #################################################################
// CompoList中に登録されているtypeとkeyに対するIDを返す。発見できない場合はzero
int FBUtility::findIDfromCmp(const CompoList* cmp, const int Nmax, const std::string key, const int type)
{
  std::string str = key;
  
  for (int i=1; i<=Nmax; i++)
  {
    if ( !strcasecmp(str.c_str(), cmp[i].medium.c_str()) && cmp[i].getType()==type ) return i;
  }
  
  return 0;
}


// #################################################################
/**
 * @brief S3D配列の初期化 (REAL_TYPE)
 * @param [out]    dst   出力
 * @param [in]     size  配列サイズ
 * @param [in]     guide ガイドセルサイズ
 * @param [in]     init  定数
 */
void FBUtility::initS3D(REAL_TYPE* dst, const int* size, const int guide, const REAL_TYPE init)
{
  REAL_TYPE s = init;
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd, s) schedule(static) collapse(2)
  for (int k=1-gd; k<=kx+gd; k++) {
    for (int j=1-gd; j<=jx+gd; j++) {
      for (int i=1-gd; i<=ix+gd; i++) {

        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        dst[m] = s;
      }
    }
  }
  
}

// #################################################################
/**
 * @brief S3D配列の初期化 (Int)
 * @param [out]    dst   出力
 * @param [in]     size  配列サイズ
 * @param [in]     guide ガイドセルサイズ
 * @param [in]     init  定数
 */
void FBUtility::initS3D(int* dst, const int* size, const int guide, const int init)
{
  int s = init;
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd, s) schedule(static) collapse(2)
  for (int k=1-gd; k<=kx+gd; k++) {
    for (int j=1-gd; j<=jx+gd; j++) {
      for (int i=1-gd; i<=ix+gd; i++) {
        
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        dst[m] = s;
      }
    }
  }
  
}


// #################################################################
// S4DEX配列の初期化
void FBUtility::initS4DEX(REAL_TYPE* dst, const int* size, const int guide, const REAL_TYPE init)
{
  REAL_TYPE s = init;
  size_t nx = (size[0]+2*guide) * (size[1]+2*guide) * (size[2]+2*guide) * 6;
  
#pragma omp parallel for firstprivate(nx, s) schedule(static)
  for (int m=0; m<nx; m++)
  {
    dst[m] = s;
  }
}


// #################################################################
// メモリ使用量を表示する
void FBUtility::MemoryRequirement(const char* mode, const double Memory, const double l_memory, FILE* fp)
{
  const double mem = Memory;
  const double lmem= l_memory;
  const double KB = 1024.0;
  const double MB = 1024.0*KB;
  const double GB = 1024.0*MB;
  const double TB = 1024.0*GB;
  const double PB = 1024.0*TB;
  const double factor = 1.05; // estimate 5% for addtional

  fprintf (fp,"\t>> Memory required for %s : ", mode);

  // Global memory
  fprintf (fp," Global=");
  if ( mem > PB ) {
    fprintf (fp,"%6.2f (PB)", mem / PB *factor);
  }
  else if ( mem > TB ) {
    fprintf (fp,"%6.2f (TB)", mem / TB *factor);
  }
  else if ( mem > GB ) {
    fprintf (fp,"%6.2f (GB)", mem / GB *factor);
  }
  else if ( mem > MB ) {
    fprintf (fp,"%6.2f (MB)", mem / MB *factor);
  }
  else if ( mem > KB ) {
    fprintf (fp,"%6.2f (KB)", mem / KB *factor);
  }
  else if ( mem <= KB ){
    fprintf (fp,"%6.2f (B)", mem *factor);
  }
  else {
    fprintf (fp,"Caution! Memory required : %d (Byte)", (int)(mem *factor) );
  }
  
  // Local memory
  fprintf (fp," : Local=");
  if ( lmem > PB ) {
    fprintf (fp,"%6.2f (PB)\n", lmem / PB *factor);
  }
  else if ( lmem > TB ) {
    fprintf (fp,"%6.2f (TB)\n", lmem / TB *factor);
  }
  else if ( lmem > GB ) {
    fprintf (fp,"%6.2f (GB)\n", lmem / GB *factor);
  }
  else if ( lmem > MB ) {
    fprintf (fp,"%6.2f (MB)\n", lmem / MB *factor);
  }
  else if ( lmem > KB ) {
    fprintf (fp,"%6.2f (KB)\n", lmem / KB *factor);
  }
  else if ( lmem <= KB ){
    fprintf (fp,"%6.2f (B)\n", lmem *factor);
  }
  else {
    fprintf (fp,"Caution! Memory required : %d (Byte)\n", (int)(lmem *factor) );
  }

  fflush(fp);
}


// #################################################################
/**
 * @brief 階層ディレクトリの作成
 * @param [in] path ディレクトリパス
 */
int FBUtility::mkdirs(string path)
{
  int len = path.size() + 4;
  char* buf = new char[len];
  
  if ( !buf )
  {
    printf("Error: create buffer(%d) %s\n", errno, strerror(errno));
    return(-1);
  }
  strcpy(buf, path.c_str());
  
  // 階層的にディレクトリを作成
  char *p = NULL;
  int ret = 0;
  
  for(p=strchr(buf+1, '/'); p; p=strchr(p+1, '/'))
  {
    *p = '\0';
    ret = c_mkdir(buf);
    if (ret != 1)
    {
      delete [] buf;
      return(-1);
    }
    *p = '/';
  }
  
  if (buf)
  {
    delete [] buf;
  }
  return(1);
}
