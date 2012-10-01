// #################################################################
//
// CAERU Library
//
// Copyright (c) 2012 All right reserved.
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
//
// #################################################################

/** 
 * @file   FBUtility.C
 * @brief  FlowBase FBUtility class
 * @author kero
 */


#include "FBUtility.h"


// #################################################################
// ディレクトリがなければ作成、既存なら何もしない
int FBUtility::c_mkdir(string path)
{
  // 標準ライブラリ呼び出し
  // パーミッションはumaskで指定
  umask(022);
  
  int ret = mkdir(path.c_str(), 0777); // rwx
  
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
// 圧力値を有次元から無次元へ変換し，scale倍
void FBUtility::prs_array_D2ND(REAL_TYPE* dst, const int* size, const int guide, const REAL_TYPE Base_prs, const REAL_TYPE Ref_rho, const REAL_TYPE Ref_v, const REAL_TYPE scale, double& flop)
{
  REAL_TYPE dp = scale / (Ref_rho * Ref_v * Ref_v);
  REAL_TYPE bp = Base_prs;

  flop += (double)size[0] * (double)size[1] * (double)size[2] * 2.0 + 10.0;
  
  size_t n = (size_t)(size[0]+2*guide) * (size_t)(size[1]+2*guide) * (size_t)(size[2]+2*guide);
  
#pragma omp parallel for firstprivate(n, dp, bp) schedule(static)
  for (size_t i=0; i<n; i++) {
    dst[i] = ( dst[i] - bp ) * dp;
  }
}


// #################################################################
// 圧力値を無次元から有次元へ変換し，scale倍
void FBUtility::prs_array_ND2D(REAL_TYPE* dst, const int* size, const int guide, const REAL_TYPE* src, const REAL_TYPE Base_prs, const REAL_TYPE Ref_rho, const REAL_TYPE Ref_v, const REAL_TYPE scale, double& flop)
{
  REAL_TYPE dp = scale * (Ref_rho * Ref_v * Ref_v);
  REAL_TYPE bp = Base_prs;
  
  flop += (double)size[0] * (double)size[1] * (double)size[2] * 2.0 + 3.0;
  
  size_t n = (size_t)(size[0]+2*guide) * (size_t)(size[1]+2*guide) * (size_t)(size[2]+2*guide);
  
#pragma omp parallel for firstprivate(n, dp, bp) schedule(static)
  for (size_t i=0; i<n; i++) {
    dst[i] = src[i] * dp + bp;
  }
}


// #################################################################
// 温度値を有次元から無次元へ変換し，scale倍
void FBUtility::tmp_array_D2ND(REAL_TYPE* dst, const int* size, const int guide, const REAL_TYPE Base_tmp, const REAL_TYPE Diff_tmp, const REAL_TYPE klv, const REAL_TYPE scale, double& flop)
{
  REAL_TYPE dp = scale / fabs(Diff_tmp);
  REAL_TYPE bt = klv - Base_tmp;
  
  size_t n = (size_t)(size[0]+2*guide) * (size_t)(size[1]+2*guide) * (size_t)(size[2]+2*guide);
  
  flop += (double)size[0] * (double)size[1] * (double)size[2] * 2.0 + 9.0;
  
#pragma omp parallel for firstprivate(n, dp, bt) schedule(static)
  for (size_t i=0; i<n; i++) {
    dst[i] = ( dst[i] + bt ) * dp;
  }
}


// #################################################################
// 温度値を無次元から有次元へ変換し，scale倍
void FBUtility::tmp_array_ND2D(REAL_TYPE* dst, const int* size, const int guide, const REAL_TYPE* src, const REAL_TYPE Base_tmp, const REAL_TYPE Diff_tmp, const REAL_TYPE klv, const REAL_TYPE scale, double& flop)
{
  REAL_TYPE dp = scale * fabs(Diff_tmp);
  REAL_TYPE bt = Base_tmp - klv;

  size_t n = (size_t)(size[0]+2*guide) * (size_t)(size[1]+2*guide) * (size_t)(size[2]+2*guide);
  
  flop += (double)size[0] * (double)size[1] * (double)size[2] * 2.0 + 3.0;
  
#pragma omp parallel for firstprivate(n, dp, bt) schedule(static)
  for (size_t i=0; i<n; i++) {
    dst[i] = src[i] * dp + bt;
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

  if ( !strcasecmp(mode,"prep") ) {
    fprintf (fp,"\t>> Memory required for Preprocessor : ");
  }
  else if ( !strcasecmp(mode,"solver") ) {
    fprintf (fp,"\t>> Memory required for Solver : ");
  }
  else if ( !strcasecmp(mode,"polygon") ) {
    fprintf (fp,"\t>> Memory required for Polygon : ");
  }
  else if ( !strcasecmp(mode,"cut") ) {
    fprintf (fp,"\t>> Memory required for Cut : ");
  }
  else if ( !strcasecmp(mode,"component") ) {
    fprintf (fp,"\t>> Memory required for Component : ");
  }
  else {
    Exit(0);
  }

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
// スカラー倍してコピー
void FBUtility::xcopy(REAL_TYPE* dst, const int* size, const int guide, const REAL_TYPE* src, const REAL_TYPE scale, const int mode, double& flop)
{
  REAL_TYPE s = scale;
  size_t n = (size_t)(size[0]+2*guide) * (size_t)(size[1]+2*guide) * (size_t)(size[2]+2*guide);
  
  if ( mode == kind_vector ) n *= 3;
  
  flop += (double)n;
  
#pragma omp parallel for firstprivate(n, s) schedule(static)
  for (size_t i=0; i<n; i++) {
    dst[i] = s * src[i];
  }
}



// #################################################################
// 初期化
void FBUtility::xset(REAL_TYPE* dst, const int* size, const int guide, const REAL_TYPE init, const int mode)
{
  REAL_TYPE s = init;
  size_t n = (size_t)(size[0]+2*guide) * (size_t)(size[1]+2*guide) * (size_t)(size[2]+2*guide);
  
  if ( mode == kind_vector ) n *= 3;
  
#pragma omp parallel for firstprivate(n, s) schedule(static)
  for (size_t i=0; i<n; i++) {
    dst[i] = s;
  }
}
