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
template<typename T>
void FBUtility::xcopy(T* dst, const int* size, const int guide, const T* src, const T scale, const int mode, double& flop)
{
  T s = scale;
  size_t n = (size_t)(size[0]+2*guide) * (size_t)(size[1]+2*guide) * (size_t)(size[2]+2*guide);
  
  if ( mode == kind_vector ) size_t n *= 3;
  
  flop += (double)n;
  
#pragma omp parallel for firstprivate(n, s) schedule(static)
  for (int i=0; i<n; i++) {
    dst[i] = s * src[i];
  }
}


// #################################################################
// 初期化
template<typename T>
void FBUtility::xset(T* dst, const int* size, const int guide, const T init, const int mode)
{
  T s = init;
  size_t n = (size_t)(size[0]+2*guide) * (size_t)(size[1]+2*guide) * (size_t)(size[2]+2*guide);
  
  if ( mode == kind_vector ) size_t n *= 3;
  
#pragma omp parallel for firstprivate(n, s) schedule(static)
  for (int i=0; i<n; i++) {
    dst[i] = s;
  }
}


// #################################################################
// ベクトルの初期化（内部のみ）
template<typename T>
void FBUtility::xsetv(T* dst, const int* size, const int guide, const T* init)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  T[3] s = {init[0], init[1], init[2]};
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd, s) schedule(static)
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        dst[_F_IDX_S3D(i, j, k, ix, jx, kx, gd)] = s;
      }
    }
  }
}
