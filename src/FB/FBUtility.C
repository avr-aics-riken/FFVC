/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file FBUtility.C
//@brief FlowBase FBUtility class Header
//@author keno, FSI Team, VCAD, RIKEN

#include "FBUtility.h"


/**
 @fn string FBUtility::getDirection(const unsigned dir)
 @brief dirの方向ラベルを返す
 @param dir 方向
 */
string FBUtility::getDirection(const unsigned dir)
{
  string face;
  if      (dir == X_MINUS) face = "X-";
  else if (dir == X_PLUS)  face = "X+";
  else if (dir == Y_MINUS) face = "Y-";
  else if (dir == Y_PLUS)  face = "Y+";
  else if (dir == Z_MINUS) face = "Z-";
  else if (dir == Z_PLUS)  face = "Z+";
  return face;
}

/**
 @fn void FBUtility::displayMemory(const char* mode, const unsigned long Memory, const unsigned long l_memory, FILE* fp, FILE* mp)
 @brief メモリ使用量を表示する
 @param mode 処理モード
 @param Memory 必要メモリ量
 @param l_memory local
 */
void FBUtility::displayMemory(const char* mode, const unsigned long Memory, const unsigned long l_memory, FILE* fp, FILE* mp)
{
  MemoryRequirement(mode, Memory, l_memory, mp);
  MemoryRequirement(mode, Memory, l_memory, fp);
}

/**
 @fn void FBUtility::MemoryRequirement(const char* mode, const unsigned long Memory, const unsigned long l_memory, FILE* fp)
 @brief メモリ使用量を表示する
 @param mode 処理モード
 @param Memory 必要メモリ量
 @param l_memory local
 */
void FBUtility::MemoryRequirement(const char* mode, const unsigned long Memory, const unsigned long l_memory, FILE* fp)
{
  const SKL_REAL mem = (SKL_REAL)Memory;
  const SKL_REAL lmem= (SKL_REAL)l_memory;
  const SKL_REAL KB = 1024.0;
  const SKL_REAL MB = 1024.0*KB;
  const SKL_REAL GB = 1024.0*MB;
  const SKL_REAL TB = 1024.0*GB;
  const SKL_REAL PB = 1024.0*TB;
  const SKL_REAL factor = 1.05; // estimate 5% for addtional

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
  else {
    assert(0);
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
    fprintf (fp,"Caution! Memory required : %d (Byte)", mem *factor);
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
    fprintf (fp,"Caution! Memory required : %d (Byte)\n", lmem *factor);
  }

  fflush(fp);
}

/**
 @fn void FBUtility::printVersion(FILE* fp, const char* str, const unsigned ver)
 @brief バージョン情報の表示
 */
void FBUtility::printVersion(FILE* fp, const char* str, const unsigned ver)
{
  unsigned a, b, c;
  a = b = c = 0;

  a = ver / 100;
  b = (ver - a*100) / 10;
  c = ver - a*100 - b*10;

  fprintf(fp,"\n");
  fprintf(fp,"\t%s \tVersion %d.%d.%d\n", str, a, b, c);
}
