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
 * @file   FileIO_SPH_sph_read.C
 * @brief  FileIO_SPH Class
 * @author aics
 */

#include "FileIO_sph.h"


/**
 * ファイルヘッダー情報を取得する。(unsigned型)
 * @param[in] fname     ファイル名
 * @param[out] real_type   データ型タイプ
 * @param[out] data_type   データ種別
 * @param[out] voxsize     ボクセルサイズ
 * @return      true=success, false=fault
 */
bool FileIO_SPH::GetInfo(const char* fname,
                         RealType* real_type,
                         DataType* data_type,
                         unsigned int voxsize[3])
{
  if( !fname || !real_type || !data_type || !voxsize ) return false;

  EMatchType eType;
  eType = isMatchEndian(fname, 8);
  if( eType == UnKnown ){ return false; }

  FILE* fp;
  if( !(fp = fopen(fname, "rb")) ) {
    fprintf(stderr, "Can't open file.(%s)\n", fname);
    return false;
  }

  int dmy;
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(dmy);
  if( dmy != 8 ) { fclose(fp); return false; }
  if( fread(data_type, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(*data_type);
  if( fread(real_type, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(*real_type);
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(dmy);
  if( dmy != 8 ) { fclose(fp); return false; }

  if( *real_type == _FLOAT ){
    if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
    if( eType == UnMatch ) BSWAP32(dmy);
    if( dmy != 12 ) { fclose(fp); return false; }
    if( fread(voxsize, sizeof(int), 3, fp) != 3 ) { fclose(fp); return false; }
    if( eType == UnMatch ){
      BSWAP32(voxsize[0]);
      BSWAP32(voxsize[1]);
      BSWAP32(voxsize[2]);
    }
    if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
    if( eType == UnMatch ) BSWAP32(dmy);
    if( dmy != 12 ) { fclose(fp); return false; }
  }else if( *real_type == _DOUBLE ){
    if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
    if( eType == UnMatch ) BSWAP32(dmy);
    if( dmy != 24 ) { fclose(fp); return false; }
    unsigned long long tmp[3];
    if( fread(tmp, sizeof(long long), 3, fp) != 3 ){ fclose(fp); return false; }
    if( eType == UnMatch ){BSWAP64(tmp[0]); BSWAP64(tmp[1]); BSWAP64(tmp[2]);}
    if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
    if( eType == UnMatch ) BSWAP32(dmy);
    if( dmy != 24 ) { fclose(fp); return false; }
    voxsize[0] = (unsigned int)tmp[0];
    voxsize[1] = (unsigned int)tmp[1];
    voxsize[2] = (unsigned int)tmp[2];
    fprintf(stderr, "WARNNING : The data described in this SPH file(\"%s\") is double precision.", fname);
  }else{
    fclose(fp);
    return false;
  }

  fclose(fp);
  return true;
}

/**
 * ファイルヘッダー情報を取得する。(long long型)
 * @param[in] fname     ファイル名
 * @param[out] real_type   データ型タイプ
 * @param[out] data_type   データ種別
 * @param[out] voxsize     ボクセルサイズ
 * @return      true=success, false=fault
 */
bool FileIO_SPH::GetInfo(const char* fname,
                         RealType* real_type,
                         DataType* data_type,
                         unsigned long long voxsize[3])
{
  if( !fname || !real_type || !data_type || !voxsize ) return false;

  EMatchType eType;
  eType = isMatchEndian(fname, 8);
  if( eType == UnKnown ){ return false; }

  FILE* fp;
  if( !(fp = fopen(fname, "rb")) ) {
    fprintf(stderr, "Can't open file.(%s)\n", fname);
    return false;
  }

  int dmy;
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(dmy);
  if( dmy != 8 ) { fclose(fp); return false; }
  if( fread(data_type, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(*data_type);
  if( fread(real_type, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(*real_type);
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(dmy);
  if( dmy != 8 ) { fclose(fp); return false; }

  if( *real_type == _FLOAT ){
    if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
    if( eType == UnMatch ) BSWAP32(dmy);
    if( dmy != 12 ) { fclose(fp); return false; }
    unsigned int tmp[3];
    if( fread(tmp, sizeof(int), 3, fp) != 3 ) { fclose(fp); return false; }
    if( eType == UnMatch ){
      BSWAP32(tmp[0]);
      BSWAP32(tmp[1]);
      BSWAP32(tmp[2]);
    }
    if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
    if( eType == UnMatch ) BSWAP32(dmy);
    if( dmy != 12 ) { fclose(fp); return false; }
    voxsize[0] = (unsigned long long)tmp[0];
    voxsize[1] = (unsigned long long)tmp[1];
    voxsize[2] = (unsigned long long)tmp[2];
    fprintf(stderr, "WARNNING : The data described in this SPH file(\"%s\") is single precision.", fname);
  }else if( *real_type == _DOUBLE ){
    if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
    if( eType == UnMatch ) BSWAP32(dmy);
    if( dmy != 24 ) { fclose(fp); return false; }
    if( fread(voxsize, sizeof(long long), 3, fp) != 3 ){ fclose(fp); return false; }
    if( eType == UnMatch ){
      BSWAP64(voxsize[0]);
      BSWAP64(voxsize[1]);
      BSWAP64(voxsize[2]);
    }
    if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
    if( eType == UnMatch ) BSWAP32(dmy);
    if( dmy != 24 ) { fclose(fp); return false; }
  }else{
    fclose(fp);
    return false;
  }

  if( *real_type != _FLOAT ){
    fclose(fp);
    return false;
  }

  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(dmy);
  if( dmy != 24 ) { fclose(fp); return false; }
  if( fread(voxsize, sizeof(long long), 3, fp) != 3 ){fclose(fp);return false;}
  if( eType == UnMatch ){
    BSWAP64(voxsize[0]);
    BSWAP64(voxsize[1]);
    BSWAP64(voxsize[2]);
  }
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(dmy);
  if( dmy != 24 ) { fclose(fp); return false; }

  fclose(fp);
  return true;
}


/**
 * ファイルからSPHヘッダデータを読み込む。(出力の型に合わせて読込む）
 * @param[in] fname     ファイル名
 * @param[out] voxorg      原点座標
 * @param[out] voxpit      ピッチ
 * @return      true=success, false=fault
 */
bool FileIO_SPH::GetHeader(const char* fname,
                           double voxorg[3],
                           double voxpit[3])
{
  
  if( !fname || !voxorg || !voxpit ) return false;
  
  EMatchType eType;
  eType = isMatchEndian(fname, 8);
  if( eType == UnKnown ){ return false; }
  
  FILE* fp;
  if( !(fp = fopen(fname, "rb")) ) {
    fprintf(stderr, "Can't open file.(%s)\n", fname);
    return false;
  }
  
  int dmy;
  DataType data_type;
  RealType real_type;
  
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(dmy);
  if( dmy != 8 ) { fclose(fp); return false; }
  if( fread(&data_type, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(data_type);
  if( fread(&real_type, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(real_type);
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(dmy);
  if( dmy != 8 ) { fclose(fp); return false; }
  
  unsigned int voxsize[3];
  int dsize;
  if( real_type == _FLOAT ) {
    dsize = 12;
  } else if( real_type == _DOUBLE ) {
    dsize = 24;
  } else {
    fclose(fp);
    return false;
  }
  
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(dmy);
  if( dmy != dsize ) { fclose(fp); return false; }
  if( real_type == _FLOAT ) {
    if( fread(voxsize, sizeof(int), 3, fp) != 3 ) { fclose(fp); return false; }
    if( eType == UnMatch ){BSWAP32(voxsize[0]); BSWAP32(voxsize[1]); BSWAP32(voxsize[2]);}
  } else {
    unsigned long long tmp[3];
    if( fread(tmp, sizeof(long long), 3, fp) != 3 ){ fclose(fp); return false; }
    if( eType == UnMatch ){BSWAP64(tmp[0]); BSWAP64(tmp[1]); BSWAP64(tmp[2]);}
    voxsize[0] = (unsigned int)tmp[0];
    voxsize[1] = (unsigned int)tmp[1];
    voxsize[2] = (unsigned int)tmp[2];
  }
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(dmy);
  if( dmy != dsize ) { fclose(fp); return false; }
  
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(dmy);
  if( dmy != dsize ) { fclose(fp); return false; }
  if( real_type == _FLOAT ) {
    float tmp[3];
    if( fread(tmp, sizeof(float), 3, fp) != 3 ) { fclose(fp); return false; }
    if( eType == UnMatch ){BSWAP32(tmp[0]); BSWAP32(tmp[1]); BSWAP32(tmp[2]);}
    voxorg[0] = (double)tmp[0];
    voxorg[1] = (double)tmp[1];
    voxorg[2] = (double)tmp[2];
  } else {
    if( fread(voxorg, sizeof(double), 3, fp) != 3 ) { fclose(fp); return false; }
    if( eType == UnMatch ){BSWAP64(voxorg[0]); BSWAP64(voxorg[1]); BSWAP64(voxorg[2]);}
  }
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(dmy);
  if( dmy != dsize ) { fclose(fp); return false; }
  
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(dmy);
  if( dmy != dsize ) { fclose(fp); return false; }
  if( real_type == _FLOAT ) {
    float tmp[3];
    if( fread(tmp, sizeof(float), 3, fp) != 3 ) { fclose(fp); return false; }
    if( eType == UnMatch ){BSWAP32(tmp[0]); BSWAP32(tmp[1]); BSWAP32(tmp[2]);}
    voxpit[0] = (double)tmp[0];
    voxpit[1] = (double)tmp[1];
    voxpit[2] = (double)tmp[2];
  } else {
    if( fread(voxpit, sizeof(double), 3, fp) != 3 ) { fclose(fp); return false; }
    if( eType == UnMatch ){BSWAP64(voxpit[0]); BSWAP64(voxpit[1]); BSWAP64(voxpit[2]);}
  }
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(dmy);
  if( dmy != dsize ) { fclose(fp); return false; }
  
  if( real_type == _FLOAT ) {
    dsize = 8;
  } else {
    dsize = 16;
  }
  
  unsigned int step;
  double time;
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(dmy);
  if( dmy != dsize ) { fclose(fp); return false; }
  if( real_type == _FLOAT ) {
    if( fread(&step, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
    if( eType == UnMatch ){BSWAP32(step);}
    float tmp;
    if( fread(&tmp, sizeof(float), 1, fp) != 1 ) { fclose(fp); return false; }
    if( eType == UnMatch ){BSWAP32(tmp);}
    time = (double)tmp;
  } else {
    unsigned long long tmp;
    if( fread(&tmp, sizeof(long long), 1, fp) != 1 ) { fclose(fp); return false; }
    if( eType == UnMatch ){BSWAP64(tmp);}
    step = (unsigned int)tmp;
    if( fread(&time, sizeof(double), 1, fp) != 1 ) { fclose(fp); return false; }
    if( eType == UnMatch ){BSWAP64(time);}
  }
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(dmy);
  if( dmy != dsize ) { fclose(fp); return false; }
  
  fclose(fp);
  return true;
}



/**
 * ファイルからSPHデータを読み込む。(float型)
 * @param[in] fname     ファイル名
 * @param[in] dsize     データサイズ
 * @param[out] data_type   データ種別
 * @param[out] voxsize     ボクセルサイズ
 * @param[out] voxorg      原点座標
 * @param[out] voxpit      ピッチ
 * @param[out] step      ステップ数
 * @param[out] time      時間
 * @param[out] data     SPHデータ
 * @return      true=success, false=fault
 */
bool FileIO_SPH::GetData(const char* fname,
                         unsigned int dsize,
                         int* data_type,
                         unsigned int voxsize[3],
                         float voxorg[3],
                         float voxpit[3],
                         int* step,
                         float* time,
                         float* data)
{
  if( !fname || !voxsize || !voxorg || !voxpit
             || !data_type || !step || !time || !data ) return false;

  DataType dt;
  RealType real_type;
  if( !GetInfo(fname, &real_type, &dt, voxsize) ) return false;
  if( real_type != _FLOAT ) return false;

  EMatchType eType;
  eType = isMatchEndian(fname, 8);
  if( eType == UnKnown ){ return false; }

  FILE* fp;
  if( !(fp = fopen(fname, "rb")) ) {
    fprintf(stderr, "Can't open file.(%s)\n", fname);
    return false;
  }

  unsigned int dmy;
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(dmy);
  if( dmy != 8 ) { fclose(fp); return false; }
  if( fread(data_type, sizeof(int), 1, fp) != 1 ) {fclose(fp);return false;}
  if( eType == UnMatch ) BSWAP32(*data_type);
  if( fread(&real_type, sizeof(int), 1, fp) != 1 ) {fclose(fp);return false;}
  if( eType == UnMatch ) BSWAP32(real_type);
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(dmy);
  if( dmy != 8 ) { fclose(fp); return false; }
  if( real_type != _FLOAT ){ fclose(fp); return false; }

  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(dmy);
  if( dmy != 12 ) { fclose(fp); return false; }
  if( fread(voxsize, sizeof(int), 3, fp) != 3 ) { fclose(fp); return false; }
  if( eType == UnMatch ){
    BSWAP32(voxsize[0]);
    BSWAP32(voxsize[1]);
    BSWAP32(voxsize[2]);
  }
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(dmy);
  if( dmy != 12 ) { fclose(fp); return false; }

  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(dmy);
  if( dmy != 12 ) { fclose(fp); return false; }
  if( fread(voxorg, sizeof(float), 3, fp) != 3 ) { fclose(fp); return false; }
  if( eType == UnMatch ){
    BSWAP32(voxorg[0]);
    BSWAP32(voxorg[1]);
    BSWAP32(voxorg[2]);
  }
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(dmy);
  if( dmy != 12 ) { fclose(fp); return false; }

  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(dmy);
  if( dmy != 12 ) { fclose(fp); return false; }
  if( fread(voxpit, sizeof(float), 3, fp) != 3 ) { fclose(fp); return false; }
  if( eType == UnMatch ){
    BSWAP32(voxpit[0]);
    BSWAP32(voxpit[1]);
    BSWAP32(voxpit[2]);
  }
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(dmy);
  if( dmy != 12 ) { fclose(fp); return false; }

  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(dmy);
  if( dmy != 8 ) { fclose(fp); return false; }
  if( fread(step, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(*step);
  if( fread(time, sizeof(float), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(*time);
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(dmy);
  if( dmy != 8 ) { fclose(fp); return false; }

  unsigned int dim;
  if( *data_type == _SCALAR ) dim = 1;
  else if( *data_type == _VECTOR ) dim = 3;
  else { fclose(fp); return false; }
  unsigned int aryLen = voxsize[0]*voxsize[1]*voxsize[2]*dim;
  if( aryLen < 1 ) { fclose(fp); return false; }
  if(dsize < aryLen) { fclose(fp); return false; }
  memset(data, 0, sizeof(float)*aryLen);

  if( fread(&dmy, sizeof(int), 1, fp) != 1 ){
    fclose(fp); /*delete [] *data; *data = NULL;*/ return false;
  }
  if( eType == UnMatch ) BSWAP32(dmy);
  if( dmy != aryLen*sizeof(float) ){
    fclose(fp); /*delete [] *data; *data = NULL;*/ return false;
  }
  if( fread(data, sizeof(float), aryLen, fp) != aryLen ){
    fclose(fp); /*delete [] *data; *data = NULL;*/ return false;
  }
  if( eType == UnMatch ) BSWAPVEC(data, aryLen);
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ){
    fclose(fp); /*delete [] *data; *data = NULL;*/ return false;
  }
  if( eType == UnMatch ) BSWAP32(dmy);
  if( dmy != aryLen*sizeof(float) ){
    fclose(fp); /*delete [] *data; *data = NULL;*/ return false;
  }

  fclose(fp);
  return true;
}

/**
 * ファイルからSPHデータを読み込む。(double型)
 * @param[in] fname     ファイル名
 * @param[in] dsize     データサイズ
 * @param[out] data_type   データ種別
 * @param[out] voxsize     ボクセルサイズ
 * @param[out] voxorg      原点座標
 * @param[out] voxpit      ピッチ
 * @param[out] step      ステップ数
 * @param[out] time      時間
 * @param[out] data      SPHデータ
 * @return      true=success, false=fault
 */
bool FileIO_SPH::GetData(const char* fname,
                         unsigned long long dsize,
                         int* data_type,
                         unsigned long long voxsize[3],
                         double voxorg[3],
                         double voxpit[3],
                         long long* step,
                         double* time,
                         double* data)
{
  if( !fname || !voxsize || !voxorg || !voxpit
             || !data_type || !step || !time || !data ) return false;

  DataType dt;
  RealType real_type;
  if( !GetInfo(fname, &real_type, &dt, voxsize) ) return false;
  if( real_type != _DOUBLE ) return false;

  EMatchType eType;
  eType = isMatchEndian(fname, 8);
  if( eType == UnKnown ){ return false; }

  FILE* fp;
  if( !(fp = fopen(fname, "rb")) ) {
    fprintf(stderr, "Can't open file.(%s)\n", fname);
    return false;
  }

  unsigned int dmy;
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(dmy);
  if( dmy != 8 ) { fclose(fp); return false; }
  if( fread(data_type, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(*data_type);
  if( fread(&real_type, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(real_type);
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(dmy);
  if( dmy != 8 ) { fclose(fp); return false; }
  if( real_type != _DOUBLE ){ fclose(fp); return false; }

  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(dmy);
  if( dmy != 24 ) { fclose(fp); return false; }
  if( fread(voxsize, sizeof(long long), 3, fp) != 3 ){fclose(fp);return false;}
  if( eType == UnMatch ){
    BSWAP64(voxsize[0]);
    BSWAP64(voxsize[1]);
    BSWAP64(voxsize[2]);
  }
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(dmy);
  if( dmy != 24 ) { fclose(fp); return false; }

  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(dmy);
  if( dmy != 24 ) { fclose(fp); return false; }
  if( fread(voxorg, sizeof(double), 3, fp) != 3 ){fclose(fp);return false;}
  if( eType == UnMatch ){
    BSWAP64(voxorg[0]);
    BSWAP64(voxorg[1]);
    BSWAP64(voxorg[2]);
  }
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(dmy);
  if( dmy != 24 ) { fclose(fp); return false; }

  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(dmy);
  if( dmy != 24 ) { fclose(fp); return false; }
  if( fread(voxpit, sizeof(double), 3, fp) != 3 ){fclose(fp);return false;}
  if( eType == UnMatch ){
    BSWAP64(voxpit[0]);
    BSWAP64(voxpit[1]);
    BSWAP64(voxpit[2]);
  }
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(dmy);
  if( dmy != 24 ) { fclose(fp); return false; }

  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(dmy);
  if( dmy != 16 ) { fclose(fp); return false; }
  if( fread(step, sizeof(long long), 1, fp) != 1 ) {fclose(fp); return false;}
  if( eType == UnMatch ) BSWAP64(*step);
  if( fread(time, sizeof(double), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP64(*time);
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( eType == UnMatch ) BSWAP32(dmy);
  if( dmy != 16 ) { fclose(fp); return false; }

  unsigned int dim;
  if( *data_type == _SCALAR ) dim = 1;
  else if( *data_type == _VECTOR ) dim = 3;
  else { fclose(fp); return false; }

  unsigned long long aryLen = voxsize[0]*voxsize[1]*voxsize[2]*(unsigned long long)dim;
  if( aryLen < 1 ) { fclose(fp); return false; }
  if(dsize < aryLen) { fclose(fp); return false; }
  memset(data, 0, sizeof(double)*aryLen);

  if( fread(&dmy, sizeof(int), 1, fp) != 1 ){
    fclose(fp); /*delete [] *data; *data = NULL;*/ return false;
  }
  if( eType == UnMatch ) BSWAP32(dmy);
  if( dmy != aryLen*sizeof(double) ){
    fclose(fp); /*delete [] *data; *data = NULL;*/ return false;
  }
  if( fread(data, sizeof(double), aryLen, fp) != aryLen ){
    fclose(fp); /*delete [] *data; *data = NULL;*/ return false;
  }
  if( eType == UnMatch ) DBSWAPVEC(data, aryLen);

  if( fread(&dmy, sizeof(int), 1, fp) != 1 ){
    fclose(fp); /*delete [] *data; *data = NULL;*/ return false;
  }
  if( eType == UnMatch ) BSWAP32(dmy);
  if( dmy != aryLen*sizeof(double) ){
    fclose(fp); /*delete [] *data; *data = NULL;*/ return false;
  }

  fclose(fp);
  return true;
}
