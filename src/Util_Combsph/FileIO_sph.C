// #################################################################
//
// FileIO_SPH
//
// Copyright (c) All right reserved. 2012
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
//
// #################################################################

/** 
 * @file   FileIO_SPH_sph.C
 * @brief  FileIO_SPH Class
 * @author kero
 */

#include "FileIO_sph.h"


/**
 * ファイル出力情報を設定する。(float型)
 * @param dtype     データ型
 * @param voxsize   ボクセルサイズ
 * @param voxorg    原点座標
 * @param voxpit    ピッチ
 * @param step    ステップ数
 * @param time      時間
 * @param data      SPHデータ
 * @return      true=success, false=fault
 */
bool FileIO_SPH::Init(int dtype,
                      unsigned int voxsize[3],
                      float voxorg[3],
                      float voxpit[3],
                      int* step,
                      float* time,
                      float** data)
{
  if( !voxsize || !voxorg || !voxpit
      || !step || !time || !data || !(*data) ) return false;
  if( !((dtype==1)||(dtype==2)) ) return false;

  m_outType = _OUT_METHOD1;
  m_rType = _FLOAT;
  m_dType = (DataType)dtype;
  m_iSize[0] = voxsize[0];
  m_iSize[1] = voxsize[1];
  m_iSize[2] = voxsize[2];
  m_fOrg[0] = voxorg[0];
  m_fOrg[1] = voxorg[1];
  m_fOrg[2] = voxorg[2];
  m_fPit[0] = voxpit[0];
  m_fPit[1] = voxpit[1];
  m_fPit[2] = voxpit[2];
  m_fRefStep = step;
  m_fRefTime = time;
  m_fData = data;

  return true;
}

/**
 * ファイル出力情報を設定する。(double型)
 * @param dtype     データ型
 * @param voxsize   ボクセルサイズ
 * @param voxorg    原点座標
 * @param voxpit    ピッチ
 * @param step    ステップ数
 * @param time      時間
 * @param data      SPHデータ
 * @return      true=success, false=fault
 */
bool FileIO_SPH::Init(int dtype,
                      unsigned long long voxsize[3],
                      double voxorg[3],
                      double voxpit[3],
                      long long* step,
                      double* time,
                      double** data)
{
  if( !voxsize || !voxorg || !voxpit
      || !step || !time || !data || !(*data) ) return false;
  if( !((dtype==1)||(dtype==2)) ) return false;

  m_outType = _OUT_METHOD1;
  m_rType = _DOUBLE;
  m_dType = (DataType)dtype;
  m_lSize[0] = voxsize[0];
  m_lSize[1] = voxsize[1];
  m_lSize[2] = voxsize[2];
  m_dOrg[0] = voxorg[0];
  m_dOrg[1] = voxorg[1];
  m_dOrg[2] = voxorg[2];
  m_dPit[0] = voxpit[0];
  m_dPit[1] = voxpit[1];
  m_dPit[2] = voxpit[2];
  m_dRefStep = step;
  m_dRefTime = time;
  m_dData = data;

  return true;
}

/**
 * ファイル出力情報を設定する。(ヘッダーレコード:float型)
 * @param dtype     データ型
 * @param voxsize   ボクセルサイズ
 * @param voxorg    原点座標
 * @param voxpit    ピッチ
 * @return      true=success, false=fault
 */
bool FileIO_SPH::Init(int dtype,
                      unsigned int voxsize[3],
                      float voxorg[3],
                      float voxpit[3])
{
  if( !voxsize || !voxorg || !voxpit ) return false;
  if( !((dtype==1)||(dtype==2)) ) return false;

  m_outType = _OUT_METHOD2;
  m_rType = _FLOAT;
  m_dType = (DataType)dtype;
  m_iSize[0] = voxsize[0];
  m_iSize[1] = voxsize[1];
  m_iSize[2] = voxsize[2];
  m_fOrg[0] = voxorg[0];
  m_fOrg[1] = voxorg[1];
  m_fOrg[2] = voxorg[2];
  m_fPit[0] = voxpit[0];
  m_fPit[1] = voxpit[1];
  m_fPit[2] = voxpit[2];

  return true;
}


/**
 * ファイル出力情報を設定する。(ヘッダーレコード:double型)
 * @param dtype     データ型
 * @param voxsize   ボクセルサイズ
 * @param voxorg    原点座標
 * @param voxpit    ピッチ
 * @return      true=success, false=fault
 */
bool FileIO_SPH::Init(int dtype,
                      unsigned long long voxsize[3],
                      double voxorg[3],
                      double voxpit[3])
{
  if( !voxsize || !voxorg || !voxpit ) return false;
  if( !((dtype==1)||(dtype==2)) ) return false;

  m_outType = _OUT_METHOD2;
  m_rType = _DOUBLE;
  m_dType = (DataType)dtype;
  m_lSize[0] = voxsize[0];
  m_lSize[1] = voxsize[1];
  m_lSize[2] = voxsize[2];
  m_dOrg[0] = voxorg[0];
  m_dOrg[1] = voxorg[1];
  m_dOrg[2] = voxorg[2];
  m_dPit[0] = voxpit[0];
  m_dPit[1] = voxpit[1];
  m_dPit[2] = voxpit[2];

  return true;
}

/**
 * SPHファイルに出力を行う。
 * @param step    ステップ数
 * @return      true=success, false=fault
 */
bool FileIO_SPH::OutputBinaryFile(unsigned int step)
{
  if( m_outType != _OUT_METHOD1 ) return false;

  char* fname;
  if( !(fname = MakeOutFileName(step)) ) return false;

  FILE* fp;
  if( (fp = fopen(fname, "wb")) == NULL ) {
    fprintf(stderr, "Can't open file.(%s)\n", fname);
    return false;
  }
  if( !WriteData(fp) ){ fclose(fp); return false; }
  fclose(fp);

  delete [] fname;
  return true;
}


/**
 * SPHファイルに出力を行う。(float型)
 * @param step    ステップ数
 * @param time      時間
 * @param data      SPHデータ
 * @return      true=success, false=fault
 */
bool FileIO_SPH::OutputBinaryFile(int step, float time, const float* data)
{
  if( !data ) return false;
  if( m_outType != _OUT_METHOD2 ) return false;
  if( m_rType != _FLOAT ) return false;

  char* fname;
  if( !(fname = MakeOutFileName(step)) ) return false;

  m_fStep = step;
  m_fTime = time;
  m_fRefStep = &m_fStep;
  m_fRefTime = &m_fTime;
  m_fData = &data;

  FILE* fp;
  if( (fp = fopen(fname, "wb")) == NULL ) {
    fprintf(stderr, "Can't open file.(%s)\n", fname);
    return false;
  }
  if( !WriteData(fp) ){ fclose(fp); return false; }
  fclose(fp);

  delete [] fname;
  return true;
}

/**
 * SPHファイルに出力を行う。(double型)
 * @param step    ステップ数
 * @param time      時間
 * @param data      SPHデータ
 * @return      true=success, false=fault
 */
bool FileIO_SPH::OutputBinaryFile(long long step,
                                  double time,
                                  const double* data)
{
  if( !data ) return false;
  if( m_outType != _OUT_METHOD2 ) return false;
  if( m_rType != _DOUBLE ) return false;

  char* fname;
  if( !(fname = MakeOutFileName(step)) ) return false;

  m_dStep = step;
  m_dTime = time;
  m_dRefStep = &m_dStep;
  m_dRefTime = &m_dTime;
  m_dData = &data;

  FILE* fp;
  if( (fp = fopen(fname, "wb")) == NULL ) {
    fprintf(stderr, "Can't open file.(%s)\n", fname);
    return false;
  }
  if( !WriteData(fp) ){ fclose(fp); return false; }
  fclose(fp);

  delete [] fname;
  return true;
}

/**
 * 出力ファイル名を作成する。
 * @param num   ステップ数
 * @return    出力ファイル名
 */
char* FileIO_SPH::MakeOutFileName(unsigned int num)
{
  if( !m_fBase ) return NULL;

  unsigned nameLen = strlen(m_fBase) + m_fSeqNumLen + 1;
  if( m_fSuffix ) nameLen = nameLen + strlen(m_fSuffix) + 1;
  char* fname = new char[nameLen];
  if( !fname ) return NULL;

  strcpy(fname, m_fBase);
  char str1[32], str2[32];
  sprintf(str1, "%%0%dd", (int)m_fSeqNumLen);
  sprintf(str2, str1, num);
  strcat(fname, str2);
  if( m_fSuffix ){
    strcat(fname, ".");
    strcat(fname, m_fSuffix);
  }

  return fname;
}

/**
 * SPHファイルに書込を行う。
 * @param fp    ファイルポインタ
 * @return      true=success, false=fault
 */
bool FileIO_SPH::WriteData(FILE* fp)
{
  if( !fp ) return false;
  if( (m_rType == _REAL_UNKNOWN) ) return false;
  if( (m_dType == _DATA_UNKNOWN) ) return false;
  if( (m_rType==_FLOAT)&&(!m_fRefStep||!m_fRefTime) ) return false;
  if( (m_rType==_DOUBLE)&&(!m_dRefStep||!m_dRefTime) ) return false;
  if( (m_rType==_FLOAT)&&(!(*m_fData)) ) return false;
  if( (m_rType==_DOUBLE)&&(!(*m_dData)) ) return false;

  unsigned int dmy;
  dmy = 8;
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
  if( fwrite(&m_dType, sizeof(int), 1, fp) != 1 ) return false;
  if( fwrite(&m_rType, sizeof(int), 1, fp) != 1 ) return false;
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;

  if( m_rType == _FLOAT ){
    dmy = 12;
    if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
    if( fwrite(m_iSize, sizeof(int), 3, fp) != 3 ) return false;
    if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;

    if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
    if( fwrite(m_fOrg, sizeof(float), 3, fp) != 3 ) return false;
    if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;

    if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
    if( fwrite(m_fPit, sizeof(float), 3, fp) != 3 ) return false;
    if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
  } else { // _DOUBLE
    dmy = 24;
    if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
    if( fwrite(m_lSize, sizeof(long long), 3, fp) != 3 ) return false;
    if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;

    if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
    if( fwrite(m_dOrg, sizeof(double), 3, fp) != 3 ) return false;
    if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;

    if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
    if( fwrite(m_dPit, sizeof(double), 3, fp) != 3 ) return false;
    if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
  }

  if( m_rType == _FLOAT ) dmy = 8;
  else dmy = 16;
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
  if( m_rType == _FLOAT ){
    if( fwrite(m_fRefStep, sizeof(int), 1, fp) != 1 ) return false;
    if( fwrite(m_fRefTime, sizeof(float), 1, fp) != 1 ) return false;
  }else{
    if( fwrite(m_dRefStep, sizeof(long long), 1, fp) != 1 ) return false;
    if( fwrite(m_dRefTime, sizeof(double), 1, fp) != 1 ) return false;
  }
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;

  size_t dLen;
  if( m_rType == _FLOAT )
    dLen = (size_t)(m_iSize[0] * m_iSize[1] * m_iSize[2]);
  else
    dLen = (size_t)(m_lSize[0] * m_lSize[1] * m_lSize[2]);

  if( m_dType == _VECTOR ) dLen *= 3;
  if( m_rType == _FLOAT ) dmy = dLen * sizeof(float);
  else dmy = dLen * sizeof(double);
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
  if( m_rType == _FLOAT ){
    if( fwrite(*m_fData, sizeof(float), dLen, fp) != dLen ) return false;
  }else{
    if( fwrite(*m_dData, sizeof(double), dLen, fp) != dLen ) return false;
  }
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;

  return true;
}

