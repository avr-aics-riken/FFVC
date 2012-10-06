// #################################################################
//
// FFV : Frontflow / violet
//
// Copyright (c) 2012 All right reserved.
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
//
// #################################################################

/** 
 * @file   PLOT3D.C
 * @brief  FlowBase FileIO_PLOT3D class Header
 * @author kero
 */

#include "PLOT3D.h"


//// CPMクラスポインタのコピー
//void FileIO_PLOT3D::importCPM(cpm_ParaManager* m_paraMngr)
//{
//  //if ( !m_paraMngr ) Exit(0);//Exit0)が未定義でコンパイラに怒られる
//  paraMngr = m_paraMngr;
//}

/**
 * @brief MoveGridを設定する
 * @retval 設定の成否
 * @param is フラグ
 */


void FileIO_PLOT3D::setMoveGrid(const int is){
  switch (is) {
    case GRID_NOT_MOVE:
    case GRID_MOVE:
      P3Op.MoveGrid=is;
      return;
      break;
    default:
      P3Op.MoveGrid=-1;
      return;
      break;
  }
  return;
}

/**
 * @brief Steadyを設定する
 * @retval 設定の成否
 * @param is フラグ
 */
void FileIO_PLOT3D::setSteady(const int is){
  switch (is) {
    case FB_STEADY:
    case FB_UNSTEADY:
      P3Op.Steady=is;
      return;
      break;
    default:
      P3Op.Steady=-1;
      return;
      break;
  }
  return;
}

/**
 * @brief IBlankFlagを設定する
 * @retval 設定の成否
 * @param is フラグ
 */
void FileIO_PLOT3D::setIBlankFlag(const int is){

  switch (is) {
    case NOT_SET_IBLANK:
    case SET_IBLANK:
      P3Op.IBlankFlag=is;
      return;
      break;
    default:
      P3Op.IBlankFlag=-1;
      return;
      break;
  }
  return;
}

/**
 * @brief PLOT3Dのファイルフォーマットを設定する
 * @retval 設定の成否
 * @param is フラグ
 */
void FileIO_PLOT3D::setFormat(const int is)
{
  switch (is) {
    case UNFORMATTED:
    case FORMATTED:
    case C_BINARY:
      P3Op.Format=is; 
      return;
      break;
    default:
      P3Op.Format=-1; 
      return;
      break;
  }
  return;
}


/**
 * @brief Gridデータのセット
 */
void FileIO_PLOT3D::setGridData(
  int m_id,
  int m_jd,
  int m_kd,
  int m_ngrid)
{
  id = m_id;
  jd = m_jd;
  kd = m_kd;
  ngrid = m_ngrid;
}

/**
　* @brief Gridデータのセット
 */
void FileIO_PLOT3D::setGridData2D(
  int m_id,
  int m_jd,
  int m_ngrid)
{
  id = m_id;
  jd = m_jd;
  ngrid = m_ngrid;
}

/**
 * @brief 形状データのセット
 */
void FileIO_PLOT3D::setXYZData(
  float* m_x,
  float* m_y,
  float* m_z,
  int* m_iblank)
{
  x = m_x;
  y = m_y;
  z = m_z;
  iblank = m_iblank;
}

/**
 * @brief 形状データのセット
 */
void FileIO_PLOT3D::setXYZData(
  double* m_x,
  double* m_y,
  double* m_z,
  int* m_iblank)
{
  dx = m_x;
  dy = m_y;
  dz = m_z;
  iblank = m_iblank;
}

/**
 * @brief 形状データのセット
 */
void FileIO_PLOT3D::setXYZData2D(
  float* m_x,
  float* m_y,
  int* m_iblank)
{
  x = m_x;
  y = m_y;
  iblank = m_iblank;
}

/**
 * @brief 形状データのセット
 */
void FileIO_PLOT3D::setXYZData2D(
  double* m_x,
  double* m_y,
  int* m_iblank)
{
  dx = m_x;
  dy = m_y;
  iblank = m_iblank;
}

/**
 * @brief 結果Qデータのセット
 */
void FileIO_PLOT3D::setQData(
  float m_fsmach,
  float m_alpha,
  float m_re,
  float m_time,
  float* m_q)
{
  fsmach  = m_fsmach;
  alpha   = m_alpha;
  re      = m_re;
  time    = m_time;
  q       = m_q;
}

/**
 * @brief 結果Qデータのセット
 */
void FileIO_PLOT3D::setQData(
  double m_fsmach,
  double m_alpha,
  double m_re,
  double m_time,
  double* m_q)
{
  dfsmach  = m_fsmach;
  dalpha   = m_alpha;
  dre      = m_re;
  dtime    = m_time;
  dq       = m_q;
}

/**
 * @brief 結果Funcデータ項目数のセット
 */
void FileIO_PLOT3D::setFuncDataNum(
  int m_nvar)
{
  nvar  = m_nvar;
}

/**
 * @brief 結果Funcデータのセット
 */
void FileIO_PLOT3D::setFuncData(
  float* m_d)
{
  d  = m_d;
}

/**
 * @brief 結果Funcデータのセット
 */
void FileIO_PLOT3D::setFuncData(
  double* m_d)
{
  dd  = m_d;
}


/**
 * @brief PLOT3Dの出力ファイルオープン
 */
bool FileIO_PLOT3D::OpenFile()
{
  int ierror=0;
  switch (P3Op.Format) {
    case UNFORMATTED:
    case FORMATTED:

      char tmp[FB_FILE_PATH_LENGTH];
      memset(tmp, 0, sizeof(char)*FB_FILE_PATH_LENGTH);
      strcpy(tmp, fname.c_str());

      int len = strlen(tmp);
      int flag=P3Op.Format;
      int fnsize=FB_FILE_PATH_LENGTH;

      //open_plot3d_file_(&flag,tmp,&ifl,&fnsize,&ierror);
      open_plot3d_file_(&flag,tmp,&ifl,&len,&ierror);

      break;

    case C_BINARY:
      if( (fp = fopen(fname.c_str(), "wb")) == NULL ) {
        fprintf(stderr, "Can't open file.(%s)\n", fname.c_str());
        return false;
      }
      ierror = true;
      break;
  }

  return ierror;
}

/**
 * @brief PLOT3Dのファイルクローズ
 */
void FileIO_PLOT3D::CloseFile()
{
  switch (P3Op.Format) {
    case UNFORMATTED:
    case FORMATTED:
      close_plot3d_file_(&ifl);
      break;
    case C_BINARY:
      fclose(fp);
      break;
  }

}
