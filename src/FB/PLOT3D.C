// #################################################################
//
// FFV : Frontflow / violet
//
// Copyright (c) All right reserved. 2012
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
//
// #################################################################

/** 
 * @file PLOT3D.C
 * @brief FlowBase FileIO_PLOT3D class Header
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
bool FileIO_PLOT3D::setMoveGrid(const int is){
  switch (is) {
    case GRID_NOT_MOVE:
    case GRID_MOVE:
      P3Op.MoveGrid=is;
      return true;
      break;
    default:
      return false;
      break;
  }
  return false;
}

/**
 * @brief Steadyを設定する
 * @retval 設定の成否
 * @param is フラグ
 */
bool FileIO_PLOT3D::setSteady(const int is){
  switch (is) {
    case FB_STEADY:
    case FB_UNSTEADY:
      P3Op.Steady=is;
      return true;
      break;
    default:
      return false;
      break;
  }
  return false;
}

/**
 * @brief IBlankFlagを設定する
 * @retval 設定の成否
 * @param is フラグ
 */
bool FileIO_PLOT3D::setIBlankFlag(const int is){

  switch (is) {
    case NOT_SET_IBLANK:
    case SET_IBLANK:
      P3Op.IBlankFlag=is;
      return true;
      break;
    default:
      return false;
      break;
  }
  return false;
}

/**
 * @brief PLOT3Dのファイルフォーマットを設定する
 * @retval 設定の成否
 * @param is フラグ
 */
bool FileIO_PLOT3D::setFormat(const int is)
{
  switch (is) {
    case UNFORMATTED:
    case FORMATTED:
    case UNFORMATTED_SPECIAL:
      P3Op.Format=is; 
      return true;
      break;
    default:
      return false;
      break;
  }
  return false;
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
  REAL_TYPE* m_x,
  REAL_TYPE* m_y,
  REAL_TYPE* m_z,
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
void FileIO_PLOT3D::setXYZData2D(
  REAL_TYPE* m_x,
  REAL_TYPE* m_y,
  int* m_iblank)
{
  x = m_x;
  y = m_y;
  iblank = m_iblank;
}

/**
 * @brief 結果Qデータのセット
 */
void FileIO_PLOT3D::setQData(
  REAL_TYPE m_fsmach,
  REAL_TYPE m_alpha,
  REAL_TYPE m_re,
  REAL_TYPE m_time,
  REAL_TYPE* m_q)
{
  fsmach  = m_fsmach;
  alpha   = m_alpha;
  re      = m_re;
  time    = m_time;
  q       = m_q;
}

/**
 * @brief 結果Funcデータのセット
 */
void FileIO_PLOT3D::setFuncData(
  int m_nvar,
  REAL_TYPE* m_d)
{
  nvar  = m_nvar;
  d  = m_d;
}

/**
 * @brief PLOT3Dの出力ファイルオープン
 */
bool FileIO_PLOT3D::OpenFile()
{
  //set file name
  char tmp[FB_FILE_PATH_LENGTH];
  memset(tmp, 0, sizeof(char)*FB_FILE_PATH_LENGTH);
  strcpy(tmp, fname.c_str());

  int len = strlen(tmp);
  //std::cout << "len = " << len << std::endl;

  //open file
  int flag=P3Op.Format;
  int fnsize=FB_FILE_PATH_LENGTH;
  int ierror=0;
  //open_plot3d_file_(&flag,tmp,&ifl,&fnsize,&ierror);
  open_plot3d_file_(&flag,tmp,&ifl,&len,&ierror);

  return ierror;
}

/**
 * @brief PLOT3Dのファイルクローズ
 */
void FileIO_PLOT3D::CloseFile()
{
  close_plot3d_file_(&ifl);
}
