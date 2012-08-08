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
 * @file   PLOT3D_read.C
 * @brief  FlowBase FileIO_PLOT3D_READ class Header
 * @author kero
 */

#include "PLOT3D_read.h"


/**
 * @brief PLOT3Dの入力ファイルオープン
 */
bool FileIO_PLOT3D_READ::OpenFile()
{
  //set file name
  char tmp[FB_FILE_PATH_LENGTH];
  memset(tmp, 0, sizeof(char)*FB_FILE_PATH_LENGTH);
  strcpy(tmp, fname.c_str());

  int len = strlen(tmp);

  //open file
  int flag=P3Op.Format;
  int fnsize=FB_FILE_PATH_LENGTH;
  int ierror=0;
  //open_plot3d_inputfile_(&flag,tmp,&ifl,&fnsize,&ierror);
  open_plot3d_inputfile_(&flag,tmp,&ifl,&len,&ierror);

  return ierror;
}


/**
 * @brief グリッド数の読み込み
 */
void FileIO_PLOT3D_READ::ReadNgrid(int* ngrid)
{
  //single grid
  if( P3Op.GridKind == SINGLE_GRID )
  {
    *ngrid=1;
    return;
  }

  //multigrid
  switch (P3Op.Format) {
    case UNFORMATTED:
      read_ngrid_data_(ngrid, &ifl);
      break;
    case FORMATTED:
      read_ngrid_data_formatted_(ngrid, &ifl);
      break;
    case UNFORMATTED_SPECIAL:
      break;
    default:
      break;
  }
}

/**
 * @brief ブロックデータの読み込み
 */
void FileIO_PLOT3D_READ::ReadBlockData(int* id, int* jd, int* kd)
{
  switch (P3Op.Format) {
    case UNFORMATTED:
      read_block_data_(id, jd, kd, &ifl);
      break;
	case FORMATTED:
      read_block_data_formatted_(id, jd, kd, &ifl);
      break;
    case UNFORMATTED_SPECIAL:
      break;
    default:
      break;
  }
}

/**
 * @brief ブロックデータの読み込み（2D）
 */
void FileIO_PLOT3D_READ::ReadBlockData2D(int* id, int* jd)
{
  switch (P3Op.Format) {
    case UNFORMATTED:
      read_block_data_2d_(id, jd, &ifl);
      break;
    case FORMATTED:
      read_block_data_2d_formatted_(id, jd, &ifl);
      break;
    case UNFORMATTED_SPECIAL:
      break;
    default:
      break;
  }
}

/**
 * @brief 計算結果Functionファイル、ブロックデータの読み込み
 */
void FileIO_PLOT3D_READ::ReadFuncBlockData(int* id, int* jd, int* kd, int* nvar, int ngrid)
{
  switch (P3Op.DimIs) {
    case DIMENSION_2D:
      switch (P3Op.Format) {
        case UNFORMATTED:
          read_func_block_data_2d_(id, jd, nvar, &ngrid, &ifl);
          break;
        case FORMATTED:
          read_func_block_data_2d_formatted_(id, jd, nvar, &ngrid, &ifl);
          break;
        case UNFORMATTED_SPECIAL:
          break;
        default:
          break;
      }
    case DIMENSION_3D:
      switch (P3Op.Format) {
        case UNFORMATTED:
          read_func_block_data_(id, jd, kd, nvar, &ngrid, &ifl);
          break;
        case FORMATTED:
          read_func_block_data_formatted_(id, jd, kd, nvar, &ngrid, &ifl);
          break;
        case UNFORMATTED_SPECIAL:
          break;
        default:
          break;
      }
      break;
  }
}


/**
 * @brief グリッド座標ファイルの読み込み（UNFORMATTED）
 */
void FileIO_PLOT3D_READ::ReadXYZ_UNFORMATTED()
{
  int d_type = (sizeof(REAL_TYPE) == 4) ? 1 : 2;  // 1-float / 2-double
  if(d_type==1)
  {
    switch (P3Op.DimIs) {
      case DIMENSION_2D:
        read_xyz_2d_(&id, &jd, x, y, &ifl);
        break;
      case DIMENSION_3D:
        if(P3Op.IBlankFlag == NOT_SET_IBLANK){
          read_xyz_3d_(&id, &jd, &kd, x, y, z, &ifl);
        }
        else if(P3Op.IBlankFlag == SET_IBLANK){
          read_xyz_3d_iblank_(&id, &jd, &kd, x, y, z, iblank, &ifl);
        }
        break;
    }
  }
  else if(d_type==2)
  {
    switch (P3Op.DimIs) {
      case DIMENSION_2D:
        dread_xyz_2d_(&id, &jd, x, y, &ifl);
        break;
      case DIMENSION_3D:
        if(P3Op.IBlankFlag == NOT_SET_IBLANK){
          dread_xyz_3d_(&id, &jd, &kd, x, y, z, &ifl);
        }
        else if(P3Op.IBlankFlag == SET_IBLANK){
          dread_xyz_3d_iblank_(&id, &jd, &kd, x, y, z, iblank, &ifl);
        }
        break;
    }
  }
}

/**
 * @brief グリッド座標ファイルの読み込み（FORMATTED）
 */
void FileIO_PLOT3D_READ::ReadXYZ_FORMATTED()
{
  int d_type = (sizeof(REAL_TYPE) == 4) ? 1 : 2;  // 1-float / 2-double
  if(d_type==1)
  {
    switch (P3Op.DimIs) {
      case DIMENSION_2D:
        read_xyz_2d_formatted_(&id, &jd, x, y, &ifl);
        break;
      case DIMENSION_3D:
        if(P3Op.IBlankFlag == NOT_SET_IBLANK){
          read_xyz_3d_formatted_(&id, &jd, &kd, x, y, z, &ifl);
        }
        else if(P3Op.IBlankFlag == SET_IBLANK){
          read_xyz_3d_iblank_formatted_(&id, &jd, &kd, x, y, z, iblank, &ifl);
        }
        break;
    }
  }
  else if(d_type==2)
  {
    switch (P3Op.DimIs) {
      case DIMENSION_2D:
        read_xyz_2d_formatted_(&id, &jd, x, y, &ifl);
        break;
      case DIMENSION_3D:
        if(P3Op.IBlankFlag == NOT_SET_IBLANK){
          dread_xyz_3d_formatted_(&id, &jd, &kd, x, y, z, &ifl);
        }
        else if(P3Op.IBlankFlag == SET_IBLANK){
          dread_xyz_3d_iblank_formatted_(&id, &jd, &kd, x, y, z, iblank, &ifl);
        }
        break;
    }
  }
}


/**
 * @brief グリッド座標ファイルの読み込み（UNFORMATTED）特別仕様フォーマット
 */
void FileIO_PLOT3D_READ::ReadXYZ_UNFORMATTED_SPECIAL()
{

}

/**
 * @brief 計算結果Qファイルの読み込み（UNFORMATTED）
 */
void FileIO_PLOT3D_READ::ReadQ_UNFORMATTED()
{
  int d_type = (sizeof(REAL_TYPE) == 4) ? 1 : 2;  // 1-float / 2-double
  if(d_type==1)
  {
    switch (P3Op.DimIs) {
      case DIMENSION_2D:
        read_q_2d_(&id, &jd, &fsmach, &alpha, &re, &time, q, &ifl);
        break;
      case DIMENSION_3D:
        read_q_3d_(&id, &jd, &kd, &fsmach, &alpha, &re, &time, q, &ifl);
        break;
    }
  }
  else if(d_type==2)
  {
    switch (P3Op.DimIs) {
      case DIMENSION_2D:
        dread_q_2d_(&id, &jd, &fsmach, &alpha, &re, &time, q, &ifl);
        break;
      case DIMENSION_3D:
        dread_q_3d_(&id, &jd, &kd, &fsmach, &alpha, &re, &time, q, &ifl);
        break;
    }
  }
}

/**
 * @brief 計算結果Qファイルの読み込み（FORMATTED）
 */
void FileIO_PLOT3D_READ::ReadQ_FORMATTED()
{
  int d_type = (sizeof(REAL_TYPE) == 4) ? 1 : 2;  // 1-float / 2-double
  if(d_type==1)
  {
    switch (P3Op.DimIs) {
      case DIMENSION_2D:
        read_q_2d_formatted_(&id, &jd, &fsmach, &alpha, &re, &time, q, &ifl);
        break;
      case DIMENSION_3D:
        read_q_3d_formatted_(&id, &jd, &kd, &fsmach, &alpha, &re, &time, q, &ifl);
        break;
    }
  }
  else if(d_type==2)
  {
    switch (P3Op.DimIs) {
      case DIMENSION_2D:
        dread_q_2d_formatted_(&id, &jd, &fsmach, &alpha, &re, &time, q, &ifl);
        break;
      case DIMENSION_3D:
        dread_q_3d_formatted_(&id, &jd, &kd, &fsmach, &alpha, &re, &time, q, &ifl);
        break;
    }
  }
}

/**
 * @brief 計算結果Qファイルの読み込み（UNFORMATTED）特別仕様フォーマット
 */
void FileIO_PLOT3D_READ::ReadQ_UNFORMATTED_SPECIAL()
{

}


/**
 * @brief 計算結果Functionファイルの読み込み（UNFORMATTED）
 */
void FileIO_PLOT3D_READ::ReadFunc_UNFORMATTED()
{
  int d_type = (sizeof(REAL_TYPE) == 4) ? 1 : 2;  // 1-float / 2-double
  if(d_type==1)
  {
    switch (P3Op.DimIs) {
      case DIMENSION_2D:
        read_func_2d_(&id, &jd, &nvar, d, &ifl);
        break;
      case DIMENSION_3D:
        read_func_3d_(&id, &jd, &kd, &nvar, d, &ifl);
        break;
    }
  }
  else if(d_type==2)
  {
    switch (P3Op.DimIs) {
      case DIMENSION_2D:
        dread_func_2d_(&id, &jd, &nvar, d, &ifl);
        break;
      case DIMENSION_3D:
        dread_func_3d_(&id, &jd, &kd, &nvar, d, &ifl);
        break;
    }
  }
}

/**
 * @brief 計算結果Functionファイルの読み込み（FORMATTED）
 */
void FileIO_PLOT3D_READ::ReadFunc_FORMATTED()
{
  int d_type = (sizeof(REAL_TYPE) == 4) ? 1 : 2;  // 1-float / 2-double
  if(d_type==1)
  {
    switch (P3Op.DimIs) {
      case DIMENSION_2D:
        read_func_2d_formatted_(&id, &jd, &nvar, d, &ifl);
        break;
      case DIMENSION_3D:
        read_func_3d_formatted_(&id, &jd, &kd, &nvar, d, &ifl);
        break;
    }
  }
  else if(d_type==2)
  {
    switch (P3Op.DimIs) {
      case DIMENSION_2D:
        dread_func_2d_formatted_(&id, &jd, &nvar, d, &ifl);
        break;
      case DIMENSION_3D:
        dread_func_3d_formatted_(&id, &jd, &kd, &nvar, d, &ifl);
        break;
    }
  }
}

/**
 * @brief 計算結果Functionファイルの読み込み（UNFORMATTED）特別仕様フォーマット
 */
void FileIO_PLOT3D_READ::ReadFunc_UNFORMATTED_SPECIAL()
{

}

/**
 * @brief グリッド座標ファイルの読み込み
 */
bool FileIO_PLOT3D_READ::ReadXYZData()
{
  switch (P3Op.Format) {
    case UNFORMATTED:
      ReadXYZ_UNFORMATTED();
      break;
    case FORMATTED:
      ReadXYZ_FORMATTED();
      break;
    case UNFORMATTED_SPECIAL:
      ReadXYZ_UNFORMATTED_SPECIAL();
      break;
    default:
      return false;
      break;
  }
  return true;
}


/**
 * @brief 計算結果Qファイルの読み込み
 */
bool FileIO_PLOT3D_READ::ReadQData(
  REAL_TYPE* m_fsmach,
  REAL_TYPE* m_alpha,
  REAL_TYPE* m_re,
  REAL_TYPE* m_time,
  REAL_TYPE* m_q)
{
  //set pointer
  q=m_q;

  switch (P3Op.Format) {
    case UNFORMATTED:
      ReadQ_UNFORMATTED();
      break;
    case FORMATTED:
      ReadQ_FORMATTED();
      break;
    case UNFORMATTED_SPECIAL:
      ReadQ_UNFORMATTED_SPECIAL();
      break;
    default:
      return false;
      break;
  }

  *m_fsmach=fsmach;
  *m_alpha=alpha;
  *m_re=re;
  *m_time=time;

  return true;
}

/**
 * @brief 計算結果Functionファイルの読み込み
 */
bool FileIO_PLOT3D_READ::ReadFuncData()
{
  switch (P3Op.Format) {
    case UNFORMATTED:
      ReadFunc_UNFORMATTED();
      break;
    case FORMATTED:
      ReadFunc_FORMATTED();
      break;
    case UNFORMATTED_SPECIAL:
      ReadFunc_UNFORMATTED_SPECIAL();
      break;
    default:
      return false;
      break;
  }
  return true;
}

/**
 * @brief FunctionNameファイルの読み込み
 */
void FileIO_PLOT3D_READ::ReadFunctionName(string* fn)
{
  int iend=0;
  int len=FB_BUFF_LENGTH;
  char tmp[FB_BUFF_LENGTH];

  memset(tmp, '\0', sizeof(char)*FB_BUFF_LENGTH);
  read_line_(tmp,&ifl,&len,&iend);
  *fn=tmp;
  return;
}


/**
 * @brief 文字列の読み込み
 */
void FileIO_PLOT3D_READ::ReadSTRING(
  string* m_buff)
{
  int iend=0;
  int len=FB_BUFF_LENGTH;
  char tmp[FB_BUFF_LENGTH];

  memset(tmp, '\0', sizeof(char)*FB_BUFF_LENGTH);
  read_line_(tmp,&ifl,&len,&iend);
  *m_buff=tmp;
  return;
}

/**
 * @brief FVBNDファイルの読み込み
 */
void FileIO_PLOT3D_READ::ReadFVBND(
  int* m_type, int* m_gridnum,
  int* m_Imin,int* m_Imax,
  int* m_Jmin,int* m_Jmax,
  int* m_Kmin,int* m_Kmax,
  string* m_ResultFlag,
  int* m_dir)
{
  int type,gridnum;
  int Imin,Imax,Jmin,Jmax,Kmin,Kmax;
  int dir;

  char tmp[FB_BUFF_LENGTH];
  string buff;

//read boundaries
  memset(tmp, '\0', sizeof(char)*FB_BUFF_LENGTH);
  read_fvbnd_boundary_(
    &type,&gridnum,&Imin,&Imax,&Jmin,&Jmax,&Kmin,&Kmax,tmp,&dir,&ifl);

  *m_type=type;
  *m_gridnum=gridnum;
  *m_Imin=Imin;
  *m_Imax=Imax;
  *m_Jmin=Jmin;
  *m_Jmax=Jmax;
  *m_Kmin=Kmin;
  *m_Kmax=Kmax;
  *m_ResultFlag=tmp;
  *m_dir=dir;
}

///**
// * @brief FVBNDファイルの読み込み
// */
//void FileIO_PLOT3D_READ::ReadFVBND(
//  const int nbname, const int nb,
//  string* boundary_name,
//  int* type, int* gridnum,
//  int* Imin,int* Imax,
//  int* Jmin,int* Jmax,
//  int* Kmin,int* Kmax,
//  string* ResultFlag,
//  int* dir)
//{
//  int iend=0;
//  int len=FB_BUFF_LENGTH;
//  char tmp[FB_BUFF_LENGTH];
//  string buff;
//
////read header "FVBND 1 4"
//  read_line_(tmp,&ifl,&len,&iend);
//  if(iend) return;
//
////read boundary name
//  for(int i=0;i<nbname;i++){
//    memset(tmp, '\0', sizeof(char)*FB_BUFF_LENGTH);
//    read_line_(tmp,&ifl,&len,&iend);
//    boundary_name[i]=tmp;
//    if(iend) return;
//  }
//
////read boundaries
//  memset(tmp, '\0', sizeof(char)*FB_BUFF_LENGTH);
//  read_line_(tmp,&ifl,&len,&iend);//read "BOUNDARIES"
//  for(int i=0;i<nb;i++){
//    memset(tmp, '\0', sizeof(char)*FB_BUFF_LENGTH);
//    read_fvbnd_boundary_(
//      &type[i],&gridnum[i],&Imin[i],&Imax[i],&Jmin[i],&Jmax[i],&Kmin[i],&Kmax[i],tmp,&dir[i],&ifl);
//    ResultFlag[i]=tmp;
//  }
//
//}
