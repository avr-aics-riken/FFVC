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
 * @file   PLOT3D_write.C
 * @brief  FlowBase FileIO_PLOT3D_WRITE class Header
 * @author kero
 */

#include "PLOT3D_write.h"


/**
 * @brief PLOT3Dの出力ファイルオープン
 */
bool FileIO_PLOT3D_WRITE::OpenFile()
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
  //open_plot3d_outputfile_(&flag,tmp,&ifl,&fnsize,&ierror);
  open_plot3d_outputfile_(&flag,tmp,&ifl,&len,&ierror);

  return ierror;
}

/**
 * @brief グリッド数の書き出し
 */
void FileIO_PLOT3D_WRITE::WriteNgrid(int ngrid)
{
  switch (P3Op.Format) {
    case UNFORMATTED:
      write_ngrid_data_(&ngrid, &ifl);
      break;
    case FORMATTED:
      write_ngrid_data_formatted_(&ngrid, &ifl);
      break;
    case UNFORMATTED_SPECIAL:
      break;
    default:
      break;
  }
}

/**
 * @brief ブロックデータの書き出し
 */
void FileIO_PLOT3D_WRITE::WriteBlockData(int id, int jd, int kd)
{
  switch (P3Op.Format) {
    case UNFORMATTED:
      write_block_data_(&id, &jd, &kd, &ifl);
      break;
    case FORMATTED:
      write_block_data_formatted_(&id, &jd, &kd, &ifl);
      break;
    case UNFORMATTED_SPECIAL:
      break;
    default:
      break;
  }
}

/**
 * @brief ブロックデータの書き出し（2D）
 */
void FileIO_PLOT3D_WRITE::WriteBlockData2D(int id, int jd)
{
  switch (P3Op.Format) {
    case UNFORMATTED:
      write_block_data_2d_(&id, &jd, &ifl);
      break;
    case FORMATTED:
      write_block_data_2d_formatted_(&id, &jd, &ifl);
      break;
    case UNFORMATTED_SPECIAL:
      break;
    default:
      break;
  }
}

/**
 * @brief 計算結果Functionファイル、ブロックデータの書き出し
 */
void FileIO_PLOT3D_WRITE::WriteFuncBlockData(int* id, int* jd, int* kd, int* nvar, int ngrid)
{
  switch (P3Op.DimIs) {
    case DIMENSION_2D:
      switch (P3Op.Format) {
        case UNFORMATTED:
          write_func_block_data_2d_(id, jd, nvar, &ngrid, &ifl);
          break;
        case FORMATTED:
          write_func_block_data_2d_formatted_(id, jd, nvar, &ngrid, &ifl);
          break;
        case UNFORMATTED_SPECIAL:
          break;
        default:
          break;
      }
    case DIMENSION_3D:
      switch (P3Op.Format) {
        case UNFORMATTED:
          write_func_block_data_(id, jd, kd, nvar, &ngrid, &ifl);
          break;
        case FORMATTED:
          write_func_block_data_formatted_(id, jd, kd, nvar, &ngrid, &ifl);
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
 * @brief グリッド座標ファイルの書き出し（UNFORMATTED）
 */
void FileIO_PLOT3D_WRITE::WriteXYZ_UNFORMATTED()
{
  int d_type = (sizeof(REAL_TYPE) == 4) ? 1 : 2;  // 1-float / 2-double
  if(d_type==1)
  {
    switch (P3Op.DimIs) {
      case DIMENSION_2D:
        write_xyz_2d_(&id, &jd, x, y, &ifl);
        break;
      case DIMENSION_3D:
        if(P3Op.IBlankFlag == NOT_SET_IBLANK){
          write_xyz_3d_(&id, &jd, &kd, x, y, z, &ifl);
        }
        else if(P3Op.IBlankFlag == SET_IBLANK){
          write_xyz_3d_iblank_(&id, &jd, &kd, x, y, z, iblank, &ifl);
        }
        break;
    }
  }
  else if(d_type==2)
  {
    switch (P3Op.DimIs) {
      case DIMENSION_2D:
        dwrite_xyz_2d_(&id, &jd, x, y, &ifl);
        break;
      case DIMENSION_3D:
        if(P3Op.IBlankFlag == NOT_SET_IBLANK){
          dwrite_xyz_3d_(&id, &jd, &kd, x, y, z, &ifl);
        }
        else if(P3Op.IBlankFlag == SET_IBLANK){
          dwrite_xyz_3d_iblank_(&id, &jd, &kd, x, y, z, iblank, &ifl);
        }
        break;
    }
  }
}

/**
 * @brief グリッド座標ファイルの書き出し（FORMATTED）
 */
void FileIO_PLOT3D_WRITE::WriteXYZ_FORMATTED()
{
  int d_type = (sizeof(REAL_TYPE) == 4) ? 1 : 2;  // 1-float / 2-double
  if(d_type==1)
  {
    switch (P3Op.DimIs) {
      case DIMENSION_2D:
        write_xyz_2d_formatted_(&id, &jd, x, y, &ifl);
        break;
      case DIMENSION_3D:
        if(P3Op.IBlankFlag == NOT_SET_IBLANK){
          write_xyz_3d_formatted_(&id, &jd, &kd, x, y, z, &ifl);
        }
        else if(P3Op.IBlankFlag == SET_IBLANK){
          write_xyz_3d_iblank_formatted_(&id, &jd, &kd, x, y, z, iblank, &ifl);
        }
        break;
    }
  }
  else if(d_type==2)
  {
    switch (P3Op.DimIs) {
      case DIMENSION_2D:
        write_xyz_2d_formatted_(&id, &jd, x, y, &ifl);
        break;
      case DIMENSION_3D:
        if(P3Op.IBlankFlag == NOT_SET_IBLANK){
          dwrite_xyz_3d_formatted_(&id, &jd, &kd, x, y, z, &ifl);
        }
        else if(P3Op.IBlankFlag == SET_IBLANK){
          dwrite_xyz_3d_iblank_formatted_(&id, &jd, &kd, x, y, z, iblank, &ifl);
        }
        break;
    }
  }
}


/**
 * @brief グリッド座標ファイルの書き出し（UNFORMATTED）特別仕様フォーマット
 */
void FileIO_PLOT3D_WRITE::WriteXYZ_UNFORMATTED_SPECIAL()
{

}

/**
 * @brief 計算結果Qファイルの書き出し（UNFORMATTED）
 */
void FileIO_PLOT3D_WRITE::WriteQ_UNFORMATTED()
{
  int d_type = (sizeof(REAL_TYPE) == 4) ? 1 : 2;  // 1-float / 2-double
  if(d_type==1)
  {
    switch (P3Op.DimIs) {
      case DIMENSION_2D:
        write_q_2d_(&id, &jd, &fsmach, &alpha, &re, &time, q, &ifl);
        break;
      case DIMENSION_3D:
        write_q_3d_(&id, &jd, &kd, &fsmach, &alpha, &re, &time, q, &ifl);
        break;
    }
  }
  else if(d_type==2)
  {
    switch (P3Op.DimIs) {
      case DIMENSION_2D:
        dwrite_q_2d_(&id, &jd, &fsmach, &alpha, &re, &time, q, &ifl);
        break;
      case DIMENSION_3D:
        dwrite_q_3d_(&id, &jd, &kd, &fsmach, &alpha, &re, &time, q, &ifl);
        break;
    }
  }
}

/**
 * @brief 計算結果Qファイルの書き出し（FORMATTED）
 */
void FileIO_PLOT3D_WRITE::WriteQ_FORMATTED()
{
  int d_type = (sizeof(REAL_TYPE) == 4) ? 1 : 2;  // 1-float / 2-double
  if(d_type==1)
  {
    switch (P3Op.DimIs) {
      case DIMENSION_2D:
        write_q_2d_formatted_(&id, &jd, &fsmach, &alpha, &re, &time, q, &ifl);
        break;
      case DIMENSION_3D:
        write_q_3d_formatted_(&id, &jd, &kd, &fsmach, &alpha, &re, &time, q, &ifl);
        break;
    }
  }
  else if(d_type==2)
  {
    switch (P3Op.DimIs) {
      case DIMENSION_2D:
        dwrite_q_2d_formatted_(&id, &jd, &fsmach, &alpha, &re, &time, q, &ifl);
        break;
      case DIMENSION_3D:
        dwrite_q_3d_formatted_(&id, &jd, &kd, &fsmach, &alpha, &re, &time, q, &ifl);
        break;
    }
  }
}

/**
 * @brief 計算結果Qファイルの書き出し（UNFORMATTED）特別仕様フォーマット
 */
void FileIO_PLOT3D_WRITE::WriteQ_UNFORMATTED_SPECIAL()
{

}


/**
 * @brief 計算結果Functionファイルの書き出し（UNFORMATTED）
 */
void FileIO_PLOT3D_WRITE::WriteFunc_UNFORMATTED()
{
  int d_type = (sizeof(REAL_TYPE) == 4) ? 1 : 2;  // 1-float / 2-double
  if(d_type==1)
  {
    switch (P3Op.DimIs) {
      case DIMENSION_2D:
        write_func_2d_(&id, &jd, &nvar, d, &ifl);
        break;
      case DIMENSION_3D:
        write_func_3d_(&id, &jd, &kd, &nvar, d, &ifl);
        break;
    }
  }
  else if(d_type==2)
  {
    switch (P3Op.DimIs) {
      case DIMENSION_2D:
        dwrite_func_2d_(&id, &jd, &nvar, d, &ifl);
        break;
      case DIMENSION_3D:
        dwrite_func_3d_(&id, &jd, &kd, &nvar, d, &ifl);
        break;
    }
  }
}

/**
 * @brief 計算結果Functionファイルの書き出し（FORMATTED）
 */
void FileIO_PLOT3D_WRITE::WriteFunc_FORMATTED()
{
  int d_type = (sizeof(REAL_TYPE) == 4) ? 1 : 2;  // 1-float / 2-double
  if(d_type==1)
  {
    switch (P3Op.DimIs) {
      case DIMENSION_2D:
        write_func_2d_formatted_(&id, &jd, &nvar, d, &ifl);
        break;
      case DIMENSION_3D:
        write_func_3d_formatted_(&id, &jd, &kd, &nvar, d, &ifl);
        break;
    }
  }
  else if(d_type==2)
  {
    switch (P3Op.DimIs) {
      case DIMENSION_2D:
        dwrite_func_2d_formatted_(&id, &jd, &nvar, d, &ifl);
        break;
      case DIMENSION_3D:
        dwrite_func_3d_formatted_(&id, &jd, &kd, &nvar, d, &ifl);
        break;
    }
  }
}

/**
 * @brief 計算結果Functionファイルの書き出し（UNFORMATTED）特別仕様フォーマット
 */
void FileIO_PLOT3D_WRITE::WriteFunc_UNFORMATTED_SPECIAL()
{

}


/**
 * @brief グリッド座標ファイルの書き出し
 */
bool FileIO_PLOT3D_WRITE::WriteXYZData()
{
  switch (P3Op.Format) {
    case UNFORMATTED:
      WriteXYZ_UNFORMATTED();
      break;
    case FORMATTED:
      WriteXYZ_FORMATTED();
      break;
    case UNFORMATTED_SPECIAL:
      WriteXYZ_UNFORMATTED_SPECIAL();
      break;
    default:
      return false;
      break;
  }
  return true;
}

/**
 * @brief 計算結果Qファイルの書き出し
 */
bool FileIO_PLOT3D_WRITE::WriteQData()
{
  switch (P3Op.Format) {
    case UNFORMATTED:
      WriteQ_UNFORMATTED();
      break;
    case FORMATTED:
      WriteQ_FORMATTED();
      break;
    case UNFORMATTED_SPECIAL:
      WriteQ_UNFORMATTED_SPECIAL();
      break;
    default:
      return false;
      break;
  }
  return true;
}

/**
 * @brief 計算結果Functionファイルの書き出し
 */
bool FileIO_PLOT3D_WRITE::WriteFuncData()
{
  switch (P3Op.Format) {
    case UNFORMATTED:
      WriteFunc_UNFORMATTED();
      break;
    case FORMATTED:
      WriteFunc_FORMATTED();
      break;
    case UNFORMATTED_SPECIAL:
      WriteFunc_UNFORMATTED_SPECIAL();
      break;
    default:
      return false;
      break;
  }
  return true;
}

/**
 * @brief FunctionNameファイルの書き出し
 */
void FileIO_PLOT3D_WRITE::WriteFunctionName(const char* fn)
{
  funcname=fn;

  char tmp[FB_BUFF_LENGTH];
  memset(tmp, '\0', sizeof(char)*FB_BUFF_LENGTH);
  strcpy(tmp, funcname.c_str());

  int len = strlen(tmp);
  write_line_(tmp,&ifl,&len);
}


/**
 * @brief FVBNDファイルHEADERの書き出し1
 */
void FileIO_PLOT3D_WRITE::WriteFVBNDHEAD1()
{
  int len;
  string buff;
  char tmp[FB_BUFF_LENGTH];

  buff = "FVBND 1 4";
  memset(tmp, '\0', sizeof(char)*FB_BUFF_LENGTH);
  strcpy(tmp, buff.c_str());
  len = strlen(tmp);
  write_line_(tmp,&ifl,&len);
}


/**
 * @brief 文字列の書き出し
 */
void FileIO_PLOT3D_WRITE::WriteSTRING(const char* buff)
{
  string dum=buff;

  char tmp[FB_BUFF_LENGTH];
  memset(tmp, '\0', sizeof(char)*FB_BUFF_LENGTH);
  strcpy(tmp, dum.c_str());

  int len = strlen(tmp);
  write_line_(tmp,&ifl,&len);
}

/**
 * @brief FVBNDファイルHEADERの書き出し2
 */
void FileIO_PLOT3D_WRITE::WriteFVBNDHEAD2()
{
  int len;
  string buff;
  char tmp[FB_BUFF_LENGTH];

  buff = "BOUNDARIES";
  memset(tmp, '\0', sizeof(char)*FB_BUFF_LENGTH);
  strcpy(tmp, buff.c_str());
  len = strlen(tmp);
  write_line_(tmp,&ifl,&len);
}

/**
 * @brief FVBNDファイルの書き出し
 */
void FileIO_PLOT3D_WRITE::WriteFVBND(
  int type, int gridnum,
  int Imin, int Imax,
  int Jmin, int Jmax,
  int Kmin, int Kmax,
  string ResultFlag,
  int dir)
{
  string buff;
  char tmp[FB_BUFF_LENGTH];

//write boundaries
  memset(tmp, '\0', sizeof(char)*FB_BUFF_LENGTH);
  strcpy(tmp, ResultFlag.c_str());
  write_fvbnd_boundary_(
    &type,&gridnum,&Imin,&Imax,&Jmin,&Jmax,&Kmin,&Kmax,tmp,&dir,&ifl);
}

///**
// * @brief FVBNDファイルの書き出し
// */
//void FileIO_PLOT3D_WRITE::WriteFVBND(
//  const int nbname, const int nb,
//  string* boundary_name,
//  int* type, int* gridnum,
//  int* Imin,int* Imax,
//  int* Jmin,int* Jmax,
//  int* Kmin,int* Kmax,
//  string* ResultFlag,
//  int* dir)
//{
//  int len;
//  string buff;
//  char tmp[FB_BUFF_LENGTH];
//
////write header
//  buff = "FVBND 1 4";
//  memset(tmp, '\0', sizeof(char)*FB_BUFF_LENGTH);
//  strcpy(tmp, buff.c_str());
//  len = strlen(tmp);
//  write_line_(tmp,&ifl,&len);
//
////write boundary name
//  for(int i=0;i<nbname;i++){
//    buff=boundary_name[i];
//    strcpy(tmp, buff.c_str());
//    len = strlen(tmp);
//    write_line_(tmp,&ifl,&len);
//  }
//
////write boundaries
//  buff = "BOUNDARIES";
//  memset(tmp, '\0', sizeof(char)*FB_BUFF_LENGTH);
//  strcpy(tmp, buff.c_str());
//  len = strlen(tmp);
//  write_line_(tmp,&ifl,&len);
//  for(int i=0;i<nb;i++){
//    buff=ResultFlag[i];
//    memset(tmp, '\0', sizeof(char)*FB_BUFF_LENGTH);
//    strcpy(tmp, buff.c_str());
//    len = strlen(tmp);
//    write_fvbnd_boundary_(
//      &type[i],&gridnum[i],&Imin[i],&Imax[i],&Jmin[i],&Jmax[i],&Kmin[i],&Kmax[i],tmp,&dir[i],&ifl);
//  }
//}
