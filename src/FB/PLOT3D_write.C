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

      //open_plot3d_outputfile_(&flag,tmp,&ifl,&fnsize,&ierror);
      open_plot3d_outputfile_(&flag,tmp,&ifl,&len,&ierror);
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
    case C_BINARY:
      fwrite(&ngrid, sizeof(int), 1, fp);
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
    case C_BINARY:
      fwrite(&id, sizeof(int), 1, fp);
      fwrite(&jd, sizeof(int), 1, fp);
      fwrite(&kd, sizeof(int), 1, fp);
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
    case C_BINARY:
      fwrite(&id, sizeof(int), 1, fp);
      fwrite(&jd, sizeof(int), 1, fp);
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
        case C_BINARY:
          for(int igrid=0;igrid<ngrid;igrid++){
            fwrite(&id[igrid], sizeof(int), 1, fp);
            fwrite(&jd[igrid], sizeof(int), 1, fp);
            fwrite(&nvar[igrid], sizeof(int), 1, fp);
          }
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
        case C_BINARY:
          for(int igrid=0;igrid<ngrid;igrid++){
            fwrite(&id[igrid], sizeof(int), 1, fp);
            fwrite(&jd[igrid], sizeof(int), 1, fp);
            fwrite(&kd[igrid], sizeof(int), 1, fp);
            fwrite(&nvar[igrid], sizeof(int), 1, fp);
          }
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
  //int d_type = (sizeof(REAL_TYPE) == 4) ? 1 : 2;  // 1-float / 2-double
  //if(d_type==1)
  if(P3Op.realtype==1)
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
  //else if(d_type==2)
  else if(P3Op.realtype==2)
  {
    switch (P3Op.DimIs) {
      case DIMENSION_2D:
        dwrite_xyz_2d_(&id, &jd, dx, dy, &ifl);
        break;
      case DIMENSION_3D:
        if(P3Op.IBlankFlag == NOT_SET_IBLANK){
          dwrite_xyz_3d_(&id, &jd, &kd, dx, dy, dz, &ifl);
        }
        else if(P3Op.IBlankFlag == SET_IBLANK){
          dwrite_xyz_3d_iblank_(&id, &jd, &kd, dx, dy, dz, iblank, &ifl);
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
  //int d_type = (sizeof(REAL_TYPE) == 4) ? 1 : 2;  // 1-float / 2-double
  //if(d_type==1)
  if(P3Op.realtype==1)
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
  //else if(d_type==2)
  else if(P3Op.realtype==2)
  {
    switch (P3Op.DimIs) {
      case DIMENSION_2D:
        dwrite_xyz_2d_formatted_(&id, &jd, dx, dy, &ifl);
        break;
      case DIMENSION_3D:
        if(P3Op.IBlankFlag == NOT_SET_IBLANK){
          dwrite_xyz_3d_formatted_(&id, &jd, &kd, dx, dy, dz, &ifl);
        }
        else if(P3Op.IBlankFlag == SET_IBLANK){
          dwrite_xyz_3d_iblank_formatted_(&id, &jd, &kd, dx, dy, dz, iblank, &ifl);
        }
        break;
    }
  }
}



/**
 * @brief グリッド座標ファイルの書き出し（C_BINARY）
 */
void FileIO_PLOT3D_WRITE::WriteXYZ_C_BINARY()
{
  //int d_type = (sizeof(REAL_TYPE) == 4) ? 1 : 2;  // 1-float / 2-double
  //if(d_type==1)
  if(P3Op.realtype==1)
  {
    switch (P3Op.DimIs) {
      case DIMENSION_2D:
        fwrite(x, sizeof(float), id*jd, fp);
        fwrite(y, sizeof(float), id*jd, fp);
        break;
      case DIMENSION_3D:
        fwrite(x, sizeof(float), id*jd*kd, fp);
        fwrite(y, sizeof(float), id*jd*kd, fp);
        fwrite(z, sizeof(float), id*jd*kd, fp);
        if(P3Op.IBlankFlag == SET_IBLANK){
          fwrite(iblank, sizeof(int), id*jd*kd, fp);
        }
        break;
    }
  }
  //else if(d_type==2)
  else if(P3Op.realtype==2)
  {
    switch (P3Op.DimIs) {
      case DIMENSION_2D:
        fwrite(dx, sizeof(double), id*jd, fp);
        fwrite(dy, sizeof(double), id*jd, fp);
        break;
      case DIMENSION_3D:
        fwrite(dx, sizeof(double), id*jd*kd, fp);
        fwrite(dy, sizeof(double), id*jd*kd, fp);
        fwrite(dz, sizeof(double), id*jd*kd, fp);
        if(P3Op.IBlankFlag == SET_IBLANK){
          fwrite(iblank, sizeof(int), id*jd*kd, fp);
        }
        break;
    }
  }
}


/**
 * @brief 計算結果Qファイルの書き出し（UNFORMATTED）
 */
void FileIO_PLOT3D_WRITE::WriteQ_UNFORMATTED()
{
  //int d_type = (sizeof(REAL_TYPE) == 4) ? 1 : 2;  // 1-float / 2-double
  //if(d_type==1)
  if(P3Op.realtype==1)
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
  //else if(d_type==2)
  else if(P3Op.realtype==2)
  {
    switch (P3Op.DimIs) {
      case DIMENSION_2D:
        dwrite_q_2d_(&id, &jd, &dfsmach, &dalpha, &dre, &dtime, dq, &ifl);
        break;
      case DIMENSION_3D:
        dwrite_q_3d_(&id, &jd, &kd, &dfsmach, &dalpha, &dre, &dtime, dq, &ifl);
        break;
    }
  }
}

/**
 * @brief 計算結果Qファイルの書き出し（FORMATTED）
 */
void FileIO_PLOT3D_WRITE::WriteQ_FORMATTED()
{
  //int d_type = (sizeof(REAL_TYPE) == 4) ? 1 : 2;  // 1-float / 2-double
  //if(d_type==1)
  if(P3Op.realtype==1)
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
  //else if(d_type==2)
  else if(P3Op.realtype==2)
  {
    switch (P3Op.DimIs) {
      case DIMENSION_2D:
        dwrite_q_2d_formatted_(&id, &jd, &dfsmach, &dalpha, &dre, &dtime, dq, &ifl);
        break;
      case DIMENSION_3D:
        dwrite_q_3d_formatted_(&id, &jd, &kd, &dfsmach, &dalpha, &dre, &dtime, dq, &ifl);
        break;
    }
  }
}


/**
 * @brief 計算結果Qファイルの書き出し（C_BINARY）
 */
void FileIO_PLOT3D_WRITE::WriteQ_C_BINARY()
{

}

/**
 * @brief 計算結果Functionファイルの書き出し（UNFORMATTED）
 * @note  出力項目すべてを一度に書き出し ---> nvar:出力項目数
 */
void FileIO_PLOT3D_WRITE::WriteFunc_UNFORMATTED()
{
  //int d_type = (sizeof(REAL_TYPE) == 4) ? 1 : 2;  // 1-float / 2-double
  //if(d_type==1)
  if(P3Op.realtype==1)
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
  //else if(d_type==2)
  else if(P3Op.realtype==2)
  {
    switch (P3Op.DimIs) {
      case DIMENSION_2D:
        dwrite_func_2d_(&id, &jd, &nvar, dd, &ifl);
        break;
      case DIMENSION_3D:
        dwrite_func_3d_(&id, &jd, &kd, &nvar, dd, &ifl);
        break;
    }
  }
}

/**
 * @brief 計算結果Functionファイルの書き出し（FORMATTED）
 * @note  出力項目すべてを一度に書き出し ---> nvar:出力項目数
 */
void FileIO_PLOT3D_WRITE::WriteFunc_FORMATTED()
{
  //int d_type = (sizeof(REAL_TYPE) == 4) ? 1 : 2;  // 1-float / 2-double
  //if(d_type==1)
  if(P3Op.realtype==1)
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
  //else if(d_type==2)
  else if(P3Op.realtype==2)
  {
    switch (P3Op.DimIs) {
      case DIMENSION_2D:
        dwrite_func_2d_formatted_(&id, &jd, &nvar, dd, &ifl);
        break;
      case DIMENSION_3D:
        dwrite_func_3d_formatted_(&id, &jd, &kd, &nvar, dd, &ifl);
        break;
    }
  }
}


/**
 * @brief 計算結果Functionファイルの書き出し（C_BINARY）
 * @note  C_BINARYでの出力は項目ごとに書き出し
 */
void FileIO_PLOT3D_WRITE::WriteFunc_C_BINARY()
{

#if 1 //項目ごと

  //int d_type = (sizeof(REAL_TYPE) == 4) ? 1 : 2;  // 1-float / 2-double
  //if(d_type==1)
  if(P3Op.realtype==1)
  {
    switch (P3Op.DimIs) {
      case DIMENSION_2D:
        fwrite(d, sizeof(float), id*jd, fp);
        break;
      case DIMENSION_3D:
        fwrite(d, sizeof(float), id*jd*kd, fp);
        break;
    }
  }
  //else if(d_type==2)
  else if(P3Op.realtype==2)
  {
    switch (P3Op.DimIs) {
      case DIMENSION_2D:
        fwrite(dd, sizeof(double), id*jd, fp);
        break;
      case DIMENSION_3D:
        fwrite(dd, sizeof(double), id*jd*kd, fp);
        break;
    }
  }

#else //全項目

  //int d_type = (sizeof(REAL_TYPE) == 4) ? 1 : 2;  // 1-float / 2-double
  //if(d_type==1)
  if(P3Op.realtype==1)
  {
    switch (P3Op.DimIs) {
      case DIMENSION_2D:
        fwrite(d, sizeof(float), id*jd*nvar, fp);
        break;
      case DIMENSION_3D:
        fwrite(d, sizeof(float), id*jd*kd*nvar, fp);
        break;
    }
  }
  //else if(d_type==2)
  else if(P3Op.realtype==2)
  {
    switch (P3Op.DimIs) {
      case DIMENSION_2D:
        fwrite(dd, sizeof(double), id*jd*nvar, fp);
        break;
      case DIMENSION_3D:
        fwrite(dd, sizeof(double), id*jd*kd*nvar, fp);
        break;
    }
  }

#endif

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
    case C_BINARY:
      WriteXYZ_C_BINARY();
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
    case C_BINARY:
      WriteQ_C_BINARY();
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
    case C_BINARY:
      WriteFunc_C_BINARY();
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

/**
 * @brief C_BINARY用ヘッダorフッタの書き出し
 */
void FileIO_PLOT3D_WRITE::WriteMarker(int dmy)
{
  fwrite(&dmy, sizeof(int), 1, fp);
}