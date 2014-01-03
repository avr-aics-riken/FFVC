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
  char tmp[FB_FILE_PATH_LENGTH];
  int len;
  int flag;
  int fnsize;
  
  int ierror=0;
  switch (P3Op.Format) {
    case UNFORMATTED:
    case FORMATTED:
      
      memset(tmp, 0, sizeof(char)*FB_FILE_PATH_LENGTH);
      strcpy(tmp, fname.c_str());
      
      len=strlen(tmp);
      flag=P3Op.Format;
      fnsize=FB_FILE_PATH_LENGTH;
      
      open_plot3d_inputfile_(&flag,tmp,&ifl,&len,&ierror);
      break;
      
    case C_BINARY:
      if( (fp = fopen(fname.c_str(), "rb")) == NULL ) {
        fprintf(stderr, "Can't open file.(%s)\n", fname.c_str());
        return false;
      }
      ierror = true;
      break;
  }
  
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
    case C_BINARY:
      fread(ngrid, sizeof(int), 1, fp);
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
    case C_BINARY:
      fread(id, sizeof(int), 1, fp);
      fread(jd, sizeof(int), 1, fp);
      fread(kd, sizeof(int), 1, fp);
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
    case C_BINARY:
      fread(id, sizeof(int), 1, fp);
      fread(jd, sizeof(int), 1, fp);
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
        case C_BINARY:
          for(int igrid=0;igrid<ngrid;igrid++){
            fread(id, sizeof(int), 1, fp);
            fread(jd, sizeof(int), 1, fp);
            fread(nvar, sizeof(int), 1, fp);
          }
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
        case C_BINARY:
          for(int igrid=0;igrid<ngrid;igrid++){
            fread(id, sizeof(int), 1, fp);
            fread(jd, sizeof(int), 1, fp);
            fread(kd, sizeof(int), 1, fp);
            fread(nvar, sizeof(int), 1, fp);
          }
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
  if(P3Op.realtype==1)
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
  else if(P3Op.realtype==2)
  {
    switch (P3Op.DimIs) {
      case DIMENSION_2D:
        dread_xyz_2d_(&id, &jd, dx, dy, &ifl);
        break;
      case DIMENSION_3D:
        if(P3Op.IBlankFlag == NOT_SET_IBLANK){
          dread_xyz_3d_(&id, &jd, &kd, dx, dy, dz, &ifl);
        }
        else if(P3Op.IBlankFlag == SET_IBLANK){
          dread_xyz_3d_iblank_(&id, &jd, &kd, dx, dy, dz, iblank, &ifl);
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
  if(P3Op.realtype==1)
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
  else if(P3Op.realtype==2)
  {
    switch (P3Op.DimIs) {
      case DIMENSION_2D:
        dread_xyz_2d_formatted_(&id, &jd, dx, dy, &ifl);
        break;
      case DIMENSION_3D:
        if(P3Op.IBlankFlag == NOT_SET_IBLANK){
          dread_xyz_3d_formatted_(&id, &jd, &kd, dx, dy, dz, &ifl);
        }
        else if(P3Op.IBlankFlag == SET_IBLANK){
          dread_xyz_3d_iblank_formatted_(&id, &jd, &kd, dx, dy, dz, iblank, &ifl);
        }
        break;
    }
  }
}


/**
 * @brief グリッド座標ファイルの読み込み（C_BINARY）
 */
void FileIO_PLOT3D_READ::ReadXYZ_C_BINARY()
{
  size_t s1 = (size_t)id * (size_t)jd;
  size_t s2 = s1 * (size_t)kd;
  
  if(P3Op.realtype==1)
  {
    switch (P3Op.DimIs) {
      case DIMENSION_2D:
        fread(x, sizeof(float), s1, fp);
        fread(y, sizeof(float), s1, fp);
        break;
      case DIMENSION_3D:
        fread(x, sizeof(float), s2, fp);
        fread(y, sizeof(float), s2, fp);
        fread(z, sizeof(float), s2, fp);
        if(P3Op.IBlankFlag == SET_IBLANK){
          fread(iblank, sizeof(int), s2, fp);
        }
        break;
    }
  }
  else if(P3Op.realtype==2)
  {
    switch (P3Op.DimIs) {
      case DIMENSION_2D:
        fread(dx, sizeof(double), s1, fp);
        fread(dy, sizeof(double), s1, fp);
        break;
      case DIMENSION_3D:
        fread(dx, sizeof(double), s2, fp);
        fread(dy, sizeof(double), s2, fp);
        fread(dz, sizeof(double), s2, fp);
        if(P3Op.IBlankFlag == SET_IBLANK){
          fread(iblank, sizeof(int), s2, fp);
        }
        break;
    }
  }
}

/**
 * @brief 計算結果Qファイルの読み込み（UNFORMATTED）
 */
void FileIO_PLOT3D_READ::ReadQ_UNFORMATTED()
{
  
  if(P3Op.realtype==1)
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
  else if(P3Op.realtype==2)
  {
    switch (P3Op.DimIs) {
      case DIMENSION_2D:
        dread_q_2d_(&id, &jd, &dfsmach, &dalpha, &dre, &dtime, dq, &ifl);
        break;
      case DIMENSION_3D:
        dread_q_3d_(&id, &jd, &kd, &dfsmach, &dalpha, &dre, &dtime, dq, &ifl);
        break;
    }
  }
}

/**
 * @brief 計算結果Qファイルの読み込み（FORMATTED）
 */
void FileIO_PLOT3D_READ::ReadQ_FORMATTED()
{
  if(P3Op.realtype==1)
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
  else if(P3Op.realtype==2)
  {
    switch (P3Op.DimIs) {
      case DIMENSION_2D:
        dread_q_2d_formatted_(&id, &jd, &dfsmach, &dalpha, &dre, &dtime, dq, &ifl);
        break;
      case DIMENSION_3D:
        dread_q_3d_formatted_(&id, &jd, &kd, &dfsmach, &dalpha, &dre, &dtime, dq, &ifl);
        break;
    }
  }
}


/**
 * @brief 計算結果Qファイルの読み込み（C_BINARY）
 */
void FileIO_PLOT3D_READ::ReadQ_C_BINARY()
{
  
}

/**
 * @brief 計算結果Functionファイルの読み込み（UNFORMATTED）
 */
void FileIO_PLOT3D_READ::ReadFunc_UNFORMATTED()
{
  if(P3Op.realtype==1)
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
  else if(P3Op.realtype==2)
  {
    switch (P3Op.DimIs) {
      case DIMENSION_2D:
        dread_func_2d_(&id, &jd, &nvar, dd, &ifl);
        break;
      case DIMENSION_3D:
        dread_func_3d_(&id, &jd, &kd, &nvar, dd, &ifl);
        break;
    }
  }
}

/**
 * @brief 計算結果Functionファイルの読み込み（FORMATTED）
 */
void FileIO_PLOT3D_READ::ReadFunc_FORMATTED()
{
  if(P3Op.realtype==1)
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
  else if(P3Op.realtype==2)
  {
    switch (P3Op.DimIs) {
      case DIMENSION_2D:
        dread_func_2d_formatted_(&id, &jd, &nvar, dd, &ifl);
        break;
      case DIMENSION_3D:
        dread_func_3d_formatted_(&id, &jd, &kd, &nvar, dd, &ifl);
        break;
    }
  }
}


/**
 * @brief 計算結果Functionファイルの読み込み（C_BINARY）
 */
void FileIO_PLOT3D_READ::ReadFunc_C_BINARY()
{
  size_t s1 = (size_t)id * (size_t)jd;
  size_t s2 = s1 * (size_t)kd;
  size_t s3 = s1 * (size_t)nvar;
  size_t s4 = s2 * (size_t)nvar;
  
  if(P3Op.realtype==1)
  {
    switch (P3Op.DimIs) {
      case DIMENSION_2D:
        fread(d, sizeof(float), s1, fp);
        break;
      case DIMENSION_3D:
        fread(d, sizeof(float), s2, fp);
        break;
    }
  }
  else if(P3Op.realtype==2)
  {
    switch (P3Op.DimIs) {
      case DIMENSION_2D:
        fread(dd, sizeof(double), s1, fp);
        break;
      case DIMENSION_3D:
        fread(dd, sizeof(double), s2, fp);
        break;
    }
  }
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
    case C_BINARY:
      ReadXYZ_C_BINARY();
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
                                   float* m_fsmach,
                                   float* m_alpha,
                                   float* m_re,
                                   float* m_time,
                                   float* m_q)
{
  q=m_q;
  
  switch (P3Op.Format) {
    case UNFORMATTED:
      ReadQ_UNFORMATTED();
      break;
    case FORMATTED:
      ReadQ_FORMATTED();
      break;
    case C_BINARY:
      ReadQ_C_BINARY();
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
 * @brief 計算結果Qファイルの読み込み
 */
bool FileIO_PLOT3D_READ::ReadQData(
                                   double* m_fsmach,
                                   double* m_alpha,
                                   double* m_re,
                                   double* m_time,
                                   double* m_q)
{
  dq=m_q;
  
  switch (P3Op.Format) {
    case UNFORMATTED:
      ReadQ_UNFORMATTED();
      break;
    case FORMATTED:
      ReadQ_FORMATTED();
      break;
    case C_BINARY:
      ReadQ_C_BINARY();
      break;
    default:
      return false;
      break;
  }
  
  *m_fsmach=dfsmach;
  *m_alpha=dalpha;
  *m_re=dre;
  *m_time=dtime;
  
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
    case C_BINARY:
      ReadFunc_C_BINARY();
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
                                   int* m_Imin, int* m_Imax,
                                   int* m_Jmin, int* m_Jmax,
                                   int* m_Kmin, int* m_Kmax,
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
