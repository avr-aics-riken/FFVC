#ifndef _PLT3D_PLOT3D_WRITE_H_
#define _PLT3D_PLOT3D_WRITE_H_

//##################################################################################
//
// FFV-C : Frontflow / violet Cartesian
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2013 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2013 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

/** 
 * @file   PLOT3D_write.h
 * @brief  FlowBase FileIO_PLOT3D_WRITE class Header
 * @author kero
 */

#include "PLOT3D.h"


#ifdef _WIN32

// PLT3D_write.f90
#define write_ngrid_data_                     WRITE_NGIRD_DATA
#define write_ngrid_data_formatted_           WRITE_NGIRD_DATA_FORMATTED
#define write_block_data_                     WRITE_BLOCK_DATA
#define write_block_data_formatted_           WRITE_BLOCK_DATA_FORMATTED
#define write_block_data_2d_                  WRITE_BLOCK_DATA_2D
#define write_block_data_2d_formatted_        WRITE_BLOCK_DATA_2D_FORMATTED
#define write_func_block_data_                WRITE_FUNC_BLOCK_DATA
#define write_func_block_data_formatted_      WRITE_FUNC_BLOCK_DATA_FORMATTED
#define write_func_block_data_2d_             WRITE_FUNC_BLOCK_DATA_2D
#define write_func_block_data_2d_formatted_   WRITE_FUNC_BLOCK_DATA_2D_FORMATTED
#define write_xyz_3d_                         WRITE_XYZ_3D
#define write_xyz_3d_formatted_               WRITE_XYZ_3D_FORMATTED
#define write_xyz_3d_iblank_                  WRITE_XYZ_3D_IBLANK
#define write_xyz_3d_iblank_formatted_        WRITE_XYZ_3D_IBLANK_FORMATTED
#define write_xyz_2d_                         WRITE_XYZ_2D
#define write_xyz_2d_formatted_               WRITE_XYZ_2D_FORMATTED
#define write_q_3d_                           WRITE_Q_3D
#define write_q_3d_formatted_                 WRITE_Q_3D_FORMATTED
#define write_q_2d_                           WRITE_Q_2D
#define write_q_2d_formatted_                 WRITE_Q_2D_FORMATTED
#define write_func_3d_                        WRITE_FUNC_3D
#define write_func_3d_formatted_              WRITE_FUNC_3D_FORMATTED
#define write_func_2d_                        WRITE_FUNC_2D
#define write_func_2d_formatted_              WRITE_FUNC_2D_FORMATTED
#define dwrite_xyz_3d_                        DWRITE_XYZ_3D
#define dwrite_xyz_3d_formatted_              DWRITE_XYZ_3D_FORMATTED
#define dwrite_xyz_3d_iblank_                 DWRITE_XYZ_3D_IBLANK
#define dwrite_xyz_3d_iblank_formatted_       DWRITE_XYZ_3D_IBLANK_FORMATTED
#define dwrite_xyz_2d_                        DWRITE_XYZ_2D
#define dwrite_xyz_2d_formatted_              DWRITE_XYZ_2D_FORMATTED
#define dwrite_q_3d_                          DWRITE_Q_3D
#define dwrite_q_3d_formatted_                DWRITE_Q_3D_FORMATTED
#define dwrite_q_2d_                          DWRITE_Q_2D
#define dwrite_q_2d_formatted_                DWRITE_Q_2D_FORMATTED
#define dwrite_func_3d_                       DWRITE_FUNC_3D
#define dwrite_func_3d_formatted_             DWRITE_FUNC_3D_FORMATTED
#define dwrite_func_2d_                       DWRITE_FUNC_2D
#define dwrite_func_2d_formatted_             DWRITE_FUNC_2D_FORMATTED

#endif // _WIN32


extern "C" {

// PLT3D_write.f90
  void write_ngrid_data_ (int* ngrid, int* ifl);
  void write_ngrid_data_formatted_  (int* ngrid, int* ifl);
  void write_block_data_ (
      int* id, int* jd, int* kd, int* ifl);
  void write_block_data_formatted_  (
      int* id, int* jd, int* kd, int* ifl);
  void write_block_data_2d_ (
      int* id, int* jd, int* ifl);
  void write_block_data_2d_formatted_  (
      int* id, int* jd, int* ifl);
  void write_func_block_data_ (
                               int* id, int* jd, int* kd, int* nvar, int* ngrid, int* ifl);
  void write_func_block_data_formatted_  (
                                          int* id, int* jd, int* kd, int* nvar, int* ngrid, int* ifl);
  void write_func_block_data_2d_ (
                                  int* id, int* jd, int* nvar, int* ngrid, int* ifl);
  void write_func_block_data_2d_formatted_  (
                                             int* id, int* jd, int* nvar, int* ngrid, int* ifl);
  void write_xyz_3d_ (
      int* id, int* jd, int* kd, float* x, float* y, float* z, int* ifl);
  void write_xyz_3d_formatted_ (
      int* id, int* jd, int* kd, float* x, float* y, float* z, int* ifl);
  void write_xyz_3d_iblank_ (
      int* id, int* jd, int* kd, float* x, float* y, float* z, int* iblank, int* ifl);
  void write_xyz_3d_iblank_formatted_ (
      int* id, int* jd, int* kd, float* x, float* y, float* z, int* iblank, int* ifl);
  void write_xyz_2d_ (
      int* id, int* jd, float* x, float* y, int* ifl);
  void write_xyz_2d_formatted_  (
      int* id, int* jd, float* x, float* y, int* ifl);
  void write_q_3d_ (
      int* id, int* jd, int* kd, float* fsmach, float* alpha, float* re, float* time, float* q, int* ifl);
  void write_q_3d_formatted_ (
      int* id, int* jd, int* kd, float* fsmach, float* alpha, float* re, float* time, float* q, int* ifl);
  void write_q_2d_ (
      int* id, int* jd, float* fsmach, float* alpha, float* re, float* time, float* q, int* ifl);
  void write_q_2d_formatted_ (
      int* id, int* jd, float* fsmach, float* alpha, float* re, float* time, float* q, int* ifl);
  void write_func_3d_ (
      int* id, int* jd, int* kd, int* nvar, float* d, int* ifl);
  void write_func_3d_formatted_ (
      int* id, int* jd, int* kd, int* nvar, float* d, int* ifl);
  void write_func_2d_ (
      int* id, int* jd, int* nvar, float* d, int* ifl);
  void write_func_2d_formatted_ (
      int* id, int* jd, int* nvar, float* d, int* ifl);
  void dwrite_xyz_3d_ (
      int* id, int* jd, int* kd, double* x, double* y, double* z, int* ifl);
  void dwrite_xyz_3d_formatted_ (
      int* id, int* jd, int* kd, double* x, double* y, double* z, int* ifl);
  void dwrite_xyz_3d_iblank_ (
      int* id, int* jd, int* kd, double* x, double* y, double* z, int* iblank, int* ifl);
  void dwrite_xyz_3d_iblank_formatted_ (
      int* id, int* jd, int* kd, double* x, double* y, double* z, int* iblank, int* ifl);
  void dwrite_xyz_2d_ (
      int* id, int* jd, double* x, double* y, int* ifl);
  void dwrite_xyz_2d_formatted_  (
      int* id, int* jd, double* x, double* y, int* ifl);
  void dwrite_q_3d_ (
      int* id, int* jd, int* kd, double* fsmach, double* alpha, double* re, double* time, double* q, int* ifl);
  void dwrite_q_3d_formatted_ (
      int* id, int* jd, int* kd, double* fsmach, double* alpha, double* re, double* time, double* q, int* ifl);
  void dwrite_q_2d_ (
      int* id, int* jd, double* fsmach, double* alpha, double* re, double* time, double* q, int* ifl);
  void dwrite_q_2d_formatted_ (
      int* id, int* jd, double* fsmach, double* alpha, double* re, double* time, double* q, int* ifl);
  void dwrite_func_3d_ (
      int* id, int* jd, int* kd, int* nvar, double* d, int* ifl);
  void dwrite_func_3d_formatted_ (
      int* id, int* jd, int* kd, int* nvar, double* d, int* ifl);
  void dwrite_func_2d_ (
      int* id, int* jd, int* nvar, double* d, int* ifl);
  void dwrite_func_2d_formatted_ (
      int* id, int* jd, int* nvar, double* d, int* ifl);

}

using namespace std;

class FileIO_PLOT3D_WRITE : public FileIO_PLOT3D {
  
public:
  /** コンストラクタ */
  FileIO_PLOT3D_WRITE() {}

  /**　デストラクタ */
  ~FileIO_PLOT3D_WRITE() {}
  
private:

protected:

  /**
   * @brief グリッド座標ファイルの書き出し（UNFORMATTED）
   */
  void WriteXYZ_UNFORMATTED();

  /**
   * @brief グリッド座標ファイルの書き出し（FORMATTED）
   */
  void WriteXYZ_FORMATTED();


  /**
   * @brief グリッド座標ファイルの書き出し（C_BINARY）
   */
  void WriteXYZ_C_BINARY();

  /**
   * @brief 計算結果Qファイルの書き出し（UNFORMATTED）
   */
  void WriteQ_UNFORMATTED();

  /**
   * @brief 計算結果Qファイルの書き出し（FORMATTED）
   */
  void WriteQ_FORMATTED();


  /**
   * @brief 計算結果Qファイルの書き出し（C_BINARY）
   */
  void WriteQ_C_BINARY();

  /**
   * @brief 計算結果Functionファイルの書き出し（UNFORMATTED）
   * @note  出力項目すべてを一度に書き出し ---> nvar:出力項目数
   */
  void WriteFunc_UNFORMATTED();

  /**
   * @brief 計算結果Functionファイルの書き出し（FORMATTED）
   * @note  出力項目すべてを一度に書き出し ---> nvar:出力項目数
   */
  void WriteFunc_FORMATTED();

  /**
   * @brief 計算結果Functionファイルの書き出し（C_BINARY）
   * @note  C_BINARYでの出力は項目ごとに書き出し
   */
  void WriteFunc_C_BINARY();

public:

  /**
   * @brief PLOT3Dの出力ファイルオープン
   */
  bool OpenFile();

  /**
   * @brief グリッド数の書き出し
   */
  void WriteNgrid(int ngrid);

  /**
   * @brief ブロックデータの書き出し
   */
  void WriteBlockData(int id, int jd, int kd);

  /**
   * @brief ブロックデータの書き出し（2D）
   */
  void WriteBlockData2D(int id, int jd);

  /**
   * @brief 計算結果Functionファイル、ブロックデータの書き出し
   */
  void WriteFuncBlockData(int* id, int* jd, int* kd, int* nvar, int ngrid);

  /**
   * @brief グリッド座標ファイルの書き出し
   */
  bool WriteXYZData();

  /**
   * @brief 計算結果Qファイルの書き出し
   */
  bool WriteQData();

  /**
   * @brief 計算結果Functionファイルの書き出し
   */
  bool WriteFuncData();

  /**
   * @brief FunctionNameファイルの書き出し
   */
  void WriteFunctionName(const char* fn);

  /**
   * @brief FVBNDファイルHEADERの書き出し1
   */
  void WriteFVBNDHEAD1();

  /**
   * @brief 文字列の書き出し
   */
  void WriteSTRING(const char* buff);

  /**
   * @brief FVBNDファイルHEADERの書き出し2
   */
  void WriteFVBNDHEAD2();

  /**
   * @brief FVBNDファイルの書き出し
   */
  void WriteFVBND(
    int type, int gridnum,
    int Imin, int Imax, int Jmin, int Jmax,
    int Kmin, int Kmax, string ResultFlag,int dir);

  /**
   * @brief C_BINARY用ヘッダorフッタの書き出し
   */
  void WriteMarker(int dmy);

};

#endif // _PLT3D_PLOT3D_WRITE_H_
