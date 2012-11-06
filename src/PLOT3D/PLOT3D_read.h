#ifndef _PLT3D_PLOT3D_READ_H_
#define _PLT3D_PLOT3D_READ_H_

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
 * @file   PLOT3D_read.h
 * @brief  FlowBase FileIO_PLOT3D_READ class Header
 * @author kero
 */

#include "PLOT3D.h"

#ifdef _WIN32

// PLT3D_read.f90
#define read_ngrid_data_                      READ_NGIRD_DATA
#define read_ngrid_data_formatted_            READ_NGIRD_DATA_FORMATTED
#define read_block_data_                      READ_BLOCK_DATA
#define read_block_data_formatted_            READ_BLOCK_DATA_FORMATTED
#define read_block_data_2d_                   READ_BLOCK_DATA_2D
#define read_block_data_2d_formatted_         READ_BLOCK_DATA_2D_FORMATTED
#define read_func_block_data_                 READ_FUNC_BLOCK_DATA
#define read_func_block_data_formatted_       READ_FUNC_BLOCK_DATA_FORMATTED
#define read_func_block_data_2d_              READ_FUNC_BLOCK_DATA_2D
#define read_func_block_data_2d_formatted_    READ_FUNC_BLOCK_DATA_2D_FORMATTED
#define read_xyz_3d_                          READ_XYZ_3D
#define read_xyz_3d_formatted_                READ_XYZ_3D_FORMATTED
#define read_xyz_3d_iblank_                   READ_XYZ_3D_IBLANK
#define read_xyz_3d_iblank_formatted_         READ_XYZ_3D_IBLANK_FORMATTED
#define read_xyz_2d_                          READ_XYZ_2D
#define read_xyz_2d_formatted_                READ_XYZ_2D_FORMATTED
#define read_q_3d_                            READ_Q_3D
#define read_q_3d_formatted_                  READ_Q_3D_FORMATTED
#define read_q_2d_                            READ_Q_2D
#define read_q_2d_formatted_                  READ_Q_2D_FORMATTED
#define read_func_3d_                         READ_FUNC_3D
#define read_func_3d_formatted_               READ_FUNC_3D_FORMATTED
#define read_func_2d_                         READ_FUNC_2D
#define read_func_2d_formatted_               READ_FUNC_2D_FORMATTED
#define dread_xyz_3d_                         DREAD_XYZ_3D
#define dread_xyz_3d_formatted_               DREAD_XYZ_3D_FORMATTED
#define dread_xyz_3d_iblank_                  DREAD_XYZ_3D_IBLANK
#define dread_xyz_3d_iblank_formatted_        DREAD_XYZ_3D_IBLANK_FORMATTED
#define dread_xyz_2d_                         DREAD_XYZ_2D
#define dread_xyz_2d_formatted_               DREAD_XYZ_2D_FORMATTED
#define dread_q_3d_                           DREAD_Q_3D
#define dread_q_3d_formatted_                 DREAD_Q_3D_FORMATTED
#define dread_q_2d_                           DREAD_Q_2D
#define dread_q_2d_formatted_                 DREAD_Q_2D_FORMATTED
#define dread_func_3d_                        DREAD_FUNC_3D
#define dread_func_3d_formatted_              DREAD_FUNC_3D_FORMATTED
#define dread_func_2d_                        DREAD_FUNC_2D
#define dread_func_2d_formatted_              DREAD_FUNC_2D_FORMATTED

#endif // _WIN32


extern "C" {

// PLT3D_read.f90
  void read_ngrid_data_ (int* ngrid, int* ifl);
  void read_ngrid_data_formatted_  (int* ngrid, int* ifl);
  void read_block_data_ (
      int* id, int* jd, int* kd, int* ifl);
  void read_block_data_formatted_  (
      int* id, int* jd, int* kd, int* ifl);
  void read_block_data_2d_ (
      int* id, int* jd, int* ifl);
  void read_block_data_2d_formatted_  (
      int* id, int* jd, int* ifl);
  void read_func_block_data_ (
      int* id, int* jd, int* kd, int* nvar, int* ngrid, int* ifl);
  void read_func_block_data_formatted_  (
      int* id, int* jd, int* kd, int* nvar, int* ngrid, int* ifl);
  void read_func_block_data_2d_ (
      int* id, int* jd, int* nvar, int* ngrid, int* ifl);
  void read_func_block_data_2d_formatted_  (
      int* id, int* jd, int* nvar, int* ngrid, int* ifl);
  void read_xyz_3d_ (
      int* id, int* jd, int* kd, float* x, float* y, float* z, int* ifl);
  void read_xyz_3d_formatted_ (
      int* id, int* jd, int* kd, float* x, float* y, float* z, int* ifl);
  void read_xyz_3d_iblank_ (
      int* id, int* jd, int* kd, float* x, float* y, float* z, int* iblank, int* ifl);
  void read_xyz_3d_iblank_formatted_ (
      int* id, int* jd, int* kd, float* x, float* y, float* z, int* iblank, int* ifl);
  void read_xyz_2d_ (
      int* id, int* jd, float* x, float* y, int* ifl);
  void read_xyz_2d_formatted_  (
      int* id, int* jd, float* x, float* y, int* ifl);
  void read_q_3d_ (
      int* id, int* jd, int* kd, float* fsmach, float* alpha, float* re, float* time, float* q, int* ifl);
  void read_q_3d_formatted_ (
      int* id, int* jd, int* kd, float* fsmach, float* alpha, float* re, float* time, float* q, int* ifl);
  void read_q_2d_ (
      int* id, int* jd, float* fsmach, float* alpha, float* re, float* time, float* q, int* ifl);
  void read_q_2d_formatted_ (
      int* id, int* jd, float* fsmach, float* alpha, float* re, float* time, float* q, int* ifl);
  void read_func_3d_ (
      int* id, int* jd, int* kd, int* nvar, float* d, int* ifl);
  void read_func_3d_formatted_ (
      int* id, int* jd, int* kd, int* nvar, float* d, int* ifl);
  void read_func_2d_ (
      int* id, int* jd, int* nvar, float* d, int* ifl);
  void read_func_2d_formatted_ (
      int* id, int* jd, int* nvar, float* d, int* ifl);
  void dread_xyz_3d_ (
      int* id, int* jd, int* kd, double* x, double* y, double* z, int* ifl);
  void dread_xyz_3d_formatted_ (
      int* id, int* jd, int* kd, double* x, double* y, double* z, int* ifl);
  void dread_xyz_3d_iblank_ (
      int* id, int* jd, int* kd, double* x, double* y, double* z, int* iblank, int* ifl);
  void dread_xyz_3d_iblank_formatted_ (
      int* id, int* jd, int* kd, double* x, double* y, double* z, int* iblank, int* ifl);
  void dread_xyz_2d_ (
      int* id, int* jd, double* x, double* y, int* ifl);
  void dread_xyz_2d_formatted_  (
      int* id, int* jd, double* x, double* y, int* ifl);
  void dread_q_3d_ (
      int* id, int* jd, int* kd, double* fsmach, double* alpha, double* re, double* time, double* q, int* ifl);
  void dread_q_3d_formatted_ (
      int* id, int* jd, int* kd, double* fsmach, double* alpha, double* re, double* time, double* q, int* ifl);
  void dread_q_2d_ (
      int* id, int* jd, double* fsmach, double* alpha, double* re, double* time, double* q, int* ifl);
  void dread_q_2d_formatted_ (
      int* id, int* jd, double* fsmach, double* alpha, double* re, double* time, double* q, int* ifl);
  void dread_func_3d_ (
      int* id, int* jd, int* kd, int* nvar, double* d, int* ifl);
  void dread_func_3d_formatted_ (
      int* id, int* jd, int* kd, int* nvar, double* d, int* ifl);
  void dread_func_2d_ (
      int* id, int* jd, int* nvar, double* d, int* ifl);
  void dread_func_2d_formatted_ (
      int* id, int* jd, int* nvar, double* d, int* ifl);
}

using namespace std;

class FileIO_PLOT3D_READ : public FileIO_PLOT3D {
  
public:
  /** コンストラクタ */
  FileIO_PLOT3D_READ() {}

  /**　デストラクタ */
  ~FileIO_PLOT3D_READ() {}
  
private:

protected:

  /**
   * @brief グリッド座標ファイルの読み込み（UNFORMATTED）
   */
  void ReadXYZ_UNFORMATTED();

  /**
   * @brief グリッド座標ファイルの読み込み（FORMATTED）
   */
  void ReadXYZ_FORMATTED();

  /**
   * @brief グリッド座標ファイルの読み込み（C_BINARY）
   */
  void ReadXYZ_C_BINARY();

  /**
   * @brief 計算結果Qファイルの読み込み（UNFORMATTED）
   */
  void ReadQ_UNFORMATTED();

  /**
   * @brief 計算結果Qファイルの読み込み（FORMATTED）
   */
  void ReadQ_FORMATTED();

  /**
   * @brief 計算結果Qファイルの読み込み（C_BINARY）
   */
  void ReadQ_C_BINARY();

  /**
   * @brief 計算結果Functionファイルの読み込み（UNFORMATTED）
   */
  void ReadFunc_UNFORMATTED();

  /**
   * @brief 計算結果Functionファイルの読み込み（FORMATTED）
   */
  void ReadFunc_FORMATTED();

  /**
   * @brief 計算結果Functionファイルの読み込み（C_BINARY）
   */
  void ReadFunc_C_BINARY();

public:

  /**
   * @brief PLOT3Dの入力ファイルオープン
   */
  bool OpenFile();

  /**
   * @brief グリッド数の読み込み
   */
  void ReadNgrid(int* ngrid);

  /**
   * @brief ブロックデータの読み込み
   */
  void ReadBlockData(int* id, int* jd, int* kd);

  /**
   * @brief ブロックデータの読み込み（2D）
   */
  void ReadBlockData2D(int* id, int* jd);

  /**
   * @brief 計算結果Functionファイル、ブロックデータの読み込み
   */
  void ReadFuncBlockData(int* id, int* jd, int* kd, int* nvar, int ngrid);

  /**
   * @brief グリッド座標ファイルの読み込み
   */
  bool ReadXYZData();

  /**
   * @brief 計算結果Qファイルの読み込み
   */
  bool ReadQData(
    float* m_fsmach,
    float* m_alpha,
    float* m_re,
    float* m_time,
    float* m_q);

  /**
   * @brief 計算結果Qファイルの読み込み
   */
  bool ReadQData(
    double* m_fsmach,
    double* m_alpha,
    double* m_re,
    double* m_time,
    double* m_q);

  /**
   * @brief 計算結果Functionファイルの読み込み
   */
  bool ReadFuncData();

  /**
   * @brief FunctionNameファイルの読み込み
   */
  void ReadFunctionName(string* fn);

  /**
   * @brief 文字列の読み込み
   */
  void ReadSTRING(string* m_buff);

  /**
   * @brief FVBNDファイルの読み込み
   */
  void ReadFVBND(
    int* m_type, int* m_gridnum,
    int* m_Imin,int* m_Imax,
    int* m_Jmin,int* m_Jmax,
    int* m_Kmin,int* m_Kmax,
    string* m_ResultFlag,
    int* m_dir);

};

#endif // _PLT3D_PLOT3D_READ_H_
