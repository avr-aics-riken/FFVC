#ifndef _FileIO_SPH_H_
#define _FileIO_SPH_H_

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
 * @file   FileIO_SPH.h
 * @brief  FileIO_SPH Class Header
 * @author kero
 */

#include <stdio.h>
#include <string.h>

#include "endianUtil.h"
#include "FB_Define.h"

/**
 * SPHファイル入出力クラス
 */
class FileIO_SPH { // : public Base {

public:

  // output method (_OUT_METHOD1, _OUT_METHOD2)
  typedef enum {_OUT_UNKNOWN=0, _OUT_METHOD1, _OUT_METHOD2} OutType;   ///< 出力タイプ
  // data type (float or double)
  typedef enum {_REAL_UNKNOWN=0, _FLOAT, _DOUBLE} RealType;       ///< データ型フラグ
  // data dims (scalar or vector)
  typedef enum {_DATA_UNKNOWN=0, _SCALAR, _VECTOR} DataType;       ///< データ種別フラグ

private:
  OutType m_outType;
  RealType m_rType;
  DataType m_dType;

  char* m_fBase;
  unsigned int m_fSeqNumLen;
  char* m_fSuffix;

  unsigned int m_iSize[3];
  float m_fOrg[3], m_fPit[3];
  int m_fStep;
  long long m_dStep;
  float m_fTime;
  double m_dTime;

  unsigned long long  m_lSize[3];
  double m_dOrg[3], m_dPit[3];
  int* m_fRefStep;
  long long* m_dRefStep;
  float* m_fRefTime;
  double* m_dRefTime;

  const float* const * m_fData;
  const double* const * m_dData;

  /**
   * 出力ファイル名を作成する。
   * @param num   ステップ数
   * @return    出力ファイル名
   */
  char* MakeOutFileName(unsigned int stepNum);

  /**
   * SPHファイルに書込を行う。
   * @param fp    ファイルポインタ
   * @return      true=success, false=fault
   */
  bool WriteData(FILE* fp);

  /**
   * 変数初期化を行う。
   */
  void Initialize(){
    m_fBase = NULL;
    m_fSeqNumLen = 0;
    m_fSuffix = NULL;
    m_outType = _OUT_UNKNOWN;
    m_rType = _REAL_UNKNOWN;
    m_dType = _DATA_UNKNOWN;
    m_iSize[0] = m_iSize[1] = m_iSize[2] = 0;
    m_fOrg[0] = m_fOrg[1] = m_fOrg[2] = 0.0f;
    m_fPit[0] = m_fPit[1] = m_fPit[2] = 0.0f;
    m_lSize[0] = m_lSize[1] = m_lSize[2] = 0;
    m_dOrg[0] = m_dOrg[1] = m_dOrg[2] = 0.0;
    m_dPit[0] = m_dPit[1] = m_dPit[2] = 0.0;
    m_fStep = 0;
    m_dStep = 0;
    m_fTime = 0.0f;
    m_dTime = 0.0;
    m_fRefStep = NULL;
    m_dRefStep = NULL;
    m_fRefTime = NULL;
    m_dRefTime = NULL;
    m_fData = NULL;
    m_dData = NULL;
  };

public:
  /**
   * コンストラクタ
   */
  FileIO_SPH(){
    Initialize();
  };

  /**
   * デストラクタ
   */
  virtual ~FileIO_SPH(){
    if( m_fBase ) delete [] m_fBase;
    if( m_fSuffix ) delete [] m_fSuffix;
    Initialize();
  };

  /**
   * ファイル出力情報を設定する。(float型)
   * @param dtype     データ型
   * @param voxsize   ボクセルサイズ
   * @param voxorg    原点座標
   * @param voxpit    ピッチ
   * @param step      ステップ数
   * @param time      時間
   * @param data      SPHデータ
   * @return      true=success, false=fault
   */
  bool Init(int dtype,
            unsigned int voxsize[3],
            float voxorg[3],
            float voxpit[3],
            int* step,
            float* time,
            float** data);

  /**
   * ファイル出力情報を設定する。(double型)
   * @param dtype     データ型
   * @param voxsize   ボクセルサイズ
   * @param voxorg    原点座標
   * @param voxpit    ピッチ
   * @param step      ステップ数
   * @param time      時間
   * @param data      SPHデータ
   * @return      true=success, false=fault
   */
  bool Init(int dtype,
            unsigned long long voxsize[3],
            double voxorg[3],
            double voxpit[3],
            long long* step,
            double* time,
            double** data);

  /**
   * ファイル出力情報を設定する。(ヘッダーレコード:float型)
   * @param dtype     データ型
   * @param voxsize   ボクセルサイズ
   * @param voxorg    原点座標
   * @param voxpit    ピッチ
   * @return      true=success, false=fault
   */
  bool Init(int dtype,
            unsigned int voxsize[3], float voxorg[3], float voxpit[3]);

  /**
   * ファイル出力情報を設定する。(ヘッダーレコード:double型)
   * @param dtype     データ型
   * @param voxsize   ボクセルサイズ
   * @param voxorg    原点座標
   * @param voxpit    ピッチ
   * @return      true=success, false=fault
   */
  bool Init(int dtype,
            unsigned long long voxsize[3], double voxorg[3], double voxpit[3]);

  /**
   * SPHファイルに出力を行う。
   * @param step    ステップ数
   * @return      true=success, false=fault
   */
  bool OutputBinaryFile(unsigned int step);

  /**
   * SPHファイルに出力を行う。(float型)
   * @param step      ステップ数
   * @param time      時間
   * @param data      SPHデータ
   * @return      true=success, false=fault
   */
  bool OutputBinaryFile(int step, float time, const float* data);

  /**
   * SPHファイルに出力を行う。(double型)
   * @param step      ステップ数
   * @param time      時間
   * @param data      SPHデータ
   * @return      true=success, false=fault
   */
  bool OutputBinaryFile(long long step, double time, const double* data);

  /**
   * ファイルヘッダー情報を取得する。(unsigned型)
   * @param[in] fname        ファイル名
   * @param[out] real_type   データ型タイプ
   * @param[out] data_type   データ種別
   * @param[out] voxsize     ボクセルサイズ
   * @return      true=success, false=fault
   */
  static bool GetInfo(const char* fname,
                      RealType* real_type,
                      DataType* data_type,
                      unsigned int voxsize[3]);

  /**
   * ファイルヘッダー情報を取得する。(long long型)
   * @param[in] fname        ファイル名
   * @param[out] real_type   データ型タイプ
   * @param[out] data_type   データ種別
   * @param[out] voxsize     ボクセルサイズ
   * @return      true=success, false=fault
   */
  static bool GetInfo(const char* fname,
                      RealType* real_type,
                      DataType* data_type,
                      unsigned long long voxsize[3]);

  /**
   * ファイルからSPHデータを読み込む。(float型)
   * @param[in] fname        ファイル名
   * @param[out] data_type   データ種別
   * @param[out] voxsize     ボクセルサイズ
   * @param[out] voxorg      原点座標
   * @param[out] voxpit      ピッチ
   * @param[out] step        ステップ数
   * @param[out] time        時間
   * @param[out] data        SPHデータ
   * @return      true=success, false=fault
   */
  static bool ReadFile(const char* fname,
                       DataType* data_type,
                       unsigned int voxsize[3],
                       float voxorg[3],
                       float voxpit[3],
                       int* step,
                       float* time,
                       float** data);

  /**
   * ファイルからSPHデータを読み込む。(double型)
   * @param[in] fname        ファイル名
   * @param[out] data_type   データ種別
   * @param[out] voxsize     ボクセルサイズ
   * @param[out] voxorg      原点座標
   * @param[out] voxpit      ピッチ
   * @param[out] step        ステップ数
   * @param[out] time        時間
   * @param[out] data        SPHデータ
   * @return      true=success, false=fault
   */
  static bool ReadFile(const char* fname,
                       DataType* data_type,
                       unsigned long long voxsize[3],
                       double voxorg[3],
                       double voxpit[3],
                       long long* step,
                       double* time,
                       double** data);


  /**
   * ファイルからSPHヘッダデータを読み込む。(float型)
   * @param[in] fname        ファイル名
   * @param[out] data_type   データ種別
   * @param[out] voxsize     ボクセルサイズ
   * @param[out] voxorg      原点座標
   * @param[out] voxpit      ピッチ
   * @param[out] step        ステップ数
   * @param[out] time        時間
   * @return      true=success, false=fault
   */
  static bool GetHeader(const char* fname,
                        //DataType* data_type,
                        int* data_type,
                        unsigned int voxsize[3],
                        float voxorg[3],
                        float voxpit[3],
                        int* step,
                        float* time);

  /**
   * ファイルからSPHヘッダデータを読み込む。(double型)
   * @param[in] fname        ファイル名
   * @param[out] data_type   データ種別
   * @param[out] voxsize     ボクセルサイズ
   * @param[out] voxorg      原点座標
   * @param[out] voxpit      ピッチ
   * @param[out] step        ステップ数
   * @param[out] time        時間
   * @return      true=success, false=fault
   */
  static bool GetHeader(const char* fname,
                        //DataType* data_type,
                        int* data_type,
                        unsigned long long voxsize[3],
                        double voxorg[3],
                        double voxpit[3],
                        long long* step,
                        double* time);

  /**
   * データ型を取得する。
   * @param[in] fname        ファイル名
   * @param[out] real_type   データ型タイプ
   * @param[out] data_type   データ種別
   * @param[out] voxsize     ボクセルサイズ
   * @return      true=success, false=fault
   */
  static bool GetDataType(const char* fname,
                          //RealType* real_type,
                          //DataType* data_type);
                          int* real_type,
                          int* data_type);

  /**
   * ファイルからSPHデータを読み込む。(float型) 
   * @param[in] fname        ファイル名
   * @param[in] dsize        データサイズ
   * @param[out] data_type   データ種別
   * @param[out] voxsize     ボクセルサイズ
   * @param[out] voxorg      原点座標
   * @param[out] voxpit      ピッチ
   * @param[out] step        ステップ数
   * @param[out] time        時間
   * @param[out] data        SPHデータ
   * @return      true=success, false=fault
   */
  static bool GetData(const char* fname,
                      unsigned int dsize,
                      int* data_type,
                      unsigned int voxsize[3],
                      float voxorg[3],
                      float voxpit[3],
                      int* step,
                      float* time,
                      float* data);

  /**
   * ファイルからSPHデータを読み込む。(double型)
   * @param[in] fname        ファイル名
   * @param[in] dsize        データサイズ
   * @param[out] data_type   データ種別
   * @param[out] voxsize     ボクセルサイズ
   * @param[out] voxorg      原点座標
   * @param[out] voxpit      ピッチ
   * @param[out] step        ステップ数
   * @param[out] time        時間
   * @param[out] data        SPHデータ
   * @return      true=success, false=fault
   */
  static bool GetData(const char* fname,
                      unsigned long long dsize,
                      int* data_type,
                      unsigned long long voxsize[3],
                      double voxorg[3],
                      double voxpit[3],
                      long long* step,
                      double* time,
                      double* data);

};

#endif // _FileIO_SPH_H_
