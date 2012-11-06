#ifndef _PLT3D_FFV_PLOT3D_H_
#define _PLT3D_FFV_PLOT3D_H_
// #################################################################
//
// FFV : Frontflow / violet Cartesian
//
// Copyright (c) 2012 All right reserved.
//
// Institute of Industrial Science, University of Tokyo, Japan.
//
// #################################################################

/**
 * @file   ffv_PLOT3D.h
 * @brief  Plot3D Class
 * @author kero
 */

#include <string>
#include <mpi.h>

#include "FB_Ffunc.h"
#include "Control.h"
#include "dfi.h"
#include "PLOT3D_write.h"
#include "PLOT3D_read.h"
#include "ffv_Ffunc.h"

class Plot3D {
private:
  int myRank;        ///< ランク番号
  int size[3];       ///< 配列サイズ
  int guide;         ///< ガイドセルサイズ
  REAL_TYPE deltaX;  ///< 等間隔格子の無次元格子幅
  int dfi_plot3d;    ///< PLOT3D管理用
  
  Control* C;                 ///< Controlクラス
  FileIO_PLOT3D_WRITE* FP3DW; ///< PLOT3D format wirte class
  DFI* dfi;                   ///< Distributed File Information
  
  REAL_TYPE* d_ws;   ///< 出力データコンテナ配列
  REAL_TYPE* d_p;    ///< 圧力データ配列
  REAL_TYPE* d_wo;   ///< ベクトル出力コンテナ配列
  REAL_TYPE* d_v;    ///< 速度ベクトルデータ配列
  REAL_TYPE* d_t;    ///< 温度データ配列
  REAL_TYPE* d_p0;   ///< スカラーデータ配列
  REAL_TYPE* d_wv;   ///< ベクトルデータ配列
  int* d_bcv;        ///< BCindexV
  int* d_bcd;        ///< BCindexID

public:
  
  /** コンストラクタ */
  Plot3D() {
    myRank = 0;
    size[0] = size[1] = size[2] = 0;
    guide = 0;
    deltaX = 0.0;
    dfi_plot3d = 0;
    
    d_ws = NULL;
    d_p  = NULL;
    d_wo = NULL;
    d_v  = NULL;
    d_t  = NULL;
    d_p0 = NULL;
    d_wv = NULL;
    d_bcv= NULL;
    d_bcd= NULL;
  }
  
  /**　デストラクタ */
  ~Plot3D() {}
  
  
protected:
  
  FBUtility U;      ///< ユーティリティクラス
  
  /**
   * @brief 指定の出力ディレクトリとファイル名を結合
   * @retval フルパスファイル名
   * @param [in] path      ディレクトリパス名
   * @param [in] fname     ファイル名
   * @param [in] io_mode   ファイル出力モード
   @ @param [in] para_mode 並列モード
   */
  string directory_prefix(string path, const string fname, const int io_mode, const int para_mode);
  
  
// float
  /**
   * @brief 内部の格子点のデータを8で割る
   * @param [out]    d        格子点data
   * @param [in]     id       セル中心d x方向サイズ
   * @param [in]     jd       セル中心d y方向サイズ
   * @param [in]     kd       セル中心d z方向サイズ
   */
  void VolumeDataDivideBy8(float* d, int id, int jd, int kd);
  
  
  /**
   * @brief 面上の格子点のデータを4で割る
   * @param [out]    d        格子点data
   * @param [in]     id       セル中心d x方向サイズ
   * @param [in]     jd       セル中心d y方向サイズ
   * @param [in]     kd       セル中心d z方向サイズ
   */
  void FaceDataDivideBy4(float* d, int id, int jd, int kd);
  
  
  /**
   * @brief 辺上の格子点のデータを2で割る
   * @param [out]    d        格子点data
   * @param [in]     id       セル中心d x方向サイズ
   * @param [in]     jd       セル中心d y方向サイズ
   * @param [in]     kd       セル中心d z方向サイズ
   */
  void LineDataDivideBy2(float* d, int id, int jd, int kd);
  
  
// double
  /**
   * @brief 内部の格子点のデータを8で割る
   * @param [out]    d        格子点data
   * @param [in]     id       セル中心d x方向サイズ
   * @param [in]     jd       セル中心d y方向サイズ
   * @param [in]     kd       セル中心d z方向サイズ
   */
  void VolumeDataDivideBy8(double* d, int id, int jd, int kd);
  
  
  /**
   * @brief 面上の格子点のデータを4で割る
   * @param [out]    d        格子点data
   * @param [in]     id       セル中心d x方向サイズ
   * @param [in]     jd       セル中心d y方向サイズ
   * @param [in]     kd       セル中心d z方向サイズ
   */
  void FaceDataDivideBy4(double* d, int id, int jd, int kd);
  
  
  /**
   * @brief 辺上の格子点のデータを2で割る
   * @param [out]    d        格子点data
   * @param [in]     id       セル中心d x方向サイズ
   * @param [in]     jd       セル中心d y方向サイズ
   * @param [in]     kd       セル中心d z方向サイズ
   */
  void LineDataDivideBy2(double* d, int id, int jd, int kd);
  
  
  /**
   * @brief 計算結果ファイル（*.func）出力
   * @param [in]     CurrentStep 現在のステップ数
   * @param [in]     CurrentTime 現在の時刻
   * @param [in]     v00         格子速度
   * @param [in,out] flop        浮動小数点演算数
   */
  void OutputPlot3D_function(const unsigned CurrentStep,
                             const double CurrentTime,
                             REAL_TYPE* v00,
                             double& flop);
  
  
  /**
   * @brief 項目別計算結果ファイル（*.func）出力
   * @param [in]     CurrentStep 現在のステップ数
   * @param [in]     CurrentTime 現在の時刻
   * @param [in]     v00         格子速度
   * @param [in,out] flop    浮動小数点演算数
   */
  void OutputPlot3D_function_divide(const unsigned CurrentStep,
                                    const double CurrentTime,
                                    REAL_TYPE* v00,
                                    double& flop);
  
  
  /**
   * @brief 圧縮性流体のための計算結果ファイル（*.q）出力（未整備）
   * @param [in,out] flop    浮動小数点演算数
   */
  void OutputPlot3D_q(double& flop);
  
  
  /**
   * @brief Iblankのセット（ガイドセルの値は計算対象にいれていない）
   * @param [out]    iblank   iblank = 1 : 計算グリッド, = 0 : 非計算グリッド, = 2 : 壁面グリッド
   * @param [in]     id       iblanx x方向サイズ
   * @param [in]     jd       iblanx y方向サイズ
   * @param [in]     kd       iblanx z方向サイズ
   */
  void setIblank(int* iblank, int id, int jd, int kd);
  
  
  /**
   * @brief Iblankのセット（ガイドセルに値があることを想定しているバージョン）
   * @param [out]    iblank   iblank = 1 : 計算グリッド, = 0 : 非計算グリッド, = 2 : 壁面グリッド
   * @param [in]     id       iblanx x方向サイズ
   * @param [in]     jd       iblanx y方向サイズ
   * @param [in]     kd       iblanx z方向サイズ
   */
  void setIblankGuide(int* iblank, int id, int jd, int kd);
  
  
// flaot >> float
  /**
   * @brief Scalarの格子点での値をセット（float -> float）
   * @param [out]    d        格子点data
   * @param [in]     data     セル中心data
   * @param [in]     id       セル中心data x方向サイズ
   * @param [in]     jd       セル中心data y方向サイズ
   * @param [in]     kd       セル中心data z方向サイズ
   */
  void setScalarGridData(float* d, float* data, int id, int jd, int kd);

  
  /**
   * @brief Scalarの格子点での値をセット（float -> float）
   * @param [out]    d        格子点data
   * @param [in]     data     セル中心data
   * @param [in]     id       セル中心data x方向サイズ
   * @param [in]     jd       セル中心data y方向サイズ
   * @param [in]     kd       セル中心data z方向サイズ
   */
  void setScalarGridDataGuide(float* d, float* data, int id, int jd, int kd, int gc_out);
  
  
  /**
   * @brief Vectorの格子点での値をセット（float -> float）
   * @param [out]    d        格子点data
   * @param [in]     data     セル中心data
   * @param [in]     id       セル中心data x方向サイズ
   * @param [in]     jd       セル中心data y方向サイズ
   * @param [in]     kd       セル中心data z方向サイズ
   */
  void setVectorGridData(float* d, float* data, int id, int jd, int kd);
  
  
  /**
   * @brief Vectorの格子点での値をセット（float -> float）
   * @param [out]    d        格子点data
   * @param [in]     data     セル中心data
   * @param [in]     id       セル中心data x方向サイズ
   * @param [in]     jd       セル中心data y方向サイズ
   * @param [in]     kd       セル中心data z方向サイズ
   */
  void setVectorGridDataGuide(float* d, float* data, int id, int jd, int kd, int gc_out);
  
  
  /**
   * @brief 成分別Vectorの格子点での値をセット（float -> float）
   * @param [out]    d        格子点data
   * @param [in]     data     セル中心data
   * @param [in]     id       セル中心data x方向サイズ
   * @param [in]     jd       セル中心data y方向サイズ
   * @param [in]     kd       セル中心data z方向サイズ
   * @param [in]     ivar     ベクトル成分 =0:x =1:y =2:z
   */
  void setVectorComponentGridData(float* d, float* data, int id, int jd, int kd, int ivar);
  
  
  /**
   * @brief 成分別Vectorの格子点での値をセット（float -> float）
   * @param [out]    d        格子点data
   * @param [in]     data     セル中心data
   * @param [in]     id       セル中心data x方向サイズ
   * @param [in]     jd       セル中心data y方向サイズ
   * @param [in]     kd       セル中心data z方向サイズ
   * @param [in]     ivar     ベクトル成分 =0:x =1:y =2:z
   */
  void setVectorComponentGridDataGuide(float* d, float* data, int id, int jd, int kd, int gc_out, int ivar);
  
  
// double >> double
  /**
   * @brief Scalarの格子点での値をセット（double -> double）
   * @param [out]    d        格子点data
   * @param [in]     data     セル中心data
   * @param [in]     id       セル中心data x方向サイズ
   * @param [in]     jd       セル中心data y方向サイズ
   * @param [in]     kd       セル中心data z方向サイズ
   */
  void setScalarGridData(double* d, double* data, int id, int jd, int kd);
  
  
  /**
   * @brief Scalarの格子点での値をセット（double -> double）
   * @param [out]    d        格子点data
   * @param [in]     data     セル中心data
   * @param [in]     id       セル中心data x方向サイズ
   * @param [in]     jd       セル中心data y方向サイズ
   * @param [in]     kd       セル中心data z方向サイズ
   */
  void setScalarGridDataGuide(double* d, double* data, int id, int jd, int kd, int gc_out);
  
  
  /**
   * @brief Vectorの格子点での値をセット（double -> double）
   * @param [out]    d        格子点data
   * @param [in]     data     セル中心data
   * @param [in]     id       セル中心data x方向サイズ
   * @param [in]     jd       セル中心data y方向サイズ
   * @param [in]     kd       セル中心data z方向サイズ
   */
  void setVectorGridData(double* d, double* data, int id, int jd, int kd);
  
  
  /**
   * @brief Vectorの格子点での値をセット（double -> double）
   * @param [out]    d        格子点data
   * @param [in]     data     セル中心data
   * @param [in]     id       セル中心data x方向サイズ
   * @param [in]     jd       セル中心data y方向サイズ
   * @param [in]     kd       セル中心data z方向サイズ
   */
  void setVectorGridDataGuide(double* d, double* data, int id, int jd, int kd, int gc_out);
  
  
  /**
   * @brief 成分別Vectorの格子点での値をセット（double -> double）
   * @param [out]    d        格子点data
   * @param [in]     data     セル中心data
   * @param [in]     id       セル中心data x方向サイズ
   * @param [in]     jd       セル中心data y方向サイズ
   * @param [in]     kd       セル中心data z方向サイズ
   * @param [in]     ivar     ベクトル成分 =0:x =1:y =2:z
   */
  void setVectorComponentGridData(double* d, double* data, int id, int jd, int kd, int ivar);
  
  
  /**
   * @brief 成分別Vectorの格子点での値をセット（double -> double）
   * @param [out]    d        格子点data
   * @param [in]     data     セル中心data
   * @param [in]     id       セル中心data x方向サイズ
   * @param [in]     jd       セル中心data y方向サイズ
   * @param [in]     kd       セル中心data z方向サイズ
   * @param [in]     ivar     ベクトル成分 =0:x =1:y =2:z
   */
  void setVectorComponentGridDataGuide(double* d, double* data, int id, int jd, int kd, int gc_out, int ivar);
  
  
  
// double >> float
  /**
   * @brief Scalarの格子点での値をセット（double -> float）
   * @param [out]    d        格子点data
   * @param [in]     data     セル中心data
   * @param [in]     id       セル中心data x方向サイズ
   * @param [in]     jd       セル中心data y方向サイズ
   * @param [in]     kd       セル中心data z方向サイズ
   */
  void setScalarGridData(float* d, double* data, int id, int jd, int kd);
  
  
  /**
   * @brief Scalarの格子点での値をセット（double -> float）
   * @param [out]    d        格子点data
   * @param [in]     data     セル中心data
   * @param [in]     id       セル中心data x方向サイズ
   * @param [in]     jd       セル中心data y方向サイズ
   * @param [in]     kd       セル中心data z方向サイズ
   */
  void setScalarGridDataGuide(float* d, double* data, int id, int jd, int kd, int gc_out);
  
  
  /**
   * @brief Vectorの格子点での値をセット（double -> float）
   * @param [out]    d        格子点data
   * @param [in]     data     セル中心data
   * @param [in]     id       セル中心data x方向サイズ
   * @param [in]     jd       セル中心data y方向サイズ
   * @param [in]     kd       セル中心data z方向サイズ
   */
  void setVectorGridData(float* d, double* data, int id, int jd, int kd);
  
  
  /**
   * @brief Vectorの格子点での値をセット（double -> float）
   * @param [out]    d        格子点data
   * @param [in]     data     セル中心data
   * @param [in]     id       セル中心data x方向サイズ
   * @param [in]     jd       セル中心data y方向サイズ
   * @param [in]     kd       セル中心data z方向サイズ
   */
  void setVectorGridDataGuide(float* d, double* data, int id, int jd, int kd, int gc_out);
  
  
  /**
   * @brief 成分別Vectorの格子点での値をセット（double -> float）
   * @param [out]    d        格子点data
   * @param [in]     data     セル中心data
   * @param [in]     id       セル中心data x方向サイズ
   * @param [in]     jd       セル中心data y方向サイズ
   * @param [in]     kd       セル中心data z方向サイズ
   * @param [in]     ivar     ベクトル成分 =0:x =1:y =2:z
   */
  void setVectorComponentGridData(float* d, double* data, int id, int jd, int kd, int ivar);
  
  
  /**
   * @brief 成分別Vectorの格子点での値をセット（double -> float）
   * @param [out]    d        格子点data
   * @param [in]     data     セル中心data
   * @param [in]     id       セル中心data x方向サイズ
   * @param [in]     jd       セル中心data y方向サイズ
   * @param [in]     kd       セル中心data z方向サイズ
   * @param [in]     ivar     ベクトル成分 =0:x =1:y =2:z
   */
  void setVectorComponentGridDataGuide(float* d, double* data, int id, int jd, int kd, int gc_out, int ivar);
  
  
  
// float >> double
  /**
   * @brief Scalarの格子点での値をセット（float -> double）
   * @param [out]    d        格子点data
   * @param [in]     data     セル中心data
   * @param [in]     id       セル中心data x方向サイズ
   * @param [in]     jd       セル中心data y方向サイズ
   * @param [in]     kd       セル中心data z方向サイズ
   */
  void setScalarGridData(double* d, float* data, int id, int jd, int kd);
  
  
  /**
   * @brief Scalarの格子点での値をセット（float -> double）
   * @param [out]    d        格子点data
   * @param [in]     data     セル中心data
   * @param [in]     id       セル中心data x方向サイズ
   * @param [in]     jd       セル中心data y方向サイズ
   * @param [in]     kd       セル中心data z方向サイズ
   */
  void setScalarGridDataGuide(double* d, float* data, int id, int jd, int kd, int gc_out);
  
  

  /**
   * @brief Vectorの格子点での値をセット（float -> double）
   * @param [out]    d        格子点data
   * @param [in]     data     セル中心data
   * @param [in]     id       セル中心data x方向サイズ
   * @param [in]     jd       セル中心data y方向サイズ
   * @param [in]     kd       セル中心data z方向サイズ
   */
  void setVectorGridData(double* d, float* data, int id, int jd, int kd);
  
  
  /**
   * @brief Vectorの格子点での値をセット（float -> double）
   * @param [out]    d        格子点data
   * @param [in]     data     セル中心data
   * @param [in]     id       セル中心data x方向サイズ
   * @param [in]     jd       セル中心data y方向サイズ
   * @param [in]     kd       セル中心data z方向サイズ
   */
  void setVectorGridDataGuide(double* d, float* data, int id, int jd, int kd, int gc_out);
  
  
  /**
   * @brief 成分別Vectorの格子点での値をセット（float -> double）
   * @param [out]    d        格子点data
   * @param [in]     data     セル中心data
   * @param [in]     id       セル中心data x方向サイズ
   * @param [in]     jd       セル中心data y方向サイズ
   * @param [in]     kd       セル中心data z方向サイズ
   * @param [in]     ivar     ベクトル成分 =0:x =1:y =2:z
   */
  void setVectorComponentGridData(double* d, float* data, int id, int jd, int kd, int ivar);
  
  
  /**
   * @brief 成分別Vectorの格子点での値をセット（float -> double）
   * @param [out]    d        格子点data
   * @param [in]     data     セル中心data
   * @param [in]     id       セル中心data x方向サイズ
   * @param [in]     jd       セル中心data y方向サイズ
   * @param [in]     kd       セル中心data z方向サイズ
   * @param [in]     ivar     ベクトル成分 =0:x =1:y =2:z
   */
  void setVectorComponentGridDataGuide(double* d, float* data, int id, int jd, int kd, int gc_out, int ivar);
  
  
public:

  /**
   * @brief 計算結果ファイルの項目（*.nam）出力
   */
  void OutputPlot3D_function_name();
  
  
  /**
   * @brief 項目別計算結果ファイルの項目（*.nam）出力
   */
  void OutputPlot3D_function_name_divide();
  
  
  /**
   * @brief 境界面定義ファイル（*.fvbnd）出力（未整備）
   * @note BCMになるとりメッシュされたときに対応できないため出力することはない？
   */
  void OutputPlot3D_fvbnd();
  
  
  /**
   * @brief PLOT3Dファイルのポスト出力
   * @param [in]     CurrentStep 現在のステップ数
   * @param [in]     CurrentTime 現在の時刻
   * @param [in]     v00         格子速度
   * @param [in]     origin      基点座標
   * @param [in]     pitch       格子幅
   * @param [in,out] flop        浮動小数点演算数
   */
  void OutputPlot3D_post(const unsigned CurrentStep,
                         const double CurrentTime,
                         REAL_TYPE* v00,
                         const REAL_TYPE* origin,
                         const REAL_TYPE* pitch,
                         double& flop);
  
  
  /**
   * @brief 形状データファイル（*.xyz）の出力
   * @param [in]     CurrentStep 現在のステップ数
   * @param [in]     origin      基点座標
   * @param [in]     pitch       格子幅
   */
  void OutputPlot3D_xyz(const unsigned CurrentStep, const REAL_TYPE* origin, const REAL_TYPE* pitch);
  
  
  /**
   * @brief グリッド数、出力項目数をセット（BCMへ移行する場合要編集）
   */
  void setValuePlot3D();

  
  /**
   * @brief メンバ変数の初期化
   * @param [in]     m_size      配列サイズ
   * @param [in]     m_guide     ガイドセルサイズ
   * @param [in]     m_deltaX    格子幅
   * @param [in]     dfi_plot3d  PLOT3D管理用
   * @param [in]     m_C         Controlクラス
   * @param [in]     m_FP3DW     FileIO_PLOT3D_WRITE
   * @param [in]     m_DFI       DFI
   * @param [in]     m_d_ws      出力データコンテナ配列
   * @param [in]     m_d_p       圧力データ
   * @param [in]     m_d_wo      ベクトル出力データコンテナ配列
   * @param [in]     m_d_v       速度データ
   * @param [in]     m_d_t       温度データ
   * @param [in]     m_d_p0      スカラーデータ
   * @param [in]     m_d_wv      ベクトルデータ
   * @param [in]     m_d_bcv     BCindexV
   * @param [in]     m_d_bcd     BCindexID
   */
  void Initialize(const int* m_size,
                  const int m_guide,
                  const REAL_TYPE m_deltaX,
                  const int dfi_plot3d,
                  Control* m_C,
                  FileIO_PLOT3D_WRITE* m_FP3DW,
                  DFI* m_dfi,
                  REAL_TYPE* m_d_ws,
                  REAL_TYPE* m_d_p,
                  REAL_TYPE* m_d_wo,
                  REAL_TYPE* m_d_v,
                  REAL_TYPE* m_d_t,
                  REAL_TYPE* m_d_p0,
                  REAL_TYPE* m_d_wv,
                  int*       m_d_bcv,
                  int*       m_d_bcd);
};

#endif // _PLT3D_FFV_PLOT3D_H_
