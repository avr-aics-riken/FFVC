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

#include <mpi.h>

#include "Control.h"
#include "dfi.h"
#include "PLOT3D_write.h"
#include "PLOT3D_read.h"

class Plot3D {
private:
  int myRank;
  int size[3];
  int guide;
  
public:
  
  /** コンストラクタ */
  Plot3D() {
    myRank = 0;
    size[0] = size[1] = size[2] = 0;
    guide = 0;
  }
  
  /**　デストラクタ */
  ~Plot3D() {}
  
  
protected:
  
  FBUtility U;      ///< ユーティリティクラス
  
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
   * @param [in,out] C           Controlクラス
   * @param [in]     CurrentStep 現在のステップ数
   * @param [in]     CurrentTime 現在の時刻
   * @param [in]     FP3DW       FileIO_PLOT3D_WRITE
   * @param [in]     DFI         DFI
   * @param [in]     d_ws        出力データ配列
   * @param [in]     d_p         圧力データ配列 
   * @param [in,out] flop        浮動小数点演算数
   */
  void OutputPlot3D_function(Control* C,
                             unsigned CurrentStep,
                             double CurrentTime,
                             FileIO_PLOT3D_WRITE* FP3DW,
                             DFI* DFI,
                             REAL_TYPE* d_ws,
                             REAL_TYPE* d_p,
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
   * @param [in,out] d_bcd    BCindex
   */
  void setIblank(int* iblank, int id, int jd, int kd, int* d_bcd);
  
  
  /**
   * @brief Iblankのセット（ガイドセルに値があることを想定しているバージョン）
   * @param [out]    iblank   iblank = 1 : 計算グリッド, = 0 : 非計算グリッド, = 2 : 壁面グリッド
   * @param [in]     id       iblanx x方向サイズ
   * @param [in]     jd       iblanx y方向サイズ
   * @param [in]     kd       iblanx z方向サイズ
   * @param [in,out] d_bcd    BCindex
   */
  void setIblankGuide(int* iblank, int id, int jd, int kd, int* d_bcd);
  
  
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
   * @param [in,out] C           Controlクラス
   * @param [in]     FP3DW       FileIO_PLOT3D_WRITE
   * @param [in]     DFI         DFI
   */
  void OutputPlot3D_function_name(Control* C, FileIO_PLOT3D_WRITE* FP3DW, DFI* DFI);
  
  
  /**
   * @brief 項目別計算結果ファイルの項目（*.nam）出力
   * @param [in,out] C           Controlクラス
   * @param [in]     FP3DW       FileIO_PLOT3D_WRITE
   * @param [in]     DFI         DFI
   */
  void OutputPlot3D_function_name_divide(Control* C, FileIO_PLOT3D_WRITE* FP3DW, DFI* DFI);
  
  
  /**
   * @brief 項目別計算結果ファイル（*.func）出力
   * @param [in,out] flop    浮動小数点演算数
   */
  void OutputPlot3D_function_divide(double& flop);
  
  
  /**
   * @brief 境界面定義ファイル（*.fvbnd）出力（未整備）
   * @note BCMになるとりメッシュされたときに対応できないため出力することはない？
   */
  void OutputPlot3D_fvbnd();
  
  
  /**
   * @brief PLOT3Dファイルのポスト出力
   * @param [in]     CurrentStep 現在のステップ数
   * @param [in]     CurrentTime 現在の時刻
   * @param [in,out] C           Controlクラス
   * @param [in]     FP3DW       FileIO_PLOT3D_WRITE
   * @param [in]     DFI         DFI
   * @param [in]     d_ws        出力データ配列 
   * @param [in]     d_p         圧力データ配列 
   * @param [in,out] flop        浮動小数点演算数
   */
  void OutputPlot3D_post(unsigned CurrentStep,
                         double CurrentTime,
                         Control* C,
                         FileIO_PLOT3D_WRITE* FP3DW,
                         DFI* DFI,
                         REAL_TYPE* d_ws,
                         REAL_TYPE* d_p,
                         double& flop);
  
  
  /**
   * @brief 形状データファイル（*.xyz）の出力
   * @param [in]     CurrentStep 現在のステップ数
   * @param [in,out] C           Controlクラス
   * @param [in]     FP3DW       FileIO_PLOT3D_WRITE
   * @param [in]     DFI         DFI
   */
  void OutputPlot3D_xyz(unsigned CurrentStep, Control* C, FileIO_PLOT3D_WRITE* FP3DW, DFI* DFI);
  
  
  /**
   * @brief グリッド数、出力項目数をセット（BCMへ移行する場合要編集）
   * @param [in,out] C Controlクラス
   */
  void setValuePlot3D(Control* C);

  
  /**
   * @brief メンバ変数の初期化
   * @param [in]     m_size        配列サイズ
   * @param [in]     m_guide       ガイドセルサイズ
   */
  void Initialize(const int* m_size, const int m_guide);
};