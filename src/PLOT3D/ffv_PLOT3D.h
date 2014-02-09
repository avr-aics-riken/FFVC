#ifndef _PLT3D_FFV_PLOT3D_H_
#define _PLT3D_FFV_PLOT3D_H_
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
 * @file   ffv_PLOT3D.h
 * @brief  Plot3D Class
 * @author aics
 */

#define PLT3D_VERS 103

#include "dfi.h"
#include <string>
#include <mpi.h>

#include "../FB/FB_Ffunc.h"
#include "../FB/Control.h"
#include "../F_CORE/ffv_Ffunc.h"
#include "PLOT3D_write.h"
#include "PLOT3D_read.h"


class Plot3D {
private:
  int myRank;        ///< ランク番号
  int size[3];       ///< 配列サイズ
  int guide;         ///< ガイドセルサイズ
  REAL_TYPE deltaX;  ///< 等間隔格子の無次元格子幅
  
  Control* C;                 ///< Controlクラス
  FileIO_PLOT3D_READ*  FP3DR; ///< PLOT3D format read class
  FileIO_PLOT3D_WRITE* FP3DW; ///< PLOT3D format wirte class
  DFI* dfi;                   ///< Distributed File Information
  
  REAL_TYPE* d_ws;   ///< 出力データコンテナ配列
  REAL_TYPE* d_p;    ///< 圧力データ配列
  REAL_TYPE* d_wo;   ///< ベクトル出力コンテナ配列
  REAL_TYPE* d_v;    ///< 速度ベクトルデータ配列
  REAL_TYPE* d_ie;   ///< 内部エネルギーデータ配列
  REAL_TYPE* d_p0;   ///< スカラーデータ配列
  REAL_TYPE* d_wv;   ///< ベクトルデータ配列
  int* d_cdf;        ///< BCindex Component Directional BC
  int* d_bcd;        ///< BCindex B

public:
  
  /** コンストラクタ */
  Plot3D() {
    myRank = 0;
    size[0] = size[1] = size[2] = 0;
    guide = 0;
    deltaX = 0.0;
    
    d_ws = NULL;
    d_p  = NULL;
    d_wo = NULL;
    d_v  = NULL;
    d_ie = NULL;
    d_p0 = NULL;
    d_wv = NULL;
    d_cdf= NULL;
    d_bcd= NULL;
    FP3DR= NULL;
    FP3DW= NULL;
  }
  
  /**　デストラクタ */
  ~Plot3D() {}
  
  
protected:
  
  FBUtility U;      ///< ユーティリティクラス

  
  /**
   * @brief 計算結果ファイル（*.func）出力
   * @param [in]     CurrentStep 現在のステップ数
   * @param [in]     CurrentTime 現在の時刻
   * @param [in]     v00         格子速度
   * @param [in,out] flop        浮動小数点演算数
   */
  void function(const unsigned CurrentStep,
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
  void function_divide(const unsigned CurrentStep,
                       const double CurrentTime,
                       REAL_TYPE* v00,
                       double& flop);
  
  
  /**
   * @brief 圧縮性流体のための計算結果ファイル（*.q）出力（未整備）
   * @param [in,out] flop    浮動小数点演算数
   */
  void Q(double& flop);
  
  
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
  

  /**
   * @brief 内部の格子点のデータを8で割る
   * @param [out]    d        格子点data
   * @param [in]     id       セル中心d x方向サイズ
   * @param [in]     jd       セル中心d y方向サイズ
   * @param [in]     kd       セル中心d z方向サイズ
   */
  template<class T>
  void VolumeDataDivideBy8(T* d, int id, int jd, int kd)
  {
    int i,j,k;
    size_t ip;
    
    for (k=1; k<kd-1; k++){
      for (j=1; j<jd-1; j++){
        for (i=1; i<id-1; i++){
          ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
          d[ip]=d[ip]*0.125;
        }
      }
    }
  };


  /**
   * @brief 面上の格子点のデータを4で割る
   * @param [out]    d        格子点data
   * @param [in]     id       セル中心d x方向サイズ
   * @param [in]     jd       セル中心d y方向サイズ
   * @param [in]     kd       セル中心d z方向サイズ
   */
  template<class T>
  void FaceDataDivideBy4(T* d, int id, int jd, int kd)
  {
    int i,j,k;
    size_t ip;
    
    i=0;
    for (k=1; k<kd-1; k++){
      for (j=1; j<jd-1; j++){
        ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
        d[ip]=d[ip]*0.25;
      }
    }
    
    i=id-1;
    for (k=1; k<kd-1; k++){
      for (j=1; j<jd-1; j++){
        ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
        d[ip]=d[ip]*0.25;
      }
    }
    
    j=0;
    for (k=1; k<kd-1; k++){
      for (i=1; i<id-1; i++){
        ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
        d[ip]=d[ip]*0.25;
      }
    }
    
    j=jd-1;
    for (k=1; k<kd-1; k++){
      for (i=1; i<id-1; i++){
        ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
        d[ip]=d[ip]*0.25;
      }
    }
    
    k=0;
    for (j=1; j<jd-1; j++){
      for (i=1; i<id-1; i++){
        ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
        d[ip]=d[ip]*0.25;
      }
    }
    
    k=kd-1;
    for (j=1; j<jd-1; j++){
      for (i=1; i<id-1; i++){
        ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
        d[ip]=d[ip]*0.25;
      }
    }
  };


  /**
   * @brief 辺上の格子点のデータを2で割る
   * @param [out]    d        格子点data
   * @param [in]     id       セル中心d x方向サイズ
   * @param [in]     jd       セル中心d y方向サイズ
   * @param [in]     kd       セル中心d z方向サイズ
   */
  template<class T>
  void LineDataDivideBy2(T* d, int id, int jd, int kd)
  {
    int i,j,k;
    size_t ip;
    
    i=0; j=0;
    for (k=1; k<kd-1; k++){
      ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
      d[ip]=d[ip]*0.5;
    }
    
    i=0; j=jd-1;
    for (k=1; k<kd-1; k++){
      ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
      d[ip]=d[ip]*0.5;
    }
    
    i=0; k=0;
    for (j=1; j<jd-1; j++){
      ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
      d[ip]=d[ip]*0.5;
    }
    
    i=0; k=kd-1;
    for (j=1; j<jd-1; j++){
      ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
      d[ip]=d[ip]*0.5;
    }
    
    j=0; k=0;
    for (i=1; i<id-1; i++){
      ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
      d[ip]=d[ip]*0.5;
    }
    
    j=0; k=kd-1;
    for (i=1; i<id-1; i++){
      ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
      d[ip]=d[ip]*0.5;
    }
    
    j=jd-1; k=0;
    for (i=1; i<id-1; i++){
      ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
      d[ip]=d[ip]*0.5;
    }
    
    j=jd-1; k=kd-1;
    for (i=1; i<id-1; i++){
      ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
      d[ip]=d[ip]*0.5;
    }
    
    i=id-1; j=0;
    for (k=1; k<kd-1; k++){
      ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
      d[ip]=d[ip]*0.5;
    }
    
    i=id-1; j=jd-1;
    for (k=1; k<kd-1; k++){
      ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
      d[ip]=d[ip]*0.5;
    }
    
    i=id-1; k=0;
    for (j=1; j<jd-1; j++){
      ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
      d[ip]=d[ip]*0.5;
    }
    
    i=id-1; k=kd-1;
    for (j=1; j<jd-1; j++){
      ip = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
      d[ip]=d[ip]*0.5;
    }
  
  };


  /**
   * @brief Scalarの格子点での値をセット
   * @param [out]    d        格子点data
   * @param [in]     data     セル中心data
   * @param [in]     id       セル中心data x方向サイズ
   * @param [in]     jd       セル中心data y方向サイズ
   * @param [in]     kd       セル中心data z方向サイズ
   */
  template<class T1, class T2>
  void setScalarGridData(T1* d, T2* data, int id, int jd, int kd, int gc_out)
  {
    int ix = size[0];
    int jx = size[1];
    int kx = size[2];
    int gd = guide;
 
    size_t mip;
    float ddd;
    int i,j,k;
 
    size_t dsize = (size_t)(id*jd*kd);
    
    for (size_t l=0; l<dsize; l++) d[l]=0.0;
    
    for (int km=1-gc_out; km<=kx+gc_out; km++) {
      for (int jm=1-gc_out; jm<=jx+gc_out; jm++) {
        for (int im=1-gc_out; im<=ix+gc_out; im++) {
          mip = _F_IDX_S3D(im, jm, km, ix, jx, kx, gd);
          ddd=(T1)data[mip];
          i=im-1+gc_out;
          j=jm-1+gc_out;
          k=km-1+gc_out;
          size_t ip1 = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
          size_t ip2 = _F_IDX_S3D(i+2, j+1, k+1, id, jd, kd, 0);
          size_t ip3 = _F_IDX_S3D(i+2, j+2, k+1, id, jd, kd, 0);
          size_t ip4 = _F_IDX_S3D(i+1, j+2, k+1, id, jd, kd, 0);
          size_t ip5 = _F_IDX_S3D(i+1, j+1, k+2, id, jd, kd, 0);
          size_t ip6 = _F_IDX_S3D(i+2, j+1, k+2, id, jd, kd, 0);
          size_t ip7 = _F_IDX_S3D(i+2, j+2, k+2, id, jd, kd, 0);
          size_t ip8 = _F_IDX_S3D(i+1, j+2, k+2, id, jd, kd, 0);
          d[ip1]=d[ip1]+ddd;
          d[ip2]=d[ip2]+ddd;
          d[ip3]=d[ip3]+ddd;
          d[ip4]=d[ip4]+ddd;
          d[ip5]=d[ip5]+ddd;
          d[ip6]=d[ip6]+ddd;
          d[ip7]=d[ip7]+ddd;
          d[ip8]=d[ip8]+ddd;
        }
      }
    }
    
    //内部の格子点のデータを8で割る
    VolumeDataDivideBy8(d, id, jd, kd);
    
    //面上の格子点のデータを4で割る
    FaceDataDivideBy4(d, id, jd, kd);
    
    //辺上の格子点のデータを2で割る
    LineDataDivideBy2(d, id, jd, kd);
  
  };

  /**
   * @brief Vectorの格子点での値をセット
   * @param [out]    d        格子点data
   * @param [in]     data     セル中心data
   * @param [in]     id       セル中心data x方向サイズ
   * @param [in]     jd       セル中心data y方向サイズ
   * @param [in]     kd       セル中心data z方向サイズ
   */
  template<class T1, class T2>
  void setVectorGridData(T1* d, T2* data, int id, int jd, int kd, int gc_out)
  {
    int ix = size[0];
    int jx = size[1];
    int kx = size[2];
    int gd = guide;
 
    size_t mip;
    T1 ddd;
    int i,j,k;
 
    size_t dsize = (size_t)(id*jd*kd);
    size_t dsize3 = (size_t)(id*jd*kd*3);
 
    for (size_t l=0; l<dsize3; l++) d[l]=0.0;
    
    for (size_t ivar=0;ivar<3;ivar++){
      
      for (int km=1-gc_out; km<=kx+gc_out; km++) {
        for (int jm=1-gc_out; jm<=jx+gc_out; jm++) {
          for (int im=1-gc_out; im<=ix+gc_out; im++) {
            //mip = _F_IDX_V3D(im, jm, km, ivar, ix, jx, kx, gd); //(i,j,k,3)
            mip = _F_IDX_V3DEX(ivar, im, jm, km, ix, jx, kx, gd); //(3,i,j,k)
            ddd=(T1)data[mip];
            i=im-1+gc_out;
            j=jm-1+gc_out;
            k=km-1+gc_out;
            size_t ip1 = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
            size_t ip2 = _F_IDX_S3D(i+2, j+1, k+1, id, jd, kd, 0);
            size_t ip3 = _F_IDX_S3D(i+2, j+2, k+1, id, jd, kd, 0);
            size_t ip4 = _F_IDX_S3D(i+1, j+2, k+1, id, jd, kd, 0);
            size_t ip5 = _F_IDX_S3D(i+1, j+1, k+2, id, jd, kd, 0);
            size_t ip6 = _F_IDX_S3D(i+2, j+1, k+2, id, jd, kd, 0);
            size_t ip7 = _F_IDX_S3D(i+2, j+2, k+2, id, jd, kd, 0);
            size_t ip8 = _F_IDX_S3D(i+1, j+2, k+2, id, jd, kd, 0);
            d[ip1+dsize*ivar]=d[ip1+dsize*ivar]+ddd;
            d[ip2+dsize*ivar]=d[ip2+dsize*ivar]+ddd;
            d[ip3+dsize*ivar]=d[ip3+dsize*ivar]+ddd;
            d[ip4+dsize*ivar]=d[ip4+dsize*ivar]+ddd;
            d[ip5+dsize*ivar]=d[ip5+dsize*ivar]+ddd;
            d[ip6+dsize*ivar]=d[ip6+dsize*ivar]+ddd;
            d[ip7+dsize*ivar]=d[ip7+dsize*ivar]+ddd;
            d[ip8+dsize*ivar]=d[ip8+dsize*ivar]+ddd;
          }
        }
      }
      
      //内部の格子点のデータを8で割る
      VolumeDataDivideBy8(&d[dsize*ivar], id, jd, kd);
      
      //面上の格子点のデータを4で割る
      FaceDataDivideBy4(&d[dsize*ivar], id, jd, kd);
      
      //辺上の格子点のデータを2で割る
      LineDataDivideBy2(&d[dsize*ivar], id, jd, kd);
      
      //境界条件処理
      
      
    }//loop ivar
  
  };


  /**
   * @brief 成分別Vectorの格子点での値をセット
   * @param [out]    d        格子点data
   * @param [in]     data     セル中心data
   * @param [in]     id       セル中心data x方向サイズ
   * @param [in]     jd       セル中心data y方向サイズ
   * @param [in]     kd       セル中心data z方向サイズ
   * @param [in]     ivar     ベクトル成分 =0:x =1:y =2:z
   */
  template<class T1, class T2>
  void setVectorComponentGridData(T1* d, T2* data, int id, int jd, int kd, int ivar, int gc_out)
  {
    int ix = size[0];
    int jx = size[1];
    int kx = size[2];
    int gd = guide;
    
    size_t mip;
    T1 ddd;
    int i,j,k;
 
    size_t dsize = (size_t)(id*jd*kd);
 
    for (size_t l=0; l<dsize; l++) d[l]=0.0;
    
    for (int km=1-gc_out; km<=kx+gc_out; km++) {
      for (int jm=1-gc_out; jm<=jx+gc_out; jm++) {
        for (int im=1-gc_out; im<=ix+gc_out; im++) {
          //mip = _F_IDX_V3D(im, jm, km, ivar, ix, jx, kx, gd); //(i,j,k,3)
          mip = _F_IDX_V3DEX(ivar, im, jm, km, ix, jx, kx, gd); //(3,i,j,k)
          ddd=(T1)data[mip];
          i=im-1+gc_out;
          j=jm-1+gc_out;
          k=km-1+gc_out;
          size_t ip1 = _F_IDX_S3D(i+1, j+1, k+1, id, jd, kd, 0);
          size_t ip2 = _F_IDX_S3D(i+2, j+1, k+1, id, jd, kd, 0);
          size_t ip3 = _F_IDX_S3D(i+2, j+2, k+1, id, jd, kd, 0);
          size_t ip4 = _F_IDX_S3D(i+1, j+2, k+1, id, jd, kd, 0);
          size_t ip5 = _F_IDX_S3D(i+1, j+1, k+2, id, jd, kd, 0);
          size_t ip6 = _F_IDX_S3D(i+2, j+1, k+2, id, jd, kd, 0);
          size_t ip7 = _F_IDX_S3D(i+2, j+2, k+2, id, jd, kd, 0);
          size_t ip8 = _F_IDX_S3D(i+1, j+2, k+2, id, jd, kd, 0);
          d[ip1]=d[ip1]+ddd;
          d[ip2]=d[ip2]+ddd;
          d[ip3]=d[ip3]+ddd;
          d[ip4]=d[ip4]+ddd;
          d[ip5]=d[ip5]+ddd;
          d[ip6]=d[ip6]+ddd;
          d[ip7]=d[ip7]+ddd;
          d[ip8]=d[ip8]+ddd;
        }
      }
    }
    
    //内部の格子点のデータを8で割る
    VolumeDataDivideBy8(d, id, jd, kd);
    
    //面上の格子点のデータを4で割る
    FaceDataDivideBy4(d, id, jd, kd);
    
    //辺上の格子点のデータを2で割る
    LineDataDivideBy2(d, id, jd, kd);
  
  };
  
public:

  /**
   * @brief 計算結果ファイルの項目（*.nam）出力
   */
  void function_name();
  
  
  /**
   * @brief 項目別計算結果ファイルの項目（*.nam）出力
   */
  void function_name_divide();
  
  
  /**
   * @brief 境界面定義ファイル（*.fvbnd）出力（未整備）
   * @note BCMになるとりメッシュされたときに対応できないため出力することはない？
   */
  void fvbnd();
  
  
  /**
   * @brief PLOT3Dファイル入出力に関するパラメータ取得
   * @param [in] tpCntl  Texparserポインタ
   */
  void getParameter(TextParser* tpCntl);
  
  
  /**
   * @brief PLOT3Dファイルのポスト出力
   * @param [in]     CurrentStep 現在のステップ数
   * @param [in]     CurrentTime 現在の時刻
   * @param [in]     v00         格子速度
   * @param [in]     origin      基点座標
   * @param [in]     pitch       格子幅
   * @param [in]     dfi_mng     dfi file for plot3d counter
   * @param [in,out] flop        浮動小数点演算数
   */
  void post(const unsigned CurrentStep,
            const double CurrentTime,
            REAL_TYPE* v00,
            const REAL_TYPE* origin,
            const REAL_TYPE* pitch,
            int& dfi_mng,
            double& flop);
  
  
  /**
   * @brief PLOT3Dのパラメータを表示
   * @param [in] fp   file pointer
   */
  void printParameters(FILE* fp);
  
  
  /**
   * @brief 形状データファイル（*.xyz）の出力
   * @param [in]     CurrentStep 現在のステップ数
   * @param [in]     origin      基点座標
   * @param [in]     pitch       格子幅
   */
  void xyz(const unsigned CurrentStep, const REAL_TYPE* origin, const REAL_TYPE* pitch);
  
  
  /**
   * @brief グリッド数、出力項目数をセット（BCMへ移行する場合要編集）
   */
  void setValuePlot3D();

  
  /**
   * @brief メンバ変数の初期化
   * @param [in]     m_size      配列サイズ
   * @param [in]     m_guide     ガイドセルサイズ
   * @param [in]     m_deltaX    格子幅
   * @param [in]     m_C         Controlクラス
   * @param [in]     m_FP3DR     FileIO_PLOT3D_READ
   * @param [in]     m_FP3DW     FileIO_PLOT3D_WRITE
   * @param [in]     m_DFI       DFI
   * @param [in]     m_d_ws      出力データコンテナ配列
   * @param [in]     m_d_p       圧力データ
   * @param [in]     m_d_wo      ベクトル出力データコンテナ配列
   * @param [in]     m_d_v       速度データ
   * @param [in]     m_d_ie       温度データ
   * @param [in]     m_d_p0      スカラーデータ
   * @param [in]     m_d_wv      ベクトルデータ
   * @param [in]     m_d_cdf     BCindex C
   * @param [in]     m_d_bcd     BCindex B
   */
  void Initialize(const int* m_size,
                  const int m_guide,
                  const REAL_TYPE m_deltaX,
                  Control* m_C,
                  FileIO_PLOT3D_READ* m_FP3DR,
                  FileIO_PLOT3D_WRITE* m_FP3DW,
                  DFI* m_dfi,
                  REAL_TYPE* m_d_ws,
                  REAL_TYPE* m_d_p,
                  REAL_TYPE* m_d_wo,
                  REAL_TYPE* m_d_v,
                  REAL_TYPE* m_d_ie,
                  REAL_TYPE* m_d_p0,
                  REAL_TYPE* m_d_wv,
                  int*       m_d_cdf,
                  int*       m_d_bcd);
};

#endif // _PLT3D_FFV_PLOT3D_H_
