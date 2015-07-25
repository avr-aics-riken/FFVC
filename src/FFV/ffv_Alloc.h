#ifndef _FFV_ALLOC_H_
#define _FFV_ALLOC_H_

//##################################################################################
//
// FFV-C : Frontflow / violet Cartesian
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2015 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2015 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################
//

/**
 * @file   ffv_Alloc.h
 * @brief  FFV array allocation Class Header
 * @author aics
 */

#include "DomainInfo.h"
#include "Control.h"
#include "ffv_Define.h"


using namespace std;

class FALLOC : public DomainInfo {

private:
  double array_size;
  
  
public:
  
  // [*]印はタイムステップ間でデータを保持しておかなければならない配列
  // それ以外は，ファイル入出力のバッファに利用可能
  //
  
  // Prep
  REAL_TYPE *d_ws;  ///<     反復中に固定のソース, データ変換のワーク, file IO用のバッファには使わない
  int *d_mid;       ///<     ID配列, IBLANK配列としても利用
  int *d_bcd;       ///< [*] BCindex ID
  int *d_bcp;       ///< [*] BCindex P
  int *d_cdf;       ///< [*] BCindex C
  REAL_TYPE *d_pvf; ///< [*] セル体積率
  
  
  // Main
  REAL_TYPE *d_v;   ///< [*] セルセンター速度
  REAL_TYPE *d_vf;  ///< [*] セルフェイス速度
  REAL_TYPE *d_p;   ///< [*] 圧力
  REAL_TYPE *d_dv;  ///< [*] \sum{u}の保存
  REAL_TYPE *d_ie;  ///< [*] 内部エネルギー
  
  REAL_TYPE *d_wv;  ///<     ワーク配列, file IO用のバッファには使わない
  
  REAL_TYPE *d_io_buffer; ///< 大きなバッファ　以下と共用
  // >> file IO用のバッファと共用
  REAL_TYPE *d_vc;  ///<     セルセンター疑似速度
  REAL_TYPE *d_v0;  ///<     n-stepの速度保持
  REAL_TYPE *d_p0;  ///<     圧力（1ステップ前）
  REAL_TYPE *d_sq;  ///<     反復中に変化するソース
  REAL_TYPE *d_b;   ///<     Ax=bの右辺ベクトル
  REAL_TYPE *d_ie0; ///<     内部エネルギー（1ステップ前）
  REAL_TYPE *d_qbc; ///<     熱BC flux保持
  // << file IO用のバッファと共用
  
  
  // 渦度関連オプション
  REAL_TYPE *d_vrt; ///< [*] 渦度ベクトル
  
  // Polygon
  long long  *d_cut;    ///< [*] 距離情報
  int        *d_bid;    ///< [*] BC
  
  
  // 平均値
  REAL_TYPE *d_ap;         ///< [*] 圧力（時間平均値）
  REAL_TYPE *d_av;         ///< [*] 速度（時間平均値）
  REAL_TYPE *d_ae;         ///< [*] 内部エネルギー（時間平均値）

  
  // Adams-bashforth
  REAL_TYPE *d_abf; ///<     ワーク

  
  // 疎格子ロード時のワーク
  REAL_TYPE *d_r_v; ///<     粗格子の速度
  REAL_TYPE *d_r_p; ///<     粗格子の圧力
  REAL_TYPE *d_r_t; ///<     粗格子の温度
  
  
  // Components
  REAL_TYPE** component_array; ///< コンポーネントワーク配列のアドレス管理
  REAL_TYPE *d_cvf; ///< [*] 体積率
  
  
  // LES計算
  REAL_TYPE *d_vt;      ///< [*] 渦粘性係数
  
  // 統計
  REAL_TYPE *d_rms_v;      ///< [*] セルセンター速度の乱流強度
  REAL_TYPE *d_rms_mean_v; ///< [*] セルセンター速度の時間平均乱流強度
  REAL_TYPE *d_rms_p;      ///< [*] 圧力の変動強度
  REAL_TYPE *d_rms_mean_p; ///< [*] 圧力の時間平均変動強度
  REAL_TYPE *d_rms_t;      ///< [*] 温度の変動強度
  REAL_TYPE *d_rms_mean_t; ///< [*] 温度の時間平均変動強度

  
  // 乱流統計量の詳細
  REAL_TYPE *d_R;          ///< [*] レイノルズ応力テンソル
  REAL_TYPE *d_aR;         ///< [*] レイノルズ応力テンソル (時間平均値)
  REAL_TYPE *d_aP;         ///< [*] 生成項 (時間平均値)
  REAL_TYPE *d_aE;         ///< [*] 散逸項 (時間平均値)
  REAL_TYPE *d_aT;         ///< [*] 乱流拡散項 (時間平均値)
  REAL_TYPE *d_aPI;        ///< [*] 速度圧力勾配相関項 (時間平均値)
  
  
  // 界面計算
  REAL_TYPE *d_vof; ///< [*] VOF値
  
  
  
  // 反復法のバッファ
  int cf_sz[3];     ///<     SOR2SMAの反復のバッファサイズ
  REAL_TYPE *cf_x;  ///<     i方向のバッファ
  REAL_TYPE *cf_y;  ///<     j方向のバッファ
  REAL_TYPE *cf_z;  ///<     k方向のバッファ
  
  
  // GMRES
  REAL_TYPE * d_wg;   ///< テンポラリの配列 [size]
  REAL_TYPE * d_res;  ///< 残差 = b - Ax
  REAL_TYPE * d_vm;   ///< Kryolov subspaceの直交基底 [size*FREQ_OF_RESTART]
  REAL_TYPE * d_zm;   ///< Right-hand side vector for the residual minimization problem [size*FREQ_OF_RESTART]
  
  
  // PCG & BiCGstab
  REAL_TYPE *d_pcg_r;
  REAL_TYPE *d_pcg_p;
  
  // PCG
  REAL_TYPE *d_pcg_z;
  
  // BiCGstab
  REAL_TYPE *d_pcg_r0;
  REAL_TYPE *d_pcg_q;
  REAL_TYPE *d_pcg_s;
  REAL_TYPE *d_pcg_t;
  
  // BiCGSTAB with Preconditioning
  REAL_TYPE *d_pcg_p_;
  REAL_TYPE *d_pcg_s_;
  REAL_TYPE *d_pcg_t_;
  
  
  
  

  /** コンストラクタ */
  FALLOC() {
    
    array_size = 0.0;
    cf_sz[0] = cf_sz[1] = cf_sz[2] = 0;
    cf_x = NULL;
    cf_y = NULL;
    cf_z = NULL;
    
    d_io_buffer = NULL;
    
    d_v = NULL;
    d_vf = NULL;
    d_vc = NULL;
    d_v0 = NULL;
    d_wv = NULL;
    d_abf = NULL;
    d_av = NULL;
    d_qbc = NULL;
    
    d_mid = NULL;
    d_bcd = NULL;
    d_bcp = NULL;
    d_cdf = NULL;
    
    d_p = NULL;
    d_p0 = NULL;
    d_ws = NULL;
    d_sq = NULL;
    d_dv = NULL;
    d_b = NULL;
    d_ie = NULL;
    d_ie0 = NULL;
    d_vrt = NULL;
    
    d_vof = NULL;
    d_ap = NULL;
    d_ae = NULL;
    d_pvf = NULL;
    d_cvf = NULL;
    
    d_vt = NULL;
    d_rms_v = NULL;
    d_rms_p = NULL;
    d_rms_t = NULL;
    d_rms_mean_v = NULL;
    d_rms_mean_p = NULL;
    d_rms_mean_t = NULL;
    

    d_R   = NULL;
    d_aR  = NULL;
    d_aP  = NULL;
    d_aE  = NULL;
    d_aT  = NULL;
    d_aPI = NULL;
    
    
    d_r_v = NULL;
    d_r_p = NULL;
    d_r_t = NULL;
    
    d_wg = NULL;
    d_res = NULL;
    d_vm = NULL;
    d_zm = NULL;
    
    d_pcg_r = NULL;
    d_pcg_p = NULL;
    
    d_pcg_z = NULL;
    
    d_pcg_r0 = NULL;
    d_pcg_q = NULL;
    d_pcg_s = NULL;
    d_pcg_t = NULL;
    
    d_pcg_p_ = NULL;
    d_pcg_s_ = NULL;
    d_pcg_t_ = NULL;
    
    
    d_cut = NULL;
    d_bid = NULL;
    
    component_array = NULL;
  };
  
  /**　デストラクタ */
  ~FALLOC() {};
  
  
  
  // Adams-Bashforth法に用いる配列のアロケーション
  void allocArray_AB2(double &total);
  
  
  // 平均処理に用いる配列のアロケーション
  void allocArray_Average(double &total, Control* C);
  
  
  // 粗格子読み込みに用いる配列のアロケーション
  void allocArray_CoarseMesh(const int* r_size, double &prep, const bool isHeat);
  
  
  // コンポーネント体積率の配列のアロケーション
  void allocArray_CompoVF(double &prep, double &total);
  
  
  // カット情報の配列
  void allocArray_Cut(double &prep, double &total);
  
  
  // コンポーネントのワーク用配列のアロケート
  void allocArray_Forcing(double& m_prep, double& m_total, FILE* fp, Control* C, CompoList* cmp);
  
  
  // 熱の主計算部分に用いる配列のアロケーション
  void allocArray_Heat(double &total);
  
  
  // 体積率の配列のアロケーション
  void allocArray_Interface(double &total);
  
  
  // Krylov-subspace法に用いる配列のアロケーション
  void allocArray_Krylov(double &total);
  
  
  // LES計算に用いる配列のアロケーション
  void allocArray_LES(double &total);
  
  
  // 統計処理に用いる配列のアロケーション
  void allocArray_Statistic(double &total, Control* C);
  
  
  // 主計算部分に用いる配列のアロケーション
  void allocArray_Main(double &total, Control* C);
  
  
  // PCG法に用いる配列のアロケーション
  void allocArray_PCG(double &total);
  
  
  // BiCGSTAB法に用いる配列のアロケーション
  void allocArray_BiCGstab(double &total);
  
  
  // BiCGSTAB /w preconditionning に用いる配列のアロケーション
  void allocArray_BiCGSTABwithPreconditioning(double &total);
  
  
  // 前処理に用いる配列のアロケーション
  void allocArray_Prep(double &prep, double &total);
  
  
  // SOR2SMAのバッファ確保
  void allocate_SOR2SMA_buffer(double &total);

  
  // スカラの配列サイズを計算
  void setArraySize()
  {
    array_size = (double)(
                          (size[0]+2*guide) *
                          (size[1]+2*guide) *
                          (size[2]+2*guide) );
  }
  
};

#endif // _FFV_ALLOC_H_
