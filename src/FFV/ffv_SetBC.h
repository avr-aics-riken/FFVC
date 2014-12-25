#ifndef _FFV_SETBC_H_
#define _FFV_SETBC_H_

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
 * @file   ffv_SetBC.h
 * @brief  FFV BC Class Header
 * @author aics
 */

#include "SetBC.h"
#include <math.h>
#include "ffv_Define.h"
#include "ffv_Ffunc.h"
#include "IP_Jet.h"


class SetBC3D : public SetBC {
public:
  
  /** コンストラクタ */
  SetBC3D() {}
  
  /**　デストラクタ */
  virtual ~SetBC3D() {}
  
  
protected:
  
  // コンポーネントの角速度成分を取り出す
  void extractAngularVel(const int n, REAL_TYPE* vec, REAL_TYPE* ctr, const double tm, const REAL_TYPE* v00);
  
  // コンポーネントの速度境界条件の成分を取り出す
  REAL_TYPE extractVelLBC (const int n, REAL_TYPE* vec, const double tm, const REAL_TYPE* v00);
  
  
  // 外部境界条件の速度境界条件の成分を取り出す
  REAL_TYPE extractVelOBC (const int n, REAL_TYPE* vec, const double tm, const REAL_TYPE* v00);
  
  
  // 温度一定の境界条件
  void psIbcConstTemp (REAL_TYPE* d_ie, const int* d_bcd, const int n);
  
  
  // 内部領域の熱流束指定境界条件
  REAL_TYPE psIbcHeatflux (REAL_TYPE* d_qbc, const int* d_cdf, const int n);
  
  
  // 発熱境界条件
  REAL_TYPE psIbcHeatGen (REAL_TYPE* d_ie, const int* d_bcd, const int n, const REAL_TYPE dt);
  
  
  // 等温境界条件
  REAL_TYPE psIbcIsoThermal (REAL_TYPE* d_qbc, const int* d_cdf, const int* d_bcd, const int n, const REAL_TYPE* d_ie0);
  
  
  // 内部領域のOutflowの境界条件処理
  REAL_TYPE psIbcOutflow (REAL_TYPE* d_ws, const int* d_cdf, const int n, const REAL_TYPE* d_v, const REAL_TYPE* d_ie, const REAL_TYPE* v00);
  
  
  // 内部領域の速度と温度の指定境界条件
  REAL_TYPE psIbcSpecVH (REAL_TYPE* d_ws, const int* d_cdf, const int n, const REAL_TYPE v00, const REAL_TYPE* vec);
  
  
  // 熱伝達境界条件タイプS
  REAL_TYPE psIbcTransferS (REAL_TYPE* d_qbc, const int* d_cdf, const int* d_bcd, const int n, const REAL_TYPE* d_ie0);
  
  
  // 熱伝達境界条件タイプSF
  REAL_TYPE psIbcTransferSF (REAL_TYPE* d_qbc, const int* d_cdf, const int* d_bcd, const int n, const REAL_TYPE* d_ie0);
  
  
  // 熱伝達境界条件タイプSN
  REAL_TYPE psIbcTransferSN (REAL_TYPE* d_qbc, const int* d_cdf, const int* d_bcd, const int n, const REAL_TYPE* d_ie0);
  
  
  // 外部領域のOutflow, In_out, TractionFreeの境界条件処理
  REAL_TYPE psObcFree (REAL_TYPE* d_ws, const int* d_bcd, const int face, const REAL_TYPE* d_vf, const REAL_TYPE* d_ie, const REAL_TYPE* v00);
  
  
  // 外部領域の熱流束指定の境界条件処理
  REAL_TYPE psObcHeatflux (REAL_TYPE* d_qbc, const int face);
  
  
  // 外部領域の熱伝達境界の境界条件処理 (TypeS)
  REAL_TYPE psObcHeatTransferS (REAL_TYPE* d_qbc, const int face, const int* d_bcd, REAL_TYPE* d_ie, const REAL_TYPE* d_ie0);
  
  
  // 外部領域の熱伝達境界の境界条件処理 (Type SF)
  REAL_TYPE psObcHeatTransferSF (REAL_TYPE* d_qbc, const int face, const int* d_bcd, REAL_TYPE* d_ie, const REAL_TYPE* d_ie0);
  
  
  // 外部領域の熱伝達境界の境界条件処理 (Type SN)
  REAL_TYPE psObcHeatTransferSN (REAL_TYPE* d_qbc, const int face, const int* d_bcd, REAL_TYPE* d_ie, const REAL_TYPE* d_ie0);
  
  
  // 外部領域の等温熱流束の境界条件処理
  REAL_TYPE psObcIsoThermal (REAL_TYPE* d_qbc, const int face, const int* d_bcd, REAL_TYPE* d_ie, const REAL_TYPE* d_ie0);
  
  
  // 外部領域の速度指定の境界条件処理
  REAL_TYPE psObcSpecVH (REAL_TYPE* d_ws, const int* d_cdf, const int face, const double tm, const REAL_TYPE* v00);
  
  
  void Pibc_Prdc                (REAL_TYPE* d_p, int* st, int* ed, int* d_bcd, int odr, int dir, REAL_TYPE pv);
  
  
  // 圧力の外部周期境界条件（圧力差）
  void PobcPeriodicDirectional (REAL_TYPE* d_p, const int* ens);
  
  
  // 圧力の外部周期境界条件（単純コピー）
  void PobcPeriodicSimple (REAL_TYPE* d_p, const int* ens);
  
  
  //温度の外部周期境界条件（単純コピー）
  void TobcPeriodicSimple (REAL_TYPE* d_ie, const int* ens);
  
  
  void Vibc_Prdc (REAL_TYPE* d_v, int* st, int* ed, int* d_bd, int odr, int dir);
  
  
  // 速度の外部周期境界条件（単純コピー）
  void VobcPeriodicSimple (REAL_TYPE* d_v, const int* ens);

  
  
  
public:
  
  /**
   * @brief ドライバ指定のチェック
   * @param [in] fp
   * @note コンポーネントと外部境界で指定された，方向と位置の情報が一致するかをチェック
   */
  void checkDriver (FILE* fp);
  
  
  void InnerPBCperiodic (REAL_TYPE* d_p, int* d_bcd);
  
  
  /**
   * @brief 拡散部分に関する内部エネルギーの内部境界処理
   * @param [in,out] d_qbc  境界条件熱流束
   * @param [in]     d_cdf  BCindex C
   * @param [in]     d_bcd  BCindex B
   * @param [in,out] d_ie   n+1時刻の内部エネルギー
   * @param [in]     d_ie0  n時刻の内部エネルギー
   */
  void InnerTBCface (REAL_TYPE* d_qbc, const int* d_cdf, const int* d_bcd, REAL_TYPE* d_ie, const REAL_TYPE* d_ie0);
  
  
  /**
   * @brief セルに対する温度の内部境界をセットする
   * @param [in,out] d_ie  内部エネルギー
   * @param [in]     d_bcd BCindex B
   * @param [in]     dt    時間積分幅
   */
  void InnerTBCvol (REAL_TYPE* d_ie, const int* d_bcd, const REAL_TYPE dt);
  
  
  // 速度の内部周期境界条件
  void InnerVBCperiodic (REAL_TYPE* d_v, int* d_bd);
  
  
  /**
   * @brief 速度境界条件による速度の発散の修正ほか
   * @param [in,out] dv     \sum{u}
   * @param [in]     d_cdf  BCindex C
   * @param [in]     tm     無次元時刻
   * @param [in]     C      Controlクラス
   * @param [in]     v00    基準速度
   * @param [in,out] vf     セルフェイス速度 u^{n+1}
   * @param [in,out] v      セルセンター速度 u^{n+1}
   * @param [in]     avr    平均値計算のテンポラリ値
   * @param [in]     flop   flop count
   * @note 外部境界面のdiv(u)の修正時に領域境界の流量などのモニタ値を計算し，BoundaryOuterクラスに保持 > 反復後にDomainMonitor()で集約
   *       avr[]のインデクスに注意 (Fortran <-> C)
   */
  void modDivergence (REAL_TYPE* dv,
                      int* d_cdf,
                      double tm_d,
                      Control* C,
                      REAL_TYPE* v00,
                      REAL_TYPE* vf,
                      REAL_TYPE* v,
                      Gemini_R* avr,
                      double& flop);
  
  
  void mod_Dir_Forcing (REAL_TYPE* d_v, int* d_bd, REAL_TYPE* d_cvf, REAL_TYPE* v00, double& flop);
  
  
  /**
   * @brief 速度境界条件によるPoisosn式のソース項の修正
   * @param [in,out] dv    \sum{u^*}
   * @param [in]     d_cdf BCindex C
   * @param [in]     tm    無次元時刻
   * @param [in]     C     Control class
   * @param [in]     v00   基準速度
   * @param [in]     vf    セルフェイス速度 u^n
   * @param [in]     vc    セルセンタ疑似速度 u^*
   * @param [in]     v0    セルセンタ速度 u^n
   * @param [in]     dt    時間積分幅
   * @param [in,out] flop  flop count
   */
  void modPsrcVBC (REAL_TYPE* dv,
                   int* d_cdf,
                   const double tm,
                   Control* C,
                   REAL_TYPE* v00,
                   REAL_TYPE* vf,
                   REAL_TYPE* vc,
                   REAL_TYPE* v0,
                   REAL_TYPE dt,
                   double& flop);
  
  
  /**
   @brief 圧力損失部によるPoisosn式のソース項の修正とワーク用の速度を保持
   @param [in,out] s_1     Poisson方程式の反復ソース項
   @param [in]     v       速度ベクトル
   @param [in]     bd      BCindex B
   @param [in]     cvf     コンポーネントの体積率
   @param [in]     v00     参照速度
   @param [in]     c_array コンポーネントワーク配列の管理ポインタ
   @param [out]    flop    flop count
   */
  void mod_Psrc_Forcing (REAL_TYPE* s_1,
                         REAL_TYPE* v,
                         int* bd,
                         REAL_TYPE* cvf,
                         REAL_TYPE* v00,
                         REAL_TYPE** c_array,
                         double& flop);
  
  
  /**
   * @brief 速度境界条件による流束の修正
   * @param [in,out] wv     疑似速度ベクトル u^*
   * @param [in]     v      セルセンター速度ベクトル u^n
   * @param [in]     d_cdf  BCindex C
   * @param [in]     tm     無次元時刻
   * @param [in]     C      Control class
   * @param [in]     v00    基準速度
   * @param [in,out] flop   flop count
   */
  void modPvecFlux (REAL_TYPE* wv,
                    REAL_TYPE* v,
                    int* d_cdf,
                    const double tm,
                    Control* C,
                    REAL_TYPE* v00,
                    double& flop);
  
  void mod_Pvec_Forcing (REAL_TYPE* d_vc, REAL_TYPE* d_v, int* d_bd, REAL_TYPE* d_cvf, REAL_TYPE* v00, REAL_TYPE dt, double& flop);
  
  /**
   * @brief 圧力損失部によるセルセンタ速度の修正と速度の発散値の修正
   * @param [in,out] v          セルセンターの速度
   * @param [in]     bd         BCindex B
   * @param [in]     cvf        コンポーネントの体積率
   * @param [in]     dv         \sum{u}
   * @param [in]     dt         時間積分幅
   * @param [in]     v00        参照速度
   * @param [in]     am         平均速度と圧力損失
   * @param [in]     c_array    コンポーネントワーク配列の管理ポインタ
   * @param [out]    flop       flop count
   * @note am[]のインデクスに注意 (Fortran <-> C)
   */
  void mod_Vdiv_Forcing (REAL_TYPE* v,
                         int* bd,
                         REAL_TYPE* cvf,
                         REAL_TYPE* dv,
                         REAL_TYPE dt,
                         REAL_TYPE* v00,
                         Gemini_R* am,
                         REAL_TYPE** c_array,
                         double& flop);
  
  
  /**
   * @brief 圧力の外部境界条件
   * @param [in,out] d_p  圧力のデータクラス
   * @param [in]     ens  周期境界方向フラグ
   */
  void OuterPBC (REAL_TYPE* d_p, const int* ens);
  
  
  /**
   * @brief 拡散項計算時の温度の外部部境界処理
   * @param [in,out] d_qbc 境界条件熱流束
   * @param [out]    d_ie  n+1時刻の内部エネルギー
   * @param [in]     d_ie0 n時刻の内部エネルギー
   * @param [in]     d_bcd BCindex B
   * @param [in]     C     コントロールクラス
   * @note 断熱BCは断熱マスクで処理
   */
  void OuterTBCdiffusion (REAL_TYPE* d_qbc, REAL_TYPE* d_ie, const REAL_TYPE* d_ie0, const int* d_bcd, Control* C);
  
  
  /**
   * @brief 温度の外部周期境界条件処理
   * @param [in,out] d_ie 内部エネルギー
   * @param [in]     ens  周期境界方向フラグ
   * @see OuterTBConvection(), OuterTBCdiffusion()
   * @note OBC_SYMMETRICは，断熱マスクで処理するため，不要
   */
  void OuterTBCperiodic (REAL_TYPE* d_ie, const int* ens);
  
  
  /**
   * @brief 速度の外部境界条件処理（Div反復内で値を指定する境界条件）
   * @param [in,out] d_v   セルセンター速度ベクトル v^{n+1}
   * @param [in]     d_vf  セルフェイス速度ベクトル v^{n+1}
   * @param [in]     d_cdf BCindex C
   * @param [in]     tm    時刻
   * @param [in]     C     コントロールクラス
   * @param [in]     v00   参照速度
   * @param [in]     ens  周期境界方向フラグ
   */
  void OuterVBC (REAL_TYPE* d_v,
                 REAL_TYPE* d_vf,
                 int* d_cdf,
                 const double tm,
                 Control* C,
                 REAL_TYPE* v00,
                 const int* ens);
  
  
  /**
   * @brief 速度の外部境界条件処理（セルフェイス）の準備
   * @param [in] d_vc   疑似速度ベクトル
   * @param [in] d_v    セルセンター速度ベクトル v^{n}
   * @param [in] d_cdf  BCindex C
   * @param [in] dt     時間積分幅
   * @param [in] C      コントロールクラス
   * @param [in] ens    周期境界方向フラグ
   * @param [in] m_step セッションステップ数
   */
  void OuterVBCfacePrep (REAL_TYPE* d_vc,
                         REAL_TYPE* d_v,
                         int* d_cdf,
                         REAL_TYPE dt,
                         Control* C,
                         const int* ens,
                         const unsigned m_step);
  
  
  /**
   * @brief 周期境界の場合のインデクスの同期
   * @param [in,out] d_bx  BCindexのデータクラス
   * @param [in]     ens   周期境界方向フラグ
   */
  void setBCIperiodic (int* d_bx, const int* ens);
  
  
  /**
   * @brief 初期温度を代入
   * @param [in]     n    エントリ
   * @param [in]     d_bx BCindex B
   * @param [in,out] d_ie 内部エネルギー
   */
  void setInitialTempCompo (const int n, const int* d_bx, REAL_TYPE* d_ie);
  
  
  REAL_TYPE setDirectForcing (REAL_TYPE* d_v, int* d_bx, int n, REAL_TYPE v00);
  
  
  /**
   * @brief 対流項計算時の流束型の温度の外部境界条件処理
   * @param [in,out] d_ws  温度の増分
   * @param [in]     d_cdf BCindex C
   * @param [in]     d_vf  n+1時刻のセルフェイス速度
   * @param [in]     d_ie0 n時刻の内部エネルギー
   * @param [in]     tm    時刻
   * @param [in,out] C     コントロールクラス
   * @param [in]     v00   参照速度
   * @note 無次元熱量はvalに保存
   */
  void TBCconvection (REAL_TYPE* d_ws, const int* d_cdf, const REAL_TYPE* d_vf, const REAL_TYPE* d_ie0, const double tm, Control* C, const REAL_TYPE* v00);
};

#endif // _FFV_SETBC_H_
