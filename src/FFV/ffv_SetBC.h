#ifndef _FFV_SETBC_H_
#define _FFV_SETBC_H_

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
 * @file   ffv_SetBC.h
 * @brief  FFV BC Class Header
 * @author kero
 */

#include <math.h>

#include "SetBC.h"
#include "ffv_Define.h"
#include "ffv_Ffunc.h"
#include "../IP/IP_Jet.h"


class SetBC3D : public SetBC {
public:
  
  /** コンストラクタ */
  SetBC3D() {}
  
  /**　デストラクタ */
  virtual ~SetBC3D() {}
  
  
protected:
  
  // コンポーネントの速度境界条件の成分を取り出す
  REAL_TYPE extractVelLBC (const int n, REAL_TYPE* vec, const double tm, const REAL_TYPE* v00);
  
  
  // 外部境界条件の速度境界条件の成分を取り出す
  REAL_TYPE extractVelOBC (const int n, REAL_TYPE* vec, const double tm, const REAL_TYPE* v00);
  
  
  // 温度一定の境界条件
  void psIbcConstTemp (REAL_TYPE* d_t, const int* d_bh2, const int n);
  
  
  // 内部領域の熱流束指定境界条件
  REAL_TYPE psIbcHeatflux (REAL_TYPE* d_qbc, const int* d_bh1, const int n);
  
  
  // 発熱境界条件
  REAL_TYPE psIbcHeatGen (REAL_TYPE* d_t, const int* d_bh2, const int n, const REAL_TYPE dt);
  
  
  // 等温境界条件
  REAL_TYPE psIbcIsoThermal (REAL_TYPE* d_qbc, const int* d_bh1, const int n, const REAL_TYPE* d_t0);
  
  
  // 内部領域のOutflowの境界条件処理
  REAL_TYPE psIbcOutflow (REAL_TYPE* d_ws, const int* d_bh1, const int n, const REAL_TYPE* d_v, const REAL_TYPE* d_t, const REAL_TYPE* v00);
  
  
  // 内部領域の速度と温度の指定境界条件
  REAL_TYPE psIbcSpecVH (REAL_TYPE* d_ws, const int* d_bh1, const int n, const REAL_TYPE v00, const REAL_TYPE* vec);
  
  
  // 熱伝達境界条件タイプB
  REAL_TYPE psIbcTransferB (REAL_TYPE* d_qbc, const int* d_bh1, const int n, const REAL_TYPE* d_t0);
  
  
  // 熱伝達境界条件タイプS
  REAL_TYPE psIbcTransferS (REAL_TYPE* d_qbc, const int* d_bh1, const int n, const REAL_TYPE* d_t0);
  
  
  // 熱伝達境界条件タイプSF
  REAL_TYPE psIbcTransferSF (REAL_TYPE* d_qbc, const int* d_bh1, const int n, const REAL_TYPE* d_t0);
  
  
  // 熱伝達境界条件タイプSN
  REAL_TYPE psIbcTransferSN (REAL_TYPE* d_qbc, const int* d_bh1, const int n, const REAL_TYPE* d_t0);
  
  
  // 外部領域のOutflow, In_out, TractionFreeの境界条件処理
  REAL_TYPE psObcFree (REAL_TYPE* d_ws, const int* d_bh, const int face, const REAL_TYPE* d_vf, const REAL_TYPE* d_t, const REAL_TYPE* v00);
  
  
  // 外部領域の熱流束指定の境界条件処理
  REAL_TYPE psObcHeatflux (REAL_TYPE* d_qbc, const int face);
  
  
  // 外部領域の熱伝達境界の境界条件処理 (TypeB/S共用)
  REAL_TYPE psObcHeatTransferBS (REAL_TYPE* d_qbc, const int face, REAL_TYPE* d_t, const REAL_TYPE* d_t0);
  
  
  // 外部領域の熱伝達境界の境界条件処理 (Type SF)
  REAL_TYPE psObcHeatTransferSF (REAL_TYPE* d_qbc, const int face, REAL_TYPE* d_t, const REAL_TYPE* d_t0);
  
  
  // 外部領域の熱伝達境界の境界条件処理 (Type SN)
  REAL_TYPE psObcHeatTransferSN (REAL_TYPE* d_qbc, const int face, REAL_TYPE* d_t, const REAL_TYPE* d_t0);
  
  
  // 外部領域の等温熱流束の境界条件処理
  REAL_TYPE psObcIsoThermal (REAL_TYPE* d_qbc, const int face, REAL_TYPE* d_t, const REAL_TYPE* d_t0);
  
  
  // 外部領域の速度指定の境界条件処理
  REAL_TYPE psObcSpecVH (REAL_TYPE* d_ws, const int* d_bh1, const int face, const double tm, const REAL_TYPE* v00);
  
  
  void Pibc_Prdc                (REAL_TYPE* d_p, int* st, int* ed, int* d_bcd, int odr, int dir, REAL_TYPE pv);
  
  void PobcPeriodicDirectional (REAL_TYPE* d_p, const int face, REAL_TYPE pv, int uod);
  
  // 圧力の外部周期境界条件（単純コピー）
  void PobcPeriodicSimple (REAL_TYPE* d_p, const int face);
  
  
  //温度の外部周期境界条件（単純コピー）
  void TobcPeriodicSimple (REAL_TYPE* d_t, const int face);
  
  
  void Vibc_Prdc (REAL_TYPE* d_v, int* st, int* ed, int* d_bd, int odr, int dir);
  
  
  // 速度の外部周期境界条件（単純コピー）
  void VobcPeriodicSimple (REAL_TYPE* d_v, const int face);

  
  
  
public:
  
  /**
   * @brief 温度指定境界条件に必要な温度をセットする
   * @param [in,out] d_t  温度場
   * @param [in]     bh   BCindex H1
   * @param [in]     tm   無次元時刻
   * @param [in]     C    Control class
   */
  void assignTemp (REAL_TYPE* d_t, int* d_bh1, const double tm, const Control* C);
  
  
  /**
   * @brief 速度指定境界条件に必要な参照速度をセットする
   * @param [in,out] d_v   セルセンタ速度ベクトル (n step)
   * @param [in]     bv    BCindex V
   * @param [in]     tm    無次元時刻
   * @param [in]     v00   参照速度
   * @param [in]     clear trueのとき，出力時に速度を壁面速度にする（デフォルトfalse）
   */
  void assignVelocity (REAL_TYPE* d_v, int* d_bv, const double tm, REAL_TYPE* v00, bool clear=false);
  
  
  /**
   * @brief ドライバ指定のチェック
   * @param [in] fp
   * @note コンポーネントと外部境界で指定された，方向と位置の情報が一致するかをチェック
   */
  void checkDriver (FILE* fp);
  
  
  void InnerPBCperiodic (REAL_TYPE* d_p, int* d_bcd);
  
  
  // 拡散部分に関する温度の内部境界処理
  void InnerTBCface (REAL_TYPE* d_qbc, const int* d_bh1, REAL_TYPE* d_t, const REAL_TYPE* d_t0);
  
  
  /**
   * @brief セルに対する温度の内部境界をセットする
   * @param [in,out] d_t 温度
   * @param [in]     bh2 BCindex H2
   * @param [in]     dt  時間積分幅
   */
  void InnerTBCvol (REAL_TYPE* d_t, const int* d_bh2, const REAL_TYPE dt);
  
  
  // 速度の内部境界条件
  void InnerVBC (REAL_TYPE* d_v, int* d_bv, const double tm, REAL_TYPE* v00);
  
  
  // 速度の内部周期境界条件
  void InnerVBCperiodic (REAL_TYPE* d_v, int* d_bd);
  
  
  // 速度境界条件による速度の発散の修正ほか
  void modDivergence (REAL_TYPE* dv,
                      int* bv,
                      double tm_d,
                      REAL_TYPE* v00,
                      Gemini_R* avr,
                      REAL_TYPE* vf,
                      REAL_TYPE* v,
                      Control* C,
                      double& flop);
  
  
  void mod_Dir_Forcing (REAL_TYPE* d_v, int* d_bd, float* d_cvf, REAL_TYPE* v00, double& flop);
  
  
  /**
   * @brief 速度境界条件によるPoisosn式のソース項の修正
   * @param [in,out] s_0   \sum{u^*}
   * @param [in]     vc    セルセンタ疑似速度 u^*
   * @param [in]     v0    セルセンタ速度 u^n
   * @param [in]     vf    セルフェイス速度 u^n
   * @param [in]     bv    BCindex V
   * @param [in]     tm    無次元時刻
   * @param [in]     dt    時間積分幅
   * @param [in]     C     Control class
   * @param [in]     v00   基準速度
   * @param [in,out] flop  flop count
   */
  void modPsrcVBC (REAL_TYPE* s_0,
                   REAL_TYPE* vc,
                   REAL_TYPE* v0,
                   REAL_TYPE* vf,
                   int* bv,
                   const double tm,
                   REAL_TYPE dt,
                   Control* C,
                   REAL_TYPE* v00,
                   double& flop);
  
  
  /**
   @brief 圧力損失部によるPoisosn式のソース項の修正とワーク用の速度を保持
   @param [in,out] s_1     Poisson方程式の反復ソース項
   @param [in]     v       速度ベクトル
   @param [in]     bd      BCindex ID
   @param [in]     cvf     コンポーネントの体積率
   @param [in]     v00     参照速度
   @param [in]     c_array コンポーネントワーク配列の管理ポインタ
   @param [out]    flop    flop count
   */
  void mod_Psrc_Forcing (REAL_TYPE* s_1,
                         REAL_TYPE* v,
                         int* bd,
                         float* cvf,
                         REAL_TYPE* v00,
                         REAL_TYPE** c_array,
                         double& flop);
  
  
  /**
   * @brief 速度境界条件による流束の修正
   * @param [in,out] wv     疑似速度ベクトル u^*
   * @param [in]     v      セルセンター速度ベクトル u^n
   * @param [in]     bv     BCindex V
   * @param [in]     tm     無次元時刻
   * @param [in]     C      Control class
   * @param [in]     v_mode 粘性項のモード (0=粘性項を計算しない, 1=粘性項を計算する, 2=壁法則)
   * @param [in]     v00    基準速度
   * @param [in,out] flop   flop count
   */
  void modPvecFlux (REAL_TYPE* wv,
                    REAL_TYPE* v,
                    int* bv,
                    const double tm,
                    Control* C,
                    int v_mode,
                    REAL_TYPE* v00,
                    double& flop);
  
  void mod_Pvec_Forcing (REAL_TYPE* d_vc, REAL_TYPE* d_v, int* d_bd, float* d_cvf, REAL_TYPE* v00, REAL_TYPE dt, double& flop);
  
  /**
   * @brief 圧力損失部によるセルセンタ速度の修正と速度の発散値の修正
   * @param [in,out] v          セルセンターの速度
   * @param [in]     bd         BCindex ID
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
                         float* cvf,
                         REAL_TYPE* dv,
                         REAL_TYPE dt,
                         REAL_TYPE* v00,
                         Gemini_R* am,
                         REAL_TYPE** c_array,
                         double& flop);
  
  
  void mod_Vis_EE (REAL_TYPE* d_vc, REAL_TYPE* d_v0, REAL_TYPE cf, int* d_bx, const double tm, REAL_TYPE dt, REAL_TYPE* v00, double& flop);
  
  
  /**
   * @brief 圧力の外部境界条件
   * @param [in,out]d_p 圧力のデータクラス
   */
  void OuterPBC (REAL_TYPE* d_p);
  
  
  /**
   * @brief 対流項計算時の流束型の温度の外部境界条件処理
   * @param [in,out] d_ws  温度の増分
   * @param [in]     d_bh1 BCindex H1
   * @param [in]     d_vf  n+1時刻のセルフェイス速度
   * @param [in]     d_t0  n時刻の温度
   * @param [in]     tm    時刻
   * @param [in,out] C     コントロールクラス
   * @param [in]     v00   参照速度
   * @note 無次元熱量はvalに保存
   */
  void OuterTBCconvection (REAL_TYPE* d_ws, const int* d_bh1, const REAL_TYPE* d_vf, const REAL_TYPE* d_t0, const double tm, Control* C, const REAL_TYPE* v00);
  
  
  /**
   * @brief 拡散項計算時の温度の外部部境界処理
   * @param [in,out] d_qbc 境界条件熱流束
   * @param [out]    d_t   n+1時刻の温度場
   * @param [in]     d_t0  n時刻の温度場
   * @param [in]     C     コントロールクラス
   * @note 断熱BCは断熱マスクで処理
   */
  void OuterTBCdiffusion (REAL_TYPE* d_qbc, REAL_TYPE* d_t, const REAL_TYPE* d_t0, Control* C);
  
  
  /**
   * @brief 温度の外部周期境界条件処理
   * @param [in,out] d_t 温度のデータクラス
   * @see OuterTBConvection(), OuterTBCdiffusion()
   * @note OBC_SYMMETRICは，断熱マスクで処理するため，不要
   */
  void OuterTBCperiodic (REAL_TYPE* d_t);
  
  
  /**
   * @brief 速度の外部境界条件処理（VP反復内で値を指定する境界条件）
   * @param [in,out] d_v  セルセンター速度ベクトル v^{n+1}
   * @param [in]     d_vf セルフェイス速度ベクトル v^{n+1}
   * @param [in]     d_bv BCindex V
   * @param [in]     tm   時刻
   * @param [in]     C    コントロールクラス
   * @param [in]     v00  参照速度
   */
  void OuterVBC (REAL_TYPE* d_v, REAL_TYPE* d_vf, int* d_bv, const double tm, Control* C, REAL_TYPE* v00);
  
  
  /**
   * @brief 疑似速度の外部境界条件処理
   * @param [out]    d_vc   疑似速度ベクトル v^*
   * @param [in]     d_bv   BCindex V
   * @param [in]     C      Control class
   */
  void OuterVBCpseudo (REAL_TYPE* d_vc, int* d_bv, Control* C);
  
  
  // 周期境界の場合のインデクスの同期
  void setBCIperiodic (int* d_bx);
  
  
  /**
   * @brief 初期温度を代入
   * @param [in]     n    エントリ
   * @param [in]     d_bx BCindex ID
   * @param [in,out] d_t  温度
   */
  void setInitialTempCompo (const int n, const int* d_bx, REAL_TYPE* d_t);
  
  
  REAL_TYPE setDirectForcing (REAL_TYPE* d_v, int* d_bx, int n, REAL_TYPE v00);
  
};

#endif // _FFV_SETBC_H_
