#ifndef _FFV_SETBC_H_
#define _FFV_SETBC_H_

// #################################################################
//
// FFV : Frontflow / violet Cartesian
//
// Copyright (c) 2012-2013  All right reserved.
//
// Institute of Industrial Science, University of Tokyo, Japan. 
//
// #################################################################

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
  
  // コンポーネントから速度境界条件の成分を取り出す
  REAL_TYPE extractVel_IBC (const int n, REAL_TYPE* vec, const double tm, const REAL_TYPE* v00, double& flop);
  
  
  // 外部境界条件リストから速度境界条件の成分を取り出す
  REAL_TYPE extractVel_OBC (const int n, REAL_TYPE* vec, const double tm, const REAL_TYPE* v00, double& flop);
  
  /**
   * @brief 温度一定の境界条件
   * @param [in,out] d_t 温度場
   * @param [in]     bh2 BCindex H2
   * @param [in]     n   境界条件コンポーネントのエントリ番号
   */
  void ps_IBC_ConstTemp (REAL_TYPE* d_t, const int* d_bh2, const int n);
  
  
  REAL_TYPE ps_IBC_Heatflux       (REAL_TYPE* d_qbc, int* d_bh1, int n, double& flop);
  REAL_TYPE ps_IBC_HeatGen_SM     (REAL_TYPE* d_t,   int* d_bh2, int n, REAL_TYPE dt, double& flop);
  REAL_TYPE ps_IBC_IsoThermal_SM  (REAL_TYPE* d_qbc, int* d_bh1, int n, REAL_TYPE* d_t, REAL_TYPE* d_t0, double& flop);
  
  
  /**
   * @brief 内部領域のOutflowの境界条件処理
   * @retval 熱量[W]
   * @param [in,out] d_ws  温度増分
   * @param [in]     d_bh1 BCindex H1
   * @param [in]     n     コンポーネントリストのインデクス
   * @param [in]     d_v   速度
   * @param [in]     d_t   温度
   * @param [in]     v00   参照速度
   * @param [in,out] flop  浮動小数演算数
   * @note モニタ量の熱量va(W)は系に対する流入量なので，基準温度に対する熱量
   */
  REAL_TYPE ps_IBC_Outflow (REAL_TYPE* d_ws, const int* d_bh1, const int n, const REAL_TYPE* d_v, const REAL_TYPE* d_t, const REAL_TYPE* v00, double& flop);
  
  
  /**
   * @brief 内部領域の速度と温度の指定境界条件
   * @retval 熱量[-]
   * @param [in,out] d_ws   温度増分
   * @param [in]     d_bh1  BCindex H1
   * @param [in]     n      コンポーネントリストのインデクス
   * @param [in]     v00    参照速度
   * @param [in]     vec    指定ベクトル
   * @param [in,out] flop   浮動小数演算数
   * @note モニタ量の熱量va(無次元)は系に対する流入量なので，基準温度に対する熱量
   */
  REAL_TYPE ps_IBC_SpecVH (REAL_TYPE* d_ws, const int* d_bh1, const int n, const REAL_TYPE v00, const REAL_TYPE* vec, double& flop);
  
  
  REAL_TYPE ps_IBC_Transfer_B_SM  (REAL_TYPE* d_qbc, int* d_bh1, int n, REAL_TYPE* d_t, REAL_TYPE* d_t0, double& flop);
  REAL_TYPE ps_IBC_Transfer_S_SM  (REAL_TYPE* d_qbc, int* d_bh1, int n, REAL_TYPE* d_t, REAL_TYPE* d_t0, double& flop);
  REAL_TYPE ps_IBC_Transfer_SF_SM (REAL_TYPE* d_qbc, int* d_bh1, int n, REAL_TYPE* d_t, REAL_TYPE* d_t0, double& flop);
  
  
  /**
   * @brief 単媒質の場合の熱伝達温境界条件タイプSN（自然対流）　流体のみを解く場合
   * @retval 熱流束の和 (W/m^2)
   * @param [out]    d_qbc  境界条件熱流束
   * @param [in]     d_bh1  BCindex H1
   * @param [in]     n      境界条件コンポーネントのエントリ番号
   * @param [in,out] d_t    n+1時刻の温度場
   * @param [in]     d_t0   n時刻の温度場
   * @note
   *    - 熱流束は加算（他の条件と合成）
   *    - pセルは流体セル
   */
  REAL_TYPE ps_IBC_Transfer_SN_SM (REAL_TYPE* d_qbc, const int* d_bh1, const int n, REAL_TYPE* d_t, const REAL_TYPE* d_t0);
  

  
  
  REAL_TYPE ps_OBC_Free           (REAL_TYPE* d_ws,  int* d_bh1, const int face, REAL_TYPE* d_v, REAL_TYPE* d_t, REAL_TYPE* v00, double& flop);
  REAL_TYPE ps_OBC_Heatflux       (REAL_TYPE* d_qbc, int* d_bh1, const int face, double& flop);
  
  REAL_TYPE ps_OBC_SpecVH (REAL_TYPE* d_ws,  int* d_bh1, const int face, REAL_TYPE* d_t, const double tm, REAL_TYPE* v00, double& flop);
  
  REAL_TYPE ps_OBC_HeatTransfer_BS(REAL_TYPE* d_qbc, int* d_bh1, const int face, REAL_TYPE* d_t, REAL_TYPE* d_t0, double& flop);
  REAL_TYPE ps_OBC_HeatTransfer_SF(REAL_TYPE* d_qbc, int* d_bh1, const int face, REAL_TYPE* d_t, REAL_TYPE* d_t0, double& flop);
  REAL_TYPE ps_OBC_HeatTransfer_SN(REAL_TYPE* d_qbc, int* d_bh1, const int face, REAL_TYPE* d_t, REAL_TYPE* d_t0, double& flop);
  REAL_TYPE ps_OBC_IsoThermal     (REAL_TYPE* d_qbc, int* d_bh1, const int face, REAL_TYPE* d_t, REAL_TYPE* d_t0, double& flop);
  
  void Pibc_Prdc                (REAL_TYPE* d_p, int* st, int* ed, int* d_bcd, int odr, int dir, REAL_TYPE pv);
  void Pobc_Prdc_Directional    (REAL_TYPE* d_p, const int face, REAL_TYPE pv, int uod);
  void Pobc_Prdc_Simple         (REAL_TYPE* d_p, const int face);
  
  
  /**
   * @brief 温度の外部周期境界条件（単純コピー）
   * @param [in] t    温度のデータクラス
   * @param [in] face 面番号
   */
  void Tobc_Prdc_Simple (REAL_TYPE* d_t, const int face);
  
  
  void Vibc_Prdc        (REAL_TYPE* d_v, int* st, int* ed, int* d_bd, int odr, int dir);
  
  
  // 速度の外部周期境界条件（単純なコピー）
  void Vobc_Prdc (REAL_TYPE* d_v, const int face);

  
  
  
public:
  void assign_Temp (REAL_TYPE* d_t, int* d_bh, const double tm, Control* C);
  
  // 境界条件速度をセットする
  void assignVelocity (REAL_TYPE* d_v, int* d_bv, const double tm, REAL_TYPE* v00, bool clear=false);
  
  void checkDriver          (FILE* fp);
  void InnerPBC_Periodic    (REAL_TYPE* d_p, int* d_bcd);
  void InnerTBCface         (REAL_TYPE* d_qbc, int* d_bx, REAL_TYPE* d_t, REAL_TYPE* d_t0, double& flop);
  void InnerTBCvol          (REAL_TYPE* d_t, int* d_bx, REAL_TYPE dt, double& flop);
  
  // 速度の内部境界条件
  void InnerVBC (REAL_TYPE* d_v, int* d_bv, const double tm, REAL_TYPE* v00, double& flop);
  
  // 速度の内部周期境界条件
  void InnerVBC_Periodic (REAL_TYPE* d_v, int* d_bd);
  
  
  // 速度境界条件による速度の発散の修正ほか
  void mod_div (REAL_TYPE* dv,
                int* bv,
                double tm_d,
                REAL_TYPE* v00,
                Gemini_R* avr,
                REAL_TYPE* vf,
                REAL_TYPE* v,
                Control* C,
                double& flop);
  
  
  void mod_Dir_Forcing (REAL_TYPE* d_v, int* d_bd, float* d_cvf, REAL_TYPE* v00, double& flop);
  
  
  // 速度境界条件によるPoisosn式のソース項の修正
  void mod_Psrc_VBC (REAL_TYPE* s_0,
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
  
  
  // 速度境界条件による流束の修正
  void mod_Pvec_Flux (REAL_TYPE* wv,
                      REAL_TYPE* v,
                      int* bv,
                      const double tm,
                      Control* C,
                      int v_mode,
                      REAL_TYPE* v00,
                      double& flop);
  
  void mod_Pvec_Forcing     (REAL_TYPE* d_vc, REAL_TYPE* d_v, int* d_bd, float* d_cvf, REAL_TYPE* v00, REAL_TYPE dt, double& flop);
  
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
  
  // 圧力の外部境界条件
  void OuterPBC (REAL_TYPE* d_p);
  
  // 温度の外部境界条件
  void OuterTBC (REAL_TYPE* d_t);
  
  // 外部境界セルフェイスに対する温度条件
  void OuterTBCface (REAL_TYPE* d_qbc, int* d_bx, REAL_TYPE* d_t, REAL_TYPE* d_t0, Control* C, double& flop);
  
  // 速度の外部境界条件処理（VP反復内で値を指定する境界条件）
  void OuterVBC (REAL_TYPE* d_v, REAL_TYPE* d_vf, int* d_bv, const double tm, Control* C, REAL_TYPE* v00, double& flop);
  
  
  // 速度の外部境界処理(タイムステップに一度ガイドセルに値を設定する)
  void OuterVBC_GC (REAL_TYPE* d_v, int* d_bv, const double tm, const Control* C, const REAL_TYPE* v00, double& flop);
  
  
  // 疑似速度の外部境界条件処理
  void OuterVBC_Pseudo (REAL_TYPE* d_vc, int* d_bv, Control* C, double& flop);
  
  
  void ps_BC_Convection (REAL_TYPE* d_ws, int* d_bh1, REAL_TYPE* d_v, REAL_TYPE* d_t, const double tm, Control* C, REAL_TYPE* v00, double& flop);
  
  
  // 周期境界の場合のインデクスの同期
  void setBCIperiodic (int* d_bx);
  
  
  void setInitialTemp_Compo (int n, int* d_bx, REAL_TYPE* d_t);
  
  REAL_TYPE setDirectForcing (REAL_TYPE* d_v, int* d_bx, int n, REAL_TYPE v00);
  
};

#endif // _FFV_SETBC_H_
