#ifndef _FFV_SETBC_H_
#define _FFV_SETBC_H_

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
 * @file   ffv_SetBC.h
 * @brief  FFV BC Class Header
 * @author kero
 */

#include <math.h>

#include "SetBC.h"
#include "ffv_Define.h"
#include "ffv_Ffunc.h"


class SetBC3D : public SetBC {
public:
  
  /** コンストラクタ */
  SetBC3D() {}
  
  /**　デストラクタ */
  virtual ~SetBC3D() {}
  
protected:
  
  /**
   * @brief コンポーネントから速度境界条件の成分を取り出す
   * @param [in]     n    コンポーネントのインデクス
   * @param [out]    vec  ベクトル成分
   * @param [in]     tm   時刻
   * @param [in]     v00  格子速度
   * @param [in,out] flop 浮動小数点演算数
   */
  REAL_TYPE extractVel_IBC (const int n, REAL_TYPE* vec, const REAL_TYPE tm, const REAL_TYPE* v00, double& flop);
  
  /**
   * @brief コンポーネントから速度境界条件の成分を取り出す
   * @param [in]     n    コンポーネントのインデクス
   * @param [out]    vec  ベクトル成分
   * @param [in]     tm   時刻
   * @param [in]     v00  格子速度
   * @param [in,out] flop 浮動小数点演算数
   */
  REAL_TYPE extractVel_OBC (const int n, REAL_TYPE* vec, const REAL_TYPE tm, const REAL_TYPE* v00, double& flop);
  
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
  REAL_TYPE ps_OBC_SpecVH         (REAL_TYPE* d_ws,  int* d_bh1, const int face, REAL_TYPE* d_t, REAL_TYPE tm, REAL_TYPE* v00, double& flop);
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
  void Vobc_Prdc         (REAL_TYPE* d_v, const int face);
  void Vobc_Prdc_CF     (REAL_TYPE* d_v, const int face);
  
public:
  void assign_Temp          (REAL_TYPE* d_t, int* d_bh, REAL_TYPE tm, Control* C);
  void assign_Velocity      (REAL_TYPE* d_v, int* d_bv, REAL_TYPE tm, REAL_TYPE* v00, bool clear=false);
  void checkDriver          (FILE* fp);
  void InnerPBC_Periodic    (REAL_TYPE* d_p, int* d_bcd);
  void InnerTBCface         (REAL_TYPE* d_qbc, int* d_bx, REAL_TYPE* d_t, REAL_TYPE* d_t0, double& flop);
  void InnerTBCvol          (REAL_TYPE* d_t, int* d_bx, REAL_TYPE dt, double& flop);
  void InnerVBC             (REAL_TYPE* d_v, int* d_bv, REAL_TYPE tm, REAL_TYPE* v00, double& flop);
  void InnerVBC_Periodic    (REAL_TYPE* d_v, int* d_bd);
  
  
  // 速度境界条件による速度の発散の修正ほか
  void mod_div (REAL_TYPE* dv,
                int* bv,
                REAL_TYPE tm,
                REAL_TYPE* v00,
                Gemini_R* avr,
                REAL_TYPE* d_vf,
                double& flop);
  
  
  void mod_Dir_Forcing (REAL_TYPE* d_v, int* d_bd, float* d_cvf, REAL_TYPE* v00, double& flop);
  
  
  // 速度境界条件によるPoisosn式のソース項の修正
  void mod_Psrc_VBC (REAL_TYPE* s_0,
                     REAL_TYPE* vc,
                     REAL_TYPE* v0,
                     REAL_TYPE* vf, 
                     int* bv,
                     REAL_TYPE tm,
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
   @brief 速度境界条件による流束の修正
   @param [in,out] wv     疑似速度ベクトル u^*
   @param [in]     v      セルセンター速度ベクトル u^n
   @param [in]     vf     セルフェイス速度ベクトル u^n
   @param [in]     bv     BCindex V
   @param [in]     tm     無次元時刻
   @param [in]     C      Control class
   @param [in]     v_mode 粘性項のモード (0=粘性項を計算しない, 1=粘性項を計算する, 2=壁法則)
   @param [in]     v00    基準速度
   @param [out]    flop   flop count
   */
  void mod_Pvec_Flux (REAL_TYPE* wv,
                      REAL_TYPE* v,
                      REAL_TYPE* vf, 
                      int* bv,
                      REAL_TYPE tm,
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
  
  
  void mod_Vis_EE           (REAL_TYPE* d_vc, REAL_TYPE* d_v0, REAL_TYPE cf, int* d_bx, REAL_TYPE tm, REAL_TYPE dt, REAL_TYPE* v00, double& flop);
  void OuterPBC             (REAL_TYPE* d_p);
  void OuterTBC             (REAL_TYPE* d_t);
  void OuterTBCface         (REAL_TYPE* d_qbc, int* d_bx, REAL_TYPE* d_t, REAL_TYPE* d_t0, Control* C, double& flop);
  void OuterVBC             (REAL_TYPE* d_v, REAL_TYPE* d_vc, int* d_bv, REAL_TYPE tm, REAL_TYPE dt, Control* C, REAL_TYPE* v00, double& flop);
  void OuterVBC_Periodic    (REAL_TYPE* d_v);
  
  void OuterVBC_Pseudo      (REAL_TYPE* d_vc,
                             REAL_TYPE* d_v0,
                             REAL_TYPE tm,
                             REAL_TYPE dt,
                             Control* C,
                             int* d_bv,
                             double& flop);
  
  void ps_BC_Convection     (REAL_TYPE* d_ws, int* d_bh1, REAL_TYPE* d_v, REAL_TYPE* d_t, REAL_TYPE tm, Control* C, REAL_TYPE* v00, double& flop);
  void setBCIperiodic       (int* d_bx);
  void setInitialTemp_Compo (int n, int* d_bx, REAL_TYPE* d_t);
  void updateOuter          (REAL_TYPE* d_v, REAL_TYPE* vc);
  
  REAL_TYPE setDirectForcing (REAL_TYPE* d_v, int* d_bx, int n, REAL_TYPE v00);
  
  /**
   * @brief *.fvbndの書き出し
  void WriteBoundaryPLOT3D(FileIO_PLOT3D_WRITE* FP3DW, std::vector<std::string>& bcname);
   */
  
};

#endif // _FFV_SETBC_H_
