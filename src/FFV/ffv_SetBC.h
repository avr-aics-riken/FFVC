#ifndef _FFV_SETBC_H_
#define _FFV_SETBC_H_

// #################################################################
//
// FFV : Frontflow / violet
//
// Copyright (c) All right reserved. 2012
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
//
// #################################################################

/** 
 * @file   ffv_SetBC.h
 * @brief  FFV BC Class Header
 * @author kero
 */

#include <math.h>

#include "../FB/SetBC.h"
#include "ffv_Ffunc.h"
#include "ffv_Define.h"

class SetBC3D : public SetBC {
public:
  
  /** コンストラクタ */
  SetBC3D() {}
  
  /**　デストラクタ */
  virtual ~SetBC3D() {}
  
protected:
  REAL_TYPE extractVel_IBC        (int n, REAL_TYPE* vec, REAL_TYPE tm, REAL_TYPE* v00, double& flop);
  REAL_TYPE extractVel_OBC        (int n, REAL_TYPE* vec, REAL_TYPE tm, REAL_TYPE* v00, double& flop);
  REAL_TYPE ps_IBC_Heatflux       (REAL_TYPE* qbc, int* bh1, int n, double& flop);
  REAL_TYPE ps_IBC_HeatGen_SM     (REAL_TYPE* t,   int* bh2, int n, REAL_TYPE dt, double& flop);
  REAL_TYPE ps_IBC_IsoThermal_SM  (REAL_TYPE* qbc, int* bh1, int n, REAL_TYPE* t, REAL_TYPE* t0, double& flop);
  REAL_TYPE ps_IBC_Transfer_B_SM  (REAL_TYPE* qbc, int* bh1, int n, REAL_TYPE* t, REAL_TYPE* t0, double& flop);
  REAL_TYPE ps_IBC_Transfer_S_SM  (REAL_TYPE* qbc, int* bh1, int n, REAL_TYPE* t, REAL_TYPE* t0, double& flop);
  REAL_TYPE ps_IBC_Transfer_SF_SM (REAL_TYPE* qbc, int* bh1, int n, REAL_TYPE* t, REAL_TYPE* t0, double& flop);
  REAL_TYPE ps_IBC_Transfer_SN_SM (REAL_TYPE* qbc, int* bh1, int n, REAL_TYPE* t, REAL_TYPE* t0, double& flop);
  REAL_TYPE ps_IBC_Outflow        (REAL_TYPE* ws,  int* bh1, int n, REAL_TYPE* v, REAL_TYPE* t, REAL_TYPE* v00, double& flop);
  REAL_TYPE ps_IBC_SpecVH         (REAL_TYPE* ws,  int* bh1, int n, REAL_TYPE v00, REAL_TYPE* vec, double& flop);
  REAL_TYPE ps_OBC_Free           (REAL_TYPE* ws,  int* bh1, int face, REAL_TYPE* v, REAL_TYPE* t, REAL_TYPE* v00, double& flop);
  REAL_TYPE ps_OBC_Heatflux       (REAL_TYPE* qbc, int* bh1, int face, double& flop);
  REAL_TYPE ps_OBC_SpecVH         (REAL_TYPE* ws,  int* bh1, int face, REAL_TYPE* t, REAL_TYPE tm, REAL_TYPE* v00, double& flop);
  REAL_TYPE ps_OBC_HeatTransfer_BS(REAL_TYPE* qbc, int* bh1, int face, REAL_TYPE* t, REAL_TYPE* t0, double& flop);
  REAL_TYPE ps_OBC_HeatTransfer_SF(REAL_TYPE* qbc, int* bh1, int face, REAL_TYPE* t, REAL_TYPE* t0, double& flop);
  REAL_TYPE ps_OBC_HeatTransfer_SN(REAL_TYPE* qbc, int* bh1, int face, REAL_TYPE* t, REAL_TYPE* t0, double& flop);
  REAL_TYPE ps_OBC_IsoThermal     (REAL_TYPE* qbc, int* bh1, int face, REAL_TYPE* t, REAL_TYPE* t0, double& flop);
  
  void Pibc_Prdc                (REAL_TYPE* d_p, int* st, int* ed, int* d_bcd, int odr, int dir, REAL_TYPE pv);
  void Pobc_Prdc_Directional    (REAL_TYPE* d_p, int face, REAL_TYPE pv, int uod);
  void Pobc_Prdc_Simple         (REAL_TYPE* d_p, int face);
  void ps_IBC_ConstTemp         (REAL_TYPE* t, int* bh2, int n);
  
  
  /**
   * @brief 温度の外部周期境界条件（単純なコピー）
   * @param [in] t    温度のデータクラス
   * @param [in] face 面番号
   */
  void Tobc_Prdc_Simple         (REAL_TYPE* d_t, const int face);
  
  void Vibc_Prdc                (REAL_TYPE* d_v, int* st, int* ed, int* d_bd, int odr, int dir);
  void Vobc_Prdc                (REAL_TYPE* d_v, int face, int no_comm_face);
  void Vobc_Prdc_CF             (REAL_TYPE* d_v, int face);
  
public:
  void assign_Temp          (REAL_TYPE* t, int* bh, REAL_TYPE tm, Control* C);
  void assign_Velocity      (REAL_TYPE* v, int* bv, REAL_TYPE tm, REAL_TYPE* v00, bool clear=false);
  void checkDriver          (FILE* fp);
  void InnerPBC_Periodic    (REAL_TYPE* d_p, int* d_bcd);
  void InnerTBCface         (REAL_TYPE* qbc, int* bx, REAL_TYPE* t, REAL_TYPE* t0, double& flop);
  void InnerTBCvol          (REAL_TYPE* t, int* bx, REAL_TYPE dt, double& flop);
  void InnerVBC             (REAL_TYPE* v, int* bv, REAL_TYPE tm, REAL_TYPE* v00, double& flop);
  void InnerVBC_Periodic    (REAL_TYPE* d_v, int* d_bd);
  void mod_div              (REAL_TYPE* div, int* bv, REAL_TYPE coef, REAL_TYPE tm, REAL_TYPE* v00, REAL_TYPE* avr, 
                             double& flop);
  void mod_Dir_Forcing      (REAL_TYPE* v, int* bd, float* cvf, REAL_TYPE* v00, REAL_TYPE &flop);
  void mod_Psrc_VBC         (REAL_TYPE* div, REAL_TYPE* vc, REAL_TYPE* v0, REAL_TYPE coef, 
                             int* bv, REAL_TYPE tm, REAL_TYPE dt, Control* C, REAL_TYPE* v00, REAL_TYPE &flop);
  void mod_Psrc_Forcing     (REAL_TYPE* src, REAL_TYPE* v, int* bd, float* cvf, REAL_TYPE dh, REAL_TYPE* v00, 
                             REAL_TYPE** c_array, REAL_TYPE &flop);
  void mod_Pvec_Flux        (REAL_TYPE* wv, REAL_TYPE* v, int* bv, REAL_TYPE tm, Control* C, int v_mode, REAL_TYPE* v00, double& flop);
  void mod_Pvec_Forcing     (REAL_TYPE* vc, REAL_TYPE* v, int* bd, float* cvf, REAL_TYPE* v00, REAL_TYPE dt, REAL_TYPE &flop);
  void mod_Vdiv_Forcing     (REAL_TYPE* v, int* bd, float* cvf, REAL_TYPE* div, REAL_TYPE dt, REAL_TYPE dh, REAL_TYPE* v00, 
                             REAL_TYPE* am, REAL_TYPE** c_array, REAL_TYPE &flop);
  void mod_Vis_EE           (REAL_TYPE* vc, REAL_TYPE* v0, REAL_TYPE cf, int* bx, REAL_TYPE tm, REAL_TYPE dt, REAL_TYPE* v00, double& flop);
  void OuterPBC             (REAL_TYPE* d_p);
  void OuterTBC             (REAL_TYPE* d_t);
  void OuterTBCface         (REAL_TYPE* qbc, int* bx, REAL_TYPE* t, REAL_TYPE* t0, Control* C, double& flop);
  void OuterVBC             (REAL_TYPE* v, REAL_TYPE* vc, int* bv, REAL_TYPE tm, REAL_TYPE dt, Control* C, REAL_TYPE* v00, double& flop);
  void OuterVBC_Periodic    (REAL_TYPE* d_v);
  void OuterVBC_Pseudo      (REAL_TYPE* vc, REAL_TYPE* v0, int* bv, REAL_TYPE tm, REAL_TYPE dt, Control* C, REAL_TYPE* v00, double& flop);
  void ps_BC_Convection     (REAL_TYPE* ws, int* bh1, REAL_TYPE* v, REAL_TYPE* t, REAL_TYPE tm, Control* C, REAL_TYPE* v00, double& flop);
  void setBCIperiodic       (int* d_bx);
  void setInitialTemp_Compo (int n, int* bx, REAL_TYPE* d_t);
  void updateOuter          (REAL_TYPE* v, REAL_TYPE* vc);
  
  REAL_TYPE setDirectForcing (REAL_TYPE* v, int* bx, int n, REAL_TYPE v00);
};

#endif // _FFV_SETBC_H_
