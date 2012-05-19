#ifndef _CBC_SETBC_H_
#define _CBC_SETBC_H_

/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file CBC_SetBC.h
//@brief SetBC3D class Header
//@author keno, FSI Team, VCAD, RIKEN

#include <math.h>
#include "SetBC.h"
#include "FortranFuncCBC.h"
#include "CBC_Define.h"

class SetBC3D : public SetBC {
public:
  SetBC3D() {}
  virtual ~SetBC3D() {}
  
protected:
  REAL_TYPE extractVel_IBC        (int n, REAL_TYPE* vec, REAL_TYPE tm, REAL_TYPE* v00, REAL_TYPE& flop);
  REAL_TYPE extractVel_OBC        (int n, REAL_TYPE* vec, REAL_TYPE tm, REAL_TYPE* v00, REAL_TYPE& flop);
  REAL_TYPE ps_IBC_Heatflux       (REAL_TYPE* qbc, unsigned* bh1, int n, REAL_TYPE& flop);
  REAL_TYPE ps_IBC_HeatGen_SM     (REAL_TYPE* t,   unsigned* bh2, int n, REAL_TYPE dt, REAL_TYPE& flop);
  REAL_TYPE ps_IBC_IsoThermal_SM  (REAL_TYPE* qbc, unsigned* bh1, int n, REAL_TYPE* t, REAL_TYPE* t0, REAL_TYPE& flop);
  REAL_TYPE ps_IBC_Transfer_B_SM  (REAL_TYPE* qbc, unsigned* bh1, int n, REAL_TYPE* t, REAL_TYPE* t0, REAL_TYPE& flop);
  REAL_TYPE ps_IBC_Transfer_S_SM  (REAL_TYPE* qbc, unsigned* bh1, int n, REAL_TYPE* t, REAL_TYPE* t0, REAL_TYPE& flop);
  REAL_TYPE ps_IBC_Transfer_SF_SM (REAL_TYPE* qbc, unsigned* bh1, int n, REAL_TYPE* t, REAL_TYPE* t0, REAL_TYPE& flop);
  REAL_TYPE ps_IBC_Transfer_SN_SM (REAL_TYPE* qbc, unsigned* bh1, int n, REAL_TYPE* t, REAL_TYPE* t0, REAL_TYPE& flop);
  REAL_TYPE ps_IBC_Outflow        (REAL_TYPE* ws,  unsigned* bh1, int n, REAL_TYPE* v, REAL_TYPE* t, REAL_TYPE* v00, REAL_TYPE& flop);
  REAL_TYPE ps_IBC_SpecVH         (REAL_TYPE* ws,  unsigned* bh1, int n, REAL_TYPE v00, REAL_TYPE* vec, REAL_TYPE& flop);
  REAL_TYPE ps_OBC_Free           (REAL_TYPE* ws,  unsigned* bh1, int face, REAL_TYPE* v, REAL_TYPE* t, REAL_TYPE* v00, REAL_TYPE& flop);
  REAL_TYPE ps_OBC_Heatflux       (REAL_TYPE* qbc, unsigned* bh1, int face, REAL_TYPE& flop);
  REAL_TYPE ps_OBC_SpecVH         (REAL_TYPE* ws,  unsigned* bh1, int face, REAL_TYPE* t, REAL_TYPE tm, REAL_TYPE* v00, REAL_TYPE& flop);
  REAL_TYPE ps_OBC_HeatTransfer_BS(REAL_TYPE* qbc, unsigned* bh1, int face, REAL_TYPE* t, REAL_TYPE* t0, REAL_TYPE& flop);
  REAL_TYPE ps_OBC_HeatTransfer_SF(REAL_TYPE* qbc, unsigned* bh1, int face, REAL_TYPE* t, REAL_TYPE* t0, REAL_TYPE& flop);
  REAL_TYPE ps_OBC_HeatTransfer_SN(REAL_TYPE* qbc, unsigned* bh1, int face, REAL_TYPE* t, REAL_TYPE* t0, REAL_TYPE& flop);
  REAL_TYPE ps_OBC_IsoThermal     (REAL_TYPE* qbc, unsigned* bh1, int face, REAL_TYPE* t, REAL_TYPE* t0, REAL_TYPE& flop);
  
  void Pibc_Prdc                (SklScalar3D<REAL_TYPE>* d_p, int* st, int* ed, SklScalar3D<unsigned>* d_bcd, int odr, int dir, REAL_TYPE pv);
  void Pobc_Prdc_Directional    (SklScalar3D<REAL_TYPE>* d_p, int face, REAL_TYPE pv, unsigned uod);
  void Pobc_Prdc_Simple         (SklScalar3D<REAL_TYPE>* d_p, int face);
  void ps_IBC_ConstTemp         (REAL_TYPE* t, unsigned* bh2, int n);
  void Tobc_Prdc_Simple         (SklScalar3D<REAL_TYPE>* d_t, int face);
  void Vibc_Prdc                (SklVector3DEx<REAL_TYPE>* d_v, int* st, int* ed, SklScalar3D<unsigned>* d_bd, int odr, int dir);
  void Vobc_Prdc                (SklVector3DEx<REAL_TYPE>* d_v, int face, unsigned no_comm_face);
  void Vobc_Prdc_CF             (SklVector3DEx<REAL_TYPE>* d_v, int face);
  
public:
  void assign_Temp          (REAL_TYPE* t, unsigned* bh, REAL_TYPE tm, Control* C);
  void assign_Velocity      (REAL_TYPE* v, unsigned* bv, REAL_TYPE tm, REAL_TYPE* v00, bool clear=false);
  void checkDriver          (FILE* fp);
  void InnerPBC_Periodic    (SklScalar3D<REAL_TYPE>* d_p, SklScalar3D<unsigned>* d_bcd);
  void InnerTBCface         (REAL_TYPE* qbc, unsigned* bx, REAL_TYPE* t, REAL_TYPE* t0, REAL_TYPE& flop);
  void InnerTBCvol          (REAL_TYPE* t, unsigned* bx, REAL_TYPE dt, REAL_TYPE& flop);
  void InnerVBC             (REAL_TYPE* v, unsigned* bv, REAL_TYPE tm, REAL_TYPE* v00, REAL_TYPE& flop);
  void InnerVBC_Periodic    (SklVector3DEx<REAL_TYPE>* d_v, SklScalar3D<unsigned>* d_bd);
  void mod_div              (REAL_TYPE* div, unsigned* bv, REAL_TYPE coef, REAL_TYPE tm, REAL_TYPE* v00, REAL_TYPE* avr, 
                             REAL_TYPE& flop);
  void mod_Dir_Forcing      (REAL_TYPE* v, unsigned* bd, float* cvf, REAL_TYPE* v00, REAL_TYPE &flop);
  void mod_Psrc_VBC         (REAL_TYPE* div, REAL_TYPE* vc, REAL_TYPE* v0, REAL_TYPE coef, unsigned* bv, REAL_TYPE tm, REAL_TYPE dt, Control* C, 
                             REAL_TYPE* v00, REAL_TYPE &flop);
  void mod_Psrc_Forcing     (REAL_TYPE* src, REAL_TYPE* v, unsigned* bd, float* cvf, REAL_TYPE dh, REAL_TYPE* v00, 
                             REAL_TYPE** c_array, REAL_TYPE &flop);
  void mod_Pvec_Flux        (REAL_TYPE* wv, REAL_TYPE* v, unsigned* bv, REAL_TYPE tm, Control* C, int v_mode, REAL_TYPE* v00, REAL_TYPE& flop);
  void mod_Pvec_Forcing     (REAL_TYPE* vc, REAL_TYPE* v, unsigned* bd, float* cvf, REAL_TYPE* v00, REAL_TYPE dt, REAL_TYPE &flop);
  void mod_Vdiv_Forcing     (REAL_TYPE* v, unsigned* bd, float* cvf, REAL_TYPE* div, REAL_TYPE dt, REAL_TYPE dh, REAL_TYPE* v00, 
                             REAL_TYPE* am, REAL_TYPE** c_array, REAL_TYPE &flop);
  void mod_Vis_EE           (REAL_TYPE* vc, REAL_TYPE* v0, REAL_TYPE cf, unsigned* bx, REAL_TYPE tm, REAL_TYPE dt, REAL_TYPE* v00, REAL_TYPE& flop);
  void OuterPBC             (SklScalar3D<REAL_TYPE>* d_p);
  void OuterTBC             (SklScalar3D<REAL_TYPE>* d_t);
  void OuterTBCface         (REAL_TYPE* qbc, unsigned* bx, REAL_TYPE* t, REAL_TYPE* t0, Control* C, REAL_TYPE& flop);
  void OuterVBC             (REAL_TYPE* v, REAL_TYPE* vc, unsigned* bv, REAL_TYPE tm, REAL_TYPE dt, Control* C, REAL_TYPE* v00, REAL_TYPE& flop);
  void OuterVBC_Periodic    (SklVector3DEx<REAL_TYPE>* d_v);
  void OuterVBC_Pseudo      (REAL_TYPE* vc, REAL_TYPE* v0, unsigned* bv, REAL_TYPE tm, REAL_TYPE dt, Control* C, REAL_TYPE* v00, REAL_TYPE& flop);
  void ps_BC_Convection     (REAL_TYPE* ws, unsigned* bh1, REAL_TYPE* v, REAL_TYPE* t, REAL_TYPE tm, Control* C, REAL_TYPE* v00, REAL_TYPE& flop);
  void setBCIperiodic       (SklScalar3D<unsigned>* d_bx);
  void setInitialTemp_Compo (unsigned n, unsigned* bx, SklScalar3D<REAL_TYPE>* d_t);
  void updateOuter          (REAL_TYPE* v, REAL_TYPE* vc);
  
  REAL_TYPE setDirectForcing (REAL_TYPE* v, unsigned* bx, int n, REAL_TYPE v00);
};

#endif // _CBC_SETBC_H_
