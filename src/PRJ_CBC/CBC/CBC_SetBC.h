#ifndef _SKL_CBC_SETBC_H_
#define _SKL_CBC_SETBC_H_

/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2011
 *
 */

//@file CBC_SetBC.h
//@brief SetBC3D class Header
//@author keno, FSI Team, VCAD, RIKEN

#include <math.h>
#include "SetBC.h"
#include "FortranFuncCBC.h"
#include "SklSolverCBCDefine.h"

class SetBC3D : public SetBC {
public:
  SetBC3D() {}
  virtual ~SetBC3D() {}
  
protected:
  SKL_REAL extractVel_IBC        (int n, SKL_REAL* vec, SKL_REAL tm, SKL_REAL* v00, SKL_REAL& flop);
  SKL_REAL extractVel_OBC        (int n, SKL_REAL* vec, SKL_REAL tm, SKL_REAL* v00, SKL_REAL& flop);
  SKL_REAL ps_IBC_Heatflux       (SKL_REAL* qbc, unsigned* bh1, int n, SKL_REAL& flop);
  SKL_REAL ps_IBC_HeatGen_SM     (SKL_REAL* t,   unsigned* bh2, int n, SKL_REAL dt, SKL_REAL& flop);
  SKL_REAL ps_IBC_IsoThermal_SM  (SKL_REAL* qbc, unsigned* bh1, int n, SKL_REAL* t, SKL_REAL* t0, SKL_REAL& flop);
  SKL_REAL ps_IBC_Transfer_B_SM  (SKL_REAL* qbc, unsigned* bh1, int n, SKL_REAL* t, SKL_REAL* t0, SKL_REAL& flop);
  SKL_REAL ps_IBC_Transfer_S_SM  (SKL_REAL* qbc, unsigned* bh1, int n, SKL_REAL* t, SKL_REAL* t0, SKL_REAL& flop);
  SKL_REAL ps_IBC_Transfer_SF_SM (SKL_REAL* qbc, unsigned* bh1, int n, SKL_REAL* t, SKL_REAL* t0, SKL_REAL& flop);
  SKL_REAL ps_IBC_Transfer_SN_SM (SKL_REAL* qbc, unsigned* bh1, int n, SKL_REAL* t, SKL_REAL* t0, SKL_REAL& flop);
  SKL_REAL ps_IBC_Outflow        (SKL_REAL* ws,  unsigned* bh1, int n, SKL_REAL* v, SKL_REAL* t, SKL_REAL* v00, SKL_REAL& flop);
  SKL_REAL ps_IBC_SpecVH         (SKL_REAL* ws,  unsigned* bh1, int n, SKL_REAL v00, SKL_REAL* vec, SKL_REAL& flop);
  SKL_REAL ps_OBC_Free           (SKL_REAL* ws,  unsigned* bh1, int face, SKL_REAL* v, SKL_REAL* t, SKL_REAL* v00, SKL_REAL& flop);
  SKL_REAL ps_OBC_Heatflux       (SKL_REAL* qbc, unsigned* bh1, int face, SKL_REAL& flop);
  SKL_REAL ps_OBC_SpecVH         (SKL_REAL* ws,  unsigned* bh1, int face, SKL_REAL* t, SKL_REAL tm, SKL_REAL* v00, SKL_REAL& flop);
  SKL_REAL ps_OBC_HeatTransfer_BS(SKL_REAL* qbc, unsigned* bh1, int face, SKL_REAL* t, SKL_REAL* t0, SKL_REAL& flop);
  SKL_REAL ps_OBC_HeatTransfer_SF(SKL_REAL* qbc, unsigned* bh1, int face, SKL_REAL* t, SKL_REAL* t0, SKL_REAL& flop);
  SKL_REAL ps_OBC_HeatTransfer_SN(SKL_REAL* qbc, unsigned* bh1, int face, SKL_REAL* t, SKL_REAL* t0, SKL_REAL& flop);
  SKL_REAL ps_OBC_IsoThermal     (SKL_REAL* qbc, unsigned* bh1, int face, SKL_REAL* t, SKL_REAL* t0, SKL_REAL& flop);
  
  void flip_ObcMask             (int face, unsigned* bv, unsigned flag);
  void Pibc_Prdc                (SklScalar3D<SKL_REAL>* d_p, int* st, int* ed, SklScalar3D<unsigned>* d_bcd, int odr, int dir, SKL_REAL pv);
  void Pobc_Prdc_Directional    (SklScalar3D<SKL_REAL>* d_p, int face, SKL_REAL pv, unsigned uod);
  void Pobc_Prdc_Simple         (SklScalar3D<SKL_REAL>* d_p, int face);
  void ps_IBC_ConstTemp         (SKL_REAL* t, unsigned* bh2, int n);
  void Tobc_Prdc_Simple         (SklScalar3D<SKL_REAL>* d_t, int face);
  void Vibc_Prdc                (SklVector3D<SKL_REAL>* d_v, int* st, int* ed, SklScalar3D<unsigned>* d_bd, int odr, int dir);
  void Vobc_Prdc                (SklVector3D<SKL_REAL>* d_v, int face, unsigned no_comm_face);
  void Vobc_Prdc_CF             (SklVector3D<SKL_REAL>* d_v, int face);
  
public:
  void assign_Temp          (SKL_REAL* t, unsigned* bh, SKL_REAL tm, Control* C);
  void assign_Velocity      (SKL_REAL* v, unsigned* bv, SKL_REAL tm, SKL_REAL* v00, bool clear=false);
  void checkDriver          (FILE* fp);
  void flipDir_OBC          (unsigned* bv, Control* C);
  void InnerPBC_Periodic    (SklScalar3D<SKL_REAL>* d_p, SklScalar3D<unsigned>* d_bcd);
  void InnerTBCface         (SKL_REAL* qbc, unsigned* bx, SKL_REAL* t, SKL_REAL* t0, SKL_REAL& flop);
  void InnerTBCvol          (SKL_REAL* t, unsigned* bx, SKL_REAL dt, SKL_REAL& flop);
  void InnerVBC             (SKL_REAL* v, unsigned* bv, SKL_REAL tm, SKL_REAL* v00, SKL_REAL& flop, bool isCDS=false);
  void InnerVBC_Periodic    (SklVector3D<SKL_REAL>* d_v, SklScalar3D<unsigned>* d_bd);
  void mod_div              (SKL_REAL* div, unsigned* bv, SKL_REAL coef, SKL_REAL tm, SKL_REAL* v00, SKL_REAL& flop, bool isCDS=false);
  void mod_Psrc_Forcing     (SKL_REAL* src, SKL_REAL* v, unsigned* bd, Control* C, SKL_REAL* v00, SKL_REAL& flop);
  void mod_Psrc_VBC         (SKL_REAL* dv, SKL_REAL* vc, SKL_REAL* v0, SKL_REAL coef, unsigned* bv, SKL_REAL tm, SKL_REAL dt, Control* C, 
                             SKL_REAL* v00, SKL_REAL& flop, bool isCDS=false);
  void mod_Pvec_Flux        (SKL_REAL* wv, SKL_REAL* v, unsigned* bv, SKL_REAL tm, Control* C, int v_mode, SKL_REAL* v00, SKL_REAL& flop, bool isCDS=false);
  void mod_Pvec_Forcing     (SKL_REAL* vc, unsigned* bd, Control* C, SKL_REAL* v00, SKL_REAL& flop);
  void mod_Vcc_Forcing      (SKL_REAL* v,  unsigned* bd, Control* C, SKL_REAL dt, SKL_REAL* v00, SKL_REAL& flop);
  void mod_Vcf_Forcing      (SKL_REAL* v, unsigned* bd, Control* C, SKL_REAL dt, SKL_REAL* v00, SKL_REAL& flop);
  void mod_Vis_EE           (SKL_REAL* vc, SKL_REAL* v0, SKL_REAL cf, unsigned* bx, SKL_REAL tm, SKL_REAL dt, SKL_REAL* v00, SKL_REAL& flop);
  void OuterPBC             (SklScalar3D<SKL_REAL>* d_p);
  void OuterTBC             (SklScalar3D<SKL_REAL>* d_t);
  void OuterTBCface         (SKL_REAL* qbc, unsigned* bx, SKL_REAL* t, SKL_REAL* t0, Control* C, SKL_REAL& flop);
  void OuterVBC             (SKL_REAL* v, SKL_REAL* vc, unsigned* bv, SKL_REAL tm, SKL_REAL dt, Control* C, SKL_REAL* v00, SKL_REAL& flop);
  void OuterVBC_Periodic    (SklVector3D<SKL_REAL>* d_v);
  void OuterVBC_Pseudo      (SKL_REAL* vc, SKL_REAL* v0, unsigned* bv, SKL_REAL tm, SKL_REAL dt, Control* C, SKL_REAL* v00, SKL_REAL& flop);
  void ps_BC_Convection     (SKL_REAL* ws, unsigned* bh1, SKL_REAL* v, SKL_REAL* t, SKL_REAL tm, Control* C, SKL_REAL* v00, SKL_REAL& flop);
  void setBCIperiodic       (SklScalar3D<unsigned>* d_bx);
  void setInitialTemp_Compo (unsigned n, unsigned* bx, SklScalar3D<SKL_REAL>* d_t);
  void updateOuter          (SKL_REAL* v, SKL_REAL* vc);
  
  SKL_REAL setDirectForcing (SKL_REAL* v, unsigned* bx, int n, SKL_REAL v00);
};

#endif // _SKL_C3D_SETBC_H_
