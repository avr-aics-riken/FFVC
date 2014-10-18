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
 * @file   ffv_LS.C
 * @brief  LS Class
 * @author aics
 */

#include "ffv_LS.h"


// #################################################################
int LinearSolver::BiCGstab(REAL_TYPE* x, REAL_TYPE* b, const double b_l2, const double r0_l2)
{
  double flop_count=0.0;  /// 浮動小数点演算数
  int lc=0;               /// ループカウント
  double var[3];          /// 誤差, 残差, 解ベクトルのL2ノルム
  double rho    = 1.0;    /// \rho
  
  
  
  //TIMING_start(tm_blas_calc_r);
  flop_count=0.0;
  if ( getNaive()==OFF )
  {
    blas_calc_rk_(pcg_r, x, b, bcp, size, &guide, &flop_count); // r = b - Ax
  }
  else
  {
    blas_calc_r_naive_(pcg_r, x, b, bcp, size, &guide, pni, &flop_count); // r = b - Ax
  }
  
  //TIMING_stop(tm_blas_calc_r, flop_count);
  
  // rho = (r, r)
  rho = Fdot1(pcg_r);
  
  // 既に収束
  if( fabs(rho) < REAL_TYPE_EPSILON )
  {
    return 0;
  }
  
  
  SyncScalar(pcg_r, 1);

  
  
  // 初期値
  //TIMING_start(tm_blas_copy);
  blas_copy_(pcg_r0, pcg_r, size, &guide); // r0 = r
  //blas_copy_(pcg_p,  pcg_r, size, &guide); // p = r
  //TIMING_stop(tm_blas_copy, 0.0, 2);
  
  blas_clear_(pcg_p , size, &guide); // p=0
  blas_clear_(pcg_q , size, &guide); // q=0
  
  rho    = 1.0;
  double alpha  = 1.0;
  double omega  = 1.0;
  
  for (lc=1; lc<getMaxIteration(); lc++)
  {
    double rho_0 = rho;
    
    
    // rho = (r0, r)
    rho = Fdot2(pcg_r0, pcg_r);
    
    if( fabs(rho) < REAL_TYPE_EPSILON )
    {
      lc = 0;
      break;
    }
    
    
    double beta = (rho / rho_0) * (alpha / omega);
    
    
    // p =  r + beta * ( p - omega * q )
    //TIMING_start(tm_blas_bicg_update_p);
    flop_count=0.0;
    bicg_update_p_(pcg_p, pcg_r, &beta, &omega, pcg_q, size, &guide, &flop_count);
    //TIMING_stop(tm_blas_bicg_update_p, flop_count);
    
    
    SyncScalar(pcg_p, 1);
    
    
    // q = A * p
    //TIMING_start(tm_blas_ax);
    flop_count=0.0;
    if ( getNaive()==OFF )
    {
      blas_calc_ax_(pcg_q, pcg_p, bcp, size, &guide, &flop_count);
    }
    else
    {
      blas_calc_ax_naive_(pcg_q, pcg_p, bcp, size, &guide, pni, &flop_count);
    }
    //TIMING_stop(tm_blas_ax, flop_count);
    
    
    // alpha = rho / (r0, q)
    double tmp = Fdot2(pcg_r0, pcg_q);
    
    if( fabs(tmp) < REAL_TYPE_EPSILON )
    {
      Hostonly_ printf("BiCGstab : Fdot2(pcg_r0, pcg_q)\n");
      lc = -1;
      break;
    }
    
    alpha = rho / tmp;
    
    
    // s = r - alpha * q
    //TIMING_start(tm_blas_z_xpay);
    tmp = -alpha;
    flop_count=0.0;
    blas_z_xpay_(pcg_s, pcg_r, &tmp, pcg_q, size, &guide, &flop_count);
    //TIMING_stop(tm_blas_z_xpay, flop_count);
    
    
    SyncScalar(pcg_s, 1);
    
    
    // t = A * s
    //TIMING_start(tm_blas_ax);
    flop_count=0.0;
    if ( getNaive()==OFF )
    {
      blas_calc_ax_(pcg_t, pcg_s, bcp, size, &guide, &flop_count);
    }
    else
    {
      blas_calc_ax_naive_(pcg_t, pcg_s, bcp, size, &guide, pni, &flop_count);
    }
    //TIMING_stop(tm_blas_ax, flop_count);
    
    
    // omega = (t, s) / (t, t)
    double t_s = Fdot2(pcg_t, pcg_s);
    
    double t_t = Fdot1(pcg_t);
    
    if ( fabs(t_t) < REAL_TYPE_EPSILON )
    {
      Hostonly_ printf("BiCGstab : Fdot1(pcg_t)\n");
      lc = -1;
      break;
    }
    
    if ( t_s == 0.0 )
    {
      Hostonly_ printf("BiCGstab : omega = 0.0\n");
      lc = -1;
      break;
    }
    
    omega = t_s / t_t;
    
    
    // x = x + alpha * p + omega * s
    double x_l2 = 0.0;
    double delta_x = 0.0;
    
    //TIMING_start(tm_blas_bicg_update_x);
    flop_count=0.0;
    bicg_update_x_(x, &alpha, pcg_p, &omega, pcg_s, bcp, &x_l2, &delta_x, size, &guide, &flop_count);
    //TIMING_stop(tm_blas_bicg_update_x, flop_count);
    
    var[0] = delta_x;
    var[2] = x_l2;
    
    
    // r = s - omega * t
    //TIMING_start(tm_blas_z_xpay);
    tmp = -omega;
    flop_count=0.0;
    blas_z_xpay_(pcg_r, pcg_s, &tmp, pcg_t, size, &guide, &flop_count);
    //TIMING_stop(tm_blas_z_xpay, flop_count);
    
    var[1] = Fdot1(pcg_r);
    
    
    // 収束判定
    if ( Fcheck(var, b_l2, r0_l2) == true ) break;

  }
  
  
  // 境界条件
  //TIMING_start(tm_poi_BC);
  BC->OuterPBC(x, ensPeriodic);
  
  if ( C->EnsCompo.periodic == ON ) BC->InnerPBCperiodic(x, bcd);
  //TIMING_stop(tm_poi_BC, 0.0);
  
  SyncScalar(x, 1);
  
  return lc;
}


// #################################################################
// 収束判定　非Div反復
bool LinearSolver::Fcheck(double* var, const double b_l2, const double r0_l2)
{
  
  // 自乗量の集約
  if ( numProc > 1 )
  {
    //TIMING_start(tm_poi_res_comm);
    double tmp[3];
    tmp[0] = var[0];
    tmp[1] = var[1];
    tmp[2] = var[2];
    if ( paraMngr->Allreduce(tmp, var, 3, MPI_SUM) != CPM_SUCCESS ) Exit(0);
    //TIMING_stop(tm_poi_res_comm, 6.0*numProc*sizeof(double) ); // 双方向 x ノード数　x 3
  }
  
  // L2 ノルム
  var[0] = sqrt(var[0]); // L2 of (x^{(m+1)} - x^{(m)})
  var[1] = sqrt(var[1]); // L2 of residual
  var[2] = sqrt(var[2]); // L2 of solution vector x^{(m)}
  
  double err = var[0];
  double res = var[1];
  double x_l2= var[2];  /// 解ベクトルのL2ノルム
  
  double ErrEPS = getErrCriterion();
  
  // 残差の保存
  switch ( getResType() )
  {
    case nrm_r_x:
      setResidual( (x_l2<ErrEPS) ? res/ErrEPS : res/x_l2 );
      break;
      
    case nrm_r_b:
      setResidual( (b_l2<ErrEPS) ? res/ErrEPS : res/b_l2 );
      break;
      
    case nrm_r_r0:
      setResidual( (r0_l2<ErrEPS) ? res/ErrEPS : res/r0_l2 );
      break;
      
    default:
      printf("\tInvalid Linear Solver for Pressure\n");
      Exit(0);
      break;
  }
  
  // 誤差の保存
  switch ( getErrType() )
  {
    case nrm_dx:
      setError( err );
      break;
      
    case nrm_dx_x:
      setError( (x_l2<ErrEPS) ? err/ErrEPS : err/x_l2 );
      break;
      
    default:
      printf("\tInvalid Error Norm for Pressure\n");
      Exit(0);
      break;
  }
  
  // 収束判定　性能測定モードでないときのみ収束判定を行う　誤差または残差が収束したら抜ける
  //if ( (C->Hide.PM_Test == OFF) && ( isResConverged() || isErrConverged()) )
  if ( (C->Hide.PM_Test == OFF) && isResConverged() )
  {
    return true;
  }
  
  return false;
}


// #################################################################
double LinearSolver::Fdot1(REAL_TYPE* x)
{
  double flop_count=0.0;          /// 浮動小数点演算数
  double xy = 0.0;
  
  //TIMING_start(tm_blas_dot1);
  flop_count=0.0;
  blas_dot1_(&xy, x, bcp, size, &guide, &flop_count);
  //TIMING_stop(tm_blas_dot1, flop_count);
  
  if ( numProc > 1 )
  {
    //TIMING_start(tm_blas_comm);
    double xy_tmp = xy;
    if  ( paraMngr->Allreduce(&xy_tmp, &xy, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
    //TIMING_stop(tm_blas_comm, 2.0*numProc*sizeof(double) );
  }
  
  return xy;
}


// #################################################################
double LinearSolver::Fdot2(REAL_TYPE* x, REAL_TYPE* y)
{
  double flop_count=0.0;          /// 浮動小数点演算数
  double xy = 0.0;
  
  //TIMING_start(tm_blas_dot2);
  flop_count=0.0;
  blas_dot2_(&xy, x, y, bcp, size, &guide, &flop_count);
  //TIMING_stop(tm_blas_dot2, flop_count);
  
  if ( numProc > 1 )
  {
    //TIMING_start(tm_blas_comm);
    double xy_tmp = xy;
    if  ( paraMngr->Allreduce(&xy_tmp, &xy, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
    //TIMING_stop(tm_blas_comm, 2.0*numProc*sizeof(double) );
  }
  
  return xy;
}


// #################################################################
void LinearSolver::Initialize(Control* C,
                              SetBC3D* BC,
                              int* bcp,
                              int* bcd,
                              REAL_TYPE* pni,
                              REAL_TYPE* pcg_p,
                              REAL_TYPE* pcg_p_,
                              REAL_TYPE* pcg_r,
                              REAL_TYPE* pcg_r0,
                              REAL_TYPE* pcg_q,
                              REAL_TYPE* pcg_s,
                              REAL_TYPE* pcg_s_,
                              REAL_TYPE* pcg_t,
                              REAL_TYPE* pcg_t_,
                              const int* ensP,
                              const int* cf_sz,
                              REAL_TYPE* cf_x,
                              REAL_TYPE* cf_y,
                              REAL_TYPE* cf_z)
{
  this->C = C;
  this->BC = BC;
  this->bcp = bcp;
  this->bcd = bcd;
  this->pni = pni;
  this->pcg_p  = pcg_p;
  this->pcg_p_ = pcg_p_;
  this->pcg_r  = pcg_r;
  this->pcg_r0 = pcg_r0;
  this->pcg_q  = pcg_q;
  this->pcg_s  = pcg_s;
  this->pcg_s_ = pcg_s_;
  this->pcg_t  = pcg_t;
  this->pcg_t_ = pcg_t_;
  this->cf_x = cf_x;
  this->cf_y = cf_y;
  this->cf_z = cf_z;
  
  for (int i=0; i<3; i++)
  {
    this->ensPeriodic[i] = ensP[i];
    this->cf_sz[i] = cf_sz[i];
  }
}

// #################################################################
// PBiCBSTAB
int LinearSolver::PBiCGstab(REAL_TYPE* x, REAL_TYPE* b, const double b_l2, const double r0_l2)
{
  double var[3];          /// 誤差, 残差, 解ベクトルのL2ノルム
  REAL_TYPE res = 0.0;
  REAL_TYPE omg = getOmega();
  
  blas_clear_(pcg_r , size, &guide);
  blas_clear_(pcg_p , size, &guide);
  
  blas_clear_(pcg_r0, size, &guide);
  blas_clear_(pcg_p_, size, &guide);
  blas_clear_(pcg_q, size, &guide);
  blas_clear_(pcg_s , size, &guide);
  blas_clear_(pcg_s_, size, &guide);
  blas_clear_(pcg_t_, size, &guide);
  
  blas_calcr_(pcg_r, x, b, bcp, size, &guide);             // (1)
  SyncScalar(pcg_r, 1);
  
  blas_copy_(pcg_r0, pcg_r, size, &guide);                 // (2)
  
  REAL_TYPE rr0 = 1.0;
  REAL_TYPE alpha = 0.0;
  REAL_TYPE gamma  = 1.0;
  REAL_TYPE gamman = -gamma;
  int lc=0;                      /// ループカウント
  
  for (lc=1; lc<getMaxIteration(); lc++)
  {
    REAL_TYPE rr1 = 0.0;
    Fdot(&rr1, pcg_r, pcg_r0);                             // (4)
    
    if( fabs(rr1) < FLT_MIN )                                  // (5)
    {
      res = rr1;
      lc = 0;
      break;
    }
    
    if( lc == 1 )
    {
      blas_copy_(pcg_p, pcg_r, size, &guide);              // (7)
    }
    else
    {
      REAL_TYPE beta = rr1/rr0*alpha/gamma;                    // (9)
      blas_axpy_(pcg_p, pcg_q, &gamman, size, &guide);    // (10)
      blas_xpay_(pcg_p, pcg_r , &beta  , size, &guide);    // (10)
    }
    SyncScalar(pcg_p, 1);
    
    blas_clear_(pcg_p_, size, &guide);
    Fpreconditioner(pcg_p_, pcg_p);                    // (12)
    
    blas_calcax_(pcg_q, pcg_p_, bcp, size, &guide);     // (13)
    
    REAL_TYPE q_r0 = 0.0;
    Fdot(&q_r0, pcg_q, pcg_r0);                           // (14)
    
    REAL_TYPE alpha  = rr1/q_r0;
    REAL_TYPE alphan = -alpha;
    blas_axpyz_(pcg_s, pcg_q, pcg_r, &alphan, size, &guide); // (15)
    SyncScalar(pcg_s, 1);
    
    blas_clear_(pcg_s_, size, &guide);
    Fpreconditioner(pcg_s_, pcg_s);                    // (17)
    
    blas_calcax_(pcg_t_, pcg_s_, bcp, size, &guide);     // (18)
    
    REAL_TYPE t_s = 0.0;
    Fdot(&t_s, pcg_t_, pcg_s);                             // (19a)
    
    REAL_TYPE t_t_ = 0.0;
    Fdot(&t_t_, pcg_t_, pcg_t_);                           // (19b)
    
    gamma  = t_s/t_t_;                                         // (19)
    gamman = -gamma;
    
    blas_axpbypz_(x      , pcg_p_, pcg_s_, &alpha , &gamma, size, &guide);  // (20)
    blas_axpyz_  (pcg_r, pcg_t_, pcg_s , &gamman, size, &guide);          // (21)
    
    REAL_TYPE rr = 0.0;
    Fdot(&rr, pcg_r, pcg_r);
    
    res = sqrt(rr);
    var[1] = res;
    

    if ( Fcheck(var, b_l2, r0_l2) == true ) break;
    
    rr0 = rr1;
  }
  
  BC->OuterPBC(x, ensPeriodic);
  if ( C->EnsCompo.periodic == ON )
  {
    BC->InnerPBCperiodic(x, bcd);
  }
  
  SyncScalar(x, 1);
  
  return lc;
}



/*
// #################################################################
int LinearSolver::PBiCGstab(REAL_TYPE* x, REAL_TYPE* b, const double b_l2, const double r0_l2)
{
  double flop_count=0.0;  /// 浮動小数点演算数
  int lc=0;               /// ループカウント
  double var[3];          /// 誤差, 残差, 解ベクトルのL2ノルム
  double rho    = 1.0;    /// \rho
  
  
  
  //TIMING_start(tm_blas_calc_r);
  flop_count=0.0;
  if ( getNaive()==OFF )
  {
    blas_calc_r_(pcg_r, x, b, bcp, size, &guide, &flop_count); // r = b - Ax
  }
  else
  {
    blas_calc_r_naive_(pcg_r, x, b, bcp, size, &guide, pni, &flop_count); // r = b - Ax
  }
  
  //TIMING_stop(tm_blas_calc_r, flop_count);
  
  // rho = (r, r)
  rho = Fdot1(pcg_r);
  
  // 既に収束
  if( fabs(rho) < REAL_TYPE_EPSILON )
  {
    return 0;
  }
  
  
  SyncScalar(pcg_r, 1);
  
  
  
  // 初期値
  //TIMING_start(tm_blas_copy);
  blas_copy_(pcg_r0, pcg_r, size, &guide); // r0 = r
  //blas_copy_(pcg_p,  pcg_r, size, &guide); // p = r
  //TIMING_stop(tm_blas_copy, 0.0, 2);
  
  blas_clear_(pcg_p , size, &guide); // p=0
  blas_clear_(pcg_q , size, &guide); // q=0
  
  
  rho    = 1.0;
  double alpha  = 1.0;
  double omega  = 1.0;
  
  for (lc=0; lc<getMaxIteration(); lc++)
  {
    double rho_0 = rho;
    
    
    // rho = (r0, r)
    rho = Fdot2(pcg_r0, pcg_r);
    
    if( fabs(rho) < REAL_TYPE_EPSILON )
    {
      lc = 0;
      break;
    }
    
    
    double beta = (rho / rho_0) * (alpha / omega);
    
    
    // p =  r + beta * ( p - omega * q )
    //TIMING_start(tm_blas_bicg_update_p);
    flop_count=0.0;
    bicg_update_p_(pcg_p, pcg_r, &beta, &omega, pcg_q, size, &guide, &flop_count);
    //TIMING_stop(tm_blas_bicg_update_p, flop_count);
    
    
    SyncScalar(pcg_p, 1);
    
    
    // p^* = M * p
    blas_clear_(pcg_p_ , size, &guide);
    Preconditioner(pcg_p_, pcg_p);
    
    
    // q = A * p^*
    //TIMING_start(tm_blas_ax);
    flop_count=0.0;
    if ( getNaive()==OFF )
    {
      blas_calc_ax_(pcg_q, pcg_p, bcp, size, &guide, &flop_count);
    }
    else
    {
      blas_calc_ax_naive_(pcg_q, pcg_p, bcp, size, &guide, pni, &flop_count);
    }
    //TIMING_stop(tm_blas_ax, flop_count);
    
    
    // alpha = rho / (r0, q)
    double tmp = Fdot2(pcg_r0, pcg_q);
    
    if( fabs(tmp) < REAL_TYPE_EPSILON )
    {
      Hostonly_ printf("BiCGstab : Fdot2(pcg_r0, pcg_q)\n");
      lc = -1;
      break;
    }
    
    alpha = rho / tmp;
    
    
    // s = r - alpha * q
    //TIMING_start(tm_blas_z_xpay);
    tmp = -alpha;
    flop_count=0.0;
    blas_z_xpay_(pcg_s, pcg_r, &tmp, pcg_q, size, &guide, &flop_count);
    //TIMING_stop(tm_blas_z_xpay, flop_count);
    
    
    SyncScalar(pcg_s, 1);
    
    
    // s^* = M * s
    blas_clear_(pcg_s_ , size, &guide);
    Preconditioner(pcg_s_, pcg_s);
    
    
    // t = A * s^*
    //TIMING_start(tm_blas_ax);
    flop_count=0.0;
    if ( getNaive()==OFF )
    {
      blas_calc_ax_(pcg_t, pcg_s_, bcp, size, &guide, &flop_count);
    }
    else
    {
      blas_calc_ax_naive_(pcg_t, pcg_s_, bcp, size, &guide, pni, &flop_count);
    }
    //TIMING_stop(tm_blas_ax, flop_count);
    
    
    // omega = (t, s) / (t, t)
    double t_s = Fdot2(pcg_t, pcg_s);
    
    double t_t = Fdot1(pcg_t);
    
    if ( fabs(t_t) < REAL_TYPE_EPSILON )
    {
      Hostonly_ printf("BiCGstab : Fdot1(pcg_t)\n");
      lc = -1;
      break;
    }
    
    if ( t_s == 0.0 )
    {
      Hostonly_ printf("BiCGstab : omega = 0.0\n");
      lc = -1;
      break;
    }
    
    omega = t_s / t_t;
    
    
    // x = x + alpha * p^* + omega * s^*
    double x_l2 = 0.0;
    double delta_x = 0.0;
    
    //TIMING_start(tm_blas_bicg_update_x);
    flop_count=0.0;
    bicg_update_x_(x, &alpha, pcg_p_, &omega, pcg_s_, bcp, &x_l2, &delta_x, size, &guide, &flop_count);
    //TIMING_stop(tm_blas_bicg_update_x, flop_count);
    
    var[0] = delta_x;
    var[2] = x_l2;
    
    
    // r = s - omega * t
    //TIMING_start(tm_blas_z_xpay);
    tmp = -omega;
    flop_count=0.0;
    blas_z_xpay_(pcg_r, pcg_s, &tmp, pcg_t, size, &guide, &flop_count);
    //TIMING_stop(tm_blas_z_xpay, flop_count);
    
    var[1] = Fdot1(pcg_r);
    
    
    // 収束判定
    if ( ConvergenceTest(var, b_l2, r0_l2) == true ) break;
    
  }
  
  
  // 境界条件
  //TIMING_start(tm_poi_BC);
  BC->OuterPBC(x, ensPeriodic);
  
  if ( C->EnsCompo.periodic == ON ) BC->InnerPBCperiodic(x, bcd);
  //TIMING_stop(tm_poi_BC, 0.0);
  
  SyncScalar(x, 1);
  
  return lc;
}
*/



// #################################################################
int LinearSolver::PointSOR(REAL_TYPE* x, REAL_TYPE* b, const double b_l2, const double r0_l2)
{
  double flop_count=0.0;          /// 浮動小数点演算数
  REAL_TYPE omg = getOmega();     /// 加速係数
  double var[3];                  /// 誤差、残差、解
  int lc=0;                       /// ループカウント
  
  // x     圧力 p^{n+1}
  // b     RHS vector
  // bcp   ビットフラグ
  
  for (lc=1; lc<getMaxIteration(); lc++)
  {
    var[0] = 0.0; // 誤差
    var[1] = 0.0; // 残差
    var[2] = 0.0; // 解
    
    // 反復処理
    //TIMING_start(tm_poi_PSOR);
    flop_count = 0.0;
    psor_(x, size, &guide, &omg, var, b, bcp, &flop_count);
    //TIMING_stop(tm_poi_PSOR, flop_count);
    
    
    // 境界条件
    //TIMING_start(tm_poi_BC);
    BC->OuterPBC(x, ensPeriodic);
    if ( C->EnsCompo.periodic == ON ) BC->InnerPBCperiodic(x, bcd);
    //TIMING_stop(tm_poi_BC, 0.0);
    
    
    // 同期処理
    SyncScalar(x, 1);
    
    
    // 収束判定
    if ( Fcheck(var, b_l2, r0_l2) == true ) break;
    
  }
  
  return lc;
}


// #################################################################
void LinearSolver::Preconditioner(REAL_TYPE* x, REAL_TYPE* b)
{
  int lc_max = 3;
  double dummy = 1.0;
  
  //for (int lc=0; lc<lc_max; lc++)
  //{
  //  Smoother(x, b);
  //}
  SOR2_SMA(x, b, lc_max, dummy, dummy);
}


// #################################################################
// Smoother
void LinearSolver::Smoother(REAL_TYPE* x, REAL_TYPE* b)
{
  int ip = 0;
  REAL_TYPE omg = getOmega();     /// 加速係数
  
  if ( numProc > 1 )
  {
    ip = (head[0]+head[1]+head[2]+1) % 2;
  }
  else
  {
    ip = 0;
  }
  
  for (int color=0; color<2; color++)
  {
    blas_smoother_core_(x, b, bcp, &ip, &color, &omg, size, &guide);
    
    BC->OuterPBC(x, ensPeriodic);
    
    if ( C->EnsCompo.periodic == ON )
    {
      BC->InnerPBCperiodic(x, bcd);
    }
    
    if ( numProc > 1 )
    {
      if ( getSyncMode() == comm_sync )
      {
        if ( paraMngr->BndCommS3D(x, size[0], size[1], size[2], guide, 1) != CPM_SUCCESS ) {
          Exit(0);
        }
      }
      else
      {
        int ireq[12];
        sma_comm_     (x, size, &guide, &color, &ip, cf_sz, cf_x, cf_y, cf_z, ireq, nID);
        sma_comm_wait_(x, size, &guide, &color, &ip, cf_sz, cf_x, cf_y, cf_z, ireq);
      }
    }
  }
}


// #################################################################
int LinearSolver::SOR2_SMA(REAL_TYPE* x, REAL_TYPE* b, const int itrMax, const double b_l2, const double r0_l2)
{
  int ip;                         /// ローカルノードの基点(1,1,1)のカラーを示すインデクス
  /// ip=0 > R, ip=1 > B
  double flop_count=0.0;          /// 浮動小数点演算数
  REAL_TYPE omg = getOmega();     /// 加速係数
  double var[3];                  /// 誤差、残差、解
  int lc=0;                       /// ループカウント
  
  // x     圧力 p^{n+1}
  // b     RHS vector
  // d_bcp ビットフラグ
  
  
  for (lc=1; lc<itrMax; lc++)
  {
    // 2色のマルチカラー(Red&Black)のセットアップ
    //TIMING_start(tm_poi_setup);
    
    // ip = 0 基点(1,1,1)がRからスタート
    //    = 1 基点(1,1,1)がBからスタート
    if ( numProc > 1 )
    {
      ip = (head[0]+head[1]+head[2]+1) % 2;
    }
    else
    {
      ip = 0;
    }
    //TIMING_stop(tm_poi_setup, 0.0);
    
    
    var[0] = 0.0; // 誤差
    var[1] = 0.0; // 残差
    var[2] = 0.0; // 解
    
    // 各カラー毎の間に同期, 残差は色間で積算する
    // R - color=0 / B - color=1
    for (int color=0; color<2; color++) {
      
      //TIMING_start(tm_poi_SOR2SMA);
      flop_count = 0.0; // 色間で積算しない
      if ( getNaive()==OFF)
      {
        psor2sma_core_(x, size, &guide, &ip, &color, &omg, var, b, bcp, &flop_count);
      }
      else
      {
        psor2sma_naive_(x, size, &guide, &ip, &color, &omg, var, b, bcp, pni, &flop_count);
      }
      
      //TIMING_stop(tm_poi_SOR2SMA, flop_count);
      
      // 境界条件
      //TIMING_start(tm_poi_BC);
      BC->OuterPBC(x, ensPeriodic);
      if ( C->EnsCompo.periodic == ON ) BC->InnerPBCperiodic(x, bcd);
      //TIMING_stop(tm_poi_BC, 0.0);
      
      
      // 同期処理
      if ( numProc > 1 )
      {
        //TIMING_start(tm_poi_comm);
        
        if ( getSyncMode() == comm_sync )
        {
          if ( paraMngr->BndCommS3D(x, size[0], size[1], size[2], guide, 1) != CPM_SUCCESS ) Exit(0); // 1 layer communication
        }
        else
        {
          int ireq[12];
          sma_comm_     (x, size, &guide, &color, &ip, cf_sz, cf_x, cf_y, cf_z, ireq, nID);
          sma_comm_wait_(x, size, &guide, &color, &ip, cf_sz, cf_x, cf_y, cf_z, ireq);
        }
        //TIMING_stop(tm_poi_comm, face_comm_size*0.5);
      }
    }
    
    
    // 収束判定 varは自乗量
    if ( Fcheck(var, b_l2, r0_l2) == true ) break;
    
  }
  
  return lc;
}



// #################################################################
// 反復変数の同期処理
void LinearSolver::SyncScalar(REAL_TYPE* d_class, const int num_layer)
{
  if ( numProc > 1 )
  {
    //TIMING_start(tm_poi_comm);
    
    if ( getSyncMode() == comm_sync )
    {
      if ( paraMngr->BndCommS3D(d_class, size[0], size[1], size[2], guide, num_layer) != CPM_SUCCESS ) Exit(0);
    }
    else
    {
      MPI_Request req[12];
      for (int i=0; i<12; i++) req[i] = MPI_REQUEST_NULL;
      
      if ( paraMngr->BndCommS3D_nowait(d_class, size[0], size[1], size[2], guide, num_layer, req ) != CPM_SUCCESS ) Exit(0);
      if ( paraMngr->wait_BndCommS3D  (d_class, size[0], size[1], size[2], guide, num_layer, req ) != CPM_SUCCESS ) Exit(0);
    }
    //TIMING_stop(tm_poi_comm, face_comm_size*(double)num_layer);
  }
}


/*
// #################################################################
// Flexible Gmres(m)
void FFV::Fgmres(IterationCtl* IC, const double rhs_nrm, const double r0)
{
  const double eps_1 = 1.0e-30;
  const double eps_2 = IC->getResCriterion();
  int mode_precond = 1; // pre-conditioning
  
  // 残差収束チェック
  if ( rhs_nrm < eps_1 ) return;
  
  // d_wg : work
  // d_ws : h^2 \Psi
  // d_p  : x
  // 
  
  // >>> Gmres section
  TIMING_start(tm_gmres_sor_sct);
  
  double beta, beta1, res, eps_abs, res_abs, al;
  double r4;
  double flop=0.0;
  
  const int Iteration_Max = IC->getMaxIteration();
  int m = FREQ_OF_RESTART;
  int s_length = (size[0]+2*guide) * (size[1]+2*guide) * (size[2]+2*guide);
  
  int isfin  = 0;
  
  const int oki  = 2;
  const int step = 3;
  int nrm  = 0;
  
  
  int Nmax = m+1; // 配列確保用，ゼロはダミー
      
  double *hg = new double[(Nmax+1) * (Nmax+1)]; // Hessnberg matrix
  double *cs = new double[Nmax+1];              // Cosine for Givens rotations
  double *sn = new double[Nmax+1];              // Sine   for Givens rotations
  double *rm = new double[Nmax+1];              // residual vector for minimization problem
  
  memset(hg, 0, (Nmax+1)*(Nmax+1)*sizeof(double) );
  memset(cs, 0, (Nmax+1)*sizeof(double) );
  memset(sn, 0, (Nmax+1)*sizeof(double) );
  memset(rm, 0, (Nmax+1)*sizeof(double) );
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  int first = 1;

  
  // Outer iteration
  for (int j=1; j<=Iteration_Max; j++) {
    
    TIMING_start(tm_gmres_mvprod);
    flop = 0.0;
    matvec_p_(d_wg, size, &guide, d_p, d_bcp, &flop); // d_wg <- Ax
    TIMING_stop(tm_gmres_mvprod, flop);
    
    
    TIMING_start(tm_gmres_res_sample);
    flop = 0.0;
    residual_(d_wg, size, &guide, &res_abs, d_ws, d_sq, d_bcp, &flop); // d_wg <- b - Ax (= r); res_abs <- \sum{r^2}
    TIMING_stop(tm_gmres_res_sample, flop);
    
    
    TIMING_start(tm_gmres_comm);
    if ( numProc > 1 )
    {
      double tmp = res_abs;
      if ( paraMngr->Allreduce(&tmp, &res_abs, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
    }
    TIMING_stop(tm_gmres_comm, 2.0*numProc*sizeof(double)); // 双方向 x ノード数
    
    
    beta = sqrt(res_abs); // beta <- \|r\|_2
    
    // convergence check
    if (beta < eps_1) goto jump_2;
    
    
    
    // set initial value of residual vector
    rm[1] = beta; // rm <- {\beta, 0, 0, ..., 0}^T

    
    beta1 = 1.0 / beta;
    
    orth_basis_(d_vm, size, &guide, &Nmax, &first, &beta1, d_wg, &flop); // d_vm(1) <- v^1 
    
    
    
    // Orthogoanl basis
    for (int i=1; i<=m; i++) {
      
      // address
      int g = guide;
      size_t adrs = _F_IDX_S4D(1-g, 1-g, 1-g, i-1, size[0], size[1], size[2], guide);
      
      if ( mode_precond == 1)
      {
        // Decide the number of iteration for inner-iteration
        int n_inner = 3; //((im / oki) + 1) * step;
        
        // Inner-iteration
        for (int i_inner=1; i_inner<=n_inner; i_inner++) {
          SOR_2_SMA(IC, &d_zm[adrs], &d_vm[adrs], rhs_nrm, r0); // z^i <- K^{-1} v^i
        }
      }
      else
      {
        cp_orth_basis_(&d_zm[adrs], size, &guide, &d_vm[adrs]); // z^i <- v^i
      }
      
      TIMING_start(tm_gmres_mvprod);
      flop = 0.0;
      matvec_p_(d_wg, size, &guide, &d_zm[adrs], d_bcp, &flop); // d_wg <- Az^i
      TIMING_stop(tm_gmres_mvprod, flop);
      
      
      for (int km=1; km<=i; km++) {
        al = 0.0;
        ml_add_1_(&al, size, &guide, &Nmax, d_vm, d_wg, &km, &flop);
        
        
        TIMING_start(tm_gmres_comm);
        if ( numProc > 1 )
        {
          double tmp = al;
          if ( paraMngr->Allreduce(&tmp, &al, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
        }
        TIMING_stop(tm_gmres_comm, 2.0*numProc*sizeof(double)); // 双方向 x ノード数
        
        
        hg[_IDX2D(km, i, Nmax)] = al;
        
        double r_al = -al;
        ml_add_3_(d_wg, size, &guide, &Nmax, &r_al, d_vm, &km, &flop);
      }
      
      
      al =0.0;
      ml_add_2_(&al, size, &guide, d_wg, &flop);
      
      
      TIMING_start(tm_gmres_comm);
      if ( numProc > 1 )
      {
        double tmp = al;
        if ( paraMngr->Allreduce(&tmp, &al, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
      }
      TIMING_stop(tm_gmres_comm, 2.0*numProc*sizeof(double)); // 双方向 x ノード数
      
      
      TIMING_start(tm_gmres_others);
      flop = 0.0;
      
      hg[_IDX2D(i+1, i, Nmax)] = sqrt(al);
      flop += 20.0;
      
      if ( hg[_IDX2D(i+1, i, Nmax)] < eps_1 )
      {
        nrm = i-1;
        TIMING_stop(tm_gmres_others, flop);
        goto jump_1;
      }
      
      if (i < Nmax)
      {
        int idx = i + 1;
        al = 1.0 / hg[_IDX2D(idx, i, Nmax)];
        flop += 13.0;
        orth_basis_(d_vm, size, &guide, &Nmax, &idx, &al, d_wg, &flop);
      }
      
      rm[_IDX2D(1, i, Nmax)] = hg[_IDX2D(1, i, Nmax)];
      
      for (int km=2; km<=i; km++) {
        rm[_IDX2D(km  , i, Nmax)] = cs[km-1]*hg[_IDX2D(km, i, Nmax)] - sn[km-1]*rm[_IDX2D(km-1, i, Nmax)];
        rm[_IDX2D(km-1, i, Nmax)] = sn[km-1]*hg[_IDX2D(km, i, Nmax)] + cs[km-1]*rm[_IDX2D(km-1, i, Nmax)];
      }
      flop += (double)(i-2+1)*6.0;
      
      al = sqrt(rm[_IDX2D(i, i, Nmax)]*rm[_IDX2D(i, i, Nmax)] + hg[_IDX2D(i+1, i, Nmax)]*hg[_IDX2D(i+1, i, Nmax)]);
      flop += 23;
      
      if ( al < eps_1 )
      {
        nrm = i - 1;
        TIMING_stop(tm_gmres_others, flop);
        goto jump_1;
      }
      
      cs[i]    = rm[_IDX2D(i,  i, Nmax)] / al;
      sn[i]    = hg[_IDX2D(i+1,i, Nmax)] / al;
      
      rm[_IDX2D(i, i, Nmax)] = cs[i] * rm[_IDX2D(i, i, Nmax)] + sn[i] * hg[_IDX2D(i+1, i, Nmax)];
      
      rm[i+1]  = - sn[i]*rm[i];
      rm[i]    =   cs[i]*rm[i];
      
      al = rm[i+1] * rm[i+1];
      
      flop += 2.0*13.0 + 4.0;
      
      if (al < eps_abs)
      {
        nrm = i;
        TIMING_stop(tm_gmres_others, flop);
        goto jump_1;
      }
      
      TIMING_stop(tm_gmres_others, flop);
      
    } // loop; im
    
    nrm = Nmax;
    
  jump_1:
    
    rm[nrm] = rm[nrm] / rm[_IDX2D(nrm, nrm, Nmax)];
    
    for (int im=nrm-1; im>=1; im--) {
      al = rm[im];
      
      for (int jm=im+1; jm<=nrm; jm++) {
        al -= rm[_IDX2D(im, jm, Nmax)] * rm[jm];
      }
      
      rm[im] = al / rm[_IDX2D(im, im, Nmax)];
    }
    
    for (int im=1; im<=nrm; im++) {
      al = rm[im];
      
      ml_add_4_(d_p, size, &guide, &Nmax, &al, d_zm, &im, &flop);
    }
    
  } // loop; j
  
jump_2:
  
  SOR_2_SMA(IC, &d_zm[0], &d_vm[0], rhs_nrm, r0);
  
  
jump_4:
  
  if ( (isfin == 0) && (res < (rhs_nrm * eps_2 * eps_2)) )
  {
    res = 100.0 * rhs_nrm * eps_2 * eps_2;
  }
  
  if ( rm ) delete [] rm;
  if ( hg ) delete [] hg;
  if ( cs ) delete [] cs;
  if ( sn ) delete [] sn;
  
  TIMING_stop(tm_gmres_sor_sct, 0.0);
  // <<< Poisson Source section

}*/

// #################################################################
// Preconditioner
int LinearSolver::Fpreconditioner(REAL_TYPE* x, REAL_TYPE* b)
{
  REAL_TYPE omg = getOmega();
  
  int lc=0;                      /// ループカウント
  int lc_max = 4;
  
  // 前処理なし(コピー)
  if( lc_max == 0 )
  {
    blas_copy_(x, b, size, &guide);
    return lc;
  }
  
  for (lc=0; lc<lc_max; lc++)
  {
    Fsmoother(x, b, omg);
  }
  
  return lc;
}

// #################################################################
// Smoother
void LinearSolver::Fsmoother(REAL_TYPE* x, REAL_TYPE* b, REAL_TYPE omg)
{
  int ip = 0;
  
  if ( numProc > 1 )
  {
    ip = (head[0]+head[1]+head[2]+1) % 2;
  }
  else
  {
    ip = 0;
  }
  
  for (int color=0; color<2; color++)
  {
    blas_smoother_core_(x, b, bcp, &ip, &color, &omg, size, &guide);
    
    BC->OuterPBC(x, ensPeriodic);
    
    if ( C->EnsCompo.periodic == ON )
    {
      BC->InnerPBCperiodic(x, bcd);
    }
    
    if ( numProc > 1 )
    {
      if (getSyncMode() == comm_sync )
      {
        if ( paraMngr->BndCommS3D(x, size[0], size[1], size[2], guide, 1) != CPM_SUCCESS ) {
          Exit(0);
        }
      }
      else
      {
        int ireq[12];
        sma_comm_     (x, size, &guide, &color, &ip, cf_sz, cf_x, cf_y, cf_z, ireq, nID);
        sma_comm_wait_(x, size, &guide, &color, &ip, cf_sz, cf_x, cf_y, cf_z, ireq);
      }
    }
  }
}

// #################################################################
// Dot
void LinearSolver::Fdot(REAL_TYPE* xy, REAL_TYPE* x, REAL_TYPE* y)
{
  blas_dot_(xy, x, y, size, &guide);
  
  if ( numProc > 1 )
  {
    REAL_TYPE xy_tmp = *xy;
    if  ( paraMngr->Allreduce(&xy_tmp, xy, 1, MPI_SUM) != CPM_SUCCESS )Exit(0);
  }
}
