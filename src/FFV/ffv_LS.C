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
// 収束判定　非Div反復
bool LinearSolver::Fcheck(double* var, const double b_l2, const double r0_l2)
{
  if ( getLS() == BiCGSTAB)
  {
    ;
  }
  else
  {
    // 自乗量の集約
    if ( numProc > 1 )
    {
      TIMING_start("A_R_Convergence");
      double tmp[3];
      tmp[0] = var[0];
      tmp[1] = var[1];
      tmp[2] = var[2];
      if ( paraMngr->Allreduce(tmp, var, 3, MPI_SUM) != CPM_SUCCESS ) Exit(0);
      TIMING_stop("A_R_Convergence", 6.0*numProc*sizeof(double) ); // 双方向 x ノード数　x 3
    }
    
    // L2 ノルム
    var[0] = sqrt(var[0]); // L2 of (x^{(m+1)} - x^{(m)})
    var[1] = sqrt(var[1]); // L2 of residual
    var[2] = sqrt(var[2]); // L2 of solution vector x^{(m)}
  }

  
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
  
  // 収束判定　性能測定モードでないときのみ収束判定を行う　残差が収束したら抜ける
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
  
  TIMING_start("Dot1");
  blas_dot1_(&xy, x, bcp, size, &guide, &flop_count);
  TIMING_stop("Dot1", flop_count);
  
  if ( numProc > 1 )
  {
    TIMING_start("A_R_Dot");
    double xy_tmp = xy;
    if  ( paraMngr->Allreduce(&xy_tmp, &xy, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
    TIMING_stop("A_R_Dot", 2.0*numProc*sizeof(double) );
  }
  
  return xy;
}


// #################################################################
double LinearSolver::Fdot2(REAL_TYPE* x, REAL_TYPE* y)
{
  double flop_count=0.0;          /// 浮動小数点演算数
  double xy = 0.0;
  
  TIMING_start("Dot2");
  blas_dot2_(&xy, x, y, bcp, size, &guide, &flop_count);
  TIMING_stop("Dot2", flop_count);
  
  if ( numProc > 1 )
  {
    TIMING_start("A_R_Dot");
    double xy_tmp = xy;
    if  ( paraMngr->Allreduce(&xy_tmp, &xy, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
    TIMING_stop("A_R_Dot", 2.0*numProc*sizeof(double) );
  }
  
  return xy;
}


// #################################################################
void LinearSolver::Initialize(Control* C,
                              SetBC3D* BC,
                              int ModeT,
                              double f_comm,
                              PerfMonitor* PM,
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
  this->ModeTiming = ModeT;
  this->face_comm_size = f_comm;
  this->PM = PM;
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
int LinearSolver::PointSOR(REAL_TYPE* x, REAL_TYPE* b, const double b_l2, const double r0_l2)
{
  double flop_count=0.0;          /// 浮動小数点演算数
  REAL_TYPE omg = getOmega();     /// 加速係数
  double var[3];                  /// 誤差、残差、解
  int lc=0;                       /// ループカウント
  
  // x     圧力 p^{n+1}
  // b     RHS vector
  // bcp   ビットフラグ
  
  for (lc=1; lc<=getMaxIteration(); lc++)
  {
    var[0] = 0.0; // 誤差
    var[1] = 0.0; // 残差
    var[2] = 0.0; // 解
    
    // 反復処理
    TIMING_start("Poisson_PSOR");
    flop_count = 0.0;
    psor_(x, size, &guide, &omg, var, b, bcp, &flop_count);
    TIMING_stop("Poisson_PSOR", flop_count);
    
    
    // 境界条件
    TIMING_start("Poisson_BC");
    BC->OuterPBC(x, ensPeriodic);
    if ( C->EnsCompo.periodic == ON ) BC->InnerPBCperiodic(x, bcd);
    TIMING_stop("Poisson_BC", 0.0);
    
    
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
  int lc_max = 4;
  double dummy = 1.0;
  
  // 前処理なし(コピー)
  if ( getPrecondition() != ON )
  {
    TIMING_start("Blas_Copy");
    blas_copy_(x, b, size, &guide);
    TIMING_stop("Blas_Copy");
    return;
  }
  
  // 前処理
  SOR2_SMA(x, b, lc_max, dummy, dummy, false);
}



// #################################################################
int LinearSolver::SOR2_SMA(REAL_TYPE* x, REAL_TYPE* b, const int itrMax, const double b_l2, const double r0_l2, bool converge_check)
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
  
  
  for (lc=1; lc<=itrMax; lc++)
  {
    // 2色のマルチカラー(Red&Black)のセットアップ
    
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
    
    
    var[0] = 0.0; // 誤差
    var[1] = 0.0; // 残差
    var[2] = 0.0; // 解
    
    // 各カラー毎の間に同期, 残差は色間で積算する
    // R - color=0 / B - color=1
    for (int color=0; color<2; color++) {
      
      TIMING_start("Poisson_SOR2_SMA");
      flop_count = 0.0; // 色間で積算しない
      if ( getNaive()==OFF)
      {
        psor2sma_core_(x, size, &guide, &ip, &color, &omg, var, b, bcp, &flop_count);
      }
      else
      {
        psor2sma_naive_(x, size, &guide, &ip, &color, &omg, var, b, bcp, pni, &flop_count);
      }
      
      TIMING_stop("Poisson_SOR2_SMA", flop_count);
      
      
      // 境界条件
      TIMING_start("Poisson_BC");
      BC->OuterPBC(x, ensPeriodic);
      if ( C->EnsCompo.periodic == ON ) BC->InnerPBCperiodic(x, bcd);
      TIMING_stop("Poisson_BC", 0.0);
      
      
      // 同期処理
      if ( numProc > 1 )
      {
        TIMING_start("Sync_Poisson");
        
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
        TIMING_stop("Sync_Poisson", face_comm_size*0.5*sizeof(REAL_TYPE));
      }
    }
    
    if ( converge_check )
    {
      // 収束判定 varは自乗量
      if ( Fcheck(var, b_l2, r0_l2) == true ) break;
    }
    
  }
  
  return lc;
}



// #################################################################
// 反復変数の同期処理
void LinearSolver::SyncScalar(REAL_TYPE* d_class, const int num_layer)
{
  if ( numProc > 1 )
  {
    TIMING_start("Sync_Poisson");
    
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
    TIMING_stop("Sync_Poisson", face_comm_size*(double)num_layer*sizeof(REAL_TYPE));
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
int LinearSolver::PointSOR_4th(REAL_TYPE* x, REAL_TYPE* b, REAL_TYPE* u_sum, REAL_TYPE* w1, REAL_TYPE* w2, REAL_TYPE dt, REAL_TYPE dh, const double b_l2, const double r0_l2)
{
  double flop_count=0.0;          /// 浮動小数点演算数
  REAL_TYPE omg = getOmega();     /// 加速係数
  double var[3];                  /// 誤差、残差、解
  int lc=0;                       /// ループカウント
  int c;
  
  
  // x  Pressure
  // b  src for Poisson
  // w1 work, source of 1st and 2nd step iteration
  // w2 work, solution of 1st step iteration
  
  //blas_clear_(w2, size, &guide);
  
  SyncScalar(u_sum, 2);
  
  //src_trnc_(w1, size, &guide, &dt, &dh, u_sum, &flop_count);
  //SyncScalar(w1, 1);
  
  /*
   for (c=1; c<getMaxIteration(); c++)
   {
   var[0] = 0.0; // 誤差
   var[1] = 0.0; // 残差
   var[2] = 0.0; // 解
   
   // 反復処理
   //TIMING_start(tm_poi_PSOR);
   flop_count = 0.0;
   psor_unmask_(w2, size, &guide, &omg, var, w1, &flop_count);
   //TIMING_stop(tm_poi_PSOR, flop_count);
   
   
   // 境界条件 >> @todo w2に適した実装が必要　マスク部はゼロ
   //TIMING_start(tm_poi_BC);
   BC->OuterPBC(w2, ensPeriodic);
   if ( C->EnsCompo.periodic == ON ) BC->InnerPBCperiodic(w2, bcd);
   //TIMING_stop(tm_poi_BC, 0.0);
   
   
   // 同期処理
   SyncScalar(w2, 1);
   
   
   // 収束判定
   if ( Fcheck(var, b_l2, r0_l2) == true ) break;
   
   }
   
   lc = c;
   printf("1st : itr=%d err=%e res=%e\n", lc, var[0], var[1]);
   */
  
  src_1st_(w1, size, &guide, &dt, &dh, u_sum, b, &flop_count);
  //src_2nd_(w1, size, &guide, &dh, b, w2, &flop_count);
  
  
  for (c=1; c<getMaxIteration(); c++)
  {
    var[0] = 0.0; // 誤差
    var[1] = 0.0; // 残差
    var[2] = 0.0; // 解
    
    // 反復処理
    TIMING_start("Poisson_PSOR");
    flop_count = 0.0;
    psor_(x, size, &guide, &omg, var, w1, bcp, &flop_count);
    TIMING_stop("Poisson_PSOR", flop_count);
    
    
    // 境界条件
    TIMING_start("Poisson_BC");
    BC->OuterPBC(x, ensPeriodic);
    if ( C->EnsCompo.periodic == ON ) BC->InnerPBCperiodic(x, bcd);
    TIMING_stop("Poisson_BC", 0.0);
    
    
    // 同期処理
    SyncScalar(x, 1);
    
    
    // 収束判定
    if ( Fcheck(var, b_l2, r0_l2) == true ) break;
    
  }
  
  lc += c;
  
  return lc;
}

// #################################################################
// PBiCBSTAB 収束判定は残差
int LinearSolver::PBiCGstab(REAL_TYPE* x, REAL_TYPE* b, const double b_l2, const double r0_l2)
{
  double var[3];          /// 誤差, 残差, 解ベクトルのL2ノルム
  var[0] = var[1] = var[2] = 0.0;
  REAL_TYPE omg = getOmega();
  double flop = 0.0;
  
  TIMING_start("Blas_Clear");
  FBUtility::initS3D(pcg_r , size, guide, 0.0);
  FBUtility::initS3D(pcg_p , size, guide, 0.0);
  FBUtility::initS3D(pcg_r0, size, guide, 0.0);
  FBUtility::initS3D(pcg_p_, size, guide, 0.0);
  FBUtility::initS3D(pcg_q , size, guide, 0.0);
  FBUtility::initS3D(pcg_s , size, guide, 0.0);
  FBUtility::initS3D(pcg_s_, size, guide, 0.0);
  FBUtility::initS3D(pcg_t_, size, guide, 0.0);
  TIMING_stop("Blas_Clear", 0.0, 8);
  
  TIMING_start("Blas_Residual");
  flop = 0.0;
  blas_calc_rk_(pcg_r, x, b, bcp, size, &guide, &flop);    // (1)
  TIMING_stop("Blas_Residual", flop);
  
  SyncScalar(pcg_r, 1);
  
  TIMING_start("Blas_Copy");
  blas_copy_(pcg_r0, pcg_r, size, &guide);                 // (2)
  TIMING_stop("Blas_Copy");
  
  double rr0 = 1.0;
  double alpha = 0.0;
  double gamma  = 1.0;
  double gamman = -gamma;
  int lc=0;                      /// ループカウント
  
  for (lc=1; lc<getMaxIteration(); lc++)
  {
    double rr1 = 0.0;
    rr1 = Fdot2(pcg_r, pcg_r0);                             // (4)
    
    if( fabs(rr1) < FLT_MIN )                              // (5)
    {
      lc = 0;
      break;
    }
    
    if( lc == 1 )
    {
      TIMING_start("Blas_Copy");
      blas_copy_(pcg_p, pcg_r, size, &guide);              // (7)
      TIMING_stop("Blas_Copy");
    }
    else
    {
      double beta = rr1/rr0*alpha/gamma;                // (9)
      TIMING_start("Blas_AXPY");
      flop = 0.0;
      blas_axpy_(pcg_p, pcg_q, &gamman, size, &guide, &flop);     // (10)
      TIMING_stop("Blas_AXPY", flop);
      
      TIMING_start("Blas_XPAY");
      flop = 0.0;
      blas_xpay_(pcg_p, pcg_r, &beta, size, &guide, &flop);    // (10)
      TIMING_stop("Blas_XPAY", flop);
    }
    SyncScalar(pcg_p, 1);
    
    TIMING_start("Blas_Clear");
    FBUtility::initS3D(pcg_p_, size, guide, 0.0);
    TIMING_stop("Blas_Clear");
    
    Preconditioner(pcg_p_, pcg_p);                    // (12)
    
    TIMING_start("Blas_AX");
    flop = 0.0;
    blas_calc_ax_(pcg_q, pcg_p_, bcp, size, &guide, &flop); // (13)
    TIMING_stop("Blas_AX", flop);
    
    double q_r0 = 0.0;
    q_r0 = Fdot2(pcg_q, pcg_r0);                           // (14)
    
    alpha  = rr1/q_r0;
    double alphan = -alpha;
    TIMING_start("Blas_AXPYZ");
    flop = 0.0;
    blas_axpyz_(pcg_s, pcg_q, pcg_r, &alphan, size, &guide, &flop); // (15)
    TIMING_stop("Blas_AXPYZ", flop);
    
    SyncScalar(pcg_s, 1);
    
    TIMING_start("Blas_Clear");
    FBUtility::initS3D(pcg_s_, size, guide, 0.0);
    TIMING_stop("Blas_Clear");
    
    Preconditioner(pcg_s_, pcg_s);                    // (17)
    
    TIMING_start("Blas_AX");
    flop = 0.0;
    blas_calc_ax_(pcg_t_, pcg_s_, bcp, size, &guide, &flop);     // (18)
    TIMING_stop("Blas_AX", flop);
    
    double t_s = 0.0;
    t_s = Fdot2(pcg_t_, pcg_s);                             // (19a)
    
    double t_t_ = 0.0;
    t_t_ = Fdot1(pcg_t_);                           // (19b)
    
    gamma  = t_s/t_t_;                                     // (19)
    gamman = -gamma;
    
    TIMING_start("Blas_AXPBYPZ");
    flop = 0.0;
    blas_axpbypz_(x, pcg_p_, pcg_s_, &alpha , &gamma, size, &guide, &flop);  // (20)
    TIMING_stop("Blas_AXPBYPZ", flop);
    
    TIMING_start("Blas_AXPYZ");
    flop = 0.0;
    blas_axpyz_  (pcg_r, pcg_t_, pcg_s, &gamman, size, &guide, &flop);            // (21)
    TIMING_stop("Blas_AXPYZ", flop);
    
    double rr = 0.0;
    rr = Fdot1(pcg_r);
    
    var[1] = rr;
    
    if ( Fcheck(var, b_l2, r0_l2) == true ) break;
    
    rr0 = rr1;
  }
  
  
  TIMING_start("Poisson_BC");
  BC->OuterPBC(x, ensPeriodic);
  if ( C->EnsCompo.periodic == ON )
  {
    BC->InnerPBCperiodic(x, bcd);
  }
  TIMING_stop("Poisson_BC");
  
  
  SyncScalar(x, 1);
  
  return lc;
}

// #################################################################
// PBiCBSTAB 収束判定は残差
int LinearSolver::RC_sor(REAL_TYPE* x, REAL_TYPE* pos_rhs, REAL_TYPE* b, int* bcp, const double r0_l2)
{
  double var[3];          /// 誤差, 残差, 解ベクトルのL2ノルム
  var[0] = var[1] = var[2] = 0.0;
  REAL_TYPE res = 0.0;
  REAL_TYPE omg = getOmega();
  double flop = 0.0;
  
  int iparam[10];
  double rparam[10];
  int err[2];
  int nrc_max = 20;
  int iter_max = 100;
  int isneum = 0;
  int i_inner, n_inner;
  int mrc, i_iter;
  double al, t_eps, f_v;
  REAL_TYPE e, ee;
  
  if (iparam[7] < 0) return 0;
  
  err[0] = 0;
  iparam[0] = 0;
  iparam[1] = 0;
  iparam[2] = 0;
  //    n_iter = min(iparam[7], iter_max)
  //    nrc    = min(iparam[7], nrc_max)
  int oki = 5;
  int step = 10;
  int n_iter = iter_max;
  int nrc    = nrc_max;

  
  if (sizeof(REAL_TYPE) == 4)
  {
    e = 2.4e-7;
    t_eps = 0.995*rparam[0]*rparam[0];
  }
  else
  {
    e = 4.4e-16;
    t_eps = 0.999*rparam[0]*rparam[0];
  }

  ee = 1.0 + e;
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int lg = 1 - guide;
  int ig = ix + guide;
  int jg = jx + guide;
  int kg = kx + guide;
  int tneum = 0;
  
  if ( (isneum !=0 ) && (guide > 0) )
  {
    tneum = 1;
  }
  
  //allocate(xt(lg:ig,lg:jg,lg:kg),      yt(lg:ig,lg:jg,lg:kg), rest(lg:ig,lg:jg,lg:kg), res_fct(lg:ig,lg:jg,lg:kg))
  //allocate(xrc(ix,jx,kx,nrc), yrc(ix,jx,kx,nrc))
  //allocate(scalprod(nrc*nrc), alph(nrc), stolb(nrc), matr(nrc*nrc))
  
  int iter = 1;
  
  double b_l2 = 0.0;
  
  //rc_calc_b_(b, pos_rhs, bcp, size, &guide, &b_l2, &flop);
  
  
  
}

