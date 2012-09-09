// #################################################################
//
// FFV : Frontflow / violet Cartesian
//
// Copyright (c) 2012 All right reserved.
//
// Institute of Industrial Science, The University of Tokyo, Japan.
//
// #################################################################

/**
 * @file   ffv_Gmres.C
 * @brief  FFV Class
 * @author kero
 */


#include "ffv.h"



void FFV::Gmres_SOR(ItrCtl* IC, double res_rhs, double& flop)
{
  int err[3];
  const int m_max = 15;
  const int iter_max = 100;
  const double fct1 = 1.0-4;
  const double fct2 = 6.0;
  double t_eps, beta, beta_1, res, eps_abs, res_abs;
  double r4[6];
  
  if ( iparam[8] < 0) return;
  
  err[1] = 0;
  int isfin  = 0;
  
  // n_iter = min(iparam(8), iter_max);
  // nrc    = min(iparam(8), nrc_max);
  
  int oki  = 4;
  int step = 6;
  int m       = m_max;
  int n_iter  = iter_max;
  
  if ( C.Mode.Precision == sizeof(float) )
  {
    t_eps = 0.995 * rparam[1] * rparam[1];
  }
  else
  {
    t_eps = 0.999 * rparam[1] * rparam[1];
  }
  
  
  int ix = size[0];
  int jx = size[2];
  int kx = size[3];
  int gd = guide;
  
  int lg = 1 - gd;
  int ig = ix + gd;
  int jg = jx + gd;
  int kg = kx + gd;
  
  if ( res_rhs < 1.0e-30 ) return;
  
  t_eps = res_rhs * t_eps;
  
  res = (double)SOR_2_SMA(IC);
  
  r4[1] = sqrt(res/res_rhs);
  
  if (res < t_eps)
  {
    printf("Final : %e\n", r4[1]);
    isfin = 1;
    goto 100;
  }

  if ( res > (fct2*t_eps) ) t_eps = res/fct2;
  
  eps_abs = 0.15 * t_eps;
  
  int iter = 0;
  int nrm  = 0;
  
  for (int i_iter=1; i_iter<=n_iter; i_iter++) {
    
    res = (double)SOR_2_SMA(IC);
    
    iparam[2]++;
    
    if (res < t_eps) goto 70;
    
    
    mv_prod_(d_yt, size, &guide, d_p, d_bcp, &flop);
    iparam[3]++;
    
    res_smpl_(res_abs);
    beta = sqrt(res_abs);
    if (beta < 1.0e-30) goto 60;
    beta_1 = 1.0/beta;
    
    
    bgm[1] = beta;
        
    for (int i=2; i<=m; i++) {
      bgm[i] = 0.0;
    }
    
    if (m > m_max) m = m_max;
    
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
}