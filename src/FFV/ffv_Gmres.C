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



double FFV::Gmres_SOR(ItrCtl* IC, double res_rhs, int *iparam, double *rparam)
{
  int err[3];
  const int m_max = 15;
  const int iter_max = 100;
  const double fct1 = 1.0-4;
  const double fct2 = 6.0;
  const double eps1 = 1.0e-30;
  double t_eps, beta, beta_1, res, eps_abs, res_abs, al;
  double r4[6];
  double flop=0.0;
  
  size_t s_length = (size[0]+2*guide) * (size[1]+2*guide) * (size[2]+2*guide);
  
  if ( iparam[8] < 0 ) return;
  
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
  
  
  double *rgm = new double[m_max, m_max];
  double *hgm = new double[m_max+1, m_max];
  double *cgm = new double[m_max];
  double *sgm = new double[m_max];
  double *bgm = new double[m_max];
  double *ygm = new double[m_max];
  
  int ix = size[0];
  int jx = size[2];
  int kx = size[3];
  int gd = guide;
  
  if ( res_rhs < eps1 ) return;
  
  t_eps = res_rhs * t_eps;
  
  res = SOR_2_SMA(IC);
  
  r4[1] = sqrt(res/res_rhs);
  
  if (res < t_eps)
  {
    printf("Final : %e\n", r4[1]);
    isfin = 1;
    goto jump_4;
  }

  if ( res > (fct2*t_eps) ) t_eps = res/fct2;
  
  eps_abs = 0.15 * t_eps;
  
  int iter = 0;
  int nrm  = 0;
  
  
  for (int i_iter=1; i_iter<=n_iter; i_iter++) {
    
    res = SOR_2_SMA(IC);
    iparam[2]++;
    
    if (res < t_eps) goto jump_3;
    
    
    mv_prod_(d_yt, size, &guide, d_p, d_bcp, &flop);
    iparam[3]++;
    
    res_smpl_(d_rest, size, &guide, &res_abs, d_ws, d_yt, d_bcp, &flop);
    
    double tmp = res_abs;
    if ( paraMngr->Allreduce(&tmp, &res_abs, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
    beta = sqrt(res_abs);
    
    if (beta < eps1) goto jump_2;
    
    beta_1 = 1.0/beta;
    
    multiply_(d_vm, size, &guide, &beta_1, d_rest, &flop);
    
    bgm[1] = beta;
    
    
    for (int i=2; i<=m; i++) {
      bgm[i] = 0.0;
    }
    
    if (m > m_max) m = m_max;
    
    
    for (int im=2; im<=m; im++) {
      iter++;
      
      copy_1_(d_rest, size, &guide, d_vm, &im);
      
#pragma omp parallel for firstprivate(s_length)
      for (size_t i=0; i<s_length; i++) {
        d_xt[i] = 0.0;
      }
      
      if ( iter <= oki )
      {
        n_inner = step;
      }
      else if ( iter <= 2*oki )
      {
        n_inner = 2*step;
      }
      else if ( iter <= 3*oki )
      {
        n_inner = 3*step;
      }
      else if ( iter <= 4*oki )
      {
        n_inner = 4*step;
      }
      else
      {
        n_inner = 5*step;
      }

      if (n_inner < 1) n_inner = 1;
      
      for (int i_inner=1; i_inner<=n_inner; i_inner++) {
        
        res = SOR_2_SMA(IC);
        iparam[2]++;
        
      }
      
      copy_2_(d_zm, size, &guide, d_xt, &im); // copy range
      
      mv_prod_(d_xt, size, &guide, d_p, d_bcp, &flop);
      iparam[3]++;
      
      for (int km=1; km<=im; km++) {
        al = 0.0;
        ml_add_1_(&al, size, &guide, d_vm, d_yt, &km, &flop);
        
        hgm[km, im] = al;
        double r_al = -al;
        ml_add_3_(d_yt, size, &guide, &r_al, d_vm, &km, &flop);
      }
      
      al =0.0;
      ml_add_2_(&al, size, &guide, d_yt, &flop);
      
      
      hgm[im+1, im] = sqrt(al);
      
      if (hgm[im+1,im] < eps1)
      {
        nrm = im-1;
        goto jump_1;
      }
      
      if (im < m)
      {
        al = 1.0/hgm[im+1, im];
        multiply_(d_vm, size, &guide, &al, d_yt, &flop);
      }
      
      rgm[1,im] = hgm[1,im];
      
      for (int km=2; km<=im; km++) {
        rgm[km  ,im] = cgm[km-1]*hgm[km,im] - sgm[km-1]*rgm[km-1,im];
        rgm[km-1,im] = sgm[km-1]*hgm[km,im] + cgm[km-1]*rgm[km-1,im];
      }
      
      al = sqrt(rgm[im,im]*rgm[im,im] + hgm[im+1,im]*hgm[im+1,im]);
      
      if (al < eps1)
      {
        nrm = im - 1;
        goto jump_1;
      }
      
      cgm[im]    = rgm[im,  im]/al;
      sgm[im]    = hgm[im+1,im]/al;
      rgm[im,im] = cgm[im]*rgm[im,im] + sgm[im]*hgm[im+1,im];
      bgm[im+1]  = - sgm[im]*bgm[im];
      bgm[im]    = cgm[im]*bgm[im];
      al = bgm[im+1]*bgm[im+1];
      iparam[1]++;
      
      if (al < eps_abs)
      {
        nrm = im;
        goto jump_1;
      }
      
    } // loop; im
    
    nrm = m;
    
jump_1:
    
    ygm[nrm] = bgm[nrm]/rgm[nrm,nrm];
    
    
    for (int im=nrm-1; im>=1; -1) {
      al = bgm[im];
      
      for (int jm=im+1; jm<=nrm; jm++) {
        al -= rgm[im,jm]*ygm[jm];
      }
      
      ygm[im] = al/rgm[im,im];
    }
    
    for (int im=1, im<=nrm; i++) {
      al = ygm[im];
      
      ml_add_3_(d_p, size, &guide, &al, d_zm, &im, &flop);
    }
    
  }
  
jump_2:
  
  res = SOR_2_SMA(IC);
  
  
jump_3:
  
  
jump_4:
  
  if ( (isfin==0) && (res < (res_rhs*rparam[1]*rparam[1])) )
  {
    res = 100.0 * res_rhs * rparam[1] * rparam[1];
  }
  
  
  if ( rgm ) delete [] rgm;
  if ( hgm ) delete [] hgm;
  if ( cgm ) delete [] cgm;
  if ( sgm ) delete [] sgm;
  if ( bgm ) delete [] bgm;
  if ( ygm ) delete [] ygm;
  
}
