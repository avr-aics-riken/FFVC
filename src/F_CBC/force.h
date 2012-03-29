!   *********************************************************
!
!   SPHERE - Skeleton for PHysical and Engineering REsearch
!  
!   Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
!
!   *********************************************************

!> @file muscl_x.h
!> @brief subroutines for CBC
!> @author keno, FSI Team, VCAD, RIKEN 

! 179 flop

  idx = bd(i,j,k)
  q_p = 0.0
  q_w = 0.0
  q_e = 0.0
  q_s = 0.0
  q_n = 0.0
  q_b = 0.0
  q_t = 0.0

  if ( ibits(idx,             0, bitw_6) == odr ) q_p = 1.0
  if ( ibits(bd(i-1,j  ,k  ), 0, bitw_6) == odr ) q_w = 1.0
  if ( ibits(bd(i+1,j  ,k  ), 0, bitw_6) == odr ) q_e = 1.0
  if ( ibits(bd(i  ,j-1,k  ), 0, bitw_6) == odr ) q_s = 1.0
  if ( ibits(bd(i  ,j+1,k  ), 0, bitw_6) == odr ) q_n = 1.0
  if ( ibits(bd(i  ,j  ,k-1), 0, bitw_6) == odr ) q_b = 1.0
  if ( ibits(bd(i  ,j  ,k+1), 0, bitw_6) == odr ) q_t = 1.0

  b0 = vf(i,j,k)
  bw = 0.5*( b0 + vf(i-1,j  ,k  ) )
  be = 0.5*( b0 + vf(i+1,j  ,k  ) )
  bs = 0.5*( b0 + vf(i  ,j-1,k  ) )
  bn = 0.5*( b0 + vf(i  ,j+1,k  ) )
  bb = 0.5*( b0 + vf(i  ,j  ,k-1) )
  bt = 0.5*( b0 + vf(i  ,j  ,k+1) )

  u_p = wk(1, ii  ,jj  ,kk  ) - u_ref
  v_p = wk(2, ii  ,jj  ,kk  ) - v_ref
  w_p = wk(3, ii  ,jj  ,kk  ) - w_ref
  g_p = sqrt(u_p*u_p + v_p*v_p + w_p*w_p)
  d_p = c1*g_p*g_p + c2*g_p + c3
  if (g_p < ep) d_p = c4*g_p*g_p
  Fx_p = -sign(1.0, u_p) * d_p * nx * q_p
  Fy_p = -sign(1.0, v_p) * d_p * ny * q_p
  Fz_p = -sign(1.0, w_p) * d_p * nz * q_p

  u_w = wk(1, ii-1,jj  ,kk  ) - u_ref
  v_w = wk(2, ii-1,jj  ,kk  ) - v_ref
  w_w = wk(3, ii-1,jj  ,kk  ) - w_ref
  g_w = sqrt(u_w*u_w + v_w*v_w + w_w*w_w)
  d_w = c1*g_w*g_w + c2*g_w + c3
  if (g_w < ep) d_w = c4*g_w*g_w
  Fx_w = -sign(1.0, u_w) * d_w * nx * q_w

  u_e = wk(1, ii+1,jj  ,kk  ) - u_ref
  v_e = wk(2, ii+1,jj  ,kk  ) - v_ref
  w_e = wk(3, ii+1,jj  ,kk  ) - w_ref
  g_e = sqrt(u_e*u_e + v_e*v_e + w_e*w_e)
  d_e = c1*g_e*g_e + c2*g_e + c3
  if (g_e < ep) d_e = c4*g_e*g_e
  Fx_e = -sign(1.0, u_e) * d_e * nx * q_e

  u_s = wk(1, ii  ,jj-1,kk  ) - u_ref
  v_s = wk(2, ii  ,jj-1,kk  ) - v_ref
  w_s = wk(3, ii  ,jj-1,kk  ) - w_ref
  g_s = sqrt(u_s*u_s + v_s*v_s + w_s*w_s)
  d_s = c1*g_s*g_s + c2*g_s + c3
  if (g_s < ep) d_s = c4*g_s*g_s
  Fy_s = -sign(1.0, v_s) * d_s * ny * q_s

  u_n = wk(1, ii  ,jj+1,kk  ) - u_ref
  v_n = wk(2, ii  ,jj+1,kk  ) - v_ref
  w_n = wk(3, ii  ,jj+1,kk  ) - w_ref
  g_n = sqrt(u_n*u_n + v_n*v_n + w_n*w_n)
  d_n = c1*g_n*g_n + c2*g_n + c3
  if (g_n < ep) d_n = c4*g_n*g_n
  Fy_n = -sign(1.0, v_n) * d_n * ny * q_n

  u_b = wk(1, ii  ,jj  ,kk-1) - u_ref
  v_b = wk(2, ii  ,jj  ,kk-1) - v_ref
  w_b = wk(3, ii  ,jj  ,kk-1) - w_ref
  g_b = sqrt(u_b*u_b + v_b*v_b + w_b*w_b)
  d_b = c1*g_b*g_b + c2*g_b + c3
  if (g_b < ep) d_b = c4*g_b*g_b
  Fz_b = -sign(1.0, w_b) * d_b * nz * q_b

  u_t = wk(1, ii  ,jj  ,kk+1) - u_ref
  v_t = wk(2, ii  ,jj  ,kk+1) - v_ref
  w_t = wk(3, ii  ,jj  ,kk+1) - w_ref
  g_t = sqrt(u_t*u_t + v_t*v_t + w_t*w_t)
  d_t = c1*g_t*g_t + c2*g_t + c3
  if (g_t < ep) d_t = c4*g_t*g_t
  Fz_t = -sign(1.0, w_t) * d_t * nz * q_t

  re = 0.5*(Fx_p + Fx_e)
  rw = 0.5*(Fx_p + Fx_w)
  rn = 0.5*(Fy_p + Fy_n)
  rs = 0.5*(Fy_p + Fy_s)
  rt = 0.5*(Fz_p + Fz_t)
  rb = 0.5*(Fz_p + Fz_b)