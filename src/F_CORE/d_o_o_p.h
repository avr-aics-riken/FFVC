!###################################################################################
!
! FFV-C
! Frontflow / violet Cartesian
!
!
! Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
! All rights reserved.
!
! Copyright (c) 2011-2014 Institute of Industrial Science, The University of Tokyo.
! All rights reserved.
!
! Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
! All rights reserved.
!
!###################################################################################

!> @file   d_o_o_p.h
!! @brief  データロードの共通部分
!! @author kero
!<

!    real                                                        ::  b_w, b_e, b_s, b_n, b_b, b_t, b_p
!    real                                                        ::  w_w, w_e, w_s, w_n, w_b, w_t
!    real                                                        ::  Ue, Uw, Vn, Vs, Wt, Wb
!    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)   ::  vf
!    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bv

      ! セルの状態 (0-solid / 1-fluid)
      b_p = real( ibits(bv(i  ,j  ,k  ), State, 1) )
      b_w = real( ibits(bv(i-1,j  ,k  ), State, 1) )
      b_e = real( ibits(bv(i+1,j  ,k  ), State, 1) )
      b_s = real( ibits(bv(i  ,j-1,k  ), State, 1) )
      b_n = real( ibits(bv(i  ,j+1,k  ), State, 1) )
      b_b = real( ibits(bv(i  ,j  ,k-1), State, 1) )
      b_t = real( ibits(bv(i  ,j  ,k+1), State, 1) )

      ! セル界面のフラグ (0-wall face / 1-fluid) > real*6+ 6= 12 flops
      w_e = real(b_e * b_p)
      w_w = real(b_w * b_p)
      w_n = real(b_n * b_p)
      w_s = real(b_s * b_p)
      w_t = real(b_t * b_p)
      w_b = real(b_b * b_p)

      ! 界面速度（スタガード位置） > 6 flops
      Uw = vf(i-1, j  , k  ,1)*w_w
      Ue = vf(i  , j  , k  ,1)*w_e
      Vs = vf(i  , j-1, k  ,2)*w_s
      Vn = vf(i  , j  , k  ,2)*w_n
      Wb = vf(i  , j  , k-1,3)*w_b
      Wt = vf(i  , j  , k  ,3)*w_t

      ! real*7 real*6 + 6 + 6 = 25 flops
