!   *********************************************************
!
!   SPHERE - Skeleton for PHysical and Engineering REsearch
!  
!   Copyright (c) RIKEN, Japan. All right reserved. 2004-2011
!
!   *********************************************************

!> @file d_o_o_p.h
!> @brief subroutines for CBC
!> @author keno, FSI Team, VCAD, RIKEN

!    real                                                        ::  u_ref, v_ref, w_ref, flop
!    real                                                        ::  b_w, b_e, b_s, b_n, b_b, b_t
!    real                                                        ::  Up0, Ue0, Uw0, Vp0, Vs0, Vn0, Wp0, Wb0, Wt0
!    real                                                        ::  Ue, Uw, Vn, Vs, Wt, Wb
!    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)   ::  v0
!    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bv

      ! 各セルフェイス位置の変数ロード
      Up0 = v0(i  ,j  ,k  ,1)
      Uw0 = v0(i-1,j  ,k  ,1)
      Ue0 = v0(i+1,j  ,k  ,1)

      Vp0 = v0(i  ,j  ,k  ,2)
      Vs0 = v0(i  ,j-1,k  ,2)
      Vn0 = v0(i  ,j+1,k  ,2)

      Wp0 = v0(i  ,j  ,k  ,3)
      Wb0 = v0(i  ,j  ,k-1,3)
      Wt0 = v0(i  ,j  ,k+1,3)

      ! 壁面による修正 (0-solid / 1-fluid)
      b_w = real( ibits(bv(i-1,j  ,k  ), State, 1) )
      b_e = real( ibits(bv(i+1,j  ,k  ), State, 1) )
      b_s = real( ibits(bv(i  ,j-1,k  ), State, 1) )
      b_n = real( ibits(bv(i  ,j+1,k  ), State, 1) )
      b_b = real( ibits(bv(i  ,j  ,k-1), State, 1) )
      b_t = real( ibits(bv(i  ,j  ,k+1), State, 1) )

      Uw = 0.5*( Up0 + Uw0 )*b_w + (1.0-b_w)*u_ref
      Ue = 0.5*( Up0 + Ue0 )*b_e + (1.0-b_e)*u_ref
      Vs = 0.5*( Vp0 + Vs0 )*b_s + (1.0-b_s)*v_ref
      Vn = 0.5*( Vp0 + Vn0 )*b_n + (1.0-b_n)*v_ref
      Wb = 0.5*( Wp0 + Wb0 )*b_b + (1.0-b_b)*w_ref
      Wt = 0.5*( Wp0 + Wt0 )*b_t + (1.0-b_t)*w_ref

      flop = flop + 36.0
