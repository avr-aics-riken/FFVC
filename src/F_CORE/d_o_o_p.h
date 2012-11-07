!********************************************************************
!
!   FFV : Frontflow / violet Cartesian
!
!   Copyright (c) 2012 All right reserved.
!
!   Institute of Industrial Science, University of Tokyo, Japan. 
!
!********************************************************************

!> @file   d_o_o_p.h
!! @brief  データロードの共通部分
!! @author kero
!<

!    real                                                        ::  u_ref, v_ref, w_ref, flop
!    real                                                        ::  b_w, b_e, b_s, b_n, b_b, b_t
!    real                                                        ::  Up0, Ue0, Uw0, Vp0, Vs0, Vn0, Wp0, Wb0, Wt0
!    real                                                        ::  Ue, Uw, Vn, Vs, Wt, Wb
!    real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  v0
!    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bv

      ! 各セルフェイス位置の変数ロード
      Up0 = v0(1, i  ,j  ,k  )
      Vp0 = v0(2, i  ,j  ,k  )
      Wp0 = v0(3, i  ,j  ,k  )
      Uw0 = v0(1, i-1,j  ,k  )
      Ue0 = v0(1, i+1,j  ,k  )
      Vs0 = v0(2, i  ,j-1,k  )
      Vn0 = v0(2, i  ,j+1,k  )
      Wb0 = v0(3, i  ,j  ,k-1)
      Wt0 = v0(3, i  ,j  ,k+1)

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

      ! real*6 + 6*6 = 42 flops
