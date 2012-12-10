!********************************************************************
!
!   FFV : Frontflow / violet Cartesian
!
!   Copyright (c) 2012 All right reserved.
!
!   Institute of Industrial Science, University of Tokyo, Japan.
!
!********************************************************************

!> @file   vobc_wall.h
!! @brief  pvec_vobc_wall core
!! @author kero
!<

! real                                                      ::  Up0, Vp0, Wp0
! real                                                      ::  u_bc, v_bc, w_bc

  Up0 = v0(i,j,k,1)
  Vp0 = v0(i,j,k,2)
  Wp0 = v0(i,j,k,3)

  wv(i,j,k,1) = wv(i,j,k,1) + (u_bc - Up0) * dh2
  wv(i,j,k,2) = wv(i,j,k,2) + (v_bc - Vp0) * dh2
  wv(i,j,k,3) = wv(i,j,k,3) + (w_bc - Wp0) * dh2

! 9 flop
