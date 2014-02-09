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

!> @file   load_var_stencil5.h
!! @brief  共通のデータロード部分
!! @author aics
!<

!   real                                                        ::  Up0, Ue1, Ue2, Uw1, Uw2, Us1, Us2, Un1, Un2, Ub1, Ub2, Ut1, Ut2
!   real                                                        ::  Vp0, Ve1, Ve2, Vw1, Vw2, Vs1, Vs2, Vn1, Vn2, Vb1, Vb2, Vt1, Vt2
!   real                                                        ::  Wp0, We1, We2, Ww1, Ww2, Ws1, Ws2, Wn1, Wn2, Wb1, Wb2, Wt1, Wt2
!   real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)   ::  v

      ! 各軸方向5点の変数ロード
      Ub2 = v(i  ,j  ,k-2, 1)
      Ub1 = v(i  ,j  ,k-1, 1)
      Us2 = v(i  ,j-2,k  , 1)
      Us1 = v(i  ,j-1,k  , 1)
      Uw2 = v(i-2,j  ,k  , 1)
      Uw1 = v(i-1,j  ,k  , 1)
      Up0 = v(i  ,j  ,k  , 1)
      Ue1 = v(i+1,j  ,k  , 1)
      Ue2 = v(i+2,j  ,k  , 1)
      Un1 = v(i  ,j+1,k  , 1)
      Un2 = v(i  ,j+2,k  , 1)
      Ut1 = v(i  ,j  ,k+1, 1)
      Ut2 = v(i  ,j  ,k+2, 1)

      Vb2 = v(i  ,j  ,k-2, 2)
      Vb1 = v(i  ,j  ,k-1, 2)
      Vs2 = v(i  ,j-2,k  , 2)
      Vs1 = v(i  ,j-1,k  , 2)
      Vw2 = v(i-2,j  ,k  , 2)
      Vw1 = v(i-1,j  ,k  , 2)
      Vp0 = v(i  ,j  ,k  , 2)
      Ve1 = v(i+1,j  ,k  , 2)
      Ve2 = v(i+2,j  ,k  , 2)
      Vn1 = v(i  ,j+1,k  , 2)
      Vn2 = v(i  ,j+2,k  , 2)
      Vt1 = v(i  ,j  ,k+1, 2)
      Vt2 = v(i  ,j  ,k+2, 2)

      Wb2 = v(i  ,j  ,k-2, 3)
      Wb1 = v(i  ,j  ,k-1, 3)
      Ws2 = v(i  ,j-2,k  , 3)
      Ws1 = v(i  ,j-1,k  , 3)
      Ww2 = v(i-2,j  ,k  , 3)
      Ww1 = v(i-1,j  ,k  , 3)
      Wp0 = v(i  ,j  ,k  , 3)
      We1 = v(i+1,j  ,k  , 3)
      We2 = v(i+2,j  ,k  , 3)
      Wn1 = v(i  ,j+1,k  , 3)
      Wn2 = v(i  ,j+2,k  , 3)
      Wt1 = v(i  ,j  ,k+1, 3)
      Wt2 = v(i  ,j  ,k+2, 3)
