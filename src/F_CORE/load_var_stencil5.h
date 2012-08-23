!********************************************************************
!
!   FFV : Frontflow / violet Cartesian
!
!   Copyright (c) 2012 All right reserved.
!
!   Institute of Industrial Science, The University of Tokyo, Japan. 
!
!********************************************************************

!> @file   load_var_stencil5.h
!! @brief  共通のデータロード部分
!! @author kero
!<

!   real                                                        ::  Up0, Ue1, Ue2, Uw1, Uw2, Us1, Us2, Un1, Un2, Ub1, Ub2, Ut1, Ut2
!   real                                                        ::  Vp0, Ve1, Ve2, Vw1, Vw2, Vs1, Vs2, Vn1, Vn2, Vb1, Vb2, Vt1, Vt2
!   real                                                        ::  Wp0, We1, We2, Ww1, Ww2, Ws1, Ws2, Wn1, Wn2, Wb1, Wb2, Wt1, Wt2
!   real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  v

      ! 各軸方向5点の変数ロード
      Ub2 = v(1, i  ,j  ,k-2)
      Vb2 = v(2, i  ,j  ,k-2)
      Wb2 = v(3, i  ,j  ,k-2)

      Ub1 = v(1, i  ,j  ,k-1)
      Vb1 = v(2, i  ,j  ,k-1)
      Wb1 = v(3, i  ,j  ,k-1)

      Us2 = v(1, i  ,j-2,k  )
      Vs2 = v(2, i  ,j-2,k  )
      Ws2 = v(3, i  ,j-2,k  )

      Us1 = v(1, i  ,j-1,k  )
      Vs1 = v(2, i  ,j-1,k  )
      Ws1 = v(3, i  ,j-1,k  )

      Uw2 = v(1, i-2,j  ,k  )
      Vw2 = v(2, i-2,j  ,k  )
      Ww2 = v(3, i-2,j  ,k  )

      Uw1 = v(1, i-1,j  ,k  )
      Vw1 = v(2, i-1,j  ,k  )
      Ww1 = v(3, i-1,j  ,k  )

      Up0 = v(1, i  ,j  ,k  )
      Vp0 = v(2, i  ,j  ,k  )
      Wp0 = v(3, i  ,j  ,k  )

      Ue1 = v(1, i+1,j  ,k  )
      Ve1 = v(2, i+1,j  ,k  )
      We1 = v(3, i+1,j  ,k  )

      Ue2 = v(1, i+2,j  ,k  )
      Ve2 = v(2, i+2,j  ,k  )
      We2 = v(3, i+2,j  ,k  )

      Un1 = v(1, i  ,j+1,k  )
      Vn1 = v(2, i  ,j+1,k  )
      Wn1 = v(3, i  ,j+1,k  )

      Un2 = v(1, i  ,j+2,k  )
      Vn2 = v(2, i  ,j+2,k  )
      Wn2 = v(3, i  ,j+2,k  )

      Ut1 = v(1, i  ,j  ,k+1)
      Vt1 = v(2, i  ,j  ,k+1)
      Wt1 = v(3, i  ,j  ,k+1)

      Ut2 = v(1, i  ,j  ,k+2)
      Vt2 = v(2, i  ,j  ,k+2)
      Wt2 = v(3, i  ,j  ,k+2)
