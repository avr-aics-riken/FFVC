!********************************************************************
!
!   FFV : Frontflow / violet Cartesian
!
!   Copyright (c) 2012-2013  All right reserved.
!
!   Institute of Industrial Science, University of Tokyo, Japan. 
!
!********************************************************************

!> @file   ffv_bc_outer.f90
!! @brief  外部境界条件
!! @author kero
!<


!> ********************************************************************
!! @brief 外部速度境界条件による対流項と粘性項の流束の修正
!! @param [out] wv   疑似ベクトルの空間項の評価値
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル長
!! @param [in]  dh   格子幅
!! @param [in]  rei  Reynolds数の逆数
!! @param [in]  v0   セルセンター速度ベクトル（n-step）
!! @param [in]  bv   BCindex V
!! @param [in]  face 外部境界処理のときの面番号
!! @param [in]  cf   流出速度
!! @param [out] flop 浮動小数点演算数
!<
    subroutine vobc_pv_oflow (wv, sz, g, dh, rei, v0, bv, face, cf, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, g, face
    integer                                                   ::  ix, jx, kx
    integer, dimension(3)                                     ::  sz
    double precision                                          ::  flop
    real                                                      ::  Up0, Vp0, Wp0
    real                                                      ::  dh, dh1, dh2, rei, c, cf, m
    real                                                      ::  fu, fv, fw, EX, EY, EZ
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v0, wv
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bv
    
    ix = sz(1)
    jx = sz(2)
    kx = sz(3)

    dh1= 1.0/dh
    dh2= rei*dh1*dh1

    c = cf
    m =0.0
    
    flop = flop + 13.0d0 ! DP 15 flops
    
!$OMP PARALLEL REDUCTION(+:m) &
!$OMP FIRSTPRIVATE(ix, jx, kx, face) &
!$OMP FIRSTPRIVATE(dh1, dh2, cf, c) &
!$OMP PRIVATE(Up0, Vp0, Wp0) &
!$OMP PRIVATE(fu, fv, fw, EX, EY, EZ)
    
    FACES : select case (face)
    case (X_minus)

      if ( cf>0.0 ) c=0.0

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do j=1,jx
        if ( ibits(bv(1, j, k), bc_face_W, bitw_5) == obc_mask ) then
          Up0 = v0(1, j, k, 1)
          Vp0 = v0(1, j, k, 2)
          Wp0 = v0(1, j, k, 3)

          fu  = c * Up0
          fv  = c * Vp0
          fw  = c * Wp0
          
          EX = v0(0, j, k, 1) - Up0
          EY = v0(0, j, k, 2) - Vp0
          EZ = v0(0, j, k, 3) - Wp0

          wv(1, j, k, 1) = wv(1, j, k, 1) + ( fu*dh1 + EX*dh2 )
          wv(1, j, k, 2) = wv(1, j, k, 2) + ( fv*dh1 + EY*dh2 )
          wv(1, j, k, 3) = wv(1, j, k, 3) + ( fw*dh1 + EZ*dh2 )
          m = m + 1.0
        end if
      end do
      end do
!$OMP END DO


    case (X_plus)

      if ( cf<0.0 ) c=0.0

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do j=1,jx
        if ( ibits(bv(ix, j, k), bc_face_E, bitw_5) == obc_mask ) then
          Up0 = v0(ix, j, k, 1)
          Vp0 = v0(ix, j, k, 2)
          Wp0 = v0(ix, j, k, 3)

          fu  = c * Up0
          fv  = c * Vp0
          fw  = c * Wp0

          EX = v0(ix+1, j, k, 1) - Up0
          EY = v0(ix+1, j, k, 2) - Vp0
          EZ = v0(ix+1, j, k, 3) - Wp0

          wv(ix, j, k, 1) = wv(ix, j, k, 1) + ( -fu*dh1 + EX*dh2 )
          wv(ix, j, k, 2) = wv(ix, j, k, 2) + ( -fv*dh1 + EY*dh2 )
          wv(ix, j, k, 3) = wv(ix, j, k, 3) + ( -fw*dh1 + EZ*dh2 )
          m = m + 1.0
        end if
      end do
      end do
!$OMP END DO

      
    case (Y_minus)

      if ( cf>0.0 ) c=0.0

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do i=1,ix
        if ( ibits(bv(i, 1, k), bc_face_S, bitw_5) == obc_mask ) then
          Up0 = v0(i, 1, k, 1)
          Vp0 = v0(i, 1, k, 2)
          Wp0 = v0(i, 1, k, 3)

          fu  = c * Up0
          fv  = c * Vp0
          fw  = c * Wp0
        
          EX = v0(i, 0, k, 1) - Up0
          EY = v0(i, 0, k, 2) - Vp0
          EZ = v0(i, 0, k, 3) - Wp0

          wv(i, 1, k, 1) = wv(i, 1, k, 1) + ( fu*dh1 + EX*dh2 )
          wv(i, 1, k, 2) = wv(i, 1, k, 2) + ( fv*dh1 + EY*dh2 )
          wv(i, 1, k, 3) = wv(i, 1, k, 3) + ( fw*dh1 + EZ*dh2 )
          m = m + 1.0
        end if
      end do
      end do
!$OMP END DO

      
    case (Y_plus)

      if ( cf<0.0 ) c=0.0

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do i=1,ix
        if ( ibits(bv(i, jx, k), bc_face_N, bitw_5) == obc_mask ) then
          Up0 = v0(i, jx, k, 1)
          Vp0 = v0(i, jx, k, 2)
          Wp0 = v0(i, jx, k, 3)

          fu  = c * Up0
          fv  = c * Vp0
          fw  = c * Wp0

          EX = v0(i, jx+1, k, 1) - Up0
          EY = v0(i, jx+1, k, 2) - Vp0
          EZ = v0(i, jx+1, k, 3) - Wp0
          
          wv(i, jx, k, 1) = wv(i, jx, k, 1) + ( -fu*dh1 + EX*dh2 )
          wv(i, jx, k, 2) = wv(i, jx, k, 2) + ( -fv*dh1 + EY*dh2 )
          wv(i, jx, k, 3) = wv(i, jx, k, 3) + ( -fw*dh1 + EZ*dh2 )
          m = m + 1.0
        end if
      end do
      end do
!$OMP END DO

      
    case (Z_minus)

      if ( cf>0.0 ) c=0.0

!$OMP DO SCHEDULE(static)
      do j=1,jx
      do i=1,ix
        if ( ibits(bv(i, j, 1), bc_face_B, bitw_5) == obc_mask ) then
          Up0 = v0(i, j, 1, 1)
          Vp0 = v0(i, j, 1, 2)
          Wp0 = v0(i, j, 1, 3)

          fu  = c * Up0
          fv  = c * Vp0
          fw  = c * Wp0
          
          EX = v0(i, j, 0, 1) - Up0
          EY = v0(i, j, 0, 2) - Vp0
          EZ = v0(i, j, 0, 3) - Wp0

          wv(i, j, 1, 1) = wv(i, j, 1, 1) + ( fu*dh1 + EX*dh2 )
          wv(i, j, 1, 2) = wv(i, j, 1, 2) + ( fv*dh1 + EY*dh2 )
          wv(i, j, 1, 3) = wv(i, j, 1, 3) + ( fw*dh1 + EZ*dh2 )
          m = m + 1.0
        end if
      end do
      end do
!$OMP END DO

      
    case (Z_plus)

      if ( cf<0.0 ) c=0.0

!$OMP DO SCHEDULE(static)
      do j=1,jx
      do i=1,ix
        if ( ibits(bv(i, j, kx), bc_face_T, bitw_5) == obc_mask ) then
          Up0 = v0(i, j, kx, 1)
          Vp0 = v0(i, j, kx, 2)
          Wp0 = v0(i, j, kx, 3)

          fu  = c * Up0
          fv  = c * Vp0
          fw  = c * Wp0
          
          EX = v0(i, j, kx+1, 1) - Up0
          EY = v0(i, j, kx+1, 2) - Vp0
          EZ = v0(i, j, kx+1, 3) - Wp0
          
          wv(i, j, kx, 1) = wv(i, j, kx, 1) + ( -fu*dh1 + EX*dh2 )
          wv(i, j, kx, 2) = wv(i, j, kx, 2) + ( -fv*dh1 + EY*dh2 )
          wv(i, j, kx, 3) = wv(i, j, kx, 3) + ( -fw*dh1 + EZ*dh2 )
          m = m + 1.0
        endif
      end do
      end do
!$OMP END DO
      
    case default
    end select FACES
!$OMP END PARALLEL

    flop = flop + dble(m)*18.0d0
      
    return
    end subroutine vobc_pv_oflow
    
!> ********************************************************************
!! @brief 外部速度境界条件による対流項と粘性項の流束の修正
!! @param [out] wv   疑似ベクトルの空間項の評価値
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル長
!! @param [in]  dh   格子幅
!! @param [in]  rei  Reynolds数の逆数
!! @param [in]  v0   速度ベクトル（n-step）
!! @param [in]  bv   BCindex V
!! @param [in]  vec  指定する速度ベクトル
!! @param [in]  face 外部境界処理のときの面番号
!! @param [out] flop 浮動小数点演算数
!! @note vecには，流入条件のとき指定速度
!!  mskで部分的な速度を与える
!<
    subroutine vobc_pv_specv (wv, sz, g, dh, rei, v0, bv, vec, face, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, g, face
    integer                                                   ::  ix, jx, kx
    integer, dimension(3)                                     ::  sz
    double precision                                          ::  flop
    real                                                      ::  Up, Vp, Wp, Ur, Vr, Wr
    real                                                      ::  dh, dh1, dh2, EX, EY, EZ, rei
    real                                                      ::  fu, fv, fw, c, ac, msk
    real                                                      ::  u_bc, v_bc, w_bc, m
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v0, wv
    real, dimension(3)                                        ::  vec
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bv
    
    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    
    dh1= 1.0/dh
    dh2= rei*dh1*dh1
    
    ! u_bcは境界速度
    u_bc = vec(1)
    v_bc = vec(2)
    w_bc = vec(3)
    
    flop = flop + 13.0d0 ! DP 18 flop

    m = 0.0
    
!$OMP PARALLEL REDUCTION(+:m) &
!$OMP FIRSTPRIVATE(ix, jx, kx, u_bc, v_bc, w_bc, face) &
!$OMP FIRSTPRIVATE(dh1, dh2) &
!$OMP PRIVATE(Up, Vp, Wp, Ur, Vr, Wr) &
!$OMP PRIVATE(fu, fv, fw, EX, EY, EZ, c, ac, msk) &
!$OMP PRIVATE(i, j, k)

    FACES : select case (face)
    case (X_minus)

      i = 1
!$OMP DO SCHEDULE(static)
      do k=1,kx
      do j=1,jx
        if ( ibits(bv(i,j,k), bc_face_W, bitw_5) == obc_mask ) then
          Up = v0(i,j,k,1)
          Vp = v0(i,j,k,2)
          Wp = v0(i,j,k,3)
          
          Ur = u_bc
          Vr = v_bc
          Wr = w_bc
          c  = u_bc
          ac = abs(c)

          EX = Up - Ur
          EY = Vp - Vr
          EZ = Wp - Wr

          fu = 0.5*(c*(Up+Ur) - ac*EX)
          fv = 0.5*(c*(Vp+Vr) - ac*EY)
          fw = 0.5*(c*(Wp+Wr) - ac*EZ)

          msk = real(ibits(bv(0,j,k), State, 1))
          
          wv(i,j,k,1) = wv(i,j,k,1) + ( fu*dh1 - EX*dh2 ) * msk
          wv(i,j,k,2) = wv(i,j,k,2) + ( fv*dh1 - EY*dh2 ) * msk
          wv(i,j,k,3) = wv(i,j,k,3) + ( fw*dh1 - EZ*dh2 ) * msk
          m = m + 1.0
        endif
      end do
      end do
!$OMP END DO

      
    case (X_plus)

      i = ix
!$OMP DO SCHEDULE(static)
      do k=1,kx
      do j=1,jx
        if ( ibits(bv(i,j,k), bc_face_E, bitw_5) == obc_mask ) then
          Up = v0(i,j,k,1)
          Vp = v0(i,j,k,2)
          Wp = v0(i,j,k,3)
          
          Ur = u_bc
          Vr = v_bc
          Wr = w_bc
          c  = u_bc
          ac = abs(c)

          EX = Ur - Up
          EY = Vr - Vp
          EZ = Wr - Wp

          fu = 0.5*(c*(Ur+Up) - ac*EX)
          fv = 0.5*(c*(Vr+Vp) - ac*EY)
          fw = 0.5*(c*(Wr+Wp) - ac*EZ)

          msk = real(ibits(bv(ix+1,j,k), State, 1))

          wv(i,j,k,1) = wv(i,j,k,1) + ( -fu*dh1 + EX*dh2 ) * msk
          wv(i,j,k,2) = wv(i,j,k,2) + ( -fv*dh1 + EY*dh2 ) * msk
          wv(i,j,k,3) = wv(i,j,k,3) + ( -fw*dh1 + EZ*dh2 ) * msk
          m = m + 1.0
        endif
      end do
      end do
!$OMP END DO

      
    case (Y_minus)

      j = 1
!$OMP DO SCHEDULE(static)
      do k=1,kx
      do i=1,ix
        if ( ibits(bv(i,j,k), bc_face_S, bitw_5) == obc_mask ) then
          Up = v0(i,j,k,1)
          Vp = v0(i,j,k,2)
          Wp = v0(i,j,k,3)
          
          Ur = u_bc
          Vr = v_bc
          Wr = w_bc
          c  = v_bc
          ac = abs(c)

          EX = Up - Ur
          EY = Vp - Vr
          EZ = Wp - Wr

          fu = 0.5*(c*(Up+Ur) - ac*EX)
          fv = 0.5*(c*(Vp+Vr) - ac*EY)
          fw = 0.5*(c*(Wp+Wr) - ac*EZ)

          msk = real(ibits(bv(i,0,k), State, 1))
          
          wv(i,j,k,1) = wv(i,j,k,1) + ( fu*dh1 - EX*dh2 ) * msk
          wv(i,j,k,2) = wv(i,j,k,2) + ( fv*dh1 - EY*dh2 ) * msk
          wv(i,j,k,3) = wv(i,j,k,3) + ( fw*dh1 - EZ*dh2 ) * msk
          m = m + 1.0
        endif
      end do
      end do
!$OMP END DO

      
    case (Y_plus)

      j = jx
!$OMP DO SCHEDULE(static)
      do k=1,kx
      do i=1,ix
        if ( ibits(bv(i,j,k), bc_face_N, bitw_5) == obc_mask ) then
          Up = v0(i,j,k,1)
          Vp = v0(i,j,k,2)
          Wp = v0(i,j,k,3)
          
          Ur = u_bc
          Vr = v_bc
          Wr = w_bc
          c  = v_bc
          ac = abs(c)

          EX = Ur - Up
          EY = Vr - Vp
          EZ = Wr - Wp

          fu = 0.5*(c*(Ur+Up) - ac*EX)
          fv = 0.5*(c*(Vr+Vp) - ac*EY)
          fw = 0.5*(c*(Wr+Wp) - ac*EZ)

          msk = real(ibits(bv(i,jx+1,k), State, 1))
          
          wv(i,j,k,1) = wv(i,j,k,1) + ( -fu*dh1 + EX*dh2 ) * msk
          wv(i,j,k,2) = wv(i,j,k,2) + ( -fv*dh1 + EY*dh2 ) * msk
          wv(i,j,k,3) = wv(i,j,k,3) + ( -fw*dh1 + EZ*dh2 ) * msk
          m = m + 1.0
        endif
      end do
      end do
!$OMP END DO

      
    case (Z_minus)

      k = 1
!$OMP DO SCHEDULE(static)
      do j=1,jx
      do i=1,ix
        if ( ibits(bv(i,j,k), bc_face_B, bitw_5) == obc_mask ) then
          Up = v0(i,j,k,1)
          Vp = v0(i,j,k,2)
          Wp = v0(i,j,k,3)
          
          Ur = u_bc
          Vr = v_bc
          Wr = w_bc
          c  = w_bc
          ac = abs(c)

          EX = Up - Ur
          EY = Vp - Vr
          EZ = Wp - Wr

          fu = 0.5*(c*(Up+Ur) - ac*EX)
          fv = 0.5*(c*(Vp+Vr) - ac*EY)
          fw = 0.5*(c*(Wp+Wr) - ac*EZ)

          msk = real(ibits(bv(i,j,0), State, 1))
          
          wv(i,j,k,1) = wv(i,j,k,1) + ( fu*dh1 - EX*dh2 ) * msk
          wv(i,j,k,2) = wv(i,j,k,2) + ( fv*dh1 - EY*dh2 ) * msk
          wv(i,j,k,3) = wv(i,j,k,3) + ( fw*dh1 - EZ*dh2 ) * msk
          m = m + 1.0
        endif
      end do
      end do
!$OMP END DO

      
    case (Z_plus)

      k = kx
!$OMP DO SCHEDULE(static)
      do j=1,jx
      do i=1,ix
        if ( ibits(bv(i,j,k), bc_face_T, bitw_5) == obc_mask ) then
          Up = v0(i,j,k,1)
          Vp = v0(i,j,k,2)
          Wp = v0(i,j,k,3)
          
          Ur = u_bc
          Vr = v_bc
          Wr = w_bc
          c  = w_bc
          ac = abs(c)

          EX = Ur - Up
          EY = Vr - Vp
          EZ = Wr - Wp

          fu = 0.5*(c*(Ur+Up) - ac*EX)
          fv = 0.5*(c*(Vr+Vp) - ac*EY)
          fw = 0.5*(c*(Wr+Wp) - ac*EZ)

          msk = real(ibits(bv(i,j,kx+1), State, 1))

          wv(i,j,k,1) = wv(i,j,k,1) + ( -fu*dh1 + EX*dh2 ) * msk
          wv(i,j,k,2) = wv(i,j,k,2) + ( -fv*dh1 + EY*dh2 ) * msk
          wv(i,j,k,3) = wv(i,j,k,3) + ( -fw*dh1 + EZ*dh2 ) * msk
          m = m + 1.0
        endif
      end do
      end do
!$OMP END DO
      
    case default
    end select FACES
!$OMP END PARALLEL

    flop = flop + dble(m)*28.0d0
    
    return
    end subroutine vobc_pv_specv
    
    
!> ********************************************************************
!! @brief 外部速度境界条件による対流項と粘性項の流束の修正
!! @param [out] wv   疑似ベクトルの空間項の評価値
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル長
!! @param [in]  dh   格子幅
!! @param [in]  rei  Reynolds数の逆数
!! @param [in]  v0   速度ベクトル（n-step）
!! @param [in]  vec  指定する速度ベクトル
!! @param [in]  face 外部境界処理のときの面番号
!! @param [out] flop 浮動小数点演算数
!! @note 対流流束はゼロ，壁面法線方向の1階微分もゼロ
!<
    subroutine vobc_pv_wall (wv, sz, g, dh, rei, v0, vec, face, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, g, face
    integer                                                   ::  ix, jx, kx
    integer, dimension(3)                                     ::  sz
    double precision                                          ::  flop, rix, rjx, rkx
    real                                                      ::  u_bc, v_bc, w_bc
    real                                                      ::  dh, dh1, dh2, rei
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v0, wv
    real, dimension(3)                                        ::  vec
    
    ix = sz(1)
    jx = sz(2)
    kx = sz(3)

    dh1= 1.0/dh
    dh2= rei*dh1*dh1*2.0
    
    u_bc = vec(1)
    v_bc = vec(2)
    w_bc = vec(3)


!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, u_bc, v_bc, w_bc, dh2, face) &
!$OMP PRIVATE(i, j, k)

    FACES : select case (face)
    
    case (X_minus)

      i = 1
!$OMP DO SCHEDULE(static)
      do k=1,kx
      do j=1,jx
      ! wv(i,j,k,1) = wv(i,j,k,1) + (u_bc - v0(i,j,k,1)) * dh2 ! du/dx = 0
        wv(i,j,k,2) = wv(i,j,k,2) + (v_bc - v0(i,j,k,2)) * dh2 ! dv/dx
        wv(i,j,k,3) = wv(i,j,k,3) + (w_bc - v0(i,j,k,3)) * dh2 ! dw/dx
      end do
      end do
!$OMP END DO

      
    case (X_plus)

      i = ix
!$OMP DO SCHEDULE(static)
      do k=1,kx
      do j=1,jx
      ! wv(i,j,k,1) = wv(i,j,k,1) + (u_bc - v0(i,j,k,1)) * dh2 ! du/dx = 0
        wv(i,j,k,2) = wv(i,j,k,2) + (v_bc - v0(i,j,k,2)) * dh2 ! dv/dx
        wv(i,j,k,3) = wv(i,j,k,3) + (w_bc - v0(i,j,k,3)) * dh2 ! dw/dx
      end do
      end do
!$OMP END DO

      
    case (Y_minus)

      j = 1
!$OMP DO SCHEDULE(static)
      do k=1,kx
      do i=1,ix
        wv(i,j,k,1) = wv(i,j,k,1) + (u_bc - v0(i,j,k,1)) * dh2 ! du/dy
      ! wv(i,j,k,2) = wv(i,j,k,2) + (v_bc - v0(i,j,k,2)) * dh2 ! dv/dy = 0
        wv(i,j,k,3) = wv(i,j,k,3) + (w_bc - v0(i,j,k,3)) * dh2 ! dw/dy
      end do
      end do
!$OMP END DO

      
    case (Y_plus)

      j = jx
!$OMP DO SCHEDULE(static)
      do k=1,kx
      do i=1,ix
        wv(i,j,k,1) = wv(i,j,k,1) + (u_bc - v0(i,j,k,1)) * dh2 ! du/dy
      ! wv(i,j,k,2) = wv(i,j,k,2) + (v_bc - v0(i,j,k,2)) * dh2 ! dv/dy = 0
        wv(i,j,k,3) = wv(i,j,k,3) + (w_bc - v0(i,j,k,3)) * dh2 ! dw/dy
      end do
      end do
!$OMP END DO

      
    case (Z_minus)

      k = 1
!$OMP DO SCHEDULE(static)
      do j=1,jx
      do i=1,ix
        wv(i,j,k,1) = wv(i,j,k,1) + (u_bc - v0(i,j,k,1)) * dh2 ! du/dz
        wv(i,j,k,2) = wv(i,j,k,2) + (v_bc - v0(i,j,k,2)) * dh2 ! dv/dz
      ! wv(i,j,k,3) = wv(i,j,k,3) + (w_bc - v0(i,j,k,3)) * dh2 ! dw/dz = 0
      end do
      end do
!$OMP END DO

      
    case (Z_plus)

      k = kx
!$OMP DO SCHEDULE(static)
      do j=1,jx
      do i=1,ix
        wv(i,j,k,1) = wv(i,j,k,1) + (u_bc - v0(i,j,k,1)) * dh2 ! du/dz
        wv(i,j,k,2) = wv(i,j,k,2) + (v_bc - v0(i,j,k,2)) * dh2 ! dv/dz
      ! wv(i,j,k,3) = wv(i,j,k,3) + (w_bc - v0(i,j,k,3)) * dh2 ! dw/dz = 0
      end do
      end do
!$OMP END DO

    case default
    end select FACES
!$OMP END PARALLEL


    rix = dble(jx)*dble(kx)
    rjx = dble(ix)*dble(kx)
    rkx = dble(ix)*dble(jx)

    flop = flop + 14.0d0


    FACES2 : select case (face)

    case (X_minus)
      flop = flop + rix*9.0d0

    case (X_plus)
      flop = flop + rix*9.0d0

    case (Y_minus)
      flop = flop + rjx*9.0d0

    case (Y_plus)
      flop = flop + rjx*9.0d0

    case (Z_minus)
      flop = flop + rkx*9.0d0

    case (Z_plus)
      flop = flop + rkx*9.0d0

    case default
    end select FACES2

    return
    end subroutine vobc_pv_wall


!> ********************************************************************
!! @brief ガイドセルの速度指定境界条件を設定するために必要な参照値をセットする
!! @param [out] v    セルセンタ速度
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル長
!! @param [in]  bv   BCindex V
!! @param [in]  face 外部境界の面番号
!! @param [in]  vec  指定する速度ベクトル
!! @note 流束型の境界条件を用いるので，内点の計算に使う参照点に値があればよい（1層）
!<
    subroutine vobc_drchlt (v, sz, g, bv, face, vec)
    implicit none
    include 'ffv_f_params.h'
    integer                                                     ::  i, j, k, g, face, ix, jx, kx
    integer, dimension(3)                                       ::  sz
    real                                                        ::  u_bc, v_bc, w_bc
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)   ::  v
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bv
    real, dimension(3)                                          ::  vec

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    
    ! u_bcは境界速度
    u_bc = vec(1)
    v_bc = vec(2)
    w_bc = vec(3)
    
!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, u_bc, v_bc, w_bc, face) &
!$OMP PRIVATE(i, j, k)

    FACES : select case (face)
    
    case (X_minus)

      i = 1
!$OMP DO SCHEDULE(static)
      do k=1,kx
      do j=1,jx
        if ( ibits(bv(i,j,k), bc_face_W, bitw_5) == obc_mask ) then
          v(i-1, j, k, 1) = u_bc
          v(i-1, j, k, 2) = v_bc
          v(i-1, j, k, 3) = w_bc
        endif
      end do
      end do
!$OMP END DO
      
      
    case (X_plus)

      i = ix
!$OMP DO SCHEDULE(static)
      do k=1,kx
      do j=1,jx
        if ( ibits(bv(i,j,k), bc_face_E, bitw_5) == obc_mask ) then
          v(i+1, j, k, 1) = u_bc
          v(i+1, j, k, 2) = v_bc
          v(i+1, j, k, 3) = w_bc
        endif
      end do
      end do
!$OMP END DO
      
      
    case (Y_minus)

      j = 1
!$OMP DO SCHEDULE(static)
      do k=1,kx
      do i=1,ix
        if ( ibits(bv(i,j,k), bc_face_S, bitw_5) == obc_mask ) then
          v(i, j-1, k, 1) = u_bc
          v(i, j-1, k, 2) = v_bc
          v(i, j-1, k, 3) = w_bc
        endif
      end do
      end do
!$OMP END DO
      
      
    case (Y_plus)

      j = jx
!$OMP DO SCHEDULE(static)
      do k=1,kx
      do i=1,ix
        if ( ibits(bv(i,j,k), bc_face_N, bitw_5) == obc_mask ) then
          v(i, j+1, k, 1) = u_bc
          v(i, j+1, k, 2) = v_bc
          v(i, j+1, k, 3) = w_bc
        endif
      end do
      end do
!$OMP END DO
      
      
    case (Z_minus)

      k = 1
!$OMP DO SCHEDULE(static)
      do j=1,jx
      do i=1,ix
        if ( ibits(bv(i,j,k), bc_face_B, bitw_5) == obc_mask ) then
          v(i, j, k-1, 1) = u_bc
          v(i, j, k-1, 2) = v_bc
          v(i, j, k-1, 3) = w_bc
        endif
      end do
      end do
!$OMP END DO
      
      
    case (Z_plus)

      k = kx
!$OMP DO SCHEDULE(static)
      do j=1,jx
      do i=1,ix
        if ( ibits(bv(i,j,k), bc_face_T, bitw_5) == obc_mask ) then
          v(i, j, k+1, 1) = u_bc
          v(i, j, k+1, 2) = v_bc
          v(i, j, k+1, 3) = w_bc
        endif
      end do
      end do
!$OMP END DO
      
    case default
    end select FACES
!$OMP END PARALLEL
    
    return
    end subroutine vobc_drchlt

!> ********************************************************************
!! @brief 遠方境界の近似
!! @param [out] v    速度ベクトル
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル長
!! @param [in]  face 外部境界面の番号
!<
    subroutine vobc_neumann (v, sz, g, face)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, face, g
    integer, dimension(3)                                     ::  sz
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, g, face)

    FACES : select case (face)
    case (X_minus)

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do j=1,jx
      do i=1-g,0
        v(i, j, k, 1) = v(1-i, j, k, 1)
        v(i, j, k, 2) = v(1-i, j, k, 2)
        v(i, j, k, 3) = v(1-i, j, k, 3)
      end do
      end do
      end do
!$OMP END DO
      

    case (X_plus)

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do j=1,jx
      do i=ix+1,ix+g
        v(i, j, k, 1) = v(2*ix+1-i, j, k, 1)
        v(i, j, k, 2) = v(2*ix+1-i, j, k, 2)
        v(i, j, k, 3) = v(2*ix+1-i, j, k, 3)
      end do
      end do
      end do
!$OMP END DO
      

    case (Y_minus)

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do j=1-g,0
      do i=1,ix
        v(i, j, k, 1) = v(i, 1-j, k, 1)
        v(i, j, k, 2) = v(i, 1-j, k, 2)
        v(i, j, k, 3) = v(i, 1-j, k, 3)
      end do
      end do
      end do
!$OMP END DO
      

    case (Y_plus)

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do j=jx+1,jx+g
      do i=1,ix
        v(i, j, k, 1) = v(i, 2*jx+1-j, k, 1)
        v(i, j, k, 2) = v(i, 2*jx+1-j, k, 2)
        v(i, j, k, 3) = v(i, 2*jx+1-j, k, 3)
      end do
      end do
      end do
!$OMP END DO
      

    case (Z_minus)

!$OMP DO SCHEDULE(static)
      do j=1,jx
      do k=1-g,0
      do i=1,ix
        v(i, j, k, 1) = v(i, j, 1-k, 1)
        v(i, j, k, 2) = v(i, j, 1-k, 2)
        v(i, j, k, 3) = v(i, j, 1-k, 3)
      end do
      end do
      end do
!$OMP END DO
      

    case (Z_plus)

!$OMP DO SCHEDULE(static)
      do j=1,jx
      do k=kx+1,kx+g
      do i=1,ix
        v(i, j, k, 1) = v(i, j, 2*kx+1-k, 1)
        v(i, j, k, 2) = v(i, j, 2*kx+1-k, 2)
        v(i, j, k, 3) = v(i, j, 2*kx+1-k, 3)
      end do
      end do
      end do
!$OMP END DO


    case default
    end select FACES
!$OMP END PARALLEL

    return 
    end subroutine vobc_neumann

!> ********************************************************************
!! @brief 外部流出境界で，次ステップの流出速度を対流流出条件で予測し，ガイドセルに参照値として代入する
!! @param [out]    v    速度 u^*
!! @param [in]     sz   配列長
!! @param [in]     g    ガイドセル長
!! @param [in]     cf   uc*dt/dh
!! @param [in]     face 外部境界の面番号
!! @param [in]     v0   セルセンタ速度 u^n
!! @param [in]     bv   BCindex V
!! @param [in,out] flop 浮動小数点演算数
!<
    subroutine vobc_outflow (v, sz, g, cf, face, v0, bv, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, g, face, ix, jx, kx
    integer, dimension(3)                                     ::  sz
    double precision                                          ::  flop
    real                                                      ::  Ue, Uw, Un, Us, Ut, Ub
    real                                                      ::  Ve, Vw, Vn, Vs, Vt, Vb
    real                                                      ::  We, Ww, Wn, Ws, Wt, Wb
    real                                                      ::  cf, c, m
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v, v0
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bv

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)

    c = cf
    m = 0.0
    
!$OMP PARALLEL REDUCTION(+:m) &
!$OMP FIRSTPRIVATE(ix, jx, kx, face, c, cf) &
!$OMP PRIVATE(Ue, Uw, Un, Us, Ut, Ub) &
!$OMP PRIVATE(Ve, Vw, Vn, Vs, Vt, Vb) &
!$OMP PRIVATE(We, Ww, Wn, Ws, Wt, Wb)
    
    FACES : select case (face)
    case (X_minus)

      if ( cf>0.0 ) c=0.0
      
!$OMP DO SCHEDULE(static)
      do k=1,kx
      do j=1,jx
        if ( ibits(bv(1, j, k), bc_face_W, bitw_5) == obc_mask ) then
          Uw = v0(0, j, k, 1)
          Ue = v0(1, j, k, 1)
          Vw = v0(0, j, k, 2)
          Ve = v0(1, j, k, 2)
          Ww = v0(0, j, k, 3)
          We = v0(1, j, k, 3)

          v(0, j ,k, 1) = Uw - c*(Ue-Uw)
          v(0, j ,k, 2) = Vw - c*(Ve-Vw)
          v(0, j ,k, 3) = Ww - c*(We-Ww)
          m = m + 1.0
        endif
      end do
      end do
!$OMP END DO
      
      
    case (X_plus)

      if ( cf<0.0 ) c=0.0
      
!$OMP DO SCHEDULE(static)
      do k=1,kx
      do j=1,jx
        if ( ibits(bv(ix, j, k), bc_face_E, bitw_5) == obc_mask ) then
          Uw = v0(ix  , j, k, 1)
          Ue = v0(ix+1, j, k, 1)
          Vw = v0(ix  , j, k, 2)
          Ve = v0(ix+1, j, k, 2)
          Ww = v0(ix  , j, k, 3)
          We = v0(ix+1, j, k, 3)

          v(ix+1, j ,k, 1) = Ue - c*(Ue-Uw)
          v(ix+1, j ,k, 2) = Ve - c*(Ve-Vw)
          v(ix+1, j ,k, 3) = We - c*(We-Ww)
          m = m + 1.0
        endif
      end do
      end do
!$OMP END DO
      
      
    case (Y_minus)

      if ( cf>0.0 ) c=0.0

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do i=1,ix
        if ( ibits(bv(i, 1, k), bc_face_S, bitw_5) == obc_mask ) then
          Us = v0(i, 0, k, 1)
          Un = v0(i, 1, k, 1)
          Vs = v0(i, 0, k, 2)
          Vn = v0(i, 1, k, 2)
          Ws = v0(i, 0, k, 3)
          Wn = v0(i, 1, k, 3)

          v(i, 0, k, 1) = Us - c*(Un-Us)
          v(i, 0, k, 2) = Vs - c*(Vn-Vs)
          v(i, 0, k, 3) = Ws - c*(Wn-Ws)
          m = m + 1.0
        endif
      end do
      end do
!$OMP END DO
      
      
    case (Y_plus)

      if ( cf<0.0 ) c=0.0
      
!$OMP DO SCHEDULE(static)
      do k=1,kx
      do i=1,ix
        if ( ibits(bv(i, jx, k), bc_face_N, bitw_5) == obc_mask ) then
          Us = v0(i, jx  , k, 1)
          Un = v0(i, jx+1, k, 1)
          Vs = v0(i, jx  , k, 2)
          Vn = v0(i, jx+1, k, 2)
          Ws = v0(i, jx  , k, 3)
          Wn = v0(i, jx+1, k, 3)

          v(i, jx+1, k, 1) = Un - c*(Un-Us)
          v(i, jx+1, k, 2) = Vn - c*(Vn-Vs)
          v(i, jx+1, k, 3) = Wn - c*(Wn-Ws)
          m = m + 1.0
        endif
      end do
      end do
!$OMP END DO
      
      
    case (Z_minus)

      if ( cf>0.0 ) c=0.0

!$OMP DO SCHEDULE(static)
      do j=1,jx
      do i=1,ix
        if ( ibits(bv(i, j, 1), bc_face_B, bitw_5) == obc_mask ) then
          Ub = v0(i, j, 0, 1)
          Ut = v0(i, j, 1, 1)
          Vb = v0(i, j, 0, 2)
          Vt = v0(i, j, 1, 2)
          Wb = v0(i, j, 0, 3)
          Wt = v0(i, j, 1, 3)

          v(i, j, 0, 1) = Ub - c*(Ut-Ub)
          v(i, j, 0, 2) = Vb - c*(Vt-Vb)
          v(i, j, 0, 3) = Wb - c*(Wt-Wb)
          m = m + 1.0
        endif
      end do
      end do
!$OMP END DO
      
      
    case (Z_plus)

      if ( cf<0.0 ) c=0.0
      
!$OMP DO SCHEDULE(static)
      do j=1,jx
      do i=1,ix
        if ( ibits(bv(i, j, kx), bc_face_T, bitw_5) == obc_mask ) then
          Ub = v0(i, j, kx  , 1)
          Ut = v0(i, j, kx+1, 1)
          Vb = v0(i, j, kx  , 2)
          Vt = v0(i, j, kx+1, 2)
          Wb = v0(i, j, kx  , 3)
          Wt = v0(i, j, kx+1, 3)

          v(i, j, kx+1, 1) = Ut - c*(Ut-Ub)
          v(i, j, kx+1, 2) = Vt - c*(Vt-Vb)
          v(i, j, kx+1, 3) = Wt - c*(Wt-Wb)
          m = m + 1.0
        endif
      end do
      end do
!$OMP END DO
      
    case default
    end select FACES
!$OMP END PARALLEL

    flop = flop + dble(m)*9.0d0

    return
    end subroutine vobc_outflow

!> ********************************************************************
!! @brief 速度の外部境界：　トラクションフリー
!! @param v 速度ベクトル
!! @param sz 配列長
!! @param g ガイドセル長
!! @param face 外部境界面の番号
!! @param flop 浮動小数演算数
!! @note トラクションフリー面は全て流体のこと
!<
    subroutine vobc_tfree (v, sz, g, face, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, face, g, ii, jj, kk
    integer, dimension(3)                                     ::  sz
    double precision                                          ::  flop, rix, rjx, rkx
    real                                                      ::  v1, v2, v3, v4
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    rix = dble(jx)*dble(kx)
    rjx = dble(ix)*dble(kx)
    rkx = dble(ix)*dble(jx)

!$OMP PARALLEL REDUCTION(+:flop) &
!$OMP FIRSTPRIVATE(ix, jx, kx, rix, rjx, rkx, g, face) &
!$OMP PRIVATE(i, j, k, ii, jj, kk, v1, v2, v3, v4)

    FACES : select case (face)
    case (X_minus)
      i = 1

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do j=1,jx

        v1 = 0.5 * (v(i-1, j+1, k  , 1) + v(i, j+1, k  , 1))
        v2 = 0.5 * (v(i-1, j-1, k  , 1) + v(i, j-1, k  , 1))
        v3 = 0.5 * (v(i-1, j  , k+1, 1) + v(i, j  , k+1, 1))
        v4 = 0.5 * (v(i-1, j  , k-1, 1) + v(i, j  , k-1, 1))

        v(i-1, j, k, 1) = v(i, j, k, 1)
        v(i-1, j, k, 2) = v(i, j, k, 2) + 0.5 * (v1 - v2)
        v(i-1, j, k, 3) = v(i, j, k, 3) + 0.5 * (v3 - v4)

      end do
      end do
!$OMP END DO
      
      flop = flop + rix*12.0d0
      

    case (X_plus)
      i = ix

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do j=1,jx

        v1 = 0.5 * (v(i+1, j+1, k  , 1) + v(i, j+1, k  , 1))
        v2 = 0.5 * (v(i+1, j-1, k  , 1) + v(i, j-1, k  , 1))
        v3 = 0.5 * (v(i+1, j  , k+1, 1) + v(i, j  , k+1, 1))
        v4 = 0.5 * (v(i+1, j  , k-1, 1) + v(i, j  , k-1, 1))
        
        v(i+1, j, k, 1) = v(i, j, k, 1)
        v(i+1, j, k, 2) = v(i, j, k, 2) - 0.5 * (v1 - v2)
        v(i+1, j, k, 3) = v(i, j, k, 3) - 0.5 * (v3 - v4)

      end do
      end do
!$OMP END DO
      
      flop = flop + rix*12.0d0
      

    case (Y_minus)
      j = 1

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do i=1,ix

        v1 = 0.5 * (v(i+1, j-1, k  , 2) + v(i+1, j, k  , 2))
        v2 = 0.5 * (v(i-1, j-1, k  , 2) + v(i-1, j, k  , 2))
        v3 = 0.5 * (v(i  , j-1, k+1, 2) + v(i  , j, k+1, 2))
        v4 = 0.5 * (v(i  , j-1, k-1, 2) + v(i  , j, k-1, 2))
                
        v(i, j-1, k, 1) = v(i, j, k, 1) + 0.5 * (v1 - v2)
        v(i, j-1, k, 2) = v(i, j, k, 2)
        v(i, j-1, k, 3) = v(i, j, k, 3) + 0.5 * (v3 - v4)

      end do
      end do
!$OMP END DO
      
      flop = flop + rjx*12.0d0
      

    case (Y_plus)
      j = jx

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do i=1,ix

        v1 = 0.5 * (v(i+1, j+1, k  , 2) + v(i+1, j, k  , 2))
        v2 = 0.5 * (v(i-1, j+1, k  , 2) + v(i-1, j, k  , 2))
        v3 = 0.5 * (v(i  , j+1, k+1, 2) + v(i  , j, k+1, 2))
        v4 = 0.5 * (v(i  , j+1, k-1, 2) + v(i  , j, k-1, 2))
                
        v(i, j+1, k, 1) = v(i, j, k, 1) - 0.5 * (v1 - v2)
        v(i, j+1, k, 2) = v(i, j, k, 2)
        v(i, j+1, k, 3) = v(i, j, k, 3) - 0.5 * (v3 - v4)

      end do
      end do
!$OMP END DO

      flop = flop + rjx*12.0d0
      

    case (Z_minus)
      k=1
      
!$OMP DO SCHEDULE(static)
      do j=1,jx
      do i=1,ix

        v1 = 0.5 * (v(i+1, j  , k-1, 3) + v(i+1, j  , k, 3))
        v2 = 0.5 * (v(i-1, j  , k-1, 3) + v(i-1, j  , k, 3))
        v3 = 0.5 * (v(i  , j+1, k-1, 3) + v(i  , j+1, k, 3))
        v4 = 0.5 * (v(i  , j-1, k-1, 3) + v(i  , j-1, k, 3))
                
        v(i, j, k-1, 1) = v(i, j, k, 1) + 0.5 * (v1 - v2)
        v(i, j, k-1, 2) = v(i, j, k, 2) + 0.5 * (v3 - v4)
        v(i, j, k-1, 3) = v(i, j, k, 3)

      end do
      end do
!$OMP END DO

      flop = flop + rkx*12.0d0
      

    case (Z_plus)
      k = kx

!$OMP DO SCHEDULE(static)
      do j=1,jx
      do i=1,ix

        v1= 0.5 * (v(i+1, j  , k+1, 3) + v(i+1, j  , k, 3))
        v2= 0.5 * (v(i-1, j  , k+1, 3) + v(i-1, j  , k, 3))
        v3= 0.5 * (v(i  , j+1, k+1, 3) + v(i  , j+1, k, 3))
        v4= 0.5 * (v(i  , j-1, k+1, 3) + v(i  , j-1, k, 3))
                
        v(i, j, k+1, 1) = v(i, j, k, 1) - 0.5 * (v1 - v2)
        v(i, j, k+1, 2) = v(i, j, k, 2) - 0.5 * (v3 - v4)
        v(i, j, k+1, 3) = v(i, j, k, 3)

      end do
      end do
!$OMP END DO

      flop = flop + rkx*12.0d0
      

    case default
    end select FACES

!$OMP END PARALLEL

    return 
    end subroutine vobc_tfree

!> ********************************************************************
!! @brief 疑似速度から次ステップ速度へ参照する速度をコピーする
!! @param [out] v    速度ベクトル（セルセンタ）
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル長
!! @param [in]  vc   セルセンタ疑似速度 u^*
!! @param [in]  face 面番号
!<
    subroutine vobc_update (v, sz, g, vc, face)
    implicit none
    include 'ffv_f_params.h'
    integer                                                     ::  i, j, k, g, ix, jx, kx, face
    integer, dimension(3)                                       ::  sz
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)   ::  v, vc

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    
!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, face)
    
    FACES : select case (face)
    case (X_minus)

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do j=1,jx
        v(0, j, k, 1) = vc(0, j, k, 1)
        v(0, j, k, 2) = vc(0, j, k, 2)
        v(0, j, k, 3) = vc(0, j, k, 3)
      end do
      end do
!$OMP END DO
      
    case (X_plus)

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do j=1,jx
        v(ix+1, j, k, 1) = vc(ix+1, j, k, 1)
        v(ix+1, j, k, 2) = vc(ix+1, j, k, 2)
        v(ix+1, j, k, 3) = vc(ix+1, j, k, 3)
      end do
      end do
!$OMP END DO
      
    case (Y_minus)

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do i=1,ix
        v(i, 0, k, 1) = vc(i, 0, k, 1)
        v(i, 0, k, 2) = vc(i, 0, k, 2)
        v(i, 0, k, 3) = vc(i, 0, k, 3)
      end do
      end do
!$OMP END DO
      
    case (Y_plus)

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do i=1,ix
        v(i, jx+1, k, 1) = vc(i, jx+1, k, 1)
        v(i, jx+1, k, 2) = vc(i, jx+1, k, 2)
        v(i, jx+1, k, 3) = vc(i, jx+1, k, 3)
      end do
      end do
!$OMP END DO
      
    case (Z_minus)

!$OMP DO SCHEDULE(static)
      do j=1,jx
      do i=1,ix
        v(i, j, 0, 1) = vc(i, j, 0, 1)
        v(i, j, 0, 2) = vc(i, j, 0, 2)
        v(i, j, 0, 3) = vc(i, j, 0, 3)
      end do
      end do
!$OMP END DO
      
    case (Z_plus)

!$OMP DO SCHEDULE(static)
      do j=1,jx
      do i=1,ix
        v(i, j, kx+1, 1) = vc(i, j, kx+1, 1)
        v(i, j, kx+1, 2) = vc(i, j, kx+1, 2)
        v(i, j, kx+1, 3) = vc(i, j, kx+1, 3)
      end do
      end do
!$OMP END DO
      
    case default
    end select FACES

!$OMP END PARALLEL
    
    return
    end subroutine vobc_update


!> ********************************************************************
!! @brief 外部指定境界条件による速度の発散の修正
!! @param [in,out] div   速度の発散
!! @param [in]     sz    配列長
!! @param [in]     g     ガイドセル長
!! @param [in]     face  面番号
!! @param [in]     bv    BCindex V
!! @param [in]     vec   指定する速度ベクトル
!! @param [in,out] flop  flop count
!! @note 固体部分は対象外とするのでループ中に判定あり
!! 部分的な境界条件の実装のため、ガイドセル部のマスク情報を利用
!<
    subroutine vobc_div_drchlt (div, sz, g, face, bv, vec, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, g, ix, jx, kx, face, bvx
    integer, dimension(3)                                     ::  sz
    double precision                                          ::  flop, rix, rjx, rkx
    real                                                      ::  u_bc, v_bc, w_bc
    real, dimension(3)                                        ::  vec
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  div
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bv

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    
    u_bc = vec(1)
    v_bc = vec(2)
    w_bc = vec(3)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, u_bc, v_bc, w_bc, face)

    FACES : select case (face)
    case (X_minus)

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do j=1,jx
        bvx = bv(1, j, k)
        if ( ibits(bvx, bc_face_W, bitw_5) == obc_mask ) then
          div(1, j, k) = div(1, j, k) - u_bc * real(ibits(bvx, State, 1)) * real(ibits(bv(0, j, k), State, 1))
        endif
      end do
      end do
!$OMP END DO
      
      
    case (X_plus)

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do j=1,jx
        bvx = bv(ix, j, k)
        if ( ibits(bvx, bc_face_E, bitw_5) == obc_mask ) then
          div(ix, j, k) = div(ix, j, k) + u_bc * real(ibits(bvx, State, 1)) * real(ibits(bv(ix+1, j, k), State, 1))
        endif
      end do
      end do
!$OMP END DO
      
      
    case (Y_minus)

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do i=1,ix
        bvx = bv(i, 1, k)
        if ( ibits(bvx, bc_face_S, bitw_5) == obc_mask ) then
          div(i, 1, k) = div(i, 1, k) - v_bc * real(ibits(bvx, State, 1)) * real(ibits(bv(i, 0, k), State, 1))
        endif
      end do
      end do
!$OMP END DO
      
      
    case (Y_plus)

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do i=1,ix
        bvx = bv(i, jx, k)
        if ( ibits(bvx, bc_face_N, bitw_5) == obc_mask ) then
          div(i, jx, k) = div(i, jx, k) + v_bc * real(ibits(bvx, State, 1)) * real(ibits(bv(i, jx+1, k), State, 1))
        endif
      end do
      end do
!$OMP END DO
    
    
    case (Z_minus)

!$OMP DO SCHEDULE(static)
      do j=1,jx
      do i=1,ix
        bvx = bv(i, j, 1)
        if ( ibits(bvx, bc_face_B, bitw_5) == obc_mask ) then
          div(i, j, 1) = div(i, j, 1) - w_bc * real(ibits(bvx, State, 1)) * real(ibits(bv(i, j, 0), State, 1))
        endif
      end do
      end do
!$OMP END DO
      
      
    case (Z_plus)

!$OMP DO SCHEDULE(static)
      do j=1,jx
      do i=1,ix
        bvx = bv(i, j, kx)
        if ( ibits(bvx, bc_face_T, bitw_5) == obc_mask ) then
          div(i, j, kx) = div(i, j, kx) + w_bc * real(ibits(bvx, State, 1)) * real(ibits(bv(i, j, kx+1), State, 1))
        endif
      end do
      end do
!$OMP END DO
    
    case default
    end select FACES

!$OMP END PARALLEL


    rix = dble(jx)*dble(kx)
    rjx = dble(ix)*dble(kx)
    rkx = dble(ix)*dble(jx)

    flop = flop + 9.0d0

    FACES2 : select case (face)
    case (X_minus)
      flop = flop + rix*3.0d0 ! 2+ real*1

    case (X_plus)
      flop = flop + rix*3.0d0 ! 2+ real*1

    case (Y_minus)
      flop = flop + rjx*3.0d0 ! 2+ real*1

    case (Y_plus)
      flop = flop + rjx*3.0d0 ! 2+ real*1

    case (Z_minus)
      flop = flop + rkx*3.0d0 ! 2+ real*1

    case (Z_plus)
      flop = flop + rkx*3.0d0 ! 2+ real*1

    case default
    end select FACES2


    return
    end subroutine vobc_div_drchlt

!> ********************************************************************
!! @brief 外部流出境界条件による疑似速度ベクトルの発散の修正
!! @param [in,out] div   速度の発散
!! @param [in]     sz    配列長
!! @param [in]     g     ガイドセル長
!! @param [in]     face  面番号
!! @param [in]     cf    uc*dt/dh
!! @param [in]     bv    BCindex V
!! @param [in]     vf    セルフェイス速度ベクトル（n-step）
!! @param [out]    flop  flop count
!<
    subroutine vobc_div_pv_oflow (div, sz, g, face, cf, bv, vf, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, g, ix, jx, kx, face
    integer, dimension(3)                                     ::  sz
    double precision                                          ::  flop
    real                                                      ::  cf, c, m
    real                                                      ::  b_w, b_e, b_s, b_n, b_b, b_t, b_p
    real                                                      ::  w_e, w_w, w_n, w_s, w_t, w_b
    real                                                      ::  Ue, Uw, Vn, Vs, Wt, Wb
    real                                                      ::  Ue_t, Uw_t, Vn_t, Vs_t, Wt_t, Wb_t
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  div
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bv
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  vf

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)

    c = cf
    m =0.0

!$OMP PARALLEL REDUCTION(+:m) &
!$OMP FIRSTPRIVATE(ix, jx, kx, face, cf, c) &
!$OMP PRIVATE(Ue, Uw, Vn, Vs, Wt, Wb) &
!$OMP PRIVATE(Ue_t, Uw_t, Vn_t, Vs_t, Wt_t, Wb_t) &
!$OMP PRIVATE(w_e, w_w, w_n, w_s, w_t, w_b) &
!$OMP PRIVATE(b_w, b_e, b_s, b_n, b_b, b_t, b_p) &
!$OMP PRIVATE(i, j, k)


    FACES : select case (face)
    
    case (X_minus)

      if ( cf>0.0 ) c=0.0

      i = 1
!$OMP DO SCHEDULE(static)
      do k=1,kx
      do j=1,jx
        if ( ibits(bv(i,j,k), bc_face_W, bitw_5) == obc_mask ) then
          include 'd_o_o_p.h' ! 25 flops

          Uw_t = Uw - c*(Ue-Uw)
          div(i,j,k) = div(i,j,k) - Uw_t
          m = m + 1.0
        end if
      end do
      end do
!$OMP END DO

      
    case (X_plus)

      if ( cf<0.0 ) c=0.0

      i = ix
!$OMP DO SCHEDULE(static)
      do k=1,kx
      do j=1,jx
        if ( ibits(bv(i,j,k), bc_face_E, bitw_5) == obc_mask ) then
          include 'd_o_o_p.h'

          Ue_t = Ue - c*(Ue-Uw)
          div(i,j,k) = div(i,j,k) + Ue_t
          m = m + 1.0
        end if
      end do
      end do
!$OMP END DO


    case (Y_minus)

      if ( cf>0.0 ) c=0.0

      j = 1
!$OMP DO SCHEDULE(static)
      do k=1,kx
      do i=1,ix
        if ( ibits(bv(i,j,k), bc_face_S, bitw_5) == obc_mask ) then
          include 'd_o_o_p.h'

          Vs_t = Vs - c*(Vn-Vs)
          div(i,j,k) = div(i,j,k) - Vs_t
          m = m + 1.0
        end if
      end do
      end do
!$OMP END DO

      
    case (Y_plus)

      if ( cf<0.0 ) c=0.0

      j = jx
!$OMP DO SCHEDULE(static)
      do k=1,kx
      do i=1,ix
        if ( ibits(bv(i,j,k), bc_face_N, bitw_5) == obc_mask ) then
          include 'd_o_o_p.h'

          Vn_t = Vn - c*(Vn-Vs)
          div(i,j,k) = div(i,j,k) + Vn_t
          m = m + 1.0
        end if
      end do
      end do
!$OMP END DO

    
    case (Z_minus)

      if ( cf>0.0 ) c=0.0

      k = 1
!$OMP DO SCHEDULE(static)
      do j=1,jx
      do i=1,ix
        if ( ibits(bv(i,j,k), bc_face_B, bitw_5) == obc_mask ) then
          include 'd_o_o_p.h'

          Wb_t = Wb - c*(Wt-Wb)
          div(i,j,k) = div(i,j,k) - Wb_t
          m = m + 1.0
        end if
      end do
      end do
!$OMP END DO

      
    case (Z_plus)

      if ( cf<0.0 ) c=0.0

      k = kx
!$OMP DO SCHEDULE(static)
      do j=1,jx
      do i=1,ix
        if ( ibits(bv(i,j,k), bc_face_T, bitw_5) == obc_mask ) then
          include 'd_o_o_p.h'

          Wt_t = Wt - c*(Wt-Wb)
          div(i,j,k) = div(i,j,k) + Wt_t
          m = m + 1.0
        end if
      end do
      end do
!$OMP END DO

    case default
    end select FACES
!$OMP END PARALLEL

    flop = flop + dble(m)*29.0d0

    return
    end subroutine vobc_div_pv_oflow

!> ********************************************************************
!! @brief 外部流出境界条件による速度ベクトルの発散の修正
!! @param [in,out] div  \sum{u}
!! @param [in]     sz   配列長
!! @param [in]     g    ガイドセル長
!! @param [in]     face 面番号
!! @param [out]    aa   領域境界の積算値
!! @param [out]    vf   セルフェイス速度 n+1
!! @param [in]     bv    BCindex V
!! @param [out]    flop flop count 近似
!! @note 指定面でも固体部分は対象外とするのでループ中に判定あり
!!       div(u)=0から，内部流出境界のセルで計算されたdivが流出速度となる
!<
    subroutine vobc_div_oflow (div, sz, g, face, aa, vf, bv, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, g, ix, jx, kx, face
    integer, dimension(3)                                     ::  sz
    double precision                                          ::  flop
    real                                                      ::  dv, a1, a2, a3, m
    real, dimension(3)                                        ::  aa
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  div
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  vf
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bv

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)

    a1 = 0.0   ! sum
    a2 = 1.0e6 ! min
    a3 =-1.0e6 ! max

    m = 0.0
    
    
!$OMP PARALLEL &
!$OMP REDUCTION(+:a1) &
!$OMP REDUCTION(min:a2) &
!$OMP REDUCTION(max:a3) &
!$OMP REDUCTION(+:m) &
!$OMP FIRSTPRIVATE(ix, jx, kx, face) &
!$OMP PRIVATE(dv, i, j, k)

    FACES : select case (face)
    
    case (X_minus)

      i = 1
!$OMP DO SCHEDULE(static)
      do k=1,kx
      do j=1,jx
        if ( ibits(bv(i,j,k), bc_face_W, bitw_5) == obc_mask ) then
          dv = div(i,j,k)
          vf(i-1,j,k,1) = dv
          a1 = a1 + dv
          a2 = min(a2, dv)
          a3 = max(a3, dv)
          div(i,j,k) = 0.0 ! 対象セルは発散をゼロにする
          m = m + 1.0
        end if
      end do
      end do
!$OMP END DO
      
      
    case (X_plus)

      i = ix
!$OMP DO SCHEDULE(static)
      do k=1,kx
      do j=1,jx
        if ( ibits(bv(i,j,k), bc_face_E, bitw_5) == obc_mask ) then
          dv = -div(i,j,k)
          vf(i,j,k,1) = dv
          a1 = a1 + dv
          a2 = min(a2, dv)
          a3 = max(a3, dv)
          div(i,j,k) = 0.0
          m = m + 1.0
        end if
      end do
      end do
!$OMP END DO
      
      
    case (Y_minus)

      j = 1
!$OMP DO SCHEDULE(static)
      do k=1,kx
      do i=1,ix
        if ( ibits(bv(i,j,k), bc_face_S, bitw_5) == obc_mask ) then
          dv = div(i,j,k)
          vf(i,j-1,k,2) = dv
          a1 = a1 + dv
          a2 = min(a2, dv)
          a3 = max(a3, dv)
          div(i,j,k) = 0.0
          m = m + 1.0
        end if
      end do
      end do
!$OMP END DO
      

    case (Y_plus)

      j = jx
!$OMP DO SCHEDULE(static)
      do k=1,kx
      do i=1,ix
        if ( ibits(bv(i,j,k), bc_face_N, bitw_5) == obc_mask ) then
          dv = -div(i,j,k)
          vf(i,j,k,2) = dv
          a1 = a1 + dv
          a2 = min(a2, dv)
          a3 = max(a3, dv)
          div(i,j,k) = 0.0
          m = m + 1.0
        end if
      end do
      end do
!$OMP END DO
    
    
    case (Z_minus)

      k = 1
!$OMP DO SCHEDULE(static)
      do j=1,jx
      do i=1,ix
        if ( ibits(bv(i,j,k), bc_face_B, bitw_5) == obc_mask ) then
          dv = div(i,j,k)
          vf(i,j,k-1,3) = dv
          a1 = a1 + dv
          a2 = min(a2, dv)
          a3 = max(a3, dv)
          div(i,j,k) = 0.0
          m = m + 1.0
        end if
      end do
      end do
!$OMP END DO
      

    case (Z_plus)

      k = kx
!$OMP DO SCHEDULE(static)
      do j=1,jx
      do i=1,ix
        if ( ibits(bv(i,j,k), bc_face_T, bitw_5) == obc_mask ) then
          dv = -div(i,j,k)
          vf(i,j,k,3) = dv
          a1 = a1 + dv
          a2 = min(a2, dv)
          a3 = max(a3, dv)
          div(i,j,k) = 0.0
          m = m + 1.0
        end if
      end do
      end do
!$OMP END DO
    
    case default
    end select FACES
!$OMP END PARALLEL
    
    aa(1) = a1 ! sum
    aa(2) = a2 ! min
    aa(3) = a3 ! max

    flop = flop + dble(m)*3.0d0

    return
    end subroutine vobc_div_oflow


!> ********************************************************************
!! @brief 外部流出境界条件による速度ベクトルの発散の修正
!! @param [in]     sz   配列長
!! @param [in]     g    ガイドセル長
!! @param [in]     face 面番号
!! @param [out]    sum  領域境界の積算値
!! @param [in]     vf   セルフェイス速度 n+1
!! @param [in]     bv   BCindex V
!! @param [out]    flop flop count 近似
!! @note 指定面でも固体部分は対象外とするのでループ中に判定あり
!!       div(u)=0から，内部流出境界のセルで計算されたdivが流出速度となる
!<
    subroutine vobc_div_vec (sz, g, face, sum, vf, bv, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, g, ix, jx, kx, face
    integer, dimension(3)                                     ::  sz
    double precision                                          ::  flop, rix, rjx, rkx
    real                                                      ::  sum, a
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  vf
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bv

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)

    a = 0.0   ! sum


!$OMP PARALLEL &
!$OMP REDUCTION(+:a) &
!$OMP FIRSTPRIVATE(ix, jx, kx, face) &
!$OMP PRIVATE(i, j, k)

    FACES : select case (face)

    case (X_minus)

      i = 1
!$OMP DO SCHEDULE(static)
      do k=1,kx
      do j=1,jx
        a = a + vf(i-1,j,k,1)
      end do
      end do
!$OMP END DO


    case (X_plus)

      i = ix
!$OMP DO SCHEDULE(static)
      do k=1,kx
      do j=1,jx
        a = a + vf(i,j,k,1)
      end do
      end do
!$OMP END DO


    case (Y_minus)

      j = 1
!$OMP DO SCHEDULE(static)
      do k=1,kx
      do i=1,ix
        a = a + vf(i,j-1,k,2)
      end do
      end do
!$OMP END DO


    case (Y_plus)

      j = jx
!$OMP DO SCHEDULE(static)
      do k=1,kx
      do i=1,ix
        a = a + vf(i,j,k,2)
      end do
      end do
!$OMP END DO


    case (Z_minus)

      k = 1
!$OMP DO SCHEDULE(static)
      do j=1,jx
      do i=1,ix
        a = a + vf(i,j,k-1,3)
      end do
      end do
!$OMP END DO


    case (Z_plus)

      k = kx
!$OMP DO SCHEDULE(static)
      do j=1,jx
      do i=1,ix
        a = a + vf(i,j,k,3)
      end do
      end do
!$OMP END DO

    case default
    end select FACES
!$OMP END PARALLEL

    sum = a


    rix = dble(jx)*dble(kx)
    rjx = dble(ix)*dble(kx)
    rkx = dble(ix)*dble(jx)

    FACES2 : select case (face)

    case (X_minus)
      flop = flop + rix*1.0d0

    case (X_plus)
      flop = flop + rix*1.0d0

    case (Y_minus)
      flop = flop + rjx*1.0d0

    case (Y_plus)
      flop = flop + rjx*1.0d0

    case (Z_minus)
      flop = flop + rkx*1.0d0

    case (Z_plus)
      flop = flop + rkx*1.0d0

    case default
    end select FACES2

    return
    end subroutine vobc_div_vec


!> ********************************************************************
!! @brief セルフェイスの値をセットする
!! @param [out] vf   セルフェイス速度
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル長
!! @param [in]  bv   BCindex V
!! @param [in]  face 外部境界の面番号
!! @param [in]  vec  指定する速度ベクトル
!! @note 部分的な境界条件の実装のため、ガイドセル部のマスク情報を利用
!<
    subroutine vobc_face_drchlt (vf, sz, g, bv, face, vec)
    implicit none
    include 'ffv_f_params.h'
    integer                                                     ::  i, j, k, g, face, ix, jx, kx
    integer, dimension(3)                                       ::  sz
    real                                                        ::  u_bc, v_bc, w_bc
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)   ::  vf
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bv
    real, dimension(3)                                          ::  vec

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)

    ! u_bcは境界速度
    u_bc = vec(1)
    v_bc = vec(2)
    w_bc = vec(3)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, u_bc, v_bc, w_bc, face)

    FACES : select case (face)

    case (X_minus)

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do j=1,jx
        if ( ibits(bv(1, j, k), bc_face_W, bitw_5) == obc_mask ) then
          vf(0, j, k, 1) = u_bc * real(ibits(bv(0, j, k), State, 1))
        endif
      end do
      end do
!$OMP END DO


    case (X_plus)

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do j=1,jx
        if ( ibits(bv(ix, j, k), bc_face_E, bitw_5) == obc_mask ) then
          vf(ix, j, k, 1) = u_bc * real(ibits(bv(ix+1, j, k), State, 1))
        endif
      end do
      end do
!$OMP END DO


    case (Y_minus)

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do i=1,ix
        if ( ibits(bv(i, 1, k), bc_face_S, bitw_5) == obc_mask ) then
          vf(i, 0, k, 2) = v_bc * real(ibits(bv(i, 0, k), State, 1))
        endif
      end do
      end do
!$OMP END DO


    case (Y_plus)

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do i=1,ix
        if ( ibits(bv(i, jx, k), bc_face_N, bitw_5) == obc_mask ) then
          vf(i, jx, k, 2) = v_bc * real(ibits(bv(i, jx+1, k), State, 1))
        endif
      end do
      end do
!$OMP END DO


    case (Z_minus)

!$OMP DO SCHEDULE(static)
      do j=1,jx
      do i=1,ix
        if ( ibits(bv(i, j ,1), bc_face_B, bitw_5) == obc_mask ) then
          vf(i, j, 0, 3) = w_bc * real(ibits(bv(i, j, 0), State, 1))
        endif
      end do
      end do
!$OMP END DO


    case (Z_plus)

!$OMP DO SCHEDULE(static)
      do j=1,jx
      do i=1,ix
        if ( ibits(bv(i, j, kx), bc_face_T, bitw_5) == obc_mask ) then
          vf(i, j, kx, 3) = w_bc * real(ibits(bv(i, j, kx+1), State, 1))
        endif
      end do
      end do
!$OMP END DO

    case default
    end select FACES
!$OMP END PARALLEL

    return
    end subroutine vobc_face_drchlt
