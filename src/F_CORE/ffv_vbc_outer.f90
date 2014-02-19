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

!> @file   ffv_bc_outer.f90
!! @brief  外部境界条件
!! @author aics
!<


!> ********************************************************************
!! @brief 外部速度境界条件による対流項と粘性項の流束の修正
!! @param [out] wv     疑似ベクトルの空間項の評価値
!! @param [in]  sz     配列長
!! @param [in]  g      ガイドセル長
!! @param [in]  m_face 外部境界処理のときの面番号
!! @param [in]  dh     格子幅
!! @param [in]  rei    Reynolds数の逆数
!! @param [in]  v0     速度ベクトル（n-step）
!! @param [in]  bv     BCindex C
!! @param [in]  vec    指定する速度ベクトル
!! @param [in]  nID    隣接ランク番号（nID[]<0の時外部境界面）
!! @param [out] flop   浮動小数点演算数
!! @note vecには，流入条件のとき指定速度
!!  mskで部分的な速度を与える
!<
    subroutine vobc_pv_specv (wv, sz, g, m_face, dh, rei, v0, bv, vec, nID, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, g, face, m_face
    integer                                                   ::  ix, jx, kx
    integer, dimension(3)                                     ::  sz
    double precision                                          ::  flop
    real                                                      ::  Up, Vp, Wp, Ur, Vr, Wr
    real                                                      ::  dh, dh1, dh2, EX, EY, EZ, rei
    real                                                      ::  fu, fv, fw, c, ac, msk
    real                                                      ::  u_bc, v_bc, w_bc, m
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v0, wv
    real, dimension(3)                                        ::  vec
    integer, dimension(0:5)                                   ::  nID
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bv

    if ( nID(m_face) >= 0 ) return

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    face = m_face
    
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
!! @param [out] wv     疑似ベクトルの空間項の評価値
!! @param [in]  sz     配列長
!! @param [in]  g      ガイドセル長
!! @param [in]  m_face 外部境界処理のときの面番号
!! @param [in]  dh     格子幅
!! @param [in]  rei    Reynolds数の逆数
!! @param [in]  v0     速度ベクトル（n-step）
!! @param [in]  vec    指定する速度ベクトル
!! @param [in]  nID    隣接ランク番号（nID[]<0の時外部境界面）
!! @param [out] flop   浮動小数点演算数
!! @note 対流流束はゼロ，壁面法線方向の1階微分もゼロ
!<
    subroutine vobc_pv_wall (wv, sz, g, m_face, dh, rei, v0, vec, nID, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, g, face, m_face
    integer                                                   ::  ix, jx, kx
    integer, dimension(3)                                     ::  sz
    double precision                                          ::  flop, rix, rjx, rkx
    real                                                      ::  u_bc, v_bc, w_bc
    real                                                      ::  dh, dh1, dh2, rei
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v0, wv
    real, dimension(3)                                        ::  vec
    integer, dimension(0:5)                                   ::  nID

    if ( nID(m_face) >= 0 ) return

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    face = m_face

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

    FACES2 : select case (face)

    case (X_minus)
      flop = flop + rix*6.0d0

    case (X_plus)
      flop = flop + rix*6.0d0

    case (Y_minus)
      flop = flop + rjx*6.0d0

    case (Y_plus)
      flop = flop + rjx*6.0d0

    case (Z_minus)
      flop = flop + rkx*6.0d0

    case (Z_plus)
      flop = flop + rkx*6.0d0

    case default
    end select FACES2

    return
    end subroutine vobc_pv_wall


!> ********************************************************************
!! @brief ガイドセルの速度指定境界条件を設定するために必要な参照値をセットする
!! @param [out] v      セルセンタ速度
!! @param [in]  sz     配列長
!! @param [in]  g      ガイドセル長
!! @param [in]  m_face 外部境界の面番号
!! @param [in]  bv     BCindex C
!! @param [in]  vec    指定する速度ベクトル
!! @param [in]  nID    隣接ランク番号（nID[]<0の時外部境界面）
!! @note 流束型の境界条件を用いるので，内点の計算に使う参照点に値があればよい（1層）
!<
    subroutine vobc_drchlt (v, sz, g, m_face, bv, vec, nID)
    implicit none
    include 'ffv_f_params.h'
    integer                                                     ::  i, j, k, g, face, ix, jx, kx, gc, m_face
    integer, dimension(3)                                       ::  sz
    real                                                        ::  u_bc, v_bc, w_bc
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)   ::  v
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bv
    real, dimension(3)                                          ::  vec
    integer, dimension(0:5)                                     ::  nID

    if ( nID(m_face) >= 0 ) return

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    gc = g
    face = m_face
    
    ! u_bcは境界速度
    u_bc = vec(1)
    v_bc = vec(2)
    w_bc = vec(3)
    
!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, u_bc, v_bc, w_bc, face, gc) &
!$OMP PRIVATE(i, j, k)

    FACES : select case (face)
    
    case (X_minus)

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do j=1,jx
      do i=1-gc, 0
        v(i, j, k, 1) = u_bc
        v(i, j, k, 2) = v_bc
        v(i, j, k, 3) = w_bc
      end do
      end do
      end do
!$OMP END DO

      
    case (X_plus)

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do j=1,jx
      do i=ix+1, ix+gc
        v(i, j, k, 1) = u_bc
        v(i, j, k, 2) = v_bc
        v(i, j, k, 3) = w_bc
      end do
      end do
      end do
!$OMP END DO

      
    case (Y_minus)

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do i=1,ix
      do j=1-gc, 0
        v(i, j, k, 1) = u_bc
        v(i, j, k, 2) = v_bc
        v(i, j, k, 3) = w_bc
      end do
      end do
      end do
!$OMP END DO

      
    case (Y_plus)

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do i=1,ix
      do j=jx+1, jx+gc
        v(i, j, k, 1) = u_bc
        v(i, j, k, 2) = v_bc
        v(i, j, k, 3) = w_bc
      end do
      end do
      end do
!$OMP END DO

      
    case (Z_minus)

!$OMP DO SCHEDULE(static)
      do j=1,jx
      do i=1,ix
      do k=1-gc, 0
        v(i, j, k, 1) = u_bc
        v(i, j, k, 2) = v_bc
        v(i, j, k, 3) = w_bc
      end do
      end do
      end do
!$OMP END DO

      
    case (Z_plus)

!$OMP DO SCHEDULE(static)
      do j=1,jx
      do i=1,ix
      do k=kx+1, kx+gc
        v(i, j, k, 1) = u_bc
        v(i, j, k, 2) = v_bc
        v(i, j, k, 3) = w_bc
      end do
      end do
      end do
!$OMP END DO

    case default
    end select FACES
!$OMP END PARALLEL
    
    return
    end subroutine vobc_drchlt


!> ********************************************************************
!! @brief ノイマン条件
!! @param [in,out] v      速度ベクトル
!! @param [in]     sz     配列長
!! @param [in]     g      ガイドセル長
!! @param [in]     m_face 外部境界面の番号
!! @param [out]    aa     積算値 \sum{v}
!! @param [in]     nID    隣接ランク番号（nID[]<0の時外部境界面）
!<
    subroutine vobc_neumann (v, sz, g, m_face, aa, nID)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, face, g, gc, m_face
    real                                                      ::  ux, uy, uz, aa
    integer, dimension(3)                                     ::  sz
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v
    integer, dimension(0:5)                                   ::  nID

    if ( nID(m_face) >= 0 ) return

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    gc = g
    face = m_face

    aa = 0.0

!$OMP PARALLEL REDUCTION(+:aa) &
!$OMP FIRSTPRIVATE(ix, jx, kx, gc, face) &
!$OMP PRIVATE(i, j, k, ux, uy, uz)

    FACES : select case (face)
    case (X_minus)

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do j=1,jx
        ux = v(1, j, k, 1)
        uy = v(1, j, k, 2)
        uz = v(1, j, k, 3)
        aa = aa + ux

        do i=1-gc, 0
          v(i, j, k, 1) = ux
          v(i, j, k, 2) = uy
          v(i, j, k, 3) = uz
        end do
      end do
      end do
!$OMP END DO


    case (X_plus)

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do j=1,jx
        ux = v(ix, j, k, 1)
        uy = v(ix, j, k, 2)
        uz = v(ix, j, k, 3)
        aa = aa + ux

        do i=ix+1, ix+gc
          v(i, j, k, 1) = ux
          v(i, j, k, 2) = uy
          v(i, j, k, 3) = uz
        end do
      end do
      end do
!$OMP END DO


    case (Y_minus)

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do i=1,ix
        ux = v(i, 1, k, 1)
        uy = v(i, 1, k, 2)
        uz = v(i, 1, k, 3)
        aa = aa + uy

        do j=1-gc, 0
          v(i, j, k, 1) = ux
          v(i, j, k, 2) = uy
          v(i, j, k, 3) = uz
        end do
      end do
      end do
!$OMP END DO


    case (Y_plus)

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do i=1,ix
        ux = v(i, jx, k, 1)
        uy = v(i, jx, k, 2)
        uz = v(i, jx, k, 3)
        aa = aa + uy

        do j=jx+1, jx+gc
          v(i, j, k, 1) = ux
          v(i, j, k, 2) = uy
          v(i, j, k, 3) = uz
        end do
      end do
      end do
!$OMP END DO


    case (Z_minus)

!$OMP DO SCHEDULE(static)
      do j=1,jx
      do i=1,ix
        ux = v(i, j, 1, 1)
        uy = v(i, j, 1, 2)
        uz = v(i, j, 1, 3)
        aa = aa + uz

        do k=1-gc, 0
          v(i, j, k, 1) = ux
          v(i, j, k, 2) = uy
          v(i, j, k, 3) = uz
        end do
      end do
      end do
!$OMP END DO


    case (Z_plus)

!$OMP DO SCHEDULE(static)
      do j=1,jx
      do i=1,ix
        ux = v(i, j, kx, 1)
        uy = v(i, j, kx, 2)
        uz = v(i, j, kx, 3)
        aa = aa + uz

        do k=kx+1, kx+gc
          v(i, j, k, 1) = ux
          v(i, j, k, 2) = uy
          v(i, j, k, 3) = uz
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
!! @brief 速度の外部境界：　トラクションフリー
!! @param [in,out] v      セルセンタ速度ベクトル
!! @param [in]     sz     配列長
!! @param [in]     g      ガイドセル長
!! @param [in]     m_face 外部境界面の番号
!! @param [in,out] vf     セルフェイス速度
!! @param [out]    sum    積算値 \sum{v}
!! @param [in]     nID    隣接ランク番号（nID[]<0の時外部境界面）
!! @note 今のところ、トラクションフリー面は全て流体
!! 後半のループは1~ix，flopはtfree1,2の合計
!<
    subroutine vobc_tfree2 (v, sz, g, m_face, vf, sum, nID)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, face, g, gc, m_face
    integer, dimension(3)                                     ::  sz
    real                                                      ::  ut, vt, wt, aa, sum
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v, vf
    integer, dimension(0:5)                                   ::  nID

    if ( nID(m_face) >= 0 ) return

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    gc = g
    face = m_face

    aa = 0.0

!$OMP PARALLEL REDUCTION(+:aa) &
!$OMP FIRSTPRIVATE(ix, jx, kx, gc, face) &
!$OMP PRIVATE(i, j, k, ut, vt, wt)

    FACES : select case (face)
    case (X_minus)

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do j=1,jx
        ut = v(1, j, k, 1)
        vt = 0.5 * (vf(0, j, k, 2) + vf(0, j-1, k  , 2))
        wt = 0.5 * (vf(0, j, k, 3) + vf(0, j  , k-1, 3))

        aa = aa + vf(0, j, k, 1)

        do i=1-gc, 0
          v(i, j, k, 1) = ut
          v(i, j, k, 2) = vt
          v(i, j, k, 3) = wt
        end do

      end do
      end do
!$OMP END DO
      

    case (X_plus)

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do j=1,jx
        ut = v(ix, j, k, 1)
        vt = 0.5 * (vf(ix+1, j, k, 2) + vf(ix+1, j-1, k  , 2))
        wt = 0.5 * (vf(ix+1, j, k, 3) + vf(ix+1, j  , k-1, 3))

        aa = aa + vf(ix, j, k, 1)

        do i=ix+1, ix+gc
          v(i, j, k, 1) = ut
          v(i, j, k, 2) = vt
          v(i, j, k, 3) = wt
        end do

      end do
      end do
!$OMP END DO
      

    case (Y_minus)

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do i=1,ix
        ut = 0.5 * (vf(i, 0, k, 1) + vf(i-1, 0, k  , 1))
        vt = v(i, 1, k, 2)
        wt = 0.5 * (vf(i, 0, k, 3) + vf(i  , 0, k-1, 3))

        aa = aa + vf(i, 0, k, 2)

        do j=1-gc, 0
          v(i, j, k, 1) = ut
          v(i, j, k, 2) = vt
          v(i, j, k, 3) = wt
        end do

      end do
      end do
!$OMP END DO
      

    case (Y_plus)

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do i=1,ix
        ut = 0.5 * (vf(i, jx+1, k, 1) + vf(i-1, jx+1, k  , 1))
        vt = v(i, jx, k, 2)
        wt = 0.5 * (vf(i, jx+1, k, 3) + vf(i  , jx+1, k-1, 3))

        aa = aa + vf(i, jx, k, 2)

        do j=jx+1, jx+gc
          v(i, j, k, 1) = ut
          v(i, j, k, 2) = vt
          v(i, j, k, 3) = wt
        end do

      end do
      end do
!$OMP END DO
      

    case (Z_minus)

!$OMP DO SCHEDULE(static)
      do j=1,jx
      do i=1,ix
        ut = 0.5 * (vf(i, j, 0, 1) + vf(i-1, j  , 0, 1))
        vt = 0.5 * (vf(i, j, 0, 2) + vf(i  , j-1, 0, 2))
        wt = v(i, j, 1, 3)

        aa = aa + vf(i, j, 0, 3)

        do k=1-gc, 0
          v(i, j, k, 1) = ut
          v(i, j, k, 2) = vt
          v(i, j, k, 3) = wt
        end do

      end do
      end do
!$OMP END DO
      

    case (Z_plus)

!$OMP DO SCHEDULE(static)
      do j=1,jx
      do i=1,ix
        ut = 0.5 * (vf(i, j, kx+1, 1) + vf(i-1, j  , kx+1, 1))
        vt = 0.5 * (vf(i, j, kx+1, 2) + vf(i  , j-1, kx+1, 2))
        wt = v(i, j, kx, 3)

        aa = aa + vf(i, j, kx, 3)

        do k=kx+1, kx+gc
          v(i, j, k, 1) = ut
          v(i, j, k, 2) = vt
          v(i, j, k, 3) = wt
        end do

      end do
      end do
!$OMP END DO
      
    case default
    end select FACES

!$OMP END PARALLEL


    sum = aa

    return
    end subroutine vobc_tfree2


!> ********************************************************************
!! @brief 速度の外部境界：　トラクションフリーの前半処理
!! @param [in,out] vf     セルフェイス速度
!! @param [in]     sz     配列長
!! @param [in]     g      ガイドセル長
!! @param [in]     m_face 外部境界面の番号
!! @param [in]     nID    隣接ランク番号（nID[]<0の時外部境界面）
!! ループは0~ixで回す．0から始めるのは，後半処理でスタガード位置の変数からセンターへ内挿するときに
!! 参照する端点を含めるため．
!! 処理の間に同期が必要なので，サブルーチンを分割
!<
  subroutine vobc_tfree1 (vf, sz, g, m_face, nID)
  implicit none
  include 'ffv_f_params.h'
  integer                                                   ::  i, j, k, ix, jx, kx, face, g, m_face
  integer, dimension(3)                                     ::  sz
  integer, dimension(0:5)                                   ::  nID
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  vf

  if ( nID(m_face) >= 0 ) return

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
  face = m_face

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, face)

  FACES : select case (face)
  case (X_minus)

!$OMP DO SCHEDULE(static) PRIVATE(j, k)
    do k=0,kx
    do j=0,jx
      vf(0, j, k, 2) = vf(1, j, k, 2) + vf(0, j+1, k  , 1) - vf(0, j, k, 1)
      vf(0, j, k, 3) = vf(1, j, k, 3) + vf(0, j  , k+1, 1) - vf(0, j, k, 1)
    end do
    end do
!$OMP END DO


  case (X_plus)

!$OMP DO SCHEDULE(static) PRIVATE(j, k)
    do k=0,kx
    do j=0,jx
      vf(ix+1, j, k, 2) = vf(ix, j, k, 2) - vf(ix, j+1, k  , 1) + vf(ix, j, k, 1)
      vf(ix+1, j, k, 3) = vf(ix, j, k, 3) - vf(ix, j  , k+1, 1) + vf(ix, j, k, 1)
    end do
    end do
!$OMP END DO


  case (Y_minus)

!$OMP DO SCHEDULE(static) PRIVATE(i, k)
    do k=0,kx
    do i=0,ix
      vf(i, 0, k, 1) = vf(i, 1, k, 1) + vf(i+1, 0, k  , 2) - vf(i, 0, k, 2)
      vf(i, 0, k, 3) = vf(i, 1, k, 3) + vf(i  , 0, k+1, 2) - vf(i, 0, k, 2)
    end do
    end do
!$OMP END DO


  case (Y_plus)

!$OMP DO SCHEDULE(static) PRIVATE(i, k)
    do k=0,kx
    do i=0,ix
      vf(i, jx+1, k, 1) = vf(i, jx, k, 1) - vf(i+1, jx, k  , 2) + vf(i, jx, k, 2)
      vf(i, jx+1, k, 3) = vf(i, jx, k, 3) - vf(i  , jx, k+1, 2) + vf(i, jx, k, 2)
    end do
    end do
!$OMP END DO


  case (Z_minus)

!$OMP DO SCHEDULE(static) PRIVATE(i, j)
    do j=0,jx
    do i=0,ix
      vf(i, j, 0, 1) = vf(i, j, 1, 1) + vf(i+1, j  , 0, 3) - vf(i, j, 0, 3)
      vf(i, j, 0, 2) = vf(i, j, 1, 2) + vf(i  , j+1, 0, 3) - vf(i, j, 0, 3)
    end do
    end do
!$OMP END DO


  case (Z_plus)

!$OMP DO SCHEDULE(static) PRIVATE(i, j)
    do j=0,jx
    do i=0,ix
      vf(i, j, kx+1, 1) = vf(i, j, kx, 1) - vf(i+1, j  , kx, 3) + vf(i, j, kx, 3)
      vf(i, j, kx+1, 2) = vf(i, j, kx, 2) - vf(i  , j+1, kx, 3) + vf(i, j, kx, 3)
    end do
    end do
!$OMP END DO


  case default
  end select FACES

!$OMP END PARALLEL

  return
  end subroutine vobc_tfree1


!> ********************************************************************
!! @brief 疑似速度から次ステップ速度へ参照する速度をコピーする
!! @param [out] v      速度ベクトル（セルセンタ）
!! @param [in]  sz     配列長
!! @param [in]  g      ガイドセル長
!! @param [in]  m_face 面番号
!! @param [in]  vc     セルセンタ疑似速度 u^*
!! @param [in]  nID    隣接ランク番号（nID[]<0の時外部境界面）
!<
    subroutine vobc_update (v, sz, g, m_face, vc, nID)
    implicit none
    include 'ffv_f_params.h'
    integer                                                     ::  i, j, k, g, ix, jx, kx, face, gc, m_face
    integer, dimension(3)                                       ::  sz
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)   ::  v, vc
    integer, dimension(0:5)                                     ::  nID

    if ( nID(m_face) >= 0 ) return

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    gc = g
    face = m_face

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, face, gc) PRIVATE(i, j, k)

    FACES : select case (face)
    case (X_minus)

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do j=1,jx
      do i=1-gc, 0
        v(i, j, k, 1) = vc(i, j, k, 1)
        v(i, j, k, 2) = vc(i, j, k, 2)
        v(i, j, k, 3) = vc(i, j, k, 3)
      end do
      end do
      end do
!$OMP END DO

    case (X_plus)

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do j=1,jx
      do i=ix+1, ix+gc
        v(i, j, k, 1) = vc(i, j, k, 1)
        v(i, j, k, 2) = vc(i, j, k, 2)
        v(i, j, k, 3) = vc(i, j, k, 3)
      end do
      end do
      end do
!$OMP END DO
      
    case (Y_minus)

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do i=1,ix
      do j=1-gc, 0
        v(i, j, k, 1) = vc(i, j, k, 1)
        v(i, j, k, 2) = vc(i, j, k, 2)
        v(i, j, k, 3) = vc(i, j, k, 3)
      end do
      end do
      end do
!$OMP END DO
      
    case (Y_plus)

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do i=1,ix
      do j=jx+1, jx+gc
        v(i, j, k, 1) = vc(i, j, k, 1)
        v(i, j, k, 2) = vc(i, j, k, 2)
        v(i, j, k, 3) = vc(i, j, k, 3)
      end do
      end do
      end do
!$OMP END DO
      
    case (Z_minus)

!$OMP DO SCHEDULE(static)
      do j=1,jx
      do i=1,ix
      do k=1-gc, 0
        v(i, j, k, 1) = vc(i, j, k, 1)
        v(i, j, k, 2) = vc(i, j, k, 2)
        v(i, j, k, 3) = vc(i, j, k, 3)
      end do
      end do
      end do
!$OMP END DO
      
    case (Z_plus)

!$OMP DO SCHEDULE(static)
      do j=1,jx
      do i=1,ix
      do k=kx+1, kx+gc
        v(i, j, k, 1) = vc(i, j, k, 1)
        v(i, j, k, 2) = vc(i, j, k, 2)
        v(i, j, k, 3) = vc(i, j, k, 3)
      end do
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
!! @param [in,out] div     速度の発散
!! @param [in]     sz      配列長
!! @param [in]     g       ガイドセル長
!! @param [in]     m_face  面番号
!! @param [in]     bv      BCindex C
!! @param [in]     vec     指定する速度ベクトル
!! @param [out]    vsum    \sum{v}
!! @param [in]     nID     隣接ランク番号（nID[]<0の時外部境界面）
!! @note 固体部分は対象外とするのでループ中に判定あり
!!       部分的な境界条件の実装のため、ガイドセル部のマスク情報を利用
!<
    subroutine vobc_div_drchlt (div, sz, g, m_face, bv, vec, vsum, nID)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, g, ix, jx, kx, face, bvx, m_face
    integer, dimension(3)                                     ::  sz
    real                                                      ::  u_bc, v_bc, w_bc, vsum, b, a
    real, dimension(3)                                        ::  vec
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  div
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bv
    integer, dimension(0:5)                                   ::  nID

    if ( nID(m_face) >= 0 ) return

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    face = m_face
    
    u_bc = vec(1)
    v_bc = vec(2)
    w_bc = vec(3)

    b = 0.0

!$OMP PARALLEL &
!$OMP REDUCTION(+:b) &
!$OMP FIRSTPRIVATE(ix, jx, kx, u_bc, v_bc, w_bc, face) &
!$OMP PRIVATE(a)

    FACES : select case (face)
    case (X_minus)

!$OMP DO SCHEDULE(static) PRIVATE(j, k)
      do k=1,kx
      do j=1,jx
        bvx = bv(1, j, k)
        if ( ibits(bvx, bc_face_W, bitw_5) == obc_mask ) then
          a = u_bc * real(ibits(bvx, State, 1)) * real(ibits(bv(0, j, k), State, 1))
          div(1, j, k) = div(1, j, k) - a
          b = b + a
        endif
      end do
      end do
!$OMP END DO
      
      
    case (X_plus)

!$OMP DO SCHEDULE(static) PRIVATE(j, k)
      do k=1,kx
      do j=1,jx
        bvx = bv(ix, j, k)
        if ( ibits(bvx, bc_face_E, bitw_5) == obc_mask ) then
          a = u_bc * real(ibits(bvx, State, 1)) * real(ibits(bv(ix+1, j, k), State, 1))
          div(ix, j, k) = div(ix, j, k) + a
          b = b + a
        endif
      end do
      end do
!$OMP END DO
      
      
    case (Y_minus)

!$OMP DO SCHEDULE(static) PRIVATE(i, k)
      do k=1,kx
      do i=1,ix
        bvx = bv(i, 1, k)
        if ( ibits(bvx, bc_face_S, bitw_5) == obc_mask ) then
          a = v_bc * real(ibits(bvx, State, 1)) * real(ibits(bv(i, 0, k), State, 1))
          div(i, 1, k) = div(i, 1, k) - a
          b = b + a
        endif
      end do
      end do
!$OMP END DO
      
      
    case (Y_plus)

!$OMP DO SCHEDULE(static) PRIVATE(i, k)
      do k=1,kx
      do i=1,ix
        bvx = bv(i, jx, k)
        if ( ibits(bvx, bc_face_N, bitw_5) == obc_mask ) then
          a = v_bc * real(ibits(bvx, State, 1)) * real(ibits(bv(i, jx+1, k), State, 1))
          div(i, jx, k) = div(i, jx, k) + a
          b = b + a
        endif
      end do
      end do
!$OMP END DO
    
    
    case (Z_minus)

!$OMP DO SCHEDULE(static) PRIVATE(i, j)
      do j=1,jx
      do i=1,ix
        bvx = bv(i, j, 1)
        if ( ibits(bvx, bc_face_B, bitw_5) == obc_mask ) then
          a = w_bc * real(ibits(bvx, State, 1)) * real(ibits(bv(i, j, 0), State, 1))
          div(i, j, 1) = div(i, j, 1) - a
          b = b + a
        endif
      end do
      end do
!$OMP END DO
      
      
    case (Z_plus)

!$OMP DO SCHEDULE(static) PRIVATE(i, j)
      do j=1,jx
      do i=1,ix
        bvx = bv(i, j, kx)
        if ( ibits(bvx, bc_face_T, bitw_5) == obc_mask ) then
          a = w_bc * real(ibits(bvx, State, 1)) * real(ibits(bv(i, j, kx+1), State, 1))
          div(i, j, kx) = div(i, j, kx) + a
          b = b + a
        endif
      end do
      end do
!$OMP END DO
    
    case default
    end select FACES

!$OMP END PARALLEL


    return
    end subroutine vobc_div_drchlt


!> ********************************************************************
!! @brief 外部境界面の流入出量を求める
!! @param [out]    sum    領域境界の流速の積算値 \sum{vf}
!! @param [in]     sz     配列長
!! @param [in]     g      ガイドセル長
!! @param [in]     m_face 面番号
!! @param [in]     v      セルセンター速度 n+1
!! @param [in]     bv     BCindex C
!! @param [in]     nID    隣接ランク番号（nID[]<0の時外部境界面）
!! @note 有効セルのマスクを掛けて、流量を積算
!<
    subroutine vobc_get_massflow (sum, sz, g, m_face, v, bv, nID)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, g, ix, jx, kx, face, m_face
    integer, dimension(3)                                     ::  sz
    real                                                      ::  sum, a, s
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bv
    integer, dimension(0:5)                                   ::  nID

    if ( nID(m_face) >= 0 ) return

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    face = m_face

    a = 0.0   ! sum


!$OMP PARALLEL &
!$OMP REDUCTION(+:a) &
!$OMP FIRSTPRIVATE(ix, jx, kx, face) &
!$OMP PRIVATE(s)

    FACES : select case (face)

    case (X_minus)

!$OMP DO SCHEDULE(static) PRIVATE(j, k)
      do k=1,kx
      do j=1,jx
        s = real(ibits(bv(0,j,k), State, 1)) * real(ibits(bv(1,j,k), State, 1))
        a = a + 0.5*( v(0,j,k,1) + v(1,j,k,1)) * s
      end do
      end do
!$OMP END DO


    case (X_plus)

!$OMP DO SCHEDULE(static) PRIVATE(j, k)
      do k=1,kx
      do j=1,jx
        s = real(ibits(bv(ix,j,k), State, 1)) * real(ibits(bv(ix+1,j,k), State, 1))
        a = a + 0.5*( v(ix,j,k,1) + v(ix+1,j,k,1)) * s
      end do
      end do
!$OMP END DO


    case (Y_minus)

!$OMP DO SCHEDULE(static) PRIVATE(i, k)
      do k=1,kx
      do i=1,ix
        s = real(ibits(bv(i,0,k), State, 1)) * real(ibits(bv(i,1,k), State, 1))
        a = a + 0.5*( v(i,0,k,2) + v(i,1,k,2)) * s
      end do
      end do
!$OMP END DO


    case (Y_plus)

!$OMP DO SCHEDULE(static) PRIVATE(i, k)
      do k=1,kx
      do i=1,ix
        s = real(ibits(bv(i,jx,k), State, 1)) * real(ibits(bv(i,jx+1,k), State, 1))
        a = a + 0.5*( v(i,jx,k,2) + v(i,jx+1,k,2)) * s
      end do
      end do
!$OMP END DO


    case (Z_minus)

!$OMP DO SCHEDULE(static) PRIVATE(i, j)
      do j=1,jx
      do i=1,ix
        s = real(ibits(bv(i,j,0), State, 1)) * real(ibits(bv(i,j,1), State, 1))
        a = a + 0.5*( v(i,j,0,3) + v(i,j,1,3)) * s
      end do
      end do
!$OMP END DO


    case (Z_plus)

!$OMP DO SCHEDULE(static) PRIVATE(i, j)
      do j=1,jx
      do i=1,ix
        s = real(ibits(bv(i,j,kx), State, 1)) * real(ibits(bv(i,j,kx+1), State, 1))
        a = a + 0.5*( v(i,j,kx,3) + v(i,j,kx+1,3)) * s
      end do
      end do
!$OMP END DO

    case default
    end select FACES
!$OMP END PARALLEL

    sum = a

    return
    end subroutine vobc_get_massflow


!> ********************************************************************
!! @brief セルフェイスの値をセットする
!! @param [out] vf     セルフェイス速度
!! @param [in]  sz     配列長
!! @param [in]  g      ガイドセル長
!! @param [in]  m_face 外部境界の面番号
!! @param [in]  bv     BCindex C
!! @param [in]  vec    指定する速度ベクトル
!! @param [out] vsum   \sum{vf}
!! @param [in]  nID    隣接ランク番号（nID[]<0の時外部境界面）
!! @note 部分的な境界条件の実装のため、ガイドセル部のマスク情報を利用
!<
    subroutine vobc_face_drchlt (vf, sz, g, m_face, bv, vec, vsum, nID)
    implicit none
    include 'ffv_f_params.h'
    integer                                                     ::  i, j, k, g, face, ix, jx, kx, m_face
    integer, dimension(3)                                       ::  sz
    real                                                        ::  u_bc, v_bc, w_bc, vsum, a, b
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)   ::  vf
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bv
    real, dimension(3)                                          ::  vec
    integer, dimension(0:5)                                     ::  nID

    if ( nID(m_face) >= 0 ) return

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    face = m_face

    ! u_bcは境界速度
    u_bc = vec(1)
    v_bc = vec(2)
    w_bc = vec(3)

    b = 0.0

!$OMP PARALLEL &
!$OMP REDUCTION(+:b) &
!$OMP FIRSTPRIVATE(ix, jx, kx, u_bc, v_bc, w_bc, face) &
!$OMP PRIVATE(a)

    FACES : select case (face)

    case (X_minus)

!$OMP DO SCHEDULE(static) PRIVATE(j, k)
      do k=1,kx
      do j=1,jx
        if ( ibits(bv(1, j, k), bc_face_W, bitw_5) == obc_mask ) then
          a = u_bc * real(ibits(bv(0,j,k), State, 1)) * real(ibits(bv(1,j,k), State, 1))
          vf(0, j, k, 1) = a
          b = b + a
        endif
      end do
      end do
!$OMP END DO


    case (X_plus)

!$OMP DO SCHEDULE(static) PRIVATE(j, k)
      do k=1,kx
      do j=1,jx
        if ( ibits(bv(ix, j, k), bc_face_E, bitw_5) == obc_mask ) then
          a = u_bc * real(ibits(bv(ix,j,k), State, 1)) * real(ibits(bv(ix+1,j,k), State, 1))
          vf(ix, j, k, 1) = a
          b = b + a
        endif
      end do
      end do
!$OMP END DO


    case (Y_minus)

!$OMP DO SCHEDULE(static) PRIVATE(i, k)
      do k=1,kx
      do i=1,ix
        if ( ibits(bv(i, 1, k), bc_face_S, bitw_5) == obc_mask ) then
          a = v_bc * real(ibits(bv(i,0,k), State, 1)) * real(ibits(bv(i,1,k), State, 1))
          vf(i, 0, k, 2) = a
          b = b + a
        endif
      end do
      end do
!$OMP END DO


    case (Y_plus)

!$OMP DO SCHEDULE(static) PRIVATE(i, k)
      do k=1,kx
      do i=1,ix
        if ( ibits(bv(i, jx, k), bc_face_N, bitw_5) == obc_mask ) then
          a = v_bc * real(ibits(bv(i,jx,k), State, 1)) * real(ibits(bv(i,jx+1,k), State, 1))
          vf(i, jx, k, 2) = a
          b = b + a
        endif
      end do
      end do
!$OMP END DO


    case (Z_minus)

!$OMP DO SCHEDULE(static) PRIVATE(i, j)
      do j=1,jx
      do i=1,ix
        if ( ibits(bv(i, j ,1), bc_face_B, bitw_5) == obc_mask ) then
          a = w_bc * real(ibits(bv(i,j,0), State, 1)) * real(ibits(bv(i,j,1), State, 1))
          vf(i, j, 0, 3) = a
          b = b + a
        endif
      end do
      end do
!$OMP END DO


    case (Z_plus)

!$OMP DO SCHEDULE(static) PRIVATE(i, j)
      do j=1,jx
      do i=1,ix
        if ( ibits(bv(i, j, kx), bc_face_T, bitw_5) == obc_mask ) then
          a = w_bc * real(ibits(bv(i,j,kx), State, 1)) * real(ibits(bv(i,j,kx+1), State, 1))
          vf(i, j, kx, 3) = a
          b = b + a
        endif
      end do
      end do
!$OMP END DO

    case default
    end select FACES
!$OMP END PARALLEL

    vsum = b

    return
    end subroutine vobc_face_drchlt


