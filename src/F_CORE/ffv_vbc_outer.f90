!********************************************************************
!
!   FFV : Frontflow / violet Cartesian
!
!   Copyright (c) 2012 All right reserved.
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
    subroutine pvec_vobc_oflow (wv, sz, g, dh, rei, v0, bv, face, cf, flop)
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
      do i=1,1
      do k=1,kx
      do j=1,jx
        if ( ibits(bv(i,j,k), bc_face_W, bitw_5) == obc_mask ) then
          Up0 = v0(i,j,k,1)
          Vp0 = v0(i,j,k,2)
          Wp0 = v0(i,j,k,3)

          fu  = c * Up0
          fv  = c * Vp0
          fw  = c * Wp0
          
          EX = v0(i-1, j, k, 1) - Up0
          EY = v0(i-1, j, k, 2) - Vp0
          EZ = v0(i-1, j, k, 3) - Wp0

          wv(i,j,k,1) = wv(i,j,k,1) + ( -fu*dh1 + EX*dh2 )
          wv(i,j,k,2) = wv(i,j,k,2) + ( -fv*dh1 + EY*dh2 )
          wv(i,j,k,3) = wv(i,j,k,3) + ( -fw*dh1 + EZ*dh2 )
          m = m + 1.0
        end if
      end do
      end do
      end do
!$OMP END DO


    case (X_plus)

      if ( cf<0.0 ) c=0.0

!$OMP DO SCHEDULE(static)
      do i=ix,ix
      do k=1,kx
      do j=1,jx
        if ( ibits(bv(i,j,k), bc_face_E, bitw_5) == obc_mask ) then
          Up0 = v0(i,j,k,1)
          Vp0 = v0(i,j,k,2)
          Wp0 = v0(i,j,k,3)

          fu  = c * Up0
          fv  = c * Vp0
          fw  = c * Wp0

          EX = v0(i+1, j, k, 1) - Up0
          EY = v0(i+1, j, k, 2) - Vp0
          EZ = v0(i+1, j, k, 3) - Wp0

          wv(i,j,k,1) = wv(i,j,k,1) + ( -fu*dh1 + EX*dh2 )
          wv(i,j,k,2) = wv(i,j,k,2) + ( -fv*dh1 + EY*dh2 )
          wv(i,j,k,3) = wv(i,j,k,3) + ( -fw*dh1 + EZ*dh2 )
          m = m + 1.0
        end if
      end do
      end do
      end do
!$OMP END DO

      
    case (Y_minus)

      if ( cf>0.0 ) c=0.0

!$OMP DO SCHEDULE(static)
      do j=1,1
      do k=1,kx
      do i=1,ix
        if ( ibits(bv(i,j,k), bc_face_S, bitw_5) == obc_mask ) then
          Up0 = v0(i,j,k,1)
          Vp0 = v0(i,j,k,2)
          Wp0 = v0(i,j,k,3)

          fu  = c * Up0
          fv  = c * Vp0
          fw  = c * Wp0
        
          EX = v0(i, j-1, k, 1) - Up0
          EY = v0(i, j-1, k, 2) - Vp0
          EZ = v0(i, j-1, k, 3) - Wp0
          
          wv(i,j,k,1) = wv(i,j,k,1) + ( fu*dh1 + EX*dh2 )
          wv(i,j,k,2) = wv(i,j,k,2) + ( fv*dh1 + EY*dh2 )
          wv(i,j,k,3) = wv(i,j,k,3) + ( fw*dh1 + EZ*dh2 )
          m = m + 1.0
        end if
      end do
      end do
      end do
!$OMP END DO

      
    case (Y_plus)

      if ( cf<0.0 ) c=0.0

!$OMP DO SCHEDULE(static)
      do j=jx,jx
      do k=1,kx
      do i=1,ix
        if ( ibits(bv(i,j,k), bc_face_N, bitw_5) == obc_mask ) then
          Up0 = v0(i,j,k,1)
          Vp0 = v0(i,j,k,2)
          Wp0 = v0(i,j,k,3)

          fu  = c * Up0
          fv  = c * Vp0
          fw  = c * Wp0

          EX = v0(i, j+1, k, 1) - Up0
          EY = v0(i, j+1, k, 2) - Vp0
          EZ = v0(i, j+1, k, 3) - Wp0
          
          wv(i,j,k,1) = wv(i,j,k,1) + ( -fu*dh1 + EX*dh2 )
          wv(i,j,k,2) = wv(i,j,k,2) + ( -fv*dh1 + EY*dh2 )
          wv(i,j,k,3) = wv(i,j,k,3) + ( -fw*dh1 + EZ*dh2 )
          m = m + 1.0
        end if
      end do
      end do
      end do
!$OMP END DO

      
    case (Z_minus)

      if ( cf>0.0 ) c=0.0

!$OMP DO SCHEDULE(static)
      do k=1,1
      do j=1,jx
      do i=1,ix
        if ( ibits(bv(i,j,k), bc_face_B, bitw_5) == obc_mask ) then
          Up0 = v0(i,j,k,1)
          Vp0 = v0(i,j,k,2)
          Wp0 = v0(i,j,k,3)

          fu  = c * Up0
          fv  = c * Vp0
          fw  = c * Wp0
          
          EX = v0(i, j, k-1, 1) - Up0
          EY = v0(i, j, k-1, 2) - Vp0
          EZ = v0(i, j, k-1, 3) - Wp0
          
          wv(i,j,k,1) = wv(i,j,k,1) + ( fu*dh1 + EX*dh2 )
          wv(i,j,k,2) = wv(i,j,k,2) + ( fv*dh1 + EY*dh2 )
          wv(i,j,k,3) = wv(i,j,k,3) + ( fw*dh1 + EZ*dh2 )
          m = m + 1.0
        end if
      end do
      end do
      end do
!$OMP END DO

      
    case (Z_plus)

      if ( cf<0.0 ) c=0.0

!$OMP DO SCHEDULE(static)
      do k=kx,kx
      do j=1,jx
      do i=1,ix
        if ( ibits(bv(i,j,k), bc_face_T, bitw_5) == obc_mask ) then
          Up0 = v0(i,j,k,1)
          Vp0 = v0(i,j,k,2)
          Wp0 = v0(i,j,k,3)

          fu  = c * Up0
          fv  = c * Vp0
          fw  = c * Wp0
          
          EX = v0(i, j, k+1, 1) - Up0
          EY = v0(i, j, k+1, 2) - Vp0
          EZ = v0(i, j, k+1, 3) - Wp0
          
          wv(i,j,k,1) = wv(i,j,k,1) + ( -fu*dh1 + EX*dh2 )
          wv(i,j,k,2) = wv(i,j,k,2) + ( -fv*dh1 + EY*dh2 )
          wv(i,j,k,3) = wv(i,j,k,3) + ( -fw*dh1 + EZ*dh2 )
          m = m + 1.0
        endif
      end do
      end do
      end do
!$OMP END DO
      
    case default

    end select FACES
    
!$OMP END PARALLEL

    flop = flop + dble(m)*18.0d0
      
    return
    end subroutine pvec_vobc_oflow
    
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
!<
    subroutine pvec_vobc_specv (wv, sz, g, dh, rei, v0, bv, vec, face, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, g, face
    integer                                                   ::  ix, jx, kx
    integer, dimension(3)                                     ::  sz
    double precision                                          ::  flop
    real                                                      ::  Up0, Ue1, Uw1, Us1, Un1, Ub1, Ut1
    real                                                      ::  Vp0, Ve1, Vw1, Vs1, Vn1, Vb1, Vt1
    real                                                      ::  Wp0, We1, Ww1, Ws1, Wn1, Wb1, Wt1
    real                                                      ::  dh, dh1, dh2, EX, EY, EZ, rei
    real                                                      ::  fu, fv, fw, c, ac
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
!$OMP PRIVATE(Up0, Ue1, Uw1, Us1, Un1, Ub1, Ut1) &
!$OMP PRIVATE(Vp0, Ve1, Vw1, Vs1, Vn1, Vb1, Vt1) &
!$OMP PRIVATE(Wp0, We1, Ww1, Ws1, Wn1, Wb1, Wt1) &
!$OMP PRIVATE(fu, fv, fw, EX, EY, EZ, c, ac)

    FACES : select case (face)
    case (X_minus)
      
!$OMP DO SCHEDULE(static)
      do i=1,1
      do k=1,kx
      do j=1,jx
        if ( ibits(bv(i,j,k), bc_face_W, bitw_5) == obc_mask ) then
          Up0 = v0(i,j,k,1)
          Vp0 = v0(i,j,k,2)
          Wp0 = v0(i,j,k,3)
          
          Uw1 = u_bc
          Vw1 = v_bc
          Ww1 = w_bc
          c   = u_bc
          ac  = abs(c)
          fu  = 0.5*(c*(Up0+Uw1) - ac*(Up0-Uw1))
          fv  = 0.5*(c*(Vp0+Vw1) - ac*(Vp0-Vw1))
          fw  = 0.5*(c*(Wp0+Ww1) - ac*(Wp0-Ww1))
          
          EX = Uw1 - Up0
          EY = Vw1 - Vp0
          EZ = Ww1 - Wp0

          wv(i,j,k,1) = wv(i,j,k,1) + ( fu*dh1 + EX*dh2 )
          wv(i,j,k,2) = wv(i,j,k,2) + ( fv*dh1 + EY*dh2 )
          wv(i,j,k,3) = wv(i,j,k,3) + ( fw*dh1 + EZ*dh2 )
          m = m + 1.0
        endif
      end do
      end do
      end do
!$OMP END DO

      
    case (X_plus)
      
!$OMP DO SCHEDULE(static)
      do i=ix,ix
      do k=1,kx
      do j=1,jx
        if ( ibits(bv(i,j,k), bc_face_E, bitw_5) == obc_mask ) then
          Up0 = v0(i,j,k,1)
          Vp0 = v0(i,j,k,2)
          Wp0 = v0(i,j,k,3)
          
          Ue1 = u_bc
          Ve1 = v_bc
          We1 = w_bc
          c   = u_bc
          ac  = abs(c)
          fu  = 0.5*(c*(Ue1+Up0) - ac*(Ue1-Up0))
          fv  = 0.5*(c*(Ve1+Vp0) - ac*(Ve1-Vp0))
          fw  = 0.5*(c*(We1+Wp0) - ac*(We1-Wp0))

          EX = Ue1 - Up0
          EY = Ve1 - Vp0
          EZ = We1 - Wp0

          wv(i,j,k,1) = wv(i,j,k,1) + ( -fu*dh1 + EX*dh2 )
          wv(i,j,k,2) = wv(i,j,k,2) + ( -fv*dh1 + EY*dh2 )
          wv(i,j,k,3) = wv(i,j,k,3) + ( -fw*dh1 + EZ*dh2 )
          m = m + 1.0
        endif
      end do
      end do
      end do
!$OMP END DO

      
    case (Y_minus)
      
!$OMP DO SCHEDULE(static)
      do j=1,1
      do k=1,kx
      do i=1,ix
        if ( ibits(bv(i,j,k), bc_face_S, bitw_5) == obc_mask ) then
          Up0 = v0(i,j,k,1)
          Vp0 = v0(i,j,k,2)
          Wp0 = v0(i,j,k,3)
          
          Us1 = u_bc
          Vs1 = v_bc
          Ws1 = w_bc
          c   = v_bc
          ac  = abs(c)
          fu  = 0.5*(c*(Up0+Us1) - ac*(Up0-Us1))
          fv  = 0.5*(c*(Vp0+Vs1) - ac*(Vp0-Vs1))
          fw  = 0.5*(c*(Wp0+Ws1) - ac*(Wp0-Ws1))
        
          EX = Us1 - Up0
          EY = Vs1 - Vp0
          EZ = Ws1 - Wp0
          
          wv(i,j,k,1) = wv(i,j,k,1) + ( fu*dh1 + EX*dh2 )
          wv(i,j,k,2) = wv(i,j,k,2) + ( fv*dh1 + EY*dh2 )
          wv(i,j,k,3) = wv(i,j,k,3) + ( fw*dh1 + EZ*dh2 )
          m = m + 1.0
        endif
      end do
      end do
      end do
!$OMP END DO

      
    case (Y_plus)
      
!$OMP DO SCHEDULE(static)
      do j=jx,jx
      do k=1,kx
      do i=1,ix
        if ( ibits(bv(i,j,k), bc_face_N, bitw_5) == obc_mask ) then
          Up0 = v0(i,j,k,1)
          Vp0 = v0(i,j,k,2)
          Wp0 = v0(i,j,k,3)
          
          Un1 = u_bc
          Vn1 = v_bc
          Wn1 = w_bc
          c   = v_bc
          ac  = abs(c)
          fu  = 0.5*(c*(Un1+Up0) - ac*(Un1-Up0))
          fv  = 0.5*(c*(Vn1+Vp0) - ac*(Vn1-Vp0))
          fw  = 0.5*(c*(Wn1+Wp0) - ac*(Wn1-Wp0))

          EX = Un1 - Up0
          EY = Vn1 - Vp0
          EZ = Wn1 - Wp0
          
          wv(i,j,k,1) = wv(i,j,k,1) + ( -fu*dh1 + EX*dh2 )
          wv(i,j,k,2) = wv(i,j,k,2) + ( -fv*dh1 + EY*dh2 )
          wv(i,j,k,3) = wv(i,j,k,3) + ( -fw*dh1 + EZ*dh2 )
          m = m + 1.0
        endif
      end do
      end do
      end do
!$OMP END DO

      
    case (Z_minus)
      
!$OMP DO SCHEDULE(static)
      do k=1,1
      do j=1,jx
      do i=1,ix
        if ( ibits(bv(i,j,k), bc_face_B, bitw_5) == obc_mask ) then
          Up0 = v0(i,j,k,1)
          Vp0 = v0(i,j,k,2)
          Wp0 = v0(i,j,k,3)
          
          Ub1 = u_bc
          Vb1 = v_bc
          Wb1 = w_bc
          c   = w_bc
          ac  = abs(c)
          fu  = 0.5*(c*(Up0+Ub1) - ac*(Up0-Ub1))
          fv  = 0.5*(c*(Vp0+Vb1) - ac*(Vp0-Vb1))
          fw  = 0.5*(c*(Wp0+Wb1) - ac*(Wp0-Wb1))
          
          EX = Ub1 - Up0
          EY = Vb1 - Vp0
          EZ = Wb1 - Wp0
          
          wv(i,j,k,1) = wv(i,j,k,1) + ( fu*dh1 + EX*dh2 )
          wv(i,j,k,2) = wv(i,j,k,2) + ( fv*dh1 + EY*dh2 )
          wv(i,j,k,3) = wv(i,j,k,3) + ( fw*dh1 + EZ*dh2 )
          m = m + 1.0
        endif
      end do
      end do
      end do
!$OMP END DO

      
    case (Z_plus)
      
!$OMP DO SCHEDULE(static)
      do k=kx,kx
      do j=1,jx
      do i=1,ix
        if ( ibits(bv(i,j,k), bc_face_T, bitw_5) == obc_mask ) then
          Up0 = v0(i,j,k,1)
          Vp0 = v0(i,j,k,2)
          Wp0 = v0(i,j,k,3)
          
          Ut1 = u_bc
          Vt1 = v_bc
          Wt1 = w_bc
          c   = w_bc
          ac  = abs(c)
          fu  = 0.5*(c*(Ut1+Up0) - ac*(Ut1-Up0))
          fv  = 0.5*(c*(Vt1+Vp0) - ac*(Vt1-Vp0))
          fw  = 0.5*(c*(Wt1+Wp0) - ac*(Wt1-Wp0))
          
          EX = Ut1 - Up0
          EY = Vt1 - Vp0
          EZ = Wt1 - Wp0
          
          wv(i,j,k,1) = wv(i,j,k,1) + ( -fu*dh1 + EX*dh2 )
          wv(i,j,k,2) = wv(i,j,k,2) + ( -fv*dh1 + EY*dh2 )
          wv(i,j,k,3) = wv(i,j,k,3) + ( -fw*dh1 + EZ*dh2 )
          m = m + 1.0
        endif
      end do
      end do
      end do
!$OMP END DO
      
    case default

    end select FACES
    
!$OMP END PARALLEL

    flop = flop + dble(m)*34.0d0
    
    return
    end subroutine pvec_vobc_specv
    
!> ********************************************************************
!! @brief 外部速度境界条件による対流項と粘性項の流束の修正
!! @param [out] wv   疑似ベクトルの空間項の評価値
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル長
!! @param [in]  dh   格子幅
!! @param [in]  rei  Reynolds数の逆数
!! @param [in]  v0   速度ベクトル（n-step）
!! @param [in]  bv   BCindex V
!! @param [in]  face 外部境界処理のときの面番号
!! @param [out] flop 浮動小数点演算数
!! @note 境界面で対流流束はゼロ，粘性流束のみ
!<
    subroutine pvec_vobc_symtrc (wv, sz, g, dh, rei, v0, bv, face, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, g, face
    integer                                                   ::  ix, jx, kx
    integer, dimension(3)                                     ::  sz
    double precision                                          ::  flop, rix, rjx, rkx
    real                                                      ::  dh, dh1, dh2, rei
    real                                                      ::  Up0, Vp0, Wp0
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v0, wv
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bv
    
    ix = sz(1)
    jx = sz(2)
    kx = sz(3)

    rix = dble(jx)*dble(kx)
    rjx = dble(ix)*dble(kx)
    rkx = dble(ix)*dble(jx)
    
    dh1= 1.0/dh
    dh2= 2.0*rei*dh1*dh1
    
    flop = flop + 14.0d0 ! DP 19 flop

!$OMP PARALLEL REDUCTION(+:flop) &
!$OMP FIRSTPRIVATE(ix, jx, kx, dh2, face) &
!$OMP PRIVATE(i, j, k, Up0, Vp0, Wp0)
    
    FACES : select case (face)
    
    case (X_minus)
      i = 1
      
!$OMP DO SCHEDULE(static)
      do k=1,kx
      do j=1,jx
        if ( ibits(bv(i,j,k), bc_face_W, bitw_5) == obc_mask ) then
          Up0 = v0(i,j,k,1)
          wv(i,j,k,1) = wv(i,j,k,1) + Up0*dh2
        endif
      end do
      end do
!$OMP END DO

      flop = flop + rix*2.0d0
      

    case (X_plus)
      i = ix
      
!$OMP DO SCHEDULE(static)
      do k=1,kx
      do j=1,jx 
        if ( ibits(bv(i,j,k), bc_face_E, bitw_5) == obc_mask ) then
          Up0 = v0(i,j,k,1)
          wv(i,j,k,1) = wv(i,j,k,1) - Up0*dh2
        endif
      end do
      end do
!$OMP END DO
    
      flop = flop + rix*2.0d0
      

    case (Y_minus)
      j = 1
      
!$OMP DO SCHEDULE(static)
      do k=1,kx
      do i=1,ix
        if ( ibits(bv(i,j,k), bc_face_S, bitw_5) == obc_mask ) then
          Vp0 = v0(i,j,k,2)
          wv(i,j,k,2) = wv(i,j,k,2) + Vp0*dh2
        endif
      end do
      end do
!$OMP END DO
      
      flop = flop + rjx*2.0d0
      

    case (Y_plus)
      j = jx
      
!$OMP DO SCHEDULE(static)
      do k=1,kx
      do i=1,ix
        if ( ibits(bv(i,j,k), bc_face_N, bitw_5) == obc_mask ) then
          Vp0 = v0(i,j,k,2)
          wv(i,j,k,2) = wv(i,j,k,2) - Vp0*dh2
        endif
      end do
      end do
!$OMP END DO
      
      flop = flop + rjx*2.0d0
      

    case (Z_minus)
      k = 1
      
!$OMP DO SCHEDULE(static)
      do j=1,jx
      do i=1,ix
        if ( ibits(bv(i,j,k), bc_face_B, bitw_5) == obc_mask ) then
          Wp0 = v0(i,j,k,3)
          wv(i,j,k,3) = wv(i,j,k,3) + Wp0*dh2
        endif
      end do
      end do
!$OMP END DO
      
      flop = flop + rkx*2.0d0
      

    case (Z_plus)
      k = kx
      
!$OMP DO SCHEDULE(static)
      do j=1,jx
      do i=1,ix
        if ( ibits(bv(i,j,k), bc_face_T, bitw_5) == obc_mask ) then
          Wp0 = v0(i,j,k,3)
          wv(i,j,k,3) = wv(i,j,k,3) - Wp0*dh2
        endif
      end do
      end do
!$OMP END DO
      
      flop = flop + rkx*2.0d0
      

    case default
    end select FACES
    
!$OMP END PARALLEL

    return
    end subroutine pvec_vobc_symtrc
    
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
!<
    subroutine pvec_vobc_wall (wv, sz, g, dh, rei, v0, vec, face, flop)
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

    rix = dble(jx)*dble(kx)
    rjx = dble(ix)*dble(kx)
    rkx = dble(ix)*dble(jx)

    dh1= 1.0/dh
    dh2= rei*dh1*dh1*2.0
    
    u_bc = vec(1)
    v_bc = vec(2)
    w_bc = vec(3)
    
    flop = flop + 14.0d0


!$OMP PARALLEL REDUCTION(+:flop) &
!$OMP FIRSTPRIVATE(ix, jx, kx, u_bc, v_bc, w_bc, dh2, face, rix, rjx, rkx)

    FACES : select case (face)
    
    case (X_minus)
      
!$OMP DO SCHEDULE(static)
      do i=1,1
      do k=1,kx
      do j=1,jx
      ! wv(i,j,k,1) = wv(i,j,k,1) + (u_bc - v0(i,j,k,1)) * dh2 ! du/dx = 0
        wv(i,j,k,2) = wv(i,j,k,2) + (v_bc - v0(i,j,k,2)) * dh2 ! dv/dx
        wv(i,j,k,3) = wv(i,j,k,3) + (w_bc - v0(i,j,k,3)) * dh2 ! dw/dx
      end do
      end do
      end do
!$OMP END DO

      flop = flop + rix*9.0d0

      
    case (X_plus)

!$OMP DO SCHEDULE(static)
      do i=ix,ix
      do k=1,kx
      do j=1,jx
      ! wv(i,j,k,1) = wv(i,j,k,1) + (u_bc - v0(i,j,k,1)) * dh2 ! du/dx = 0
        wv(i,j,k,2) = wv(i,j,k,2) + (v_bc - v0(i,j,k,2)) * dh2 ! dv/dx
        wv(i,j,k,3) = wv(i,j,k,3) + (w_bc - v0(i,j,k,3)) * dh2 ! dw/dx
      end do
      end do
      end do
!$OMP END DO

      flop = flop + rix*9.0d0

      
    case (Y_minus)

!$OMP DO SCHEDULE(static)
      do j=1,1
      do k=1,kx
      do i=1,ix
        wv(i,j,k,1) = wv(i,j,k,1) + (u_bc - v0(i,j,k,1)) * dh2 ! du/dy
      ! wv(i,j,k,2) = wv(i,j,k,2) + (v_bc - v0(i,j,k,2)) * dh2 ! dv/dy = 0
        wv(i,j,k,3) = wv(i,j,k,3) + (w_bc - v0(i,j,k,3)) * dh2 ! dw/dy
      end do
      end do
      end do
!$OMP END DO

      flop = flop + rjx*6.0d0

      
    case (Y_plus)

!$OMP DO SCHEDULE(static)
      do j=jx,jx
      do k=1,kx
      do i=1,ix
        wv(i,j,k,1) = wv(i,j,k,1) + (u_bc - v0(i,j,k,1)) * dh2 ! du/dy
      ! wv(i,j,k,2) = wv(i,j,k,2) + (v_bc - v0(i,j,k,2)) * dh2 ! dv/dy = 0
        wv(i,j,k,3) = wv(i,j,k,3) + (w_bc - v0(i,j,k,3)) * dh2 ! dw/dy
      end do
      end do
      end do
!$OMP END DO

      flop = flop + rjx*6.0d0

      
    case (Z_minus)

!$OMP DO SCHEDULE(static)
      do k=1,1
      do j=1,jx
      do i=1,ix
        wv(i,j,k,1) = wv(i,j,k,1) + (u_bc - v0(i,j,k,1)) * dh2 ! du/dz
        wv(i,j,k,2) = wv(i,j,k,2) + (v_bc - v0(i,j,k,2)) * dh2 ! dv/dz
      ! wv(i,j,k,3) = wv(i,j,k,3) + (w_bc - v0(i,j,k,3)) * dh2 ! dw/dz = 0
      end do
      end do
      end do
!$OMP END DO

      flop = flop + rkx*6.0d0

      
    case (Z_plus)

!$OMP DO SCHEDULE(static)
      do k=kx,kx
      do j=1,jx
      do i=1,ix
        wv(i,j,k,1) = wv(i,j,k,1) + (u_bc - v0(i,j,k,1)) * dh2 ! du/dz
        wv(i,j,k,2) = wv(i,j,k,2) + (v_bc - v0(i,j,k,2)) * dh2 ! dv/dz
      ! wv(i,j,k,3) = wv(i,j,k,3) + (w_bc - v0(i,j,k,3)) * dh2 ! dw/dz = 0
      end do
      end do
      end do
!$OMP END DO

      flop = flop + rkx*6.0d0

    case default

    end select FACES
    
!$OMP END PARALLEL
    
    return
    end subroutine pvec_vobc_wall


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
!$OMP FIRSTPRIVATE(ix, jx, kx, u_bc, v_bc, w_bc, face)

    FACES : select case (face)
    
    case (X_minus)

!$OMP DO SCHEDULE(static)
      do i=1,1
      do k=1,kx
      do j=1,jx
        if ( ibits(bv(i,j,k), bc_face_W, bitw_5) == obc_mask ) then
          v(i-1, j, k, 1) = u_bc
          v(i-1, j, k, 2) = v_bc
          v(i-1, j, k, 3) = w_bc
        endif
      end do
      end do
      end do
!$OMP END DO
      
      
    case (X_plus)

!$OMP DO SCHEDULE(static)
      do i=ix,ix
      do k=1,kx
      do j=1,jx
        if ( ibits(bv(i,j,k), bc_face_E, bitw_5) == obc_mask ) then
          v(i+1, j, k, 1) = u_bc
          v(i+1, j, k, 2) = v_bc
          v(i+1, j, k, 3) = w_bc
        endif
      end do
      end do
      end do
!$OMP END DO
      
      
    case (Y_minus)

!$OMP DO SCHEDULE(static)
      do j=1,1
      do k=1,kx
      do i=1,ix
        if ( ibits(bv(i,j,k), bc_face_S, bitw_5) == obc_mask ) then
          v(i, j-1, k, 1) = u_bc
          v(i, j-1, k, 2) = v_bc
          v(i, j-1, k, 3) = w_bc
        endif
      end do
      end do
      end do
!$OMP END DO
      
      
    case (Y_plus)

!$OMP DO SCHEDULE(static)
      do j=jx,jx
      do k=1,kx
      do i=1,ix
        if ( ibits(bv(i,j,k), bc_face_N, bitw_5) == obc_mask ) then
          v(i, j+1, k, 1) = u_bc
          v(i, j+1, k, 2) = v_bc
          v(i, j+1, k, 3) = w_bc
        endif
      end do
      end do
      end do
!$OMP END DO
      
      
    case (Z_minus)

!$OMP DO SCHEDULE(static)
      do k=1,1
      do j=1,jx
      do i=1,ix
        if ( ibits(bv(i,j,k), bc_face_B, bitw_5) == obc_mask ) then
          v(i, j, k-1, 1) = u_bc
          v(i, j, k-1, 2) = v_bc
          v(i, j, k-1, 3) = w_bc
        endif
      end do
      end do
      end do
!$OMP END DO
      
      
    case (Z_plus)

!$OMP DO SCHEDULE(static)
      do k=kx,kx
      do j=1,jx
      do i=1,ix
        if ( ibits(bv(i,j,k), bc_face_T, bitw_5) == obc_mask ) then
          v(i, j, k+1, 1) = u_bc
          v(i, j, k+1, 2) = v_bc
          v(i, j, k+1, 3) = w_bc
        endif
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
!$OMP FIRSTPRIVATE(ix, jx, kx, face)

    FACES : select case (face)
    case (X_minus)

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do j=1,jx
        v(0, j, k, 1) = v(1, j, k, 1)
        v(0, j, k, 2) = v(1, j, k, 2)
        v(0, j, k, 3) = v(1, j, k, 3)
      end do
      end do
!$OMP END DO
      

    case (X_plus)

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do j=1,jx
        v(ix+1, j, k, 1) = v(ix, j, k, 1)
        v(ix+1, j, k, 2) = v(ix, j, k, 2)
        v(ix+1, j, k, 3) = v(ix, j, k, 3)
      end do
      end do
!$OMP END DO
      

    case (Y_minus)

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do i=1,ix
        v(i, 0, k, 1) = v(i, 1, k, 1)
        v(i, 0, k, 2) = v(i, 1, k, 2)
        v(i, 0, k, 3) = v(i, 1, k, 3)
      end do
      end do
!$OMP END DO
      

    case (Y_plus)

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do i=1,ix      
        v(i, jx+1, k, 1) = v(i, jx, k, 1)
        v(i, jx+1, k, 2) = v(i, jx, k, 2)
        v(i, jx+1, k, 3) = v(i, jx, k, 3)
      end do
      end do
!$OMP END DO
      

    case (Z_minus)

!$OMP DO SCHEDULE(static)
      do j=1,jx
      do i=1,ix
        v(i, j, 0, 1) = v(i, j, 1, 1)
        v(i, j, 0, 2) = v(i, j, 1, 2)
        v(i, j, 0, 3) = v(i, j, 1, 3)
      end do
      end do
!$OMP END DO
      

    case (Z_plus)

!$OMP DO SCHEDULE(static)
      do j=1,jx
      do i=1,ix
        v(i, j, kx+1, 1) = v(i, j, kx, 1)
        v(i, j, kx+1, 2) = v(i, j, kx, 2)
        v(i, j, kx+1, 3) = v(i, j, kx, 3)
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
      do i=1,1
      do k=1,kx
      do j=1,jx
        if ( ibits(bv(i,j,k), bc_face_W, bitw_5) == obc_mask ) then
          Uw = v0(i-1,j  ,k  , 1)
          Ue = v0(i  ,j  ,k  , 1)
          Vw = v0(i-1,j  ,k  , 2)
          Ve = v0(i  ,j  ,k  , 2)
          Ww = v0(i-1,j  ,k  , 3)
          We = v0(i  ,j  ,k  , 3)

          v(i-1, j  ,k  , 1) = Uw - c*(Ue-Uw)
          v(i-1, j  ,k  , 2) = Vw - c*(Ve-Vw)
          v(i-1, j  ,k  , 3) = Ww - c*(We-Ww)
          m = m + 1.0
        endif
      end do
      end do
      end do
!$OMP END DO
      
      
    case (X_plus)

      if ( cf<0.0 ) c=0.0
      
!$OMP DO SCHEDULE(static)
      do i=ix,ix
      do k=1,kx
      do j=1,jx
        if ( ibits(bv(i,j,k), bc_face_E, bitw_5) == obc_mask ) then
          Uw = v0(i  ,j  ,k  , 1)
          Ue = v0(i+1,j  ,k  , 1)
          Vw = v0(i  ,j  ,k  , 2)
          Ve = v0(i+1,j  ,k  , 2)
          Ww = v0(i  ,j  ,k  , 3)
          We = v0(i+1,j  ,k  , 3)

          v(i+1, j  ,k  , 1) = Ue - c*(Ue-Uw)
          v(i+1, j  ,k  , 2) = Ve - c*(Ve-Vw)
          v(i+1, j  ,k  , 3) = We - c*(We-Ww)
          m = m + 1.0
        endif
      end do
      end do
      end do
!$OMP END DO
      
      
    case (Y_minus)

      if ( cf>0.0 ) c=0.0

!$OMP DO SCHEDULE(static)
      do j=1,1
      do k=1,kx
      do i=1,ix
        if ( ibits(bv(i,j,k), bc_face_S, bitw_5) == obc_mask ) then
          Us = v0(i  ,j-1,k  , 1)
          Un = v0(i  ,j  ,k  , 1)
          Vs = v0(i  ,j-1,k  , 2)
          Vn = v0(i  ,j  ,k  , 2)
          Ws = v0(i  ,j-1,k  , 3)
          Wn = v0(i  ,j  ,k  , 3)

          v(i  ,j-1, k  , 1) = Us - c*(Un-Us)
          v(i  ,j-1, k  , 2) = Vs - c*(Vn-Vs)
          v(i  ,j-1, k  , 3) = Ws - c*(Wn-Ws)
          m = m + 1.0
        endif
      end do
      end do
      end do
!$OMP END DO
      
      
    case (Y_plus)

      if ( cf<0.0 ) c=0.0
      
!$OMP DO SCHEDULE(static)
      do j=jx,jx
      do k=1,kx
      do i=1,ix
        if ( ibits(bv(i,j,k), bc_face_N, bitw_5) == obc_mask ) then
          Us = v0(i  ,j  ,k  , 1)
          Un = v0(i  ,j+1,k  , 1)
          Vs = v0(i  ,j  ,k  , 2)
          Vn = v0(i  ,j+1,k  , 2)
          Ws = v0(i  ,j  ,k  , 3)
          Wn = v0(i  ,j+1,k  , 3)

          v(i  ,j+1, k  , 1) = Un - c*(Un-Us)
          v(i  ,j+1, k  , 2) = Vn - c*(Vn-Vs)
          v(i  ,j+1, k  , 3) = Wn - c*(Wn-Ws)
          m = m + 1.0
        endif
      end do
      end do
      end do
!$OMP END DO
      
      
    case (Z_minus)

      if ( cf>0.0 ) c=0.0

!$OMP DO SCHEDULE(static)
      do k=1,1
      do j=1,jx
      do i=1,ix
        if ( ibits(bv(i,j,k), bc_face_B, bitw_5) == obc_mask ) then
          Ub = v0(i  ,j  ,k-1, 1)
          Ut = v0(i  ,j  ,k  , 1)
          Vb = v0(i  ,j  ,k-1, 2)
          Vt = v0(i  ,j  ,k  , 2)
          Wb = v0(i  ,j  ,k-1, 3)
          Wt = v0(i  ,j  ,k  , 3)

          v(i  ,j  , k-1, 1) = Ub - c*(Ut-Ub)
          v(i  ,j  , k-1, 2) = Vb - c*(Vt-Vb)
          v(i  ,j  , k-1, 3) = Wb - c*(Wt-Wb)
          m = m + 1.0
        endif
      end do
      end do
      end do
!$OMP END DO
      
      
    case (Z_plus)

      if ( cf<0.0 ) c=0.0
      
!$OMP DO SCHEDULE(static)
      do k=kx,kx
      do j=1,jx
      do i=1,ix
        if ( ibits(bv(i,j,k), bc_face_T, bitw_5) == obc_mask ) then
          Ub = v0(i  ,j  ,k  , 1)
          Ut = v0(i  ,j  ,k+1, 1)
          Vb = v0(i  ,j  ,k  , 2)
          Vt = v0(i  ,j  ,k+1, 2)
          Wb = v0(i  ,j  ,k  , 3)
          Wt = v0(i  ,j  ,k+1, 3)

          v(i  ,j  , k+1, 1) = Ut - c*(Ut-Ub)
          v(i  ,j  , k+1, 2) = Vt - c*(Vt-Vb)
          v(i  ,j  , k+1, 3) = Wt - c*(Wt-Wb)
          m = m + 1.0
        endif
      end do
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
      do i=0,0
      do k=1,kx
      do j=1,jx
        v(i,j,k,1) = vc(i,j,k,1)
        v(i,j,k,2) = vc(i,j,k,2)
        v(i,j,k,3) = vc(i,j,k,3)
      end do
      end do
      end do
!$OMP END DO
      
    case (X_plus)

!$OMP DO SCHEDULE(static)
      do i=ix+1,ix+1
      do k=1,kx
      do j=1,jx
        v(i,j,k,1) = vc(i,j,k,1)
        v(i,j,k,2) = vc(i,j,k,2)
        v(i,j,k,3) = vc(i,j,k,3)
      end do
      end do
      end do
!$OMP END DO
      
    case (Y_minus)

!$OMP DO SCHEDULE(static)
      do j=0,0
      do k=1,kx
      do i=1,ix
        v(i,j,k,1) = vc(i,j,k,1)
        v(i,j,k,2) = vc(i,j,k,2)
        v(i,j,k,3) = vc(i,j,k,3)
      end do
      end do
      end do
!$OMP END DO
      
    case (Y_plus)

!$OMP DO SCHEDULE(static)
      do j=jx+1,jx+1
      do k=1,kx
      do i=1,ix
        v(i,j,k,1) = vc(i,j,k,1)
        v(i,j,k,2) = vc(i,j,k,2)
        v(i,j,k,3) = vc(i,j,k,3)
      end do
      end do
      end do
!$OMP END DO
      
    case (Z_minus)

!$OMP DO SCHEDULE(static)
      do k=0,0
      do j=1,jx
      do i=1,ix
        v(i,j,k,1) = vc(i,j,k,1)
        v(i,j,k,2) = vc(i,j,k,2)
        v(i,j,k,3) = vc(i,j,k,3)
      end do
      end do
      end do
!$OMP END DO
      
    case (Z_plus)

!$OMP DO SCHEDULE(static)
      do k=kx+1,kx+1
      do j=1,jx
      do i=1,ix
        v(i,j,k,1) = vc(i,j,k,1)
        v(i,j,k,2) = vc(i,j,k,2)
        v(i,j,k,3) = vc(i,j,k,3)
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
!! @param [in,out] div   速度の発散
!! @param [in]     sz    配列長
!! @param [in]     g     ガイドセル長
!! @param [in]     face  面番号
!! @param [in]     bv    BCindex V
!! @param [in]     vec   指定する速度ベクトル
!! @param [in,out] flop  flop count
!! @note 指定面でも固体部分は対象外とするのでループ中に判定あり
!<
    subroutine div_obc_drchlt (div, sz, g, face, bv, vec, flop)
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

    rix = dble(jx)*dble(kx)
    rjx = dble(ix)*dble(kx)
    rkx = dble(ix)*dble(jx)
    
    u_bc = vec(1)
    v_bc = vec(2)
    w_bc = vec(3)
    
    flop = flop + 3.0d0

!$OMP PARALLEL REDUCTION(+:flop) &
!$OMP FIRSTPRIVATE(ix, jx, kx, u_bc, v_bc, w_bc, face, rix, rjx, rkx) &
!$OMP PRIVATE(bvx)

    FACES : select case (face)
    case (X_minus)
      
!$OMP DO SCHEDULE(static)
      do i=1,1
      do k=1,kx
      do j=1,jx
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_W, bitw_5) == obc_mask ) then
          div(i,j,k) = div(i,j,k) - u_bc * real(ibits(bvx, State, 1))
        endif
      end do
      end do
      end do
!$OMP END DO

      flop = flop + rix*3.0d0 ! 2+ real*1
      
      
    case (X_plus)
      
!$OMP DO SCHEDULE(static)
      do i=ix,ix
      do k=1,kx
      do j=1,jx
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_E, bitw_5) == obc_mask ) then
          div(i,j,k) = div(i,j,k) + u_bc * real(ibits(bvx, State, 1))
        endif
      end do
      end do
      end do
!$OMP END DO
      
      flop = flop + rix*3.0d0 ! 2+ real*1
      
      
    case (Y_minus)

!$OMP DO SCHEDULE(static)
      do j=1,1
      do k=1,kx
      do i=1,ix
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_S, bitw_5) == obc_mask ) then
          div(i,j,k) = div(i,j,k) - v_bc * real(ibits(bvx, State, 1))
        endif
      end do
      end do
      end do
!$OMP END DO
      
      flop = flop + rjx*3.0d0 ! 2+ real*1
      
      
    case (Y_plus)

!$OMP DO SCHEDULE(static)
      do j=jx,jx
      do k=1,kx
      do i=1,ix
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_N, bitw_5) == obc_mask ) then
          div(i,j,k) = div(i,j,k) + v_bc * real(ibits(bvx, State, 1))
        endif
      end do
      end do
      end do
!$OMP END DO
      
      flop = flop + rjx*3.0d0 ! 2+ real*1
    
    
    case (Z_minus)

!$OMP DO SCHEDULE(static)
      do k=1,1
      do j=1,jx
      do i=1,ix
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_B, bitw_5) == obc_mask ) then
          div(i,j,k) = div(i,j,k) - w_bc * real(ibits(bvx, State, 1))
        endif
      end do
      end do
      end do
!$OMP END DO
      
      flop = flop + rkx*3.0d0 ! 2+ real*1
      
      
    case (Z_plus)

!$OMP DO SCHEDULE(static)
      do k=kx,kx
      do j=1,jx
      do i=1,ix
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_T, bitw_5) == obc_mask ) then
          div(i,j,k) = div(i,j,k) + w_bc * real(ibits(bvx, State, 1))
        endif
      end do
      end do
      end do
!$OMP END DO
      
      flop = flop + rkx*3.0d0 ! 2+ real*1
    
    case default

    end select FACES

!$OMP END PARALLEL

    return
    end subroutine div_obc_drchlt

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
    subroutine div_obc_oflow_pvec (div, sz, g, face, cf, bv, vf, flop)
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
!$OMP PRIVATE(b_w, b_e, b_s, b_n, b_b, b_t, b_p)


    FACES : select case (face)
    
    case (X_minus)

      if ( cf>0.0 ) c=0.0

!$OMP DO SCHEDULE(static)
      do i=1,1
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
      end do
!$OMP END DO

      
    case (X_plus)

      if ( cf<0.0 ) c=0.0

!$OMP DO SCHEDULE(static)
      do i=ix,ix
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
      end do
!$OMP END DO


    case (Y_minus)

      if ( cf>0.0 ) c=0.0

!$OMP DO SCHEDULE(static)
      do j=1,1
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
      end do
!$OMP END DO

      
    case (Y_plus)

      if ( cf<0.0 ) c=0.0

!$OMP DO SCHEDULE(static)
      do j=jx,jx
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
      end do
!$OMP END DO

    
    case (Z_minus)

      if ( cf>0.0 ) c=0.0

!$OMP DO SCHEDULE(static)
      do k=1,1
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
      end do
!$OMP END DO

      
    case (Z_plus)

      if ( cf<0.0 ) c=0.0

!$OMP DO SCHEDULE(static)
      do k=kx,kx
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
      end do
!$OMP END DO

    case default

    end select FACES

!$OMP END PARALLEL

    flop = flop + dble(m)*29.0d0

    return
    end subroutine div_obc_oflow_pvec

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
    subroutine div_obc_oflow_vec (div, sz, g, face, aa, vf, bv, flop)
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
!$OMP PRIVATE(dv)

    FACES : select case (face)
    
    case (X_minus)

!$OMP DO SCHEDULE(static)
      do i=1,1
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
      end do
!$OMP END DO
      
      
    case (X_plus)

!$OMP DO SCHEDULE(static)
      do i=ix,ix
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
      end do
!$OMP END DO
      
      
    case (Y_minus)

!$OMP DO SCHEDULE(static)
      do j=1,1
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
      end do
!$OMP END DO
      

    case (Y_plus)

!$OMP DO SCHEDULE(static)
      do j=jx,jx
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
      end do
!$OMP END DO
    
    
    case (Z_minus)

!$OMP DO SCHEDULE(static)
      do k=1,1
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
      end do
!$OMP END DO
      

    case (Z_plus)

!$OMP DO SCHEDULE(static)
      do k=kx,kx
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
    end subroutine div_obc_oflow_vec


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
    subroutine div_obc_vec (sz, g, face, sum, vf, bv, flop)
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

    rix = dble(jx)*dble(kx)
    rjx = dble(ix)*dble(kx)
    rkx = dble(ix)*dble(jx)

    a = 0.0   ! sum

!$OMP PARALLEL &
!$OMP REDUCTION(+:a) &
!$OMP REDUCTION(+:flop) &
!$OMP FIRSTPRIVATE(ix, jx, kx, rix, rjx, rkx, face)

    FACES : select case (face)

    case (X_minus)

!$OMP DO SCHEDULE(static)
      do i=1,1
      do k=1,kx
      do j=1,jx
        if ( ibits(bv(i,j,k), bc_face_W, bitw_5) == obc_mask ) then
          a = a + vf(i-1,j,k,1)
        end if
      end do
      end do
      end do
!$OMP END DO

      flop = flop + rix*1.0d0


    case (X_plus)

!$OMP DO SCHEDULE(static)
      do i=ix,ix
      do k=1,kx
      do j=1,jx
        if ( ibits(bv(i,j,k), bc_face_E, bitw_5) == obc_mask ) then
          a = a + vf(i,j,k,1)
        end if
      end do
      end do
      end do
!$OMP END DO

      flop = flop + rix*1.0d0


    case (Y_minus)

!$OMP DO SCHEDULE(static)
      do j=1,1
      do k=1,kx
      do i=1,ix
        if ( ibits(bv(i,j,k), bc_face_S, bitw_5) == obc_mask ) then
          a = a + vf(i,j-1,k,2)
        end if
      end do
      end do
      end do
!$OMP END DO

      flop = flop + rjx*1.0d0


    case (Y_plus)

!$OMP DO SCHEDULE(static)
      do j=jx,jx
      do k=1,kx
      do i=1,ix
        if ( ibits(bv(i,j,k), bc_face_N, bitw_5) == obc_mask ) then
          a = a + vf(i,j,k,2)
        end if
      end do
      end do
      end do
!$OMP END DO

      flop = flop + rjx*1.0d0


    case (Z_minus)

!$OMP DO SCHEDULE(static)
      do k=1,1
      do j=1,jx
      do i=1,ix
        if ( ibits(bv(i,j,k), bc_face_B, bitw_5) == obc_mask ) then
          a = a + vf(i,j,k-1,3)
        end if
      end do
      end do
      end do
!$OMP END DO

      flop = flop + rkx*1.0d0


    case (Z_plus)

!$OMP DO SCHEDULE(static)
      do k=kx,kx
      do j=1,jx
      do i=1,ix
        if ( ibits(bv(i,j,k), bc_face_T, bitw_5) == obc_mask ) then
          a = a + vf(i,j,k,3)
        end if
      end do
      end do
      end do
!$OMP END DO

      flop = flop + rkx*1.0d0

    case default

    end select FACES

!$OMP END PARALLEL

    sum = a

    return
    end subroutine div_obc_vec


!> ********************************************************************
!! @brief セルフェイスの値をセットする
!! @param [out] v    セルフェイス速度
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル長
!! @param [in]  bv   BCindex V
!! @param [in]  face 外部境界の面番号
!! @param [in]  vec  指定する速度ベクトル
!<
    subroutine vobc_drchlt_vf (v, sz, g, bv, face, vec)
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
!$OMP FIRSTPRIVATE(ix, jx, kx, u_bc, v_bc, w_bc, face)

    FACES : select case (face)

    case (X_minus)

!$OMP DO SCHEDULE(static)
      do i=1,1
      do k=1,kx
      do j=1,jx
        if ( ibits(bv(i,j,k), bc_face_W, bitw_5) == obc_mask ) then
          v(i-1, j, k, 1) = u_bc
        endif
      end do
      end do
      end do
!$OMP END DO


    case (X_plus)

!$OMP DO SCHEDULE(static)
      do i=ix,ix
      do k=1,kx
      do j=1,jx
        if ( ibits(bv(i,j,k), bc_face_E, bitw_5) == obc_mask ) then
          v(i, j, k, 1) = u_bc
        endif
      end do
      end do
      end do
!$OMP END DO


    case (Y_minus)

!$OMP DO SCHEDULE(static)
      do j=1,1
      do k=1,kx
      do i=1,ix
        if ( ibits(bv(i,j,k), bc_face_S, bitw_5) == obc_mask ) then
          v(i, j-1, k, 2) = v_bc
        endif
      end do
      end do
      end do
!$OMP END DO


    case (Y_plus)

!$OMP DO SCHEDULE(static)
      do j=jx,jx
      do k=1,kx
      do i=1,ix
        if ( ibits(bv(i,j,k), bc_face_N, bitw_5) == obc_mask ) then
          v(i, j, k, 2) = v_bc
        endif
      end do
      end do
      end do
!$OMP END DO


    case (Z_minus)

!$OMP DO SCHEDULE(static)
      do k=1,1
      do j=1,jx
      do i=1,ix
        if ( ibits(bv(i,j,k), bc_face_B, bitw_5) == obc_mask ) then
          v(i, j, k-1, 3) = w_bc
        endif
      end do
      end do
      end do
!$OMP END DO


    case (Z_plus)

!$OMP DO SCHEDULE(static)
      do k=kx,kx
      do j=1,jx
      do i=1,ix
        if ( ibits(bv(i,j,k), bc_face_T, bitw_5) == obc_mask ) then
          v(i, j, k, 3) = w_bc
        endif
      end do
      end do
      end do
!$OMP END DO

    case default
    end select FACES

!$OMP END PARALLEL

    return
    end subroutine vobc_drchlt_vf
