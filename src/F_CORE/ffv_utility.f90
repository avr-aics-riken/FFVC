!###################################################################################
!
! FFV-C
! Frontflow / violet Cartesian
!
!
! Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
! All rights reserved.
!
! Copyright (c) 2011-2015 Institute of Industrial Science, The University of Tokyo.
! All rights reserved.
!
! Copyright (c) 2012-2015 Advanced Institute for Computational Science, RIKEN.
! All rights reserved.
!
!###################################################################################

!> @file   ffv_utility.f90
!! @brief  Utilitiy functions
!! @author aics
!<

!> ********************************************************************
!! @brief 有効セルに対する発散の自乗和を計算
!! @param [out] ds   残差の自乗和
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル長
!! @param [in]  div  発散値のベース
!! @param [in]  bp   BCindex P
!! @param [out] flop flop count
!<
    subroutine norm_v_div_l2 (ds, sz, g, div, bp, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, g
    integer, dimension(3)                                     ::  sz
    double precision                                          ::  flop
    real                                                      ::  ds, r
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  div
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    ds = 0.0

    flop = flop + dble(ix)*dble(jx)*dble(kx)*3.0d0

!$OMP PARALLEL &
!$OMP PRIVATE(r) &
!$OMP REDUCTION(+:ds) &
!$OMP FIRSTPRIVATE(ix, jx, kx)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
    do k=1,kx
    do j=1,jx
    do i=1,ix
      r = div(i,j,k) * real(ibits(bp(i,j,k), vld_cnvg, 1)) ! 有効セルの場合 1.0
      ds = ds + r*r
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    return
    end subroutine norm_v_div_l2

!> ********************************************************************
!! @brief 発散量の最大値を計算する
!! @param [out] ds   最大値
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル長
!! @param [in]  div  発散値
!! @param [in]  bp   BCindex P
!! @param [out] flop flop count
!<
    subroutine norm_v_div_max (ds, sz, g, div, bp, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, g
    integer, dimension(3)                                     ::  sz
    double precision                                          ::  flop
    real                                                      ::  ds, r
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  div
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    ds = 0.0

    flop = flop + dble(ix)*dble(jx)*dble(kx)*2.0d0

!$OMP PARALLEL &
!$OMP REDUCTION(max:ds) &
!$OMP PRIVATE(r) &
!$OMP FIRSTPRIVATE(ix, jx, kx)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
    do k=1,kx
    do j=1,jx
    do i=1,ix
      r = div(i,j,k) * real(ibits(bp(i,j,k), vld_cnvg, 1)) ! 有効セルの場合 1.0
      ds = max(ds, abs(r) )
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    return
    end subroutine norm_v_div_max


!> ********************************************************************
!! @brief 速度成分の最大値を計算する
!! @param [out] v_max 最大値
!! @param [in]  sz    配列長
!! @param [in]  g     ガイドセル長
!! @param [in]  v00   参照速度
!! @param [in]  v     速度ベクトル
!! @param [out] flop  flop count
!<
    subroutine find_vmax (v_max, sz, g, v00, v, flop)
    implicit none
    integer                                                   ::  i, j, k, ix, jx, kx, g
    integer, dimension(3)                                     ::  sz
    double precision                                          ::  flop
    real                                                      ::  vm1, vm2, vm3, v_max, vx, vy, vz
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v
    real, dimension(0:3)                                      ::  v00

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    vm1 = 0.0
    vm2 = 0.0
    vm3 = 0.0
    vx = v00(1)
    vy = v00(2)
    vz = v00(3)
    flop = flop + dble(ix)*dble(jx)*dble(kx)*9.0d0 + 2.0d0

!$OMP PARALLEL &
!$OMP REDUCTION(max:vm1) &
!$OMP REDUCTION(max:vm2) &
!$OMP REDUCTION(max:vm3) &
!$OMP FIRSTPRIVATE(ix, jx, kx, vx, vy, vz)

!$OMP DO SCHEDULE(static) COLLAPSE(2)

    do k=1,kx
    do j=1,jx
    do i=1,ix
      vm1 = max(vm1, abs(v(i,j,k,1)-vx ) )
      vm2 = max(vm2, abs(v(i,j,k,2)-vy ) )
      vm3 = max(vm3, abs(v(i,j,k,3)-vz ) )
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL
    
    v_max = max(vm1, vm2, vm3) ! maxss %xmm0, %xmm1, x 2 times > 2 flop

    return
    end subroutine find_vmax


!> ********************************************************************
!! @brief 速度勾配テンソルの第２不変量の計算
!! @param [out] q    速度勾配テンソルの第２不変量
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル長
!! @param [in]  dh   格子幅
!! @param [in]  v    セルセンター速度ベクトル
!! @param [in]  bv   BCindex C
!! @param [in]  v00  参照速度
!! @param [out] flop 浮動小数点演算数
!<
    subroutine i2vgt (q, sz, g, dh, v, bv, v00, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, g, idx
    integer, dimension(3)                                     ::  sz
    double precision                                          ::  flop
    integer                                                   ::  b_e1, b_w1, b_n1, b_s1, b_t1, b_b1
    real                                                      ::  u_ref, v_ref, w_ref, actv, rx, ry, rz
    real                                                      ::  Up0, Ue1, Uw1, Us1, Un1, Ub1, Ut1
    real                                                      ::  Vp0, Ve1, Vw1, Vs1, Vn1, Vb1, Vt1
    real                                                      ::  Wp0, We1, Ww1, Ws1, Wn1, Wb1, Wt1
    real                                                      ::  w_e, w_w, w_n, w_s, w_t, w_b, uq, vq, wq
    real                                                      ::  d11, d22, d33, d12, d13, d23
    real                                                      ::  w12, w13, w23
    real                                                      ::  q11, q22, q33, q12, q13, q23
    real, dimension(3)                                        ::  dh
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  q
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bv
    real, dimension(0:3)                                      ::  v00

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    rx = 0.25 / dh(1)
    ry = 0.25 / dh(2)
    rz = 0.25 / dh(3)

    ! 参照座標系の速度
    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)

    flop = flop + dble(ix)*dble(jx)*dble(kx)*78.0d0 + 24.0d0


!$OMP PARALLEL &
!$OMP PRIVATE(idx, actv, uq, vq, wq) &
!$OMP PRIVATE(Up0, Ue1, Uw1, Us1, Un1, Ub1, Ut1) &
!$OMP PRIVATE(Vp0, Ve1, Vw1, Vs1, Vn1, Vb1, Vt1) &
!$OMP PRIVATE(Wp0, We1, Ww1, Ws1, Wn1, Wb1, Wt1) &
!$OMP PRIVATE(b_e1, b_w1, b_n1, b_s1, b_t1, b_b1) &
!$OMP PRIVATE(w_e, w_w, w_n, w_s, w_t, w_b) &
!$OMP PRIVATE(d11, d22, d33, d12, d13, d23) &
!$OMP PRIVATE(w12, w13, w23) &
!$OMP PRIVATE(q11, q22, q33, q12, q13, q23) &
!$OMP FIRSTPRIVATE(ix, jx, kx, u_ref, v_ref, w_ref, rx, ry, rz)

!$OMP DO SCHEDULE(static) COLLAPSE(2)

    do k=1,kx
    do j=1,jx
    do i=1,ix
      idx = bv(i,j,k)

      Ub1 = v(i  ,j  ,k-1, 1)
      Us1 = v(i  ,j-1,k  , 1)
      Uw1 = v(i-1,j  ,k  , 1)
      Up0 = v(i  ,j  ,k  , 1)
      Ue1 = v(i+1,j  ,k  , 1)
      Un1 = v(i  ,j+1,k  , 1)
      Ut1 = v(i  ,j  ,k+1, 1)

      Vb1 = v(i  ,j  ,k-1, 2)
      Vs1 = v(i  ,j-1,k  , 2)
      Vw1 = v(i-1,j  ,k  , 2)
      Vp0 = v(i  ,j  ,k  , 2)
      Ve1 = v(i+1,j  ,k  , 2)
      Vn1 = v(i  ,j+1,k  , 2)
      Vt1 = v(i  ,j  ,k+1, 2)

      Wb1 = v(i  ,j  ,k-1, 3)
      Ws1 = v(i  ,j-1,k  , 3)
      Ww1 = v(i-1,j  ,k  , 3)
      Wp0 = v(i  ,j  ,k  , 3)
      We1 = v(i+1,j  ,k  , 3)
      Wn1 = v(i  ,j+1,k  , 3)
      Wt1 = v(i  ,j  ,k+1, 3)
      
      ! セル状態 (0-solid / 1-fluid) > 1 flop
      actv= real(ibits(idx, State, 1))
      
      b_w1= ibits(bv(i-1,j  ,k  ), State, 1)
      b_e1= ibits(bv(i+1,j  ,k  ), State, 1)
      b_s1= ibits(bv(i  ,j-1,k  ), State, 1)
      b_n1= ibits(bv(i  ,j+1,k  ), State, 1)
      b_b1= ibits(bv(i  ,j  ,k-1), State, 1)
      b_t1= ibits(bv(i  ,j  ,k+1), State, 1)
      
      ! 隣接セルの状態で置換フラグをセット セル状態 (0-solid / 1-fluid) > 12 flop
      w_e = real(b_e1) * actv
      w_w = real(b_w1) * actv
      w_n = real(b_n1) * actv
      w_s = real(b_s1) * actv
      w_t = real(b_t1) * actv
      w_b = real(b_b1) * actv
      
      ! セルセンターからの壁面修正速度 > 6 flop
      uq = 2.0*u_ref - Up0
      vq = 2.0*v_ref - Vp0
      wq = 2.0*w_ref - Wp0
      
      ! 壁面の場合の参照速度の修正
      if ( b_e1 == 0 ) then
        Ue1 = uq
        Ve1 = vq
        We1 = wq
      endif
      
      if ( b_w1 == 0 ) then
        Uw1 = uq
        Vw1 = vq
        Ww1 = wq
      end if
      
      if ( b_n1 == 0 ) then
        Un1 = uq
        Vn1 = vq
        Wn1 = wq
      end if
      
      if ( b_s1 == 0 ) then
        Us1 = uq
        Vs1 = vq
        Ws1 = wq
      end if
      
      if ( b_t1 == 0 ) then
        Ut1 = uq
        Vt1 = vq
        Wt1 = wq
      end if
      
      if ( b_b1 == 0 ) then
        Ub1 = uq
        Vb1 = vq
        Wb1 = wq
      end if

      d11 = rx * ( 2.0*(Ue1-Uw1) )
      d22 = ry * ( 2.0*(Vn1-Vs1) )
      d33 = rz * ( 2.0*(Wt1-Wb1) )

      d12 = ry * (Un1-Us1) + rx * (Ve1-Vw1)
      d13 = rz * (Ut1-Ub1) + rx * (We1-Ww1)
      d23 = rz * (Vt1-Vb1) + ry * (Wn1-Ws1)
      
      w12 = ry * (Un1-Us1) - rx * (Ve1-Vw1)
      w13 = rz * (Ut1-Ub1) - rx * (We1-Ww1)
      w23 = rz * (Vt1-Vb1) - ry * (Wn1-Ws1)
      
      q11 = -d11*d11
      q22 = -d22*d22
      q33 = -d33*d33
      q12 = w12*w12 - d12*d12
      q13 = w13*w13 - d13*d13
      q23 = w23*w23 - d23*d23
      
      q(i,j,k) = 0.5*( q11+q22+q33 + 2.0*(q12+q13+q23) ) * actv
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    return
    end subroutine i2vgt

!> ********************************************************************
!! @brief 渦度ベクトルの計算
!! @param [out] rot  渦度ベクトル
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル長
!! @param [in]  dh   格子幅
!! @param [in]  v    セルセンター速度ベクトル
!! @param [in]  bv   BCindex C
!! @param [in]  v00  参照速度
!! @param [out] flop 浮動小数点演算数
!<
    subroutine rot_v (rot, sz, g, dh, v, bv, v00, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, g, idx
    integer, dimension(3)                                     ::  sz
    double precision                                          ::  flop
    integer                                                   ::  b_e1, b_w1, b_n1, b_s1, b_t1, b_b1
    real                                                      ::  u_ref, v_ref, w_ref, actv, rx, ry, rz
    real                                                      ::  Up0, Ue1, Uw1, Us1, Un1, Ub1, Ut1
    real                                                      ::  Vp0, Ve1, Vw1, Vs1, Vn1, Vb1, Vt1
    real                                                      ::  Wp0, We1, Ww1, Ws1, Wn1, Wb1, Wt1
    real                                                      ::  w_e, w_w, w_n, w_s, w_t, w_b, uq, vq, wq
    real                                                      ::  r1, r2, r3
    real, dimension(3)                                        ::  dh
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v, rot
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bv
    real, dimension(0:3)                                      ::  v00

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    rx = 0.5 / dh(1)
    ry = 0.5 / dh(2)
    rz = 0.5 / dh(3)
    
    ! 参照座標系の速度
    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)

    flop = flop + dble(ix)*dble(jx)*dble(kx)*37.0d0 + 24.0d0


!$OMP PARALLEL &
!$OMP PRIVATE(idx, actv, uq, vq, wq) &
!$OMP PRIVATE(Up0, Ue1, Uw1, Us1, Un1, Ub1, Ut1) &
!$OMP PRIVATE(Vp0, Ve1, Vw1, Vs1, Vn1, Vb1, Vt1) &
!$OMP PRIVATE(Wp0, We1, Ww1, Ws1, Wn1, Wb1, Wt1) &
!$OMP PRIVATE(b_e1, b_w1, b_n1, b_s1, b_t1, b_b1) &
!$OMP PRIVATE(w_e, w_w, w_n, w_s, w_t, w_b) &
!$OMP PRIVATE(r1, r2, r3) &
!$OMP FIRSTPRIVATE(ix, jx, kx, u_ref, v_ref, w_ref, rx, ry, rz)

!$OMP DO SCHEDULE(static) COLLAPSE(2)

    do k=1,kx
    do j=1,jx
    do i=1,ix
      idx = bv(i,j,k)

      Ub1 = v(i  ,j  ,k-1, 1)
      Us1 = v(i  ,j-1,k  , 1)
      Uw1 = v(i-1,j  ,k  , 1)
      Up0 = v(i  ,j  ,k  , 1)
      Ue1 = v(i+1,j  ,k  , 1)
      Un1 = v(i  ,j+1,k  , 1)
      Ut1 = v(i  ,j  ,k+1, 1)

      Vb1 = v(i  ,j  ,k-1, 2)
      Vs1 = v(i  ,j-1,k  , 2)
      Vw1 = v(i-1,j  ,k  , 2)
      Vp0 = v(i  ,j  ,k  , 2)
      Ve1 = v(i+1,j  ,k  , 2)
      Vn1 = v(i  ,j+1,k  , 2)
      Vt1 = v(i  ,j  ,k+1, 2)

      Wb1 = v(i  ,j  ,k-1, 3)
      Ws1 = v(i  ,j-1,k  , 3)
      Ww1 = v(i-1,j  ,k  , 3)
      Wp0 = v(i  ,j  ,k  , 3)
      We1 = v(i+1,j  ,k  , 3)
      Wn1 = v(i  ,j+1,k  , 3)
      Wt1 = v(i  ,j  ,k+1, 3)
      
      ! セル状態 (0-solid / 1-fluid) > 1 flop
      actv= real(ibits(idx, State, 1))
      
      b_w1= ibits(bv(i-1,j  ,k  ), State, 1)
      b_e1= ibits(bv(i+1,j  ,k  ), State, 1)
      b_s1= ibits(bv(i  ,j-1,k  ), State, 1)
      b_n1= ibits(bv(i  ,j+1,k  ), State, 1)
      b_b1= ibits(bv(i  ,j  ,k-1), State, 1)
      b_t1= ibits(bv(i  ,j  ,k+1), State, 1)
      
      ! 隣接セルの状態で置換フラグをセット セル状態 (0-solid / 1-fluid) > 12 flop
      w_e = real(b_e1) * actv
      w_w = real(b_w1) * actv
      w_n = real(b_n1) * actv
      w_s = real(b_s1) * actv
      w_t = real(b_t1) * actv
      w_b = real(b_b1) * actv
      
      ! セルセンターからの壁面修正速度 > 6 flop
      uq = 2.0*u_ref - Up0
      vq = 2.0*v_ref - Vp0
      wq = 2.0*w_ref - Wp0
      
      ! 壁面の場合の参照速度の修正
      if ( b_e1 == 0 ) then
        Ue1 = uq
        Ve1 = vq
        We1 = wq
      endif
      
      if ( b_w1 == 0 ) then
        Uw1 = uq
        Vw1 = vq
        Ww1 = wq
      end if
      
      if ( b_n1 == 0 ) then
        Un1 = uq
        Vn1 = vq
        Wn1 = wq
      end if
      
      if ( b_s1 == 0 ) then
        Us1 = uq
        Vs1 = vq
        Ws1 = wq
      end if
      
      if ( b_t1 == 0 ) then
        Ut1 = uq
        Vt1 = vq
        Wt1 = wq
      end if
      
      if ( b_b1 == 0 ) then
        Ub1 = uq
        Vb1 = vq
        Wb1 = wq
      end if
      
      r1 = ry * (Wn1-Ws1) - rz * (Vt1-Vb1) ! 12 flop
      r2 = rz * (Ut1-Ub1) - rx * (We1-Ww1)
      r3 = rx * (Ve1-Vw1) - ry * (Un1-Us1)
      
      rot(i,j,k,1) = r1 * actv ! 3 flop
      rot(i,j,k,2) = r2 * actv
      rot(i,j,k,3) = r3 * actv
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    return
    end subroutine rot_v

!> ********************************************************************
!! @brief ヘリシティの計算
!! @param [out] ht   ヘリシティ
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル長
!! @param [in]  dh   格子幅
!! @param [in]  v    セルセンター速度ベクトル
!! @param [in]  bv   BCindex C
!! @param [in]  v00  参照速度
!! @param [out] flop flop count
!<
    subroutine helicity (ht, sz, g, dh, v, bv, v00, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, g, idx
    integer, dimension(3)                                     ::  sz
    double precision                                          ::  flop
    integer                                                   ::  b_e1, b_w1, b_n1, b_s1, b_t1, b_b1
    real                                                      ::  u_ref, v_ref, w_ref, actv, rx, ry, rz
    real                                                      ::  Up0, Ue1, Uw1, Us1, Un1, Ub1, Ut1
    real                                                      ::  Vp0, Ve1, Vw1, Vs1, Vn1, Vb1, Vt1
    real                                                      ::  Wp0, We1, Ww1, Ws1, Wn1, Wb1, Wt1
    real                                                      ::  w_e, w_w, w_n, w_s, w_t, w_b, uq, vq, wq
    real                                                      ::  r1, r2, r3, u1, u2, u3
    real, dimension(3)                                        ::  dh
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  ht
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bv
    real, dimension(0:3)                                      ::  v00

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)

    rx = 0.5 / dh(1)
    ry = 0.5 / dh(2)
    rz = 0.5 / dh(3)

    ! 参照座標系の速度
    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)

    flop = flop + dble(ix)*dble(jx)*dble(kx)*43.0d0 + 24.0d0


!$OMP PARALLEL &
!$OMP PRIVATE(idx, actv, uq, vq, wq) &
!$OMP PRIVATE(Up0, Ue1, Uw1, Us1, Un1, Ub1, Ut1) &
!$OMP PRIVATE(Vp0, Ve1, Vw1, Vs1, Vn1, Vb1, Vt1) &
!$OMP PRIVATE(Wp0, We1, Ww1, Ws1, Wn1, Wb1, Wt1) &
!$OMP PRIVATE(b_e1, b_w1, b_n1, b_s1, b_t1, b_b1) &
!$OMP PRIVATE(w_e, w_w, w_n, w_s, w_t, w_b) &
!$OMP PRIVATE(r1, r2, r3, u1, u2, u3) &
!$OMP FIRSTPRIVATE(ix, jx, kx, u_ref, v_ref, w_ref, rx, ry, rz)

!$OMP DO SCHEDULE(static) COLLAPSE(2)

    do k=1,kx
    do j=1,jx
    do i=1,ix
      idx = bv(i,j,k)

      Ub1 = v(i  ,j  ,k-1, 1)
      Us1 = v(i  ,j-1,k  , 1)
      Uw1 = v(i-1,j  ,k  , 1)
      Up0 = v(i  ,j  ,k  , 1)
      Ue1 = v(i+1,j  ,k  , 1)
      Un1 = v(i  ,j+1,k  , 1)
      Ut1 = v(i  ,j  ,k+1, 1)

      Vb1 = v(i  ,j  ,k-1, 2)
      Vs1 = v(i  ,j-1,k  , 2)
      Vw1 = v(i-1,j  ,k  , 2)
      Vp0 = v(i  ,j  ,k  , 2)
      Ve1 = v(i+1,j  ,k  , 2)
      Vn1 = v(i  ,j+1,k  , 2)
      Vt1 = v(i  ,j  ,k+1, 2)

      Wb1 = v(i  ,j  ,k-1, 3)
      Ws1 = v(i  ,j-1,k  , 3)
      Ww1 = v(i-1,j  ,k  , 3)
      Wp0 = v(i  ,j  ,k  , 3)
      We1 = v(i+1,j  ,k  , 3)
      Wn1 = v(i  ,j+1,k  , 3)
      Wt1 = v(i  ,j  ,k+1, 3)

      ! セル状態 (0-solid / 1-fluid) > 1 flop
      actv= real(ibits(idx, State, 1))

      b_w1= ibits(bv(i-1,j  ,k  ), State, 1)
      b_e1= ibits(bv(i+1,j  ,k  ), State, 1)
      b_s1= ibits(bv(i  ,j-1,k  ), State, 1)
      b_n1= ibits(bv(i  ,j+1,k  ), State, 1)
      b_b1= ibits(bv(i  ,j  ,k-1), State, 1)
      b_t1= ibits(bv(i  ,j  ,k+1), State, 1)
      
      ! 隣接セルの状態で置換フラグをセット セル状態 (0-solid / 1-fluid) > 12 flop
      w_e = real(b_e1) * actv
      w_w = real(b_w1) * actv
      w_n = real(b_n1) * actv
      w_s = real(b_s1) * actv
      w_t = real(b_t1) * actv
      w_b = real(b_b1) * actv
      
      ! セルセンターからの壁面修正速度 > 6 flop
      uq = 2.0*u_ref - Up0
      vq = 2.0*v_ref - Vp0
      wq = 2.0*w_ref - Wp0
      
      ! 壁面の場合の参照速度の修正
      if ( b_e1 == 0 ) then
        Ue1 = uq
        Ve1 = vq
        We1 = wq
      endif
      
      if ( b_w1 == 0 ) then
        Uw1 = uq
        Vw1 = vq
        Ww1 = wq
      end if
      
      if ( b_n1 == 0 ) then
        Un1 = uq
        Vn1 = vq
        Wn1 = wq
      end if
      
      if ( b_s1 == 0 ) then
        Us1 = uq
        Vs1 = vq
        Ws1 = wq
      end if
      
      if ( b_t1 == 0 ) then
        Ut1 = uq
        Vt1 = vq
        Wt1 = wq
      end if
      
      if ( b_b1 == 0 ) then
        Ub1 = uq
        Vb1 = vq
        Wb1 = wq
      end if
      
      r1 = ry * (Wn1-Ws1) - rz * (Vt1-Vb1) ! 12 flop
      r2 = rz * (Ut1-Ub1) - rx * (We1-Ww1)
      r3 = rx * (Ve1-Vw1) - ry * (Un1-Us1)
      
      u1 = Up0 - u_ref ! 3 flop
      u2 = Vp0 - v_ref
      u3 = Wp0 - w_ref
      
      ht(i,j,k) = (u1*r1 + u2*r2 + u3*r3) * actv ! 6 flop
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    return
    end subroutine helicity

!> ********************************************************************
!! @brief 指定面の平均圧力を求める
!! @param p 圧力
!! @param sz 配列長
!! @param g ガイドセル長
!! @param face 外部境界面の番号
!! @param avr 平均値
!<
    subroutine face_avr_sampling (p, sz, g, face, avr)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, face, g
    integer, dimension(3)                                     ::  sz
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  p
    real                                                      ::  avr, rix, rjx, rkx

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    
    rix = real(jx)*real(kx)
    rjx = real(ix)*real(kx)
    rkx = real(ix)*real(jx)
    
    avr = 0.0

!$OMP PARALLEL REDUCTION(+:avr) &
!$OMP FIRSTPRIVATE(ix, jx, kx, face, rix, rjx, rkx)

    FACES : select case (face)
    
    case (X_minus)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
      do k=1,kx
      do j=1,jx
        avr = avr + p(1,j,k)
      end do
      end do
!$OMP END DO

      avr = avr / rix
      
      
    case (X_plus)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
      do k=1,kx
      do j=1,jx
        avr = avr + p(ix,j,k)
      end do
      end do
!$OMP END DO

      avr = avr / rix
      
      
    case (Y_minus)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
      do k=1,kx
      do i=1,ix
        avr = avr + p(i,1,k)
      end do
      end do
!$OMP END DO

      avr = avr / rjx
      
      
    case (Y_plus)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
      do k=1,kx
      do i=1,ix
        avr = avr + p(i,jx,k)
      end do
      end do
!$OMP END DO

      avr = avr / rjx
      
      
    case (Z_minus)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
      do j=1,jx
      do i=1,ix
        avr = avr + p(i,j,1)
      end do
      end do
!$OMP END DO

      avr = avr / rkx
      
    
    case (Z_plus)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
      do j=1,jx
      do i=1,ix
        avr = avr + p(i,j,kx)
      end do
      end do
!$OMP END DO
      
      avr = avr / rkx
      
    case default
    end select FACES

!$OMP END PARALLEL

    return
    end subroutine face_avr_sampling

!> ********************************************************************
!! @brief 圧力を指定値だけシフトする
!! @param p 圧力
!! @param sz 配列長
!! @param g ガイドセル長
!! @param avr 平均値
!<
    subroutine shift_pressure (p, sz, g, avr)
    implicit none
    integer                                                   ::  i, j, k, ix, jx, kx, g
    integer, dimension(3)                                     ::  sz
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  p
    real                                                      ::  avr

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, avr)

!$OMP DO SCHEDULE(static) COLLAPSE(2)

    do k=1,kx
    do j=1,jx
    do i=1,ix
      p(i,j,k) = p(i,j,k) - avr
    end do
    end do
    end do
!$OMP END DO

!$OMP END PARALLEL

    return
    end subroutine shift_pressure

!> ********************************************************************
!! @brief 物体表面の力の成分を計算する
!! @param [out] frc  力の成分
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル長
!! @param [in]  tgt  コンポーネントのエントリ番号
!! @param [in]  p    圧力
!! @param [in]  bid  交点id配列
!! @param [in]  dh   無次元格子幅
!! @param [in]  st   開始インデクス
!! @param [in]  ed   終了インデクス
!! @param [out] flop flop count
!!
!!
!!  力の符号は、軸の正方向の力をプラスの符号とする
!<
  subroutine force_compo (frc, sz, g, tgt, p, bid, dh, st, ed, flop)
  implicit none
  include 'ffv_f_params.h'
  integer                                                   ::  i, j, k, g, bd, is, ie, js, je, ks, ke
  integer, dimension(3)                                     ::  sz, st, ed
  double precision                                          ::  flop
  integer                                                   ::  idw, ide, ids, idn, idb, idt, tgt, tg
  real                                                      ::  fx, fy, fz, pp
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  p
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bid
  real, dimension(3)                                        ::  frc, dh

  is = st(1)
  js = st(2)
  ks = st(3)
  ie = ed(1)
  je = ed(2)
  ke = ed(3)

  tg = tgt

  fx = 0.0
  fy = 0.0
  fz = 0.0

  flop = flop + dble(ie-is+1)*dble(je-js+1)*dble(ke-ks+1)*7.0d0 + 5.0d0

!$OMP PARALLEL &
!$OMP REDUCTION(+:fx) &
!$OMP REDUCTION(+:fy) &
!$OMP REDUCTION(+:fz) &
!$OMP FIRSTPRIVATE(is, ie, js, je, ks, ke, tg) &
!$OMP PRIVATE(pp, bd) &
!$OMP PRIVATE(idw, ide, ids, idn, idb, idt)

!$OMP DO SCHEDULE(static) COLLAPSE(3)

  do k=ks,ke
  do j=js,je
  do i=is,ie

    bd = bid(i,j,k)

    ! 検査方向面のIDを取り出す >> CompoListのエントリ番号
    idw = ibits(bd, bc_face_W, bitw_5)
    ide = ibits(bd, bc_face_E, bitw_5)
    ids = ibits(bd, bc_face_S, bitw_5)
    idn = ibits(bd, bc_face_N, bitw_5)
    idb = ibits(bd, bc_face_B, bitw_5)
    idt = ibits(bd, bc_face_T, bitw_5)

    ! セルマスク　Fluid -> 1.0 / Wall -> 0.0
    pp = p(i,j,k) * real(ibits(bd, State, 1))

    ! 力の積算、軸方向を正にとる
    if ( idw == tg ) fx = fx - pp
    if ( ide == tg ) fx = fx + pp
    if ( ids == tg ) fy = fy - pp
    if ( idn == tg ) fy = fy + pp
    if ( idb == tg ) fz = fz - pp
    if ( idt == tg ) fz = fz + pp

  end do
  end do
  end do

!$OMP END DO
!$OMP END PARALLEL

  frc(1) = fx * dh(2)*dh(3)
  frc(2) = fy * dh(1)*dh(3)
  frc(3) = fz * dh(1)*dh(2)

  return
  end subroutine force_compo



!> ********************************************************************
!! @brief 乱流の統計情報
!! @param [out]    rms     瞬間の変動速度
!! @param [in,out] rmsmean 変動速度の二乗和の平方根（実効値）
!! @param [in]     sz      配列長
!! @param [in]     g       ガイドセル長
!! @param [in]     v       速度
!! @param [in]     av      速度の時間平均値
!! @param [in]     accum   積算ステップ数
!! @param [out]    flop    flop count
!!
!!  u          ; 瞬時値
!!  \var{u}    ; 時間平均値
!!  u^{\prime} ; 変動値 = u - \var{u}
!!  \sigma     ; 標準偏差 RMS, \sigma = \sqrt{ \var{ {u^{\prime}}^2 } }
!! Turbulent Intensity = \frac{\sigma}{\var{U}}
!<
subroutine calc_rms_v(rms, rmsmean, sz, g, v, av, accum, flop)
implicit none
include 'ffv_f_params.h'
integer                                                   :: ix, jx, kx, i, j, k, g
double precision                                          :: flop
integer, dimension(3)                                     :: sz
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) :: v, av, rms, rmsmean
real                                                      :: val1, val2, u1, u2, u3, accum

ix = sz(1)
jx = sz(2)
kx = sz(3)

val2 = 1.0/accum
val1 = 1.0 - val2

flop = flop + dble(ix+2*g) * dble(jx+2*g) * dble(kx+2*g) * 45.0d0 + 9.0d0

!$OMP PARALLEL &
!$OMP PRIVATE(u1, u2, u3) &
!$OMP FIRSTPRIVATE(ix, jx, kx, g, val1, val2)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k = 1-g, kx+g
do j = 1-g, jx+g
do i = 1-g, ix+g
u1 = v(i, j, k, 1) - av(i, j, k, 1)
u2 = v(i, j, k, 2) - av(i, j, k, 2)
u3 = v(i, j, k, 3) - av(i, j, k, 3)
u1 = sqrt( u1*u1 )
u2 = sqrt( u2*u2 )
u3 = sqrt( u3*u3 )

! 瞬間の標準偏差
rms(i, j, k, 1) = u1
rms(i, j, k, 2) = u2
rms(i, j, k, 3) = u3

! 時間平均 >> RMS
rmsmean(i, j, k, 1) = val1 * rmsmean(i, j, k, 1) + val2 * u1
rmsmean(i, j, k, 2) = val1 * rmsmean(i, j, k, 2) + val2 * u2
rmsmean(i, j, k, 3) = val1 * rmsmean(i, j, k, 3) + val2 * u3
end do
end do
end do
!$OMP END PARALLEL

return
end subroutine calc_rms_v


!> ********************************************************************
!! @brief スカラ変数の変動統計情報
!! @param [out]    rms     瞬間の変動値
!! @param [in,out] rmsmean 変動値の二乗和の平方根（実効値）
!! @param [in]     sz      配列長
!! @param [in]     g       ガイドセル長
!! @param [in]     s       スカラ変数
!! @param [in]     as      スカラ変数の時間平均値
!! @param [in]     accum   積算ステップ数
!! @param [out]    flop    flop count
!!
!!  u          ; 瞬時値
!!  \var{u}    ; 時間平均値
!!  u^{\prime} ; 変動値 = u - \var{u}
!!  \sigma     ; 標準偏差 RMS, \sigma = \sqrt{ \var{ {u^{\prime}}^2 } }
!! Turbulent Intensity = \frac{\sigma}{\var{U}}
!<
subroutine calc_rms_s(rms, rmsmean, sz, g, s, as, accum, flop)
implicit none
include 'ffv_f_params.h'
integer                                                   :: ix, jx, kx, i, j, k, g
double precision                                          :: flop
integer, dimension(3)                                     :: sz
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    :: s, as, rms, rmsmean
real                                                      :: val1, val2, u, u2, accum

ix = sz(1)
jx = sz(2)
kx = sz(3)

val2 = 1.0/accum
val1 = 1.0 - val2

flop = flop + dble(ix+2*g) * dble(jx+2*g) * dble(kx+2*g) * 15.0d0 + 9.0d0

!$OMP PARALLEL &
!$OMP PRIVATE(u, u2) &
!$OMP FIRSTPRIVATE(ix, jx, kx, g, val1, val2)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k = 1-g, kx+g
do j = 1-g, jx+g
do i = 1-g, ix+g
u = s(i, j, k) - as(i, j, k)
u2 = sqrt( u*u )

! 瞬間の標準偏差
rms(i, j, k) = u2

! 時間平均 >> RMS
rmsmean(i, j, k) = val1 * rmsmean(i, j, k) + val2 * u2
end do
end do
end do
!$OMP END PARALLEL

return
end subroutine calc_rms_s


!> ********************************************************************
!! @brief スカラ変数のガイドセルを調整してパックする
!! @param [out]  iblk  IBLANK
!! @param [in]   sz    配列長
!! @param [in]   g     アロケート時のガイドセル長
!! @param [in]   bcd   BCibdex D
!! @param [in]   w     出力するガイドセル数
!<
subroutine generate_iblank (iblk, sz, g, bcd, w)
implicit none
include 'ffv_f_params.h'
integer                                                   ::  i, j, k, ix, jx, kx, g, w, bd
integer, dimension(3)                                     ::  sz
integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bcd
integer, dimension(1-w:sz(1)+w, 1-w:sz(2)+w, 1-w:sz(3)+w) ::  iblk

ix = sz(1)
jx = sz(2)
kx = sz(3)

!$OMP PARALLEL &
!$OMP PRIVATE(i, j, k) &
!$OMP FIRSTPRIVATE(ix, jx, kx, w)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1-w, kx+w
do j=1-w, jx+w
do i=1-w, ix+w
  iblk(i,j,k) = 0
end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL


!$OMP PARALLEL &
!$OMP PRIVATE(i, j, k, bd) &
!$OMP FIRSTPRIVATE(ix, jx, kx, w)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1-w, kx+w
do j=1-w, jx+w
do i=1-w, ix+w
  bd = bcd(i,j,k)
  if ( ibits(bd, Active, 1) == 1 ) then

    if ( ibits(bd, State, 1) == 1 ) then
       iblk(i,j,k) = 1  ! Fluid >> 1
    else
       iblk(i,j,k) = 2  ! Solid >> 2
    end if

  else
    iblk(i,j,k) = 0  ! Inactive >> 0
  end if
end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL

return
end subroutine generate_iblank


!> ********************************************************************
!! @brief スカラ変数のガイドセルを調整してパックする
!! @param [out]  dst   出力先
!! @param [in]   sz    配列長
!! @param [in]   g     アロケート時のガイドセル長
!! @param [in]   src   入力
!! @param [in]   w     出力するガイドセル数
!<
subroutine pack_scalar (dst, sz, g, src, w)
implicit none
integer                                                   ::  i, j, k, ix, jx, kx, g, w
integer, dimension(3)                                     ::  sz
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  src
real, dimension(1-w:sz(1)+w, 1-w:sz(2)+w, 1-w:sz(3)+w)    ::  dst

ix = sz(1)
jx = sz(2)
kx = sz(3)

!$OMP PARALLEL &
!$OMP PRIVATE(i, j, k) &
!$OMP FIRSTPRIVATE(ix, jx, kx, w)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1-w, kx+w
do j=1-w, jx+w
do i=1-w, ix+w
  dst(i,j,k) = src(i,j,k)
end do
end do
end do

!$OMP END DO
!$OMP END PARALLEL

return
end subroutine pack_scalar


!> ********************************************************************
!! @brief ベクタ変数のガイドセルを調整してパックする
!! @param [out]  dst   出力先
!! @param [in]   sz    配列長
!! @param [in]   g     アロケート時のガイドセル長
!! @param [in]   src   入力
!! @param [in]   w     出力するガイドセル数
!<
subroutine pack_vector (dst, sz, g, src, w)
implicit none
integer                                                   ::  i, j, k, ix, jx, kx, g, w
integer, dimension(3)                                     ::  sz
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  src
real, dimension(1-w:sz(1)+w, 1-w:sz(2)+w, 1-w:sz(3)+w, 3) ::  dst

ix = sz(1)
jx = sz(2)
kx = sz(3)

!$OMP PARALLEL &
!$OMP PRIVATE(i, j, k) &
!$OMP FIRSTPRIVATE(ix, jx, kx, w)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1-w, kx+w
do j=1-w, jx+w
do i=1-w, ix+w
dst(i,j,k,1) = src(i,j,k,1)
dst(i,j,k,2) = src(i,j,k,2)
dst(i,j,k,3) = src(i,j,k,3)
end do
end do
end do

!$OMP END DO
!$OMP END PARALLEL

return
end subroutine pack_vector


!> ********************************************************************
!! @brief スカラ変数のガイドセルを調整してアンパックする
!! @param [out]  dst   出力先
!! @param [in]   sz    配列長
!! @param [in]   g     アロケート時のガイドセル長
!! @param [in]   src   入力
!! @param [in]   w     リスタートファイルに記録されたガイドセル数
!<
subroutine unpack_scalar (dst, sz, g, src, w)
implicit none
integer                                                   ::  i, j, k, ix, jx, kx, g, w
integer, dimension(3)                                     ::  sz
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  dst
real, dimension(1-w:sz(1)+w, 1-w:sz(2)+w, 1-w:sz(3)+w)    ::  src

ix = sz(1)
jx = sz(2)
kx = sz(3)

!$OMP PARALLEL &
!$OMP PRIVATE(i, j, k) &
!$OMP FIRSTPRIVATE(ix, jx, kx, w)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1-w, kx+w
do j=1-w, jx+w
do i=1-w, ix+w
  dst(i,j,k) = src(i,j,k)
end do
end do
end do

!$OMP END DO
!$OMP END PARALLEL

return
end subroutine unpack_scalar


!> ********************************************************************
!! @brief ベクタ変数のガイドセルを調整してアンパックする
!! @param [out]  dst   出力先
!! @param [in]   sz    配列長
!! @param [in]   g     アロケート時のガイドセル長
!! @param [in]   src   入力
!! @param [in]   w     リスタートファイルに記録されたガイドセル数
!<
subroutine unpack_vector (dst, sz, g, src, w)
implicit none
integer                                                   ::  i, j, k, ix, jx, kx, g, w
integer, dimension(3)                                     ::  sz
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  dst
real, dimension(1-w:sz(1)+w, 1-w:sz(2)+w, 1-w:sz(3)+w, 3) ::  src

ix = sz(1)
jx = sz(2)
kx = sz(3)

!$OMP PARALLEL &
!$OMP PRIVATE(i, j, k) &
!$OMP FIRSTPRIVATE(ix, jx, kx, w)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1-w, kx+w
do j=1-w, jx+w
do i=1-w, ix+w
dst(i,j,k,1) = src(i,j,k,1)
dst(i,j,k,2) = src(i,j,k,2)
dst(i,j,k,3) = src(i,j,k,3)
end do
end do
end do

!$OMP END DO
!$OMP END PARALLEL

return
end subroutine unpack_vector



!> ********************************************************************
!! @brief PLOT3D function format出力
!! @param [out]  buf   出力バッファ
!! @param [in]   sz    配列長
!! @param [in]   g     出力ガイドセル長
!! @param [in]   nvar  ブロック数
!! @param [in]   fname filename
!<
subroutine read_plot3d (buf, sz, g, nvar, fname)
implicit none
integer                                                  ::  i, j, k, ix, jx, kx, g, n, nx, nvar
integer, dimension(3)                                    ::  sz
real, dimension(sz(1)+2*g, sz(2)+2*g, sz(3)+2*g, nvar)   ::  buf
character*256                                            ::  fname

open(unit=30, file=fname, form='unformatted')
read(30) ix, jx, kx, nx
if ( (ix /= sz(1)+2*g) .or. (jx /= sz(2)+2*g) .or. (kx /= sz(3)+2*g) .or. (nx /= nvar) ) then
write(*,*) 'Read header error'
stop
end if
read(30) ((((buf(i,j,k,n),i=1,ix),j=1,jx),k=1,kx),n=1,nx)
close(unit=30)

return
end subroutine read_plot3d


!> ********************************************************************
!! @brief PLOT3D function format出力
!! @param [out]  buf   出力バッファ
!! @param [in]   sz    配列長
!! @param [in]   g     出力ガイドセル長
!! @param [in]   nx    ブロック数
!! @param [in]   fname filename
!<
subroutine write_plot3d (buf, sz, g, nx, fname)
implicit none
integer                                                  ::  i, j, k, ix, jx, kx, g, n, nx
integer, dimension(3)                                    ::  sz
real, dimension(sz(1)+2*g, sz(2)+2*g, sz(3)+2*g, nx)     ::  buf
character*256                                            ::  fname

ix = sz(1)+2*g
jx = sz(2)+2*g
kx = sz(3)+2*g

open(unit=30, file=fname, form='unformatted')
write(30) ix, jx, kx, nx
write(30) ((((buf(i,j,k,n),i=1,ix),j=1,jx),k=1,kx),n=1,nx)
close(unit=30)

return
end subroutine write_plot3d


!********************************************************************
! [in]  v         セルセンター速度ベクトル（n-step）
! [in]  vmean     セルセンター時間平均速度ベクトル（n-step）
subroutine output_vtk(step, G_origin, G_division, G_size, myRank, sz, dh, g, v, p)
implicit none
integer :: step
real, dimension(3) ::  G_origin
integer, dimension(3) :: G_division, G_size, sz
integer :: g
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  p
integer :: ix, jx, kx, i, j, k, myRank, ip, jp, kp, ii, jj, kk
real,dimension(3) :: dh

character fname*128

ix = sz(1)
jx = sz(2)
kx = sz(3)

if ( step < 10 ) then
if ( myRank == 0 ) then
write(fname, "('u_p0_step', i1.0, '.vtk')"), step
open(10, file = fname)
else if ( myRank < 10 ) then
write(fname, "('u_p', i1.0, '_step', i1.0, '.vtk')"), myRank, step
open(10, file = fname)
else if ( myRank < 100 ) then
write(fname, "('u_p', i2.0, '_step', i1.0, '.vtk')"), myRank, step
open(10, file = fname)
end if
else if ( step < 100 ) then
if ( myRank == 0 ) then
write(fname, "('u_p0_step', i2.0, '.vtk')"), step
open(10, file = fname)
else if ( myRank < 10 ) then
write(fname, "('u_p', i1.0, '_step', i2.0, '.vtk')"), myRank, step
open(10, file = fname)
else if ( myRank < 100 ) then
write(fname, "('u_p', i2.0, '_step', i2.0, '.vtk')"), myRank, step
open(10, file = fname)
end if
else if ( step < 1000 ) then
if ( myRank == 0 ) then
write(fname, "('u_p0_step', i3.0, '.vtk')"), step
open(10, file = fname)
else if ( myRank < 10 ) then
write(fname, "('u_p', i1.0, '_step', i3.0, '.vtk')"), myRank, step
open(10, file = fname)
else if ( myRank < 100 ) then
write(fname, "('u_p', i2.0, '_step', i3.0, '.vtk')"), myRank, step
open(10, file = fname)
end if
else if ( step < 10000 ) then
if ( myRank == 0 ) then
write(fname, "('u_p0_step', i4.0, '.vtk')"), step
open(10, file = fname)
else if ( myRank < 10 ) then
write(fname, "('u_p', i1.0, '_step', i4.0, '.vtk')"), myRank, step
open(10, file = fname)
else if ( myRank < 100 ) then
write(fname, "('u_p', i2.0, '_step', i4.0, '.vtk')"), myRank, step
open(10, file = fname)
end if
else if ( step < 100000 ) then
if ( myRank == 0 ) then
write(fname, "('u_p0_step', i5.0, '.vtk')"), step
open(10, file = fname)
else if ( myRank < 10 ) then
write(fname, "('u_p', i1.0, '_step', i5.0, '.vtk')"), myRank, step
open(10, file = fname)
else if ( myRank < 100 ) then
write(fname, "('u_p', i2.0, '_step', i5.0, '.vtk')"), myRank, step
open(10, file = fname)
end if
else if ( step < 1000000 ) then
if ( myRank == 0 ) then
write(fname, "('u_p0_step', i6.0, '.vtk')"), step
open(10, file = fname)
else if ( myRank < 10 ) then
write(fname, "('u_p', i1.0, '_step', i6.0, '.vtk')"), myRank, step
open(10, file = fname)
else if ( myRank < 100 ) then
write(fname, "('u_p', i2.0, '_step', i6.0, '.vtk')"), myRank, step
open(10, file = fname)
end if
end if



!----- writing vtk data ---
write(10, "('# vtk DataFile Version 3.0')")
write(10, "('v_vmean_p.vtk')")
write(10, "('ASCII')")

write(10, "('DATASET STRUCTURED_GRID')")
write(10, "('DIMENSIONS ', 3(1x, i8))") ix+2*g, jx+2*g, kx+2*g
write(10, "('POINTS ', i20, ' float')") (ix+2*g)*(jx+2*g)*(kx+2*g)

do kp = 1-g, kx+g
do jp = 1-g, jx+g
do ip = 1-g, ix+g
kk = int(myRank/(G_division(1)*G_division(2)))
!jj = int((myRank - kk*G_division(1)*G_division(3))/G_division(1))
!ii = mod(myRank, G_division(1))
jj = 1
ii = 1

k = kp + kk * int(G_size(3)/G_division(3))
j = jp + jj * int(G_size(2)/G_division(2))
i = ip + ii * int(G_size(1)/G_division(1))

write(10, *) G_origin(1) + (i - 0.5d0)*dh(1),  &
G_origin(2) + (j - 0.5d0)*dh(2),  &
G_origin(3) + (k - 0.5d0)*dh(3)
end do
end do
end do


write(10, "('POINT_DATA ', i9)") (ix+2*g)*(jx+2*g)*(kx+2*g)

! velocity
write(10, "('VECTORS velocity float')")
do k = 1-g, kx+g
do j = 1-g, jx+g
do i = 1-g, ix+g
write(10, *) v(i, j, k, 1), v(i, j, k, 2), v(i, j, k, 3)
end do
end do
end do

!       ! umean
!       write(10, "('SCALARS umean float')")
!       write(10, "('LOOKUP_TABLE default')")
!       do k = 0, kx
!          do j = 0, jx
!             do i = 0, ix
!                write(10, *) vmean(i, j, k, 1)
!             end do
!          end do
!       end do

! pressure
write(10, "('SCALARS pressure float')")
write(10, "('LOOKUP_TABLE default')")
do k = 1-g, kx+g
do j = 1-g, jx+g
do i = 1-g, ix+g
write(10, *) p(i, j, k)
end do
end do
end do

close(10)

return
end subroutine output_vtk




!********************************************************************
subroutine output_mean(step, G_origin, G_region, G_division, G_size, myRank, sz, dh, g, vmean, rms, rmsmean)
implicit none
integer                                                   :: step
real, dimension(3)                                        :: G_origin, G_region, dh
integer, dimension(3)                                     :: G_division, G_size, sz
integer                                                   :: g
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) :: vmean, rms, rmsmean
real, dimension(1:sz(2), 3)                               :: global_vmean, global_rmsmean
integer                                                   :: ix, jx, kx, i, j, k, ip, jp, kp, ii, jj, kk, myRank
real                                                      :: xc, yc, zc, ymean

character fname1*128, fname2*128, fname3*128, fname4*128, fname5*128, fname6*128, fname7*128, fname8*128

ix = sz(1)
jx = sz(2)
kx = sz(3)

if ( step < 10 ) then
write(fname1, "('umean_step',    i1.0, '.dat')"), step
write(fname2, "('urmsmean_step', i1.0, '.dat')"), step
write(fname3, "('vrmsmean_step', i1.0, '.dat')"), step
write(fname4, "('wrmsmean_step', i1.0, '.dat')"), step
write(fname5, "('global_umean_step', i1.0, '.dat')"), step
write(fname6, "('global_urmsmean_step', i1.0, '.dat')"), step
write(fname7, "('global_vrmsmean_step', i1.0, '.dat')"), step
write(fname8, "('global_wrmsmean_step', i1.0, '.dat')"), step
else if ( step < 100 ) then
write(fname1, "('umean_step',    i2.0, '.dat')"), step
write(fname2, "('urmsmean_step', i2.0, '.dat')"), step
write(fname3, "('vrmsmean_step', i2.0, '.dat')"), step
write(fname4, "('wrmsmean_step', i2.0, '.dat')"), step
write(fname5, "('global_umean_step', i2.0, '.dat')"), step
write(fname6, "('global_urmsmean_step', i2.0, '.dat')"), step
write(fname7, "('global_vrmsmean_step', i2.0, '.dat')"), step
write(fname8, "('global_wrmsmean_step', i2.0, '.dat')"), step
else if ( step < 1000 ) then
write(fname1, "('umean_step',    i3.0, '.dat')"), step
write(fname2, "('urmsmean_step', i3.0, '.dat')"), step
write(fname3, "('vrmsmean_step', i3.0, '.dat')"), step
write(fname4, "('wrmsmean_step', i3.0, '.dat')"), step
write(fname5, "('global_umean_step', i3.0, '.dat')"), step
write(fname6, "('global_urmsmean_step', i3.0, '.dat')"), step
write(fname7, "('global_vrmsmean_step', i3.0, '.dat')"), step
write(fname8, "('global_wrmsmean_step', i3.0, '.dat')"), step
else if ( step < 10000 ) then
write(fname1, "('umean_step',    i4.0, '.dat')"), step
write(fname2, "('urmsmean_step', i4.0, '.dat')"), step
write(fname3, "('vrmsmean_step', i4.0, '.dat')"), step
write(fname4, "('wrmsmean_step', i4.0, '.dat')"), step
write(fname5, "('global_umean_step', i4.0, '.dat')"), step
write(fname6, "('global_urmsmean_step', i4.0, '.dat')"), step
write(fname7, "('global_vrmsmean_step', i4.0, '.dat')"), step
write(fname8, "('global_wrmsmean_step', i4.0, '.dat')"), step
else if ( step < 100000 ) then
write(fname1, "('umean_step',    i5.0, '.dat')"), step
write(fname2, "('urmsmean_step', i5.0, '.dat')"), step
write(fname3, "('vrmsmean_step', i5.0, '.dat')"), step
write(fname4, "('wrmsmean_step', i5.0, '.dat')"), step
write(fname5, "('global_umean_step', i5.0, '.dat')"), step
write(fname6, "('global_urmsmean_step', i5.0, '.dat')"), step
write(fname7, "('global_vrmsmean_step', i5.0, '.dat')"), step
write(fname8, "('global_wrmsmean_step', i5.0, '.dat')"), step
else if ( step < 1000000 ) then
write(fname1, "('umean_step',    i6.0, '.dat')"), step
write(fname2, "('urmsmean_step', i6.0, '.dat')"), step
write(fname3, "('vrmsmean_step', i6.0, '.dat')"), step
write(fname4, "('wrmsmean_step', i6.0, '.dat')"), step
write(fname5, "('global_umean_step', i6.0, '.dat')"), step
write(fname6, "('global_urmsmean_step', i6.0, '.dat')"), step
write(fname7, "('global_vrmsmean_step', i6.0, '.dat')"), step
write(fname8, "('global_wrmsmean_step', i6.0, '.dat')"), step
end if

do j = 1, jx
global_vmean(j, 1)   = 0.0d0
global_vmean(j, 2)   = 0.0d0
global_vmean(j, 3)   = 0.0d0
global_rmsmean(j, 1) = 0.0d0
global_rmsmean(j, 2) = 0.0d0
global_rmsmean(j, 3) = 0.0d0
end do

do j = 1, jx
do k = 1, kx
do i = 1, ix
global_vmean(j, 1)   = global_vmean(j, 1) + vmean(i, j, k, 1)
global_vmean(j, 2)   = global_vmean(j, 2) + vmean(i, j, k, 2)
global_vmean(j, 3)   = global_vmean(j, 3) + vmean(i, j, k, 3)
global_rmsmean(j, 1) = global_rmsmean(j, 1) + rmsmean(i, j, k, 1)
global_rmsmean(j, 2) = global_rmsmean(j, 2) + rmsmean(i, j, k, 2)
global_rmsmean(j, 3) = global_rmsmean(j, 3) + rmsmean(i, j, k, 3)
end do
end do
global_vmean(j, 1)   = global_vmean(j, 1)/(ix*kx)
global_vmean(j, 2)   = global_vmean(j, 2)/(ix*kx)
global_vmean(j, 3)   = global_vmean(j, 3)/(ix*kx)
global_rmsmean(j, 1) = global_rmsmean(j, 1)/(ix*kx)
global_rmsmean(j, 2) = global_rmsmean(j, 2)/(ix*kx)
global_rmsmean(j, 3) = global_rmsmean(j, 3)/(ix*kx)
end do


do kp = 1, kx
do jp = 1, jx
do ip = 1, ix
kk = int(myRank/(G_division(1)*G_division(2)))
jj = int((myRank - kk*G_division(1)*G_division(3))/G_division(1))
ii = mod(myRank, G_division(1))
k = kp + kk * int(G_size(3)/G_division(3))
j = jp + jj * int(G_size(2)/G_division(2))
i = ip + ii * int(G_size(1)/G_division(1))
xc = G_origin(1) + (i - 0.5d0)*dh(1)
yc = G_origin(2) + (j - 0.5d0)*dh(2)
zc = G_origin(3) + (k - 0.5d0)*dh(3)

if ( (abs(xc - G_region(1)/2.0d0) < dh(1)).and.(abs(yc - G_region(2)/2.0d0) < dh(2)) &
.and.(abs(zc - G_region(3)/2.0d0) < dh(3)).and.(xc < G_region(1)/2.0d0) &
.and.(yc < G_region(2)/2.0d0).and.(zc < G_region(3)/2.0d0) ) then
open(10, file = fname1)
open(11, file = fname2)
open(12, file = fname3)
open(13, file = fname4)
open(20, file = "umean_latest.dat")
open(21, file = "urmsmean_latest.dat")
open(22, file = "vrmsmean_latest.dat")
open(23, file = "wrmsmean_latest.dat")
do jj = 1, jx
ymean = G_origin(2) + (jj - 0.5d0)*dh(2)
write(10, *) ymean, vmean(ip, jj, kp, 1)
write(11, *) ymean, rmsmean(ip, jj, kp, 1)
write(12, *) ymean, rmsmean(ip, jj, kp, 2)
write(13, *) ymean, rmsmean(ip, jj, kp, 3)
write(20, *) ymean, vmean(ip, jj, kp, 1)
write(21, *) ymean, rmsmean(ip, jj, kp, 1)
write(22, *) ymean, rmsmean(ip, jj, kp, 2)
write(23, *) ymean, rmsmean(ip, jj, kp, 3)
end do
close(10)
close(11)
close(12)
close(13)
close(20)
close(21)
close(22)
close(23)
end if

if ( (i == 1).and.(k == 1) ) then
open(14, file = fname5)
open(15, file = fname6)
open(16, file = fname7)
open(17, file = fname8)
open(24, file = "global_umean_latest.dat")
open(25, file = "global_urmsmean_latest.dat")
open(26, file = "global_vrmsmean_latest.dat")
open(27, file = "global_wrmsmean_latest.dat")
do jj = 1, jx
ymean = G_origin(2) + (jj - 0.5d0)*dh(2)
write(14, *) ymean, global_vmean(jj, 1)
write(15, *) ymean, global_rmsmean(jj, 1)
write(16, *) ymean, global_rmsmean(jj, 2)
write(17, *) ymean, global_rmsmean(jj, 3)
write(24, *) ymean, global_vmean(jj, 1)
write(25, *) ymean, global_rmsmean(jj, 1)
write(26, *) ymean, global_rmsmean(jj, 2)
write(27, *) ymean, global_rmsmean(jj, 3)
end do
close(14)
close(15)
close(16)
close(17)
close(24)
close(25)
close(26)
close(27)
end if

end do
end do
end do


return
end subroutine output_mean
