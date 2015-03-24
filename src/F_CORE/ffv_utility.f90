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
subroutine output_vtk(step, G_origin, G_division, G_size, myRank, sz, dx, dy, dz, g, v, p)
implicit none
integer :: step
real, dimension(3) ::  G_origin
integer, dimension(3) :: G_division, G_size, sz
integer :: g
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  p
integer :: ix, jx, kx, i, j, k, myRank, ip, jp, kp, ii, jj, kk
real :: dx, dy, dz
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
write(10, "('DIMENSIONS ', 3(1x, i8))") ix, jx, kx
write(10, "('POINTS ', i20, ' float')") ix*jx*kx

do kp = 1, kx
do jp = 1, jx
do ip = 1, ix

ii = 0
jj = 0

k = kp + myRank * int(G_size(3)/G_division(3))
j = jp + jj * int(G_size(2)/G_division(2))
i = ip + ii * int(G_size(1)/G_division(1))

write(10, *) G_origin(1) + (i - 0.5d0)*dx,  &
G_origin(2) + (j - 0.5d0)*dy,  &
G_origin(3) + (k - 0.5d0)*dz
end do
end do
end do


write(10, "('POINT_DATA ', i9)") ix*jx*kx

! velocity
write(10, "('VECTORS velocity float')")
do k = 1, kx
do j = 1, jx
do i = 1, ix
write(10, *) v(i, j, k, 1), v(i, j, k, 2), v(i, j, k, 3)
end do
end do
end do

! pressure
write(10, "('SCALARS pressure float')")
write(10, "('LOOKUP_TABLE default')")
do k = 1, kx
do j = 1, jx
do i = 1, ix
write(10, *) p(i, j, k)
end do
end do
end do

close(10)

return
end subroutine output_vtk





!=====uzawa
!> ********************************************************************
!! @param [out]      R       レイノルズ応力テンソル
!! @param [in, out]  R_ave   レイノルズ応力テンソル (時間平均値)
!! @param [in]       sz      配列長
!! @param [in]       g       ガイドセル長
!! @param [in]       v       セルセンター速度ベクトル
!! @param [in]       v_ave   セルセンター時間平均速度ベクトル
!! @param [in]       nadd    加算回数
!! @param [out]      flop    flop count
!<
subroutine calc_reynolds_stress(R, R_ave, sz, g, v, v_ave, nadd, flop)
implicit none
include 'ffv_f_params.h'
integer                                                   :: ix, jx, kx, i, j, k, g
integer, dimension(3)                                     :: sz
real, dimension(6, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) :: R, R_ave
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) :: v, v_ave
real                                                      :: v1p, v2p, v3p
real                                                      :: nadd, val1, val2
double precision                                          :: flop

ix = sz(1)
jx = sz(2)
kx = sz(3)

! 9 flop
val2 = 1.0/nadd
val1 = 1.0 - val2

flop = flop + dble(ix+2*g) * dble(jx+2*g) * dble(kx+2*g) * 27.0d0 + 9.0d0

!$OMP PARALLEL &
!$OMP PRIVATE(v1p, v2p, v3p) &
!$OMP FIRSTPRIVATE(ix, jx, kx, g, val1, val2)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k = 1-g, kx+g
do j = 1-g, jx+g
do i = 1-g, ix+g

      ! 変動速度ベクトルのロード
      ! 3 flop
      v1p = v(i, j, k, 1) - v_ave(i, j, k, 1)
      v2p = v(i, j, k, 2) - v_ave(i, j, k, 2)
      v3p = v(i, j, k, 3) - v_ave(i, j, k, 3)

      ! レイノルズ応力テンソル
      ! 6 flop
      R(1, i, j, k) = v1p * v1p
      R(2, i, j, k) = v1p * v2p
      R(3, i, j, k) = v1p * v3p
      R(4, i, j, k) = v2p * v2p
      R(5, i, j, k) = v2p * v3p
      R(6, i, j, k) = v3p * v3p

      ! レイノルズ応力テンソル (時間平均値)
      ! 18 flop
      R_ave(1, i, j, k) = val1 * R_ave(1, i, j, k) + val2 * R(1, i, j, k)
      R_ave(2, i, j, k) = val1 * R_ave(2, i, j, k) + val2 * R(2, i, j, k)
      R_ave(3, i, j, k) = val1 * R_ave(3, i, j, k) + val2 * R(3, i, j, k)
      R_ave(4, i, j, k) = val1 * R_ave(4, i, j, k) + val2 * R(4, i, j, k)
      R_ave(5, i, j, k) = val1 * R_ave(5, i, j, k) + val2 * R(5, i, j, k)
      R_ave(6, i, j, k) = val1 * R_ave(6, i, j, k) + val2 * R(6, i, j, k)

end do
end do
end do
!$OMP END PARALLEL

return
end subroutine calc_reynolds_stress
!> ********************************************************************





!> ********************************************************************
!! @brief 生成項の計算
!! @param [in, out]  P_ave  生成項 (時間平均値)
!! @param [in]       sz     配列長
!! @param [in]       dh     格子幅
!! @param [in]       g      ガイドセル長
!! @param [in]       v_ave  セルセンター時間平均速度ベクトル
!! @param [in]       R      レイノルズ応力テンソル
!! @param [in]       bv     BCindex C
!! @param [in]       nadd   加算回数
!! @param [out]      flop   flop count
!<
subroutine calc_production_rate (P_ave, sz, dh, g, v_ave, R, bv, nadd, flop)
implicit none
include 'ffv_f_params.h'
integer                                                   ::  i, j, k, ix, jx, kx, g, idx
integer, dimension(3)                                     ::  sz
real, dimension(3)                                        ::  dh
real                                                      ::  R11, R12, R13, R21, R22, R23, R31, R32, R33
integer                                                   ::  b_e1, b_w1, b_n1, b_s1, b_t1, b_b1
real                                                      ::  actv, rx, ry, rz
real                                                      ::  Up0, Ue1, Uw1, Us1, Un1, Ub1, Ut1
real                                                      ::  Vp0, Ve1, Vw1, Vs1, Vn1, Vb1, Vt1
real                                                      ::  Wp0, We1, Ww1, Ws1, Wn1, Wb1, Wt1
real                                                      ::  w_e, w_w, w_n, w_s, w_t, w_b, uq, vq, wq
real                                                      ::  gradU11, gradU12, gradU13
real                                                      ::  gradU21, gradU22, gradU23
real                                                      ::  gradU31, gradU32, gradU33
real                                                      ::  P11_1, P12_1, P13_1, P22_1, P23_1, P33_1
real                                                      ::  P21_1, P31_1, P32_1
real                                                      ::  P11, P12, P13, P22, P23, P33
real, dimension(6, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  P_ave, R
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v_ave
real                                                      ::  nadd, val1, val2
integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bv
double precision                                          ::  flop

ix = sz(1)
jx = sz(2)
kx = sz(3)

! 27 flop
rx = 1.0 / (2.0*dh(1))
ry = 1.0 / (2.0*dh(2))
rz = 1.0 / (2.0*dh(3))

! 9 flop
val2 = 1.0/nadd
val1 = 1.0 - val2

flop = flop + dble(ix)*dble(jx)*dble(kx)*133.0d0 + 36.0d0

!$OMP PARALLEL &
!$OMP PRIVATE(idx, actv, uq, vq, wq) &
!$OMP PRIVATE(Up0, Ue1, Uw1, Us1, Un1, Ub1, Ut1) &
!$OMP PRIVATE(Vp0, Ve1, Vw1, Vs1, Vn1, Vb1, Vt1) &
!$OMP PRIVATE(Wp0, We1, Ww1, Ws1, Wn1, Wb1, Wt1) &
!$OMP PRIVATE(b_e1, b_w1, b_n1, b_s1, b_t1, b_b1) &
!$OMP PRIVATE(w_e, w_w, w_n, w_s, w_t, w_b) &
!$OMP PRIVATE(R11, R12, R13, R22, R23, R33) &
!$OMP PRIVATE(gradU11, gradU12, gradU13) &
!$OMP PRIVATE(gradU21, gradU22, gradU23) &
!$OMP PRIVATE(gradU31, gradU32, gradU33) &
!$OMP PRIVATE(P11_1, P12_1, P13_1, P22_1, P23_1, P33_1) &
!$OMP PRIVATE(P21_1, P31_1, P32_1) &
!$OMP PRIVATE(P11, P12, P13, P22, P23, P33) &
!$OMP FIRSTPRIVATE(ix, jx, kx, rx, ry, rz)
!$OMP DO SCHEDULE(static)

do k = 1, kx
do j = 1, jx
do i = 1, ix

      ! レイノルズ応力テンソル
      R11 = R(1, i, j, k)
      R12 = R(2, i, j, k)
      R13 = R(3, i, j, k)
      R21 = R12
      R22 = R(4, i, j, k)
      R23 = R(5, i, j, k)
      R31 = R13
      R32 = R23
      R33 = R(6, i, j, k)

      idx = bv(i,j,k)

      Ub1 = v_ave(i  ,j  ,k-1, 1)
      Us1 = v_ave(i  ,j-1,k  , 1)
      Uw1 = v_ave(i-1,j  ,k  , 1)
      Up0 = v_ave(i  ,j  ,k  , 1)
      Ue1 = v_ave(i+1,j  ,k  , 1)
      Un1 = v_ave(i  ,j+1,k  , 1)
      Ut1 = v_ave(i  ,j  ,k+1, 1)

      Vb1 = v_ave(i  ,j  ,k-1, 2)
      Vs1 = v_ave(i  ,j-1,k  , 2)
      Vw1 = v_ave(i-1,j  ,k  , 2)
      Vp0 = v_ave(i  ,j  ,k  , 2)
      Ve1 = v_ave(i+1,j  ,k  , 2)
      Vn1 = v_ave(i  ,j+1,k  , 2)
      Vt1 = v_ave(i  ,j  ,k+1, 2)

      Wb1 = v_ave(i  ,j  ,k-1, 3)
      Ws1 = v_ave(i  ,j-1,k  , 3)
      Ww1 = v_ave(i-1,j  ,k  , 3)
      Wp0 = v_ave(i  ,j  ,k  , 3)
      We1 = v_ave(i+1,j  ,k  , 3)
      Wn1 = v_ave(i  ,j+1,k  , 3)
      Wt1 = v_ave(i  ,j  ,k+1, 3)

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
      uq = -Up0
      vq = -Vp0
      wq = -Wp0

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

      ! 平均速度勾配テンソル 
      ! 27 flop
      gradU11 = rx * ( Ue1 - Uw1 ) * actv
      gradU12 = rx * ( Ve1 - Vw1 ) * actv
      gradU13 = rx * ( We1 - Ww1 ) * actv
      gradU21 = ry * ( Un1 - Us1 ) * actv
      gradU22 = ry * ( Vn1 - Vs1 ) * actv
      gradU23 = ry * ( Wn1 - Ws1 ) * actv
      gradU31 = rz * ( Ut1 - Ub1 ) * actv
      gradU32 = rz * ( Vt1 - Vb1 ) * actv
      gradU33 = rz * ( Wt1 - Wb1 ) * actv

      ! 生成項第一項
      ! 45 flop
      P11_1 = R11*gradU11 + R12*gradU21 + R13*gradU31
      P12_1 = R11*gradU12 + R12*gradU22 + R13*gradU32
      P13_1 = R11*gradU13 + R12*gradU23 + R13*gradU33
      P21_1 = R21*gradU11 + R22*gradU21 + R23*gradU31
      P22_1 = R21*gradU12 + R22*gradU22 + R23*gradU32
      P23_1 = R21*gradU13 + R22*gradU23 + R23*gradU33
      P31_1 = R31*gradU11 + R32*gradU21 + R33*gradU31
      P32_1 = R31*gradU12 + R32*gradU22 + R33*gradU32
      P33_1 = R31*gradU13 + R32*gradU23 + R33*gradU33

      ! 生成項 (生成項第一項 + 第一項の転置項)
      ! 24 flop
      P11 = -( P11_1 + P11_1 ) * actv
      P12 = -( P12_1 + P21_1 ) * actv
      P13 = -( P13_1 + P31_1 ) * actv
      P22 = -( P22_1 + P22_1 ) * actv
      P23 = -( P23_1 + P32_1 ) * actv
      P33 = -( P33_1 + P33_1 ) * actv

      ! 生成項 (時間平均値)
      ! 18 flop
      P_ave(1, i, j, k) = val1 * P_ave(1, i, j, k) + val2 * P11
      P_ave(2, i, j, k) = val1 * P_ave(2, i, j, k) + val2 * P12
      P_ave(3, i, j, k) = val1 * P_ave(3, i, j, k) + val2 * P13
      P_ave(4, i, j, k) = val1 * P_ave(4, i, j, k) + val2 * P22
      P_ave(5, i, j, k) = val1 * P_ave(5, i, j, k) + val2 * P23
      P_ave(6, i, j, k) = val1 * P_ave(6, i, j, k) + val2 * P33

end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL

return
end subroutine calc_production_rate
!> ********************************************************************





!> ********************************************************************
!! @brief 散逸項の計算
!! @param [out]  E_ave  散逸項 (時間平均値)
!! @param [in]   sz     配列長
!! @param [in]   dh     格子幅
!! @param [in]   g      ガイドセル長
!! @param [in]   nu     動粘性係数
!! @param [in]   v      セルセンター速度ベクトル
!! @param [in]   v_ave  セルセンター時間平均速度ベクトル
!! @param [in]   bv     BCindex C
!! @param [in]   nadd   加算回数
!! @param [out]  flop   flop count
!<
subroutine calc_dissipation_rate (E_ave, sz, dh, g, nu, v, v_ave, bv, nadd, flop)
implicit none
include 'ffv_f_params.h'
integer                                                   ::  i, j, k, ix, jx, kx, g, idx
integer, dimension(3)                                     ::  sz
real, dimension(3)                                        ::  dh
integer                                                   ::  b_e1, b_w1, b_n1, b_s1, b_t1, b_b1
real                                                      ::  actv, rx, ry, rz
real                                                      ::  Up0, Ue1, Uw1, Us1, Un1, Ub1, Ut1
real                                                      ::  Vp0, Ve1, Vw1, Vs1, Vn1, Vb1, Vt1
real                                                      ::  Wp0, We1, Ww1, Ws1, Wn1, Wb1, Wt1
real                                                      ::  w_e, w_w, w_n, w_s, w_t, w_b, uq, vq, wq
real                                                      ::  gradU11, gradU12, gradU13
real                                                      ::  gradU21, gradU22, gradU23
real                                                      ::  gradU31, gradU32, gradU33
real                                                      ::  gradUt11, gradUt12, gradUt13
real                                                      ::  gradUt21, gradUt22, gradUt23
real                                                      ::  gradUt31, gradUt32, gradUt33
real                                                      ::  E11, E12, E13, E22, E23, E33
real, dimension(6, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  E_ave
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v, v_ave
real                                                      ::  nadd, val1, val2, nu
integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bv
double precision                                          ::  flop

ix = sz(1)
jx = sz(2)
kx = sz(3)

! 27 flop
rx = 1.0 / (2.0*dh(1))
ry = 1.0 / (2.0*dh(2))
rz = 1.0 / (2.0*dh(3))

! 9 flop
val2 = 1.0/nadd
val1 = 1.0 - val2

flop = flop + dble(ix)*dble(jx)*dble(kx)*131.0d0 + 36.0d0


!$OMP PARALLEL &
!$OMP PRIVATE(idx, actv, uq, vq, wq) &
!$OMP PRIVATE(Up0, Ue1, Uw1, Us1, Un1, Ub1, Ut1) &
!$OMP PRIVATE(Vp0, Ve1, Vw1, Vs1, Vn1, Vb1, Vt1) &
!$OMP PRIVATE(Wp0, We1, Ww1, Ws1, Wn1, Wb1, Wt1) &
!$OMP PRIVATE(b_e1, b_w1, b_n1, b_s1, b_t1, b_b1) &
!$OMP PRIVATE(w_e, w_w, w_n, w_s, w_t, w_b) &
!$OMP PRIVATE(gradU11, gradU12, gradU13) &
!$OMP PRIVATE(gradU21, gradU22, gradU23) &
!$OMP PRIVATE(gradU31, gradU32, gradU33) &
!$OMP PRIVATE(gradUt11, gradUt12, gradUt13) &
!$OMP PRIVATE(gradUt21, gradUt22, gradUt23) &
!$OMP PRIVATE(gradUt31, gradUt32, gradUt33) &
!$OMP PRIVATE(E11, E12, E13, E22, E23, E33) &
!$OMP FIRSTPRIVATE(ix, jx, kx, rx, ry, rz, nu)
!$OMP DO SCHEDULE(static)

do k = 1, kx
do j = 1, jx
do i = 1, ix
      idx = bv(i,j,k)

      ! 21 flop
      Ub1 = v(i  ,j  ,k-1, 1) - v_ave(i  ,j  ,k-1, 1)
      Us1 = v(i  ,j-1,k  , 1) - v_ave(i  ,j-1,k  , 1)
      Uw1 = v(i-1,j  ,k  , 1) - v_ave(i-1,j  ,k  , 1)
      Up0 = v(i  ,j  ,k  , 1) - v_ave(i  ,j  ,k  , 1)
      Ue1 = v(i+1,j  ,k  , 1) - v_ave(i+1,j  ,k  , 1)
      Un1 = v(i  ,j+1,k  , 1) - v_ave(i  ,j+1,k  , 1)
      Ut1 = v(i  ,j  ,k+1, 1) - v_ave(i  ,j  ,k+1, 1)

      Vb1 = v(i  ,j  ,k-1, 2) - v_ave(i  ,j  ,k-1, 2)
      Vs1 = v(i  ,j-1,k  , 2) - v_ave(i  ,j-1,k  , 2)
      Vw1 = v(i-1,j  ,k  , 2) - v_ave(i-1,j  ,k  , 2)
      Vp0 = v(i  ,j  ,k  , 2) - v_ave(i  ,j  ,k  , 2)
      Ve1 = v(i+1,j  ,k  , 2) - v_ave(i+1,j  ,k  , 2)
      Vn1 = v(i  ,j+1,k  , 2) - v_ave(i  ,j+1,k  , 2)
      Vt1 = v(i  ,j  ,k+1, 2) - v_ave(i  ,j  ,k+1, 2)

      Wb1 = v(i  ,j  ,k-1, 3) - v_ave(i  ,j  ,k-1, 3)
      Ws1 = v(i  ,j-1,k  , 3) - v_ave(i  ,j-1,k  , 3)
      Ww1 = v(i-1,j  ,k  , 3) - v_ave(i-1,j  ,k  , 3)
      Wp0 = v(i  ,j  ,k  , 3) - v_ave(i  ,j  ,k  , 3)
      We1 = v(i+1,j  ,k  , 3) - v_ave(i+1,j  ,k  , 3)
      Wn1 = v(i  ,j+1,k  , 3) - v_ave(i  ,j+1,k  , 3)
      Wt1 = v(i  ,j  ,k+1, 3) - v_ave(i  ,j  ,k+1, 3)

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
      uq = -Up0
      vq = -Vp0
      wq = -Wp0

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

      ! 速度変動勾配テンソル 
      ! 27 flop
      gradU11 = rx * ( Ue1 - Uw1 ) * actv
      gradU12 = rx * ( Ve1 - Vw1 ) * actv
      gradU13 = rx * ( We1 - Ww1 ) * actv
      gradU21 = ry * ( Un1 - Us1 ) * actv
      gradU22 = ry * ( Vn1 - Vs1 ) * actv
      gradU23 = ry * ( Wn1 - Ws1 ) * actv
      gradU31 = rz * ( Ut1 - Ub1 ) * actv
      gradU32 = rz * ( Vt1 - Vb1 ) * actv
      gradU33 = rz * ( Wt1 - Wb1 ) * actv

      ! 速度変動勾配テンソルの転置
      gradUt11 = gradU11
      gradUt12 = gradU21
      gradUt13 = gradU31
      gradUt21 = gradU12
      gradUt22 = gradU22
      gradUt23 = gradU32
      gradUt31 = gradU13
      gradUt32 = gradU23
      gradUt33 = gradU33

      ! 散逸項
      ! 48 flop
      E11 = 2 * nu * ( gradUt11*gradU11 + gradUt12*gradU21 + gradUt13*gradU31 ) * actv
      E12 = 2 * nu * ( gradUt11*gradU12 + gradUt12*gradU22 + gradUt13*gradU32 ) * actv
      E13 = 2 * nu * ( gradUt11*gradU13 + gradUt12*gradU23 + gradUt13*gradU33 ) * actv
      E22 = 2 * nu * ( gradUt21*gradU12 + gradUt22*gradU22 + gradUt23*gradU32 ) * actv
      E23 = 2 * nu * ( gradUt21*gradU13 + gradUt22*gradU23 + gradUt23*gradU33 ) * actv
      E33 = 2 * nu * ( gradUt31*gradU13 + gradUt32*gradU23 + gradUt33*gradU33 ) * actv

      ! 散逸項 (時間平均値)
      ! 18 flop
      E_ave(1, i, j, k) = val1 * E_ave(1, i, j, k) + val2 * E11
      E_ave(2, i, j, k) = val1 * E_ave(2, i, j, k) + val2 * E12
      E_ave(3, i, j, k) = val1 * E_ave(3, i, j, k) + val2 * E13
      E_ave(4, i, j, k) = val1 * E_ave(4, i, j, k) + val2 * E22
      E_ave(5, i, j, k) = val1 * E_ave(5, i, j, k) + val2 * E23
      E_ave(6, i, j, k) = val1 * E_ave(6, i, j, k) + val2 * E33

end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL


return
end subroutine calc_dissipation_rate
!> ********************************************************************





!> ********************************************************************
!! @brief 乱流拡散項の計算
!! @param [in, out]  T_ave  乱流拡散項 (時間平均値)
!! @param [in]       sz     配列長
!! @param [in]       dh     格子幅
!! @param [in]       g      ガイドセル長
!! @param [in]       v      セルセンター速度ベクトル
!! @param [in]       v_ave  セルセンター時間平均速度ベクトル
!! @param [in]       R      レイノルズ応力テンソル
!! @param [in]       bv     BCindex C
!! @param [in]       nadd   加算回数
!! @param [out]      flop   flop count
!<
subroutine calc_turb_transport_rate (T_ave, sz, dh, g, v, v_ave, R, bv, nadd, flop)
implicit none
include 'ffv_f_params.h'
integer                                                   ::  i, j, k, ix, jx, kx, g, idx
real, dimension(3)                                        ::  dh
double precision                                          ::  flop
integer, dimension(3)                                     ::  sz
integer                                                   ::  b_e1, b_w1, b_n1, b_s1, b_t1, b_b1
real                                                      ::  actv, rx, ry, rz
real                                                      ::  v1p, v2p, v3p
real                                                      ::  R11_p0, R12_p0, R13_p0, R22_p0, R23_p0, R33_p0
real                                                      ::  R11_w1, R11_e1, R11_s1, R11_n1, R11_b1, R11_t1
real                                                      ::  R12_w1, R12_e1, R12_s1, R12_n1, R12_b1, R12_t1
real                                                      ::  R13_w1, R13_e1, R13_s1, R13_n1, R13_b1, R13_t1
real                                                      ::  R22_w1, R22_e1, R22_s1, R22_n1, R22_b1, R22_t1
real                                                      ::  R23_w1, R23_e1, R23_s1, R23_n1, R23_b1, R23_t1
real                                                      ::  R33_w1, R33_e1, R33_s1, R33_n1, R33_b1, R33_t1
real                                                      ::  T11, T12, T13, T22, T23, T33
real                                                      ::  w_e, w_w, w_n, w_s, w_t, w_b, uq, vq, wq
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v, v_ave
real, dimension(6, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  R, T_ave
integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bv
real                                                      ::  nadd, val1, val2

ix = sz(1)
jx = sz(2)
kx = sz(3)

! 27 flop
rx = 1.0 / (2.0*dh(1))
ry = 1.0 / (2.0*dh(2))
rz = 1.0 / (2.0*dh(3))

! 9 flop
val2 = 1.0/nadd
val1 = 1.0 - val2

flop = flop + dble(ix)*dble(jx)*dble(kx)*118.0d0 + 36.0d0

!$OMP PARALLEL &
!$OMP PRIVATE(idx, actv, uq, vq, wq) &
!$OMP PRIVATE(v1p, v2p, v3p) &
!$OMP PRIVATE(R11_w1, R11_e1, R11_s1, R11_n1, R11_b1, R11_t1) &
!$OMP PRIVATE(R12_w1, R12_e1, R12_s1, R12_n1, R12_b1, R12_t1) &
!$OMP PRIVATE(R13_w1, R13_e1, R13_s1, R13_n1, R13_b1, R13_t1) &
!$OMP PRIVATE(R22_w1, R22_e1, R22_s1, R22_n1, R22_b1, R22_t1) &
!$OMP PRIVATE(R23_w1, R23_e1, R23_s1, R23_n1, R23_b1, R23_t1) &
!$OMP PRIVATE(R33_w1, R33_e1, R33_s1, R33_n1, R33_b1, R33_t1) &
!$OMP PRIVATE(b_e1, b_w1, b_n1, b_s1, b_t1, b_b1) &
!$OMP PRIVATE(w_e, w_w, w_n, w_s, w_t, w_b) &
!$OMP FIRSTPRIVATE(ix, jx, kx, rx, ry, rz, val1, val2)

!$OMP DO SCHEDULE(static) COLLAPSE(2)

do k=1,kx
do j=1,jx
do i=1,ix

      idx = bv(i,j,k)

      ! 変動速度ベクトル
      ! 3 flop
      v1p = v(i, j, k, 1) - v_ave(i, j, k, 1)
      v2p = v(i, j, k, 2) - v_ave(i, j, k, 2)
      v3p = v(i, j, k, 3) - v_ave(i, j, k, 3)

      ! レイノルズ応力テンソル
      R11_b1 = R(1, i,   j,   k-1)
      R11_s1 = R(1, i,   j-1, k  )
      R11_w1 = R(1, i-1, j,   k  )
      R11_p0 = R(1, i,   j,   k  )
      R11_e1 = R(1, i+1, j,   k  )
      R11_n1 = R(1, i,   j+1, k  )
      R11_t1 = R(1, i,   j,   k+1)

      R12_b1 = R(2, i,   j,   k-1)
      R12_s1 = R(2, i,   j-1, k  )
      R12_w1 = R(2, i-1, j,   k  )
      R12_p0 = R(2, i,   j,   k  )
      R12_e1 = R(2, i+1, j,   k  )
      R12_n1 = R(2, i,   j+1, k  )
      R12_t1 = R(2, i,   j,   k+1)

      R13_b1 = R(3, i,   j,   k-1)
      R13_s1 = R(3, i,   j-1, k  )
      R13_w1 = R(3, i-1, j,   k  )
      R13_p0 = R(3, i,   j,   k  )
      R13_e1 = R(3, i+1, j,   k  )
      R13_n1 = R(3, i,   j+1, k  )
      R13_t1 = R(3, i,   j,   k+1)

      R22_b1 = R(4, i,   j,   k-1)
      R22_s1 = R(4, i,   j-1, k  )
      R22_w1 = R(4, i-1, j,   k  )
      R22_p0 = R(4, i,   j,   k  )
      R22_e1 = R(4, i+1, j,   k  )
      R22_n1 = R(4, i,   j+1, k  )
      R22_t1 = R(4, i,   j,   k+1)

      R23_b1 = R(5, i,   j,   k-1)
      R23_s1 = R(5, i,   j-1, k  )
      R23_w1 = R(5, i-1, j,   k  )
      R23_p0 = R(5, i,   j,   k  )
      R23_e1 = R(5, i+1, j,   k  )
      R23_n1 = R(5, i,   j+1, k  )
      R23_t1 = R(5, i,   j,   k+1)

      R33_b1 = R(6, i,   j,   k-1)
      R33_s1 = R(6, i,   j-1, k  )
      R33_w1 = R(6, i-1, j,   k  )
      R33_p0 = R(6, i,   j,   k  )
      R33_e1 = R(6, i+1, j,   k  )
      R33_n1 = R(6, i,   j+1, k  )
      R33_t1 = R(6, i,   j,   k+1)

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
      
      ! 壁面の場合の参照速度の修正
      if ( b_e1 == 0 ) then
        R11_e1 = -R11_p0
        R12_e1 = -R12_p0
        R13_e1 = -R13_p0
        R22_e1 = -R22_p0
        R23_e1 = -R23_p0
        R33_e1 = -R33_p0
      endif
      
      if ( b_w1 == 0 ) then
        R11_w1 = -R11_p0
        R12_w1 = -R12_p0
        R13_w1 = -R13_p0
        R22_w1 = -R22_p0
        R23_w1 = -R23_p0
        R33_w1 = -R33_p0
      end if
      
      if ( b_n1 == 0 ) then
        R11_n1 = -R11_p0
        R12_n1 = -R12_p0
        R13_n1 = -R13_p0
        R22_n1 = -R22_p0
        R23_n1 = -R23_p0
        R33_n1 = -R33_p0
      end if
      
      if ( b_s1 == 0 ) then
        R11_s1 = -R11_p0
        R12_s1 = -R12_p0
        R13_s1 = -R13_p0
        R22_s1 = -R22_p0
        R23_s1 = -R23_p0
        R33_s1 = -R33_p0
      end if
      
      if ( b_t1 == 0 ) then
        R11_t1 = -R11_p0
        R12_t1 = -R12_p0
        R13_t1 = -R13_p0
        R22_t1 = -R22_p0
        R23_t1 = -R23_p0
        R33_t1 = -R33_p0
      end if
      
      if ( b_b1 == 0 ) then
        R11_b1 = -R11_p0
        R12_b1 = -R12_p0
        R13_b1 = -R13_p0
        R22_b1 = -R22_p0
        R23_b1 = -R23_p0
        R33_b1 = -R33_p0
      end if
      
      ! 乱流拡散項
      ! 84 flop
      T11 = - ( v1p * rx * (R11_e1 - R11_w1) &
              + v2p * ry * (R11_n1 - R11_s1) &
              + v3p * rz * (R11_t1 - R11_b1) &
              ) * actv

      T12 = - ( v1p * rx * (R12_e1 - R12_w1) &
              + v2p * ry * (R12_n1 - R12_s1) &
              + v3p * rz * (R12_t1 - R12_b1) &
              ) * actv

      T13 = - ( v1p * rx * (R13_e1 - R13_w1) &
              + v2p * ry * (R13_n1 - R13_s1) &
              + v3p * rz * (R13_t1 - R13_b1) &
              ) * actv

      T22 = - ( v1p * rx * (R22_e1 - R22_w1) &
              + v2p * ry * (R22_n1 - R22_s1) &
              + v3p * rz * (R22_t1 - R22_b1) &
              ) * actv

      T23 = - ( v1p * rx * (R23_e1 - R23_w1) &
              + v2p * ry * (R23_n1 - R23_s1) &
              + v3p * rz * (R23_t1 - R23_b1) &
              ) * actv

      T33 = - ( v1p * rx * (R33_e1 - R33_w1) &
              + v2p * ry * (R33_n1 - R33_s1) &
              + v3p * rz * (R33_t1 - R33_b1) &
              ) * actv

      ! 乱流拡散項 (時間平均値)
      ! 18 flop
      T_ave(1, i, j, k) = val1 * T_ave(1, i, j, k) + val2 * T11
      T_ave(2, i, j, k) = val1 * T_ave(2, i, j, k) + val2 * T12
      T_ave(3, i, j, k) = val1 * T_ave(3, i, j, k) + val2 * T13
      T_ave(4, i, j, k) = val1 * T_ave(4, i, j, k) + val2 * T22
      T_ave(5, i, j, k) = val1 * T_ave(5, i, j, k) + val2 * T23
      T_ave(6, i, j, k) = val1 * T_ave(6, i, j, k) + val2 * T33

end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL

return
end subroutine calc_turb_transport_rate
!> ********************************************************************





!> ********************************************************************
!! @brief 速度圧力勾配相関項の計算
!! @param [out]  PI_ave  散逸項 (時間平均値)
!! @param [in]   sz      配列長
!! @param [in]   dh      格子幅
!! @param [in]   g       ガイドセル長
!! @param [in]   v       セルセンター速度ベクトル
!! @param [in]   v_ave   セルセンター時間平均速度ベクトル
!! @param [in]   p       セルセンター圧力
!! @param [in]   p_ave   セルセンター時間平均圧力
!! @param [in]   bp      BCindex P
!! @param [in]   nadd    加算回数
!! @param [out]  flop    flop count
!<
subroutine calc_vel_pregrad_term (PI_ave, sz, dh, g, v, v_ave, p, p_ave, bp, nadd, flop)
implicit none
include 'ffv_f_params.h'
integer                                                   ::  i, j, k, ix, jx, kx, g, idx
integer, dimension(3)                                     ::  sz
real, dimension(3)                                        ::  dh
real                                                      ::  actv, rx, ry, rz
real                                                      ::  v1p, v2p, v3p
real                                                      ::  p0, pw1, pe1, ps1, pn1, pb1, pt1
real                                                      ::  PI11, PI12, PI13, PI22, PI23, PI33
real                                                      ::  ugradp11, ugradp12, ugradp13
real                                                      ::  ugradp21, ugradp22, ugradp23
real                                                      ::  ugradp31, ugradp32, ugradp33
real                                                      ::  ugradpt11, ugradpt12, ugradpt13
real                                                      ::  ugradpt21, ugradpt22, ugradpt23
real                                                      ::  ugradpt31, ugradpt32, ugradpt33
real, dimension(6, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  PI_ave
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v, v_ave
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  p, p_ave
real                                                      ::  nadd, val1, val2
integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp
double precision                                          ::  flop

ix = sz(1)
jx = sz(2)
kx = sz(3)

! 27 flop
! rx = 1.0 / (2.0*dh(1))
! ry = 1.0 / (2.0*dh(2))
! rz = 1.0 / (2.0*dh(3))
rx = 1.0 / (1.0*dh(1))
ry = 1.0 / (1.0*dh(2))
rz = 1.0 / (1.0*dh(3))

! 9 flop
val2 = 1.0/nadd
val1 = 1.0 - val2

flop = flop + dble(ix)*dble(jx)*dble(kx)*89.0d0 + 36.0d0

!$OMP PARALLEL &
!$OMP PRIVATE(idx, actv) &
!$OMP PRIVATE(v1p, v2p, v3p) &
!$OMP PRIVATE(p0, pw1, pe1, ps1, pn1, pb1, pt1) &
!$OMP PRIVATE(PI11, PI12, PI13, PI22, PI23, PI33) &
!$OMP PRIVATE(ugradp11, ugradp12, ugradp13) &
!$OMP PRIVATE(ugradp21, ugradp22, ugradp23) &
!$OMP PRIVATE(ugradp31, ugradp32, ugradp33) &
!$OMP PRIVATE(ugradpt11, ugradpt12, ugradpt13) &
!$OMP PRIVATE(ugradpt21, ugradpt22, ugradpt23) &
!$OMP PRIVATE(ugradpt31, ugradpt32, ugradpt33) &
!$OMP FIRSTPRIVATE(ix, jx, kx, rx, ry, rz)
!$OMP DO SCHEDULE(static)

do k = 1, kx
do j = 1, jx
do i = 1, ix

      ! セル状態 (0-solid / 1-fluid) > 1 flop
      actv= real(ibits(bp(i,j,k), State, 1))

      ! 変動速度ベクトル
      ! 3 flop
      v1p = v(i, j, k, 1) - v_ave(i, j, k, 1)
      v2p = v(i, j, k, 2) - v_ave(i, j, k, 2)
      v3p = v(i, j, k, 3) - v_ave(i, j, k, 3)

      ! 変動圧力
      ! 7 flop
      p0  = p(i,   j,   k  ) - p_ave(i,   j,   k  )
      pw1 = p(i-1, j,   k  ) - p_ave(i-1, j,   k  )
      pe1 = p(i+1, j,   k  ) - p_ave(i+1, j,   k  )
      ps1 = p(i,   j-1, k  ) - p_ave(i,   j-1, k  )
      pn1 = p(i,   j+1, k  ) - p_ave(i,   j+1, k  )
      pb1 = p(i,   j,   k-1) - p_ave(i,   j,   k-1)
      pt1 = p(i,   j,   k+1) - p_ave(i,   j,   k+1)

      ! 壁面の場合の参照速度の修正 (Neumann 条件)
      if ( ibits(bp(i, j, k), bc_n_W, 1) == 0 ) pw1 = p0
      if ( ibits(bp(i, j, k), bc_n_E, 1) == 0 ) pe1 = p0
      if ( ibits(bp(i, j, k), bc_n_S, 1) == 0 ) ps1 = p0
      if ( ibits(bp(i, j, k), bc_n_N, 1) == 0 ) pn1 = p0
      if ( ibits(bp(i, j, k), bc_n_B, 1) == 0 ) pb1 = p0
      if ( ibits(bp(i, j, k), bc_n_T, 1) == 0 ) pt1 = p0

      ! 壁面の場合の参照速度の修正 (Dirichlet 条件)
      if ( ibits(bp(i, j, k), bc_d_W, 1) == 0 ) pw1 = -p0
      if ( ibits(bp(i, j, k), bc_d_E, 1) == 0 ) pe1 = -p0
      if ( ibits(bp(i, j, k), bc_d_S, 1) == 0 ) ps1 = -p0
      if ( ibits(bp(i, j, k), bc_d_N, 1) == 0 ) pn1 = -p0
      if ( ibits(bp(i, j, k), bc_d_B, 1) == 0 ) pb1 = -p0
      if ( ibits(bp(i, j, k), bc_d_T, 1) == 0 ) pt1 = -p0

      ! 速度圧力勾配相関項 (第一項)
      ! 36 flop
!      ugradp11 = v1p * rx * ( pe1 - pw1 ) * actv
!      ugradp12 = v1p * ry * ( pn1 - ps1 ) * actv
!      ugradp13 = v1p * rz * ( pt1 - pb1 ) * actv
!      ugradp21 = v2p * rx * ( pe1 - pw1 ) * actv
!      ugradp22 = v2p * ry * ( pn1 - ps1 ) * actv
!      ugradp23 = v2p * rz * ( pt1 - pb1 ) * actv
!      ugradp31 = v3p * rx * ( pe1 - pw1 ) * actv
!      ugradp32 = v3p * ry * ( pn1 - ps1 ) * actv
!      ugradp33 = v3p * rz * ( pt1 - pb1 ) * actv
      ugradp11 = v1p * rx * ( pe1 - p0 ) * actv
      ugradp12 = v1p * ry * ( pn1 - p0 ) * actv
      ugradp13 = v1p * rz * ( pt1 - p0 ) * actv
      ugradp21 = v2p * rx * ( pe1 - p0 ) * actv
      ugradp22 = v2p * ry * ( pn1 - p0 ) * actv
      ugradp23 = v2p * rz * ( pt1 - p0 ) * actv
      ugradp31 = v3p * rx * ( pe1 - p0 ) * actv
      ugradp32 = v3p * ry * ( pn1 - p0 ) * actv
      ugradp33 = v3p * rz * ( pt1 - p0 ) * actv

      ! 第一項の転置
      ugradpt11 = ugradp11
      ugradpt12 = ugradp21
      ugradpt13 = ugradp31
      ugradpt21 = ugradp12
      ugradpt22 = ugradp22
      ugradpt23 = ugradp32
      ugradpt31 = ugradp13
      ugradpt32 = ugradp23
      ugradpt33 = ugradp33

      ! 速度圧力勾配相関項
      ! 24 flop
      PI11 = -( ugradp11 + ugradpt11 ) * actv
      PI12 = -( ugradp12 + ugradpt12 ) * actv
      PI13 = -( ugradp13 + ugradpt13 ) * actv
      PI22 = -( ugradp22 + ugradpt22 ) * actv
      PI23 = -( ugradp23 + ugradpt23 ) * actv
      PI33 = -( ugradp33 + ugradpt33 ) * actv

      ! 速度圧力勾配相関項 (時間平均値)
      ! 18 flop
      PI_ave(1, i, j, k) = val1 * PI_ave(1, i, j, k) + val2 * PI11
      PI_ave(2, i, j, k) = val1 * PI_ave(2, i, j, k) + val2 * PI12
      PI_ave(3, i, j, k) = val1 * PI_ave(3, i, j, k) + val2 * PI13
      PI_ave(4, i, j, k) = val1 * PI_ave(4, i, j, k) + val2 * PI22
      PI_ave(5, i, j, k) = val1 * PI_ave(5, i, j, k) + val2 * PI23
      PI_ave(6, i, j, k) = val1 * PI_ave(6, i, j, k) + val2 * PI33

end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL

return
end subroutine calc_vel_pregrad_term
!> ********************************************************************






