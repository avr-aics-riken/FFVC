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

!> @file   ffv_utility.f90
!! @brief  Utilitiy functions
!! @author aics
!<

!> ********************************************************************
!! @brief 有効セルに対する発散の最大値を計算，絶対値の最大値の位置を返す
!! @param [out] ds   残差の絶対値
!! @param [out] idx  ノルムの最大値の位置
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル長
!! @param [in]  div  発散値のベース
!! @param [in]  coef 係数
!! @param [in]  bp   BCindex P
!! @param [out] flop flop count
!<
    subroutine norm_v_div_dbg (ds, idx, sz, g, div, coef, bp, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, g, i0, j0, k0
    integer, dimension(3)                                     ::  sz, idx
    double precision                                          ::  flop, ds, r, d
    real                                                      ::  coef
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  div
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    ds = 0.0
    i0 = 0
    j0 = 0
    k0 = 0
    ds = 0.0

    flop = flop + dble(ix)*dble(jx)*dble(kx)*5.0d0

!$OMP PARALLEL &
!$OMP PRIVATE(r, d) &
!$OMP FIRSTPRIVATE(ix, jx, kx, coef) &
!$OMP SHARED(ds, i0, j0, k0)

!$OMP DO SCHEDULE(static)

    do k=1,kx
    do j=1,jx
    do i=1,ix
      r = dble(div(i,j,k) * coef) * dble(ibits(bp(i,j,k), vld_cnvg, 1)) ! 有効セルの場合 1.0
      d = abs(r)
      
      if ( d > ds ) then
!$OMP CRITICAL
        i0 = i
        j0 = j
        k0 = k
        ds = d
!$OMP END CRITICAL
      endif
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL
    
    idx(1) = i0
    idx(2) = j0
    idx(3) = k0

    return
    end subroutine norm_v_div_dbg

!> ********************************************************************
!! @brief 速度成分の最大値を計算する
!! @param [out] ds   最大値
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル長
!! @param [in]  div  セルに対する速度の和
!! @param [in]  coef 係数
!! @param [in]  bp   BCindex P
!! @param [out] flop flop count
!<
    subroutine norm_v_div_max (ds, sz, g, div, coef, bp, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, g
    integer, dimension(3)                                     ::  sz
    double precision                                          ::  flop, ds, r
    real                                                      ::  coef
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  div
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    ds = 0.0

    flop = flop + dble(ix)*dble(jx)*dble(kx)*6.0d0

!$OMP PARALLEL &
!$OMP REDUCTION(max:ds) &
!$OMP PRIVATE(r) &
!$OMP FIRSTPRIVATE(ix, jx, kx, coef)

!$OMP DO SCHEDULE(static)
    do k=1,kx
    do j=1,jx
    do i=1,ix
      r = dble(div(i,j,k) * coef) * dble(ibits(bp(i,j,k), vld_cnvg, 1)) ! 有効セルの場合 1.0
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

!$OMP DO SCHEDULE(static)

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
    real                                                      ::  u_ref, v_ref, w_ref, dh, h, actv
    real                                                      ::  Up0, Ue1, Uw1, Us1, Un1, Ub1, Ut1
    real                                                      ::  Vp0, Ve1, Vw1, Vs1, Vn1, Vb1, Vt1
    real                                                      ::  Wp0, We1, Ww1, Ws1, Wn1, Wb1, Wt1
    real                                                      ::  w_e, w_w, w_n, w_s, w_t, w_b, uq, vq, wq
    real                                                      ::  d11, d22, d33, d12, d13, d23
    real                                                      ::  w12, w13, w23
    real                                                      ::  q11, q22, q33, q12, q13, q23
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  q
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bv
    real, dimension(0:3)                                      ::  v00

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    h = 0.25 / dh
    
    ! 参照座標系の速度
    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)

    flop = flop + dble(ix)*dble(jx)*dble(kx)*72.0d0 + 8.0d0
    ! flop = flop + dble(ix)*dble(jx)*dble(kx)*72.0d0 + 13.0d0 ! DP

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
!$OMP FIRSTPRIVATE(ix, jx, kx, h, u_ref, v_ref, w_ref)

!$OMP DO SCHEDULE(static)

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

      d11 = h*( 2.0*(Ue1-Uw1) )
      d22 = h*( 2.0*(Vn1-Vs1) )
      d33 = h*( 2.0*(Wt1-Wb1) )
      d12 = h*( (Un1-Us1) + (Ve1-Vw1) )
      d13 = h*( (Ut1-Ub1) + (We1-Ww1) )
      d23 = h*( (Vt1-Vb1) + (Wn1-Ws1) )
      
      w12 = h*( (Un1-Us1) - (Ve1-Vw1) )
      w13 = h*( (Ut1-Ub1) - (We1-Ww1) )
      w23 = h*( (Vt1-Vb1) - (Wn1-Ws1) )
      
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
    real                                                      ::  u_ref, v_ref, w_ref, dh, h, actv
    real                                                      ::  Up0, Ue1, Uw1, Us1, Un1, Ub1, Ut1
    real                                                      ::  Vp0, Ve1, Vw1, Vs1, Vn1, Vb1, Vt1
    real                                                      ::  Wp0, We1, Ww1, Ws1, Wn1, Wb1, Wt1
    real                                                      ::  w_e, w_w, w_n, w_s, w_t, w_b, uq, vq, wq
    real                                                      ::  r1, r2, r3
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v, rot
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bv
    real, dimension(0:3)                                      ::  v00

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    h = 0.5 / dh
    
    ! 参照座標系の速度
    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)

    flop = flop + dble(ix)*dble(jx)*dble(kx)*34.0d0 + 8.0d0
    ! flop = flop + dble(ix)*dble(jx)*dble(kx)*34.0d0 + 13.0d0 ! DP

!$OMP PARALLEL &
!$OMP PRIVATE(idx, actv, uq, vq, wq) &
!$OMP PRIVATE(Up0, Ue1, Uw1, Us1, Un1, Ub1, Ut1) &
!$OMP PRIVATE(Vp0, Ve1, Vw1, Vs1, Vn1, Vb1, Vt1) &
!$OMP PRIVATE(Wp0, We1, Ww1, Ws1, Wn1, Wb1, Wt1) &
!$OMP PRIVATE(b_e1, b_w1, b_n1, b_s1, b_t1, b_b1) &
!$OMP PRIVATE(w_e, w_w, w_n, w_s, w_t, w_b) &
!$OMP PRIVATE(r1, r2, r3) &
!$OMP FIRSTPRIVATE(ix, jx, kx, h, u_ref, v_ref, w_ref)

!$OMP DO SCHEDULE(static)

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
      
      r1 = h*( (Wn1-Ws1) - (Vt1-Vb1) ) ! 12 flop
      r2 = h*( (Ut1-Ub1) - (We1-Ww1) )
      r3 = h*( (Ve1-Vw1) - (Un1-Us1) )
      
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
!! @param[out] ht ヘリシティ
!! @param sz 配列長
!! @param g ガイドセル長
!! @param dh 格子幅
!! @param v セルセンター速度ベクトル
!! @param bv BCindex C
!! @param v00 参照速度
!! @param[out] flop
!<
    subroutine helicity (ht, sz, g, dh, v, bv, v00, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, g, idx
    integer, dimension(3)                                     ::  sz
    double precision                                          ::  flop
    integer                                                   ::  b_e1, b_w1, b_n1, b_s1, b_t1, b_b1
    real                                                      ::  u_ref, v_ref, w_ref, dh, h, actv
    real                                                      ::  Up0, Ue1, Uw1, Us1, Un1, Ub1, Ut1
    real                                                      ::  Vp0, Ve1, Vw1, Vs1, Vn1, Vb1, Vt1
    real                                                      ::  Wp0, We1, Ww1, Ws1, Wn1, Wb1, Wt1
    real                                                      ::  w_e, w_w, w_n, w_s, w_t, w_b, uq, vq, wq
    real                                                      ::  r1, r2, r3, u1, u2, u3
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  ht
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bv
    real, dimension(0:3)                                      ::  v00

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    h = 0.5 / dh
    
    ! 参照座標系の速度
    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)

    flop = flop + dble(ix)*dble(jx)*dble(kx)*40.0d0 + 8.0d0
    ! flop = flop + dble(ix)*dble(jx)*dble(kx)*40.0d0 + 13.0d0 ! DP

!$OMP PARALLEL &
!$OMP PRIVATE(idx, actv, uq, vq, wq) &
!$OMP PRIVATE(Up0, Ue1, Uw1, Us1, Un1, Ub1, Ut1) &
!$OMP PRIVATE(Vp0, Ve1, Vw1, Vs1, Vn1, Vb1, Vt1) &
!$OMP PRIVATE(Wp0, We1, Ww1, Ws1, Wn1, Wb1, Wt1) &
!$OMP PRIVATE(b_e1, b_w1, b_n1, b_s1, b_t1, b_b1) &
!$OMP PRIVATE(w_e, w_w, w_n, w_s, w_t, w_b) &
!$OMP PRIVATE(r1, r2, r3, u1, u2, u3) &
!$OMP FIRSTPRIVATE(ix, jx, kx, h, u_ref, v_ref, w_ref)

!$OMP DO SCHEDULE(static)

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
      
      r1 = h*( (Wn1-Ws1) - (Vt1-Vb1) ) ! 12 flop
      r2 = h*( (Ut1-Ub1) - (We1-Ww1) )
      r3 = h*( (Ve1-Vw1) - (Un1-Us1) )
      
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

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do j=1,jx
        avr = avr + p(1,j,k)
      end do
      end do
!$OMP END DO

      avr = avr / rix
      
      
    case (X_plus)

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do j=1,jx
        avr = avr + p(ix,j,k)
      end do
      end do
!$OMP END DO

      avr = avr / rix
      
      
    case (Y_minus)

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do i=1,ix
        avr = avr + p(i,1,k)
      end do
      end do
!$OMP END DO

      avr = avr / rjx
      
      
    case (Y_plus)

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do i=1,ix
        avr = avr + p(i,jx,k)
      end do
      end do
!$OMP END DO

      avr = avr / rjx
      
      
    case (Z_minus)

!$OMP DO SCHEDULE(static)
      do j=1,jx
      do i=1,ix
        avr = avr + p(i,j,1)
      end do
      end do
!$OMP END DO

      avr = avr / rkx
      
    
    case (Z_plus)

!$OMP DO SCHEDULE(static)
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

!$OMP DO SCHEDULE(static)

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
  real                                                      ::  fx, fy, fz
  real                                                      ::  pp, dh, cf
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  p
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bid
  real, dimension(3)                                        ::  frc

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

!$OMP DO SCHEDULE(static)

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

  cf = 2.0 * dh * dh
  frc(1) = fx * cf
  frc(2) = fy * cf
  frc(3) = fz * cf

  return
  end subroutine force_compo
