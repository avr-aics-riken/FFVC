!********************************************************************
!
!   FFV : Frontflow / violet Cartesian
!
!   Copyright (c) 2012 All right reserved.
!
!   Institute of Industrial Science, The University of Tokyo, Japan. 
!
!********************************************************************

!> @file   ffv_utility.f90
!! @brief  Utilitiy functions
!! @author kero
!<

!> ********************************************************************
!! @brief 有効セルに対する発散の最大値と自乗和を計算，絶対値の最大値の位置を返す
!! @param ds 残差の絶対値
!! @param rm 残差の自乗和
!! @param idx ノルムの最大値の位置
!! @param sz 配列長
!! @param g ガイドセル長
!! @param div 発散値のベース
!! @param coef 係数
!! @param bp BCindex P
!! @param flop
!<
    subroutine norm_v_div_dbg (ds, rm, idx, sz, g, div, coef, bp, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, g, i0, j0, k0
    integer, dimension(3)                                     ::  sz, idx
    double precision                                          ::  flop
    real                                                      ::  ds, r, coef, rm, d
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
    rm = 0.0

    flop = flop + dble(ix)*dble(jx)*dble(kx)*5.0d0

!$OMP PARALLEL &
!$OMP PRIVATE(r, d) &
!$OMP FIRSTPRIVATE(ix, jx, kx, coef) &
!$OMP SHARED(ds, i0, j0, k0)

#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1) &
#elif defined _STATIC
!$OMP DO SCHEDULE(static) &
#else
!$OMP DO SCHEDULE(hoge)
#endif
!$OMP REDUCTION(+:rm)
    do k=1,kx
    do j=1,jx
    do i=1,ix
      r = div(i,j,k) * coef * real(ibits(bp(i,j,k), vld_cnvg, 1)) ! 有効セルの場合 1.0
      d = abs(r)
      rm = rm + r*r
      
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
!! @param ds 最大値
!! @param sz 配列長
!! @param g ガイドセル長
!! @param div 速度の発散
!! @param coef 係数
!! @param bp BCindex P
!! @param flop
!<
    subroutine norm_v_div_l2 (ds, sz, g, div, coef, bp, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, g
    integer, dimension(3)                                     ::  sz
    double precision                                          ::  flop
    real                                                      ::  ds, r, coef
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  div
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    ds = 0.0

    flop = flop + dble(ix)*dble(jx)*dble(kx)*5.0d0

!$OMP PARALLEL &
!$OMP PRIVATE(r) &
!$OMP FIRSTPRIVATE(ix, jx, kx, coef)

#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1) &
#elif defined _STATIC
!$OMP DO SCHEDULE(static) &
#else
!$OMP DO SCHEDULE(hoge)
#endif
!$OMP REDUCTION(+:ds) 
    do k=1,kx
    do j=1,jx
    do i=1,ix
      r = div(i,j,k) * coef * real(ibits(bp(i,j,k), vld_cnvg, 1)) ! 有効セルの場合 1.0
      ds = ds + r*r
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    return
    end subroutine norm_v_div_l2

!> ********************************************************************
!! @brief 速度成分の最大値を計算する
!! @param ds 最大値
!! @param sz 配列長
!! @param g ガイドセル長
!! @param div 速度の発散
!! @param coef 係数
!! @param bp BCindex P
!! @param flop
!<
    subroutine norm_v_div_max (ds, sz, g, div, coef, bp, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, g
    integer, dimension(3)                                     ::  sz
    double precision                                          ::  flop
    real                                                      ::  ds, r, coef
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  div
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    ds = 0.0

    flop = flop + dble(ix)*dble(jx)*dble(kx)*5.0d0

!$OMP PARALLEL &
!$OMP PRIVATE(r) &
!$OMP FIRSTPRIVATE(ix, jx, kx, coef)

#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1) &
#elif defined _STATIC
!$OMP DO SCHEDULE(static) &
#else
!$OMP DO SCHEDULE(hoge)
#endif
!$OMP REDUCTION(max:ds)
    do k=1,kx
    do j=1,jx
    do i=1,ix
      r = div(i,j,k) * coef * real(ibits(bp(i,j,k), vld_cnvg, 1)) ! 有効セルの場合 1.0
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
!! @param v_max 最大値
!! @param sz 配列長
!! @param g ガイドセル長
!! @param v00 参照速度
!! @param v 速度ベクトル
!! @param flop
!<
    subroutine find_vmax (v_max, sz, g, v00, v, flop)
    implicit none
    integer                                                   ::  i, j, k, ix, jx, kx, g
    integer, dimension(3)                                     ::  sz
    double precision                                          ::  flop
    real                                                      ::  vm1, vm2, vm3, v_max, vx, vy, vz
    real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  v
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
!$OMP FIRSTPRIVATE(ix, jx, kx, vx, vy, vz)

#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1) &
#elif defined _STATIC
!$OMP DO SCHEDULE(static) &
#else
!$OMP DO SCHEDULE(hoge)
#endif
!$OMP REDUCTION(max:vm1) &
!$OMP REDUCTION(max:vm2) &
!$OMP REDUCTION(max:vm3)
    do k=1,kx
    do j=1,jx
    do i=1,ix
      vm1 = max(vm1, abs(v(1,i,j,k)-vx ) )
      vm2 = max(vm2, abs(v(2,i,j,k)-vy ) )
      vm3 = max(vm3, abs(v(3,i,j,k)-vz ) )
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
!! @param[out] q 速度勾配テンソルの第２不変量
!! @param sz 配列長
!! @param g ガイドセル長
!! @param dh 格子幅
!! @param v セルセンター速度ベクトル
!! @param bv BCindex V
!! @param v00 参照速度
!! @param[out] flop
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
    real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  v
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

#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
    do k=1,kx
    do j=1,jx
    do i=1,ix
      idx = bv(i,j,k)

      Ub1 = v(1, i  ,j  ,k-1)
      Vb1 = v(2, i  ,j  ,k-1)
      Wb1 = v(3, i  ,j  ,k-1)

      Us1 = v(1, i  ,j-1,k  )
      Vs1 = v(2, i  ,j-1,k  )
      Ws1 = v(3, i  ,j-1,k  )

      Uw1 = v(1, i-1,j  ,k  )
      Vw1 = v(2, i-1,j  ,k  )
      Ww1 = v(3, i-1,j  ,k  )

      Up0 = v(1, i  ,j  ,k  )
      Vp0 = v(2, i  ,j  ,k  )
      Wp0 = v(3, i  ,j  ,k  )
      
      Ue1 = v(1, i+1,j  ,k  )
      Ve1 = v(2, i+1,j  ,k  )
      We1 = v(3, i+1,j  ,k  )
      
      Un1 = v(1, i  ,j+1,k  )
      Vn1 = v(2, i  ,j+1,k  )
      Wn1 = v(3, i  ,j+1,k  )

      Ut1 = v(1, i  ,j  ,k+1)
      Vt1 = v(2, i  ,j  ,k+1)
      Wt1 = v(3, i  ,j  ,k+1)
      
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
!! @brief 速度勾配テンソルの第２不変量の計算
!! @param[out] rot 渦度ベクトル
!! @param sz 配列長
!! @param g ガイドセル長
!! @param dh 格子幅
!! @param v セルセンター速度ベクトル
!! @param bv BCindex V
!! @param v00 参照速度
!! @param[out] flop
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
    real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  v, rot
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

#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
    do k=1,kx
    do j=1,jx
    do i=1,ix
      idx = bv(i,j,k)

      Ub1 = v(1, i  ,j  ,k-1)
      Vb1 = v(2, i  ,j  ,k-1)
      Wb1 = v(3, i  ,j  ,k-1)

      Us1 = v(1, i  ,j-1,k  )
      Vs1 = v(2, i  ,j-1,k  )
      Ws1 = v(3, i  ,j-1,k  )

      Uw1 = v(1, i-1,j  ,k  )
      Vw1 = v(2, i-1,j  ,k  )
      Ww1 = v(3, i-1,j  ,k  )

      Up0 = v(1, i  ,j  ,k  )
      Vp0 = v(2, i  ,j  ,k  )
      Wp0 = v(3, i  ,j  ,k  )

      Ue1 = v(1, i+1,j  ,k  )
      Ve1 = v(2, i+1,j  ,k  )
      We1 = v(3, i+1,j  ,k  )

      Un1 = v(1, i  ,j+1,k  )
      Vn1 = v(2, i  ,j+1,k  )
      Wn1 = v(3, i  ,j+1,k  )

      Ut1 = v(1, i  ,j  ,k+1)
      Vt1 = v(2, i  ,j  ,k+1)
      Wt1 = v(3, i  ,j  ,k+1)
      
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
      
      rot(1,i,j,k) = r1 * actv ! 3 flop
      rot(2,i,j,k) = r2 * actv
      rot(3,i,j,k) = r3 * actv
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
!! @param bv BCindex V
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
    real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  v
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

#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
    do k=1,kx
    do j=1,jx
    do i=1,ix
      idx = bv(i,j,k)

      Ub1 = v(1, i  ,j  ,k-1)
      Vb1 = v(2, i  ,j  ,k-1)
      Wb1 = v(3, i  ,j  ,k-1)

      Us1 = v(1, i  ,j-1,k  )
      Vs1 = v(2, i  ,j-1,k  )
      Ws1 = v(3, i  ,j-1,k  )

      Uw1 = v(1, i-1,j  ,k  )
      Vw1 = v(2, i-1,j  ,k  )
      Ww1 = v(3, i-1,j  ,k  )

      Up0 = v(1, i  ,j  ,k  )
      Vp0 = v(2, i  ,j  ,k  )
      Wp0 = v(3, i  ,j  ,k  )

      Ue1 = v(1, i+1,j  ,k  )
      Ve1 = v(2, i+1,j  ,k  )
      We1 = v(3, i+1,j  ,k  )

      Un1 = v(1, i  ,j+1,k  )
      Vn1 = v(2, i  ,j+1,k  )
      Wn1 = v(3, i  ,j+1,k  )

      Ut1 = v(1, i  ,j  ,k+1)
      Vt1 = v(2, i  ,j  ,k+1)
      Wt1 = v(3, i  ,j  ,k+1)
      
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
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do k=1,kx
      do j=1,jx
        avr = avr + p(1,j,k)
      end do
      end do
!$OMP END DO

      avr = avr / rix
      
      
    case (X_plus)
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do k=1,kx
      do j=1,jx
        avr = avr + p(ix,j,k)
      end do
      end do
!$OMP END DO

      avr = avr / rix
      
      
    case (Y_minus)
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do k=1,kx
      do i=1,ix
        avr = avr + p(i,1,k)
      end do
      end do
!$OMP END DO

      avr = avr / rjx
      
      
    case (Y_plus)
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do k=1,kx
      do i=1,ix
        avr = avr + p(i,jx,k)
      end do
      end do
!$OMP END DO

      avr = avr / rjx
      
      
    case (Z_minus)
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do j=1,jx
      do i=1,ix
        avr = avr + p(i,j,1)
      end do
      end do
!$OMP END DO

      avr = avr / rkx
      
    
    case (Z_plus)
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
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

#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif

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
!! @brief 物体表面の力を計算する
!! @param[out] frc 力の成分
!! @param sz 配列長
!! @param g ガイドセル長
!! @param p 圧力
!! @param bp BCindex P
!! @param dh 無次元格子幅
!! @param[out] flop flop count
!<
  subroutine force (frc, sz, g, p, bp, dh, flop)
  implicit none
  include 'ffv_f_params.h'
  integer                                                   ::  i, j, k, ix, jx, kx, g
  integer                                                   ::  bw, be, bs, bn, bb, bt
  integer, dimension(3)                                     ::  sz
  double precision                                          ::  flop
  real                                                      ::  fx, fy, fz
  real                                                      ::  qw, qe, qs, qn, qb, qt
  real                                                      ::  actv, pp, dh, cf
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  p
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp
  real, dimension(3)                                        ::  frc

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)

  flop = flop + dble(ix)*dble(jx)*dble(kx)*12.0d0 + 5.0d0

  fx = 0.0
  fy = 0.0
  fz = 0.0

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx) &
!$OMP PRIVATE(actv, pp) &
!$OMP PRIVATE(bw, be, bs, bn, bb, bt) &
!$OMP PRIVATE(qw, qe, qs, qn, qb, qt)

#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1) &
#elif defined _STATIC
!$OMP DO SCHEDULE(static) &
#else
!$OMP DO SCHEDULE(hoge)
#endif
!$OMP REDUCTION(+:fx) &
!$OMP REDUCTION(+:fy) &
!$OMP REDUCTION(+:fz)

  do k=1,kx
  do j=1,jx
  do i=1,ix

    ! Fluid -> 1.0
    actv= real(ibits(bp(i,j,k), State, 1))
    bw = ibits(bp(i-1, j  , k  ), State, 1)
    be = ibits(bp(i+1, j  , k  ), State, 1)
    bs = ibits(bp(i  , j-1, k  ), State, 1)
    bn = ibits(bp(i  , j+1, k  ), State, 1)
    bb = ibits(bp(i  , j  , k-1), State, 1)
    bt = ibits(bp(i  , j  , k+1), State, 1)


    ! if not cut -> q=0.0
    qw = 0.0
    qe = 0.0
    qs = 0.0
    qn = 0.0
    qb = 0.0
    qt = 0.0

    ! 参照先が固体である場合のみ
    if ( bw == 0 ) qw = 1.0
    if ( be == 0 ) qe = 1.0
    if ( bs == 0 ) qs = 1.0
    if ( bn == 0 ) qn = 1.0
    if ( bb == 0 ) qb = 1.0
    if ( bt == 0 ) qt = 1.0

    pp = p(i,j,k) * actv

    ! 各方向に壁がある場合、かつ流体セルのみ力を積算
    fx = fx + pp * ( qe - qw )
    fy = fy + pp * ( qn - qs )
    fz = fz + pp * ( qt - qb )

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
  end subroutine force
