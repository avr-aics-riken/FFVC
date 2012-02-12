!   *********************************************************
!
!   SPHERE - Skeleton for PHysical and Engineering REsearch
!
!   Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
!
!   *********************************************************
!
!> @file cbc_utility.f90
!> @brief Utilities
!> @author keno, FSI Team, VCAD, RIKEN
!<

!  ***************************************************************
!> @subroutine cbc_norm_v_div_max (ds, sz, g, div, coef, bp, flop)
!! @brief 速度成分の最大値を計算する
!! @param ds 最大値
!! @param sz 配列長
!! @param g ガイドセル長
!! @param div 速度の発散
!! @param coef 係数
!! @param bp BCindex P
!! @param flop
!<
    subroutine cbc_norm_v_div_max (ds, sz, g, div, coef, bp, flop)
    implicit none
    include 'cbc_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, g
    integer, dimension(3)                                     ::  sz
    real                                                      ::  ds, flop, r, coef
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  div
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    ds = 0.0

    flop = flop + real(ix)*real(jx)*real(kx)*5.0

    do k=1,kx
    do j=1,jx
    do i=1,ix
      r = div(i,j,k) * coef * real(ibits(bp(i,j,k), vld_cnvg, 1)) ! 有効セルの場合 1.0
      ds = max(ds, abs(r) )
    end do
    end do
    end do

    return
    end subroutine cbc_norm_v_div_max
    
!  *************************************************
!> @subroutine cbc_vmax (v_max, sz, g, v00, v, flop)
!! @brief 速度成分の最大値を計算する
!! @param v_max 最大値
!! @param sz 配列長
!! @param g ガイドセル長
!! @param v00 参照速度
!! @param v 速度ベクトル
!! @param flop
!<
    subroutine cbc_vmax (v_max, sz, g, v00, v, flop)
    implicit none
    integer                                                   ::  i, j, k, ix, jx, kx, g
    integer, dimension(3)                                     ::  sz
    real                                                      ::  vm1, vm2, vm3, v_max, flop
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v
    real, dimension(0:3)                                      ::  v00

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    vm1 = 0.0
    vm2 = 0.0
    vm3 = 0.0
    flop = flop + real(ix)*real(jx)*real(kx)*9.0 + 2.0

    do k=1,kx
    do j=1,jx
    do i=1,ix
      vm1 = max(vm1, abs(v(i,j,k,1)-v00(1) ) )
      vm2 = max(vm2, abs(v(i,j,k,2)-v00(2) ) )
      vm3 = max(vm3, abs(v(i,j,k,3)-v00(3) ) )
    end do
    end do
    end do
    
    v_max = max(vm1, vm2, vm3) ! maxss %xmm0, %xmm1, x 2 times > 2 flop

    return
    end subroutine cbc_vmax
    
!  ******************************************************
!> @subroutine cbc_i2vgt (q, sz, g, dh, v, bv, v00, flop)
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
    subroutine cbc_i2vgt (q, sz, g, dh, v, bv, v00, flop)
    implicit none
    include 'cbc_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, g, idx
    integer, dimension(3)                                     ::  sz
    integer                                                   ::  b_e1, b_w1, b_n1, b_s1, b_t1, b_b1
    real                                                      ::  u_ref, v_ref, w_ref, dh, h, flop, actv
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

    flop = flop + real(ix)*real(jx)*real(kx)*72.0 + 8.0
    ! flop = flop + real(ix)*real(jx)*real(kx)*72.0 + 13.0 ! DP

    do k=1,kx
    do j=1,jx
    do i=1,ix
      idx = bv(i,j,k)

      Up0 = v(i  ,j  ,k  ,1)
      Uw1 = v(i-1,j  ,k  ,1)
      Ue1 = v(i+1,j  ,k  ,1)
      Us1 = v(i  ,j-1,k  ,1)
      Un1 = v(i  ,j+1,k  ,1)
      Ub1 = v(i  ,j  ,k-1,1)
      Ut1 = v(i  ,j  ,k+1,1)
			
      Vp0 = v(i  ,j  ,k  ,2)
      Vw1 = v(i-1,j  ,k  ,2)
      Ve1 = v(i+1,j  ,k  ,2)
      Vs1 = v(i  ,j-1,k  ,2)
      Vn1 = v(i  ,j+1,k  ,2)
      Vb1 = v(i  ,j  ,k-1,2)
      Vt1 = v(i  ,j  ,k+1,2)

      Wp0 = v(i  ,j  ,k  ,3)
      Ww1 = v(i-1,j  ,k  ,3)
      We1 = v(i+1,j  ,k  ,3)
      Ws1 = v(i  ,j-1,k  ,3)
      Wn1 = v(i  ,j+1,k  ,3)
      Wb1 = v(i  ,j  ,k-1,3)
      Wt1 = v(i  ,j  ,k+1,3)
      
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

    return
    end subroutine cbc_i2vgt

!  ********************************************************
!> @subroutine cbc_rot_v (rot, sz, g, dh, v, bv, v00, flop)
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
    subroutine cbc_rot_v (rot, sz, g, dh, v, bv, v00, flop)
    implicit none
    include 'cbc_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, g, idx
    integer, dimension(3)                                     ::  sz
    integer                                                   ::  b_e1, b_w1, b_n1, b_s1, b_t1, b_b1
    real                                                      ::  u_ref, v_ref, w_ref, dh, h, flop, actv
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

    flop = flop + real(ix)*real(jx)*real(kx)*34.0 + 8.0
    ! flop = flop + real(ix)*real(jx)*real(kx)*34.0 + 13.0 ! DP

    do k=1,kx
    do j=1,jx
    do i=1,ix
      idx = bv(i,j,k)

      Up0 = v(i  ,j  ,k  ,1)
      Uw1 = v(i-1,j  ,k  ,1)
      Ue1 = v(i+1,j  ,k  ,1)
      Us1 = v(i  ,j-1,k  ,1)
      Un1 = v(i  ,j+1,k  ,1)
      Ub1 = v(i  ,j  ,k-1,1)
      Ut1 = v(i  ,j  ,k+1,1)
			
      Vp0 = v(i  ,j  ,k  ,2)
      Vw1 = v(i-1,j  ,k  ,2)
      Ve1 = v(i+1,j  ,k  ,2)
      Vs1 = v(i  ,j-1,k  ,2)
      Vn1 = v(i  ,j+1,k  ,2)
      Vb1 = v(i  ,j  ,k-1,2)
      Vt1 = v(i  ,j  ,k+1,2)

      Wp0 = v(i  ,j  ,k  ,3)
      Ww1 = v(i-1,j  ,k  ,3)
      We1 = v(i+1,j  ,k  ,3)
      Ws1 = v(i  ,j-1,k  ,3)
      Wn1 = v(i  ,j+1,k  ,3)
      Wb1 = v(i  ,j  ,k-1,3)
      Wt1 = v(i  ,j  ,k+1,3)
      
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

    return
    end subroutine cbc_rot_v

!  **********************************************************
!> @subroutine cbc_helicity (ht, sz, g, dh, v, bv, v00, flop)
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
    subroutine cbc_helicity (ht, sz, g, dh, v, bv, v00, flop)
    implicit none
    include 'cbc_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, g, idx
    integer, dimension(3)                                     ::  sz
    integer                                                   ::  b_e1, b_w1, b_n1, b_s1, b_t1, b_b1
    real                                                      ::  u_ref, v_ref, w_ref, dh, h, flop, actv
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

    flop = flop + real(ix)*real(jx)*real(kx)*40.0 + 8.0
    ! flop = flop + real(ix)*real(jx)*real(kx)*40.0 + 13.0 ! DP

    do k=1,kx
    do j=1,jx
    do i=1,ix
      idx = bv(i,j,k)

      Up0 = v(i  ,j  ,k  ,1)
      Uw1 = v(i-1,j  ,k  ,1)
      Ue1 = v(i+1,j  ,k  ,1)
      Us1 = v(i  ,j-1,k  ,1)
      Un1 = v(i  ,j+1,k  ,1)
      Ub1 = v(i  ,j  ,k-1,1)
      Ut1 = v(i  ,j  ,k+1,1)
			
      Vp0 = v(i  ,j  ,k  ,2)
      Vw1 = v(i-1,j  ,k  ,2)
      Ve1 = v(i+1,j  ,k  ,2)
      Vs1 = v(i  ,j-1,k  ,2)
      Vn1 = v(i  ,j+1,k  ,2)
      Vb1 = v(i  ,j  ,k-1,2)
      Vt1 = v(i  ,j  ,k+1,2)

      Wp0 = v(i  ,j  ,k  ,3)
      Ww1 = v(i-1,j  ,k  ,3)
      We1 = v(i+1,j  ,k  ,3)
      Ws1 = v(i  ,j-1,k  ,3)
      Wn1 = v(i  ,j+1,k  ,3)
      Wb1 = v(i  ,j  ,k-1,3)
      Wt1 = v(i  ,j  ,k+1,3)
      
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

    return
    end subroutine cbc_helicity
