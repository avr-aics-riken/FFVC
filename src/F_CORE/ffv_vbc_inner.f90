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

!> @file   ffv_bc_inner.f90
!! @brief  計算領域内部の境界条件
!! @author aics
!<


!> ********************************************************************
!! @brief 乱流の初期擾乱 Y方向の壁面
!! @param [out] v       速度場
!! @param [in]  sz      配列長
!! @param [in]  g       ガイドセル長
!! @param [in]  dh      格子幅
!! @param [in]  origin  サブドメインの領域基点
!! @param [in]  width   代表流路幅
!! @param [in]  Re_tau  レイノルズ数
!! @param [in]  Ubar    Bulk velocity
!! @param [in]  nu      動粘性係数
!! @param [in]  mode    cell center(1) or face(2)
!!<
subroutine perturb_u_y (v, sz, g, dh, origin, width, Re_tau, Ubar, nu, mode)
implicit none
include 'ffv_f_params.h'
integer                                                   :: ix, jx, kx, i, j, k, g, mode
integer, dimension(3)                                     :: sz
real                                                      :: duplus, deviation, eps, nu
real                                                      :: width, Re_tau, u_tau, Ubar, tmp1, tmp2, rhalf
real                                                      :: half, dh, xc, yc, zc, xplus, yplus, zplus, perturb
real, dimension(3)                                        :: origin
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) :: v
real, parameter :: sigmaC     = 0.00055d0,                           & !
                   duplusC    = 0.25d0,                              & !
                   epsilonC   = 0.05d0,                              & !
                   deviationC = 0.2d0,                               & !
                   betaPlusC  = 2.0d0 * 3.14d0 * (1.0d0 / 200.0d0),  & !
                   alphaPlusC = 2.0d0 * 3.14d0 * (1.0d0 / 500.0d0)     !


half      = width * 0.5d0          ! channel half width
u_tau     = Re_tau * nu / half     ! friction velocity
duplus    = Ubar * duplusC / u_tau
eps       = Ubar * epsilonC
rhalf     = Re_tau / half

ix = sz(1)
jx = sz(2)
kx = sz(3)


!$OMP PARALLEL &
!$OMP PRIVATE(xc, yc, zc, i, j, k, tmp1, tmp2, xplus, yplus, zplus) &
!$OMP PRIVATE(perturb, deviation) &
!$OMP FIRSTPRIVATE(ix, jx, kx, g, dh, half, rhalf, width, u_tau, duplus, eps, mode)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k = 1-g, kx+g
do j = 1-g, jx+g
do i = 1-g, ix+g

  if ( mode == 1 ) then
    xc = origin(1) + (i - 0.5d0)*dh
  else
    xc = origin(1) + i * dh
  end if

  yc = origin(2) + (j - 0.5d0)*dh
  zc = origin(3) + (k - 0.5d0)*dh

! 壁面からの距離
  if ( yc .gt. half ) then
    yc = width - yc
  end if

  xplus = xc * rhalf
  yplus = yc * rhalf
  zplus = zc * rhalf

  perturb = 0.0
  call RANDOM_NUMBER(perturb)

  deviation = 1.0d0 + deviationC * perturb

  tmp1 = -sigmaC * yplus * yplus
  tmp2 = yc / half

  v(i, j, k, 1) = 3.0d0 * Ubar * (tmp2 - 0.5d0 * tmp2*tmp2) &
                + u_tau * duplus * 0.5d0 * (yplus/40.0d0)   &
                * exp(tmp1 + 0.5d0) * cos(betaPlusC * zplus) * deviation
  v(i, j, k, 2) = 0.0d0
  v(i, j, k, 3) = eps * sin(alphaPlusC * xplus) * yplus * exp(tmp1) * deviation
end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL

return

end subroutine perturb_u_y


!> ********************************************************************
!! @brief 乱流の初期擾乱 Z方向
!! @param [out] v       セルセンター速度場
!! @param [in]  sz      配列長
!! @param [in]  g       ガイドセル長
!! @param [in]  dh      格子幅
!! @param [in]  origin  サブドメインの領域基点
!! @param [in]  width   代表流路幅
!! @param [in]  Re_tau  レイノルズ数
!! @param [in]  Ubar    Bulk velocity
!! @param [in]  nu      動粘性係数
!! @param [in]  mode    cell center(1) or face(2)
!!<
subroutine perturb_u_z (v, sz, g, dh, origin, width, Re_tau, Ubar, nu, mode)
implicit none
include 'ffv_f_params.h'
integer                                                   :: ix, jx, kx, i, j, k, g, mode
integer, dimension(3)                                     :: sz
real                                                      :: duplus, deviation, eps, nu
real                                                      :: width, Re_tau, u_tau, Ubar, tmp1, tmp2, rhalf
real                                                      :: half, dh, xc, yc, zc, xplus, yplus, zplus, perturb
real, dimension(3)                                        :: origin
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) :: v
real, parameter :: sigmaC     = 0.00055d0,                           &
                   duplusC    = 0.25d0,                              &
                   epsilonC   = 0.05d0,                              &
                   deviationC = 0.2d0,                               &
                   betaPlusC  = 2.0d0 * 3.14d0 * (1.0d0 / 200.0d0),  &
                   alphaPlusC = 2.0d0 * 3.14d0 * (1.0d0 / 500.0d0)


half      = width * 0.5d0          ! channel half width
u_tau     = Re_tau * nu / half
duplus    = Ubar * duplusC / u_tau
eps       = Ubar * epsilonC
rhalf     = Re_tau / half

ix = sz(1)
jx = sz(2)
kx = sz(3)


!$OMP PARALLEL &
!$OMP PRIVATE(xc, yc, zc, i, j, k, tmp1, tmp2, xplus, yplus, zplus) &
!$OMP PRIVATE(perturb, deviation) &
!$OMP FIRSTPRIVATE(ix, jx, kx, g, dh, half, rhalf, width, u_tau, duplus, eps, mode)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k = 1-g, kx+g
do j = 1-g, jx+g
do i = 1-g, ix+g

  if ( mode == 1 ) then
    xc = origin(1) + (i - 0.5d0)*dh
  else
    xc = origin(1) + i * dh
  end if

  yc = origin(2) + (j - 0.5d0)*dh
  zc = origin(3) + (k - 0.5d0)*dh

  ! 壁面からの距離
  if ( zc .gt. half ) then
    zc = width - zc
  end if

  xplus = xc * rhalf
  yplus = yc * rhalf
  zplus = zc * rhalf

  perturb = 0.0
  call RANDOM_NUMBER(perturb)

  deviation = 1.0d0 + deviationC * perturb

  tmp1 = -sigmaC * zplus * zplus
  tmp2 = zc / half

  v(i, j, k, 1) = 3.0d0 * Ubar * (tmp2 - 0.5d0 * tmp2*tmp2) &
                + u_tau * duplus * 0.5d0 * (zplus/40.0d0)   &
                * exp(tmp1 + 0.5d0) * cos(betaPlusC * yplus) * deviation
  v(i, j, k, 2) = eps * sin(alphaPlusC * xplus) * zplus * exp(tmp1) * deviation
  v(i, j, k, 3) = 0.0d0
end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL

return

end subroutine perturb_u_z



!> ********************************************************************
!! @brief 内部速度境界条件による対流項と粘性項の流束の修正
!! @param [out] wv   疑似ベクトルの空間項の評価値
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル長
!! @param [in]  st   ループの開始インデクス
!! @param [in]  ed   ループの終了インデクス
!! @param [in]  dh   格子幅
!! @param [in]  rei  Reynolds数の逆数
!! @param [in]  v    セルセンター速度ベクトル（u^n）
!! @param [in]  bv   BCindex C
!! @param [in]  odr  内部境界処理時の速度境界条件のエントリ
!! @param [in]  vec  指定する速度ベクトル
!! @param [out] flop 浮動小数点演算数
!! @note vecには，流出条件のとき対流流出速度
!! @todo 内部と外部の分離 do loopの内側に条件分岐を入れているので修正
!! @todo 流出境界はローカルの流束となるように変更する（外部境界参照）
!<
    subroutine pvec_ibc_oflow (wv, sz, g, st, ed, dh, rei, v, bv, odr, vec, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, g, bvx, odr, is, ie, js, je, ks, ke
    integer, dimension(3)                                     ::  sz, st, ed
    double precision                                          ::  flop
    real                                                      ::  Up0, Ue1, Uw1, Us1, Un1, Ub1, Ut1
    real                                                      ::  Vp0, Ve1, Vw1, Vs1, Vn1, Vb1, Vt1
    real                                                      ::  Wp0, We1, Ww1, Ws1, Wn1, Wb1, Wt1
    real                                                      ::  c_e, c_w, c_n, c_s, c_t, c_b
    real                                                      ::  dh, dh1, dh2, EX, EY, EZ, rei
    real                                                      ::  cnv_u, cnv_v, cnv_w, cr, cl, m
    real                                                      ::  u_bc, v_bc, w_bc
    real                                                      ::  fu_r, fu_l, fv_r, fv_l, fw_r, fw_l
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v, wv
    real, dimension(3)                                        ::  vec
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bv
    
    dh1= 1.0/dh
    dh2= rei*dh1*dh1
    
    ! u_bcは境界速度
    u_bc = vec(1)
    v_bc = vec(2)
    w_bc = vec(3)
    
    is = st(1)
    ie = ed(1)
    js = st(2)
    je = ed(2)
    ks = st(3)
    ke = ed(3)
    
    flop = flop + 10.0d0 ! DP 15 flops
    
    m = 0.0

!$OMP PARALLEL REDUCTION(+:m) &
!$OMP FIRSTPRIVATE(is, ie, js, je, ks, ke, u_bc, v_bc, w_bc, odr) &
!$OMP FIRSTPRIVATE(dh1, dh2) &
!$OMP PRIVATE(bvx, cnv_u, cnv_v, cnv_w, cr, cl) &
!$OMP PRIVATE(Up0, Ue1, Uw1, Us1, Un1, Ub1, Ut1) &
!$OMP PRIVATE(Vp0, Ve1, Vw1, Vs1, Vn1, Vb1, Vt1) &
!$OMP PRIVATE(Wp0, We1, Ww1, Ws1, Wn1, Wb1, Wt1) &
!$OMP PRIVATE(c_e, c_w, c_n, c_s, c_t, c_b) &
!$OMP PRIVATE(fu_r, fu_l, fv_r, fv_l, fw_r, fw_l) &
!$OMP PRIVATE(EX, EY, EZ)

!$OMP DO SCHEDULE(static)

    do k=ks,ke
    do j=js,je
    do i=is,ie
      bvx = bv(i,j,k)

      if ( 0 /= iand(bvx, bc_mask30) ) then ! 6面のうちのどれか速度境界フラグが立っている場合
        cnv_u = 0.0
        cnv_v = 0.0
        cnv_w = 0.0
        
        ! 変数のロード
        Up0 = v(i,j,k,1)
        Vp0 = v(i,j,k,2)
        Wp0 = v(i,j,k,3)
      
        Uw1 = 0.0
        Ue1 = 0.0
        Us1 = 0.0
        Un1 = 0.0
        Ub1 = 0.0
        Ut1 = 0.0
        Vw1 = 0.0
        Ve1 = 0.0
        Vs1 = 0.0
        Vn1 = 0.0
        Vb1 = 0.0
        Vt1 = 0.0
        Ww1 = 0.0
        We1 = 0.0
        Ws1 = 0.0
        Wn1 = 0.0
        Wb1 = 0.0
        Wt1 = 0.0
      
        ! 内部境界のときの各面のBCフラグ ibits() = 0(Normal) / others(BC) >> c_e = 0.0(Normal) / 1.0(BC) 
        c_e = 0.0
        c_w = 0.0
        c_n = 0.0
        c_s = 0.0
        c_t = 0.0
        c_b = 0.0
        if ( ibits(bvx, bc_face_E, bitw_5) == odr ) c_e = 1.0
        if ( ibits(bvx, bc_face_W, bitw_5) == odr ) c_w = 1.0
        if ( ibits(bvx, bc_face_N, bitw_5) == odr ) c_n = 1.0
        if ( ibits(bvx, bc_face_S, bitw_5) == odr ) c_s = 1.0
        if ( ibits(bvx, bc_face_T, bitw_5) == odr ) c_t = 1.0
        if ( ibits(bvx, bc_face_B, bitw_5) == odr ) c_b = 1.0
			
        ! X方向 ---------------------------------------
        if ( c_w == 1.0 ) then
          Uw1  = Up0
          Vw1  = Vp0
          Ww1  = Wp0
          cl   = u_bc
          if ( cl>0.0 ) cl=0.0
          fu_l = cl*Up0
          fv_l = cl*Vp0
          fw_l = cl*Wp0 ! 3 flops
        end if
        
        if ( c_e == 1.0 ) then
          Ue1  = Up0
          Ve1  = Vp0
          We1  = Wp0
          cr   = u_bc
          if ( cr<0.0 ) cr=0.0
          fu_r = cr*Up0
          fv_r = cr*Vp0
          fw_r = cr*Wp0 ! 3 flops
        end if
        
        cnv_u = cnv_u + fu_r*c_e - fu_l*c_w
        cnv_v = cnv_v + fv_r*c_e - fv_l*c_w
        cnv_w = cnv_w + fw_r*c_e - fw_l*c_w ! 12 flops
			
        ! Y方向 ---------------------------------------
        if ( c_s == 1.0 ) then
          Us1  = Up0
          Vs1  = Vp0
          Ws1  = Wp0
          cl   = v_bc
          if ( cl>0.0 ) cl=0.0
          fu_l = cl*Up0
          fv_l = cl*Vp0
          fw_l = cl*Wp0
        end if
        
        if ( c_n == 1.0 ) then
          Un1  = Up0
          Vn1  = Vp0
          Wn1  = Wp0
          cr   = v_bc
          if ( cr<0.0 ) cr=0.0
          fu_r = cr*Up0
          fv_r = cr*Vp0
          fw_r = cr*Wp0
        end if
        
        cnv_u = cnv_u + fu_r*c_n - fu_l*c_s
        cnv_v = cnv_v + fv_r*c_n - fv_l*c_s
        cnv_w = cnv_w + fw_r*c_n - fw_l*c_s
			
        ! Z方向 ---------------------------------------
        if ( c_b == 1.0 ) then
          Ub1  = Up0
          Vb1  = Vp0
          Wb1  = Wp0
          cl   = w_bc
          if ( cl>0.0 ) cl=0.0
          fu_l = cl*Up0
          fv_l = cl*Vp0
          fw_l = cl*Wp0
        end if

        if ( c_t == 1.0 ) then
          Ut1  = Up0
          Vt1  = Vp0
          Wt1  = Wp0
          cr   = w_bc
          if ( cr<0.0 ) cr=0.0
          fu_r = cr*Up0
          fv_r = cr*Vp0
          fw_r = cr*Wp0
        end if
        
        cnv_u = cnv_u + fu_r*c_t - fu_l*c_b
        cnv_v = cnv_v + fv_r*c_t - fv_l*c_b
        cnv_w = cnv_w + fw_r*c_t - fw_l*c_b
      
        ! 粘性項
        EX = ( Ue1 - Up0 ) * c_e &
           + ( Uw1 - Up0 ) * c_w &
           + ( Un1 - Up0 ) * c_n &
           + ( Us1 - Up0 ) * c_s &
           + ( Ut1 - Up0 ) * c_t &
           + ( Ub1 - Up0 ) * c_b ! 17 flops
        EY = ( Ve1 - Vp0 ) * c_e &
           + ( Vw1 - Vp0 ) * c_w &
           + ( Vn1 - Vp0 ) * c_n &
           + ( Vs1 - Vp0 ) * c_s &
           + ( Vt1 - Vp0 ) * c_t &
           + ( Vb1 - Vp0 ) * c_b
        EZ = ( We1 - Wp0 ) * c_e &
           + ( Ww1 - Wp0 ) * c_w &
           + ( Wn1 - Wp0 ) * c_n &
           + ( Ws1 - Wp0 ) * c_s &
           + ( Wt1 - Wp0 ) * c_t &
           + ( Wb1 - Wp0 ) * c_b

        wv(i,j,k,1) = wv(i,j,k,1) + ( -cnv_u*dh1 + EX*dh2 )
        wv(i,j,k,2) = wv(i,j,k,2) + ( -cnv_v*dh1 + EY*dh2 )
        wv(i,j,k,3) = wv(i,j,k,3) + ( -cnv_w*dh1 + EZ*dh2 ) ! 4*3 = 12 flops
        m = m + 1.0
      endif
    end do
    end do
    end do
    
!$OMP END DO
!$OMP END PARALLEL

    ! loop :  (3+3+12)*3dir + 17*3 + 12 = 117 flops

    flop = flop + m*117.0
    
    return
    end subroutine pvec_ibc_oflow

!> ********************************************************************
!! @brief 内部速度境界条件による対流項と粘性項の流束の修正（有限体積型のスキーム）
!! @param [out] wv   疑似ベクトルの空間項の評価値
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル長
!! @param [in]  st   ループの開始インデクス
!! @param [in]  ed   ループの終了インデクス
!! @param [in]  dh   格子幅
!! @param [in]  v00  参照速度
!! @param [in]  rei  Reynolds数の逆数
!! @param [in]  v    セルセンター速度ベクトル（u^n）
!! @param [in]  bv   BCindex C
!! @param [in]  odr  内部境界処理時の速度境界条件のエントリ
!! @param [in]  vec  指定する速度ベクトル
!! @param [out] flop 浮動小数点演算数
!! @note vecには，流入条件のとき指定速度，流出条件のとき対流流出速度，カット位置に関わらず指定速度で流束を計算
!! @todo 流出境界はローカルの流束となるように変更する（外部境界参照）
!<
    subroutine pvec_ibc_specv_fvm (wv, sz, g, st, ed, dh, v00, rei, v, bv, odr, vec, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                     ::  i, j, k, g, bvx, odr, is, ie, js, je, ks, ke
    integer, dimension(3)                                       ::  sz, st, ed
    double precision                                            ::  flop
    real                                                        ::  Up0, Ue1, Uw1, Us1, Un1, Ub1, Ut1
    real                                                        ::  Vp0, Ve1, Vw1, Vs1, Vn1, Vb1, Vt1
    real                                                        ::  Wp0, We1, Ww1, Ws1, Wn1, Wb1, Wt1
    real                                                        ::  c_e, c_w, c_n, c_s, c_t, c_b
    real                                                        ::  dh, dh1, dh2, EX, EY, EZ, rei
    real                                                        ::  u_ref, v_ref, w_ref, m1, m2
    real                                                        ::  cnv_u, cnv_v, cnv_w, cr, cl, acr, acl
    real                                                        ::  u_bc, v_bc, w_bc, u_bc_ref, v_bc_ref, w_bc_ref
    real                                                        ::  fu_r, fu_l, fv_r, fv_l, fw_r, fw_l
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)   ::  v, wv
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bv
    real, dimension(0:3)                                        ::  v00
    real, dimension(3)                                          ::  vec
    
    is = st(1)
    ie = ed(1)
    js = st(2)
    je = ed(2)
    ks = st(3)
    ke = ed(3)

    dh1= 1.0/dh
    dh2= rei*dh1*dh1

    ! 参照座標系の速度
    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)
    
    ! u_bcは境界速度
    u_bc = vec(1)
    v_bc = vec(2)
    w_bc = vec(3)
    
    ! u_bc_refは参照座標系での境界速度
    u_bc_ref = u_bc + u_ref
    v_bc_ref = v_bc + v_ref
    w_bc_ref = w_bc + w_ref
    
    m1 = 0.0
    m2 = 0.0
    
    flop = flop + 13.0d0 ! DP 18 flop

!$OMP PARALLEL &
!$OMP REDUCTION(+:m1) &
!$OMP REDUCTION(+:m2) &
!$OMP FIRSTPRIVATE(is, ie, js, je, ks, ke, u_bc, v_bc, w_bc, u_bc_ref, v_bc_ref, w_bc_ref) &
!$OMP FIRSTPRIVATE(dh1, dh2, odr) &
!$OMP PRIVATE(bvx, cnv_u, cnv_v, cnv_w, EX, EY, EZ, cr, cl, acr, acl) &
!$OMP PRIVATE(Up0, Ue1, Uw1, Us1, Un1, Ub1, Ut1) &
!$OMP PRIVATE(Vp0, Ve1, Vw1, Vs1, Vn1, Vb1, Vt1) &
!$OMP PRIVATE(Wp0, We1, Ww1, Ws1, Wn1, Wb1, Wt1) &
!$OMP PRIVATE(c_e, c_w, c_n, c_s, c_t, c_b) &
!$OMP PRIVATE(fu_r, fu_l, fv_r, fv_l, fw_r, fw_l)

!$OMP DO SCHEDULE(static)
    
    do k=ks,ke
    do j=js,je
    do i=is,ie
      bvx = bv(i,j,k)

      if ( 0 /= iand(bvx, bc_mask30) ) then ! 6面のうちのどれか速度境界フラグが立っている場合
        cnv_u = 0.0
        cnv_v = 0.0
        cnv_w = 0.0
        
        ! 変数のロード
        Up0 = v(i,j,k,1)
        Vp0 = v(i,j,k,2)
        Wp0 = v(i,j,k,3)

        Uw1 = 0.0
        Ue1 = 0.0
        Us1 = 0.0
        Un1 = 0.0
        Ub1 = 0.0
        Ut1 = 0.0
        Vw1 = 0.0
        Ve1 = 0.0
        Vs1 = 0.0
        Vn1 = 0.0
        Vb1 = 0.0
        Vt1 = 0.0
        Ww1 = 0.0
        We1 = 0.0
        Ws1 = 0.0
        Wn1 = 0.0
        Wb1 = 0.0
        Wt1 = 0.0
      
        ! 内部境界のときの各面のBCフラグ ibits() = 0(Normal) / others(BC) >> c_e = 0.0(Normal) / 1.0(BC) 
        c_e = 0.0
        c_w = 0.0
        c_n = 0.0
        c_s = 0.0
        c_t = 0.0
        c_b = 0.0
        if ( ibits(bvx, bc_face_E, bitw_5) == odr ) c_e = 1.0
        if ( ibits(bvx, bc_face_W, bitw_5) == odr ) c_w = 1.0
        if ( ibits(bvx, bc_face_N, bitw_5) == odr ) c_n = 1.0
        if ( ibits(bvx, bc_face_S, bitw_5) == odr ) c_s = 1.0
        if ( ibits(bvx, bc_face_T, bitw_5) == odr ) c_t = 1.0
        if ( ibits(bvx, bc_face_B, bitw_5) == odr ) c_b = 1.0
      
        ! X方向 ---------------------------------------
        if ( c_w == 1.0 ) then
          Uw1  = u_bc_ref
          Vw1  = v_bc_ref
          Ww1  = w_bc_ref
          cl   = u_bc
          acl  = abs(cl)
          fu_l = 0.5*(cl*(Up0+Uw1) - acl*(Up0-Uw1))
          fv_l = 0.5*(cl*(Vp0+Vw1) - acl*(Vp0-Vw1))
          fw_l = 0.5*(cl*(Wp0+Ww1) - acl*(Wp0-Ww1)) ! 18+1 flops
          m2 = m2 + 1.0
        end if
        
        if ( c_e == 1.0 ) then
          Ue1  = u_bc_ref
          Ve1  = v_bc_ref
          We1  = w_bc_ref
          cr   = u_bc
          acr  = abs(cr)
          fu_r = 0.5*(cr*(Ue1+Up0) - acr*(Ue1-Up0))
          fv_r = 0.5*(cr*(Ve1+Vp0) - acr*(Ve1-Vp0))
          fw_r = 0.5*(cr*(We1+Wp0) - acr*(We1-Wp0)) ! 18+1 flops
          m2 = m2 + 1.0
        end if
        
        cnv_u = cnv_u + fu_r*c_e - fu_l*c_w
        cnv_v = cnv_v + fv_r*c_e - fv_l*c_w
        cnv_w = cnv_w + fw_r*c_e - fw_l*c_w ! 12 flops
			
        ! Y方向 ---------------------------------------
        if ( c_s == 1.0 ) then
          Us1  = u_bc_ref
          Vs1  = v_bc_ref
          Ws1  = w_bc_ref
          cl   = v_bc
          acl  = abs(cl)
          fu_l = 0.5*(cl*(Up0+Us1) - acl*(Up0-Us1))
          fv_l = 0.5*(cl*(Vp0+Vs1) - acl*(Vp0-Vs1))
          fw_l = 0.5*(cl*(Wp0+Ws1) - acl*(Wp0-Ws1))
          m2 = m2 + 1.0
        end if
        
        if ( c_n == 1.0 ) then
          Un1  = u_bc_ref
          Vn1  = v_bc_ref
          Wn1  = w_bc_ref
          cr   = v_bc
          acr  = abs(cr)
          fu_r = 0.5*(cr*(Un1+Up0) - acr*(Un1-Up0))
          fv_r = 0.5*(cr*(Vn1+Vp0) - acr*(Vn1-Vp0))
          fw_r = 0.5*(cr*(Wn1+Wp0) - acr*(Wn1-Wp0))
          m2 = m2 + 1.0
        end if
        
        cnv_u = cnv_u + fu_r*c_n - fu_l*c_s
        cnv_v = cnv_v + fv_r*c_n - fv_l*c_s
        cnv_w = cnv_w + fw_r*c_n - fw_l*c_s
			
        ! Z方向 ---------------------------------------
        if ( c_b == 1.0 ) then
          Ub1  = u_bc_ref
          Vb1  = v_bc_ref
          Wb1  = w_bc_ref
          cl   = w_bc
          acl  = abs(cl)
          fu_l = 0.5*(cl*(Up0+Ub1) - acl*(Up0-Ub1))
          fv_l = 0.5*(cl*(Vp0+Vb1) - acl*(Vp0-Vb1))
          fw_l = 0.5*(cl*(Wp0+Wb1) - acl*(Wp0-Wb1))
          m2 = m2 + 1.0
        end if

        if ( c_t == 1.0 ) then
          Ut1  = u_bc_ref
          Vt1  = v_bc_ref
          Wt1  = w_bc_ref
          cr   = w_bc
          acr  = abs(cr)
          fu_r = 0.5*(cr*(Ut1+Up0) - acr*(Ut1-Up0))
          fv_r = 0.5*(cr*(Vt1+Vp0) - acr*(Vt1-Vp0))
          fw_r = 0.5*(cr*(Wt1+Wp0) - acr*(Wt1-Wp0))
          m2 = m2 + 1.0
        end if
        
        cnv_u = cnv_u + fu_r*c_t - fu_l*c_b
        cnv_v = cnv_v + fv_r*c_t - fv_l*c_b
        cnv_w = cnv_w + fw_r*c_t - fw_l*c_b
      
        ! 粘性項
        EX = ( Ue1 - Up0 ) * c_e &
           + ( Uw1 - Up0 ) * c_w &
           + ( Un1 - Up0 ) * c_n &
           + ( Us1 - Up0 ) * c_s &
           + ( Ut1 - Up0 ) * c_t &
           + ( Ub1 - Up0 ) * c_b ! 17 flops
        EY = ( Ve1 - Vp0 ) * c_e &
           + ( Vw1 - Vp0 ) * c_w &
           + ( Vn1 - Vp0 ) * c_n &
           + ( Vs1 - Vp0 ) * c_s &
           + ( Vt1 - Vp0 ) * c_t &
           + ( Vb1 - Vp0 ) * c_b
        EZ = ( We1 - Wp0 ) * c_e &
           + ( Ww1 - Wp0 ) * c_w &
           + ( Wn1 - Wp0 ) * c_n &
           + ( Ws1 - Wp0 ) * c_s &
           + ( Wt1 - Wp0 ) * c_t &
           + ( Wb1 - Wp0 ) * c_b

        wv(i,j,k,1) = wv(i,j,k,1) + ( -cnv_u*dh1 + EX*dh2 )
        wv(i,j,k,2) = wv(i,j,k,2) + ( -cnv_v*dh1 + EY*dh2 )
        wv(i,j,k,3) = wv(i,j,k,3) + ( -cnv_w*dh1 + EZ*dh2 ) ! 4*3 = 12 flops
        m1 = m1 + 1.0
        
      endif
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    ! loop :  (12*3 + 17*3 + 12)*m1 + 19*m2 = 99*m1 + 19*m2

    flop = flop + m1*99.0 + m2*19.0

    return
    end subroutine pvec_ibc_specv_fvm



!> ********************************************************************
!! @brief 内部速度境界条件による対流項と粘性項の流束の修正（有限差分型のスキーム）
!! @param [out] wv   疑似ベクトルの空間項の評価値
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル長
!! @param [in]  st   ループの開始インデクス
!! @param [in]  ed   ループの終了インデクス
!! @param [in]  dh   格子幅
!! @param [in]  v00  参照速度
!! @param [in]  rei  Reynolds数の逆数
!! @param [in]  v    セルセンター速度ベクトル（u^n）
!! @param [in]  bv   BCindex C
!! @param [in]  odr  内部境界処理時の速度境界条件のエントリ
!! @param [in]  vec  指定する速度ベクトル
!! @param [out] flop 浮動小数点演算数
!! @note vecには，流入条件のとき指定速度，流出条件のとき対流流出速度，カット位置に関わらず指定速度で流束を計算
!! @todo 流出境界はローカルの流束となるように変更する（外部境界参照）
!!
!! wv(i,j,k) = -frac{1}{2} ( u {\frac{\partial u}{\partial x}}_R + u {\frac{\partial u}{\partial x}}_L )
!!             + ( {\frac{\partial u}{\partial x}}_R - {\frac{\partial u}{\partial x}}_L ) \frac{1}{Re h}
!!
!<
subroutine pvec_ibc_specv_fdm (wv, sz, g, st, ed, dh, v00, rei, v, bv, odr, vec, flop)
implicit none
include 'ffv_f_params.h'
integer                                                     ::  i, j, k, g, bvx, odr, is, ie, js, je, ks, ke
integer, dimension(3)                                       ::  sz, st, ed
double precision                                            ::  flop
real                                                        ::  c_e, c_w, c_n, c_s, c_t, c_b
real                                                        ::  dh, dh1, dh2, rei, Up, Vp, Wp
real                                                        ::  u_ref, v_ref, w_ref, m2, gu, gv, gw
real                                                        ::  u_bc, v_bc, w_bc, u_bc_ref, v_bc_ref, w_bc_ref
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)   ::  v, wv
integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bv
real, dimension(0:3)                                        ::  v00
real, dimension(3)                                          ::  vec

is = st(1)
ie = ed(1)
js = st(2)
je = ed(2)
ks = st(3)
ke = ed(3)

dh1= 1.0/dh
dh2= rei*dh1

! 参照座標系の速度
u_ref = v00(1)
v_ref = v00(2)
w_ref = v00(3)

! u_bcは境界速度
u_bc = vec(1)
v_bc = vec(2)
w_bc = vec(3)

! u_bc_refは参照座標系での境界速度
u_bc_ref = u_bc + u_ref
v_bc_ref = v_bc + v_ref
w_bc_ref = w_bc + w_ref

m2 = 0.0

flop = flop + 13.0d0 ! DP 18 flop

!$OMP PARALLEL &
!$OMP REDUCTION(+:m2) &
!$OMP FIRSTPRIVATE(is, ie, js, je, ks, ke, u_bc, v_bc, w_bc, u_bc_ref, v_bc_ref, w_bc_ref) &
!$OMP FIRSTPRIVATE(dh1, dh2, odr) &
!$OMP PRIVATE(bvx, Up, Vp, Wp, gu, gv, gw, i, j, k) &
!$OMP PRIVATE(c_e, c_w, c_n, c_s, c_t, c_b)

!$OMP DO SCHEDULE(static) COLLAPSE(3)

do k=ks,ke
do j=js,je
do i=is,ie
bvx = bv(i,j,k)

if ( 0 /= iand(bvx, bc_mask30) ) then ! 6面のうちのどれか速度境界フラグが立っている場合

! 変数のロード
Up = v(i,j,k,1)
Vp = v(i,j,k,2)
Wp = v(i,j,k,3)

! 内部境界のときの各面のBCフラグ ibits() = 0(Normal) / others(BC) >> c_e = 0.0(Normal) / 1.0(BC)
c_e = 0.0
c_w = 0.0
c_n = 0.0
c_s = 0.0
c_t = 0.0
c_b = 0.0
if ( ibits(bvx, bc_face_E, bitw_5) == odr ) c_e = 1.0
if ( ibits(bvx, bc_face_W, bitw_5) == odr ) c_w = 1.0
if ( ibits(bvx, bc_face_N, bitw_5) == odr ) c_n = 1.0
if ( ibits(bvx, bc_face_S, bitw_5) == odr ) c_s = 1.0
if ( ibits(bvx, bc_face_T, bitw_5) == odr ) c_t = 1.0
if ( ibits(bvx, bc_face_B, bitw_5) == odr ) c_b = 1.0


! X方向 ---------------------------------------
if ( c_w == 1.0 ) then
gu = 2.0 * (Up - u_bc_ref) * dh1
gv = 2.0 * (Vp - v_bc_ref) * dh1
gw = 2.0 * (Wp - w_bc_ref) * dh1
wv(i,j,k,1) = wv(i,j,k,1) - ( 0.5 * u_bc * gu + gu * dh2 )
wv(i,j,k,2) = wv(i,j,k,2) - ( 0.5 * u_bc * gv + gv * dh2 )
wv(i,j,k,3) = wv(i,j,k,3) - ( 0.5 * u_bc * gw + gw * dh2 )
m2 = m2 + 1.0
end if

if ( c_e == 1.0 ) then
gu = 2.0 * (u_bc_ref - Up) * dh1
gv = 2.0 * (v_bc_ref - Vp) * dh1
gw = 2.0 * (w_bc_ref - Wp) * dh1
wv(i,j,k,1) = wv(i,j,k,1) - ( 0.5 * u_bc * gu - gu * dh2 )
wv(i,j,k,2) = wv(i,j,k,2) - ( 0.5 * u_bc * gv - gv * dh2 )
wv(i,j,k,3) = wv(i,j,k,3) - ( 0.5 * u_bc * gw - gw * dh2 )
m2 = m2 + 1.0
end if


! Y方向 ---------------------------------------
if ( c_s == 1.0 ) then
gu = 2.0 * (Up - u_bc_ref) * dh1
gv = 2.0 * (Vp - v_bc_ref) * dh1
gw = 2.0 * (Wp - w_bc_ref) * dh1
wv(i,j,k,1) = wv(i,j,k,1) - ( 0.5 * v_bc * gu + gu * dh2 )
wv(i,j,k,2) = wv(i,j,k,2) - ( 0.5 * v_bc * gv + gv * dh2 )
wv(i,j,k,3) = wv(i,j,k,3) - ( 0.5 * v_bc * gw + gw * dh2 )
m2 = m2 + 1.0
end if

if ( c_n == 1.0 ) then
gu = 2.0 * (u_bc_ref - Up) * dh1
gv = 2.0 * (v_bc_ref - Vp) * dh1
gw = 2.0 * (w_bc_ref - Wp) * dh1
wv(i,j,k,1) = wv(i,j,k,1) - ( 0.5 * v_bc * gu - gu * dh2 )
wv(i,j,k,2) = wv(i,j,k,2) - ( 0.5 * v_bc * gv - gv * dh2 )
wv(i,j,k,3) = wv(i,j,k,3) - ( 0.5 * v_bc * gw - gw * dh2 )
m2 = m2 + 1.0
end if


! Z方向 ---------------------------------------
if ( c_b == 1.0 ) then
gu = 2.0 * (Up - u_bc_ref) * dh1
gv = 2.0 * (Vp - v_bc_ref) * dh1
gw = 2.0 * (Wp - w_bc_ref) * dh1
wv(i,j,k,1) = wv(i,j,k,1) - ( 0.5 * w_bc * gu + gu * dh2 )
wv(i,j,k,2) = wv(i,j,k,2) - ( 0.5 * w_bc * gv + gv * dh2 )
wv(i,j,k,3) = wv(i,j,k,3) - ( 0.5 * w_bc * gw + gw * dh2 )
m2 = m2 + 1.0
end if

if ( c_t == 1.0 ) then
gu = 2.0 * (u_bc_ref - Up) * dh1
gv = 2.0 * (v_bc_ref - Vp) * dh1
gw = 2.0 * (w_bc_ref - Wp) * dh1
wv(i,j,k,1) = wv(i,j,k,1) - ( 0.5 * w_bc * gu - gu * dh2 )
wv(i,j,k,2) = wv(i,j,k,2) - ( 0.5 * w_bc * gv - gv * dh2 )
wv(i,j,k,3) = wv(i,j,k,3) - ( 0.5 * w_bc * gw - gw * dh2 )
m2 = m2 + 1.0
end if


endif
end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL

flop = flop + m2*24.0

return
end subroutine pvec_ibc_specv_fdm


!> ********************************************************************
!! @brief 内部回転境界条件による対流項と粘性項の流束の修正（有限体積型のスキーム）
!! @param [out] wv   疑似ベクトルの空間項の評価値
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル長
!! @param [in]  st   ループの開始インデクス
!! @param [in]  ed   ループの終了インデクス
!! @param [in]  dh   格子幅
!! @param [in]  rei  Reynolds数の逆数
!! @param [in]  v    セルセンター速度ベクトル（u^n）
!! @param [in]  bv   BCindex C
!! @param [in]  odr  内部境界処理時の速度境界条件のエントリ
!! @param [in]  vec  角速度ベクトル
!! @param [in]  org  領域基点座標
!! @param [in]  ctr  回転中心座標
!! @param [out] flop 浮動小数点演算数
!! @note カット位置に関わらず指定速度で流束を計算
!<
subroutine pvec_ibc_sldrev_fvm (wv, sz, g, st, ed, dh, rei, v, bv, odr, vec, ctr, org, flop)
implicit none
include 'ffv_f_params.h'
integer                                                     ::  i, j, k, g, bvx, odr, is, ie, js, je, ks, ke
integer, dimension(3)                                       ::  sz, st, ed
double precision                                            ::  flop
real                                                        ::  Up0, Ue1, Uw1, Us1, Un1, Ub1, Ut1
real                                                        ::  Vp0, Ve1, Vw1, Vs1, Vn1, Vb1, Vt1
real                                                        ::  Wp0, We1, Ww1, Ws1, Wn1, Wb1, Wt1
real                                                        ::  c_e, c_w, c_n, c_s, c_t, c_b
real                                                        ::  dh, dh1, dh2, EX, EY, EZ, rei
real                                                        ::  m1, m2, cx, cy, cz
real                                                        ::  cnv_u, cnv_v, cnv_w, cr, cl, acr, acl
real                                                        ::  u_bc, v_bc, w_bc, dx, dy, dz
real                                                        ::  fu_r, fu_l, fv_r, fv_l, fw_r, fw_l
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)   ::  v, wv
integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bv
real, dimension(3)                                          ::  vec, org, ctr
real                                                        ::  r_u, r_v, r_w, omg_u, omg_v, omg_w, ox, oy, oz

is = st(1)
ie = ed(1)
js = st(2)
je = ed(2)
ks = st(3)
ke = ed(3)

dh1= 1.0/dh
dh2= rei*dh1*dh1

ox = org(1)
oy = org(2)
oz = org(3)

cx = ctr(1)
cy = ctr(2)
cz = ctr(3)

omg_u = vec(1)
omg_v = vec(2)
omg_w = vec(3)

dx = dh
dy = dh
dz = dh

m1 = 0.0
m2 = 0.0

flop = flop + 13.0d0 ! DP 18 flop

!$OMP PARALLEL &
!$OMP REDUCTION(+:m1) &
!$OMP REDUCTION(+:m2) &
!$OMP FIRSTPRIVATE(is, ie, js, je, ks, ke, ox, oy, oz) &
!$OMP FIRSTPRIVATE(dh1, dh2, odr, dx, dy, dz, cx, cy, cz, omg_u, omg_v, omg_w) &
!$OMP PRIVATE(bvx, cnv_u, cnv_v, cnv_w, EX, EY, EZ, cr, cl, acr, acl) &
!$OMP PRIVATE(Up0, Ue1, Uw1, Us1, Un1, Ub1, Ut1) &
!$OMP PRIVATE(Vp0, Ve1, Vw1, Vs1, Vn1, Vb1, Vt1) &
!$OMP PRIVATE(Wp0, We1, Ww1, Ws1, Wn1, Wb1, Wt1) &
!$OMP PRIVATE(c_e, c_w, c_n, c_s, c_t, c_b) &
!$OMP PRIVATE(fu_r, fu_l, fv_r, fv_l, fw_r, fw_l) &
!$OMP PRIVATE(r_u, r_v, r_w, u_bc, v_bc, w_bc)

!$OMP DO SCHEDULE(static)

do k=ks,ke
do j=js,je
do i=is,ie
bvx = bv(i,j,k)

if ( 0 /= iand(bvx, bc_mask30) ) then ! 6面のうちのどれか速度境界フラグが立っている場合
cnv_u = 0.0
cnv_v = 0.0
cnv_w = 0.0

! 変数のロード
Up0 = v(i,j,k,1)
Vp0 = v(i,j,k,2)
Wp0 = v(i,j,k,3)

Uw1 = 0.0
Ue1 = 0.0
Us1 = 0.0
Un1 = 0.0
Ub1 = 0.0
Ut1 = 0.0
Vw1 = 0.0
Ve1 = 0.0
Vs1 = 0.0
Vn1 = 0.0
Vb1 = 0.0
Vt1 = 0.0
Ww1 = 0.0
We1 = 0.0
Ws1 = 0.0
Wn1 = 0.0
Wb1 = 0.0
Wt1 = 0.0

! 内部境界のときの各面のBCフラグ ibits() = 0(Normal) / others(BC) >> c_e = 0.0(Normal) / 1.0(BC)
c_e = 0.0
c_w = 0.0
c_n = 0.0
c_s = 0.0
c_t = 0.0
c_b = 0.0
if ( ibits(bvx, bc_face_E, bitw_5) == odr ) c_e = 1.0
if ( ibits(bvx, bc_face_W, bitw_5) == odr ) c_w = 1.0
if ( ibits(bvx, bc_face_N, bitw_5) == odr ) c_n = 1.0
if ( ibits(bvx, bc_face_S, bitw_5) == odr ) c_s = 1.0
if ( ibits(bvx, bc_face_T, bitw_5) == odr ) c_t = 1.0
if ( ibits(bvx, bc_face_B, bitw_5) == odr ) c_b = 1.0


! X方向 ---------------------------------------
if ( c_w == 1.0 ) then
r_u = ox + (real(i)-1.0)*dx - cx
r_v = oy + (real(j)-0.5)*dy - cy
r_w = oz + (real(k)-0.5)*dz - cz
u_bc = omg_v * r_w - omg_w * r_v
v_bc = omg_w * r_u - omg_u * r_w
w_bc = omg_u * r_v - omg_v * r_u
Uw1  = u_bc
Vw1  = v_bc
Ww1  = w_bc
cl   = u_bc
acl  = abs(cl)
fu_l = 0.5*(cl*(Up0+Uw1) - acl*(Up0-Uw1))
fv_l = 0.5*(cl*(Vp0+Vw1) - acl*(Vp0-Vw1))
fw_l = 0.5*(cl*(Wp0+Ww1) - acl*(Wp0-Ww1)) ! 40+1 flops
m2 = m2 + 1.0
end if

if ( c_e == 1.0 ) then
r_u = ox +  real(i)     *dx - cx
r_v = oy + (real(j)-0.5)*dy - cy
r_w = oz + (real(k)-0.5)*dz - cz
u_bc = omg_v * r_w - omg_w * r_v
v_bc = omg_w * r_u - omg_u * r_w
w_bc = omg_u * r_v - omg_v * r_u
Ue1  = u_bc
Ve1  = v_bc
We1  = w_bc
cr   = u_bc
acr  = abs(cr)
fu_r = 0.5*(cr*(Ue1+Up0) - acr*(Ue1-Up0))
fv_r = 0.5*(cr*(Ve1+Vp0) - acr*(Ve1-Vp0))
fw_r = 0.5*(cr*(We1+Wp0) - acr*(We1-Wp0))
m2 = m2 + 1.0
end if

cnv_u = cnv_u + fu_r*c_e - fu_l*c_w
cnv_v = cnv_v + fv_r*c_e - fv_l*c_w
cnv_w = cnv_w + fw_r*c_e - fw_l*c_w ! 12 flops

! Y方向 ---------------------------------------
if ( c_s == 1.0 ) then
r_u = ox + (real(i)-0.5)*dx - cx
r_v = oy + (real(j)-1.0)*dy - cy
r_w = oz + (real(k)-0.5)*dz - cz
u_bc = omg_v * r_w - omg_w * r_v
v_bc = omg_w * r_u - omg_u * r_w
w_bc = omg_u * r_v - omg_v * r_u
Us1  = u_bc
Vs1  = v_bc
Ws1  = w_bc
cl   = v_bc
acl  = abs(cl)
fu_l = 0.5*(cl*(Up0+Us1) - acl*(Up0-Us1))
fv_l = 0.5*(cl*(Vp0+Vs1) - acl*(Vp0-Vs1))
fw_l = 0.5*(cl*(Wp0+Ws1) - acl*(Wp0-Ws1))
m2 = m2 + 1.0
end if

if ( c_n == 1.0 ) then
r_u = ox + (real(i)-0.5)*dx - cx
r_v = oy +  real(j)     *dy - cy
r_w = oz + (real(k)-0.5)*dz - cz
u_bc = omg_v * r_w - omg_w * r_v
v_bc = omg_w * r_u - omg_u * r_w
w_bc = omg_u * r_v - omg_v * r_u
Un1  = u_bc
Vn1  = v_bc
Wn1  = w_bc
cr   = v_bc
acr  = abs(cr)
fu_r = 0.5*(cr*(Un1+Up0) - acr*(Un1-Up0))
fv_r = 0.5*(cr*(Vn1+Vp0) - acr*(Vn1-Vp0))
fw_r = 0.5*(cr*(Wn1+Wp0) - acr*(Wn1-Wp0))
m2 = m2 + 1.0
end if

cnv_u = cnv_u + fu_r*c_n - fu_l*c_s
cnv_v = cnv_v + fv_r*c_n - fv_l*c_s
cnv_w = cnv_w + fw_r*c_n - fw_l*c_s

! Z方向 ---------------------------------------
if ( c_b == 1.0 ) then
r_u = ox + (real(i)-0.5)*dx - cx
r_v = oy + (real(j)-0.5)*dy - cy
r_w = oz + (real(k)-1.0)*dz - cz
u_bc = omg_v * r_w - omg_w * r_v
v_bc = omg_w * r_u - omg_u * r_w
w_bc = omg_u * r_v - omg_v * r_u
Ub1  = u_bc
Vb1  = v_bc
Wb1  = w_bc
cl   = w_bc
acl  = abs(cl)
fu_l = 0.5*(cl*(Up0+Ub1) - acl*(Up0-Ub1))
fv_l = 0.5*(cl*(Vp0+Vb1) - acl*(Vp0-Vb1))
fw_l = 0.5*(cl*(Wp0+Wb1) - acl*(Wp0-Wb1))
m2 = m2 + 1.0
end if

if ( c_t == 1.0 ) then
r_u = ox + (real(i)-0.5)*dx - cx
r_v = oy + (real(j)-0.5)*dy - cy
r_w = oz +  real(k)     *dz - cz
u_bc = omg_v * r_w - omg_w * r_v
v_bc = omg_w * r_u - omg_u * r_w
w_bc = omg_u * r_v - omg_v * r_u
Ut1  = u_bc
Vt1  = v_bc
Wt1  = w_bc
cr   = w_bc
acr  = abs(cr)
fu_r = 0.5*(cr*(Ut1+Up0) - acr*(Ut1-Up0))
fv_r = 0.5*(cr*(Vt1+Vp0) - acr*(Vt1-Vp0))
fw_r = 0.5*(cr*(Wt1+Wp0) - acr*(Wt1-Wp0))
m2 = m2 + 1.0
end if

cnv_u = cnv_u + fu_r*c_t - fu_l*c_b
cnv_v = cnv_v + fv_r*c_t - fv_l*c_b
cnv_w = cnv_w + fw_r*c_t - fw_l*c_b

! 粘性項
EX = ( Ue1 - Up0 ) * c_e &
   + ( Uw1 - Up0 ) * c_w &
   + ( Un1 - Up0 ) * c_n &
   + ( Us1 - Up0 ) * c_s &
   + ( Ut1 - Up0 ) * c_t &
   + ( Ub1 - Up0 ) * c_b ! 17 flops
EY = ( Ve1 - Vp0 ) * c_e &
   + ( Vw1 - Vp0 ) * c_w &
   + ( Vn1 - Vp0 ) * c_n &
   + ( Vs1 - Vp0 ) * c_s &
   + ( Vt1 - Vp0 ) * c_t &
   + ( Vb1 - Vp0 ) * c_b
EZ = ( We1 - Wp0 ) * c_e &
   + ( Ww1 - Wp0 ) * c_w &
   + ( Wn1 - Wp0 ) * c_n &
   + ( Ws1 - Wp0 ) * c_s &
   + ( Wt1 - Wp0 ) * c_t &
   + ( Wb1 - Wp0 ) * c_b

wv(i,j,k,1) = wv(i,j,k,1) + ( -cnv_u*dh1 + EX*dh2 )
wv(i,j,k,2) = wv(i,j,k,2) + ( -cnv_v*dh1 + EY*dh2 )
wv(i,j,k,3) = wv(i,j,k,3) + ( -cnv_w*dh1 + EZ*dh2 ) ! 4*3 = 12 flops
m1 = m1 + 1.0

endif
end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL

! loop :  (12*3 + 17*3 + 12)*m1 + 41*m2 = 99*m1 + 41*m2

flop = flop + m1*99.0 + m2*41.0

return
end subroutine pvec_ibc_sldrev_fvm



!> ********************************************************************
!! @brief 内部速度指定境界条件による疑似速度の\sum{u}の修正
!! @param [in,out] div  \sum{u_j}
!! @param [in]     sz   配列長
!! @param [in]     g    ガイドセル長
!! @param [in]     st   ループの開始インデクス
!! @param [in]     ed   ループの終了インデクス
!! @param [in]     v00  参照速度
!! @param [in]     bv   BCindex C
!! @param [in]     odr  速度境界条件のエントリ
!! @param [in]     vec  指定する速度ベクトル
!! @param [in,out] flop flop count 近似
!<
    subroutine div_ibc_drchlt (div, sz, g, st, ed, v00, bv, odr, vec, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, g, bvx, odr, is, ie, js, je, ks, ke
    integer, dimension(3)                                     ::  sz, st, ed
    double precision                                          ::  flop, m
    real                                                      ::  Ue_t, Uw_t, Vn_t, Vs_t, Wt_t, Wb_t
    real                                                      ::  u_bc_ref, v_bc_ref, w_bc_ref
    real, dimension(3)                                        ::  vec
    real, dimension(0:3)                                      ::  v00
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  div
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bv

    ! u_bc_refは参照座標系での境界速度
    u_bc_ref = vec(1) + v00(1)
    v_bc_ref = vec(2) + v00(2)
    w_bc_ref = vec(3) + v00(3)

    is = st(1)
    ie = ed(1)
    js = st(2)
    je = ed(2)
    ks = st(3)
    ke = ed(3)

    m = dble( (ie-is+1)*(je-js+1)*(ke-ks+1) )
    flop = flop + m*6.0d0

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(is, ie, js, je, ks, ke, u_bc_ref, v_bc_ref, w_bc_ref, odr) &
!$OMP PRIVATE(bvx, Ue_t, Uw_t, Vn_t, Vs_t, Wt_t, Wb_t)

!$OMP DO SCHEDULE(static)

    do k=ks,ke
    do j=js,je
    do i=is,ie
      bvx = bv(i,j,k)
      if ( 0 /= iand(bvx, bc_mask30) ) then ! 6面のうちのどれか速度境界フラグが立っている場合
        Ue_t = 0.0
        Uw_t = 0.0
        Vn_t = 0.0
        Vs_t = 0.0
        Wt_t = 0.0
        Wb_t = 0.0
        
        if ( ibits(bvx, bc_face_W, bitw_5) == odr ) Uw_t = u_bc_ref
        if ( ibits(bvx, bc_face_E, bitw_5) == odr ) Ue_t = u_bc_ref
        if ( ibits(bvx, bc_face_S, bitw_5) == odr ) Vs_t = v_bc_ref
        if ( ibits(bvx, bc_face_N, bitw_5) == odr ) Vn_t = v_bc_ref
        if ( ibits(bvx, bc_face_B, bitw_5) == odr ) Wb_t = w_bc_ref
        if ( ibits(bvx, bc_face_T, bitw_5) == odr ) Wt_t = w_bc_ref

        ! VBCの面だけUe_tなどは値をもつ  対象セルは流体なのでマスク不要
        div(i,j,k) = div(i,j,k) + ( Ue_t - Uw_t + Vn_t - Vs_t + Wt_t - Wb_t )
      end if
    end do
    end do
    end do
    
!$OMP END DO
!$OMP END PARALLEL

    return
    end subroutine div_ibc_drchlt



!> ********************************************************************
!! @brief 回転速度境界条件による疑似速度の\sum{u}の修正
!! @param [in,out] div  \sum{u_j}
!! @param [in]     sz   配列長
!! @param [in]     g    ガイドセル長
!! @param [in]     st   ループの開始インデクス
!! @param [in]     ed   ループの終了インデクス
!! @param [in]     pch  格子幅
!! @param [in]     v00  参照速度
!! @param [in]     bv   BCindex C
!! @param [in]     odr  速度境界条件のエントリ
!! @param [in]     vec  角速度ベクトル
!! @param [in]     org  領域基点座標
!! @param [in]     ctr  回転中心座標
!! @param [in,out] flop flop count 近似
!<
subroutine div_ibc_sldrev (div, sz, g, st, ed, pch, v00, bv, odr, vec, ctr, org, flop)
implicit none
include 'ffv_f_params.h'
integer                                                   ::  i, j, k, g, bvx, odr, is, ie, js, je, ks, ke
integer, dimension(3)                                     ::  sz, st, ed
double precision                                          ::  flop, m
real                                                      ::  Ue_t, Uw_t, Vn_t, Vs_t, Wt_t, Wb_t
real                                                      ::  u_bc, v_bc, w_bc, cx, cy, cz, dx, dy, dz
real, dimension(3)                                        ::  vec, ctr, org, pch
real, dimension(0:3)                                      ::  v00
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  div
integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bv
real                                                      ::  r_u, r_v, r_w, omg_u, omg_v, omg_w, ox, oy, oz

ox = org(1)
oy = org(2)
oz = org(3)

cx = ctr(1)
cy = ctr(2)
cz = ctr(3)

omg_u = vec(1)
omg_v = vec(2)
omg_w = vec(3)

dx = pch(1)
dy = pch(2)
dz = pch(3)

is = st(1)
ie = ed(1)
js = st(2)
je = ed(2)
ks = st(3)
ke = ed(3)

m = dble( (ie-is+1)*(je-js+1)*(ke-ks+1) )
flop = flop + m*6.0d0

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(is, ie, js, je, ks, ke) &
!$OMP FIRSTPRIVATE(dx, dy, dz, cx, cy, cz, omg_u, omg_v, omg_w, ox, oy, oz) &
!$OMP PRIVATE(bvx, Ue_t, Uw_t, Vn_t, Vs_t, Wt_t, Wb_t) &
!$OMP PRIVATE(r_u, r_v, r_w, u_bc, v_bc, w_bc)

!$OMP DO SCHEDULE(static)

do k=ks,ke
do j=js,je
do i=is,ie
bvx = bv(i,j,k)
if ( 0 /= iand(bvx, bc_mask30) ) then ! 6面のうちのどれか速度境界フラグが立っている場合
Ue_t = 0.0
Uw_t = 0.0
Vn_t = 0.0
Vs_t = 0.0
Wt_t = 0.0
Wb_t = 0.0

if ( ibits(bvx, bc_face_W, bitw_5) == odr ) then
!r_u = ox + (real(i)-1.0)*dx - cx
r_v = oy + (real(j)-0.5)*dy - cy
r_w = oz + (real(k)-0.5)*dz - cz
u_bc = omg_v * r_w - omg_w * r_v
!v_bc = omg_w * r_u - omg_u * r_w
!w_bc = omg_u * r_v - omg_v * r_u
Uw_t = u_bc
endif

if ( ibits(bvx, bc_face_E, bitw_5) == odr ) then
!r_u = ox +  real(i)     *dx - cx
r_v = oy + (real(j)-0.5)*dy - cy
r_w = oz + (real(k)-0.5)*dz - cz
u_bc = omg_v * r_w - omg_w * r_v
!v_bc = omg_w * r_u - omg_u * r_w
!w_bc = omg_u * r_v - omg_v * r_u
Ue_t = u_bc
endif

if ( ibits(bvx, bc_face_S, bitw_5) == odr ) then
r_u = ox + (real(i)-0.5)*dx - cx
!r_v = oy + (real(j)-1.0)*dy - cy
r_w = oz + (real(k)-0.5)*dz - cz
!u_bc = omg_v * r_w - omg_w * r_v
v_bc = omg_w * r_u - omg_u * r_w
!w_bc = omg_u * r_v - omg_v * r_u
Vs_t = v_bc
endif

if ( ibits(bvx, bc_face_N, bitw_5) == odr ) then
r_u = ox + (real(i)-0.5)*dx - cx
!r_v = oy +  real(j)     *dy - cy
r_w = oz + (real(k)-0.5)*dz - cz
!u_bc = omg_v * r_w - omg_w * r_v
v_bc = omg_w * r_u - omg_u * r_w
!w_bc = omg_u * r_v - omg_v * r_u
Vn_t = v_bc
endif

if ( ibits(bvx, bc_face_B, bitw_5) == odr ) then
r_u = ox + (real(i)-0.5)*dx - cx
r_v = oy + (real(j)-0.5)*dy - cy
!r_w = oz + (real(k)-1.0)*dz - cz
!u_bc = omg_v * r_w - omg_w * r_v
!v_bc = omg_w * r_u - omg_u * r_w
w_bc = omg_u * r_v - omg_v * r_u
Wb_t = w_bc
endif

if ( ibits(bvx, bc_face_T, bitw_5) == odr ) then
r_u = ox + (real(i)-0.5)*dx - cx
r_v = oy + (real(j)-0.5)*dy - cy
!r_w = oz +  real(k)     *dz - cz
!u_bc = omg_v * r_w - omg_w * r_v
!v_bc = omg_w * r_u - omg_u * r_w
w_bc = omg_u * r_v - omg_v * r_u
Wt_t = w_bc
endif

! VBCの面だけUe_tなどは値をもつ  対象セルは流体なのでマスク不要
div(i,j,k) = div(i,j,k) + ( Ue_t - Uw_t + Vn_t - Vs_t + Wt_t - Wb_t )
end if
end do
end do
end do

!$OMP END DO
!$OMP END PARALLEL

return
end subroutine div_ibc_sldrev



!> ********************************************************************
!! @brief 内部流出境界条件による疑似速度ベクトルの発散の修正
!! @param [in,out] div  速度の発散
!! @param [in]     sz   配列長
!! @param [in]     g    ガイドセル長
!! @param [in]     st   ループの開始インデクス
!! @param [in]     ed   ループの終了インデクス
!! @param [in]     v00  参照速度
!! @param [in]     cf   u_out*dt/dh
!! @param [in]     bv   BCindex C
!! @param [in]     odr  速度境界条件のエントリ
!! @param [in]     v0   セルセンター速度　u^n
!! @param [in]     vf   セルフェイス速度ベクトル（n-step）
!! @param [in,out] flop flop count
!! @note 流出境界面ではu_e^{n+1}=u_e^n-cf*(u_e^n-u_w^n)を予測値としてdivの寄与として加算
!! @note flop countはコスト軽減のため近似
!<
    subroutine div_ibc_oflow_pvec (div, sz, g, st, ed, v00, cf, bv, odr, v0, vf, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, g, bvx, odr, is, ie, js, je, ks, ke
    integer, dimension(3)                                     ::  sz, st, ed
    double precision                                          ::  flop
    real                                                      ::  m
    real                                                      ::  b_w, b_e, b_s, b_n, b_b, b_t, b_p
    real                                                      ::  w_e, w_w, w_n, w_s, w_t, w_b
    real                                                      ::  Ue, Uw, Vn, Vs, Wt, Wb
    real                                                      ::  Ue_t, Uw_t, Vn_t, Vs_t, Wt_t, Wb_t
    real                                                      ::  cf, u_ref, v_ref, w_ref
    real, dimension(0:3)                                      ::  v00
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  div
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v0, vf
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bv

    ! 参照速度
    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)

    is = st(1)
    ie = ed(1)
    js = st(2)
    je = ed(2)
    ks = st(3)
    ke = ed(3)

    m = 0.0
    
!$OMP PARALLEL REDUCTION(+:m) &
!$OMP FIRSTPRIVATE(is, ie, js, je, ks, ke, u_ref, v_ref, w_ref, odr, cf) &
!$OMP PRIVATE(bvx) &
!$OMP PRIVATE(b_w, b_e, b_s, b_n, b_b, b_t, b_p) &
!$OMP PRIVATE(w_e, w_w, w_n, w_s, w_t, w_b) &
!$OMP PRIVATE(Ue, Uw, Vn, Vs, Wt, Wb) &
!$OMP PRIVATE(Ue_t, Uw_t, Vn_t, Vs_t, Wt_t, Wb_t)

!$OMP DO SCHEDULE(static)

    do k=ks,ke
    do j=js,je
    do i=is,ie
      bvx = bv(i,j,k)
      if ( 0 /= iand(bvx, bc_mask30) ) then ! 6面のうちのどれか速度境界フラグが立っている場合

        ! セルの状態 (0-solid / 1-fluid)
        b_p = real( ibits(bv(i  ,j  ,k  ), State, 1) )
        b_w = real( ibits(bv(i-1,j  ,k  ), State, 1) )
        b_e = real( ibits(bv(i+1,j  ,k  ), State, 1) )
        b_s = real( ibits(bv(i  ,j-1,k  ), State, 1) )
        b_n = real( ibits(bv(i  ,j+1,k  ), State, 1) )
        b_b = real( ibits(bv(i  ,j  ,k-1), State, 1) )
        b_t = real( ibits(bv(i  ,j  ,k+1), State, 1) )

        ! セル界面のフラグ (0-wall face / 1-fluid) > real*6+ 6= 12 flops
        w_e = real(b_e * b_p)
        w_w = real(b_w * b_p)
        w_n = real(b_n * b_p)
        w_s = real(b_s * b_p)
        w_t = real(b_t * b_p)
        w_b = real(b_b * b_p)

        ! 界面速度（スタガード位置） > 6 flops
        Uw = vf(i-1, j  , k  ,1)*w_w
        Ue = vf(i  , j  , k  ,1)*w_e
        Vs = vf(i  , j-1, k  ,2)*w_s
        Vn = vf(i  , j  , k  ,2)*w_n
        Wb = vf(i  , j  , k-1,3)*w_b
        Wt = vf(i  , j  , k  ,3)*w_t ! real*7 real*6 + 6 + 6 = 25 flopss

        ! 寄与をゼロにしておく
        Ue_t = 0.0
        Uw_t = 0.0
        Vn_t = 0.0
        Vs_t = 0.0
        Wt_t = 0.0
        Wb_t = 0.0
        
        ! X方向 ---------------------------------------
        if ( ibits(bvx, bc_face_W, bitw_5) == odr ) then
          !Uw = Ue + (Vn - Vs + Wt - Wb) ! 連続の式から流出面の速度を推定，これは移動座標系上の速度成分
          if ( cf>0.0 ) cf=0.0
          Uw_t = Uw - cf*(Ue-Uw)
        endif
        
        if ( ibits(bvx, bc_face_E, bitw_5) == odr ) then
          !Ue = Uw - (Vn - Vs + Wt - Wb)
          if ( cf<0.0 ) cf=0.0
          Ue_t = Ue - cf*(Ue-Uw)
        endif
        
        ! Y方向 ---------------------------------------
        if ( ibits(bvx, bc_face_S, bitw_5) == odr ) then
          !Vs = Vn + (Ue - Uw + Wt - Wb)
          if ( cf>0.0 ) cf=0.0
          Vs_t = Vs - cf*(Vn-Vs)
        endif
        
        if ( ibits(bvx, bc_face_N, bitw_5) == odr ) then
          !Vn = Vs - (Ue - Uw + Wt - Wb)
          if ( cf<0.0 ) cf=0.0
          Vn_t = Vn - cf*(Vn-Vs)
        endif
        
        ! Z方向 ---------------------------------------
        if ( ibits(bvx, bc_face_B, bitw_5) == odr ) then
          !Wb = Wt + (Ue - Uw + Vn - Vs)
          if ( cf>0.0 ) cf=0.0
          Wb_t = Wb - cf*(Wt-Wb)
        endif
        
        if ( ibits(bvx, bc_face_T, bitw_5) == odr ) then
          !Wt = Wb - (Ue - Uw + Vn - Vs)
          if ( cf<0.0 ) cf=0.0
          Wt_t = Wt - cf*(Wt-Wb)
        endif

        ! VBCの面だけUe_tなどは値をもつ
        div(i,j,k) = div(i,j,k) + ( Ue_t - Uw_t + Vn_t - Vs_t + Wt_t - Wb_t ) ! 対象セルは流体なのでマスク不要
        m = m + 1.0
      end if
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    flop = flop + dble(m)*66.0d0

    return
    end subroutine div_ibc_oflow_pvec
    
!> ********************************************************************
!! @brief 内部流出境界条件によるn+1時刻の速度の発散の修正と流量の積算
!! @param [in,out] div  速度の発散
!! @param [in]     sz   配列長
!! @param [in]     g    ガイドセル長
!! @param [in]     st   ループの開始インデクス
!! @param [in]     ed   ループの終了インデクス
!! @param [in]     bv   BCindex C
!! @param [in]     odr  速度境界条件のエントリ
!! @param [out]    av   積算速度と積算数
!! @param [out]    flop flop count
!! @note div(u)=0から，内部流出境界のセルで計算されたdivの値が流出速度となる
!! @note 1つのセルに複数の速度境界条件がある場合にはだめ
!! @note flop countはコスト軽減のため近似
!<
    subroutine div_ibc_oflow_vec (div, sz, g, st, ed, bv, odr, av, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, g, bvx, odr, is, ie, js, je, ks, ke
    integer, dimension(3)                                     ::  sz, st, ed
    double precision                                          ::  flop
    real                                                      ::  dv, a1, m
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  div
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bv
    real, dimension(2)                                        ::  av
    
    a1 = 0.0
    m = 0.0

    is = st(1)
    ie = ed(1)
    js = st(2)
    je = ed(2)
    ks = st(3)
    ke = ed(3)


!$OMP PARALLEL &
!$OMP REDUCTION(+:a1) &
!$OMP REDUCTION(+:m) &
!$OMP FIRSTPRIVATE(is, ie, js, je, ks, ke, odr) &
!$OMP PRIVATE(bvx, dv)

!$OMP DO SCHEDULE(static)
    do k=ks,ke
    do j=js,je
    do i=is,ie
      bvx = bv(i,j,k)
      if ( 0 /= iand(bvx, bc_mask30) ) then ! 6面のうちのどれか速度境界フラグが立っている場合
        dv = div(i,j,k)
        
        if ( ibits(bvx, bc_face_W, bitw_5) == odr ) then ! u_w
          a1 = a1 + dv
        endif

        if ( ibits(bvx, bc_face_E, bitw_5) == odr ) then ! u_e
          a1 = a1 - dv
        endif
        
        if ( ibits(bvx, bc_face_S, bitw_5) == odr ) then
          a1 = a1 + dv
        endif
        
        if ( ibits(bvx, bc_face_N, bitw_5) == odr ) then
          a1 = a1 - dv
        endif
        
        if ( ibits(bvx, bc_face_B, bitw_5) == odr ) then
          a1 = a1 + dv
        endif
        
        if ( ibits(bvx, bc_face_T, bitw_5) == odr ) then
          a1 = a1 - dv
        endif

        div(i,j,k) = 0.0 ! 対象セルは発散をゼロにする
        m = m + 1.0
      end if
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    av(1) = a1
    av(2) = m
    flop = flop + dble(m)*2.0d0

    return
    end subroutine div_ibc_oflow_vec
