!********************************************************************
!
!   FFV : Frontflow / violet Cartesian
!
!   Copyright (c) 2012-2013  All right reserved.
!
!   Institute of Industrial Science, University of Tokyo, Japan.
!
!********************************************************************

!> @file   IP_boundary.f90
!! @brief  組み込み例題の境界条件
!! @author kero
!<


!> ********************************************************************
!! @brief JETの流入条件：対流項と粘性項の流束の修正
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
subroutine vobc_pv_jet (wv, sz, g, dh, rei, v0, bv, vec, face, &
  pat1, r1i, r1o, omg1, q1, a1, pat2, r2i, r2o, omg2, q2, a2, &
  flop)
implicit none
include '../F_CORE/ffv_f_params.h'
integer                                                   ::  i, j, k, g, face
integer                                                   ::  ix, jx, kx, pat1, pat2
integer, dimension(3)                                     ::  sz
double precision                                          ::  flop
real                                                      ::  Up, Vp, Wp, Ur, Vr, Wr
real                                                      ::  dh, dh1, dh2, EX, EY, EZ, rei
real                                                      ::  fu, fv, fw, c, ac, msk
real                                                      ::  u_bc, v_bc, w_bc, m
real                                                      ::  r1i, r1o, omg1, q1, a1
real                                                      ::  r2i, r2o, omg2, q2, a2
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
i = 1



! X_MINUS


!$OMP PARALLEL REDUCTION(+:m) &
!$OMP FIRSTPRIVATE(i, ix, jx, kx, u_bc, v_bc, w_bc, face) &
!$OMP FIRSTPRIVATE(dh1, dh2) &
!$OMP PRIVATE(Up, Vp, Wp, Ur, Vr, Wr) &
!$OMP PRIVATE(fu, fv, fw, EX, EY, EZ, c, ac, msk) &
!$OMP PRIVATE(j, k)
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

flop = flop + dble(m)*28.0d0

return
end subroutine vobc_pv_jet