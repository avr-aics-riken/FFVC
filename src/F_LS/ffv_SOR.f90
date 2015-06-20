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

!> @file   ffv_SOR.f90
!! @brief  SOR series routine
!! @author aics
!<


!> ********************************************************************
!! @brief point SOR法
!! @param [in,out] p    圧力
!! @param [in]     sz   配列長
!! @param [in]     g    ガイドセル長
!! @param [in]     dh   格子幅
!! @param [in]     omg  加速係数
!! @param [out]    cnv  収束判定値　修正量の自乗和と残差の自乗和、解ベクトルの自乗和
!! @param [in]     b    RHS vector
!! @param [in]     bp   BCindex P
!! @param [out]    flop flop count
!! @note Activeマスクの位置は，固体中のラプラス式を解くように，更新式にはかけず残差のみにする
!<
  subroutine psor (p, sz, g, dh, omg, cnv, b, bp, flop)
  implicit none
  include 'ffv_f_params.h'
  integer                                                   ::  i, j, k, ix, jx, kx, g, idx
  integer, dimension(3)                                     ::  sz
  double precision                                          ::  flop, res, err, aa, xl2
  real                                                      ::  omg, dd, ss, dp, pp, bb, de, pn, dsw
  real                                                      ::  c_w, c_e, c_s, c_n, c_b, c_t
  real                                                      ::  d_w, d_e, d_s, d_n, d_b, d_t
  real                                                      ::  r_xx, r_xy, r_xz, r_x2, r_y2, r_z2
  real, dimension(3)                                        ::  dh
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  p, b
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp
  double precision, dimension(3)                            ::  cnv

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
  res = 0.0
  err = 0.0
  xl2 = 0.0

  r_xx = 1.0
  r_xy = dh(1) / dh(2)
  r_xz = dh(1) / dh(3)
  r_x2 = r_xx * r_xx
  r_y2 = r_xy * r_xy
  r_z2 = r_xz * r_xz

  flop = flop + dble(ix)*dble(jx)*dble(kx)*56.0d0 + 19.0d0


!$OMP PARALLEL &
!$OMP REDUCTION(+:res) &
!$OMP REDUCTION(+:err) &
!$OMP REDUCTION(+:xl2) &
!$OMP PRIVATE(c_w, c_e, c_s, c_n, c_b, c_t) &
!$OMP PRIVATE(d_w, d_e, d_s, d_n, d_b, d_t) &
!$OMP PRIVATE(dd, ss, dp, idx, aa, pp, bb, de, pn, dsw) &
!$OMP FIRSTPRIVATE(ix, jx, kx, omg) &
!$OMP FIRSTPRIVATE(r_x2, r_y2, r_z2)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
  do k=1,kx
  do j=1,jx
  do i=1,ix
    idx = bp(i,j,k)

    c_w = real(ibits(idx, bc_ndag_W, 1))  ! w
    c_e = real(ibits(idx, bc_ndag_E, 1))  ! e
    c_s = real(ibits(idx, bc_ndag_S, 1))  ! s
    c_n = real(ibits(idx, bc_ndag_N, 1))  ! n
    c_b = real(ibits(idx, bc_ndag_B, 1))  ! b
    c_t = real(ibits(idx, bc_ndag_T, 1))  ! t

    d_w = real(ibits(idx, bc_dn_W, 1))
    d_e = real(ibits(idx, bc_dn_E, 1))
    d_s = real(ibits(idx, bc_dn_S, 1))
    d_n = real(ibits(idx, bc_dn_N, 1))
    d_b = real(ibits(idx, bc_dn_B, 1))
    d_t = real(ibits(idx, bc_dn_T, 1))

    dsw = real(ibits(idx, bc_diag, 1))

    dd = r_x2 * (c_w + c_e) &
       + r_y2 * (c_s + c_n) &
       + r_z2 * (c_b + c_t) &
       + 2.0                &
       *(r_x2 * (d_w + d_e) &
       + r_y2 * (d_s + d_n) &
       + r_z2 * (d_b + d_t) )

    dd = dsw * dd + 1.0 - dsw

    aa = dble(ibits(idx, Active, 1))

    pp = p(i,j,k)
    bb = b(i,j,k)

    ss = r_x2 * ( c_e * p(i+1,j  ,k  ) + c_w * p(i-1,j  ,k  ) ) &
       + r_y2 * ( c_n * p(i  ,j+1,k  ) + c_s * p(i  ,j-1,k  ) ) &
       + r_z2 * ( c_t * p(i  ,j  ,k+1) + c_b * p(i  ,j  ,k-1) )

    dp = ( (ss - bb)/dd - pp ) * omg
    pn = pp + dp
    p(i,j,k) = pn

    de  = bb - (ss - pn * dd)
    res = res + dble(de*de) * aa
    xl2 = xl2 + dble(pn*pn) * aa
    err = err + dble(dp*dp) * aa
  end do
  end do
  end do
!$OMP END DO
!$OMP END PARALLEL

    cnv(1) = err
    cnv(2) = res
    cnv(3) = xl2

    return
    end subroutine psor


!> ********************************************************************
!! @brief 2-colored SOR法 stride memory access
!! @param [in,out] p     圧力
!! @param [in]     sz    配列長
!! @param [in]     g     ガイドセル長
!! @param [in]     dh    格子幅
!! @param [in]     ip    開始点インデクス
!! @param [in]     color グループ番号
!! @param [in]     omg   加速係数
!! @param [in,out] cnv   収束判定値　修正量の自乗和と残差の自乗和、解ベクトルの自乗和
!! @param [in]     b     RHS vector
!! @param [in]     bp    BCindex P
!! @param [out]    flop  浮動小数演算数
!! @note resは積算
!<
subroutine psor2sma_core (p, sz, g, dh, ip, color, omg, cnv, b, bp, flop)
implicit none
include 'ffv_f_params.h'
integer                                                   ::  i, j, k, ix, jx, kx, g, idx
integer, dimension(3)                                     ::  sz
double precision                                          ::  flop, res, err, xl2, aa
real                                                      ::  omg, dd, ss, dp, pp, bb, de, pn, dsw
real                                                      ::  c_w, c_e, c_s, c_n, c_b, c_t
real                                                      ::  d_w, d_e, d_s, d_n, d_b, d_t
real                                                      ::  r_xx, r_xy, r_xz, r_x2, r_y2, r_z2
real, dimension(3)                                        ::  dh
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  p, b
integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp
integer                                                   ::  ip, color
double precision, dimension(3)                            ::  cnv

ix = sz(1)
jx = sz(2)
kx = sz(3)

err = 0.0
res = 0.0
xl2 = 0.0

r_xx = 1.0
r_xy = dh(1) / dh(2)
r_xz = dh(1) / dh(3)
r_x2 = r_xx * r_xx
r_y2 = r_xy * r_xy
r_z2 = r_xz * r_xz

flop = flop + (dble(ix)*dble(jx)*dble(kx) * 56.0d0 +19.0d0 ) * 0.5d0


!$OMP PARALLEL &
!$OMP REDUCTION(+:res) &
!$OMP REDUCTION(+:err) &
!$OMP REDUCTION(+:xl2) &
!$OMP PRIVATE(c_w, c_e, c_s, c_n, c_b, c_t) &
!$OMP PRIVATE(d_w, d_e, d_s, d_n, d_b, d_t) &
!$OMP PRIVATE(idx, dsw, aa, dd, pp, bb, ss, dp, de, pn) &
!$OMP FIRSTPRIVATE(ix, jx, kx, color, ip, omg) &
!$OMP FIRSTPRIVATE(r_x2, r_y2, r_z2)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do j=1,jx
do i=1+mod(k+j+color+ip,2), ix, 2
  idx = bp(i,j,k)

  c_w = real(ibits(idx, bc_ndag_W, 1))  ! w
  c_e = real(ibits(idx, bc_ndag_E, 1))  ! e
  c_s = real(ibits(idx, bc_ndag_S, 1))  ! s
  c_n = real(ibits(idx, bc_ndag_N, 1))  ! n
  c_b = real(ibits(idx, bc_ndag_B, 1))  ! b
  c_t = real(ibits(idx, bc_ndag_T, 1))  ! t

  d_w = real(ibits(idx, bc_dn_W, 1))
  d_e = real(ibits(idx, bc_dn_E, 1))
  d_s = real(ibits(idx, bc_dn_S, 1))
  d_n = real(ibits(idx, bc_dn_N, 1))
  d_b = real(ibits(idx, bc_dn_B, 1))
  d_t = real(ibits(idx, bc_dn_T, 1))

  dsw = real(ibits(idx, bc_diag, 1))

  dd = r_x2 * (c_w + c_e) &
     + r_y2 * (c_s + c_n) &
     + r_z2 * (c_b + c_t) &
     + 2.0                &
     *(r_x2 * (d_w + d_e) &
     + r_y2 * (d_s + d_n) &
     + r_z2 * (d_b + d_t) )

  dd = dsw * dd + 1.0 - dsw

  aa = dble(ibits(idx, Active, 1))

  pp = p(i,j,k)
  bb = b(i,j,k)

  ss = r_x2 * ( c_e * p(i+1,j  ,k  ) + c_w * p(i-1,j  ,k  ) ) &
     + r_y2 * ( c_n * p(i  ,j+1,k  ) + c_s * p(i  ,j-1,k  ) ) &
     + r_z2 * ( c_t * p(i  ,j  ,k+1) + c_b * p(i  ,j  ,k-1) )

  dp = ( (ss - bb)/dd - pp ) * omg
  pn = pp + dp
  p(i,j,k) = pn

  de  = bb - (ss - pn * dd)
  res = res + dble(de*de) * aa
  xl2 = xl2 + dble(pn*pn) * aa
  err = err + dble(dp*dp) * aa
end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL

cnv(1) = cnv(1) + err
cnv(2) = cnv(2) + res
cnv(3) = cnv(3) + xl2

return
end subroutine psor2sma_core



!> ***********************************************************************************
!! @brief SOR2SMAの非同期通信処理
!! @param p 圧力
!! @param sz 配列長
!! @param g ガイドセル長
!! @param col オーダリングカラーの番号
!! @param ip オーダリングカラー0の最初のインデクス
!! @param cf_sz バッファサイズ
!! @param cf_x x方向のバッファ
!! @param cf_y y方向のバッファ
!! @param cf_z z方向のバッファ
!! @param key 送信ID
!<
  subroutine sma_comm(p, sz, g, col, ip, cf_sz, cf_x, cf_y, cf_z, key, nID)
  implicit none
  include 'cpm_fparam.fi'
  integer                                                ::  ix, jx, kx, g
  integer                                                ::  i, j, k, ic, icnt, ierr, iret
  integer                                                ::  col ! color No. 0 or 1
  integer                                                ::  ip  ! top index type of color0
          !  0 : color 0 start is (1,1,1)
          !  1 : color 0 start is (2,1,1)
  integer, dimension(3)                                  ::  sz, cf_sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  p 
  real, dimension(cf_sz(1), 4)                           ::  cf_x
  real, dimension(cf_sz(2), 4)                           ::  cf_y
  real, dimension(cf_sz(3), 4)                           ::  cf_z
  integer, dimension(0:5, 2)                             ::  key
  integer, dimension(0:5)                                ::  nID

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
  ic = mod(col+ip,2)
  iret = 0
  ierr = 0

  do i=X_MINUS,Z_PLUS ! (0:5)
    key(i,1) = -1 !send
    key(i,2) = -1 !recv
  end do


! X_MINUS
! send
  if( nID(X_MINUS).ge.0 ) then
    icnt = 1
    i = 1
    do k=1,kx
    do j=1+mod(k+ic+1,2),jx,2
      cf_x(icnt,1) = p(i,j,k)
      icnt = icnt+1
    end do
    end do

    call cpm_Isend(cf_x(1,1), cf_sz(1), CPM_REAL, nID(X_MINUS), 0, key(X_MINUS,1), ierr)
  endif

! recv
  if( nID(X_PLUS).ge.0 ) then
    call cpm_Irecv(cf_x(1,3), cf_sz(1), CPM_REAL, nID(X_PLUS), 0, key(X_MINUS,2), ierr)
  endif

! X_PLUS
! send
  if( nID(X_PLUS).ge.0 ) then
    icnt = 1
    i = ix
    do k=1,kx
    do j=1+mod(k+ic+ix,2),jx,2
      cf_x(icnt,2) = p(i,j,k)
      icnt = icnt+1
    end do
    end do

    call cpm_Isend(cf_x(1,2), cf_sz(1), CPM_REAL, nID(X_PLUS), 0, key(X_PLUS,1), ierr)
  endif

! recv
  if( nID(X_MINUS).ge.0 ) then
    call cpm_Irecv(cf_x(1,4), cf_sz(1), CPM_REAL, nID(X_MINUS), 0, key(X_PLUS,2), ierr)
  endif

! Y_MINUS
! send
  if( nID(Y_MINUS).ge.0 ) then
    icnt = 1
    j = 1
    do k=1,kx
    do i=1+mod(k+ic+1,2),ix,2
      cf_y(icnt,1) = p(i,j,k)
      icnt = icnt+1
    end do
    end do

    call cpm_Isend(cf_y(1,1), cf_sz(2), CPM_REAL, nID(Y_MINUS), 0, key(Y_MINUS,1), ierr)
  endif

! recv
  if( nID(Y_PLUS).ge.0 ) then
    call cpm_Irecv(cf_y(1,3), cf_sz(2), CPM_REAL, nID(Y_PLUS), 0, key(Y_MINUS,2), ierr)
  endif

! Y_PLUS
! send
  if( nID(Y_PLUS).ge.0 ) then
    icnt = 1
    j = jx
    do k=1,kx
    do i=1+mod(k+ic+jx,2),ix,2
      cf_y(icnt,2) = p(i,j,k)
      icnt = icnt+1
    end do
    end do

    call cpm_Isend(cf_y(1,2), cf_sz(2), CPM_REAL, nID(Y_PLUS), 0, key(Y_PLUS,1), ierr)
  endif

! recv
  if( nID(Y_MINUS).ge.0 ) then
    call cpm_Irecv(cf_y(1,4), cf_sz(2), CPM_REAL, nID(Y_MINUS), 0, key(Y_PLUS,2), ierr)
  endif

! Z_MINUS
! send
  if( nID(Z_MINUS).ge.0 ) then
    icnt = 1
    k = 1
    do j=1,jx
    do i=1+mod(j+ic+1,2),ix,2
      cf_z(icnt,1) = p(i,j,k)
      icnt = icnt+1
    end do
    end do

    call cpm_Isend(cf_z(1,1), cf_sz(3), CPM_REAL, nID(Z_MINUS), 0, key(Z_MINUS,1), ierr)
  endif

! recv
  if( nID(Z_PLUS).ge.0 ) then
    call cpm_Irecv(cf_z(1,3), cf_sz(3), CPM_REAL, nID(Z_PLUS), 0, key(Z_MINUS,2), ierr)
  endif

! Z_PLUS
! send
  if( nID(Z_PLUS).ge.0 ) then
    icnt = 1
    k = kx
    do j=1,jx
    do i=1+mod(j+ic+kx,2),ix,2
      cf_z(icnt,2) = p(i,j,k)
      icnt = icnt+1
    end do
    end do

    call cpm_Isend(cf_z(1,2), cf_sz(3), CPM_REAL, nID(Z_PLUS), 0, key(Z_PLUS,1), ierr)
  endif

! recv
  if( nID(Z_MINUS).ge.0 ) then
    call cpm_Irecv(cf_z(1,4), cf_sz(3), CPM_REAL, nID(Z_MINUS), 0, key(Z_PLUS,2), ierr)
  endif

  end subroutine sma_comm

 
!> ******************************************************************************
!! @brief SOR2の非同期通信処理
!! @param p 圧力
!! @param sz 配列長
!! @param g ガイドセル長
!! @param col オーダリングカラーの番号
!! @param ip オーダリングカラー0の最初のインデクス
!! @param cf_sz バッファサイズ
!! @param cf_x x方向のバッファ
!! @param cf_y y方向のバッファ
!! @param cf_z z方向のバッファ
!! @param key 送信ID
!<
  subroutine sma_comm_wait(p, sz, g, col, ip, cf_sz, cf_x, cf_y, cf_z, key)
  implicit none
  include 'cpm_fparam.fi'
  integer                                                ::  ix, jx, kx, g
  integer                                                ::  i, j, k, ic, icnt, ierr
  integer                                                ::  col ! color No. 0 or 1
  integer                                                ::  ip  ! top index type of color0
          !  0 : color 0 start is (1,1,1)
          !  1 : color 0 start is (2,1,1)
  integer, dimension(3)                                  ::  sz, cf_sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  p 
  real, dimension(cf_sz(1), 4)                           ::  cf_x
  real, dimension(cf_sz(2), 4)                           ::  cf_y
  real, dimension(cf_sz(3), 4)                           ::  cf_z
  integer, dimension(0:5, 2)                             ::  key

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
  ic = mod(col+ip,2)
  ierr = 0

! wait for recv
! from X_MINUS
  if( key(X_MINUS,2).ge.0 ) then
    call cpm_Wait(key(X_MINUS,2), ierr)
    icnt = 1
    i = ix+1
    do k=1,kx
    do j=1+mod(k+ic+ix+1,2),jx,2
      p(i,j,k) = cf_x(icnt,3)
      icnt = icnt+1
    end do
    end do
  endif

! from X_PLUS
  if( key(X_PLUS,2).ge.0 ) then
    call cpm_Wait(key(X_PLUS,2), ierr)
    icnt = 1
    i = 0
    do k=1,kx
    do j=1+mod(k+ic,2),jx,2
      p(i,j,k) = cf_x(icnt,4)
      icnt = icnt+1
    end do
    end do
  endif

! from Y_MINUS
  if( key(Y_MINUS,2).ge.0 ) then
    call cpm_Wait(key(Y_MINUS,2), ierr)
    icnt = 1
    j = jx+1
    do k=1,kx
    do i=1+mod(k+ic+jx+1,2),ix,2
      p(i,j,k) = cf_y(icnt,3)
      icnt = icnt+1
    end do
    end do
  endif

! from Y_PLUS
  if( key(Y_PLUS,2).ge.0 ) then
    call cpm_Wait(key(Y_PLUS,2), ierr)
    icnt = 1
    j = 0
    do k=1,kx
    do i=1+mod(k+ic,2),ix,2
      p(i,j,k) = cf_y(icnt,4)
      icnt = icnt+1
    end do
    end do
  endif

! from Z_MINUS
  if( key(Z_MINUS,2).ge.0 ) then
    call cpm_Wait(key(Z_MINUS,2), ierr)
    icnt = 1
    k = kx+1
    do j=1,jx
    do i=1+mod(j+ic+kx+1,2),ix,2
      p(i,j,k) = cf_z(icnt,3)
      icnt = icnt+1
    end do
    end do
  endif

! from Z_PLUS
  if( key(Z_PLUS,2).ge.0 ) then
    call cpm_Wait(key(Z_PLUS,2), ierr)
    icnt = 1
    k = 0
    do j=1,jx
    do i=1+mod(j+ic,2),ix,2
      p(i,j,k) = cf_z(icnt,4)
      icnt = icnt+1
    end do
    end do
  endif

! wait for send

  if( key(X_MINUS,1).ge.0) then
    call cpm_Wait(key(X_MINUS,1), ierr)
  endif

  if( key(X_PLUS,1).ge.0) then
    call cpm_Wait(key(X_PLUS,1), ierr)
  endif

  if( key(Y_MINUS,1).ge.0) then
    call cpm_Wait(key(Y_MINUS,1), ierr)
  endif

  if( key(Y_PLUS,1).ge.0) then
    call cpm_Wait(key(Y_PLUS,1), ierr)
  endif

  if( key(Z_MINUS,1).ge.0) then
    call cpm_Wait(key(Z_MINUS,1), ierr)
  endif

  if( key(Z_PLUS,1).ge.0) then
    call cpm_Wait(key(Z_PLUS,1), ierr)
  endif

end subroutine sma_comm_wait
