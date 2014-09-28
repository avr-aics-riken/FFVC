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

!> @file   ffv_poisson.f90
!! @brief  Poisson routine
!! @author aics
!<


!> ********************************************************************
!! @brief 圧力Poissonの定数項の計算
!! @param [out] rhs  右辺ベクトルbの自乗和
!! @param [out] b    RHS vector b
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル長
!! @param [in]  s_0  \sum {u^*}
!! @param [in]  s_1  \sum {\beta F}
!! @param [in]  bp   BCindex P
!! @param [in]  dh   格子幅
!! @param [in]  dt   時間積分幅
!! @param [out] flop flop count
!<
subroutine poi_rhs (rhs, b, sz, g, s_0, s_1, bp, dh, dt, flop)
implicit none
include 'ffv_f_params.h'
integer                                                     ::  i, j, k, ix, jx, kx, g, idx
integer, dimension(3)                                       ::  sz
double precision                                            ::  flop, rhs
real                                                        ::  dh, dt, c1, dv, dd, d0, d1, d2
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  s_0, s_1, b
integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bp

ix = sz(1)
jx = sz(2)
kx = sz(3)
rhs = 0.0
c1 = dh / dt

flop = flop + dble(ix)*dble(jx)*dble(kx)*19.0d0 + 8.0d0
! flop = flop + dble(ix)*dble(jx)*dble(kx)*20.0d0

!$OMP PARALLEL &
!$OMP REDUCTION(+:rhs) &
!$OMP PRIVATE(dv, dd, d0, d1, d2, idx) &
!$OMP FIRSTPRIVATE(ix, jx, kx, c1, dh)

!$OMP DO SCHEDULE(static)

do k=1,kx
do j=1,jx
do i=1,ix
  idx = bp(i,j,k)
  d0 = real(ibits(idx, bc_diag + 0, 1)) ! 京のコンパイラが３ビットデコードをSIMD化できれば，ibits(idx, bc_diag, 3)
  d1 = real(ibits(idx, bc_diag + 1, 1))
  d2 = real(ibits(idx, bc_diag + 2, 1))
  dd = 1.0 / (d2*4.0 + d1*2.0 + d0)  ! diagonal
  dv = dd * real( c1 * s_0(i,j,k) + dh * s_1(i,j,k) ) * real(ibits(idx, Active, 1))
  b(i,j,k) = dv ! \frac{h^2}{\Delta t} \nabla u^*
  rhs = rhs + dble(dv*dv)
end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL

return
end subroutine poi_rhs


!> ********************************************************************
!! @brief 残差計算
!! @param [out] res  残差
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル長
!! @param [in]  p    圧力
!! @param [in]  b    RHS vector
!! @param [in]  bp   BCindex P
!! @param [out] flop flop count
!<
subroutine poi_residual (res, sz, g, p, b, bp, flop)
implicit none
include 'ffv_f_params.h'
integer                                                   ::  i, j, k, ix, jx, kx, g, idx
integer, dimension(3)                                     ::  sz
double precision                                          ::  flop, res
real                                                      ::  ndag_e, ndag_w, ndag_n, ndag_s, ndag_t, ndag_b
real                                                      ::  dd, ss, dp, d0, d1, d2
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  p, b
integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp

ix = sz(1)
jx = sz(2)
kx = sz(3)
res = 0.0

flop = flop + dble(ix)*dble(jx)*dble(kx)*29.0d0
! flop = flop + dble(ix)*dble(jx)*dble(kx)*39.0d0 ! DP

!$OMP PARALLEL &
!$OMP REDUCTION(+:res) &
!$OMP PRIVATE(ndag_w, ndag_e, ndag_s, ndag_n, ndag_b, ndag_t, dd, ss, dp, idx, d0, d1, d2) &
!$OMP FIRSTPRIVATE(ix, jx, kx)

!$OMP DO SCHEDULE(static)

do k=1,kx
do j=1,jx
do i=1,ix
idx = bp(i,j,k)
ndag_e = real(ibits(idx, bc_ndag_E, 1))  ! e, non-diagonal
ndag_w = real(ibits(idx, bc_ndag_W, 1))  ! w
ndag_n = real(ibits(idx, bc_ndag_N, 1))  ! n
ndag_s = real(ibits(idx, bc_ndag_S, 1))  ! s
ndag_t = real(ibits(idx, bc_ndag_T, 1))  ! t
ndag_b = real(ibits(idx, bc_ndag_B, 1))  ! b

d0 = real(ibits(idx, bc_diag + 0, 1))
d1 = real(ibits(idx, bc_diag + 1, 1))
d2 = real(ibits(idx, bc_diag + 2, 1))
dd = 1.0 / (d2*4.0 + d1*2.0 + d0)  ! diagonal
! dd = 1.0 / real(ibits(idx, bc_diag, 3))  iff, K compiler is improved

ss =  ndag_e * p(i+1,j  ,k  ) &
    + ndag_w * p(i-1,j  ,k  ) &
    + ndag_n * p(i  ,j+1,k  ) &
    + ndag_s * p(i  ,j-1,k  ) &
    + ndag_t * p(i  ,j  ,k+1) &
    + ndag_b * p(i  ,j  ,k-1)
dp = ( b(i,j,k) + dd*ss - p(i,j,k) ) * real(ibits(idx, Active, 1))
res = res + dble(dp*dp)
end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL

return
end subroutine poi_residual



!> ********************************************************************
!! @brief point SOR法
!! @param [in,out] p    圧力
!! @param [in]     sz   配列長
!! @param [in]     g    ガイドセル長
!! @param [in]     omg  加速係数
!! @param [out]    cnv  収束判定値　修正量の自乗和と残差の自乗和、解ベクトルの自乗和
!! @param [in]     b    RHS vector
!! @param [in]     bp   BCindex P
!! @param [out]    flop flop count
!! @note Activeマスクの位置は，固体中のラプラス式を解くように，更新式にはかけず残差のみにする
!<
    subroutine psor (p, sz, g, omg, cnv, b, bp, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, g, idx
    integer, dimension(3)                                     ::  sz
    double precision                                          ::  flop, res, err, aa, xl2
    real                                                      ::  ndag_e, ndag_w, ndag_n, ndag_s, ndag_t, ndag_b
    real                                                      ::  omg, dd, ss, dp, pp, d0, d1, d2, bb, de
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  p, b
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp
    double precision, dimension(3)                            ::  cnv

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
		res = 0.0
    err = 0.0
    xl2 = 0.0

    flop = flop + dble(ix)*dble(jx)*dble(kx)*39.0d0
    ! flop = flop + dble(ix)*dble(jx)*dble(kx)*41.0d0 ! DP

!$OMP PARALLEL &
!$OMP REDUCTION(+:res) &
!$OMP REDUCTION(+:err) &
!$OMP REDUCTION(+:xl2) &
!$OMP PRIVATE(ndag_w, ndag_e, ndag_s, ndag_n, ndag_b, ndag_t) &
!$OMP PRIVATE(idx, aa, d0, d1, d2, dd, pp, bb, ss, dp, de) &
!$OMP FIRSTPRIVATE(ix, jx, kx, omg)

!$OMP DO SCHEDULE(static)

    do k=1,kx
    do j=1,jx
    do i=1,ix
      idx = bp(i,j,k)
      ndag_e = real(ibits(idx, bc_ndag_E, 1))  ! e, non-diagonal
      ndag_w = real(ibits(idx, bc_ndag_W, 1))  ! w 
      ndag_n = real(ibits(idx, bc_ndag_N, 1))  ! n
      ndag_s = real(ibits(idx, bc_ndag_S, 1))  ! s
      ndag_t = real(ibits(idx, bc_ndag_T, 1))  ! t
      ndag_b = real(ibits(idx, bc_ndag_B, 1))  ! b

      aa = dble(ibits(idx, Active, 1))

      d0 = real(ibits(idx, bc_diag + 0, 1))
      d1 = real(ibits(idx, bc_diag + 1, 1))
      d2 = real(ibits(idx, bc_diag + 2, 1))
      dd = d2*4.0 + d1*2.0 + d0  ! diagonal

      ! dd = real(ibits(idx, bc_diag, 3))  Kのコンパイラがよくなれば
      pp = p(i,j,k)
      bb = b(i,j,k)
      ss = ndag_e * p(i+1,j  ,k  ) &
         + ndag_w * p(i-1,j  ,k  ) &
         + ndag_n * p(i  ,j+1,k  ) &
         + ndag_s * p(i  ,j-1,k  ) &
         + ndag_t * p(i  ,j  ,k+1) &
         + ndag_b * p(i  ,j  ,k-1)
      dp = ( (ss - bb)/dd - pp ) * omg
      p(i,j,k) = pp + dp

      de  = bb - (ss - pp * dd)
      res = res + dble(de*de) * aa
      xl2 = xl2 + dble(pp*pp) * aa
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
!! @param [in]     ip    開始点インデクス
!! @param [in]     color グループ番号
!! @param [in]     omg   加速係数
!! @param [in,out] cnv   収束判定値　修正量の自乗和と残差の自乗和、解ベクトルの自乗和
!! @param [in]     b     RHS vector
!! @param [in]     bp    BCindex P
!! @param [out]    flop  浮動小数演算数
!! @note resは積算
!<
subroutine psor2sma_core (p, sz, g, ip, color, omg, cnv, b, bp, flop)
implicit none
include 'ffv_f_params.h'
integer                                                   ::  i, j, k, ix, jx, kx, g, idx
integer, dimension(3)                                     ::  sz
double precision                                          ::  flop, res, err, xl2, aa
real                                                      ::  ndag_e, ndag_w, ndag_n, ndag_s, ndag_t, ndag_b
real                                                      ::  omg, dd, ss, dp, pp, d0, d1, d2, bb, de
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

flop = flop + dble(ix)*dble(jx)*dble(kx) * 39.0d0 * 0.5d0
! flop = flop + dble(ix)*dble(jx)*dble(kx) * 41.0d0 * 0.5d0 ! DP

!$OMP PARALLEL &
!$OMP REDUCTION(+:res) &
!$OMP REDUCTION(+:err) &
!$OMP REDUCTION(+:xl2) &
!$OMP PRIVATE(ndag_w, ndag_e, ndag_s, ndag_n, ndag_b, ndag_t) &
!$OMP PRIVATE(idx, aa, d0, d1, d2, dd, pp, bb, ss, dp, de) &
!$OMP FIRSTPRIVATE(ix, jx, kx, color, ip, omg)

!$OMP DO SCHEDULE(static)

do k=1,kx
do j=1,jx
do i=1+mod(k+j+color+ip,2), ix, 2
  idx = bp(i,j,k)

  ndag_e = real(ibits(idx, bc_ndag_E, 1))  ! e, non-diagonal
  ndag_w = real(ibits(idx, bc_ndag_W, 1))  ! w
  ndag_n = real(ibits(idx, bc_ndag_N, 1))  ! n
  ndag_s = real(ibits(idx, bc_ndag_S, 1))  ! s
  ndag_t = real(ibits(idx, bc_ndag_T, 1))  ! t
  ndag_b = real(ibits(idx, bc_ndag_B, 1))  ! b

  aa = dble(ibits(idx, Active, 1))

  d0 = real(ibits(idx, bc_diag + 0, 1))
  d1 = real(ibits(idx, bc_diag + 1, 1))
  d2 = real(ibits(idx, bc_diag + 2, 1))
  dd = d2*4.0 + d1*2.0 + d0  ! diagonal
  ! dd = real(ibits(idx, bc_diag, 3))  Kのコンパイラがよくなれば

  pp = p(i,j,k)
  bb = b(i,j,k)
  ss = ndag_e * p(i+1,j  ,k  ) &
     + ndag_w * p(i-1,j  ,k  ) &
     + ndag_n * p(i  ,j+1,k  ) &
     + ndag_s * p(i  ,j-1,k  ) &
     + ndag_t * p(i  ,j  ,k+1) &
     + ndag_b * p(i  ,j  ,k-1)
  dp = ( (ss - bb)/dd - pp ) * omg
  p(i,j,k) = pp + dp

  de  = bb - (ss - pp * dd)
  res = res + dble(de*de) * aa
  xl2 = xl2 + dble(pp*pp) * aa
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



!> ********************************************************************
!! @brief 2-colored SOR法 stride memory access Naive Imprementation
!! @param [in,out] p     圧力
!! @param [in]     sz    配列長
!! @param [in]     g     ガイドセル長
!! @param [in]     ip    開始点インデクス
!! @param [in]     color グループ番号
!! @param [in]     omg   加速係数
!! @param [in,out] cnv   収束判定値　修正量の自乗和と残差の自乗和、解ベクトルの自乗和
!! @param [in]     b     RHS vector
!! @param [in]     bp    BCindex P
!! @param [in]     pn    係数行列
!! @param [out]    flop  浮動小数演算数
!! @note resは積算
!<
  subroutine psor2sma_naive (p, sz, g, ip, color, omg, cnv, b, bp, pn, flop)
  implicit none
  include 'ffv_f_params.h'
  integer                                                   ::  i, j, k, ix, jx, kx, g
  integer, dimension(3)                                     ::  sz
  double precision                                          ::  flop, res, err, xl2, aa
  real                                                      ::  ndag_e, ndag_w, ndag_n, ndag_s, ndag_t, ndag_b
  real                                                      ::  omg, dd, ss, dp, pp, bb, de
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  p, b
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 7) ::  pn
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp
  integer                                                   ::  ip, color
  double precision, dimension(3)                            ::  cnv

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)

  err = 0.0
  res = 0.0
  xl2 = 0.0

  flop = flop + dble(ix)*dble(jx)*dble(kx) * 34.0d0 * 0.5d0
! flop = flop + dble(ix)*dble(jx)*dble(kx) * 41.0d0 * 0.5d0 ! DP

!$OMP PARALLEL &
!$OMP REDUCTION(+:res) &
!$OMP REDUCTION(+:err) &
!$OMP REDUCTION(+:xl2) &
!$OMP PRIVATE(ndag_w, ndag_e, ndag_s, ndag_n, ndag_b, ndag_t) &
!$OMP PRIVATE(dd, aa, pp, bb, ss, dp, de) &
!$OMP FIRSTPRIVATE(ix, jx, kx, color, ip, omg)

!$OMP DO SCHEDULE(static)

  do k=1,kx
  do j=1,jx
  do i=1+mod(k+j+color+ip,2), ix, 2

    ndag_w = pn(i,j,k,1)  ! w
    ndag_e = pn(i,j,k,2)  ! e, non-diagonal
    ndag_s = pn(i,j,k,3)  ! s
    ndag_n = pn(i,j,k,4)  ! n
    ndag_b = pn(i,j,k,5)  ! b
    ndag_t = pn(i,j,k,6)  ! t
    dd     = pn(i,j,k,7)

    aa = dble(ibits(bp(i,j,k), Active, 1))

    pp = p(i,j,k)
    bb = b(i,j,k)
    ss = ndag_e * p(i+1,j  ,k  ) &
       + ndag_w * p(i-1,j  ,k  ) &
       + ndag_n * p(i  ,j+1,k  ) &
       + ndag_s * p(i  ,j-1,k  ) &
       + ndag_t * p(i  ,j  ,k+1) &
       + ndag_b * p(i  ,j  ,k-1)
    dp = ( (ss - bb)/dd - pp ) * omg
    p(i,j,k) = pp + dp

    de  = bb - (ss - pp * dd)
    res = res + dble(de*de) * aa
    xl2 = xl2 + dble(pp*pp) * aa
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
  end subroutine psor2sma_naive



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
