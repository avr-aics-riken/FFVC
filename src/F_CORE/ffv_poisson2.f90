!********************************************************************
!
!   FFV : Frontflow / violet Cartesian
!
!   Copyright (c) 2012 All right reserved.
!
!   Institute of Industrial Science, University of Tokyo, Japan.
!
!********************************************************************

!> @file   ffv_poisson.f90
!! @brief  Poisson routine
!! @author kero
!<

!> ********************************************************************
!! @brief Matrix Vector product
!! @param [out] ax   Ax
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル長
!! @param [in]  p    圧力
!! @param [in]  bp   BCindex P
!! @param [out] flop flop count
!<
  subroutine matvec_p (ax, sz, g, p, bp, flop)
  implicit none
  include 'ffv_f_params.h'
  integer                                                   ::  i, j, k, ix, jx, kx, g, idx
  integer, dimension(3)                                     ::  sz
  double precision                                          ::  flop
  real                                                      ::  ndag_e, ndag_w, ndag_n, ndag_s, ndag_t, ndag_b
  real                                                      ::  dd, ss
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  p, ax
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)

  flop = flop + dble(ix)*dble(jx)*dble(kx)*30.0d0

!$OMP PARALLEL &
!$OMP PRIVATE(ndag_w, ndag_e, ndag_s, ndag_n, ndag_b, ndag_t, dd, ss, idx) &
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

    dd = 1.0 / real(ibits(idx, bc_diag, 3))  ! diagonal

    ss =  ndag_e * p(i+1,j  ,k  ) &
        + ndag_w * p(i-1,j  ,k  ) &
        + ndag_n * p(i  ,j+1,k  ) &
        + ndag_s * p(i  ,j-1,k  ) &
        + ndag_t * p(i  ,j  ,k+1) &
        + ndag_b * p(i  ,j  ,k-1)

    ax(i,j,k) = ( p(i,j,k) - dd*ss ) * real(ibits(idx, Active, 1))
  end do
  end do
  end do
!$OMP END DO
!$OMP END PARALLEL

  return
  end subroutine matvec_p



!> ********************************************************************
!! @brief 残差の計算
!! @param [in,out] rs    Axの値 / 残差 rs=b-Ax
!! @param [in]     sz    配列長
!! @param [in]     g     ガイドセル長
!! @param [out]    res_a 残差の自乗和
!! @param [in]     div   div(u) dh/dt
!! @param [in]     src   div(f)
!! @param [in]     bp    BCindex P
!! @param [out]    flop  flop count
!<
  subroutine residual (rs, sz, g, res_a, div, src, bp, flop)
  implicit none
  include 'ffv_f_params.h'
  integer                                                   ::  i, j, k, ix, jx, kx, g, idx
  integer, dimension(3)                                     ::  sz
  double precision                                          ::  flop, res_a, res
  real                                                      ::  al, dd
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  div, rs, src
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
  res = 0.0

  flop = flop + dble(ix)*dble(jx)*dble(kx)*17.0d0

!$OMP PARALLEL &
!$OMP REDUCTION(+:res) &
!$OMP PRIVATE(al, idx, dd) &
!$OMP FIRSTPRIVATE(ix, jx, kx)

!$OMP DO SCHEDULE(static)

  do k=1,kx
  do j=1,jx
  do i=1,ix
    idx = bp(i,j,k)
    dd = 1.0 / real(ibits(idx, bc_diag, 3))  ! diagonal
    al = ( -dd * (div(i,j,k) + src(i,j,k)) - rs(i,j,k) ) * real(ibits(idx, Active, 1))
    rs(i,j,k) = al
    res = res + dble(al*al)
  end do
  end do
  end do
!$OMP END DO
!$OMP END PARALLEL

  res_a = res

  return
  end subroutine residual


!> ********************************************************************
!! @brief 直交基底
!! @param [out] v    直交基底
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル長
!! @param [in]  nc   サイズ
!! @param [in]  l    インデクス
!! @param [in]  s    スカラー値
!! @param [in]  res  残差
!! @param [out] flop  flop count
!<
  subroutine orth_basis (v, sz, g, nc, l, s, res, flop)
  implicit none
  include 'ffv_f_params.h'
  integer                                                     ::  i, j, k, ix, jx, kx, g, l, nc
  integer, dimension(3)                                       ::  sz
  double precision                                            ::  flop, s
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, nc)  ::  v
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  res

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)

  flop = flop + dble(ix)*dble(jx)*dble(kx)*1.0d0

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, s, l)

!$OMP DO SCHEDULE(static)

  do k=1,kx
  do j=1,jx
  do i=1,ix
    v(i, j, k, l) = real(s) * res(i, j, k)
  end do
  end do
  end do
!$OMP END DO
!$OMP END PARALLEL

  return
  end subroutine orth_basis


!> ********************************************************************
!! @brief コピー
!! @param [out] z    前処理された直交基底?  
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル長
!! @param [in]  v    直交基底
!<
  subroutine cp_orth_basis (z, sz, g, v)
  implicit none
  include 'ffv_f_params.h'
  integer                                                     ::  i, j, k, ix, jx, kx, g
  integer, dimension(3)                                       ::  sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  z, v

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx)

!$OMP DO SCHEDULE(static)

  do k=1,kx
  do j=1,jx
  do i=1,ix
    z(i, j, k) = v(i, j, k)
  end do
  end do
  end do
!$OMP END DO
!$OMP END PARALLEL

  return
  end subroutine cp_orth_basis


!> ********************************************************************
!! @brief コピー
!! @param [out] dst  
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル長
!! @param [in]  nc   サイズ
!! @param [in]  src  被積分配列
!! @param [in]  im   列番号
!<
  subroutine copy_1 (dst, sz, g, nc, src, im)
  implicit none
  include 'ffv_f_params.h'
  integer                                                     ::  i, j, k, ix, jx, kx, g, im, nc
  integer, dimension(3)                                       ::  sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, nc)  ::  src
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  dst

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, im)

!$OMP DO SCHEDULE(static)

  do k=1,kx
  do j=1,jx
  do i=1,ix
    dst(i, j, k) = src(i, j, k, im)
  end do
  end do
  end do
!$OMP END DO
!$OMP END PARALLEL

  return
  end subroutine copy_1


!> ********************************************************************
!! @brief コピー
!! @param [out] dst
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル長
!! @param [in]  nc   サイズ
!! @param [in]  src  被積分配列
!! @param [in]  im   列番号
!<
  subroutine copy_2 (dst, sz, g, nc, src, im)
  implicit none
  include 'ffv_f_params.h'
  integer                                                     ::  i, j, k, ix, jx, kx, g, im, nc
  integer, dimension(3)                                       ::  sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, nc)  ::  dst
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  src

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, g, im)

!$OMP DO SCHEDULE(static)

  do k=1-g,kx+g
  do j=1-g,jx+g
  do i=1-g,ix+g
    dst(i, j, k, im) = src(i, j, k)
  end do
  end do
  end do
!$OMP END DO
!$OMP END PARALLEL

  return
  end subroutine copy_2


!> ********************************************************************
!! @brief 積和関数
!! @param [out] ac   和
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル長
!! @param [in]  nc   サイズ
!! @param [in]  s4   被積分配列（4次元）
!! @param [in]  s3   被積分配列（3次元）
!! @param [in]  lm   列番号
!! @param [out] flop flop count
!<
  subroutine ml_add_1 (ac, sz, g, nc, s4, s3, lm, flop)
  implicit none
  include 'ffv_f_params.h'
  integer                                                     ::  i, j, k, ix, jx, kx, g, lm, nc
  integer, dimension(3)                                       ::  sz
  double precision                                            ::  flop, ac
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, nc)  ::  s4
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  s3

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)

  flop = flop + dble(ix)*dble(jx)*dble(kx)*2.0d0

!$OMP PARALLEL &
!$OMP REDUCTION(+:ac) &
!$OMP FIRSTPRIVATE(ix, jx, kx, lm)

!$OMP DO SCHEDULE(static)

  do k=1,kx
  do j=1,jx
  do i=1,ix
    ac = ac + dble(s4(i, j, k, lm) * s3(i, j, k))
  end do
  end do
  end do
!$OMP END DO
!$OMP END PARALLEL

  return
  end subroutine ml_add_1


!> ********************************************************************
!! @brief 積和関数
!! @param [out] ac   和
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル長
!! @param [in]  s3   被積分配列（3次元）
!! @param [out] flop flop count
!<
  subroutine ml_add_2 (ac, sz, g, s3, flop)
  implicit none
  include 'ffv_f_params.h'
  integer                                                     ::  i, j, k, ix, jx, kx, g
  integer, dimension(3)                                       ::  sz
  double precision                                            ::  flop, ac
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  s3

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)

  flop = flop + dble(ix)*dble(jx)*dble(kx)*2.0d0

!$OMP PARALLEL &
!$OMP REDUCTION(+:ac) &
!$OMP FIRSTPRIVATE(ix, jx, kx)

!$OMP DO SCHEDULE(static)

  do k=1,kx
  do j=1,jx
  do i=1,ix
    ac = ac + dble(s3(i, j, k) * s3(i, j, k))
  end do
  end do
  end do
!$OMP END DO
!$OMP END PARALLEL

  return
  end subroutine ml_add_2


!> ********************************************************************
!! @brief 積和関数
!! @param [out] s3   出力
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル長
!! @param [in]  nc   サイズ
!! @param [in]  s    スカラー
!! @param [in]  s4   被積分配列（4次元）
!! @param [in]  lm   列番号
!! @param [out] flop flop count
!<
  subroutine ml_add_3 (s3, sz, g, nc, s, s4, lm, flop)
  implicit none
  include 'ffv_f_params.h'
  integer                                                     ::  i, j, k, ix, jx, kx, g, lm, nc
  integer, dimension(3)                                       ::  sz
  double precision                                            ::  flop, s
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, nc)  ::  s4
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  s3

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)

  flop = flop + dble(ix)*dble(jx)*dble(kx)*2.0d0

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, lm, s)

!$OMP DO SCHEDULE(static)

  do k=1,kx
  do j=1,jx
  do i=1,ix
    s3(i,j,k) = s3(i,j,k) + real(s) * s4(i, j, k, lm)
  end do
  end do
  end do
!$OMP END DO
!$OMP END PARALLEL

  return
  end subroutine ml_add_3


!> ********************************************************************
!! @brief 積和関数
!! @param [out] s3   出力
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル長
!! @param [in]  nc   サイズ
!! @param [in]  s    スカラー
!! @param [in]  s4   被積分配列（4次元）
!! @param [in]  lm   列番号
!! @param [out] flop flop count
!<
  subroutine ml_add_4 (s3, sz, g, nc, s, s4, lm, flop)
  implicit none
  include 'ffv_f_params.h'
  integer                                                     ::  i, j, k, ix, jx, kx, g, lm, nc
  integer, dimension(3)                                       ::  sz
  double precision                                            ::  flop, s
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, nc)  ::  s4
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  s3

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)

  flop = flop + dble(ix)*dble(jx)*dble(kx)*2.0d0

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, lm, g, s)

!$OMP DO SCHEDULE(static)

  do k=1-g,kx+g
  do j=1-g,jx+g
  do i=1-g,ix+g
    s3(i,j,k) = s3(i,j,k) + real(s) * s4(i, j, k, lm)
  end do
  end do
  end do
!$OMP END DO
!$OMP END PARALLEL

  return
  end subroutine ml_add_4
