!********************************************************************
!
!   FFV : Frontflow / violet Cartesian
!
!   Copyright (c) 2012 All right reserved.
!
!   Institute of Industrial Science, The University of Tokyo, Japan.
!
!********************************************************************

!> @file   ffv_poisson.f90
!! @brief  Poisson routine
!! @author kero
!<

!> ********************************************************************
!! @brief Matrix x Vector product
!! @param [out] ax   Ax
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル長
!! @param [in]  p    圧力
!! @param [in]  bp   BCindex P
!! @param [out] flop flop count
!<
  subroutine mv_prod (ax, sz, g, p, bp, flop)
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

  flop = flop + dble(ix)*dble(jx)*dble(kx)*15.0d0

!$OMP PARALLEL &
!$OMP PRIVATE(ndag_w, ndag_e, ndag_s, ndag_n, ndag_b, ndag_t, dd, ss, idx) &
!$OMP FIRSTPRIVATE(ix, jx, kx)

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
  end subroutine mv_prod



!> ********************************************************************
!! @brief 残差の計算
!! @param [out] rest  残差 b-Ax
!! @param [in]  sz    配列長
!! @param [in]  g     ガイドセル長
!! @param [out] res_a 残差のノルム
!! @param [in]  div   rhs
!! @param [in]  ax    Ax
!! @param [in]  bp    BCindex P
!! @param [out] flop  flop count
!<
  subroutine res_smpl (rest, sz, g, res_a, div, ax, bp, flop)
  implicit none
  include 'ffv_f_params.h'
  integer                                                   ::  i, j, k, ix, jx, kx, g, idx
  integer, dimension(3)                                     ::  sz
  double precision                                          ::  flop, res_a, res
  real                                                      ::  al, dd
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  div, ax, rest
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
  res = 0.0

  flop = flop + dble(ix)*dble(jx)*dble(kx)*7.0d0

!$OMP PARALLEL &
!$OMP PRIVATE(al, idx, dd) &
!$OMP FIRSTPRIVATE(ix, jx, kx)

#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1) &
#elif defined _STATIC
!$OMP DO SCHEDULE(static) &
#else
!$OMP DO SCHEDULE(hoge)
#endif
!$OMP REDUCTION(+:res)

  do k=1,kx
  do j=1,jx
  do i=1,ix
    idx = bp(i,j,k)
    dd = 1.0 / real(ibits(idx, bc_diag, 3))  ! diagonal
    al = ( -dd * div(i,j,k) - ax(i,j,k) ) * real(ibits(idx, Active, 1))
    rest(i,j,k) = al
    res = res + dble(al*al)
  end do
  end do
  end do
!$OMP END DO
!$OMP END PARALLEL

  res_a = res

  return
  end subroutine res_smpl


!> ********************************************************************
!! @brief 乗算関数
!! @param [out] dst  目的の積
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル長
!! @param [in]  nc   サイズ
!! @param [in]  s    スカラー値
!! @param [in]  src  被積分配列
!! @param [out] flop  flop count
!<
  subroutine multiply (dst, sz, g, s, src, flop)
  implicit none
  include 'ffv_f_params.h'
  integer                                                     ::  i, j, k, ix, jx, kx, g
  integer, dimension(3)                                       ::  sz
  double precision                                            ::  flop, s
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 15)  ::  dst
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  src

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)

  flop = flop + dble(ix)*dble(jx)*dble(kx)*1.0d0

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, s)

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
    dst(i, j, k, 1) = real(s) * src(i,j,k)
  end do
  end do
  end do
!$OMP END DO
!$OMP END PARALLEL

  return
  end subroutine multiply


!> ********************************************************************
!! @brief コピー
!! @param [out] dst  
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル長
!! @param [in]  src  被積分配列
!! @param [in]  im   列番号
!<
  subroutine copy_1 (dst, sz, g, src, im)
  implicit none
  include 'ffv_f_params.h'
  integer                                                     ::  i, j, k, ix, jx, kx, g, im
  integer, dimension(3)                                       ::  sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 15)  ::  src
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  dst

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, im)

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
!! @param [in]  src  被積分配列
!! @param [in]  im   列番号
!<
  subroutine copy_2 (dst, sz, g, src, im)
  implicit none
  include 'ffv_f_params.h'
  integer                                                     ::  i, j, k, ix, jx, kx, g, im
  integer, dimension(3)                                       ::  sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 15)  ::  dst
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  src

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, im)

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
!! @param [in]  s4   被積分配列（4次元）
!! @param [in]  s3   被積分配列（3次元）
!! @param [in]  lm   列番号
!! @param [out] flop flop count
!<
  subroutine ml_add_1 (ac, sz, g, s4, s3, lm, flop)
  implicit none
  include 'ffv_f_params.h'
  integer                                                     ::  i, j, k, ix, jx, kx, g, lm
  integer, dimension(3)                                       ::  sz
  double precision                                            ::  flop, ac
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 15)  ::  s4
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  s3

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)

  flop = flop + dble(ix)*dble(jx)*dble(kx)*2.0d0

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, lm)

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
!! @param [out] s3   出力
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル長
!! @param [in]  s    スカラー
!! @param [in]  s4   被積分配列（4次元）
!! @param [in]  lm   列番号
!! @param [out] flop flop count
!<
  subroutine ml_add_3 (s3, sz, g, s, s4, lm, flop)
  implicit none
  include 'ffv_f_params.h'
  integer                                                     ::  i, j, k, ix, jx, kx, g, lm
  integer, dimension(3)                                       ::  sz
  double precision                                            ::  flop, s
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 15)  ::  s4
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  s3

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)

  flop = flop + dble(ix)*dble(jx)*dble(kx)*2.0d0

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, lm)

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
!$OMP FIRSTPRIVATE(ix, jx, kx)

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
    ac = ac + dble(s3(i, j, k) * s3(i, j, k))
  end do
  end do
  end do
!$OMP END DO
!$OMP END PARALLEL

  return
  end subroutine ml_add_2
