!********************************************************************
!
!   FFV : Frontflow / violet
!
!   Copyright (c) All right reserved. 2012
!
!   Institute of Industrial Science, The University of Tokyo, Japan. 
!
!********************************************************************

!> @file   ffv_poisson.f90
!! @brief  Poisson routine
!! @author kero
!<

!> ********************************************************************
!! @brief 連立一次方程式の定数項の計算
!! @param div 速度の発散値
!! @param sz 配列長
!! @param g ガイドセル長
!! @param b2 定数ベクトルの自乗和
!! @param bp BCindex P
!! @param[out] flop flop count
!<
    subroutine div_cnst (div, sz, g, b2, bp, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                     ::  i, j, k, ix, jx, kx, g
    integer, dimension(3)                                       ::  sz
    double precision                                            ::  flop
    real                                                        ::  b2, dv
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  div
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bp

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    b2 = 0.0
    
    flop = flop + dble(ix)*dble(jx)*dble(kx)*4.0d0

!$OMP PARALLEL &
!$OMP PRIVATE(dv) &
!$OMP FIRSTPRIVATE(ix, jx, kx)

#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1) &
#elif defined _STATIC
!$OMP DO SCHEDULE(static) &
#else
!$OMP DO SCHEDULE(hoge)
#endif
!$OMP REDUCTION(+:b2)
    do k=1,kx
    do j=1,jx
    do i=1,ix
      dv = div(i,j,k)*real(ibits(bp(i,j,k), Active, 1))
      b2 = b2 + dv*dv
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    return
    end subroutine div_cnst
    
!> ********************************************************************
!! @brief point SOR法
!! @param[in/out] p 圧力
!! @param sz 配列長
!! @param g ガイドセル長
!! @param omg 加速係数
!! @param[out] res 絶対残差と相対残差
!! @param src0 固定ソース項
!! @param src1 反復毎に変化するソース項
!! @param bp BCindex P
!! @param[out] flop
!<
    subroutine psor (p, sz, g, omg, res, src0, src1, bp, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, g, idx
    integer, dimension(3)                                     ::  sz
    double precision                                          ::  flop
    real                                                      ::  ndag_e, ndag_w, ndag_n, ndag_s, ndag_t, ndag_b
    real                                                      ::  omg, dd, ss, dp, res
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  p, src0, src1
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
		res = 0.0 ! absolute
    
    flop = flop + dble(ix)*dble(jx)*dble(kx)*36.0d0
    ! flop = flop + dble(ix)*dble(jx)*dble(kx)*41.0d0 ! DP

!$OMP PARALLEL &
!$OMP PRIVATE(ndag_w, ndag_e, ndag_s, ndag_n, ndag_b, ndag_t, dd, ss, dp, idx) &
!$OMP FIRSTPRIVATE(ix, jx, kx, omg)

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
      ndag_e = real(ibits(idx, bc_ndag_E, 1))  ! e, non-diagonal
      ndag_w = real(ibits(idx, bc_ndag_W, 1))  ! w 
      ndag_n = real(ibits(idx, bc_ndag_N, 1))  ! n
      ndag_s = real(ibits(idx, bc_ndag_S, 1))  ! s
      ndag_t = real(ibits(idx, bc_ndag_T, 1))  ! t
      ndag_b = real(ibits(idx, bc_ndag_B, 1))  ! b
      
      dd = 1.0 / real(ibits(idx, bc_diag, 3))  ! diagonal
      
      ss = ndag_e * p(i+1,j  ,k  ) &
         + ndag_w * p(i-1,j  ,k  ) &
         + ndag_n * p(i  ,j+1,k  ) &
         + ndag_s * p(i  ,j-1,k  ) &
         + ndag_t * p(i  ,j  ,k+1) &
         + ndag_b * p(i  ,j  ,k-1) &
         - src0(i,j,k)             &
         + src1(i,j,k)
      dp = ( dd*ss - p(i,j,k) )
      p(i,j,k) = p(i,j,k) + omg*dp
      res = res + dp*dp*real(ibits(idx, Active, 1))
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    return
    end subroutine psor

!> ********************************************************************
!! @brief 2-colored SOR法 stride memory access
!! @param[out] p 圧力
!! @param sz 配列長
!! @param g ガイドセル長
!! @param ip 開始点インデクス
!! @param color グループ番号
!! @param omg 加速係数
!! @param res 絶対残差と相対残差
!! @param src0 固定ソース項
!! @param src1 反復毎に変化するソース項
!! @param bp BCindex P
!! @param[out] flop 浮動小数演算数
!<
    subroutine psor2sma_core (p, sz, g, ip, color, omg, res, src0, src1, bp, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, g, idx
    integer, dimension(3)                                     ::  sz
    double precision                                          ::  flop
    real                                                      ::  ndag_e, ndag_w, ndag_n, ndag_s, ndag_t, ndag_b
    real                                                      ::  omg, dd, ss, dp, res, pp
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  p, src0, src1
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp
    integer                                                   ::  ip, color

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    flop = flop + dble(ix)*dble(jx)*dble(kx) * 36.0d0 * 0.5d0
    ! flop = flop + dble(ix)*dble(jx)*dble(kx) * 41.0d0 * 0.5d0 ! DP

!$OMP PARALLEL &
!$OMP PRIVATE(ndag_w, ndag_e, ndag_s, ndag_n, ndag_b, ndag_t, dd, pp, ss, dp, idx) &
!$OMP FIRSTPRIVATE(ix, jx, kx, color, ip, omg)

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
    do i=1+mod(k+j+color+ip,2), ix, 2
      idx = bp(i,j,k)
      
      ndag_e = real(ibits(idx, bc_ndag_E, 1))  ! e, non-diagonal
      ndag_w = real(ibits(idx, bc_ndag_W, 1))  ! w 
      ndag_n = real(ibits(idx, bc_ndag_N, 1))  ! n
      ndag_s = real(ibits(idx, bc_ndag_S, 1))  ! s
      ndag_t = real(ibits(idx, bc_ndag_T, 1))  ! t
      ndag_b = real(ibits(idx, bc_ndag_B, 1))  ! b
      
      dd = 1.0 / real(ibits(idx, bc_diag, 3))  ! diagonal
      pp = p(i,j,k)
      ss = ndag_e * p(i+1,j  ,k  ) &
         + ndag_w * p(i-1,j  ,k  ) &
         + ndag_n * p(i  ,j+1,k  ) &
         + ndag_s * p(i  ,j-1,k  ) &
         + ndag_t * p(i  ,j  ,k+1) &
         + ndag_b * p(i  ,j  ,k-1) &
         - src0(i,j,k)             &
         + src1(i,j,k)
      dp = (dd*ss - pp)
      p(i,j,k) = pp + omg*dp
      res = res + dp*dp*real(ibits(idx, Active, 1))
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    return
    end subroutine psor2sma_core
