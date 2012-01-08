!   *********************************************************
!
!   SPHERE - Skeleton for PHysical and Engineering REsearch
!
!   Copyright (c) RIKEN, Japan. All right reserved. 2004-2011
!
!   *********************************************************
!
!> @file cbc_poisson.f90
!> @brief Poisson calculation
!> @author keno, FSI Team, VCAD, RIKEN
!<
!  ***************************************************
!> @subroutine cbc_div_cnst (div, sz, g, b2, bp, flop)
!! @brief 連立一次方程式の定数項の計算
!! @param div 速度の発散値
!! @param sz 配列長
!! @param g ガイドセル長
!! @param b2 定数ベクトルの自乗和
!! @param bp BCindex P
!! @param[out] flop flop count
!<
    subroutine cbc_div_cnst (div, sz, g, b2, bp, flop)
    implicit none
    include 'cbc_f_params.h'
    integer                                                     ::  i, j, k, ix, jx, kx, g
    integer, dimension(3)                                       ::  sz
    real                                                        ::  b2, dv, flop
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  div
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bp

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    b2 = 0.0
    flop = flop + real(ix*jx*kx*2) ! 20100706

    do k=1,kx
    do j=1,jx
    do i=1,ix
      dv = div(i,j,k)*real(ibits(bp(i,j,k), Active, 1))
      b2 = b2 + dv*dv
    end do
    end do
    end do

    return
    end subroutine cbc_div_cnst
    
!  ***************************************************************
!> @subroutine cbc_psor (p, sz, g, omg, res, src0, src1, bp, flop)
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
    subroutine cbc_psor (p, sz, g, omg, res, src0, src1, bp, flop)
    implicit none
    include 'cbc_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, g, idx
    integer, dimension(3)                                     ::  sz
    real                                                      ::  ndag_e, ndag_w, ndag_n, ndag_s, ndag_t, ndag_b
    real                                                      ::  omg, dd, ss, dp, res, flop
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  p, src0, src1
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
		res = 0.0 ! absolute
    
    flop = flop + real(ix*jx*kx*22)

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

    return
    end subroutine cbc_psor

!  ******************************************************************
!> @subroutine cbc_psor_if (p, sz, g, omg, res, src0, src1, bp, flop)
!! @brief point SOR 非計算セルをif-statementでスキップする
!! @param[out] p 圧力
!! @param sz 配列長
!! @param g ガイドセル長
!! @param omg 加速係数
!! @param[out] res 絶対残差と相対残差
!! @param src0 固定ソース項
!! @param src1 反復毎に変化するソース項
!! @param bp BCindex P
!! @param[out] flop
!<
    subroutine cbc_psor_if (p, sz, g, omg, res, src0, src1, bp, flop)
    implicit none
    include 'cbc_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, g, idx, c
    integer, dimension(3)                                     ::  sz
    real                                                      ::  ndag_e, ndag_w, ndag_n, ndag_s, ndag_t, ndag_b
    real                                                      ::  omg, dd, ss, dp, res, flop, f0
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  p, src0, src1
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
		res = 0.0 ! absolute
    c = 0

    do k=1,kx
    do j=1,jx
    do i=1,ix
      idx = bp(i,j,k)
      f0 = real(ibits(idx, Active, 1))
      if ( f0 == 1.0 ) then
        c = c+1
      
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
        dp = (dd*ss - p(i,j,k))
        p(i,j,k) = p(i,j,k) + omg*dp
        res = res + dp*dp*f0
      end if
    end do
    end do
    end do
    
    flop = flop + real(c*22)

    return
    end subroutine cbc_psor_if
    
!  *************************************************************************************
!> @subroutine cbc_psor_index3 (p, sz, g, omg, res, src0, src1, bp, index, idx_sz, flop)
!! @brief SOR　インデクスによる非計算セルのスキップ．インデクスはinteger*3
!! @param[out] p 圧力
!! @param sz 配列長
!! @param g ガイドセル長
!! @param omg 加速係数
!! @param res 絶対残差と相対残差
!! @param src0 固定ソース項
!! @param src1 反復毎に変化するソース項
!! @param bp BCindex P
!! @param index インデクスベクトル
!! @param idx_sz C.Fcell
!! @param[out] flop
!<
    subroutine cbc_psor_index3 (p, sz, g, omg, res, src0, src1, bp, index, idx_sz, flop)
    implicit none
    include 'cbc_f_params.h'
    integer                                                   ::  i, j, k, m, g, idx
    integer, dimension(3)                                     ::  sz
    real                                                      ::  ndag_e, ndag_w, ndag_n, ndag_s, ndag_t, ndag_b
    real                                                      ::  omg, dd, ss, dp, res, flop
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  p, src0, src1
    integer                                                   ::  idx_sz
    integer, dimension(3, idx_sz)                             ::  index
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp

		res = 0.0 ! absolute
    
    flop = flop + real(idx_sz*22)

    do m=1, idx_sz
      i = index(1, m)
      j = index(2, m)
      k = index(3, m)
      
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
      dp = (dd*ss - p(i,j,k))
      p(i,j,k) = p(i,j,k) + omg*dp
      res = res + dp*dp*real(ibits(idx, Active, 1))
    end do

    return
    end subroutine cbc_psor_index3

!  ************************************************************************************
!> @subroutine cbc_psor_index (p, sz, g, omg, res, src0, src1, bp, index, idx_sz, flop)
!! @brief SOR　インデクスによる非計算セルのスキップ．インデクスはi,j,kをビットエンコード処理
!! @param[out] p 圧力
!! @param sz 配列長
!! @param g ガイドセル長
!! @param omg 加速係数
!! @param res 絶対残差と相対残差
!! @param src0 固定ソース項
!! @param src1 反復毎に変化するソース項
!! @param bp BCindex P
!! @param index インデクスベクトル
!! @param idx_sz C.Fcell
!! @param[out] flop
!<
    subroutine cbc_psor_index (p, sz, g, omg, res, src0, src1, bp, index, idx_sz, flop)
    implicit none
    include 'cbc_f_params.h'
    integer                                                   ::  i, j, k, m, g, idx, idx_sz, ldx
    integer, dimension(3)                                     ::  sz
    real                                                      ::  ndag_e, ndag_w, ndag_n, ndag_s, ndag_t, ndag_b
    real                                                      ::  omg, dd, ss, dp, res, flop
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  p, src0, src1
    integer, dimension(idx_sz)                                ::  index
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp

		res = 0.0 ! absolute
    
    flop = flop + real(idx_sz*22)

    do m=1, idx_sz
      ldx = index(m)
      i = ibits(ldx,  0, 10)
      j = ibits(ldx, 10, 10)
      k = ibits(ldx, 20, 10)
      
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
      dp = (dd*ss - p(i,j,k))
      p(i,j,k) = p(i,j,k) + omg*dp
      res = res + dp*dp*real(ibits(idx, Active, 1))
    end do

    return
    end subroutine cbc_psor_index
    
!  **********************************************************************
!> @subroutine cbc_jacobi (p, sz, g, omg, res, src0, src1, bp, wk2, flop)
!! @brief 緩和Jacobi法
!! @param[out] p 圧力
!! @param sz 配列長
!! @param g ガイドセル長
!! @param omg 加速係数
!! @param res 絶対残差と相対残差
!! @param src0 固定ソース項
!! @param src1 反復毎に変化するソース項
!! @param wk2 ワーク用配列
!! @param bp BCindex P
!! @param[out] flop flop count
!<
    subroutine cbc_jacobi (p, sz, g, omg, res, src0, src1, bp, wk2, flop)
    implicit none
    include 'cbc_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, g, idx
    integer, dimension(3)                                     ::  sz
    real                                                      ::  ndag_e, ndag_w, ndag_n, ndag_s, ndag_t, ndag_b
    real                                                      ::  omg, dd, ss, dp, res, flop
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  p, wk2, src0, src1
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
		res = 0.0 ! absolute
    
    flop = flop + real(ix*jx*kx*22)

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
      dp = (dd*ss - p(i,j,k))
      wk2(i,j,k) = p(i,j,k) + omg*dp
      res = res + dp*dp*real(ibits(idx, Active, 1))
    end do
    end do
    end do

    do k=1,kx
    do j=1,jx
    do i=1,ix
      p(i,j,k)=wk2(i,j,k)
    end do
    end do
    end do

    return
    end subroutine cbc_jacobi

!  *************************************************************************
!> @subroutine cbc_jacobi_if (p, sz, g, omg, res, src0, src1, bp, wk2, flop)
!! @brief 緩和Jacobi法，非計算セルのスキップ
!! @param[out] p 圧力
!! @param sz 配列長
!! @param g ガイドセル長
!! @param omg 加速係数
!! @param res 絶対残差と相対残差
!! @param wrk ソース項
!! @param wk2 ワーク用配列
!! @param bp BCindex P
!! @param[out] flop flop count
!<
    subroutine cbc_jacobi_if (p, sz, g, omg, res, src0, src1, bp, wk2, flop)
    implicit none
    include 'cbc_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, g, c, idx
    integer, dimension(3)                                     ::  sz
    real                                                      ::  ndag_e, ndag_w, ndag_n, ndag_s, ndag_t, ndag_b
    real                                                      ::  omg, dd, ss, dp, res, flop, f0
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  p, wk2, src0, src1
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
		res = 0.0 ! absolute
    c = 0

    do k=1,kx
    do j=1,jx
    do i=1,ix
      idx = bp(i,j,k)
      f0 = real(ibits(idx, Active, 1))
      if ( f0 == 1.0 ) then
        c = c+1
      
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
        dp = (dd*ss - p(i,j,k))
        wk2(i,j,k) = p(i,j,k) + omg*dp
        res = res + dp*dp*f0
      end if
    end do
    end do
    end do

    do k=1,kx
    do j=1,jx
    do i=1,ix
      p(i,j,k)=wk2(i,j,k)
    end do
    end do
    end do
    
    flop = flop + real(c*22)

    return
    end subroutine cbc_jacobi_if

!  ***********************************************************************************
!> @subroutine cbc_psor2sma_core (p, sz, g, ip, color, omg, res, src0, src1, bp, flop)
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
    subroutine cbc_psor2sma_core (p, sz, g, ip, color, omg, res, src0, src1, bp, flop)
    implicit none
    include 'cbc_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, g, idx, w
    integer, dimension(3)                                     ::  sz
    real                                                      ::  ndag_e, ndag_w, ndag_n, ndag_s, ndag_t, ndag_b
    real                                                      ::  omg, dd, ss, dp, res, flop, pp
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  p, src0, src1
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp
    integer                                                   ::  ip, color

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    w = 0

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
      w = w+1
    end do
    end do
    end do
    
    flop = flop + real(w*22)

    return
    end subroutine cbc_psor2sma_core

!  **************************************************************************************
!> @subroutine cbc_psor2sma_core_if (p, sz, g, ip, color, omg, res, src0, src1, bp, flop)
!! @brief 2-colored SOR法，stride memory access 非計算セルのスキップ
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
!! @param[out] flop
!<
    subroutine cbc_psor2sma_core_if (p, sz, g, ip, color, omg, res, src0, src1, bp, flop)
    implicit none
    include 'cbc_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, g, idx, w
    integer, dimension(3)                                     ::  sz
    real                                                      ::  ndag_e, ndag_w, ndag_n, ndag_s, ndag_t, ndag_b
    real                                                      ::  omg, dd, ss, dp, res, flop, f0, pp
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  p, src0, src1
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp
    integer                                                   ::  ip, color

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
		res = 0.0

    do k=1,kx
    do j=1,jx
    do i=1+mod(k+j+color+ip,2), ix, 2
      idx = bp(i,j,k)
      f0 = real(ibits(idx, Active, 1))
      if ( f0 == 1.0 ) then
        w = w+1

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
        res = res + dp*dp*f0
      end if
    end do
    end do
    end do
    
    flop = flop + real(w*22)

    return
    end subroutine cbc_psor2sma_core_if
    
!  ***********************************************************************************
!> @subroutine cbc_sma_comm(p, sz, g, col, ip, cf_sz, cf_x, cf_y, cf_z, key, para_key)
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
!! @param para_key parallel managerの識別ID
!<
    subroutine cbc_sma_comm(p, sz, g, col, ip, cf_sz, cf_x, cf_y, cf_z, key, para_key)
    implicit none
    include 'sklparaf.h'
    integer                                                ::  ix, jx, kx, g
    integer                                                ::  i, j, k, ic, icnt, ierr, para_key, iret
    integer                                                ::  col ! color No. 0 or 1
    integer                                                ::  ip  ! top index type of color0
                                                                   !  0 : color 0 start is (1,1,1)
                                                                   !  1 : color 0 start is (2,1,1)
    integer, dimension(3)                                  ::  sz, cf_sz
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  p 
    real, dimension(cf_sz(1), 4)                           ::  cf_x
    real, dimension(cf_sz(2), 4)                           ::  cf_y
    real, dimension(cf_sz(3), 4)                           ::  cf_z
    integer, dimension(6, 2)                               ::  key
    integer, dimension(6)                                  ::  nID
!
    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    ic = mod(col+ip,2)
    iret = 0
    ierr = 0

    do i=1,6
      key(i,1) = -1 !send
      key(i,2) = -1 !recv
    end do

    ! get neighbor domain ID
    call SklGetCommID(para_key, -1,  0,  0, nID(1), SKL_DEFAULT_GROUP, ierr) !X_MINUS
    call SklGetCommID(para_key,  1,  0,  0, nID(2), SKL_DEFAULT_GROUP, ierr) !X_PLUS
    call SklGetCommID(para_key,  0, -1,  0, nID(3), SKL_DEFAULT_GROUP, ierr) !Y_MINUS
    call SklGetCommID(para_key,  0,  1,  0, nID(4), SKL_DEFAULT_GROUP, ierr) !Y_PLUS
    call SklGetCommID(para_key,  0,  0, -1, nID(5), SKL_DEFAULT_GROUP, ierr) !Z_MINUS
    call SklGetCommID(para_key,  0,  0,  1, nID(6), SKL_DEFAULT_GROUP, ierr) !Z_PLUS

    call SklIsParallel(iret)
    
! X_MINUS
    ! send
    if( nID(1).ge.0 ) then
      icnt = 1
      i = 1
      do k=1,kx
      do j=1+mod(k+ic+1,2),jx,2
        cf_x(icnt,1) = p(i,j,k)
        icnt = icnt+1
      end do
      end do
      
      if ( iret == 1 ) then
        call SklImmediateSend(cf_x(1,1), cf_sz(1), SKL_REAL, nID(1), SKL_DEFAULT_GROUP, key(1,1), ierr)
      endif
    endif

    ! recv
    if( nID(2).ge.0 ) then
      if ( iret == 1 ) then
        call SklImmediateRecv(cf_x(1,3), cf_sz(1), SKL_REAL, nID(2), SKL_DEFAULT_GROUP, key(1,2), ierr)
      end if
    endif

! X_PLUS
    ! send
    if( nID(2).ge.0 ) then
      icnt = 1
      i = ix
      do k=1,kx
      do j=1+mod(k+ic+ix,2),jx,2
        cf_x(icnt,2) = p(i,j,k)
        icnt = icnt+1
      end do
      end do
      
      if ( iret == 1 ) then
        call SklImmediateSend(cf_x(1,2), cf_sz(1), SKL_REAL, nID(2), SKL_DEFAULT_GROUP, key(2,1), ierr)
      end if
    endif

    ! recv
    if( nID(1).ge.0 ) then
      if ( iret == 1 ) then
        call SklImmediateRecv(cf_x(1,4), cf_sz(1), SKL_REAL, nID(1), SKL_DEFAULT_GROUP, key(2,2), ierr)
      endif
    endif

! Y_MINUS
    ! send
    if( nID(3).ge.0 ) then
      icnt = 1
      j = 1
      do k=1,kx
      do i=1+mod(k+ic+1,2),ix,2
        cf_y(icnt,1) = p(i,j,k)
        icnt = icnt+1
      end do
      end do
      
      if ( iret == 1 ) then
        call SklImmediateSend(cf_y(1,1), cf_sz(2), SKL_REAL, nID(3), SKL_DEFAULT_GROUP, key(3,1), ierr)
      endif
    endif

    ! recv
    if( nID(4).ge.0 ) then
      if ( iret == 1 ) then
        call SklImmediateRecv(cf_y(1,3), cf_sz(2), SKL_REAL, nID(4), SKL_DEFAULT_GROUP, key(3,2), ierr)
      endif
    endif

! Y_PLUS
    ! send
    if( nID(4).ge.0 ) then
      icnt = 1
      j = jx
      do k=1,kx
      do i=1+mod(k+ic+jx,2),ix,2
        cf_y(icnt,2) = p(i,j,k)
        icnt = icnt+1
      end do
      end do
      
      if ( iret == 1 ) then
        call SklImmediateSend(cf_y(1,2), cf_sz(2), SKL_REAL, nID(4), SKL_DEFAULT_GROUP, key(4,1), ierr)
      endif
    endif

    ! recv
    if( nID(3).ge.0 ) then
      if ( iret == 1 ) then
        call SklImmediateRecv(cf_y(1,4), cf_sz(2), SKL_REAL, nID(3), SKL_DEFAULT_GROUP, key(4,2), ierr)
      endif
    endif

! Z_MINUS
    ! send
    if( nID(5).ge.0 ) then
      icnt = 1
      k = 1
      do j=1,jx
      do i=1+mod(j+ic+1,2),ix,2
        cf_z(icnt,1) = p(i,j,k)
        icnt = icnt+1
      end do
      end do
      
      if ( iret == 1 ) then
        call SklImmediateSend(cf_z(1,1), cf_sz(3), SKL_REAL, nID(5), SKL_DEFAULT_GROUP, key(5,1), ierr)
      endif
    endif

    ! recv
    if( nID(6).ge.0 ) then
      if ( iret == 1 ) then
        call SklImmediateRecv(cf_z(1,3), cf_sz(3), SKL_REAL, nID(6), SKL_DEFAULT_GROUP, key(5,2), ierr)
      endif
    endif

! Z_PLUS
    ! send
    if( nID(6).ge.0 ) then
      icnt = 1
      k = kx
      do j=1,jx
      do i=1+mod(j+ic+kx,2),ix,2
        cf_z(icnt,2) = p(i,j,k)
        icnt = icnt+1
      end do
      end do
      
      if ( iret == 1 ) then
        call SklImmediateSend(cf_z(1,2), cf_sz(3), SKL_REAL, nID(6), SKL_DEFAULT_GROUP, key(6,1), ierr)
      endif
    endif

    ! recv
    if( nID(5).ge.0 ) then
      if ( iret == 1 ) then
        call SklImmediateRecv(cf_z(1,4), cf_sz(3), SKL_REAL, nID(5), SKL_DEFAULT_GROUP, key(6,2), ierr)
      endif
    endif

    end subroutine cbc_sma_comm

!  ******************************************************************************
!> @subroutine cbc_sma_comm_wait(p, sz, g, col, ip, cf_sz, cf_x, cf_y, cf_z, key)
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
    subroutine cbc_sma_comm_wait(p, sz, g, col, ip, cf_sz, cf_x, cf_y, cf_z, key)
    implicit none
    include 'sklparaf.h'
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
    integer, dimension(6, 2)                               ::  key

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    ic = mod(col+ip,2)
    ierr = 0

! wait for recv
    ! from X_MINUS
    if( key(1,2).ge.0 ) then
      call SklWait(key(1,2), ierr)
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
    if( key(2,2).ge.0 ) then
      call SklWait(key(2,2), ierr)
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
    if( key(3,2).ge.0 ) then
      call SklWait(key(3,2), ierr)
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
    if( key(4,2).ge.0 ) then
      call SklWait(key(4,2), ierr)
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
    if( key(5,2).ge.0 ) then
      call SklWait(key(5,2), ierr)
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
    if( key(6,2).ge.0 ) then
      call SklWait(key(6,2), ierr)
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
    if( key(1,1).ge.0) then
      call SklWait(key(1,1), ierr)
    endif
    if( key(2,1).ge.0) then
      call SklWait(key(2,1), ierr)
    endif
    if( key(3,1).ge.0) then
      call SklWait(key(3,1), ierr)
    endif
    if( key(4,1).ge.0) then
      call SklWait(key(4,1), ierr)
    endif
    if( key(5,1).ge.0) then
      call SklWait(key(5,1), ierr)
    endif
    if( key(6,1).ge.0) then
      call SklWait(key(6,1), ierr)
    endif

    end subroutine cbc_sma_comm_wait