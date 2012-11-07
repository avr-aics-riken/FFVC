!********************************************************************
!
!   FFV : Frontflow / violet Cartesian
!
!   Copyright (c) 2012 All right reserved.
!
!   Institute of Industrial Science, University of Tokyo, Japan. 
!
!********************************************************************

!> @file   unused.f90
!! @brief  予備の関数群
!! @author kero
!<

!> ********************************************************************
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
    subroutine ffv_psor2sma_core_if (p, sz, g, ip, color, omg, res, src0, src1, bp, flop)
    implicit none
    include '../FB/ffv_f_params.h'
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
    
    flop = flop + real(w)*36.0
    ! flop = flop + real(w)*41.0 ! DP

    return
    end subroutine ffv_psor2sma_core_if
    
!> ********************************************************************
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
    subroutine ffv_jacobi (p, sz, g, omg, res, src0, src1, bp, wk2, flop)
    implicit none
    include '../FB/ffv_f_params.h'
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
    
    flop = flop + real(ix)*real(jx)*real(kx)*36.0
    ! flop = flop + real(ix)*real(jx)*real(kx)*36.0 ! DP

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
    end subroutine ffv_jacobi

!> ********************************************************************
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
    subroutine ffv_psor_index (p, sz, g, omg, res, src0, src1, bp, index, idx_sz, flop)
    implicit none
    include '../FB/ffv_f_params.h'
    integer                                                   ::  i, j, k, m, g, idx, idx_sz, ldx
    integer, dimension(3)                                     ::  sz
    real                                                      ::  ndag_e, ndag_w, ndag_n, ndag_s, ndag_t, ndag_b
    real                                                      ::  omg, dd, ss, dp, res, flop
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  p, src0, src1
    integer, dimension(idx_sz)                                ::  index
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp

		res = 0.0 ! absolute
    
    flop = flop + real(idx_sz)*36.0
    ! flop = flop + real(idx_sz)*36.0 ! DP

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
    end subroutine ffv_psor_index

!> ********************************************************************
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
subroutine ffv_psor_index3 (p, sz, g, omg, res, src0, src1, bp, index, idx_sz, flop)
implicit none
include '../FB/ffv_f_params.h'
integer                                                   ::  i, j, k, m, g, idx
integer, dimension(3)                                     ::  sz
real                                                      ::  ndag_e, ndag_w, ndag_n, ndag_s, ndag_t, ndag_b
real                                                      ::  omg, dd, ss, dp, res, flop
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  p, src0, src1
integer                                                   ::  idx_sz
integer, dimension(3, idx_sz)                             ::  index
integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp

res = 0.0 ! absolute

flop = flop + real(idx_sz)*36.0
! flop = flop + real(idx_sz)*36.0 ! DP

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
end subroutine ffv_psor_index3

!> ********************************************************************
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
subroutine ffv_psor_if (p, sz, g, omg, res, src0, src1, bp, flop)
implicit none
include '../FB/ffv_f_params.h'
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

flop = flop + real(c)*36.0
! flop = flop + real(c)*41.0 ! DP

return
end subroutine ffv_psor_if

!> ********************************************************************
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
subroutine ffv_jacobi_if (p, sz, g, omg, res, src0, src1, bp, wk2, flop)
implicit none
include '../FB/ffv_f_params.h'
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

flop = flop + real(c)*36.0
! flop = flop + real(c)*41.0 ! DP

return
end subroutine ffv_jacobi_if