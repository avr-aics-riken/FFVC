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

!> @file   unused.f90
!! @brief  予備の関数群
!! @author aics
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

!> ********************************************************************
!! @brief 壁関数での摩擦速度の計算
!! @param ut 摩擦速度
!! @param sz 配列長
!! @param g ガイドセル長
!! @param dh 格子幅
!! @param re レイノルズ数
!! @param v 速度ベクトル（nステップ）
!! @param bp BCindex P
!! @param range_Yp Y+の最小値と最大値
!! @param range_Ut 摩擦速度（無次元）の最小値と最大値
!! @param v00 参照速度
!! @param flop
!! @note NOCHECK
!<
subroutine friction_velocity (ut, sz, g, dh, re, v, bp, range_Yp, range_Ut, v00, flop)
implicit none
include 'ffv_f_params.h'

integer                                                   ::  i, j, k, ix, jx, kx, g, bpx, itr, itrMax, iret, ierr
integer, dimension(3)                                     ::  sz
double precision                                          ::  flop
real, dimension(2)                                        ::  range_Yp, range_Ut, tmp_Max, tmp_Min, tmp
real                                                      ::  dh, dd, dis, eps1, eps2, re
real                                                      ::  T_w, U_t, U_t0, Yp, dut
real                                                      ::  u_ref, v_ref, w_ref, u1, u2, u3, ug
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  ut
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v
real, dimension(0:3)                                      ::  v00
integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp

ix = sz(1)
jx = sz(2)
kx = sz(3)

range_Yp(1) =  1.0e6
range_Yp(2) = -1.0e6
range_Ut(1) =  1.0e6
range_Ut(2) = -1.0e6
tmp(1) = 0.0
tmp(2) = 0.0

dd    = 2.0/(dh*re)
dis   = 0.5*dh
u_ref = v00(1)
v_ref = v00(2)
w_ref = v00(3)

eps1 = 1.0e-3 ! 速度がeps1よりも小さい値のときには摩擦速度をゼロにする
eps2 = 1.0e-3 ! 反復の収束判定閾値
itrMax = 10   ! Newton反復の最大反復数
flop = flop + 4.0d0

do k=1,kx
do j=1,jx
do i=1,ix
bpx = bp(i,j,k)

u1 = v(i,j,k,1) - u_ref
u2 = v(i,j,k,2) - v_ref
u3 = v(i,j,k,3) - w_ref
ug = sqrt(u1*u1 + u2*u2 + u3*u3)
U_t0 = 0.0

flop = flop + 23.0d0

! 6面のうちのどれか隣接壁への方向フラグが立っている場合，かつ速度がeps1以上
if ( (0 /= ibits(bpx, facing_W, 6)) .and. (ug >= eps1) ) then

T_w = ug*dd
U_t0= sqrt(T_w)
Yp  = dis*U_t0*re
flop = flop + 18.0d0

do itr=1,itrMax
dut = ( ug - (5.75*log10(Yp)+5.5)*U_t0 ) / ( ug/U_t0 + 1.0/0.4 )
U_t = U_t0 + dut
if ( U_t <= eps1 ) U_t=eps1

U_t0 = U_t
Yp   = dis*U_t0*re
if ( abs(dut)/abs(U_t0) < eps2 ) exit
end do

range_Yp(1) = min(range_Yp(1), Yp)
range_Yp(2) = max(range_Yp(2), Yp)
range_Ut(1) = min(range_Ut(1), U_t0)
range_Ut(2) = max(range_Ut(2), U_t0)

flop = flop + dble(itr)*40.0d0
endif

ut(i,j,k) = U_t0

end do
end do
end do

!    call SklIsParallel(iret)
if ( iret == 1 ) then
tmp_Min(1) = range_Yp(1)
tmp_Min(2) = range_Ut(1)
!      call SklAllreduce(tmp_Min, tmp, 2, SKL_REAL, SKL_MIN, SKL_DEFAULT_GROUP, ierr)
range_Yp(1) = tmp(1)
range_Ut(1) = tmp(2)

tmp_Max(1) = range_Yp(2)
tmp_Max(2) = range_Ut(2)
!      call SklAllreduce(tmp_Max, tmp, 2, SKL_REAL, SKL_MAX, SKL_DEFAULT_GROUP, ierr)
range_Yp(2) = tmp(1)
range_Ut(2) = tmp(2)
end if

return
end subroutine friction_velocity