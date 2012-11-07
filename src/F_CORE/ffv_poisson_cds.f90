!********************************************************************
!
!   FFV : Frontflow / violet Cartesian
!
!   Copyright (c) 2012 All right reserved.
!
!   Institute of Industrial Science, University of Tokyo, Japan.
!
!********************************************************************
!
!> @file   ffv_poisson_cds.f90
!! @brief  速度計算のルーチン群（バイナリモデル）
!! @author kero
!<
    
!  *********************************************************************
!> @subroutine cds_psor (p, sz, g, omg, res, div, bnd, cut, e, para_key)
!! @brief Poissonから導かれた連立一次方程式をSOR法で解く
!! @param p 圧力
!! @param sz 配列長
!! @param g ガイドセル長
!! @param omg 加速係数
!! @param res 残差
!! @param div ソース項
!! @param bnd マスク
!! @param cut カット情報
!! @param e マシンイプシロン
!! @param parakey parallel managerの識別ID
!! @todo
!!    - 2階微分の距離の考慮
!<
    subroutine cds_psor (p, sz, g, omg, res, div, bnd, cut, e, para_key)
    implicit none
    include 'sklparaf.h'
    integer                                                   ::  i, j, k, ix, jx, kx, g
    integer, dimension(3)                                     ::  sz
    real                                                      ::  dp, r1, e, ee, omg, cpd, res, s0, tmp
    real                                                      ::  m1, m2, m3, m4, m5, m6
    real                                                      ::  d1, d2, d3, d4, d5, d6
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  p, div, bnd
    real, dimension(6, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  cut
    integer                                                   ::  ierr, iret, para_key

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    ee = 1.0+e
    iret = 0
    ierr = 0
    
    r1  = 0.0
    res = 0.0

    do k=1,kx
    do j=1,jx
    do i=1,ix
      d1 = cut(1,i,j,k) ! Cut(i)   X+
      d2 = cut(2,i,j,k) ! Cut(i)   X-
      d3 = cut(3,i,j,k) ! Cut(j)   Y+
      d4 = cut(4,i,j,k) ! Cut(j)   Y-
      d5 = cut(5,i,j,k) ! Cut(k)   Z+
      d6 = cut(6,i,j,k) ! Cut(k)   Z-
      
      ! セル内にカットがあればゼロ，それ以外は1
      m1 = aint(d1 + 0.5)
      m2 = aint(d2 + 0.5)
      m3 = aint(d3 + 0.5)
      m4 = aint(d4 + 0.5)
      m5 = aint(d5 + 0.5)
      m6 = aint(d6 + 0.5)
      
      cpd=ee/(e+m1+m2+m3+m4+m5+m6)
      s0 = m1*p(i+1,j  ,k  )+m2*p(i-1,j  ,k  )   &
         + m3*p(i  ,j+1,k  )+m4*p(i  ,j-1,k  )   &
         + m5*p(i  ,j  ,k+1)+m6*p(i  ,j  ,k-1) + div(i,j,k)
      dp = (cpd*s0-p(i,j,k))*bnd(i,j,k)
      p(i,j,k)=p(i,j,k) + omg*dp
      r1 = r1 + dp*dp
    end do
    end do
    end do
    
    res = r1
    call SklIsParallel(iret)
    if ( iret == 1 ) then
      tmp = r1
      call SklAllreduce(tmp, res, 1, SKL_REAL, SKL_SUM, SKL_DEFAULT_GROUP, ierr)
    end if

    return
    end subroutine cds_psor
