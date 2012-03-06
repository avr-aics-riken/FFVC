!   *********************************************************
!
!   SPHERE - Skeleton for PHysical and Engineering REsearch
!  
!   Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
!
!   *********************************************************
!
!> @file BCprs_cc.f90
!> @brief Boundary conditons for pressure
!> @author keno, FSI Team, VCAD, RIKEN
!<

!  ************************************************
!> @subroutine cbc_pobc_drchlt (p, sz, g, face, pv)
!! @brief 圧力の外部ノイマン境界
!! @param p 圧力
!! @param sz 配列長
!! @param g ガイドセル長
!! @param face 外部境界面の番号
!! @param pv 値
!<
    subroutine cbc_pobc_drchlt (p, sz, g, face, pv)
    implicit none
    include '../FB/cbc_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, face, g
    integer, dimension(3)                                     ::  sz
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  p
    real                                                      ::  pv

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    
    FACES : select case (face)
    case (X_minus)
include '../FB/omp_head.h'
!$OMP   FIRSTPRIVATE(jx,kx,pv)
      do k=1,kx
      do j=1,jx
        p(0,j,k) = pv
      end do
      end do
include '../FB/omp_tail.h'
      
    case (X_plus)
include '../FB/omp_head.h'
!$OMP   FIRSTPRIVATE(jx,kx,pv)
      do k=1,kx
      do j=1,jx
        p(ix+1,j,k) = pv
      end do
      end do
include '../FB/omp_tail.h'
      
    case (Y_minus)
include '../FB/omp_head.h'
!$OMP   FIRSTPRIVATE(ix,kx,pv)
      do k=1,kx
      do i=1,ix
        p(i,0,k) = pv
      end do
      end do
include '../FB/omp_tail.h'
      
    case (Y_plus)
include '../FB/omp_head.h'
!$OMP   FIRSTPRIVATE(ix,kx,pv)
      do k=1,kx
      do i=1,ix
        p(i,jx+1,k) = pv
      end do
      end do
include '../FB/omp_tail.h'
      
    case (Z_minus)
include '../FB/omp_head.h'
!$OMP   FIRSTPRIVATE(ix,jx,pv)
      do j=1,jx
      do i=1,ix
        p(i,j,0) = pv
      end do
      end do
include '../FB/omp_tail.h'
    
    case (Z_plus)
include '../FB/omp_head.h'
!$OMP   FIRSTPRIVATE(ix,jx,pv)
      do j=1,jx
      do i=1,ix
        p(i,j,kx+1) = pv
      end do
      end do
include '../FB/omp_tail.h'
      
    case default
    end select FACES
    
    return
    end subroutine cbc_pobc_drchlt
    
!  *********************************************
!> @subroutine cbc_pobc_neumann (p, sz, g, face)
!! @brief 圧力の外部ノイマン境界
!! @param p 圧力
!! @param sz 配列長
!! @param g ガイドセル長
!! @param face 外部境界面の番号
!<
    subroutine cbc_pobc_neumann (p, sz, g, face)
    implicit none
    include '../FB/cbc_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, face, g
    integer, dimension(3)                                     ::  sz
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  p

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    
    FACES : select case (face)
    case (X_minus)
include '../FB/omp_head.h'
!$OMP   FIRSTPRIVATE(jx,kx)
      do k=1,kx
      do j=1,jx
        p(0,j,k) = p(1,j,k)
      end do
      end do
include '../FB/omp_tail.h'
      
    case (X_plus)
include '../FB/omp_head.h'
!$OMP   FIRSTPRIVATE(jx,kx)
      do k=1,kx
      do j=1,jx
        p(ix+1,j,k) = p(ix,j,k)
      end do
      end do
include '../FB/omp_tail.h'
      
    case (Y_minus)
include '../FB/omp_head.h'
!$OMP   FIRSTPRIVATE(ix,kx)
      do k=1,kx
      do i=1,ix
        p(i,0,k) = p(i,1,k)
      end do
      end do
include '../FB/omp_tail.h'
      
    case (Y_plus)
include '../FB/omp_head.h'
!$OMP   FIRSTPRIVATE(ix,kx)
      do k=1,kx
      do i=1,ix
        p(i,jx+1,k) = p(i,jx,k)
      end do
      end do
include '../FB/omp_tail.h'
      
    case (Z_minus)
include '../FB/omp_head.h'
!$OMP   FIRSTPRIVATE(ix,jx)
      do j=1,jx
      do i=1,ix
        p(i,j,0) = p(i,j,1)
      end do
      end do
include '../FB/omp_tail.h'
    
    case (Z_plus)
include '../FB/omp_head.h'
!$OMP   FIRSTPRIVATE(ix,jx)
      do j=1,jx
      do i=1,ix
        p(i,j,kx+1) = p(i,j,kx)
      end do
      end do
include '../FB/omp_tail.h'
      
    case default
    end select FACES
    
    return
    end subroutine cbc_pobc_neumann
