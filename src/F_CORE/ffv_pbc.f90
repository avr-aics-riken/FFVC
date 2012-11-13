!********************************************************************
!
!   FFV : Frontflow / violet Cartesian
!
!   Copyright (c) 2012 All right reserved.
!
!   Institute of Industrial Science, University of Tokyo, Japan. 
!
!********************************************************************

!> @file   ffv_pbc.f90
!! @brief  圧力境界条件
!! @author kero
!<
    

!> ********************************************************************
!! @brief 圧力の外部ディリクレ境界
!! @param p 圧力
!! @param sz 配列長
!! @param g ガイドセル長
!! @param face 外部境界面の番号
!! @param pv 値
!<
    subroutine pobc_drchlt (p, sz, g, face, pv)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, face, g
    integer, dimension(3)                                     ::  sz
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  p
    real                                                      ::  pv

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, pv)

    FACES : select case (face)
    case (X_minus)

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do j=1,jx
        p(0,j,k) = pv
      end do
      end do
!$OMP END DO
      
    case (X_plus)

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do j=1,jx
        p(ix+1,j,k) = pv
      end do
      end do
!$OMP END DO
      
    case (Y_minus)

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do i=1,ix
        p(i,0,k) = pv
      end do
      end do
!$OMP END DO
      
    case (Y_plus)

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do i=1,ix
        p(i,jx+1,k) = pv
      end do
      end do
!$OMP END DO
      
    case (Z_minus)

!$OMP DO SCHEDULE(static)
      do j=1,jx
      do i=1,ix
        p(i,j,0) = pv
      end do
      end do
!$OMP END DO
    
    case (Z_plus)

!$OMP DO SCHEDULE(static)
      do j=1,jx
      do i=1,ix
        p(i,j,kx+1) = pv
      end do
      end do
!$OMP END DO
      
    case default
    end select FACES

!$OMP END PARALLEL

    return
    end subroutine pobc_drchlt
    
    
!> ********************************************************************
!! @brief 圧力の外部ノイマン境界
!! @param p 圧力
!! @param sz 配列長
!! @param g ガイドセル長
!! @param face 外部境界面の番号
!<
    subroutine pobc_neumann (p, sz, g, face)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, face, g
    integer, dimension(3)                                     ::  sz
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  p

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx)

    FACES : select case (face)
    case (X_minus)

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do j=1,jx
        p(0,j,k) = p(1,j,k)
      end do
      end do
!$OMP END DO
      
    case (X_plus)

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do j=1,jx
        p(ix+1,j,k) = p(ix,j,k)
      end do
      end do
!$OMP END DO
      
    case (Y_minus)

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do i=1,ix
        p(i,0,k) = p(i,1,k)
      end do
      end do
!$OMP END DO
      
    case (Y_plus)

!$OMP DO SCHEDULE(static)
      do k=1,kx
      do i=1,ix
        p(i,jx+1,k) = p(i,jx,k)
      end do
      end do
!$OMP END DO
      
    case (Z_minus)

!$OMP DO SCHEDULE(static)
      do j=1,jx
      do i=1,ix
        p(i,j,0) = p(i,j,1)
      end do
      end do
!$OMP END DO
    
    case (Z_plus)

!$OMP DO SCHEDULE(static)
      do j=1,jx
      do i=1,ix
        p(i,j,kx+1) = p(i,j,kx)
      end do
      end do
!$OMP END DO
      
    case default
    end select FACES

!$OMP END PARALLEL

    return
    end subroutine pobc_neumann
