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

!> @file   ffv_pbc.f90
!! @brief  圧力境界条件
!! @author kero
!<
    

!> ********************************************************************
!! @brief 圧力の外部ディリクレ境界
!! @param [in,out] p      圧力
!! @param [in]     sz     配列長
!! @param [in]     g      ガイドセル長
!! @param [in]     m_face 外部境界面の番号
!! @param [in]     pv     値
!! @param [in]     nID    隣接ランク番号（nID[]<0の時外部境界面）
!<
    subroutine pobc_drchlt (p, sz, g, m_face, pv, nID)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, face, g, m_face
    integer, dimension(3)                                     ::  sz
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  p
    real                                                      ::  pv
    integer, dimension(0:5)                                   ::  nID

    if ( nID(m_face) >= 0 ) return

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    face = m_face

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, pv, face)

    FACES : select case (face)
    case (X_minus)

!$OMP DO SCHEDULE(static) PRIVATE(j, k)
      do k=1,kx
      do j=1,jx
        p(0,j,k) = pv
      end do
      end do
!$OMP END DO
      
    case (X_plus)

!$OMP DO SCHEDULE(static) PRIVATE(j, k)
      do k=1,kx
      do j=1,jx
        p(ix+1,j,k) = pv
      end do
      end do
!$OMP END DO
      
    case (Y_minus)

!$OMP DO SCHEDULE(static) PRIVATE(i, k)
      do k=1,kx
      do i=1,ix
        p(i,0,k) = pv
      end do
      end do
!$OMP END DO
      
    case (Y_plus)

!$OMP DO SCHEDULE(static) PRIVATE(i, k)
      do k=1,kx
      do i=1,ix
        p(i,jx+1,k) = pv
      end do
      end do
!$OMP END DO
      
    case (Z_minus)

!$OMP DO SCHEDULE(static) PRIVATE(i, j)
      do j=1,jx
      do i=1,ix
        p(i,j,0) = pv
      end do
      end do
!$OMP END DO
    
    case (Z_plus)

!$OMP DO SCHEDULE(static) PRIVATE(i, j)
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
!! @param [in,out] p      圧力
!! @param [in]     sz     配列長
!! @param [in]     g      ガイドセル長
!! @param [in]     m_face 外部境界面の番号
!! @param [in]     nID    隣接ランク番号（nID[]<0の時外部境界面）
!<
    subroutine pobc_neumann (p, sz, g, m_face, nID)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, face, g, m_face
    integer, dimension(3)                                     ::  sz
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  p
    integer, dimension(0:5)                                   ::  nID

    if ( nID(m_face) >= 0 ) return

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    face = m_face

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, face)

    FACES : select case (face)
    case (X_minus)

!$OMP DO SCHEDULE(static) PRIVATE(j, k)
      do k=1,kx
      do j=1,jx
        p(0,j,k) = p(1,j,k)
      end do
      end do
!$OMP END DO
      
    case (X_plus)

!$OMP DO SCHEDULE(static) PRIVATE(j, k)
      do k=1,kx
      do j=1,jx
        p(ix+1,j,k) = p(ix,j,k)
      end do
      end do
!$OMP END DO
      
    case (Y_minus)

!$OMP DO SCHEDULE(static) PRIVATE(i, k)
      do k=1,kx
      do i=1,ix
        p(i,0,k) = p(i,1,k)
      end do
      end do
!$OMP END DO
      
    case (Y_plus)

!$OMP DO SCHEDULE(static) PRIVATE(i, k)
      do k=1,kx
      do i=1,ix
        p(i,jx+1,k) = p(i,jx,k)
      end do
      end do
!$OMP END DO
      
    case (Z_minus)

!$OMP DO SCHEDULE(static) PRIVATE(i, j)
      do j=1,jx
      do i=1,ix
        p(i,j,0) = p(i,j,1)
      end do
      end do
!$OMP END DO
    
    case (Z_plus)

!$OMP DO SCHEDULE(static) PRIVATE(i, j)
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
