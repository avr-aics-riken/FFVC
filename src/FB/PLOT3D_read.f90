!********************************************************************
!
!   FFV : Frontflow / violet
!
!   Copyright (c) All right reserved. 2012
!
!   Institute of Industrial Science, The University of Tokyo, Japan. 
!
!********************************************************************

!> @file   PLOT3D_read.f90
!! @brief  FlowBase PLOT3D related code
!! @author kero
!<

!  ***************************************************
!> @subroutine read_ngrid_data(ngrid, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine read_ngrid_data(ngrid, ifl)
    implicit none
    integer :: ngrid, ifl
    
    read(ifl) ngrid
    
    return
    end subroutine read_ngrid_data
 
!  ***************************************************
!> @subroutine read_ngrid_data_formatted(ngrid, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine read_ngrid_data_formatted(ngrid, ifl)
    implicit none
    integer :: ngrid, ifl
    
    read(ifl,'(i5)') ngrid
    
    return
    end subroutine read_ngrid_data_formatted

!  ***************************************************
!> @subroutine read_block_data(id, jd, kd, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine read_block_data(id, jd, kd, ifl)
    implicit none
    integer :: id, jd, kd, ifl
    
    read(ifl) id,jd,kd
    
    return
    end subroutine read_block_data

!  ***************************************************
!> @subroutine read_block_data_formatted(id, jd, kd, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine read_block_data_formatted(id, jd, kd, ifl)
    implicit none
    integer :: id, jd, kd, ifl
    
    read(ifl,'(3i5)') id,jd,kd
    
    return
    end subroutine read_block_data_formatted

!  ***************************************************
!> @subroutine read_block_data_2d(id, jd, kd, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine read_block_data_2d(id, jd, ifl)
    implicit none
    integer :: id, jd, ifl
    
    read(ifl) id,jd
    
    return
    end subroutine read_block_data_2d

!  ***************************************************
!> @subroutine read_block_data_2d_formatted(id, jd, kd, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine read_block_data_2d_formatted(id, jd, ifl)
    implicit none
    integer :: id, jd, ifl
    
    read(ifl,'(2i5)') id,jd
    
    return
    end subroutine read_block_data_2d_formatted

!  ***************************************************
!> @subroutine read_func_block_data(id, jd, kd, nvar, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine read_func_block_data(id, jd, kd, nvar, ifl)
    implicit none
    integer :: id, jd, kd, nvar, ifl
    
    read(ifl) id,jd,kd,nvar
    
    return
    end subroutine read_func_block_data

!  ***************************************************
!> @subroutine read_func_block_data_formatted(id, jd, kd, nvar, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine read_func_block_data_formatted(id, jd, kd, nvar, ifl)
    implicit none
    integer :: id, jd, kd, nvar, ifl
    
    read(ifl,'(4i5)') id,jd,kd,nvar
    
    return
    end subroutine read_func_block_data_formatted

!  ***************************************************
!> @subroutine read_xyz_3d(id, jd, kd, x, y, z, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine read_xyz_3d(id, jd, kd, x, y, z, ifl)
    implicit none
    integer :: id, jd, kd, ifl
    real, dimension(id, jd, kd) :: x, y, z
    integer :: i, j, k
    
    read(ifl) (((x(i,j,k),i=1,id),j=1,jd),k=1,kd)
    read(ifl) (((y(i,j,k),i=1,id),j=1,jd),k=1,kd)
    read(ifl) (((z(i,j,k),i=1,id),j=1,jd),k=1,kd)

    return
    end subroutine read_xyz_3d

!  ***************************************************
!> @subroutine read_xyz_3d_formatted(id, jd, kd, x, y, z, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine read_xyz_3d_formatted(id, jd, kd, x, y, z, ifl)
    implicit none
    integer :: id, jd, kd, ifl
    real, dimension(id, jd, kd) :: x, y, z
    integer :: i, j, k
    
    do k=1,kd
      read(ifl,'(10e15.6)') ((x(i,j,k),i=1,id),j=1,jd)
    end do
    do k=1,kd
      read(ifl,'(10e15.6)') ((y(i,j,k),i=1,id),j=1,jd)
    end do
    do k=1,kd
      read(ifl,'(10e15.6)') ((z(i,j,k),i=1,id),j=1,jd)
    end do

    return
    end subroutine read_xyz_3d_formatted

!  ***************************************************
!> @subroutine read_xyz_3d_iblank(id, jd, kd, x, y, z, iblank, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine read_xyz_3d_iblank(id, jd, kd, x, y, z, iblank, ifl)
    implicit none
    integer :: id, jd, kd, ifl
    real, dimension(id, jd, kd) :: x, y, z
    integer, dimension(id, jd, kd) :: iblank
    integer :: i, j, k
    
    read(ifl) (((x(i,j,k),i=1,id),j=1,jd),k=1,kd)
    read(ifl) (((y(i,j,k),i=1,id),j=1,jd),k=1,kd)
    read(ifl) (((z(i,j,k),i=1,id),j=1,jd),k=1,kd)
    read(ifl) (((iblank(i,j,k),i=1,id),j=1,jd),k=1,kd)

    return
    end subroutine read_xyz_3d_iblank

!  ***************************************************
!> @subroutine read_xyz_3d_iblank_formatted(id, jd, kd, x, y, z, iblank, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine read_xyz_3d_iblank_formatted(id, jd, kd, x, y, z, iblank, ifl)
    implicit none
    integer :: id, jd, kd, ifl
    real, dimension(id, jd, kd) :: x, y, z
!    real, dimension(id*jd*kd) :: x, y, z
    integer, dimension(id, jd, kd) :: iblank
    integer :: i, j, k

    do k=1,kd
      read(ifl,'(10e15.6)') ((x(i,j,k),i=1,id),j=1,jd)
    end do
    do k=1,kd
      read(ifl,'(10e15.6)') ((y(i,j,k),i=1,id),j=1,jd)
    end do
    do k=1,kd
      read(ifl,'(10e15.6)') ((z(i,j,k),i=1,id),j=1,jd)
    end do
    do k=1,kd
      read(ifl,'(10i2)') ((iblank(i,j,k),i=1,id),j=1,jd)
    end do

    return
    end subroutine read_xyz_3d_iblank_formatted

!  ***************************************************
!> @subroutine read_xyz_2d(id, jd, x, y, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine read_xyz_2d(id, jd, x, y, ifl)
    implicit none
    integer :: id, jd, ifl
    real, dimension(id, jd) :: x, y
    integer :: i, j
    
    read(ifl) ((x(i,j),i=1,id),j=1,jd)
    read(ifl) ((y(i,j),i=1,id),j=1,jd)
    
    return
    end subroutine read_xyz_2d

!  ***************************************************
!> @subroutine read_xyz_2d_formatted(id, jd, x, y, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine read_xyz_2d_formatted(id, jd, x, y, ifl)
    implicit none
    integer :: id, jd, ifl
    real, dimension(id, jd) :: x, y
    integer :: i, j

    read(ifl,'(10e15.6)') ((x(i,j),i=1,id),j=1,jd)
    read(ifl,'(10e15.6)') ((y(i,j),i=1,id),j=1,jd)
    
    return
    end subroutine read_xyz_2d_formatted

!  ***************************************************
!> @subroutine read_q_3d(id, jd, kd, fsmach, alpha, re, time, q, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine read_q_3d(id, jd, kd, fsmach, alpha, re, time, q, ifl)
    implicit none
    integer :: id, jd, kd, ifl
    real :: fsmach, alpha, re, time
    real, dimension(id, jd, kd, 5) :: q
    integer :: i, j, k, nx
    
    read(ifl) fsmach,alpha,re,time
    read(ifl) ((((q(i,j,k,nx),i=1,id),j=1,jd),k=1,kd),nx=1,5)

    return
    end subroutine read_q_3d

!  ***************************************************
!> @subroutine read_q_3d_formatted(id, jd, kd, fsmach, alpha, re, time, q, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine read_q_3d_formatted(id, jd, kd, fsmach, alpha, re, time, q, ifl)
    implicit none
    integer :: id, jd, kd, ifl
    real :: fsmach, alpha, re, time
    real, dimension(id, jd, kd, 5) :: q
    integer :: i, j, k, nx

    read(ifl,'(4f8.4)') fsmach,alpha,re,time
    do nx=1,5
    do k=1,kd
      read(ifl,'(10e15.6)') ((q(i,j,k,nx),i=1,id),j=1,jd)
    end do
    end do

    return
    end subroutine read_q_3d_formatted

!  ***************************************************
!> @subroutine read_q_2d(id, jd, fsmach, alpha, re, time, q, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine read_q_2d(id, jd, fsmach, alpha, re, time, q, ifl)
    implicit none
    integer :: id, jd, ifl
    real :: fsmach, alpha, re, time
    real, dimension(id, jd, 5) :: q
    integer :: i, j, nx
    
    read(ifl) fsmach,alpha,re,time
    read(ifl) (((q(i,j,nx),i=1,id),j=1,jd),nx=1,5)
    
    return
    end subroutine read_q_2d

!  ***************************************************
!> @subroutine read_q_2d_formatted(id, jd, fsmach, alpha, re, time, q, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine read_q_2d_formatted(id, jd, fsmach, alpha, re, time, q, ifl)
    implicit none
    integer :: id, jd, ifl
    real :: fsmach, alpha, re, time
    real, dimension(id, jd, 5) :: q
    integer :: i, j, nx

    read(ifl,'(4f8.4)') fsmach,alpha,re,time
    do nx=1,5
      read(ifl,'(10e15.6)') ((q(i,j,nx),i=1,id),j=1,jd)
    end do
    
    return
    end subroutine read_q_2d_formatted

!  ***************************************************
!> @subroutine read_func_3d(id, jd, kd, nvar, d, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine read_func_3d(id, jd, kd, nvar, d, ifl)
    implicit none
    integer :: id, jd, kd, nvar, ifl
    real, dimension(id, jd, kd, nvar) :: d
    integer :: i, j, k, nx
    
    read(ifl) ((((d(i,j,k,nx),i=1,id),j=1,jd),k=1,kd),nx=1,nvar)

    return
    end subroutine read_func_3d

!  ***************************************************
!> @subroutine read_func_3d_formatted(id, jd, kd, nvar, d, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine read_func_3d_formatted(id, jd, kd, nvar, d, ifl)
    implicit none
    integer :: id, jd, kd, nvar, ifl
    real, dimension(id, jd, kd, nvar) :: d
    integer :: i, j, k, nx

    do nx=1,nvar
    do k=1,kd
      read(ifl,'(10e15.6)') ((d(i,j,k,nx),i=1,id),j=1,jd)
    end do
    end do

    return
    end subroutine read_func_3d_formatted

!  ***************************************************
!> @subroutine read_func_2d(id, jd, nvar, d, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine read_func_2d(id, jd, nvar, d, ifl)
    implicit none
    integer :: id, jd, nvar, ifl
    real, dimension(id, jd, nvar) :: d
    integer :: i, j, nx
    
    read(ifl) (((d(i,j,nx),i=1,id),j=1,jd),nx=1,nvar)
    
    return
    end subroutine read_func_2d

!  ***************************************************
!> @subroutine read_func_2d_formatted(id, jd, nvar, d, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine read_func_2d_formatted(id, jd, nvar, d, ifl)
    implicit none
    integer :: id, jd, nvar, ifl
    real, dimension(id, jd, nvar) :: d
    integer :: i, j, nx

    do nx=1,nvar
    read(ifl,'(10e15.6)') ((d(i,j,nx),i=1,id),j=1,jd)
    end do

    return
    end subroutine read_func_2d_formatted

!  ****************************************************************************
!  ****************************************************************************
!  ****************************************************************************
!  ****************************************************************************
!  ****************************************************************************
!  ****************************************************************************
!  ****************************************************************************
!  ****************************************************************************
!  ****************************************************************************
!  ****************************************************************************
 
 !  ***************************************************
!> @subroutine dread_xyz_3d(id, jd, kd, x, y, z, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine dread_xyz_3d(id, jd, kd, x, y, z, ifl)
    implicit none
    integer :: id, jd, kd, ifl
    real*8, dimension(id, jd, kd) :: x, y, z
    integer :: i, j, k
    
    read(ifl) (((x(i,j,k),i=1,id),j=1,jd),k=1,kd)
    read(ifl) (((y(i,j,k),i=1,id),j=1,jd),k=1,kd)
    read(ifl) (((z(i,j,k),i=1,id),j=1,jd),k=1,kd)

    return
    end subroutine dread_xyz_3d

!  ***************************************************
!> @subroutine dread_xyz_3d_formatted(id, jd, kd, x, y, z, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine dread_xyz_3d_formatted(id, jd, kd, x, y, z, ifl)
    implicit none
    integer :: id, jd, kd, ifl
    real*8, dimension(id, jd, kd) :: x, y, z
    integer :: i, j, k

    do k=1,kd
      read(ifl,'(10e15.6)') ((x(i,j,k),i=1,id),j=1,jd)
    end do
    do k=1,kd
      read(ifl,'(10e15.6)') ((y(i,j,k),i=1,id),j=1,jd)
    end do
    do k=1,kd
      read(ifl,'(10e15.6)') ((z(i,j,k),i=1,id),j=1,jd)
    end do

    return
    end subroutine dread_xyz_3d_formatted

!  ***************************************************
!> @subroutine dread_xyz_3d_iblank(id, jd, kd, x, y, z, iblank, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine dread_xyz_3d_iblank(id, jd, kd, x, y, z, iblank, ifl)
    implicit none
    integer :: id, jd, kd, ifl
    real*8, dimension(id, jd, kd) :: x, y, z
    integer, dimension(id, jd, kd) :: iblank
    integer :: i, j, k
    
    read(ifl) (((x(i,j,k),i=1,id),j=1,jd),k=1,kd)
    read(ifl) (((y(i,j,k),i=1,id),j=1,jd),k=1,kd)
    read(ifl) (((z(i,j,k),i=1,id),j=1,jd),k=1,kd)
    read(ifl) (((iblank(i,j,k),i=1,id),j=1,jd),k=1,kd)

    return
    end subroutine dread_xyz_3d_iblank

!  ***************************************************
!> @subroutine dread_xyz_3d_iblank_formatted(id, jd, kd, x, y, z, iblank, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine dread_xyz_3d_iblank_formatted(id, jd, kd, x, y, z, iblank, ifl)
    implicit none
    integer :: id, jd, kd, ifl
    real*8, dimension(id, jd, kd) :: x, y, z
!    real*8, dimension(id*jd*kd) :: x, y, z
    integer, dimension(id, jd, kd) :: iblank
    integer :: i, j, k
    
    do k=1,kd
      read(ifl,'(10e15.6)') ((x(i,j,k),i=1,id),j=1,jd)
    end do
    do k=1,kd
      read(ifl,'(10e15.6)') ((y(i,j,k),i=1,id),j=1,jd)
    end do
    do k=1,kd
      read(ifl,'(10e15.6)') ((z(i,j,k),i=1,id),j=1,jd)
    end do
    do k=1,kd
      read(ifl,'(10i2)') ((iblank(i,j,k),i=1,id),j=1,jd)
    end do
    
    return
    end subroutine dread_xyz_3d_iblank_formatted

!  ***************************************************
!> @subroutine dread_xyz_2d(id, jd, x, y, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine dread_xyz_2d(id, jd, x, y, ifl)
    implicit none
    integer :: id, jd, ifl
    real*8, dimension(id, jd) :: x, y
    integer :: i, j
    
    read(ifl) ((x(i,j),i=1,id),j=1,jd)
    read(ifl) ((y(i,j),i=1,id),j=1,jd)
    
    return
    end subroutine dread_xyz_2d

!  ***************************************************
!> @subroutine dread_xyz_2d_formatted(id, jd, x, y, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine dread_xyz_2d_formatted(id, jd, x, y, ifl)
    implicit none
    integer :: id, jd, ifl
    real*8, dimension(id, jd) :: x, y
    integer :: i, j

    read(ifl,'(10e15.6)') ((x(i,j),i=1,id),j=1,jd)
    read(ifl,'(10e15.6)') ((y(i,j),i=1,id),j=1,jd)
    
    return
    end subroutine dread_xyz_2d_formatted

!  ***************************************************
!> @subroutine dread_q_3d(id, jd, kd, fsmach, alpha, re, time, q, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine dread_q_3d(id, jd, kd, fsmach, alpha, re, time, q, ifl)
    implicit none
    integer :: id, jd, kd, ifl
    real*8 :: fsmach, alpha, re, time
    real*8, dimension(id, jd, kd, 5) :: q
    integer :: i, j, k, nx
    
    read(ifl) fsmach,alpha,re,time
    read(ifl) ((((q(i,j,k,nx),i=1,id),j=1,jd),k=1,kd),nx=1,5)

    return
    end subroutine dread_q_3d

!  ***************************************************
!> @subroutine dread_q_3d_formatted(id, jd, kd, fsmach, alpha, re, time, q, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine dread_q_3d_formatted(id, jd, kd, fsmach, alpha, re, time, q, ifl)
    implicit none
    integer :: id, jd, kd, ifl
    real*8 :: fsmach, alpha, re, time
    real*8, dimension(id, jd, kd, 5) :: q
    integer :: i, j, k, nx

    read(ifl,'(4f8.4)') fsmach,alpha,re,time
    do nx=1,5
    do k=1,kd
      read(ifl,'(10e15.6)') ((q(i,j,k,nx),i=1,id),j=1,jd)
    end do
    end do

    return
    end subroutine dread_q_3d_formatted

!  ***************************************************
!> @subroutine dread_q_2d(id, jd, fsmach, alpha, re, time, q, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine dread_q_2d(id, jd, fsmach, alpha, re, time, q, ifl)
    implicit none
    integer :: id, jd, ifl
    real*8 :: fsmach, alpha, re, time
    real*8, dimension(id, jd, 5) :: q
    integer :: i, j, nx
    
    read(ifl) fsmach,alpha,re,time
    read(ifl) (((q(i,j,nx),i=1,id),j=1,jd),nx=1,5)
    
    return
    end subroutine dread_q_2d

!  ***************************************************
!> @subroutine dread_q_2d_formatted(id, jd, fsmach, alpha, re, time, q, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine dread_q_2d_formatted(id, jd, fsmach, alpha, re, time, q, ifl)
    implicit none
    integer :: id, jd, ifl
    real*8 :: fsmach, alpha, re, time
    real*8, dimension(id, jd, 5) :: q
    integer :: i, j, nx

    read(ifl,'(4f8.4)') fsmach,alpha,re,time
    do nx=1,5
      read(ifl,'(10e15.6)') ((q(i,j,nx),i=1,id),j=1,jd)
    end do
    
    return
    end subroutine dread_q_2d_formatted

!  ***************************************************
!> @subroutine dread_func_3d(id, jd, kd, nvar, d, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine dread_func_3d(id, jd, kd, nvar, d, ifl)
    implicit none
    integer :: id, jd, kd, nvar, ifl
    real*8, dimension(id, jd, kd, nvar) :: d
    integer :: i, j, k, nx
    
    read(ifl) ((((d(i,j,k,nx),i=1,id),j=1,jd),k=1,kd),nx=1,nvar)

    return
    end subroutine dread_func_3d

!  ***************************************************
!> @subroutine dread_func_3d_formatted(id, jd, kd, nvar, data, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine dread_func_3d_formatted(id, jd, kd, nvar, d, ifl)
    implicit none
    integer :: id, jd, kd, nvar, ifl
    real*8, dimension(id, jd, kd, nvar) :: d
    integer :: i, j, k, nx

    do nx=1,nvar
    do k=1,kd
      read(ifl,'(10e15.6)') ((d(i,j,k,nx),i=1,id),j=1,jd)
    end do
    end do

    return
    end subroutine dread_func_3d_formatted

!  ***************************************************
!> @subroutine dread_func_2d(id, jd, fsmach, alpha, re, time, q, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine dread_func_2d(id, jd, nvar, d, ifl)
    implicit none
    integer :: id, jd, nvar, ifl
    real*8, dimension(id, jd, nvar) :: d
    integer :: i, j, nx
    
    read(ifl) (((d(i,j,nx),i=1,id),j=1,jd),nx=1,nvar)
    
    return
    end subroutine dread_func_2d

!  ***************************************************
!> @subroutine dread_func_2d_formatted(id, jd, nvar, d, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine dread_func_2d_formatted(id, jd, nvar, d, ifl)
    implicit none
    integer :: id, jd, nvar, ifl
    real*8, dimension(id, jd, nvar) :: d
    integer :: i, j, nx

    do nx=1,nvar
    read(ifl,'(10e15.6)') ((d(i,j,nx),i=1,id),j=1,jd)
    end do

    return
    end subroutine dread_func_2d_formatted







