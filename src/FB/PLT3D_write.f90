!********************************************************************
!
!   FFV : Frontflow / violet
!
!   Copyright (c) 2012 All right reserved.
!
!   Institute of Industrial Science, The University of Tokyo, Japan. 
!
!********************************************************************

!> @file   PLT3D_write.f90
!! @brief  FlowBase PLOT3D related code
!! @author kero
!<

!  ***************************************************
!> @subroutine write_ngrid_data(ngrid, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine write_ngrid_data(ngrid, ifl)
    implicit none
    integer :: ngrid, ifl
    
    write(ifl) ngrid
    
    return
    end subroutine write_ngrid_data

 !  ***************************************************
!> @subroutine write_ngrid_data_formatted(ngrid, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine write_ngrid_data_formatted(ngrid, ifl)
    implicit none
    integer :: ngrid, ifl
    
    write(ifl,'(i5)') ngrid
    
    return
    end subroutine write_ngrid_data_formatted

!  ***************************************************
!> @subroutine write_block_data(id, jd, kd, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine write_block_data(id, jd, kd, ifl)
    implicit none
    integer :: id, jd, kd, ifl
    
    write(ifl) id,jd,kd
    
    return
    end subroutine write_block_data

!  ***************************************************
!> @subroutine write_block_data_formatted(id, jd, kd, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine write_block_data_formatted(id, jd, kd, ifl)
    implicit none
    integer :: id, jd, kd, ifl
    
    write(ifl,'(3i5)') id,jd,kd
    
    return
    end subroutine write_block_data_formatted

!  ***************************************************
!> @subroutine write_block_data_2d(id, jd, kd, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine write_block_data_2d(id, jd, ifl)
    implicit none
    integer :: id, jd, ifl
    
    write(ifl) id,jd
    
    return
    end subroutine write_block_data_2d

!  ***************************************************
!> @subroutine write_block_data_2d_formatted(id, jd, kd, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine write_block_data_2d_formatted(id, jd, ifl)
    implicit none
    integer :: id, jd, ifl
    
    write(ifl,'(2i5)') id,jd
    
    return
    end subroutine write_block_data_2d_formatted

!  ***************************************************
!> @subroutine write_func_block_data(id, jd, kd, nvar, ngrid, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine write_func_block_data(id, jd, kd, nvar, ngrid, ifl)
    implicit none
    integer :: ngrid, ifl
    integer, dimension(ngrid) :: id, jd, kd, nvar
    integer :: igrid

    write(ifl) (id(igrid),jd(igrid),kd(igrid),nvar(igrid),igrid=1,ngrid)
    
    return
    end subroutine write_func_block_data

!  ***************************************************
!> @subroutine write_func_block_data_formatted(id, jd, kd, nvar, ngrid, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine write_func_block_data_formatted(id, jd, kd, nvar, ngrid, ifl)
    implicit none
    integer :: ngrid, ifl
    integer, dimension(ngrid) :: id, jd, kd, nvar
    integer :: igrid

    do igrid=1,ngrid
      write(ifl,'(4i5)') id(igrid),jd(igrid),kd(igrid),nvar(igrid)
    end do
    
    return
    end subroutine write_func_block_data_formatted

!  ***************************************************
!> @subroutine write_func_block_data_2d(id, jd, nvar, ngrid, ifl)
!! @brief
!! @param
!! @param[out]
!<
    subroutine write_func_block_data_2d(id, jd, nvar, ngrid, ifl)
    implicit none
    integer :: ngrid, ifl
    integer, dimension(ngrid) :: id, jd, nvar
    integer :: igrid

    write(ifl) (id(igrid),jd(igrid),nvar(igrid),igrid=1,ngrid)

    return
    end subroutine write_func_block_data_2d

!  ***************************************************
!> @subroutine write_func_block_data_2d_formatted(id, jd, nvar, ngrid, ifl)
!! @brief
!! @param
!! @param[out]
!<
    subroutine write_func_block_data_2d_formatted(id, jd, nvar, ngrid, ifl)
    implicit none
    integer :: ngrid, ifl
    integer, dimension(ngrid) :: id, jd, nvar
    integer :: igrid

    do igrid=1,ngrid
      write(ifl,'(4i5)') id(igrid),jd(igrid),nvar(igrid)
    end do

    return
    end subroutine write_func_block_data_2d_formatted

!  ***************************************************
!> @subroutine write_xyz_3d(id, jd, kd, x, y, z, ifl)
!! @brief
!! @param
!! @param[out]
!<
    subroutine write_xyz_3d(id, jd, kd, x, y, z, ifl)
    implicit none
    integer :: id, jd, kd, ifl
    real*4, dimension(id, jd, kd) :: x, y, z
    integer :: i, j, k

    write(ifl) (((x(i,j,k),i=1,id),j=1,jd),k=1,kd), &
               (((y(i,j,k),i=1,id),j=1,jd),k=1,kd), &
               (((z(i,j,k),i=1,id),j=1,jd),k=1,kd)

    return
    end subroutine write_xyz_3d

!  ***************************************************
!> @subroutine write_xyz_3d_formatted(id, jd, kd, x, y, z, ifl)
!! @brief
!! @param
!! @param[out]
!<
    subroutine write_xyz_3d_formatted(id, jd, kd, x, y, z, ifl)
    implicit none
    integer :: id, jd, kd, ifl
    real*4, dimension(id, jd, kd) :: x, y, z
    integer :: i, j, k
    
    do k=1,kd
      write(ifl,'(10e15.6)') ((x(i,j,k),i=1,id),j=1,jd)
    end do
    do k=1,kd
      write(ifl,'(10e15.6)') ((y(i,j,k),i=1,id),j=1,jd)
    end do
    do k=1,kd
      write(ifl,'(10e15.6)') ((z(i,j,k),i=1,id),j=1,jd)
    end do

    return
    end subroutine write_xyz_3d_formatted

!  ***************************************************
!> @subroutine write_xyz_3d_iblank(id, jd, kd, x, y, z, iblank, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine write_xyz_3d_iblank(id, jd, kd, x, y, z, iblank, ifl)
    implicit none
    integer :: id, jd, kd, ifl
    real*4, dimension(id, jd, kd) :: x, y, z
    integer, dimension(id, jd, kd) :: iblank
    integer :: i, j, k
    
    write(ifl) (((x(i,j,k),i=1,id),j=1,jd),k=1,kd), &
               (((y(i,j,k),i=1,id),j=1,jd),k=1,kd), &
               (((z(i,j,k),i=1,id),j=1,jd),k=1,kd), &
               (((iblank(i,j,k),i=1,id),j=1,jd),k=1,kd)

    return
    end subroutine write_xyz_3d_iblank

!  ***************************************************
!> @subroutine write_xyz_3d_iblank_formatted(id, jd, kd, x, y, z, iblank, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine write_xyz_3d_iblank_formatted(id, jd, kd, x, y, z, iblank, ifl)
    implicit none
    integer :: id, jd, kd, ifl
    real*4, dimension(id, jd, kd) :: x, y, z
    integer, dimension(id, jd, kd) :: iblank
    integer :: i, j, k

    do k=1,kd
      write(ifl,'(10e15.6)') ((x(i,j,k),i=1,id),j=1,jd)
    end do
    do k=1,kd
      write(ifl,'(10e15.6)') ((y(i,j,k),i=1,id),j=1,jd)
    end do
    do k=1,kd
      write(ifl,'(10e15.6)') ((z(i,j,k),i=1,id),j=1,jd)
    end do
    do k=1,kd
      write(ifl,'(10i2)') ((iblank(i,j,k),i=1,id),j=1,jd)
    end do

    return
    end subroutine write_xyz_3d_iblank_formatted

!  ***************************************************
!> @subroutine write_xyz_2d(id, jd, x, y, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine write_xyz_2d(id, jd, x, y, ifl)
    implicit none
    integer :: id, jd, ifl
    real*4, dimension(id, jd) :: x, y
    integer :: i, j
    
    write(ifl) ((x(i,j),i=1,id),j=1,jd), &
               ((y(i,j),i=1,id),j=1,jd)
    
    return
    end subroutine write_xyz_2d

!  ***************************************************
!> @subroutine write_xyz_2d_formatted(id, jd, x, y, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine write_xyz_2d_formatted(id, jd, x, y, ifl)
    implicit none
    integer :: id, jd, ifl
    real*4, dimension(id, jd) :: x, y
    integer :: i, j

    write(ifl,'(10e15.6)') ((x(i,j),i=1,id),j=1,jd)
    write(ifl,'(10e15.6)') ((y(i,j),i=1,id),j=1,jd)
    
    return
    end subroutine write_xyz_2d_formatted

!  ***************************************************
!> @subroutine write_q_3d(id, jd, kd, fsmach, alpha, re, time, q, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine write_q_3d(id, jd, kd, fsmach, alpha, re, time, q, ifl)
    implicit none
    integer :: id, jd, kd, ifl
    real*4 :: fsmach, alpha, re, time
    real*4, dimension(id, jd, kd, 5) :: q
    integer :: i, j, k, nx
    
    write(ifl) fsmach,alpha,re,time
    write(ifl) ((((q(i,j,k,nx),i=1,id),j=1,jd),k=1,kd),nx=1,5)

    return
    end subroutine write_q_3d

!  ***************************************************
!> @subroutine write_q_3d_formatted(id, jd, kd, fsmach, alpha, re, time, q, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine write_q_3d_formatted(id, jd, kd, fsmach, alpha, re, time, q, ifl)
    implicit none
    integer :: id, jd, kd, ifl
    real*4 :: fsmach, alpha, re, time
    real*4, dimension(id, jd, kd, 5) :: q
    integer :: i, j, k, nx

    write(ifl,'(4f8.4)') fsmach,alpha,re,time
    do nx=1,5
    do k=1,kd
      write(ifl,'(10e15.6)') ((q(i,j,k,nx),i=1,id),j=1,jd)
    end do
    end do

    return
    end subroutine write_q_3d_formatted

!  ***************************************************
!> @subroutine write_q_2d(id, jd, fsmach, alpha, re, time, q, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine write_q_2d(id, jd, fsmach, alpha, re, time, q, ifl)
    implicit none
    integer :: id, jd, ifl
    real*4 :: fsmach, alpha, re, time
    real*4, dimension(id, jd, 5) :: q
    integer :: i, j, nx
    
    write(ifl) fsmach,alpha,re,time
    write(ifl) (((q(i,j,nx),i=1,id),j=1,jd),nx=1,5)
    
    return
    end subroutine write_q_2d

!  ***************************************************
!> @subroutine write_q_2d_formatted(id, jd, fsmach, alpha, re, time, q, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine write_q_2d_formatted(id, jd, fsmach, alpha, re, time, q, ifl)
    implicit none
    integer :: id, jd, ifl
    real*4 :: fsmach, alpha, re, time
    real*4, dimension(id, jd, 5) :: q
    integer :: i, j, nx

    write(ifl,'(4f8.4)') fsmach,alpha,re,time
    do nx=1,5
      write(ifl,'(10e15.6)') ((q(i,j,nx),i=1,id),j=1,jd)
    end do
    
    return
    end subroutine write_q_2d_formatted

!  ***************************************************
!> @subroutine write_func_3d(id, jd, kd, nvar, d, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine write_func_3d(id, jd, kd, nvar, d, ifl)
    implicit none
    integer :: id, jd, kd, nvar, ifl
    real*4, dimension(id, jd, kd, nvar) :: d
    integer :: i, j, k, nx
    
    write(ifl) ((((d(i,j,k,nx),i=1,id),j=1,jd),k=1,kd),nx=1,nvar)

    return
    end subroutine write_func_3d

!  ***************************************************
!> @subroutine write_func_3d_formatted(id, jd, kd, nvar, d, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine write_func_3d_formatted(id, jd, kd, nvar, d, ifl)
    implicit none
    integer :: id, jd, kd, nvar, ifl
    real*4, dimension(id, jd, kd, nvar) :: d
    integer :: i, j, k, nx

    do nx=1,nvar
    do k=1,kd
      write(ifl,'(10e15.6)') ((d(i,j,k,nx),i=1,id),j=1,jd)
    end do
    end do

    return
    end subroutine write_func_3d_formatted

!  ***************************************************
!> @subroutine write_func_2d(id, jd, nvar, d, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine write_func_2d(id, jd, nvar, d, ifl)
    implicit none
    integer :: id, jd, nvar, ifl
    real*4, dimension(id, jd, nvar) :: d
    integer :: i, j, nx
    
    write(ifl) (((d(i,j,nx),i=1,id),j=1,jd),nx=1,nvar)
    
    return
    end subroutine write_func_2d

!  ***************************************************
!> @subroutine write_func_2d_formatted(id, jd, nvar, d, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine write_func_2d_formatted(id, jd, nvar, d, ifl)
    implicit none
    integer :: id, jd, nvar, ifl
    real*4, dimension(id, jd, nvar) :: d
    integer :: i, j, nx

    do nx=1,nvar
    write(ifl,'(10e15.6)') ((d(i,j,nx),i=1,id),j=1,jd)
    end do

    return
    end subroutine write_func_2d_formatted

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
!> @subroutine dwrite_xyz_3d(id, jd, kd, x, y, z, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine dwrite_xyz_3d(id, jd, kd, x, y, z, ifl)
    implicit none
    integer :: id, jd, kd, ifl
    real*8, dimension(id, jd, kd) :: x, y, z
    integer :: i, j, k
    
    write(ifl) (((x(i,j,k),i=1,id),j=1,jd),k=1,kd), &
               (((y(i,j,k),i=1,id),j=1,jd),k=1,kd), &
               (((z(i,j,k),i=1,id),j=1,jd),k=1,kd)

    return
    end subroutine dwrite_xyz_3d

!  ***************************************************
!> @subroutine dwrite_xyz_3d_formatted(id, jd, kd, x, y, z, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine dwrite_xyz_3d_formatted(id, jd, kd, x, y, z, ifl)
    implicit none
    integer :: id, jd, kd, ifl
    real*8, dimension(id, jd, kd) :: x, y, z
    integer :: i, j, k

    do k=1,kd
      write(ifl,'(10e15.6)') ((x(i,j,k),i=1,id),j=1,jd)
    end do
    do k=1,kd
      write(ifl,'(10e15.6)') ((y(i,j,k),i=1,id),j=1,jd)
    end do
    do k=1,kd
      write(ifl,'(10e15.6)') ((z(i,j,k),i=1,id),j=1,jd)
    end do

    return
    end subroutine dwrite_xyz_3d_formatted

!  ***************************************************
!> @subroutine dwrite_xyz_3d_iblank(id, jd, kd, x, y, z, iblank, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine dwrite_xyz_3d_iblank(id, jd, kd, x, y, z, iblank, ifl)
    implicit none
    integer :: id, jd, kd, ifl
    real*8, dimension(id, jd, kd) :: x, y, z
    integer, dimension(id, jd, kd) :: iblank
    integer :: i, j, k
    
    write(ifl) (((x(i,j,k),i=1,id),j=1,jd),k=1,kd), &
               (((y(i,j,k),i=1,id),j=1,jd),k=1,kd), &
               (((z(i,j,k),i=1,id),j=1,jd),k=1,kd), &
               (((iblank(i,j,k),i=1,id),j=1,jd),k=1,kd)

    return
    end subroutine dwrite_xyz_3d_iblank

!  ***************************************************
!> @subroutine dwrite_xyz_3d_iblank_formatted(id, jd, kd, x, y, z, iblank, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine dwrite_xyz_3d_iblank_formatted(id, jd, kd, x, y, z, iblank, ifl)
    implicit none
    integer :: id, jd, kd, ifl
    real*8, dimension(id, jd, kd) :: x, y, z
    integer, dimension(id, jd, kd) :: iblank
    integer :: i, j, k
    
    do k=1,kd
      write(ifl,'(10e15.6)') ((x(i,j,k),i=1,id),j=1,jd)
    end do
    do k=1,kd
      write(ifl,'(10e15.6)') ((y(i,j,k),i=1,id),j=1,jd)
    end do
    do k=1,kd
      write(ifl,'(10e15.6)') ((z(i,j,k),i=1,id),j=1,jd)
    end do
    do k=1,kd
      write(ifl,'(10i2)') ((iblank(i,j,k),i=1,id),j=1,jd)
    end do
    
    return
    end subroutine dwrite_xyz_3d_iblank_formatted

!  ***************************************************
!> @subroutine dwrite_xyz_2d(id, jd, x, y, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine dwrite_xyz_2d(id, jd, x, y, ifl)
    implicit none
    integer :: id, jd, ifl
    real*8, dimension(id, jd) :: x, y
    integer :: i, j
    
    write(ifl) ((x(i,j),i=1,id),j=1,jd), &
               ((y(i,j),i=1,id),j=1,jd)
    
    return
    end subroutine dwrite_xyz_2d

!  ***************************************************
!> @subroutine dwrite_xyz_2d_formatted(id, jd, x, y, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine dwrite_xyz_2d_formatted(id, jd, x, y, ifl)
    implicit none
    integer :: id, jd, ifl
    real*8, dimension(id, jd) :: x, y
    integer :: i, j

    write(ifl,'(10e15.6)') ((x(i,j),i=1,id),j=1,jd)
    write(ifl,'(10e15.6)') ((y(i,j),i=1,id),j=1,jd)
    
    return
    end subroutine dwrite_xyz_2d_formatted

!  ***************************************************
!> @subroutine dwrite_q_3d(id, jd, kd, fsmach, alpha, re, time, q, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine dwrite_q_3d(id, jd, kd, fsmach, alpha, re, time, q, ifl)
    implicit none
    integer :: id, jd, kd, ifl
    real*8 :: fsmach, alpha, re, time
    real*8, dimension(id, jd, kd, 5) :: q
    integer :: i, j, k, nx
    
    write(ifl) fsmach,alpha,re,time
    write(ifl) ((((q(i,j,k,nx),i=1,id),j=1,jd),k=1,kd),nx=1,5)

    return
    end subroutine dwrite_q_3d

!  ***************************************************
!> @subroutine dwrite_q_3d_formatted(id, jd, kd, fsmach, alpha, re, time, q, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine dwrite_q_3d_formatted(id, jd, kd, fsmach, alpha, re, time, q, ifl)
    implicit none
    integer :: id, jd, kd, ifl
    real*8 :: fsmach, alpha, re, time
    real*8, dimension(id, jd, kd, 5) :: q
    integer :: i, j, k, nx

    write(ifl,'(4f8.4)') fsmach,alpha,re,time
    do nx=1,5
    do k=1,kd
      write(ifl,'(10e15.6)') ((q(i,j,k,nx),i=1,id),j=1,jd)
    end do
    end do

    return
    end subroutine dwrite_q_3d_formatted

!  ***************************************************
!> @subroutine dwrite_q_2d(id, jd, fsmach, alpha, re, time, q, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine dwrite_q_2d(id, jd, fsmach, alpha, re, time, q, ifl)
    implicit none
    integer :: id, jd, ifl
    real*8 :: fsmach, alpha, re, time
    real*8, dimension(id, jd, 5) :: q
    integer :: i, j, nx
    
    write(ifl) fsmach,alpha,re,time
    write(ifl) (((q(i,j,nx),i=1,id),j=1,jd),nx=1,5)
    
    return
    end subroutine dwrite_q_2d

!  ***************************************************
!> @subroutine dwrite_q_2d_formatted(id, jd, fsmach, alpha, re, time, q, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine dwrite_q_2d_formatted(id, jd, fsmach, alpha, re, time, q, ifl)
    implicit none
    integer :: id, jd, ifl
    real*8 :: fsmach, alpha, re, time
    real*8, dimension(id, jd, 5) :: q
    integer :: i, j, nx

    write(ifl,'(4f8.4)') fsmach,alpha,re,time
    do nx=1,5
      write(ifl,'(10e15.6)') ((q(i,j,nx),i=1,id),j=1,jd)
    end do
    
    return
    end subroutine dwrite_q_2d_formatted

!  ***************************************************
!> @subroutine dwrite_func_3d(id, jd, kd, nvar, d, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine dwrite_func_3d(id, jd, kd, nvar, d, ifl)
    implicit none
    integer :: id, jd, kd, nvar, ifl
    real*8, dimension(id, jd, kd, nvar) :: d
    integer :: i, j, k, nx
    
    write(ifl) ((((d(i,j,k,nx),i=1,id),j=1,jd),k=1,kd),nx=1,nvar)

    return
    end subroutine dwrite_func_3d

!  ***************************************************
!> @subroutine dwrite_func_3d_formatted(id, jd, kd, nvar, data, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine dwrite_func_3d_formatted(id, jd, kd, nvar, d, ifl)
    implicit none
    integer :: id, jd, kd, nvar, ifl
    real*8, dimension(id, jd, kd, nvar) :: d
    integer :: i, j, k, nx

    do nx=1,nvar
    do k=1,kd
      write(ifl,'(10e15.6)') ((d(i,j,k,nx),i=1,id),j=1,jd)
    end do
    end do

    return
    end subroutine dwrite_func_3d_formatted

!  ***************************************************
!> @subroutine dwrite_func_2d(id, jd, fsmach, alpha, re, time, q, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine dwrite_func_2d(id, jd, nvar, d, ifl)
    implicit none
    integer :: id, jd, nvar, ifl
    real*8, dimension(id, jd, nvar) :: d
    integer :: i, j, nx
    
    write(ifl) (((d(i,j,nx),i=1,id),j=1,jd),nx=1,nvar)
    
    return
    end subroutine dwrite_func_2d

!  ***************************************************
!> @subroutine dwrite_func_2d_formatted(id, jd, nvar, d, ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine dwrite_func_2d_formatted(id, jd, nvar, d, ifl)
    implicit none
    integer :: id, jd, nvar, ifl
    real*8, dimension(id, jd, nvar) :: d
    integer :: i, j, nx

    do nx=1,nvar
    write(ifl,'(10e15.6)') ((d(i,j,nx),i=1,id),j=1,jd)
    end do

    return
    end subroutine dwrite_func_2d_formatted
