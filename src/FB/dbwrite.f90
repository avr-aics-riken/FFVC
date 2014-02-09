!###################################################################################
!
! Flow Base class
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

!> @file   dbwrite.f90
!! @brief  FlowBase utilities for CIOlib
!! @author aics
!<

!> ********************************************************************
!! @brief S3D型ベクトルデータのデバッグ出力
!! @param [in] ID  ファイルディスクリプタ
!! @param [in] s   出力データ
!! @param [in] sz  配列長
!! @param [in] g   ガイドセル長
!<
    subroutine s3dwrite(ID, s, sz, g)

    implicit none
    integer                                                 :: i, j, k, ix, jx, kx, ID, g
    integer,dimension(3)                                    :: sz
    real,dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   :: s

    !write(*,*) "fort.",ID

    !write(ID,*) "sz : ",sz(1),sz(2),sz(3),g

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)

    do k=1,kx
    do j=1,jx
    do i=1,ix
      write(ID,100) i,j,k,s(i,j,k)
    enddo
    enddo
    enddo

100 format(1h ,3i4,e12.5)

    return
    end subroutine s3dwrite


!> ********************************************************************
!! @brief V3D型ベクトルデータのデバッグ出力
!! @param [in] ID  ファイルディスクリプタ
!! @param [in] v   出力データ
!! @param [in] sz  配列長
!! @param [in] g   ガイドセル長
!<
    subroutine v3dwrite(ID, v, sz, g)

    implicit none
    integer                                                 :: i, j, k, ix, jx, kx, ID, g
    integer,dimension(3)                                    :: sz
    real,dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g,3) :: v
    

    !write(*,*) "fort.",ID

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)

    do k=1,kx
    do j=1,jx
    do i=1,ix
      write(ID, 100) i,j,k,v(i,j,k,1),v(i,j,k,2),v(i,j,k,3)
    enddo
    enddo
    enddo

 100 format(1h ,3i4,3e12.5)

    return
    end subroutine v3dwrite


!> ********************************************************************
!! @brief V3DEX型ベクトルデータのデバッグ出力
!! @param [in] ID  ファイルディスクリプタ
!! @param [in] v   出力データ
!! @param [in] sz  配列長
!! @param [in] g   ガイドセル長
!<
    subroutine nv3dwrite(ID,v,sz,g)

    implicit none
    integer                                                 :: i, j, k, ix, jx, kx, ID, g
    integer,dimension(3)                                    :: sz
    real,dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g,3) :: v

    !write(*,*) "fort.",ID

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)

    do k=1,kx
    do j=1,jx
    do i=1,ix
      write(ID,100) i,j,k,v(1,i,j,k),v(2,i,j,k),v(3,i,j,k)
    enddo
    enddo
    enddo

 100 format(1h ,3i4,3e12.5)

    return
    end subroutine nv3dwrite
