! #################################################################
!
! CIOlib - Cartesian Input / Output library
!
! Copyright (c) 2013 Advanced Institute for Computational Science, RIKEN.
! All right reserved.
!
! #################################################################

!> @file   cio_DFI_util.f90
!! @brief  
!! @author kero
!<

!> ********************************************************************
!! @brief データコピー
!! @param [out] dst   コピー先のデータ
!! @param [in]  sz    サイズ
!! @param [in]  dg    コピー先の仮想セル
!! @param [in]  src   コピー元のデータ 
!! @param [in]  sg    コピー元の仮想セル
!! @param [in]  ncomp 成分数
!<

subroutine cio_dfi_copy_ijkn(dst,sz,dg,src,sg,ncomp)

  implicit none 

  integer, dimension(3) :: sz 
  integer               :: dg,sg,ncomp 
  real,dimension(1-dg:sz(1)+dg, 1-dg:sz(2)+dg, 1-dg:sz(3)+dg,ncomp) :: dst
  real,dimension(1-sg:sz(1)+sg, 1-sg:sz(2)+sg, 1-sg:sz(3)+sg,ncomp) :: src

  integer               :: i,j,k,n,gc

  gc = min(dg,sg)
  do n=1,ncomp
  do k=1-gc,sz(3)+gc
  do j=1-gc,sz(2)+gc
  do i=1-gc,sz(1)+gc
    dst(i,j,k,n)=src(i,j,k,n)
  enddo
  enddo
  enddo
  enddo

  return

end subroutine cio_dfi_copy_ijkn

!> ********************************************************************
!! @brief データコピー
!! @param [out] dst   コピー先のデータ
!! @param [in]  sz    サイズ
!! @param [in]  dg    コピー先の仮想セル
!! @param [in]  src   コピー元のデータ
!! @param [in]  sg    コピー元の仮想セル
!! @param [in]  ncomp 成分数
!<

subroutine cio_dfi_copy_nijk(dst,sz,dg,src,sg,ncomp)

  implicit none

  integer, dimension(3) :: sz
  integer               :: dg,sg,ncomp
  real,dimension(ncomp,1-dg:sz(1)+dg, 1-dg:sz(2)+dg, 1-dg:sz(3)+dg) :: dst
  real,dimension(ncomp,1-sg:sz(1)+sg, 1-sg:sz(2)+sg, 1-sg:sz(3)+sg) :: src

  integer               :: i,j,k,n,gc

  gc = min(dg,sg)
  do k=1-gc,sz(3)+gc
  do j=1-gc,sz(2)+gc
  do i=1-gc,sz(1)+gc
  do n=1,ncomp
    dst(n,i,j,k)=src(n,i,j,k)
  enddo
  enddo
  enddo
  enddo

  return

end subroutine cio_dfi_copy_nijk
