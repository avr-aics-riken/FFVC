! #################################################################
!
! CIOlib - Cartesian Input / Output library
!
! Copyright (c) AICS, RIKEN. All right reserved. 2013
!
! #################################################################

!> @file   cio_DFI_COARSE.f90
!! @brief  
!! @author kero
!<

!> ********************************************************************
!! @brief 粗い格子から密な格子への補間（ゼロ次）
!! @param [out] dst   密な格子系
!! @param [in]  dhead dstの開始インデックス（グローバル）
!! @param [in]  dtail dstの終了インデックス（グローバル）
!! @param [in]  dg    dstの仮想セル
!! @param [in]  src   粗い格子系
!! @param [in]  shead srcの開始インデックス（グローバル）
!! @param [in]  stail srcの終了インデックス（グローバル）
!! @param [in]  sg    srcの仮想セル
!! @param [in]  ncomp 成分数
!! @param [in]  ista  補間の開始インデクス
!! @param [in]  iend  補間の終了インデクス
!<

subroutine cio_dfi_coarse_ijkn(dst,dhead,dtail,dg,src,shead,stail,sg,ncomp,ista,iend,irank)

  implicit none 

  integer               :: irank

  integer, dimension(3) :: dhead,dtail,shead,stail,ista,iend
  integer               :: dg,sg,ncomp 
  real,dimension(dhead(1)-dg:dtail(1)+dg, dhead(2)-dg:dtail(2)+dg, dhead(3)-dg:dtail(3)+dg,ncomp) :: dst
  real,dimension(shead(1)-sg:stail(1)+sg, shead(2)-sg:stail(2)+sg, shead(3)-sg:stail(3)+sg,ncomp) :: src

  integer               :: i,j,k,n,ii,jj,kk
  real                  :: q

  integer               :: start_i,start_j,start_k

  write(*,*) "rank : ",irank," **** cio_dfi_coarse_ijkn "

  start_i=ista(1)+1
  start_j=ista(2)+1
  start_k=ista(3)+1
  if( start_i.le.0 ) start_i=1
  if( start_j.le.0 ) start_j=1
  if( start_k.le.0 ) start_k=1

  if( ista(3).lt.0 ) then
    do n=1,ncomp
    do k=0,ista(3)+1,-1
      kk=(k-1)/2
    do j=start_j,iend(2)+1
      jj=(j-1)/2+1
    do i=start_i,iend(1)+1
      ii=(i-1)/2+1
      dst(i,j,k,n)=src(ii,jj,kk,n)
    enddo
    enddo
    enddo
    enddo
  endif

  if( ista(2).lt.0 ) then
    do n=1,ncomp
    do k=start_k,iend(3)+1
      kk=(k-1)/2+1
    do j=0,ista(2)+1,-1
      jj=(j-1)/2
    do i=start_i,iend(1)+1
      ii=(i-1)/2+1
      dst(i,j,k,n)=src(ii,jj,kk,n)
    enddo
    enddo
    enddo
    enddo
  endif

  if( ista(1).lt.0 ) then
    do n=1,ncomp
    do k=start_k,iend(3)+1
      kk=(k-1)/2+1
    do j=start_j,iend(2)+1
      jj=(j-1)/2+1
    do i=0,ista(1)+1,-1
      ii=(i-1)/2
      dst(i,j,k,n)=src(ii,jj,kk,n)
    enddo
    enddo
    enddo
    enddo
  endif

  do n=1,ncomp
  do k=0,ista(3)+1,-1
    kk=(k-1)/2
  do j=0,ista(2)+1,-1
    jj=(j-1)/2
  do i=0,ista(1)+1,-1
    ii=(i-1)/2
    dst(i,j,k,n)=src(ii,jj,kk,n)
  enddo
  enddo
  enddo
  enddo

  do n=1,ncomp
  do k=start_k,iend(3)+1
    kk=(k-1)/2+1
  do j=start_j,iend(2)+1
    jj=(j-1)/2+1
  do i=start_i,iend(1)+1
    ii=(i-1)/2+1
    dst(i,j,k,n)=src(ii,jj,kk,n)
  enddo
  enddo
  enddo
  enddo

!  do n=1,ncomp
!  do k=ista(3)+1,iend(3)+1
!    kk=(k-1)/2+1
!  do j=ista(2)+1,iend(2)+1 
!    jj=(j-1)/2+1
!  do i=ista(1)+1,iend(1)+1
!    ii=(i-1)/2+1
!    dst(i,j,k,n)=src(ii,jj,kk,n)
!  enddo
!  enddo
!  enddo
!  enddo

  return

end subroutine cio_dfi_coarse_ijkn
!
!
!> ********************************************************************
!! @brief 粗い格子から密な格子への補間（ゼロ次）
!! @param [out] dst   密な格子系
!! @param [in]  dhead dstの開始インデックス（グローバル）
!! @param [in]  dtail dstの終了インデックス（グローバル）
!! @param [in]  dg    dstの仮想セル
!! @param [in]  src   粗い格子系
!! @param [in]  shead srcの開始インデックス（グローバル）
!! @param [in]  stail srcの終了インデックス（グローバル）
!! @param [in]  sg    srcの仮想セル
!! @param [in]  ncomp 成分数
!! @param [in]  ista  補間の開始インデクス
!! @param [in]  iend  補間の終了インデクス
!<

subroutine cio_dfi_coarse_nijk(dst,dhead,dtail,dg,src,shead,stail,sg,ncomp,ista,iend,irank)

  implicit none

  integer               :: irank

  integer, dimension(3) :: dhead,dtail,shead,stail,ista,iend
  integer               :: dg,sg,ncomp
  real,dimension(ncomp,dhead(1)-dg:dtail(1)+dg, dhead(2)-dg:dtail(2)+dg,dhead(3)-dg:dtail(3)+dg) :: dst
  real,dimension(ncomp,shead(1)-sg:stail(1)+sg, shead(2)-sg:stail(2)+sg,shead(3)-sg:stail(3)+sg) :: src

  integer               :: i,j,k,n,ii,jj,kk
  real                  :: q

  integer               :: start_i,start_j,start_k


  write(*,*) "rank : ",irank," **** cio_dfi_coarse_nijk "

  write(*,*) "rank : ",irank,"ista : ",ista
  write(*,*) "rank : ",irank,"iend : ",iend

  start_i=ista(1)+1
  start_j=ista(2)+1
  start_k=ista(3)+1
  if( start_i.le.0 ) start_i=1 
  if( start_j.le.0 ) start_j=1 
  if( start_k.le.0 ) start_k=1 

  if( ista(3).lt.0 ) then
    do k=0,ista(3)+1,-1
      kk=(k-1)/2
    do j=start_j,iend(2)+1
      jj=(j-1)/2+1
    do i=start_i,iend(1)+1
      ii=(i-1)/2+1
    do n=1,ncomp
      dst(n,i,j,k)=src(n,ii,jj,kk)
    enddo
    enddo
    enddo
    enddo
  endif

  if( ista(2).lt.0 ) then
    do k=start_k,iend(3)+1
      kk=(k-1)/2+1
    do j=0,ista(2)+1,-1
      jj=(j-1)/2
    do i=start_i,iend(1)+1
      ii=(i-1)/2+1
    do n=1,ncomp
      dst(n,i,j,k)=src(n,ii,jj,kk)
    enddo
    enddo
    enddo
    enddo
  endif

  if( ista(1).lt.0 ) then
    do k=start_k,iend(3)+1
      kk=(k-1)/2+1
    do j=start_j,iend(2)+1
      jj=(j-1)/2+1
    do i=0,ista(1)+1,-1
      ii=(i-1)/2
    do n=1,ncomp
      dst(n,i,j,k)=src(n,ii,jj,kk)
    enddo
    enddo
    enddo
    enddo
  endif

  do k=0,ista(3)+1,-1
    kk=(k-1)/2
  do j=0,ista(2)+1,-1
    jj=(j-1)/2
  do i=0,ista(1)+1,-1
    ii=(i-1)/2
  do n=1,ncomp
    dst(n,i,j,k)=src(n,ii,jj,kk)
  enddo
  enddo
  enddo
  enddo

  do k=start_k,iend(3)+1
    kk=(k-1)/2+1
  do j=start_j,iend(2)+1
    jj=(j-1)/2+1
  do i=start_i,iend(1)+1
    ii=(i-1)/2+1
  do n=1,ncomp
    dst(n,i,j,k)=src(n,ii,jj,kk)
  enddo
  enddo
  enddo
  enddo

  return

end subroutine cio_dfi_coarse_nijk
