!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!add.s CIO dbwrite
    subroutine v3dwrite(ID,v,sz,g)

    implicit none

    integer :: ID,g
    integer,dimension(3) :: sz
    real,dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g,3) :: v

    integer :: i,j,k,ix,jx,kx

    !write(*,*) "fort.",ID

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)

    do k=1,kx
    do j=1,jx
    do i=1,ix
      write(ID,100) i,j,k,v(i,j,k,1),v(i,j,k,2),v(i,j,k,3)
    enddo
    enddo
    enddo

 100 format(1h ,3i4,3e12.5)

    return
    end subroutine v3dwrite

    subroutine nv3dwrite(ID,v,sz,g)

    implicit none

    integer :: ID,g
    integer,dimension(3) :: sz
    real,dimension(3,1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) :: v

    integer :: i,j,k,ix,jx,kx

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

!!add.e
