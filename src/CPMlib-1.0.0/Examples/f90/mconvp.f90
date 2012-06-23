!   *********************************************************
!   ***      Sample code of porting
!   ***      2010 keno@VCAD, RIKEN
!
!   *********************************************************
    module array_def
    real, dimension(:,:,:), allocatable, save :: p, q, w
    end module array_def

!   *********************************************************
    program diffusion

    use array_def
    implicit none
    integer       :: nx, ny, nz
    integer       :: n, scheme, itrMax, stp, stpMax
    real          :: dx, dy, dz
    real          :: eps, er, omg, dt, cf

!   Initialize
    write (*,*) 'Dimension size'
    write (*,*) 'Input dimension x size ='
    read  (*,*)  nx
    write (*,*) 'Input dimension y size ='
    read  (*,*)  ny
    write (*,*) 'Input dimension z size ='
    read  (*,*)  nz

    write (*,*) 'SELECT time integration method'
    write (*,*) '        1 --- Euler Explicit'
    write (*,*) '        2 --- Euler Implicit'
    read  (*,*) scheme
    
    write (*,*) 'Time step max='
    read  (*,*) stpMax
    
    write (*,*) 'Time increment='
    read  (*,*) dt
    
    write (*,*) 'Coef. of diffusion='
    read  (*,*) cf
      
    if ( scheme == 2 ) then
      write (*,*) 'Iteration max='
      read  (*,*) itrMax
      write (*,*) 'Epsilon ='
      read  (*,*) eps
      write (*,*) 'Relazation coef. ='
      read  (*,*) omg
    end if

    dx = 1.0/real(nx)
    dy = 1.0/real(ny)
    dz = 1.0/real(nz)

    write (*,*) 'nx, ny, nz = ', nx, ny, nz 
    write (*,*) 'dx, dy, dz = ', dx, dy, dz 
    write (*,*) 'dt         = ', dt
    if ( scheme == 1 ) then
      write(*,*) 'Euler Explicit ftcs'
    else
      write(*,*) 'Euler Implicit jacobi'
    end if
    
    call allocate_array (nx, ny, nz)
    call pbc (nx, ny, nz, p, dx, dy, dz)
    
!   time step loop
    do stp=1, stpMax
      if ( scheme == 1 ) then
        call ftcs (nx, ny, nz, p, q, cf, dx, dy, dz, dt, er)
        call pbc (nx, ny, nz, p, dx, dy, dz)
        write (*,*) stp, er
      else
        er = 0.0
        w(0:nx+1, 0:ny+1, 0:nz+1) = p(0:nx+1, 0:ny+1, 0:nz+1)
        ! iteration
        do n=1, itrMax
          call jacobi (nx, ny, nz, p, q, w, cf, dx, dy, dz, dt, omg, er)
          call pbc (nx, ny, nz, p, dx, dy, dz)
          if (er<eps) exit
        end do
        write (*,*) stp, n, er
      end if
    end do

!   Post
    call fileout (nx, ny, nz, p, w, dx, dy, dz)


    stop
    end program diffusion

!   ******************************
    subroutine allocate_array (nx, ny, nz)
    use array_def
    implicit none
    integer :: status
    integer :: nx, ny, nz
    
    allocate ( p(0:nx+1, 0:ny+1, 0:nz+1) , stat=status )
    allocate ( w(0:nx+1, 0:ny+1, 0:nz+1) , stat=status )
    allocate ( q(0:nx+1, 0:ny+1, 0:nz+1) , stat=status )
    p(0:nx+1, 0:ny+1, 0:nz+1) = 0.d0

    return
    end subroutine allocate_array

!   ****************************************************
    subroutine jacobi (nx, ny, nz, p, q, w, cf, dx, dy, dz, dt, omg, er)
    implicit none
    integer                                 :: nx, ny, nz 
    real, dimension (0:nx+1,0:ny+1,0:nz+1)  :: p, q, w
    real                                    :: cf, dx, dy, dz, dt, omg, er

    integer                                 :: i, j, k
    real                                    :: ddx, ddy, ddz, cpd
    real                                    :: s0, ss, sx, sy, sz
    
    er = 0.0
    ddx = cf*dt/(dx*dx)
    ddy = cf*dt/(dy*dy)
    ddz = cf*dt/(dz*dz)
    cpd=1.0/(1.0+2.0*ddx+2.0*ddy+2.0*ddz)
    
    do k=1,nz
    do j=1,ny
    do i=1,nx
      sx= p(i+1,j  ,k  ) + p(i-1,j  ,k  )
      sy= p(i  ,j+1,k  ) + p(i  ,j-1,k  )
      sz= p(i  ,j  ,k+1) + p(i  ,j  ,k-1)
      s0=ddx*sx+ddy*sy+ddz*sz
      ss=(cpd*(s0+w(i,j,k))-p(i,j,k))*omg
      q(i,j,k)=p(i,j,k)+ss
      er = er + ss*ss
    end do
    end do
    end do

    p(1:nx, 1:ny, 1:nz) = q(1:nx, 1:ny, 1:nz)
    er = sqrt(er)
    
    return
    end subroutine jacobi
    
!   ***************************************
    subroutine ftcs (nx, ny, nz, p, q, cf, dx, dy, dz, dt, er)
    implicit none

    integer                                 :: nx, ny, nz
    real, dimension (0:nx+1,0:ny+1,0:nz+1)  :: p, q
    real                                    :: cf, dx, dy, dz, dt, er

    integer                                 :: i, j, k
    real                                    :: ddx, ddy, ddz
    real                                    :: ss, sx, sy, sz
    
    er = 0.0
    ddx=cf*dt/(dx*dx)
    ddy=cf*dt/(dy*dy)
    ddz=cf*dt/(dz*dz)
    
    do k=1,nz
    do j=1,ny
    do i=1,nx
      sx= p(i+1,j  ,k  ) + p(i-1,j  ,k  ) &
        - p(i  ,j  ,k  )*2.0
      sy= p(i  ,j+1,k  ) + p(i  ,j-1,k  ) &
        - p(i  ,j  ,k  )*2.0
      sz= p(i  ,j  ,k+1) + p(i  ,j  ,k-1) &
        - p(i  ,j  ,k  )*2.0
      ss=ddx*sx+ddy*sy+ddz*sz
      q(i,j,k)=p(i,j,k)+ss
      er = er + ss*ss
    end do
    end do
    end do

    p(1:nx, 1:ny, 1:nz) = q(1:nx, 1:ny, 1:nz)

    er = sqrt(er)
    
    return
    end subroutine ftcs
    
!   **************************
    subroutine pbc (nx, ny, nz, p, dx, dy, dz)
    implicit none
    integer                               :: nx, ny, nz
    real, dimension(0:nx+1,0:ny+1,0:nz+1) :: p
    real                                  :: dx, dy, dz

    integer                               :: i, j, k
    real                                  :: pi, x, y

    pi = 2.0*asin(1.0)
!kero left lower coner of cell(1,1,1) is defined as an original point O(0.0, 0.0, 0.0)

    do j=1,ny
    do i=1,nx
      x= dx*(i-1)+dx*0.5
      y= dy*(j-1)+dy*0.5
      p(i,j,0)   =2.0*sin(pi*x)*sin(pi*y)-p(i,j,1)
      p(i,j,nz+1)=2.0*sin(pi*x)*sin(pi*y)-p(i,j,nz)
    end do
    end do
    
    do k=1,nz
    do j=1,ny
      p(0,   j,k) = p(1, j,k)
      p(nx+1,j,k) = p(nx,j,k)
    end do
    end do

    do k=1,nz
    do i=1,nx
      p(i,0,   k) = p(i,1, k)
      p(i,ny+1,k) = p(i,ny,k)
    end do
    end do
    
    return
    end subroutine pbc

!   ******************************
    subroutine fileout (nx, ny, nz, p, w, dx, dy, dz)
    implicit none
    integer                                 :: nx, ny, nz
    real  , dimension(0:nx+1,0:ny+1,0:nz+1) :: p
    real*4, dimension(0:nx+1,0:ny+1,0:nz+1) :: w
    real                                    :: dx, dy, dz

    integer                                 :: i,j,k
    integer                                 :: step
    real*4                                  :: org(3), pch(3), time
    character*256                           :: ename

!kero left lower coner of cell(1,1,1) is defined as an original point O(0.0, 0.0, 0.0)
!kero the original point is modified according to the above defined point

    ename='e.sph'

    step = 0
    time = 0.0
    pch(1) = real(dx)
    pch(2) = real(dy)
    pch(3) = real(dz)
    org    = -pch
    do k=0,nz
    do j=0,ny
    do i=0,nx
      w(i,j,k) = real(p(i,j,k))
    enddo
    enddo
    enddo

    open (unit=22,file=ename,form='unformatted')
    write (22) 1, 1
    write (22) nx+2,ny+2,nz+2
    write (22) org
    write (22) pch
    write (22) step, time
    write (22) w
    close (unit=22)

    return
    end subroutine fileout

