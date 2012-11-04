!********************************************************************
!
!   FFV : Frontflow / violet
!
!   Copyright (c) 2012 All right reserved.
!
!   Institute of Industrial Science, University of Tokyo, Japan. 
!
!********************************************************************

!> @file   PLT3D.f90
!! @brief  FlowBase PLOT3D related code
!! @author kero
!<

!  ***************************************************
!> @subroutine open_plot3d_file(iflag,fname,ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine open_plot3d_file(iflag,fname,ifl,fnsize,ierror)
    implicit none
    integer :: iflag,ifl,fnsize,ierror
    character(fnsize) :: fname
    
    if(iflag.eq.1) then !UNFORMATTED
       open(unit=ifl, err=10, file=fname, status='new', access='sequential', form='unformatted')
    else if(iflag.eq.2) then !FORMATTED
      open(unit=ifl, err=10, file=fname, status='new', access='sequential', form='formatted')
    end if
    ierror=1 !true
    return
    
 10 ierror=0 !false
    return
    
    end subroutine open_plot3d_file

!  ***************************************************
!> @subroutine open_plot3d_outputfile(iflag,fname,ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine open_plot3d_outputfile(iflag,fname,ifl,fnsize,ierror)
    implicit none
    integer :: iflag,ifl,fnsize,ierror
    character(fnsize) :: fname
    
    if(iflag.eq.1) then !UNFORMATTED
       open(unit=ifl, err=10, file=fname, status='unknown', access='sequential', form='unformatted')
    else if(iflag.eq.2) then !FORMATTED
      open(unit=ifl, err=10, file=fname, status='unknown', access='sequential', form='formatted')
    end if
    ierror=1 !true
    return
    
 10 ierror=0 !false
    return
    
    end subroutine open_plot3d_outputfile
 
!  ***************************************************
!> @subroutine open_plot3d_inputfile(iflag,fname,ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine open_plot3d_inputfile(iflag,fname,ifl,fnsize,ierror)
    implicit none
    integer :: iflag,ifl,fnsize,ierror
    character(fnsize) :: fname
    
    if(iflag.eq.1) then !UNFORMATTED
       open(unit=ifl, err=10, file=fname, status='old', access='sequential', form='unformatted')
    else if(iflag.eq.2) then !FORMATTED
      open(unit=ifl, err=10, file=fname, status='old', access='sequential', form='formatted')
    end if
    ierror=1 !true
    return
    
 10 ierror=0 !false
    return
    
    end subroutine open_plot3d_inputfile
 
 !  ***************************************************
!> @subroutine close_plot3d_file(ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine close_plot3d_file(ifl)
    implicit none
    integer :: ifl
    
    close(ifl)
    
    return
    end subroutine close_plot3d_file

!  ***************************************************
!> @subroutine write_line(buff,ifl,fnsize)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine write_line(buff,ifl,fnsize)
    implicit none
    integer :: ifl,fnsize
    character(fnsize) :: buff

    write(ifl,'(a)') buff
    
    end subroutine write_line
 
 !  ***************************************************
!> @subroutine read_line(buff,ifl,fnsize,iend)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine read_line(buff,ifl,fnsize,iend)
    implicit none
    integer :: ifl,fnsize,iend
    character(fnsize) :: buff

    iend=0
    read(ifl,'(a)',end=10) buff
    return
    
 10 continue
    iend=1
    return   
    
    end subroutine read_line

!  ***************************************************
!> @subroutine write_fvbnd_boundary(tp,gn,imin,imsx,jmin,jmax,kmin,kmax,flag,dir,ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine write_fvbnd_boundary(tp,gn,imin,imsx,jmin,jmax,kmin,kmax,flag,dir,ifl)
    implicit none
    integer :: tp,gn,imin,imsx,jmin,jmax,kmin,kmax,dir,ifl
    character(1) :: flag

    write(ifl,'(8i8,1x,1a,1x,i8)') tp,gn,imin,imsx,jmin,jmax,kmin,kmax,flag,dir
    
    end subroutine write_fvbnd_boundary

 !  ***************************************************
!> @subroutine read_fvbnd_boundary(tp,gn,imin,imsx,jmin,jmax,kmin,kmax,flag,dir,ifl)
!! @brief 
!! @param 
!! @param[out] 
!<
    subroutine read_fvbnd_boundary(tp,gn,imin,imsx,jmin,jmax,kmin,kmax,flag,dir,ifl)
    implicit none
    integer :: tp,gn,imin,imsx,jmin,jmax,kmin,kmax,dir,ifl
    character(1) :: flag

    read(ifl,'(8i8,1x,1a,1x,i8)') tp,gn,imin,imsx,jmin,jmax,kmin,kmax,flag,dir
    
    end subroutine read_fvbnd_boundary
