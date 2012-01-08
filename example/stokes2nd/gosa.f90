!1234567891*********2*********3*********4*********5*********6*********72
!   exact.f90 2011.5.19-2011.11.17(stokes2nd)
!***********************************************************************

	real pi,sqrt2
	real u, y(0:59), exact(0:59,4)
	integer j, m, n

 	pi=acos(-1.0)
	sqrt2=sqrt(2.0)

	open(1, file='sample_x=0_y.log', status='old')
	open(10, file='stokes_exact.ascii', status='unknown')
	open(11, file='stokes2nd_gosa_t07.ascii', status='new')
	open(12, file='stokes2nd_gosa_t14.ascii', status='new')
	open(13, file='stokes2nd_gosa_t21.ascii', status='new')
	open(14, file='stokes2nd_gosa_t28.ascii', status='new')
	open(15, file='stokes2nd_gosa_t35.ascii', status='new')
	open(16, file='stokes2nd_gosa_t42.ascii', status='new')
	open(17, file='stokes2nd_gosa_t49.ascii', status='new')
	open(18, file='stokes2nd_gosa_t56.ascii', status='new')

	read(1,*)

	do n=0,59 
	  read(1,*) x, y(n)
	  do m=1,4
	    exact(n,m)=exp(-y(n)*7.16*0.5*sqrt2)*cos(0.5*m*pi-y(n)*7.16*0.5*sqrt2)
	  end do
	  write(10,*) y(n)*7.16, (exact(n,m), m=1,4)
	end do

	do 
	  read(1,*) j
	  if(j.eq.4500) then	!t=7, 0.5*pi
	    exit
	  end if
	end do
	do n=0,59
	  read(1,*) u
	  write(11,*) u-exact(n,1), y(n)*7.16
	end do
	close(11)

	do 
	  read(1,*) j
	  if(j.eq.9001) then	!t=14, 1.0*pi
	    exit
	  end if
	end do
 	do n=0,59
	  read(1,*) u
	  write(12,*) u-exact(n,2), y(n)*7.16
 	end do
	close(12)

	do 
	  read(1,*) j
	  if(j.eq.13502) then	!t=21, 1.5*pi
	    exit
	  end if
	end do
 	do n=0,59
	  read(1,*) u
	  write(13,*) u-exact(n,3), y(n)*7.16
 	end do
	close(13)

	do 
	  read(1,*) j
	  if(j.eq.18003) then	!t=28, 2.0*pi
	    exit
	  end if
	end do
 	do n=0,59
	  read(1,*) u
	  write(14,*) u-exact(n,4), y(n)*7.16
 	end do
	close(14)

	do 
	  read(1,*) j
	  if(j.eq.22502) then	!t=35	0.5*pi
	    exit
	  end if
	end do
 	do n=0,59
	  read(1,*) u
	  write(15,*) u-exact(n,1), y(n)*7.16
 	end do
	close(15)

	do 
	  read(1,*) j
	  if(j.eq.26998) then	!t=42	1.0*pi
	    exit
	  end if
	end do
 	do n=0,59
	  read(1,*) u
	  write(16,*) u-exact(n,2), y(n)*7.16
 	end do
	close(16)

	do 
	  read(1,*) j
	  if(j.eq.31493) then	!t=49	1.5*pi
	    exit
	  end if
	end do
 	do n=0,59
	  read(1,*) u
	  write(17,*) u-exact(n,3), y(n)*7.16
 	end do
	close(17)

	do 
	  read(1,*) j
	  if(j.eq.35989) then	!t=56	2.0*pi
	    exit
	  end if
	end do
 	do n=0,59
	  read(1,*) u
	  write(18,*) u-exact(n,4), y(n)*7.16
 	end do
	close(18)

	end
