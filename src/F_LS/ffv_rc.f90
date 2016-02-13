!###################################################################################
!
! FFV-C
! Frontflow / violet Cartesian
!
!
! Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
! All rights reserved.
!
! Copyright (c) 2011-2015 Institute of Industrial Science, The University of Tokyo.
! All rights reserved.
!
! Copyright (c) 2012-2016 Advanced Institute for Computational Science, RIKEN.
! All rights reserved.
!
!###################################################################################

!> @file   ffv_SOR.f90
!! @brief  SOR series routine
!! @author aics
!<


!> ********************************************************************
!! @brief point SOR法
!! @param [in,out] p    圧力
!! @param [in]     sz   配列長
!! @param [in]     g    ガイドセル長
!! @param [in]     omg  加速係数
!! @param [out]    cnv  収束判定値　修正量の自乗和と残差の自乗和、解ベクトルの自乗和
!! @param [in]     b    RHS vector
!! @param [in]     bp   BCindex P
!! @param [out]    flop flop count
!! @note Activeマスクの位置は，固体中のラプラス式を解くように，更新式にはかけず残差のみにする
!<
subroutine rc_sor (para_key, p, sz, g, omg, res, wrk, bnd, commmode, do_time, commtime, dtype, iparam, rparam)
implicit none
include 'ffv_f_params.h'

integer                                                   ::  para_key, g, dtype
integer                                                   ::  iparam(*)
integer                                                   ::  commmode ! 0:CommBndCell, 1:Commface(NoHide), 2:Commface(Hide)
integer                                                   ::  do_time, commtime
integer, dimension(3)                                     ::  sz
real                                                      ::  omg, res
real*8                                                    ::  rparam(*)
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  wrk, p, bnd

real*8,dimension(:),allocatable                           ::  scalprod, matr, alph, stolb
real,dimension(:,:,:),allocatable                         ::  xt, yt, rest, res_fct
real,dimension(:,:,:,:),allocatable                       ::  xrc, yrc
integer                                                   ::  i, j, k, ix, jx, kx, it, ierr
integer                                                   ::  lg, ig, jg, kg, tneum
integer                                                   ::  mrc, nrc, iter, i_iter, n_iter
integer, parameter                                        ::  nrc_max = 20, iter_max = 100
integer, parameter                                        ::  isneum = 0
integer                                                   ::  i_inner, n_inner, oki,step
real                                                      ::  e, ee
real                                                      ::  fc, res_abs, res_0
real*8                                                    ::  prev, TuneGetTime
real*8                                                    ::  al, t_eps, f_v

integer, dimension(2)                                     ::  err


if(iparam(8).lt.0) return
err(1) = 0
iparam(1) = 0
iparam(2) = 0
iparam(3) = 0
!    n_iter = min(iparam(8), iter_max)
!    nrc    = min(iparam(8), nrc_max)
oki = 5; step = 10
n_iter = iter_max
nrc    = nrc_max

if (dtype == 4) then
  e = 2.4e-7
  t_eps = 0.995*rparam(1)*rparam(1)
else
  e = 4.4e-16
  t_eps = 0.999*rparam(1)*rparam(1)
endif
ee = 1.0 + e

ix = sz(1)
jx = sz(2)
kx = sz(3)
lg = 1 - g
ig = ix + g
jg = jx + g
kg = kx + g
tneum = 0

if((isneum.ne.0).and.(g.gt.0)) then; tneum = 1; endif
! Distributed like "wrk", "p", "bnd"
allocate(xt(lg:ig,lg:jg,lg:kg),      yt(lg:ig,lg:jg,lg:kg), rest(lg:ig,lg:jg,lg:kg), res_fct(lg:ig,lg:jg,lg:kg))
allocate(xrc(ix,jx,kx,nrc), yrc(ix,jx,kx,nrc))
! Allocated only on main processor
allocate(scalprod(nrc*nrc), alph(nrc), stolb(nrc), matr(nrc*nrc))

iter = 1

!$OMP PARALLEL SHARED  (ix,jx,kx,p,bnd,wrk,omg,e) &
!              PRIVATE (pl,pc,pr,ytmp,k,j,i,it,fc)
f_v = 0.0
!$OMP BARRIER
!$OMP DO
do k=1,kx;
do j=1,jx;
do i=1,ix
fc = e + (bnd(i+1,j  ,k  ) + bnd(i-1,j  ,k  ) &
+ bnd(i  ,j+1,k  ) + bnd(i  ,j-1,k  ) &
+ bnd(i  ,j  ,k+1) + bnd(i  ,j  ,k-1))*bnd(i,j,k)
res_fct(i,j,k) = bnd(i,j,k)*ee/fc
al = res_fct(i,j,k)*wrk(i,j,k);
rest(i,j,k) = al
f_v = f_v + al*al
end do;
end do;
end do
!$OMP END DO

t_eps = f_v*t_eps
res_0 = f_v
if(f_v.lt.1.0e-30) return

call cbs3d_mvprod (p, sz, g, e, bnd, yt, res_fct)

call cbs3d_ressmpl (sz, g, res_abs, rest, yt)

iparam(3) = iparam(3) + 1
res_0 = sqrt(res_abs/f_v)
!    if(res_abs.lt.t_eps) goto 10


do i_iter=1,n_iter
!$OMP BARRIER
!$OMP DO
do k=lg,kg;
do j=lg,jg;
do i=lg,ig
  xt(i,j,k) = 0.0
end do;
end do;
end do
!$OMP END DO


if(    i_iter.le.  oki) then; n_inner =   step
elseif(i_iter.le.2*oki) then; n_inner = 2*step
elseif(i_iter.le.3*oki) then; n_inner = 3*step
elseif(i_iter.le.4*oki) then; n_inner = 4*step
else;                         n_inner = 5*step
endif

if(n_inner.lt.1) n_inner = 1
do i_inner = 1, n_inner

call cbs3d_psor_m1 (para_key, xt, sz, g, omg, e, res, rest, bnd, commmode, do_time, commtime, dtype)
iparam(2) = iparam(2) + 1
enddo

if(res.lt.t_eps) then
!$OMP DO
do k=lg,kg;
do j=lg,jg;
do i=lg,ig
  p(i,j,k) = p(i,j,k) + xt(i,j,k)
end do;
end do;
end do
!$OMP END DO

call cbs3d_psor (para_key, p, sz, g, omg, res, wrk, bnd, commmode, do_time, commtime, dtype)
iparam(2) = iparam(2) + 1

if(res.lt.t_eps) then

print*,'Final',i_iter,res,sngl(sqrt(res/f_v))
exit

else

do k=1,kx
do j=1,jx
do i=1,ix
rest(i,j,k) = res_fct(i,j,k)*wrk(i,j,k)
end do
end do
end do

call cbs3d_mvprod (p, sz, g, e, bnd, yt, res_fct)

call cbs3d_ressmpl (sz, g, res_abs, rest, yt)

cycle

endif
endif

call cbs3d_mvprod (xt, sz, g, e, bnd, yt, res_fct)

iparam(3) = iparam(3) + 1

call rescut_real_struct1(sz, g, xt, yt, rest, xrc, yrc, scalprod, alph, stolb, matr, iter, nrc, err)
iparam(1) = iparam(1) + 1

if(err(1).ne.0) return

if(tneum.eq.1) then
  call cbs3d_set_neum (xt, sz, g)
endif
!$OMP BARRIER

!$OMP DO
do k=lg,kg
do j=lg,jg
do i=lg,ig
  p(i,j,k) = p(i,j,k) + xt(i,j,k)
end do
end do
end do
!$OMP END DO

call cbs3d_ressmpl (sz, g, res_abs, rest, yt)

end do
!$OMP END PARALLEL
10  continue

if(allocated(matr))	deallocate(matr)
if(allocated(scalprod))	deallocate(scalprod)
if(allocated(alph))	deallocate(alph)
if(allocated(stolb))	deallocate(stolb)

if(allocated(xt))	deallocate(xt)
if(allocated(yt))	deallocate(yt)
if(allocated(rest))	deallocate(rest)
if(allocated(res_fct))	deallocate(res_fct)
if(allocated(xrc))	deallocate(xrc)
if(allocated(yrc))	deallocate(yrc)

return
end subroutine rc_sor


!   *******************************************************************************************
subroutine cbs3d_psor_m1 (para_key, p, sz, g, omg, e, res, wrk_m, bnd, cm_mode, do_time, cm_time, dtype)
implicit none
include 'ffv_f_params.h'
integer                                                   ::  para_key, g, dtype
integer, dimension(3)                                     ::  sz
real                                                      ::  omg, e, res
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  wrk_m, p, bnd
integer                                                   ::  cm_mode ! 0:CommBndCell, 1:CommBndCell2
integer                                                   ::  do_time
real*8                                                    ::  cm_time

integer                                                   ::  ierr, wait_num, req(12)
integer                                                   ::  i, j, k, ix, jx, kx, iret
real                                                      ::  cpd, s0, ee
real                                                      ::  f1, f2, f3, f4, f5, f6, f0
real                                                      ::  tmp, r1, dp
real*8                                                    ::  prev, TuneGetTime

ix = sz(1)
jx = sz(2)
kx = sz(3)
ee= 1.0+e
r1 = 0.0

!$OMP PARALLEL SHARED  (ix,jx,kx,p,bnd,wrk,omg,e,ee) &
!              PRIVATE (k,j,i,s0,ss,f0,f1,f2,f3,f4,f5,f6,cpd)
!$OMP BARRIER
!$OMP DO
do k=1,kx
do j=1,jx
do i=1,ix
f0=bnd(i,j,k)
f1=bnd(i+1,j  ,k  )
f2=bnd(i-1,j  ,k  )
f3=bnd(i  ,j+1,k  )
f4=bnd(i  ,j-1,k  )
f5=bnd(i  ,j  ,k+1)
f6=bnd(i  ,j  ,k-1)
cpd=ee/(e+(f1+f2+f3+f4+f5+f6)*f0)

s0=  (f1*p(i+1,j  ,k  )+f2*p(i-1,j  ,k  )  &
+ f3*p(i  ,j+1,k  )+f4*p(i  ,j-1,k  )  &
+ f5*p(i  ,j  ,k+1)+f6*p(i  ,j  ,k-1))*f0
dp = wrk_m(i,j,k) + (cpd*s0-p(i,j,k))*f0
p(i,j,k)=p(i,j,k) + omg*dp
r1 = r1 + dp*dp
end do
end do
end do
!$OMP END DO NOWAIT
!$OMP END PARALLEL

res = r1
call SklIsParallel(iret)
if ( iret .eq. 1 ) then
tmp = r1
if ( dtype == 4 ) then
call SklAllreduce(tmp, res, 1, SKL_REAL4, SKL_SUM, SKL_DEFAULT_GROUP, ierr)
else
call SklAllreduce(tmp, res, 1, SKL_REAL8, SKL_SUM, SKL_DEFAULT_GROUP, ierr)
endif
end if

if( do_time.eq.1 ) then
call SklBarrier(SKL_DEFAULT_GROUP, ierr)
prev = TuneGetTime()
endif

if( cm_mode.eq.0 ) then    ! $BF14|DL?.(B
call SklCommBndCell(p, 1, ierr)
else                        ! $BHsF14|DL?.(B
call SklCommBndCell2(p, 1, wait_num, req, ierr)
call SklWaitAll(wait_num, req, ierr);
endif

if( do_time.eq.1 ) then
cm_time = cm_time + (TuneGetTime()-prev)
endif

return
end subroutine cbs3d_psor_m1


!   *******************************************************************************************
subroutine cbs3d_mvprod (p, sz, g, e, bnd, bt, res_fct)
implicit none
include 'sklparaf.h'
integer                                                   ::  g
integer, dimension(3)                                     ::  sz
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  p, bnd, bt, res_fct
real                                                      ::  e

integer                                                   ::  i, j, k, ix, jx, kx
real                                                      ::  f0, f1, f2, f3, f4,f5,f6,fc,ee1

ix = sz(1)
jx = sz(2)
kx = sz(3)
ee1 = 1.0/(1.0 + e)
!$OMP PARALLEL SHARED  (ix,jx,kx,p,bnd,bt,e,ee1) &
!              PRIVATE (k,j,i,f1,f2,f3,f4,f5,f6,fc)
!$OMP BARRIER
!$OMP DO
do k=1,kx
do j=1,jx
do i=1,ix
f0=bnd(i  ,j  ,k  )
f1=bnd(i+1,j  ,k  )
f2=bnd(i-1,j  ,k  )
f3=bnd(i  ,j+1,k  )
f4=bnd(i  ,j-1,k  )
f5=bnd(i  ,j  ,k+1)
f6=bnd(i  ,j  ,k-1)
fc = ee1*(e+(f1+f2+f3+f4+f5+f6)*f0)
bt(i,j,k) =              (    fc*p(i  ,j  ,k  )  &
- f0*(f1*p(i+1,j  ,k  ) + f2*p(i-1,j  ,k  )  &
+ f3*p(i  ,j+1,k  ) + f4*p(i  ,j-1,k  )  &
+ f5*p(i  ,j  ,k+1) + f6*p(i  ,j  ,k-1)))*res_fct(i,j,k)
end do
end do
end do
!$OMP END DO NOWAIT
!$OMP END PARALLEL
return
end subroutine cbs3d_mvprod


!   *******************************************************************************************
subroutine cbs3d_ressmpl (sz, g, res, zansa, bt)
implicit none
include 'sklparaf.h'
integer                                                   ::  g
integer, dimension(3)                                     ::  sz
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  zansa, bt
real                                                      ::  res

integer                                                   ::  i, j, k, ix, jx, kx
integer                                                   ::  ierr
real                                                      ::  tmp_res
real                                                      ::  al

ix = sz(1)
jx = sz(2)
kx = sz(3)
res = 0.0

!$OMP PARALLEL SHARED  (ix,jx,kx,zansa,bt,res_fct) &
!              PRIVATE (k,j,i,fc,al)
!$OMP BARRIER
!$OMP DO
do k=1,kx
do j=1,jx
do i=1,ix
al = zansa(i,j,k) - bt(i,j,k); zansa(i,j,k) = al
res = res + al*al
end do
end do
end do
!$OMP END DO NOWAIT
!$OMP END PARALLEL

tmp_res = res
call SklAllreduce(tmp_res, res, 1, SKL_REAL4, SKL_SUM, SKL_DEFAULT_GROUP, ierr)
res = res
return
end subroutine cbs3d_ressmpl


!   *******************************************************************************************
subroutine rescut_real_struct1(sz, g, x, y, res, xrc, yrc, scalprod, alph, stolb, matr, iter, nrc, err)
implicit none
include 'sklparaf.h'
integer                                               ::  g
integer                                               ::  iter, nrc
integer, dimension(2)                                 ::  err
integer, dimension(3)                                 ::  sz
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  x, y, res
real, dimension(sz(1), sz(2), sz(3), nrc)             ::  xrc, yrc
real*8, dimension(nrc*nrc)                            ::  scalprod, matr
real*8, dimension(nrc)                                ::  alph, stolb

real*8                                                ::  al, tmp_al, bl, a1
integer                                               ::  i, j, k, ix, jx, kx
integer                                               ::  irc, mrc, l, chan, ierr

ix = sz(1)
jx = sz(2)
kx = sz(3)

if(iter.le.nrc) then; mrc = iter; chan = iter
else
mrc = nrc
chan = modulo(iter, nrc); if(chan.eq.0) chan = nrc
endif
!$OMP PARALLEL SHARED  (ix,jx,kx,x,y,res,xrc,yrc) PRIVATE (k,j,i,a1)
!	print*, iter, mrc, chan
al = 0.0
!$OMP BARRIER
!$OMP DO
do k=1,kx
do j=1,jx
do i=1,ix
al = al + y(i,j,k)*y(i,j,k)
end do
end do
end do
!$OMP END DO NOWAIT
!	  print*,'Alpha calculated',al,mrc
!	tmp_al = al
!	call SklAllreduce(tmp_al, al, 1, SKL_REAL4, SKL_SUM, SKL_DEFAULT_GROUP, ierr)

if(.not.((al.gt.(1.0d-300)).and.(al.lt.(1.0d300)))) then
print*,'ResCut error',al
err(1) = 98435; return
endif
al = 1.0/dsqrt(al)

!$OMP BARRIER
!$OMP DO
do k=1,kx
do j=1,jx
do i=1,ix
x(i,j,k) = al*x(i,j,k); xrc(i,j,k,chan) = x(i,j,k)
y(i,j,k) = al*y(i,j,k); yrc(i,j,k,chan) = y(i,j,k)
end do
end do
end do
!$OMP END DO NOWAIT

do irc = 1, mrc
al = 0.0; bl = 0.0
!$OMP BARRIER
!$OMP DO
do k=1,kx
do j=1,jx
do i=1,ix
a1 = yrc(i,j,k,irc)
al = al + y(i,j,k)*a1; bl = bl + res(i,j,k)*a1
end do
end do
end do
!$OMP END DO NOWAIT
!	  tmp_al = al
!	  call SklAllreduce(tmp_al, al, 1, SKL_REAL4, SKL_SUM, SKL_DEFAULT_GROUP, ierr)
!	  tmp_al = bl
!	  call SklAllreduce(tmp_al, bl, 1, SKL_REAL4, SKL_SUM, SKL_DEFAULT_GROUP, ierr)

scalprod((irc-1)*nrc+chan) = al; scalprod((chan-1)*nrc+irc) = al
stolb(irc) = bl; alph(irc) = bl
!	  if(i.eq.mrc) print*, real(al), real(bl/al)
enddo

!	call print_matrix(scalprod, mrc, nrc)

k = 1
do i = 1, mrc
l = (i-1)*nrc + 1
do j = 1, mrc; matr(k) = scalprod(l); k = k+1; l = l+1; enddo
enddo
call decomp_cholu_r8(matr, mrc, err)
if(err(1).ne.0) return
!	call dpptrf('L', mrc, matr, chan)
!	if(chan.ne.0) then; err(1) = abs(chan); goto 10; endif
call solve_cholu_r8(matr, alph, mrc)
!	al = 0.0
!	do i = 1, mrc
!	  l = (i-1)*nrc + 1; bl = stolb(i)
!	  do j = 1, mrc; bl = bl - alph(j)*scalprod(l); l = l+1; enddo
!	  al = al + bl*bl
!	enddo
!	print*,'ResCut zansa',dsqrt(al)
!$OMP BARRIER
!$OMP DO
do k=1,kx
do j=1,jx
do i=1,ix
x(i,j,k) = 0.0
y(i,j,k) = 0.0
end do
end do
end do
!$OMP END DO NOWAIT

do irc = 1, mrc
al = alph(irc)
!	print*,'Alpha ResCut',irc,al,stolb(irc)
!$OMP DO
do k=1,kx
do j=1,jx
do i=1,ix
x(i,j,k) = x(i,j,k) + al*xrc(i,j,k,irc)
y(i,j,k) = y(i,j,k) + al*yrc(i,j,k,irc)
end do
end do
end do
!$OMP END DO NOWAIT
enddo
!$OMP END PARALLEL
10	iter = iter + 1
end subroutine rescut_real_struct1

!************************************************************************
subroutine cbs3d_psor (para_key, p, sz, g, omg, res, wrk, bnd, commmode, do_time, commtime, ec)
implicit none
include '../../../include/sklparaf.h'
integer                                                   ::  i, j, k, ix, jx, kx, g, ec
integer, dimension(3)                                     ::  sz
real                                                      ::  f1, f2, f3, f4, f5, f6, f0, ds, dp, r, r1, r2, rs, b
real, dimension(2)                                        ::  tmp, res
real                                                      ::  omg, cpd, s0, e, ee
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  wrk, p, bnd
integer                                                   ::  ierr, wait_num, req(12)
integer                                                   ::  commmode ! 0:CommBndCell, 1:CommBndCell2
integer                                                   ::  do_time
real*8                                                    ::  commtime, prev, TuneGetTime
integer                                                   ::  para_key

ix = sz(1)
jx = sz(2)
kx = sz(3)
e = 1.0e-6
ee= 1.0+e
r1 = 0.0 ! Residual
r2 = 0.0 ! delta P

!$OMP PARALLEL SHARED  (ix,jx,kx,p,bnd,wrk,omg,e,ee) &
!              PRIVATE (k,j,i,s0,ss,f0,f1,f2,f3,f4,f5,f6,cpd)
!$OMP BARRIER
!$OMP DO
do k=1,kx
do j=1,jx
do i=1,ix
f0=bnd(i,j,k)
f1=bnd(i+1,j  ,k  )*f0
f2=bnd(i-1,j  ,k  )*f0
f3=bnd(i  ,j+1,k  )*f0
f4=bnd(i  ,j-1,k  )*f0
f5=bnd(i  ,j  ,k+1)*f0
f6=bnd(i  ,j  ,k-1)*f0
b = wrk(i,j,k)
cpd=ee/(e+f1+f2+f3+f4+f5+f6)
s0 = f1*p(i+1,j  ,k  )+f2*p(i-1,j  ,k  )   &
+ f3*p(i  ,j+1,k  )+f4*p(i  ,j-1,k  )   &
+ f5*p(i  ,j  ,k+1)+f6*p(i  ,j  ,k-1) + b
dp = (cpd*s0-p(i,j,k))*f0
ds = omg*dp
p(i,j,k)=p(i,j,k) + ds
r2 = r2 + ds*ds
r  = (1.0/cpd-omg*(f2*p(i-1,j  ,k  )+f4*p(i  ,j-1,k  )+f6*p(i  ,j  ,k-1)))/omg * dp
rs = abs( r ) / (abs(b)+e)
r1 = r1 + rs*rs
end do
end do
end do
!$OMP END DO NOWAIT
!$OMP END PARALLEL

tmp(1) = r1
tmp(2) = r2
res(1) = r1
res(2) = r2
call SklAllreduce(tmp, res, 2, SKL_REAL4, SKL_SUM, SKL_DEFAULT_GROUP, ierr)
res(1) = sqrt(res(1))/real(ec)
res(2) = sqrt(res(2))/real(ec)

if( do_time.eq.1 ) then
call SklBarrier(SKL_DEFAULT_GROUP,ierr)
prev = TuneGetTime()
endif
if( commmode.eq.0 ) then
call SklCommBndCell(p, 1, ierr)
else
call SklCommBndCell2(p, 1, wait_num, req, ierr)
call SklWaitAll(wait_num, req, ierr);
endif
if( do_time.eq.1 ) then
commtime = commtime + (TuneGetTime()-prev)
endif

return
end subroutine cbs3d_psor


!> ********************************************************************
!! @brief AX
!! @param [out] b       AX=bのRHS
!! @param [in]  pos_rhs Poissonの右辺
!! @param [in]  bp      BCindexP
!! @param [in]  sz      配列長
!! @param [in]  g       ガイドセル
!! @param [out] b_l2    L2 norm
!! @param [out] flop    flop count
!<
subroutine rc_calc_b(b, pos_rhs, bp, sz, g, b_l2, flop)
implicit none
include 'ffv_f_params.h'
integer                                                   ::  i, j, k, ix, jx, kx, g, idx
integer, dimension(3)                                     ::  sz
real                                                      ::  dd, ss, d0, d1, d2
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  b, pos_rhs
integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp
double precision                                          ::  flop, b_l2

ix = sz(1)
jx = sz(2)
kx = sz(3)

b_l2 = 0.0d0;

flop = flop + dble(ix)*dble(jx)*dble(kx)*18.0d0

!$OMP PARALLEL &
!$OMP REDUCTION(+:b_l2) &
!$OMP PRIVATE(dd, ss, idx, d0, d1, d2) &
!$OMP FIRSTPRIVATE(ix, jx, kx)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do j=1,jx
do i=1,ix
idx = bp(i,j,k)
d0 = real(ibits(idx, bc_diag + 0, 1))
d1 = real(ibits(idx, bc_diag + 1, 1))
d2 = real(ibits(idx, bc_diag + 2, 1))
dd = d2*4.0 + d1*2.0 + d0  ! diagonal
!dd = real(ibits(idx, bc_diag, 3))  ! diagonal

ss = pos_rhs(i, j, k) / dd * real(ibits(idx, Active, 1))
bb(i, j, k) = ss
b_l2 = b_l2 + dble(ss) * dble(ss)
end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL

return
end subroutine rc_calc_b
