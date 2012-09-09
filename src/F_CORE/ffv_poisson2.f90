!********************************************************************
!
!   FFV : Frontflow / violet Cartesian
!
!   Copyright (c) 2012 All right reserved.
!
!   Institute of Industrial Science, The University of Tokyo, Japan.
!
!********************************************************************

!> @file   ffv_poisson.f90
!! @brief  Poisson routine
!! @author kero
!<

!> ********************************************************************
!! @brief 残差の計算
!! @param [out] bt   残差
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル長
!! @param [in]  p    圧力
!! @param [in]  bp   BCindex P
!! @param [out] flop flop count
!<
  subroutine mv_prod (bt, sz, g, p, bp, flop)
  implicit none
  include 'ffv_f_params.h'
  integer                                                   ::  i, j, k, ix, jx, kx, g, idx
  integer, dimension(3)                                     ::  sz
  double precision                                          ::  flop
  real                                                      ::  ndag_e, ndag_w, ndag_n, ndag_s, ndag_t, ndag_b
  real                                                      ::  dd, ss
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  p, bt
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)

  flop = flop + dble(ix)*dble(jx)*dble(kx)*15.0d0

!$OMP PARALLEL &
!$OMP PRIVATE(ndag_w, ndag_e, ndag_s, ndag_n, ndag_b, ndag_t, dd, ss, idx) &
!$OMP FIRSTPRIVATE(ix, jx, kx)

#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1) &
#elif defined _STATIC
!$OMP DO SCHEDULE(static) &
#else
!$OMP DO SCHEDULE(hoge)
#endif
!$OMP REDUCTION(+:res)

  do k=1,kx
  do j=1,jx
  do i=1,ix
    idx = bp(i,j,k)
    ndag_e = real(ibits(idx, bc_ndag_E, 1))  ! e, non-diagonal
    ndag_w = real(ibits(idx, bc_ndag_W, 1))  ! w
    ndag_n = real(ibits(idx, bc_ndag_N, 1))  ! n
    ndag_s = real(ibits(idx, bc_ndag_S, 1))  ! s
    ndag_t = real(ibits(idx, bc_ndag_T, 1))  ! t
    ndag_b = real(ibits(idx, bc_ndag_B, 1))  ! b

    dd = 1.0 / real(ibits(idx, bc_diag, 3))  ! diagonal

    ss =  ndag_e * p(i+1,j  ,k  ) &
        + ndag_w * p(i-1,j  ,k  ) &
        + ndag_n * p(i  ,j+1,k  ) &
        + ndag_s * p(i  ,j-1,k  ) &
        + ndag_t * p(i  ,j  ,k+1) &
        + ndag_b * p(i  ,j  ,k-1)

    bt(i,j,k) = ( p(i,j,k) - dd*ss ) * real(ibits(idx, Active, 1))
  end do
  end do
  end do
!$OMP END DO
!$OMP END PARALLEL

  return
  end subroutine mv_prod



!> ********************************************************************
!! @brief 残差の計算
!! @param [out] bt   残差
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル長
!! @param [in]  p    圧力
!! @param [in]  bp   BCindex P
!! @param [out] flop flop count
!<
  subroutine res_smpl (res, sz, g, zansa, bt, flop)
  implicit none
  include 'ffv_f_params.h'
  integer                                                   ::  i, j, k, ix, jx, kx, g
  integer, dimension(3)                                     ::  sz
  double precision                                          ::  flop
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  zansa, bt
  real                                                      ::  res, al

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
  res = 0.0

!$OMP PARALLEL &
!$OMP PRIVATE(al) &
!$OMP FIRSTPRIVATE(ix, jx, kx)

#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1) &
#elif defined _STATIC
!$OMP DO SCHEDULE(static) &
#else
!$OMP DO SCHEDULE(hoge)
#endif
!$OMP REDUCTION(+:res)

  do k=1,kx
  do j=1,jx
  do i=1,ix
    al = zansa(i,j,k) - bt(i,j,k)
    zansa(i,j,k) = al
    res = res + al*al
  end do
  end do
  end do
!$OMP END DO
!$OMP END PARALLEL

  return
  end subroutine res_smpl


!> ********************************************************************
!! @brief 乗算関数
!! @param [out] dst  目的の積
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル長
!! @param [in]  nc   サイズ
!! @param [in]  s    スカラー値
!! @param [in]  src  被積分配列
!! @param [out] flop  flop count
!<
  subroutine multiply (dst, sz, g, s, src, flop)
  implicit none
  include 'ffv_f_params.h'
  integer                                                     ::  i, j, k, ix, jx, kx, g
  integer, dimension(3)                                       ::  sz
  double precision                                            ::  flop
  real                                                        ::  s
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 15)  ::  dst
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  src

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, s)

#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1) &
#elif defined _STATIC
!$OMP DO SCHEDULE(static) &
#else
!$OMP DO SCHEDULE(hoge)
#endif
!$OMP REDUCTION(+:res)
  do k=1,kx
  do j=1,jx
  do i=1,ix
    dst(i, j, k, 1) = s * src(i,j,k)
  end do
  end do
  end do
!$OMP END DO
!$OMP END PARALLEL

  return
  end subroutine multiply




subroutine rescut_real_struct1(sz, g, x, y, res, xrc, yrc, &
scalprod, alph, stolb, matr, iter, nrc, err)
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

if(iter.le.nrc) then
mrc = iter; chan = iter
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
err(1) = 98435;
return
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
al = 0.0
bl = 0.0

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

scalprod((irc-1)*nrc+chan) = al
scalprod((chan-1)*nrc+irc) = al
stolb(irc) = bl; alph(irc) = bl
!	  if(i.eq.mrc) print*, real(al), real(bl/al)
enddo

!	call print_matrix(scalprod, mrc, nrc)

k = 1
do i = 1, mrc
l = (i-1)*nrc + 1
do j = 1, mrc;
matr(k) = scalprod(l);
k = k+1;
l = l+1;
enddo
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


!   ******************************************************
subroutine set_neum (t, sz, g)
implicit none
integer                                                   ::  g
integer, dimension(3)                                     ::  sz
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  t
integer                                                   ::  i, j, k, ix, jx, kx
integer                                                   ::  ix21, jx21, kx21, lg

ix = sz(1); ix21 = 2*ix+1
jx = sz(2); jx21 = 2*jx+1
kx = sz(3); kx21 = 2*kx+1
lg = 1 - g
! X-
do k=1,kx
do j=1,jx
do i=lg, 0
t(i,j,k) = t(1-i,j,k)
end do
end do
end do

! X+
do k=1,kx
do j=1,jx
do i=ix+1, ix+g
t(i,j,k) = t(ix21-i,j,k)
end do
end do
end do

!  Y-
do k=1,kx
do j=lg, 0
do i=1,ix
t(i,j,k) = t(i,1-j,k)
end do
end do
end do

!  Y+
do k=1,kx
do j=jx+1, jx+g
do i=1,ix
t(i,j,k) = t(i,jx21-j,k)
end do
end do
end do

!  Z-
do k=lg, 0
do j=1,jx
do i=1,ix
t(i,j,k) = t(i,j,1-k)
end do
end do
end do

!  Z+
do k=kx+1, kx+g
do j=1,jx
do i=1,ix
t(i,j,k) = t(i,j,kx21-k)
end do
end do
end do

return
end subroutine set_neum

!   ******************************************************
subroutine set_neum8 (t, sz, g)
implicit none
integer                                                   ::  g
integer, dimension(3)                                     ::  sz
real*8, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  ::  t
integer                                                   ::  i, j, k, ix, jx, kx
integer                                                   ::  ix21, jx21, kx21, lg

ix = sz(1); ix21 = 2*ix+1
jx = sz(2); jx21 = 2*jx+1
kx = sz(3); kx21 = 2*kx+1
lg = 1 - g
! X-
do k=1,kx
do j=1,jx
do i=lg, 0
t(i,j,k) = t(1-i,j,k)
end do
end do
end do

! X+
do k=1,kx
do j=1,jx
do i=ix+1, ix+g
t(i,j,k) = t(ix21-i,j,k)
end do
end do
end do

!  Y-
do k=1,kx
do j=lg, 0
do i=1,ix
t(i,j,k) = t(i,1-j,k)
end do
end do
end do

!  Y+
do k=1,kx
do j=jx+1, jx+g
do i=1,ix
t(i,j,k) = t(i,jx21-j,k)
end do
end do
end do

!  Z-
do k=lg, 0
do j=1,jx
do i=1,ix
t(i,j,k) = t(i,j,1-k)
end do
end do
end do

!  Z+
do k=kx+1, kx+g
do j=1,jx
do i=1,ix
t(i,j,k) = t(i,j,kx21-k)
end do
end do
end do

return
end subroutine set_neum8


!> ********************************************************************
!! @brief Residula-Cut SOR
!! @param[in,out] p 圧力
!! @param sz 配列長
!! @param g ガイドセル長
!! @param omg 加速係数
!! @param[out] res 絶対残差と相対残差
!! @param src0 固定ソース項
!! @param src1 反復毎に変化するソース項
!! @param bp BCindex P
!! @param[out] flop
!<
subroutine psor_m1 (para_key, p, sz, g, omg, e, res, wrk_m, bnd, cm_mode, do_time, cm_time, dtype)
implicit none
include 'ffv_f_params.h'
integer                                                   ::  para_key, g, dtype
integer, dimension(3)                                     ::  sz
real                                                      ::  omg, e, res
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  wrk_m, p, bnd
integer                                                   ::  cm_mode ! 0:CommBndCell, 1:CommBndCell2
integer                                                   ::  do_time
double precision                                          ::  cm_time

integer                                                   ::  ierr, wait_num, req(12)
integer                                                   ::  i, j, k, ix, jx, kx, iret
real                                                      ::  cpd, s0, ee
real                                                      ::  f1, f2, f3, f4, f5, f6, f0
real                                                      ::  tmp, r1, dp
double precision                                          ::  prev, TuneGetTime

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
end subroutine psor_m1

!   ***************************************************
subroutine cbs3d_gmressor (para_key, p, sz, g, omg, res, wrk, bnd, commmode, do_time, commtime, dtype, iparam, rparam)
implicit none
include 'sklparaf.h'

integer                                                   ::  para_key, g, dtype
integer                                                   ::  iparam(*)
integer                                                   ::  commmode ! 0:CommBndCell, 1:Commface(NoHide), 2:Commface(Hide)
integer                                                   ::  do_time, commtime
integer, dimension(3)                                     ::  sz
real                                                      ::  omg, res
real*8                                                    ::  rparam(*)
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  wrk, p, bnd

real,dimension(:,:,:),allocatable                         ::  xt, yt, rest, res_fct
real,dimension(:,:,:,:),allocatable                       ::  vm, zm
integer                                                   ::  i, j, k, ix, jx, kx, it, m
integer                                                   ::  ig, jg, kg, lg, im, jm, km, lm
integer                                                   ::  i_iter, n_iter, tneum
integer                                                   ::  iter, ierr, nrm, isfin
integer, parameter                                        ::  m_max = 15, iter_max = 100
integer, parameter                                        ::  isneum = 0
integer                                                   ::  i_inner, n_inner, oki, step
real                                                      ::  e, ee, fc, eps_abs, res_abs
real*8                                                    ::  rgm(m_max, m_max)
real*8                                                    ::  hgm(m_max+1, m_max)
real*8                                                    ::  cgm(m_max), sgm(m_max)
real*8                                                    ::  bgm(m_max), ygm(m_max)
real*8                                                    ::  prev, TuneGetTime
real*8                                                    ::  al, t_eps, f_v, beta, beta_1
real*8, parameter                                         ::  fct1 = 1.0d-4, fct2 = 6.0d0
integer, dimension(2)                                     ::  err
real*4, dimension(5)                                      ::  r4

if(iparam(8).lt.0) return
err(1) = 0
isfin = 0
!    n_iter = min(iparam(8), iter_max)
!    nrc    = min(iparam(8), nrc_max)
oki = 4; step = 6
m       = m_max
n_iter  = iter_max

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
if((isneum.ne.0).and.(g.gt.0)) then
tneum = 1
endif

! Distributed like "wrk", "p", "bnd"
allocate(xt(lg:ig,lg:jg,lg:kg),      yt(lg:ig,lg:jg,lg:kg), rest(lg:ig,lg:jg,lg:kg), res_fct(lg:ig,lg:jg,lg:kg))
allocate(vm(ix,jx,kx,m_max), zm(lg:ig,lg:jg,lg:kg,m_max))

!$OMP PARALLEL SHARED  (ix,jx,kx,p,bnd,wrk,omg,e) &
!              PRIVATE (pl,pc,pr,ytmp,k,j,i,it,fc)
f_v = 0.0
!$OMP BARRIER
!$OMP DO
do k=1,kx
do j=1,jx
do i=1,ix
fc = e + (bnd(i+1,j  ,k  ) + bnd(i-1,j  ,k  ) &
+ bnd(i  ,j+1,k  ) + bnd(i  ,j-1,k  ) &
+ bnd(i  ,j  ,k+1) + bnd(i  ,j  ,k-1))*bnd(i,j,k)
res_fct(i,j,k) = bnd(i,j,k)*ee/fc
al = res_fct(i,j,k)*wrk(i,j,k)
f_v = f_v + al*al
end do
end do
end do
!$OMP END DO
!$OMP BARRIER

if (f_v.lt.1.0e-30) return
t_eps   = f_v*t_eps

call cbs3d_psor (para_key, p, sz, g, omg, res, wrk, bnd, commmode, do_time, commtime, dtype)
r4(1) = sqrt(res/f_v)

if(res.lt.t_eps) then
print*,'Final',r4(1)
isfin = 1
goto 100
endif

if(tneum.eq.1) then
call cbs3d_set_neum (p, sz, g)
endif

!    if((res.gt.(fct1*f_v)).and.(res.gt.(fct2*t_eps))) then

if(res.gt.(fct2*t_eps)) then
t_eps = res/fct2
!    else
!      t_eps = sqrt(res*t_eps)
endif

eps_abs = 0.15*t_eps
iter = 0
nrm = 0

do i_iter=1,n_iter
if(tneum.eq.1) then
call cbs3d_psor_neumann (para_key, p, sz, g, omg, res, wrk, bnd, commmode, do_time, commtime, dtype)
call cbs3d_set_neum (p, sz, g)
else
call cbs3d_psor (para_key, p, sz, g, omg, res, wrk, bnd, commmode, do_time, commtime, dtype)
endif

iparam(2) = iparam(2) + 1

if(res.lt.t_eps) goto 70
!      if(res.lt.t_eps.and.i_iter.gt.1) goto 70
!    print*,'Intermid',i_iter,res/f_v

do k=1,kx
do j=1,jx
do i=1,ix
rest(i,j,k) = res_fct(i,j,k)*wrk(i,j,k)
end do
end do
end do

call cbs3d_mvprod (p, sz, g, e, bnd, yt, res_fct)
iparam(3) = iparam(3) + 1

call cbs3d_ressmpl (sz, g, res_abs, rest, yt)
beta = sqrt(res_abs)

if(beta.lt.1.0e-30) goto 60
beta_1 = 1.0/beta

!$OMP DO
do k=1,kx
do j=1,jx
do i=1,ix
vm(i,j,k,1) = beta_1*rest(i,j,k)
end do
end do
end do

bgm(1) = beta

do lm=2,m
bgm(lm) = 0.0;
enddo

if(m.gt.m_max) m = m_max

!$OMP END DO
do im=1,m
iter = iter + 1
!$OMP BARRIER

do k=1,kx
do j=1,jx
do i=1,ix
rest(i,j,k) = vm(i,j,k,im)
end do
end do
end do

!$OMP DO
do k=lg,kg
do j=lg,jg
do i=lg,ig
xt(i,j,k) = 0.0
end do
end do
end do
!$OMP END DO

if(    iter.le.  oki) then; n_inner =   step
elseif(iter.le.2*oki) then; n_inner = 2*step
elseif(iter.le.3*oki) then; n_inner = 3*step
elseif(iter.le.4*oki) then; n_inner = 4*step
else;                       n_inner = 5*step
endif

if(n_inner.lt.1) n_inner = 1

do i_inner = 1, n_inner
if(tneum.eq.1) then
call cbs3d_psor_neumann_m1(para_key, xt, sz, g, omg, e, res, rest, bnd, commmode, do_time, commtime, dtype)
call cbs3d_set_neum (xt, sz, g)
else
call cbs3d_psor_m1 (para_key, xt, sz, g, omg, e, res, rest, bnd, commmode, do_time, commtime, dtype)
endif

iparam(2) = iparam(2) + 1
enddo


do k=lg,kg
do j=lg,jg
do i=lg,ig
zm(i,j,k,im) = xt(i,j,k)
end do
end do
end do

call cbs3d_mvprod (xt, sz, g, e, bnd, yt, res_fct)
iparam(3) = iparam(3) + 1

do km=1,im
al = 0.0

do k=1,kx
do j=1,jx
do i=1,ix
al = al + vm(i,j,k,km)*yt(i,j,k)
end do
end do
end do

hgm(km,im) = al

do k=1,kx
do j=1,jx
do i=1,ix
yt(i,j,k) = yt(i,j,k) - al*vm(i,j,k,km)
end do
end do
end do
end do

al = 0.0

do k=1,kx
do j=1,jx
do i=1,ix
al = al + yt(i,j,k)*yt(i,j,k)
end do
end do
end do

hgm(im+1,im) = sqrt(al)

if(hgm(im+1,im).lt.1.0e-30) then
nrm = im-1
goto 50
endif

if(im.lt.m) then
al = 1.0/hgm(im+1,im)

do k=1,kx
do j=1,jx
do i=1,ix
vm(i,j,k,im+1) = al*yt(i,j,k)
end do
end do
end do
endif

rgm(1,im) = hgm(1,im)

do km=2,im
rgm(km  ,im) = cgm(km-1)*hgm(km,im) - sgm(km-1)*rgm(km-1,im)
rgm(km-1,im) = sgm(km-1)*hgm(km,im) + cgm(km-1)*rgm(km-1,im)
end do

al = sqrt(rgm(im,im)*rgm(im,im) + hgm(im+1,im)*hgm(im+1,im))

if(al.lt.1.0e-30) then
nrm = im-1
goto 50
endif

cgm(im) = rgm(im,im)/al
sgm(im) = hgm(im+1,im)/al
rgm(im,im) = cgm(im)*rgm(im,im) + sgm(im)*hgm(im+1,im)
bgm(im+1) = - sgm(im)*bgm(im)
bgm(im) = cgm(im)*bgm(im)
al = bgm(im+1)*bgm(im+1)
iparam(1) = iparam(1) + 1

if(al.lt.eps_abs) then
nrm = im
goto 50
endif
!      end "im"
end do

nrm = m

50  ygm(nrm) = bgm(nrm)/rgm(nrm,nrm)

do im = nrm-1, 1, -1
al = bgm(im)

do jm = im+1, nrm
al = al - rgm(im,jm)*ygm(jm)
enddo

ygm(im) = al/rgm(im,im)
end do

do im = 1, nrm
al = ygm(im)

do k=lg,kg
do j=lg,jg
do i=lg,ig
p(i,j,k) = p(i,j,k) + al*zm(i,j,k,im)
end do
end do
end do
end do
!      end "i_iter"
end do

60  continue

if(tneum.eq.1) then
call cbs3d_psor_neumann (para_key, p, sz, g, omg, res, wrk, bnd, commmode, do_time, commtime, dtype)
call cbs3d_set_neum (p, sz, g)
else
call cbs3d_psor (para_key, p, sz, g, omg, res, wrk, bnd, commmode, do_time, commtime, dtype)
endif

70 continue
!    call cbs3d_mvprod (p, sz, g, e, bnd, yt, res_fct)
!    do k=1,kx; do j=1,jx; do i=1,ix
!      rest(i,j,k) = res_fct(i,j,k)*wrk(i,j,k)
!    end do; end do; end do
!    call cbs3d_ressmpl (sz, g, res_abs, rest, yt)
!    r4(2) = sqrt(res/f_v); r4(3) = sqrt(res_abs/f_v)
!    print*,'Intermid',i_iter-1,nrm,r4(1),r4(2),r4(3)
!    r4(2) = sqrt(res/f_v)
!    print*,'Intermid',i_iter-1,nrm,r4(1),r4(2)
!$OMP END PARALLEL
100 continue

if((isfin.eq.0).and.(res.lt.(f_v*rparam(1)*rparam(1)))) res = 100.0*f_v*rparam(1)*rparam(1)

if(allocated(xt))	deallocate(xt)
if(allocated(yt))	deallocate(yt)
if(allocated(rest))	deallocate(rest)
if(allocated(res_fct))	deallocate(res_fct)
if(allocated(vm))	deallocate(vm)
if(allocated(zm))	deallocate(zm)

return
end subroutine cbs3d_gmressor


!   ***************************************************
subroutine cbs3d_gmresadi (para_key, p, sz, g, res, wrk, bnd,         &
commmode, do_time, commtime, dtype, iparam, rparam)
implicit none
include 'sklparaf.h'

integer                                                   ::  para_key, g, dtype
integer                                                   ::  iparam(*)
integer                                                   ::  commmode ! 0:CommBndCell, 1:Commface(NoHide), 2:Commface(Hide)
integer                                                   ::  do_time, commtime
integer, dimension(3)                                     ::  sz
real                                                      ::  res
real*8                                                    ::  rparam(*)
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  wrk, p, bnd

real,dimension(:),allocatable                             ::  ytmp, pl, pc, pr
real,dimension(:,:,:),allocatable                         ::  xt, yt, rest, res_fct
real,dimension(:,:,:,:),allocatable                       ::  vm, zm
integer                                                   ::  i, j, k, ix, jx, kx, it, m
integer                                                   ::  ig, jg, kg, lg, im, jm, km, lm
integer                                                   ::  i_iter, n_iter, tneum
integer                                                   ::  iter, ierr, nrm, dir, la
integer, parameter                                        ::  m_max = 15, iter_max = 100
integer, parameter                                        ::  isneum = 0, is_fast = 1
integer                                                   ::  i_inner, n_inner, oki, step
real                                                      ::  e, ee, fc, tomg
real                                                      ::  eps_abs, res_abs
real*8                                                    ::  rgm(m_max, m_max)
real*8                                                    ::  hgm(m_max+1, m_max)
real*8                                                    ::  cgm(m_max), sgm(m_max)
real*8                                                    ::  bgm(m_max), ygm(m_max)
real*8                                                    ::  prev, TuneGetTime
real*8                                                    ::  al, t_eps, f_v, t_tiny
real*8                                                    ::  beta, beta_1

integer, dimension(2)                                     ::  err

if(iparam(8).lt.0) return
err(1) = 0
iparam(1) = 0
iparam(2) = 0
iparam(3) = 0
!    n_iter = min(iparam(8), iter_max)
!    nrc    = min(iparam(8), nrc_max)
oki = 3; step = 3
m       = m_max
n_iter  = iter_max
if (dtype == 4) then
e = 2.4e-7
t_eps = 0.995*rparam(1)*rparam(1)
else
e = 4.4e-16
t_eps = 0.999*rparam(1)*rparam(1)
endif
eps_abs = 0.15*rparam(1)*rparam(1)
tomg    = 1.2
ee = 1.0 + e
t_tiny = 1.0e-36

ix = sz(1)
jx = sz(2)
kx = sz(3)
lg = 1 - g; ig = ix + g; jg = jx + g; kg = kx + g
la = max(ix, jx, kx); tneum = 0
if((isneum.ne.0).and.(g.gt.0)) then; la = la + 2; tneum = 1; endif
! Distributed like "wrk", "p", "bnd"
allocate(xt(lg:ig,lg:jg,lg:kg),      yt(lg:ig,lg:jg,lg:kg), &
rest(lg:ig,lg:jg,lg:kg), res_fct(lg:ig,lg:jg,lg:kg))
allocate(vm(ix,jx,kx,m_max), zm(lg:ig,lg:jg,lg:kg,m_max))
! Allocate on each processor
allocate(ytmp(la), pl(la), pc(la), pr(la))

!$OMP PARALLEL SHARED  (ix,jx,kx,p,bnd,wrk,tomg,e) &
!              PRIVATE (pl,pc,pr,ytmp,k,j,i,it,fc)
f_v = 0.0
!$OMP BARRIER
!$OMP DO
do k=1,kx; do j=1,jx; do i=1,ix
fc = e + (bnd(i+1,j  ,k  ) + bnd(i-1,j  ,k  ) &
+ bnd(i  ,j+1,k  ) + bnd(i  ,j-1,k  ) &
+ bnd(i  ,j  ,k+1) + bnd(i  ,j  ,k-1))*bnd(i,j,k)
res_fct(i,j,k) = bnd(i,j,k)*ee/fc
al = res_fct(i,j,k)*wrk(i,j,k)
f_v = f_v + al*al
end do; end do; end do
!$OMP END DO
!$OMP BARRIER
if(f_v.lt.1.0e-30) return
t_eps   = f_v*t_eps
eps_abs = f_v*eps_abs
iter = 0
dir = 1
nrm = 0

do i_iter=1,n_iter
call cbs3d_psor (para_key, p, sz, g, tomg, res, wrk, bnd, commmode, do_time, commtime, dtype)
iparam(2) = iparam(2) + 1
if(res.lt.t_eps) goto 70
!    print*,'Intermid',i_iter,res
do k=1,kx; do j=1,jx; do i=1,ix
rest(i,j,k) = res_fct(i,j,k)*wrk(i,j,k)
end do; end do; end do
call cbs3d_mvprod (p, sz, g, e, bnd, yt, res_fct)
call cbs3d_ressmpl (sz, g, res_abs, rest, yt)
iparam(3) = iparam(3) + 1
beta = sqrt(res_abs)
if(beta.lt.1.0e-30) goto 70
beta_1 = 1.0/beta
!$OMP DO
do k=1,kx; do j=1,jx; do i=1,ix
vm(i,j,k,1) = beta_1*rest(i,j,k)
end do; end do; end do
!$OMP END DO
bgm(1) = beta; do lm=2,m; bgm(lm) = 0.0; enddo
if(m.gt.m_max) m = m_max

do im=1,m
iter = iter + 1
!$OMP BARRIER
!$OMP DO
do k=1,kx; do j=1,jx; do i=1,ix
if(abs(res_fct(i,j,k)).gt.e) then
rest(i,j,k) = vm(i,j,k,im)/res_fct(i,j,k)
else
rest(i,j,k) = vm(i,j,k,im)/e
endif
end do; end do; end do
!$OMP END DO
!$OMP DO
do k=lg,kg; do j=lg,jg; do i=lg,ig
xt(i,j,k) = 0.0
end do; end do; end do
!$OMP END DO

if(    iter.le.  oki) then; n_inner =   step
elseif(iter.le.2*oki) then; n_inner = 2*step
elseif(iter.le.3*oki) then; n_inner = 3*step
elseif(iter.le.4*oki) then; n_inner = 4*step
else;                       n_inner = 5*step
endif

if(n_inner.lt.1) n_inner = 1
do i_inner = 1, n_inner
if(is_fast.eq.1) then
call cbs3d_adiiter_x1 (xt, sz, g, rest, bnd, tneum, t_tiny, e,    &
ytmp, pl, pc, pr, err,     &
commmode, do_time, commtime)
else
if(dir.eq.1) then
call cbs3d_adiiter_x (xt, sz, g, rest, bnd, tneum, t_tiny, e,    &
ytmp, pl, pc, pr, err,     &
commmode, do_time, commtime)
if(err(1).ne.0) return
dir = 2
elseif(dir.eq.2) then
call cbs3d_adiiter_y (xt, sz, g, rest, bnd, tneum, t_tiny, e,    &
ytmp, pl, pc, pr, err,     &
commmode, do_time, commtime)
if(err(1).ne.0) return
dir = 3
else
call cbs3d_adiiter_z (xt, sz, g, rest, bnd, tneum, t_tiny, e,    &
ytmp, pl, pc, pr, err,     &
commmode, do_time, commtime)
if(err(1).ne.0) return
dir = 1
endif
endif
iparam(2) = iparam(2) + 1
enddo
if(tneum.eq.1) then
call cbs3d_set_neum (xt, sz, g)
endif

do k=lg,kg; do j=lg,jg; do i=lg,ig
zm(i,j,k,im) = xt(i,j,k)
end do; end do; end do
call cbs3d_mvprod (xt, sz, g, e, bnd, yt, res_fct)
iparam(3) = iparam(3) + 1
do km=1,im
al = 0.0
do k=1,kx; do j=1,jx; do i=1,ix
al = al + vm(i,j,k,km)*yt(i,j,k)
end do; end do; end do
hgm(km,im) = al
do k=1,kx; do j=1,jx; do i=1,ix
yt(i,j,k) = yt(i,j,k) - al*vm(i,j,k,km)
end do; end do; end do
end do
al = 0.0
do k=1,kx; do j=1,jx; do i=1,ix
al = al + yt(i,j,k)*yt(i,j,k)
end do; end do; end do
hgm(im+1,im) = sqrt(al)
if(hgm(im+1,im).lt.1.0e-30) then; nrm = im-1; goto 50; endif
if(im.lt.m) then
al = 1.0/hgm(im+1,im)
do k=1,kx; do j=1,jx; do i=1,ix
vm(i,j,k,im+1) = al*yt(i,j,k)
end do; end do; end do
endif
rgm(1,im) = hgm(1,im)
do km=2,im
rgm(km  ,im) = cgm(km-1)*hgm(km,im) - sgm(km-1)*rgm(km-1,im)
rgm(km-1,im) = sgm(km-1)*hgm(km,im) + cgm(km-1)*rgm(km-1,im)
end do
al = sqrt(rgm(im,im)*rgm(im,im) + hgm(im+1,im)*hgm(im+1,im))
if(al.lt.1.0e-30) then; nrm = im-1; goto 50; endif
cgm(im) = rgm(im,im)/al; sgm(im) = hgm(im+1,im)/al
rgm(im,im) = cgm(im)*rgm(im,im) + sgm(im)*hgm(im+1,im)
bgm(im+1) = - sgm(im)*bgm(im); bgm(im) = cgm(im)*bgm(im)
al = bgm(im+1)*bgm(im+1)
iparam(1) = iparam(1) + 1
if(al.lt.eps_abs) then; nrm = im; goto 50; endif
!      end "im"
end do
nrm = m
50  ygm(nrm) = bgm(nrm)/rgm(nrm,nrm)
do im = nrm-1, 1, -1
al = bgm(im)
do jm = im+1, nrm; al = al - rgm(im,jm)*ygm(jm); enddo
ygm(im) = al/rgm(im,im)
end do
do im = 1, nrm
al = ygm(im)
do k=lg,kg; do j=lg,jg; do i=lg,ig
p(i,j,k) = p(i,j,k) + al*zm(i,j,k,im)
end do; end do; end do
end do
!      end "i_iter"
end do

60  continue
call cbs3d_psor (para_key, p, sz, g, tomg, res, wrk, bnd, commmode, do_time, commtime, dtype)
70  continue
!    call cbs3d_mvprod (p, sz, g, e, bnd, yt, res_fct)
!    do k=1,kx; do j=1,jx; do i=1,ix
!      rest(i,j,k) = res_fct(i,j,k)*wrk(i,j,k)
!    end do; end do; end do
!    call cbs3d_ressmpl (sz, g, res_abs, rest, yt)
!    print*,i_iter-1,nrm,res,res_abs,sngl(sqrt(res_abs/f_v))
!$OMP END PARALLEL
100 continue

if(allocated(ytmp))	deallocate(ytmp)
if(allocated(pl))	deallocate(pl)
if(allocated(pc))	deallocate(pc)
if(allocated(pr))	deallocate(pr)

if(allocated(xt))	deallocate(xt)
if(allocated(yt))	deallocate(yt)
if(allocated(rest))	deallocate(rest)
if(allocated(res_fct))	deallocate(res_fct)
if(allocated(vm))	deallocate(vm)
if(allocated(zm))	deallocate(zm)

return
end subroutine cbs3d_gmresadi
