!   *********************************************************
!
!   SPHERE - Skeleton for PHysical and Engineering REsearch
!  
!   Copyright (c) RIKEN, Japan. All right reserved. 2004-2010
!
!   *********************************************************
!
!> @file c3d_vof.f90
!> @brief VOF routines
!> @author keno, FSI Team, VCAD, RIKEN
!<
!  ***********************************************************
!> @subroutine vof_uwd (f, sz, g, v00, dt, dh, v, q, bx, flop)
!! @brief VOFの移流を一次風上で計算
!! @param[out] f VOF
!! @param sz 配列長
!! @param g ガイドセル長
!! @param v00 参照速度
!! @param dt 時間積分幅
!! @param dh 格子幅
!! @param v セルフェイスの速度ベクトル
!! @param q ワーク用配列
!! @param bx BC index for V
!! @param[out] flop flop count
!<
    subroutine vof_uwd (f, sz, g, v00, dt, dh, v, q, bx, flop)
    implicit none
    include '../F_CBC/cbc_f_params.h'
    integer                                                   :: i, j, k, ix, jx, kx, g, idx
    integer, dimension(3)                                     :: sz
    real                                                      :: ur, ul, vr, vl, wr, wl, dt, dh, dth
    real                                                      :: fc, fe, fw, fn, fs, ft, fb
    real                                                      :: bc, be, bw, bn, bs, bt, bb
    real                                                      :: re, rw, rn, rs, rt, rb, flop
    real                                                      :: u_ref, v_ref, w_ref
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) :: v
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    :: f, q
    real, dimension(0:3)                                      :: v00
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) :: bx

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)
    dth= dt/dh
    flop = real(ix*jx*kx*75 + 1)

    do k=1,kx
    do j=1,jx
    do i=1,ix
      ur = v(i  ,j  ,k  ,1) - u_ref
      ul = v(i-1,j  ,k  ,1) - u_ref
      vr = v(i  ,j  ,k  ,2) - v_ref
      vl = v(i  ,j-1,k  ,2) - v_ref
      wr = v(i  ,j  ,k  ,3) - w_ref
      wl = v(i  ,j  ,k-1,3) - w_ref
      
      idx = bx(i,j,k)
      bc = real(ibits(idx, Active,    1))
      bw = real(ibits(idx, bc_face_W, 1))  ! w 
      be = real(ibits(idx, bc_face_E, 1))  ! e 
      bs = real(ibits(idx, bc_face_S, 1))  ! s
      bn = real(ibits(idx, bc_face_N, 1))  ! n
      bb = real(ibits(idx, bc_face_B, 1))  ! b
      bt = real(ibits(idx, bc_face_T, 1))  ! t
      
      fc = f(i  ,j  ,k  )
      fe = f(i+1,j  ,k  )
      fw = f(i-1,j  ,k  )
      fn = f(i  ,j+1,k  )
      fs = f(i  ,j-1,k  )
      ft = f(i  ,j  ,k+1)
      fb = f(i  ,j  ,k-1)
      
      re = 0.5*(ur*(fe+fc)-abs(ur)*(fe-fc))
      rw = 0.5*(ul*(fc+fw)-abs(ul)*(fc-fw))
      rn = 0.5*(vr*(fn+fc)-abs(vr)*(fn-fc))
      rs = 0.5*(vl*(fc+fs)-abs(vl)*(fc-fs))
      rt = 0.5*(wr*(ft+fc)-abs(wr)*(ft-fc))
      rb = 0.5*(wl*(fc+fb)-abs(wl)*(fc-fb))
      
      q(i,j,k) = fc - dth*(re*be - rw*bw + rn*bn - rs*bs + rt*bt - rb*bb) * bc
    end do
    end do
    end do
    
    do k=1,kx
    do j=1,jx
    do i=1,ix
      ur = v(i  ,j  ,k  ,1)
      ul = v(i-1,j  ,k  ,1)
      vr = v(i  ,j  ,k  ,2)
      vl = v(i  ,j-1,k  ,2)
      wr = v(i  ,j  ,k  ,3)
      wl = v(i  ,j  ,k-1,3)
      f(i,j,k) = q(i,j,k) + dth*f(i,j,k)*(ur-ul+vr-vl+wr-wl)
    end do
    end do
    end do
    
    return
    end subroutine vof_uwd

!  *******************************************************************
!> @subroutine c3d_vof_muscl (f, sz, g, v00, dt, dh, v, fn, bnd, flop)
!! @brief 温度のパッシブスカラ方程式の移流部分をMUSCL法で解く
!! @param fx 対流流束
!! @param sz 配列長
!! @param g ガイドセル長
!! @param v00 参照速度
!! @param v 速度ベクトル
!! @param t 温度
!! @param bnd マスク
!! @note
!!    - リミター付き
!<
    subroutine c3d_vof_muscl (f, sz, g, v00, dt, dh, v, fn, bnd, flop)
    implicit none
    integer                                                   ::  i, j, k, ix, jx, kx, g
    integer, dimension(3)                                     ::  sz
    real                                                      ::  dh, ck, b, dt, flop
    real                                                      ::  di1, di2, di3, si1, si2, si3, si4, gi1, gi2, gi3, gi4, tr1, tl1, c1
    real                                                      ::  dj1, dj2, dj3, sj1, sj2, sj3, sj4, gj1, gj2, gj3, gj4, tr2, tl2, c2
    real                                                      ::  dk1, dk2, dk3, sk1, sk2, sk3, sk4, gk1, gk2, gk3, gk4, tr3, tl3, c3
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  bnd, f, fn
    real, dimension(0:3)                                      ::  v00

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    ck=1.0/3.0
    b=(3.0-ck)/(1.0-ck)

    do k=0,kx
    do j=0,jx
    do i=0,ix
         c1 =v(i  ,j  ,k  ,1)-v00(1)
         si1=bnd(i  ,j,k)*bnd(i-1,j,k)
         si2=bnd(i+1,j,k)*bnd(i  ,j,k)
         si3=bnd(i+2,j,k)*bnd(i+1,j,k)
         di1=(f(i  ,j  ,k  )-f(i-1,j  ,k  ))*si1
         di2=(f(i+1,j  ,k  )-f(i  ,j  ,k  ))*si2
         di3=(f(i+2,j  ,k  )-f(i+1,j  ,k  ))*si3
         si1=sign(1.0,di1)
         gi1=si1*max(0.0,min(abs(di1),b*di2*si1))
         si2=sign(1.0,di2)
         gi2=si2*max(0.0,min(abs(di2),b*di1*si2))
         si3=sign(1.0,di3)
         gi3=si3*max(0.0,min(abs(di3),b*di2*si3))
         si4=sign(1.0,di2)
         gi4=si4*max(0.0,min(abs(di2),b*di3*si4))
         tl1=f(i  ,j  ,k  )+0.25*si2*((1.0-ck)*gi1+(1.0+ck)*gi2)
         tr1=f(i+1,j  ,k  )-0.25*si2*((1.0-ck)*gi3+(1.0+ck)*gi4)
         fn(i,j,k)=0.5*(c1*(tr1+tl1)-abs(c1)*(tr1-tl1))

         c2 =v(i  ,j  ,k  ,2)-v00(2)
         sj1=bnd(i,j  ,k)*bnd(i,j-1,k)
         sj2=bnd(i,j+1,k)*bnd(i,j  ,k)
         sj3=bnd(i,j+2,k)*bnd(i,j+1,k)
         dj1=(f(i  ,j  ,k  )-f(i  ,j-1,k  ))*sj1
         dj2=(f(i  ,j+1,k  )-f(i  ,j  ,k  ))*sj2
         dj3=(f(i  ,j+2,k  )-f(i  ,j+1,k  ))*sj3
         sj1=sign(1.0,dj1)
         gj1=sj1*max(0.0,min(abs(dj1),b*dj2*sj1))
         sj2=sign(1.0,dj2)
         gj2=sj2*max(0.0,min(abs(dj2),b*dj1*sj2))
         sj3=sign(1.0,dj3)
         gj3=sj3*max(0.0,min(abs(dj3),b*dj2*sj3))
         sj4=sign(1.0,dj2)
         gj4=sj4*max(0.0,min(abs(dj2),b*dj3*sj4))
         tl2=f(i  ,j  ,k  )+0.25*sj2*((1.0-ck)*gj1+(1.0+ck)*gj2)
         tr2=f(i  ,j+1,k  )-0.25*sj2*((1.0-ck)*gj3+(1.0+ck)*gj4)
         fn(i,j,k)=0.5*(c2*(tr2+tl2)-abs(c2)*(tr2-tl2))

         c3 =v(i  ,j  ,k  ,3)-v00(3)
         sk1=bnd(i,j,k  )*bnd(i,j,k-1)
         sk2=bnd(i,j,k+1)*bnd(i,j,k  )
         sk3=bnd(i,j,k+2)*bnd(i,j,k+1)
         dk1=(f(i  ,j  ,k  )-f(i  ,j  ,k-1))*sk1
         dk2=(f(i  ,j  ,k+1)-f(i  ,j  ,k  ))*sk2
         dk3=(f(i  ,j  ,k+2)-f(i  ,j  ,k+1))*sk3
         sk1=sign(1.0,dk1)
         gk1=sk1*max(0.0,min(abs(dk1),b*dk2*sk1))
         sk2=sign(1.0,dk2)
         gk2=sk2*max(0.0,min(abs(dk2),b*dk1*sk2))
         sk3=sign(1.0,dk3)
         gk3=sk3*max(0.0,min(abs(dk3),b*dk2*sk3))
         sk4=sign(1.0,dk2)
         gk4=sk4*max(0.0,min(abs(dk2),b*dk3*sk4))
         tl3=f(i  ,j  ,k  )+0.25*sk2*((1.0-ck)*gk1+(1.0+ck)*gk2)
         tr3=f(i  ,j  ,k+1)-0.25*sk2*((1.0-ck)*gk3+(1.0+ck)*gk4)
         fn(i,j,k)=0.5*(c3*(tr3+tl3)-abs(c3)*(tr3-tl3))
    end do
    end do
    end do

    return
    end subroutine c3d_vof_muscl