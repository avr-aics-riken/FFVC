!   *********************************************************
!
!   SPHERE - Skeleton for PHysical and Engineering REsearch
!  
!   Copyright (c) RIKEN, Japan. All right reserved. 2004-2011
!
!   *********************************************************

!> @file BCvec.f90
!> @brief Boundary conditons for collocated velocity 
!> @author keno, FSI Team, VCAD, RIKEN

!  ************************************************
!> @subroutine cbc_vobc_update (v, sz, g, vc, face)
!! @brief 疑似速度から次ステップ速度へ参照する速度をコピーする
!! @param[out] v 速度ベクトル（セルセンタ）
!! @param sz 配列長
!! @param g ガイドセル長
!! @param vc セルセンタ疑似速度 u^*
!! @param face 面番号
!<
    subroutine cbc_vobc_update (v, sz, g, vc, face)
    implicit none
    include 'cbc_f_params.h'
    integer                                                     ::  i, j, k, g, ix, jx, kx, face
    integer, dimension(3)                                       ::  sz
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)   ::  v, vc

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    
    FACES : select case (face)
    case (X_minus)
      i = 0
      do k=1,kx
      do j=1,jx
        v(i,j,k,1) = vc(i,j,k,1)
        v(i,j,k,2) = vc(i,j,k,2)
        v(i,j,k,3) = vc(i,j,k,3)
      end do
      end do
      
    case (X_plus)
      i = ix+1
      do k=1,kx
      do j=1,jx
        v(i,j,k,1) = vc(i,j,k,1)
        v(i,j,k,2) = vc(i,j,k,2)
        v(i,j,k,3) = vc(i,j,k,3)
      end do
      end do
      
    case (Y_minus)
      j = 0
      do k=1,kx
      do i=1,ix
        v(i,j,k,1) = vc(i,j,k,1)
        v(i,j,k,2) = vc(i,j,k,2)
        v(i,j,k,3) = vc(i,j,k,3)
      end do
      end do
      
    case (Y_plus)
      j = jx+1
      do k=1,kx
      do i=1,ix
        v(i,j,k,1) = vc(i,j,k,1)
        v(i,j,k,2) = vc(i,j,k,2)
        v(i,j,k,3) = vc(i,j,k,3)
      end do
      end do
      
    case (Z_minus)
      k = 0
      do j=1,jx
      do i=1,ix
        v(i,j,k,1) = vc(i,j,k,1)
        v(i,j,k,2) = vc(i,j,k,2)
        v(i,j,k,3) = vc(i,j,k,3)
      end do
      end do
      
    case (Z_plus)
      k = kx+1
      do j=1,jx
      do i=1,ix
        v(i,j,k,1) = vc(i,j,k,1)
        v(i,j,k,2) = vc(i,j,k,2)
        v(i,j,k,3) = vc(i,j,k,3)
      end do
      end do
      
    case default
    end select FACES
    
    return
    end subroutine cbc_vobc_update
    
!  **************************************************************
!> @subroutine cbc_vobc_outflow (v, sz, g, c, bv, face, v0, flop)
!! @brief 外部流出境界で，次ステップの流出速度を対流流出条件で予測し，ガイドセルに参照値として代入する
!! @param[out] v 速度 u^*
!! @param sz 配列長
!! @param g ガイドセル長
!! @param c uc*dt/dh
!! @param bv BCindex V
!! @param face 外部境界の面番号
!! @param v0 セルセンタ速度 u^n
!! @param[out] flop
!<
    subroutine cbc_vobc_outflow (v, sz, g, c, bv, face, v0, flop)
    implicit none
    include 'cbc_f_params.h'
    integer                                                     ::  i, j, k, g, idx, face, ix, jx, kx, m
    integer, dimension(3)                                       ::  sz
    real                                                        ::  Ue, Uw, Un, Us, Ut, Ub
    real                                                        ::  Ve, Vw, Vn, Vs, Vt, Vb
    real                                                        ::  We, Ww, Wn, Ws, Wt, Wb
    real                                                        ::  flop, c
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)   ::  v, v0
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bv

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    
    m = 0
    
    FACES : select case (face)
    case (X_minus)
      if ( c>0.0 ) c=0.0
      i = 1
      do k=1,kx
      do j=1,jx
        idx = bv(i,j,k)
        if ( ibits(idx, bc_face_W, bitw_5) == obc_mask ) then
          Ue = v0(i  ,j  ,k  ,1)
          Uw = v0(i-1,j  ,k  ,1)
          Ve = v0(i  ,j  ,k  ,2)
          Vw = v0(i-1,j  ,k  ,2)
          We = v0(i  ,j  ,k  ,3)
          Ww = v0(i-1,j  ,k  ,3)
          v(i-1,j  ,k  ,1) = Uw - c*(Ue-Uw)
          v(i-1,j  ,k  ,2) = Vw - c*(Ve-Vw)
          v(i-1,j  ,k  ,3) = Ww - c*(We-Ww)
          m = m+1
        endif
      end do
      end do
      
    case (X_plus)
      if ( c<0.0 ) c=0.0
      i = ix
      do k=1,kx
      do j=1,jx
        idx = bv(i,j,k)
        if ( ibits(idx, bc_face_E, bitw_5) == obc_mask ) then
          Ue = v0(i+1,j  ,k  ,1)
          Uw = v0(i  ,j  ,k  ,1)
          Ve = v0(i+1,j  ,k  ,2)
          Vw = v0(i  ,j  ,k  ,2)
          We = v0(i+1,j  ,k  ,3)
          Ww = v0(i  ,j  ,k  ,3)
          v(i+1,j  ,k  ,1) = Ue - c*(Ue-Uw)
          v(i+1,j  ,k  ,2) = Ve - c*(Ve-Vw)
          v(i+1,j  ,k  ,3) = We - c*(We-Ww)
          m = m+1
        endif
      end do
      end do
      
    case (Y_minus)
    if ( c>0.0 ) c=0.0
      j = 1
      do k=1,kx
      do i=1,ix
        idx = bv(i,j,k)
        if ( ibits(idx, bc_face_S, bitw_5) == obc_mask ) then
          Un = v0(i  ,j  ,k  ,1)
          Us = v0(i  ,j-1,k  ,1)
          Vn = v0(i  ,j  ,k  ,2)
          Vs = v0(i  ,j-1,k  ,2)
          Wn = v0(i  ,j  ,k  ,3)
          Ws = v0(i  ,j-1,k  ,3)
          v(i  ,j-1,k  ,1) = Us - c*(Un-Us)
          v(i  ,j-1,k  ,2) = Vs - c*(Vn-Vs)
          v(i  ,j-1,k  ,3) = Ws - c*(Wn-Ws)
          m = m+1
        endif
      end do
      end do
      
    case (Y_plus)
      if ( c<0.0 ) c=0.0
      j = jx
      do k=1,kx
      do i=1,ix
        idx = bv(i,j,k)
        if ( ibits(idx, bc_face_N, bitw_5) == obc_mask ) then
          Un = v0(i  ,j+1,k  ,1)
          Us = v0(i  ,j  ,k  ,1)
          Vn = v0(i  ,j+1,k  ,2)
          Vs = v0(i  ,j  ,k  ,2)
          Wn = v0(i  ,j+1,k  ,3)
          Ws = v0(i  ,j  ,k  ,3)
          v(i  ,j+1,k  ,1) = Un - c*(Un-Us)
          v(i  ,j+1,k  ,2) = Vn - c*(Vn-Vs)
          v(i  ,j+1,k  ,3) = Wn - c*(Wn-Ws)
          m = m+1
        endif
      end do
      end do
      
    case (Z_minus)
    if ( c>0.0 ) c=0.0
      k = 1
      do j=1,jx
      do i=1,ix
        idx = bv(i,j,k)
        if ( ibits(idx, bc_face_B, bitw_5) == obc_mask ) then
          Ut = v0(i  ,j  ,k  ,1)
          Ub = v0(i  ,j  ,k-1,1)
          Vt = v0(i  ,j  ,k  ,2)
          Vb = v0(i  ,j  ,k-1,2)
          Wt = v0(i  ,j  ,k  ,3)
          Wb = v0(i  ,j  ,k-1,3)
          v(i  ,j  ,k-1,1) = Ub - c*(Ut-Ub)
          v(i  ,j  ,k-1,2) = Vb - c*(Vt-Vb)
          v(i  ,j  ,k-1,3) = Wb - c*(Wt-Wb)
          m = m+1
        endif
      end do
      end do
      
    case (Z_plus)
      if ( c<0.0 ) c=0.0
      k = kx
      do j=1,jx
      do i=1,ix
        idx = bv(i,j,k)
        if ( ibits(idx, bc_face_T, bitw_5) == obc_mask ) then
          Ut = v0(i  ,j  ,k+1,1)
          Ub = v0(i  ,j  ,k  ,1)
          Vt = v0(i  ,j  ,k+1,2)
          Vb = v0(i  ,j  ,k  ,2)
          Wt = v0(i  ,j  ,k+1,3)
          Wb = v0(i  ,j  ,k  ,3)
          v(i  ,j  ,k+1,1) = Ut - c*(Ut-Ub)
          v(i  ,j  ,k+1,2) = Vt - c*(Vt-Vb)
          v(i  ,j  ,k+1,3) = Wt - c*(Wt-Wb)
          m = m+1
        endif
      end do
      end do
      
    case default
    end select FACES
    
    flop = flop + real(m*9)
    
    return
    end subroutine cbc_vobc_outflow

!  **********************************************************
!> @subroutine cbc_vobc_tfree (v, sz, g, face, bv, v00, flop)
!! @brief 速度の外部境界：　トラクションフリー
!! @param v 速度ベクトル
!! @param sz 配列長
!! @param g ガイドセル長
!! @param face 外部境界面の番号
!! @param bv BCindex V
!! @param v00 参照速度
!! @param flop 浮動小数演算数
!<
    subroutine cbc_vobc_tfree (v, sz, g, face, bv, v00, flop)
    implicit none
    include 'cbc_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, face, g
    integer, dimension(3)                                     ::  sz
    real                                                      ::  b0, b1, b2, b3, b4, b5
    real, dimension(0:3)                                      ::  v00
    real                                                      ::  c1, c2, c3, v0, v1, v2, v3, v4
    real                                                      ::  uu, vv, ww
    real                                                      ::  u_ref, v_ref, w_ref
    real                                                      ::  d1, d2, d3, d4, r, flop
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bv

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)
    
    FACES : select case (face)
    case (X_minus)
      do k=1,kx
      do j=1,jx
        c1 = v(1,j,k,1)
        c2 = v(1,j,k,2)
        c3 = v(1,j,k,3)
        
        v0=      v(0,j  ,k  ,1)
        v1= 0.5*(v(0,j+1,k  ,1)+v(1,j+1,k  ,1))
        v2= 0.5*(v(0,j-1,k  ,1)+v(1,j-1,k  ,1))
        v3= 0.5*(v(0,j  ,k+1,1)+v(1,j  ,k+1,1))
        v4= 0.5*(v(0,j  ,k-1,1)+v(1,j  ,k-1,1))
        
        b0 = real( ibits(bv(0,j  ,k  ), State, 1) * ibits(bv(1,j  ,k  ), State, 1) )
        b1 = real( ibits(bv(0,j+1,k  ), State, 1) * ibits(bv(1,j+1,k  ), State, 1) )
        b2 = real( ibits(bv(0,j-1,k  ), State, 1) * ibits(bv(1,j-1,k  ), State, 1) )
        b3 = real( ibits(bv(0,j  ,k+1), State, 1) * ibits(bv(1,j  ,k+1), State, 1) )
        b4 = real( ibits(bv(0,j  ,k-1), State, 1) * ibits(bv(1,j  ,k-1), State, 1) )
        
        r  = 2.0*u_ref - 0.5*(v0+c1)
        d1 = v1*b1 + (1.0-b1)*r
        d2 = v2*b2 + (1.0-b2)*r
        d3 = v3*b3 + (1.0-b3)*r
        d4 = v4*b4 + (1.0-b4)*r
        
        b5 = 1.0-b0
        uu =  c1               *b0 + b5*u_ref
        vv = (c2 + 0.5*(d1-d2))*b0 + b5*v_ref
        ww = (c3 + 0.5*(d3-d4))*b0 + b5*w_ref

        do i=1-g, 0
          v(i,j,k,1)= uu
          v(i,j,k,2)= vv
          v(i,j,k,3)= ww
        end do
      end do
      end do
      
      flop = flop + real(jx*kx*49)
      
    case (X_plus)
      do k=1,kx
      do j=1,jx
        c1 = v(ix,j,k,1)
        c2 = v(ix,j,k,2)
        c3 = v(ix,j,k,3)
        
        v0=      v(ix+1,j  ,k  ,1)
        v1= 0.5*(v(ix+1,j+1,k  ,1)+v(ix,j+1,k  ,1))
        v2= 0.5*(v(ix+1,j-1,k  ,1)+v(ix,j-1,k  ,1))
        v3= 0.5*(v(ix+1,j  ,k+1,1)+v(ix,j  ,k+1,1))
        v4= 0.5*(v(ix+1,j  ,k-1,1)+v(ix,j  ,k-1,1))

        b0 = real( ibits(bv(ix+1,j  ,k  ), State, 1) * ibits(bv(ix,j  ,k  ), State, 1) )
        b1 = real( ibits(bv(ix+1,j+1,k  ), State, 1) * ibits(bv(ix,j+1,k  ), State, 1) )
        b2 = real( ibits(bv(ix+1,j-1,k  ), State, 1) * ibits(bv(ix,j-1,k  ), State, 1) )
        b3 = real( ibits(bv(ix+1,j  ,k+1), State, 1) * ibits(bv(ix,j  ,k+1), State, 1) )
        b4 = real( ibits(bv(ix+1,j  ,k-1), State, 1) * ibits(bv(ix,j  ,k-1), State, 1) )
        
        r  = 2.0*u_ref - 0.5*(v0+c1)
        d1 = v1*b1 + (1.0-b1)*r
        d2 = v2*b2 + (1.0-b2)*r
        d3 = v3*b3 + (1.0-b3)*r
        d4 = v4*b4 + (1.0-b4)*r

        b5 = 1.0-b0
        uu =  c1               *b0 + b5*u_ref
        vv = (c2 + 0.5*(d1-d2))*b0 + b5*v_ref
        ww = (c3 + 0.5*(d3-d4))*b0 + b5*w_ref
        
        do i=ix+1, ix+g
          v(i,j,k,1)= uu
          v(i,j,k,2)= vv
          v(i,j,k,3)= ww
        end do
      end do
      end do
      
      flop = flop + real(jx*kx*49)
      
    case (Y_minus)
      do k=1,kx
      do i=1,ix
        c1 = v(i,1,k,1)
        c2 = v(i,1,k,2)
        c3 = v(i,1,k,3)

        v0=      v(i  ,0,k  ,2)
        v1= 0.5*(v(i+1,0,k  ,2)+v(i+1,1,k  ,2))
        v2= 0.5*(v(i-1,0,k  ,2)+v(i-1,1,k  ,2))
        v3= 0.5*(v(i  ,0,k+1,2)+v(i  ,1,k+1,2))
        v4= 0.5*(v(i  ,0,k-1,2)+v(i  ,1,k-1,2))

        b0 = real( ibits(bv(i  ,0,k  ), State, 1) * ibits(bv(i  ,1,k  ), State, 1) )
        b1 = real( ibits(bv(i+1,0,k  ), State, 1) * ibits(bv(i+1,1,k  ), State, 1) )
        b2 = real( ibits(bv(i-1,0,k  ), State, 1) * ibits(bv(i-1,1,k  ), State, 1) )
        b3 = real( ibits(bv(i  ,0,k+1), State, 1) * ibits(bv(i  ,1,k+1), State, 1) )
        b4 = real( ibits(bv(i  ,0,k-1), State, 1) * ibits(bv(i  ,1,k-1), State, 1) )
        
        r  = 2.0*v_ref - 0.5*(v0+c2)
        d1 = v1*b1 + (1.0-b1)*r
        d2 = v2*b2 + (1.0-b2)*r
        d3 = v3*b3 + (1.0-b3)*r
        d4 = v4*b4 + (1.0-b4)*r

        b5 = 1.0-b0
        uu = (c1 + 0.5*(d1-d2))*b0 + b5*u_ref
        vv =  c2               *b0 + b5*v_ref
        ww = (c3 + 0.5*(d3-d4))*b0 + b5*w_ref
                
        do j=1-g, 0
          v(i,j,k,1)= uu
          v(i,j,k,2)= vv
          v(i,j,k,3)= ww
        end do
      end do
      end do
      
      flop = flop + real(ix*kx*49)
      
    case (Y_plus)
      do k=1,kx
      do i=1,ix
        c1 = v(i,jx,k,1)
        c2 = v(i,jx,k,2)
        c3 = v(i,jx,k,3)
        
        v0=      v(i  ,jx+1,k  ,2)
        v1= 0.5*(v(i+1,jx+1,k  ,2)+v(i+1,jx,k  ,2))
        v2= 0.5*(v(i-1,jx+1,k  ,2)+v(i-1,jx,k  ,2))
        v3= 0.5*(v(i  ,jx+1,k+1,2)+v(i  ,jx,k+1,2))
        v4= 0.5*(v(i  ,jx+1,k-1,2)+v(i  ,jx,k-1,2))

        b0 = real( ibits(bv(i  ,jx+1,k  ), State, 1) * ibits(bv(i  ,jx,k  ), State, 1) )
        b1 = real( ibits(bv(i+1,jx+1,k  ), State, 1) * ibits(bv(i+1,jx,k  ), State, 1) )
        b2 = real( ibits(bv(i-1,jx+1,k  ), State, 1) * ibits(bv(i-1,jx,k  ), State, 1) )
        b3 = real( ibits(bv(i  ,jx+1,k+1), State, 1) * ibits(bv(i  ,jx,k+1), State, 1) )
        b4 = real( ibits(bv(i  ,jx+1,k-1), State, 1) * ibits(bv(i  ,jx,k-1), State, 1) )
        
        r  = 2.0*v_ref - 0.5*(v0+c2)
        d1 = v1*b1 + (1.0-b1)*r
        d2 = v2*b2 + (1.0-b2)*r
        d3 = v3*b3 + (1.0-b3)*r
        d4 = v4*b4 + (1.0-b4)*r

        b5 = 1.0-b0
        uu = (c1 + 0.5*(d1-d2))*b0 + b5*u_ref
        vv =  c2               *b0 + b5*v_ref
        ww = (c3 + 0.5*(d3-d4))*b0 + b5*w_ref
                
        do j=jx+1, jx+g
          v(i,j,k,1)= uu
          v(i,j,k,2)= vv
          v(i,j,k,3)= ww
        end do
      end do
      end do
      
      flop = flop + real(ix*kx*49)
      
    case (Z_minus)
      do j=1,jx
      do i=1,ix
        c1 = v(i,j,1,1)
        c2 = v(i,j,1,2)
        c3 = v(i,j,1,3)

        v0=      v(i  ,j  ,0,3)
        v1= 0.5*(v(i+1,j  ,0,3)+v(i+1,j  ,1,3))
        v2= 0.5*(v(i-1,j  ,0,3)+v(i-1,j  ,1,3))
        v3= 0.5*(v(i  ,j+1,0,3)+v(i  ,j+1,1,3))
        v4= 0.5*(v(i  ,j-1,0,3)+v(i  ,j-1,1,3))
        
        b0 = real( ibits(bv(i  ,j  ,0), State, 1) * ibits(bv(i  ,j  ,1), State, 1) )
        b1 = real( ibits(bv(i+1,j  ,0), State, 1) * ibits(bv(i+1,j  ,1), State, 1) )
        b2 = real( ibits(bv(i-1,j  ,0), State, 1) * ibits(bv(i-1,j  ,1), State, 1) )
        b3 = real( ibits(bv(i  ,j+1,0), State, 1) * ibits(bv(i  ,j+1,1), State, 1) )
        b4 = real( ibits(bv(i  ,j-1,0), State, 1) * ibits(bv(i  ,j-1,1), State, 1) )
        
        r  = 2.0*w_ref - 0.5*(v0+c3)
        d1 = v1*b1 + (1.0-b1)*r
        d2 = v2*b2 + (1.0-b2)*r
        d3 = v3*b3 + (1.0-b3)*r
        d4 = v4*b4 + (1.0-b4)*r

        b5 = 1.0-b0
        uu = (c1 + 0.5*(d1-d2))*b0 + b5*u_ref
        vv = (c2 + 0.5*(d3-d4))*b0 + b5*v_ref
        ww =  c3               *b0 + b5*w_ref
                
        do k=1-g, 0
          v(i,j,k,1)= uu
          v(i,j,k,2)= vv
          v(i,j,k,3)= ww
        end do
      end do
      end do
      
      flop = flop + real(ix*jx*49)
      
    case (Z_plus)
      do j=1,jx
      do i=1,ix
        c1 = v(i,j,kx,1)
        c2 = v(i,j,kx,2)
        c3 = v(i,j,kx,3)
        
        v0=      v(i  ,j  ,kx+1,3)
        v1= 0.5*(v(i+1,j  ,kx+1,3)+v(i+1,j  ,kx,3))
        v2= 0.5*(v(i-1,j  ,kx+1,3)+v(i-1,j  ,kx,3))
        v3= 0.5*(v(i  ,j+1,kx+1,3)+v(i  ,j+1,kx,3))
        v4= 0.5*(v(i  ,j-1,kx+1,3)+v(i  ,j-1,kx,3))
        
        b0 = real( ibits(bv(i  ,j  ,kx+1), State, 1) * ibits(bv(i  ,j  ,kx), State, 1) )
        b1 = real( ibits(bv(i+1,j  ,kx+1), State, 1) * ibits(bv(i+1,j  ,kx), State, 1) )
        b2 = real( ibits(bv(i-1,j  ,kx+1), State, 1) * ibits(bv(i-1,j  ,kx), State, 1) )
        b3 = real( ibits(bv(i  ,j+1,kx+1), State, 1) * ibits(bv(i  ,j+1,kx), State, 1) )
        b4 = real( ibits(bv(i  ,j-1,kx+1), State, 1) * ibits(bv(i  ,j-1,kx), State, 1) )
        
        r  = 2.0*w_ref - 0.5*(v0+c3)
        d1 = v1*b1 + (1.0-b1)*r
        d2 = v2*b2 + (1.0-b2)*r
        d3 = v3*b3 + (1.0-b3)*r
        d4 = v4*b4 + (1.0-b4)*r

        b5 = 1.0-b0
        uu = (c1 + 0.5*(d1-d2))*b0 + b5*u_ref
        vv = (c2 + 0.5*(d3-d4))*b0 + b5*v_ref
        ww =  c3               *b0 + b5*w_ref
                
        do k=kx+1, kx+g
          v(i,j,k,1)= uu
          v(i,j,k,2)= vv
          v(i,j,k,3)= ww
        end do
      end do
      end do
      
      flop = flop + real(ix*jx*49)
      
    case default
    end select FACES

    return 
    end subroutine cbc_vobc_tfree

!  *****************************************************************
!> @subroutine cbc_vibc_drchlt (v, sz, g, st, ed, v00, bv, odr, vec)
!! @brief 計算領域内部の速度指定境界条件を設定するために必要な参照値をセットする
!! @param[out] v 速度ベクトル（セルセンタ）
!! @param sz 配列長
!! @param g ガイドセル長
!! @param st ループの開始インデクス
!! @param ed ループの終了インデクス
!! @param v00 参照速度
!! @param bv BCindex V
!! @param odr 内部境界処理時の速度境界条件のエントリ
!! @param vec 指定する速度ベクトル
!<
    subroutine cbc_vibc_drchlt (v, sz, g, st, ed, v00, bv, odr, vec)
    implicit none
    include 'cbc_f_params.h'
    integer                                                     ::  i, j, k, g, idx, odr
    integer, dimension(3)                                       ::  sz, st, ed
    real                                                        ::  u_bc_ref, v_bc_ref, w_bc_ref
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)   ::  v
    real, dimension(0:3)                                        ::  v00
    real, dimension(3)                                          ::  vec
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bv
    
    ! u_bc_refは参照座標系での境界速度
    u_bc_ref = vec(1) + v00(1)
    v_bc_ref = vec(2) + v00(2)
    w_bc_ref = vec(3) + v00(3)
    
    do k=st(3),ed(3)
    do j=st(2),ed(2)
    do i=st(1),ed(1)
      idx = bv(i,j,k)

      if ( 0 /= iand(idx, bc_mask30) ) then ! 6面のうちのどれか速度境界フラグが立っている場合

        if ( ibits(idx, bc_face_W, bitw_5) == odr ) then
          v(i-1,j,k,1) = u_bc_ref
          v(i-1,j,k,2) = v_bc_ref
          v(i-1,j,k,3) = w_bc_ref
        end if

        if ( ibits(idx, bc_face_E, bitw_5) == odr ) then
          v(i+1,j,k,1) = u_bc_ref
          v(i+1,j,k,2) = v_bc_ref
          v(i+1,j,k,3) = w_bc_ref
        end if

        if ( ibits(idx, bc_face_S, bitw_5) == odr ) then
          v(i,j-1,k,1) = u_bc_ref
          v(i,j-1,k,2) = v_bc_ref
          v(i,j-1,k,3) = w_bc_ref
        end if
        
        if ( ibits(idx, bc_face_N, bitw_5) == odr ) then
          v(i,j+1,k,1) = u_bc_ref
          v(i,j+1,k,2) = v_bc_ref
          v(i,j+1,k,3) = w_bc_ref
        end if
        
        if ( ibits(idx, bc_face_B, bitw_5) == odr ) then
          v(i,j,k-1,1) = u_bc_ref
          v(i,j,k-1,2) = v_bc_ref
          v(i,j,k-1,3) = w_bc_ref
        end if
			
        if ( ibits(idx, bc_face_T, bitw_5) == odr ) then
          v(i,j,k+1,1) = u_bc_ref
          v(i,j,k+1,2) = v_bc_ref
          v(i,j,k+1,3) = w_bc_ref
        end if

      endif
    end do
    end do
    end do
    
    return
    end subroutine cbc_vibc_drchlt

!  **********************************************************
!> @subroutine cbc_vobc_drchlt (v, sz, g, v00, bv, face, vec)
!! @brief ガイドセルの速度指定境界条件を設定するために必要な参照値をセットする
!! @param[out] v セルセンタ速度
!! @param sz 配列長
!! @param g ガイドセル長
!! @param v00 参照速度
!! @param bv BCindex V
!! @param face 外部境界の面番号
!! @param vec 指定する速度ベクトル
!! @note 流束型の境界条件を用いるので，内点の計算に使う参照点に値があればよい（1層）
!<
    subroutine cbc_vobc_drchlt (v, sz, g, v00, bv, face, vec)
    implicit none
    include 'cbc_f_params.h'
    integer                                                     ::  i, j, k, g, face, ix, jx, kx, bvx
    integer, dimension(3)                                       ::  sz
    real                                                        ::  u_bc_ref, v_bc_ref, w_bc_ref
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)   ::  v
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bv
    real, dimension(3)                                          ::  vec
    real, dimension(0:3)                                        ::  v00

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    
    ! u_bc_refは参照座標系での境界速度
    u_bc_ref = vec(1) + v00(1)
    v_bc_ref = vec(2) + v00(2)
    w_bc_ref = vec(3) + v00(3)
    
    FACES : select case (face)
    case (X_minus)
      i = 1
      do k=1,kx
      do j=1,jx
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_W, bitw_5) == obc_mask ) then
          v(i-1,j,k,1) = u_bc_ref
          v(i-1,j,k,2) = v_bc_ref
          v(i-1,j,k,3) = w_bc_ref
        endif
      end do
      end do
      
    case (X_plus)
      i = ix
      do k=1,kx
      do j=1,jx
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_E, bitw_5) == obc_mask ) then
          v(i+1,j,k,1) = u_bc_ref
          v(i+1,j,k,2) = v_bc_ref
          v(i+1,j,k,3) = w_bc_ref
        endif
      end do
      end do
      
    case (Y_minus)
      j = 1
      do k=1,kx
      do i=1,ix
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_S, bitw_5) == obc_mask ) then
          v(i,j-1,k,1) = u_bc_ref
          v(i,j-1,k,2) = v_bc_ref
          v(i,j-1,k,3) = w_bc_ref
        endif
      end do
      end do
      
    case (Y_plus)
      j = jx
      do k=1,kx
      do i=1,ix
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_N, bitw_5) == obc_mask ) then
          v(i,j+1,k,1) = u_bc_ref
          v(i,j+1,k,2) = v_bc_ref
          v(i,j+1,k,3) = w_bc_ref
        endif
      end do
      end do
      
    case (Z_minus)
      k = 1
      do j=1,jx
      do i=1,ix
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_B, bitw_5) == obc_mask ) then
          v(i,j,k-1,1) = u_bc_ref
          v(i,j,k-1,2) = v_bc_ref
          v(i,j,k-1,3) = w_bc_ref
        endif
      end do
      end do
      
    case (Z_plus)
      k = kx
      do j=1,jx
      do i=1,ix
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_T, bitw_5) == obc_mask ) then
          v(i,j,k+1,1) = u_bc_ref
          v(i,j,k+1,2) = v_bc_ref
          v(i,j,k+1,3) = w_bc_ref
        endif
      end do
      end do
      
    case default
    end select FACES
    
    return
    end subroutine cbc_vobc_drchlt
    
!  ************************************************************
!> @subroutine cbc_vibc_outflow (v, sz, g, st, ed, cf, bv, odr)
!! @brief 速度指定境界条件を設定するために必要な参照値を，流出下流のセル（固体セル）にセットする
!! @param[out] v 速度 u^{n+1}
!! @param sz 配列長
!! @param g ガイドセル長
!! @param st ループの開始インデクス
!! @param ed ループの終了インデクス
!! @param bv BCindex V
!! @param odr 内部境界処理時の速度境界条件のエントリ
!! @note 参照される頻度の少ない（逆流時の）近似値として，ゼロ時近似を与える．発散しないようにする消極的な代入．
!<
    subroutine cbc_vibc_outflow (v, sz, g, st, ed, bv, odr)
    implicit none
    include 'cbc_f_params.h'
    integer                                                     ::  i, j, k, g, idx, odr
    integer, dimension(3)                                       ::  sz, st, ed
    real                                                        ::  Up, Vp, Wp
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)   ::  v
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bv
    
    do k=st(3),ed(3)
    do j=st(2),ed(2)
    do i=st(1),ed(1)
    
      idx = bv(i,j,k)
      if ( 0 /= iand(idx, bc_mask30) ) then ! 6面のうちのどれか速度境界フラグが立っている場合

        Up = v(i,j,k,1)
        Vp = v(i,j,k,2)
        Wp = v(i,j,k,3)
        
        ! X-
        if ( ibits(idx, bc_face_W, bitw_5) == odr ) then
          v(i-1,j,k,1) = Up
          v(i-1,j,k,2) = Vp
          v(i-1,j,k,3) = Wp
        end if
        
        ! X+
        if ( ibits(idx, bc_face_E, bitw_5) == odr ) then
          v(i+1,j,k,1) = Up
          v(i+1,j,k,2) = Vp
          v(i+1,j,k,3) = Wp
        end if
        
        ! Y-
        if ( ibits(idx, bc_face_S, bitw_5) == odr ) then
          v(i,j-1,k,1) = Up
          v(i,j-1,k,2) = Vp
          v(i,j-1,k,3) = Wp
        end if
        
        ! Y+
        if ( ibits(idx, bc_face_N, bitw_5) == odr ) then
          v(i,j+1,k,1) = Up
          v(i,j+1,k,2) = Vp
          v(i,j+1,k,3) = Wp
        end if
        
        ! Z-
        if ( ibits(idx, bc_face_B, bitw_5) == odr ) then
          v(i,j,k-1,1) = Up
          v(i,j,k-1,2) = Vp
          v(i,j,k-1,3) = Wp
        end if
			  
        ! Z+
        if ( ibits(idx, bc_face_T, bitw_5) == odr ) then
          v(i,j,k+1,1) = Up
          v(i,j,k+1,2) = Vp
          v(i,j,k+1,3) = Wp
        end if
        
      endif
    end do
    end do
    end do
    
    return
    end subroutine cbc_vibc_outflow

!  *****************************************************************************************
!> @subroutine cbc_div_ibc_oflow_pvec (div, sz, g, st, ed, v00, cf, coef, bv, odr, v0, flop)
!! @brief 内部流出境界条件による疑似速度ベクトルの発散の修正
!! @param[out] div 速度の発散
!! @param sz 配列長
!! @param g ガイドセル長
!! @param st ループの開始インデクス
!! @param ed ループの終了インデクス
!! @param v00 参照速度
!! @param cf u_out*dt/dh
!! @param coef 係数
!! @param bv BCindex V
!! @param odr 速度境界条件のエントリ
!! @param v0 セルセンター速度　u^n
!! @param[out] flop flop count
!! @note 流出境界面ではu_e^{n+1}=u_e^n-cf*(u_e^n-u_w^n)を予測値としてdivの寄与として加算．u_e^nの値は連続の式から計算する．
!<
    subroutine cbc_div_ibc_oflow_pvec (div, sz, g, st, ed, v00, cf, coef, bv, odr, v0, flop)
    implicit none
    include 'cbc_f_params.h'
    integer                                                     ::  i, j, k, g, bvx, odr
    integer, dimension(3)                                       ::  sz, st, ed
    real                                                        ::  flop, coef
    real                                                        ::  b_w, b_e, b_s, b_n, b_b, b_t
    real                                                        ::  Ue, Uw, Vn, Vs, Wt, Wb
    real                                                        ::  Ue_t, Uw_t, Vn_t, Vs_t, Wt_t, Wb_t
    real                                                        ::  cf, u_ref, v_ref, w_ref
    real                                                        ::  Up0, Ue0, Uw0, Vp0, Vs0, Vn0, Wp0, Wb0, Wt0
    real, dimension(0:3)                                        ::  v00
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  div
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)   ::  vc, v0
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bv

    ! 参照速度
    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)
    
    do k=st(3),ed(3)
    do j=st(2),ed(2)
    do i=st(1),ed(1)
      bvx = bv(i,j,k)
      if ( 0 /= iand(bvx, bc_mask30) ) then ! 6面のうちのどれか速度境界フラグが立っている場合

        include 'd_o_o_p.h'
        
        ! 寄与をゼロにしておく
        Ue_t = 0.0
        Uw_t = 0.0
        Vn_t = 0.0
        Vs_t = 0.0
        Wt_t = 0.0
        Wb_t = 0.0
        
        ! X方向 ---------------------------------------
        if ( ibits(bvx, bc_face_W, bitw_5) == odr ) then
          Uw = Ue + (Vn - Vs + Wt - Wb) ! 連続の式から流出面の速度を推定，これは移動座標系上の速度成分
          if ( cf>0.0 ) cf=0.0
          Uw_t = Uw - cf*(Ue-Uw)
          flop = flop + 7.0
        endif
        
        if ( ibits(bvx, bc_face_E, bitw_5) == odr ) then
          Ue = Uw - (Vn - Vs + Wt - Wb)
          if ( cf<0.0 ) cf=0.0
          Ue_t = Ue - cf*(Ue-Uw)
          flop = flop + 7.0
        endif
        
        ! Y方向 ---------------------------------------
        if ( ibits(bvx, bc_face_S, bitw_5) == odr ) then
          Vs = Vn + (Ue - Uw + Wt - Wb)
          if ( cf>0.0 ) cf=0.0
          Vs_t = Vs - cf*(Vn-Vs)
          flop = flop + 7.0
        endif
        
        if ( ibits(bvx, bc_face_N, bitw_5) == odr ) then
          Vn = Vs - (Ue - Uw + Wt - Wb)
          if ( cf<0.0 ) cf=0.0
          Vn_t = Vn - cf*(Vn-Vs)
          flop = flop + 7.0
        endif
        
        ! Z方向 ---------------------------------------
        if ( ibits(bvx, bc_face_B, bitw_5) == odr ) then
          Wb = Wt + (Ue - Uw + Vn - Vs)
          if ( cf>0.0 ) cf=0.0
          Wb_t = Wb - cf*(Wt-Wb)
          flop = flop + 7.0
        endif
        
        if ( ibits(bvx, bc_face_T, bitw_5) == odr ) then
          Wt = Wb - (Ue - Uw + Vn - Vs)
          if ( cf<0.0 ) cf=0.0
          Wt_t = Wt - cf*(Wt-Wb)
          flop = flop + 7.0
        endif

        ! VBCの面だけUe_tなどは値をもつ
        div(i,j,k) = div(i,j,k) + ( Ue_t - Uw_t + Vn_t - Vs_t + Wt_t - Wb_t ) * coef ! 対象セルは流体なのでマスク不要
        flop = flop + 7.0
      end if
    end do
    end do
    end do

    return
    end subroutine cbc_div_ibc_oflow_pvec
    
!  *************************************************************************************
!> @subroutine cbc_div_ibc_oflow_vec (div, sz, g, st, ed, v00, coef, bv, odr, avr, flop)
!! @brief 内部流出境界条件によるn+1時刻の速度の発散の修正と流量の積算
!! @param[in/out] div 速度の発散
!! @param sz 配列長
!! @param g ガイドセル長
!! @param st ループの開始インデクス
!! @param ed ループの終了インデクス
!! @param v00 参照速度
!! @param coef 係数 h/dt
!! @param bv BCindex V
!! @param odr 速度境界条件のエントリ
!! @param[out] avr 平均流出速度（面直方向のみ）
!! @param[out] flop flop count
!! @note div(u)=0から，内部流出境界のセルで計算されたdivの値が流出速度となる
!<
    subroutine cbc_div_ibc_oflow_vec (div, sz, g, st, ed, v00, coef, bv, odr, avr, flop)
    implicit none
    include 'cbc_f_params.h'
    include 'sklparaf.h'
    integer                                                     ::  i, j, k, g, bvx, m, odr, iret, ierr
    integer, dimension(3)                                       ::  sz, st, ed
    real                                                        ::  flop, coef, rc
    real                                                        ::  u_ref, v_ref, w_ref, dv, a1, avr
    real, dimension(0:3)                                        ::  v00
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  div
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bv
    real, dimension(2)                                          ::  av, tmp
    
    ! 参照速度
    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)
    a1 = 0.0
    m = 0
    rc = 1.0/coef
    
    do k=st(3),ed(3)
    do j=st(2),ed(2)
    do i=st(1),ed(1)
      bvx = bv(i,j,k)
      if ( 0 /= iand(bvx, bc_mask30) ) then ! 6面のうちのどれか速度境界フラグが立っている場合
        dv = div(i,j,k) * rc
        flop = flop + 1.0
        
        if ( ibits(bvx, bc_face_W, bitw_5) == odr ) then ! u_w
          a1 = a1 + dv - u_ref
          m = m+1
        endif
        
        if ( ibits(bvx, bc_face_E, bitw_5) == odr ) then ! u_e
          a1 = a1 - dv - u_ref
          m = m+1
        endif
        
        if ( ibits(bvx, bc_face_S, bitw_5) == odr ) then
          a1 = a1 + dv - v_ref
          m = m+1
        endif
        
        if ( ibits(bvx, bc_face_N, bitw_5) == odr ) then
          a1 = a1 - dv - v_ref
          m = m+1
        endif
        
        if ( ibits(bvx, bc_face_B, bitw_5) == odr ) then
          a1 = a1 + dv - w_ref
          m = m+1
        endif
        
        if ( ibits(bvx, bc_face_T, bitw_5) == odr ) then
          a1 = a1 - dv - w_ref
          m = m+1
        endif

        div(i,j,k) = 0.0 ! 対象セルは発散をゼロにする
      end if
    end do
    end do
    end do

    av(1) = a1
    av(2) = real(m)
    flop = flop + real(m*2)

    call SklIsParallel(iret)
    if ( iret == 1 ) then
      tmp(1) = av(1)
      tmp(2) = av(2)
      call SklAllreduce(tmp, av, 2, SKL_REAL, SKL_SUM, SKL_DEFAULT_GROUP, ierr)
    end if

    avr = av(1)/av(2)

    return
    end subroutine cbc_div_ibc_oflow_vec

!  **********************************************************************************
!> @subroutine cbc_div_ibc_drchlt (div, sz, g, st, ed, v00, coef, bv, odr, vec, flop)
!! @brief 内部速度指定境界条件による疑似速度の発散の修正
!! @param[out] div 速度の発散
!! @param sz 配列長
!! @param g ガイドセル長
!! @param st ループの開始インデクス
!! @param ed ループの終了インデクス
!! @param v00 参照速度
!! @param coef 係数
!! @param bv BCindex V
!! @param odr 速度境界条件のエントリ
!! @param vec 指定する速度ベクトル
!! @param[out] flop flop count
!<
    subroutine cbc_div_ibc_drchlt (div, sz, g, st, ed, v00, coef, bv, odr, vec, flop)
    implicit none
    include 'cbc_f_params.h'
    integer                                                     ::  i, j, k, g, bvx, m, odr
    integer, dimension(3)                                       ::  sz, st, ed
    real                                                        ::  flop, coef
    real                                                        ::  Ue_t, Uw_t, Vn_t, Vs_t, Wt_t, Wb_t
    real                                                        ::  u_bc_ref, v_bc_ref, w_bc_ref
    real, dimension(3)                                          ::  vec
    real, dimension(0:3)                                        ::  v00
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  div
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bv
    
    ! u_bc_refは参照座標系での境界速度
    u_bc_ref = vec(1) + v00(1)
    v_bc_ref = vec(2) + v00(2)
    w_bc_ref = vec(3) + v00(3)
    
    m = 0
    
    do k=st(3),ed(3)
    do j=st(2),ed(2)
    do i=st(1),ed(1)
      bvx = bv(i,j,k)
      if ( 0 /= iand(bvx, bc_mask30) ) then ! 6面のうちのどれか速度境界フラグが立っている場合
        Ue_t = 0.0
        Uw_t = 0.0
        Vn_t = 0.0
        Vs_t = 0.0
        Wt_t = 0.0
        Wb_t = 0.0
        
        if ( ibits(bvx, bc_face_W, bitw_5) == odr ) Uw_t = u_bc_ref
        if ( ibits(bvx, bc_face_E, bitw_5) == odr ) Ue_t = u_bc_ref
        if ( ibits(bvx, bc_face_S, bitw_5) == odr ) Vs_t = v_bc_ref
        if ( ibits(bvx, bc_face_N, bitw_5) == odr ) Vn_t = v_bc_ref
        if ( ibits(bvx, bc_face_B, bitw_5) == odr ) Wb_t = w_bc_ref
        if ( ibits(bvx, bc_face_T, bitw_5) == odr ) Wt_t = w_bc_ref

        ! VBCの面だけUe_tなどは値をもつ
        div(i,j,k) = div(i,j,k) + ( Ue_t - Uw_t + Vn_t - Vs_t + Wt_t - Wb_t ) * coef ! 対象セルは流体なのでマスク不要
        m = m+1
      end if
    end do
    end do
    end do
    flop = flop + real(m*7)

    return
    end subroutine cbc_div_ibc_drchlt

!  ***************************************************************************
!> @subroutine cbc_div_obc_drchlt (div, sz, g, face, v00, coef, bv, vec, flop)
!! @brief 外部指定境界条件による速度の発散の修正
!! @param[out] div 速度の発散
!! @param sz 配列長
!! @param g ガイドセル長
!! @param face 面番号
!! @param v00 参照速度
!! @param coef 係数
!! @param bv BCindex V
!! @param vec 指定する速度ベクトル
!! @param[out] flop flop count
!! @note 指定面でも固体部分は対象外とするのでループ中に判定あり
!<
    subroutine cbc_div_obc_drchlt (div, sz, g, face, v00, coef, bv, vec, flop)
    implicit none
    include 'cbc_f_params.h'
    integer                                                     ::  i, j, k, g, ix, jx, kx, face, bvx
    integer, dimension(3)                                       ::  sz
    real                                                        ::  flop, coef
    real                                                        ::  u_bc_ref, v_bc_ref, w_bc_ref
    real, dimension(0:3)                                        ::  v00
    real, dimension(3)                                          ::  vec
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  div
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bv

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    
    ! 参照座標系の速度に係数をかけておく
    u_bc_ref = (vec(1) + v00(1)) * coef
    v_bc_ref = (vec(2) + v00(2)) * coef
    w_bc_ref = (vec(3) + v00(3)) * coef
    
    flop = flop + 6.0

    FACES : select case (face)
    case (X_minus)
      i = 1
      do k=1,kx
      do j=1,jx
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_W, bitw_5) == obc_mask ) then
          div(i,j,k) = div(i,j,k) - u_bc_ref * real(ibits(bvx, State, 1))
          flop = flop + 2.0
        endif
      end do
      end do
      
    case (X_plus)
      i = ix
      do k=1,kx
      do j=1,jx
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_E, bitw_5) == obc_mask ) then
          div(i,j,k) = div(i,j,k) + u_bc_ref * real(ibits(bvx, State, 1))
          flop = flop + 2.0
        endif
      end do
      end do
      
    case (Y_minus)
      j = 1
      do k=1,kx
      do i=1,ix
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_S, bitw_5) == obc_mask ) then
          div(i,j,k) = div(i,j,k) - v_bc_ref * real(ibits(bvx, State, 1))
          flop = flop + 2.0
        endif
      end do
      end do
      
    case (Y_plus)
      j = jx
      do k=1,kx
      do i=1,ix
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_N, bitw_5) == obc_mask ) then
          div(i,j,k) = div(i,j,k) + v_bc_ref * real(ibits(bvx, State, 1))
          flop = flop + 2.0
        endif
      end do
      end do
    
    case (Z_minus)
      k = 1
      do j=1,jx
      do i=1,ix
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_B, bitw_5) == obc_mask ) then
          div(i,j,k) = div(i,j,k) - w_bc_ref * real(ibits(bvx, State, 1))
          flop = flop + 2.0
        endif
      end do
      end do
      
    case (Z_plus)
      k = kx
      do j=1,jx
      do i=1,ix
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_T, bitw_5) == obc_mask ) then
          div(i,j,k) = div(i,j,k) + w_bc_ref * real(ibits(bvx, State, 1))
          flop = flop + 2.0
        endif
      end do
      end do
    
    case default
    end select FACES

    return
    end subroutine cbc_div_obc_drchlt

!  *************************************************************************************
!> @subroutine cbc_div_obc_oflow_pvec (div, sz, g, face, v00, v_out, coef, bv, v0, flop)
!! @brief 外部流出境界条件による疑似速度ベクトルの発散の修正
!! @param[out] div 速度の発散
!! @param sz 配列長
!! @param g ガイドセル長
!! @param face 面番号
!! @param v00 参照速度
!! @param v_out u_out*dt/dh
!! @param coef 係数 dh/dt
!! @param bv BCindex V
!! @param v0 速度ベクトル u^n
!! @param[out] flop flop count
!! @note 指定面でも固体部分は対象外とするのでループ中に判定あり
!<
    subroutine cbc_div_obc_oflow_pvec (div, sz, g, face, v00, v_out, coef, bv, v0, flop)
    implicit none
    include 'cbc_f_params.h'
    integer                                                     ::  i, j, k, g, ix, jx, kx, face, bvx
    integer, dimension(3)                                       ::  sz
    real                                                        ::  flop, coef, v_out
    real                                                        ::  b_w, b_e, b_s, b_n, b_b, b_t
    real                                                        ::  Up0, Ue0, Uw0, Vp0, Vs0, Vn0, Wp0, Wb0, Wt0
    real                                                        ::  Ue, Uw, Vn, Vs, Wt, Wb
    real                                                        ::  Ue_t, Uw_t, Vn_t, Vs_t, Wt_t, Wb_t
    real                                                        ::  u_ref, v_ref, w_ref
    real, dimension(0:3)                                        ::  v00
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  div
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bv
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)   ::  v0

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    
    ! 参照速度
    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)

    FACES : select case (face)
    case (X_minus)
      if ( v_out>0.0 ) v_out=0.0
      i = 1
      do k=1,kx
      do j=1,jx
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_W, bitw_5) == obc_mask ) then
        
          include 'd_o_o_p.h'
          
          Uw = Ue + (Vn - Vs + Wt - Wb) ! 連続の式から流出面の速度を推定，これは移動座標系上の速度成分
          Uw_t = Uw - v_out*(Ue-Uw)
          div(i,j,k) = div(i,j,k) - Uw_t * coef
          flop = flop + 9.0
        endif
      end do
      end do
      
    case (X_plus)
      if ( v_out<0.0 ) v_out=0.0
      i = ix
      do k=1,kx
      do j=1,jx
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_E, bitw_5) == obc_mask ) then
        
          include 'd_o_o_p.h'
          
          Ue = Uw - (Vn - Vs + Wt - Wb)
          Ue_t = Ue - v_out*(Ue-Uw)
          div(i,j,k) = div(i,j,k) + Ue_t * coef
          flop = flop + 9.0
        endif
      end do
      end do
      
    case (Y_minus)
      if ( v_out>0.0 ) v_out=0.0
      j = 1
      do k=1,kx
      do i=1,ix
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_S, bitw_5) == obc_mask ) then
        
          include 'd_o_o_p.h'
        
          Vs = Vn + (Ue - Uw + Wt - Wb)
          Vs_t = Vs - v_out*(Vn-Vs)
          div(i,j,k) = div(i,j,k) - Vs_t * coef
          flop = flop + 9.0
        endif
      end do
      end do
      
    case (Y_plus)
      if ( v_out<0.0 ) v_out=0.0
      j = jx
      do k=1,kx
      do i=1,ix
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_N, bitw_5) == obc_mask ) then
        
          include 'd_o_o_p.h'
        
          Vn = Vs - (Ue - Uw + Wt - Wb)
          Vn_t = Vn - v_out*(Vn-Vs)
          div(i,j,k) = div(i,j,k) + Vn_t * coef
          flop = flop + 9.0
        endif
      end do
      end do
    
    case (Z_minus)
      if ( v_out>0.0 ) v_out=0.0
      k = 1
      do j=1,jx
      do i=1,ix
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_B, bitw_5) == obc_mask ) then
        
          include 'd_o_o_p.h'
        
          Wb = Wt + (Ue - Uw + Vn - Vs)
          Wb_t = Wb - v_out*(Wt-Wb)
          div(i,j,k) = div(i,j,k) - Wb_t * coef
          flop = flop + 9.0
        endif
      end do
      end do
      
    case (Z_plus)
      if ( v_out<0.0 ) v_out=0.0
      k = kx
      do j=1,jx
      do i=1,ix
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_T, bitw_5) == obc_mask ) then
        
          include 'd_o_o_p.h'
        
          Wt = Wb - (Ue - Uw + Vn - Vs)
          Wt_t = Wt - v_out*(Wt-Wb)
          div(i,j,k) = div(i,j,k) + Wt_t * coef
          flop = flop + 9.0
        endif
      end do
      end do
    
    case default
    end select FACES

    return
    end subroutine cbc_div_obc_oflow_pvec

!  *****************************************************************************
!> @subroutine cbc_div_obc_oflow_vec (div, sz, g, face, v00, coef, bv, aa, flop)
!! @brief 外部流出境界条件による速度ベクトルの発散の修正
!! @param[out] div 速度の発散
!! @param sz 配列長
!! @param g ガイドセル長
!! @param face 面番号
!! @param v00 参照速度
!! @param coef 係数
!! @param bv BCindex V
!! @param[out] aa 領域境界の積算値
!! @param[out] flop flop count
!! @note 指定面でも固体部分は対象外とするのでループ中に判定あり
!!       div(u)=0から，内部流出境界のセルで計算されたdivが流出速度となる
!<
    subroutine cbc_div_obc_oflow_vec (div, sz, g, face, v00, coef, bv, aa, flop)
    implicit none
    include 'cbc_f_params.h'
    integer                                                     ::  i, j, k, g, ix, jx, kx, face, bvx
    integer, dimension(3)                                       ::  sz
    real                                                        ::  flop, coef, dv, a1, a2, a3, u_ref, v_ref, w_ref, rc
    real, dimension(0:3)                                        ::  v00
    real, dimension(3)                                          ::  aa
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  div
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bv

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)

    a1 = 0.0   ! sum
    a2 = 1.0e6 ! min
    a3 =-1.0e6 ! max
    
    rc = 1.0/coef
    
    ! 参照速度
    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)
    
    FACES : select case (face)
    case (X_minus)
      i = 1
      do k=1,kx
      do j=1,jx
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_W, bitw_5) == obc_mask ) then
          dv = div(i,j,k) * rc - u_ref
          a1 = a1 + dv
          a2 = min(a2, dv)
          a3 = max(a3, dv)
          div(i,j,k) = 0.0 ! 対象セルは発散をゼロにする
          flop = flop + 5.0
        endif
      end do
      end do
      
    case (X_plus)
      i = ix
      do k=1,kx
      do j=1,jx
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_E, bitw_5) == obc_mask ) then
          dv = -div(i,j,k) * rc  - u_ref
          a1 = a1 + dv
          a2 = min(a2, dv)
          a3 = max(a3, dv)
          div(i,j,k) = 0.0
          flop = flop + 5.0
        endif
      end do
      end do
      
    case (Y_minus)
      j = 1
      do k=1,kx
      do i=1,ix
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_S, bitw_5) == obc_mask ) then
          dv = div(i,j,k) * rc - v_ref
          a1 = a1 + dv
          a2 = min(a2, dv)
          a3 = max(a3, dv)
          div(i,j,k) = 0.0
          flop = flop + 5.0
        endif
      end do
      end do
      
    case (Y_plus)
      j = jx
      do k=1,kx
      do i=1,ix
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_N, bitw_5) == obc_mask ) then
          dv = -div(i,j,k) * rc - v_ref
          a1 = a1 + dv
          a2 = min(a2, dv)
          a3 = max(a3, dv)
          div(i,j,k) = 0.0
          flop = flop + 5.0
        endif
      end do
      end do
    
    case (Z_minus)
      k = 1
      do j=1,jx
      do i=1,ix
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_B, bitw_5) == obc_mask ) then
          dv = div(i,j,k) * rc - w_ref
          a1 = a1 + dv
          a2 = min(a2, dv)
          a3 = max(a3, dv)
          div(i,j,k) = 0.0
          flop = flop + 5.0
        endif
      end do
      end do
      
    case (Z_plus)
      k = kx
      do j=1,jx
      do i=1,ix
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_T, bitw_5) == obc_mask ) then
          dv = -div(i,j,k) * rc - w_ref
          a1 = a1 + dv
          a2 = min(a2, dv)
          a3 = max(a3, dv)
          div(i,j,k) = 0.0
          flop = flop + 5.0
        endif
      end do
      end do
    
    case default
    end select FACES
    
    aa(1) = a1
    aa(2) = a2
    aa(3) = a3

    return
    end subroutine cbc_div_obc_oflow_vec
