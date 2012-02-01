!   *********************************************************
!
!   SPHERE - Skeleton for PHysical and Engineering REsearch
!  
!   Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
!
!   *********************************************************

!> @file cbc_3d.f90
!> @brief subroutines for CBC
!> @author keno, FSI Team, VCAD, RIKEN

!  ******************************************************************************************************
!> @subroutine cbc_pvec_muscl (wv, sz, g, dh, c_scheme, v00, rei, v, bv, bp, v_mode, ut, wall_type, flop)
!! @brief 対流項と粘性項の計算
!! @param[out] wv 疑似ベクトルの空間項の評価値
!! @param sz 配列長
!! @param g ガイドセル長
!! @param dh 格子幅
!! @param c_scheme 対流項スキームのモード（1-UWD, 2-center, 3-MUSCL）
!! @param v00 参照速度
!! @param rei レイノルズ数の逆数
!! @param v 速度ベクトル（n-step, collocated）
!! @param bv BCindex V
!! @param bp BCindex P
!! @param v_mode 粘性項のモード (0=対流項のみ, 1=対流項と粘性項，2=粘性項は壁法則)
!! @param ut 摩擦速度
!! @param wall_type 壁面条件 (0=no_slip, 1=slip)
!! @param[out] flop
!<
    subroutine cbc_pvec_muscl (wv, sz, g, dh, c_scheme, v00, rei, v, bv, bp, v_mode, ut, wall_type, flop)
    implicit none
    include 'cbc_f_params.h'
    integer                                                     ::  i, j, k, ix, jx, kx, g, c_scheme, bvx, v_mode, bpx, wall_type
    integer                                                     ::  b_e1, b_w1, b_n1, b_s1, b_t1, b_b1, b_e2, b_w2, b_n2, b_s2, b_t2, b_b2, b_p
    integer, dimension(3)                                       ::  sz
    real                                                        ::  UPe, UPw, VPn, VPs, WPt, WPb, u1, u2, u3, ug, e1, e2, e3, u_tau
    real                                                        ::  Up0, Ue1, Ue2, Uw1, Uw2, Us1, Us2, Un1, Un2, Ub1, Ub2, Ut1, Ut2
    real                                                        ::  Vp0, Ve1, Ve2, Vw1, Vw2, Vs1, Vs2, Vn1, Vn2, Vb1, Vb2, Vt1, Vt2
    real                                                        ::  Wp0, We1, We2, Ww1, Ww2, Ws1, Ws2, Wn1, Wn2, Wb1, Wb2, Wt1, Wt2
    real                                                        ::  ck, u_ref, v_ref, w_ref, dh, dh1, dh2, flop, vcs
    real                                                        ::  c_e, c_w, c_n, c_s, c_t, c_b, wls, wm1, wm2
    real                                                        ::  w_e, w_w, w_n, w_s, w_t, w_b
    real                                                        ::  uu_e, uu_w, uu_s, uu_n, uu_b, uu_t
    real                                                        ::  vv_e, vv_w, vv_s, vv_n, vv_b, vv_t
    real                                                        ::  ww_e, ww_w, ww_s, ww_n, ww_b, ww_t
    real                                                        ::  d1, d2, d3, d4, g1, g2, g3, g4, g5, g6, s1, s2, s3, s4, b
    real                                                        ::  Urr, Url, Ulr, Ull, Vrr, Vrl, Vlr, Vll, Wrr, Wrl, Wlr, Wll
    real                                                        ::  cr, cl, acr, acl, cnv_u, cnv_v, cnv_w, EX, EY, EZ, rei, beta, qtz
    real                                                        ::  fu_r, fu_l, fv_r, fv_l, fw_r, fw_l, uq, vq, wq, ss
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)   ::  v, wv
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  ut
    real, dimension(0:3)                                        ::  v00
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bv, bp
    
    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    
    dh1= 1.0/dh
    dh2= rei*dh1
    qtz = 1.0/255.0
    
    ! 粘性項のマスク
    if ( v_mode == 0 ) then ! 粘性項は計算しない
      vcs = 0.0
    else
      vcs = 1.0
    endif
    
    ! 参照座標系の速度
    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)
    
    ck = 0.0
    b  = 0.0
    ss = 1.0
            
    if ( c_scheme == 1 ) then      !     1st order upwind
      ss = 0.0
    else if ( c_scheme == 2 ) then !     2nd order central 
      ck = 1.0     
    else if ( c_scheme == 3 ) then !     3rd order MUSCL
      ck = 1.0/3.0
      b  = (3.0-ck)/(1.0-ck)
    else
      write(*,*) 'out of scheme selection'
      stop
    endif
    
    ! 壁面条件
    if ( wall_type == 0 ) then ! no slip
      wls = 0.0
    else if ( wall_type == 1 ) then ! slip
      wls = 1.0
    endif
    
    wm1 = 2.0*(1.0-wls)
    wm2 = 2.0*wls-1.0
    
    flop = flop + 6.0
    
    do k=1,kx
    do j=1,jx
    do i=1,ix
      cnv_u = 0.0
      cnv_v = 0.0
      cnv_w = 0.0

      ! 変数のロード
      include 'load_var_stencil5.h'
      
      bvx = bv(i,j,k)
      bpx = bp(i,j,k)
      
      ! セル状態 (0-solid / 1-fluid)
      b_p = ibits(bvx,             State, 1)
      b_w2= ibits(bv(i-2,j  ,k  ), State, 1)
      b_w1= ibits(bv(i-1,j  ,k  ), State, 1)
      b_e1= ibits(bv(i+1,j  ,k  ), State, 1)
      b_e2= ibits(bv(i+2,j  ,k  ), State, 1)
      b_s2= ibits(bv(i  ,j-2,k  ), State, 1)
      b_s1= ibits(bv(i  ,j-1,k  ), State, 1)
      b_n1= ibits(bv(i  ,j+1,k  ), State, 1)
      b_n2= ibits(bv(i  ,j+2,k  ), State, 1)
      b_b2= ibits(bv(i  ,j  ,k-2), State, 1)
      b_b1= ibits(bv(i  ,j  ,k-1), State, 1)
      b_t1= ibits(bv(i  ,j  ,k+1), State, 1)
      b_t2= ibits(bv(i  ,j  ,k+2), State, 1)
      
      ! セル界面のフラグ (0-wall face / 1-fluid)
      w_e = real(b_e1 * b_p)
      w_w = real(b_w1 * b_p)
      w_n = real(b_n1 * b_p)
      w_s = real(b_s1 * b_p)
      w_t = real(b_t1 * b_p)
      w_b = real(b_b1 * b_p)
      
      ! 各面のVBCフラグ ibits() = 0(Normal) / others(BC) >> c_e = 1.0(Normal) / 0.0(BC) 
      c_e = 1.0
      c_w = 1.0
      c_n = 1.0
      c_s = 1.0
      c_t = 1.0
      c_b = 1.0
      if ( ibits(bvx, bc_face_E, bitw_5) /= 0 ) c_e = 0.0
      if ( ibits(bvx, bc_face_W, bitw_5) /= 0 ) c_w = 0.0
      if ( ibits(bvx, bc_face_N, bitw_5) /= 0 ) c_n = 0.0
      if ( ibits(bvx, bc_face_S, bitw_5) /= 0 ) c_s = 0.0
      if ( ibits(bvx, bc_face_T, bitw_5) /= 0 ) c_t = 0.0
      if ( ibits(bvx, bc_face_B, bitw_5) /= 0 ) c_b = 0.0
      
      ! 界面速度（スタガード位置）
      UPe = 0.5*(Up0+Ue1)*w_e + u_ref*(1.0-w_e)
      UPw = 0.5*(Up0+Uw1)*w_w + u_ref*(1.0-w_w)
      VPn = 0.5*(Vp0+Vn1)*w_n + v_ref*(1.0-w_n)
      VPs = 0.5*(Vp0+Vs1)*w_s + v_ref*(1.0-w_s)
      WPt = 0.5*(Wp0+Wt1)*w_t + w_ref*(1.0-w_t)
      WPb = 0.5*(Wp0+Wb1)*w_b + w_ref*(1.0-w_b)
      
      ! セルセンターからの壁面修正速度
      uq = 2.0*u_ref - Up0
      vq = 2.0*v_ref - Vp0
      wq = 2.0*w_ref - Wp0
      flop = flop + 48.0
			
      ! X方向 ---------------------------------------
      
      ! 速度指定の場合にMUSCLスキームの参照先として，固体内にテンポラリに与えた値を使う
      if ( (b_e2 == 0)  ) then
        Ue2 = 2.0*u_ref - v(i+1,j  ,k  ,1)
        Ve2 = wm1*v_ref + v(i+1,j  ,k  ,2)*wm2
        We2 = wm1*w_ref + v(i+1,j  ,k  ,3)*wm2
        flop = flop + 8.0
        !Ve2 = 2.0*v_ref - v(i+1,j  ,k  ,2)
        !We2 = 2.0*w_ref - v(i+1,j  ,k  ,3)
      endif
      
      ! 壁面の場合の参照速度の修正
      if ( b_e1 == 0 ) then
        Ue1 = uq
        Ve1 = wm1*v_ref + Vp0*wm2
        We1 = wm1*w_ref + Wp0*wm2
        flop = flop + 6.0
        !Ve1 = vq
        !We1 = wq
      endif
      
      if ( b_w1 == 0 ) then
        Uw1 = uq
        Vw1 = wm1*v_ref + Vp0*wm2
        Ww1 = wm1*w_ref + Wp0*wm2
        flop = flop + 6.0
        !Vw1 = vq
        !Ww1 = wq
      end if
      
      if ( (b_w2 == 0)  ) then
        Uw2 = 2.0*u_ref - v(i-1,j  ,k  ,1)
        Vw2 = wm1*v_ref + v(i-1,j  ,k  ,2)*wm2
        Ww2 = wm1*w_ref + v(i-1,j  ,k  ,3)*wm2
        flop = flop + 8.0
        !Vw2 = 2.0*v_ref - v(i-1,j  ,k  ,2)
        !Ww2 = 2.0*w_ref - v(i-1,j  ,k  ,3)
      end if
      
      ! 流束　流体のみと固体壁の影響を含み，隣接セルが固体の場合にはマスクする
      cr  = UPe - u_ref
      cl  = UPw - u_ref
      acr = abs(cr)
      acl = abs(cl)
      
      d4 = Ue2-Ue1
      d3 = Ue1-Up0
      d2 = Up0-Uw1
      d1 = Uw1-Uw2
      flop = flop + 10.0
      
      include 'muscl.h'
      
      Urr = Ue1 - 0.25*((1-ck)*g6+(1+ck)*g5)*ss
      Url = Up0 + 0.25*((1-ck)*g3+(1+ck)*g4)*ss
      Ulr = Up0 - 0.25*((1-ck)*g4+(1+ck)*g3)*ss
      Ull = Uw1 + 0.25*((1-ck)*g1+(1+ck)*g2)*ss
      fu_r = 0.5*(cr*(Urr+Url) - acr*(Urr-Url)) * w_e
      fu_l = 0.5*(cl*(Ulr+Ull) - acl*(Ulr-Ull)) * w_w

      d4 = Ve2-Ve1
      d3 = Ve1-Vp0
      d2 = Vp0-Vw1
      d1 = Vw1-Vw2
      flop = flop + 50.0
      
      include 'muscl.h'
      
      Vrr = Ve1 - 0.25*((1-ck)*g6+(1+ck)*g5)*ss
      Vrl = Vp0 + 0.25*((1-ck)*g3+(1+ck)*g4)*ss
      Vlr = Vp0 - 0.25*((1-ck)*g4+(1+ck)*g3)*ss
      Vll = Vw1 + 0.25*((1-ck)*g1+(1+ck)*g2)*ss
      fv_r = 0.5*(cr*(Vrr+Vrl) - acr*(Vrr-Vrl)) * w_e
      fv_l = 0.5*(cl*(Vlr+Vll) - acl*(Vlr-Vll)) * w_w

      d4 = We2-We1
      d3 = We1-Wp0
      d2 = Wp0-Ww1
      d1 = Ww1-Ww2
      flop = flop + 50.0
      
      include 'muscl.h'
      
      Wrr = We1 - 0.25*((1-ck)*g6+(1+ck)*g5)*ss
      Wrl = Wp0 + 0.25*((1-ck)*g3+(1+ck)*g4)*ss
      Wlr = Wp0 - 0.25*((1-ck)*g4+(1+ck)*g3)*ss
      Wll = Ww1 + 0.25*((1-ck)*g1+(1+ck)*g2)*ss
      fw_r = 0.5*(cr*(Wrr+Wrl) - acr*(Wrr-Wrl)) * w_e
      fw_l = 0.5*(cl*(Wlr+Wll) - acl*(Wlr-Wll)) * w_w
      
      ! 流束の加算　VBCでない面の寄与のみを評価する
      cnv_u = cnv_u + fu_r*c_e - fu_l*c_w
      cnv_v = cnv_v + fv_r*c_e - fv_l*c_w
      cnv_w = cnv_w + fw_r*c_e - fw_l*c_w
      flop = flop + 46.0 + 12.0
			
      ! Y方向 ---------------------------------------
      
      if ( (b_n2 == 0)  ) then
        Un2 = wm1*u_ref + v(i  ,j+1,k  ,1)*wm2
        Vn2 = 2.0*v_ref - v(i  ,j+1,k  ,2)
        Wn2 = wm1*w_ref + v(i  ,j+1,k  ,3)*wm2
        flop = flop + 8.0
        !Un2 = 2.0*u_ref - v(i  ,j+1,k  ,1)
        !Wn2 = 2.0*w_ref - v(i  ,j+1,k  ,3)
      endif
      
      if ( b_n1 == 0 ) then
        Un1 = wm1*u_ref + Up0*wm2
        Vn1 = vq
        Wn1 = wm1*w_ref + Wp0*wm2
        flop = flop + 6.0
        !Un1 = uq
        !Wn1 = wq
      endif
      
      if ( b_s1 == 0 ) then
        Us1 = wm1*u_ref + Up0*wm2
        Vs1 = vq
        Ws1 = wm1*w_ref + Wp0*wm2
        flop = flop + 6.0
        !Us1 = uq
        !Ws1 = wq
      endif
      
      if ( (b_s2 == 0)  ) then
        Us2 = wm1*u_ref + v(i  ,j-1,k  ,1)*wm2
        Vs2 = 2.0*v_ref - v(i  ,j-1,k  ,2)
        Ws2 = wm1*w_ref + v(i  ,j-1,k  ,3)*wm2
        flop = flop + 8.0
        !Us2 = 2.0*u_ref - v(i  ,j-1,k  ,1)
        !Ws2 = 2.0*w_ref - v(i  ,j-1,k  ,3)
      endif
      
      ! 流束　流体のみと固体壁の影響を含み，隣接セルが固体の場合にはマスクする > 216 flops
      cr  = VPn - v_ref
      cl  = VPs - v_ref
      acr = abs(cr)
      acl = abs(cl)
      
      d4 = Un2-Un1
      d3 = Un1-Up0
      d2 = Up0-Us1
      d1 = Us1-Us2
      flop = flop + 10.0
      
      include 'muscl.h'
      
      Urr = Un1 - 0.25*((1-ck)*g6+(1+ck)*g5)*ss
      Url = Up0 + 0.25*((1-ck)*g3+(1+ck)*g4)*ss
      Ulr = Up0 - 0.25*((1-ck)*g4+(1+ck)*g3)*ss
      Ull = Us1 + 0.25*((1-ck)*g1+(1+ck)*g2)*ss
      fu_r = 0.5*(cr*(Urr+Url) - acr*(Urr-Url)) * w_n
      fu_l = 0.5*(cl*(Ulr+Ull) - acl*(Ulr-Ull)) * w_s

      d4 = Vn2-Vn1
      d3 = Vn1-Vp0
      d2 = Vp0-Vs1
      d1 = Vs1-Vs2
      flop = flop + 50.0
      
      include 'muscl.h'
      
      Vrr = Vn1 - 0.25*((1-ck)*g6+(1+ck)*g5)*ss
      Vrl = Vp0 + 0.25*((1-ck)*g3+(1+ck)*g4)*ss
      Vlr = Vp0 - 0.25*((1-ck)*g4+(1+ck)*g3)*ss
      Vll = Vs1 + 0.25*((1-ck)*g1+(1+ck)*g2)*ss
      fv_r = 0.5*(cr*(Vrr+Vrl) - acr*(Vrr-Vrl)) * w_n
      fv_l = 0.5*(cl*(Vlr+Vll) - acl*(Vlr-Vll)) * w_s

      d4 = Wn2-Wn1
      d3 = Wn1-Wp0
      d2 = Wp0-Ws1
      d1 = Ws1-Ws2
      flop = flop + 50.0
      
      include 'muscl.h'
      
      Wrr = Wn1 - 0.25*((1-ck)*g6+(1+ck)*g5)*ss
      Wrl = Wp0 + 0.25*((1-ck)*g3+(1+ck)*g4)*ss
      Wlr = Wp0 - 0.25*((1-ck)*g4+(1+ck)*g3)*ss
      Wll = Ws1 + 0.25*((1-ck)*g1+(1+ck)*g2)*ss
      fw_r = 0.5*(cr*(Wrr+Wrl) - acr*(Wrr-Wrl)) * w_n
      fw_l = 0.5*(cl*(Wlr+Wll) - acl*(Wlr-Wll)) * w_s
      
      ! 流束の加算　VBCでない面の寄与のみを評価する
      cnv_u = cnv_u + fu_r*c_n - fu_l*c_s
      cnv_v = cnv_v + fv_r*c_n - fv_l*c_s
      cnv_w = cnv_w + fw_r*c_n - fw_l*c_s
      flop = flop + 46.0+12.0
			
      ! Z方向 ---------------------------------------
      
      ! 壁面の場合の参照速度の修正
      if ( (b_t2 == 0)  ) then
        Ut2 = wm1*u_ref + v(i  ,j  ,k+1,1)*wm2
        Vt2 = wm1*v_ref + v(i  ,j  ,k+1,2)*wm2
        Wt2 = 2.0*w_ref - v(i  ,j  ,k+1,3)
        flop = flop + 8.0
        !Ut2 = 2.0*u_ref - v(i  ,j  ,k+1,1)
        !Vt2 = 2.0*v_ref - v(i  ,j  ,k+1,2)
      end if
      
      if ( b_t1 == 0 ) then
        Ut1 = wm1*u_ref + Up0*wm2
        Vt1 = wm1*v_ref + Vp0*wm2
        Wt1 = wq
        flop = flop + 6.0
        !Ut1 = uq
        !Vt1 = vq
      end if
      
      if ( b_b1 == 0 ) then
        Ub1 = wm1*u_ref + Up0*wm2
        Vb1 = wm1*v_ref + Vp0*wm2
        Wb1 = wq
        flop = flop + 6.0
        !Ub1 = uq
        !Vb1 = vq
      end if
      
      if ( (b_b2 == 0)  ) then
        Ub2 = wm1*u_ref + v(i  ,j  ,k-1,1)*wm2
        Vb2 = wm1*v_ref + v(i  ,j  ,k-1,2)*wm2
        Wb2 = 2.0*w_ref - v(i  ,j  ,k-1,3)
        flop = flop + 8.0
        !Ub2 = 2.0*u_ref - v(i  ,j  ,k-1,1)
        !Vb2 = 2.0*v_ref - v(i  ,j  ,k-1,2)
      end if
      
      ! 流束　流体のみと固体壁の影響を含み，隣接セルが固体の場合にはマスクする > 216 flops
      cr  = WPt - w_ref
      cl  = WPb - w_ref
      acr = abs(cr)
      acl = abs(cl)

      d4 = Ut2-Ut1
      d3 = Ut1-Up0
      d2 = Up0-Ub1
      d1 = Ub1-Ub2
      flop = flop + 10.0
      
      include 'muscl.h'
      
      Urr = Ut1 - 0.25*((1-ck)*g6+(1+ck)*g5)*ss
      Url = Up0 + 0.25*((1-ck)*g3+(1+ck)*g4)*ss
      Ulr = Up0 - 0.25*((1-ck)*g4+(1+ck)*g3)*ss
      Ull = Ub1 + 0.25*((1-ck)*g1+(1+ck)*g2)*ss
      fu_r = 0.5*(cr*(Urr+Url) - acr*(Urr-Url)) * w_t
      fu_l = 0.5*(cl*(Ulr+Ull) - acl*(Ulr-Ull)) * w_b

      d4 = Vt2-Vt1
      d3 = Vt1-Vp0
      d2 = Vp0-Vb1
      d1 = Vb1-Vb2
      flop = flop + 50.0
      
      include 'muscl.h'
            
      Vrr = Vt1 - 0.25*((1-ck)*g6+(1+ck)*g5)*ss
      Vrl = Vp0 + 0.25*((1-ck)*g3+(1+ck)*g4)*ss
      Vlr = Vp0 - 0.25*((1-ck)*g4+(1+ck)*g3)*ss
      Vll = Vb1 + 0.25*((1-ck)*g1+(1+ck)*g2)*ss
      fv_r = 0.5*(cr*(Vrr+Vrl) - acr*(Vrr-Vrl)) * w_t
      fv_l = 0.5*(cl*(Vlr+Vll) - acl*(Vlr-Vll)) * w_b

      d4 = Wt2-Wt1
      d3 = Wt1-Wp0
      d2 = Wp0-Wb1
      d1 = Wb1-Wb2
      flop = flop + 50.0
      
      include 'muscl.h'
            
      Wrr = Wt1 - 0.25*((1-ck)*g6+(1+ck)*g5)*ss
      Wrl = Wp0 + 0.25*((1-ck)*g3+(1+ck)*g4)*ss
      Wlr = Wp0 - 0.25*((1-ck)*g4+(1+ck)*g3)*ss
      Wll = Wb1 + 0.25*((1-ck)*g1+(1+ck)*g2)*ss
      fw_r = 0.5*(cr*(Wrr+Wrl) - acr*(Wrr-Wrl)) * w_t
      fw_l = 0.5*(cl*(Wlr+Wll) - acl*(Wlr-Wll)) * w_b
      
      ! 流束の加算　VBCでない面の寄与のみを評価する
      cnv_u = cnv_u + fu_r*c_t - fu_l*c_b
      cnv_v = cnv_v + fv_r*c_t - fv_l*c_b
      cnv_w = cnv_w + fw_r*c_t - fw_l*c_b
      flop = flop + 46.0+12.0
      
      ! 粘性項の計算　セル界面の剪断力を計算し，必要に応じて置換する
      uu_e = ( Ue1 - Up0 ) * dh2 ! 1/Re du/dx, e
      uu_w = ( Up0 - Uw1 ) * dh2
      uu_n = ( Un1 - Up0 ) * dh2
      uu_s = ( Up0 - Us1 ) * dh2
      uu_t = ( Ut1 - Up0 ) * dh2
      uu_b = ( Up0 - Ub1 ) * dh2
      
      vv_e = ( Ve1 - Vp0 ) * dh2
      vv_w = ( Vp0 - Vw1 ) * dh2
      vv_n = ( Vn1 - Vp0 ) * dh2
      vv_s = ( Vp0 - Vs1 ) * dh2
      vv_t = ( Vt1 - Vp0 ) * dh2
      vv_b = ( Vp0 - Vb1 ) * dh2
      
      ww_e = ( We1 - Wp0 ) * dh2
      ww_w = ( Wp0 - Ww1 ) * dh2
      ww_n = ( Wn1 - Wp0 ) * dh2
      ww_s = ( Wp0 - Ws1 ) * dh2
      ww_t = ( Wt1 - Wp0 ) * dh2
      ww_b = ( Wp0 - Wb1 ) * dh2
      flop = flop + 36.0
      
      if ( v_mode == 2 ) then ! 壁法則の場合の壁面摩擦による剪断応力の置換
        if ( ibits(bpx, facing_W, 6) /= 0 ) then ! 6面のうちのどれか方向フラグが立っている場合，つまり壁面隣接セル
          u1 = Up0 - u_ref
          u2 = Vp0 - v_ref
          u3 = Wp0 - w_ref
          ug = sqrt(u1*u1 + u2*u2 + u3*u3)
          flop = flop + 10.0
          
          if ( ug > 0.0 ) then
            e1 = abs(u1/ug)
            e2 = abs(u2/ug)
            e3 = abs(u3/ug)
            flop = flop + 6.0
          else
            e1 = 0.0
            e2 = 0.0
            e3 = 0.0
          endif
          
          u_tau = ut(i,j,k)*ut(i,j,k)
          
          if ( ibits(bpx, facing_E, 1) == 1 ) then ! ビットON -> 壁面
          ! uu_eは法線方向で，既に壁面条件で計算されている
            vv_e = u_tau * e2*e2
            ww_e = u_tau * e3*e3
            flop = flop + 4.0
          endif
          
          if ( ibits(bpx, facing_W, 1) == 1 ) then
            vv_w = u_tau * e2*e2
            ww_w = u_tau * e3*e3
            flop = flop + 4.0
          endif
          
          if ( ibits(bpx, facing_N, 1) == 1 ) then
            uu_n = u_tau * e1*e1
            ww_n = u_tau * e3*e3
            flop = flop + 4.0
          endif
          
          if ( ibits(bpx, facing_S, 1) == 1 ) then
            uu_s = u_tau * e1*e1
            ww_s = u_tau * e3*e3
            flop = flop + 4.0
          endif
          
          if ( ibits(bpx, facing_T, 1) == 1 ) then
            uu_t = u_tau * e1*e1
            vv_t = u_tau * e2*e2
            flop = flop + 4.0
          endif
          
          if ( ibits(bpx, facing_B, 1) == 1 ) then
            uu_b = u_tau * e1*e1
            vv_b = u_tau * e2*e2
            flop = flop + 4.0
          endif
          
        endif
      endif
      
      beta = 1.0;
      !if (ibits(bdx, forcing_bit, 1) == 1) then ! 圧力損失コンポの場合
      !  beta = 1.0 - real(ibits( bdx, top_vf, bitw_vf )) * qtz ! 1-体積率
      !endif
      
      EX = ( uu_e * c_e &
           - uu_w * c_w &
           + uu_n * c_n &
           - uu_s * c_s &
           + uu_t * c_t &
           - uu_b * c_b ) * dh1
           
      EY = ( vv_e * c_e &
           - vv_w * c_w &
           + vv_n * c_n &
           - vv_s * c_s &
           + vv_t * c_t &
           - vv_b * c_b ) * dh1
           
      EZ = ( ww_e * c_e &
           - ww_w * c_w &
           + ww_n * c_n &
           - ww_s * c_s &
           + ww_t * c_t &
           - ww_b * c_b ) * dh1
			
      ! 対流項と粘性項の和
      wv(i,j,k,1) = -cnv_u*dh1 + beta*EX*vcs
      wv(i,j,k,2) = -cnv_v*dh1 + beta*EY*vcs
      wv(i,j,k,3) = -cnv_w*dh1 + beta*EZ*vcs
      flop = flop + 72.0
    end do
    end do
    end do

    return
    end subroutine cbc_pvec_muscl

!  ***********************************************************************************************
!> @subroutine cbc_pvec_vibc (wv, sz, g, st, ed, dh, v00, rei, v, bv, odr, vec, v_mode, ofi, flop)
!! @brief 内部速度境界条件による対流項と粘性項の流束の修正
!! @param[out] wv 疑似ベクトルの空間項の評価値
!! @param sz 配列長
!! @param g ガイドセル長
!! @param st ループの開始インデクス
!! @param ed ループの終了インデクス
!! @param dh 格子幅
!! @param v00 参照速度
!! @param rei Reynolds数の逆数
!! @param v 速度ベクトル（u^n, セルセンタ）
!! @param bv BCindex V
!! @param odr 内部境界処理時の速度境界条件のエントリ
!! @param vec 指定する速度ベクトル
!! @param v_mode 粘性項のモード (0=対流項のみ, 1=対流項と粘性項，2=粘性項は壁法則)
!! @param ofi 識別子 (0=specvel, 1=outlfow, 2-wall, 3-symmetric, 4-periodic)
!! @param[out] flop
!! @note vecには，流入条件のとき指定速度，流出条件のとき対流流出速度
!! @todo 内部と外部の分離 do loopの内側に条件分岐を入れているので修正
!! @todo 流出境界はローカルの流束となるように変更する（外部境界参照）
!<
    subroutine cbc_pvec_vibc (wv, sz, g, st, ed, dh, v00, rei, v, bv, odr, vec, v_mode, ofi, flop)
    implicit none
    include 'cbc_f_params.h'
    integer                                                     ::  i, j, k, g, bvx, ofi, odr, v_mode
    integer                                                     ::  b_e1, b_w1, b_n1, b_s1, b_t1, b_b1, b_p
    integer, dimension(3)                                       ::  sz, st, ed
    real                                                        ::  Up0, Ue1, Uw1, Us1, Un1, Ub1, Ut1
    real                                                        ::  Vp0, Ve1, Vw1, Vs1, Vn1, Vb1, Vt1
    real                                                        ::  Wp0, We1, Ww1, Ws1, Wn1, Wb1, Wt1
    real                                                        ::  c_e, c_w, c_n, c_s, c_t, c_b
    real                                                        ::  w_e, w_w, w_n, w_s, w_t, w_b
    real                                                        ::  dh, dh1, dh2, flop, vcs, EX, EY, EZ, rei
    real                                                        ::  u_ref, v_ref, w_ref
    real                                                        ::  cnv_u, cnv_v, cnv_w, cr, cl, acr, acl
    real                                                        ::  u_bc, v_bc, w_bc, u_bc_ref, v_bc_ref, w_bc_ref
    real                                                        ::  fu_r, fu_l, fv_r, fv_l, fw_r, fw_l
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)   ::  v, wv
    real, dimension(0:3)                                        ::  v00
    real, dimension(3)                                          ::  vec
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bv
    
    if      ( v_mode == 0 ) then 
        vcs=0.0
    else if ( v_mode == 1 ) then 
        vcs=1.0
    else if ( v_mode == 2 ) then 
        vcs=1.0
    end if
    
    dh1= 1.0/dh
    dh2= rei*dh1*dh1*vcs

    ! 参照座標系の速度
    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)
    
    ! u_bcは境界速度
    u_bc = vec(1)
    v_bc = vec(2)
    w_bc = vec(3)
    
    ! u_bc_refは参照座標系での境界速度
    u_bc_ref = u_bc + u_ref
    v_bc_ref = v_bc + v_ref
    w_bc_ref = w_bc + w_ref
    
    flop = flop + 6.0
    
    do k=st(3),ed(3)
    do j=st(2),ed(2)
    do i=st(1),ed(1)
      bvx = bv(i,j,k)

      if ( 0 /= iand(bvx, bc_mask30) ) then ! 6面のうちのどれか速度境界フラグが立っている場合
        cnv_u = 0.0
        cnv_v = 0.0
        cnv_w = 0.0
        
        ! セル状態 (0-solid / 1-fluid)
        b_p = ibits(bvx,             State, 1)
        b_w1= ibits(bv(i-1,j  ,k  ), State, 1)
        b_e1= ibits(bv(i+1,j  ,k  ), State, 1)
        b_s1= ibits(bv(i  ,j-1,k  ), State, 1)
        b_n1= ibits(bv(i  ,j+1,k  ), State, 1)
        b_b1= ibits(bv(i  ,j  ,k-1), State, 1)
        b_t1= ibits(bv(i  ,j  ,k+1), State, 1)
      
        ! セル界面のフラグ (0-wall face / 1-fluid)
        w_e = real(b_e1 * b_p)
        w_w = real(b_w1 * b_p)
        w_n = real(b_n1 * b_p)
        w_s = real(b_s1 * b_p)
        w_t = real(b_t1 * b_p)
        w_b = real(b_b1 * b_p)
        
        ! 変数のロード
        Up0 = v(i  ,j  ,k  ,1)
        Vp0 = v(i  ,j  ,k  ,2)
        Wp0 = v(i  ,j  ,k  ,3)
      
        Uw1 = 0.0
        Ue1 = 0.0
        Us1 = 0.0
        Un1 = 0.0
        Ub1 = 0.0
        Ut1 = 0.0
        Vw1 = 0.0
        Ve1 = 0.0
        Vs1 = 0.0
        Vn1 = 0.0
        Vb1 = 0.0
        Vt1 = 0.0
        Ww1 = 0.0
        We1 = 0.0
        Ws1 = 0.0
        Wn1 = 0.0
        Wb1 = 0.0
        Wt1 = 0.0
      
        ! 内部境界のときの各面のBCフラグ ibits() = 0(Normal) / others(BC) >> c_e = 0.0(Normal) / 1.0(BC) 
        c_e = 0.0
        c_w = 0.0
        c_n = 0.0
        c_s = 0.0
        c_t = 0.0
        c_b = 0.0
        if ( ibits(bvx, bc_face_E, bitw_5) == odr ) c_e = 1.0
        if ( ibits(bvx, bc_face_W, bitw_5) == odr ) c_w = 1.0
        if ( ibits(bvx, bc_face_N, bitw_5) == odr ) c_n = 1.0
        if ( ibits(bvx, bc_face_S, bitw_5) == odr ) c_s = 1.0
        if ( ibits(bvx, bc_face_T, bitw_5) == odr ) c_t = 1.0
        if ( ibits(bvx, bc_face_B, bitw_5) == odr ) c_b = 1.0
			
        ! X方向 ---------------------------------------
        if ( c_w == 1.0 ) then
          if ( ofi == id_specvel ) then
            Uw1  = u_bc_ref
            Vw1  = v_bc_ref
            Ww1  = w_bc_ref
            cl   = u_bc
            acl  = abs(cl)
            fu_l = 0.5*(cl*(Up0+Uw1) - acl*(Up0-Uw1))
            fv_l = 0.5*(cl*(Vp0+Vw1) - acl*(Vp0-Vw1))
            fw_l = 0.5*(cl*(Wp0+Ww1) - acl*(Wp0-Ww1))
            flop = flop + 20.0
            
          else if ( ofi == id_outflow ) then
            Uw1  = Up0
            Vw1  = Vp0
            Ww1  = Wp0
            cl   = u_bc
            if ( cl>0.0 ) cl=0.0
            fu_l = cl*Up0
            fv_l = cl*Vp0
            fw_l = cl*Wp0
            flop = flop + 3.0
          endif
        end if
        
        if ( c_e == 1.0 ) then
          if ( ofi == id_specvel ) then
            Ue1  = u_bc_ref
            Ve1  = v_bc_ref
            We1  = w_bc_ref
            cr   = u_bc
            acr  = abs(cr)
            fu_r = 0.5*(cr*(Ue1+Up0) - acr*(Ue1-Up0))
            fv_r = 0.5*(cr*(Ve1+Vp0) - acr*(Ve1-Vp0))
            fw_r = 0.5*(cr*(We1+Wp0) - acr*(We1-Wp0))
            flop = flop + 20.0
            
          else if ( ofi == id_outflow ) then
            Ue1  = Up0
            Ve1  = Vp0
            We1  = Wp0
            cr   = u_bc
            if ( cr<0.0 ) cr=0.0
            fu_r = cr*Up0
            fv_r = cr*Vp0
            fw_r = cr*Wp0
            flop = flop + 3.0
          endif
        end if
        
        cnv_u = cnv_u + fu_r*c_e - fu_l*c_w
        cnv_v = cnv_v + fv_r*c_e - fv_l*c_w
        cnv_w = cnv_w + fw_r*c_e - fw_l*c_w
			
        ! Y方向 ---------------------------------------
        if ( c_s == 1.0 ) then
          if ( ofi == id_specvel ) then
            Us1  = u_bc_ref
            Vs1  = v_bc_ref
            Ws1  = w_bc_ref
            cl   = v_bc
            acl  = abs(cl)
            fu_l = 0.5*(cl*(Up0+Us1) - acl*(Up0-Us1))
            fv_l = 0.5*(cl*(Vp0+Vs1) - acl*(Vp0-Vs1))
            fw_l = 0.5*(cl*(Wp0+Ws1) - acl*(Wp0-Ws1))
            flop = flop + 20.0
            
          else if ( ofi == id_outflow ) then
            Us1  = Up0
            Vs1  = Vp0
            Ws1  = Wp0
            cl   = v_bc
            if ( cl>0.0 ) cl=0.0
            fu_l = cl*Up0
            fv_l = cl*Vp0
            fw_l = cl*Wp0
            flop = flop + 3.0
          endif
        end if
        
        if ( c_n == 1.0 ) then
          if ( ofi == id_specvel ) then
            Un1  = u_bc_ref
            Vn1  = v_bc_ref
            Wn1  = w_bc_ref
            cr   = v_bc
            acr  = abs(cr)
            fu_r = 0.5*(cr*(Un1+Up0) - acr*(Un1-Up0))
            fv_r = 0.5*(cr*(Vn1+Vp0) - acr*(Vn1-Vp0))
            fw_r = 0.5*(cr*(Wn1+Wp0) - acr*(Wn1-Wp0))
            flop = flop + 20.0
            
          else if ( ofi == id_outflow ) then
            Un1  = Up0
            Vn1  = Vp0
            Wn1  = Wp0
            cr   = v_bc
            if ( cr<0.0 ) cr=0.0
            fu_r = cr*Up0
            fv_r = cr*Vp0
            fw_r = cr*Wp0
            flop = flop + 3.0
          endif
        end if
        
        cnv_u = cnv_u + fu_r*c_n - fu_l*c_s
        cnv_v = cnv_v + fv_r*c_n - fv_l*c_s
        cnv_w = cnv_w + fw_r*c_n - fw_l*c_s
			
        ! Z方向 ---------------------------------------
        if ( c_b == 1.0 ) then
          if ( ofi == id_specvel ) then
            Ub1  = u_bc_ref
            Vb1  = v_bc_ref
            Wb1  = w_bc_ref
            cl   = w_bc
            acl  = abs(cl)
            fu_l = 0.5*(cl*(Up0+Ub1) - acl*(Up0-Ub1))
            fv_l = 0.5*(cl*(Vp0+Vb1) - acl*(Vp0-Vb1))
            fw_l = 0.5*(cl*(Wp0+Wb1) - acl*(Wp0-Wb1))
            flop = flop + 20.0
            
          else if ( ofi == id_outflow ) then
            Ub1  = Up0
            Vb1  = Vp0
            Wb1  = Wp0
            cl   = w_bc
            if ( cl>0.0 ) cl=0.0
            fu_l = cl*Up0
            fv_l = cl*Vp0
            fw_l = cl*Wp0
            flop = flop + 3.0
          endif
        end if

        if ( c_t == 1.0 ) then
          if ( ofi == id_specvel ) then
            Ut1  = u_bc_ref
            Vt1  = v_bc_ref
            Wt1  = w_bc_ref
            cr   = w_bc
            acr  = abs(cr)
            fu_r = 0.5*(cr*(Ut1+Up0) - acr*(Ut1-Up0))
            fv_r = 0.5*(cr*(Vt1+Vp0) - acr*(Vt1-Vp0))
            fw_r = 0.5*(cr*(Wt1+Wp0) - acr*(Wt1-Wp0))
            flop = flop + 20.0
            
          else if ( ofi == id_outflow ) then
            Ut1  = Up0
            Vt1  = Vp0
            Wt1  = Wp0
            cr   = w_bc
            if ( cr<0.0 ) cr=0.0
            fu_r = cr*Up0
            fv_r = cr*Vp0
            fw_r = cr*Wp0
            flop = flop + 3.0
          endif
        end if
        
        cnv_u = cnv_u + fu_r*c_t - fu_l*c_b
        cnv_v = cnv_v + fv_r*c_t - fv_l*c_b
        cnv_w = cnv_w + fw_r*c_t - fw_l*c_b
      
        ! 粘性項
        EX = ( Ue1 - Up0 ) * c_e &
           + ( Uw1 - Up0 ) * c_w &
           + ( Un1 - Up0 ) * c_n &
           + ( Us1 - Up0 ) * c_s &
           + ( Ut1 - Up0 ) * c_t &
           + ( Ub1 - Up0 ) * c_b
        EY = ( Ve1 - Vp0 ) * c_e &
           + ( Vw1 - Vp0 ) * c_w &
           + ( Vn1 - Vp0 ) * c_n &
           + ( Vs1 - Vp0 ) * c_s &
           + ( Vt1 - Vp0 ) * c_t &
           + ( Vb1 - Vp0 ) * c_b
        EZ = ( We1 - Wp0 ) * c_e &
           + ( Ww1 - Wp0 ) * c_w &
           + ( Wn1 - Wp0 ) * c_n &
           + ( Ws1 - Wp0 ) * c_s &
           + ( Wt1 - Wp0 ) * c_t &
           + ( Wb1 - Wp0 ) * c_b

        wv(i,j,k,1) = wv(i,j,k,1) + ( -cnv_u*dh1 + EX*dh2 )
        wv(i,j,k,2) = wv(i,j,k,2) + ( -cnv_v*dh1 + EY*dh2 )
        wv(i,j,k,3) = wv(i,j,k,3) + ( -cnv_w*dh1 + EZ*dh2 )

        flop = flop + 78.0
      endif
    end do
    end do
    end do
    
    return
    end subroutine cbc_pvec_vibc

!  ****************************************************************************************
!> @subroutine cbc_pvec_vobc (wv, sz, g, dh, v00, rei, v0, bv, vec, v_mode, ofi, face, flop)
!! @brief 外部速度境界条件による対流項と粘性項の流束の修正
!! @param[out] wv 疑似ベクトルの空間項の評価値
!! @param sz 配列長
!! @param g ガイドセル長
!! @param dh 格子幅
!! @param v00 参照速度
!! @param rei Reynolds数の逆数
!! @param v0 速度ベクトル（n-step）
!! @param bv BCindex V
!! @param vec 指定する速度ベクトル
!! @param v_mode 粘性項のモード (0=対流項のみ, 1=対流項と粘性項，2=粘性項は壁法則)
!! @param ofi 識別子 (0=specvel, 1=outlfow, 2-wall, 3-symmetric)
!! @param face 外部境界処理のときの面番号
!! @param[out] flop
!! @note vecには，流入条件のとき指定速度，流出境界の流束はローカルのセルフェイス速度を使うこと
!! @todo 内部と外部の分離 do loopの内側に条件分岐を入れているので修正
!<
    subroutine cbc_pvec_vobc (wv, sz, g, dh, v00, rei, v0, bv, vec, v_mode, ofi, face, flop)
    implicit none
    include 'cbc_f_params.h'
    integer                                                     ::  i, j, k, g, bvx, ofi, face, v_mode
    integer                                                     ::  ix, jx, kx
    integer, dimension(3)                                       ::  sz
    real                                                        ::  Up0, Ue1, Uw1, Us1, Un1, Ub1, Ut1
    real                                                        ::  Vp0, Ve1, Vw1, Vs1, Vn1, Vb1, Vt1
    real                                                        ::  Wp0, We1, Ww1, Ws1, Wn1, Wb1, Wt1
    real                                                        ::  dh, dh1, dh2, flop, vcs, EX, EY, EZ, rei
    real                                                        ::  u_ref, v_ref, w_ref
    real                                                        ::  cnv_u, cnv_v, cnv_w, cr, cl, acr, acl
    real                                                        ::  u_bc, v_bc, w_bc, u_bc_ref, v_bc_ref, w_bc_ref
    real                                                        ::  fu_r, fu_l, fv_r, fv_l, fw_r, fw_l
    real                                                        ::  Ue0, Uw0, Vs0, Vn0, Wb0, Wt0
    real                                                        ::  Ue, Uw, Vn, Vs, Wt, Wb
    real                                                        ::  b_w, b_e, b_s, b_n, b_b, b_t
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)   ::  v0, wv
    real, dimension(0:3)                                        ::  v00
    real, dimension(3)                                          ::  vec
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bv
    
    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    
    if      ( v_mode == 0 ) then 
        vcs=0.0
    else if ( v_mode == 1 ) then 
        vcs=1.0
    else if ( v_mode == 2 ) then 
        vcs=1.0
    end if
    
    dh1= 1.0/dh
    dh2= rei*dh1*dh1*vcs

    ! 参照座標系の速度
    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)
    
    ! u_bcは境界速度
    u_bc = vec(1)
    v_bc = vec(2)
    w_bc = vec(3)
    
    ! u_bc_refは参照座標系での境界速度
    u_bc_ref = u_bc + u_ref
    v_bc_ref = v_bc + v_ref
    w_bc_ref = w_bc + w_ref
    
    flop = flop + 6.0
    
    FACES : select case (face)
    case (X_minus)
      i = 1
      do k=1,kx
      do j=1,jx
        bvx = bv(i,j,k)
        
        if ( ibits(bvx, bc_face_W, bitw_5) == obc_mask ) then
          Up0 = v0(i  ,j  ,k  ,1)
          Vp0 = v0(i  ,j  ,k  ,2)
          Wp0 = v0(i  ,j  ,k  ,3)
          Uw1 = 0.0
          Vw1 = 0.0
          Ww1 = 0.0
          
          if ( ofi == id_specvel ) then
            Uw1  = u_bc_ref
            Vw1  = v_bc_ref
            Ww1  = w_bc_ref
            !cl   = u_bc
            !acl  = abs(cl)
            !fu_l = 0.5*(cl*(Up0+Uw1) - acl*(Up0-Uw1))
            !fv_l = 0.5*(cl*(Vp0+Vw1) - acl*(Vp0-Vw1))
            !fw_l = 0.5*(cl*(Wp0+Ww1) - acl*(Wp0-Ww1))
            fu_l = u_bc * u_bc_ref
            fv_l = v_bc * v_bc_ref
            fw_l = w_bc * w_bc_ref
            flop = flop + 3.0
            
          else if ( ofi == id_outflow ) then
            include 'd_o_o_p.h'
            Uw1  = v0(i-1,j  ,k  ,1)
            Vw1  = v0(i-1,j  ,k  ,2)
            Ww1  = v0(i-1,j  ,k  ,3)
            cl   = Ue + (Vn - Vs + Wt - Wb) - u_ref
            if ( cl>0.0 ) cl=0.0
            fu_l = cl * Up0
            fv_l = cl * Vp0
            fw_l = cl * Wp0
            flop = flop + 8.0
            
          else if ( ofi == id_wall ) then
            Uw1  = 2.0*u_bc_ref - Up0
            Vw1  = 2.0*v_bc_ref - Vp0
            Ww1  = 2.0*w_bc_ref - Wp0
            fu_l = 0.0
            fv_l = 0.0
            fw_l = 0.0
            flop = flop + 6.0
            
          else if ( ofi == id_symmetric ) then
            Uw1  = -Up0
            Vw1  =  Vp0
            Ww1  =  Wp0
            fu_l = 0.0
            fv_l = 0.0
            fw_l = 0.0
          endif
          
          cnv_u = - fu_l
          cnv_v = - fv_l
          cnv_w = - fw_l
          
          EX = Uw1 - Up0
          EY = Vw1 - Vp0
          EZ = Ww1 - Wp0

          wv(i,j,k,1) = wv(i,j,k,1) + ( -cnv_u*dh1 + EX*dh2 )
          wv(i,j,k,2) = wv(i,j,k,2) + ( -cnv_v*dh1 + EY*dh2 )
          wv(i,j,k,3) = wv(i,j,k,3) + ( -cnv_w*dh1 + EZ*dh2 )
          flop = flop + 21.0
        endif
      end do
      end do
      
    case (X_plus)
      i = ix
      do k=1,kx
      do j=1,jx
        bvx = bv(i,j,k)
        
        if ( ibits(bvx, bc_face_E, bitw_5) == obc_mask ) then
          Up0 = v0(i  ,j  ,k  ,1)
          Vp0 = v0(i  ,j  ,k  ,2)
          Wp0 = v0(i  ,j  ,k  ,3)
          Ue1 = 0.0
          Ve1 = 0.0
          We1 = 0.0
          
          if ( ofi == id_specvel ) then
            Ue1  = u_bc_ref
            Ve1  = v_bc_ref
            We1  = w_bc_ref
            cr   = u_bc
            acr  = abs(cr)
            fu_r = 0.5*(cr*(Ue1+Up0) - acr*(Ue1-Up0))
            fv_r = 0.5*(cr*(Ve1+Vp0) - acr*(Ve1-Vp0))
            fw_r = 0.5*(cr*(We1+Wp0) - acr*(We1-Wp0))
            flop = flop + 20.0
            
          else if ( ofi == id_outflow ) then
            include 'd_o_o_p.h'
            Ue1  = v0(i+1,j  ,k  ,1)
            Ve1  = v0(i+1,j  ,k  ,2)
            We1  = v0(i+1,j  ,k  ,3)
            cr   = Uw - (Vn - Vs + Wt - Wb) - u_ref
            if ( cr<0.0 ) cr=0.0
            fu_r = cr * Up0
            fv_r = cr * Vp0
            fw_r = cr * Wp0
            flop = flop + 8.0
            
          else if ( ofi == id_wall ) then
            Ue1  = 2.0*u_bc_ref - Up0
            Ve1  = 2.0*v_bc_ref - Vp0
            We1  = 2.0*w_bc_ref - Wp0
            fu_r = 0.0
            fv_r = 0.0
            fw_r = 0.0
            flop = flop + 6.0
            
          else if ( ofi == id_symmetric ) then
            Ue1  = -Up0
            Ve1  =  Vp0
            We1  =  Wp0
            fu_r = 0.0
            fv_r = 0.0
            fw_r = 0.0
          endif
        
          cnv_u = fu_r
          cnv_v = fv_r
          cnv_w = fw_r

          EX = Ue1 - Up0
          EY = Ve1 - Vp0
          EZ = We1 - Wp0

          wv(i,j,k,1) = wv(i,j,k,1) + ( -cnv_u*dh1 + EX*dh2 )
          wv(i,j,k,2) = wv(i,j,k,2) + ( -cnv_v*dh1 + EY*dh2 )
          wv(i,j,k,3) = wv(i,j,k,3) + ( -cnv_w*dh1 + EZ*dh2 )
        endif
        flop = flop + 21.0
      end do
      end do
      
    case (Y_minus)
      j = 1
      do k=1,kx
      do i=1,ix
        bvx = bv(i,j,k)
        
        if ( ibits(bvx, bc_face_S, bitw_5) == obc_mask ) then
          Up0 = v0(i  ,j  ,k  ,1)
          Vp0 = v0(i  ,j  ,k  ,2)
          Wp0 = v0(i  ,j  ,k  ,3)
          Us1 = 0.0
          Vs1 = 0.0
          Ws1 = 0.0
          
          if ( ofi == id_specvel ) then
            Us1  = u_bc_ref
            Vs1  = v_bc_ref
            Ws1  = w_bc_ref
            cl   = v_bc
            acl  = abs(cl)
            fu_l = 0.5*(cl*(Up0+Us1) - acl*(Up0-Us1))
            fv_l = 0.5*(cl*(Vp0+Vs1) - acl*(Vp0-Vs1))
            fw_l = 0.5*(cl*(Wp0+Ws1) - acl*(Wp0-Ws1))
            flop = flop + 20.0
            
          else if ( ofi == id_outflow ) then
            include 'd_o_o_p.h'
            Us1  = v0(i  ,j-1,k  ,1)
            Vs1  = v0(i  ,j-1,k  ,2)
            Ws1  = v0(i  ,j-1,k  ,3)
            cl   = Vn + (Ue - Uw + Wt - Wb) - v_ref
            if ( cl>0.0 ) cl=0.0
            fu_l = cl * Up0
            fv_l = cl * Vp0
            fw_l = cl * Wp0
            flop = flop + 8.0
            
          else if ( ofi == id_wall ) then
            Us1  = 2.0*u_bc_ref - Up0
            Vs1  = 2.0*v_bc_ref - Vp0
            Ws1  = 2.0*w_bc_ref - Wp0
            fu_l = 0.0
            fv_l = 0.0
            fw_l = 0.0
            flop = flop + 6.0
            
          else if ( ofi == id_symmetric ) then
            Us1  =  Up0
            Vs1  = -Vp0
            Ws1  =  Wp0
            fu_l = 0.0
            fv_l = 0.0
            fw_l = 0.0
          endif

          cnv_u = - fu_l
          cnv_v = - fv_l
          cnv_w = - fw_l
        
          EX = Us1 - Up0
          EY = Vs1 - Vp0
          EZ = Ws1 - Wp0
          
          wv(i,j,k,1) = wv(i,j,k,1) + ( -cnv_u*dh1 + EX*dh2 )
          wv(i,j,k,2) = wv(i,j,k,2) + ( -cnv_v*dh1 + EY*dh2 )
          wv(i,j,k,3) = wv(i,j,k,3) + ( -cnv_w*dh1 + EZ*dh2 )
          flop = flop + 21.0
        endif
      end do
      end do
      
    case (Y_plus)
      j = jx
      do k=1,kx
      do i=1,ix
        bvx = bv(i,j,k)
        
        if ( ibits(bvx, bc_face_N, bitw_5) == obc_mask ) then
          Up0 = v0(i  ,j  ,k  ,1)
          Vp0 = v0(i  ,j  ,k  ,2)
          Wp0 = v0(i  ,j  ,k  ,3)
          Un1 = 0.0
          Vn1 = 0.0
          Wn1 = 0.0
          
          if ( ofi == id_specvel ) then
            Un1  = u_bc_ref
            Vn1  = v_bc_ref
            Wn1  = w_bc_ref
            cr   = v_bc
            acr  = abs(cr)
            fu_r = 0.5*(cr*(Un1+Up0) - acr*(Un1-Up0))
            fv_r = 0.5*(cr*(Vn1+Vp0) - acr*(Vn1-Vp0))
            fw_r = 0.5*(cr*(Wn1+Wp0) - acr*(Wn1-Wp0))
            flop = flop + 20.0
            
          else if ( ofi == id_outflow ) then
            include 'd_o_o_p.h'
            Un1  = v0(i  ,j+1,k  ,1)
            Vn1  = v0(i  ,j+1,k  ,2)
            Wn1  = v0(i  ,j+1,k  ,3)
            cr   = Vs - (Ue - Uw + Wt - Wb) - v_ref
            if ( cr<0.0 ) cr=0.0
            fu_r = cr * Up0
            fv_r = cr * Vp0
            fw_r = cr * Wp0
            flop = flop + 8.0
            
          else if ( ofi == id_wall ) then
            Un1  = 2.0*u_bc_ref - Up0
            Vn1  = 2.0*v_bc_ref - Vp0
            Wn1  = 2.0*w_bc_ref - Wp0
            fu_r = 0.0
            fv_r = 0.0
            fw_r = 0.0
            flop = flop + 6.0
            
          else if ( ofi == id_symmetric ) then
            Un1  =  Up0
            Vn1  = -Vp0
            Wn1  =  Wp0
            fu_r = 0.0
            fv_r = 0.0
            fw_r = 0.0
          endif

          cnv_u = fu_r
          cnv_v = fv_r
          cnv_w = fw_r

          EX = Un1 - Up0
          EY = Vn1 - Vp0
          EZ = Wn1 - Wp0
          
          wv(i,j,k,1) = wv(i,j,k,1) + ( -cnv_u*dh1 + EX*dh2 )
          wv(i,j,k,2) = wv(i,j,k,2) + ( -cnv_v*dh1 + EY*dh2 )
          wv(i,j,k,3) = wv(i,j,k,3) + ( -cnv_w*dh1 + EZ*dh2 )
          flop = flop + 21.0
        endif
      end do
      end do
      
    case (Z_minus)
      k = 1
      do j=1,jx
      do i=1,ix
        bvx = bv(i,j,k)
        
        if ( ibits(bvx, bc_face_B, bitw_5) == obc_mask ) then
          Up0 = v0(i  ,j  ,k  ,1)
          Vp0 = v0(i  ,j  ,k  ,2)
          Wp0 = v0(i  ,j  ,k  ,3)
          Ub1 = 0.0
          Vb1 = 0.0
          Wb1 = 0.0
          
          if ( ofi == id_specvel ) then
            Ub1  = u_bc_ref
            Vb1  = v_bc_ref
            Wb1  = w_bc_ref
            cl   = w_bc
            acl  = abs(cl)
            fu_l = 0.5*(cl*(Up0+Ub1) - acl*(Up0-Ub1))
            fv_l = 0.5*(cl*(Vp0+Vb1) - acl*(Vp0-Vb1))
            fw_l = 0.5*(cl*(Wp0+Wb1) - acl*(Wp0-Wb1))
            flop = flop + 20.0
            
          else if ( ofi == id_outflow ) then
            include 'd_o_o_p.h'
            Ub1  = v0(i  ,j  ,k-1,1)
            Vb1  = v0(i  ,j  ,k-1,2)
            Wb1  = v0(i  ,j  ,k-1,3)
            cl   = Wt + (Ue - Uw + Vn - Vs) - w_ref
            if ( cl>0.0 ) cl=0.0
            fu_l = cl * Up0
            fv_l = cl * Vp0
            fw_l = cl * Wp0
            flop = flop + 8.0
            
          else if ( ofi == id_wall ) then
            Ub1  = 2.0*u_bc_ref - Up0
            Vb1  = 2.0*v_bc_ref - Vp0
            Wb1  = 2.0*w_bc_ref - Wp0
            fu_l = 0.0
            fv_l = 0.0
            fw_l = 0.0
            flop = flop + 6.0
            
          else if ( ofi == id_symmetric ) then
            Ub1  =  Up0
            Vb1  =  Vp0
            Wb1  = -Wp0
            fu_l = 0.0
            fv_l = 0.0
            fw_l = 0.0
          endif

          cnv_u = - fu_l
          cnv_v = - fv_l
          cnv_w = - fw_l
          
          EX = Ub1 - Up0
          EY = Vb1 - Vp0
          EZ = Wb1 - Wp0
          
          wv(i,j,k,1) = wv(i,j,k,1) + ( -cnv_u*dh1 + EX*dh2 )
          wv(i,j,k,2) = wv(i,j,k,2) + ( -cnv_v*dh1 + EY*dh2 )
          wv(i,j,k,3) = wv(i,j,k,3) + ( -cnv_w*dh1 + EZ*dh2 )
          flop = flop + 21.0
        endif
      end do
      end do
      
    case (Z_plus)
      k = kx
      do j=1,jx
      do i=1,ix
        bvx = bv(i,j,k)
        
        if ( ibits(bvx, bc_face_T, bitw_5) == obc_mask ) then
          Up0 = v0(i  ,j  ,k  ,1)
          Vp0 = v0(i  ,j  ,k  ,2)
          Wp0 = v0(i  ,j  ,k  ,3)
          Ut1 = 0.0
          Vt1 = 0.0
          Wt1 = 0.0
          
          if ( ofi == id_specvel ) then
            Ut1  = u_bc_ref
            Vt1  = v_bc_ref
            Wt1  = w_bc_ref
            cr   = w_bc
            acr  = abs(cr)
            fu_r = 0.5*(cr*(Ut1+Up0) - acr*(Ut1-Up0))
            fv_r = 0.5*(cr*(Vt1+Vp0) - acr*(Vt1-Vp0))
            fw_r = 0.5*(cr*(Wt1+Wp0) - acr*(Wt1-Wp0))
            flop = flop + 20.0
            
          else if ( ofi == id_outflow ) then
            include 'd_o_o_p.h'
            Ut1  = v0(i  ,j  ,k+1,1)
            Vt1  = v0(i  ,j  ,k+1,2)
            Wt1  = v0(i  ,j  ,k+1,3)
            cr   = Wb - (Ue - Uw + Vn - Vs) - w_ref
            if ( cr<0.0 ) cr=0.0
            fu_r = cr * Up0
            fv_r = cr * Vp0
            fw_r = cr * Wp0
            flop = flop + 8.0
            
          else if ( ofi == id_wall ) then
            Ut1  = 2.0*u_bc_ref - Up0
            Vt1  = 2.0*v_bc_ref - Vp0
            Wt1  = 2.0*w_bc_ref - Wp0
            fu_r = 0.0
            fv_r = 0.0
            fw_r = 0.0
            flop = flop + 6.0
            
          else if ( ofi == id_symmetric ) then
            Ut1  =  Up0
            Vt1  =  Vp0
            Wt1  = -Wp0
            fu_r = 0.0
            fv_r = 0.0
            fw_r = 0.0
          endif

          cnv_u = fu_r
          cnv_v = fv_r
          cnv_w = fw_r
          
          EX = Ut1 - Up0
          EY = Vt1 - Vp0
          EZ = Wt1 - Wp0
          
          wv(i,j,k,1) = wv(i,j,k,1) + ( -cnv_u*dh1 + EX*dh2 )
          wv(i,j,k,2) = wv(i,j,k,2) + ( -cnv_v*dh1 + EY*dh2 )
          wv(i,j,k,3) = wv(i,j,k,3) + ( -cnv_w*dh1 + EZ*dh2 )
          flop = flop + 21.0
        endif
      end do
      end do
      
    case default
    end select FACES
    
    return
    end subroutine cbc_pvec_vobc

!  **********************************************************************************
!> @subroutine cbc_update_vec (v, div, sz, g, dt, dh, vc, p, bp, bv, v00, coef, flop)
!! @brief 次ステップのセルセンターの速度を更新
!! @param[out] v n+1時刻の速度ベクトル
!! @param[out] div 速度の発散値
!! @param sz 配列長
!! @param g ガイドセル長
!! @param dt 時間積分幅
!! @param dh 格子幅
!! @param vc 疑似速度ベクトル
!! @param p 圧力
!! @param bp BCindex P
!! @param bv BCindex V
!! @param v00 参照速度
!! @param coef 係数
!! @param[out] flop flop count
!! @note 
!!    - actvのマスクはSPEC_VEL/OUTFLOWの参照セルをマスクしないようにbvを使う
!!    - VBC(OUTFLOW, SPEC_VEL)の参照セルでは不定値となるが，InnerVBC()で上書き
!<
    subroutine cbc_update_vec (v, div, sz, g, dt, dh, vc, p, bp, bv, v00, coef, flop)
    implicit none
    include 'cbc_f_params.h'
    integer                                                     ::  i, j, k, ix, jx, kx, g, bpx, bvx
    integer, dimension(3)                                       ::  sz
    real                                                        ::  dh, dt, dd, flop, coef, actv, r_actv
    real                                                        ::  pc, px, py, pz, pxw, pxe, pys, pyn, pzb, pzt
    real                                                        ::  u_ref, v_ref, w_ref
    real                                                        ::  Ue0, Uw0, Vn0, Vs0, Wt0, Wb0, Up0, Vp0, Wp0
    real                                                        ::  Ue, Uw, Vn, Vs, Wt, Wb
    real                                                        ::  c1, c2, c3, c4, c5, c6
    real                                                        ::  N_e, N_w, N_n, N_s, N_t, N_b
    real                                                        ::  b_w, b_e, b_s, b_n, b_b, b_t
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)   ::  v, vc
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  div, p
    real, dimension(0:3)                                        ::  v00
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bp, bv
    
    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    dd = dt/dh
    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)
    flop = flop + real(ix*jx*kx)*95.0 + 2.0

    do k=1,kx
    do j=1,jx
    do i=1,ix
      bpx = bp(i,j,k)
      bvx = bv(i,j,k)
      actv = real(ibits(bvx, State,  1))
      r_actv = 1.0 - actv
      
      c1 = 1.0
      c2 = 1.0
      c3 = 1.0
      c4 = 1.0
      c5 = 1.0
      c6 = 1.0
      
      ! Neumann条件のとき，0.0
      N_w = real(ibits(bpx, bc_n_W, 1))  ! w
      N_e = real(ibits(bpx, bc_n_E, 1))  ! e
      N_s = real(ibits(bpx, bc_n_S, 1))  ! s
      N_n = real(ibits(bpx, bc_n_N, 1))  ! n
      N_b = real(ibits(bpx, bc_n_B, 1))  ! b
      N_t = real(ibits(bpx, bc_n_T, 1))  ! t
      
      ! 疑似ベクトル
      Up0 = vc(i  ,j  ,k  ,1)
      Vp0 = vc(i  ,j  ,k  ,2)
      Wp0 = vc(i  ,j  ,k  ,3)
      Uw0 = vc(i-1,j  ,k  ,1)
      Ue0 = vc(i+1,j  ,k  ,1)
      Vs0 = vc(i  ,j-1,k  ,2)
      Vn0 = vc(i  ,j+1,k  ,2)
      Wb0 = vc(i  ,j  ,k-1,3)
      Wt0 = vc(i  ,j  ,k+1,3)      
      
      ! 壁面による修正 (0-solid / 1-fluid)
      b_w = real( ibits(bv(i-1,j  ,k  ), State, 1) )
      b_e = real( ibits(bv(i+1,j  ,k  ), State, 1) )
      b_s = real( ibits(bv(i  ,j-1,k  ), State, 1) )
      b_n = real( ibits(bv(i  ,j+1,k  ), State, 1) )
      b_b = real( ibits(bv(i  ,j  ,k-1), State, 1) )
      b_t = real( ibits(bv(i  ,j  ,k+1), State, 1) )

      Uw = 0.5*( Up0 + Uw0 )*b_w + (1.0-b_w)*u_ref ! 36 flop
      Ue = 0.5*( Up0 + Ue0 )*b_e + (1.0-b_e)*u_ref
      Vs = 0.5*( Vp0 + Vs0 )*b_s + (1.0-b_s)*v_ref
      Vn = 0.5*( Vp0 + Vn0 )*b_n + (1.0-b_n)*v_ref
      Wb = 0.5*( Wp0 + Wb0 )*b_b + (1.0-b_b)*w_ref
      Wt = 0.5*( Wp0 + Wt0 )*b_t + (1.0-b_t)*w_ref
      
      ! c=0.0(VBC), 1.0(Fluid); VBCは内部と外部の両方
      if ( ibits(bvx, bc_face_W, bitw_5) /= 0 ) c1 = 0.0
      if ( ibits(bvx, bc_face_E, bitw_5) /= 0 ) c2 = 0.0
      if ( ibits(bvx, bc_face_S, bitw_5) /= 0 ) c3 = 0.0
      if ( ibits(bvx, bc_face_N, bitw_5) /= 0 ) c4 = 0.0
      if ( ibits(bvx, bc_face_B, bitw_5) /= 0 ) c5 = 0.0
      if ( ibits(bvx, bc_face_T, bitw_5) /= 0 ) c6 = 0.0
      
      ! 圧力勾配
      pc  = p(i,  j,  k  )
      pxw = (pc - p(i-1,j  ,k  )) * N_w
      pxe = (p(i+1,j  ,k  ) - pc) * N_e
      pys = (pc - p(i  ,j-1,k  )) * N_s
      pyn = (p(i  ,j+1,k  ) - pc) * N_n
      pzb = (pc - p(i  ,j  ,k-1)) * N_b
      pzt = (p(i  ,j  ,k+1) - pc) * N_t
      px = 0.5*(pxe + pxw)
      py = 0.5*(pyn + pys)
      pz = 0.5*(pzt + pzb)

      ! 発散値 VBCの寄与と壁面の影響は除外
      div(i,j,k) = ((Ue - dd * pxe) * c2 &
                  - (Uw - dd * pxw) * c1 &
                  + (Vn - dd * pyn) * c4 &
                  - (Vs - dd * pys) * c3 &
                  + (Wt - dd * pzt) * c6 &
                  - (Wb - dd * pzb) * c5 ) * coef * actv
      
      ! セルセンタの速度更新
      v(i,j,k,1) = ( Up0-px*dd )*actv + r_actv*u_ref
      v(i,j,k,2) = ( Vp0-py*dd )*actv + r_actv*v_ref
      v(i,j,k,3) = ( Wp0-pz*dd )*actv + r_actv*w_ref
    end do
    end do
    end do
    
    return
    end subroutine cbc_update_vec

!  **********************************************************
!> @subroutine cbc_div (div, sz, g, coef, v0, bv, v00, flopc)
!! @brief 速度の発散を計算する
!! @param div 速度の発散値
!! @param sz 配列長
!! @param g ガイドセル長
!! @param coef 係数
!! @param v0 疑似ベクトル
!! @param bv BCindex V
!! @param v00 参照速度
!! @param[out] flop flop count
!! @note divの値は計算対象のFluidセルでのみ有効であればよいので，bvを使ってよい
!<
    subroutine cbc_div (div, sz, g, coef, v0, bv, v00, flopc)
    implicit none
    include 'cbc_f_params.h'
    integer                                                     ::  i, j, k, ix, jx, kx, g, bvx
    integer, dimension(3)                                       ::  sz
    real                                                        ::  coef, flop, flopc
    real                                                        ::  Ue, Uw, Vn, Vs, Wt, Wb, actv
    real                                                        ::  c_e, c_w, c_n, c_s, c_t, c_b
    real                                                        ::  u_ref, v_ref, w_ref
    real                                                        ::  Ue0, Uw0, Vn0, Vs0, Wt0, Wb0, Up0, Vp0, Wp0
    real                                                        ::  b_w, b_e, b_s, b_n, b_b, b_t
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)   ::  v0
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  div
    real, dimension(0:3)                                        ::  v00
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bv

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)
    flopc = flopc + real(ix*jx*kx)*49.0 ! flopはd_o_o_p.hで使われている

    do k=1,kx
    do j=1,jx
    do i=1,ix
      bvx = bv(i,j,k)
      actv= real(ibits(bvx, State, 1))
      
      include 'd_o_o_p.h'
      
      ! 各面のVBCフラグ ibits() = 0(Normal) / others(BC) >> c_e = 1.0(Normal) / 0.0(BC) 
      c_e = 1.0
      c_w = 1.0
      c_n = 1.0
      c_s = 1.0
      c_t = 1.0
      c_b = 1.0
      if ( ibits(bvx, bc_face_E, bitw_5) /= 0 ) c_e = 0.0
      if ( ibits(bvx, bc_face_W, bitw_5) /= 0 ) c_w = 0.0
      if ( ibits(bvx, bc_face_N, bitw_5) /= 0 ) c_n = 0.0
      if ( ibits(bvx, bc_face_S, bitw_5) /= 0 ) c_s = 0.0
      if ( ibits(bvx, bc_face_T, bitw_5) /= 0 ) c_t = 0.0
      if ( ibits(bvx, bc_face_B, bitw_5) /= 0 ) c_b = 0.0
      
      ! VBC面の影響をフラグで無効化 >> OBC_SPEC_VEL, OBC_WALL, OBC_SYMMETRIC, OBC_OUTFLOW
      div(i,j,k) =( Ue*c_e - Uw*c_w + Vn*c_n - Vs*c_s + Wt*c_t - Wb*c_b ) * coef * actv
    end do
    end do
    end do

    return
    end subroutine cbc_div
    
!  ***********************************************
!> @subroutine cbc_ee (vc, sz, g, dt, v, bd, flop)
!! @brief 疑似ベクトルの時間積分（Euler陽解法）
!! @param[in/out] vc 疑似ベクトル
!! @param sz 配列長
!! @param g ガイドセル長
!! @param dt 時間積分幅
!! @param v 速度ベクトル（n-step, collocated）
!! @param bd BCindex ID
!! @param[out] flop
!! @note ここのマスクはIDのこと，VSPEC, OUTFLOWの増分をキャンセルするため
!<
    subroutine cbc_ee (vc, sz, g, dt, v, bd, flop)
    implicit none
    include 'cbc_f_params.h'
    integer                                                     ::  i, j, k, ix, jx, kx, g
    integer, dimension(3)                                       ::  sz
    real                                                        ::  flop, actv, dt
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)   ::  vc, v
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bd
    
    ix = sz(1)
    jx = sz(2)
    kx = sz(3)

    flop = flop + real(ix*jx*kx)*9.0
    
    do k=1,kx
    do j=1,jx
    do i=1,ix
      actv = real(ibits(bd(i,j,k), State, 1))

      vc(i,j,k,1) = v(i,j,k,1) + dt*vc(i,j,k,1)* actv
      vc(i,j,k,2) = v(i,j,k,2) + dt*vc(i,j,k,2)* actv
      vc(i,j,k,3) = v(i,j,k,3) + dt*vc(i,j,k,3)* actv
    end do
    end do
    end do
    
    return
    end subroutine cbc_ee

!  *************************************************************************
!> @subroutine cbc_vis_ee (vc, sz, g, dh, dt, v00, rei, wv, v, bx, cf, flop)
!! @brief 粘性項の計算（Euler陽解法，壁面境界の影響のみ考慮する）
!! @param[out] vc 疑似ベクトル
!! @param sz 配列長
!! @param g ガイドセル長
!! @param dh 格子幅
!! @param dt 時間積分幅
!! @param v00 参照速度
!! @param rei 1/Re
!! @param wv 疑似ベクトルの対流項分
!! @param v 速度ベクトル（n-step, collocated）
!! @param bx BC index for v
!! @param cf 粘性項の係数 (0.0=only convection, 1.0=convection+viscous)
!! @param[out] flop
!<
    subroutine cbc_vis_ee (vc, sz, g, dh, dt, v00, rei, wv, v, bx, cf, flop)
    implicit none
    include 'cbc_f_params.h'
    integer                                                     ::  i, j, k, ix, jx, kx, g, bvx
    integer                                                     ::  b_e1, b_w1, b_n1, b_s1, b_t1, b_b1
    integer, dimension(3)                                       ::  sz
    real                                                        ::  Up0, Ue1, Uw1, Us1, Un1, Ub1, Ut1
    real                                                        ::  Vp0, Ve1, Vw1, Vs1, Vn1, Vb1, Vt1
    real                                                        ::  Wp0, We1, Ww1, Ws1, Wn1, Wb1, Wt1
    real                                                        ::  u_ref, v_ref, w_ref, dh, dh1, dh2, flop
    real                                                        ::  EX, EY, EZ, rei, cf, dt, uq, vq, wq, actv
    real                                                        ::  c_e, c_w, c_n, c_s, c_t, c_b
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)   ::  v, wv, vc
    real, dimension(0:3)                                        ::  v00
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bx
    
    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    dh1= 1.0/dh
    dh2= dt*rei*dh1*dh1 * cf
    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)
    flop = flop + real(ix*jx*kx)*48.0 + 5.0
    
    do k=1,kx
    do j=1,jx
    do i=1,ix
      Up0 = v(i  ,j  ,k  ,1)
      Uw1 = v(i-1,j  ,k  ,1)
      Ue1 = v(i+1,j  ,k  ,1)
      Us1 = v(i  ,j-1,k  ,1)
      Un1 = v(i  ,j+1,k  ,1)
      Ub1 = v(i  ,j  ,k-1,1)
      Ut1 = v(i  ,j  ,k+1,1)
			
      Vp0 = v(i  ,j  ,k  ,2)
      Vw1 = v(i-1,j  ,k  ,2)
      Ve1 = v(i+1,j  ,k  ,2)
      Vs1 = v(i  ,j-1,k  ,2)
      Vn1 = v(i  ,j+1,k  ,2)
      Vb1 = v(i  ,j  ,k-1,2)
      Vt1 = v(i  ,j  ,k+1,2)

      Wp0 = v(i  ,j  ,k  ,3)
      Ww1 = v(i-1,j  ,k  ,3)
      We1 = v(i+1,j  ,k  ,3)
      Ws1 = v(i  ,j-1,k  ,3)
      Wn1 = v(i  ,j+1,k  ,3)
      Wb1 = v(i  ,j  ,k-1,3)
      Wt1 = v(i  ,j  ,k+1,3)
      
      uq = 2.0*u_ref - Up0
      vq = 2.0*v_ref - Vp0
      wq = 2.0*w_ref - Wp0
        
      ! セルフェイスのマスク関数
      bvx = bx(i,j,k)
      actv= real(ibits(bvx,        State, 1))
      b_w1= ibits(bx(i-1,j  ,k  ), State, 1)
      b_e1= ibits(bx(i+1,j  ,k  ), State, 1)
      b_s1= ibits(bx(i  ,j-1,k  ), State, 1)
      b_n1= ibits(bx(i  ,j+1,k  ), State, 1)
      b_b1= ibits(bx(i  ,j  ,k-1), State, 1)
      b_t1= ibits(bx(i  ,j  ,k+1), State, 1)
      
      ! 各面のVBCフラグ ibits() = 0(Normal) / others(BC) >> c_e = 1.0(Normal) / 0.0(BC) 
      c_e = 1.0
      c_w = 1.0
      c_n = 1.0
      c_s = 1.0
      c_t = 1.0
      c_b = 1.0
      if ( ibits(bvx, bc_face_E, bitw_5) /= 0 ) c_e = 0.0
      if ( ibits(bvx, bc_face_W, bitw_5) /= 0 ) c_w = 0.0
      if ( ibits(bvx, bc_face_N, bitw_5) /= 0 ) c_n = 0.0
      if ( ibits(bvx, bc_face_S, bitw_5) /= 0 ) c_s = 0.0
      if ( ibits(bvx, bc_face_T, bitw_5) /= 0 ) c_t = 0.0
      if ( ibits(bvx, bc_face_B, bitw_5) /= 0 ) c_b = 0.0
			
      if ( b_e1 == 0 ) then
        Ue1 = uq
        Ve1 = vq
        We1 = wq
      endif
      
      if ( b_w1 == 0 ) then
        Uw1 = uq
        Vw1 = vq
        Ww1 = wq
      end if
			
      if ( b_n1 == 0 ) then
        Un1 = uq
        Vn1 = vq
        Wn1 = wq
      end if
      
      if ( b_s1 == 0 ) then
        Us1 = uq
        Vs1 = vq
        Ws1 = wq
      end if
			
      if ( b_t1 == 0 ) then
        Ut1 = uq
        Vt1 = vq
        Wt1 = wq
      end if
      
      if ( b_b1 == 0 ) then
        Ub1 = uq
        Vb1 = vq
        Wb1 = wq
      end if
      
      ! 粘性項　Euler陽解法
      EX = ( Ue1 - Up0 ) * c_e &
         + ( Uw1 - Up0 ) * c_w &
         + ( Un1 - Up0 ) * c_n &
         + ( Us1 - Up0 ) * c_s &
         + ( Ut1 - Up0 ) * c_t &
         + ( Ub1 - Up0 ) * c_b
      EY = ( Ve1 - Vp0 ) * c_e &
         + ( Vw1 - Vp0 ) * c_w &
         + ( Vn1 - Vp0 ) * c_n &
         + ( Vs1 - Vp0 ) * c_s &
         + ( Vt1 - Vp0 ) * c_t &
         + ( Vb1 - Vp0 ) * c_b
      EZ = ( We1 - Wp0 ) * c_e &
         + ( Ww1 - Wp0 ) * c_w &
         + ( Wn1 - Wp0 ) * c_n &
         + ( Ws1 - Wp0 ) * c_s &
         + ( Wt1 - Wp0 ) * c_t &
         + ( Wb1 - Wp0 ) * c_b
			
      vc(i,j,k,1) = ( wv(i,j,k,1) + EX*dh2 )
      vc(i,j,k,2) = ( wv(i,j,k,2) + EY*dh2 )
      vc(i,j,k,3) = ( wv(i,j,k,3) + EZ*dh2 )
    end do
    end do
    end do

    return
    end subroutine cbc_vis_ee

!  *******************************************************************************************
!> @subroutine cbc_vis_ee_vbc (vc, sz, g, st, ed, dh, dt, v00, rei, v, bx, odr, cf, vec, flop)
!! @brief 速度境界条件による粘性項の修正（Euler陽解法）
!! @param[out] vc 疑似ベクトル
!! @param sz 配列長
!! @param g ガイドセル長
!! @param st ループの開始インデクス
!! @param ed ループの終了インデクス
!! @param dh 格子幅
!! @param dt 時間積分幅
!! @param v00 参照速度
!! @param rei 1/Re
!! @param v 速度ベクトル（n-step, collocated）
!! @param bx BC index for v
!! @param odr 内部境界処理時の速度境界条件のエントリ
!! @param cf 粘性項の係数
!! @param vec 指定する速度ベクトル
!! @param[out] flop
!<
    subroutine cbc_vis_ee_vbc (vc, sz, g, st, ed, dh, dt, v00, rei, v, bx, odr, cf, vec, flop)
    implicit none
    include 'cbc_f_params.h'
    integer                                                     ::  i, j, k, g, bvx, m, odr
    integer, dimension(3)                                       ::  sz, st, ed
    real                                                        ::  Up0, Ue1, Uw1, Us1, Un1, Ub1, Ut1
    real                                                        ::  Vp0, Ve1, Vw1, Vs1, Vn1, Vb1, Vt1
    real                                                        ::  Wp0, We1, Ww1, Ws1, Wn1, Wb1, Wt1
    real                                                        ::  u_ref, v_ref, w_ref, dh, dh1, dh2, flop
    real                                                        ::  EX, EY, EZ, rei, cf, dt
    real                                                        ::  u_bc_ref, v_bc_ref, w_bc_ref
    real                                                        ::  c_e, c_w, c_n, c_s, c_t, c_b
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)   ::  vc, v
    real, dimension(0:3)                                        ::  v00
    real, dimension(3)                                          ::  vec
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bx
    
    dh1= 1.0/dh
    dh2= dt*rei*dh1*dh1 * cf
    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)
    u_bc_ref = vec(1) + u_ref
    v_bc_ref = vec(2) + v_ref
    w_bc_ref = vec(3) + w_ref
    m = 0
    
    Uw1 = 0.0
    Ue1 = 0.0
    Us1 = 0.0
    Un1 = 0.0
    Ub1 = 0.0
    Ut1 = 0.0
    Vw1 = 0.0
    Ve1 = 0.0
    Vs1 = 0.0
    Vn1 = 0.0
    Vb1 = 0.0
    Vt1 = 0.0
    Ww1 = 0.0
    We1 = 0.0
    Ws1 = 0.0
    Wn1 = 0.0
    Wb1 = 0.0
    Wt1 = 0.0
    
    do k=st(3),ed(3)
    do j=st(2),ed(2)
    do i=st(1),ed(1)
      bvx = bx(i,j,k)
      if ( 0 /= iand(bvx, bc_mask30) ) then ! 6面のうちのどれか速度境界フラグが立っている場合
        
        Up0 = v(i  ,j  ,k  ,1)
        Vp0 = v(i  ,j  ,k  ,2)
        Wp0 = v(i  ,j  ,k  ,3)
        
        ! 内部境界のときの各面のVBCフラグ ibits() = 0(Normal) / others(BC) >> c_e = 0.0(Normal) / 1.0(BC) 
        c_e = 0.0
        c_w = 0.0
        c_n = 0.0
        c_s = 0.0
        c_t = 0.0
        c_b = 0.0
        
        if ( ibits(bvx, bc_face_E, bitw_5) == odr ) then
          c_e = 1.0
          Ue1 = u_bc_ref
          Ve1 = v_bc_ref
          We1 = w_bc_ref
        endif
        
        if ( ibits(bvx, bc_face_W, bitw_5) == odr ) then
          c_w = 1.0
          Uw1 = u_bc_ref
          Vw1 = v_bc_ref
          Ww1 = w_bc_ref
        endif
        
        if ( ibits(bvx, bc_face_N, bitw_5) == odr ) then
          c_n = 1.0
          Un1 = u_bc_ref
          Vn1 = v_bc_ref
          Wn1 = w_bc_ref
        endif
        
        if ( ibits(bvx, bc_face_S, bitw_5) == odr ) then
          c_s = 1.0
          Us1 = u_bc_ref
          Vs1 = v_bc_ref
          Ws1 = w_bc_ref
        endif
        
        if ( ibits(bvx, bc_face_T, bitw_5) == odr ) then
          c_t = 1.0
          Ut1 = u_bc_ref
          Vt1 = v_bc_ref
          Wt1 = w_bc_ref
        endif
        
        if ( ibits(bvx, bc_face_B, bitw_5) == odr ) then
          c_b = 1.0
          Ub1 = u_bc_ref
          Vb1 = v_bc_ref
          Wb1 = w_bc_ref
        endif
        
        EX = ( Ue1 - Up0 ) * c_e &
           + ( Uw1 - Up0 ) * c_w &
           + ( Un1 - Up0 ) * c_n &
           + ( Us1 - Up0 ) * c_s &
           + ( Ut1 - Up0 ) * c_t &
           + ( Ub1 - Up0 ) * c_b
        EY = ( Ve1 - Vp0 ) * c_e &
           + ( Vw1 - Vp0 ) * c_w &
           + ( Vn1 - Vp0 ) * c_n &
           + ( Vs1 - Vp0 ) * c_s &
           + ( Vt1 - Vp0 ) * c_t &
           + ( Vb1 - Vp0 ) * c_b
        EZ = ( We1 - Wp0 ) * c_e &
           + ( Ww1 - Wp0 ) * c_w &
           + ( Wn1 - Wp0 ) * c_n &
           + ( Ws1 - Wp0 ) * c_s &
           + ( Wt1 - Wp0 ) * c_t &
           + ( Wb1 - Wp0 ) * c_b
			
        vc(i,j,k,1) = vc(i,j,k,1) + EX*dh2
        vc(i,j,k,2) = vc(i,j,k,2) + EY*dh2
        vc(i,j,k,3) = vc(i,j,k,3) + EZ*dh2
        
        m = m+1
      end if
    end do
    end do
    end do
    
    flop = flop + real(8 + 42*m) ! 20100628

    return
    end subroutine cbc_vis_ee_vbc
    
!  ******************************************************************************
!> @subroutine cbc_vis_cn_sor (v, sz, g, dh, v00, rei, omg, vc, bx, cf, dv, flop)
!! @brief SOR法による粘性項の計算
!! @param[out] v 疑似ベクトル
!! @param sz 配列長
!! @param g ガイドセル長
!! @param dh 格子幅
!! @param dt 時間積分幅
!! @param v00 参照速度
!! @param rei 1/Re
!! @param omg SOR法の加速係数（1<omg<2）
!! @param vc 疑似ベクトル（n-step, collocated）
!! @param bx BC index for V
!! @param cf 粘性項の係数 (0.5=CN, 1.0=Euler Implicit)
!! @param dv 残差
!! @param[out] flop
!<
    subroutine cbc_vis_cn_sor (v, sz, g, dh, dt, v00, rei, omg, vc, bx, cf, dv, flop)
    implicit none
    include 'cbc_f_params.h'
    integer                                                     ::  i, j, k, ix, jx, kx, g, bvx
    integer                                                     ::  b_e1, b_w1, b_n1, b_s1, b_t1, b_b1
    integer, dimension(3)                                       ::  sz
    real                                                        ::  Up0, Ue1, Uw1, Us1, Un1, Ub1, Ut1
    real                                                        ::  Vp0, Ve1, Vw1, Vs1, Vn1, Vb1, Vt1
    real                                                        ::  Wp0, We1, Ww1, Ws1, Wn1, Wb1, Wt1
    real                                                        ::  u_ref, v_ref, w_ref, dh, dh1, dh2, flop
    real                                                        ::  rei, dd, dv, dv1, dv2, dv3, uq, vq, wq, ddv
    real                                                        ::  omg, s1, s2, s3, cf, dt, actv, b_vbc
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)   ::  v, vc
    real, dimension(0:3)                                        ::  v00
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bx
    
    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    dh1= 1.0/dh
    dh2= dt*rei*dh1*dh1*cf
    dd = 1.0 / (1.0 + 6.0*dh2)
    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)
    flop = flop + real(ix*jx*kx)*52.0 + 8.0
    
    do k=1,kx
    do j=1,jx
    do i=1,ix
      Up0 = v(i  ,j  ,k  ,1)
      Uw1 = v(i-1,j  ,k  ,1)
      Ue1 = v(i+1,j  ,k  ,1)
      Us1 = v(i  ,j-1,k  ,1)
      Un1 = v(i  ,j+1,k  ,1)
      Ub1 = v(i  ,j  ,k-1,1)
      Ut1 = v(i  ,j  ,k+1,1)
			
      Vp0 = v(i  ,j  ,k  ,2)
      Vw1 = v(i-1,j  ,k  ,2)
      Ve1 = v(i+1,j  ,k  ,2)
      Vs1 = v(i  ,j-1,k  ,2)
      Vn1 = v(i  ,j+1,k  ,2)
      Vb1 = v(i  ,j  ,k-1,2)
      Vt1 = v(i  ,j  ,k+1,2)

      Wp0 = v(i  ,j  ,k  ,3)
      Ww1 = v(i-1,j  ,k  ,3)
      We1 = v(i+1,j  ,k  ,3)
      Ws1 = v(i  ,j-1,k  ,3)
      Wn1 = v(i  ,j+1,k  ,3)
      Wb1 = v(i  ,j  ,k-1,3)
      Wt1 = v(i  ,j  ,k+1,3)
      
      uq = 2.0*u_ref - Up0
      vq = 2.0*v_ref - Vp0
      wq = 2.0*w_ref - Wp0
      
      ! セルフェイスのマスク関数
      bvx = bx(i,j,k)
      actv= real(ibits(bvx,        State, 1))
      b_w1= ibits(bx(i-1,j  ,k  ), State, 1)
      b_e1= ibits(bx(i+1,j  ,k  ), State, 1)
      b_s1= ibits(bx(i  ,j-1,k  ), State, 1)
      b_n1= ibits(bx(i  ,j+1,k  ), State, 1)
      b_b1= ibits(bx(i  ,j  ,k-1), State, 1)
      b_t1= ibits(bx(i  ,j  ,k+1), State, 1)
      
      b_vbc = 1.0 ! dummy
			
      if ( b_e1 == 0 ) then
        Ue1 = uq
        Ve1 = vq
        We1 = wq
      endif
      
      if ( b_w1 == 0 ) then
        Uw1 = uq
        Vw1 = vq
        Ww1 = wq
      end if
			
      if ( b_n1 == 0 ) then
        Un1 = uq
        Vn1 = vq
        Wn1 = wq
      end if
      
      if ( b_s1 == 0 ) then
        Us1 = uq
        Vs1 = vq
        Ws1 = wq
      end if
			
      if ( b_t1 == 0 ) then
        Ut1 = uq
        Vt1 = vq
        Wt1 = wq
      end if
      
      if ( b_b1 == 0 ) then
        Ub1 = uq
        Vb1 = vq
        Wb1 = wq
      end if
      
      s1 = (Ue1 + Uw1 + Un1 + Us1 + Ut1 + Ub1) * dh2 + vc(i,j,k,1)
      s2 = (Ve1 + Vw1 + Vn1 + Vs1 + Vt1 + Vb1) * dh2 + vc(i,j,k,2)
      s3 = (We1 + Ww1 + Wn1 + Ws1 + Wt1 + Wb1) * dh2 + vc(i,j,k,3)
      dv1= (dd*s1-Up0)*actv * b_vbc ! to exclude vbc cell
      dv2= (dd*s2-Vp0)*actv * b_vbc
      dv3= (dd*s3-Wp0)*actv * b_vbc
      v(i,j,k,1) = Up0 + omg*dv1
      v(i,j,k,2) = Vp0 + omg*dv2
      v(i,j,k,3) = Wp0 + omg*dv3
      ddv = dv1*dv1 + dv2*dv2 + dv3*dv3
      dv = max( dv, sqrt(ddv) )
    end do
    end do
    end do

    return
    end subroutine cbc_vis_cn_sor

!  ***************************************************************************************************
!> @subroutine cbc_vis_cn_mod_sor (v, sz, g, st, ed, dh, dt, v00, rei, omg, vc, bx, cf, dv, vec, flop)
!! @brief SOR法による粘性項計算の境界条件による修正
!! @param[in/out] v 疑似ベクトル
!! @param sz 配列長
!! @param g ガイドセル長
!! @param st ループの開始インデクス
!! @param ed ループの終了インデクス
!! @param dh 格子幅
!! @param dt 時間積分幅
!! @param v00 参照速度
!! @param rei 1/Re
!! @param omg SOR法の加速係数（1<omg<2）
!! @param vc 疑似ベクトル（n-step, collocated）
!! @param bx BC index for V
!! @param cf 粘性項の係数 (0.5=CN, 1.0=Euler Implicit)
!! @param[out] dv 残差
!! @param vec 指定する速度ベクトル
!! @param[out] flop
!<
    subroutine cbc_vis_cn_mod_sor (v, sz, g, st, ed, dh, dt, v00, rei, omg, vc, bx, cf, dv, vec, flop)
    implicit none
    include 'cbc_f_params.h'
    integer                                                     ::  i, j, k, g, bvx, m
    integer                                                     ::  b_e1, b_w1, b_n1, b_s1, b_t1, b_b1
    integer, dimension(3)                                       ::  sz, st, ed
    real                                                        ::  Up0, Ue1, Uw1, Us1, Un1, Ub1, Ut1
    real                                                        ::  Vp0, Ve1, Vw1, Vs1, Vn1, Vb1, Vt1
    real                                                        ::  Wp0, We1, Ww1, Ws1, Wn1, Wb1, Wt1
    real                                                        ::  u_ref, v_ref, w_ref, dh, dh1, dh2, flop
    real                                                        ::  rei, dd, dv, dv1, dv2, dv3, ddv
    real                                                        ::  u_bc, v_bc, w_bc, uq, vq, wq
    real                                                        ::  omg, s1, s2, s3, cf, dt, actv
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)   ::  v, vc
    real, dimension(0:3)                                        ::  v00
    real, dimension(3)                                          ::  vec
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bx
    
    dh1= 1.0/dh
    dh2= dt*rei*dh1*dh1*cf
    dd = 1.0 / (1.0 + 6.0*dh2)
    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)
    u_bc = vec(1) + u_ref
    v_bc = vec(2) + v_ref
    w_bc = vec(3) + w_ref
    m = 0

    do k=st(3),ed(3)
    do j=st(2),ed(2)
    do i=st(1),ed(1)
      bvx = bx(i,j,k)
      if ( 0 /= iand(bvx, bc_mask30) ) then ! 6面のうちのどれか速度境界フラグが立っている場合
        m = m+1
        Up0 = v(i  ,j  ,k  ,1)
        Uw1 = v(i-1,j  ,k  ,1)
        Ue1 = v(i+1,j  ,k  ,1)
        Us1 = v(i  ,j-1,k  ,1)
        Un1 = v(i  ,j+1,k  ,1)
        Ub1 = v(i  ,j  ,k-1,1)
        Ut1 = v(i  ,j  ,k+1,1)
			
        Vp0 = v(i  ,j  ,k  ,2)
        Vw1 = v(i-1,j  ,k  ,2)
        Ve1 = v(i+1,j  ,k  ,2)
        Vs1 = v(i  ,j-1,k  ,2)
        Vn1 = v(i  ,j+1,k  ,2)
        Vb1 = v(i  ,j  ,k-1,2)
        Vt1 = v(i  ,j  ,k+1,2)

        Wp0 = v(i  ,j  ,k  ,3)
        Ww1 = v(i-1,j  ,k  ,3)
        We1 = v(i+1,j  ,k  ,3)
        Ws1 = v(i  ,j-1,k  ,3)
        Wn1 = v(i  ,j+1,k  ,3)
        Wb1 = v(i  ,j  ,k-1,3)
        Wt1 = v(i  ,j  ,k+1,3)
        
        uq = 2.0*u_ref - Up0
        vq = 2.0*v_ref - Vp0
        wq = 2.0*w_ref - Wp0
      
        actv= real(ibits(bvx,        State, 1))
        b_w1= ibits(bx(i-1,j  ,k  ), State, 1)
        b_e1= ibits(bx(i+1,j  ,k  ), State, 1)
        b_s1= ibits(bx(i  ,j-1,k  ), State, 1)
        b_n1= ibits(bx(i  ,j+1,k  ), State, 1)
        b_b1= ibits(bx(i  ,j  ,k-1), State, 1)
        b_t1= ibits(bx(i  ,j  ,k+1), State, 1)
			
        if ( b_e1 == 0 ) then
          Ue1 = uq
          Ve1 = vq
          We1 = wq
        endif
      
        if ( b_w1 == 0 ) then
          Uw1 = uq
          Vw1 = vq
          Ww1 = wq
        end if
			
        if ( b_n1 == 0 ) then
          Un1 = uq
          Vn1 = vq
          Wn1 = wq
        end if
      
        if ( b_s1 == 0 ) then
          Us1 = uq
          Vs1 = vq
          Ws1 = wq
        end if
			
        if ( b_t1 == 0 ) then
          Ut1 = uq
          Vt1 = vq
          Wt1 = wq
        end if
      
        if ( b_b1 == 0 ) then
          Ub1 = uq
          Vb1 = vq
          Wb1 = wq
        end if
        
        if ( ibits(bvx, bc_face_E, 1) == 0 ) then
          Ue1 = u_bc
          Ve1 = v_bc
          We1 = w_bc
        end if
        
        if ( ibits(bvx, bc_face_W, 1) == 0 ) then
          Uw1 = u_bc
          Vw1 = v_bc
          Ww1 = w_bc
        end if
        
        if ( ibits(bvx, bc_face_N, 1) == 0 ) then
          Un1 = u_bc
          Vn1 = v_bc
          Wn1 = w_bc
        end if
        
        if ( ibits(bvx, bc_face_S, 1) == 0 ) then
          Us1 = u_bc
          Vs1 = v_bc
          Ws1 = w_bc
        end if
        
        if ( ibits(bvx, bc_face_T, 1) == 0 ) then
          Ut1 = u_bc
          Vt1 = v_bc
          Wt1 = w_bc
        end if
        
        if ( ibits(bvx, bc_face_B, 1) == 0 ) then
          Ub1 = u_bc
          Vb1 = v_bc
          Wb1 = w_bc
        end if
      
        s1 = (Ue1 + Uw1 + Un1 + Us1 + Ut1 + Ub1) * dh2 + vc(i,j,k,1)
        s2 = (Ve1 + Vw1 + Vn1 + Vs1 + Vt1 + Vb1) * dh2 + vc(i,j,k,2)
        s3 = (We1 + Ww1 + Wn1 + Ws1 + Wt1 + Wb1) * dh2 + vc(i,j,k,3)
        dv1= (dd*s1-Up0)*actv
        dv2= (dd*s2-Vp0)*actv
        dv3= (dd*s3-Wp0)*actv
        v(i,j,k,1) = Up0 + omg*dv1
        v(i,j,k,2) = Vp0 + omg*dv2
        v(i,j,k,3) = Wp0 + omg*dv3

        ddv = dv1*dv1 + dv2*dv2 + dv3*dv3
        dv = max( dv, sqrt(ddv) )
      end if
    end do
    end do
    end do
    
    flop = flop + real(m)*48.0 + 8.0

    return
    end subroutine cbc_vis_cn_mod_sor

!  **************************************************************************************
!> @subroutine cbc_eddy_viscosity (vt, sz, g, dh, re, cs, v, bx, vt_range, yp_range, v00)
!! @brief Smagorinsky Modelによる乱流渦粘性係数の計算，減衰関数を併用
!! @param vt 乱流渦粘性係数
!! @param sz 配列長
!! @param g ガイドセル長
!! @param dh 格子幅
!! @param re レイノルズ数
!! @param cs 須磨五輪スキー定数
!! @param v 速度ベクトル
!! @param bx BC index for V
!! @param vt_range 渦粘性係数の最小値と最大値
!! @param yp_range Y+の最小値と最大値
!! @param v00 参照速度
!! @note
!!    - vtmin, vtmax > vt_range(2)
!!    - ypmin, ypmax > yp_range(2)
!! @todo
!!    - 境界条件は必要か？
!<
    subroutine cbc_eddy_viscosity (vt, sz, g, dh, re, cs, v, bx, vt_range, yp_range, v00)
    implicit none
    include 'cbc_f_params.h'
    integer                                                     ::  i, j, k, ix, jx, kx, g, m
    integer, dimension(3)                                       ::  sz
    real, dimension(2)                                          ::  vt_range, yp_range
    real                                                        ::  dh, re, cs, yp, Ut, Ut0
    real                                                        ::  dudx, dudy, dudz
    real                                                        ::  dvdx, dvdy, dvdz
    real                                                        ::  dwdx, dwdy, dwdz
    real                                                        ::  d1, d2, ddd, c1, c2, c3, delta, dd
    real                                                        ::  fs, aaa, Vmag, dis, tw, up1
    real                                                        ::  u_ref, v_ref, w_ref, u1, u2, u3
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  vt
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)   ::  v
    real, dimension(0:3)                                        ::  v00
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bx

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    vt_range(1) =  1.0e6
    vt_range(2) = -1.0e6
    yp_range(1) =  1.0e6
    yp_range(2) = -1.0e6
    Delta = (dh*dh*dh)**(0.333)
    dd    = 0.5/dh
    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)

    do k=1,kx
    do j=1,jx
    do i=1,ix
      DUDX=(v(i+1,j,k,1)-v(i-1,j,k,1)) * dd
      DUDY=(v(i,j+1,k,1)-v(i,j-1,k,1)) * dd
      DUDZ=(v(i,j,k+1,1)-v(i,j,k-1,1)) * dd

      DVDX=(v(i+1,j,k,2)-v(i-1,j,k,2)) * dd
      DVDY=(v(i,j+1,k,2)-v(i,j-1,k,2)) * dd
      DVDZ=(v(i,j,k+1,2)-v(i,j,k-1,2)) * dd

      DWDX=(v(i+1,j,k,3)-v(i-1,j,k,3)) * dd
      DWDY=(v(i,j+1,k,3)-v(i,j-1,k,3)) * dd
      DWDZ=(v(i,j,k+1,3)-v(i,j,k-1,3)) * dd

      D1  = DUDX*DUDX + DVDY*DVDY + DWDZ*DWDZ
      c1  = DUDY+DVDX
      c2  = DVDZ+DWDY
      c3  = DWDX+DUDZ
      D2  = c1*c1 + c2*c2 + c3*c3
      DDD = SQRT(2.0D0*D1 + D2)

!     ------- Y+ & fs--
      fs =1.0
      !AAA= bx(i,j,k) * bx(i,j,k)   &
      !    *bx(i,j,k) * bx(i,j,k)   &
      !    *bx(i,j,k) * bx(i,j,k)
      u1 = v(i,j,k,1) - u_ref
      u2 = v(i,j,k,2) - v_ref
      u3 = v(i,j,k,3) - w_ref
      Vmag = SQRT(u1*u1 + u2*u2 + u3*u3)

      IF ( (AAA < 1.0) .AND. (Vmag >= 0.001) ) THEN
        !DIS = AMIN1(bnd(i,j,k), bnd(i,j,k),    &
        !            bnd(i,j,k), bnd(i,j,k),    &
        !            bnd(i,j,k), bnd(i,j,k) ) * dh
        Tw = Vmag/(DIS*RE)
        Ut0= SQRT(Tw)
        YP = DIS*Ut0*RE
        Up1= Vmag
        
        DO M=1,5
          Ut = Ut0 + (Up1-Ut0*(ALOG(YP)/0.4+5.5)) / (Up1/Ut0+1.0/0.4)
          IF (ABS(Ut-Ut0) <= 1.0E-6) exit
          IF (Ut <= 0.010) Ut=0.0100
          Ut0=Ut
          YP =DIS*Ut0*RE
        end do
        
        fs = 1.0-exp(-YP/25.0)
        yp_range(1) = amin1(yp_range(1), YP)
        yp_range(2) = amax1(yp_range(2), YP)
      ENDIF

      VT(i,j,k) = ((Cs*fs*Delta)**2)*DDD*RE
      vt_range(1) = MIN(vt_range(1), VT(i,j,k))
      vt_range(2) = MAX(vt_range(2), VT(i,j,k))
    end do
    end do
    end do

    do k=1,kx
    do i=1,ix
      VT(i,0   ,k)=VT(i,1 ,k)
      VT(i,jx+1,k)=VT(i,jx,k)
    end do
    end do

    do k=1,kx
    do j=1,jx
      VT(0   ,j,k)=VT(1 ,j,k)
      VT(ix+1,j,k)=VT(ix,j,k)
    end do
    end do
 
    do j=1,jx
    do i=1,ix
      VT(i,j,0   )=VT(i,j,1 )
      VT(i,j,kx+1)=VT(i,j,kx)
    end do
    end do

    return
    end subroutine cbc_eddy_viscosity

!  *********************************************************
!> @subroutine cbc_ab2 (vc, sz, g, dt, v, ab, bd, v00, flop)
!! @brief Adams-Bashforth法による疑似ベクトルの時間積分
!! @param[out] vc 疑似ベクトル
!! @param sz 配列長
!! @param g ガイドセル長
!! @param dt 時間積分幅
!! @param v 速度ベクトル（n-step, collocated）
!! @param ab 前ステップの対流項（＋粘性項）の計算値
!! @param bd BCindex ID
!! @param v00 参照速度
!! @param[out] flop
!<
    subroutine cbc_ab2 (vc, sz, g, dt, v, ab, bd, v00, flop)
    implicit none
    include 'cbc_f_params.h'
    integer                                                     ::  i, j, k, ix, jx, kx, g
    integer, dimension(3)                                       ::  sz
    real                                                        ::  flop, actv, dt, ab_u, ab_v, ab_w, u_ref, v_ref, w_ref
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)   ::  vc, v, ab
    real, dimension(0:3)                                        ::  v00
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bd
    
    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)
    flop = flop + real(27*ix*jx*kx) ! 20100706
    
    do k=1,kx
    do j=1,jx
    do i=1,ix
      actv = real(ibits(bd(i,j,k), State, 1))
      
      ab_u = ab(i,j,k,1)
      ab_v = ab(i,j,k,2)
      ab_w = ab(i,j,k,3)
      
      ab(i,j,k,1) = vc(i,j,k,1)
      ab(i,j,k,2) = vc(i,j,k,2)
      ab(i,j,k,3) = vc(i,j,k,3)

      vc(i,j,k,1) = ( v(i,j,k,1) + dt*0.5*( 3.0*vc(i,j,k,1)-ab_u ) )* actv + (1.0-actv)*u_ref
      vc(i,j,k,2) = ( v(i,j,k,2) + dt*0.5*( 3.0*vc(i,j,k,2)-ab_v ) )* actv + (1.0-actv)*v_ref
      vc(i,j,k,3) = ( v(i,j,k,3) + dt*0.5*( 3.0*vc(i,j,k,3)-ab_w ) )* actv + (1.0-actv)*w_ref
    end do
    end do
    end do

    return
    end subroutine cbc_ab2
    
!  **********************************************************************************
!> @subroutine cbc_vis_cn_jcb (v, sz, g, dh, v00, rei, omg, vc, bx, wk, cf, dv, flop)
!! @brief 緩和Jacobi法による粘性項の計算
!! @param[out] v 疑似ベクトル
!! @param sz 配列長
!! @param g ガイドセル長
!! @param dh 格子幅
!! @param dt 時間積分幅
!! @param v00 参照速度
!! @param rei 1/Re
!! @param omg 緩和Jacobi法の加速係数（0<omg<1）
!! @param vc 疑似ベクトル（n-step, collocated）
!! @param bx BC index for V
!! @param wk ワーク用配列
!! @param cf 粘性項の係数 (0.5=CN, 1.0=Euler Implicit)
!! @param[out] dv 残差
!! @param[out] flop
!<
    subroutine cbc_vis_cn_jcb (v, sz, g, dh, dt, v00, rei, omg, vc, bx, wk, cf, dv, flop)
    implicit none
    include 'cbc_f_params.h'
    integer                                                     ::  i, j, k, ix, jx, kx, g, bvx
    integer                                                     ::  b_e1, b_w1, b_n1, b_s1, b_t1, b_b1
    integer, dimension(3)                                       ::  sz
    real                                                        ::  Up0, Ue1, Uw1, Us1, Un1, Ub1, Ut1
    real                                                        ::  Vp0, Ve1, Vw1, Vs1, Vn1, Vb1, Vt1
    real                                                        ::  Wp0, We1, Ww1, Ws1, Wn1, Wb1, Wt1
    real                                                        ::  u_ref, v_ref, w_ref, dh, dh1, dh2, flop
    real                                                        ::  rei, dd, dv, dv1, dv2, dv3, uq, vq, wq, ddv
    real                                                        ::  omg, s1, s2, s3, cf, dt, actv, b_vbc
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)   ::  v, vc, wk
    real, dimension(0:3)                                        ::  v00
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bx
    
    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    dh1= 1.0/dh
    dh2= dt*rei*dh1*dh1*cf
    dd = 1.0 / (1.0 + 6.0*dh2)
    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)
    flop = flop + real(8 + 52*ix*jx*kx)
    
    do k=1,kx
    do j=1,jx
    do i=1,ix
      Up0 = v(i  ,j  ,k  ,1)
      Uw1 = v(i-1,j  ,k  ,1)
      Ue1 = v(i+1,j  ,k  ,1)
      Us1 = v(i  ,j-1,k  ,1)
      Un1 = v(i  ,j+1,k  ,1)
      Ub1 = v(i  ,j  ,k-1,1)
      Ut1 = v(i  ,j  ,k+1,1)
			
      Vp0 = v(i  ,j  ,k  ,2)
      Vw1 = v(i-1,j  ,k  ,2)
      Ve1 = v(i+1,j  ,k  ,2)
      Vs1 = v(i  ,j-1,k  ,2)
      Vn1 = v(i  ,j+1,k  ,2)
      Vb1 = v(i  ,j  ,k-1,2)
      Vt1 = v(i  ,j  ,k+1,2)

      Wp0 = v(i  ,j  ,k  ,3)
      Ww1 = v(i-1,j  ,k  ,3)
      We1 = v(i+1,j  ,k  ,3)
      Ws1 = v(i  ,j-1,k  ,3)
      Wn1 = v(i  ,j+1,k  ,3)
      Wb1 = v(i  ,j  ,k-1,3)
      Wt1 = v(i  ,j  ,k+1,3)
      
      uq = 2.0*u_ref - Up0
      vq = 2.0*v_ref - Vp0
      wq = 2.0*w_ref - Wp0
      
      ! セルフェイスのマスク関数
      bvx = bx(i,j,k)
      actv= real(ibits(bvx,        State, 1))
      b_w1= ibits(bx(i-1,j  ,k  ), State, 1)
      b_e1= ibits(bx(i+1,j  ,k  ), State, 1)
      b_s1= ibits(bx(i  ,j-1,k  ), State, 1)
      b_n1= ibits(bx(i  ,j+1,k  ), State, 1)
      b_b1= ibits(bx(i  ,j  ,k-1), State, 1)
      b_t1= ibits(bx(i  ,j  ,k+1), State, 1)
      
      b_vbc = 1.0 ! dummy
			
      if ( b_e1 == 0 ) then
        Ue1 = uq
        Ve1 = vq
        We1 = wq
      endif
      
      if ( b_w1 == 0 ) then
        Uw1 = uq
        Vw1 = vq
        Ww1 = wq
      end if
			
      if ( b_n1 == 0 ) then
        Un1 = uq
        Vn1 = vq
        Wn1 = wq
      end if
      
      if ( b_s1 == 0 ) then
        Us1 = uq
        Vs1 = vq
        Ws1 = wq
      end if
			
      if ( b_t1 == 0 ) then
        Ut1 = uq
        Vt1 = vq
        Wt1 = wq
      end if
      
      if ( b_b1 == 0 ) then
        Ub1 = uq
        Vb1 = vq
        Wb1 = wq
      end if
      
      s1 = (Ue1 + Uw1 + Un1 + Us1 + Ut1 + Ub1) * dh2 + vc(i,j,k,1)
      s2 = (Ve1 + Vw1 + Vn1 + Vs1 + Vt1 + Vb1) * dh2 + vc(i,j,k,2)
      s3 = (We1 + Ww1 + Wn1 + Ws1 + Wt1 + Wb1) * dh2 + vc(i,j,k,3)
      dv1= (dd*s1-Up0)*actv * b_vbc ! to exclude vbc cell
      dv2= (dd*s2-Vp0)*actv * b_vbc
      dv3= (dd*s3-Wp0)*actv * b_vbc
      wk(i,j,k,1) = Up0 + omg*dv1
      wk(i,j,k,2) = Vp0 + omg*dv2
      wk(i,j,k,3) = Wp0 + omg*dv3
      ddv = dv1*dv1 + dv2*dv2 + dv3*dv3
      dv = max( dv, sqrt(ddv) )
    end do
    end do
    end do
    
! update should be done after modification

    return
    end subroutine cbc_vis_cn_jcb

!  *******************************************************************************************************
!> @subroutine cbc_vis_cn_mod_jcb (v, sz, g, st, ed, dh, dt, v00, rei, omg, vc, bx, wk, cf, dv, vec, flop)
!! @brief 緩和Jacobi法による粘性項の計算
!! @param[in/out] v 疑似ベクトル
!! @param sz 配列長
!! @param g ガイドセル長
!! @param st ループの開始インデクス
!! @param ed ループの終了インデクス
!! @param dh 格子幅
!! @param dt 時間積分幅
!! @param v00 参照速度
!! @param rei 1/Re
!! @param omg 緩和Jacobi法の加速係数（0<omg<1）
!! @param vc 疑似ベクトル（n-step, collocated）
!! @param bx BC index for V
!! @param wk ワーク用配列
!! @param cf 粘性項の係数 (0.5=CN, 1.0=Euler Implicit)
!! @param[out] dv 残差
!! @param vec 指定する速度ベクトル
!! @param[out] flop
!<
    subroutine cbc_vis_cn_mod_jcb (v, sz, g, st, ed, dh, dt, v00, rei, omg, vc, bx, wk, cf, dv, vec, flop)
    implicit none
    include 'cbc_f_params.h'
    integer                                                     ::  i, j, k, g, bvx, m
    integer                                                     ::  b_e1, b_w1, b_n1, b_s1, b_t1, b_b1
    integer, dimension(3)                                       ::  sz, st, ed
    real                                                        ::  Up0, Ue1, Uw1, Us1, Un1, Ub1, Ut1
    real                                                        ::  Vp0, Ve1, Vw1, Vs1, Vn1, Vb1, Vt1
    real                                                        ::  Wp0, We1, Ww1, Ws1, Wn1, Wb1, Wt1
    real                                                        ::  u_ref, v_ref, w_ref, dh, dh1, dh2, flop
    real                                                        ::  rei, dd, dv, dv1, dv2, dv3, ddv
    real                                                        ::  u_bc, v_bc, w_bc, uq, vq, wq
    real                                                        ::  omg, s1, s2, s3, cf, dt, actv
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)   ::  v, vc, wk
    real, dimension(0:3)                                        ::  v00
    real, dimension(3)                                          ::  vec
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bx
    
    dh1= 1.0/dh
    dh2= dt*rei*dh1*dh1*cf
    dd = 1.0 / (1.0 + 6.0*dh2)
    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)
    u_bc = vec(1) + u_ref
    v_bc = vec(2) + v_ref
    w_bc = vec(3) + w_ref
    m = 0

    do k=st(3),ed(3)
    do j=st(2),ed(2)
    do i=st(1),ed(1)
      bvx = bx(i,j,k)
      if ( 0 /= iand(bvx, bc_mask30) ) then ! 6面のうちのどれか速度境界フラグが立っている場合
        m = m+1
        Up0 = v(i  ,j  ,k  ,1)
        Uw1 = v(i-1,j  ,k  ,1)
        Ue1 = v(i+1,j  ,k  ,1)
        Us1 = v(i  ,j-1,k  ,1)
        Un1 = v(i  ,j+1,k  ,1)
        Ub1 = v(i  ,j  ,k-1,1)
        Ut1 = v(i  ,j  ,k+1,1)
			
        Vp0 = v(i  ,j  ,k  ,2)
        Vw1 = v(i-1,j  ,k  ,2)
        Ve1 = v(i+1,j  ,k  ,2)
        Vs1 = v(i  ,j-1,k  ,2)
        Vn1 = v(i  ,j+1,k  ,2)
        Vb1 = v(i  ,j  ,k-1,2)
        Vt1 = v(i  ,j  ,k+1,2)

        Wp0 = v(i  ,j  ,k  ,3)
        Ww1 = v(i-1,j  ,k  ,3)
        We1 = v(i+1,j  ,k  ,3)
        Ws1 = v(i  ,j-1,k  ,3)
        Wn1 = v(i  ,j+1,k  ,3)
        Wb1 = v(i  ,j  ,k-1,3)
        Wt1 = v(i  ,j  ,k+1,3)
        
        uq = 2.0*u_ref - Up0
        vq = 2.0*v_ref - Vp0
        wq = 2.0*w_ref - Wp0
      
        actv= real(ibits(bvx,        State, 1))
        b_w1= ibits(bx(i-1,j  ,k  ), State, 1)
        b_e1= ibits(bx(i+1,j  ,k  ), State, 1)
        b_s1= ibits(bx(i  ,j-1,k  ), State, 1)
        b_n1= ibits(bx(i  ,j+1,k  ), State, 1)
        b_b1= ibits(bx(i  ,j  ,k-1), State, 1)
        b_t1= ibits(bx(i  ,j  ,k+1), State, 1)
			
        if ( b_e1 == 0 ) then
          Ue1 = uq
          Ve1 = vq
          We1 = wq
        endif
      
        if ( b_w1 == 0 ) then
          Uw1 = uq
          Vw1 = vq
          Ww1 = wq
        end if
			
        if ( b_n1 == 0 ) then
          Un1 = uq
          Vn1 = vq
          Wn1 = wq
        end if
      
        if ( b_s1 == 0 ) then
          Us1 = uq
          Vs1 = vq
          Ws1 = wq
        end if
			
        if ( b_t1 == 0 ) then
          Ut1 = uq
          Vt1 = vq
          Wt1 = wq
        end if
      
        if ( b_b1 == 0 ) then
          Ub1 = uq
          Vb1 = vq
          Wb1 = wq
        end if
        
        if ( ibits(bvx, bc_face_E, 1) == 0 ) then
          Ue1 = u_bc
          Ve1 = v_bc
          We1 = w_bc
        end if
        
        if ( ibits(bvx, bc_face_W, 1) == 0 ) then
          Uw1 = u_bc
          Vw1 = v_bc
          Ww1 = w_bc
        end if
        
        if ( ibits(bvx, bc_face_N, 1) == 0 ) then
          Un1 = u_bc
          Vn1 = v_bc
          Wn1 = w_bc
        end if
        
        if ( ibits(bvx, bc_face_S, 1) == 0 ) then
          Us1 = u_bc
          Vs1 = v_bc
          Ws1 = w_bc
        end if
        
        if ( ibits(bvx, bc_face_T, 1) == 0 ) then
          Ut1 = u_bc
          Vt1 = v_bc
          Wt1 = w_bc
        end if
        
        if ( ibits(bvx, bc_face_B, 1) == 0 ) then
          Ub1 = u_bc
          Vb1 = v_bc
          Wb1 = w_bc
        end if
      
        s1 = (Ue1 + Uw1 + Un1 + Us1 + Ut1 + Ub1) * dh2 + vc(i,j,k,1)
        s2 = (Ve1 + Vw1 + Vn1 + Vs1 + Vt1 + Vb1) * dh2 + vc(i,j,k,2)
        s3 = (We1 + Ww1 + Wn1 + Ws1 + Wt1 + Wb1) * dh2 + vc(i,j,k,3)
        dv1= (dd*s1-Up0)*actv
        dv2= (dd*s2-Vp0)*actv
        dv3= (dd*s3-Wp0)*actv
        wk(i,j,k,1) = Up0 + omg*dv1
        wk(i,j,k,2) = Vp0 + omg*dv2
        wk(i,j,k,3) = Wp0 + omg*dv3

        ddv = dv1*dv1 + dv2*dv2 + dv3*dv3
        dv = max( dv, sqrt(ddv) )
      end if
    end do
    end do
    end do
    
    flop = flop + real(8 + 48*m)

    return
    end subroutine cbc_vis_cn_mod_jcb

!  *******************************************************************************************
!> @subroutine cbc_friction_velocity (ut, sz, g, dh, re, v, bp, range_Yp, range_Ut, v00, flop)
!! @brief 壁関数での摩擦速度の計算
!! @param ut 摩擦速度
!! @param sz 配列長
!! @param g ガイドセル長
!! @param dh 格子幅
!! @param re レイノルズ数
!! @param v 速度ベクトル（nステップ）
!! @param bp BCindex P
!! @param range_Yp Y+の最小値と最大値
!! @param range_Ut 摩擦速度（無次元）の最小値と最大値
!! @param v00 参照速度
!! @param flop
!<
    subroutine cbc_friction_velocity (ut, sz, g, dh, re, v, bp, range_Yp, range_Ut, v00, flop)
    implicit none
    include 'cbc_f_params.h'
    include 'sklparaf.h'
    integer                                                     ::  i, j, k, ix, jx, kx, g, bpx, itr, itrMax, iret, ierr
    integer, dimension(3)                                       ::  sz
    real, dimension(2)                                          ::  range_Yp, range_Ut, tmp_Max, tmp_Min, tmp
    real                                                        ::  dh, dd, dis, eps1, eps2, re, flop
    real                                                        ::  T_w, U_t, U_t0, Yp, dut
    real                                                        ::  u_ref, v_ref, w_ref, u1, u2, u3, ug
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  ut
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)   ::  v
    real, dimension(0:3)                                        ::  v00
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bp

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)

    range_Yp(1) =  1.0e6
    range_Yp(2) = -1.0e6
    range_Ut(1) =  1.0e6
    range_Ut(2) = -1.0e6
    tmp(1) = 0.0
    tmp(2) = 0.0

    dd    = 2.0/(dh*re)
    dis   = 0.5*dh
    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)
    
    eps1 = 1.0e-3 ! 速度がeps1よりも小さい値のときには摩擦速度をゼロにする
    eps2 = 1.0e-3 ! 反復の収束判定閾値
    itrMax = 10   ! Newton反復の最大反復数
    flop = flop + 4.0

    do k=1,kx
    do j=1,jx
    do i=1,ix
      bpx = bp(i,j,k)

      u1 = v(i,j,k,1) - u_ref
      u2 = v(i,j,k,2) - v_ref
      u3 = v(i,j,k,3) - w_ref
      ug = sqrt(u1*u1 + u2*u2 + u3*u3)
      U_t0 = 0.0
      
      flop = flop + 23.0
      
      if ( (0 /= ibits(bpx, facing_W, 6)) .and. (ug >= eps1) ) then ! 6面のうちのどれか隣接壁への方向フラグが立っている場合，かつ速度がeps1以上
        T_w = ug*dd
        U_t0= sqrt(T_w)
        Yp  = dis*U_t0*re
        flop = flop + 18.0
        
        do itr=1,itrMax
          dut = ( ug - (5.75*log10(Yp)+5.5)*U_t0 ) / ( ug/U_t0 + 1.0/0.4 )
          U_t = U_t0 + dut
          if ( U_t <= eps1 ) U_t=eps1
          
          U_t0 = U_t
          Yp   = dis*U_t0*re
          if ( abs(dut)/abs(U_t0) < eps2 ) exit
        end do
        
        range_Yp(1) = min(range_Yp(1), Yp)
        range_Yp(2) = max(range_Yp(2), Yp)
        range_Ut(1) = min(range_Ut(1), U_t0)
        range_Ut(2) = max(range_Ut(2), U_t0)
        
        flop = flop + real(itr)*40.0
      endif
      
      ut(i,j,k) = U_t0
      
    end do
    end do
    end do
    
    call SklIsParallel(iret)
		if ( iret == 1 ) then
			tmp_Min(1) = range_Yp(1)
      tmp_Min(2) = range_Ut(1)
      call SklAllreduce(tmp_Min, tmp, 2, SKL_REAL, SKL_MIN, SKL_DEFAULT_GROUP, ierr)
      range_Yp(1) = tmp(1)
      range_Ut(1) = tmp(2)
      
      tmp_Max(1) = range_Yp(2)
      tmp_Max(2) = range_Ut(2)
      call SklAllreduce(tmp_Max, tmp, 2, SKL_REAL, SKL_MAX, SKL_DEFAULT_GROUP, ierr)
      range_Yp(2) = tmp(1)
      range_Ut(2) = tmp(2)
		end if

    return
    end subroutine cbc_friction_velocity
