!###################################################################################
!
! FFV-C
! Frontflow / violet Cartesian
!
!
! Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
! All rights reserved.
!
! Copyright (c) 2011-2014 Institute of Industrial Science, The University of Tokyo.
! All rights reserved.
!
! Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
! All rights reserved.
!
!###################################################################################

!> @file   ffv_velocity_binary.f90
!! @brief  速度計算のルーチン群（バイナリモデル）
!! @author aics
!<

!> ********************************************************************
!! @brief 対流項と粘性項の計算
!! @param [out] wv        疑似ベクトルの空間項
!! @param [in]  sz        配列長
!! @param [in]  g         ガイドセル長
!! @param [in]  dh        格子幅
!! @param [in]  c_scheme  対流項スキームのモード（1-UWD, 2-center, 3-MUSCL）
!! @param [in]  v00       参照速度
!! @param [in]  rei       レイノルズ数の逆数
!! @param [in]  v         セルセンター速度ベクトル（n-step）
!! @param [in]  vf        セルフェイス速度ベクトル（n-step）
!! @param [in]  bv        BCindex C
!! @param [in]  bp        BCindex P
!! @param [in]  v_mode    粘性項のモード (0=対流項のみ, 1=対流項と粘性項，2=粘性項は壁法則)
!! @param [in]  ut        摩擦速度
!! @param [in]  wall_type 壁面条件 (0=no_slip, 1=slip)
!! @param [in]  bd        BCindex B
!! @param [in]  cvf       コンポーネント体積率
!! @param [in]  vcs_coef  粘性項の係数
!! @param [out] flop      浮動小数点演算数
!<
    subroutine pvec_muscl (wv, sz, g, dh, c_scheme, v00, rei, v, vf, bv, bp, v_mode, ut, wall_type, bd, cvf, vcs_coef, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, g, c_scheme, bvx, v_mode, bpx, wall_type, bdx
    integer, dimension(3)                                     ::  sz
    double precision                                          ::  flop, vflop
    real                                                      ::  b_e1, b_w1, b_n1, b_s1, b_t1, b_b1
    real                                                      ::  b_e2, b_w2, b_n2, b_s2, b_t2, b_b2, b_p
    real                                                      ::  UPe, UPw, VPn, VPs, WPt, WPb, u1, u2, u3, ug, e1, e2, e3, u_tau
    real                                                      ::  Up0, Ue1, Ue2, Uw1, Uw2, Us1, Us2, Un1, Un2, Ub1, Ub2, Ut1, Ut2
    real                                                      ::  Vp0, Ve1, Ve2, Vw1, Vw2, Vs1, Vs2, Vn1, Vn2, Vb1, Vb2, Vt1, Vt2
    real                                                      ::  Wp0, We1, We2, Ww1, Ww2, Ws1, Ws2, Wn1, Wn2, Wb1, Wb2, Wt1, Wt2
    real                                                      ::  ck, dh, dh1, dh2, vcs, vcs_coef, tmp1, tmp2
    real                                                      ::  u_ref, v_ref, w_ref, u_ref2, v_ref2, w_ref2
    real                                                      ::  c_e, c_w, c_n, c_s, c_t, c_b, wls, wm1, wm2, cm1, cm2, ss_4
    real                                                      ::  uu_e, uu_w, uu_s, uu_n, uu_b, uu_t
    real                                                      ::  vv_e, vv_w, vv_s, vv_n, vv_b, vv_t
    real                                                      ::  ww_e, ww_w, ww_s, ww_n, ww_b, ww_t
    real                                                      ::  dv1, dv2, dv3, dv4, g1, g2, g3, g4, g5, g6, s1, s2, s3, s4, b
    real                                                      ::  Urr, Url, Ulr, Ull, Vrr, Vrl, Vlr, Vll, Wrr, Wrl, Wlr, Wll
    real                                                      ::  cr, cl, acr, acl, cnv_u, cnv_v, cnv_w, EX, EY, EZ, rei, beta, qtz
    real                                                      ::  fu_r, fu_l, fv_r, fv_l, fw_r, fw_l, uq, vq, wq, ss
    real                                                      ::  lmt_w, lmt_e, lmt_s, lmt_n, lmt_b, lmt_t
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v, wv, vf
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  ut
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  cvf
    real, dimension(0:3)                                      ::  v00
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bv, bp, bd
    
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

    ! vcs = 1.0 (Euler Explicit) / 0.5 (CN)
    vcs = vcs * vcs_coef
    
    ! 参照座標系の速度
    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)
    u_ref2 = 2.0*u_ref
    v_ref2 = 2.0*v_ref
    w_ref2 = 2.0*w_ref
    
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
    
    ss_4 = 0.25*ss
    
    cm1 = 1.0 - ck
    cm2 = 1.0 + ck
    
    ! 壁面条件
    if ( wall_type == 0 ) then ! no slip
      wls = 0.0
    else if ( wall_type == 1 ) then ! slip
      wls = 1.0
    endif
    
    wm1 = 2.0*(1.0-wls)
    wm2 = 2.0*wls-1.0
    
! /4 + 12 = 8*4+13  = 45 flops
! /4 + 12 = 13*4+13 = 65 flops ! DP
! ループ内セットアップ   6 + 36 + 2 = 44 flops
! 対流項の各方向　( 7+6+7 + 78*3 + 12 ) * 3dir = 798
! 粘性項 36 + 36 + 12 = 84
! 粘性項置換部　67 ! DP 92
! Total : 45 + 44 + 798 + 84 = 971 ! DP 65 + 44 + 798 + 84 = 991

    flop = flop + dble(ix)*dble(jx)*dble(kx)*798.0d0
    ! flop = flop + dble(ix)*dble(jx)*dble(kx)*991.0d0 ! DP
    
    vflop = 0.0d0

!$OMP PARALLEL &
!$OMP REDUCTION(+:vflop) &
!$OMP FIRSTPRIVATE(ix, jx, kx, dh1, dh2, qtz, vcs, b, ck, ss_4, ss, cm1, cm2, wls) &
!$OMP FIRSTPRIVATE(u_ref, v_ref, w_ref, u_ref2, v_ref2, w_ref2, wm1, wm2, rei, v_mode) &
!$OMP PRIVATE(cnv_u, cnv_v, cnv_w, bvx, bpx, uq, vq, wq, tmp1, tmp2, bdx) &
!$OMP PRIVATE(Up0, Ue1, Ue2, Uw1, Uw2, Us1, Us2, Un1, Un2, Ub1, Ub2, Ut1, Ut2) &
!$OMP PRIVATE(Vp0, Ve1, Ve2, Vw1, Vw2, Vs1, Vs2, Vn1, Vn2, Vb1, Vb2, Vt1, Vt2) &
!$OMP PRIVATE(Wp0, We1, We2, Ww1, Ww2, Ws1, Ws2, Wn1, Wn2, Wb1, Wb2, Wt1, Wt2) &
!$OMP PRIVATE(b_e1, b_w1, b_n1, b_s1, b_t1, b_b1, b_e2, b_w2, b_n2, b_s2, b_t2, b_b2, b_p) &
!$OMP PRIVATE(c_e, c_w, c_n, c_s, c_t, c_b) &
!$OMP PRIVATE(UPe, UPw, VPn, VPs, WPt, WPb) &
!$OMP PRIVATE(cr, cl, acr, acl, beta) &
!$OMP PRIVATE(dv1, dv2, dv3, dv4, g1, g2, g3, g4, g5, g6, s1, s2, s3, s4) &
!$OMP PRIVATE(Urr, Url, Ulr, Ull, Vrr, Vrl, Vlr, Vll, Wrr, Wrl, Wlr, Wll) &
!$OMP PRIVATE(fu_r, fu_l, fv_r, fv_l, fw_r, fw_l) &
!$OMP PRIVATE(uu_e, uu_w, uu_s, uu_n, uu_b, uu_t) &
!$OMP PRIVATE(vv_e, vv_w, vv_s, vv_n, vv_b, vv_t) &
!$OMP PRIVATE(ww_e, ww_w, ww_s, ww_n, ww_b, ww_t) &
!$OMP PRIVATE(u1, u2, u3, ug, e1, e2, e3, u_tau) &
!$OMP PRIVATE(lmt_w, lmt_e, lmt_s, lmt_n, lmt_b, lmt_t)

!$OMP DO SCHEDULE(static)

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
      
      ! (i,j,k)からみたセル状態 (0-solid / 1-fluid)
      b_p = real(ibits(bvx, State, 1))

      ! 各方向に物体があれば，マスクはゼロ
      ! セル界面のフラグとしても利用 (0-wall face / 1-fluid)
      b_w1 = real(ibits(bpx, bc_n_W, 1))
      b_e1 = real(ibits(bpx, bc_n_E, 1))
      b_s1 = real(ibits(bpx, bc_n_S, 1))
      b_n1 = real(ibits(bpx, bc_n_N, 1))
      b_b1 = real(ibits(bpx, bc_n_B, 1))
      b_t1 = real(ibits(bpx, bc_n_T, 1))

      ! (i,j,k)を基準にした遠い方向なので，隣接セルで判断
      b_w2 = real(ibits(bp(i-1,j  ,k  ), bc_n_W, 1))
      b_e2 = real(ibits(bp(i+1,j  ,k  ), bc_n_E, 1))
      b_s2 = real(ibits(bp(i  ,j-1,k  ), bc_n_S, 1))
      b_n2 = real(ibits(bp(i  ,j+1,k  ), bc_n_N, 1))
      b_b2 = real(ibits(bp(i  ,j  ,k-1), bc_n_B, 1))
      b_t2 = real(ibits(bp(i  ,j  ,k+1), bc_n_T, 1))


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
      
      ! ステンシルの参照先がvspec, outflowである場合のスキームの破綻を回避，１次精度におとす
      lmt_w = 1.0
      lmt_e = 1.0
      lmt_s = 1.0
      lmt_n = 1.0
      lmt_b = 1.0
      lmt_t = 1.0

      if ( ibits(bp(i-1, j  , k  ), vbc_uwd, 1) == 1) lmt_w = 0.0
      if ( ibits(bp(i+1, j  , k  ), vbc_uwd, 1) == 1) lmt_e = 0.0
      if ( ibits(bp(i  , j-1, k  ), vbc_uwd, 1) == 1) lmt_s = 0.0
      if ( ibits(bp(i  , j+1, k  ), vbc_uwd, 1) == 1) lmt_n = 0.0
      if ( ibits(bp(i  , j  , k-1), vbc_uwd, 1) == 1) lmt_b = 0.0
      if ( ibits(bp(i  , j  , k+1), vbc_uwd, 1) == 1) lmt_t = 0.0

!if ( (ibits(bv(i-1, j  , k  ), bc_face_W, bitw_5) /= 0) .and. (ibits(bp(i-1, j  , k  ), vbc_uwd, 1) == 1) ) lmt_w = 0.0
!if ( (ibits(bv(i+1, j  , k  ), bc_face_E, bitw_5) /= 0) .and. (ibits(bp(i+1, j  , k  ), vbc_uwd, 1) == 1) ) lmt_e = 0.0
!if ( (ibits(bv(i  , j-1, k  ), bc_face_S, bitw_5) /= 0) .and. (ibits(bp(i  , j-1, k  ), vbc_uwd, 1) == 1) ) lmt_s = 0.0
!if ( (ibits(bv(i  , j+1, k  ), bc_face_N, bitw_5) /= 0) .and. (ibits(bp(i  , j+1, k  ), vbc_uwd, 1) == 1) ) lmt_n = 0.0
!if ( (ibits(bv(i  , j  , k-1), bc_face_B, bitw_5) /= 0) .and. (ibits(bp(i  , j  , k-1), vbc_uwd, 1) == 1) ) lmt_b = 0.0
!if ( (ibits(bv(i  , j  , k+1), bc_face_T, bitw_5) /= 0) .and. (ibits(bp(i  , j  , k+1), vbc_uwd, 1) == 1) ) lmt_t = 0.0


      ! 界面速度（スタガード位置） > 24 flops
      UPe = vf(i  , j  , k  ,1)*b_e1 + u_ref*(1.0-b_e1)
      UPw = vf(i-1, j  , k  ,1)*b_w1 + u_ref*(1.0-b_w1)
      VPn = vf(i  , j  , k  ,2)*b_n1 + v_ref*(1.0-b_n1)
      VPs = vf(i  , j-1, k  ,2)*b_s1 + v_ref*(1.0-b_s1)
      WPt = vf(i  , j  , k  ,3)*b_t1 + w_ref*(1.0-b_t1)
      WPb = vf(i  , j  , k-1,3)*b_b1 + w_ref*(1.0-b_b1)


      ! セルセンターからの壁面修正速度 > 2 flops
      uq = u_ref2 - Up0
      vq = v_ref2 - Vp0
      wq = w_ref2 - Wp0
			
      ! X方向 ---------------------------------------
      
      ! 速度指定の場合にMUSCLスキームの参照先として，固体内にテンポラリに与えた値を使う
      if ( (b_e2 == 0.0)  ) then  ! 7 flops
        Ue2 = u_ref2    - v(i+1,j  ,k  , 1)
        Ve2 = wm1*v_ref + v(i+1,j  ,k  , 2)*wm2
        We2 = wm1*w_ref + v(i+1,j  ,k  , 3)*wm2
      endif
      
      ! 壁面の場合の参照速度の修正
      tmp1 = wm1*v_ref + Vp0*wm2 ! 6 flops
      tmp2 = wm1*w_ref + Wp0*wm2
      
      if ( b_e1 == 0.0 ) then
        Ue1 = uq
        Ve1 = tmp1
        We1 = tmp2
      endif
      
      if ( b_w1 == 0.0 ) then
        Uw1 = uq
        Vw1 = tmp1
        Ww1 = tmp2
      end if
      
      if ( (b_w2 == 0.0)  ) then ! 7 flops
        Uw2 = u_ref2    - v(i-1,j  ,k  , 1)
        Vw2 = wm1*v_ref + v(i-1,j  ,k  , 2)*wm2
        Ww2 = wm1*w_ref + v(i-1,j  ,k  , 3)*wm2
      end if
      
      ! 流束　流体のみと固体壁の影響を含み，隣接セルが固体の場合にはマスクする
      cr  = UPe - u_ref
      cl  = UPw - u_ref
      acr = abs(cr)
      acl = abs(cl)
      
      dv4 = Ue2-Ue1
      dv3 = Ue1-Up0
      dv2 = Up0-Uw1
      dv1 = Uw1-Uw2
      
      include 'muscl.h'
      
      Urr = Ue1 - (cm1*g6+cm2*g5)*ss_4 * lmt_e
      Url = Up0 + (cm1*g3+cm2*g4)*ss_4
      Ulr = Up0 - (cm1*g4+cm2*g3)*ss_4
      Ull = Uw1 + (cm1*g1+cm2*g2)*ss_4 * lmt_w
      fu_r = 0.5*(cr*(Urr+Url) - acr*(Urr-Url)) * b_e1
      fu_l = 0.5*(cl*(Ulr+Ull) - acl*(Ulr-Ull)) * b_w1 ! > 4 + 4 + 36 + 5*4+7*2 = 78 flops

      dv4 = Ve2-Ve1
      dv3 = Ve1-Vp0
      dv2 = Vp0-Vw1
      dv1 = Vw1-Vw2
      
      include 'muscl.h'
      
      Vrr = Ve1 - (cm1*g6+cm2*g5)*ss_4 * lmt_e
      Vrl = Vp0 + (cm1*g3+cm2*g4)*ss_4
      Vlr = Vp0 - (cm1*g4+cm2*g3)*ss_4
      Vll = Vw1 + (cm1*g1+cm2*g2)*ss_4 * lmt_w
      fv_r = 0.5*(cr*(Vrr+Vrl) - acr*(Vrr-Vrl)) * b_e1
      fv_l = 0.5*(cl*(Vlr+Vll) - acl*(Vlr-Vll)) * b_w1

      dv4 = We2-We1
      dv3 = We1-Wp0
      dv2 = Wp0-Ww1
      dv1 = Ww1-Ww2
      
      include 'muscl.h'
      
      Wrr = We1 - (cm1*g6+cm2*g5)*ss_4 * lmt_e
      Wrl = Wp0 + (cm1*g3+cm2*g4)*ss_4
      Wlr = Wp0 - (cm1*g4+cm2*g3)*ss_4
      Wll = Ww1 + (cm1*g1+cm2*g2)*ss_4 * lmt_w
      fw_r = 0.5*(cr*(Wrr+Wrl) - acr*(Wrr-Wrl)) * b_e1
      fw_l = 0.5*(cl*(Wlr+Wll) - acl*(Wlr-Wll)) * b_w1
      
      ! 流束の加算　VBCでない面の寄与のみを評価する
      cnv_u = cnv_u + fu_r*c_e - fu_l*c_w
      cnv_v = cnv_v + fv_r*c_e - fv_l*c_w
      cnv_w = cnv_w + fw_r*c_e - fw_l*c_w ! > 4*3 = 12 flops
			
      ! Y方向 ---------------------------------------
      
      if ( (b_n2 == 0.0)  ) then
        Un2 = wm1*u_ref + v(i  ,j+1,k  , 1)*wm2
        Vn2 = v_ref2    - v(i  ,j+1,k  , 2)
        Wn2 = wm1*w_ref + v(i  ,j+1,k  , 3)*wm2
      endif
      
      tmp1 = wm1*u_ref + Up0*wm2
      tmp2 = wm1*w_ref + Wp0*wm2
      
      if ( b_n1 == 0.0 ) then
        Un1 = tmp1
        Vn1 = vq
        Wn1 = tmp2
      endif
      
      if ( b_s1 == 0.0 ) then
        Us1 = tmp1
        Vs1 = vq
        Ws1 = tmp2
      endif
      
      if ( (b_s2 == 0.0)  ) then
        Us2 = wm1*u_ref + v(i  ,j-1,k  , 1)*wm2
        Vs2 = v_ref2    - v(i  ,j-1,k  , 2)
        Ws2 = wm1*w_ref + v(i  ,j-1,k  , 3)*wm2
      endif
      
      ! 流束　流体のみと固体壁の影響を含み，隣接セルが固体の場合にはマスクする
      cr  = VPn - v_ref
      cl  = VPs - v_ref
      acr = abs(cr)
      acl = abs(cl)
      
      dv4 = Un2-Un1
      dv3 = Un1-Up0
      dv2 = Up0-Us1
      dv1 = Us1-Us2
      
      include 'muscl.h'
      
      Urr = Un1 - (cm1*g6+cm2*g5)*ss_4 * lmt_n
      Url = Up0 + (cm1*g3+cm2*g4)*ss_4
      Ulr = Up0 - (cm1*g4+cm2*g3)*ss_4
      Ull = Us1 + (cm1*g1+cm2*g2)*ss_4 * lmt_s
      fu_r = 0.5*(cr*(Urr+Url) - acr*(Urr-Url)) * b_n1
      fu_l = 0.5*(cl*(Ulr+Ull) - acl*(Ulr-Ull)) * b_s1

      dv4 = Vn2-Vn1
      dv3 = Vn1-Vp0
      dv2 = Vp0-Vs1
      dv1 = Vs1-Vs2
      
      include 'muscl.h'
      
      Vrr = Vn1 - (cm1*g6+cm2*g5)*ss_4 * lmt_n
      Vrl = Vp0 + (cm1*g3+cm2*g4)*ss_4
      Vlr = Vp0 - (cm1*g4+cm2*g3)*ss_4
      Vll = Vs1 + (cm1*g1+cm2*g2)*ss_4 * lmt_s
      fv_r = 0.5*(cr*(Vrr+Vrl) - acr*(Vrr-Vrl)) * b_n1
      fv_l = 0.5*(cl*(Vlr+Vll) - acl*(Vlr-Vll)) * b_s1

      dv4 = Wn2-Wn1
      dv3 = Wn1-Wp0
      dv2 = Wp0-Ws1
      dv1 = Ws1-Ws2
      
      include 'muscl.h'
      
      Wrr = Wn1 - (cm1*g6+cm2*g5)*ss_4 * lmt_n
      Wrl = Wp0 + (cm1*g3+cm2*g4)*ss_4
      Wlr = Wp0 - (cm1*g4+cm2*g3)*ss_4
      Wll = Ws1 + (cm1*g1+cm2*g2)*ss_4 * lmt_s
      fw_r = 0.5*(cr*(Wrr+Wrl) - acr*(Wrr-Wrl)) * b_n1
      fw_l = 0.5*(cl*(Wlr+Wll) - acl*(Wlr-Wll)) * b_s1
      
      ! 流束の加算　VBCでない面の寄与のみを評価する
      cnv_u = cnv_u + fu_r*c_n - fu_l*c_s
      cnv_v = cnv_v + fv_r*c_n - fv_l*c_s
      cnv_w = cnv_w + fw_r*c_n - fw_l*c_s
			
      ! Z方向 ---------------------------------------
      
      ! 壁面の場合の参照速度の修正
      if ( (b_t2 == 0.0)  ) then
        Ut2 = wm1*u_ref + v(i  ,j  ,k+1, 1)*wm2
        Vt2 = wm1*v_ref + v(i  ,j  ,k+1, 2)*wm2
        Wt2 = w_ref2    - v(i  ,j  ,k+1, 3)
      end if
      
      tmp1 = wm1*u_ref + Up0*wm2
      tmp2 = wm1*v_ref + Vp0*wm2
      
      if ( b_t1 == 0.0 ) then
        Ut1 = tmp1
        Vt1 = tmp2
        Wt1 = wq
      end if
      
      if ( b_b1 == 0.0 ) then
        Ub1 = tmp1
        Vb1 = tmp2
        Wb1 = wq
      end if
      
      if ( (b_b2 == 0.0)  ) then
        Ub2 = wm1*u_ref + v(i  ,j  ,k-1, 1)*wm2
        Vb2 = wm1*v_ref + v(i  ,j  ,k-1, 2)*wm2
        Wb2 = w_ref2    - v(i  ,j  ,k-1, 3)
      end if
      
      ! 流束　流体のみと固体壁の影響を含み，隣接セルが固体の場合にはマスクする
      cr  = WPt - w_ref
      cl  = WPb - w_ref
      acr = abs(cr)
      acl = abs(cl)

      dv4 = Ut2-Ut1
      dv3 = Ut1-Up0
      dv2 = Up0-Ub1
      dv1 = Ub1-Ub2
      
      include 'muscl.h'
      
      Urr = Ut1 - (cm1*g6+cm2*g5)*ss_4 * lmt_t
      Url = Up0 + (cm1*g3+cm2*g4)*ss_4
      Ulr = Up0 - (cm1*g4+cm2*g3)*ss_4
      Ull = Ub1 + (cm1*g1+cm2*g2)*ss_4 * lmt_b
      fu_r = 0.5*(cr*(Urr+Url) - acr*(Urr-Url)) * b_t1
      fu_l = 0.5*(cl*(Ulr+Ull) - acl*(Ulr-Ull)) * b_b1

      dv4 = Vt2-Vt1
      dv3 = Vt1-Vp0
      dv2 = Vp0-Vb1
      dv1 = Vb1-Vb2
      
      include 'muscl.h'
            
      Vrr = Vt1 - (cm1*g6+cm2*g5)*ss_4 * lmt_t
      Vrl = Vp0 + (cm1*g3+cm2*g4)*ss_4
      Vlr = Vp0 - (cm1*g4+cm2*g3)*ss_4
      Vll = Vb1 + (cm1*g1+cm2*g2)*ss_4 * lmt_b
      fv_r = 0.5*(cr*(Vrr+Vrl) - acr*(Vrr-Vrl)) * b_t1
      fv_l = 0.5*(cl*(Vlr+Vll) - acl*(Vlr-Vll)) * b_b1

      dv4 = Wt2-Wt1
      dv3 = Wt1-Wp0
      dv2 = Wp0-Wb1
      dv1 = Wb1-Wb2
      
      include 'muscl.h'
            
      Wrr = Wt1 - (cm1*g6+cm2*g5)*ss_4 * lmt_t
      Wrl = Wp0 + (cm1*g3+cm2*g4)*ss_4
      Wlr = Wp0 - (cm1*g4+cm2*g3)*ss_4
      Wll = Wb1 + (cm1*g1+cm2*g2)*ss_4 * lmt_b
      fw_r = 0.5*(cr*(Wrr+Wrl) - acr*(Wrr-Wrl)) * b_t1
      fw_l = 0.5*(cl*(Wlr+Wll) - acl*(Wlr-Wll)) * b_b1
      
      ! 流束の加算　VBCでない面の寄与のみを評価する
      cnv_u = cnv_u + fu_r*c_t - fu_l*c_b
      cnv_v = cnv_v + fv_r*c_t - fv_l*c_b
      cnv_w = cnv_w + fw_r*c_t - fw_l*c_b

      
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
      ww_b = ( Wp0 - Wb1 ) * dh2 ! > (2*6)*3 = 36 flops
      
      if ( v_mode == 2 ) then ! 壁法則の場合の壁面摩擦による剪断応力の置換
        if ( ibits(bpx, facing_W, 6) /= 0 ) then ! 6面のうちのどれか方向フラグが立っている場合，つまり壁面隣接セル
          u1 = Up0 - u_ref
          u2 = Vp0 - v_ref
          u3 = Wp0 - w_ref
          ug = sqrt(u1*u1 + u2*u2 + u3*u3)
          
          if ( ug > 0.0 ) then
            e1 = abs(u1/ug)
            e2 = abs(u2/ug)
            e3 = abs(u3/ug)
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
          endif
          
          if ( ibits(bpx, facing_W, 1) == 1 ) then
            vv_w = u_tau * e2*e2
            ww_w = u_tau * e3*e3
          endif
          
          if ( ibits(bpx, facing_N, 1) == 1 ) then
            uu_n = u_tau * e1*e1
            ww_n = u_tau * e3*e3
          endif
          
          if ( ibits(bpx, facing_S, 1) == 1 ) then
            uu_s = u_tau * e1*e1
            ww_s = u_tau * e3*e3
          endif
          
          if ( ibits(bpx, facing_T, 1) == 1 ) then
            uu_t = u_tau * e1*e1
            vv_t = u_tau * e2*e2
          endif
          
          if ( ibits(bpx, facing_B, 1) == 1 ) then
            uu_b = u_tau * e1*e1
            vv_b = u_tau * e2*e2
          endif
          
          vflop = vflop + 67.0d0
          ! > 9 + sqrt*1 + /3 + 4*6 = 9+10+8*3+24 = 67 ! DP 9+20+13*3+24 = 92
        endif
      endif 
      
      beta = 1.0
      bdx = bd(i,j,k)
      if (ibits(bdx, forcing_bit, 1) == 1) then ! 圧力損失コンポの場合
        ! beta = 1.0 - real(ibits( bdx, top_vf, bitw_8 )) * qtz ! 1-体積率
        beta = 1.0 - cvf(i,j,k)
      endif
      
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
           - ww_b * c_b ) * dh1 ! > 2*6*3 = 36 flops
			
      ! 対流項と粘性項の和 > 4*3 = 12 flops
      wv(i,j,k,1) = -cnv_u*dh1 + beta*EX*vcs
      wv(i,j,k,2) = -cnv_v*dh1 + beta*EY*vcs
      wv(i,j,k,3) = -cnv_w*dh1 + beta*EZ*vcs
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    flop = flop + vflop

    return
    end subroutine pvec_muscl

!> ********************************************************************
!! @brief 次ステップのセルセンターの速度を更新
!! @param [out] v    n+1時刻の速度ベクトル
!! @param [out] div  \sum {u^{n+1}}
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル長
!! @param [in]  dt   時間積分幅
!! @param [in]  dh   格子幅
!! @param [in]  vc   セルセンター疑似速度ベクトル
!! @param [in]  vf   セルフェイス速度ベクトル
!! @param [in]  p    圧力
!! @param [in]  bp   BCindex P
!! @param [in]  bv   BCindex C
!! @param [out] flop 浮動小数点演算数
!! @note 
!!    - actvのマスクはSPEC_VEL/OUTFLOWの参照セルをマスクしないようにbvを使う
!<
    subroutine update_vec (v, div, sz, g, dt, dh, vc, vf, p, bp, bv, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, g, bpx, bvx
    integer, dimension(3)                                     ::  sz
    double precision                                          ::  flop
    real                                                      ::  dh, dt, dd, actv
    real                                                      ::  pc, px, py, pz, pxw, pxe, pys, pyn, pzb, pzt
    real                                                      ::  Ue0, Uw0, Vn0, Vs0, Wt0, Wb0, Up0, Vp0, Wp0
    real                                                      ::  Ue, Uw, Vn, Vs, Wt, Wb
    real                                                      ::  Uef, Uwf, Vnf, Vsf, Wtf, Wbf
    real                                                      ::  c1, c2, c3, c4, c5, c6
    real                                                      ::  N_e, N_w, N_n, N_s, N_t, N_b
    real                                                      ::  D_e, D_w, D_n, D_s, D_t, D_b
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v, vc, vf
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  div, p
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp, bv
    
    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    dd = dt/dh

    flop = flop + dble(ix)*dble(jx)*dble(kx)*106.0 + 8.0d0


!$OMP PARALLEL &
!$OMP PRIVATE(bpx, bvx, actv) &
!$OMP PRIVATE(c1, c2, c3, c4, c5, c6) &
!$OMP PRIVATE(N_e, N_w, N_n, N_s, N_t, N_b) &
!$OMP PRIVATE(D_e, D_w, D_n, D_s, D_t, D_b) &
!$OMP PRIVATE(Ue0, Uw0, Vn0, Vs0, Wt0, Wb0, Up0, Vp0, Wp0) &
!$OMP PRIVATE(Ue, Uw, Vn, Vs, Wt, Wb) &
!$OMP PRIVATE(Uef, Uwf, Vnf, Vsf, Wtf, Wbf) &
!$OMP PRIVATE(pc, px, py, pz, pxw, pxe, pys, pyn, pzb, pzt) &
!$OMP FIRSTPRIVATE(ix, jx, kx, dd)

!$OMP DO SCHEDULE(static)
    do k=1,kx
    do j=1,jx
    do i=1,ix
      bpx = bp(i,j,k)
      bvx = bv(i,j,k)
      actv = real(ibits(bvx, State,  1))

      c1 = 1.0
      c2 = 1.0
      c3 = 1.0
      c4 = 1.0
      c5 = 1.0
      c6 = 1.0
      
      ! Neumann条件のとき，0.0 > 6 flop
      ! 物体があればノイマン条件なので，セルフェイスマスクとしても利用
      N_w = real(ibits(bpx, bc_n_W, 1))  ! w
      N_e = real(ibits(bpx, bc_n_E, 1))  ! e
      N_s = real(ibits(bpx, bc_n_S, 1))  ! s
      N_n = real(ibits(bpx, bc_n_N, 1))  ! n
      N_b = real(ibits(bpx, bc_n_B, 1))  ! b
      N_t = real(ibits(bpx, bc_n_T, 1))  ! t

      ! Dirichlet(p=0)条件のとき 0.0 > 6 flop
      D_w = real(ibits(bpx, bc_d_W, 1))  ! w
      D_e = real(ibits(bpx, bc_d_E, 1))  ! e
      D_s = real(ibits(bpx, bc_d_S, 1))  ! s
      D_n = real(ibits(bpx, bc_d_N, 1))  ! n
      D_b = real(ibits(bpx, bc_d_B, 1))  ! b
      D_t = real(ibits(bpx, bc_d_T, 1))  ! t


      ! 疑似ベクトル
      Uw0 = vc(i-1,j  ,k  , 1)
      Up0 = vc(i  ,j  ,k  , 1)
      Ue0 = vc(i+1,j  ,k  , 1)

      Vs0 = vc(i  ,j-1,k  , 2)
      Vp0 = vc(i  ,j  ,k  , 2)
      Vn0 = vc(i  ,j+1,k  , 2)

      Wb0 = vc(i  ,j  ,k-1, 3)
      Wp0 = vc(i  ,j  ,k  , 3)
      Wt0 = vc(i  ,j  ,k+1, 3)

      Uw = 0.5*( Up0 + Uw0 )*N_w ! 18 flop
      Ue = 0.5*( Up0 + Ue0 )*N_e
      Vs = 0.5*( Vp0 + Vs0 )*N_s
      Vn = 0.5*( Vp0 + Vn0 )*N_n
      Wb = 0.5*( Wp0 + Wb0 )*N_b
      Wt = 0.5*( Wp0 + Wt0 )*N_t
      
      ! c=0.0(VBC), 1.0(Fluid); VBCは内部と外部の両方
      if ( ibits(bvx, bc_face_W, bitw_5) /= 0 ) c1 = 0.0
      if ( ibits(bvx, bc_face_E, bitw_5) /= 0 ) c2 = 0.0
      if ( ibits(bvx, bc_face_S, bitw_5) /= 0 ) c3 = 0.0
      if ( ibits(bvx, bc_face_N, bitw_5) /= 0 ) c4 = 0.0
      if ( ibits(bvx, bc_face_B, bitw_5) /= 0 ) c5 = 0.0
      if ( ibits(bvx, bc_face_T, bitw_5) /= 0 ) c6 = 0.0
      
      ! 圧力勾配 18flop
      pc  = p(i,  j,  k  )
      pxw = (-p(i-1,j  ,k  )*D_w + pc) * N_w
      pxe = ( p(i+1,j  ,k  )*D_e - pc) * N_e
      pys = (-p(i  ,j-1,k  )*D_s + pc) * N_s
      pyn = ( p(i  ,j+1,k  )*D_n - pc) * N_n
      pzb = (-p(i  ,j  ,k-1)*D_b + pc) * N_b
      pzt = ( p(i  ,j  ,k+1)*D_t - pc) * N_t
      px = 0.5*(pxe + pxw)
      py = 0.5*(pyn + pys)
      pz = 0.5*(pzt + pzb)

      ! セルフェイス VBCの寄与と壁面の影響は除外 18flop
      Uwf = (Uw - dd * pxw) * c1
      Uef = (Ue - dd * pxe) * c2
      Vsf = (Vs - dd * pys) * c3
      Vnf = (Vn - dd * pyn) * c4
      Wbf = (Wb - dd * pzb) * c5
      Wtf = (Wt - dd * pzt) * c6

      ! i=1...ix >> vfは0...ixの範囲をカバーするので，通信不要 6flop
      vf(i-1,j  ,k  ,1) = Uwf * N_w
      vf(i  ,j  ,k  ,1) = Uef * N_e
      vf(i  ,j-1,k  ,2) = Vsf * N_s
      vf(i  ,j  ,k  ,2) = Vnf * N_n
      vf(i  ,j  ,k-1,3) = Wbf * N_b
      vf(i  ,j  ,k  ,3) = Wtf * N_t

      div(i,j,k) = (Uef - Uwf + Vnf - Vsf + Wtf - Wbf) * actv ! 6flop
      
      ! セルセンタの速度更新 9flop
      v(i,j,k,1) = ( Up0-px*dd ) * actv
      v(i,j,k,2) = ( Vp0-py*dd ) * actv
      v(i,j,k,3) = ( Wp0-pz*dd ) * actv

    end do
    end do
    end do
!$OMP END DO

!$OMP END PARALLEL

    return
    end subroutine update_vec

!> ********************************************************************
!! @brief 速度の発散に使う\sum{u^*}を計算する
!! @param [out]    div  速度の発散値
!! @param [in]     sz   配列長
!! @param [in]     g    ガイドセル長
!! @param [in]     v0   セルセンター疑似ベクトル
!! @param [in]     bv   BCindex C
!! @param [in]     bp   BCindex P
!! @param [in,out] flop 浮動小数点演算数
!<
    subroutine divergence (div, sz, g, v0, bv, bp, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, g, bvx, bpx
    integer, dimension(3)                                     ::  sz
    double precision                                          ::  flop
    real                                                      ::  Ue, Uw, Vn, Vs, Wt, Wb, actv
    real                                                      ::  Ue0, Uw0, Up0, Vn0, Vs0, Vp0, Wb0, Wt0, Wp0
    real                                                      ::  c_e, c_w, c_n, c_s, c_t, c_b
    real                                                      ::  b_w, b_e, b_s, b_n, b_b, b_t
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v0
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  div
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bv, bp

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    
    ! loop : 1 + 42 + 13
    flop  = flop + dble(ix)*dble(jx)*dble(kx)*37.0d0

!$OMP PARALLEL &
!$OMP PRIVATE(bvx, actv, bpx) &
!$OMP PRIVATE(b_w, b_e, b_s, b_n, b_b, b_t) &
!$OMP PRIVATE(Ue, Uw, Vn, Vs, Wt, Wb) &
!$OMP PRIVATE(Ue0, Uw0, Up0, Vn0, Vs0, Vp0, Wb0, Wt0, Wp0) &
!$OMP PRIVATE(c_e, c_w, c_n, c_s, c_t, c_b) &
!$OMP FIRSTPRIVATE(ix, jx, kx)

!$OMP DO SCHEDULE(static)

    do k=1,kx
    do j=1,jx
    do i=1,ix
      bvx = bv(i,j,k)
      bpx = bp(i,j,k)
      actv= real(ibits(bvx, State, 1))
      
      ! 各セルセンター位置の変数ロード
      Uw0 = v0(i-1,j  ,k  , 1)
      Up0 = v0(i  ,j  ,k  , 1)
      Ue0 = v0(i+1,j  ,k  , 1)

      Vs0 = v0(i  ,j-1,k  , 2)
      Vp0 = v0(i  ,j  ,k  , 2)
      Vn0 = v0(i  ,j+1,k  , 2)

      Wb0 = v0(i  ,j  ,k-1, 3)
      Wp0 = v0(i  ,j  ,k  , 3)
      Wt0 = v0(i  ,j  ,k+1, 3)

      ! 物体があればノイマン条件なので，セルフェイスマスクとしても利用 (0-solid / 1-fluid)
      b_w = real(ibits(bpx, bc_n_W, 1))
      b_e = real(ibits(bpx, bc_n_E, 1))
      b_s = real(ibits(bpx, bc_n_S, 1))
      b_n = real(ibits(bpx, bc_n_N, 1))
      b_b = real(ibits(bpx, bc_n_B, 1))
      b_t = real(ibits(bpx, bc_n_T, 1))

      Uw = 0.5*( Up0 + Uw0 )*b_w
      Ue = 0.5*( Up0 + Ue0 )*b_e
      Vs = 0.5*( Vp0 + Vs0 )*b_s
      Vn = 0.5*( Vp0 + Vn0 )*b_n
      Wb = 0.5*( Wp0 + Wb0 )*b_b
      Wt = 0.5*( Wp0 + Wt0 )*b_t

      ! 各面のVBCフラグ ibits() = 0(Normal) / others(BC) >> c_e = 1.0(Normal) / 0.0(BC)
      c_w = 1.0
      c_e = 1.0
      c_s = 1.0
      c_n = 1.0
      c_b = 1.0
      c_t = 1.0

      if ( ibits(bvx, bc_face_W, bitw_5) /= 0 ) c_w = 0.0
      if ( ibits(bvx, bc_face_E, bitw_5) /= 0 ) c_e = 0.0
      if ( ibits(bvx, bc_face_S, bitw_5) /= 0 ) c_s = 0.0
      if ( ibits(bvx, bc_face_N, bitw_5) /= 0 ) c_n = 0.0
      if ( ibits(bvx, bc_face_B, bitw_5) /= 0 ) c_b = 0.0
      if ( ibits(bvx, bc_face_T, bitw_5) /= 0 ) c_t = 0.0
      
      ! VBC面の影響をフラグで無効化 >> OBC_SPEC_VEL, OBC_WALL  13flops
      div(i,j,k) = ( Ue*c_e - Uw*c_w + Vn*c_n - Vs*c_s + Wt*c_t - Wb*c_b ) * actv
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    return
    end subroutine divergence
    
!> ********************************************************************
!! @brief 疑似ベクトルの時間積分（Euler陽解法）
!! @param [in,out] vc   疑似ベクトル
!! @param [in]     sz   配列長
!! @param [in]     g    ガイドセル長
!! @param [in]     dt   時間積分幅
!! @param [in]     v    速度ベクトル（n-step, collocated）
!! @param [in]     bd   BCindex B
!! @param [in,out] flop 浮動小数点演算数
!! @note ここのマスクはIDのこと，VSPEC, OUTFLOWの増分をキャンセルするため
!<
    subroutine euler_explicit (vc, sz, g, dt, v, bd, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, g
    integer, dimension(3)                                     ::  sz
    double precision                                          ::  flop
    real                                                      ::  actv, dt
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  vc, v
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bd
    
    ix = sz(1)
    jx = sz(2)
    kx = sz(3)

    flop = flop + dble(ix)*dble(jx)*dble(kx)*8.0d0

!$OMP PARALLEL &
!$OMP PRIVATE(actv) &
!$OMP FIRSTPRIVATE(ix, jx, kx, dt)

!$OMP DO SCHEDULE(static)
    do k=1,kx
    do j=1,jx
    do i=1,ix
      actv = dt * real(ibits(bd(i,j,k), State, 1))

      vc(i,j,k,1) = v(i,j,k,1) + vc(i,j,k,1)* actv
      vc(i,j,k,2) = v(i,j,k,2) + vc(i,j,k,2)* actv
      vc(i,j,k,3) = v(i,j,k,3) + vc(i,j,k,3)* actv
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL
    
    return
    end subroutine euler_explicit


!> ********************************************************************
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
    subroutine vis_ee (vc, sz, g, dh, dt, v00, rei, wv, v, bx, cf, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, g, bvx
    integer                                                   ::  b_e1, b_w1, b_n1, b_s1, b_t1, b_b1
    integer, dimension(3)                                     ::  sz
    double precision                                          ::  flop
    real                                                      ::  Up0, Ue1, Uw1, Us1, Un1, Ub1, Ut1
    real                                                      ::  Vp0, Ve1, Vw1, Vs1, Vn1, Vb1, Vt1
    real                                                      ::  Wp0, We1, Ww1, Ws1, Wn1, Wb1, Wt1
    real                                                      ::  u_ref, v_ref, w_ref, dh, dh1, dh2
    real                                                      ::  EX, EY, EZ, rei, cf, dt, uq, vq, wq, actv
    real                                                      ::  c_e, c_w, c_n, c_s, c_t, c_b
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v, wv, vc
    real, dimension(0:3)                                      ::  v00
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bx
    
    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    dh1= 1.0/dh
    dh2= dt*rei*dh1*dh1 * cf
    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)
    flop = flop + dble(ix*jx*kx)*48.0d0 + 5.0d0
    
    do k=1,kx
    do j=1,jx
    do i=1,ix
      Up0 = v(i  ,j  ,k  , 1)
      Uw1 = v(i-1,j  ,k  , 1)
      Ue1 = v(i+1,j  ,k  , 1)
      Us1 = v(i  ,j-1,k  , 1)
      Un1 = v(i  ,j+1,k  , 1)
      Ub1 = v(i  ,j  ,k-1, 1)
      Ut1 = v(i  ,j  ,k+1, 1)
			
      Vp0 = v(i  ,j  ,k  , 2)
      Vw1 = v(i-1,j  ,k  , 2)
      Ve1 = v(i+1,j  ,k  , 2)
      Vs1 = v(i  ,j-1,k  , 2)
      Vn1 = v(i  ,j+1,k  , 2)
      Vb1 = v(i  ,j  ,k-1, 2)
      Vt1 = v(i  ,j  ,k+1, 2)

      Wp0 = v(i  ,j  ,k  , 3)
      Ww1 = v(i-1,j  ,k  , 3)
      We1 = v(i+1,j  ,k  , 3)
      Ws1 = v(i  ,j-1,k  , 3)
      Wn1 = v(i  ,j+1,k  , 3)
      Wb1 = v(i  ,j  ,k-1, 3)
      Wt1 = v(i  ,j  ,k+1, 3)
      
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
    end subroutine vis_ee

!> ********************************************************************
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
!! @note NOCHECK
!<
    subroutine vis_ee_vbc (vc, sz, g, st, ed, dh, dt, v00, rei, v, bx, odr, cf, vec, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, g, bvx, m, odr
    integer, dimension(3)                                     ::  sz, st, ed
    double precision                                          ::  flop
    real                                                      ::  Up0, Ue1, Uw1, Us1, Un1, Ub1, Ut1
    real                                                      ::  Vp0, Ve1, Vw1, Vs1, Vn1, Vb1, Vt1
    real                                                      ::  Wp0, We1, Ww1, Ws1, Wn1, Wb1, Wt1
    real                                                      ::  u_ref, v_ref, w_ref, dh, dh1, dh2
    real                                                      ::  EX, EY, EZ, rei, cf, dt
    real                                                      ::  u_bc_ref, v_bc_ref, w_bc_ref
    real                                                      ::  c_e, c_w, c_n, c_s, c_t, c_b
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  vc, v
    real, dimension(0:3)                                      ::  v00
    real, dimension(3)                                        ::  vec
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bx
    
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
    
    flop = flop + dble(m)*42d0 + 8.0d0

    return
    end subroutine vis_ee_vbc
    
!> ********************************************************************
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
!! @note NOCHECK
!<
    subroutine vis_cn_sor (v, sz, g, dh, dt, v00, rei, omg, vc, bx, cf, dv, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, g, bvx
    integer                                                   ::  b_e1, b_w1, b_n1, b_s1, b_t1, b_b1
    integer, dimension(3)                                     ::  sz
    double precision                                          ::  flop
    real                                                      ::  Up0, Ue1, Uw1, Us1, Un1, Ub1, Ut1
    real                                                      ::  Vp0, Ve1, Vw1, Vs1, Vn1, Vb1, Vt1
    real                                                      ::  Wp0, We1, Ww1, Ws1, Wn1, Wb1, Wt1
    real                                                      ::  u_ref, v_ref, w_ref, dh, dh1, dh2
    real                                                      ::  rei, dd, dv, dv1, dv2, dv3, uq, vq, wq, ddv
    real                                                      ::  omg, s1, s2, s3, cf, dt, actv, b_vbc
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v, vc
    real, dimension(0:3)                                      ::  v00
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bx
    
    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    dh1= 1.0/dh
    dh2= dt*rei*dh1*dh1*cf
    dd = 1.0 / (1.0 + 6.0*dh2)
    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)
    flop = flop + dble(ix*jx*kx)*52.0d0 + 8.0d0
    
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
    end subroutine vis_cn_sor

!> ********************************************************************
!! @brief SOR法による粘性項計算の境界条件による修正
!! @param[in,out] v 疑似ベクトル
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
!! @note NOCHECK
!<
    subroutine vis_cn_mod_sor (v, sz, g, st, ed, dh, dt, v00, rei, omg, vc, bx, cf, dv, vec, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, g, bvx, m
    integer                                                   ::  b_e1, b_w1, b_n1, b_s1, b_t1, b_b1
    integer, dimension(3)                                     ::  sz, st, ed
    double precision                                          ::  flop
    real                                                      ::  Up0, Ue1, Uw1, Us1, Un1, Ub1, Ut1
    real                                                      ::  Vp0, Ve1, Vw1, Vs1, Vn1, Vb1, Vt1
    real                                                      ::  Wp0, We1, Ww1, Ws1, Wn1, Wb1, Wt1
    real                                                      ::  u_ref, v_ref, w_ref, dh, dh1, dh2
    real                                                      ::  rei, dd, dv, dv1, dv2, dv3, ddv
    real                                                      ::  u_bc, v_bc, w_bc, uq, vq, wq
    real                                                      ::  omg, s1, s2, s3, cf, dt, actv
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v, vc
    real, dimension(0:3)                                      ::  v00
    real, dimension(3)                                        ::  vec
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bx
    
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
    
    flop = flop + dble(m)*48.0d0 + 8.0d0

    return
    end subroutine vis_cn_mod_sor

!> ********************************************************************
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
!! @note NOCHECK
!<
    subroutine eddy_viscosity (vt, sz, g, dh, re, cs, v, bx, vt_range, yp_range, v00)
    implicit none
    include 'ffv_f_params.h'
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
      DUDX=(v(i+1,j  ,k  , 1)-v(i-1,j  ,k  , 1)) * dd
      DUDY=(v(i  ,j+1,k  , 1)-v(i  ,j-1,k  , 1)) * dd
      DUDZ=(v(i  ,j  ,k+1, 1)-v(i  ,j  ,k-1, 1)) * dd

      DVDX=(v(i+1,j  ,k  , 2)-v(i-1,j  ,k  , 2)) * dd
      DVDY=(v(i  ,j+1,k  , 2)-v(i  ,j-1,k  , 2)) * dd
      DVDZ=(v(i  ,j  ,k+1, 2)-v(i  ,j  ,k-1, 2)) * dd

      DWDX=(v(i+1,j  ,k  , 3)-v(i-1,j  ,k  , 3)) * dd
      DWDY=(v(i  ,j+1,k  , 3)-v(i  ,j-1,k  , 3)) * dd
      DWDZ=(v(i  ,j  ,k+1, 3)-v(i  ,j  ,k-1, 3)) * dd

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
    end subroutine eddy_viscosity

!> ********************************************************************
!! @brief Adams-Bashforth法による疑似ベクトルの時間積分
!! @param[out] vc 疑似ベクトル
!! @param sz 配列長
!! @param g ガイドセル長
!! @param dt 時間積分幅
!! @param v 速度ベクトル（n-step, collocated）
!! @param ab 前ステップの対流項（＋粘性項）の計算値
!! @param bd BCindex B
!! @param v00 参照速度
!! @param[out] flop
!! @note NOCHECK
!<
    subroutine ab2 (vc, sz, g, dt, v, ab, bd, v00, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, g
    integer, dimension(3)                                     ::  sz
    double precision                                          ::  flop
    real                                                      ::  actv, dt, ab_u, ab_v, ab_w, u_ref, v_ref, w_ref
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  vc, v, ab
    real, dimension(0:3)                                      ::  v00
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bd
    
    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)
    flop = flop + dble(ix*jx*kx) * 27.0d0
    
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
    end subroutine ab2
    
!> ********************************************************************
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
!! @note NOCHECK
!<
    subroutine vis_cn_jcb (v, sz, g, dh, dt, v00, rei, omg, vc, bx, wk, cf, dv, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, g, bvx
    integer                                                   ::  b_e1, b_w1, b_n1, b_s1, b_t1, b_b1
    integer, dimension(3)                                     ::  sz
    double precision                                          ::  flop
    real                                                      ::  Up0, Ue1, Uw1, Us1, Un1, Ub1, Ut1
    real                                                      ::  Vp0, Ve1, Vw1, Vs1, Vn1, Vb1, Vt1
    real                                                      ::  Wp0, We1, Ww1, Ws1, Wn1, Wb1, Wt1
    real                                                      ::  u_ref, v_ref, w_ref, dh, dh1, dh2
    real                                                      ::  rei, dd, dv, dv1, dv2, dv3, uq, vq, wq, ddv
    real                                                      ::  omg, s1, s2, s3, cf, dt, actv, b_vbc
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v, vc, wk
    real, dimension(0:3)                                      ::  v00
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bx
    
    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    dh1= 1.0/dh
    dh2= dt*rei*dh1*dh1*cf
    dd = 1.0 / (1.0 + 6.0*dh2)
    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)
    flop = flop + dble(ix*jx*kx)*52.0d0 + 8.0d0
    
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
    end subroutine vis_cn_jcb

!> ********************************************************************
!! @brief 緩和Jacobi法による粘性項の計算
!! @param[in,out] v 疑似ベクトル
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
!! @note NOCHECK
!<
    subroutine vis_cn_mod_jcb (v, sz, g, st, ed, dh, dt, v00, rei, omg, vc, bx, wk, cf, dv, vec, flop)
    implicit none
    include 'ffv_f_params.h'
    integer                                                   ::  i, j, k, g, bvx, m
    integer                                                   ::  b_e1, b_w1, b_n1, b_s1, b_t1, b_b1
    integer, dimension(3)                                     ::  sz, st, ed
    double precision                                          ::  flop
    real                                                      ::  Up0, Ue1, Uw1, Us1, Un1, Ub1, Ut1
    real                                                      ::  Vp0, Ve1, Vw1, Vs1, Vn1, Vb1, Vt1
    real                                                      ::  Wp0, We1, Ww1, Ws1, Wn1, Wb1, Wt1
    real                                                      ::  u_ref, v_ref, w_ref, dh, dh1, dh2
    real                                                      ::  rei, dd, dv, dv1, dv2, dv3, ddv
    real                                                      ::  u_bc, v_bc, w_bc, uq, vq, wq
    real                                                      ::  omg, s1, s2, s3, cf, dt, actv
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v, vc, wk
    real, dimension(0:3)                                      ::  v00
    real, dimension(3)                                        ::  vec
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bx
    
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
    
    flop = flop + dble(m)*48.0d0 + 8.0d0

    return
    end subroutine vis_cn_mod_jcb

!> ********************************************************************
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
!! @note NOCHECK
!<
    subroutine friction_velocity (ut, sz, g, dh, re, v, bp, range_Yp, range_Ut, v00, flop)
    implicit none
    include 'ffv_f_params.h'

    integer                                                   ::  i, j, k, ix, jx, kx, g, bpx, itr, itrMax, iret, ierr
    integer, dimension(3)                                     ::  sz
    double precision                                          ::  flop
    real, dimension(2)                                        ::  range_Yp, range_Ut, tmp_Max, tmp_Min, tmp
    real                                                      ::  dh, dd, dis, eps1, eps2, re
    real                                                      ::  T_w, U_t, U_t0, Yp, dut
    real                                                      ::  u_ref, v_ref, w_ref, u1, u2, u3, ug
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  ut
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v
    real, dimension(0:3)                                      ::  v00
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp

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
    flop = flop + 4.0d0

    do k=1,kx
    do j=1,jx
    do i=1,ix
      bpx = bp(i,j,k)

      u1 = v(i,j,k,1) - u_ref
      u2 = v(i,j,k,2) - v_ref
      u3 = v(i,j,k,3) - w_ref
      ug = sqrt(u1*u1 + u2*u2 + u3*u3)
      U_t0 = 0.0
      
      flop = flop + 23.0d0

      ! 6面のうちのどれか隣接壁への方向フラグが立っている場合，かつ速度がeps1以上
      if ( (0 /= ibits(bpx, facing_W, 6)) .and. (ug >= eps1) ) then

        T_w = ug*dd
        U_t0= sqrt(T_w)
        Yp  = dis*U_t0*re
        flop = flop + 18.0d0
        
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
        
        flop = flop + dble(itr)*40.0d0
      endif
      
      ut(i,j,k) = U_t0
      
    end do
    end do
    end do
    
!    call SklIsParallel(iret)
		if ( iret == 1 ) then
			tmp_Min(1) = range_Yp(1)
      tmp_Min(2) = range_Ut(1)
!      call SklAllreduce(tmp_Min, tmp, 2, SKL_REAL, SKL_MIN, SKL_DEFAULT_GROUP, ierr)
      range_Yp(1) = tmp(1)
      range_Ut(1) = tmp(2)
      
      tmp_Max(1) = range_Yp(2)
      tmp_Max(2) = range_Ut(2)
!      call SklAllreduce(tmp_Max, tmp, 2, SKL_REAL, SKL_MAX, SKL_DEFAULT_GROUP, ierr)
      range_Yp(2) = tmp(1)
      range_Ut(2) = tmp(2)
		end if

    return
    end subroutine friction_velocity
