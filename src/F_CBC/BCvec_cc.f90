!   *********************************************************
!
!   SPHERE - Skeleton for PHysical and Engineering REsearch
!  
!   Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
!
!   *********************************************************

!> @file BCvec.f90
!> @brief Boundary conditons for collocated velocity 
!> @author keno, FSI Team, VCAD, RIKEN

!  ***********************************************************************************
!> @subroutine cbc_pvec_vibc_oflow (wv, sz, g, st, ed, dh, rei, v, bv, odr, vec, flop)
!! @brief 内部速度境界条件による対流項と粘性項の流束の修正
!! @param[out] wv 疑似ベクトルの空間項の評価値
!! @param sz 配列長
!! @param g ガイドセル長
!! @param st ループの開始インデクス
!! @param ed ループの終了インデクス
!! @param dh 格子幅
!! @param rei Reynolds数の逆数
!! @param v 速度ベクトル（u^n, セルセンタ）
!! @param bv BCindex V
!! @param odr 内部境界処理時の速度境界条件のエントリ
!! @param vec 指定する速度ベクトル
!! @param[out] flop
!! @note vecには，流入条件のとき指定速度，流出条件のとき対流流出速度
!! @todo 内部と外部の分離 do loopの内側に条件分岐を入れているので修正
!! @todo 流出境界はローカルの流束となるように変更する（外部境界参照）
!<
    subroutine cbc_pvec_vibc_oflow (wv, sz, g, st, ed, dh, rei, v, bv, odr, vec, flop)
    implicit none
    include '../FB/cbc_f_params.h'
    integer                                                     ::  i, j, k, g, bvx, odr, is, ie, js, je, ks, ke
    integer, dimension(3)                                       ::  sz, st, ed
    real                                                        ::  Up0, Ue1, Uw1, Us1, Un1, Ub1, Ut1
    real                                                        ::  Vp0, Ve1, Vw1, Vs1, Vn1, Vb1, Vt1
    real                                                        ::  Wp0, We1, Ww1, Ws1, Wn1, Wb1, Wt1
    real                                                        ::  c_e, c_w, c_n, c_s, c_t, c_b
    real                                                        ::  dh, dh1, dh2, flop, EX, EY, EZ, rei
    real                                                        ::  cnv_u, cnv_v, cnv_w, cr, cl, m
    real                                                        ::  u_bc, v_bc, w_bc
    real                                                        ::  fu_r, fu_l, fv_r, fv_l, fw_r, fw_l
    real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  v, wv
    real, dimension(3)                                          ::  vec
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bv
    
    dh1= 1.0/dh
    dh2= rei*dh1*dh1
    
    ! u_bcは境界速度
    u_bc = vec(1)
    v_bc = vec(2)
    w_bc = vec(3)
    
    is = st(1)
    ie = ed(1)
    js = st(2)
    je = ed(2)
    ks = st(3)
    ke = ed(3)
    
    flop = flop + 10.0 ! DP 15 flops
    
    m = 0.0

!$OMP PARALLEL REDUCTION(+:m) &
!$OMP FIRSTPRIVATE(is, ie, js, je, ks, ke, u_bc, v_bc, w_bc, odr) &
!$OMP FIRSTPRIVATE(dh1, dh2) &
!$OMP PRIVATE(bvx, cnv_u, cnv_v, cnv_w, cr, cl) &
!$OMP PRIVATE(Up0, Ue1, Uw1, Us1, Un1, Ub1, Ut1) &
!$OMP PRIVATE(Vp0, Ve1, Vw1, Vs1, Vn1, Vb1, Vt1) &
!$OMP PRIVATE(Wp0, We1, Ww1, Ws1, Wn1, Wb1, Wt1) &
!$OMP PRIVATE(c_e, c_w, c_n, c_s, c_t, c_b) &
!$OMP PRIVATE(fu_r, fu_l, fv_r, fv_l, fw_r, fw_l) &
!$OMP PRIVATE(EX, EY, EZ)

#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif

    do k=ks,ke
    do j=js,je
    do i=is,ie
      bvx = bv(i,j,k)

      if ( 0 /= iand(bvx, bc_mask30) ) then ! 6面のうちのどれか速度境界フラグが立っている場合
        cnv_u = 0.0
        cnv_v = 0.0
        cnv_w = 0.0
        
        ! 変数のロード
        Up0 = v(1,i,j,k)
        Vp0 = v(2,i,j,k)
        Wp0 = v(3,i,j,k)
      
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
          Uw1  = Up0
          Vw1  = Vp0
          Ww1  = Wp0
          cl   = u_bc
          if ( cl>0.0 ) cl=0.0
          fu_l = cl*Up0
          fv_l = cl*Vp0
          fw_l = cl*Wp0 ! 3 flops
        end if
        
        if ( c_e == 1.0 ) then
          Ue1  = Up0
          Ve1  = Vp0
          We1  = Wp0
          cr   = u_bc
          if ( cr<0.0 ) cr=0.0
          fu_r = cr*Up0
          fv_r = cr*Vp0
          fw_r = cr*Wp0 ! 3 flops
        end if
        
        cnv_u = cnv_u + fu_r*c_e - fu_l*c_w
        cnv_v = cnv_v + fv_r*c_e - fv_l*c_w
        cnv_w = cnv_w + fw_r*c_e - fw_l*c_w ! 12 flops
			
        ! Y方向 ---------------------------------------
        if ( c_s == 1.0 ) then
          Us1  = Up0
          Vs1  = Vp0
          Ws1  = Wp0
          cl   = v_bc
          if ( cl>0.0 ) cl=0.0
          fu_l = cl*Up0
          fv_l = cl*Vp0
          fw_l = cl*Wp0
        end if
        
        if ( c_n == 1.0 ) then
          Un1  = Up0
          Vn1  = Vp0
          Wn1  = Wp0
          cr   = v_bc
          if ( cr<0.0 ) cr=0.0
          fu_r = cr*Up0
          fv_r = cr*Vp0
          fw_r = cr*Wp0
        end if
        
        cnv_u = cnv_u + fu_r*c_n - fu_l*c_s
        cnv_v = cnv_v + fv_r*c_n - fv_l*c_s
        cnv_w = cnv_w + fw_r*c_n - fw_l*c_s
			
        ! Z方向 ---------------------------------------
        if ( c_b == 1.0 ) then
          Ub1  = Up0
          Vb1  = Vp0
          Wb1  = Wp0
          cl   = w_bc
          if ( cl>0.0 ) cl=0.0
          fu_l = cl*Up0
          fv_l = cl*Vp0
          fw_l = cl*Wp0
        end if

        if ( c_t == 1.0 ) then
          Ut1  = Up0
          Vt1  = Vp0
          Wt1  = Wp0
          cr   = w_bc
          if ( cr<0.0 ) cr=0.0
          fu_r = cr*Up0
          fv_r = cr*Vp0
          fw_r = cr*Wp0
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
           + ( Ub1 - Up0 ) * c_b ! 17 flops
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

        wv(1,i,j,k) = wv(1,i,j,k) + ( -cnv_u*dh1 + EX*dh2 )
        wv(2,i,j,k) = wv(2,i,j,k) + ( -cnv_v*dh1 + EY*dh2 )
        wv(3,i,j,k) = wv(3,i,j,k) + ( -cnv_w*dh1 + EZ*dh2 ) ! 4*3 = 12 flops
        m = m + 1.0
      endif
    end do
    end do
    end do
    
!$OMP END DO
!$OMP END PARALLEL

    ! loop :  (3+3+12)*3dir + 17*3 + 12 = 117 flops

    flop = flop + m*117.0
    
    return
    end subroutine cbc_pvec_vibc_oflow

!  ****************************************************************************************
!> @subroutine cbc_pvec_vibc_specv (wv, sz, g, st, ed, dh, v00, rei, v, bv, odr, vec, flop)
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
!! @param[out] flop
!! @note vecには，流入条件のとき指定速度，流出条件のとき対流流出速度，カット位置に関わらず指定速度で流束を計算
!! @todo 流出境界はローカルの流束となるように変更する（外部境界参照）
!<
    subroutine cbc_pvec_vibc_specv (wv, sz, g, st, ed, dh, v00, rei, v, bv, odr, vec, flop)
    implicit none
    include '../FB/cbc_f_params.h'
    integer                                                     ::  i, j, k, g, bvx, odr, is, ie, js, je, ks, ke
    integer, dimension(3)                                       ::  sz, st, ed
    real                                                        ::  Up0, Ue1, Uw1, Us1, Un1, Ub1, Ut1
    real                                                        ::  Vp0, Ve1, Vw1, Vs1, Vn1, Vb1, Vt1
    real                                                        ::  Wp0, We1, Ww1, Ws1, Wn1, Wb1, Wt1
    real                                                        ::  c_e, c_w, c_n, c_s, c_t, c_b
    real                                                        ::  dh, dh1, dh2, flop, EX, EY, EZ, rei
    real                                                        ::  u_ref, v_ref, w_ref, m1, m2
    real                                                        ::  cnv_u, cnv_v, cnv_w, cr, cl, acr, acl
    real                                                        ::  u_bc, v_bc, w_bc, u_bc_ref, v_bc_ref, w_bc_ref
    real                                                        ::  fu_r, fu_l, fv_r, fv_l, fw_r, fw_l
    real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  v, wv
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bv
    real, dimension(0:3)                                        ::  v00
    real, dimension(3)                                          ::  vec
    
    is = st(1)
    ie = ed(1)
    js = st(2)
    je = ed(2)
    ks = st(3)
    ke = ed(3)
    
    dh1= 1.0/dh
    dh2= rei*dh1*dh1

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
    
    m1 = 0.0
    m2 = 0.0
    
    flop = flop + 13.0 ! DP 18 flop

!$OMP PARALLEL &
!$OMP REDUCTION(+:m1) &
!$OMP REDUCTION(+:m2) &
!$OMP FIRSTPRIVATE(is, ie, js, je, ks, ke, u_bc, v_bc, w_bc, u_bc_ref, v_bc_ref, w_bc_ref) &
!$OMP FIRSTPRIVATE(dh1, dh2, odr) &
!$OMP PRIVATE(bvx, cnv_u, cnv_v, cnv_w, EX, EY, EZ, cr, cl, acr, acl) &
!$OMP PRIVATE(Up0, Ue1, Uw1, Us1, Un1, Ub1, Ut1) &
!$OMP PRIVATE(Vp0, Ve1, Vw1, Vs1, Vn1, Vb1, Vt1) &
!$OMP PRIVATE(Wp0, We1, Ww1, Ws1, Wn1, Wb1, Wt1) &
!$OMP PRIVATE(c_e, c_w, c_n, c_s, c_t, c_b) &
!$OMP PRIVATE(fu_r, fu_l, fv_r, fv_l, fw_r, fw_l)
    
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
    
    do k=ks,ke
    do j=js,je
    do i=is,ie
      bvx = bv(i,j,k)

      if ( 0 /= iand(bvx, bc_mask30) ) then ! 6面のうちのどれか速度境界フラグが立っている場合
        cnv_u = 0.0
        cnv_v = 0.0
        cnv_w = 0.0
        
        ! 変数のロード
        Up0 = v(1,i,j,k)
        Vp0 = v(2,i,j,k)
        Wp0 = v(3,i,j,k)
      
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
          Uw1  = u_bc_ref
          Vw1  = v_bc_ref
          Ww1  = w_bc_ref
          cl   = u_bc
          acl  = abs(cl)
          fu_l = 0.5*(cl*(Up0+Uw1) - acl*(Up0-Uw1))
          fv_l = 0.5*(cl*(Vp0+Vw1) - acl*(Vp0-Vw1))
          fw_l = 0.5*(cl*(Wp0+Ww1) - acl*(Wp0-Ww1)) ! 18+1 flops
          m2 = m2 + 1.0
        end if
        
        if ( c_e == 1.0 ) then
          Ue1  = u_bc_ref
          Ve1  = v_bc_ref
          We1  = w_bc_ref
          cr   = u_bc
          acr  = abs(cr)
          fu_r = 0.5*(cr*(Ue1+Up0) - acr*(Ue1-Up0))
          fv_r = 0.5*(cr*(Ve1+Vp0) - acr*(Ve1-Vp0))
          fw_r = 0.5*(cr*(We1+Wp0) - acr*(We1-Wp0)) ! 18+1 flops
          m2 = m2 + 1.0
        end if
        
        cnv_u = cnv_u + fu_r*c_e - fu_l*c_w
        cnv_v = cnv_v + fv_r*c_e - fv_l*c_w
        cnv_w = cnv_w + fw_r*c_e - fw_l*c_w ! 12 flops
			
        ! Y方向 ---------------------------------------
        if ( c_s == 1.0 ) then
          Us1  = u_bc_ref
          Vs1  = v_bc_ref
          Ws1  = w_bc_ref
          cl   = v_bc
          acl  = abs(cl)
          fu_l = 0.5*(cl*(Up0+Us1) - acl*(Up0-Us1))
          fv_l = 0.5*(cl*(Vp0+Vs1) - acl*(Vp0-Vs1))
          fw_l = 0.5*(cl*(Wp0+Ws1) - acl*(Wp0-Ws1))
          m2 = m2 + 1.0
        end if
        
        if ( c_n == 1.0 ) then
          Un1  = u_bc_ref
          Vn1  = v_bc_ref
          Wn1  = w_bc_ref
          cr   = v_bc
          acr  = abs(cr)
          fu_r = 0.5*(cr*(Un1+Up0) - acr*(Un1-Up0))
          fv_r = 0.5*(cr*(Vn1+Vp0) - acr*(Vn1-Vp0))
          fw_r = 0.5*(cr*(Wn1+Wp0) - acr*(Wn1-Wp0))
          m2 = m2 + 1.0
        end if
        
        cnv_u = cnv_u + fu_r*c_n - fu_l*c_s
        cnv_v = cnv_v + fv_r*c_n - fv_l*c_s
        cnv_w = cnv_w + fw_r*c_n - fw_l*c_s
			
        ! Z方向 ---------------------------------------
        if ( c_b == 1.0 ) then
          Ub1  = u_bc_ref
          Vb1  = v_bc_ref
          Wb1  = w_bc_ref
          cl   = w_bc
          acl  = abs(cl)
          fu_l = 0.5*(cl*(Up0+Ub1) - acl*(Up0-Ub1))
          fv_l = 0.5*(cl*(Vp0+Vb1) - acl*(Vp0-Vb1))
          fw_l = 0.5*(cl*(Wp0+Wb1) - acl*(Wp0-Wb1))
          m2 = m2 + 1.0
        end if

        if ( c_t == 1.0 ) then
          Ut1  = u_bc_ref
          Vt1  = v_bc_ref
          Wt1  = w_bc_ref
          cr   = w_bc
          acr  = abs(cr)
          fu_r = 0.5*(cr*(Ut1+Up0) - acr*(Ut1-Up0))
          fv_r = 0.5*(cr*(Vt1+Vp0) - acr*(Vt1-Vp0))
          fw_r = 0.5*(cr*(Wt1+Wp0) - acr*(Wt1-Wp0))
          m2 = m2 + 1.0
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
           + ( Ub1 - Up0 ) * c_b ! 17 flops
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

        wv(1,i,j,k) = wv(1,i,j,k) + ( -cnv_u*dh1 + EX*dh2 )
        wv(2,i,j,k) = wv(2,i,j,k) + ( -cnv_v*dh1 + EY*dh2 )
        wv(3,i,j,k) = wv(3,i,j,k) + ( -cnv_w*dh1 + EZ*dh2 ) ! 4*3 = 12 flops
        m1 = m1 + 1.0
        
      endif
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    ! loop :  (12*3 + 17*3 + 12)*m1 + 19*m2 = 99*m1 + 19*m2

    flop = flop + m1*99.0 + m2*19.0

    return
    end subroutine cbc_pvec_vibc_specv


!  **********************************************************************************
!> @subroutine cbc_pvec_vobc_oflow (wv, sz, g, dh, v00, rei, v0, bv, vec, face, flop)
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
!! @param face 外部境界処理のときの面番号
!! @param[out] flop
!! @note vecには，流入条件のとき指定速度，流出境界の流束はローカルのセルフェイス速度を使うこと
!<
    subroutine cbc_pvec_vobc_oflow (wv, sz, g, dh, v00, rei, v0, bv, vec, face, flop)
    implicit none
    include '../FB/cbc_f_params.h'
    integer                                                     ::  i, j, k, g, bvx, face
    integer                                                     ::  ix, jx, kx
    integer, dimension(3)                                       ::  sz
    real                                                        ::  Up0, Ue1, Uw1, Us1, Un1, Ub1, Ut1
    real                                                        ::  Vp0, Ve1, Vw1, Vs1, Vn1, Vb1, Vt1
    real                                                        ::  Wp0, We1, Ww1, Ws1, Wn1, Wb1, Wt1
    real                                                        ::  dh, dh1, dh2, flop, rei
    real                                                        ::  u_ref, v_ref, w_ref, m
    real                                                        ::  u_bc, v_bc, w_bc
    real                                                        ::  fu, fv, fw, c, EX, EY, EZ
    real                                                        ::  Ue0, Uw0, Vs0, Vn0, Wb0, Wt0
    real                                                        ::  Ue, Uw, Vn, Vs, Wt, Wb
    real                                                        ::  b_w, b_e, b_s, b_n, b_b, b_t
    real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  v0, wv
    real, dimension(0:3)                                        ::  v00
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bv
    real, dimension(3)                                          ::  vec
    
    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    
    dh1= 1.0/dh
    dh2= rei*dh1*dh1

    ! 参照座標系の速度
    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)
    
    ! u_bcは境界速度
    u_bc = vec(1)
    v_bc = vec(2)
    w_bc = vec(3)
    
    flop = flop + 10.0 ! DP 15 flops

    m = 0.0
    
!$OMP PARALLEL REDUCTION(+:m) &
!$OMP FIRSTPRIVATE(ix, jx, kx, face) &
!$OMP FIRSTPRIVATE(u_ref, v_ref, w_ref, dh1, dh2) &
!$OMP PRIVATE(i, j, k, bvx) &
!$OMP PRIVATE(Up0, Ue1, Uw1, Us1, Un1, Ub1, Ut1) &
!$OMP PRIVATE(Vp0, Ve1, Vw1, Vs1, Vn1, Vb1, Vt1) &
!$OMP PRIVATE(Wp0, We1, Ww1, Ws1, Wn1, Wb1, Wt1) &
!$OMP PRIVATE(b_w, b_e, b_s, b_n, b_b, b_t) &
!$OMP PRIVATE(Ue, Uw, Vn, Vs, Wt, Wb) &
!$OMP PRIVATE(Ue0, Uw0, Vs0, Vn0, Wb0, Wt0) &
!$OMP PRIVATE(fu, fv, fw, c, EX, EY, EZ)
    
    FACES : select case (face)
    case (X_minus)
      i = 1
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do k=1,kx
      do j=1,jx
        bvx = bv(i,j,k)
        
        if ( ibits(bvx, bc_face_W, bitw_5) == obc_mask ) then
          Up0 = v0(1,i,j,k)
          Vp0 = v0(2,i,j,k)
          Wp0 = v0(3,i,j,k)
          
          include 'd_o_o_p.h' ! 42 flop
          
          Uw1 = v0(i-1,j  ,k  ,1)
          Vw1 = v0(i-1,j  ,k  ,2)
          Ww1 = v0(i-1,j  ,k  ,3)
          c   = Ue + (Vn - Vs + Wt - Wb) - u_ref
          if ( c>0.0 ) c=0.0
          fu  = c * Up0
          fv  = c * Vp0
          fw  = c * Wp0
          
          EX = Uw1 - Up0
          EY = Vw1 - Vp0
          EZ = Ww1 - Wp0

          wv(1,i,j,k) = wv(1,i,j,k) + ( -fu*dh1 + EX*dh2 )
          wv(2,i,j,k) = wv(2,i,j,k) + ( -fv*dh1 + EY*dh2 )
          wv(3,i,j,k) = wv(3,i,j,k) + ( -fw*dh1 + EZ*dh2 )
          m = m + 1.0
        endif
      end do
      end do
!$OMP END DO

      
    case (X_plus)
      i = ix
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do k=1,kx
      do j=1,jx
        bvx = bv(i,j,k)
        
        if ( ibits(bvx, bc_face_E, bitw_5) == obc_mask ) then
          Up0 = v0(1,i,j,k)
          Vp0 = v0(2,i,j,k)
          Wp0 = v0(3,i,j,k)
          
          include 'd_o_o_p.h'
          
          Ue1 = v0(i+1,j  ,k  ,1)
          Ve1 = v0(i+1,j  ,k  ,2)
          We1 = v0(i+1,j  ,k  ,3)
          c   = Uw - (Vn - Vs + Wt - Wb) - u_ref
          if ( c<0.0 ) c=0.0
          fu  = c * Up0
          fv  = c * Vp0
          fw  = c * Wp0

          EX = Ue1 - Up0
          EY = Ve1 - Vp0
          EZ = We1 - Wp0

          wv(1,i,j,k) = wv(1,i,j,k) + ( -fu*dh1 + EX*dh2 )
          wv(2,i,j,k) = wv(2,i,j,k) + ( -fv*dh1 + EY*dh2 )
          wv(3,i,j,k) = wv(3,i,j,k) + ( -fw*dh1 + EZ*dh2 )
          m = m + 1.0
        endif
      end do
      end do
!$OMP END DO

      
    case (Y_minus)
      j = 1
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do k=1,kx
      do i=1,ix
        bvx = bv(i,j,k)
        
        if ( ibits(bvx, bc_face_S, bitw_5) == obc_mask ) then
          Up0 = v0(1,i,j,k)
          Vp0 = v0(2,i,j,k)
          Wp0 = v0(3,i,j,k)
            
          include 'd_o_o_p.h'
          
          Us1 = v0(i  ,j-1,k  ,1)
          Vs1 = v0(i  ,j-1,k  ,2)
          Ws1 = v0(i  ,j-1,k  ,3)
          c   = Vn + (Ue - Uw + Wt - Wb) - v_ref
          if ( c>0.0 ) c=0.0
          fu  = c * Up0
          fv  = c * Vp0
          fw  = c * Wp0
        
          EX = Us1 - Up0
          EY = Vs1 - Vp0
          EZ = Ws1 - Wp0
          
          wv(1,i,j,k) = wv(1,i,j,k) + ( fu*dh1 + EX*dh2 )
          wv(2,i,j,k) = wv(2,i,j,k) + ( fv*dh1 + EY*dh2 )
          wv(3,i,j,k) = wv(3,i,j,k) + ( fw*dh1 + EZ*dh2 )
          m = m + 1.0
        endif
      end do
      end do
!$OMP END DO

      
    case (Y_plus)
      j = jx
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do k=1,kx
      do i=1,ix
        bvx = bv(i,j,k)
        
        if ( ibits(bvx, bc_face_N, bitw_5) == obc_mask ) then
          Up0 = v0(1,i,j,k)
          Vp0 = v0(2,i,j,k)
          Wp0 = v0(3,i,j,k)
            
          include 'd_o_o_p.h'
            
          Un1 = v0(i  ,j+1,k  ,1)
          Vn1 = v0(i  ,j+1,k  ,2)
          Wn1 = v0(i  ,j+1,k  ,3)
          c   = Vs - (Ue - Uw + Wt - Wb) - v_ref
          if ( c<0.0 ) c=0.0
          fu  = c * Up0
          fv  = c * Vp0
          fw  = c * Wp0

          EX = Un1 - Up0
          EY = Vn1 - Vp0
          EZ = Wn1 - Wp0
          
          wv(1,i,j,k) = wv(1,i,j,k) + ( -fu*dh1 + EX*dh2 )
          wv(2,i,j,k) = wv(2,i,j,k) + ( -fv*dh1 + EY*dh2 )
          wv(3,i,j,k) = wv(3,i,j,k) + ( -fw*dh1 + EZ*dh2 )
          m = m + 1.0
        endif
      end do
      end do
!$OMP END DO

      
    case (Z_minus)
      k = 1
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do j=1,jx
      do i=1,ix
        bvx = bv(i,j,k)
        
        if ( ibits(bvx, bc_face_B, bitw_5) == obc_mask ) then
          Up0 = v0(1,i,j,k)
          Vp0 = v0(2,i,j,k)
          Wp0 = v0(3,i,j,k)

          include 'd_o_o_p.h'
          
          Ub1 = v0(i  ,j  ,k-1,1)
          Vb1 = v0(i  ,j  ,k-1,2)
          Wb1 = v0(i  ,j  ,k-1,3)
          c   = Wt + (Ue - Uw + Vn - Vs) - w_ref
          if ( c>0.0 ) c=0.0
          fu  = c * Up0
          fv  = c * Vp0
          fw  = c * Wp0
          
          EX = Ub1 - Up0
          EY = Vb1 - Vp0
          EZ = Wb1 - Wp0
          
          wv(1,i,j,k) = wv(1,i,j,k) + ( fu*dh1 + EX*dh2 )
          wv(2,i,j,k) = wv(2,i,j,k) + ( fv*dh1 + EY*dh2 )
          wv(3,i,j,k) = wv(3,i,j,k) + ( fw*dh1 + EZ*dh2 )
          m = m + 1.0
        endif
      end do
      end do
!$OMP END DO

      
    case (Z_plus)
      k = kx
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do j=1,jx
      do i=1,ix
        bvx = bv(i,j,k)
        
        if ( ibits(bvx, bc_face_T, bitw_5) == obc_mask ) then
          Up0 = v0(1,i,j,k)
          Vp0 = v0(2,i,j,k)
          Wp0 = v0(3,i,j,k)

          include 'd_o_o_p.h'
          
          Ut1 = v0(i  ,j  ,k+1,1)
          Vt1 = v0(i  ,j  ,k+1,2)
          Wt1 = v0(i  ,j  ,k+1,3)
          c   = Wb - (Ue - Uw + Vn - Vs) - w_ref
          if ( c<0.0 ) c=0.0
          fu  = c * Up0
          fv  = c * Vp0
          fw  = c * Wp0
          
          EX = Ut1 - Up0
          EY = Vt1 - Vp0
          EZ = Wt1 - Wp0
          
          wv(1,i,j,k) = wv(1,i,j,k) + ( -fu*dh1 + EX*dh2 )
          wv(2,i,j,k) = wv(2,i,j,k) + ( -fv*dh1 + EY*dh2 )
          wv(3,i,j,k) = wv(3,i,j,k) + ( -fw*dh1 + EZ*dh2 )
          m = m + 1.0
        endif
      end do
      end do
!$OMP END DO

      
    case default
    end select FACES
    
!$OMP END PARALLEL

    flop = flop + m*68.0
      
    return
    end subroutine cbc_pvec_vobc_oflow
    
!  **********************************************************************************
!> @subroutine cbc_pvec_vobc_specv (wv, sz, g, dh, v00, rei, v0, bv, vec, face, flop)
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
!! @param face 外部境界処理のときの面番号
!! @param[out] flop
!! @note vecには，流入条件のとき指定速度，流出境界の流束はローカルのセルフェイス速度を使うこと
!<
    subroutine cbc_pvec_vobc_specv (wv, sz, g, dh, v00, rei, v0, bv, vec, face, flop)
    implicit none
    include '../FB/cbc_f_params.h'
    integer                                                     ::  i, j, k, g, bvx, face
    integer                                                     ::  ix, jx, kx
    integer, dimension(3)                                       ::  sz
    real                                                        ::  Up0, Ue1, Uw1, Us1, Un1, Ub1, Ut1
    real                                                        ::  Vp0, Ve1, Vw1, Vs1, Vn1, Vb1, Vt1
    real                                                        ::  Wp0, We1, Ww1, Ws1, Wn1, Wb1, Wt1
    real                                                        ::  dh, dh1, dh2, flop, EX, EY, EZ, rei
    real                                                        ::  u_ref, v_ref, w_ref, m
    real                                                        ::  fu, fv, fw, c, ac
    real                                                        ::  u_bc, v_bc, w_bc, u_bc_ref, v_bc_ref, w_bc_ref
    real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  v0, wv
    real, dimension(0:3)                                        ::  v00
    real, dimension(3)                                          ::  vec
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bv
    
    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    
    dh1= 1.0/dh
    dh2= rei*dh1*dh1

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
    
    flop = flop + 13.0 ! DP 18 flop

    m = 0.0
    
!$OMP PARALLEL REDUCTION(+:m) &
!$OMP FIRSTPRIVATE(ix, jx, kx, u_bc_ref, v_bc_ref, w_bc_ref, face) &
!$OMP FIRSTPRIVATE(dh1, dh2, u_bc, v_bc, w_bc) &
!$OMP PRIVATE(i, j, k, bvx) &
!$OMP PRIVATE(Up0, Ue1, Uw1, Us1, Un1, Ub1, Ut1) &
!$OMP PRIVATE(Vp0, Ve1, Vw1, Vs1, Vn1, Vb1, Vt1) &
!$OMP PRIVATE(Wp0, We1, Ww1, Ws1, Wn1, Wb1, Wt1) &
!$OMP PRIVATE(fu, fv, fw, EX, EY, EZ, c, ac)

    FACES : select case (face)
    case (X_minus)
      i = 1
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do k=1,kx
      do j=1,jx
        bvx = bv(i,j,k)
        
        if ( ibits(bvx, bc_face_W, bitw_5) == obc_mask ) then
          Up0 = v0(1,i,j,k)
          Vp0 = v0(2,i,j,k)
          Wp0 = v0(3,i,j,k)
          
          Uw1 = u_bc_ref
          Vw1 = v_bc_ref
          Ww1 = w_bc_ref
          c   = u_bc
          ac  = abs(c)
          fu  = 0.5*(c*(Up0+Uw1) - ac*(Up0-Uw1))
          fv  = 0.5*(c*(Vp0+Vw1) - ac*(Vp0-Vw1))
          fw  = 0.5*(c*(Wp0+Ww1) - ac*(Wp0-Ww1))
          
          EX = Uw1 - Up0
          EY = Vw1 - Vp0
          EZ = Ww1 - Wp0

          wv(1,i,j,k) = wv(1,i,j,k) + ( fu*dh1 + EX*dh2 )
          wv(2,i,j,k) = wv(2,i,j,k) + ( fv*dh1 + EY*dh2 )
          wv(3,i,j,k) = wv(3,i,j,k) + ( fw*dh1 + EZ*dh2 )
          m = m + 1.0
        endif
      end do
      end do
!$OMP END DO

      
    case (X_plus)
      i = ix
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do k=1,kx
      do j=1,jx
        bvx = bv(i,j,k)
        
        if ( ibits(bvx, bc_face_E, bitw_5) == obc_mask ) then ! 方向によって実装が異なるのでチェック
          Up0 = v0(1,i,j,k)
          Vp0 = v0(2,i,j,k)
          Wp0 = v0(3,i,j,k)
          
          Ue1 = u_bc_ref
          Ve1 = v_bc_ref
          We1 = w_bc_ref
          c   = u_bc
          ac  = abs(c)
          fu  = 0.5*(c*(Ue1+Up0) - ac*(Ue1-Up0))
          fv  = 0.5*(c*(Ve1+Vp0) - ac*(Ve1-Vp0))
          fw  = 0.5*(c*(We1+Wp0) - ac*(We1-Wp0))

          EX = Ue1 - Up0
          EY = Ve1 - Vp0
          EZ = We1 - Wp0

          wv(1,i,j,k) = wv(1,i,j,k) + ( -fu*dh1 + EX*dh2 )
          wv(2,i,j,k) = wv(2,i,j,k) + ( -fv*dh1 + EY*dh2 )
          wv(3,i,j,k) = wv(3,i,j,k) + ( -fw*dh1 + EZ*dh2 )
          m = m + 1.0
        endif
      end do
      end do
!$OMP END DO

      
    case (Y_minus)
      j = 1
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do k=1,kx
      do i=1,ix
        bvx = bv(i,j,k)
        
        if ( ibits(bvx, bc_face_S, bitw_5) == obc_mask ) then
          Up0 = v0(1,i,j,k)
          Vp0 = v0(2,i,j,k)
          Wp0 = v0(3,i,j,k)
          
          Us1 = u_bc_ref
          Vs1 = v_bc_ref
          Ws1 = w_bc_ref
          c   = v_bc
          ac  = abs(c)
          fu  = 0.5*(c*(Up0+Us1) - ac*(Up0-Us1))
          fv  = 0.5*(c*(Vp0+Vs1) - ac*(Vp0-Vs1))
          fw  = 0.5*(c*(Wp0+Ws1) - ac*(Wp0-Ws1))
        
          EX = Us1 - Up0
          EY = Vs1 - Vp0
          EZ = Ws1 - Wp0
          
          wv(1,i,j,k) = wv(1,i,j,k) + ( fu*dh1 + EX*dh2 )
          wv(2,i,j,k) = wv(2,i,j,k) + ( fv*dh1 + EY*dh2 )
          wv(3,i,j,k) = wv(3,i,j,k) + ( fw*dh1 + EZ*dh2 )
          m = m + 1.0
        endif
      end do
      end do
!$OMP END DO

      
    case (Y_plus)
      j = jx
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do k=1,kx
      do i=1,ix
        bvx = bv(i,j,k)
        
        if ( ibits(bvx, bc_face_N, bitw_5) == obc_mask ) then
          Up0 = v0(1,i,j,k)
          Vp0 = v0(2,i,j,k)
          Wp0 = v0(3,i,j,k)
          
          Un1 = u_bc_ref
          Vn1 = v_bc_ref
          Wn1 = w_bc_ref
          c   = v_bc
          ac  = abs(c)
          fu  = 0.5*(c*(Un1+Up0) - ac*(Un1-Up0))
          fv  = 0.5*(c*(Vn1+Vp0) - ac*(Vn1-Vp0))
          fw  = 0.5*(c*(Wn1+Wp0) - ac*(Wn1-Wp0))

          EX = Un1 - Up0
          EY = Vn1 - Vp0
          EZ = Wn1 - Wp0
          
          wv(1,i,j,k) = wv(1,i,j,k) + ( -fu*dh1 + EX*dh2 )
          wv(2,i,j,k) = wv(2,i,j,k) + ( -fv*dh1 + EY*dh2 )
          wv(3,i,j,k) = wv(3,i,j,k) + ( -fw*dh1 + EZ*dh2 )
          m = m + 1.0
        endif
      end do
      end do
!$OMP END DO

      
    case (Z_minus)
      k = 1
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do j=1,jx
      do i=1,ix
        bvx = bv(i,j,k)
        
        if ( ibits(bvx, bc_face_B, bitw_5) == obc_mask ) then
          Up0 = v0(1,i,j,k)
          Vp0 = v0(2,i,j,k)
          Wp0 = v0(3,i,j,k)
          
          Ub1 = u_bc_ref
          Vb1 = v_bc_ref
          Wb1 = w_bc_ref
          c   = w_bc
          ac  = abs(c)
          fu  = 0.5*(c*(Up0+Ub1) - ac*(Up0-Ub1))
          fv  = 0.5*(c*(Vp0+Vb1) - ac*(Vp0-Vb1))
          fw  = 0.5*(c*(Wp0+Wb1) - ac*(Wp0-Wb1))
          
          EX = Ub1 - Up0
          EY = Vb1 - Vp0
          EZ = Wb1 - Wp0
          
          wv(1,i,j,k) = wv(1,i,j,k) + ( fu*dh1 + EX*dh2 )
          wv(2,i,j,k) = wv(2,i,j,k) + ( fv*dh1 + EY*dh2 )
          wv(3,i,j,k) = wv(3,i,j,k) + ( fw*dh1 + EZ*dh2 )
          m = m + 1.0
        endif
      end do
      end do
!$OMP END DO
      
    case (Z_plus)
      k = kx
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do j=1,jx
      do i=1,ix
        bvx = bv(i,j,k)
        
        if ( ibits(bvx, bc_face_T, bitw_5) == obc_mask ) then
          Up0 = v0(1,i,j,k)
          Vp0 = v0(2,i,j,k)
          Wp0 = v0(3,i,j,k)
          
          Ut1 = u_bc_ref
          Vt1 = v_bc_ref
          Wt1 = w_bc_ref
          c   = w_bc
          ac  = abs(c)
          fu  = 0.5*(c*(Ut1+Up0) - ac*(Ut1-Up0))
          fv  = 0.5*(c*(Vt1+Vp0) - ac*(Vt1-Vp0))
          fw  = 0.5*(c*(Wt1+Wp0) - ac*(Wt1-Wp0))
          
          EX = Ut1 - Up0
          EY = Vt1 - Vp0
          EZ = Wt1 - Wp0
          
          wv(1,i,j,k) = wv(1,i,j,k) + ( -fu*dh1 + EX*dh2 )
          wv(2,i,j,k) = wv(2,i,j,k) + ( -fv*dh1 + EY*dh2 )
          wv(3,i,j,k) = wv(3,i,j,k) + ( -fw*dh1 + EZ*dh2 )
          m = m + 1.0
        endif
      end do
      end do
!$OMP END DO
      
    case default
    end select FACES
    
!$OMP END PARALLEL

    flop = flop + m*34.0
    
    return
    end subroutine cbc_pvec_vobc_specv
    
!  *************************************************************************
!> @subroutine cbc_pvec_vobc_symtrc (wv, sz, g, dh, rei, v0, bv, face, flop)
!! @brief 外部速度境界条件による対流項と粘性項の流束の修正
!! @param[out] wv 疑似ベクトルの空間項の評価値
!! @param sz 配列長
!! @param g ガイドセル長
!! @param dh 格子幅
!! @param rei Reynolds数の逆数
!! @param v0 速度ベクトル（n-step）
!! @param bv BCindex V
!! @param face 外部境界処理のときの面番号
!! @param[out] flop
!! @note 境界面で対流流束はゼロ，粘性流束のみ
!<
    subroutine cbc_pvec_vobc_symtrc (wv, sz, g, dh, rei, v0, bv, face, flop)
    implicit none
    include '../FB/cbc_f_params.h'
    integer                                                     ::  i, j, k, g, bvx, face
    integer                                                     ::  ix, jx, kx
    integer, dimension(3)                                       ::  sz
    real                                                        ::  dh, dh1, dh2, flop, rei
    real                                                        ::  rix, rjx, rkx, Up0, Vp0, Wp0
    real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  v0, wv
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bv
    
    ix = sz(1)
    jx = sz(2)
    kx = sz(3)

    rix = real(jx)*real(kx)
    rjx = real(ix)*real(kx)
    rkx = real(ix)*real(jx)
    
    dh1= 1.0/dh
    dh2= 2.0*rei*dh1*dh1
    
    flop = flop + 14.0 ! DP 19 flop

!$OMP PARALLEL REDUCTION(+:flop) &
!$OMP FIRSTPRIVATE(ix, jx, kx, dh2, face) &
!$OMP PRIVATE(i, j, k, bvx, Up0, Vp0, Wp0)
    
    FACES : select case (face)
    
    case (X_minus)
      i = 1
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do k=1,kx
      do j=1,jx
        bvx = bv(i,j,k)
        
        if ( ibits(bvx, bc_face_W, bitw_5) == obc_mask ) then
          Up0 = v0(1,i,j,k)
          wv(1,i,j,k) = wv(1,i,j,k) + Up0*dh2
        endif
      end do
      end do
!$OMP END DO

      flop = flop + rix*2.0
      

    case (X_plus)
      i = ix
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do k=1,kx
      do j=1,jx
        bvx = bv(i,j,k)
        
        if ( ibits(bvx, bc_face_E, bitw_5) == obc_mask ) then
          Up0 = v0(1,i,j,k)
          wv(1,i,j,k) = wv(1,i,j,k) - Up0*dh2
        endif
      end do
      end do
!$OMP END DO
    
      flop = flop + rix*2.0
      

    case (Y_minus)
      j = 1
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do k=1,kx
      do i=1,ix
        bvx = bv(i,j,k)
        
        if ( ibits(bvx, bc_face_S, bitw_5) == obc_mask ) then
          Vp0 = v0(2,i,j,k)
          wv(2,i,j,k) = wv(2,i,j,k) + Vp0*dh2
        endif
      end do
      end do
!$OMP END DO
      
      flop = flop + rjx*2.0
      

    case (Y_plus)
      j = jx
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do k=1,kx
      do i=1,ix
        bvx = bv(i,j,k)
        
        if ( ibits(bvx, bc_face_N, bitw_5) == obc_mask ) then
          Vp0 = v0(2,i,j,k)
          wv(2,i,j,k) = wv(2,i,j,k) - Vp0*dh2
        endif
      end do
      end do
!$OMP END DO
      
      flop = flop + rjx*2.0
      

    case (Z_minus)
      k = 1
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do j=1,jx
      do i=1,ix
        bvx = bv(i,j,k)
        
        if ( ibits(bvx, bc_face_B, bitw_5) == obc_mask ) then
          Wp0 = v0(3,i,j,k)
          wv(3,i,j,k) = wv(3,i,j,k) + Wp0*dh2
        endif
      end do
      end do
!$OMP END DO
      
      flop = flop + rkx*2.0
      

    case (Z_plus)
      k = kx
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do j=1,jx
      do i=1,ix
        bvx = bv(i,j,k)
        
        if ( ibits(bvx, bc_face_T, bitw_5) == obc_mask ) then
          Wp0 = v0(3,i,j,k)
          wv(3,i,j,k) = wv(3,i,j,k) - Wp0*dh2
        endif
      end do
      end do
!$OMP END DO
      
      flop = flop + rkx*2.0
      

    case default
    end select FACES
    
!$OMP END PARALLEL

    return
    end subroutine cbc_pvec_vobc_symtrc
    
!  *********************************************************************************
!> @subroutine cbc_pvec_vobc_wall (wv, sz, g, dh, v00, rei, v0, bv, vec, face, flop)
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
!! @param face 外部境界処理のときの面番号
!! @param[out] flop
!! @note vecには，流入条件のとき指定速度，流出境界の流束はローカルのセルフェイス速度を使うこと
!! @todo 内部と外部の分離 do loopの内側に条件分岐を入れているので修正
!<
    subroutine cbc_pvec_vobc_wall  (wv, sz, g, dh, v00, rei, v0, bv, vec, face, flop)
    implicit none
    include '../FB/cbc_f_params.h'
    integer                                                     ::  i, j, k, g, bvx, face
    integer                                                     ::  ix, jx, kx
    integer, dimension(3)                                       ::  sz
    real                                                        ::  Up0, Ue1, Uw1, Us1, Un1, Ub1, Ut1
    real                                                        ::  Vp0, Ve1, Vw1, Vs1, Vn1, Vb1, Vt1
    real                                                        ::  Wp0, We1, Ww1, Ws1, Wn1, Wb1, Wt1
    real                                                        ::  dh, dh1, dh2, flop, EX, EY, EZ, rei
    real                                                        ::  u_ref, v_ref, w_ref, m
    real                                                        ::  u_bc, v_bc, w_bc, u_bc_ref, v_bc_ref, w_bc_ref
    real                                                        ::  u_bc_ref2, v_bc_ref2, w_bc_ref2
    real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  v0, wv
    real, dimension(0:3)                                        ::  v00
    real, dimension(3)                                          ::  vec
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bv
    
    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    
    dh1= 1.0/dh
    dh2= rei*dh1*dh1

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
    
    u_bc_ref2 = 2.0*u_bc_ref
    v_bc_ref2 = 2.0*v_bc_ref
    w_bc_ref2 = 2.0*w_bc_ref
    
    flop = flop + 16.0 ! DP 21 flops

    m = 0.0

!$OMP PARALLEL REDUCTION(+:m) &
!$OMP FIRSTPRIVATE(ix, jx, kx, u_bc_ref2, v_bc_ref2, w_bc_ref2, dh2, face) &
!$OMP PRIVATE(i, j, k, bvx, Up0, Vp0, Wp0, EX, EY, EZ) &
!$OMP PRIVATE(Ue1, Uw1, Us1, Un1, Ub1, Ut1) &
!$OMP PRIVATE(Ve1, Vw1, Vs1, Vn1, Vb1, Vt1) &
!$OMP PRIVATE(We1, Ww1, Ws1, Wn1, Wb1, Wt1)

    FACES : select case (face)
    
    case (X_minus)
      i = 1
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do k=1,kx
      do j=1,jx
        bvx = bv(i,j,k)
        
        if ( ibits(bvx, bc_face_W, bitw_5) == obc_mask ) then
          Up0 = v0(1,i,j,k)
          Vp0 = v0(2,i,j,k)
          Wp0 = v0(3,i,j,k)
          
          Uw1  = u_bc_ref2 - Up0
          Vw1  = v_bc_ref2 - Vp0
          Ww1  = w_bc_ref2 - Wp0
          
          EX = Uw1 - Up0
          EY = Vw1 - Vp0
          EZ = Ww1 - Wp0

          wv(1,i,j,k) = wv(1,i,j,k) + EX*dh2
          wv(2,i,j,k) = wv(2,i,j,k) + EY*dh2
          wv(3,i,j,k) = wv(3,i,j,k) + EZ*dh2
          m = m + 1.0
        endif
      end do
      end do
!$OMP END DO
      
      
    case (X_plus)
      i = ix
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do k=1,kx
      do j=1,jx
        bvx = bv(i,j,k)
        
        if ( ibits(bvx, bc_face_E, bitw_5) == obc_mask ) then
          Up0 = v0(1,i,j,k)
          Vp0 = v0(2,i,j,k)
          Wp0 = v0(3,i,j,k)

          Ue1  = u_bc_ref2 - Up0
          Ve1  = v_bc_ref2 - Vp0
          We1  = w_bc_ref2 - Wp0

          EX = Ue1 - Up0
          EY = Ve1 - Vp0
          EZ = We1 - Wp0

          wv(1,i,j,k) = wv(1,i,j,k) + EX*dh2
          wv(2,i,j,k) = wv(2,i,j,k) + EY*dh2
          wv(3,i,j,k) = wv(3,i,j,k) + EZ*dh2
          m = m + 1.0
        endif
      end do
      end do
!$OMP END DO
      
      
    case (Y_minus)
      j = 1
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do k=1,kx
      do i=1,ix
        bvx = bv(i,j,k)
        
        if ( ibits(bvx, bc_face_S, bitw_5) == obc_mask ) then
          Up0 = v0(1,i,j,k)
          Vp0 = v0(2,i,j,k)
          Wp0 = v0(3,i,j,k)

          Us1  = u_bc_ref2 - Up0
          Vs1  = v_bc_ref2 - Vp0
          Ws1  = w_bc_ref2 - Wp0
        
          EX = Us1 - Up0
          EY = Vs1 - Vp0
          EZ = Ws1 - Wp0

          wv(1,i,j,k) = wv(1,i,j,k) + EX*dh2
          wv(2,i,j,k) = wv(2,i,j,k) + EY*dh2
          wv(3,i,j,k) = wv(3,i,j,k) + EZ*dh2
          m = m + 1.0
        endif
      end do
      end do
!$OMP END DO
      
      
    case (Y_plus)
      j = jx
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do k=1,kx
      do i=1,ix
        bvx = bv(i,j,k)
        
        if ( ibits(bvx, bc_face_N, bitw_5) == obc_mask ) then
          Up0 = v0(1,i,j,k)
          Vp0 = v0(2,i,j,k)
          Wp0 = v0(3,i,j,k)

          Un1  = u_bc_ref2 - Up0
          Vn1  = v_bc_ref2 - Vp0
          Wn1  = w_bc_ref2 - Wp0

          EX = Un1 - Up0
          EY = Vn1 - Vp0
          EZ = Wn1 - Wp0
          
          wv(1,i,j,k) = wv(1,i,j,k) + EX*dh2
          wv(2,i,j,k) = wv(2,i,j,k) + EY*dh2
          wv(3,i,j,k) = wv(3,i,j,k) + EZ*dh2
          m = m + 1.0
        endif
      end do
      end do
!$OMP END DO
      
      
    case (Z_minus)
      k = 1
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do j=1,jx
      do i=1,ix
        bvx = bv(i,j,k)
        
        if ( ibits(bvx, bc_face_B, bitw_5) == obc_mask ) then
          Up0 = v0(1,i,j,k)
          Vp0 = v0(2,i,j,k)
          Wp0 = v0(3,i,j,k)

          Ub1  = u_bc_ref2 - Up0
          Vb1  = v_bc_ref2 - Vp0
          Wb1  = w_bc_ref2 - Wp0
          
          EX = Ub1 - Up0
          EY = Vb1 - Vp0
          EZ = Wb1 - Wp0
          
          wv(1,i,j,k) = wv(1,i,j,k) + EX*dh2
          wv(2,i,j,k) = wv(2,i,j,k) + EY*dh2
          wv(3,i,j,k) = wv(3,i,j,k) + EZ*dh2
          m = m + 1.0
        endif
      end do
      end do
!$OMP END DO
      
      
    case (Z_plus)
      k = kx
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do j=1,jx
      do i=1,ix
        bvx = bv(i,j,k)
        
        if ( ibits(bvx, bc_face_T, bitw_5) == obc_mask ) then
          Up0 = v0(1,i,j,k)
          Vp0 = v0(2,i,j,k)
          Wp0 = v0(3,i,j,k)

          Ut1  = u_bc_ref2 - Up0
          Vt1  = v_bc_ref2 - Vp0
          Wt1  = w_bc_ref2 - Wp0
          
          EX = Ut1 - Up0
          EY = Vt1 - Vp0
          EZ = Wt1 - Wp0
          
          wv(1,i,j,k) = wv(1,i,j,k) + EX*dh2
          wv(2,i,j,k) = wv(2,i,j,k) + EY*dh2
          wv(3,i,j,k) = wv(3,i,j,k) + EZ*dh2
          m = m + 1.0
        endif
      end do
      end do
!$OMP END DO
  
    case default
    end select FACES
    
!$OMP END PARALLEL

    flop = flop + m*12.0
    
    return
    end subroutine cbc_pvec_vobc_wall


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
    include '../FB/cbc_f_params.h'
    integer                                                     ::  i, j, k, g, idx, odr, is, ie, js, je, ks, ke
    integer, dimension(3)                                       ::  sz, st, ed
    real                                                        ::  u_bc_ref, v_bc_ref, w_bc_ref
    real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  v
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bv
    real, dimension(0:3)                                        ::  v00
    real, dimension(3)                                          ::  vec
    
    ! u_bc_refは参照座標系での境界速度
    u_bc_ref = vec(1) + v00(1)
    v_bc_ref = vec(2) + v00(2)
    w_bc_ref = vec(3) + v00(3)
    
    is = st(1)
    ie = ed(1)
    js = st(2)
    je = ed(2)
    ks = st(3)
    ke = ed(3)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(is, ie, js, je, ks, ke, u_bc_ref, v_bc_ref, w_bc_ref, odr) &
!$OMP PRIVATE(idx)

#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
    do k=ks,ke
    do j=js,je
    do i=is,ie
      idx = bv(i,j,k)

      if ( 0 /= iand(idx, bc_mask30) ) then ! 6面のうちのどれか速度境界フラグが立っている場合

        if ( ibits(idx, bc_face_W, bitw_5) == odr ) then
          v(1, i-1,j,k) = u_bc_ref
          v(2, i-1,j,k) = v_bc_ref
          v(3, i-1,j,k) = w_bc_ref
        end if

        if ( ibits(idx, bc_face_E, bitw_5) == odr ) then
          v(1, i+1,j,k) = u_bc_ref
          v(2, i+1,j,k) = v_bc_ref
          v(3, i+1,j,k) = w_bc_ref
        end if

        if ( ibits(idx, bc_face_S, bitw_5) == odr ) then
          v(1, i,j-1,k) = u_bc_ref
          v(2, i,j-1,k) = v_bc_ref
          v(3, i,j-1,k) = w_bc_ref
        end if
        
        if ( ibits(idx, bc_face_N, bitw_5) == odr ) then
          v(1, i,j+1,k) = u_bc_ref
          v(2, i,j+1,k) = v_bc_ref
          v(3, i,j+1,k) = w_bc_ref
        end if
        
        if ( ibits(idx, bc_face_B, bitw_5) == odr ) then
          v(1, i,j,k-1) = u_bc_ref
          v(2, i,j,k-1) = v_bc_ref
          v(3, i,j,k-1) = w_bc_ref
        end if
			
        if ( ibits(idx, bc_face_T, bitw_5) == odr ) then
          v(1, i,j,k+1) = u_bc_ref
          v(2, i,j,k+1) = v_bc_ref
          v(3, i,j,k+1) = w_bc_ref
        end if

      endif
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL
    
    return
    end subroutine cbc_vibc_drchlt

!  ********************************************************
!> @subroutine cbc_vibc_outflow (v, sz, g, st, ed, bv, odr)
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
    include '../FB/cbc_f_params.h'
    integer                                                     ::  i, j, k, g, idx, odr, is, ie, js, je, ks, ke
    integer, dimension(3)                                       ::  sz, st, ed
    real                                                        ::  Up, Vp, Wp
    real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  v
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bv
    
    is = st(1)
    ie = ed(1)
    js = st(2)
    je = ed(2)
    ks = st(3)
    ke = ed(3)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(is, ie, js, je, ks, ke, odr) &
!$OMP PRIVATE(idx, Up, Vp, Wp)

#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
    do k=ks,ke
    do j=js,je
    do i=is,ie
    
      idx = bv(i,j,k)
      if ( 0 /= iand(idx, bc_mask30) ) then ! 6面のうちのどれか速度境界フラグが立っている場合

        Up = v(1,i,j,k)
        Vp = v(2,i,j,k)
        Wp = v(3,i,j,k)
        
        ! X-
        if ( ibits(idx, bc_face_W, bitw_5) == odr ) then
          v(1, i-1,j,k) = Up
          v(2, i-1,j,k) = Vp
          v(3, i-1,j,k) = Wp
        end if
        
        ! X+
        if ( ibits(idx, bc_face_E, bitw_5) == odr ) then
          v(1, i+1,j,k) = Up
          v(2, i+1,j,k) = Vp
          v(3, i+1,j,k) = Wp
        end if
        
        ! Y-
        if ( ibits(idx, bc_face_S, bitw_5) == odr ) then
          v(1, i,j-1,k) = Up
          v(2, i,j-1,k) = Vp
          v(3, i,j-1,k) = Wp
        end if
        
        ! Y+
        if ( ibits(idx, bc_face_N, bitw_5) == odr ) then
          v(1, i,j+1,k) = Up
          v(2, i,j+1,k) = Vp
          v(3, i,j+1,k) = Wp
        end if
        
        ! Z-
        if ( ibits(idx, bc_face_B, bitw_5) == odr ) then
          v(1, i,j,k-1) = Up
          v(2, i,j,k-1) = Vp
          v(3, i,j,k-1) = Wp
        end if
			  
        ! Z+
        if ( ibits(idx, bc_face_T, bitw_5) == odr ) then
          v(1, i,j,k+1) = Up
          v(2, i,j,k+1) = Vp
          v(3, i,j,k+1) = Wp
        end if
        
      endif
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL
    
    return
    end subroutine cbc_vibc_outflow


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
    include '../FB/cbc_f_params.h'
    integer                                                     ::  i, j, k, g, face, ix, jx, kx, bvx
    integer, dimension(3)                                       ::  sz
    real                                                        ::  u_bc_ref, v_bc_ref, w_bc_ref
    real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  v
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
    
!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, u_bc_ref, v_bc_ref, w_bc_ref, face) &
!$OMP PRIVATE(i, j, k, bvx)

    FACES : select case (face)
    case (X_minus)
      i = 1
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do k=1,kx
      do j=1,jx
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_W, bitw_5) == obc_mask ) then
          v(1, i-1,j,k) = u_bc_ref
          v(2, i-1,j,k) = v_bc_ref
          v(3, i-1,j,k) = w_bc_ref
        endif
      end do
      end do
!$OMP END DO
      
      
    case (X_plus)
      i = ix
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do k=1,kx
      do j=1,jx
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_E, bitw_5) == obc_mask ) then
          v(1, i-1,j,k) = u_bc_ref
          v(2, i-1,j,k) = v_bc_ref
          v(3, i-1,j,k) = w_bc_ref
        endif
      end do
      end do
!$OMP END DO
      
      
    case (Y_minus)
      j = 1
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do k=1,kx
      do i=1,ix
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_S, bitw_5) == obc_mask ) then
          v(1, i-1,j,k) = u_bc_ref
          v(2, i-1,j,k) = v_bc_ref
          v(3, i-1,j,k) = w_bc_ref
        endif
      end do
      end do
!$OMP END DO
      
      
    case (Y_plus)
      j = jx
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do k=1,kx
      do i=1,ix
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_N, bitw_5) == obc_mask ) then
          v(1, i-1,j,k) = u_bc_ref
          v(2, i-1,j,k) = v_bc_ref
          v(3, i-1,j,k) = w_bc_ref
        endif
      end do
      end do
!$OMP END DO
      
      
    case (Z_minus)
      k = 1
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do j=1,jx
      do i=1,ix
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_B, bitw_5) == obc_mask ) then
          v(1, i-1,j,k) = u_bc_ref
          v(2, i-1,j,k) = v_bc_ref
          v(3, i-1,j,k) = w_bc_ref
        endif
      end do
      end do
!$OMP END DO
      
      
    case (Z_plus)
      k = kx
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do j=1,jx
      do i=1,ix
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_T, bitw_5) == obc_mask ) then
          v(1, i-1,j,k) = u_bc_ref
          v(2, i-1,j,k) = v_bc_ref
          v(3, i-1,j,k) = w_bc_ref
        endif
      end do
      end do
!$OMP END DO
      
    case default
    end select FACES

!$OMP END PARALLEL
    
    return
    end subroutine cbc_vobc_drchlt

!  *********************************************
!> @subroutine cbc_vobc_neumann (v, sz, g, face)
!! @brief 遠方境界の近似
!! @param v 速度ベクトル
!! @param sz 配列長
!! @param g ガイドセル長
!! @param face 外部境界面の番号
!<
    subroutine cbc_vobc_neumann (v, sz, g, face)
    implicit none
    include '../FB/cbc_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, face, g
    integer, dimension(3)                                     ::  sz
    real                                                      ::  uu, vv, ww
    real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  v

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, g, face) &
!$OMP PRIVATE(i, j, k, uu, vv, ww)

    FACES : select case (face)
    case (X_minus)

#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif

      do k=1,kx
      do j=1,jx

        uu = v(1, 1, j, k)
        vv = v(2, 1, j, k)
        ww = v(3, 1, j, k)

        do i=1-g, 0
          v(1, i, j, k) = uu
          v(2, i, j, k) = vv
          v(3, i, j, k) = ww
        end do

      end do
      end do
!$OMP END DO
      

    case (X_plus)

#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif

      do k=1,kx
      do j=1,jx

        uu = v(1, ix, j, k)
        vv = v(2, ix, j, k)
        ww = v(3, ix, j, k)
        
        do i=ix+1, ix+g
          v(1, i, j, k) = uu
          v(2, i, j, k) = vv
          v(3, i, j, k) = ww
        end do

      end do
      end do
!$OMP END DO
      

    case (Y_minus)

#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif

      do k=1,kx
      do i=1,ix
        
        uu = v(1, i, 1, k)
        vv = v(2, i, 1, k)
        ww = v(3, i, 1, k)
                
        do j=1-g, 0
          v(1, i, j, k) = uu
          v(2, i, j, k) = vv
          v(3, i, j, k) = ww
        end do

      end do
      end do
!$OMP END DO
      

    case (Y_plus)

#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif

      do k=1,kx
      do i=1,ix

        uu = v(1, i, jx, k)
        vv = v(2, i, jx, k)
        ww = v(3, i, jx, k)
                
        do j=jx+1, jx+g
          v(1, i, j, k) = uu
          v(2, i, j, k) = vv
          v(3, i, j, k) = ww
        end do

      end do
      end do
!$OMP END DO
      

    case (Z_minus)

#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif

      do j=1,jx
      do i=1,ix

        uu = v(1, i, j, 1)
        vv = v(2, i, j, 1)
        ww = v(3, i, j, 1)
                
        do k=1-g, 0
          v(1, i, j, k) = uu
          v(2, i, j, k) = vv
          v(3, i, j, k) = ww
        end do

      end do
      end do
!$OMP END DO
      

    case (Z_plus)

#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif

      do j=1,jx
      do i=1,ix

        uu = v(1, i, j, kx)
        vv = v(2, i, j, kx)
        ww = v(3, i, j, kx)
                
        do k=kx+1, kx+g
          v(1, i, j, k) = uu
          v(2, i, j, k) = vv
          v(3, i, j, k) = ww
        end do

      end do
      end do
!$OMP END DO
      

    case default
    end select FACES

!$OMP END PARALLEL

    return 
    end subroutine cbc_vobc_neumann

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
    include '../FB/cbc_f_params.h'
    integer                                                     ::  i, j, k, g, idx, face, ix, jx, kx
    integer, dimension(3)                                       ::  sz
    real                                                        ::  Ue, Uw, Un, Us, Ut, Ub
    real                                                        ::  Ve, Vw, Vn, Vs, Vt, Vb
    real                                                        ::  We, Ww, Wn, Ws, Wt, Wb
    real                                                        ::  flop, c, rix, rjx, rkx
    real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  v, v0
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bv

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    rix = real(jx)*real(kx)
    rjx = real(ix)*real(kx)
    rkx = real(ix)*real(jx)
    
!$OMP PARALLEL REDUCTION(+:flop) &
!$OMP FIRSTPRIVATE(ix, jx, kx, face, c) &
!$OMP PRIVATE(i, j, k, idx) &
!$OMP PRIVATE(Ue, Uw, Un, Us, Ut, Ub) &
!$OMP PRIVATE(Ve, Vw, Vn, Vs, Vt, Vb) &
!$OMP PRIVATE(We, Ww, Wn, Ws, Wt, Wb)
    
    FACES : select case (face)
    case (X_minus)
      if ( c>0.0 ) c=0.0
      i = 1
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      
      do k=1,kx
      do j=1,jx
        idx = bv(i,j,k)
        if ( ibits(idx, bc_face_W, bitw_5) == obc_mask ) then
          Uw = v0(1, i-1,j  ,k  )
          Vw = v0(2, i-1,j  ,k  )
          Ww = v0(3, i-1,j  ,k  )
          Ue = v0(1, i  ,j  ,k  )
          Ve = v0(2, i  ,j  ,k  )
          We = v0(3, i  ,j  ,k  )

          v(1, i-1,j  ,k  ) = Uw - c*(Ue-Uw)
          v(2, i-1,j  ,k  ) = Vw - c*(Ve-Vw)
          v(3, i-1,j  ,k  ) = Ww - c*(We-Ww)
        endif
      end do
      end do
!$OMP END DO
      
      flop = flop + rix*9.0
      
      
    case (X_plus)
      if ( c<0.0 ) c=0.0
      i = ix
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      
      do k=1,kx
      do j=1,jx
        idx = bv(i,j,k)
        if ( ibits(idx, bc_face_E, bitw_5) == obc_mask ) then
          Uw = v0(1, i  ,j  ,k  )
          Vw = v0(2, i  ,j  ,k  )
          Ww = v0(3, i  ,j  ,k  )
          Ue = v0(1, i+1,j  ,k  )
          Ve = v0(2, i+1,j  ,k  )
          We = v0(3, i+1,j  ,k  )
          
          v(i+1,j  ,k  ,1) = Ue - c*(Ue-Uw)
          v(i+1,j  ,k  ,2) = Ve - c*(Ve-Vw)
          v(i+1,j  ,k  ,3) = We - c*(We-Ww)
        endif
      end do
      end do
!$OMP END DO
      
      flop = flop + rix*9.0
      
      
    case (Y_minus)
    if ( c>0.0 ) c=0.0
      j = 1
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      
      do k=1,kx
      do i=1,ix
        idx = bv(i,j,k)
        if ( ibits(idx, bc_face_S, bitw_5) == obc_mask ) then
          Us = v0(1, i  ,j-1,k  )
          Vs = v0(2, i  ,j-1,k  )
          Ws = v0(3, i  ,j-1,k  )
          Un = v0(1, i  ,j  ,k  )
          Vn = v0(2, i  ,j  ,k  )
          Wn = v0(3, i  ,j  ,k  )

          v(i  ,j-1,k  ,1) = Us - c*(Un-Us)
          v(i  ,j-1,k  ,2) = Vs - c*(Vn-Vs)
          v(i  ,j-1,k  ,3) = Ws - c*(Wn-Ws)
        endif
      end do
      end do
!$OMP END DO
      
      flop = flop + rjx*9.0
      
      
    case (Y_plus)
      if ( c<0.0 ) c=0.0
      j = jx
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      
      do k=1,kx
      do i=1,ix
        idx = bv(i,j,k)
        if ( ibits(idx, bc_face_N, bitw_5) == obc_mask ) then
          Us = v0(1, i  ,j  ,k  )
          Vs = v0(2, i  ,j  ,k  )
          Ws = v0(3, i  ,j  ,k  )
          Un = v0(1, i  ,j+1,k  )
          Vn = v0(2, i  ,j+1,k  )
          Wn = v0(3, i  ,j+1,k  )

          v(i  ,j+1,k  ,1) = Un - c*(Un-Us)
          v(i  ,j+1,k  ,2) = Vn - c*(Vn-Vs)
          v(i  ,j+1,k  ,3) = Wn - c*(Wn-Ws)
        endif
      end do
      end do
!$OMP END DO
      
      flop = flop + rjx*9.0
      
      
    case (Z_minus)
    if ( c>0.0 ) c=0.0
      k = 1
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      
      do j=1,jx
      do i=1,ix
        idx = bv(i,j,k)
        if ( ibits(idx, bc_face_B, bitw_5) == obc_mask ) then
          Ub = v0(1, i  ,j  ,k-1)
          Vb = v0(2, i  ,j  ,k-1)
          Wb = v0(3, i  ,j  ,k-1)
          Ut = v0(1, i  ,j  ,k  )
          Vt = v0(2, i  ,j  ,k  )
          Wt = v0(3, i  ,j  ,k  )

          v(i  ,j  ,k-1,1) = Ub - c*(Ut-Ub)
          v(i  ,j  ,k-1,2) = Vb - c*(Vt-Vb)
          v(i  ,j  ,k-1,3) = Wb - c*(Wt-Wb)
        endif
      end do
      end do
!$OMP END DO
      
      flop = flop + rkx*9.0
      
      
    case (Z_plus)
      if ( c<0.0 ) c=0.0
      k = kx
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      
      do j=1,jx
      do i=1,ix
        idx = bv(i,j,k)
        if ( ibits(idx, bc_face_T, bitw_5) == obc_mask ) then
          Ub = v0(1, i  ,j  ,k  )
          Vb = v0(2, i  ,j  ,k  )
          Wb = v0(3, i  ,j  ,k  )
          Ut = v0(1, i  ,j  ,k+1)
          Vt = v0(2, i  ,j  ,k+1)
          Wt = v0(3, i  ,j  ,k+1)

          v(i  ,j  ,k+1,1) = Ut - c*(Ut-Ub)
          v(i  ,j  ,k+1,2) = Vt - c*(Vt-Vb)
          v(i  ,j  ,k+1,3) = Wt - c*(Wt-Wb)
        endif
      end do
      end do
!$OMP END DO
      
      flop = flop + rkx*9.0
      
    case default
    end select FACES
    
!$OMP END PARALLEL
    
    return
    end subroutine cbc_vobc_outflow

!  *************************************************
!> @subroutine cbc_vobc_tfree (v, sz, g, face, flop)
!! @brief 速度の外部境界：　トラクションフリー
!! @param v 速度ベクトル
!! @param sz 配列長
!! @param g ガイドセル長
!! @param face 外部境界面の番号
!! @param flop 浮動小数演算数
!! @note トラクションフリー面は全て流体のこと
!<
    subroutine cbc_vobc_tfree (v, sz, g, face, flop)
    implicit none
    include '../FB/cbc_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, face, g, ii, jj, kk
    integer, dimension(3)                                     ::  sz
    real                                                      ::  v1, v2, v3, v4
    real                                                      ::  flop, rix, rjx, rkx
    real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  v

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    rix = real(jx)*real(kx)
    rjx = real(ix)*real(kx)
    rkx = real(ix)*real(jx)

!$OMP PARALLEL REDUCTION(+:flop) &
!$OMP FIRSTPRIVATE(ix, jx, kx, rix, rjx, rkx, g, face) &
!$OMP PRIVATE(i, j, k, ii, jj, kk, v1, v2, v3, v4)

    FACES : select case (face)
    case (X_minus)
      i = 1

#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do k=1,kx
      do j=1,jx

        v1 = 0.5 * (v(1, i-1, j+1, k  ) + v(1, i, j+1, k  ))
        v2 = 0.5 * (v(1, i-1, j-1, k  ) + v(1, i, j-1, k  ))
        v3 = 0.5 * (v(1, i-1, j  , k+1) + v(1, i, j  , k+1))
        v4 = 0.5 * (v(1, i-1, j  , k-1) + v(1, i, j  , k-1))

        v(1, i-1, j, k) = v(1, i, j, k)
        v(2, i-1, j, k) = v(2, i, j, k) + 0.5 * (v1 - v2)
        v(3, i-1, j, k) = v(3, i, j, k) + 0.5 * (v3 - v4)

      end do
      end do
!$OMP END DO
      
      flop = flop + rix*12.0
      

    case (X_plus)
      i = ix

#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do k=1,kx
      do j=1,jx

        v1 = 0.5 * (v(1, i+1, j+1, k  ) + v(1, i, j+1, k  ))
        v2 = 0.5 * (v(1, i+1, j-1, k  ) + v(1, i, j-1, k  ))
        v3 = 0.5 * (v(1, i+1, j  , k+1) + v(1, i, j  , k+1))
        v4 = 0.5 * (v(1, i+1, j  , k-1) + v(1, i, j  , k-1))
        
        v(1, i+1, j, k) = v(1, i, j, k)
        v(2, i+1, j, k) = v(2, i, j, k) - 0.5 * (v1 - v2)
        v(3, i+1, j, k) = v(3, i, j, k) - 0.5 * (v3 - v4)

      end do
      end do
!$OMP END DO
      
      flop = flop + rix*12.0
      

    case (Y_minus)
      j = 1

#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do k=1,kx
      do i=1,ix

        v1 = 0.5 * (v(2, i+1, j-1, k  ) + v(2, i+1, j, k  ))
        v2 = 0.5 * (v(2, i-1, j-1, k  ) + v(2, i-1, j, k  ))
        v3 = 0.5 * (v(2, i  , j-1, k+1) + v(2, i  , j, k+1))
        v4 = 0.5 * (v(2, i  , j-1, k-1) + v(2, i  , j, k-1))
                
        v(1, i, j-1, k) = v(1, i, j, k) + 0.5 * (v1 - v2)
        v(2, i, j-1, k) = v(2, i, j, k)
        v(3, i, j-1, k) = v(3, i, j, k) + 0.5 * (v3 - v4)

      end do
      end do
!$OMP END DO
      
      flop = flop + rjx*12.0
      

    case (Y_plus)
      j = jx

#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do k=1,kx
      do i=1,ix

        v1 = 0.5 * (v(2, i+1, j+1, k  ) + v(2, i+1, j, k  ))
        v2 = 0.5 * (v(2, i-1, j+1, k  ) + v(2, i-1, j, k  ))
        v3 = 0.5 * (v(2, i  , j+1, k+1) + v(2, i  , j, k+1))
        v4 = 0.5 * (v(2, i  , j+1, k-1) + v(2, i  , j, k-1))
                
        v(1, i, j+1, k) = v(1, i, j, k) - 0.5 * (v1 - v2)
        v(2, i, j+1, k) = v(2, i, j, k)
        v(3, i, j+1, k) = v(3, i, j, k) - 0.5 * (v3 - v4)

      end do
      end do
!$OMP END DO

      flop = flop + rjx*12.0
      

    case (Z_minus)
      k=1
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do j=1,jx
      do i=1,ix

        v1 = 0.5 * (v(3, i+1, j  , k-1) + v(3, i+1, j  , k))
        v2 = 0.5 * (v(3, i-1, j  , k-1) + v(3, i-1, j  , k))
        v3 = 0.5 * (v(3, i  , j+1, k-1) + v(3, i  , j+1, k))
        v4 = 0.5 * (v(3, i  , j-1, k-1) + v(3, i  , j-1, k))
                
        v(1, i, j, k-1) = v(1, i, j, k) + 0.5 * (v1 - v2)
        v(2, i, j, k-1) = v(2, i, j, k) + 0.5 * (v3 - v4)
        v(3, i, j, k-1) = v(3, i, j, k)

      end do
      end do
!$OMP END DO

      flop = flop + rkx*12.0
      

    case (Z_plus)
      k = kx

#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do j=1,jx
      do i=1,ix

        v1= 0.5 * (v(3, i+1, j  , k+1) + v(3, i+1, j  , k))
        v2= 0.5 * (v(3, i-1, j  , k+1) + v(3, i-1, j  , k))
        v3= 0.5 * (v(3, i  , j+1, k+1) + v(3, i  , j+1, k))
        v4= 0.5 * (v(3, i  , j-1, k+1) + v(3, i  , j-1, k))
                
        v(1, i, j, k+1) = v(1, i, j, k) - 0.5 * (v1 - v2)
        v(2, i, j, k+1) = v(2, i, j, k) - 0.5 * (v3 - v4)
        v(3, i, j, k+1) = v(3, i, j, k)

      end do
      end do
!$OMP END DO

      flop = flop + rkx*12.0
      

    case default
    end select FACES

!$OMP END PARALLEL

    return 
    end subroutine cbc_vobc_tfree

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
    include '../FB/cbc_f_params.h'
    integer                                                     ::  i, j, k, g, ix, jx, kx, face
    integer, dimension(3)                                       ::  sz
    real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  v, vc

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    
!$OMP PARALLEL &
!$OMP PRIVATE(i, j, k) &
!$OMP FIRSTPRIVATE(ix, jx, kx, face)
    
    FACES : select case (face)
    case (X_minus)
      i = 0
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do k=1,kx
      do j=1,jx
        v(1,i,j,k) = vc(1,i,j,k)
        v(2,i,j,k) = vc(2,i,j,k)
        v(3,i,j,k) = vc(3,i,j,k)
      end do
      end do
!$OMP END DO
      
    case (X_plus)
      i = ix+1
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do k=1,kx
      do j=1,jx
        v(1,i,j,k) = vc(1,i,j,k)
        v(2,i,j,k) = vc(2,i,j,k)
        v(3,i,j,k) = vc(3,i,j,k)
      end do
      end do
!$OMP END DO
      
    case (Y_minus)
      j = 0
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do k=1,kx
      do i=1,ix
        v(1,i,j,k) = vc(1,i,j,k)
        v(2,i,j,k) = vc(2,i,j,k)
        v(3,i,j,k) = vc(3,i,j,k)
      end do
      end do
!$OMP END DO
      
    case (Y_plus)
      j = jx+1
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do k=1,kx
      do i=1,ix
        v(1,i,j,k) = vc(1,i,j,k)
        v(2,i,j,k) = vc(2,i,j,k)
        v(3,i,j,k) = vc(3,i,j,k)
      end do
      end do
!$OMP END DO
      
    case (Z_minus)
      k = 0
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do j=1,jx
      do i=1,ix
        v(1,i,j,k) = vc(1,i,j,k)
        v(2,i,j,k) = vc(2,i,j,k)
        v(3,i,j,k) = vc(3,i,j,k)
      end do
      end do
!$OMP END DO
      
    case (Z_plus)
      k = kx+1
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do j=1,jx
      do i=1,ix
        v(1,i,j,k) = vc(1,i,j,k)
        v(2,i,j,k) = vc(2,i,j,k)
        v(3,i,j,k) = vc(3,i,j,k)
      end do
      end do
!$OMP END DO
      
    case default
    end select FACES

!$OMP END PARALLEL
    
    return
    end subroutine cbc_vobc_update


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
!! @param[out] flop flop count 近似
!<
    subroutine cbc_div_ibc_drchlt (div, sz, g, st, ed, v00, coef, bv, odr, vec, flop)
    implicit none
    include '../FB/cbc_f_params.h'
    integer                                                     ::  i, j, k, g, bvx, odr, is, ie, js, je, ks, ke
    integer, dimension(3)                                       ::  sz, st, ed
    real                                                        ::  flop, coef
    real                                                        ::  Ue_t, Uw_t, Vn_t, Vs_t, Wt_t, Wb_t
    real                                                        ::  u_bc_ref, v_bc_ref, w_bc_ref, m
    real, dimension(3)                                          ::  vec
    real, dimension(0:3)                                        ::  v00
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  div
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bv
    
    ! u_bc_refは参照座標系での境界速度
    u_bc_ref = vec(1) + v00(1)
    v_bc_ref = vec(2) + v00(2)
    w_bc_ref = vec(3) + v00(3)

    is = st(1)
    ie = ed(1)
    js = st(2)
    je = ed(2)
    ks = st(3)
    ke = ed(3)
    
    m = real( (ie-is+1)*(je-js+1)*(ke-ks+1) )
    flop = flop + m*7.0 + 3.0

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(is, ie, js, je, ks, ke, u_bc_ref, v_bc_ref, w_bc_ref, odr, coef) &
!$OMP PRIVATE(bvx, Ue_t, Uw_t, Vn_t, Vs_t, Wt_t, Wb_t)

#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
    do k=ks,ke
    do j=js,je
    do i=is,ie
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

        ! VBCの面だけUe_tなどは値をもつ  対象セルは流体なのでマスク不要
        div(i,j,k) = div(i,j,k) + ( Ue_t - Uw_t + Vn_t - Vs_t + Wt_t - Wb_t ) * coef
      end if
    end do
    end do
    end do
    
!$OMP END DO
!$OMP END PARALLEL

    return
    end subroutine cbc_div_ibc_drchlt

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
!! @note flop countはコスト軽減のため近似
!<
    subroutine cbc_div_ibc_oflow_pvec (div, sz, g, st, ed, v00, cf, coef, bv, odr, v0, flop)
    implicit none
    include '../FB/cbc_f_params.h'
    integer                                                     ::  i, j, k, g, bvx, odr, is, ie, js, je, ks, ke
    integer, dimension(3)                                       ::  sz, st, ed
    real                                                        ::  flop, coef, m
    real                                                        ::  b_w, b_e, b_s, b_n, b_b, b_t
    real                                                        ::  Ue, Uw, Vn, Vs, Wt, Wb
    real                                                        ::  Ue_t, Uw_t, Vn_t, Vs_t, Wt_t, Wb_t
    real                                                        ::  cf, u_ref, v_ref, w_ref
    real                                                        ::  Up0, Ue0, Uw0, Vp0, Vs0, Vn0, Wp0, Wb0, Wt0
    real, dimension(0:3)                                        ::  v00
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  div
    real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  v0
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bv

    ! 参照速度
    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)

    is = st(1)
    ie = ed(1)
    js = st(2)
    je = ed(2)
    ks = st(3)
    ke = ed(3)

    m = 0.0
    
!$OMP PARALLEL REDUCTION(+:m) &
!$OMP FIRSTPRIVATE(is, ie, js, je, ks, ke, u_ref, v_ref, w_ref, odr, coef, cf) &
!$OMP PRIVATE(bvx) &
!$OMP PRIVATE(b_w, b_e, b_s, b_n, b_b, b_t) &
!$OMP PRIVATE(Ue, Uw, Vn, Vs, Wt, Wb) &
!$OMP PRIVATE(Ue_t, Uw_t, Vn_t, Vs_t, Wt_t, Wb_t) &
!$OMP PRIVATE(Up0, Ue0, Uw0, Vp0, Vs0, Vn0, Wp0, Wb0, Wt0)

#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif

    do k=ks,ke
    do j=js,je
    do i=is,ie
      bvx = bv(i,j,k)
      if ( 0 /= iand(bvx, bc_mask30) ) then ! 6面のうちのどれか速度境界フラグが立っている場合

        include 'd_o_o_p.h' ! 42 flops
        
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
        endif
        
        if ( ibits(bvx, bc_face_E, bitw_5) == odr ) then
          Ue = Uw - (Vn - Vs + Wt - Wb)
          if ( cf<0.0 ) cf=0.0
          Ue_t = Ue - cf*(Ue-Uw)
        endif
        
        ! Y方向 ---------------------------------------
        if ( ibits(bvx, bc_face_S, bitw_5) == odr ) then
          Vs = Vn + (Ue - Uw + Wt - Wb)
          if ( cf>0.0 ) cf=0.0
          Vs_t = Vs - cf*(Vn-Vs)
        endif
        
        if ( ibits(bvx, bc_face_N, bitw_5) == odr ) then
          Vn = Vs - (Ue - Uw + Wt - Wb)
          if ( cf<0.0 ) cf=0.0
          Vn_t = Vn - cf*(Vn-Vs)
        endif
        
        ! Z方向 ---------------------------------------
        if ( ibits(bvx, bc_face_B, bitw_5) == odr ) then
          Wb = Wt + (Ue - Uw + Vn - Vs)
          if ( cf>0.0 ) cf=0.0
          Wb_t = Wb - cf*(Wt-Wb)
        endif
        
        if ( ibits(bvx, bc_face_T, bitw_5) == odr ) then
          Wt = Wb - (Ue - Uw + Vn - Vs)
          if ( cf<0.0 ) cf=0.0
          Wt_t = Wt - cf*(Wt-Wb)
        endif

        ! VBCの面だけUe_tなどは値をもつ
        div(i,j,k) = div(i,j,k) + ( Ue_t - Uw_t + Vn_t - Vs_t + Wt_t - Wb_t ) * coef ! 対象セルは流体なのでマスク不要
        m = m + 1.0
      end if
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    flop = flop + m*91.0

    return
    end subroutine cbc_div_ibc_oflow_pvec
    
!  ************************************************************************************
!> @subroutine cbc_div_ibc_oflow_vec (div, sz, g, st, ed, v00, coef, bv, odr, av, flop)
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
!! @param[out] av 積算流量と面積
!! @param[out] flop flop count
!! @note div(u)=0から，内部流出境界のセルで計算されたdivの値が流出速度となる
!! @note flop countはコスト軽減のため近似
!<
    subroutine cbc_div_ibc_oflow_vec (div, sz, g, st, ed, v00, coef, bv, odr, av, flop)
    implicit none
    include '../FB/cbc_f_params.h'
    integer                                                     ::  i, j, k, g, bvx, odr, is, ie, js, je, ks, ke
    integer, dimension(3)                                       ::  sz, st, ed
    real                                                        ::  flop, coef, rc
    real                                                        ::  u_ref, v_ref, w_ref, dv, a1, m
    real, dimension(0:3)                                        ::  v00
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  div
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bv
    real, dimension(2)                                          ::  av
    
    ! 参照速度
    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)
    a1 = 0.0
    m = 0.0
    rc = 1.0/coef

    is = st(1)
    ie = ed(1)
    js = st(2)
    je = ed(2)
    ks = st(3)
    ke = ed(3)

    flop = flop + 8.0 ! DP 13 flop

!$OMP PARALLEL &
!$OMP REDUCTION(+:a1) &
!$OMP REDUCTION(+:m) &
!$OMP FIRSTPRIVATE(is, ie, js, je, ks, ke, u_ref, v_ref, w_ref, rc, odr) &
!$OMP PRIVATE(bvx, dv)

#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif

    do k=ks,ke
    do j=js,je
    do i=is,ie
      bvx = bv(i,j,k)
      if ( 0 /= iand(bvx, bc_mask30) ) then ! 6面のうちのどれか速度境界フラグが立っている場合
        dv = div(i,j,k) * rc
        
        if ( ibits(bvx, bc_face_W, bitw_5) == odr ) then ! u_w
          a1 = a1 + dv - u_ref
        endif
        
        if ( ibits(bvx, bc_face_E, bitw_5) == odr ) then ! u_e
          a1 = a1 - dv - u_ref
        endif
        
        if ( ibits(bvx, bc_face_S, bitw_5) == odr ) then
          a1 = a1 + dv - v_ref
        endif
        
        if ( ibits(bvx, bc_face_N, bitw_5) == odr ) then
          a1 = a1 - dv - v_ref
        endif
        
        if ( ibits(bvx, bc_face_B, bitw_5) == odr ) then
          a1 = a1 + dv - w_ref
        endif
        
        if ( ibits(bvx, bc_face_T, bitw_5) == odr ) then
          a1 = a1 - dv - w_ref
        endif

        div(i,j,k) = 0.0 ! 対象セルは発散をゼロにする
        m = m + 1.0
      end if
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    av(1) = a1
    av(2) = m
    flop = flop + m*13.0

    return
    end subroutine cbc_div_ibc_oflow_vec

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
    include '../FB/cbc_f_params.h'
    integer                                                     ::  i, j, k, g, ix, jx, kx, face, bvx
    integer, dimension(3)                                       ::  sz
    real                                                        ::  flop, coef, rix, rjx, rkx
    real                                                        ::  u_bc_ref, v_bc_ref, w_bc_ref
    real, dimension(0:3)                                        ::  v00
    real, dimension(3)                                          ::  vec
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  div
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bv

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)

    rix = real(jx)*real(kx)
    rjx = real(ix)*real(kx)
    rkx = real(ix)*real(jx)
    
    ! 参照座標系の速度に係数をかけておく
    u_bc_ref = (vec(1) + v00(1)) * coef
    v_bc_ref = (vec(2) + v00(2)) * coef
    w_bc_ref = (vec(3) + v00(3)) * coef
    
    flop = flop + 6.0

!$OMP PARALLEL REDUCTION(+:flop) &
!$OMP FIRSTPRIVATE(ix, jx, kx, u_bc_ref, v_bc_ref, w_bc_ref, face) &
!$OMP PRIVATE(i, j, k, bvx)

    FACES : select case (face)
    case (X_minus)
      i = 1
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do k=1,kx
      do j=1,jx
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_W, bitw_5) == obc_mask ) then
          div(i,j,k) = div(i,j,k) - u_bc_ref * real(ibits(bvx, State, 1))
        endif
      end do
      end do
!$OMP END DO

      flop = flop + rix*3.0 ! 2+ real*1
      
      
    case (X_plus)
      i = ix
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do k=1,kx
      do j=1,jx
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_E, bitw_5) == obc_mask ) then
          div(i,j,k) = div(i,j,k) + u_bc_ref * real(ibits(bvx, State, 1))
        endif
      end do
      end do
!$OMP END DO
      
      flop = flop + rix*3.0 ! 2+ real*1
      
      
    case (Y_minus)
      j = 1
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do k=1,kx
      do i=1,ix
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_S, bitw_5) == obc_mask ) then
          div(i,j,k) = div(i,j,k) - v_bc_ref * real(ibits(bvx, State, 1))
        endif
      end do
      end do
!$OMP END DO
      
      flop = flop + rjx*3.0 ! 2+ real*1
      
      
    case (Y_plus)
      j = jx
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do k=1,kx
      do i=1,ix
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_N, bitw_5) == obc_mask ) then
          div(i,j,k) = div(i,j,k) + v_bc_ref * real(ibits(bvx, State, 1))
        endif
      end do
      end do
!$OMP END DO
      
      flop = flop + rjx*3.0 ! 2+ real*1
    
    
    case (Z_minus)
      k = 1
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do j=1,jx
      do i=1,ix
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_B, bitw_5) == obc_mask ) then
          div(i,j,k) = div(i,j,k) - w_bc_ref * real(ibits(bvx, State, 1))
        endif
      end do
      end do
!$OMP END DO
      
      flop = flop + rkx*3.0 ! 2+ real*1
      
      
    case (Z_plus)
      k = kx
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do j=1,jx
      do i=1,ix
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_T, bitw_5) == obc_mask ) then
          div(i,j,k) = div(i,j,k) + w_bc_ref * real(ibits(bvx, State, 1))
        endif
      end do
      end do
!$OMP END DO
      
      flop = flop + rkx*3.0 ! 2+ real*1
    
    case default
    end select FACES

!$OMP END PARALLEL

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
    include '../FB/cbc_f_params.h'
    integer                                                     ::  i, j, k, g, ix, jx, kx, face, bvx
    integer, dimension(3)                                       ::  sz
    real                                                        ::  flop, coef, v_out
    real                                                        ::  b_w, b_e, b_s, b_n, b_b, b_t
    real                                                        ::  Up0, Ue0, Uw0, Vp0, Vs0, Vn0, Wp0, Wb0, Wt0
    real                                                        ::  Ue, Uw, Vn, Vs, Wt, Wb
    real                                                        ::  Ue_t, Uw_t, Vn_t, Vs_t, Wt_t, Wb_t
    real                                                        ::  u_ref, v_ref, w_ref, m
    real, dimension(0:3)                                        ::  v00
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  div
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bv
    real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  v0

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    
    ! 参照速度
    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)

    m = 0.0

!$OMP PARALLEL REDUCTION(+:m) &
!$OMP FIRSTPRIVATE(ix, jx, kx, u_ref, v_ref, w_ref, face, coef, v_out) &
!$OMP PRIVATE(i, j, k, bvx) &
!$OMP PRIVATE(Up0, Ue0, Uw0, Vp0, Vs0, Vn0, Wp0, Wb0, Wt0) &
!$OMP PRIVATE(Ue, Uw, Vn, Vs, Wt, Wb) &
!$OMP PRIVATE(Ue_t, Uw_t, Vn_t, Vs_t, Wt_t, Wb_t) &
!$OMP PRIVATE(b_w, b_e, b_s, b_n, b_b, b_t)


    FACES : select case (face)
    
    case (X_minus)
      if ( v_out>0.0 ) v_out=0.0
      i = 1
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do k=1,kx
      do j=1,jx
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_W, bitw_5) == obc_mask ) then
        
          include 'd_o_o_p.h' ! 42 flops
          
          Uw = Ue + (Vn - Vs + Wt - Wb) ! 連続の式から流出面の速度を推定，これは移動座標系上の速度成分
          Uw_t = Uw - v_out*(Ue-Uw)
          div(i,j,k) = div(i,j,k) - Uw_t * coef
          m = m + 1.0
        endif
      end do
      end do
!$OMP END DO
      
      
    case (X_plus)
      if ( v_out<0.0 ) v_out=0.0
      i = ix
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do k=1,kx
      do j=1,jx
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_E, bitw_5) == obc_mask ) then
        
          include 'd_o_o_p.h'
          
          Ue = Uw - (Vn - Vs + Wt - Wb)
          Ue_t = Ue - v_out*(Ue-Uw)
          div(i,j,k) = div(i,j,k) + Ue_t * coef
          m = m + 1.0
        endif
      end do
      end do
!$OMP END DO


    case (Y_minus)
      if ( v_out>0.0 ) v_out=0.0
      j = 1
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do k=1,kx
      do i=1,ix
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_S, bitw_5) == obc_mask ) then
        
          include 'd_o_o_p.h'
        
          Vs = Vn + (Ue - Uw + Wt - Wb)
          Vs_t = Vs - v_out*(Vn-Vs)
          div(i,j,k) = div(i,j,k) - Vs_t * coef
          m = m + 1.0
        endif
      end do
      end do
!$OMP END DO
      
      
    case (Y_plus)
      if ( v_out<0.0 ) v_out=0.0
      j = jx
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do k=1,kx
      do i=1,ix
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_N, bitw_5) == obc_mask ) then
        
          include 'd_o_o_p.h'
        
          Vn = Vs - (Ue - Uw + Wt - Wb)
          Vn_t = Vn - v_out*(Vn-Vs)
          div(i,j,k) = div(i,j,k) + Vn_t * coef
          m = m + 1.0
        endif
      end do
      end do
!$OMP END DO
    
    
    case (Z_minus)
      if ( v_out>0.0 ) v_out=0.0
      k = 1
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do j=1,jx
      do i=1,ix
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_B, bitw_5) == obc_mask ) then
        
          include 'd_o_o_p.h'
        
          Wb = Wt + (Ue - Uw + Vn - Vs)
          Wb_t = Wb - v_out*(Wt-Wb)
          div(i,j,k) = div(i,j,k) - Wb_t * coef
          m = m + 1.0
        endif
      end do
      end do
!$OMP END DO
      
      
    case (Z_plus)
      if ( v_out<0.0 ) v_out=0.0
      k = kx
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do j=1,jx
      do i=1,ix
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_T, bitw_5) == obc_mask ) then
        
          include 'd_o_o_p.h'
        
          Wt = Wb - (Ue - Uw + Vn - Vs)
          Wt_t = Wt - v_out*(Wt-Wb)
          div(i,j,k) = div(i,j,k) + Wt_t * coef
          m = m + 1.0
        endif
      end do
      end do
!$OMP END DO
    
    case default
    end select FACES

!$OMP END PARALLEL
    
    flop = flop + m*51.0

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
!! @param[out] flop flop count 近似
!! @note 指定面でも固体部分は対象外とするのでループ中に判定あり
!!       div(u)=0から，内部流出境界のセルで計算されたdivが流出速度となる
!<
    subroutine cbc_div_obc_oflow_vec (div, sz, g, face, v00, coef, bv, aa, flop)
    implicit none
    include '../FB/cbc_f_params.h'
    integer                                                     ::  i, j, k, g, ix, jx, kx, face, bvx
    integer, dimension(3)                                       ::  sz
    real                                                        ::  flop, coef, dv, a1, a2, a3, u_ref, v_ref, w_ref, rc
    real                                                        ::  rix, rjx, rkx
    real, dimension(0:3)                                        ::  v00
    real, dimension(3)                                          ::  aa
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  div
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bv

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    rix = real(jx)*real(kx)
    rjx = real(ix)*real(kx)
    rkx = real(ix)*real(jx)

    a1 = 0.0   ! sum
    a2 = 1.0e6 ! min
    a3 =-1.0e6 ! max
    
    rc = 1.0/coef

    flop = flop + 8.0 ! DP 13 flop
    
    ! 参照速度
    u_ref = v00(1)
    v_ref = v00(2)
    w_ref = v00(3)
    
!$OMP PARALLEL &
!$OMP REDUCTION(+:a1) &
!$OMP REDUCTION(min:a2) &
!$OMP REDUCTION(max:a3) &
!$OMP REDUCTION(+:flop) &
!$OMP FIRSTPRIVATE(ix, jx, kx, u_ref, v_ref, w_ref, rix, rjx, rkx, rc, face) &
!$OMP PRIVATE(i, j, k, bvx, dv)

    FACES : select case (face)
    
    case (X_minus)
      i = 1
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do k=1,kx
      do j=1,jx
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_W, bitw_5) == obc_mask ) then
          dv = div(i,j,k) * rc - u_ref
          a1 = a1 + dv
          a2 = min(a2, dv)
          a3 = max(a3, dv)
          div(i,j,k) = 0.0 ! 対象セルは発散をゼロにする
        endif
      end do
      end do
!$OMP END DO
      
      flop = flop + rix*5.0
      
      
    case (X_plus)
      i = ix
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do k=1,kx
      do j=1,jx
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_E, bitw_5) == obc_mask ) then
          dv = -div(i,j,k) * rc  - u_ref
          a1 = a1 + dv
          a2 = min(a2, dv)
          a3 = max(a3, dv)
          div(i,j,k) = 0.0
        endif
      end do
      end do
!$OMP END DO

      flop = flop + rix*5.0
      
      
    case (Y_minus)
      j = 1
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do k=1,kx
      do i=1,ix
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_S, bitw_5) == obc_mask ) then
          dv = div(i,j,k) * rc - v_ref
          a1 = a1 + dv
          a2 = min(a2, dv)
          a3 = max(a3, dv)
          div(i,j,k) = 0.0
        endif
      end do
      end do
!$OMP END DO

      flop = flop + rjx*5.0
      
      
    case (Y_plus)
      j = jx
      
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do k=1,kx
      do i=1,ix
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_N, bitw_5) == obc_mask ) then
          dv = -div(i,j,k) * rc - v_ref
          a1 = a1 + dv
          a2 = min(a2, dv)
          a3 = max(a3, dv)
          div(i,j,k) = 0.0
        endif
      end do
      end do
!$OMP END DO

      flop = flop + rjx*5.0
    
    
    case (Z_minus)
      k = 1
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do j=1,jx
      do i=1,ix
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_B, bitw_5) == obc_mask ) then
          dv = div(i,j,k) * rc - w_ref
          a1 = a1 + dv
          a2 = min(a2, dv)
          a3 = max(a3, dv)
          div(i,j,k) = 0.0
        endif
      end do
      end do
!$OMP END DO

      flop = flop + rkx*5.0
      
      
    case (Z_plus)
      k = kx
#ifdef _DYNAMIC
!$OMP DO SCHEDULE(dynamic,1)
#elif defined _STATIC
!$OMP DO SCHEDULE(static)
#else
!$OMP DO SCHEDULE(hoge)
#endif
      do j=1,jx
      do i=1,ix
        bvx = bv(i,j,k)
        if ( ibits(bvx, bc_face_T, bitw_5) == obc_mask ) then
          dv = -div(i,j,k) * rc - w_ref
          a1 = a1 + dv
          a2 = min(a2, dv)
          a3 = max(a3, dv)
          div(i,j,k) = 0.0
        endif
      end do
      end do
!$OMP END DO

      flop = flop + rkx*5.0
    
    case default
    end select FACES

!$OMP END PARALLEL
    
    aa(1) = a1
    aa(2) = a2
    aa(3) = a3

    return
    end subroutine cbc_div_obc_oflow_vec
    