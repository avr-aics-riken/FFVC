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
! Copyright (c) 2012-2015 Advanced Institute for Computational Science, RIKEN.
! All rights reserved.
!
!###################################################################################
    
!> @file   ffv_pscalar.f90
!! @brief  パッシブスカラー計算のルーチン群
!! @author aics
!<

!> ********************************************************************
!! @brief パッシブスカラの移流部分の積分
!! @param [in,out] ws     内部エネルギーの対流項による増分
!! @param [in]     sz     配列長
!! @param [in]     g      ガイドセル長
!! @param [in]     dh     格子幅
!! @param [in]     scheme 対流項スキームのモード（1-UWD, 3-MUSCL）
!! @param [in]     v00    参照速度
!! @param [in]     vf     セルフェイス速度ベクトル（n+1 step）
!! @param [in]     ie     内部エネルギー
!! @param [in]     bid    Cut ID
!! @param [in]     cdf    BCindex C
!! @param [in]     bcd    BCindex B
!! @param [in]     swt    固定壁の扱い（0-断熱，1-共役熱移動）
!! @param [out]    flop   浮動小数点演算数
!<
subroutine ps_muscl (ws, sz, g, dh, scheme, v00, vf, ie, bid, cdf, bcd, swt, flop)
implicit none
include '../FB/ffv_f_params.h'
integer                                                     ::  i, j, k, ix, jx, kx, g, scheme, idx, swt, hdx, bix
integer                                                     ::  b_e1, b_w1, b_n1, b_s1, b_t1, b_b1
integer                                                     ::  b_e2, b_w2, b_n2, b_s2, b_t2, b_b2, b_p
integer, dimension(3)                                       ::  sz
double precision                                            ::  flop
real                                                        ::  UPe, UPw, VPn, VPs, WPt, WPb
real                                                        ::  Fp0, Fe1, Fe2, Fw1, Fw2, Fs1, Fs2, Fn1, Fn2, Fb1, Fb2, Ft1, Ft2
real                                                        ::  ck, u_ref, v_ref, w_ref, actv, rx, ry, rz
real                                                        ::  c_e1, c_w1, c_n1, c_s1, c_t1, c_b1
real                                                        ::  c_e2, c_w2, c_n2, c_s2, c_t2, c_b2
real                                                        ::  a_e, a_w, a_n, a_s, a_t, a_b
real                                                        ::  dv1, dv2, dv3, dv4, g1, g2, g3, g4, g5, g6, s1, s2, s3, s4
real                                                        ::  Fr_r, Fr_l, Fl_r, Fl_l
real                                                        ::  cr, cl, acr, acl, cnv, ss, b, cm1, cm2, ss_4
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)   ::  vf
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  ie, ws
real, dimension(0:3)                                        ::  v00
real, dimension(3)                                          ::  dh
integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  cdf, bcd, bid

ix = sz(1)
jx = sz(2)
kx = sz(3)

rx = 1.0/dh(1)
ry = 1.0/dh(2)
rz = 1.0/dh(3)

ss = 1.0
ck = 0.0
b = 0.0

! 参照座標系の速度
u_ref = v00(1)
v_ref = v00(2)
w_ref = v00(3)

if ( scheme == 1 ) then      !     1st order upwind
ss = 0.0
else if ( scheme == 3 ) then !     3rd order MUSCL
ck = 1.0/3.0
b  = (3.0-ck)/(1.0-ck)
endif

ss_4 = 0.25*ss

cm1 = 1.0 - ck
cm2 = 1.0 + ck

! 24 + 3*84 + 3
flop = flop + dble(ix)*dble(jx)*dble(kx)*279.0d0 + 45.0d0


!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, rx, ry, rz, ss, ck, b, u_ref, v_ref, w_ref, swt, cm1, cm2, ss_4) &
!$OMP PRIVATE(cnv, idx, hdx, bix, actv) &
!$OMP PRIVATE(b_e1, b_w1, b_n1, b_s1, b_t1, b_b1, b_e2, b_w2, b_n2, b_s2, b_t2, b_b2, b_p) &
!$OMP PRIVATE(Fp0, Fe1, Fe2, Fw1, Fw2, Fs1, Fs2, Fn1, Fn2, Fb1, Fb2, Ft1, Ft2) &
!$OMP PRIVATE(c_e1, c_w1, c_n1, c_s1, c_t1, c_b1) &
!$OMP PRIVATE(c_e2, c_w2, c_n2, c_s2, c_t2, c_b2) &
!$OMP PRIVATE(a_e, a_w, a_n, a_s, a_t, a_b) &
!$OMP PRIVATE(UPe, UPw, VPn, VPs, WPt, WPb) &
!$OMP PRIVATE(dv1, dv2, dv3, dv4, g1, g2, g3, g4, g5, g6, s1, s2, s3, s4) &
!$OMP PRIVATE(Fr_r, Fr_l, Fl_r, Fl_l) &
!$OMP PRIVATE(cr, cl, acr, acl)

!$OMP DO SCHEDULE(static)

do k=1,kx
do j=1,jx
do i=1,ix

cnv = 0.0

! 変数のロード
Fp0 = ie(i  ,j  ,k  )
Fw2 = ie(i-2,j  ,k  )
Fw1 = ie(i-1,j  ,k  )
Fe1 = ie(i+1,j  ,k  )
Fe2 = ie(i+2,j  ,k  )
Fs2 = ie(i  ,j-2,k  )
Fs1 = ie(i  ,j-1,k  )
Fn1 = ie(i  ,j+1,k  )
Fn2 = ie(i  ,j+2,k  )
Fb2 = ie(i  ,j  ,k-2)
Fb1 = ie(i  ,j  ,k-1)
Ft1 = ie(i  ,j  ,k+1)
Ft2 = ie(i  ,j  ,k+2)

idx = cdf(i,j,k)
hdx = bcd(i,j,k)
bix = bid(i,j,k)
actv= real(ibits(hdx, State, 1)) ! 対流マスクなのでbcdの状態を参照

! 各方向に物体があれば，マスクはゼロ
! セル界面のフラグ b_?1={0.0-wall face / 1.0-fluid} <<= bix={0-fluid, ID-wall}
b_w1 = 1.0
b_e1 = 1.0
b_s1 = 1.0
b_n1 = 1.0
b_b1 = 1.0
b_t1 = 1.0
if ( ibits(bix, bc_face_W, bitw_5) /= 0 ) b_w1 = 0.0
if ( ibits(bix, bc_face_E, bitw_5) /= 0 ) b_e1 = 0.0
if ( ibits(bix, bc_face_S, bitw_5) /= 0 ) b_s1 = 0.0
if ( ibits(bix, bc_face_N, bitw_5) /= 0 ) b_n1 = 0.0
if ( ibits(bix, bc_face_B, bitw_5) /= 0 ) b_b1 = 0.0
if ( ibits(bix, bc_face_T, bitw_5) /= 0 ) b_t1 = 0.0


! (i,j,k)を基準にした遠い方向なので，隣接セルで判断
b_w2 = 1.0
b_e2 = 1.0
b_s2 = 1.0
b_n2 = 1.0
b_b2 = 1.0
b_t2 = 1.0
if ( ibits(bid(i-1,j  ,k  ), bc_face_W, bitw_5) /= 0 ) b_w2 = 0.0
if ( ibits(bid(i+1,j  ,k  ), bc_face_E, bitw_5) /= 0 ) b_e2 = 0.0
if ( ibits(bid(i  ,j-1,k  ), bc_face_S, bitw_5) /= 0 ) b_s2 = 0.0
if ( ibits(bid(i  ,j+1,k  ), bc_face_N, bitw_5) /= 0 ) b_n2 = 0.0
if ( ibits(bid(i  ,j  ,k-1), bc_face_B, bitw_5) /= 0 ) b_b2 = 0.0
if ( ibits(bid(i  ,j  ,k+1), bc_face_T, bitw_5) /= 0 ) b_t2 = 0.0


! 各面のVBCフラグの有無 => flux mask
! ibits() = 1(Normal) / 0(VBC)
c_w1 = real( ibits(hdx, bc_d_W, 1) )
c_e1 = real( ibits(hdx, bc_d_E, 1) )
c_s1 = real( ibits(hdx, bc_d_S, 1) )
c_n1 = real( ibits(hdx, bc_d_N, 1) )
c_b1 = real( ibits(hdx, bc_d_B, 1) )
c_t1 = real( ibits(hdx, bc_d_T, 1) )


! ステンシルの参照先がvspec, outflowである場合のスキームの破綻を回避，１次精度におとすSW
c_w2 = real( ibits(bcd(i-1, j  , k  ), bc_d_W, 1) )
c_e2 = real( ibits(bcd(i+1, j  , k  ), bc_d_E, 1) )
c_s2 = real( ibits(bcd(i  , j-1, k  ), bc_d_S, 1) )
c_n2 = real( ibits(bcd(i  , j+1, k  ), bc_d_N, 1) )
c_b2 = real( ibits(bcd(i  , j  , k-1), bc_d_B, 1) )
c_t2 = real( ibits(bcd(i  , j  , k+1), bc_d_T, 1) )



! 断熱マスク
a_e = real(ibits(hdx, adbtc_E, 1))
a_w = real(ibits(hdx, adbtc_W, 1))
a_n = real(ibits(hdx, adbtc_N, 1))
a_s = real(ibits(hdx, adbtc_S, 1))
a_t = real(ibits(hdx, adbtc_T, 1))
a_b = real(ibits(hdx, adbtc_B, 1))



! 界面速度（スタガード位置） > 24 flops
UPe = vf(i  , j  , k  ,1)*b_e1 + u_ref*(1.0-b_e1)
UPw = vf(i-1, j  , k  ,1)*b_w1 + u_ref*(1.0-b_w1)
VPn = vf(i  , j  , k  ,2)*b_n1 + v_ref*(1.0-b_n1)
VPs = vf(i  , j-1, k  ,2)*b_s1 + v_ref*(1.0-b_s1)
WPt = vf(i  , j  , k  ,3)*b_t1 + w_ref*(1.0-b_t1)
WPb = vf(i  , j  , k-1,3)*b_b1 + w_ref*(1.0-b_b1)


! X方向 --------------------------------------- >> 4 + 36 + 22 + 4 + 18 = 84
! 壁面の扱い　共役熱移動の場合は，固体セル内部の値をそのまま使う
if (swt == 0) then ! (i,j,k)からみて断熱壁の処理
  if ( b_e2 == 0 ) Fe2 = Fe1
  if ( b_e1 == 0 ) Fe1 = Fp0
  if ( b_w1 == 0 ) Fw1 = Fp0
  if ( b_w2 == 0 ) Fw2 = Fw1
endif

dv4 = Fe2 - Fe1
dv3 = Fe1 - Fp0
dv2 = Fp0 - Fw1
dv1 = Fw1 - Fw2 ! 4 flop

s4 = sign(1.0, dv4) ! sign is zero flop
s3 = sign(1.0, dv3)
s2 = sign(1.0, dv2)
s1 = sign(1.0, dv1)

g6 = s4 * max(0.0, min( abs(dv4), s4 * b * dv3))
g5 = s3 * max(0.0, min( abs(dv3), s3 * b * dv4))
g4 = s3 * max(0.0, min( abs(dv3), s3 * b * dv2))
g3 = s2 * max(0.0, min( abs(dv2), s2 * b * dv3))
g2 = s2 * max(0.0, min( abs(dv2), s2 * b * dv1))
g1 = s1 * max(0.0, min( abs(dv1), s1 * b * dv2)) ! 36 flop

Fr_r = Fe1 - (cm1*g6+cm2*g5)*ss_4 * c_e2
Fr_l = Fp0 + (cm1*g3+cm2*g4)*ss_4
Fl_r = Fp0 - (cm1*g4+cm2*g3)*ss_4
Fl_l = Fw1 + (cm1*g1+cm2*g2)*ss_4 * c_w2 ! 22 flop

! 流束　壁面上で速度ゼロ->対流熱流束がゼロになる >  4 flop
cr  = UPe - u_ref
cl  = UPw - u_ref
acr = abs(cr)
acl = abs(cl)

! 境界条件の面の寄与はスキップ > 18 flop
cnv = (0.5*(cr*(Fr_r+Fr_l) - acr*(Fr_r-Fr_l)) * c_e1 * a_e &
    -  0.5*(cl*(Fl_r+Fl_l) - acl*(Fl_r-Fl_l)) * c_w1 * a_w ) * rx



! Y方向 ---------------------------------------
! 壁面の扱い
if (swt == 0) then ! (i,j,k)からみて断熱壁の処理
  if ( b_n2 == 0 ) Fn2 = Fn1
  if ( b_n1 == 0 ) Fn1 = Fp0
  if ( b_s1 == 0 ) Fs1 = Fp0
  if ( b_s2 == 0 ) Fs2 = Fs1
endif

dv4 = Fn2 - Fn1
dv3 = Fn1 - Fp0
dv2 = Fp0 - Fs1
dv1 = Fs1 - Fs2

s4 = sign(1.0, dv4) ! sign is zero flop
s3 = sign(1.0, dv3)
s2 = sign(1.0, dv2)
s1 = sign(1.0, dv1)

g6 = s4 * max(0.0, min( abs(dv4), s4 * b * dv3))
g5 = s3 * max(0.0, min( abs(dv3), s3 * b * dv4))
g4 = s3 * max(0.0, min( abs(dv3), s3 * b * dv2))
g3 = s2 * max(0.0, min( abs(dv2), s2 * b * dv3))
g2 = s2 * max(0.0, min( abs(dv2), s2 * b * dv1))
g1 = s1 * max(0.0, min( abs(dv1), s1 * b * dv2))

Fr_r = Fn1 - (cm1*g6+cm2*g5)*ss_4 * c_n2
Fr_l = Fp0 + (cm1*g3+cm2*g4)*ss_4
Fl_r = Fp0 - (cm1*g4+cm2*g3)*ss_4
Fl_l = Fs1 + (cm1*g1+cm2*g2)*ss_4 * c_s2

cr  = VPn - v_ref
cl  = VPs - v_ref
acr = abs(cr)
acl = abs(cl)

cnv = (0.5*(cr*(Fr_r+Fr_l) - acr*(Fr_r-Fr_l)) * c_n1 * a_n &
    -  0.5*(cl*(Fl_r+Fl_l) - acl*(Fl_r-Fl_l)) * c_s1 * a_s ) * ry + cnv



! Z方向 ---------------------------------------
! 壁面の扱い
if (swt == 0) then ! (i,j,k)からみて断熱壁の処理
  if ( b_t2 == 0 ) Ft2 = Ft1
  if ( b_t1 == 0 ) Ft1 = Fp0
  if ( b_b1 == 0 ) Fb1 = Fp0
  if ( b_b2 == 0 ) Fb2 = Fb1
endif

dv4 = Ft2 - Ft1
dv3 = Ft1 - Fp0
dv2 = Fp0 - Fb1
dv1 = Fb1 - Fb2

s4 = sign(1.0, dv4) ! sign is zero flop
s3 = sign(1.0, dv3)
s2 = sign(1.0, dv2)
s1 = sign(1.0, dv1)

g6 = s4 * max(0.0, min( abs(dv4), s4 * b * dv3))
g5 = s3 * max(0.0, min( abs(dv3), s3 * b * dv4))
g4 = s3 * max(0.0, min( abs(dv3), s3 * b * dv2))
g3 = s2 * max(0.0, min( abs(dv2), s2 * b * dv3))
g2 = s2 * max(0.0, min( abs(dv2), s2 * b * dv1))
g1 = s1 * max(0.0, min( abs(dv1), s1 * b * dv2))

Fr_r = Ft1 - (cm1*g6+cm2*g5)*ss_4 * c_t2
Fr_l = Fp0 + (cm1*g3+cm2*g4)*ss_4
Fl_r = Fp0 - (cm1*g4+cm2*g3)*ss_4
Fl_l = Fb1 + (cm1*g1+cm2*g2)*ss_4 * c_b2

cr  = WPt - w_ref
cl  = WPb - w_ref
acr = abs(cr)
acl = abs(cl)

cnv = (0.5*(cr*(Fr_r+Fr_l) - acr*(Fr_r-Fr_l)) * c_t1 * a_t &
    -  0.5*(cl*(Fl_r+Fl_l) - acl*(Fl_r-Fl_l)) * c_b1 * a_b ) * rz + cnv

! ---------------------------------------
ws(i,j,k) = -cnv * actv

end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL

return
end subroutine ps_muscl


!> ********************************************************************
!! @brief 浮力項の計算（単相流）
!! @param [in,out] v     速度ベクトル
!! @param [in]     sz    配列長
!! @param [in]     g     ガイドセル長
!! @param [in]     dgr   係数
!! @param [in]     ie    内部エネルギー
!! @param [in]     bd    BCindex B
!! @param [in]     ncompo コンポーネント数
!! @param [in]     mtbl   コンポーネントの物性値
!! @param [in,out] flop  浮動小数点演算数
!! @todo 対象セルは流体だけでよい？　active flag?
!<
    subroutine ps_buoyancy (v, sz, g, dgr, ie, bd, ncompo, mtbl, flop)
    implicit none
    include '../FB/ffv_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, g, ncompo, l, idx
    integer, dimension(3)                                     ::  sz
    double precision                                          ::  flop
    real                                                      ::  dgr, r, rcp, t
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  ie
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bd
    real*8, dimension(3, 0:ncompo)                            ::  mtbl

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    r = dgr
    
    flop = flop + dble(ix)*dble(jx)*dble(kx)*12.0d0

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, r) &
!$OMP PRIVATE(idx, l, rcp, t)

!$OMP DO SCHEDULE(static)
    do k=1,kx
    do j=1,jx
    do i=1,ix
      idx = bd(i,j,k)
      l = ibits(idx, 0, bitw_5)
      rcp = mtbl(1, l) * mtbl(2, l)
      t = ie(i, j, k) / rcp
      v(i,j,k,3) = v(i,j,k,3) + r * t * real(ibits(idx, State, 1))
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    return
    end subroutine ps_buoyancy



!> ********************************************************************
!! @brief 移流項のEuler陽解法による時間積分
!! @param [in,out] ie    内部エネルギー
!! @param [in]     sz    配列長
!! @param [in]     g     ガイドセル長
!! @param [in]     dt    時間積分幅
!! @param [in]     bd    BCindex B
!! @param [in]     ie0   内部エネルギー n-step
!! @param [in,out] flop  浮動小数点演算数
!<
subroutine ps_convection_ee(ie, sz, g, dt, bd, ie0, flop)
implicit none
include '../FB/ffv_f_params.h'
integer                                                   ::  i, j, k, ix, jx, kx, g
integer, dimension(3)                                     ::  sz
double precision                                          ::  flop
real                                                      ::  dt
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  ie, ie0
integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bd

ix = sz(1)
jx = sz(2)
kx = sz(3)

flop = flop + dble(ix*jx*kx)*3.0d0

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, dt)

!$OMP DO SCHEDULE(static)
do k=1,kx
do j=1,jx
do i=1,ix
  ie(i,j,k) = ie0(i,j,k) + dt * ie(i,j,k) * real(ibits(bd(i,j,k), State, 1))
end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL

return
end subroutine ps_convection_ee



!> ********************************************************************
!! @brief 温度の拡散項の半陰的時間積分（対流項を積分した結果を用いて粘性項を計算）
!! @param [in,out] ie     内部エネルギー
!! @param [in]     sz     配列長
!! @param [in]     g      ガイドセル長
!! @param [out]    res_l2 残差の自乗和
!! @param [in]     dh     格子幅
!! @param [in]     dt     時間積分幅
!! @param [in]     qbc    境界条件の熱流束
!! @param [in]     bcd    BCindex B
!! @param [in]     ws     部分段階の温度
!! @param [in]     ncompo コンポーネント数
!! @param [in]     mtbl   コンポーネントの物性値
!! @param [in]     h_mode mode(0-conjugate heat transfer / 1-othres)
!! @param [in,out] flop   浮動小数演算数
!<
subroutine ps_diff_ee (ie, sz, g, res_l2, dh, dt, qbc, bcd, ws, ncompo, mtbl, h_mode, flop)
implicit none
include '../FB/ffv_f_params.h'
integer                                                   ::  i, j, k, ix, jx, kx, g, idx, ncompo, h_mode, hm
integer                                                   ::  l_p, l_w, l_e, l_s, l_n, l_b, l_t
integer, dimension(3)                                     ::  sz
double precision                                          ::  flop, res, res_l2
real                                                      ::  dt, delta, sw, rx, ry, rz, rx2, ry2, rz2
real                                                      ::  t_p, t_w, t_e, t_s, t_n, t_b, t_t
real                                                      ::  g_w, g_e, g_s, g_n, g_b, g_t
real                                                      ::  a_w, a_e, a_s, a_n, a_b, a_t
real                                                      ::  lmd_p, lmd_w, lmd_e, lmd_s, lmd_n, lmd_b, lmd_t
real                                                      ::  rcp_p, rcp_w, rcp_e, rcp_s, rcp_n, rcp_b, rcp_t
real                                                      ::  tc_w, tc_e, tc_s, tc_n, tc_b, tc_t
real                                                      ::  c_w, c_e, c_s, c_n, c_b, c_t
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  ie, ws
real, dimension(6, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  qbc
integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bcd
real*8, dimension(3, 0:ncompo)                            ::  mtbl
real, dimension(3)                                        ::  dh

ix = sz(1)
jx = sz(2)
kx = sz(3)

rx = 1.0/dh(1)
ry = 1.0/dh(2)
rz = 1.0/dh(3)
rx2= rx*rx
ry2= ry*ry
rz2= rz*rz

res  = 0.0

hm = 1
if ( h_mode == 0) hm = 0

! 6 + 7 + 56 + 66 + 54 + 4 = 193
flop = flop + dble(ix)*dble(jx)*dble(kx)*193.0d0 + 27.0d0


!$OMP PARALLEL &
!$OMP REDUCTION(+:res) &
!$OMP FIRSTPRIVATE(ix, jx, kx, dt, rx, ry, rz, rx2, ry2, rz2, hm) &
!$OMP PRIVATE(idx, delta, sw) &
!$OMP PRIVATE(t_p, t_w, t_e, t_s, t_n, t_b, t_t) &
!$OMP PRIVATE(g_w, g_e, g_s, g_n, g_b, g_t) &
!$OMP PRIVATE(lmd_p, lmd_w, lmd_e, lmd_s, lmd_n, lmd_b, lmd_t) &
!$OMP PRIVATE(rcp_p, rcp_w, rcp_e, rcp_s, rcp_n, rcp_b, rcp_t) &
!$OMP PRIVATE(tc_w, tc_e, tc_s, tc_n, tc_b, tc_t) &
!$OMP PRIVATE(c_w, c_e, c_s, c_n, c_b, c_t) &
!$OMP PRIVATE(l_p, l_w, l_e, l_s, l_n, l_b, l_t) &
!$OMP PRIVATE(a_w, a_e, a_s, a_n, a_b, a_t)

!$OMP DO SCHEDULE(static)

do k=1,kx
do j=1,jx
do i=1,ix

idx = bcd(i,j,k)

sw = real(ibits(idx, Active,  1))
if (hm == 0) sw = 1.0 ! conjugate heat transferのときマスクなし

a_w = real(ibits(idx, adbtc_W, 1))
a_e = real(ibits(idx, adbtc_E, 1))
a_s = real(ibits(idx, adbtc_S, 1))
a_n = real(ibits(idx, adbtc_N, 1))
a_b = real(ibits(idx, adbtc_B, 1))
a_t = real(ibits(idx, adbtc_T, 1))

g_w = real(ibits(idx, gma_W, 1))
g_e = real(ibits(idx, gma_E, 1))
g_s = real(ibits(idx, gma_S, 1))
g_n = real(ibits(idx, gma_N, 1))
g_b = real(ibits(idx, gma_B, 1))
g_t = real(ibits(idx, gma_T, 1))

c_w = g_w * a_w
c_e = g_e * a_e
c_s = g_s * a_s
c_n = g_n * a_n
c_b = g_b * a_b
c_t = g_t * a_t ! 6

! 媒質ID
l_p = ibits(idx,                0, bitw_5)
l_w = ibits(bcd(i-1, j  , k  ), 0, bitw_5)
l_e = ibits(bcd(i+1, j  , k  ), 0, bitw_5)
l_s = ibits(bcd(i  , j-1, k  ), 0, bitw_5)
l_n = ibits(bcd(i  , j+1, k  ), 0, bitw_5)
l_b = ibits(bcd(i  , j  , k-1), 0, bitw_5)
l_t = ibits(bcd(i  , j  , k+1), 0, bitw_5)

rcp_p = mtbl(1, l_p) * mtbl(2, l_p)
rcp_w = mtbl(1, l_w) * mtbl(2, l_w)
rcp_e = mtbl(1, l_e) * mtbl(2, l_e)
rcp_s = mtbl(1, l_s) * mtbl(2, l_s)
rcp_n = mtbl(1, l_n) * mtbl(2, l_n)
rcp_b = mtbl(1, l_b) * mtbl(2, l_b)
rcp_t = mtbl(1, l_t) * mtbl(2, l_t) ! 7

lmd_p = mtbl(3, l_p)
lmd_w = mtbl(3, l_w)
lmd_e = mtbl(3, l_e)
lmd_s = mtbl(3, l_s)
lmd_n = mtbl(3, l_n)
lmd_b = mtbl(3, l_b)
lmd_t = mtbl(3, l_t)

t_p = ie(i  , j  , k  ) / rcp_p
t_w = ie(i-1, j  , k  ) / rcp_w
t_e = ie(i+1, j  , k  ) / rcp_e
t_s = ie(i  , j-1, k  ) / rcp_s
t_n = ie(i  , j+1, k  ) / rcp_n
t_b = ie(i  , j  , k-1) / rcp_b
t_t = ie(i  , j  , k+1) / rcp_t ! 7*8 = 56

tc_w = lmd_p * lmd_w / (lmd_p + lmd_w) * 2.0
tc_e = lmd_p * lmd_e / (lmd_p + lmd_e) * 2.0
tc_s = lmd_p * lmd_s / (lmd_p + lmd_s) * 2.0
tc_n = lmd_p * lmd_n / (lmd_p + lmd_n) * 2.0
tc_b = lmd_p * lmd_b / (lmd_p + lmd_b) * 2.0
tc_t = lmd_p * lmd_t / (lmd_p + lmd_t) * 2.0 ! (3+8)*6 = 66

delta =(rx2*( c_w * tc_w * (t_w - t_p)  & ! west
            + c_e * tc_e * (t_e - t_p)) & ! east
       +ry2*( c_s * tc_s * (t_s - t_p)  & ! south
            + c_n * tc_n * (t_n - t_p)) & ! north
       +rz2*( c_b * tc_b * (t_b - t_p)  & ! bottom
            + c_t * tc_t * (t_t - t_p)) & ! top
       +rx*( (1.0 - g_w) * a_w * qbc(1, i, j, k)  & ! west   gamma
            -(1.0 - g_e) * a_e * qbc(2, i, j, k)) & ! east   gamma
       +ry*( (1.0 - g_s) * a_s * qbc(3, i, j, k)  & ! south  gamma
            -(1.0 - g_n) * a_n * qbc(4, i, j, k)) & ! north  gamma
       +rz*( (1.0 - g_b) * a_b * qbc(5, i, j, k)  & ! bottom gamma
            +(1.0 - g_t) * a_t * qbc(6, i, j, k)) & ! top    gamma
       ) * sw ! 26 + 27 + 1 = 54
ie(i,j,k) = ws(i,j,k) + delta * dt
res = res + dble(delta*delta)
end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL

res_l2 = res

return
end subroutine ps_diff_ee
