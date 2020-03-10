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
!! @param [in]  c_scheme  対流項スキームのモード（1-UWD, 3-MUSCL）
!! @param [in]  rei       レイノルズ数の逆数
!! @param [in]  v         セルセンター速度ベクトル（n-step）
!! @param [in]  vf        セルフェイス速度ベクトル（n-step）
!! @param [in]  bid       Cut ID
!! @param [in]  bcd       BCindex B
!! @param [in]  vcs_coef  粘性項の係数（粘性項を計算しない場合には0.0）
!! @param [out] flop      浮動小数点演算数
!<
subroutine pvec_muscl (wv, sz, g, dh, c_scheme, rei, v, vf, bid, bcd, vcs_coef, flop)
implicit none
include 'ffv_f_params.h'
integer                                                   ::  i, j, k, ix, jx, kx, g, c_scheme, bix, bdx
integer, dimension(3)                                     ::  sz
double precision                                          ::  flop
real                                                      ::  b_e1, b_w1, b_n1, b_s1, b_t1, b_b1
real                                                      ::  b_e2, b_w2, b_n2, b_s2, b_t2, b_b2, b_p
real                                                      ::  UPe, UPw, VPn, VPs, WPt, WPb
real                                                      ::  Up0, Ue1, Ue2, Uw1, Uw2, Us1, Us2, Un1, Un2, Ub1, Ub2, Ut1, Ut2
real                                                      ::  Vp0, Ve1, Ve2, Vw1, Vw2, Vs1, Vs2, Vn1, Vn2, Vb1, Vb2, Vt1, Vt2
real                                                      ::  Wp0, We1, We2, Ww1, Ww2, Ws1, Ws2, Wn1, Wn2, Wb1, Wb2, Wt1, Wt2
real                                                      ::  ck, vcs, vcs_coef, rx, ry, rz, rx2, ry2 ,rz2
real                                                      ::  c_e1, c_w1, c_n1, c_s1, c_t1, c_b1, cm1, cm2, ss_4
real                                                      ::  c_e2, c_w2, c_n2, c_s2, c_t2, c_b2
real                                                      ::  dv1, dv2, dv3, dv4, g1, g2, g3, g4, g5, g6, s1, s2, s3, s4, b
real                                                      ::  Urr, Url, Ulr, Ull, Vrr, Vrl, Vlr, Vll, Wrr, Wrl, Wlr, Wll
real                                                      ::  cr, cl, acr, acl, cnv_u, cnv_v, cnv_w, EX, EY, EZ, rei
real                                                      ::  fu_r, fu_l, fv_r, fv_l, fw_r, fw_l, uq, vq, wq, ss
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v, wv, vf
real, dimension(3)                                        ::  dh
integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bid, bcd

ix = sz(1)
jx = sz(2)
kx = sz(3)

rx = 1.0/dh(1)
ry = 1.0/dh(2)
rz = 1.0/dh(3)

rx2 = rei * rx * rx
ry2 = rei * ry * ry
rz2 = rei * rz * rz

! vcs = 1.0 (Euler Explicit) / 0.5 (CN) / 0.0(No)
vcs = vcs_coef


ck = 0.0
b  = 0.0
ss = 1.0

if ( c_scheme == 1 ) then      !     1st order upwind
ss = 0.0
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


! Total : 36 + 24 + 3 + (14 + 78 * 3 + 12) + 69 + 12 = 888

flop = flop + dble(ix)*dble(jx)*dble(kx)*888.0d0 + 36.0d0


!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, rx, ry, rz, rx2, ry2 ,rz2, vcs, b, ck, ss_4, ss, cm1, cm2) &
!$OMP FIRSTPRIVATE(rei) &
!$OMP PRIVATE(cnv_u, cnv_v, cnv_w, bix, bdx, uq, vq, wq) &
!$OMP PRIVATE(Up0, Ue1, Ue2, Uw1, Uw2, Us1, Us2, Un1, Un2, Ub1, Ub2, Ut1, Ut2) &
!$OMP PRIVATE(Vp0, Ve1, Ve2, Vw1, Vw2, Vs1, Vs2, Vn1, Vn2, Vb1, Vb2, Vt1, Vt2) &
!$OMP PRIVATE(Wp0, We1, We2, Ww1, Ww2, Ws1, Ws2, Wn1, Wn2, Wb1, Wb2, Wt1, Wt2) &
!$OMP PRIVATE(b_e1, b_w1, b_n1, b_s1, b_t1, b_b1, b_e2, b_w2, b_n2, b_s2, b_t2, b_b2, b_p) &
!$OMP PRIVATE(c_e1, c_w1, c_n1, c_s1, c_t1, c_b1, c_e2, c_w2, c_n2, c_s2, c_t2, c_b2) &
!$OMP PRIVATE(UPe, UPw, VPn, VPs, WPt, WPb) &
!$OMP PRIVATE(cr, cl, acr, acl) &
!$OMP PRIVATE(dv1, dv2, dv3, dv4, g1, g2, g3, g4, g5, g6, s1, s2, s3, s4) &
!$OMP PRIVATE(Urr, Url, Ulr, Ull, Vrr, Vrl, Vlr, Vll, Wrr, Wrl, Wlr, Wll) &
!$OMP PRIVATE(fu_r, fu_l, fv_r, fv_l, fw_r, fw_l) &
!$OMP PRIVATE(EX, EY, EZ)

!$OMP DO SCHEDULE(static) COLLAPSE(2)

do k=1,kx
do j=1,jx
do i=1,ix
cnv_u = 0.0
cnv_v = 0.0
cnv_w = 0.0

! 各軸方向5点の変数ロード
Ub2 = v(i  ,j  ,k-2, 1)
Ub1 = v(i  ,j  ,k-1, 1)
Us2 = v(i  ,j-2,k  , 1)
Us1 = v(i  ,j-1,k  , 1)
Uw2 = v(i-2,j  ,k  , 1)
Uw1 = v(i-1,j  ,k  , 1)
Up0 = v(i  ,j  ,k  , 1)
Ue1 = v(i+1,j  ,k  , 1)
Ue2 = v(i+2,j  ,k  , 1)
Un1 = v(i  ,j+1,k  , 1)
Un2 = v(i  ,j+2,k  , 1)
Ut1 = v(i  ,j  ,k+1, 1)
Ut2 = v(i  ,j  ,k+2, 1)

Vb2 = v(i  ,j  ,k-2, 2)
Vb1 = v(i  ,j  ,k-1, 2)
Vs2 = v(i  ,j-2,k  , 2)
Vs1 = v(i  ,j-1,k  , 2)
Vw2 = v(i-2,j  ,k  , 2)
Vw1 = v(i-1,j  ,k  , 2)
Vp0 = v(i  ,j  ,k  , 2)
Ve1 = v(i+1,j  ,k  , 2)
Ve2 = v(i+2,j  ,k  , 2)
Vn1 = v(i  ,j+1,k  , 2)
Vn2 = v(i  ,j+2,k  , 2)
Vt1 = v(i  ,j  ,k+1, 2)
Vt2 = v(i  ,j  ,k+2, 2)

Wb2 = v(i  ,j  ,k-2, 3)
Wb1 = v(i  ,j  ,k-1, 3)
Ws2 = v(i  ,j-2,k  , 3)
Ws1 = v(i  ,j-1,k  , 3)
Ww2 = v(i-2,j  ,k  , 3)
Ww1 = v(i-1,j  ,k  , 3)
Wp0 = v(i  ,j  ,k  , 3)
We1 = v(i+1,j  ,k  , 3)
We2 = v(i+2,j  ,k  , 3)
Wn1 = v(i  ,j+1,k  , 3)
Wn2 = v(i  ,j+2,k  , 3)
Wt1 = v(i  ,j  ,k+1, 3)
Wt2 = v(i  ,j  ,k+2, 3)

bix = bid(i,j,k)
bdx = bcd(i,j,k)

! (i,j,k)からみたセル状態 (0-solid / 1-fluid)
b_p = real(ibits(bdx, State, 1))

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
c_w1 = real( ibits(bdx, bc_d_W, 1) )
c_e1 = real( ibits(bdx, bc_d_E, 1) )
c_s1 = real( ibits(bdx, bc_d_S, 1) )
c_n1 = real( ibits(bdx, bc_d_N, 1) )
c_b1 = real( ibits(bdx, bc_d_B, 1) )
c_t1 = real( ibits(bdx, bc_d_T, 1) )


! ステンシルの参照先がvspec, outflowである場合のスキームの破綻を回避，１次精度におとすSW
c_w2 = real( ibits(bcd(i-1, j  , k  ), bc_d_W, 1) )
c_e2 = real( ibits(bcd(i+1, j  , k  ), bc_d_E, 1) )
c_s2 = real( ibits(bcd(i  , j-1, k  ), bc_d_S, 1) )
c_n2 = real( ibits(bcd(i  , j+1, k  ), bc_d_N, 1) )
c_b2 = real( ibits(bcd(i  , j  , k-1), bc_d_B, 1) )
c_t2 = real( ibits(bcd(i  , j  , k+1), bc_d_T, 1) )



! 界面速度（スタガード位置） > 24 flops
UPe = vf(i  , j  , k  ,1)*b_e1
UPw = vf(i-1, j  , k  ,1)*b_w1
VPn = vf(i  , j  , k  ,2)*b_n1
VPs = vf(i  , j-1, k  ,2)*b_s1
WPt = vf(i  , j  , k  ,3)*b_t1
WPb = vf(i  , j  , k-1,3)*b_b1


! セルセンターからの壁面修正速度 > 3 flops
uq =  - Up0
vq =  - Vp0
wq =  - Wp0

! X方向 ---------------------------------------

! 速度指定の場合にMUSCLスキームの参照先として，固体内にテンポラリに与えた値を使う
if ( (b_e2 == 0.0)  ) then  ! 7 flops
Ue2 =  - v(i+1,j  ,k  , 1)
Ve2 =  - v(i+1,j  ,k  , 2)
We2 =  - v(i+1,j  ,k  , 3)
endif

if ( b_e1 == 0.0 ) then
Ue1 = uq
Ve1 = vq
We1 = wq
endif

if ( b_w1 == 0.0 ) then
Uw1 = uq
Vw1 = vq
Ww1 = wq
end if

if ( (b_w2 == 0.0)  ) then ! 7 flops
Uw2 = - v(i-1,j  ,k  , 1)
Vw2 = - v(i-1,j  ,k  , 2)
Ww2 = - v(i-1,j  ,k  , 3)
end if

! 流束　流体のみと固体壁の影響を含み，隣接セルが固体の場合にはマスクする
cr  = UPe
cl  = UPw
acr = abs(cr)
acl = abs(cl)

dv4 = Ue2-Ue1
dv3 = Ue1-Up0
dv2 = Up0-Uw1
dv1 = Uw1-Uw2

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

Urr = Ue1 - (cm1*g6+cm2*g5)*ss_4 * c_e2
Url = Up0 + (cm1*g3+cm2*g4)*ss_4
Ulr = Up0 - (cm1*g4+cm2*g3)*ss_4
Ull = Uw1 + (cm1*g1+cm2*g2)*ss_4 * c_w2
fu_r = 0.5*(cr*(Urr+Url) - acr*(Urr-Url)) * b_e1
fu_l = 0.5*(cl*(Ulr+Ull) - acl*(Ulr-Ull)) * b_w1 ! > 4 + 4 + 36 + 5*4+7*2 = 78 flops

dv4 = Ve2-Ve1
dv3 = Ve1-Vp0
dv2 = Vp0-Vw1
dv1 = Vw1-Vw2

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

Vrr = Ve1 - (cm1*g6+cm2*g5)*ss_4 * c_e2
Vrl = Vp0 + (cm1*g3+cm2*g4)*ss_4
Vlr = Vp0 - (cm1*g4+cm2*g3)*ss_4
Vll = Vw1 + (cm1*g1+cm2*g2)*ss_4 * c_w2
fv_r = 0.5*(cr*(Vrr+Vrl) - acr*(Vrr-Vrl)) * b_e1
fv_l = 0.5*(cl*(Vlr+Vll) - acl*(Vlr-Vll)) * b_w1

dv4 = We2-We1
dv3 = We1-Wp0
dv2 = Wp0-Ww1
dv1 = Ww1-Ww2

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

Wrr = We1 - (cm1*g6+cm2*g5)*ss_4 * c_e2
Wrl = Wp0 + (cm1*g3+cm2*g4)*ss_4
Wlr = Wp0 - (cm1*g4+cm2*g3)*ss_4
Wll = Ww1 + (cm1*g1+cm2*g2)*ss_4 * c_w2
fw_r = 0.5*(cr*(Wrr+Wrl) - acr*(Wrr-Wrl)) * b_e1
fw_l = 0.5*(cl*(Wlr+Wll) - acl*(Wlr-Wll)) * b_w1

! 流束の加算　VBCでない面の寄与のみを評価する
cnv_u = cnv_u + fu_r*c_e1 - fu_l*c_w1
cnv_v = cnv_v + fv_r*c_e1 - fv_l*c_w1
cnv_w = cnv_w + fw_r*c_e1 - fw_l*c_w1 ! > 4*3 = 12 flops



! Y方向 ---------------------------------------

if ( (b_n2 == 0.0)  ) then
Un2 = - v(i  ,j+1,k  , 1)
Vn2 = - v(i  ,j+1,k  , 2)
Wn2 = - v(i  ,j+1,k  , 3)
endif

if ( b_n1 == 0.0 ) then
Un1 = uq
Vn1 = vq
Wn1 = wq
endif

if ( b_s1 == 0.0 ) then
Us1 = uq
Vs1 = vq
Ws1 = wq
endif

if ( (b_s2 == 0.0)  ) then
Us2 = - v(i  ,j-1,k  , 1)
Vs2 = - v(i  ,j-1,k  , 2)
Ws2 = - v(i  ,j-1,k  , 3)
endif

! 流束　流体のみと固体壁の影響を含み，隣接セルが固体の場合にはマスクする
cr  = VPn
cl  = VPs
acr = abs(cr)
acl = abs(cl)

dv4 = Un2-Un1
dv3 = Un1-Up0
dv2 = Up0-Us1
dv1 = Us1-Us2

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

Urr = Un1 - (cm1*g6+cm2*g5)*ss_4 * c_n2
Url = Up0 + (cm1*g3+cm2*g4)*ss_4
Ulr = Up0 - (cm1*g4+cm2*g3)*ss_4
Ull = Us1 + (cm1*g1+cm2*g2)*ss_4 * c_s2
fu_r = 0.5*(cr*(Urr+Url) - acr*(Urr-Url)) * b_n1
fu_l = 0.5*(cl*(Ulr+Ull) - acl*(Ulr-Ull)) * b_s1

dv4 = Vn2-Vn1
dv3 = Vn1-Vp0
dv2 = Vp0-Vs1
dv1 = Vs1-Vs2

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

Vrr = Vn1 - (cm1*g6+cm2*g5)*ss_4 * c_n2
Vrl = Vp0 + (cm1*g3+cm2*g4)*ss_4
Vlr = Vp0 - (cm1*g4+cm2*g3)*ss_4
Vll = Vs1 + (cm1*g1+cm2*g2)*ss_4 * c_s2
fv_r = 0.5*(cr*(Vrr+Vrl) - acr*(Vrr-Vrl)) * b_n1
fv_l = 0.5*(cl*(Vlr+Vll) - acl*(Vlr-Vll)) * b_s1

dv4 = Wn2-Wn1
dv3 = Wn1-Wp0
dv2 = Wp0-Ws1
dv1 = Ws1-Ws2

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

Wrr = Wn1 - (cm1*g6+cm2*g5)*ss_4 * c_n2
Wrl = Wp0 + (cm1*g3+cm2*g4)*ss_4
Wlr = Wp0 - (cm1*g4+cm2*g3)*ss_4
Wll = Ws1 + (cm1*g1+cm2*g2)*ss_4 * c_s2
fw_r = 0.5*(cr*(Wrr+Wrl) - acr*(Wrr-Wrl)) * b_n1
fw_l = 0.5*(cl*(Wlr+Wll) - acl*(Wlr-Wll)) * b_s1

! 流束の加算　VBCでない面の寄与のみを評価する
cnv_u = cnv_u + fu_r*c_n1 - fu_l*c_s1
cnv_v = cnv_v + fv_r*c_n1 - fv_l*c_s1
cnv_w = cnv_w + fw_r*c_n1 - fw_l*c_s1



! Z方向 ---------------------------------------

! 壁面の場合の参照速度の修正
if ( (b_t2 == 0.0)  ) then
Ut2 = - v(i  ,j  ,k+1, 1)
Vt2 = - v(i  ,j  ,k+1, 2)
Wt2 = - v(i  ,j  ,k+1, 3)
end if

if ( b_t1 == 0.0 ) then
Ut1 = uq
Vt1 = vq
Wt1 = wq
end if

if ( b_b1 == 0.0 ) then
Ub1 = uq
Vb1 = vq
Wb1 = wq
end if

if ( (b_b2 == 0.0)  ) then
Ub2 = - v(i  ,j  ,k-1, 1)
Vb2 = - v(i  ,j  ,k-1, 2)
Wb2 = - v(i  ,j  ,k-1, 3)
end if

! 流束　流体のみと固体壁の影響を含み，隣接セルが固体の場合にはマスクする
cr  = WPt
cl  = WPb
acr = abs(cr)
acl = abs(cl)

dv4 = Ut2-Ut1
dv3 = Ut1-Up0
dv2 = Up0-Ub1
dv1 = Ub1-Ub2

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

Urr = Ut1 - (cm1*g6+cm2*g5)*ss_4 * c_t2
Url = Up0 + (cm1*g3+cm2*g4)*ss_4
Ulr = Up0 - (cm1*g4+cm2*g3)*ss_4
Ull = Ub1 + (cm1*g1+cm2*g2)*ss_4 * c_b2
fu_r = 0.5*(cr*(Urr+Url) - acr*(Urr-Url)) * b_t1
fu_l = 0.5*(cl*(Ulr+Ull) - acl*(Ulr-Ull)) * b_b1

dv4 = Vt2-Vt1
dv3 = Vt1-Vp0
dv2 = Vp0-Vb1
dv1 = Vb1-Vb2

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

Vrr = Vt1 - (cm1*g6+cm2*g5)*ss_4 * c_t2
Vrl = Vp0 + (cm1*g3+cm2*g4)*ss_4
Vlr = Vp0 - (cm1*g4+cm2*g3)*ss_4
Vll = Vb1 + (cm1*g1+cm2*g2)*ss_4 * c_b2
fv_r = 0.5*(cr*(Vrr+Vrl) - acr*(Vrr-Vrl)) * b_t1
fv_l = 0.5*(cl*(Vlr+Vll) - acl*(Vlr-Vll)) * b_b1

dv4 = Wt2-Wt1
dv3 = Wt1-Wp0
dv2 = Wp0-Wb1
dv1 = Wb1-Wb2

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

Wrr = Wt1 - (cm1*g6+cm2*g5)*ss_4 * c_t2
Wrl = Wp0 + (cm1*g3+cm2*g4)*ss_4
Wlr = Wp0 - (cm1*g4+cm2*g3)*ss_4
Wll = Wb1 + (cm1*g1+cm2*g2)*ss_4 * c_b2
fw_r = 0.5*(cr*(Wrr+Wrl) - acr*(Wrr-Wrl)) * b_t1
fw_l = 0.5*(cl*(Wlr+Wll) - acl*(Wlr-Wll)) * b_b1

! 流束の加算　VBCでない面の寄与のみを評価する
cnv_u = cnv_u + fu_r*c_t1 - fu_l*c_b1
cnv_v = cnv_v + fv_r*c_t1 - fv_l*c_b1
cnv_w = cnv_w + fw_r*c_t1 - fw_l*c_b1



! 粘性項の計算　セル界面の剪断力を計算し，必要に応じて置換する 23*3 = 69 flops
EX =  ( Ue1 - Up0 ) * c_e1 * rx2 &
    - ( Up0 - Uw1 ) * c_w1 * rx2 &
    + ( Un1 - Up0 ) * c_n1 * ry2 &
    - ( Up0 - Us1 ) * c_s1 * ry2 &
    + ( Ut1 - Up0 ) * c_t1 * rz2 &
    - ( Up0 - Ub1 ) * c_b1 * rz2

EY =  ( Ve1 - Vp0 ) * c_e1 * rx2 &
    - ( Vp0 - Vw1 ) * c_w1 * rx2 &
    + ( Vn1 - Vp0 ) * c_n1 * ry2 &
    - ( Vp0 - Vs1 ) * c_s1 * ry2 &
    + ( Vt1 - Vp0 ) * c_t1 * rz2 &
    - ( Vp0 - Vb1 ) * c_b1 * rz2

EZ =  ( We1 - Wp0 ) * c_e1 * rx2 &
    - ( Wp0 - Ww1 ) * c_w1 * rx2 &
    + ( Wn1 - Wp0 ) * c_n1 * ry2 &
    - ( Wp0 - Ws1 ) * c_s1 * ry2 &
    + ( Wt1 - Wp0 ) * c_t1 * rz2 &
    - ( Wp0 - Wb1 ) * c_b1 * rz2

! 対流項と粘性項の和 > 4*3 = 12 flops
wv(i,j,k,1) = -cnv_u * rx + EX*vcs
wv(i,j,k,2) = -cnv_v * ry + EY*vcs
wv(i,j,k,3) = -cnv_w * rz + EZ*vcs
end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL

return
end subroutine pvec_muscl


!> ********************************************************************
!! @brief 対流項と粘性項の計算
!! @param [out] wv        疑似ベクトルの空間項
!! @param [in]  sz        配列長
!! @param [in]  g         ガイドセル長
!! @param [in]  dh        格子幅
!! @param [in]  c_scheme  対流項スキームのモード（1-UWD, 3-MUSCL）
!! @param [in]  rei       レイノルズ数の逆数
!! @param [in]  v         セルセンター速度ベクトル（n-step）
!! @param [in]  vf        セルフェイス速度ベクトル（n-step）
!! @param [in]  bid       Cut ID
!! @param [in]  bcd       BCindex B
!! @param [in]  vcs_coef  粘性項の係数（粘性項を計算しない場合には0.0）
!! @param [in]  Cs        定数CS
!! @param [in]  imodel    乱流モデル
!! @param [in]  m_nu      動粘性係数
!! @param [in]  rho       密度
!! @param [out] flop      浮動小数点演算数
!<
subroutine pvec_muscl_les (wv, sz, g, dh, c_scheme, rei, v, vf, bid, bcd, vcs_coef, Cs, imodel, m_nu, rho, flop)
implicit none
include 'ffv_f_params.h'
integer                                                   ::  i, j, k, ix, jx, kx, g, c_scheme, bix, bdx
integer, dimension(3)                                     ::  sz
double precision                                          ::  flop
real                                                      ::  b_e1, b_w1, b_n1, b_s1, b_t1, b_b1
real                                                      ::  b_e2, b_w2, b_n2, b_s2, b_t2, b_b2, b_p
real                                                      ::  UPe, UPw, VPn, VPs, WPt, WPb
real                                                      ::  Up0, Ue1, Ue2, Uw1, Uw2, Us1, Us2, Un1, Un2, Ub1, Ub2, Ut1, Ut2
real                                                      ::  Vp0, Ve1, Ve2, Vw1, Vw2, Vs1, Vs2, Vn1, Vn2, Vb1, Vb2, Vt1, Vt2
real                                                      ::  Wp0, We1, We2, Ww1, Ww2, Ws1, Ws2, Wn1, Wn2, Wb1, Wb2, Wt1, Wt2
real                                                      ::  ck, vcs, vcs_coef, rx, ry, rz, rx2, ry2 ,rz2, dx, dy, dz
real                                                      ::  c_e1, c_w1, c_n1, c_s1, c_t1, c_b1, cm1, cm2, ss_4
real                                                      ::  c_e2, c_w2, c_n2, c_s2, c_t2, c_b2
real                                                      ::  dv1, dv2, dv3, dv4, g1, g2, g3, g4, g5, g6, s1, s2, s3, s4, b
real                                                      ::  Urr, Url, Ulr, Ull, Vrr, Vrl, Vlr, Vll, Wrr, Wrl, Wlr, Wll
real                                                      ::  cr, cl, acr, acl, cnv_u, cnv_v, cnv_w, EX, EY, EZ, rei
real                                                      ::  fu_r, fu_l, fv_r, fv_l, fw_r, fw_l, uq, vq, wq, ss
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v, wv, vf
real, dimension(3)                                        ::  dh
integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bid, bcd
integer                                                   ::  imodel
real                                                      ::  Cs, nu, rho, Cw, m_nu, r_nu
real                                                      ::  DUDX, DUDY, DUDZ
real                                                      ::  DVDX, DVDY, DVDZ
real                                                      ::  DWDX, DWDY, DWDZ
real                                                      ::  fs, DUDY_w, tauw, utau, yc, yp, up, min_h
real                                                      ::  S11, S12, S13, S21, S22, S23, S31, S32, S33, SSS
real                                                      ::  W11, W12, W13, W21, W22, W23, W31, W32, W33, WWW
real                                                      ::  S_1, S_2, S_3, W_1, W_2, W_3
real                                                      ::  S11d, S12d, S13d, S21d, S22d, S23d, S31d, S32d, S33d
real                                                      ::  Fcs, E_csm, Q_csm, Sijd2, nut
double precision                                          ::  EPS


ix = sz(1)
jx = sz(2)
kx = sz(3)

dx = dh(1)
dy = dh(2)
dz = dh(3)

rx = 1.0/dh(1)
ry = 1.0/dh(2)
rz = 1.0/dh(3)

rx2= rei * rx * rx
ry2= rei * ry * ry
rz2= rei * rz * rz

EPS = 1.0d-10
Cw  = 0.325d0
nu  = m_nu
r_nu= 1.0d0/nu

! vcs = 1.0 (Euler Explicit) / 0.5 (CN) / 0.0(No)
vcs = vcs_coef

ck = 0.0
b  = 0.0
ss = 1.0

if ( c_scheme == 1 ) then      !     1st order upwind
ss = 0.0
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


flop = flop + 36.0

!$OMP PARALLEL &
!$OMP REDUCTION(+:flop) &
!$OMP FIRSTPRIVATE(ix, jx, kx, rx, ry, rz, rx2, ry2 ,rz2, vcs, b, ck, ss_4, ss, cm1, cm2) &
!$OMP FIRSTPRIVATE(rei, dx, dy, dz) &
!$OMP PRIVATE(cnv_u, cnv_v, cnv_w, bix, bdx, uq, vq, wq) &
!$OMP PRIVATE(Up0, Ue1, Ue2, Uw1, Uw2, Us1, Us2, Un1, Un2, Ub1, Ub2, Ut1, Ut2) &
!$OMP PRIVATE(Vp0, Ve1, Ve2, Vw1, Vw2, Vs1, Vs2, Vn1, Vn2, Vb1, Vb2, Vt1, Vt2) &
!$OMP PRIVATE(Wp0, We1, We2, Ww1, Ww2, Ws1, Ws2, Wn1, Wn2, Wb1, Wb2, Wt1, Wt2) &
!$OMP PRIVATE(b_e1, b_w1, b_n1, b_s1, b_t1, b_b1, b_e2, b_w2, b_n2, b_s2, b_t2, b_b2, b_p) &
!$OMP PRIVATE(c_e1, c_w1, c_n1, c_s1, c_t1, c_b1, c_e2, c_w2, c_n2, c_s2, c_t2, c_b2) &
!$OMP PRIVATE(UPe, UPw, VPn, VPs, WPt, WPb) &
!$OMP PRIVATE(cr, cl, acr, acl) &
!$OMP PRIVATE(dv1, dv2, dv3, dv4, g1, g2, g3, g4, g5, g6, s1, s2, s3, s4) &
!$OMP PRIVATE(Urr, Url, Ulr, Ull, Vrr, Vrl, Vlr, Vll, Wrr, Wrl, Wlr, Wll) &
!$OMP PRIVATE(fu_r, fu_l, fv_r, fv_l, fw_r, fw_l) &
!$OMP PRIVATE(EX, EY, EZ) &
!$OMP FIRSTPRIVATE(imodel, EPS, Cs, nu, rho, Cw, r_nu) &
!$OMP PRIVATE(DUDX, DUDY, DUDZ, DVDX, DVDY, DVDZ, DWDX, DWDY, DWDZ) &
!$OMP PRIVATE(S11, S12, S13, S21, S22, S23, S31, S32, S33, SSS) &
!$OMP PRIVATE(W11, W12, W13, W21, W22, W23, W31, W32, W33, WWW) &
!$OMP PRIVATE(S_1, S_2, S_3, W_1, W_2, W_3) &
!$OMP PRIVATE(S11d, S12d, S13d, S21d, S22d, S23d, S31d, S32d, S33d) &
!$OMP PRIVATE(Fcs, E_csm, Q_csm, Sijd2, nut) &
!$OMP PRIVATE(fs, DUDY_w, tauw, utau, yc, yp, up, min_h)

!$OMP DO SCHEDULE(static) COLLAPSE(2)

do k=1,kx
do j=1,jx
do i=1,ix
cnv_u = 0.0
cnv_v = 0.0
cnv_w = 0.0

! 各軸方向5点の変数ロード
Ub2 = v(i  ,j  ,k-2, 1)
Ub1 = v(i  ,j  ,k-1, 1)
Us2 = v(i  ,j-2,k  , 1)
Us1 = v(i  ,j-1,k  , 1)
Uw2 = v(i-2,j  ,k  , 1)
Uw1 = v(i-1,j  ,k  , 1)
Up0 = v(i  ,j  ,k  , 1)
Ue1 = v(i+1,j  ,k  , 1)
Ue2 = v(i+2,j  ,k  , 1)
Un1 = v(i  ,j+1,k  , 1)
Un2 = v(i  ,j+2,k  , 1)
Ut1 = v(i  ,j  ,k+1, 1)
Ut2 = v(i  ,j  ,k+2, 1)

Vb2 = v(i  ,j  ,k-2, 2)
Vb1 = v(i  ,j  ,k-1, 2)
Vs2 = v(i  ,j-2,k  , 2)
Vs1 = v(i  ,j-1,k  , 2)
Vw2 = v(i-2,j  ,k  , 2)
Vw1 = v(i-1,j  ,k  , 2)
Vp0 = v(i  ,j  ,k  , 2)
Ve1 = v(i+1,j  ,k  , 2)
Ve2 = v(i+2,j  ,k  , 2)
Vn1 = v(i  ,j+1,k  , 2)
Vn2 = v(i  ,j+2,k  , 2)
Vt1 = v(i  ,j  ,k+1, 2)
Vt2 = v(i  ,j  ,k+2, 2)

Wb2 = v(i  ,j  ,k-2, 3)
Wb1 = v(i  ,j  ,k-1, 3)
Ws2 = v(i  ,j-2,k  , 3)
Ws1 = v(i  ,j-1,k  , 3)
Ww2 = v(i-2,j  ,k  , 3)
Ww1 = v(i-1,j  ,k  , 3)
Wp0 = v(i  ,j  ,k  , 3)
We1 = v(i+1,j  ,k  , 3)
We2 = v(i+2,j  ,k  , 3)
Wn1 = v(i  ,j+1,k  , 3)
Wn2 = v(i  ,j+2,k  , 3)
Wt1 = v(i  ,j  ,k+1, 3)
Wt2 = v(i  ,j  ,k+2, 3)

bix = bid(i,j,k)
bdx = bcd(i,j,k)

! (i,j,k)からみたセル状態 (0-solid / 1-fluid)
b_p = real(ibits(bdx, State, 1))

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
c_w1 = real( ibits(bdx, bc_d_W, 1) )
c_e1 = real( ibits(bdx, bc_d_E, 1) )
c_s1 = real( ibits(bdx, bc_d_S, 1) )
c_n1 = real( ibits(bdx, bc_d_N, 1) )
c_b1 = real( ibits(bdx, bc_d_B, 1) )
c_t1 = real( ibits(bdx, bc_d_T, 1) )


! ステンシルの参照先がvspec, outflowである場合のスキームの破綻を回避，１次精度におとすSW
c_w2 = real( ibits(bcd(i-1, j  , k  ), bc_d_W, 1) )
c_e2 = real( ibits(bcd(i+1, j  , k  ), bc_d_E, 1) )
c_s2 = real( ibits(bcd(i  , j-1, k  ), bc_d_S, 1) )
c_n2 = real( ibits(bcd(i  , j+1, k  ), bc_d_N, 1) )
c_b2 = real( ibits(bcd(i  , j  , k-1), bc_d_B, 1) )
c_t2 = real( ibits(bcd(i  , j  , k+1), bc_d_T, 1) )


! 界面速度（スタガード位置） > 24 flops
UPe = vf(i  , j  , k  ,1)*b_e1
UPw = vf(i-1, j  , k  ,1)*b_w1
VPn = vf(i  , j  , k  ,2)*b_n1
VPs = vf(i  , j-1, k  ,2)*b_s1
WPt = vf(i  , j  , k  ,3)*b_t1
WPb = vf(i  , j  , k-1,3)*b_b1


! セルセンターからの壁面修正速度 > 3 flops
uq = - Up0
vq = - Vp0
wq = - Wp0

! X方向 ---------------------------------------

! 速度指定の場合にMUSCLスキームの参照先として，固体内にテンポラリに与えた値を使う
if ( (b_e2 == 0.0)  ) then  ! 7 flops
Ue2 = - v(i+1,j  ,k  , 1)
Ve2 = - v(i+1,j  ,k  , 2)
We2 = - v(i+1,j  ,k  , 3)
endif

if ( b_e1 == 0.0 ) then
Ue1 = uq
Ve1 = vq
We1 = wq
endif

if ( b_w1 == 0.0 ) then
Uw1 = uq
Vw1 = vq
Ww1 = wq
end if

if ( (b_w2 == 0.0)  ) then ! 7 flops
Uw2 = - v(i-1,j  ,k  , 1)
Vw2 = - v(i-1,j  ,k  , 2)
Ww2 = - v(i-1,j  ,k  , 3)
end if

! 流束　流体のみと固体壁の影響を含み，隣接セルが固体の場合にはマスクする
cr  = UPe
cl  = UPw
acr = abs(cr)
acl = abs(cl)

dv4 = Ue2-Ue1
dv3 = Ue1-Up0
dv2 = Up0-Uw1
dv1 = Uw1-Uw2

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

Urr = Ue1 - (cm1*g6+cm2*g5)*ss_4 * c_e2
Url = Up0 + (cm1*g3+cm2*g4)*ss_4
Ulr = Up0 - (cm1*g4+cm2*g3)*ss_4
Ull = Uw1 + (cm1*g1+cm2*g2)*ss_4 * c_w2
fu_r = 0.5*(cr*(Urr+Url) - acr*(Urr-Url)) * b_e1
fu_l = 0.5*(cl*(Ulr+Ull) - acl*(Ulr-Ull)) * b_w1 ! > 4 + 4 + 36 + 5*4+7*2 = 78 flops

dv4 = Ve2-Ve1
dv3 = Ve1-Vp0
dv2 = Vp0-Vw1
dv1 = Vw1-Vw2

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

Vrr = Ve1 - (cm1*g6+cm2*g5)*ss_4 * c_e2
Vrl = Vp0 + (cm1*g3+cm2*g4)*ss_4
Vlr = Vp0 - (cm1*g4+cm2*g3)*ss_4
Vll = Vw1 + (cm1*g1+cm2*g2)*ss_4 * c_w2
fv_r = 0.5*(cr*(Vrr+Vrl) - acr*(Vrr-Vrl)) * b_e1
fv_l = 0.5*(cl*(Vlr+Vll) - acl*(Vlr-Vll)) * b_w1

dv4 = We2-We1
dv3 = We1-Wp0
dv2 = Wp0-Ww1
dv1 = Ww1-Ww2

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

Wrr = We1 - (cm1*g6+cm2*g5)*ss_4 * c_e2
Wrl = Wp0 + (cm1*g3+cm2*g4)*ss_4
Wlr = Wp0 - (cm1*g4+cm2*g3)*ss_4
Wll = Ww1 + (cm1*g1+cm2*g2)*ss_4 * c_w2
fw_r = 0.5*(cr*(Wrr+Wrl) - acr*(Wrr-Wrl)) * b_e1
fw_l = 0.5*(cl*(Wlr+Wll) - acl*(Wlr-Wll)) * b_w1

! 流束の加算　VBCでない面の寄与のみを評価する
cnv_u = cnv_u + fu_r*c_e1 - fu_l*c_w1
cnv_v = cnv_v + fv_r*c_e1 - fv_l*c_w1
cnv_w = cnv_w + fw_r*c_e1 - fw_l*c_w1 ! > 4*3 = 12 flops

! Y方向 ---------------------------------------

if ( (b_n2 == 0.0)  ) then
Un2 = - v(i  ,j+1,k  , 1)
Vn2 = - v(i  ,j+1,k  , 2)
Wn2 = - v(i  ,j+1,k  , 3)
endif

if ( b_n1 == 0.0 ) then
Un1 = uq
Vn1 = vq
Wn1 = wq
endif

if ( b_s1 == 0.0 ) then
Us1 = uq
Vs1 = vq
Ws1 = wq
endif

if ( (b_s2 == 0.0)  ) then
Us2 = - v(i  ,j-1,k  , 1)
Vs2 = - v(i  ,j-1,k  , 2)
Ws2 = - v(i  ,j-1,k  , 3)
endif

! 流束　流体のみと固体壁の影響を含み，隣接セルが固体の場合にはマスクする
cr  = VPn
cl  = VPs
acr = abs(cr)
acl = abs(cl)

dv4 = Un2-Un1
dv3 = Un1-Up0
dv2 = Up0-Us1
dv1 = Us1-Us2

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

Urr = Un1 - (cm1*g6+cm2*g5)*ss_4 * c_n2
Url = Up0 + (cm1*g3+cm2*g4)*ss_4
Ulr = Up0 - (cm1*g4+cm2*g3)*ss_4
Ull = Us1 + (cm1*g1+cm2*g2)*ss_4 * c_s2
fu_r = 0.5*(cr*(Urr+Url) - acr*(Urr-Url)) * b_n1
fu_l = 0.5*(cl*(Ulr+Ull) - acl*(Ulr-Ull)) * b_s1

dv4 = Vn2-Vn1
dv3 = Vn1-Vp0
dv2 = Vp0-Vs1
dv1 = Vs1-Vs2

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

Vrr = Vn1 - (cm1*g6+cm2*g5)*ss_4 * c_n2
Vrl = Vp0 + (cm1*g3+cm2*g4)*ss_4
Vlr = Vp0 - (cm1*g4+cm2*g3)*ss_4
Vll = Vs1 + (cm1*g1+cm2*g2)*ss_4 * c_s2
fv_r = 0.5*(cr*(Vrr+Vrl) - acr*(Vrr-Vrl)) * b_n1
fv_l = 0.5*(cl*(Vlr+Vll) - acl*(Vlr-Vll)) * b_s1

dv4 = Wn2-Wn1
dv3 = Wn1-Wp0
dv2 = Wp0-Ws1
dv1 = Ws1-Ws2

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

Wrr = Wn1 - (cm1*g6+cm2*g5)*ss_4 * c_n2
Wrl = Wp0 + (cm1*g3+cm2*g4)*ss_4
Wlr = Wp0 - (cm1*g4+cm2*g3)*ss_4
Wll = Ws1 + (cm1*g1+cm2*g2)*ss_4 * c_s2
fw_r = 0.5*(cr*(Wrr+Wrl) - acr*(Wrr-Wrl)) * b_n1
fw_l = 0.5*(cl*(Wlr+Wll) - acl*(Wlr-Wll)) * b_s1

! 流束の加算　VBCでない面の寄与のみを評価する
cnv_u = cnv_u + fu_r*c_n1 - fu_l*c_s1
cnv_v = cnv_v + fv_r*c_n1 - fv_l*c_s1
cnv_w = cnv_w + fw_r*c_n1 - fw_l*c_s1

! Z方向 ---------------------------------------

! 壁面の場合の参照速度の修正
if ( (b_t2 == 0.0)  ) then
Ut2 = - v(i  ,j  ,k+1, 1)
Vt2 = - v(i  ,j  ,k+1, 2)
Wt2 = - v(i  ,j  ,k+1, 3)
end if

if ( b_t1 == 0.0 ) then
Ut1 = uq
Vt1 = vq
Wt1 = wq
end if

if ( b_b1 == 0.0 ) then
Ub1 = uq
Vb1 = vq
Wb1 = wq
end if

if ( (b_b2 == 0.0)  ) then
Ub2 = - v(i  ,j  ,k-1, 1)
Vb2 = - v(i  ,j  ,k-1, 2)
Wb2 = - v(i  ,j  ,k-1, 3)
end if

! 流束　流体のみと固体壁の影響を含み，隣接セルが固体の場合にはマスクする
cr  = WPt
cl  = WPb
acr = abs(cr)
acl = abs(cl)

dv4 = Ut2-Ut1
dv3 = Ut1-Up0
dv2 = Up0-Ub1
dv1 = Ub1-Ub2

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

Urr = Ut1 - (cm1*g6+cm2*g5)*ss_4 * c_t2
Url = Up0 + (cm1*g3+cm2*g4)*ss_4
Ulr = Up0 - (cm1*g4+cm2*g3)*ss_4
Ull = Ub1 + (cm1*g1+cm2*g2)*ss_4 * c_b2
fu_r = 0.5*(cr*(Urr+Url) - acr*(Urr-Url)) * b_t1
fu_l = 0.5*(cl*(Ulr+Ull) - acl*(Ulr-Ull)) * b_b1

dv4 = Vt2-Vt1
dv3 = Vt1-Vp0
dv2 = Vp0-Vb1
dv1 = Vb1-Vb2

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

Vrr = Vt1 - (cm1*g6+cm2*g5)*ss_4 * c_t2
Vrl = Vp0 + (cm1*g3+cm2*g4)*ss_4
Vlr = Vp0 - (cm1*g4+cm2*g3)*ss_4
Vll = Vb1 + (cm1*g1+cm2*g2)*ss_4 * c_b2
fv_r = 0.5*(cr*(Vrr+Vrl) - acr*(Vrr-Vrl)) * b_t1
fv_l = 0.5*(cl*(Vlr+Vll) - acl*(Vlr-Vll)) * b_b1

dv4 = Wt2-Wt1
dv3 = Wt1-Wp0
dv2 = Wp0-Wb1
dv1 = Wb1-Wb2

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

Wrr = Wt1 - (cm1*g6+cm2*g5)*ss_4 * c_t2
Wrl = Wp0 + (cm1*g3+cm2*g4)*ss_4
Wlr = Wp0 - (cm1*g4+cm2*g3)*ss_4
Wll = Wb1 + (cm1*g1+cm2*g2)*ss_4 * c_b2
fw_r = 0.5*(cr*(Wrr+Wrl) - acr*(Wrr-Wrl)) * b_t1
fw_l = 0.5*(cl*(Wlr+Wll) - acl*(Wlr-Wll)) * b_b1

! 流束の加算　VBCでない面の寄与のみを評価する
cnv_u = cnv_u + fu_r*c_t1 - fu_l*c_b1
cnv_v = cnv_v + fv_r*c_t1 - fv_l*c_b1
cnv_w = cnv_w + fw_r*c_t1 - fw_l*c_b1



! 粘性項の評価   27 flops
DUDX = ( Ue1 - Uw1 ) * rx * 0.5d0
DUDY = ( Un1 - Us1 ) * ry * 0.5d0
DUDZ = ( Ut1 - Ub1 ) * rz * 0.5d0

DVDX = ( Ve1 - Vw1 ) * rx * 0.5d0
DVDY = ( Vn1 - Vs1 ) * ry * 0.5d0
DVDZ = ( Vt1 - Vb1 ) * rz * 0.5d0

DWDX = ( We1 - Ww1 ) * rx * 0.5d0
DWDY = ( We1 - Ws1 ) * ry * 0.5d0
DWDZ = ( We1 - Wb1 ) * rz * 0.5d0


!----------- Sij (strain rate tensor)  6 + 20 + 13 = 28
S11 = DUDX
S12 = (DVDX + DUDY) * 0.5d0
S13 = (DWDX + DUDZ) * 0.5d0

S21 = S12
S22 = DVDY
S23 = (DWDY + DVDZ) * 0.5d0

S31 = S13
S32 = S23
S33 = DWDZ

SSS = DSQRT( 2.0d0 * (S11*S11 + S22*S22 + S33*S33) &
           + 4.0d0 * (S12*S12 + S23*S23 + S31*S31) )


!----------- Wij (vorticity tensor)  12 + 20 + 12 = 34
W11 = 0.0d0
W12 = (DVDX - DUDY) * 0.5d0
W13 = (DWDX - DUDZ) * 0.5d0

W21 = (DUDY - DVDX) * 0.5d0
W22 = 0.0d0
W23 = (DWDY - DVDZ) * 0.5d0

W31 = (DUDZ - DWDX) * 0.5d0
W32 = (DVDZ - DWDY) * 0.5d0
W33 = 0.0d0

WWW = DSQRT( 2.0d0 * (W12*W12 + W13*W13 + W21*W21 &
                    + W23*W23 + W31*W31 + W32*W32 ) )



!----------- Smagorinsky model
if (imodel == 1) then
yc     = (j - 0.5d0)*dy
min_h  = min(yc, dy*jx - yc)
DUDY_w = v(i, jx, k, 1) * (ry*2.0d0)
tauw   = (nu * rho) * abs(DUDY_w)
utau   = sqrt(tauw/rho)
yp     = min_h * utau * r_nu
up     = Up0 / utau
fs     = 1.0d0 - exp(-yp/26.0d0)
nut = (Cs * fs * dy) * (Cs * fs * dy) * SSS
flop = flop + 83.0d0
end if


!----------- CSM model
if (imodel == 2) then
Q_csm   = (WWW * WWW - SSS * SSS) * 0.25d0
E_csm   = (WWW * WWW + SSS * SSS) * 0.25d0
Fcs     = Q_csm / (E_csm + EPS)

nut = (1.0d0/22.0d0) * abs(Fcs) * sqrt(abs(Fcs)) * (1.0d0 - Fcs) * (dx*dy*dz)**(2.0d0/3.0d0) * SSS
flop = flop + 45.0d0
end if


!----------- WALE model
if (imodel == 3) then

S11d = 0.5d0*(DUDX*DUDX + DUDX*DUDX) - (1.0d0/3.0d0)*(DUDX*DUDX + DVDY*DVDY + DWDZ*DWDZ)
S12d = 0.5d0*(DUDY*DUDY + DVDX*DVDX)
S13d = 0.5d0*(DUDZ*DUDZ + DWDX*DWDX)
S21d = S12d
S22d = 0.5d0*(DVDY*DVDY + DVDY*DVDY) - (1.0d0/3.0d0)*(DUDX*DUDX + DVDY*DVDY + DWDZ*DWDZ)
S23d = 0.5d0*(DVDZ*DVDZ + DWDY*DWDY)
S31d = S13d
S32d = S23d
S33d = 0.5d0*(DWDZ*DWDZ + DWDZ*DWDZ) - (1.0d0/3.0d0)*(DUDX*DUDX + DVDY*DVDY + DWDZ*DWDZ)

Sijd2 =  S11d*S11d + S21d*S21d + S31d*S31d &
       + S12d*S12d + S22d*S22d + S32d*S32d &
       + S13d*S13d + S23d*S23d + S33d*S33d

nut = (Cw * (dx*dy*dz)**(1.0d0/3.0d0))* (Cw * (dx*dy*dz)**(1.0d0/3.0d0)) * (Sijd2)**(3.0d0/2.0d0)  &
    / ( (SSS * SSS)**(5.0d0/2.0d0) + (Sijd2)**(5.0d0/4.0d0) )
flop = flop + 213.0d0
end if
!===============

!----------- DNS
if (imodel == 0) then
  nut = 0.0
end if


! 粘性項の計算　セル界面の剪断力を計算し，必要に応じて置換する   23*3 = 69
EX =   ( Ue1 - Up0 ) * c_e1 * rx2 &
     - ( Up0 - Uw1 ) * c_w1 * rx2 &
     + ( Un1 - Up0 ) * c_n1 * ry2 &
     - ( Up0 - Us1 ) * c_s1 * ry2 &
     + ( Ut1 - Up0 ) * c_t1 * rz2 &
     - ( Up0 - Ub1 ) * c_b1 * rz2

EY =   ( Ve1 - Vp0 ) * c_e1 * rx2 &
     - ( Vp0 - Vw1 ) * c_w1 * rx2 &
     + ( Vn1 - Vp0 ) * c_n1 * ry2 &
     - ( Vp0 - Vs1 ) * c_s1 * ry2 &
     + ( Vt1 - Vp0 ) * c_t1 * rz2 &
     - ( Vp0 - Vb1 ) * c_b1 * rz2

EZ =   ( We1 - Wp0 ) * c_e1 * rx2 &
     - ( Wp0 - Ww1 ) * c_w1 * rx2 &
     + ( Wn1 - Wp0 ) * c_n1 * ry2 &
     - ( Wp0 - Ws1 ) * c_s1 * ry2 &
     + ( Wt1 - Wp0 ) * c_t1 * rz2 &
     - ( Wp0 - Wb1 ) * c_b1 * rz2


! 対流項と粘性項の和 > 6*3 + 1 = 19 flops
nut = nut * r_nu
wv(i, j, k, 1) = -cnv_u * rx + EX * ( 1.0d0 + nut ) * vcs
wv(i, j, k, 2) = -cnv_v * ry + EY * ( 1.0d0 + nut ) * vcs
wv(i, j, k, 3) = -cnv_w * rz + EZ * ( 1.0d0 + nut ) * vcs

! 24 + 3 + 3 * (14 + 78 * 3 + 12) + 27 + 28 + 34 + 69 + 19 = 984
flop = flop + 984.0d0

end do
end do
end do

!$OMP END DO
!$OMP END PARALLEL


return
end subroutine pvec_muscl_les


!> ********************************************************************
!! @brief 次ステップのセルセンター，フェイスの速度と発散値を更新
!! @param [out] v    n+1時刻のセルセンター速度ベクトル
!! @param [out] vf   n+1時刻のセルフェイス速度ベクトル
!! @param [out] div  div {u^{n+1}}
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル長
!! @param [in]  dt   時間積分幅
!! @param [in]  dh   格子幅
!! @param [in]  vc   セルセンター疑似速度ベクトル
!! @param [in]  p    圧力
!! @param [in]  bp   BCindex P
!! @param [in]  bcd  BCindex B
!! @param [out] flop 浮動小数点演算数
!<
subroutine update_vec (v, vf, div, sz, g, dt, dh, vc, p, bp, bcd, flop)
implicit none
include 'ffv_f_params.h'
integer                                                   ::  i, j, k, ix, jx, kx, g, bpx, bdx
integer, dimension(3)                                     ::  sz
double precision                                          ::  flop
real                                                      ::  dt, actv, rx, ry, rz
real                                                      ::  pc, px, py, pz, pxw, pxe, pys, pyn, pzb, pzt
real                                                      ::  Ue0, Uw0, Vn0, Vs0, Wt0, Wb0, Up0, Vp0, Wp0
real                                                      ::  Ue, Uw, Vn, Vs, Wt, Wb
real                                                      ::  Uef, Uwf, Vnf, Vsf, Wtf, Wbf
real                                                      ::  c1, c2, c3, c4, c5, c6
real                                                      ::  N_e, N_w, N_n, N_s, N_t, N_b
real                                                      ::  c_w, c_e, c_s, c_n, c_b, c_t
real                                                      ::  d_w, d_e, d_s, d_n, d_b, d_t
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v, vc, vf
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  div, p
integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp, bcd
real, dimension(3)                                        ::  dh

ix = sz(1)
jx = sz(2)
kx = sz(3)

rx = 1.0 / dh(1)
ry = 1.0 / dh(2)
rz = 1.0 / dh(3)


flop = flop + dble(ix)*dble(jx)*dble(kx)*84.0 + 24.0d0


!$OMP PARALLEL &
!$OMP PRIVATE(bpx, actv, bdx) &
!$OMP PRIVATE(c1, c2, c3, c4, c5, c6) &
!$OMP PRIVATE(N_e, N_w, N_n, N_s, N_t, N_b) &
!$OMP PRIVATE(c_w, c_e, c_s, c_n, c_b, c_t) &
!$OMP PRIVATE(d_w, d_e, d_s, d_n, d_b, d_t) &
!$OMP PRIVATE(Ue0, Uw0, Vn0, Vs0, Wt0, Wb0, Up0, Vp0, Wp0) &
!$OMP PRIVATE(Ue, Uw, Vn, Vs, Wt, Wb) &
!$OMP PRIVATE(Uef, Uwf, Vnf, Vsf, Wtf, Wbf) &
!$OMP PRIVATE(pc, px, py, pz, pxw, pxe, pys, pyn, pzb, pzt) &
!$OMP FIRSTPRIVATE(ix, jx, kx, dt, rx, ry, rz)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do j=1,jx
do i=1,ix
bpx = bp(i,j,k)
bdx = bcd(i,j,k)
actv = real(ibits(bdx, State,  1))


! Neumann条件のとき，0.0
! 物体があればノイマン条件なので，セルフェイスマスクとしても利用
N_w = real(ibits(bpx, bc_n_W, 1))  ! w
N_e = real(ibits(bpx, bc_n_E, 1))  ! e
N_s = real(ibits(bpx, bc_n_S, 1))  ! s
N_n = real(ibits(bpx, bc_n_N, 1))  ! n
N_b = real(ibits(bpx, bc_n_B, 1))  ! b
N_t = real(ibits(bpx, bc_n_T, 1))  ! t

! \phi^D \phi^N
c_w = real(ibits(bpx, bc_ndag_W, 1))  ! w
c_e = real(ibits(bpx, bc_ndag_E, 1))  ! e
c_s = real(ibits(bpx, bc_ndag_S, 1))  ! s
c_n = real(ibits(bpx, bc_ndag_N, 1))  ! n
c_b = real(ibits(bpx, bc_ndag_B, 1))  ! b
c_t = real(ibits(bpx, bc_ndag_T, 1))  ! t

! (1 - \phi^D) \phi^N
d_w = real(ibits(bpx, bc_dn_W, 1))
d_e = real(ibits(bpx, bc_dn_E, 1))
d_s = real(ibits(bpx, bc_dn_S, 1))
d_n = real(ibits(bpx, bc_dn_N, 1))
d_b = real(ibits(bpx, bc_dn_B, 1))
d_t = real(ibits(bpx, bc_dn_T, 1))


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

Uw = 0.5 * ( Up0 + Uw0 ) * N_w ! 18 flop
Ue = 0.5 * ( Up0 + Ue0 ) * N_e
Vs = 0.5 * ( Vp0 + Vs0 ) * N_s
Vn = 0.5 * ( Vp0 + Vn0 ) * N_n
Wb = 0.5 * ( Wp0 + Wb0 ) * N_b
Wt = 0.5 * ( Wp0 + Wt0 ) * N_t


! 各面のVBCフラグの有無 => flux mask
! ibits() = 1(Normal) / 0(VBC)
c1 = real( ibits(bdx, bc_d_W, 1) )
c2 = real( ibits(bdx, bc_d_E, 1) )
c3 = real( ibits(bdx, bc_d_S, 1) )
c4 = real( ibits(bdx, bc_d_N, 1) )
c5 = real( ibits(bdx, bc_d_B, 1) )
c6 = real( ibits(bdx, bc_d_T, 1) )


! 圧力勾配 24flop >> DirichletとNeumannの値を0としている
pc  = p(i, j, k)
pxw = rx * (-p(i-1,j  ,k  )*c_w + (c_w + 2.0*d_w) * pc )
pxe = rx * ( p(i+1,j  ,k  )*c_e - (c_e + 2.0*d_e) * pc )
pys = ry * (-p(i  ,j-1,k  )*c_s + (c_s + 2.0*d_s) * pc )
pyn = ry * ( p(i  ,j+1,k  )*c_n - (c_n + 2.0*d_n) * pc )
pzb = rz * (-p(i  ,j  ,k-1)*c_b + (c_b + 2.0*d_b) * pc )
pzt = rz * ( p(i  ,j  ,k+1)*c_t - (c_t + 2.0*d_t) * pc )
px = 0.5 * (pxe + pxw)
py = 0.5 * (pyn + pys)
pz = 0.5 * (pzt + pzb)

! セルフェイス VBCの寄与と壁面の影響は除外 24flop
Uwf = (Uw - dt * pxw) * c1 * N_w
Uef = (Ue - dt * pxe) * c2 * N_e
Vsf = (Vs - dt * pys) * c3 * N_s
Vnf = (Vn - dt * pyn) * c4 * N_n
Wbf = (Wb - dt * pzb) * c5 * N_b
Wtf = (Wt - dt * pzt) * c6 * N_t

! i=1...ix >> vfは0...ixの範囲をカバーするので，通信不要
vf(i-1,j  ,k  ,1) = Uwf
vf(i  ,j  ,k  ,1) = Uef
vf(i  ,j-1,k  ,2) = Vsf
vf(i  ,j  ,k  ,2) = Vnf
vf(i  ,j  ,k-1,3) = Wbf
vf(i  ,j  ,k  ,3) = Wtf

div(i,j,k) = ((Uef - Uwf) * rx + (Vnf - Vsf) * ry + (Wtf - Wbf) * rz) * actv ! 9flop

! セルセンタの速度更新 9flop
v(i,j,k,1) = ( Up0 - dt * px ) * actv
v(i,j,k,2) = ( Vp0 - dt * py ) * actv
v(i,j,k,3) = ( Wp0 - dt * pz ) * actv

end do
end do
end do
!$OMP END DO

!$OMP END PARALLEL

return
end subroutine update_vec



!> ********************************************************************
!! @brief 次ステップのセルセンターの速度を更新
!! @param [out] v        n+1時刻の速度ベクトル
!! @param [out] vf       n+1時刻のセルフェイス速度ベクトル
!! @param [out] div      \sum {u^{n+1}}
!! @param [in]  sz       配列長
!! @param [in]  g        ガイドセル長
!! @param [in]  dt       時間積分幅
!! @param [in]  dh       格子幅
!! @param [in]  vc       セルセンター疑似速度ベクトル
!! @param [in]  p        圧力
!! @param [in]  bp       BCindex P
!! @param [in]  bid      Cut ID
!! @param [in]  bcd      BCindex B
!! @param [out] flop     浮動小数点演算数
!! @param [in   c_scheme 対流項スキーム
!<
subroutine update_vec4 (v, vf, div, sz, g, dt, dh, vc, p, bp, bid, bcd, flop, c_scheme)
implicit none
include 'ffv_f_params.h'
integer                                                   ::  i, j, k, ix, jx, kx, g, bpx, c_scheme, bdx
integer, dimension(3)                                     ::  sz
double precision                                          ::  flop
real                                                      ::  c_w1, c_e1, c_s1, c_n1, c_b1, c_t1
real                                                      ::  c_w2, c_e2, c_s2, c_n2, c_b2, c_t2
real                                                      ::  dt, actv, ss, c1, c2, sw1, sw2
real                                                      ::  N_e, N_w, N_n, N_s, N_t, N_b
real                                                      ::  D_e, D_w, D_n, D_s, D_t, D_b
real                                                      ::  N_e2, N_w2, N_n2, N_s2, N_t2, N_b2
real                                                      ::  D_e2, D_w2, D_n2, D_s2, D_t2, D_b2
real                                                      ::  Ue2, Uw2, Vn2, Vs2, Wt2, Wb2
real                                                      ::  Ue1, Uw1, Vn1, Vs1, Wt1, Wb1
real                                                      ::  Up0, Vp0, Wp0, drx, dry, drz, rx, ry, rz
real                                                      ::  Ue, Uw, Vn, Vs, Wt, Wb
real                                                      ::  pc, pe1, pe2, pw1, pw2, pn1, pn2, ps1, ps2, pt1, pt2, pb1, pb2
real                                                      ::  dp_e4, dp_w4, dp_n4, dp_s4, dp_t4, dp_b4
real                                                      ::  dp_e2, dp_w2, dp_n2, dp_s2, dp_t2, dp_b2
real                                                      ::  dpx, dpy, dpz, dp_e, dp_w, dp_n, dp_s, dp_t, dp_b
real                                                      ::  Uef, Uwf, Vnf, Vsf, Wtf, Wbf
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v, vc, vf
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  div, p
integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bp, bid, bcd
real, dimension(3)                                        ::  dh


c1 = 9.0/8.0
c2 = 1.0/24.0

if ( c_scheme == 2 ) then      ! 2nd order accuracy
ss = 0.0
else if ( c_scheme == 4 ) then ! 4th order accuracy
ss = 1.0
else
write(*, *) 'out of scheme selection'
stop
end if

ix = sz(1)
jx = sz(2)
kx = sz(3)

rx = 1.0 / dh(1)
ry = 1.0 / dh(2)
rz = 1.0 / dh(3)

drx = dt * rx
dry = dt * ry
drz = dt * rz

flop = flop + dble(ix)*dble(jx)*dble(kx)*159.0 + 27.0d0

!$OMP PARALLEL &
!$OMP PRIVATE(bpx, actv, bdx) &
!$OMP PRIVATE(c_w1, c_e1, c_s1, c_n1, c_b1, c_t1) &
!$OMP PRIVATE(c_w2, c_e2, c_s2, c_n2, c_b2, c_t2) &
!$OMP PRIVATE(N_e, N_w, N_n, N_s, N_t, N_b) &
!$OMP PRIVATE(D_e, D_w, D_n, D_s, D_t, D_b) &
!$OMP PRIVATE(N_e2, N_w2, N_n2, N_s2, N_t2, N_b2) &
!$OMP PRIVATE(D_e2, D_w2, D_n2, D_s2, D_t2, D_b2) &
!$OMP PRIVATE(Uw2, Uw1, Up0, Ue1, Ue2) &
!$OMP PRIVATE(Vs2, Vs1, Vp0, Vn1, Vn2) &
!$OMP PRIVATE(Wb2, Wb1, Wp0, Wt1, Wt2) &
!$OMP PRIVATE(Ue, Uw, Vn, Vs, Wt, Wb) &
!$OMP PRIVATE(sw1, sw2) &
!$OMP PRIVATE(pc, pe1, pe2, pw1, pw2, pn1, pn2, ps1, ps2, pt1, pt2, pb1, pb2) &
!$OMP PRIVATE(dp_e4, dp_w4, dp_n4, dp_s4, dp_t4, dp_b4) &
!$OMP PRIVATE(dp_e2, dp_w2, dp_n2, dp_s2, dp_t2, dp_b2) &
!$OMP PRIVATE(dpx, dpy, dpz, dp_e, dp_w, dp_n, dp_s, dp_t, dp_b) &
!$OMP PRIVATE(Uef, Uwf, Vnf, Vsf, Wtf, Wbf) &
!$OMP FIRSTPRIVATE(ix, jx, kx, ss, c1, c2, drx, dry, drz, rx, ry, rz)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do j=1,jx
do i=1,ix
bpx = bp(i,j,k)
bdx = bcd(i,j,k)
actv = real(ibits(bdx, State,  1))


! Neumann条件のとき，0.0 > 6 flop
! 物体があればノイマン条件なので，セルフェイスマスクとしても利用
N_w = real(ibits(bpx, bc_n_W, 1))  ! w
N_e = real(ibits(bpx, bc_n_E, 1))  ! e
N_s = real(ibits(bpx, bc_n_S, 1))  ! s
N_n = real(ibits(bpx, bc_n_N, 1))  ! n
N_b = real(ibits(bpx, bc_n_B, 1))  ! b
N_t = real(ibits(bpx, bc_n_T, 1))  ! t

N_w2 = real( ibits(bp(i-1, j,   k  ), bc_n_W, 1) )
N_e2 = real( ibits(bp(i+1, j,   k  ), bc_n_E, 1) )
N_s2 = real( ibits(bp(i,   j-1, k  ), bc_n_S, 1) )
N_n2 = real( ibits(bp(i,   j+1, k  ), bc_n_N, 1) )
N_b2 = real( ibits(bp(i,   j,   k-1), bc_n_B, 1) )
N_t2 = real( ibits(bp(i,   j,   k+1), bc_n_T, 1) )


! Dirichlet(p=0)条件のとき 0.0
D_w = real(ibits(bpx, bc_d_W, 1))  ! w
D_e = real(ibits(bpx, bc_d_E, 1))  ! e
D_s = real(ibits(bpx, bc_d_S, 1))  ! s
D_n = real(ibits(bpx, bc_d_N, 1))  ! n
D_b = real(ibits(bpx, bc_d_B, 1))  ! b
D_t = real(ibits(bpx, bc_d_T, 1))  ! t

D_w2 = real( ibits(bp(i-1, j,   k  ), bc_d_W, 1) )
D_e2 = real( ibits(bp(i+1, j,   k  ), bc_d_E, 1) )
D_s2 = real( ibits(bp(i,   j-1, k  ), bc_d_S, 1) )
D_n2 = real( ibits(bp(i,   j+1, k  ), bc_d_N, 1) )
D_b2 = real( ibits(bp(i,   j,   k-1), bc_d_B, 1) )
D_t2 = real( ibits(bp(i,   j,   k+1), bc_d_T, 1) )


! 疑似ベクトル
Uw2 = vc(i-2,j  ,k  , 1)
Uw1 = vc(i-1,j  ,k  , 1)
Up0 = vc(i  ,j  ,k  , 1)
Ue1 = vc(i+1,j  ,k  , 1)
Ue2 = vc(i+2,j  ,k  , 1)

Vs2 = vc(i  ,j-2,k  , 2)
Vs1 = vc(i  ,j-1,k  , 2)
Vp0 = vc(i  ,j  ,k  , 2)
Vn1 = vc(i  ,j+1,k  , 2)
Vn2 = vc(i  ,j+2,k  , 2)

Wb2 = vc(i  ,j  ,k-2, 3)
Wb1 = vc(i  ,j  ,k-1, 3)
Wp0 = vc(i  ,j  ,k  , 3)
Wt1 = vc(i  ,j  ,k+1, 3)
Wt2 = vc(i  ,j  ,k+2, 3)

Uw = 0.5*( Up0 + Uw1 )*N_w ! 18 flop
Ue = 0.5*( Up0 + Ue1 )*N_e
Vs = 0.5*( Vp0 + Vs1 )*N_s
Vn = 0.5*( Vp0 + Vn1 )*N_n
Wb = 0.5*( Wp0 + Wb1 )*N_b
Wt = 0.5*( Wp0 + Wt1 )*N_t


! 各面のVBCフラグの有無 => flux mask
! ibits() = 1(Normal) / 0(VBC)
c_w1 = real( ibits(bdx, bc_d_W, 1) )
c_e1 = real( ibits(bdx, bc_d_E, 1) )
c_s1 = real( ibits(bdx, bc_d_S, 1) )
c_n1 = real( ibits(bdx, bc_d_N, 1) )
c_b1 = real( ibits(bdx, bc_d_B, 1) )
c_t1 = real( ibits(bdx, bc_d_T, 1) )


! ステンシルの参照先がvspec, outflowである場合のスキームの破綻を回避，１次精度におとすSW
c_w2 = real( ibits(bcd(i-1, j  , k  ), bc_d_W, 1) )
c_e2 = real( ibits(bcd(i+1, j  , k  ), bc_d_E, 1) )
c_s2 = real( ibits(bcd(i  , j-1, k  ), bc_d_S, 1) )
c_n2 = real( ibits(bcd(i  , j+1, k  ), bc_d_N, 1) )
c_b2 = real( ibits(bcd(i  , j  , k-1), bc_d_B, 1) )
c_t2 = real( ibits(bcd(i  , j  , k+1), bc_d_T, 1) )


!---- target cell >> 24 flops
pc = p(i, j, k)

!---- x-direction
pe1 = p(i+1, j, k) * D_e  * N_e
pe2 = p(i+2, j, k) * D_e2 * N_e2
pw1 = p(i-1, j, k) * D_w  * N_w
pw2 = p(i-2, j, k) * D_w2 * N_w2

!---- y-direction
pn1 = p(i, j+1, k) * D_n  * N_n
pn2 = p(i, j+2, k) * D_n2 * N_n2
ps1 = p(i, j-1, k) * D_s  * N_s
ps2 = p(i, j-2, k) * D_s2 * N_s2

!---- z-direction
pt1 = p(i, j, k+1) * D_t  * N_t
pt2 = p(i, j, k+2) * D_t2 * N_t2
pb1 = p(i, j, k-1) * D_b  * N_b
pb2 = p(i, j, k-2) * D_b2 * N_b2

!---- pressure difference for each direction >> 26*3 flops
! N_e = N_w = 0: wall face, 1: fluid face
sw1   = N_e * c_e2 * ss
sw2   = N_w * c_w2 * ss
dp_e2 = pe1 - pc
dp_w2 = pc  - pw1
dp_e4 = c1 * dp_e2 - c2 * (pe2 - pw1)
dp_w4 = c1 * dp_w2 - c2 * (pe1 - pw2)
dp_e  = dp_e4 * sw1 + dp_e2 * (1.0 - sw1)
dp_w  = dp_w4 * sw2 + dp_w2 * (1.0 - sw2)
dpx   = 0.5d0 * ( dp_e * c_e1 + dp_w * c_w1 )

! N_n = N_s = 0: wall face, 1: fluid face
sw1   = N_n * c_n2 * ss
sw2   = N_s * c_s2 * ss
dp_n2 = pn1 - pc
dp_s2 = pc  - ps1
dp_n4 = c1 * dp_n2 - c2 * (pn2 - ps1)
dp_s4 = c1 * dp_s2 - c2 * (pn1 - ps2)
dp_n  = dp_n4 * sw1 + dp_n2 * (1.0 - sw1)
dp_s  = dp_s4 * sw2 + dp_s2 * (1.0 - sw2)
dpy   = 0.5d0 * ( dp_n * c_n1 + dp_s * c_s1 )

! N_t = N_b = 0: wall face, 1: fluid face
sw1   = N_t * c_t2 * ss
sw2   = N_b * c_b2 * ss
dp_t2 = pt1 - pc
dp_b2 = pc  - pb1
dp_t4 = c1 * dp_t2 - c2 * (pt2 - pb1)
dp_b4 = c1 * dp_b2 - c2 * (pt1 - pb2)
dp_t  = dp_t4 * sw1 + dp_t2 * (1.0 - sw1)
dp_b  = dp_b4 * sw2 + dp_b2 * (1.0 - sw2)
dpz   = 0.5d0 * ( dp_t * c_t1 + dp_b * c_b1 )


! セルフェイス VBCの寄与と壁面の影響は除外 24flop
Uwf = (Uw - drx * dp_w) * c_w1 * N_w
Uef = (Ue - drx * dp_e) * c_e1 * N_e
Vsf = (Vs - dry * dp_s) * c_s1 * N_s
Vnf = (Vn - dry * dp_n) * c_n1 * N_n
Wbf = (Wb - drz * dp_b) * c_b1 * N_b
Wtf = (Wt - drz * dp_t) * c_t1 * N_t

! i=1...ix >> vfは0...ixの範囲をカバーするので，通信不要
vf(i-1,j  ,k  ,1) = Uwf
vf(i  ,j  ,k  ,1) = Uef
vf(i  ,j-1,k  ,2) = Vsf
vf(i  ,j  ,k  ,2) = Vnf
vf(i  ,j  ,k-1,3) = Wbf
vf(i  ,j  ,k  ,3) = Wtf

div(i,j,k) = ( (Uef - Uwf)*rx +(Vnf - Vsf)*ry + (Wtf - Wbf)*rz ) * actv ! 6flop

! セルセンタの速度更新 9flop
v(i, j, k, 1) = ( Up0 - dpx * drx ) * actv
v(i, j, k, 2) = ( Vp0 - dpy * dry ) * actv
v(i, j, k, 3) = ( Wp0 - dpz * drz ) * actv

end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL

return
end subroutine update_vec4


!> ********************************************************************
!! @brief 速度の発散に使う div{u^*} を計算する
!! @param [out]    div  速度の和
!! @param [in]     sz   配列長
!! @param [in]     g    ガイドセル長
!! @param [in]     dh   格子幅
!! @param [in]     vc   セルセンター疑似ベクトル
!! @param [in]     bid  Cut ID
!! @param [in]     bcd  BCindex B
!! @param [in,out] flop 浮動小数点演算数
!<
subroutine divergence_cc (div, sz, g, dh, vc, bid, bcd, flop)
implicit none
include 'ffv_f_params.h'
integer                                                   ::  i, j, k, ix, jx, kx, g, bix, bdx
integer, dimension(3)                                     ::  sz
double precision                                          ::  flop
real                                                      ::  Ue, Uw, Vn, Vs, Wt, Wb, actv, rx, ry, rz
real                                                      ::  Ue0, Uw0, Up0, Vn0, Vs0, Vp0, Wb0, Wt0, Wp0
real                                                      ::  c_e, c_w, c_n, c_s, c_t, c_b
real                                                      ::  b_w, b_e, b_s, b_n, b_b, b_t
real, dimension(3)                                        ::  dh
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  vc
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  div
integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bid, bcd

ix = sz(1)
jx = sz(2)
kx = sz(3)

rx = 1.0/dh(1)
ry = 1.0/dh(2)
rz = 1.0/dh(3)

flop  = flop + dble(ix)*dble(jx)*dble(kx)*33.0d0 + 24.0d0

!$OMP PARALLEL &
!$OMP PRIVATE(actv, bix, bdx) &
!$OMP PRIVATE(b_w, b_e, b_s, b_n, b_b, b_t) &
!$OMP PRIVATE(Ue, Uw, Vn, Vs, Wt, Wb) &
!$OMP PRIVATE(Ue0, Uw0, Up0, Vn0, Vs0, Vp0, Wb0, Wt0, Wp0) &
!$OMP PRIVATE(c_e, c_w, c_n, c_s, c_t, c_b) &
!$OMP FIRSTPRIVATE(ix, jx, kx, rx, ry, rz)

!$OMP DO SCHEDULE(static) COLLAPSE(2)

do k=1,kx
do j=1,jx
do i=1,ix
bix = bid(i,j,k)
bdx = bcd(i,j,k)

actv= real(ibits(bdx, State, 1))

! 各セルセンター位置の変数ロード
Uw0 = vc(i-1,j  ,k  , 1)
Up0 = vc(i  ,j  ,k  , 1)
Ue0 = vc(i+1,j  ,k  , 1)

Vs0 = vc(i  ,j-1,k  , 2)
Vp0 = vc(i  ,j  ,k  , 2)
Vn0 = vc(i  ,j+1,k  , 2)

Wb0 = vc(i  ,j  ,k-1, 3)
Wp0 = vc(i  ,j  ,k  , 3)
Wt0 = vc(i  ,j  ,k+1, 3)

! 0-solid / 1-fluid
b_w = 1.0
b_e = 1.0
b_s = 1.0
b_n = 1.0
b_b = 1.0
b_t = 1.0
if ( ibits(bix, bc_face_W, bitw_5) /= 0 ) b_w = 0.0
if ( ibits(bix, bc_face_E, bitw_5) /= 0 ) b_e = 0.0
if ( ibits(bix, bc_face_S, bitw_5) /= 0 ) b_s = 0.0
if ( ibits(bix, bc_face_N, bitw_5) /= 0 ) b_n = 0.0
if ( ibits(bix, bc_face_B, bitw_5) /= 0 ) b_b = 0.0
if ( ibits(bix, bc_face_T, bitw_5) /= 0 ) b_t = 0.0

Uw = 0.5*( Up0 + Uw0 )*b_w ! 18 flops
Ue = 0.5*( Up0 + Ue0 )*b_e
Vs = 0.5*( Vp0 + Vs0 )*b_s
Vn = 0.5*( Vp0 + Vn0 )*b_n
Wb = 0.5*( Wp0 + Wb0 )*b_b
Wt = 0.5*( Wp0 + Wt0 )*b_t


! 各面のVBCフラグの有無 => flux mask
! ibits() = 1(Normal) / 0(VBC)
c_w = real( ibits(bdx, bc_d_W, 1) )
c_e = real( ibits(bdx, bc_d_E, 1) )
c_s = real( ibits(bdx, bc_d_S, 1) )
c_n = real( ibits(bdx, bc_d_N, 1) )
c_b = real( ibits(bdx, bc_d_B, 1) )
c_t = real( ibits(bdx, bc_d_T, 1) )


! VBC面の影響をフラグで無効化 >> OBC_SPEC_VEL, OBC_WALL  15flops
div(i,j,k) = ( rx * ( Ue * c_e - Uw * c_w ) &
             + ry * ( Vn * c_n - Vs * c_s ) &
             + rz * ( Wt * c_t - Wb * c_b ) ) * actv
end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL

return
end subroutine divergence_cc

!> ********************************************************************
!! @brief 疑似ベクトルの時間積分（Euler陽解法）
!! @param [in,out] vc   対流項と粘性項の和 > 疑似ベクトル
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

!$OMP DO SCHEDULE(static) COLLAPSE(2)
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
!! @note
!!    - vtmin, vtmax > vt_range(2)
!!    - ypmin, ypmax > yp_range(2)
!! @todo
!!    - 境界条件は必要か？
!! @note NOCHECK
!<
subroutine eddy_viscosity (vt, sz, g, dh, re, cs, v, bx, vt_range, yp_range)
implicit none
include 'ffv_f_params.h'
integer                                                     ::  i, j, k, ix, jx, kx, g, m
integer, dimension(3)                                       ::  sz
real, dimension(2)                                          ::  vt_range, yp_range
real                                                        ::  re, cs, yp, Ut, Ut0
real                                                        ::  DUDX, DUDY, DUDZ
real                                                        ::  DVDX, DVDY, DVDZ
real                                                        ::  DWDX, DWDY, DWDZ
real                                                        ::  d1, d2, ddd, c1, c2, c3, delta
real                                                        ::  fs, aaa, Vmag, dis, tw, up1
real                                                        ::  u1, u2, u3
real                                                        ::  dx, dy, dz, rx, ry, rz
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)      ::  vt
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)   ::  v
integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)   ::  bx
real, dimension(3)                                          ::  dh

ix = sz(1)
jx = sz(2)
kx = sz(3)

dx = dh(1)
dy = dh(2)
dz = dh(3)

rx = 0.5 / dx
ry = 0.5 / dy
rz = 0.5 / dz

vt_range(1) =  1.0e6
vt_range(2) = -1.0e6
yp_range(1) =  1.0e6
yp_range(2) = -1.0e6

Delta = (dx*dy*dz)**(1.0/3.0)

    do k=1,kx
    do j=1,jx
    do i=1,ix
      DUDX=(v(i+1,j  ,k  , 1)-v(i-1,j  ,k  , 1)) * rx
      DUDY=(v(i  ,j+1,k  , 1)-v(i  ,j-1,k  , 1)) * ry
      DUDZ=(v(i  ,j  ,k+1, 1)-v(i  ,j  ,k-1, 1)) * rz

      DVDX=(v(i+1,j  ,k  , 2)-v(i-1,j  ,k  , 2)) * rx
      DVDY=(v(i  ,j+1,k  , 2)-v(i  ,j-1,k  , 2)) * ry
      DVDZ=(v(i  ,j  ,k+1, 2)-v(i  ,j  ,k-1, 2)) * rz

      DWDX=(v(i+1,j  ,k  , 3)-v(i-1,j  ,k  , 3)) * rx
      DWDY=(v(i  ,j+1,k  , 3)-v(i  ,j-1,k  , 3)) * ry
      DWDZ=(v(i  ,j  ,k+1, 3)-v(i  ,j  ,k-1, 3)) * rz

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
      u1 = v(i,j,k,1)
      u2 = v(i,j,k,2)
      u3 = v(i,j,k,3)
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
!! @param [out] vc  疑似ベクトル
!! @param [in]  sz  配列長
!! @param [in]  g   ガイドセル長
!! @param [in]  dt  時間積分幅
!! @param [in]  v   速度ベクトル（n-step, collocated）
!! @param [in]  ab  前ステップの対流項（＋粘性項）の計算値
!! @param [in]  bd  BCindex B
!! @param[out] flop
!<
subroutine ab2 (vc, sz, g, dt, v, ab, bd, flop)
implicit none
include 'ffv_f_params.h'
integer                                                   ::  i, j, k, ix, jx, kx, g
integer, dimension(3)                                     ::  sz
double precision                                          ::  flop
real                                                      ::  actv, dt, ab_u, ab_v, ab_w
real                                                      ::  cf, ra
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  vc, v, ab
integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bd

ix = sz(1)
jx = sz(2)
kx = sz(3)

cf = 0.5 * dt

flop = flop + dble(ix*jx*kx) * 22.0d0


!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, cf) &
!$OMP PRIVATE(actv, ra, ab_u, ab_v, ab_w)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do j=1,jx
do i=1,ix
actv = real(ibits(bd(i,j,k), State, 1))
ra = 1.0 - actv

ab_u = ab(i,j,k,1)
ab_v = ab(i,j,k,2)
ab_w = ab(i,j,k,3)

ab(i,j,k,1) = vc(i,j,k,1)
ab(i,j,k,2) = vc(i,j,k,2)
ab(i,j,k,3) = vc(i,j,k,3)

vc(i,j,k,1) = ( v(i,j,k,1) + cf * ( 3.0 * vc(i,j,k,1) - ab_u ) ) * actv
vc(i,j,k,2) = ( v(i,j,k,2) + cf * ( 3.0 * vc(i,j,k,2) - ab_v ) ) * actv
vc(i,j,k,3) = ( v(i,j,k,3) + cf * ( 3.0 * vc(i,j,k,3) - ab_w ) ) * actv
end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL

return
end subroutine ab2


!> ********************************************************************
!! @brief 対流項と粘性項の計算
!! @param [out] wv        疑似ベクトルの空間項 u \frac{\partial u}{\partial x}
!! @param [in]  sz        配列長
!! @param [in]  g         ガイドセル長
!! @param [in]  dh        格子幅
!! @param [in]  c_scheme  対流項スキームのモード（2-Central_2nd, 4-Central_4th）
!! @param [in]  rei       レイノルズ数の逆数
!! @param [in]  v         セルセンター速度ベクトル（n-step）
!! @param [in]  vf        セルフェイス速度ベクトル（n-step）
!! @param [in]  bid       Cut ID
!! @param [in]  bcd       BCindex B
!! @param [in]  vcs_coef  粘性項の係数（粘性項を計算しない場合には0.0）
!! @param [out] flop      浮動小数点演算数
!<
subroutine pvec_central (wv, sz, g, dh, c_scheme, rei, v, vf, bid, bcd, vcs_coef, flop)
implicit none
include 'ffv_f_params.h'
integer                                                   ::  i, j, k, ix, jx, kx, g, c_scheme, bix, bdx
integer, dimension(3)                                     ::  sz
double precision                                          ::  flop
real                                                      ::  b_e1, b_w1, b_n1, b_s1, b_t1, b_b1
real                                                      ::  b_e2, b_w2, b_n2, b_s2, b_t2, b_b2, b_p
real                                                      ::  Up0, Ue1, Ue2, Uw1, Uw2, Us1, Us2, Un1, Un2, Ub1, Ub2, Ut1, Ut2
real                                                      ::  Vp0, Ve1, Ve2, Vw1, Vw2, Vs1, Vs2, Vn1, Vn2, Vb1, Vb2, Vt1, Vt2
real                                                      ::  Wp0, We1, We2, Ww1, Ww2, Ws1, Ws2, Wn1, Wn2, Wb1, Wb2, Wt1, Wt2
real                                                      ::  UPe, UPw, VPn, VPs, WPt, WPb
real                                                      ::  vcs, vcs_coef, ss, sw1, sw2, c1, c2
real                                                      ::  c_e1, c_w1, c_n1, c_s1, c_t1, c_b1
real                                                      ::  c_e2, c_w2, c_n2, c_s2, c_t2, c_b2
real                                                      ::  cnv_u, cnv_v, cnv_w, EX, EY, EZ, rei
real                                                      ::  uq, vq, wq, rx, ry, rz, rx2, ry2 ,rz2
real                                                      ::  ufr2, ufl2, vfr2, vfl2, wfr2, wfl2
real                                                      ::  ufr4, ufl4, vfr4, vfl4, wfr4, wfl4
real                                                      ::  fu_r, fu_l, fv_r, fv_l, fw_r, fw_l
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v, wv, vf
real, dimension(3)                                        ::  dh
integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bid, bcd

ix = sz(1)
jx = sz(2)
kx = sz(3)

rx = 1.0/dh(1)
ry = 1.0/dh(2)
rz = 1.0/dh(3)

rx2= rei * rx
ry2= rei * ry
rz2= rei * rz

c1 = 9.0/8.0
c2 = 1.0/24.0

! vcs = 1.0 (Euler Explicit) / 0.5 (CN) / 0.0(No)
vcs = vcs_coef


if ( c_scheme == 2 ) then      !     2nd
ss = 0.0
else if ( c_scheme == 4 ) then !     4th
ss = 1.0
else
write(*,*) 'out of scheme selection'
stop
endif


! 24 + 3 + 3 * 106 + 21 + 9 = 375
flop = flop + dble(ix)*dble(jx)*dble(kx)*375.0d0 + 46.0d0


!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, rx, ry, rz, rx2, ry2 ,rz2, vcs, ss, c1, c2) &
!$OMP FIRSTPRIVATE(rei) &
!$OMP PRIVATE(cnv_u, cnv_v, cnv_w, bix, bdx, uq, vq, wq, sw1, sw2) &
!$OMP PRIVATE(Up0, Ue1, Ue2, Uw1, Uw2, Us1, Us2, Un1, Un2, Ub1, Ub2, Ut1, Ut2) &
!$OMP PRIVATE(Vp0, Ve1, Ve2, Vw1, Vw2, Vs1, Vs2, Vn1, Vn2, Vb1, Vb2, Vt1, Vt2) &
!$OMP PRIVATE(Wp0, We1, We2, Ww1, Ww2, Ws1, Ws2, Wn1, Wn2, Wb1, Wb2, Wt1, Wt2) &
!$OMP PRIVATE(b_e1, b_w1, b_n1, b_s1, b_t1, b_b1, b_e2, b_w2, b_n2, b_s2, b_t2, b_b2, b_p) &
!$OMP PRIVATE(c_e1, c_w1, c_n1, c_s1, c_t1, c_b1, c_e2, c_w2, c_n2, c_s2, c_t2, c_b2) &
!$OMP PRIVATE(UPe, UPw, VPn, VPs, WPt, WPb) &
!$OMP PRIVATE(ufr2, ufl2, vfr2, vfl2, wfr2, wfl2) &
!$OMP PRIVATE(ufr4, ufl4, vfr4, vfl4, wfr4, wfl4) &
!$OMP PRIVATE(fu_r, fu_l, fv_r, fv_l, fw_r, fw_l) &
!$OMP PRIVATE(EX, EY, EZ)

!$OMP DO SCHEDULE(static) COLLAPSE(2)

do k=1,kx
do j=1,jx
do i=1,ix
cnv_u = 0.0
cnv_v = 0.0
cnv_w = 0.0
EX = 0.0
EY = 0.0
EZ = 0.0

! 変数のロード
! 各軸方向5点の変数ロード
Ub2 = v(i  ,j  ,k-2, 1)
Ub1 = v(i  ,j  ,k-1, 1)
Us2 = v(i  ,j-2,k  , 1)
Us1 = v(i  ,j-1,k  , 1)
Uw2 = v(i-2,j  ,k  , 1)
Uw1 = v(i-1,j  ,k  , 1)
Up0 = v(i  ,j  ,k  , 1)
Ue1 = v(i+1,j  ,k  , 1)
Ue2 = v(i+2,j  ,k  , 1)
Un1 = v(i  ,j+1,k  , 1)
Un2 = v(i  ,j+2,k  , 1)
Ut1 = v(i  ,j  ,k+1, 1)
Ut2 = v(i  ,j  ,k+2, 1)

Vb2 = v(i  ,j  ,k-2, 2)
Vb1 = v(i  ,j  ,k-1, 2)
Vs2 = v(i  ,j-2,k  , 2)
Vs1 = v(i  ,j-1,k  , 2)
Vw2 = v(i-2,j  ,k  , 2)
Vw1 = v(i-1,j  ,k  , 2)
Vp0 = v(i  ,j  ,k  , 2)
Ve1 = v(i+1,j  ,k  , 2)
Ve2 = v(i+2,j  ,k  , 2)
Vn1 = v(i  ,j+1,k  , 2)
Vn2 = v(i  ,j+2,k  , 2)
Vt1 = v(i  ,j  ,k+1, 2)
Vt2 = v(i  ,j  ,k+2, 2)

Wb2 = v(i  ,j  ,k-2, 3)
Wb1 = v(i  ,j  ,k-1, 3)
Ws2 = v(i  ,j-2,k  , 3)
Ws1 = v(i  ,j-1,k  , 3)
Ww2 = v(i-2,j  ,k  , 3)
Ww1 = v(i-1,j  ,k  , 3)
Wp0 = v(i  ,j  ,k  , 3)
We1 = v(i+1,j  ,k  , 3)
We2 = v(i+2,j  ,k  , 3)
Wn1 = v(i  ,j+1,k  , 3)
Wn2 = v(i  ,j+2,k  , 3)
Wt1 = v(i  ,j  ,k+1, 3)
Wt2 = v(i  ,j  ,k+2, 3)

bix = bid(i,j,k)
bdx = bcd(i,j,k)

! (i,j,k)からみたセル状態 (0-solid / 1-fluid)
b_p = real(ibits(bdx, State, 1))

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
c_w1 = real( ibits(bdx, bc_d_W, 1) )
c_e1 = real( ibits(bdx, bc_d_E, 1) )
c_s1 = real( ibits(bdx, bc_d_S, 1) )
c_n1 = real( ibits(bdx, bc_d_N, 1) )
c_b1 = real( ibits(bdx, bc_d_B, 1) )
c_t1 = real( ibits(bdx, bc_d_T, 1) )


! ステンシルの参照先がvspec, outflowである場合のスキームの破綻を回避，１次精度におとすSW
c_w2 = real( ibits(bcd(i-1, j  , k  ), bc_d_W, 1) )
c_e2 = real( ibits(bcd(i+1, j  , k  ), bc_d_E, 1) )
c_s2 = real( ibits(bcd(i  , j-1, k  ), bc_d_S, 1) )
c_n2 = real( ibits(bcd(i  , j+1, k  ), bc_d_N, 1) )
c_b2 = real( ibits(bcd(i  , j  , k-1), bc_d_B, 1) )
c_t2 = real( ibits(bcd(i  , j  , k+1), bc_d_T, 1) )


! 界面速度（スタガード位置） 24 flops
UPe = vf(i  , j  , k  ,1)*b_e1
UPw = vf(i-1, j  , k  ,1)*b_w1
VPn = vf(i  , j  , k  ,2)*b_n1
VPs = vf(i  , j-1, k  ,2)*b_s1
WPt = vf(i  , j  , k  ,3)*b_t1
WPb = vf(i  , j  , k-1,3)*b_b1

! セルセンターからの壁面修正速度 3 flops
uq = - Up0
vq = - Vp0
wq = - Wp0

! X方向 --------------------------------------- >> 6 + 4 + 20*3 + 21 + 15 = 106

! 速度指定の場合に参照先として，固体内にテンポラリに与えた値を使う   6flops
if ( (b_e2 == 0.0)  ) then
  Ue2 = - v(i+1,j  ,k  , 1)
  Ve2 = - v(i+1,j  ,k  , 2)
  We2 = - v(i+1,j  ,k  , 3)
endif

if ( b_e1 == 0.0 ) then
  Ue1 = uq
  Ve1 = vq
  We1 = wq
endif

if ( b_w1 == 0.0 ) then
  Uw1 = uq
  Vw1 = vq
  Ww1 = wq
end if

if ( (b_w2 == 0.0)  ) then
  Uw2 = - v(i-1,j  ,k  , 1)
  Vw2 = - v(i-1,j  ,k  , 2)
  Ww2 = - v(i-1,j  ,k  , 3)
end if

! flux u \frac{\partial u}{\partial x}
! 壁面がある場合　b_e1, b_w1 (0-wall face / 1-fluid)
! vspec/outflowを参照する場合(c_e2, c_w2)は二次精度へおとす
! 二次精度の時には、ssで強制
sw1 = b_e1 * c_e2 * ss ! 4 flops
sw2 = b_w1 * c_w2 * ss

ufr2 = Ue1-Up0 ! 20 flops
ufl2 = Up0-Uw1
ufr4 = c1 * ufr2 - c2 * (Ue2-Uw1)
ufl4 = c1 * ufl2 - c2 * (Ue1-Uw2)
fu_r = ( ufr4 * sw1 + ufr2 * (1.0 - sw1) ) * rx
fu_l = ( ufl4 * sw2 + ufl2 * (1.0 - sw2) ) * rx

vfr2 = Ve1-Vp0
vfl2 = Vp0-Vw1
vfr4 = c1 * vfr2 - c2 * (Ve2-Vw1)
vfl4 = c1 * vfl2 - c2 * (Ve1-Vw2)
fv_r = ( vfr4 * sw1 + vfr2 * (1.0 - sw1) ) * rx
fv_l = ( vfl4 * sw2 + vfl2 * (1.0 - sw2) ) * rx

wfr2 = We1-Wp0
wfl2 = Wp0-Ww1
wfr4 = c1 * wfr2 - c2 * (We2-Ww1)
wfl4 = c1 * wfl2 - c2 * (We1-Ww2)
fw_r = ( wfr4 * sw1 + wfr2 * (1.0 - sw1) ) * rx
fw_l = ( wfl4 * sw2 + wfl2 * (1.0 - sw2) ) * rx


! 流束の加算　平均勾配　VBCでない面の寄与のみを評価する 21 flops
cnv_u = cnv_u + 0.5*(Upe * fu_r * c_e1 + Upw * fu_l * c_w1)
cnv_v = cnv_v + 0.5*(Upe * fv_r * c_e1 + Upw * fv_l * c_w1)
cnv_w = cnv_w + 0.5*(Upe * fw_r * c_e1 + Upw * fw_l * c_w1)

! 粘性項の加算 15 flops
EX = EX + (fu_r * c_e1 - fu_l * c_w1) * rx2
EY = EY + (fv_r * c_e1 - fv_l * c_w1) * rx2
EZ = EZ + (fw_r * c_e1 - fw_l * c_w1) * rx2


! Y方向 ---------------------------------------

if ( (b_n2 == 0.0)  ) then
  Un2 = - v(i  ,j+1,k  , 1)
  Vn2 = - v(i  ,j+1,k  , 2)
  Wn2 = - v(i  ,j+1,k  , 3)
endif

if ( b_n1 == 0.0 ) then
  Un1 = uq
  Vn1 = vq
  Wn1 = wq
endif

if ( b_s1 == 0.0 ) then
  Us1 = uq
  Vs1 = vq
  Ws1 = wq
endif

if ( (b_s2 == 0.0)  ) then
  Us2 = - v(i  ,j-1,k  , 1)
  Vs2 = - v(i  ,j-1,k  , 2)
  Ws2 = - v(i  ,j-1,k  , 3)
endif

! flux
sw1 = b_n1 * c_n2 * ss
sw2 = b_s1 * c_s2 * ss

ufr2 = Un1-Up0
ufl2 = Up0-Us1
ufr4 = c1 * ufr2 - c2 * (Un2-Us1)
ufl4 = c1 * ufl2 - c2 * (Un1-Us2)
fu_r = ( ufr4 * sw1 + ufr2 * (1.0 - sw1) ) * ry
fu_l = ( ufl4 * sw2 + ufl2 * (1.0 - sw2) ) * ry

vfr2 = Vn1-Vp0
vfl2 = Vp0-Vs1
vfr4 = c1 * vfr2 - c2 * (Vn2-Vs1)
vfl4 = c1 * vfl2 - c2 * (Vn1-Vs2)
fv_r = ( vfr4 * sw1 + vfr2 * (1.0 - sw1) ) * ry
fv_l = ( vfl4 * sw2 + vfl2 * (1.0 - sw2) ) * ry

wfr2 = Wn1-Wp0
wfl2 = Wp0-Ws1
wfr4 = c1 * wfr2 - c2 * (Wn2-Ws1)
wfl4 = c1 * wfl2 - c2 * (Wn1-Ws2)
fw_r = ( wfr4 * sw1 + wfr2 * (1.0 - sw1) ) * ry
fw_l = ( wfl4 * sw2 + wfl2 * (1.0 - sw2) ) * ry


! 流束の加算　VBCでない面の寄与のみを評価する
cnv_u = cnv_u + 0.5*(Vpn * fu_r * c_n1 + Vps * fu_l * c_s1)
cnv_v = cnv_v + 0.5*(Vpn * fv_r * c_n1 + Vps * fv_l * c_s1)
cnv_w = cnv_w + 0.5*(Vpn * fw_r * c_n1 + Vps * fw_l * c_s1)

! 粘性項の加算
EX = EX + (fu_r * c_n1 - fu_l * c_s1) * ry2
EY = EY + (fv_r * c_n1 - fv_l * c_s1) * ry2
EZ = EZ + (fw_r * c_n1 - fw_l * c_s1) * ry2


! Z方向 ---------------------------------------

! 壁面の場合の参照速度の修正
if ( (b_t2 == 0.0)  ) then
  Ut2 = - v(i  ,j  ,k+1, 1)
  Vt2 = - v(i  ,j  ,k+1, 2)
  Wt2 = - v(i  ,j  ,k+1, 3)
end if

if ( b_t1 == 0.0 ) then
  Ut1 = uq
  Vt1 = vq
  Wt1 = wq
end if

if ( b_b1 == 0.0 ) then
  Ub1 = uq
  Vb1 = vq
  Wb1 = wq
end if

if ( (b_b2 == 0.0)  ) then
  Ub2 = - v(i  ,j  ,k-1, 1)
  Vb2 = - v(i  ,j  ,k-1, 2)
  Wb2 = - v(i  ,j  ,k-1, 3)
end if

! flux
sw1 = b_t1 * c_t2 * ss
sw2 = b_b1 * c_b2 * ss

ufr2 = Ut1-Up0
ufl2 = Up0-Ub1
ufr4 = c1 * ufr2 - c2 * (Ut2-Ub1)
ufl4 = c1 * ufl2 - c2 * (Ut1-Ub2)
fu_r = ( ufr4 * sw1 + ufr2 * (1.0 - sw1) ) * rz
fu_l = ( ufl4 * sw2 + ufl2 * (1.0 - sw2) ) * rz

vfr2 = Vt1-Vp0
vfl2 = Vp0-Vb1
vfr4 = c1 * vfr2 - c2 * (Vt2-Vb1)
vfl4 = c1 * vfl2 - c2 * (Vt1-Vb2)
fv_r = ( vfr4 * sw1 + vfr2 * (1.0 - sw1) ) * rz
fv_l = ( vfl4 * sw2 + vfl2 * (1.0 - sw2) ) * rz

wfr2 = Wt1-Wp0
wfl2 = Wp0-Wb1
wfr4 = c1 * wfr2 - c2 * (Wt2-Wb1)
wfl4 = c1 * wfl2 - c2 * (Wt1-Wb2)
fw_r = ( wfr4 * sw1 + wfr2 * (1.0 - sw1) ) * rz
fw_l = ( wfl4 * sw2 + wfl2 * (1.0 - sw2) ) * rz

! 流束の加算　VBCでない面の寄与のみを評価する  21 flops
cnv_u = cnv_u + 0.5*(Wpt * fu_r * c_t1 + Wpb * fu_l * c_b1)
cnv_v = cnv_v + 0.5*(Wpt * fv_r * c_t1 + Wpb * fv_l * c_b1)
cnv_w = cnv_w + 0.5*(Wpt * fw_r * c_t1 + Wpb * fw_l * c_b1)

! 粘性項の加算  15 flops
EX = EX + (fu_r * c_t1 - fu_l * c_b1) * rz2
EY = EY + (fv_r * c_t1 - fv_l * c_b1) * rz2
EZ = EZ + (fw_r * c_t1 - fw_l * c_b1) * rz2


! 対流項と粘性項の和 > 9 flops
wv(i,j,k,1) = -cnv_u + EX * vcs
wv(i,j,k,2) = -cnv_v + EY * vcs
wv(i,j,k,3) = -cnv_w + EZ * vcs
end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL

return
end subroutine pvec_central


!> ********************************************************************
!! @brief 対流項と粘性項の計算
!! @param [out] wv        疑似ベクトルの空間項 u \frac{\partial u}{\partial x}
!! @param [in]  sz        配列長
!! @param [in]  g         ガイドセル長
!! @param [in]  dh        格子幅
!! @param [in]  c_scheme  対流項スキームのモード（2-Central_2nd, 4-Central_4th）
!! @param [in]  rei       レイノルズ数の逆数
!! @param [in]  v         セルセンター速度ベクトル（n-step）
!! @param [in]  vf        セルフェイス速度ベクトル（n-step）
!! @param [in]  bid       Cut ID
!! @param [in]  bcd       BCindex B
!! @param [in]  vcs_coef  粘性項の係数（粘性項を計算しない場合には0.0）
!! @param [in]  Cs        定数CS
!! @param [in]  imodel    乱流モデル
!! @param [in]  nu        動粘性係数
!! @param [in]  rho       密度
!! @param [out] flop      浮動小数点演算数
!<
subroutine pvec_central_les (wv, sz, g, dh, c_scheme, rei, v, vf, bid, bcd, vcs_coef, Cs, imodel, nu, rho, flop)
implicit none
include 'ffv_f_params.h'
integer                                                   ::  i, j, k, ix, jx, kx, g, c_scheme, bix, bdx
integer, dimension(3)                                     ::  sz
double precision                                          ::  flop
real                                                      ::  b_e1, b_w1, b_n1, b_s1, b_t1, b_b1
real                                                      ::  b_e2, b_w2, b_n2, b_s2, b_t2, b_b2, b_p
real                                                      ::  Up0, Ue1, Ue2, Uw1, Uw2, Us1, Us2, Un1, Un2, Ub1, Ub2, Ut1, Ut2
real                                                      ::  Vp0, Ve1, Ve2, Vw1, Vw2, Vs1, Vs2, Vn1, Vn2, Vb1, Vb2, Vt1, Vt2
real                                                      ::  Wp0, We1, We2, Ww1, Ww2, Ws1, Ws2, Wn1, Wn2, Wb1, Wb2, Wt1, Wt2
real                                                      ::  UPe, UPw, VPn, VPs, WPt, WPb
real                                                      ::  vcs, vcs_coef, ss, sw1, sw2, c1, c2
real                                                      ::  c_e1, c_w1, c_n1, c_s1, c_t1, c_b1
real                                                      ::  c_e2, c_w2, c_n2, c_s2, c_t2, c_b2
real                                                      ::  cnv_u, cnv_v, cnv_w, EX, EY, EZ, rei
real                                                      ::  uq, vq, wq, rx, ry, rz, rx2, ry2 ,rz2, dx, dy, dz
real                                                      ::  ufr2, ufl2, vfr2, vfl2, wfr2, wfl2
real                                                      ::  ufr4, ufl4, vfr4, vfl4, wfr4, wfl4
real                                                      ::  fu_r, fu_l, fv_r, fv_l, fw_r, fw_l
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v, wv, vf
real, dimension(3)                                        ::  dh
integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bid, bcd
integer                                                   ::  imodel
real                                                      ::  Cs, nu, rho, Cw
real                                                      ::  DUDX, DUDY, DUDZ
real                                                      ::  DVDX, DVDY, DVDZ
real                                                      ::  DWDX, DWDY, DWDZ
real                                                      ::  fs, DUDY_w, tauw, utau, yc, yp, up, min_h
real                                                      ::  S11, S12, S13, S21, S22, S23, S31, S32, S33, SSS
real                                                      ::  W11, W12, W13, W21, W22, W23, W31, W32, W33, WWW
real                                                      ::  S_1, S_2, S_3, W_1, W_2, W_3
real                                                      ::  S11d, S12d, S13d, S21d, S22d, S23d, S31d, S32d, S33d
real                                                      ::  Fcs, E_csm, Q_csm, Sijd2, nut
double precision                                          ::  EPS

ix = sz(1)
jx = sz(2)
kx = sz(3)

dx = dh(1)
dy = dh(2)
dz = dh(3)

rx = 1.0/dh(1)
ry = 1.0/dh(2)
rz = 1.0/dh(3)

rx2= rei * rx
ry2= rei * ry
rz2= rei * rz

EPS = 1.0d-10
Cw  = 0.325d0

c1 = 9.0/8.0
c2 = 1.0/24.0

! vcs = 1.0 (Euler Explicit) / 0.5 (CN) / 0.0(No)
vcs = vcs_coef


if ( c_scheme == 2 ) then      !     2nd
ss = 0.0
else if ( c_scheme == 4 ) then !     4th
ss = 1.0
else
write(*,*) 'out of scheme selection'
stop
endif


flop = flop + 46.0d0


!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, rx, ry, rz, rx2, ry2 ,rz2, vcs, ss, c1, c2) &
!$OMP FIRSTPRIVATE(rei, dx, dy, dz) &
!$OMP PRIVATE(cnv_u, cnv_v, cnv_w, bix, bdx, uq, vq, wq, sw1, sw2) &
!$OMP PRIVATE(Up0, Ue1, Ue2, Uw1, Uw2, Us1, Us2, Un1, Un2, Ub1, Ub2, Ut1, Ut2) &
!$OMP PRIVATE(Vp0, Ve1, Ve2, Vw1, Vw2, Vs1, Vs2, Vn1, Vn2, Vb1, Vb2, Vt1, Vt2) &
!$OMP PRIVATE(Wp0, We1, We2, Ww1, Ww2, Ws1, Ws2, Wn1, Wn2, Wb1, Wb2, Wt1, Wt2) &
!$OMP PRIVATE(b_e1, b_w1, b_n1, b_s1, b_t1, b_b1, b_e2, b_w2, b_n2, b_s2, b_t2, b_b2, b_p) &
!$OMP PRIVATE(c_e1, c_w1, c_n1, c_s1, c_t1, c_b1, c_e2, c_w2, c_n2, c_s2, c_t2, c_b2) &
!$OMP PRIVATE(UPe, UPw, VPn, VPs, WPt, WPb) &
!$OMP PRIVATE(ufr2, ufl2, vfr2, vfl2, wfr2, wfl2) &
!$OMP PRIVATE(ufr4, ufl4, vfr4, vfl4, wfr4, wfl4) &
!$OMP PRIVATE(fu_r, fu_l, fv_r, fv_l, fw_r, fw_l) &
!$OMP PRIVATE(EX, EY, EZ) &
!$OMP FIRSTPRIVATE(imodel, EPS, Cs, nu, rho, Cw) &
!$OMP PRIVATE(DUDX, DUDY, DUDZ, DVDX, DVDY, DVDZ, DWDX, DWDY, DWDZ) &
!$OMP PRIVATE(S11, S12, S13, S21, S22, S23, S31, S32, S33, SSS) &
!$OMP PRIVATE(W11, W12, W13, W21, W22, W23, W31, W32, W33, WWW) &
!$OMP PRIVATE(S_1, S_2, S_3, W_1, W_2, W_3) &
!$OMP PRIVATE(S11d, S12d, S13d, S21d, S22d, S23d, S31d, S32d, S33d) &
!$OMP PRIVATE(Fcs, E_csm, Q_csm, Sijd2, nut) &
!$OMP PRIVATE(fs, DUDY_w, tauw, utau, yc, yp, up, min_h)

!$OMP DO SCHEDULE(static) COLLAPSE(2)

do k=1,kx
do j=1,jx
do i=1,ix
cnv_u = 0.0
cnv_v = 0.0
cnv_w = 0.0
EX = 0.0
EY = 0.0
EZ = 0.0

! 変数のロード
! 各軸方向5点の変数ロード
Ub2 = v(i  ,j  ,k-2, 1)
Ub1 = v(i  ,j  ,k-1, 1)
Us2 = v(i  ,j-2,k  , 1)
Us1 = v(i  ,j-1,k  , 1)
Uw2 = v(i-2,j  ,k  , 1)
Uw1 = v(i-1,j  ,k  , 1)
Up0 = v(i  ,j  ,k  , 1)
Ue1 = v(i+1,j  ,k  , 1)
Ue2 = v(i+2,j  ,k  , 1)
Un1 = v(i  ,j+1,k  , 1)
Un2 = v(i  ,j+2,k  , 1)
Ut1 = v(i  ,j  ,k+1, 1)
Ut2 = v(i  ,j  ,k+2, 1)

Vb2 = v(i  ,j  ,k-2, 2)
Vb1 = v(i  ,j  ,k-1, 2)
Vs2 = v(i  ,j-2,k  , 2)
Vs1 = v(i  ,j-1,k  , 2)
Vw2 = v(i-2,j  ,k  , 2)
Vw1 = v(i-1,j  ,k  , 2)
Vp0 = v(i  ,j  ,k  , 2)
Ve1 = v(i+1,j  ,k  , 2)
Ve2 = v(i+2,j  ,k  , 2)
Vn1 = v(i  ,j+1,k  , 2)
Vn2 = v(i  ,j+2,k  , 2)
Vt1 = v(i  ,j  ,k+1, 2)
Vt2 = v(i  ,j  ,k+2, 2)

Wb2 = v(i  ,j  ,k-2, 3)
Wb1 = v(i  ,j  ,k-1, 3)
Ws2 = v(i  ,j-2,k  , 3)
Ws1 = v(i  ,j-1,k  , 3)
Ww2 = v(i-2,j  ,k  , 3)
Ww1 = v(i-1,j  ,k  , 3)
Wp0 = v(i  ,j  ,k  , 3)
We1 = v(i+1,j  ,k  , 3)
We2 = v(i+2,j  ,k  , 3)
Wn1 = v(i  ,j+1,k  , 3)
Wn2 = v(i  ,j+2,k  , 3)
Wt1 = v(i  ,j  ,k+1, 3)
Wt2 = v(i  ,j  ,k+2, 3)

bix = bid(i,j,k)
bdx = bcd(i,j,k)

! (i,j,k)からみたセル状態 (0-solid / 1-fluid)
b_p = real(ibits(bdx, State, 1))

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
c_w1 = real( ibits(bdx, bc_d_W, 1) )
c_e1 = real( ibits(bdx, bc_d_E, 1) )
c_s1 = real( ibits(bdx, bc_d_S, 1) )
c_n1 = real( ibits(bdx, bc_d_N, 1) )
c_b1 = real( ibits(bdx, bc_d_B, 1) )
c_t1 = real( ibits(bdx, bc_d_T, 1) )


! ステンシルの参照先がvspec, outflowである場合のスキームの破綻を回避，１次精度におとすSW
c_w2 = real( ibits(bcd(i-1, j  , k  ), bc_d_W, 1) )
c_e2 = real( ibits(bcd(i+1, j  , k  ), bc_d_E, 1) )
c_s2 = real( ibits(bcd(i  , j-1, k  ), bc_d_S, 1) )
c_n2 = real( ibits(bcd(i  , j+1, k  ), bc_d_N, 1) )
c_b2 = real( ibits(bcd(i  , j  , k-1), bc_d_B, 1) )
c_t2 = real( ibits(bcd(i  , j  , k+1), bc_d_T, 1) )


! 界面速度（スタガード位置） 24 flops
UPe = vf(i  , j  , k  ,1)*b_e1
UPw = vf(i-1, j  , k  ,1)*b_w1
VPn = vf(i  , j  , k  ,2)*b_n1
VPs = vf(i  , j-1, k  ,2)*b_s1
WPt = vf(i  , j  , k  ,3)*b_t1
WPb = vf(i  , j  , k-1,3)*b_b1

! セルセンターからの壁面修正速度 3 flops
uq = - Up0
vq = - Vp0
wq = - Wp0

! X方向 ---------------------------------------

! 速度指定の場合に参照先として，固体内にテンポラリに与えた値を使う
if ( (b_e2 == 0.0)  ) then
Ue2 = - v(i+1,j  ,k  , 1)
Ve2 = - v(i+1,j  ,k  , 2)
We2 = - v(i+1,j  ,k  , 3)
endif

if ( b_e1 == 0.0 ) then
Ue1 = uq
Ve1 = vq
We1 = wq
endif

if ( b_w1 == 0.0 ) then
Uw1 = uq
Vw1 = vq
Ww1 = wq
end if

if ( (b_w2 == 0.0)  ) then
Uw2 = - v(i-1,j  ,k  , 1)
Vw2 = - v(i-1,j  ,k  , 2)
Ww2 = - v(i-1,j  ,k  , 3)
end if

! flux u \frac{\partial u}{\partial x}
! 壁面がある場合　b_e1, b_w1 (0-wall face / 1-fluid)
! vspec/outflowを参照する場合(c_e2, c_w2)は二次精度へおとす
! 二次精度の時には、ssで強制
sw1 = b_e1 * c_e2 * ss ! 4 flops
sw2 = b_w1 * c_w2 * ss

ufr2 = Ue1-Up0 ! 20 flops
ufl2 = Up0-Uw1
ufr4 = c1 * ufr2 - c2 * (Ue2-Uw1)
ufl4 = c1 * ufl2 - c2 * (Ue1-Uw2)
fu_r = ( ufr4 * sw1 + ufr2 * (1.0 - sw1) ) * rx
fu_l = ( ufl4 * sw2 + ufl2 * (1.0 - sw2) ) * rx
DUDX = 0.5d0 * ( fu_r * c_e1 + fu_l * c_w1 )


vfr2 = Ve1-Vp0
vfl2 = Vp0-Vw1
vfr4 = c1 * vfr2 - c2 * (Ve2-Vw1)
vfl4 = c1 * vfl2 - c2 * (Ve1-Vw2)
fv_r = ( vfr4 * sw1 + vfr2 * (1.0 - sw1) ) * rx
fv_l = ( vfl4 * sw2 + vfl2 * (1.0 - sw2) ) * rx
DVDX = 0.5d0 * ( fv_r * c_e1 + fv_l * c_w1 )


wfr2 = We1-Wp0
wfl2 = Wp0-Ww1
wfr4 = c1 * wfr2 - c2 * (We2-Ww1)
wfl4 = c1 * wfl2 - c2 * (We1-Ww2)
fw_r = ( wfr4 * sw1 + wfr2 * (1.0 - sw1) ) * rx
fw_l = ( wfl4 * sw2 + wfl2 * (1.0 - sw2) ) * rx
DWDX = 0.5d0 * ( fw_r * c_e1 + fw_l * c_w1 )


! 流束の加算　平均勾配　VBCでない面の寄与のみを評価する 21 flops
cnv_u = cnv_u + 0.5*(Upe * fu_r * c_e1 + Upw * fu_l * c_w1)
cnv_v = cnv_v + 0.5*(Upe * fv_r * c_e1 + Upw * fv_l * c_w1)
cnv_w = cnv_w + 0.5*(Upe * fw_r * c_e1 + Upw * fw_l * c_w1)

! 粘性項の加算 15 flops
EX = EX + (fu_r * c_e1 - fu_l * c_w1) * rx2
EY = EY + (fv_r * c_e1 - fv_l * c_w1) * rx2
EZ = EZ + (fw_r * c_e1 - fw_l * c_w1) * rx2


! Y方向 ---------------------------------------

if ( (b_n2 == 0.0)  ) then
Un2 = - v(i  ,j+1,k  , 1)
Vn2 = - v(i  ,j+1,k  , 2)
Wn2 = - v(i  ,j+1,k  , 3)
endif

if ( b_n1 == 0.0 ) then
Un1 = uq
Vn1 = vq
Wn1 = wq
endif

if ( b_s1 == 0.0 ) then
Us1 = uq
Vs1 = vq
Ws1 = wq
endif

if ( (b_s2 == 0.0)  ) then
Us2 = - v(i  ,j-1,k  , 1)
Vs2 = - v(i  ,j-1,k  , 2)
Ws2 = - v(i  ,j-1,k  , 3)
endif

! flux
sw1 = b_n1 * c_n2 * ss
sw2 = b_s1 * c_s2 * ss

ufr2 = Un1-Up0
ufl2 = Up0-Us1
ufr4 = c1 * ufr2 - c2 * (Un2-Us1)
ufl4 = c1 * ufl2 - c2 * (Un1-Us2)
fu_r = ( ufr4 * sw1 + ufr2 * (1.0 - sw1) ) * ry
fu_l = ( ufl4 * sw2 + ufl2 * (1.0 - sw2) ) * ry
DUDY = 0.5d0 * ( fu_r * c_n1 + fu_l * c_s1 )


vfr2 = Vn1-Vp0
vfl2 = Vp0-Vs1
vfr4 = c1 * vfr2 - c2 * (Vn2-Vs1)
vfl4 = c1 * vfl2 - c2 * (Vn1-Vs2)
fv_r = ( vfr4 * sw1 + vfr2 * (1.0 - sw1) ) * ry
fv_l = ( vfl4 * sw2 + vfl2 * (1.0 - sw2) ) * ry
DVDY = 0.5d0 * ( fv_r * c_n1 + fv_l * c_s1 )


wfr2 = Wn1-Wp0
wfl2 = Wp0-Ws1
wfr4 = c1 * wfr2 - c2 * (Wn2-Ws1)
wfl4 = c1 * wfl2 - c2 * (Wn1-Ws2)
fw_r = ( wfr4 * sw1 + wfr2 * (1.0 - sw1) ) * ry
fw_l = ( wfl4 * sw2 + wfl2 * (1.0 - sw2) ) * ry
DWDY = 0.5d0 * ( fw_r * c_n1 + fw_l * c_s1 )


! 流束の加算　VBCでない面の寄与のみを評価する
cnv_u = cnv_u + 0.5*(Vpn * fu_r * c_n1 + Vps * fu_l * c_s1)
cnv_v = cnv_v + 0.5*(Vpn * fv_r * c_n1 + Vps * fv_l * c_s1)
cnv_w = cnv_w + 0.5*(Vpn * fw_r * c_n1 + Vps * fw_l * c_s1)

! 粘性項の加算
EX = EX + (fu_r * c_n1 - fu_l * c_s1) * ry2
EY = EY + (fv_r * c_n1 - fv_l * c_s1) * ry2
EZ = EZ + (fw_r * c_n1 - fw_l * c_s1) * ry2


! Z方向 ---------------------------------------

! 壁面の場合の参照速度の修正
if ( (b_t2 == 0.0)  ) then
Ut2 = - v(i  ,j  ,k+1, 1)
Vt2 = - v(i  ,j  ,k+1, 2)
Wt2 = - v(i  ,j  ,k+1, 3)
end if

if ( b_t1 == 0.0 ) then
Ut1 = uq
Vt1 = vq
Wt1 = wq
end if

if ( b_b1 == 0.0 ) then
Ub1 = uq
Vb1 = vq
Wb1 = wq
end if

if ( (b_b2 == 0.0)  ) then
Ub2 = - v(i  ,j  ,k-1, 1)
Vb2 = - v(i  ,j  ,k-1, 2)
Wb2 = - v(i  ,j  ,k-1, 3)
end if

! flux
sw1 = b_t1 * c_t2 * ss
sw2 = b_b1 * c_b2 * ss

ufr2 = Ut1-Up0
ufl2 = Up0-Ub1
ufr4 = c1 * ufr2 - c2 * (Ut2-Ub1)
ufl4 = c1 * ufl2 - c2 * (Ut1-Ub2)
fu_r = ( ufr4 * sw1 + ufr2 * (1.0 - sw1) ) * rz
fu_l = ( ufl4 * sw2 + ufl2 * (1.0 - sw2) ) * rz
DUDZ = 0.5d0 * ( fu_r * c_t1 + fu_l * c_b1 )


vfr2 = Vt1-Vp0
vfl2 = Vp0-Vb1
vfr4 = c1 * vfr2 - c2 * (Vt2-Vb1)
vfl4 = c1 * vfl2 - c2 * (Vt1-Vb2)
fv_r = ( vfr4 * sw1 + vfr2 * (1.0 - sw1) ) * rz
fv_l = ( vfl4 * sw2 + vfl2 * (1.0 - sw2) ) * rz
DVDZ = 0.5d0 * ( fv_r * c_t1 + fv_l * c_b1 )


wfr2 = Wt1-Wp0
wfl2 = Wp0-Wb1
wfr4 = c1 * wfr2 - c2 * (Wt2-Wb1)
wfl4 = c1 * wfl2 - c2 * (Wt1-Wb2)
fw_r = ( wfr4 * sw1 + wfr2 * (1.0 - sw1) ) * rz
fw_l = ( wfl4 * sw2 + wfl2 * (1.0 - sw2) ) * rz
DWDZ = 0.5d0 * ( fw_r * c_t1 + fw_l * c_b1 )


! 流束の加算　VBCでない面の寄与のみを評価する
cnv_u = cnv_u + 0.5*(Wpt * fu_r * c_t1 + Wpb * fu_l * c_b1)
cnv_v = cnv_v + 0.5*(Wpt * fv_r * c_t1 + Wpb * fv_l * c_b1)
cnv_w = cnv_w + 0.5*(Wpt * fw_r * c_t1 + Wpb * fw_l * c_b1)

! 粘性項の加算
EX = EX + (fu_r * c_t1 - fu_l * c_b1) * rz2
EY = EY + (fv_r * c_t1 - fv_l * c_b1) * rz2
EZ = EZ + (fw_r * c_t1 - fw_l * c_b1) * rz2



!----------- Sij (strain rate tensor)  10 + 18 = 28
S11 = DUDX
S12 = (DVDX + DUDY) * 0.5d0
S13 = (DWDX + DUDZ) * 0.5d0

S21 = S12
S22 = DVDY
S23 = (DWDY + DVDZ) * 0.5d0

S31 = S13
S32 = S23
S33 = DWDZ

SSS = DSQRT( 2.0d0 * (S11*S11 + S22*S22 + S33*S33) &
           + 4.0d0 * (S12*S12 + S23*S23 + S31*S31) )


!----------- Wij (vorticity tensor)  12 + 22 = 34
W11 = 0.0d0
W12 = (DVDX - DUDY) * 0.5d0
W13 = (DWDX - DUDZ) * 0.5d0

W21 = (DUDY - DVDX) * 0.5d0
W22 = 0.0d0
W23 = (DWDY - DVDZ) * 0.5d0

W31 = (DUDZ - DWDX) * 0.5d0
W32 = (DVDZ - DWDY) * 0.5d0
W33 = 0.0d0

WWW = DSQRT( 2.0d0 * (W12*W12 + W13*W13 + W21*W21 &
                    + W23*W23 + W31*W31 + W32*W32 ) )


!----------- Smagorinsky model   83 flops
if (imodel == 1) then
yc     = (j - 0.5d0)*dy
min_h  = min(yc, dy*jx - yc)
DUDY_w = v(i, jx, k, 1) * (ry*2.0d0)
tauw   = (nu * rho) * abs(DUDY_w)
utau   = sqrt(tauw/rho)
yp     = min_h * utau / nu
up     = Up0 / utau
fs     = 1.0d0 - exp(-yp/26.0d0)
nut = (Cs * fs * dy) * (Cs * fs * dy) * SSS
flop = flop + 83.0d0
end if


!----------- CSM model 17 + 28 = 45
if (imodel == 2) then
Q_csm = (WWW * WWW - SSS * SSS) * 0.25d0
E_csm = (WWW * WWW + SSS * SSS) * 0.25d0
Fcs   = Q_csm / (E_csm + EPS)
nut = (1.0d0/22.0d0) * abs(Fcs) * sqrt(abs(Fcs)) * (1.0d0 - Fcs) * (dx*dy*dz)*(2.0d0/3.0d0) * SSS
flop = flop + 45.0d0
end if


!----------- WALE model
if (imodel == 3) then

S11d = 0.5d0*(DUDX*DUDX + DUDX*DUDX) - (1.0d0/3.0d0)*(DUDX*DUDX + DVDY*DVDY + DWDZ*DWDZ)
S12d = 0.5d0*(DUDY*DUDY + DVDX*DVDX)
S13d = 0.5d0*(DUDZ*DUDZ + DWDX*DWDX)
S21d = S12d
S22d = 0.5d0*(DVDY*DVDY + DVDY*DVDY) - (1.0d0/3.0d0)*(DUDX*DUDX + DVDY*DVDY + DWDZ*DWDZ)
S23d = 0.5d0*(DVDZ*DVDZ + DWDY*DWDY)
S31d = S13d
S32d = S23d
S33d = 0.5d0*(DWDZ*DWDZ + DWDZ*DWDZ) - (1.0d0/3.0d0)*(DUDX*DUDX + DVDY*DVDY + DWDZ*DWDZ)

Sijd2 =  S11d*S11d + S21d*S21d + S31d*S31d &
       + S12d*S12d + S22d*S22d + S32d*S32d &
       + S13d*S13d + S23d*S23d + S33d*S33d

nut = (Cw * (dx*dy*dz)**(1.0d0/3.0d0))* (Cw * (dx*dy*dz)**(1.0d0/3.0d0)) * (Sijd2)**(3.0d0/2.0d0)  &
    / ( (SSS * SSS)**(5.0d0/2.0d0) + (Sijd2)**(5.0d0/4.0d0) )
flop = flop + 213.0d0
end if

if (imodel == 0) then
nut = 0.0
end if

! 対流項と粘性項の和 > 36 flops
wv(i, j, k, 1) = -cnv_u + EX * ( 1.0d0 + nut/nu ) * vcs
wv(i, j, k, 2) = -cnv_v + EY * ( 1.0d0 + nut/nu ) * vcs
wv(i, j, k, 3) = -cnv_w + EZ * ( 1.0d0 + nut/nu ) * vcs

! 24 + 3 + 3*106 + 28 + 34 + 36 = 443
flop = flop + 443.0d0

end do
end do
end do

!$OMP END DO
!$OMP END PARALLEL


return
end subroutine pvec_central_les

!> ********************************************************************
!! @brief 4次精度の打ち切り誤差を持つポアソン反復の1ステップ目のソース項
!! @param [out]    rhs   ソース項
!! @param [in]     sz    配列長
!! @param [in]     g     ガイドセル長
!! @param [in]     dt    時間積分幅
!! @param [in]     dh    格子幅
!! @param [in]     div   divergence
!! @param [in,out] flop  浮動小数点演算数
!<
subroutine src_trnc (rhs, sz, g, dt, dh, div, flop)
implicit none
include 'ffv_f_params.h'
integer                                                   ::  i, j, k, ix, jx, kx, g
integer, dimension(3)                                     ::  sz
double precision                                          ::  flop
real                                                      ::  dt, rx, ry, rz, drx, dry, drz
real                                                      ::  bw2, bw1, be2, be1, bs2, bs1, bn2, bn1, bb2, bb1, bt2, bt1, b0
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  rhs, div
real, dimension(3)                                        ::  dh

ix = sz(1)
jx = sz(2)
kx = sz(3)

rx = 1.0 / dh(1)
ry = 1.0 / dh(2)
rz = 1.0 / dh(3)

drx = rx / dt
dry = ry / dt
drz = rz / dt

flop  = flop + dble(ix)*dble(jx)*dble(kx)*24.0d0 + 48.0d0

!$OMP PARALLEL &
!$OMP PRIVATE(bw2, bw1, be2, be1, bs2, bs1, bn2, bn1, bb2, bb1, bt2, bt1, b0) &
!$OMP FIRSTPRIVATE(ix, jx, kx, drx, dry, drz)

!$OMP DO SCHEDULE(static) COLLAPSE(2)

do k=1,kx
do j=1,jx
do i=1,ix

b0  = div(i  , j  , k  )

bw2 = div(i-2, j  , k  )
bw1 = div(i-1, j  , k  )
be1 = div(i+1, j  , k  )
be2 = div(i+2, j  , k  )

bs2 = div(i  , j-2, k  )
bs1 = div(i  , j-1, k  )
bn1 = div(i  , j+1, k  )
bn2 = div(i  , j+2, k  )

bb2 = div(i  , j  , k-2)
bb1 = div(i  , j  , k-1)
bt1 = div(i  , j  , k+1)
bt2 = div(i  , j  , k+2)

rhs(i,j,k) = ( bw2 - 4.0*bw1 + 6.0*b0 - 4.0*be1 + be2 ) * drx &
           + ( bs2 - 4.0*bs1 + 6.0*b0 - 4.0*bn1 + bn2 ) * dry &
           + ( bb2 - 4.0*bb1 + 6.0*b0 - 4.0*bt1 + bt2 ) * drz

end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL

return
end subroutine src_trnc


!> ********************************************************************
!! @brief 4次精度の打ち切り誤差を持つポアソン反復の1ステップ目のソース項
!! @param [out]    rhs  ソース項
!! @param [in]     sz   配列長
!! @param [in]     g    ガイドセル長
!! @param [in]     dt   時間積分幅
!! @param [in]     dh   格子幅
!! @param [in]     div  divergence
!! @param [in]     b    \nabla{u}
!! @param [in,out] flop 浮動小数点演算数
!<
subroutine src_1st (rhs, sz, g, dt, dh, div, b, flop)
implicit none
include 'ffv_f_params.h'
integer                                                   ::  i, j, k, ix, jx, kx, g
integer, dimension(3)                                     ::  sz
double precision                                          ::  flop
real                                                      ::  dt, rx, ry, rz
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  rhs, div, b
real, dimension(3)                                        ::  dh

ix = sz(1)
jx = sz(2)
kx = sz(3)

rx = 1.0 / dh(1)
ry = 1.0 / dh(2)
rz = 1.0 / dh(3)

flop  = flop + dble(ix)*dble(jx)*dble(kx)*9.0d0 + 9.0d0

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx)

!$OMP DO SCHEDULE(static) COLLAPSE(2)

do k=1,kx
do j=1,jx
do i=1,ix

rhs(i,j,k) = b(i,j,k) - &
( div(i+1,j  ,k  ) + div(i-1,j  ,k  ) &
+ div(i  ,j+1,k  ) + div(i  ,j-1,k  ) &
+ div(i  ,j  ,k+1) + div(i  ,j  ,k-1) &
- 6.0 * div(i,j,k) )

end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL

return
end subroutine src_1st


!> ********************************************************************
!! @brief 4次精度の打ち切り誤差を持つポアソン反復の2ステップ目のソース項
!! @param [out]    rhs  ソース項
!! @param [in]     sz   配列長
!! @param [in]     g    ガイドセル長
!! @param [in]     dh   格子幅
!! @param [in]     b    h^2/dt \nabla{u}
!! @param [in]     psi  1ステップ目の反復結果
!! @param [in,out] flop 浮動小数点演算数
!<
subroutine src_2nd (rhs, sz, g, dh, b, psi, flop)
implicit none
include 'ffv_f_params.h'
integer                                                   ::  i, j, k, ix, jx, kx, g
integer, dimension(3)                                     ::  sz
double precision                                          ::  flop
real                                                      ::  h2
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  rhs, b, psi
real, dimension(3)                                        ::  dh

ix = sz(1)
jx = sz(2)
kx = sz(3)

h2 = 0.25*dh(1)*dh(1)

flop  = flop + dble(ix)*dble(jx)*dble(kx)*2.0d0

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, h2)

!$OMP DO SCHEDULE(static) COLLAPSE(2)

do k=1,kx
do j=1,jx
do i=1,ix

rhs(i,j,k) = b(i,j,k) - h2 * psi(i,j,k)

end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL

return
end subroutine src_2nd


!> ********************************************************************
!! @brief 安定化項
!! @param [in,out] vc      疑似ベクトル
!! @param [in]     sz      配列長
!! @param [in]     g       ガイドセル長
!! @param [in]     dt      時間積分幅
!! @param [in]     v       速度ベクトル（n-step, collocated）
!! @param [in]     bcd     BCindex B
!! @param [in]     vref    参照速度の大きさ v00[0]
!! @param [in]     st      区間開始の無次元速度
!! @param [in]     ed      区間終了の無次元速度
!! @param [in]     penalty ペナルティ数
!! @param [in]     count   修正数
!! @param [in,out] flop    浮動小数点演算数
!<
subroutine stabilize (vc, sz, g, dt, v, bcd, vref, st, ed, penalty, count, flop)
implicit none
include 'ffv_f_params.h'
integer                                                   ::  i, j, k, ix, jx, kx, g, count
integer, dimension(3)                                     ::  sz
double precision                                          ::  flop
real                                                      ::  actv, dt, st, ed, penalty, pai
real                                                      ::  u1, v1, w1, uu, gma, df, vref, s1, s2
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  vc, v
integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bcd

ix = sz(1)
jx = sz(2)
kx = sz(3)

flop = flop + dble(ix)*dble(jx)*dble(kx)*42.0d0

count = 0

if ( vref == 0.0 ) vref=1.0
pai = 2.0 * asin(1.0)
df = ed - st
s1 = st*vref
s2 = ed*vref

!$OMP PARALLEL &
!$OMP REDUCTION(+:count) &
!$OMP PRIVATE(actv, uu, u1, v1, w1, gma) &
!$OMP FIRSTPRIVATE(ix, jx, kx, dt, st, ed, penalty, df, pai, vref)

!$OMP DO SCHEDULE(static) COLLAPSE(2)
do k=1,kx
do j=1,jx
do i=1,ix

u1 = v(i,j,k,1)
v1 = v(i,j,k,2)
w1 = v(i,j,k,3)
uu = sqrt(u1*u1 + v1*v1 + w1*w1) / vref

if ( uu < st ) then
  gma = 0.0
else if ( uu > ed ) then
  gma = penalty
else
  gma = 0.5 * (1.0 - cos(pai*(uu-st)/df) ) * penalty
  count = count + 1
end if

actv = dt * real(ibits(bcd(i,j,k), State, 1)) * gma

vc(i,j,k,1) = vc(i,j,k,1) - v(i,j,k,1) * actv
vc(i,j,k,2) = vc(i,j,k,2) - v(i,j,k,2) * actv
vc(i,j,k,3) = vc(i,j,k,3) - v(i,j,k,3) * actv

end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL

return
end subroutine stabilize
