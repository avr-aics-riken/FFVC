!###################################################################################
!
! Flow Base class
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
!
!> @file ffv_f_params.h
!! @brief BC Indexのマスク操作のビット定義, FB_Define.hと整合
!! @author aics
!<

integer     ::  Active, State, id_fluid, id_solid
integer     ::  bc_diag, bc_mask30, file_path_length
integer     ::  bc_d_T, bc_d_B, bc_d_N, bc_d_S, bc_d_E, bc_d_W
integer     ::  bc_n_T, bc_n_B, bc_n_N, bc_n_S, bc_n_E, bc_n_W
integer     ::  bc_ndag_T, bc_ndag_B, bc_ndag_N, bc_ndag_S, bc_ndag_E, bc_ndag_W
integer     ::  bc_face_T, bc_face_B, bc_face_N, bc_face_S, bc_face_E, bc_face_W
integer     ::  adbtc_T, adbtc_B, adbtc_N, adbtc_S, adbtc_E, adbtc_W
integer     ::  X_plus, X_minus, Y_plus, Y_minus, Z_plus, Z_minus
integer     ::  id_specvel, id_wall, id_symmetric, id_periodic
integer     ::  obc_periodic, obc_outflow, obc_mask
integer     ::  cmp_typ_outflow, cmp_typ_hex, cmp_typ_vspec, cmp_typ_solidrev
integer     ::  facing_t, facing_b, facing_n, facing_s, facing_e, facing_w, vld_cnvg
integer     ::  gma_T, gma_B, gma_N, gma_S, gma_E, gma_W, h_diag
integer     ::  top_vf, bitw_8, forcing_bit, bitw_6, bitw_5, vbc_uwd, bitw_9, bitw_10

parameter ( bc_mask30 = Z'3fffffff') ! 16進表記，VBCの6面(30bit)をまとめたマスク
parameter ( bitw_10= 10) ! FB_Define.h MASK_10 10bit幅
parameter ( bitw_9 = 9)  ! FB_Define.h MASK_9 9bit幅
parameter ( bitw_8 = 8)  ! FB_Define.h MASK_8 8bit幅
parameter ( bitw_6 = 6)  ! FB_Define.h MASK_6 6bit幅
parameter ( bitw_5 = 5)  ! FB_Define.h MASK_5 5bit幅

parameter ( X_minus = 0 ) ! FB_Define.h X_MINUS or enum
parameter ( X_plus  = 1 ) ! FB_Define.h X_PLUS
parameter ( Y_minus = 2 ) ! FB_Define.h Y_MINUS
parameter ( Y_plus  = 3 ) ! FB_Define.h Y_PLUS
parameter ( Z_minus = 4 ) ! FB_Define.h Z_MINUS
parameter ( Z_plus  = 5 ) ! FB_Define.h Z_PLUS

! 状態
parameter ( id_solid  = 0 ) ! FB_Define.h SOLID
parameter ( id_fluid  = 1 ) ! FB_Define.h FLUID

! ビットフラグ共通
parameter ( Active  = 31 ) ! FB_Define.h ACTIVE_BIT
parameter ( State   = 30 ) ! FB_Define.h STATE_BIT

! BCindex B
parameter ( forcing_bit = 28) ! FB_Define.h FORCING_BIT
parameter ( top_vf      = 20) ! FB_Define.h TOP_VF
parameter ( adbtc_T = 19) ! FB_Define.h ADIABATIC_T
parameter ( adbtc_B = 18) ! FB_Define.h ADIABATIC_B
parameter ( adbtc_N = 17) ! FB_Define.h ADIABATIC_N
parameter ( adbtc_S = 16) ! FB_Define.h ADIABATIC_S
parameter ( adbtc_E = 15) ! FB_Define.h ADIABATIC_E
parameter ( adbtc_W = 14) ! FB_Define.h ADIABATIC_W
parameter ( gma_T   = 13) ! FB_Define.h GMA_T
parameter ( gma_B   = 12) ! FB_Define.h GMA_B
parameter ( gma_N   = 11) ! FB_Define.h GMA_N
parameter ( gma_S   = 10) ! FB_Define.h GMA_S
parameter ( gma_E   =  9) ! FB_Define.h GMA_E
parameter ( gma_W   =  8) ! FB_Define.h GMA_W
parameter ( h_diag  =  5) ! FB_Define.h H_DIAG


! BCindex P
parameter ( bc_d_T    = 29 ) ! FB_Define.h BC_D_T
parameter ( bc_d_B    = 28 ) ! FB_Define.h BC_D_B
parameter ( bc_d_N    = 27 ) ! FB_Define.h BC_D_N
parameter ( bc_d_S    = 26 ) ! FB_Define.h BC_D_S
parameter ( bc_d_E    = 25 ) ! FB_Define.h BC_D_E
parameter ( bc_d_W    = 24 ) ! FB_Define.h BC_D_W
parameter ( bc_n_T    = 23 ) ! FB_Define.h BC_N_T
parameter ( bc_n_B    = 22 ) ! FB_Define.h BC_N_B
parameter ( bc_n_N    = 21 ) ! FB_Define.h BC_N_N
parameter ( bc_n_S    = 20 ) ! FB_Define.h BC_N_S
parameter ( bc_n_E    = 19 ) ! FB_Define.h BC_N_E
parameter ( bc_n_W    = 18 ) ! FB_Define.h BC_N_W
parameter ( bc_ndag_T = 17 ) ! FB_Define.h BC_NDAG_T
parameter ( bc_ndag_B = 16 ) ! FB_Define.h BC_NDAG_B
parameter ( bc_ndag_N = 15 ) ! FB_Define.h BC_NDAG_N
parameter ( bc_ndag_S = 14 ) ! FB_Define.h BC_NDAG_S
parameter ( bc_ndag_E = 13 ) ! FB_Define.h BC_NDAG_E
parameter ( bc_ndag_W = 12 ) ! FB_Define.h BC_NDAG_W
parameter ( bc_diag   =  9 ) ! FB_Define.h BC_DIAG
parameter ( facing_t  =  8 ) ! FB_Define.h FACING_T
parameter ( facing_b  =  7 ) ! FB_Define.h FACING_B
parameter ( facing_n  =  6 ) ! FB_Define.h FACING_N
parameter ( facing_s  =  5 ) ! FB_Define.h FACING_S
parameter ( facing_e  =  4 ) ! FB_Define.h FACING_E
parameter ( facing_w  =  3 ) ! FB_Define.h FACING_W
parameter ( vld_cnvg  =  2 ) ! FB_Define.h VLD_CNVG
parameter ( vbc_uwd   =  1 ) ! FB_Define.h VBC_UWD

! BCindex C
parameter ( bc_face_T = 25 ) ! FB_Define.h BC_FACE_T
parameter ( bc_face_B = 20 ) ! FB_Define.h BC_FACE_B
parameter ( bc_face_N = 15 ) ! FB_Define.h BC_FACE_N
parameter ( bc_face_S = 10 ) ! FB_Define.h BC_FACE_S
parameter ( bc_face_E = 5  ) ! FB_Define.h BC_FACE_E
parameter ( bc_face_W = 0  ) ! FB_Define.h BC_FACE_W


! 速度境界条件の種類
parameter ( id_specvel   = 0 ) ! SetBC.h outer boundary condition enum
parameter ( id_wall      = 2 ) ! SetBC.h 
parameter ( id_symmetric = 3 ) ! SetBC.h 
parameter ( id_periodic  = 4 ) ! SetBC.h 

! 内部コンポーネント番号
parameter ( cmp_typ_vspec   = 7  ) ! FB_Define.h SPEC_VEL
parameter ( cmp_typ_outflow = 8  ) ! FB_Define.h OUTFLOW
parameter ( cmp_typ_hex     = 12 ) ! FB_Define.h HEX
parameter ( cmp_typ_solidrev= 17 ) ! FB_Define.h SOLIDREV

! 外部境界条件番号
parameter ( obc_mask     = 31) ! FB_Define.h OBC_MASK 内外部境界条件の識別子
parameter ( obc_outflow  = 3 ) ! FB_Define.h OBC_OUTFLOW 
parameter ( obc_periodic = 5 ) ! FB_Define.h OBC_PERIODIC

! ファイルパスの長さ
parameter ( file_path_length = 64 ) ! #define FB_FILE_PATH_LENGTH 64


