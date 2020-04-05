!##################################################################################
!
! RIAM-COMPACT HPC version : RAinWATER
!
! Copyright (C) 2015-2018 Research Institute for Applied Mechanics(RIAM)
!                       / Research Institute for Information Technology(RIIT), Kyushu University.
! All rights reserved.
!
! Copyright (C) 2015-2018 Advanced Institute for Computational Science, RIKEN.
! All rights reserved.
!
! Copyright (c) 2016-2020 Research Institute for Information Research, Kyushu University.
! All rights reserved.
!
!##################################################################################


!> ********************************************************************
!! @brief sphフォーマットのスカラーデータ書き出し
!! @param [in]     sz          配列長
!! @param [in]     org         起点座標
!! @param [in]     pch         格子幅（ダミー）
!! @param [in]     stp         ステップ数
!! @param [in]     tm          時刻
!! @param [in]     R           tensor配列（ローカル）
!! @param [in]     c           番号
!! @param [in]     fname       ファイル名
!! @param [in]     dtype       1-float / 2-double
!! @note 倍精度の場合、整数は4バイトのままかもしれない
!<
subroutine sph_write_tensor (sz, g, org, pch, stp, tm, R, c, fname, dtype)
implicit none
integer                                                  ::  i, j, k, ix, jx, kx, g, c
integer                                                  ::  sv, dtype, stp
integer, dimension(3)                                    ::  sz
real                                                     ::  tm
real, dimension(3)                                       ::  org, pch
real, dimension(6, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)::  R
character*100                                            ::  fname

ix = sz(1)
jx = sz(2)
kx = sz(3)
sv = 1

open (unit=12, file=fname, form='unformatted', status='replace')
write(12) sv, dtype
write(12) ix, jx, kx
write(12) org(1), org(2), org(3)
write(12) pch(1), pch(2), pch(3)
write(12) stp, tm
write(12) (((R(c,i,j,k), i=1,ix), j=1,jx), k=1,kx)
close(unit=12)

return
end subroutine sph_write_tensor
