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

!> @file   FB_util.f90
!! @brief  FlowBase utilities
!! @author aics
!<

!> ********************************************************************
!! @brief 粗い格子から密な格子への補間（ゼロ次）
!! @param [out] dst 密な格子系
!! @param [in]  sz 配列長
!! @param [in]  g ガイドセル長
!! @param [in]  src 粗い格子系
!! @param [in]  st 粗い格子の開始インデクス
!! @param [in]  bk ブロック数
!<
  subroutine fb_interp_coarse0_s(dst, sz, g, src, st, bk)
  implicit none
  integer                                                      ::  i, j, k          ! 粗い格子のループインデクス
  integer                                                      ::  ii, jj, kk       ! 密な格子のループインデクス
  integer                                                      ::  ix, jx, kx, g, si, sj, sk
  integer, dimension(3)                                        ::  sz, st, bk
  real                                                         ::  q
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)       ::  dst
  real, dimension(1-g:sz(1)*bk(1)/2+g, 1-g:sz(2)*bk(2)/2+g, 1-g:sz(3)*bk(3)/2+g) ::  src

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
  si = st(1)
  sj = st(2)
  sk = st(3)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, si, sj, sk) &
!$OMP PRIVATE(q) &
!$OMP PRIVATE(ii, jj, kk)

!$OMP DO SCHEDULE(static)

  do k=sk, sk+kx/2-1
    kk = 2*k - (2*k-1)/kx * kx
  do j=sj, sj+jx/2-1
    jj = 2*j - (2*j-1)/jx * jx
  do i=si, si+ix/2-1
    ii = 2*i - (2*i-1)/ix * ix

    q = src(i  , j  , k  )

    dst(ii-1, jj-1, kk-1) = q
    dst(ii  , jj-1, kk-1) = q
    dst(ii-1, jj  , kk-1) = q
    dst(ii  , jj  , kk-1) = q
    dst(ii-1, jj-1, kk  ) = q
    dst(ii  , jj-1, kk  ) = q
    dst(ii-1, jj  , kk  ) = q
    dst(ii  , jj  , kk  ) = q

  end do
  end do
  end do

!$OMP END DO
!$OMP END PARALLEL

  return
  end subroutine fb_interp_coarse0_s


!> ********************************************************************
!! @brief 粗い格子から密な格子への補間（ゼロ次）
!! @param [out] dst 密な格子系
!! @param [in]  sz 配列長
!! @param [in]  g ガイドセル長
!! @param [in]  src 粗い格子系
!! @param [in]  st 粗い格子の開始インデクス
!! @param [in]  bk ブロック数
!<
  subroutine fb_interp_coarse0_v(dst, sz, g, src, st, bk)
  implicit none
  integer                                                         ::  i, j, k          ! 粗い格子のループインデクス
  integer                                                         ::  ii, jj, kk       ! 密な格子のループインデクス
  integer                                                         ::  ix, jx, kx, g, si, sj, sk
  integer, dimension(3)                                           ::  sz, st, bk
  real                                                            ::  q1, q2, q3
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)       ::  dst
  real, dimension(1-g:sz(1)*bk(1)/2+g, 1-g:sz(2)*bk(2)/2+g, 1-g:sz(3)*bk(3)/2+g, 3) ::  src


  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
  si = st(1)
  sj = st(2)
  sk = st(3)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, si, sj, sk) &
!$OMP PRIVATE(ii, jj, kk, q1, q2, q3)

!$OMP DO SCHEDULE(static)

  do k=sk, sk+kx/2-1
    kk = 2*k - (2*k-1)/kx * kx
  do j=sj, sj+jx/2-1
    jj = 2*j - (2*j-1)/jx * jx
  do i=si, si+ix/2-1
    ii = 2*i - (2*i-1)/ix * ix

    q1 = src(i  , j  , k  , 1)
    q2 = src(i  , j  , k  , 2)
    q3 = src(i  , j  , k  , 3)

    dst(ii-1, jj-1, kk-1, 1) = q1
    dst(ii  , jj-1, kk-1, 1) = q1
    dst(ii-1, jj  , kk-1, 1) = q1
    dst(ii  , jj  , kk-1, 1) = q1
    dst(ii-1, jj-1, kk  , 1) = q1
    dst(ii  , jj-1, kk  , 1) = q1
    dst(ii-1, jj  , kk  , 1) = q1
    dst(ii  , jj  , kk  , 1) = q1

    dst(ii-1, jj-1, kk-1, 2) = q2
    dst(ii  , jj-1, kk-1, 2) = q2
    dst(ii-1, jj  , kk-1, 2) = q2
    dst(ii  , jj  , kk-1, 2) = q2
    dst(ii-1, jj-1, kk  , 2) = q2
    dst(ii  , jj-1, kk  , 2) = q2
    dst(ii-1, jj  , kk  , 2) = q2
    dst(ii  , jj  , kk  , 2) = q2

    dst(ii-1, jj-1, kk-1, 3) = q3
    dst(ii  , jj-1, kk-1, 3) = q3
    dst(ii-1, jj  , kk-1, 3) = q3
    dst(ii  , jj  , kk-1, 3) = q3
    dst(ii-1, jj-1, kk  , 3) = q3
    dst(ii  , jj-1, kk  , 3) = q3
    dst(ii-1, jj  , kk  , 3) = q3
    dst(ii  , jj  , kk  , 3) = q3

  end do
  end do
  end do

!$OMP END DO
!$OMP END PARALLEL

  return
  end subroutine fb_interp_coarse0_v



!> ********************************************************************
!! @brief 粗い格子から密な格子への補間
!! @param [out] dst 密な格子系
!! @param [in]  sz 配列長
!! @param [in]  g ガイドセル長
!! @param [in]  src 粗い格子系
!! @param [in]  st 粗い格子の開始インデクス
!! @param [in]  bk ブロック数
!<
  subroutine fb_interp_coarse_s(dst, sz, g, src, st, bk)
  implicit none
  integer                                                      ::  i, j, k          ! 粗い格子のループインデクス
  integer                                                      ::  ii, jj, kk       ! 密な格子のループインデクス
  integer                                                      ::  ix, jx, kx, g, si, sj, sk
  integer, dimension(3)                                        ::  sz, st, bk
  real                                                         ::  g_x,  g_y,  g_z, q
  real                                                         ::  g_xx, g_yy, g_zz
  real                                                         ::  g_xy, g_yz, g_zx
  real                                                         ::         s_112,        s_121, s_122, s_123,        s_132
  real                                                         ::  s_211, s_212, s_213, s_221, s_222, s_223, s_231, s_232, s_233
  real                                                         ::         s_312,        s_321, s_322, s_323,        s_332
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)       ::  dst
  real, dimension(1-g:sz(1)*bk(1)/2+g, 1-g:sz(2)*bk(2)/2+g, 1-g:sz(3)*bk(3)/2+g) ::  src

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
  si = st(1)
  sj = st(2)
  sk = st(3)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, si, sj, sk) &
!$OMP PRIVATE(s_112, s_121, s_122, s_123, s_132) &
!$OMP PRIVATE(s_211, s_212, s_213, s_221, s_222, s_223, s_231, s_232, s_233) &
!$OMP PRIVATE(s_312, s_321, s_322, s_323, s_332) &
!$OMP PRIVATE(g_x, g_y, g_z, g_xx, g_yy, g_zz, g_xy, g_yz, g_zx, q) &
!$OMP PRIVATE(ii, jj, kk)

!$OMP DO SCHEDULE(static)

  do k=sk, sk+kx/2-1
    kk = 2*k - (2*k-1)/kx * kx
  do j=sj, sj+jx/2-1
    jj = 2*j - (2*j-1)/jx * jx
  do i=si, si+ix/2-1
    ii = 2*i - (2*i-1)/ix * ix

    s_112 = src(i-1, j-1, k  )
    s_121 = src(i-1, j  , k-1)
    s_122 = src(i-1, j  , k  )
    s_123 = src(i-1, j  , k+1)
    s_132 = src(i-1, j+1, k  )
        
    s_211 = src(i  , j-1, k-1)
    s_212 = src(i  , j-1, k  )
    s_213 = src(i  , j-1, k+1)
    s_221 = src(i  , j  , k-1)
    s_222 = src(i  , j  , k  )
    s_223 = src(i  , j  , k+1)
    s_231 = src(i  , j+1, k-1)
    s_232 = src(i  , j+1, k  )
    s_233 = src(i  , j+1, k+1)
            
    s_312 = src(i+1, j-1, k  )
    s_321 = src(i+1, j  , k-1)
    s_322 = src(i+1, j  , k  )
    s_323 = src(i+1, j  , k+1)
    s_332 = src(i+1, j+1, k  )

    g_x = ( s_322 - s_122 ) * 0.25
    g_y = ( s_232 - s_212 ) * 0.25
    g_z = ( s_223 - s_221 ) * 0.25

    g_xx = s_322 - 2.0 * s_222 + s_122
    g_yy = s_232 - 2.0 * s_222 + s_212
    g_zz = s_223 - 2.0 * s_222 + s_221

    g_xy = ( s_332 - s_312 - s_132 + s_112 ) * 0.25
    g_yz = ( s_233 - s_231 - s_213 + s_211 ) * 0.25
    g_zx = ( s_323 - s_123 - s_321 + s_121 ) * 0.25
    
    if (i == 1) then
      if (j == 1)  g_xy = s_332 - s_322 - s_232 + s_222
      if (j == jx) g_xy = s_322 - s_222 - s_312 + s_212
      
      if (k == 1)  g_zx = s_323 - s_322 - s_223 + s_222
      if (k == kx) g_zx = s_322 - s_321 - s_222 + s_221
    end if
    
    if (i == ix) then
      if (j == 1)  g_xy = s_232 - s_222 - s_132 + s_122
      if (j == jx) g_xy = s_222 - s_122 - s_212 + s_112
      
      if (k == 1)  g_zx = s_223 - s_222 - s_123 + s_122
      if (k == kx) g_zx = s_222 - s_221 - s_122 + s_121
    end if
    
    if (j == 1) then
      if (k == 1)  g_yz = s_233 - s_232 - s_223 + s_222
      if (k == kx) g_yz = s_232 - s_231 - s_222 + s_221
    end if

    if (j == jx) then
      if (k == 1)  g_yz = s_223 - s_222 - s_213 + s_212
      if (k == kx) g_yz = s_222 - s_221 - s_212 + s_211
    end if

    q = s_222 + 0.125 * ( g_xx + g_yy + g_zz )
    
    dst(ii-1, jj-1, kk-1) = q - g_x - g_y - g_z + 0.25 * (+ g_xy + g_yz + g_zx )
    dst(ii  , jj-1, kk-1) = q + g_x - g_y - g_z + 0.25 * (- g_xy + g_yz + g_zx )
    dst(ii-1, jj  , kk-1) = q - g_x + g_y - g_z + 0.25 * (- g_xy - g_yz + g_zx )
    dst(ii  , jj  , kk-1) = q + g_x + g_y - g_z + 0.25 * (+ g_xy - g_yz - g_zx )
    dst(ii-1, jj-1, kk  ) = q - g_x - g_y + g_z + 0.25 * (+ g_xy - g_yz - g_zx )
    dst(ii  , jj-1, kk  ) = q + g_x - g_y + g_z + 0.25 * (- g_xy - g_yz + g_zx )
    dst(ii-1, jj  , kk  ) = q - g_x + g_y + g_z + 0.25 * (- g_xy + g_yz - g_zx )
    dst(ii  , jj  , kk  ) = q + g_x + g_y + g_z + 0.25 * (+ g_xy + g_yz + g_zx )

  end do
  end do
  end do

!$OMP END DO
!$OMP END PARALLEL

  return
  end subroutine fb_interp_coarse_s

!> ********************************************************************
!! @brief 粗い格子から密な格子への補間
!! @param [out] dst 密な格子系
!! @param [in]  sz 配列長
!! @param [in]  g ガイドセル長
!! @param [in]  src 粗い格子系
!! @param [in]  st 粗い格子の開始インデクス
!! @param [in]  bk ブロック数
!<
  subroutine fb_interp_coarse_v(dst, sz, g, src, st, bk)
  implicit none
  integer                                                         ::  i, j, k          ! 粗い格子のループインデクス
  integer                                                         ::  ii, jj, kk       ! 密な格子のループインデクス
  integer                                                         ::  ix, jx, kx, g, si, sj, sk
  integer, dimension(3)                                           ::  sz, st, bk
  real                                                            ::  u_x, u_y, u_z, u_xx, u_yy, u_zz, u_xy, u_yz, u_zx
  real                                                            ::  v_x, v_y, v_z, v_xx, v_yy, v_zz, v_xy, v_yz, v_zx
  real                                                            ::  w_x, w_y, w_z, w_xx, w_yy, w_zz, w_xy, w_yz, w_zx
  real                                                            ::         s_112,        s_121, s_122, s_123,        s_132
  real                                                            ::  s_211, s_212, s_213, s_221, s_222, s_223, s_231, s_232, s_233
  real                                                            ::         s_312,        s_321, s_322, s_323,        s_332
  real                                                            ::  q
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3)       ::  dst
  real, dimension(1-g:sz(1)*bk(1)/2+g, 1-g:sz(2)*bk(2)/2+g, 1-g:sz(3)*bk(3)/2+g, 3) ::  src


  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
  si = st(1)
  sj = st(2)
  sk = st(3)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, si, sj, sk) &
!$OMP PRIVATE(s_112, s_121, s_122, s_123, s_132) &
!$OMP PRIVATE(s_211, s_212, s_213, s_221, s_222, s_223, s_231, s_232, s_233) &
!$OMP PRIVATE(s_312, s_321, s_322, s_323, s_332) &
!$OMP PRIVATE(u_x, u_y, u_z, u_xx, u_yy, u_zz, u_xy, u_yz, u_zx) &
!$OMP PRIVATE(v_x, v_y, v_z, v_xx, v_yy, v_zz, v_xy, v_yz, v_zx) &
!$OMP PRIVATE(w_x, w_y, w_z, w_xx, w_yy, w_zz, w_xy, w_yz, w_zx) &
!$OMP PRIVATE(ii, jj, kk, q)

!$OMP DO SCHEDULE(static)

  do k=sk, sk+kx/2-1
    kk = 2*k - (2*k-1)/kx * kx
  do j=sj, sj+jx/2-1
    jj = 2*j - (2*j-1)/jx * jx
  do i=si, si+ix/2-1
    ii = 2*i - (2*i-1)/ix * ix

!   u
    s_112 = src(i-1, j-1, k  , 1)
    s_121 = src(i-1, j  , k-1, 1)
    s_122 = src(i-1, j  , k  , 1)
    s_123 = src(i-1, j  , k+1, 1)
    s_132 = src(i-1, j+1, k  , 1)
        
    s_211 = src(i  , j-1, k-1, 1)
    s_212 = src(i  , j-1, k  , 1)
    s_213 = src(i  , j-1, k+1, 1)
    s_221 = src(i  , j  , k-1, 1)
    s_222 = src(i  , j  , k  , 1)
    s_223 = src(i  , j  , k+1, 1)
    s_231 = src(i  , j+1, k-1, 1)
    s_232 = src(i  , j+1, k  , 1)
    s_233 = src(i  , j+1, k+1, 1)
            
    s_312 = src(i+1, j-1, k  , 1)
    s_321 = src(i+1, j  , k-1, 1)
    s_322 = src(i+1, j  , k  , 1)
    s_323 = src(i+1, j  , k+1, 1)
    s_332 = src(i+1, j+1, k  , 1)

    u_x = ( s_322 - s_122 ) * 0.25
    u_y = ( s_232 - s_212 ) * 0.25
    u_z = ( s_223 - s_221 ) * 0.25

    u_xx = s_322 - 2.0 * s_222 + s_122
    u_yy = s_232 - 2.0 * s_222 + s_212
    u_zz = s_223 - 2.0 * s_222 + s_221

    u_xy = ( s_332 - s_312 - s_132 + s_112 ) * 0.25
    u_yz = ( s_233 - s_231 - s_213 + s_211 ) * 0.25
    u_zx = ( s_323 - s_123 - s_321 + s_121 ) * 0.25
    
    if (i == 1) then
      if (j == 1)  u_xy = s_332 - s_322 - s_232 + s_222
      if (j == jx) u_xy = s_322 - s_222 - s_312 + s_212
      
      if (k == 1)  u_zx = s_323 - s_322 - s_223 + s_222
      if (k == kx) u_zx = s_322 - s_321 - s_222 + s_221
    end if
    
    if (i == ix) then
      if (j == 1)  u_xy = s_232 - s_222 - s_132 + s_122
      if (j == jx) u_xy = s_222 - s_122 - s_212 + s_112
      
      if (k == 1)  u_zx = s_223 - s_222 - s_123 + s_122
      if (k == kx) u_zx = s_222 - s_221 - s_122 + s_121
    end if
    
    if (j == 1) then
      if (k == 1)  u_yz = s_233 - s_232 - s_223 + s_222
      if (k == kx) u_yz = s_232 - s_231 - s_222 + s_221
    end if

    if (j == jx) then
      if (k == 1)  u_yz = s_223 - s_222 - s_213 + s_212
      if (k == kx) u_yz = s_222 - s_221 - s_212 + s_211
    end if

    q = s_222 + 0.125 * ( u_xx + u_yy + u_zz )
    
    dst(ii-1, jj-1, kk-1, 1) = q - u_x - u_y - u_z + 0.25 * (+ u_xy + u_yz + u_zx )
    dst(ii  , jj-1, kk-1, 1) = q + u_x - u_y - u_z + 0.25 * (- u_xy + u_yz + u_zx )
    dst(ii-1, jj  , kk-1, 1) = q - u_x + u_y - u_z + 0.25 * (- u_xy - u_yz + u_zx )
    dst(ii  , jj  , kk-1, 1) = q + u_x + u_y - u_z + 0.25 * (+ u_xy - u_yz - u_zx )
    dst(ii-1, jj-1, kk  , 1) = q - u_x - u_y + u_z + 0.25 * (+ u_xy - u_yz - u_zx )
    dst(ii  , jj-1, kk  , 1) = q + u_x - u_y + u_z + 0.25 * (- u_xy - u_yz + u_zx )
    dst(ii-1, jj  , kk  , 1) = q - u_x + u_y + u_z + 0.25 * (- u_xy + u_yz - u_zx )
    dst(ii  , jj  , kk  , 1) = q + u_x + u_y + u_z + 0.25 * (+ u_xy + u_yz + u_zx )


!   v
    s_112 = src(i-1, j-1, k  , 2)
    s_121 = src(i-1, j  , k-1, 2)
    s_122 = src(i-1, j  , k  , 2)
    s_123 = src(i-1, j  , k+1, 2)
    s_132 = src(i-1, j+1, k  , 2)
        
    s_211 = src(i  , j-1, k-1, 2)
    s_212 = src(i  , j-1, k  , 2)
    s_213 = src(i  , j-1, k+1, 2)
    s_221 = src(i  , j  , k-1, 2)
    s_222 = src(i  , j  , k  , 2)
    s_223 = src(i  , j  , k+1, 2)
    s_231 = src(i  , j+1, k-1, 2)
    s_232 = src(i  , j+1, k  , 2)
    s_233 = src(i  , j+1, k+1, 2)
            
    s_312 = src(i+1, j-1, k  , 2)
    s_321 = src(i+1, j  , k-1, 2)
    s_322 = src(i+1, j  , k  , 2)
    s_323 = src(i+1, j  , k+1, 2)
    s_332 = src(i+1, j+1, k  , 2)

    v_x = ( s_322 - s_122 ) * 0.25
    v_y = ( s_232 - s_212 ) * 0.25
    v_z = ( s_223 - s_221 ) * 0.25

    v_xx = s_322 - 2.0 * s_222 + s_122
    v_yy = s_232 - 2.0 * s_222 + s_212
    v_zz = s_223 - 2.0 * s_222 + s_221

    v_xy = ( s_332 - s_312 - s_132 + s_112 ) * 0.25
    v_yz = ( s_233 - s_231 - s_213 + s_211 ) * 0.25
    v_zx = ( s_323 - s_123 - s_321 + s_121 ) * 0.25
    
    if (i == 1) then
      if (j == 1)  v_xy = s_332 - s_322 - s_232 + s_222
      if (j == jx) v_xy = s_322 - s_222 - s_312 + s_212
      
      if (k == 1)  v_zx = s_323 - s_322 - s_223 + s_222
      if (k == kx) v_zx = s_322 - s_321 - s_222 + s_221
    end if
    
    if (i == ix) then
      if (j == 1)  v_xy = s_232 - s_222 - s_132 + s_122
      if (j == jx) v_xy = s_222 - s_122 - s_212 + s_112
      
      if (k == 1)  v_zx = s_223 - s_222 - s_123 + s_122
      if (k == kx) v_zx = s_222 - s_221 - s_122 + s_121
    end if
    
    if (j == 1) then
      if (k == 1)  v_yz = s_233 - s_232 - s_223 + s_222
      if (k == kx) v_yz = s_232 - s_231 - s_222 + s_221
    end if

    if (j == jx) then
      if (k == 1)  v_yz = s_223 - s_222 - s_213 + s_212
      if (k == kx) v_yz = s_222 - s_221 - s_212 + s_211
    end if

    q = s_222 + 0.125 * ( v_xx + v_yy + v_zz )
    
    dst(ii-1, jj-1, kk-1, 2) = q - v_x - v_y - v_z + 0.25 * (+ v_xy + v_yz + v_zx )
    dst(ii  , jj-1, kk-1, 2) = q + v_x - v_y - v_z + 0.25 * (- v_xy + v_yz + v_zx )
    dst(ii-1, jj  , kk-1, 2) = q - v_x + v_y - v_z + 0.25 * (- v_xy - v_yz + v_zx )
    dst(ii  , jj  , kk-1, 2) = q + v_x + v_y - v_z + 0.25 * (+ v_xy - v_yz - v_zx )
    dst(ii-1, jj-1, kk  , 2) = q - v_x - v_y + v_z + 0.25 * (+ v_xy - v_yz - v_zx )
    dst(ii  , jj-1, kk  , 2) = q + v_x - v_y + v_z + 0.25 * (- v_xy - v_yz + v_zx )
    dst(ii-1, jj  , kk  , 2) = q - v_x + v_y + v_z + 0.25 * (- v_xy + v_yz - v_zx )
    dst(ii  , jj  , kk  , 2) = q + v_x + v_y + v_z + 0.25 * (+ v_xy + v_yz + v_zx )


!   w
    s_112 = src(i-1, j-1, k  , 3)
    s_121 = src(i-1, j  , k-1, 3)
    s_122 = src(i-1, j  , k  , 3)
    s_123 = src(i-1, j  , k+1, 3)
    s_132 = src(i-1, j+1, k  , 3)
        
    s_211 = src(i  , j-1, k-1, 3)
    s_212 = src(i  , j-1, k  , 3)
    s_213 = src(i  , j-1, k+1, 3)
    s_221 = src(i  , j  , k-1, 3)
    s_222 = src(i  , j  , k  , 3)
    s_223 = src(i  , j  , k+1, 3)
    s_231 = src(i  , j+1, k-1, 3)
    s_232 = src(i  , j+1, k  , 3)
    s_233 = src(i  , j+1, k+1, 3)
            
    s_312 = src(i+1, j-1, k  , 3)
    s_321 = src(i+1, j  , k-1, 3)
    s_322 = src(i+1, j  , k  , 3)
    s_323 = src(i+1, j  , k+1, 3)
    s_332 = src(i+1, j+1, k  , 3)

    w_x = ( s_322 - s_122 ) * 0.25
    w_y = ( s_232 - s_212 ) * 0.25
    w_z = ( s_223 - s_221 ) * 0.25

    w_xx = s_322 - 2.0 * s_222 + s_122
    w_yy = s_232 - 2.0 * s_222 + s_212
    w_zz = s_223 - 2.0 * s_222 + s_221

    w_xy = ( s_332 - s_312 - s_132 + s_112 ) * 0.25
    w_yz = ( s_233 - s_231 - s_213 + s_211 ) * 0.25
    w_zx = ( s_323 - s_123 - s_321 + s_121 ) * 0.25
    
    if (i == 1) then
      if (j == 1)  w_xy = s_332 - s_322 - s_232 + s_222
      if (j == jx) w_xy = s_322 - s_222 - s_312 + s_212
      
      if (k == 1)  w_zx = s_323 - s_322 - s_223 + s_222
      if (k == kx) w_zx = s_322 - s_321 - s_222 + s_221
    end if
    
    if (i == ix) then
      if (j == 1)  w_xy = s_232 - s_222 - s_132 + s_122
      if (j == jx) w_xy = s_222 - s_122 - s_212 + s_112
      
      if (k == 1)  w_zx = s_223 - s_222 - s_123 + s_122
      if (k == kx) w_zx = s_222 - s_221 - s_122 + s_121
    end if
    
    if (j == 1) then
      if (k == 1)  w_yz = s_233 - s_232 - s_223 + s_222
      if (k == kx) w_yz = s_232 - s_231 - s_222 + s_221
    end if

    if (j == jx) then
      if (k == 1)  w_yz = s_223 - s_222 - s_213 + s_212
      if (k == kx) w_yz = s_222 - s_221 - s_212 + s_211
    end if

    q = s_222 + 0.125 * ( w_xx + w_yy + w_zz )
    
    dst(ii-1, jj-1, kk-1, 3) = q - w_x - w_y - w_z + 0.25 * (+ w_xy + w_yz + w_zx )
    dst(ii  , jj-1, kk-1, 3) = q + w_x - w_y - w_z + 0.25 * (- w_xy + w_yz + w_zx )
    dst(ii-1, jj  , kk-1, 3) = q - w_x + w_y - w_z + 0.25 * (- w_xy - w_yz + w_zx )
    dst(ii  , jj  , kk-1, 3) = q + w_x + w_y - w_z + 0.25 * (+ w_xy - w_yz - w_zx )
    dst(ii-1, jj-1, kk  , 3) = q - w_x - w_y + w_z + 0.25 * (+ w_xy - w_yz - w_zx )
    dst(ii  , jj-1, kk  , 3) = q + w_x - w_y + w_z + 0.25 * (- w_xy - w_yz + w_zx )
    dst(ii-1, jj  , kk  , 3) = q - w_x + w_y + w_z + 0.25 * (- w_xy + w_yz - w_zx )
    dst(ii  , jj  , kk  , 3) = q + w_x + w_y + w_z + 0.25 * (+ w_xy + w_yz + w_zx )

  end do
  end do
  end do

!$OMP END DO
!$OMP END PARALLEL

  return
  end subroutine fb_interp_coarse_v

  
  
!> ********************************************************************
!! @brief 速度ベクトルの配列の添え字変換
!! @param [out]    vo    変換されたベクトル（平均場の場合は積算値）
!! @param [in]     sz    配列長
!! @param [in]     g     ガイドセル長
!! @param [in]     vi    入力速度場
!! @param [in]     v00   参照速度
!! @param [in]     refv  代表速度
!! @param [out]    flop  浮動小数演算数
!! @note dst[] = src[]/refv + v00 有次元のときrefvは次元速度，無次元のとき1.0
!<
  subroutine fb_vin_nijk (vo, sz, g, vi, v00, refv, flop)
  implicit none
  integer                                                   ::  i, j, k, ix, jx, kx, g
  integer, dimension(3)                                     ::  sz
  real                                                      ::  u_ref, v_ref, w_ref, refv, rr
  double precision                                          ::  flop
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  vo
  real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  vi
  real, dimension(0:3)                                      ::  v00

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
  
  u_ref = v00(1)
  v_ref = v00(2)
  w_ref = v00(3)
  
  rr = 1.0/refv
  
  flop = flop + dble(ix)*dble(jx)*dble(kx)*9.0d0 + 8.0d0


!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, u_ref, v_ref, w_ref, rr)

!$OMP DO SCHEDULE(static)

  do k=1,kx
  do j=1,jx
  do i=1,ix 
    vo(i,j,k,1) = vi(1,i,j,k) * rr + u_ref
    vo(i,j,k,2) = vi(2,i,j,k) * rr + v_ref
    vo(i,j,k,3) = vi(3,i,j,k) * rr + w_ref
  end do
  end do
  end do
  
!$OMP END DO
!$OMP END PARALLEL

  return
  end subroutine fb_vin_nijk
  
!> ********************************************************************
!! @brief 速度ベクトルの添え字変換をして，scale倍する
!! @param [out] vout   変換されたベクトル
!! @param [in]  vin    変換前
!! @param [in]  sz     配列長
!! @param [in]  g      ガイドセル長
!! @param [in]  v00    参照速度
!! @param [in]  unit_v 無次元のとき1.0，有次元のとき代表速度(m/s)
!! @param [out] flop   浮動小数演算数
!! @note dst[] = ( src[] * stepAvr ) - v00
!<
  subroutine fb_vout_nijk (vout, vin, sz, g, v00, unit_v, flop)
  implicit none
  integer                                                   ::  i, j, k, ix, jx, kx, g
  integer, dimension(3)                                     ::  sz
  real                                                      ::  u_ref, v_ref, w_ref, unit_v, unit
  double precision                                          ::  flop
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  vin
  real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  vout
  real, dimension(0:3)                                      ::  v00

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
  
  u_ref = v00(1)
  v_ref = v00(2)
  w_ref = v00(3)

  unit = unit_v
  
  flop = flop + dble(ix)*dble(jx)*dble(kx)*9.0d0


!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, u_ref, v_ref, w_ref, unit)

!$OMP DO SCHEDULE(static)

  do k=1,kx
  do j=1,jx
  do i=1,ix 
    vout(1,i,j,k) = ( vin(i,j,k,1) - u_ref ) * unit
    vout(2,i,j,k) = ( vin(i,j,k,2) - v_ref ) * unit
    vout(3,i,j,k) = ( vin(i,j,k,3) - w_ref ) * unit
  end do
  end do
  end do
  
!$OMP END DO
!$OMP END PARALLEL

  return
  end subroutine fb_vout_nijk

!> ********************************************************************
!! @brief 速度ベクトルをscale倍する
!! @param [out] vout   変換されたベクトル
!! @param [in]  vin    変換前
!! @param [in]  sz     配列長
!! @param [in]  g      ガイドセル長
!! @param [in]  v00    参照速度
!! @param [in]  unit_v 無次元のとき1.0，有次元のとき代表速度(m/s)
!! @param [out] flop   浮動小数演算数
!! @note dst[] = ( src[] * stepAvr ) - v00
!<
  subroutine fb_vout_ijkn (vout, vin, sz, g, v00, unit_v, flop)
  implicit none
  integer                                                   ::  i, j, k, ix, jx, kx, g
  integer, dimension(3)                                     ::  sz
  real                                                      ::  u_ref, v_ref, w_ref, unit_v, unit
  double precision                                          ::  flop
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  vin, vout
  real, dimension(0:3)                                      ::  v00

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)

  u_ref = v00(1)
  v_ref = v00(2)
  w_ref = v00(3)

  unit = unit_v

  flop = flop + dble(ix)*dble(jx)*dble(kx)*9.0d0


!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, u_ref, v_ref, w_ref, unit)

!$OMP DO SCHEDULE(static)

  do k=1,kx
  do j=1,jx
  do i=1,ix
    vout(i,j,k,1) = ( vin(i,j,k,1) - u_ref ) * unit
    vout(i,j,k,2) = ( vin(i,j,k,2) - v_ref ) * unit
    vout(i,j,k,3) = ( vin(i,j,k,3) - w_ref ) * unit
  end do
  end do
  end do

!$OMP END DO
!$OMP END PARALLEL

  return
  end subroutine fb_vout_ijkn


!> ********************************************************************
!! @brief 有効セルに対する，1タイムステップ進行時のベクトルの絶対値の変化量の和のみ
!! @param [out] d    戻り値（変化量の2乗和）
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル長
!! @param [in]  vn   ベクトル値 n+1 step
!! @param [in]  vo   ベクトル値 n step
!! @param [in]  bx   BCindex
!! @param [out] flop 浮動小数演算数
!<
  subroutine fb_delta_v (d, sz, g, vn, vo, bx, flop)
  implicit none
  include 'ffv_f_params.h'
  integer                                                   ::  i, j, k, ix, jx, kx, g
  integer, dimension(3)                                     ::  sz
  double precision                                          ::  u, v, w, x, y, z, actv
  double precision                                          ::  flop, av, rm
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  vn, vo
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bx
  double precision, dimension(2)                            ::  d

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)

  flop = flop + dble(ix)*dble(jx)*dble(kx)*58.0d0
  ! sqrt float->10, double->20

  av = 0.0
  rm = 0.0

!$OMP PARALLEL &
!$OMP REDUCTION(+:av) &
!$OMP REDUCTION(+:rm) &
!$OMP PRIVATE(actv, u, v, w, x, y, z) &
!$OMP FIRSTPRIVATE(ix, jx, kx)

!$OMP DO SCHEDULE(static)

  do k=1,kx
  do j=1,jx
  do i=1,ix

    actv = real(ibits(bx(i,j,k), State, 1))
    
    !u = dble(vn(i,j,k,1))
    !v = dble(vn(i,j,k,2))
    !w = dble(vn(i,j,k,3))
    !av = av + sqrt(u*u + v*v + w*w)*actv
    
    x = u - dble(vo(i,j,k,1))
    y = v - dble(vo(i,j,k,2))
    z = w - dble(vo(i,j,k,3))
    rm = rm + sqrt(x*x + y*y + z*z)*actv

  end do
  end do
  end do
!$OMP END DO
!$OMP END PARALLEL

  d(1) = rm
  !d(2) = av

  return
  end subroutine fb_delta_v



!> ********************************************************************
!! @brief 有効セルに対する，1タイムステップ進行時の変化量の2乗和と値の和
!! @param [out] d    戻り値（変化量の2乗和と値の和）
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル長
!! @param [in]  sn   スカラー値 n+1 step
!! @param [in]  so   スカラー値 n step
!! @param [in]  bx   BCindex
!! @param [out] flop 浮動小数演算数
!<
  subroutine fb_delta_s (d, sz, g, sn, so, bx, flop)
  implicit none
  include 'ffv_f_params.h'
  integer                                                   ::  i, j, k, ix, jx, kx, g
  integer, dimension(3)                                     ::  sz
  double precision                                          ::  a, s, av, rm, actv, flop
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  sn, so
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bx
  double precision, dimension(2)                            ::  d

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)

  flop = flop + dble(ix)*dble(jx)*dble(kx)*7.0d0

  av = 0.0
  rm = 0.0

!$OMP PARALLEL &
!$OMP REDUCTION(+:av) &
!$OMP REDUCTION(+:rm) &
!$OMP PRIVATE(actv, s, a) &
!$OMP FIRSTPRIVATE(ix, jx, kx)

!$OMP DO SCHEDULE(static)

  do k=1,kx
  do j=1,jx
  do i=1,ix
    actv = dble(ibits(bx(i,j,k), State,  1))
    
    s = dble(sn(i,j,k))
    av = av + s * actv
    
    a = ( s - dble(so(i,j,k)) )*actv
    rm = rm + a*a
  end do
  end do
  end do
!$OMP END DO
!$OMP END PARALLEL
  
  d(1) = rm
  d(2) = av


  return
  end subroutine fb_delta_s



!> ********************************************************************
!! @brief ベクトル値を加算する
!! @param [in,out] avr  平均値
!! @param [in]     sz   配列長
!! @param [in]     g    ガイドセル長
!! @param [in]     v    ベクトル値
!! @param [in]     nadd 加算回数
!! @param [out]    flop 浮動小数演算数
!<
  subroutine fb_average_v (avr, sz, g, v, nadd, flop)
  implicit none
  integer                                                   ::  i, j, k, ix, jx, kx, g
  integer, dimension(3)                                     ::  sz
  double precision                                          ::  flop
  real                                                      ::  nadd, val1, val2
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v, avr

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)

  flop = flop + dble(ix)*dble(jx)*dble(kx)*9.0d0

  val2 = 1.0/nadd
  val1 = 1.0 - val2

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx, val1, val2)

!$OMP DO SCHEDULE(static)

  do k=1,kx
  do j=1,jx
  do i=1,ix
    avr(i,j,k,1) = val1 * avr(i,j,k,1) + val2 * v(i,j,k,1)
    avr(i,j,k,2) = val1 * avr(i,j,k,2) + val2 * v(i,j,k,2)
    avr(i,j,k,3) = val1 * avr(i,j,k,3) + val2 * v(i,j,k,3)
  end do
  end do
  end do
!$OMP END DO
!$OMP END PARALLEL

  return
  end subroutine fb_average_v



!> ********************************************************************
!! @brief スカラ値を加算する
!! @param [in,out] avr  平均値
!! @param [in]     sz   配列長
!! @param [in]     g    ガイドセル長
!! @param [in]     s    スカラ値
!! @param [in]     nadd 加算回数
!! @param [out]    flop 浮動小数演算数
!<
  subroutine fb_average_s (avr, sz, g, s, nadd, flop)
  implicit none
  integer                                                   ::  i, j, k, ix, jx, kx, g
  integer, dimension(3)                                     ::  sz
  double precision                                          ::  flop
  real                                                      ::  nadd, val1, val2
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  s, avr

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)

  flop = flop + dble(ix)*dble(jx)*dble(kx)*3.0d0

  val2 = 1.0/nadd
  val1 = 1.0 - val2

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(ix, jx, kx)

!$OMP DO SCHEDULE(static)

  do k=1,kx
  do j=1,jx
  do i=1,ix
    avr(i,j,k) = val1 * avr(i,j,k) + val2 * s(i,j,k)
  end do
  end do
  end do
!$OMP END DO
!$OMP END PARALLEL

  return
  end subroutine fb_average_s



!> ********************************************************************
!! @brief 全圧を計算する
!! @param [out] tp   全圧
!! @param [in]  sz   配列長
!! @param [in]  g    ガイドセル長
!! @param [in]  v    速度ベクトル
!! @param [in]  p    圧力
!! @param [in]  v00  参照速度
!! @param [out] flop 浮動小数演算数
!<
    subroutine fb_totalp (tp, sz, g, v, p, v00, flop)
    implicit none
    integer                                                   ::  i, j, k, ix, jx, kx, g
    integer, dimension(3)                                     ::  sz
    real                                                      ::  u1, u2, u3, vx, vy, vz
    double precision                                          ::  flop
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  p, tp
    real, dimension(0:3)                                      ::  v00

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    vx = v00(1)
    vy = v00(2)
    vz = v00(3)

!$OMP PARALLEL &
!$OMP PRIVATE(u1, u2, u3) &
!$OMP FIRSTPRIVATE(ix, jx, kx, vx, vy, vz)

!$OMP DO SCHEDULE(static)

    do k=1,kx
    do j=1,jx
    do i=1,ix
      u1 = v(i,j,k,1) - vx
      u2 = v(i,j,k,2) - vy
      u3 = v(i,j,k,3) - vz
      tp(i,j,k) = 0.5*(u1*u1 + u2*u2 + u3*u3) + p(i,j,k)
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL
        
    flop = flop + dble(ix)*dble(jx)*dble(kx)*10.0d0
    
    return
    end subroutine fb_totalp



!> ********************************************************************
!! @brief 速度の最小値と最大値を計算する (shape:ijkn)
!! @param [out] v_min 最小値
!! @param [out] v_max 最大値
!! @param [in]  sz    配列長
!! @param [in]  g     ガイドセル長
!! @param [in]  v00   参照速度
!! @param [in]  v     速度ベクトル
!! @param [out] flop  浮動小数演算数
!<
    subroutine fb_minmax_v (v_min, v_max, sz, g, v00, v)
    implicit none
    integer                                                   ::  i, j, k, ix, jx, kx, g
    integer, dimension(3)                                     ::  sz
    real                                                      ::  u1, u2, u3, uu, flop, vx, vy, vz
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  v
    real, dimension(0:3)                                      ::  v00, v_min, v_max

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    vx = v00(1)
    vy = v00(2)
    vz = v00(3)
    v_min(0) =  1.0e6
    v_min(1) =  1.0e6
    v_min(2) =  1.0e6
    v_min(3) =  1.0e6
    v_max(0) = -1.0e6
    v_max(1) = -1.0e6
    v_max(2) = -1.0e6
    v_max(3) = -1.0e6
    
    ! 16 + sqrt*1 = 26 ! DP 36
    flop = flop + real(ix)*real(jx)*real(kx)*26.0d0
    ! flop = flop + real(ix)*real(jx)*real(kx)*30.0 ! DP

!$OMP PARALLEL &
!$OMP REDUCTION(min:v_min) &
!$OMP REDUCTION(max:v_max) &
!$OMP PRIVATE(u1, u2, u3, uu) &
!$OMP FIRSTPRIVATE(ix, jx, kx, vx, vy, vz)

!$OMP DO SCHEDULE(static)

    do k=1,kx
    do j=1,jx
    do i=1,ix
			u1 = v(i,j,k,1)-vx
			u2 = v(i,j,k,2)-vy
			u3 = v(i,j,k,3)-vz
			uu = sqrt( u1*u1 + u2*u2 + u3*u3 )
      v_min(0) = min(v_min(0), uu)
      v_max(0) = max(v_max(0), uu)
      v_min(1) = min(v_min(1), u1)
      v_max(1) = max(v_max(1), u1)
      v_min(2) = min(v_min(2), u2)
      v_max(2) = max(v_max(2), u2)
      v_min(3) = min(v_min(3), u3)
      v_max(3) = max(v_max(3), u3)
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    return
    end subroutine fb_minmax_v

!> ********************************************************************
!! @brief 速度の最小値と最大値を計算する (shape:nijk)
!! @param [out] v_min 最小値
!! @param [out] v_max 最大値
!! @param [in]  sz    配列長
!! @param [in]  g     ガイドセル長
!! @param [in]  v00   参照速度
!! @param [in]  v     速度ベクトル
!! @param [out] flop  浮動小数演算数
!<
  subroutine fb_minmax_vex (v_min, v_max, sz, g, v00, v)
  implicit none
  integer                                                   ::  i, j, k, ix, jx, kx, g
  integer, dimension(3)                                     ::  sz
  real                                                      ::  u1, u2, u3, uu, flop, vx, vy, vz
  real, dimension(3, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  v
  real, dimension(0:3)                                      ::  v00, v_min, v_max

  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
  vx = v00(1)
  vy = v00(2)
  vz = v00(3)
  v_min(0) =  1.0e6
  v_min(1) =  1.0e6
  v_min(2) =  1.0e6
  v_min(3) =  1.0e6
  v_max(0) = -1.0e6
  v_max(1) = -1.0e6
  v_max(2) = -1.0e6
  v_max(3) = -1.0e6

  ! 10 + sqrt*1 = 20 ! DP 30
  flop = flop + real(ix)*real(jx)*real(kx)*26.0d0
  ! flop = flop + real(ix)*real(jx)*real(kx)*30.0 ! DP

!$OMP PARALLEL &
!$OMP REDUCTION(min:v_min) &
!$OMP REDUCTION(max:v_max) &
!$OMP PRIVATE(u1, u2, u3, uu) &
!$OMP FIRSTPRIVATE(ix, jx, kx, vx, vy, vz)

!$OMP DO SCHEDULE(static)

  do k=1,kx
  do j=1,jx
  do i=1,ix
    u1 = v(1,i,j,k)-vx
    u2 = v(2,i,j,k)-vy
    u3 = v(3,i,j,k)-vz
    uu = sqrt( u1*u1 + u2*u2 + u3*u3 )
    v_min(0) = min(v_min(0), uu)
    v_max(0) = max(v_max(0), uu)
    v_min(1) = min(v_min(1), u1)
    v_max(1) = max(v_max(1), u1)
    v_min(2) = min(v_min(2), u2)
    v_max(2) = max(v_max(2), u2)
    v_min(3) = min(v_min(3), u3)
    v_max(3) = max(v_max(3), u3)
  end do
  end do
  end do
!$OMP END DO
!$OMP END PARALLEL

  return
  end subroutine fb_minmax_vex


!> ********************************************************************
!! @brief スカラ値の最小値と最大値を計算する
!! @param [out] f_min 最小値
!! @param [out] f_max 最大値
!! @param [in]  sz    配列長
!! @param [in]  g     ガイドセル長
!! @param [in]  s     スカラ値
!! @param [out] flop  浮動小数演算数
!<
    subroutine fb_minmax_s (f_min, f_max, sz, g, s, flop)
    implicit none
    integer                                                   ::  i, j, k, ix, jx, kx, g
    integer, dimension(3)                                     ::  sz
    real                                                      ::  f_min, f_max
    double precision                                          ::  flop
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  s

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    f_min =  1.0e6
    f_max = -1.0e6
    
    flop = flop + dble(ix)*dble(jx)*dble(kx)*2.0d0

!$OMP PARALLEL &
!$OMP REDUCTION(min:f_min) &
!$OMP REDUCTION(max:f_max) &
!$OMP FIRSTPRIVATE(ix, jx, kx)

!$OMP DO SCHEDULE(static)

    do k=1,kx
    do j=1,jx
    do i=1,ix
      f_min = min(f_min, s(i,j,k))
      f_max = max(f_max, s(i,j,k))
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    return
    end subroutine fb_minmax_s

!> ********************************************************************
!! @brief ベクトル値を設定する
!! @param [out] var ベクトル配列
!! @param [in]  sz  配列長
!! @param [in]  g   ガイドセル長
!! @param [in]  val ベクトル値
!! @param [in]  bv  BCindex C
!<
    subroutine fb_set_vector (var, sz, g, val, bv)
    implicit none
    include '../FB/ffv_f_params.h'
    integer                                                   ::  i, j, k, ix, jx, kx, g, bvx
    integer, dimension(3)                                     ::  sz
    real                                                      ::  u1, u2, u3
    real, dimension(3)                                        ::  val
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  var
    integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bv

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)
    u1 = val(1)
    u2 = val(2)
    u3 = val(3)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(u1, u2, u3) &
!$OMP FIRSTPRIVATE(ix, jx, kx, g) &
!$OMP PRIVATE(bvx)

!$OMP DO SCHEDULE(static)

    do k=1-g, kx+g
    do j=1-g, jx+g
    do i=1-g, ix+g
        bvx = ibits(bv(i,j,k), State, 1)
        var(i,j,k,1) = u1 * real(bvx)
        var(i,j,k,2) = u2 * real(bvx)
        var(i,j,k,3) = u3 * real(bvx)
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    return
    end subroutine fb_set_vector

!> ********************************************************************
!! @brief セルフェイスベクトル値を設定する
!! @param [out] var ベクトル配列
!! @param [in]  sz  配列長
!! @param [in]  g   ガイドセル長
!! @param [in]  val ベクトル値
!! @param [in]  bv  BCindex C
!<
subroutine fb_set_fvector (var, sz, g, val, bv)
implicit none
include '../FB/ffv_f_params.h'
integer                                                   ::  i, j, k, ix, jx, kx, g, bvx
integer, dimension(3)                                     ::  sz
real                                                      ::  u1, u2, u3
real, dimension(3)                                        ::  val
real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 3) ::  var
integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g) ::  bv

ix = sz(1)
jx = sz(2)
kx = sz(3)
u1 = val(1)
u2 = val(2)
u3 = val(3)

!$OMP PARALLEL &
!$OMP FIRSTPRIVATE(u1, u2, u3) &
!$OMP FIRSTPRIVATE(ix, jx, kx) &
!$OMP PRIVATE(bvx)

!$OMP DO SCHEDULE(static)

do k=0, kx
do j=0, jx
do i=0, ix
var(i,j,k,1) = u1 * real( ibits(bv(i,j,k), State, 1) * ibits(bv(i+1,j,k), State, 1) )
var(i,j,k,2) = u2 * real( ibits(bv(i,j,k), State, 1) * ibits(bv(i,j+1,k), State, 1) )
var(i,j,k,3) = u3 * real( ibits(bv(i,j,k), State, 1) * ibits(bv(i,j,k+1), State, 1) )
end do
end do
end do
!$OMP END DO
!$OMP END PARALLEL

return
end subroutine fb_set_fvector



!> ********************************************************************
!! @brief スカラ値の値を[0, 1]に制限する
!! @param [in,out] t  スカラ値
!! @param [in]     sz 配列長
!! @param [in]     g  ガイドセル長
!<
    subroutine fb_limit_scalar (t, sz, g)
    implicit none
    integer                                                   ::  i, j, k, ix, jx, kx, g
    integer, dimension(3)                                     ::  sz
    real                                                      ::  tmp
    real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)    ::  t

    ix = sz(1)
    jx = sz(2)
    kx = sz(3)

!$OMP PARALLEL &
!$OMP PRIVATE(tmp) &
!$OMP FIRSTPRIVATE(ix, jx, kx)

!$OMP DO SCHEDULE(static)

    do k=1,kx
    do j=1,jx
    do i=1,ix
      tmp = t(i,j,k)
      if (tmp<0.0) t(i,j,k)=0.0
      if (tmp>1.0) t(i,j,k)=1.0
    end do
    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    return
    end subroutine fb_limit_scalar
