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
! Copyright (c) 2016-2020 Research Institute for Information Research, Kyushu University.
! All rights reserved.
!
!###################################################################################

!> @file   core_psor.h
!! @brief  point sorのコア anisotropic Cartesian
!! @author aics
!<

! 57 flop

idx = bp(i,j,k)

c_w = real(ibits(idx, bc_ndag_W, 1))  ! w
c_e = real(ibits(idx, bc_ndag_E, 1))  ! e
c_s = real(ibits(idx, bc_ndag_S, 1))  ! s
c_n = real(ibits(idx, bc_ndag_N, 1))  ! n
c_b = real(ibits(idx, bc_ndag_B, 1))  ! b
c_t = real(ibits(idx, bc_ndag_T, 1))  ! t

d_w = real(ibits(idx, bc_dn_W, 1))
d_e = real(ibits(idx, bc_dn_E, 1))
d_s = real(ibits(idx, bc_dn_S, 1))
d_n = real(ibits(idx, bc_dn_N, 1))
d_b = real(ibits(idx, bc_dn_B, 1))
d_t = real(ibits(idx, bc_dn_T, 1))

dsw = real(ibits(idx, bc_diag, 1))

dd = c_w + c_e &
   + c_s + c_n &
   + c_b + c_t &
   + 2.0       &
   *(d_w + d_e &
   + d_s + d_n &
   + d_b + d_t ) &
   + cf

dd = dsw * dd + 1.0 - dsw

aa = dble(ibits(idx, Active, 1))

pp = p(i,j,k)
bb = b(i,j,k)

ss = c_e * p(i+1,j  ,k  ) + c_w * p(i-1,j  ,k  ) &
   + c_n * p(i  ,j+1,k  ) + c_s * p(i  ,j-1,k  ) &
   + c_t * p(i  ,j  ,k+1) + c_b * p(i  ,j  ,k-1)

dp = ( (ss - bb)/dd - pp ) * omg
pn = pp + dp
p(i,j,k) = pn

de  = bb - (ss - pn * dd)
res = res + dble(de*de) * aa
xl2 = xl2 + dble(pn*pn) * aa
err = err + dble(dp*dp) * aa
