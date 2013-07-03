!###################################################################################
!
! FFV-C
! Frontflow / violet Cartesian
!
!
! Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
! All rights reserved.
!
! Copyright (c) 2011-2013 Institute of Industrial Science, The University of Tokyo.
! All rights reserved.
!
! Copyright (c) 2012-2013 Advanced Institute for Computational Science, RIKEN.
! All rights reserved.
!
!###################################################################################

!> @file   muscl.h
!! @brief  MUSCL reconstruction core
!! @author kero
!<

! b  = (3.0-ck)/(1.0-ck), ck=1/3
! w_?; セル界面フラグ　(0-wall face / 1-fluid)
! 6*6 = 36 flops 

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
