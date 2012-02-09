!   *********************************************************
!
!   SPHERE - Skeleton for PHysical and Engineering REsearch
!  
!   Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
!
!   *********************************************************

!> @file muscl_x.h
!> @brief subroutines for CBC
!> @author keno, FSI Team, VCAD, RIKEN 

! b  = (3.0-ck)/(1.0-ck), ck=1/3
! w_?; セル界面フラグ　(0-wall face / 1-fluid)
! 4*sign + 18 + abs*6 + min*6 + max*6 

      s4 = sign(1.0, d4)
      s3 = sign(1.0, d3)
      s2 = sign(1.0, d2)
      s1 = sign(1.0, d1)

      g6 = s4*max(0.0, min(abs(d4), s4*b*d3))
      g5 = s3*max(0.0, min(abs(d3), s3*b*d4))
      g4 = s3*max(0.0, min(abs(d3), s3*b*d2))
      g3 = s2*max(0.0, min(abs(d2), s2*b*d3))
      g2 = s2*max(0.0, min(abs(d2), s2*b*d1))
      g1 = s1*max(0.0, min(abs(d1), s1*b*d2))
