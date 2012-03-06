!   *********************************************************
!
!   SPHERE - Skeleton for PHysical and Engineering REsearch
!  
!   Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
!
!   *********************************************************
!
!> @file omp_head.h
!! @brief OpenMPのディレクティブ
!! @author keno, FSI Team, VCAD, RIKEN
!<

!$OMP PARALLEL
!$OMP DO SCHEDULE(dynamic,1) &

!>>$OMP DO SCHEDULE(static,2) &
!>>$OMP DO &

