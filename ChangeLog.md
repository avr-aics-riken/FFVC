# FFV-C

## TODO


## Compiled

|Environment|Serial|MPI |Host|
|:--|:--:|:--:|:--|
|Intel / Intel 17.0.1 |x|ok(1)|Apple MacBook Pro 10.12.4|
|Intel / GNU 6.2.0    |x|ok(1)|Apple MacBook Pro 10.12.4|
|Fujitsu / fx10       |x|ok|Hayaka /opt/FJSVfxlang/1.2.1|


- (1) OpenMPI 1.10.4


## REVISION HISTORY

---
- 2018-06-15 Version 2.5.3
  - intel_F_TCS
  - ISNAN => isnan()

---
- 2017-07-06 Version 2.5.2
  - set(CMAKE_FIND_ROOT_PATH /opt/FJSVtclang/1.2.0) in Toolchain_K.cmake


---
- 2017-04-22 Version 2.5.1
  - fix cmake install
  - fx10 compilation


---
- 2017-04-21 Version 2.5.0
  - Cmake version
  - serial version has a problem of ambiguous MPI_SUM to compile


---
- 2016-02-04 Version 2.4.4
  - r1639 update year 2016


---
- 2015-11-28 Version 2.4.3
  - r1638 bug fix in ab2 routine


---
- 2015-11-28 Version 2.4.2
  - r1637 --with-comp becomes essential
    - change configure.ac, INSTALL, and NEWS


---
- 2015-09-13 Version 2.4.1
  - r1636 suppress check print


---
- 2015-09-12 Version 2.4.0
  - r1635 bugfix: check duplication in loadBCs()

  - r1634 version 2.4.0
    - version 2.4.x => PMlib that corresponds to process group


---
- 2015-07-27 Version 2.3.15
  - r1633 introduce process group
  - r1632 update PMlib-4.0.1


---
- 2015-07-26 Version 2.3.14
  - r1631 merge the rest of pull request #23


---
- 2015-07-25 Version 2.3.13
  - r1630 merge pull request #23 manually and bug fix
    - update user guide : Output statistics


---
- 2015-07-25 Version 2.3.12
  - r1629 modify wall parameter
    - remove getWallType()
    - add /TurbulenceModeling/VelocityProfile
    - update included user guide


---
- 2015-07-25 Version 2.3.11
  - r1628 use MPIPolylib::load_only_in_rank0() and  distribute_only_from_rank0()
    - get accurate BC area


---
- 2015-07-15 Version 2.3.10
  - r1627 update examples for latest parameters
  - included examples; 2Dcavity, 3Dcavity, 2Dcyl, Jet, LDC, PMT


---
- 2015-07-10 Version 2.3.8
  - r1626 change parameter orgnization
    - /Output
    - LimitedCompresibility
    - Preconditioner


---
- 2015-07-10 Version 2.3.7
  - r1625 change mask from d_cdf to b_bcd
  - update_vec(), update_vec4_(), divergence_cc_(), pvec_central_(), pvec_central_les_()



---
- 2015-07-09 Version 2.3.6
  - modify bicg
  - add Limited compressibility
  - non-dim param for stability
  - modify Stability Control
  - add Stability Control function
  - gather isnan flag in Loop()


---
- 2015-07-05 Version 2.3.5
  - modify Vspec implementation
  - encVbitIBC(), encVbitIBCrev(), encVbitOBC()
  - expire VBC_UWD flag


---
- 2015-06-30 Version 2.3.4
  - suppress error message in mergeSameSolidCell()
  - expire paintCutOnPoint()


---
- 2015-06-29 Version 2.3.3
  - add mergeSameSolidCell()


---
- 2015-06-28 Version 2.3.2
  - add core_psor.h, core_sor2sma.h


---
- 2015-06-28 Version 2.3.1
  - add *cmp and *mat in Geometry class
  - brush up algorithms in Geometry


---
- 2015-06-21 Version 2.2.3
- introduce ffvc-config
- introduce dryrunBC()
- modify extractVelLBC() and related parts
- expire mask for residuals in psor_() and psor2sma_core_()


---
- 2015-06-19 Version 2.2.2
  - remove mask from blas_calc_rk_() and blas_calc_ax_() to make stable
  - rm VLD_CNVG >> active_bit
  - modify encPbit(), encPbitN(), norm_v_div_l2_(), norm_v_div_max_()
  - replacing if-branch to mask combination in psor_(), psor2sma_core_()
  - arrange Geometry.C
  - move onBit()/offBit() from VoxInfo.h to FB_Define.h
  - add special BC polynomial6 for chiba-u
  - modify quantizeCut()
  - add range for Glyph


---
- 2015-06-14 Version 2.2.1
  - error: MainLoop() C.Interval[Control::tg_End].printInfo("tg_End") destroys array


---
- 2015-06-10 Version 2.2.0
  - clean package


---
- 2015-06-09 Version 2.1.9
  - This fill algorithm can successfully generate a voxel model for Nasal2015 case.


---
- 2015-06-06 Version 2.1.8
  - Change to run configure
  - Change configure.ac
  - ducky, cavity, BMW_M3, nasal(0.2, 0.4mm)
  - enable fileout of bvx format
  - remove obsolete svx output


---
- 2015-06-01 Version 2.1.7
  - remove fill_bid.h


---
- 2015-03-01 Version 2.1.5
  - define ISNAN in FBdefine.h
  - update libraries
  - set dummy values for cmp[m].medium in case of Periodic
  - skip rading filepath for intrinsic problems
  - confirm detail later
  - reference point must be dimensional, partially it was non-dimensional


---
- 2015-02-25 Version 2.1.4
  - generate polylib.tp file from boundary section
  - modify fill medium


---
- 2015-02-17 Version 2.1.3
  - 2Dcavity, 3Dcavity, 2Dcyl, Cylinder, Jet, LDC, PMT, ThermalCavity, SHC1D
  - check Cylinder
  - check 2Dcavity, 3Dcavity, 2Dcyl, Cylinder, LDC, PMT, SHC1D ,Jet, HeatedPlate, ThermalCavity


---
- 2015-02-11 Version 2.1.2
  - introduction of fill hint


---
- 2015-02-11 Version 2.1.1
  - introduction of fill hint

---
- 2015-02-07 Version 2.1.0
  - insert mpi.h at first in ParseMat.h and CompoFraction.h
  - change calculation policy from ceil to round off to get g_voxel in SD_getParameter()


-----
- 2015-02-05 Version 2.0.9
  - add printForceAvr()


---
- 2015-02-01 Version 2.0.8

---
- 2015-02-01 Version 2.0.7

---
- 2015-01-27 Version 2.0.6

---
- 2015-01-17 Version 2.0.5

---
- 2015-01-11 Version 2.0.4


---
- 2014-12-25 Version 2.0.3

---
- 2014-12-05 Version 2.0.2

---
- 2014-12-04 Version 2.0.1

---
- 2014-11-29 Version 2.0.0

---
- 2014-11-23 Version 1.9.9

---
- 2014-11-16 Version 1.9.8

---
- 2014-11-06 Version 1.9.6

---
- 2014-11-05 Version 1.9.5

---
- 2014-10-19 Version 1.8.9

---
- 2014-10-18 Version 1.8.8

---
- 2014-10-18 Version 1.8.7

---
- 2014-09-22 Version 1.8.6

---
- 2014-09-16 Version 1.8.5

---
- 2014-09-06 Version 1.8.4

---
- 2014-08-30 Version 1.8.3

---
- 2014-08-26 Version 1.8.2

---
- 2014-08-01 Version 1.8.1

---
- 2014-04-24 Version 1.8.0

---
- 2014-04-16 Version 1.7.9

---
- 2014-04-03 Version 1.7.8

---
- 2014-04-02 Version 1.7.7

---
- 2014-03-23 Version 1.7.6

---
- 2014-03-22 Version 1.7.5

---
- 2014-03-19 Version 1.7.4

---
- 2014-03-18 Version 1.7.3

---
- 2014-03-17 Version 1.7.2

---
- 2014-03-16 Version 1.7.1

---
- 2014-03-16 Version 1.7.0

---
- 2014-03-11 Version 1.6.9

---
- 2014-03-10 Version 1.6.8

---
- 2014-03-09 Version 1.6.7

---
- 2014-03-09 Version 1.6.6

---
- 2014-03-08 Version 1.6.5

---
- 2014-03-04 Version 1.6.4

---
- 2014-03-04 Version 1.6.3

---
- 2014-02-26 Version 1.6.2

---
- 2014-02-19 Version 1.6.1

---
- 2014-02-11 Version 1.6.0

---
- 2014-02-09 Version 1.5.9

---
- 2014-02-08 Version 1.5.8

---
- 2014-01-03 Version 1.5.4

---
- 2013-12-29 Version 1.5.3

---
- 2013-12-23 Version 1.5.2

---
- 2013-12-05 Version 1.5.1

---
- 2013-11-13 Version 1.5.0

---
- 2013-11-03 Version 1.4.9

---
- 2013-10-30 Version 1.4.8

---
- 2013-10-27 Version 1.4.7

---
- 2013-10-12 Version 1.4.6

---
- 2013-10-11 Version 1.4.5

---
- 2013-10-06 Version 1.4.4

---
- 2013-10-05 Version 1.4.3

---
- 2013-10-03 Version 1.4.2

---
- 2013-10-01 Version 1.4.1

---
- 2013-09-29 Version 1.4.0

---
- 2013-09-22 Version 1.3.9

---
- 2013-09-13 Version 1.3.8

---
- 2013-09-08 Version 1.3.7

---
- 2013-08-31 Version 1.3.5

---
- 2013-08-22 Version 1.3.4

---
- 2013-08-19 Version 1.3.3

---
- 2013-08-14 Version 1.3.2

---
- 2013-08-12 Version 1.3.1

---
- 2013-08-05 Version 1.3.0

---
- 2013-08-02 Version 1.2.9

---
- 2013-08-01 Version 1.2.8

---
- 2013-07-28 Version 1.2.7

---
- 2013-07-27 Version 1.2.6

---
- 2013-07-16 Version 1.2.5

---
- 2013-07-13 Version 1.2.4

---
- 2013-07-11 Version 1.2.3

---
- 2013-07-11 Version 1.2.2

---
- 2013-07-10 Version 1.2.1

---
- 2013-07-06 Version 1.2.0

---
- 2013-07-05 Version 1.1.8

---
- 2013-06-15 Version 1.1.7

---
- 2013-05-08 Version 1.1.4

---
- 2013-04-26 Version 1.1.3

---
- 2013-03-31 Version 1.0.8

---
- 2013-02-17 Version 1.0.3

---
- 2013-01-26 Version 1.0.2

---
- 2012-12-01 Version 0.9.7

---
- 2012-11-16 Version 0.9.6

---
- 2012-10-18 Version 0.8.9

---
- 2012-09-30 Version 0.8.4

---
- 2012-08-25 Version 0.7.2

---
- 2012-07-30 Version 0.6.0

---
- 2012-06-30 Version 0.2.0
