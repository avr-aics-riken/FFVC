# FFV-C   Frontflow / violet Cartesian

## OUTLINE


Frontflow/violet Cartesian is a three-dimensional unsteady incompressible thermal flow simulator on a Cartesian grid system. This solver is designed so as to assist practical design in industry fields by powerful voxel based approach and to achieve higher parallel performance on massively parallel computers.


## COPYRIGHT
- Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.  All rights reserved.
- Copyright (c) 2011-2015 Institute of Industrial Science, The University of Tokyo.  All rights reserved.
- Copyright (c) 2012-2016 Advanced Institute for Computational Science, RIKEN.  All rights reserved.
- Copyright (c) 2016-2020 Research Institute for Information Technology(RIIT), Kyushu University.  All rights reserved.


## SOFTWARE REQUIREMENT

- Cmake
- MPI library
- TextParser - Text Parser library
- PMlib - Performance Monitor library
- CPMlib - Computational space Partitioning Management library
- CDMlib - Cartesian Data Management library
- Polylib - Polygon management library


## INGREDIENTS

ChangeLog.md    History of development
CMakeLists.txt  Makefile for cmake
License.txt     License to apply
Readme.md       This document
doc/            Documents
doxygen/        Doxygen files
example/        Example files
src/		        Source files
   ASD/         Acttive Subdomain class
   F_CORE/      f90 subroutines in ffvc
   F_LS/        f90 subroutines of linear solvers
   F_VOF/       VOF routines
   FB/          flow base class
   FFV/         ffvc class
   FILE_IO/     file io utilities
   Geometry/    Geometry manipulation
   IP/          Intrinsic problems
src_util/       Combsph


## How to build

### Build

~~~
$ export FFVC_HOME=/hogehoge
$ mkdir build
$ cd build
$ cmake [options] ..
$ make
$ sudo make install
~~~


### Options

`-D INSTALL_DIR=` *Install_directory*

>  Specify the directory that this library will be installed. Built library is
   installed at `install_directory/lib` and the header files are placed at
   `install_directory/include`.

`-D enable_OPENMP=` {yes | no}

>  This option makes OpenMP directives effect. Default is yes.

`-D with_MPI=` {yes | no}

>  If you use an MPI library, specify `with_MPI=yes` (default).

`-D real_type=` {float | double}

>  Specify the type of floating point. If this option is omitted, the default is float.

`-D with_TP=` *Installed_Directory*

> Specify the directory path that TextParser is installed.

`-D with_PM=` *Installed_Directory*

> Specify the directory path that PMlib is installed.

`-D with_PL=` *Installed_Directory*

> Specify the directory path that Polylib is installed.

`-D with_CPM=` *Installed_Directory*

> Specify the directory path that CPMlib is installed.

`-D with_CDM=` *Installed_Directory*

> Specify the directory path that CDMlib is installed.

`-D with_PAPI=` *Installed_Directory*

> Specify the directory path that PAPI is installed.

## Configure Examples

`$ export FFVC_HOME=hogehoge`

In following examples, assuming that TextParser, PMlib, Polylib, CPMlib, and CDMlib are installed under the FFVC_HOME directory. If not, please specify applicable directory paths.


### INTEL/GNU compiler

~~~
$ cmake -DINSTALL_DIR=${FFVC_HOME}/FFVC \
        -Dreal_type=float \
        -Denable_OPENMP=yes \
        -Dwith_MPI=yes \
        -Dwith_TP=${FFVC_HOME}/TextParser \
        -Dwith_PM=${FFVC_HOME}/PMlib \
        -Dwith_PAPI=OFF \
        -Dwith_PL=${FFVC_HOME}/Polylib \
        -Dwith_CPM=${FFVC_HOME}/CPMlib \
        -Dwith_CDM=${FFVC_HOME}/CDMlib ..
~~~

#### NOTE
In case of some Intel compiler environment, please specify environemnt variables
`export CC=icc CXX=icpc F90=ifort FC=ifort` before compiling.


### FUJITSU compiler / FX10, FX100, K computer on login nodes (Cross compilation)

~~~
$ cmake -DINSTALL_DIR=${FFVC_HOME}/FFVC \
        -DCMAKE_TOOLCHAIN_FILE=../cmake/Toolchain_fx10.cmake \
        -Dreal_type=float \
        -Denable_OPENMP=yes \
        -Dwith_MPI=yes \
        -Dwith_TP=${FFVC_HOME}/TextParser \
        -Dwith_PM=${FFVC_HOME}/PMlib \
        -Dwith_PAPI=OFF \
        -Dwith_PL=${FFVC_HOME}/Polylib \
        -Dwith_CPM=${FFVC_HOME}/CPMlib \
        -Dwith_CDM=${FFVC_HOME}/CDMlib ..

$ cmake -DINSTALL_DIR=${FFVC_HOME}/FFVC \
        -DCMAKE_TOOLCHAIN_FILE=../cmake/Toolchain_fx100.cmake \
        -Dreal_type=float \
        -Denable_OPENMP=yes \
        -Dwith_MPI=yes \
        -Dwith_TP=${FFVC_HOME}/TextParser \
        -Dwith_PM=${FFVC_HOME}/PMlib \
        -Dwith_PAPI=OFF \
        -Dwith_PL=${FFVC_HOME}/Polylib \
        -Dwith_CPM=${FFVC_HOME}/CPMlib \
        -Dwith_CDM=${FFVC_HOME}/CDMlib ..

$ cmake -DINSTALL_DIR=${FFVC_HOME}/FFVC \
        -DCMAKE_TOOLCHAIN_FILE=../cmake/Toolchain_K.cmake \
        -Dreal_type=float \
        -Denable_OPENMP=yes \
        -Dwith_MPI=yes \
        -Dwith_TP=${FFVC_HOME}/TextParser \
        -Dwith_PM=${FFVC_HOME}/PMlib \
        -Dwith_PAPI=OFF \
        -Dwith_PL=${FFVC_HOME}/Polylib \
        -Dwith_CPM=${FFVC_HOME}/CPMlib \
        -Dwith_CDM=${FFVC_HOME}/CDMlib ..
~~~

##### Note
- On Fujitsu machines(fx10, K, fx100), confirm appropriate directrory path for compiler environment.
- Before building, execute following command to clean for sure. `$ make distclean`


## Contributors

- Kenji       Ono           keno@{riken, cc.kyushu-u.ac}.jp
- Masako      Koyama(Iwata)
- Tsuyoshi    Tamaki
- Yasuhiro    Kawashima
- Kei         Akasaka
- Soichiro    Suzuki
- Junya       Onishi
- Ken         Uzawa
- Kazuhiro    Hamaguchi


## Note of FLOP COUNT
The number of floating operations in a loop is estimated by the value measured on the K computer.

|function|float|double|
|:--|--:|--:|
|sqrt|10|20|
|exp|22|26|
|log|21|25|
|log10|24|30|
|sin|29|31|
|cos|29|29|
|min|1|1|
|max|1|1|
|abs|1|1|
|add|1|1|
|subtract|1|1|
|multiply|1|1|
|division|8|13|
|aint|6|6|
|anint|10|10|
|ceiling|2|2|
|floor|2|2|
|nint|6|6|

### Comment of precision of Fortran
For example, the `-CcdRR8` option for fortran preprocessor convert variables, functions, and constants higher precision version in source code. Thus, functions in source code is described by floating version.
