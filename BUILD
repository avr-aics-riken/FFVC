-----
Compilation order of FFV-C

Kenji Ono 	keno@riken.jp
AICS, RIKEN
July 2012
-----

インストールの詳細は，docディレクトリにあるffvc_ug.pdfをご覧ください．

1.　V-Sphere　インストール
下記のシェルを用いて，/usr/local/sphereにインストールします．

$ configure_sph_ompi_float.sh /usr/local/sphere
$ make
$ sudo make install または make install
$ make distclean

----------------------
configure_sph_ompi_float.sh
----------------------
#! /bin/sh
./configure --prefix=$1 \
            --with-comp=INTEL \
            --with-ompi=/usr/local/ompi \
            CC=icc \
            CFLAGS=-O3\
            CXX=icpc \
            CXXFLAGS=-O3\
            FC=ifort \
            FCFLAGS=-O3\
            F90=ifort \
            F90FLAGS=-O3\
            LDFLAGS=-L/opt/intel/composerxe/lib
----------------------

1.　sphPrjTool を用いたCBCのコンパイル
既に,並列ライブラリ(mpich または OpenMPI)と V-Sphere が正しくインストールされていることを前提とします．
まず，提供される project_local_settings ファイルをコピーしほかの名前にしておきます（ex. project_local_settings_old）．

つぎに環境変数SPHEREDIRを設定します．下記で，INSTALL_DIRにはV-Sphere のインストールディレクトリを指定します．
その後，src/PRJ_CBC ディレクトリで sphPrjTool を起動し,プロジェクト環境をリセットします.
localsettings オプション を指定して reset コマンドを実行すると，sphere ライブラリの config/sph-cfg.xml に記録されているコンパイル 環境情報を元にして，プロジェクト環境情報 PRJ CBC.xml を再設定します.

$ export SPHEREDIR=INSTALL_DIR
$ cd src/PRJ_CBC
$ sphPrjTool PRJ_CBC.xml
sphPrjTool> reset localsettings 
sphPrjTool> print 
sphPrjTool> save 
sphPrjTool> quit

2.　$ make
    実行モジュールが各ソルバークラスのbinの下に生成されます．

3.　各プラットホーム向けのヒントはユーザーガイドを参照してください．



	


