-----
Compilation of FFV-C

Kenji Ono 	keno@riken.jp
AICS, RIKEN
September 2012
-----


インストールの詳細は，docディレクトリにあるffvc_ug.pdfをご覧ください．

1. TextParserのインストール
textParser-x.x.tar.gzを展開し，トップディレクトリのconfig_tp.shを実行．

$ configure_tp.sh /usr/local/textparser
$ make
$ sudo make install または make install
$ make distclean


2.　CPMlibのインストール

CPMlib-x.x.x.tar.gzを展開し，トップディレクトリのconfig_cpm.shを実行．

$ configure_cpm.sh /usr/local/cpm
$ make
$ sudo make install または make install
$ make distclean

configureがうまくいかない場合は，Makefile_handを使い，コンパイルする．



3.　FFV-Cのコンパイル
以下の順序でコンパイルを行う．コンパイルをやり直す場合には，make clean, make allclean を実行．

3-1. Polylib
Makefileのマクロ変数を指定し，コンパイルする．

$ make depend
$ make

3-2. Cutlib
Makefileのマクロ変数を指定し，コンパイルする．

$ make depend
$ make

3-3. PMlib
Makefileのマクロ変数を指定し，コンパイルする．

$ make depend
$ make

3-4. FFV-C
make_settingのマクロ変数を指定し，コンパイルする．

$ make depend
$ make




Direction of compile

Please find a ffvc_ug.pdf file in document directory and see how to install in detail.



	


