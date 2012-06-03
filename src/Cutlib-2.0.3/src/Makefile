#####################################################################################
#
#  Cutlib インストール方法
#  -----------------------
#
#  利用するPolylibのバージョンにあわせてPOLYLIB_DIRを指定してください．
#
#  「環境に依存するマクロ変数」を確認・修正後
#      make
#  で，libディレクトリにlibcutlib.aがインストールされます．
#
#
#  環境に依存するマクロ変数
#  ------------------------
#
#    - CXX: C++コンパイラ
#    - AR:  アーカイバコマンド(Mac OS X では"u"を追加)
#    - POLYLIB_DIR: Polylibインストールディレクトリ
#    - CFLAGS: 最適化オプションの変更
#
#
#  Octree版関数の作成方法
#  ----------------------
#
#    1. 「with Octree lib」下の行のコメントを外す．
#    2. 「without Octree lib」下をコメントアウト.
#    3. OCTREE_DIRにOctreeライブラリのインストールディレクトリを指定.
# 
#    (*)Octreeライブラリは, MRpFTT_SGI_20080428 による動作確認をしています．
#
#
#  計算時間測定機能
#  ----------------
#
#   Cutlib内部の計算時間測定機能を有効にするには，OBJSにCutTiming.oを，
#   DEFINESに-DCUTLIB_DEBUGを追加します．
#
#####################################################################################

CUTLIB = ../lib/libcutlib.a

all: $(CUTLIB)

OBJS = Cutlib.o CutlibCore.o CutUtil.o CutTriangle.o CutBoundary.o
#OBJS = Cutlib.o CutlibCore.o CutUtil.o CutTriangle.o CutBoundary.o CutTiming.o
 
#CXX = g++
#CXX = mpiFCCpx
CXX = icpc

#AR = ar cr
AR = ar cru
RANLIB = ranlib

POLYLIB_DIR = ../../Polylib_2_0_3

POLYLIB_FLAGS = -I$(POLYLIB_DIR)/include -I/usr/include/libxml2

# FX10
#POLYLIB_FLAGS = -I$(POLYLIB_DIR)/include \
#        -I/home/c15046/LIBXML2/libxml2-2.7.8_lib/include/libxml2 \
#        -fno-strict-aliasing

# with Octree lib ############################################################
#OCTREE_DIR = ../../sgi-octree
#OCTREE_FLAGS = -DCUTLIB_OCTREE -DTD_USE_NAMESPACE -I$(OCTREE_DIR)/include
##############################################################################

# without Octree lib #########################################################
OCTREE_FLAGS = 
##############################################################################

#DEFINES = -DCUTLIB_TIMING
#DEFINES = -DCUTLIB_DEBUG

# ポリゴンの表裏が統一されていることを仮定せずに計算
DEFINES = -DIGNORE_NORMAL_DIRECTION

# STLファイル中のポリゴン法線データを信用せず，毎回法線計算
DEFINES += -DIGNORE_STL_NORMAL
#DEFINES = 

#CFLAGS = -g -Wall -I../include $(POLYLIB_FLAGS) $(OCTREE_FLAGS) $(DEFINES)
CFLAGS = -O3 -Wall -I../include $(POLYLIB_FLAGS) $(OCTREE_FLAGS) $(DEFINES)

# FX10
#CFLAGS = -Kfast,ocl,preex,simd=2 \
#        -I../include  $(POLYLIB_FLAGS) $(OCTREE_FLAGS) $(DEFINES)

$(CUTLIB): $(OBJS)
	$(AR) $@ $(OBJS)
	$(RANLIB) $@

.cpp.o:
	$(CXX) $(CFLAGS) -c $<

clean:
	rm -rf *.o $(CUTLIB)

depend: $(OBJS:.o=.cpp)
	@ rm -rf depend.inc
	@ for i in $^; do\
		$(CXX) $(CFLAGS) -MM $$i >> depend.inc;\
	done

-include depend.inc
