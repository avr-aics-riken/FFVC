#####################################################################################
#
#  PMlib コンパイル方法
#  -----------------------
#
#
#  「環境に依存するマクロ変数」を確認・修正後
#      make
#  で，libディレクトリにlibpmlib.aがインストールされます．
#
#
#  環境に依存するマクロ変数
#  ------------------------
#
#    - CXX: C++コンパイラ
#    - AR:  アーカイバコマンド(Mac OS X では"u"を追加)
#    - CFLAGS: 最適化オプションの変更
#
#
#####################################################################################

# FX10の場合，MPI_DIRはコメントアウト
MPI_DIR		= /opt/openmpi

PMLIB = ../lib/libpmlib.a

all: $(PMLIB)

OBJS = PerfWatch.o PerfMonitor.o
 
#CXX = g++
#CXX = mpiFCCpx
CXX = icpc

#AR = ar cr
AR = ar cru
RANLIB = ranlib


CFLAGS = -O3 -Wall -I../include -I$(MPI_DIR)/include
#CFLAGS = -O3 -I../include -I$(MPI_HOME)/include

$(PMLIB): $(OBJS)
	$(AR) $@ $(OBJS)
	$(RANLIB) $@

.cpp.o:
	$(CXX) $(CFLAGS) -c $<

clean:
	rm -rf *.o $(PMLIB)

depend: $(OBJS:.o=.cpp)
	@ rm -rf depend.inc
	@ for i in $^; do\
		$(CXX) $(CFLAGS) -MM $$i >> depend.inc;\
	done

-include depend.inc
