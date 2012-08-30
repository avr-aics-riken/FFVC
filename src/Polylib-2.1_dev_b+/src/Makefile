############################################################################
#
# sigle版アーカイブ、mpi版アーカイブ作成用Makefile
#													Create	2010.08.20
#													Update	2010.10.25
# sigle版アーカイブのみ作成		: make single
# mpi版アーカイブのみ作成		: make mpi
# single版、mpi版をまとめて作成	: make
# includeファイルの依存関係更新 : make depend
#
############################################################################

# FX10の場合，MPI_DIRはコメントアウト
MPI_INC_DIR	= /opt/openmpi/include
TP_INC_DIR  = /usr/local/textparser/include

#CXX			= g++
CXX			= icpc
#CXX			= mpiFCCpx

#AR			= ar cr
AR			= ar cru

RANLIB		= ranlib

TARGET_DIR	= ../lib

###########################################################################
# 以下のワーニングを抑制するために-fno-strict-aliasingを追加 2010.10.14
# warning: type-punning to incomplete type might break strict-aliasing rules
#
# デバッグ時にはCFLAGSに以下を追加のこと
# -g -DDEBUG
#
# 三角形IDファイルをテキストモードで出力するにはCFLAGSに以下を追加のこと
# -DSAVE_ID_ASCII
###########################################################################

CFLAGS		= -O3 -Wall \
            -I../include \
            -I$(MPI_INC_DIR) \
            -fno-strict-aliasing \
            -I$(TP_INC_DIR) \
#            -I/usr/include \
#            -DSAVE_ID_ASCII -DDEBUG

# FX10
#CFLAGS		= -Kfast, simd=2, ocl, preex \
#		-I../include\
#		-I/opt/FJSVfxlang/GM-1.2.0-03/include \
#		-I/opt/FJSVfxlang/GM-1.2.0-03/include/mpi/fujitsu \
#		-I/usr/include \
#		-I$(TP_INC_DIR) \
#		-fno-strict-aliasing

# K
#CFLAGS	 = -Kfast, simd=2, ocl, preex \
#   -I../include\
#   -I/opt/FJSVtclang/GM-1.2.0-07/include \
#   -I/opt/FJSVtclang/GM-1.2.0-07/include/mpi/fujitsu \
#   -I/usr/include \
#   -I$(TP_INC_DIR) \
#   -fno-strict-aliasing

#####################################
# single版マクロ
#####################################
SNG_TARGET	= $(TARGET_DIR)/libPolylib.a
GFX_OBJS	= Polylib.o \
			  polygons/VTree.o \
			  polygons/Polygons.o \
			  polygons/TriMesh.o \
			  groups/PolygonGroup.o \
			  file_io/TriMeshIO.o \
			  file_io/stl.o \
			  file_io/triangle_id.o \
			  c_lang/CPolylib.o 

#####################################
# MPI版マクロ
#####################################
MPI_TARGET	= $(TARGET_DIR)/libMPIPolylib.a
MPI_OBJS	= MPIPolylib.o \
			  c_lang/CMPIPolylib.o

#####################################
# ユーティリティ
#####################################
UTIL_OBJS	= util/time.o \

#####################################
# ヘッダファイル名作成用マクロ
#####################################
TMP = $(addsuffix .h, $(basename $(GFX_OBJS) $(MPI_OBJS) $(UTIL_OBJS)))
HDRS = $(addprefix ../include/, $(TMP))

#####################################
# single版、MPI版作成ルール
#####################################
all : depend single mpi

#####################################
# single版作成ルール
#####################################
single			: $(SNG_TARGET)
$(SNG_TARGET)	: $(GFX_OBJS) $(UTIL_OBJS)
	$(AR) $@ $^
	$(RANLIB) $(SNG_TARGET)

#####################################
# MPI版作成ルール
#####################################
mpi				: $(MPI_TARGET)
$(MPI_TARGET)	: $(GFX_OBJS) $(MPI_OBJS) $(UTIL_OBJS)
	$(AR) $@ $^
	$(RANLIB) $(MPI_TARGET)

#####################################
# サフィックスルール定義
#####################################
%.o: %.cxx
	$(CXX) $(CFLAGS) -o $@ -c $<

#####################################
# ヘッダファイル依存関係
#####################################
depend		 : $(HDRS)
	-@ $(RM) depend.inc
	-@ for i in $^; do\
		cpp $(CFLAGS) -MM $$i >> depend.inc;\
	done
-include depend.inc

clean:
	rm -f $(SNG_TARGET) $(MPI_TARGET) $(GFX_OBJS) $(MPI_OBJS) $(UTIL_OBJS)
