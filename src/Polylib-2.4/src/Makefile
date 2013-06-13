############################################################################
#
# sigle版アーカイブ、mpi版アーカイブ作成用Makefile
#													Create	2010.08.20
#													Update	2010.10.25
# sigle版アーカイブのみ作成		  : make single
# mpi版アーカイブのみ作成	      : make mpi
# single版、mpi版をまとめて作成	: make
# includeファイルの依存関係更新 : make depend
#
# 以下のワーニングを抑制するために-fno-strict-aliasingを追加 2010.10.14
# warning: type-punning to incomplete type might break strict-aliasing rules
#
# デバッグ時にはCXXFLAGSに以下を追加のこと
# -g -DDEBUG
#
# 三角形IDファイルをテキストモードで出力するにはCXXFLAGSに以下を追加のこと
# -DSAVE_ID_ASCII
###########################################################################

TARGET_DIR	= ../lib


#CXXFLAGS += -DSAVE_ID_ASCII -DDEBUG
#CXXFLAGS += -g -DDEBUG



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
single		: $(SNG_TARGET)
$(SNG_TARGET)	: $(GFX_OBJS) $(UTIL_OBJS)
	$(AR) $@ $^
	$(RANLIB) $(SNG_TARGET)


#####################################
# MPI版作成ルール
#####################################
mpi		: $(MPI_TARGET)
$(MPI_TARGET)	: $(GFX_OBJS) $(MPI_OBJS) $(UTIL_OBJS)
	$(AR) $@ $^
	$(RANLIB) $(MPI_TARGET)


#####################################
# サフィックスルール定義
#####################################
%.o: %.cxx
	$(CXX) $(CXXFLAGS) -I../include $(UDEF_INC_PATH) -o $@ -c $<


#####################################
# ヘッダファイル依存関係
#####################################
depend		 : $(HDRS)
	-@ $(RM) depend.inc
	-@ for i in $^; do\
		cpp -I../include -MM $$i >> depend.inc;\
	done
-include depend.inc

clean:
	$(RM) $(SNG_TARGET) $(MPI_TARGET) $(GFX_OBJS) $(MPI_OBJS) $(UTIL_OBJS)
