#####################################################################################
#
#  Cutlib （単体ビルド用）
#  -----------------------
#
#  環境に合わせて、make_settingマクロを編集後、make
#
#
#####################################################################################

include make_setting

all:
	( \
	cd src; \
	make \
		CXX='$(CXX)' \
		AR='$(AR)' \
		RANLIB='$(RANLIB)' \
		RM='$(RM)' \
		TP_DIR='$(TP_DIR)' \
	)

clean:
	(cd src; make clean)

allclean: clean
	(cd lib; $(RM) $(TARGET) )

depend:
	(cd src; make depend)
