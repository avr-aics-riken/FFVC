#####################################################################################
#
#  Cutlib 
#  -----------------------
#
#  At first, edit a make_setting.* file for an each environment, then make
#  If Cutlib is placed inside the directory of FFVC solver, this makefile is ignored.
#
#####################################################################################


# Default environment
MACHINE=intel


### compiler detection

ifeq ($(MACHINE),gnu)
include make_setting.gnu
endif

ifeq ($(MACHINE),fx)
include make_setting.fx
endif

ifeq ($(MACHINE),ibm)
include make_setting.ibm
endif

ifeq ($(MACHINE),intel)
include make_setting.intel
endif


all:
	( \
	cd src; \
	make \
		CXX='$(CXX)' \
		CXXFLAGS='$(CXXFLAGS)' \
		AR='$(AR)' \
		RANLIB='$(RANLIB)' \
		RM='$(RM)' \
		TP_DIR='$(TP_DIR)' \
		POLYLIB_DIR='$(POLYLIB_DIR)' \
	)

clean:
	(cd src; make clean)

allclean: clean
	(cd lib; $(RM) $(TARGET) )

depend:
	(cd src; make depend)
