# #################################################################
#
# output layoutfiles
#
# Copyright (c) All right reserved. 2012-2013
#
# Institute of Industrial Science, University of Tokyo, Japan.
#
# #################################################################

include ../FFV/make_setting

all:
	( \
	cd $(PMLIB_DIR)/src; \
	make \
		CXX='$(CXX)' \
		AR='$(AR)' \
		RANLIB='$(RANLIB)' \
		RM='$(RM)' \
		MPI_DIR='$(MPI_DIR)' \
	)

	( \
	cd ../FB; \
	make \
		CXX='$(CXX)' \
		CXXFLAGS='$(CXXFLAGS)' \
		CC='$(CC)' \
		CFLAGS='$(CFLAGS)' \
		FC='$(FC)' \
		FCFLAGS='$(FCFLAGS)' \
		F90='$(F90)' \
		F90FLAGS='$(F90FLAGS)' \
		UDEF_INC_PATH='$(UDEF_INC_PATH)' \
		AR='$(AR)' \
		RANLIB='$(RANLIB)' \
		RM='$(RM)' \
	)

	( \
	make -f MakeLAYOUT \
		CXX='$(CXX)' \
		CXXFLAGS='$(CXXFLAGS)' \
		F90='$(F90)' \
		F90FLAGS='$(F90FLAGS)' \
		LDFLAGS='$(LDFLAGS)' \
		LIBS='$(LIBS)' \
		UDEF_INC_PATH='$(UDEF_INC_PATH)' \
		UDEF_LIB_PATH='$(UDEF_LIB_PATH)' \
		UDEF_LIBS='$(UDEF_LIBS_UTIL)' \
		RM='$(RM)' \
	)

clean: 
	(make -f MakeLAYOUT clean)
  
allclean: clean
	(cd ../FB; make clean)
	(cd $(PMLIB_DIR)/src; make clean)

depend:
	(cd $(PMLIB_DIR)/src; make depend)
	(cd ../FB; make depend)
	(make -f MakeLAYOUT depend)



# DO NOT DELETE

main.o: layout.h /usr/include/stdio.h /usr/include/sys/cdefs.h
main.o: /usr/include/sys/_symbol_aliasing.h
main.o: /usr/include/sys/_posix_availability.h /usr/include/Availability.h
main.o: /usr/include/AvailabilityInternal.h /usr/include/_types.h
main.o: /usr/include/sys/_types.h /usr/include/machine/_types.h
main.o: /usr/include/i386/_types.h /usr/include/secure/_stdio.h
main.o: /usr/include/secure/_common.h /usr/include/string.h
main.o: /usr/include/strings.h /usr/include/secure/_string.h
main.o: /usr/include/errno.h /usr/include/sys/errno.h
main.o: /usr/include/sys/types.h /usr/include/sys/appleapiopts.h
main.o: /usr/include/machine/types.h /usr/include/i386/types.h
main.o: /usr/include/machine/endian.h /usr/include/i386/endian.h
main.o: /usr/include/sys/_endian.h /usr/include/libkern/_OSByteOrder.h
main.o: /usr/include/libkern/i386/_OSByteOrder.h /usr/include/sys/_structs.h
main.o: /usr/include/sys/stat.h /usr/include/sys/param.h
main.o: /usr/include/sys/syslimits.h /usr/include/machine/param.h
main.o: /usr/include/i386/param.h /usr/include/i386/_param.h
main.o: /usr/include/limits.h /usr/include/machine/limits.h
main.o: /usr/include/i386/limits.h /usr/include/i386/_limits.h
main.o: /usr/include/sys/signal.h /usr/include/machine/signal.h
main.o: /usr/include/i386/signal.h /usr/include/i386/_structs.h
main.o: /usr/include/stdlib.h /usr/include/sys/wait.h
main.o: /usr/include/sys/resource.h /usr/include/alloca.h
main.o: /usr/include/unistd.h /usr/include/sys/unistd.h
main.o: /usr/include/sys/select.h /usr/include/sys/_select.h
main.o: /usr/include/dirent.h /usr/include/sys/dirent.h ../FB/FB_Define.h
main.o: ../FB/mydebug.h ../FB/TPControl.h /usr/include/math.h ../FB/dfiinfo.h
main.o: ../FB/DomainInfo.h LAYOUT_Define.h
layout.o: layout.h /usr/include/stdio.h /usr/include/sys/cdefs.h
layout.o: /usr/include/sys/_symbol_aliasing.h
layout.o: /usr/include/sys/_posix_availability.h /usr/include/Availability.h
layout.o: /usr/include/AvailabilityInternal.h /usr/include/_types.h
layout.o: /usr/include/sys/_types.h /usr/include/machine/_types.h
layout.o: /usr/include/i386/_types.h /usr/include/secure/_stdio.h
layout.o: /usr/include/secure/_common.h /usr/include/string.h
layout.o: /usr/include/strings.h /usr/include/secure/_string.h
layout.o: /usr/include/errno.h /usr/include/sys/errno.h
layout.o: /usr/include/sys/types.h /usr/include/sys/appleapiopts.h
layout.o: /usr/include/machine/types.h /usr/include/i386/types.h
layout.o: /usr/include/machine/endian.h /usr/include/i386/endian.h
layout.o: /usr/include/sys/_endian.h /usr/include/libkern/_OSByteOrder.h
layout.o: /usr/include/libkern/i386/_OSByteOrder.h
layout.o: /usr/include/sys/_structs.h /usr/include/sys/stat.h
layout.o: /usr/include/sys/param.h /usr/include/sys/syslimits.h
layout.o: /usr/include/machine/param.h /usr/include/i386/param.h
layout.o: /usr/include/i386/_param.h /usr/include/limits.h
layout.o: /usr/include/machine/limits.h /usr/include/i386/limits.h
layout.o: /usr/include/i386/_limits.h /usr/include/sys/signal.h
layout.o: /usr/include/machine/signal.h /usr/include/i386/signal.h
layout.o: /usr/include/i386/_structs.h /usr/include/stdlib.h
layout.o: /usr/include/sys/wait.h /usr/include/sys/resource.h
layout.o: /usr/include/alloca.h /usr/include/unistd.h
layout.o: /usr/include/sys/unistd.h /usr/include/sys/select.h
layout.o: /usr/include/sys/_select.h /usr/include/dirent.h
layout.o: /usr/include/sys/dirent.h ../FB/FB_Define.h ../FB/mydebug.h
layout.o: ../FB/TPControl.h /usr/include/math.h ../FB/dfiinfo.h
layout.o: ../FB/DomainInfo.h LAYOUT_Define.h
