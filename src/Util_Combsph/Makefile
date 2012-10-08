###################################################################
#
# CAERU Library
#
# Copyright (c) All right reserved. 2012
#
# Institute of Industrial Science, The University of Tokyo, Japan. 
#
###################################################################

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
	make -f MakeCOMB \
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
	(make -f MakeCOMB clean)
  
allclean: clean
	(cd ../FB; make clean)
	(cd $(PMLIB_DIR)/src; make clean)

depend:
	(cd $(PMLIB_DIR)/src; make depend)
	(cd ../FB; make depend)
	(make -f MakeCOMB depend)

# DO NOT DELETE

main.o: comb.h /usr/include/stdio.h /usr/include/sys/cdefs.h
main.o: /usr/include/sys/_symbol_aliasing.h
main.o: /usr/include/sys/_posix_availability.h /usr/include/Availability.h
main.o: /usr/include/AvailabilityInternal.h /usr/include/_types.h
main.o: /usr/include/sys/_types.h /usr/include/machine/_types.h
main.o: /usr/include/i386/_types.h /usr/include/secure/_stdio.h
main.o: /usr/include/secure/_common.h /usr/include/stdlib.h
main.o: /usr/include/sys/wait.h /usr/include/sys/signal.h
main.o: /usr/include/sys/appleapiopts.h /usr/include/machine/signal.h
main.o: /usr/include/i386/signal.h /usr/include/i386/_structs.h
main.o: /usr/include/sys/_structs.h /usr/include/machine/_structs.h
main.o: /usr/include/sys/resource.h /usr/include/machine/endian.h
main.o: /usr/include/i386/endian.h /usr/include/sys/_endian.h
main.o: /usr/include/libkern/_OSByteOrder.h
main.o: /usr/include/libkern/i386/_OSByteOrder.h /usr/include/alloca.h
main.o: /usr/include/machine/types.h /usr/include/i386/types.h
main.o: /usr/include/string.h /usr/include/strings.h
main.o: /usr/include/secure/_string.h /usr/include/errno.h
main.o: /usr/include/sys/errno.h /usr/include/dirent.h
main.o: /usr/include/sys/dirent.h ../FB/FB_Define.h ../FB/mydebug.h
main.o: ../FB/FBUtility.h /usr/include/math.h /usr/include/sys/stat.h
main.o: ../FB/TPControl.h ../FB/dfi.h ../FB/PLOT3D_read.h ../FB/PLOT3D.h
main.o: /usr/include/time.h /usr/include/_structs.h /usr/include/fcntl.h
main.o: /usr/include/sys/fcntl.h ../FB/DomainInfo.h ../FB/FB_Ffunc.h
main.o: ../FB/PLOT3D_write.h FileIO_sph.h endianUtil.h
main.o: /usr/include/sys/types.h /usr/include/unistd.h
main.o: /usr/include/sys/unistd.h /usr/include/sys/select.h
main.o: /usr/include/sys/_select.h ../FB/dfiinfo.h /usr/include/limits.h
main.o: /usr/include/machine/limits.h /usr/include/i386/limits.h
main.o: /usr/include/i386/_limits.h /usr/include/sys/syslimits.h
main.o: COMB_Define.h
comb.o: comb.h /usr/include/stdio.h /usr/include/sys/cdefs.h
comb.o: /usr/include/sys/_symbol_aliasing.h
comb.o: /usr/include/sys/_posix_availability.h /usr/include/Availability.h
comb.o: /usr/include/AvailabilityInternal.h /usr/include/_types.h
comb.o: /usr/include/sys/_types.h /usr/include/machine/_types.h
comb.o: /usr/include/i386/_types.h /usr/include/secure/_stdio.h
comb.o: /usr/include/secure/_common.h /usr/include/stdlib.h
comb.o: /usr/include/sys/wait.h /usr/include/sys/signal.h
comb.o: /usr/include/sys/appleapiopts.h /usr/include/machine/signal.h
comb.o: /usr/include/i386/signal.h /usr/include/i386/_structs.h
comb.o: /usr/include/sys/_structs.h /usr/include/machine/_structs.h
comb.o: /usr/include/sys/resource.h /usr/include/machine/endian.h
comb.o: /usr/include/i386/endian.h /usr/include/sys/_endian.h
comb.o: /usr/include/libkern/_OSByteOrder.h
comb.o: /usr/include/libkern/i386/_OSByteOrder.h /usr/include/alloca.h
comb.o: /usr/include/machine/types.h /usr/include/i386/types.h
comb.o: /usr/include/string.h /usr/include/strings.h
comb.o: /usr/include/secure/_string.h /usr/include/errno.h
comb.o: /usr/include/sys/errno.h /usr/include/dirent.h
comb.o: /usr/include/sys/dirent.h ../FB/FB_Define.h ../FB/mydebug.h
comb.o: ../FB/FBUtility.h /usr/include/math.h /usr/include/sys/stat.h
comb.o: ../FB/TPControl.h ../FB/dfi.h ../FB/PLOT3D_read.h ../FB/PLOT3D.h
comb.o: /usr/include/time.h /usr/include/_structs.h /usr/include/fcntl.h
comb.o: /usr/include/sys/fcntl.h ../FB/DomainInfo.h ../FB/FB_Ffunc.h
comb.o: ../FB/PLOT3D_write.h FileIO_sph.h endianUtil.h
comb.o: /usr/include/sys/types.h /usr/include/unistd.h
comb.o: /usr/include/sys/unistd.h /usr/include/sys/select.h
comb.o: /usr/include/sys/_select.h ../FB/dfiinfo.h /usr/include/limits.h
comb.o: /usr/include/machine/limits.h /usr/include/i386/limits.h
comb.o: /usr/include/i386/_limits.h /usr/include/sys/syslimits.h
comb.o: COMB_Define.h
comb_sph.o: comb.h /usr/include/stdio.h /usr/include/sys/cdefs.h
comb_sph.o: /usr/include/sys/_symbol_aliasing.h
comb_sph.o: /usr/include/sys/_posix_availability.h
comb_sph.o: /usr/include/Availability.h /usr/include/AvailabilityInternal.h
comb_sph.o: /usr/include/_types.h /usr/include/sys/_types.h
comb_sph.o: /usr/include/machine/_types.h /usr/include/i386/_types.h
comb_sph.o: /usr/include/secure/_stdio.h /usr/include/secure/_common.h
comb_sph.o: /usr/include/stdlib.h /usr/include/sys/wait.h
comb_sph.o: /usr/include/sys/signal.h /usr/include/sys/appleapiopts.h
comb_sph.o: /usr/include/machine/signal.h /usr/include/i386/signal.h
comb_sph.o: /usr/include/i386/_structs.h /usr/include/sys/_structs.h
comb_sph.o: /usr/include/machine/_structs.h /usr/include/sys/resource.h
comb_sph.o: /usr/include/machine/endian.h /usr/include/i386/endian.h
comb_sph.o: /usr/include/sys/_endian.h /usr/include/libkern/_OSByteOrder.h
comb_sph.o: /usr/include/libkern/i386/_OSByteOrder.h /usr/include/alloca.h
comb_sph.o: /usr/include/machine/types.h /usr/include/i386/types.h
comb_sph.o: /usr/include/string.h /usr/include/strings.h
comb_sph.o: /usr/include/secure/_string.h /usr/include/errno.h
comb_sph.o: /usr/include/sys/errno.h /usr/include/dirent.h
comb_sph.o: /usr/include/sys/dirent.h ../FB/FB_Define.h ../FB/mydebug.h
comb_sph.o: ../FB/FBUtility.h /usr/include/math.h /usr/include/sys/stat.h
comb_sph.o: ../FB/TPControl.h ../FB/dfi.h ../FB/PLOT3D_read.h ../FB/PLOT3D.h
comb_sph.o: /usr/include/time.h /usr/include/_structs.h /usr/include/fcntl.h
comb_sph.o: /usr/include/sys/fcntl.h ../FB/DomainInfo.h ../FB/FB_Ffunc.h
comb_sph.o: ../FB/PLOT3D_write.h FileIO_sph.h endianUtil.h
comb_sph.o: /usr/include/sys/types.h /usr/include/unistd.h
comb_sph.o: /usr/include/sys/unistd.h /usr/include/sys/select.h
comb_sph.o: /usr/include/sys/_select.h ../FB/dfiinfo.h /usr/include/limits.h
comb_sph.o: /usr/include/machine/limits.h /usr/include/i386/limits.h
comb_sph.o: /usr/include/i386/_limits.h /usr/include/sys/syslimits.h
comb_sph.o: COMB_Define.h
comb_plot3d.o: comb.h /usr/include/stdio.h /usr/include/sys/cdefs.h
comb_plot3d.o: /usr/include/sys/_symbol_aliasing.h
comb_plot3d.o: /usr/include/sys/_posix_availability.h
comb_plot3d.o: /usr/include/Availability.h
comb_plot3d.o: /usr/include/AvailabilityInternal.h /usr/include/_types.h
comb_plot3d.o: /usr/include/sys/_types.h /usr/include/machine/_types.h
comb_plot3d.o: /usr/include/i386/_types.h /usr/include/secure/_stdio.h
comb_plot3d.o: /usr/include/secure/_common.h /usr/include/stdlib.h
comb_plot3d.o: /usr/include/sys/wait.h /usr/include/sys/signal.h
comb_plot3d.o: /usr/include/sys/appleapiopts.h /usr/include/machine/signal.h
comb_plot3d.o: /usr/include/i386/signal.h /usr/include/i386/_structs.h
comb_plot3d.o: /usr/include/sys/_structs.h /usr/include/machine/_structs.h
comb_plot3d.o: /usr/include/sys/resource.h /usr/include/machine/endian.h
comb_plot3d.o: /usr/include/i386/endian.h /usr/include/sys/_endian.h
comb_plot3d.o: /usr/include/libkern/_OSByteOrder.h
comb_plot3d.o: /usr/include/libkern/i386/_OSByteOrder.h /usr/include/alloca.h
comb_plot3d.o: /usr/include/machine/types.h /usr/include/i386/types.h
comb_plot3d.o: /usr/include/string.h /usr/include/strings.h
comb_plot3d.o: /usr/include/secure/_string.h /usr/include/errno.h
comb_plot3d.o: /usr/include/sys/errno.h /usr/include/dirent.h
comb_plot3d.o: /usr/include/sys/dirent.h ../FB/FB_Define.h ../FB/mydebug.h
comb_plot3d.o: ../FB/FBUtility.h /usr/include/math.h /usr/include/sys/stat.h
comb_plot3d.o: ../FB/TPControl.h ../FB/dfi.h ../FB/PLOT3D_read.h
comb_plot3d.o: ../FB/PLOT3D.h /usr/include/time.h /usr/include/_structs.h
comb_plot3d.o: /usr/include/fcntl.h /usr/include/sys/fcntl.h
comb_plot3d.o: ../FB/DomainInfo.h ../FB/FB_Ffunc.h ../FB/PLOT3D_write.h
comb_plot3d.o: FileIO_sph.h endianUtil.h /usr/include/sys/types.h
comb_plot3d.o: /usr/include/unistd.h /usr/include/sys/unistd.h
comb_plot3d.o: /usr/include/sys/select.h /usr/include/sys/_select.h
comb_plot3d.o: ../FB/dfiinfo.h /usr/include/limits.h
comb_plot3d.o: /usr/include/machine/limits.h /usr/include/i386/limits.h
comb_plot3d.o: /usr/include/i386/_limits.h /usr/include/sys/syslimits.h
comb_plot3d.o: COMB_Define.h
FileIO_sph.o: FileIO_sph.h /usr/include/stdio.h /usr/include/sys/cdefs.h
FileIO_sph.o: /usr/include/sys/_symbol_aliasing.h
FileIO_sph.o: /usr/include/sys/_posix_availability.h
FileIO_sph.o: /usr/include/Availability.h /usr/include/AvailabilityInternal.h
FileIO_sph.o: /usr/include/_types.h /usr/include/sys/_types.h
FileIO_sph.o: /usr/include/machine/_types.h /usr/include/i386/_types.h
FileIO_sph.o: /usr/include/secure/_stdio.h /usr/include/secure/_common.h
FileIO_sph.o: /usr/include/string.h /usr/include/strings.h
FileIO_sph.o: /usr/include/secure/_string.h endianUtil.h
FileIO_sph.o: /usr/include/sys/types.h /usr/include/sys/appleapiopts.h
FileIO_sph.o: /usr/include/machine/types.h /usr/include/i386/types.h
FileIO_sph.o: /usr/include/machine/endian.h /usr/include/i386/endian.h
FileIO_sph.o: /usr/include/sys/_endian.h /usr/include/libkern/_OSByteOrder.h
FileIO_sph.o: /usr/include/libkern/i386/_OSByteOrder.h
FileIO_sph.o: /usr/include/sys/_structs.h /usr/include/machine/_structs.h
FileIO_sph.o: /usr/include/i386/_structs.h /usr/include/unistd.h
FileIO_sph.o: /usr/include/sys/unistd.h /usr/include/sys/select.h
FileIO_sph.o: /usr/include/sys/_select.h ../FB/FB_Define.h ../FB/mydebug.h
FileIO_sph.o: /usr/include/stdlib.h /usr/include/sys/wait.h
FileIO_sph.o: /usr/include/sys/signal.h /usr/include/machine/signal.h
FileIO_sph.o: /usr/include/i386/signal.h /usr/include/sys/resource.h
FileIO_sph.o: /usr/include/alloca.h
FileIO_read_sph.o: FileIO_sph.h /usr/include/stdio.h /usr/include/sys/cdefs.h
FileIO_read_sph.o: /usr/include/sys/_symbol_aliasing.h
FileIO_read_sph.o: /usr/include/sys/_posix_availability.h
FileIO_read_sph.o: /usr/include/Availability.h
FileIO_read_sph.o: /usr/include/AvailabilityInternal.h /usr/include/_types.h
FileIO_read_sph.o: /usr/include/sys/_types.h /usr/include/machine/_types.h
FileIO_read_sph.o: /usr/include/i386/_types.h /usr/include/secure/_stdio.h
FileIO_read_sph.o: /usr/include/secure/_common.h /usr/include/string.h
FileIO_read_sph.o: /usr/include/strings.h /usr/include/secure/_string.h
FileIO_read_sph.o: endianUtil.h /usr/include/sys/types.h
FileIO_read_sph.o: /usr/include/sys/appleapiopts.h
FileIO_read_sph.o: /usr/include/machine/types.h /usr/include/i386/types.h
FileIO_read_sph.o: /usr/include/machine/endian.h /usr/include/i386/endian.h
FileIO_read_sph.o: /usr/include/sys/_endian.h
FileIO_read_sph.o: /usr/include/libkern/_OSByteOrder.h
FileIO_read_sph.o: /usr/include/libkern/i386/_OSByteOrder.h
FileIO_read_sph.o: /usr/include/sys/_structs.h
FileIO_read_sph.o: /usr/include/machine/_structs.h
FileIO_read_sph.o: /usr/include/i386/_structs.h /usr/include/unistd.h
FileIO_read_sph.o: /usr/include/sys/unistd.h /usr/include/sys/select.h
FileIO_read_sph.o: /usr/include/sys/_select.h ../FB/FB_Define.h
FileIO_read_sph.o: ../FB/mydebug.h /usr/include/stdlib.h
FileIO_read_sph.o: /usr/include/sys/wait.h /usr/include/sys/signal.h
FileIO_read_sph.o: /usr/include/machine/signal.h /usr/include/i386/signal.h
FileIO_read_sph.o: /usr/include/sys/resource.h /usr/include/alloca.h
