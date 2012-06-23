###################################################################
#
# CAERU Library
#
# Copyright (c) All right reserved. 2012
#
# Institute of Industrial Science, The University of Tokyo, Japan. 
#
###################################################################

TARGET = libIP.a

CSRCS = \


CXXSRCS = \
    IP_PMT.C \
    IP_Rect.C \
    IP_Polygon.C \
    IP_Cylinder.C \
    IP_Step.C \
    IP_PPLT2D.C \
    IP_Duct.C \
    IP_SHC1D.C \
    IP_Sphere.C \


F90SRCS = \


SRCS  = $(CSRCS) $(CXXSRCS) $(F90SRCS)

.SUFFIXES: .o .c .C .f90
COBJS = $(CSRCS:.c=.o)
CXXOBJS = $(CXXSRCS:.C=.o)
F90OBJS = $(F90SRCS:.f90=.o)
OBJS  = $(COBJS) $(CXXOBJS) $(F90OBJS)

$(TARGET):$(OBJS)
	$(AR) $(TARGET) $(OBJS)
	$(RANLIB) $(TARGET)

.c.o:
	$(CC) $(CFLAGS) $(UDEF_INC_PATH) -o $@ -c $<

.C.o:
	$(CXX) $(CXXFLAGS) $(UDEF_INC_PATH) -o $@ -c $<

.f90.o:
	$(F90) $(F90FLAGS) -o $@ -c $<

clean:
	$(RM) $(OBJS) $(TARGET)

depend:
	makedepend -I. -Y $(SRCS)

# DO NOT DELETE

IP_PMT.o: IP_PMT.h ../FB/Intrinsic.h ../FB/FB_Define.h ../FB/mydebug.h
IP_PMT.o: ../FB/BndOuter.h ../FB/Control.h ../FB/Medium.h ../FB/Component.h
IP_PMT.o: ../FB/FBUtility.h ../FB/Interval_Mngr.h ../FB/TPControl.h
IP_PMT.o: ../FB/vec3.h ../FB/vec3_func.h ../FB/vec3f_func.h IP_Define.h
IP_Rect.o: IP_Rect.h ../FB/Intrinsic.h ../FB/FB_Define.h ../FB/mydebug.h
IP_Rect.o: ../FB/BndOuter.h ../FB/Control.h ../FB/Medium.h ../FB/Component.h
IP_Rect.o: ../FB/FBUtility.h ../FB/Interval_Mngr.h ../FB/TPControl.h
IP_Rect.o: ../FB/vec3.h ../FB/vec3_func.h ../FB/vec3f_func.h IP_Define.h
IP_Polygon.o: IP_Polygon.h ../FB/Intrinsic.h ../FB/FB_Define.h
IP_Polygon.o: ../FB/mydebug.h ../FB/BndOuter.h ../FB/Control.h ../FB/Medium.h
IP_Polygon.o: ../FB/Component.h ../FB/FBUtility.h ../FB/Interval_Mngr.h
IP_Polygon.o: ../FB/TPControl.h ../FB/vec3.h ../FB/vec3_func.h
IP_Polygon.o: ../FB/vec3f_func.h IP_Define.h
IP_Cylinder.o: IP_Cylinder.h ../FB/Intrinsic.h ../FB/FB_Define.h
IP_Cylinder.o: ../FB/mydebug.h ../FB/BndOuter.h ../FB/Control.h
IP_Cylinder.o: ../FB/Medium.h ../FB/Component.h ../FB/FBUtility.h
IP_Cylinder.o: ../FB/Interval_Mngr.h ../FB/TPControl.h ../FB/vec3.h
IP_Cylinder.o: ../FB/vec3_func.h ../FB/vec3f_func.h IP_Define.h
IP_Step.o: IP_Step.h ../FB/Intrinsic.h ../FB/FB_Define.h ../FB/mydebug.h
IP_Step.o: ../FB/BndOuter.h ../FB/Control.h ../FB/Medium.h ../FB/Component.h
IP_Step.o: ../FB/FBUtility.h ../FB/Interval_Mngr.h ../FB/TPControl.h
IP_Step.o: ../FB/vec3.h ../FB/vec3_func.h ../FB/vec3f_func.h IP_Define.h
IP_PPLT2D.o: IP_PPLT2D.h ../FB/Intrinsic.h ../FB/FB_Define.h ../FB/mydebug.h
IP_PPLT2D.o: ../FB/BndOuter.h ../FB/Control.h ../FB/Medium.h
IP_PPLT2D.o: ../FB/Component.h ../FB/FBUtility.h ../FB/Interval_Mngr.h
IP_PPLT2D.o: ../FB/TPControl.h ../FB/vec3.h ../FB/vec3_func.h
IP_PPLT2D.o: ../FB/vec3f_func.h IP_Define.h
IP_Duct.o: IP_Duct.h ../FB/Intrinsic.h ../FB/FB_Define.h ../FB/mydebug.h
IP_Duct.o: ../FB/BndOuter.h ../FB/Control.h ../FB/Medium.h ../FB/Component.h
IP_Duct.o: ../FB/FBUtility.h ../FB/Interval_Mngr.h ../FB/TPControl.h
IP_Duct.o: ../FB/vec3.h ../FB/vec3_func.h ../FB/vec3f_func.h IP_Define.h
IP_SHC1D.o: IP_SHC1D.h ../FB/Intrinsic.h ../FB/FB_Define.h ../FB/mydebug.h
IP_SHC1D.o: ../FB/BndOuter.h ../FB/Control.h ../FB/Medium.h ../FB/Component.h
IP_SHC1D.o: ../FB/FBUtility.h ../FB/Interval_Mngr.h ../FB/TPControl.h
IP_SHC1D.o: ../FB/vec3.h ../FB/vec3_func.h ../FB/vec3f_func.h IP_Define.h
IP_Sphere.o: IP_Sphere.h ../FB/Intrinsic.h ../FB/FB_Define.h ../FB/mydebug.h
IP_Sphere.o: ../FB/BndOuter.h ../FB/Control.h ../FB/Medium.h
IP_Sphere.o: ../FB/Component.h ../FB/FBUtility.h ../FB/Interval_Mngr.h
IP_Sphere.o: ../FB/TPControl.h ../FB/vec3.h ../FB/vec3_func.h
IP_Sphere.o: ../FB/vec3f_func.h IP_Define.h
