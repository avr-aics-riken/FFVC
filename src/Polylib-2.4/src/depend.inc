VTree.o: ../include/polygons/VTree.h ../include/common/BBox.h \
  ../include/common/Vec2.h ../include/common/Vec3.h \
  ../include/common/vec3_func.h ../include/common/vec3f_func.h \
  ../include/common/axis.h ../include/common/PolylibStat.h \
  ../include/common/PolylibCommon.h
Polygons.o: ../include/polygons/Polygons.h ../include/polygons/VTree.h \
  ../include/common/BBox.h ../include/common/Vec2.h \
  ../include/common/Vec3.h ../include/common/vec3_func.h \
  ../include/common/vec3f_func.h ../include/common/axis.h \
  ../include/common/PolylibStat.h ../include/common/PolylibCommon.h \
  ../include/common/tt.h ../include/common/Vec3.h
TriMesh.o: ../include/polygons/TriMesh.h
TriMeshIO.o: ../include/file_io/TriMeshIO.h \
  ../include/common/PolylibStat.h ../include/common/PolylibCommon.h
stl.o: ../include/file_io/stl.h ../include/common/PolylibCommon.h
triangle_id.o: ../include/file_io/triangle_id.h \
  ../include/polygons/Triangle.h ../include/common/Vec3.h \
  ../include/common/vec3_func.h ../include/common/vec3f_func.h \
  ../include/common/axis.h ../include/common/PolylibStat.h \
  ../include/common/PolylibCommon.h
CPolylib.o: ../include/c_lang/CPolylib.h ../include/common/PolylibStat.h
time.o: ../include/util/time.h
