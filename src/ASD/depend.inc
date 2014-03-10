SubDomain.o: SubDomain.C SubDomain.h ../FB/FB_Define.h ../FB/mydebug.h
ASDmodule.o: ASDmodule.C ASDmodule.h ../FB/DomainInfo.h \
 /usr/local/FFV/CPMlib/include/cpm_ParaManager.h \
 /usr/local/FFV/CPMlib/include/cpm_Base.h \
 /usr/local/FFV/CPMlib/include/cpm_Define.h /opt/openmpi/include/mpi.h \
 /opt/openmpi/include/mpi_portable_platform.h \
 /opt/openmpi/include/openmpi/ompi/mpi/cxx/mpicxx.h /opt/openmpi/include/mpi.h \
 /opt/openmpi/include/openmpi/ompi/mpi/cxx/constants.h \
 /opt/openmpi/include/openmpi/ompi/mpi/cxx/functions.h \
 /opt/openmpi/include/openmpi/ompi/mpi/cxx/datatype.h \
 /opt/openmpi/include/openmpi/ompi/mpi/cxx/exception.h \
 /opt/openmpi/include/openmpi/ompi/mpi/cxx/op.h \
 /opt/openmpi/include/openmpi/ompi/mpi/cxx/status.h \
 /opt/openmpi/include/openmpi/ompi/mpi/cxx/request.h \
 /opt/openmpi/include/openmpi/ompi/mpi/cxx/group.h \
 /opt/openmpi/include/openmpi/ompi/mpi/cxx/comm.h \
 /opt/openmpi/include/openmpi/ompi/mpi/cxx/win.h \
 /opt/openmpi/include/openmpi/ompi/mpi/cxx/file.h \
 /opt/openmpi/include/openmpi/ompi/mpi/cxx/errhandler.h \
 /opt/openmpi/include/openmpi/ompi/mpi/cxx/intracomm.h \
 /opt/openmpi/include/openmpi/ompi/mpi/cxx/topology.h \
 /opt/openmpi/include/openmpi/ompi/mpi/cxx/intercomm.h \
 /opt/openmpi/include/openmpi/ompi/mpi/cxx/info.h \
 /opt/openmpi/include/openmpi/ompi/mpi/cxx/datatype_inln.h \
 /opt/openmpi/include/openmpi/ompi/mpi/cxx/functions_inln.h \
 /opt/openmpi/include/openmpi/ompi/mpi/cxx/request_inln.h \
 /opt/openmpi/include/openmpi/ompi/mpi/cxx/comm_inln.h \
 /opt/openmpi/include/openmpi/ompi/mpi/cxx/intracomm_inln.h \
 /opt/openmpi/include/openmpi/ompi/mpi/cxx/topology_inln.h \
 /opt/openmpi/include/openmpi/ompi/mpi/cxx/intercomm_inln.h \
 /opt/openmpi/include/openmpi/ompi/mpi/cxx/group_inln.h \
 /opt/openmpi/include/openmpi/ompi/mpi/cxx/op_inln.h \
 /opt/openmpi/include/openmpi/ompi/mpi/cxx/errhandler_inln.h \
 /opt/openmpi/include/openmpi/ompi/mpi/cxx/status_inln.h \
 /opt/openmpi/include/openmpi/ompi/mpi/cxx/info_inln.h \
 /opt/openmpi/include/openmpi/ompi/mpi/cxx/win_inln.h \
 /opt/openmpi/include/openmpi/ompi/mpi/cxx/file_inln.h \
 /usr/local/FFV/CPMlib/include/cpm_Version.h \
 /usr/local/FFV/CPMlib/include/cpm_DomainInfo.h \
 /usr/local/FFV/CPMlib/include/cpm_EndianUtil.h \
 /usr/local/FFV/CPMlib/include/cpm_VoxelInfo.h \
 /usr/local/FFV/CPMlib/include/cpm_ObjList.h \
 /usr/local/FFV/CPMlib/include/inline/cpm_ParaManager_inline.h \
 /usr/local/FFV/CPMlib/include/inline/cpm_ParaManager_BndComm.h \
 /usr/local/FFV/CPMlib/include/inline/cpm_ParaManager_BndCommEx.h \
 ../FB/FB_Define.h ../FB/mydebug.h SubDomain.h ../FB/VoxInfo.h ../FB/FBUtility.h \
 ../FB/Medium.h ../FB/Component.h ../FB/SetBC.h ../FB/BndOuter.h ../FB/Control.h \
 ../FB/CompoFraction.h ../FB/vec3.h ../FB/vec3_func.h ../FB/vec3f_func.h \
 ../FB/IterationControl.h /usr/local/FFV/TextParser/include/TextParser.h \
 /usr/local/FFV/TextParser/include/TextParserCommon.h \
 /usr/local/FFV/TextParser/include/tpVersion.h ../FB/IntervalManager.h \
 ../FB/Intrinsic.h /usr/local/FFV/Polylib/include/Polylib.h \
 /usr/local/FFV/Polylib/include/polygons/Polygons.h \
 /usr/local/FFV/Polylib/include/common/PolylibStat.h \
 /usr/local/FFV/Polylib/include/common/PolylibCommon.h \
 /usr/local/FFV/Polylib/include/common/PolylibDefine.h \
 /usr/local/FFV/Polylib/include/common/Vec3.h \
 /usr/local/FFV/Polylib/include/common/vec3_func.h \
 /usr/local/FFV/Polylib/include/common/vec3f_func.h \
 /usr/local/FFV/Polylib/include/common/axis.h \
 /usr/local/FFV/Polylib/include/polygons/TriMesh.h \
 /usr/local/FFV/Polylib/include/common/BBox.h \
 /usr/local/FFV/Polylib/include/common/Vec2.h \
 /usr/local/FFV/Polylib/include/polygons/Triangle.h \
 /usr/local/FFV/Polylib/include/polygons/Vertex.h \
 /usr/local/FFV/Polylib/include/polygons/DVertexTriangle.h \
 /usr/local/FFV/Polylib/include/polygons/PrivateTriangle.h \
 /usr/local/FFV/Polylib/include/groups/PolygonGroup.h \
 /usr/local/FFV/Polylib/include/groups/PolygonGroupFactory.h \
 /usr/local/FFV/Polylib/include/polyVersion.h \
 /usr/local/FFV/Polylib/include/MPIPolylib.h \
 /usr/local/FFV/Polylib/include/polygons/DVertex.h \
 /usr/local/FFV/Polylib/include/polygons/DVertexManager.h \
 /usr/local/FFV/Polylib/include/file_io/TriMeshIO.h \
 /usr/local/FFV/Polylib/include/file_io/stl.h \
 /usr/local/FFV/Polylib/include/polygons/VertexList.h \
 /usr/local/FFV/Polylib/include/file_io/obj.h \
 /usr/local/FFV/Polylib/include/common/tt.h \
 /usr/local/FFV/Polylib/include/file_io/vtk.h \
 /usr/local/FFV/Cutlib/include/Cutlib.h \
 /usr/local/FFV/Cutlib/include/RepairPolygonData/RepairPolygonData.h \
 /usr/local/FFV/Cutlib/include/CutInfo/CutInfoArray.h \
 /usr/local/FFV/Cutlib/include/CutInfo/CutInfo.h \
 /usr/local/FFV/Cutlib/include/CutInfo/CutNormalArray.h \
 /usr/local/FFV/Cutlib/include/Cutlib.h \
 /usr/local/FFV/Cutlib/include/GridAccessor/GridAccessor.h \
 /usr/local/FFV/Cutlib/include/Version.h \
 /usr/local/FFV/Cutlib/include/GridAccessor/Cell.h
