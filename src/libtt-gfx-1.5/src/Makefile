CC = icc
CXX = icpc
AR = ar cru
CFLAGS = -g -O2 -I../include
RANLIB = ranlib

TARGET = libtt-gfx.a
UTIL_SRCS = \
	util/string.cxx \
	util/message.cxx \
	util/filename.cxx \
	util/file.cxx \

GFX_SRCS = \
	core/table.cxx \
	core/mat4x4f.cxx \
	gfx/itsect_box_tri.cxx \
	gfx/shape.cxx \
	gfx/transform.cxx \
	shapes/mesh.cxx \
	shapes/mesh/mesh_edit.cxx \
	shapes/mesh/mesh_io.cxx \
	shapes/mesh/vtree.cxx \
	shapes/mesh/halfedge.cxx \
	shapes/mesh/fmt/obj.cxx \
	shapes/mesh/fmt/off.cxx \
	shapes/mesh/fmt/ply_a.cxx \
	shapes/mesh/fmt/ply_b.cxx \
	shapes/mesh/fmt/stl_a.cxx \
	shapes/mesh/fmt/stl_b.cxx \
	shapes/mesh/fmt/tm.cxx \
	vx/grid.cxx \
	vx/octree.cxx \
	vx/voxelizer_octree.cxx \
	vx/vof.cxx \
	vx/itsect_axis.cxx \
	vx/fmt/svx.cxx \
	vx/fmt/ovx.cxx \
	vx/fmt/sph.cxx \

SRCS = $(UTIL_SRCS) $(GFX_SRCS)
OBJS = $(addsuffix .o, $(basename $(SRCS)))

$(TARGET): $(OBJS)
	$(AR) $@ $^
	$(RANLIB) $(TARGET)

%.o: %.cxx
	$(CXX) $(CFLAGS) -o $@ -c $<

clean:
	rm -f $(TARGET) $(OBJS)

