#ifndef voxelizer_h
#define voxelizer_h

#include <vector>
#include <tt/gfx/bbox.h>

class TriMesh;

class VxVoxelizer
{
public:
	virtual ~VxVoxelizer() {}

	virtual void getOrigin(float* ox, float* oy, float* oz) const = 0;
	virtual void getSize(int* nx, int* ny, int* nz) const = 0;
	virtual void getPitch(float* px, float* py, float* pz) const = 0;
	virtual void setNumId(int nr) = 0;
	virtual int  getNumId() const = 0;

	virtual const BBox& getBBox() const = 0;

	virtual int  voxelizeSurface(const TriMesh* mesh, int id) = 0;
	virtual int  voxelizeSurface(const TriMesh* mesh, std::vector<int>& indices, int id) = 0;

	virtual int  getId(int x, int y, int z) const = 0;
	virtual Vec3f computeAvgNormalForVoxel(int x, int y, int z) const = 0;
};

#endif  // voxelizer_h

