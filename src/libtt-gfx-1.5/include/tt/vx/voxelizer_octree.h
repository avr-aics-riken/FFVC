#ifndef vxvoxelizer_octree_h
#define vxvoxelizer_octree_h

#include <vector>
#include <tt/gfx/bbox.h>
#include "voxelizer.h"

class VxOctree;
class VxOctreeNode;
class TriMesh;

enum
{
	VX_RATIO_BINARY = 0,
	VX_RATIO_PLANE,
};

class VxVoxelizerOctree : public VxVoxelizer
{
public:
	VxVoxelizerOctree(const BBox& bbox, int max_depth);
	~VxVoxelizerOctree();

	void getOrigin(float* ox, float* oy, float* oz) const {
		*ox = m_bbox.min[0]; *oy = m_bbox.min[1]; *oz = m_bbox.min[2];
	}
	void getSize(int* nx, int* ny, int* nz) const {
		*nx = m_nelm[0]; *ny = m_nelm[1]; *nz = m_nelm[2];
	}
	void getPitch(float* px, float* py, float* pz) const {
		*px = m_pitch[0]; *py = m_pitch[1]; *pz = m_pitch[2];
	}
	const BBox& getBBox() const { return m_bbox; }
	VxOctree* getOctree() const { return m_octree; }

	void setNumId(int nr) { m_nid = nr; }
	int  getNumId() const { return m_nid; }

	void voxelize(const TriMesh* mesh, int id);

	// access by a grid index: x, y, z
	VxOctreeNode* getLeafNode(int x, int y, int z) const;
	int  getId(int x, int y, int z) const;
	void computeRatio(int x, int y, int z, float* volume, float area[6]) const;
	Vec3f computeAvgNormalForVoxel(int x, int y, int z) const;

	// access by a node
	void computeRatio(VxOctreeNode* node, float* volume, float area[6]) const;

	void setRatioType(int type) { m_ratio_type = type; }

private:
	int  voxelizeSurface(const TriMesh* mesh, int id);
	int  voxelizeSurface(const TriMesh* mesh, std::vector<int>& indices, int id);

	int  fillOctreeByNewId(int id);

	const TriMesh* m_mesh;
	VxOctree* m_octree;
	BBox m_bbox;
	int m_nelm[3];
	float m_pitch[3];
	float m_offset[3];
	int m_nid;
	int m_fluid_id;
	int m_ratio_type;
};

#endif  // vxvoxelizer_octree_h

