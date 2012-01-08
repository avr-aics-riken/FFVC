#include <math.h>
#include <iostream>
#include <algorithm>
#include <tt/util/message.h>
#include <tt/util/time.h>
#include <tt/shapes/mesh.h>
#include "voxelizer_octree.h"
#include "grid.h"
#include "octree.h"
#include "vof.h"

using namespace std;

VxVoxelizerOctree::VxVoxelizerOctree(const BBox& bbox, int max_depth)
{
	int nx, ny, nz;
	m_bbox = bbox;
	nx = ny = nz = int(powf(2, max_depth));
	m_nelm[0] = nx;
	m_nelm[1] = ny;
	m_nelm[2] = nz;
	m_pitch[0] = bbox.length(0) / nx;
	m_pitch[1] = bbox.length(1) / ny;
	m_pitch[2] = bbox.length(2) / nz;
	m_offset[0] = m_pitch[0] * .5;
	m_offset[1] = m_pitch[1] * .5;
	m_offset[2] = m_pitch[2] * .5;
	m_octree = new VxOctree(m_bbox, max_depth);
	m_nid = 0;
	m_fluid_id = 0;
	m_ratio_type = VX_RATIO_PLANE;
}

VxVoxelizerOctree::~VxVoxelizerOctree()
{
	delete m_octree;
}

int
VxVoxelizerOctree::voxelizeSurface(const TriMesh* mesh, int id)
{
	tt_info("[vxoctree] voxelize surface...(start ID %d)\n", id);
	m_mesh = mesh;
	m_octree->setNumId(id);
	m_octree->build(mesh);
	m_octree->smoothLevel();
	id += mesh->nGroups();
	tt_info("[vxoctree] voxelize surface...done (end ID %d)\n", id);
	return id;
}

int
VxVoxelizerOctree::voxelizeSurface(const TriMesh* mesh, std::vector<int>& indices, int id)
{
	m_mesh = mesh;
	m_octree->setNumId(id);
	m_octree->build(mesh, indices);
	m_octree->smoothLevel();
	id += mesh->nGroups();
	return id;
}

#define MAX_ID 255
int
VxVoxelizerOctree::fillOctreeByNewId(int id)
{
	tt_info("[vxoctree] fill octree by new id...(start ID %d)\n", id);
	TtTime tm;
	tm.start();

	int nx = m_nelm[0];
	int ny = m_nelm[1];
	int nz = m_nelm[2];
	int min_cnt = int(log10(nx) * log10(ny) * log10(nz) + 1);

	//cout << "min_cnt: " << min_cnt << endl;

	std::vector<VxOctreeNode*> leafs;
	VxOctreeNode* node;
	m_octree->getLeafNodes(leafs);
	int sz = leafs.size();
	for (int i=0; i<sz; i++)
	{
		node = leafs[i];
		if (!node->getMark())
		{
			int cnt = m_octree->fill(node, id);
			//cout << "id, cnt: " << id << " " << cnt << endl;
			if (id < MAX_ID && cnt > min_cnt) id++;
		}
	}
	id++;

	tm.end();
	tt_info("[vxoctree] fill octree by new id... done (end ID %d, %.2fsec)\n", id, tm.getElapsedSec());
	return id;
}

void
VxVoxelizerOctree::voxelize(const TriMesh* mesh, int id)
{
	id = voxelizeSurface(mesh, id);
	m_fluid_id = id;
	id = fillOctreeByNewId(id);
	setNumId(id);
}

VxOctreeNode*
VxVoxelizerOctree::getLeafNode(int x, int y, int z) const
{
	VxOctreePedigree pedig;
	int depth = m_octree->getMaxDepth();
	pdg_xyz2pedigree(&pedig, x, y, z, depth);
	return m_octree->getNode(pedig);
}

int
VxVoxelizerOctree::getId(int x, int y, int z) const
{
	VxOctreeNode* node = getLeafNode(x, y, z);
	return node->getId();
}

// x, y, z: position of the voxel
// volume: volume ratio
// area[6]: area ratio
void
VxVoxelizerOctree::computeRatio(int x, int y, int z, float* volume, float area[6]) const
{
	VxOctreeNode* node = getLeafNode(x, y, z);
	computeRatio(node, volume, area);
}

void
VxVoxelizerOctree::computeRatio(VxOctreeNode* node, float* volume, float area[6]) const
{
	*volume = 0;
	for (int i=0; i<6; i++) area[i] = 0;
	switch (m_ratio_type)
	{
		case VX_RATIO_BINARY:
			{
				int id = node->getId();
				if (id == m_fluid_id)
				{
					*volume = 1;
					for (int i=0; i<6; i++) area[i] = 1;
				}
			}
			break;
		case VX_RATIO_PLANE:
			{
				if (node->isBoundary())
				{
					f_compute_ratio(m_mesh, m_octree, node, volume, area);
				}
				else
				{
					int id = node->getId();
					if (id == m_fluid_id)
					{
						*volume = 1;
						for (int i=0; i<6; i++) area[i] = 1;
					}
				}
			}
			break;
	}
}

Vec3f
VxVoxelizerOctree::computeAvgNormalForVoxel(int x, int y, int z) const
{
	VxOctreeNode* node = getLeafNode(x, y, z);
	if (node->isBoundary())
	{
		VxOctreeBoundaryData* data = node->getBoundaryData();
		int nelm = data->m_nelm;
		int* elms = data->m_elms;
		return computeAvgNormal(m_mesh, elms, nelm);
	}
	return Vec3f(0, 0, 0);
}

