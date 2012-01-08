#include <iostream>
#include <queue>
#include <tt/util/message.h>
#include <tt/util/time.h>
#include <tt/util/string.h>
#include <tt/gfx/itsect.h>
#include <tt/shapes/mesh.h>
#include "octree.h"

#define VX_PARALLEL_DEG 1
#define VX_OVERLAP_SCALE 1e-4

using namespace std;

// if a point is on a grid
bool vx_check_point_on_grid(const float orig[3], const float pitch[3], const float v[3])
{
	for (int i=0; i<3; i++)
	{
		float s = gmAbs(v[i]-orig[i]);
		int n = s/pitch[i];
		float d = s - n*pitch[i];
		float e = pitch[i] * VX_OVERLAP_SCALE;
		if (gmFuzEQ(d, 0, e) || gmFuzEQ(d, 1, e))
			return 1;
	}
	return 0;
}

// if a plane is on a box
bool vx_check_plane_on_box(const float min[3], const float max[3], const float nml[3], const float v[3], int* face, int* dir)
{
	double d = cos(VX_PARALLEL_DEG * M_PI/180);
	for (int i=0; i<3; i++)
	{
		if (gmAbs(nml[i]) > d)
		{
			double e = (max[i] - min[i]) * VX_OVERLAP_SCALE;
			if (gmFuzEQ(v[i], min[i], e))
				*face = 2 * i;
			else if (gmFuzEQ(v[i], max[i], e))
				*face = 2 * i + 1;
			else
				continue;

			if (nml[i] < 0)
				*dir = -1;
			else
				*dir = 1;
			return 1;
		}
	}
	return 0;
}

void vx_modify_triangle(Vec3f& p0, Vec3f& p1, Vec3f& p2, const float pitch[3], const Vec3f& nml)
{
	Vec3f t(0, 0, 0);
	int axis = gmMinAxis(pitch[0], pitch[1], pitch[2]);
	double e = pitch[axis] * VX_OVERLAP_SCALE;
	t += -nml * 2 * e;

	Vec3f center = (p0 + p1 + p2) / 3;
	Vec3f len((p0-p1).length(), (p1-p2).length(), (p2-p0).length());
	double l = gmMax3(len[0], len[1], len[2]);
	double s = 1 - 4*e/l;

	Mat4x4f M;
	M = Mat4x4f::translate(t);
	M *= Mat4x4f::translate(center);
	M *= Mat4x4f::scale(s);
	M *= Mat4x4f::translate(-center);
	p0 = M.transform(p0);
	p1 = M.transform(p1);
	p2 = M.transform(p2);
}

class VxOctreeStandardTerminator : public VxOctreeTerminator
{
public:
	VxOctreeStandardTerminator(int size, int depth)
	{
		m_size = size;
		m_depth = depth;
	}
	bool shouldTerminate(int size, int depth)
	{
		return (size <= m_size || depth >= m_depth) ? true : false;
	}
private:
	int m_size;
	int m_depth;
};

class VxOctreeGetLeafNode : public VxOctreeTravFunc
{
public:
	VxOctreeGetLeafNode(std::vector<VxOctreeNode*>& leafs)
	: m_leafs(leafs)
	{
	}
	void doLeaf(VxOctreeNode* node)
	{
		m_leafs.push_back(node);
	}
	std::vector<VxOctreeNode*>& m_leafs;
};

class VxOctreeSmoothLevel : public VxOctreeTravFunc
{
public:
	VxOctreeSmoothLevel(VxOctree* octree)
	{
		m_octree = octree;
		m_min[0] = m_min[1] = m_min[2] = 0;
		for (int i=0; i<16; i++)
		{
			m_max[i][0] = m_max[i][1] = m_max[i][2] = int(powf(2, i)) - 1;
		}
		m_update = 0;
		m_max_depth = octree->getMaxDepth();
	}

	void doLeaf(VxOctreeNode* node)
	{
		if (node->getSubdivideFlag()) return;

		int posi[3];
		int depth;
		VxOctreeNode* neighbors[4];
		int nneighbor;
		int i;

		VxOctreePosi node_posi = node->getPosi();
		pdg_posi2xyz(node_posi, &posi[0], &posi[1], &posi[2], &depth);

		if (depth > m_max_depth - 2) return;

		for (int axis=0; axis<3; axis++)
		{
			if (posi[axis] > m_min[axis])
			{
				nneighbor = m_octree->getNeighbors(node, axis, -1, neighbors);
				if (nneighbor == 4)
				{
					for (i=0; i<nneighbor; i++)
					{
						if (!neighbors[i]->isLeaf())
						{
							node->setSubdivideFlag(1);
							m_update = 1;
							return;
						}
					}
				}
			}
			if (posi[axis] < m_max[depth][axis])
			{
				nneighbor = m_octree->getNeighbors(node, axis, +1, neighbors);
				if (nneighbor == 4)
				{
					for (i=0; i<nneighbor; i++)
					{
						if (!neighbors[i]->isLeaf())
						{
							node->setSubdivideFlag(1);
							m_update = 1;
							return;
						}
					}
				}
			}
		}
	}
	VxOctree* m_octree;
	int m_min[3];
	int m_max[16][3];
	int m_update;
	int m_max_depth;
};

class VxOctreeSubdivideNode : public VxOctreeTravFunc
{
public:
	VxOctreeSubdivideNode(VxOctree* octree)
	{
		m_octree = octree;
		m_term = new VxOctreeStandardTerminator(-1, 1);
	}
	~VxOctreeSubdivideNode()
	{
		delete m_term;
	}
	void doLeaf(VxOctreeNode* node)
	{
		if (node->getSubdivideFlag())
		{
			m_octree->subdivideNode(node, m_term);
		}
	}
	VxOctree* m_octree;
	VxOctreeTerminator* m_term;
};

class VxOctreeInfo : public VxOctreeTravFunc
{
public:
	VxOctreeInfo()
	{
		m_max_size = 0;
		m_max_depth = 0;
		m_memory = 0;
		m_nelm = 0;
		m_nleaf = 0;
		m_nboundary = 0;
		m_ninterior = 0;
	}

	void doLeaf(VxOctreeNode* node)
	{
		VxOctreePosi node_posi = node->getPosi();
		int depth = pdg_get_depth(node_posi);
		m_max_depth = gmMax_i(depth, m_max_depth);
		if (node->isBoundary())
		{
			VxOctreeBoundaryData* data = node->getBoundaryData();
			int nelm = data->m_nelm;
			m_nelm += nelm;
			m_max_size = gmMax_i(nelm, m_max_size);
			m_memory += sizeof(VxOctreeNode) + sizeof(VxOctreeBoundaryData) + sizeof(int) * nelm;
			m_nboundary++;
		}
		else
		{
			m_memory += sizeof(VxOctreeNode);
			m_nleaf++;
		}
	}
	void doInterior(VxOctreeNode* node)
	{
		m_memory += sizeof(VxOctreeNode);
		m_ninterior++;
	}

	uint m_max_size;
	uint m_max_depth;
	unsigned long long m_memory;
	uint m_nelm;
	uint m_nleaf;
	uint m_nboundary;
	uint m_ninterior;
};

//=========================================================================

VxOctree::VxOctree(const BBox& bbox, int max_depth)
{
	m_bbox = bbox;
	m_max_size = 0;
	m_max_depth = max_depth;
	m_root = new VxOctreeNode();
	m_nid = 0;
	m_orig = bbox.min;
	m_pitch = bbox.size() / int(powf(2, max_depth));
}

VxOctree::~VxOctree()
{
	delete m_root;
}

void
VxOctree::build(const TriMesh* mesh)
{
	int sz = mesh->nTriangles();
	std::vector<int> indices(sz);
	for (int i=0; i<sz; i++) indices[i] = i;
	build(mesh, indices);
}

void
VxOctree::build(const TriMesh* mesh, std::vector<int>& indices)
{
	tt_info("[vxoctree] building tree... (max_size %d, max_depth %d)\n", m_max_size, m_max_depth);

	m_mesh = mesh;

	TtTime tm;
	tm.start();

	int sz = indices.size();
	int depth = 0;
	VxOctreeStandardTerminator term(m_max_size, m_max_depth);
	subdivide(m_root, m_bbox, indices, depth, &term);

	tm.end();

	VxOctreeInfo info;
	traverse(m_root, &info);

	tt_info("[vxoctree] building tree... done (%.2f sec, max_size %d, max_depth %d, memory %s, nelm %d/%d (x%.2f %s), leaf %d (%s), boundary %d (%s), interior %d (%s))\n",
			tm.getElapsedSec(),
			info.m_max_size,
			info.m_max_depth,
			tt_bytes_to_string(info.m_memory).c_str(),
			info.m_nelm,
			sz,
			float(info.m_nelm)/sz,
			tt_bytes_to_string(info.m_nelm * sizeof(int)).c_str(),
			info.m_nleaf,
			tt_bytes_to_string(info.m_nleaf * sizeof(VxOctreeNode)).c_str(),
			info.m_nboundary,
			tt_bytes_to_string(info.m_nboundary * (sizeof(VxOctreeNode) + sizeof(VxOctreeBoundaryData))).c_str(),
			info.m_ninterior,
			tt_bytes_to_string(info.m_ninterior * sizeof(VxOctreeNode)).c_str()
			);
}

void
VxOctree::subdivide(VxOctreeNode* node, const BBox& bbox, std::vector<int>& indices, int depth, VxOctreeTerminator* term)
{
	int sz = indices.size();

	if (term->shouldTerminate(sz, depth))
	{
		//
		// code for a leaf node
		//
		if (sz == 0)
		{
			VxOctreeLeafData* data = 0;
			node->setLeafData(data);
		}
		else
		{
			VxOctreeBoundaryData* data = new VxOctreeBoundaryData();
			node->setBoundaryData(data);

			int* elms = new int[sz];
			for (int i=0; i<sz; i++) elms[i] = indices[i];
			indices.clear();
			data->m_nelm = sz;
			data->m_elms = elms;

			int gid = m_mesh->getGroupID(elms[0]);
			for (int i=1; i<sz; i++)
			{
				int x = m_mesh->getGroupID(elms[i]);
				gid = std::max(gid, x);
			}
			node->setId(gid + m_nid);
			node->setMark(1);
		}

		return;
	}

	//
	// code for an interior node
	//
	VxOctreeNode* children = new VxOctreeNode[8];
	for (int i=0; i<8; i++)
	{
		VxOctreePosi posi = node->getPosi();
		pdg_push(posi, i);
		children[i].setPosi(posi);
		children[i].setSubdivideFlag(0);
	}
	node->setInterior(children);

	// bboxes for eight children
	BBox c_bbox[8];
	for (int i=0; i<8; i++)
		c_bbox[i] = getOctalBBox(bbox, i);

	// indices for eight children
	std::vector<int> c_indices[8];

	for (int i=0; i<sz; i++)
	{
		int idx = indices[i];
		Vec3f p0 = m_mesh->vertex(idx, 0);
		Vec3f p1 = m_mesh->vertex(idx, 1);
		Vec3f p2 = m_mesh->vertex(idx, 2);
		Vec3f nml = m_mesh->fnormal(idx);

		// modify the position and size of a triangle for robust intersection computation
		if (vx_check_point_on_grid(m_orig, m_pitch, p0)
		 || vx_check_point_on_grid(m_orig, m_pitch, p1)
		 || vx_check_point_on_grid(m_orig, m_pitch, p2))
		{
			vx_modify_triangle(p0, p1, p2, m_pitch, nml);
		}

		for (int j=0; j<8; j++)
		{
			Vec3f bboxcenter = c_bbox[j].center();
			Vec3f bboxhalfsize = .5 * c_bbox[j].size();
			if (f_itsect_box_tri(bboxcenter, bboxhalfsize, p0, p1, p2))
			{
				c_indices[j].push_back(idx);
			}
		}
	}

	indices.clear();

	for (int i=0; i<8; i++)
	{
		subdivide(node->childNode(i), c_bbox[i], c_indices[i], depth + 1, term);
	}
}

void
VxOctree::subdivideNode(VxOctreeNode* node, VxOctreeTerminator* term)
{
	if (!node->isLeaf()) return;

	std::vector<int> indices;
	BBox bbox = getBBox(node);
	if (node->isBoundary())
	{
		VxOctreeBoundaryData* data = node->getBoundaryData();
		int sz = data->m_nelm;
		int* elms = data->m_elms;
		indices.resize(sz);
		for (int i=0; i<sz; i++) indices[i] = elms[i];
	}
	node->deleteData();

	subdivide(node, bbox, indices, 0, term);
}

BBox
VxOctree::getBBox(const VxOctreeNode* node) const
{
	int x, y, z, depth;
	VxOctreePosi node_posi = node->getPosi();
	pdg_posi2xyz(node_posi, &x, &y, &z, &depth);
	Vec3f length = m_bbox.size();
	length /= powf(2, depth);
	Vec3f org = m_bbox.min + length * Vec3f(x, y, z);
	return BBox(org, org + length);
}

VxOctreeNode*
VxOctree::getNode(VxOctreePedigree pedig) const
{
	int oct[16];
	int i;
	int depth;

	depth = pdg_get_depth(pedig);
	for (i=0; i<depth; i++)
	{
		oct[i] = pdg_top(pedig);
		pdg_pop(pedig);
	}

	VxOctreeNode* node = m_root;
	for (i=depth-1; i>=0; i--)
	{
		if (node->isLeaf()) break;
		node = node->childNode(oct[i]);
	}
	return node;
}

VxOctreeNode*
VxOctree::getNeighbor(const VxOctreeNode* node, int axis, int d) const
{
#if 0
	int x, y, z, depth;
	int dd[3];
	VxOctreePedigree pedig2;

	dd[0] = dd[1] = dd[2] = 0;
	dd[axis] += d;
	VxOctreePosi node_posi = node->getPosi();
	pdg_posi2xyz(node_posi, &x, &y, &z, &depth);
	pdg_xyz2pedigree(&pedig2, x + dd[0], y + dd[1], z + dd[2], depth);
	return getNode(pedig2);
#else
	VxOctreePosi node_posi = node->getPosi();
	int oct = pdg_top(node_posi);
	int bro = getBrotherOctalId(oct, axis, d);
	if (bro >= 0)
	{
		return (VxOctreeNode*)(node - oct + bro);
	}

	int x, y, z, depth;
	int dd[3];
	VxOctreePedigree pedig2;

	dd[0] = dd[1] = dd[2] = 0;
	dd[axis] += d;
	pdg_posi2xyz(node_posi, &x, &y, &z, &depth);
	pdg_xyz2pedigree(&pedig2, x + dd[0], y + dd[1], z + dd[2], depth);
	return getNode(pedig2);
#endif
}

int
VxOctree::getNeighbors(const VxOctreeNode* node, int axis, int d, VxOctreeNode* neighbors[4]) const
{
	int i;
	int nneighbor = 0;
	VxOctreeNode* node2 = getNeighbor(node, axis, d);
	// axis ordered cells
	static int axis_cell[3][8] = {
		{0,2,4,6,1,3,5,7},
		{0,1,4,5,2,3,6,7},
		{0,1,2,3,4,5,6,7},
	};
	if (node2->isLeaf())
	{
		neighbors[0] = node2;
		nneighbor = 1;
	}
	else
	{
		if (d == 1)
		{
			for (i=0; i<4; i++)
			{
				neighbors[i] = node2->childNode(axis_cell[axis][i]);
			}
		}
		else if (d == -1)
		{
			for (i=0; i<4; i++)
			{
				neighbors[i] = node2->childNode(axis_cell[axis][i+4]);
			}
		}
		nneighbor = 4;
	}
	return nneighbor;
}

void
VxOctree::traverse(VxOctreeNode* node, VxOctreeTravFunc* func) const
{
	if (node->isLeaf())
	{
		func->doLeaf(node);
		return;
	}

	func->doInterior(node);
	for (int i=0; i<8; i++)
	{
		traverse(node->childNode(i), func);
	}
}

inline void pushLeafNode(std::queue<VxOctreeNode*>& q, VxOctreeNode* node, int id)
{
	node->setMark(1);
	node->setId(id);
	q.push(node);
}

int
VxOctree::fill(VxOctreeNode* node, int id)
{
	if (!node->isLeaf()) return 0;

	std::queue<VxOctreeNode*> q;
	int posi[3];
	int depth;
	int min[3];
	int max[16][3];
	int i;
	VxOctreeNode* neighbors[4];
	int nneighbor;
	int axis;
	int count = 0;

	min[0] = min[1] = min[2] = 0;
	for (i=0; i<16; i++)
	{
		max[i][0] = max[i][1] = max[i][2] = int(powf(2, i)) - 1;
	}

	if (!node->getMark())
	{
		pushLeafNode(q, node, id);
		count++;
	}
	while (q.size() > 0)
	{
		node = q.front();
		q.pop();
		VxOctreePosi node_posi = node->getPosi();
		pdg_posi2xyz(node_posi, &posi[0], &posi[1], &posi[2], &depth);
		for (axis=0; axis<3; axis++)
		{
			if (posi[axis] > min[axis])
			{
				nneighbor = getNeighbors(node, axis, -1, neighbors);
				for (i=0; i<nneighbor; i++)
				{
					if (!neighbors[i]->getMark())
					{
						pushLeafNode(q, neighbors[i], id);
						count++;
					}
				}
			}
			if (posi[axis] < max[depth][axis])
			{
				nneighbor = getNeighbors(node, axis, +1, neighbors);
				for (i=0; i<nneighbor; i++)
				{
					if (!neighbors[i]->getMark())
					{
						pushLeafNode(q, neighbors[i], id);
						count++;
					}
				}
			}
		}
	}
	return count;
}

void
VxOctree::getLeafNodes(std::vector<VxOctreeNode*>& leafs)
{
	VxOctreeGetLeafNode func(leafs);
	traverse(m_root, &func);
}

void
VxOctree::smoothLevel()
{
	tt_info("[vxoctree] smoothing level...\n");

	TtTime tm;
	tm.start();

	int count = 0;
	while (1)
	{
		VxOctreeSmoothLevel smooth(this);
		traverse(m_root, &smooth);
		if (smooth.m_update)
		{
			VxOctreeSubdivideNode subdiv(this);
			traverse(m_root, &subdiv);
			count++;
		}
		else
		{
			break;
		}
	}

	tm.end();

	VxOctreeInfo info;
	traverse(m_root, &info);

	tt_info("[vxoctree] smoothing level... done (%.2f sec, %d iterations, max_size %d, max_depth %d, memory %s, nelm %d (%s), leaf %d (%s), boundary %d (%s), interior %d (%s))\n",
			tm.getElapsedSec(),
			count,
			info.m_max_size,
			info.m_max_depth,
			tt_bytes_to_string(info.m_memory).c_str(),
			info.m_nelm,
			tt_bytes_to_string(info.m_nelm * sizeof(int)).c_str(),
			info.m_nleaf,
			tt_bytes_to_string(info.m_nleaf * sizeof(VxOctreeNode)).c_str(),
			info.m_nboundary,
			tt_bytes_to_string(info.m_nboundary * (sizeof(VxOctreeNode) + sizeof(VxOctreeBoundaryData))).c_str(),
			info.m_ninterior,
			tt_bytes_to_string(info.m_ninterior * sizeof(VxOctreeNode)).c_str()
			);
}

//=========================================================================

BBox
getOctalBBox(const BBox& bbox, int oct)
{
	Vec3f org = bbox.min;
	Vec3f half_size = .5 * bbox.size();
	if (oct & 0x1)
		org[0] += half_size[0];
	if (oct & (0x1<<1))
		org[1] += half_size[1];
	if (oct & (0x1<<2))
		org[2] += half_size[2];
	return BBox(org, org + half_size);
}

int
getOctalId(const float posf[3], const float center[3])
{
	int oct = 0;
	if (posf[0] > center[0])
		oct |= 0x1;
	if (posf[1] > center[1])
		oct |= (0x1<<1);
	if (posf[2] > center[2])
		oct |= (0x1<<2);
	return oct;
}

Vec3f
computeAvgNormal(const TriMesh* mesh, const int* elms, int size)
{
	Vec3f nml(0, 0, 0);
	for (int i=0; i<size; i++)
	{
		nml += mesh->fnormal(elms[i]);
	}
	return nml.normalize();
}

