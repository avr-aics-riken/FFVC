#ifndef vxoctree_h
#define vxoctree_h

#include <vector>
#include <tt/gfx/bbox.h>
#include <tt/util/bits.h>
#include "pedig2.h"

class TriMesh;
class VxOctreeNode;

//=========================================================================

class VxOctreeTerminator
{
public:
	virtual ~VxOctreeTerminator() {}
	virtual bool shouldTerminate(int size, int depth) = 0;
};

//=========================================================================

class VxOctreeTravFunc
{
public:
	virtual ~VxOctreeTravFunc() {}
	virtual void doLeaf(VxOctreeNode* node) {}
	virtual void doInterior(VxOctreeNode* node) {}
};

//=========================================================================

class VxOctreeLeafData
{
public:
	VxOctreeLeafData() {}
};

class VxOctreeBoundaryData : public VxOctreeLeafData
{
public:
	VxOctreeBoundaryData() : m_elms(0) {}
	~VxOctreeBoundaryData()
	{
		if (m_elms) delete [] m_elms;
	}
	int m_nelm;
	int* m_elms;
};

//=========================================================================

class VxOctreeNode
{
private:
	enum
	{
		LEAF_BIT = 4,
		BOUNDARY_BIT,
		SUBDIV_BIT,
		MARK_BIT,
		ID_BIT,
	};
public:
	VxOctreeNode() : m_children(0) {}
	~VxOctreeNode()
	{
		deleteData();
	}
	void deleteData()
	{
		if (getLeafFlag())
		{
			if (getBoundaryFlag())
			{
				if (m_bdata) delete m_bdata;
			}
			else
			{
				if (m_data) delete m_data;
			}
		}
		else
		{
			if (m_children) delete [] m_children;
		}
		m_children = 0;
	}

	void setPosi(const VxOctreePosi& posi) { m_posi = posi; }
	VxOctreePosi getPosi() const { return m_posi; }

	void setLeafFlag(int b) { tt_bits_set(m_posi.w, LEAF_BIT, b); }
	bool getLeafFlag() const { return tt_bits_get(m_posi.w, LEAF_BIT); }
	void setBoundaryFlag(int b) { tt_bits_set(m_posi.w, BOUNDARY_BIT, b); }
	bool getBoundaryFlag() const { return tt_bits_get(m_posi.w, BOUNDARY_BIT); }
	void setSubdivideFlag(int b) { tt_bits_set(m_posi.w, SUBDIV_BIT, b); }
	bool getSubdivideFlag() const { return tt_bits_get(m_posi.w, SUBDIV_BIT); }

	bool isLeaf() const { return getLeafFlag(); }
	bool isBoundary() const { return getBoundaryFlag(); }

	// methods for an interior node
	void setInterior(VxOctreeNode* children) {
		setLeafFlag(0);
		setBoundaryFlag(0);
		m_children = children;
	}
	VxOctreeNode* childNode(int i) { return &(m_children[i]); }
	const VxOctreeNode* childNode(int i) const { return &(m_children[i]); }

	// methods for a leaf node
	void setLeafData(VxOctreeLeafData* data) {
		setLeafFlag(1);
		setBoundaryFlag(0);
		m_data = data;
	}
	VxOctreeLeafData* getLeafData() const { return m_data; }
	void setId(unsigned char id) { tt_bits_setn<unsigned short>(m_posi.w, ID_BIT, 8, id); }
	unsigned char getId() const { return tt_bits_getn(m_posi.w, ID_BIT, 8); }
	void setMark(int b) { tt_bits_set(m_posi.w, MARK_BIT, b); }
	bool getMark() const { return tt_bits_get(m_posi.w, MARK_BIT); }

	// methods for a boundary leaf node
	void setBoundaryData(VxOctreeBoundaryData* data) {
		setLeafFlag(1);
		setBoundaryFlag(1);
		m_bdata = data;
	}
	VxOctreeBoundaryData* getBoundaryData() const { return m_bdata; }

private:
	VxOctreePosi m_posi;
	union {
		VxOctreeNode* m_children;
		VxOctreeLeafData* m_data;
		VxOctreeBoundaryData* m_bdata;
	};
};

//=========================================================================

class VxOctree
{
public:
	VxOctree(const BBox& bbox, int max_depth);
	~VxOctree();

	const BBox& getBBox() const { return m_bbox; }
	BBox getBBox(const VxOctreeNode* node) const;
	int getMaxDepth() const { return m_max_depth; }
	VxOctreeNode* root() const { return m_root; }
	VxOctreeNode* getNode(VxOctreePedigree pedig) const;
	VxOctreeNode* getNeighbor(const VxOctreeNode* node, int axis, int d) const;
	int getNeighbors(const VxOctreeNode* node, int axis, int d, VxOctreeNode* neighbors[4]) const;

	void build(const TriMesh* mesh);
	void build(const TriMesh* mesh, std::vector<int>& indices);
	void smoothLevel();
	int  fill(VxOctreeNode* node, int id);
	void getLeafNodes(std::vector<VxOctreeNode*>& leafs);

	const TriMesh* getTriMesh() const { return m_mesh; }

	void traverse(VxOctreeNode* node, VxOctreeTravFunc* func) const;
	void subdivideNode(VxOctreeNode* node, VxOctreeTerminator* term);

	void setNumId(int nr) { m_nid = nr; }
	int  getNumId() const { return m_nid; }

private:
	void subdivide(VxOctreeNode* node, const BBox& bbox, std::vector<int>& indices, int depth, VxOctreeTerminator* term);

	const TriMesh* m_mesh;
	BBox m_bbox;
	VxOctreeNode* m_root;
	int m_max_size;
	int m_max_depth;
	int m_nid;
	Vec3f m_orig;
	Vec3f m_pitch;
};

//=========================================================================

BBox getOctalBBox(const BBox& bbox, int oct);
int getOctalId(const float posf[3], const float center[3]);
Vec3f computeAvgNormal(const TriMesh* mesh, const int* elms, int size);

// return values
//   [0, 7]: brother's octal id
//   -1 : no brother
inline int getBrotherOctalId(int oct, int axis, int dir)
{
	// {-x, +x, -y, +y, -z, +z}
	static int bro[8][6] = {
		{-1,  1, -1,  2, -1,  4},
		{ 0, -1, -1,  3, -1,  5},
		{-1,  3,  0, -1, -1,  6},
		{ 2, -1,  1, -1, -1,  7},
		{-1,  5, -1,  6,  0, -1},
		{ 4, -1, -1,  7,  1, -1},
		{-1,  7,  4, -1,  2, -1},
		{ 6, -1,  5, -1,  3, -1},
	};
	dir = dir < 0 ? 0 : dir;
	int x = axis * 2 + dir;
	return bro[oct][x];
}

#endif  // vxoctree_h

