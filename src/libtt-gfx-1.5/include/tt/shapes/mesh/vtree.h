#ifndef vtree_h
#define vtree_h

#include <tt/gm.h>
#include <vector>
#include <list>

class TriMesh;
class VTreeQuery;

//=========================================================================

class VElement
{
public:
	// position
	Vec3f m_pos;

	// new vertex id
	int m_vid;

	// a list of triangle indices shared this vertex
	std::list<int> m_tidlist;
};

//=========================================================================

class VNode
{
public:
	VNode();
	~VNode();

	void split();

	bool isLeaf() const { return (m_left == 0); }
	int nElements() const { return m_vlist.size(); }

	VNode* m_left;
	VNode* m_right;
	BBox m_bbox;
	int m_axis;
	std::list<VElement*> m_vlist;
};

//=========================================================================

class VTree
{
public:
	VTree(TriMesh* m);
	~VTree();

	void destroy();
	void create(float sqradius);

	std::list<VElement*> m_vdata;
	std::vector< std::vector<VElement*> > m_trivtx;

private:
	VNode* search(VNode* root, Vec3f vtx) const;
	void traverse(VNode* root, VTreeQuery* q) const;

	TriMesh* m_mesh;
	VNode* m_root;
	int m_max_elements;
};

#endif  // vtree_h

