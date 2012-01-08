#include <tt/shapes/mesh.h>
#include "vtree.h"

//=========================================================================

VNode::VNode()
{
	m_left = 0;
	m_right = 0;
	m_axis = 0;
}

VNode::~VNode()
{
	delete m_left;
	delete m_right;
}

void
VNode::split()
{
	m_left = new VNode();
	m_right = new VNode();
	m_left->m_bbox = m_bbox;
	m_right->m_bbox = m_bbox;
	float x = .5 * (m_bbox.min[m_axis] + m_bbox.max[m_axis]);
	m_left->m_bbox.max[m_axis] = x;
	m_right->m_bbox.min[m_axis] = x;

	std::list<VElement*>::const_iterator itr = m_vlist.begin();
	for (; itr != m_vlist.end(); itr++)
	{
		if ((*itr)->m_pos[m_axis] < x)
		{
			m_left->m_vlist.push_back((*itr));
		}
		else
		{
			m_right->m_vlist.push_back((*itr));
		}
	}
	m_vlist.clear();

	// set the next axis to split a bounding box
	int axis = (m_axis == 2) ? 0 : m_axis + 1;
	m_left->m_axis = axis;
	m_right->m_axis = axis;
}

//=========================================================================

// find the nearest VElement around m_vtx
class VTreeQuery
{
	public:
		Vec3f m_vtx;
		float m_sqdist;

		// return the nearest VElement
		VElement* m_velement;

		// return the VNode having m_vtx
		VNode* m_vnode;

		VTreeQuery()
		{
			m_velement = 0;
			m_vnode = 0;
		}

		void query(VNode* vn)
		{
			std::list<VElement*>::const_iterator itr = vn->m_vlist.begin();
			for (; itr != vn->m_vlist.end(); itr++)
			{
				float sqdist = distanceSquared((*itr)->m_pos, m_vtx);
				if (sqdist <= m_sqdist)
				{
					m_velement = *itr;
					m_sqdist = sqdist;
				}
			}
			if (m_vnode == 0)
			{
				m_vnode = vn;
			}
		}

};

//=========================================================================

VTree::VTree(TriMesh* m)
{
	m_mesh = m;
	m_root = 0;
	m_max_elements = 15;
}

VTree::~VTree()
{
	destroy();
}

void
VTree::destroy()
{
	if (m_root)
	{
		delete m_root;
		m_root = 0;
	}

	std::list<VElement*>::iterator itr = m_vdata.begin();
	for (; itr != m_vdata.end(); itr++)
	{
		delete *itr;
	}
	m_vdata.clear();
	m_trivtx.clear();
}

void
VTree::create(float sqradius)
{
	destroy();

	m_root = new VNode();
	m_root->m_bbox = m_mesh->getBBox();
	m_root->m_axis = 0;

	VTreeQuery query;
	int newvid = 0;
	int ntri = m_mesh->nTriangles();
	m_trivtx.resize(ntri);
	for (int i=0; i<ntri; i++)
	{
		m_trivtx[i].resize(3);
		for (int j=0; j<3; j++)
		{
			// set query function
			query.m_vtx = m_mesh->vertex(i, j);
			query.m_sqdist = sqradius;
			query.m_velement = 0;
			query.m_vnode = 0;

			traverse(m_root, &query);

			VElement* elm = 0;
			if (query.m_velement == 0)
			{
				// the vtx didn't find in the tree
				// add a new vertex
				elm = new VElement();
				elm->m_pos = query.m_vtx;
				elm->m_vid = newvid;
				elm->m_tidlist.push_back(i);

				m_vdata.push_back(elm);
				newvid++;

				VNode* vnode = query.m_vnode;
				vnode->m_vlist.push_back(elm);
				if (vnode->nElements() > m_max_elements)
				{
					vnode->split();
				}
			}
			else
			{
				// we found a shared vertex.
				elm = query.m_velement;
				elm->m_tidlist.push_back(i);
			}
			m_trivtx[i][j] = elm;
		}
	}
}

// search VNode including "vtx" under the "vn"
VNode*
VTree::search(VNode* vn, Vec3f vtx) const
{
	if (vn->isLeaf())
	{
		return vn;
	}

	int axis = vn->m_axis;
	float x = vn->m_left->m_bbox.max[axis];
	if (vtx[axis] < x)
	{
		return search(vn->m_left, vtx);
	}
	else
	{
		return search(vn->m_right, vtx);
	}
}

void
VTree::traverse(VNode* vn, VTreeQuery* q) const
{
	if (vn->isLeaf())
	{
		q->query(vn);
		return;
	}

	Vec3f& vtx = q->m_vtx;
	float& sqdist = q->m_sqdist;
	int axis = vn->m_axis;
	float x = vn->m_left->m_bbox.max[axis];
	if (vtx[axis] < x)
	{
		traverse(vn->m_left, q);
		float d = x - vtx[axis];
		if (d*d < sqdist)
		{
			traverse(vn->m_right, q);
		}
	}
	else
	{
		traverse(vn->m_right, q);
		float d = vtx[axis] - x;
		if (d*d < sqdist)
		{
			traverse(vn->m_left, q);
		}
	}
}

