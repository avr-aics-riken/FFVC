/*
 * Polylib - Polygon Management Library.
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2010-
 *
 */
#include "common/PolylibCommon.h"
#include "common/Vec3.h"
#include "common/BBox.h"
#include "polygons/Triangle.h"
#include "polygons/Polygons.h"
#include "polygons/TriMesh.h"
#include "polygons/VTree.h"


namespace PolylibNS {

using namespace std;

#ifdef DEBUG_VTREE
static vector<VNode*> m_vnode;
#endif

/************************************************************************
 *  
 * VElementクラス
 *  @attention KD木構造の要素クラス
 *  
 ***********************************************************************/
// public /////////////////////////////////////////////////////////////////////
VElement::VElement(
	PrivateTriangle* tri
) {
	m_tri = tri;
	for(int i=0; i<3; i++){
		m_bbox.add(tri->get_vertex()[i]);
	}
	m_pos = m_bbox.center();
}

/************************************************************************
 *  
 * VNodeクラス
 *  @attention KD木構造のノードクラス
 *  
 ***********************************************************************/
// public /////////////////////////////////////////////////////////////////////
VNode::VNode()
{
	m_left = NULL;
	m_right = NULL;
	m_axis = AXIS_X;
	m_bbox_search.init();
#ifdef USE_DEPTH
	m_depth = 0;
#endif
}

// public /////////////////////////////////////////////////////////////////////
VNode::~VNode()
{
	vector<VElement*>::iterator itr = m_vlist.begin();
	for (; itr != m_vlist.end(); itr++) {
		delete *itr;
	}
	m_vlist.clear();
	if (m_left!=NULL) {delete m_left; m_left=NULL;}
	if (m_right!=NULL){delete m_right; m_right=NULL;}
}

// public /////////////////////////////////////////////////////////////////////
void VNode::split(const int& max_elem)
{
	m_left = new VNode();
	m_right = new VNode();

	BBox left_bbox = m_bbox;
	BBox right_bbox = m_bbox;

	float x = .5 * (m_bbox.min[m_axis] + m_bbox.max[m_axis]);
	left_bbox.max[m_axis] = x;
	right_bbox.min[m_axis] = x;

	m_left->set_bbox(left_bbox);
	m_right->set_bbox(right_bbox);

#ifdef USE_DEPTH
	m_left->m_depth = m_depth+1;
	m_right->m_depth = m_depth+1;
#endif

	vector<VElement*>::const_iterator itr = m_vlist.begin();
	for (; itr != m_vlist.end(); itr++) {
		if ((*itr)->get_pos()[m_axis] < x) {
			m_left->m_vlist.push_back((*itr));
			m_left->set_bbox_search((*itr));
		}
		else {
			m_right->m_vlist.push_back((*itr));
			m_right->set_bbox_search((*itr));
		}
	}
	m_vlist.clear();

	// set the next axis to split a bounding box
	AxisEnum axis;
	if (m_axis == AXIS_Z)			axis = AXIS_X;
	else if (m_axis == AXIS_X)		axis = AXIS_Y;
	else							axis = AXIS_Z;
	m_left->set_axis(axis);
	m_right->set_axis(axis);

	if (m_left->get_elements_num() > max_elem) {
		m_left->split(max_elem);
	}
	if (m_right->get_elements_num() > max_elem) {
		m_right->split(max_elem);
	}
}

#ifdef USE_DEPTH
// public /////////////////////////////////////////////////////////////////////
void dump_depth(int n) {
	if (is_leaf()) {
		for (int i = 0; i < n; i++) cout << " ";
		cout << "depth: " << m_depth;
		cout << " elements: " << get_elements_num() << endl;
	}
	else {
		for (int i = 0; i < n; i++) cout << " ";
		cout << "depth: " << m_depth;
		cout << " left : "	<< " max: " << m_left->get_bbox().max 
							<< " min: " << m_left->get_bbox().min << endl;
		m_left->dump_depth(n+1);

		for (int i = 0; i < n; i++) cout << " ";
		cout << "depth: " << m_depth;
		cout << " right: "	<< " max: " << m_right->get_bbox().max 
							<< " min: " << m_right->get_bbox().min << endl;
		m_right->dump_depth(n+1);
	}
}
#endif

/************************************************************************
 *
 * VTreeクラス
 * @attention リーフを三角形ポリゴンとするKD木クラス
 *
 ***********************************************************************/
// public /////////////////////////////////////////////////////////////////////
VTree::VTree(
	int							max_elem, 
	const BBox					bbox, 
	vector<PrivateTriangle*>	*tri_list
) {
	m_root = NULL;
	create(max_elem, bbox, tri_list);
}

// public /////////////////////////////////////////////////////////////////////
VTree::~VTree()
{
	destroy();
}

// public /////////////////////////////////////////////////////////////////////
void VTree::destroy()
{
	if (m_root) {
		delete m_root;
		m_root = NULL;
	}
}

// public /////////////////////////////////////////////////////////////////////
vector<PrivateTriangle*>* VTree::search(
	BBox	*bbox, 
	bool	every
) const {
#ifdef DEBUG_VTREE
	PL_DBGOSH << "VTree::search1:@@@------------------------@@@" << endl;
	Vec3f min = bbox->getPoint(0);
	Vec3f max = bbox->getPoint(7);
	PL_DBGOSH << "VTree::min(" << min << "),max(" << max << ")" << endl;
#endif

	if (m_root == 0) {
		cerr << "Polylib::vtree::Error" << endl;
		exit(1);
	}
	vector<VElement*> vlist;
	search_recursive(m_root, *bbox, every, &vlist);

	vector<VElement*>::iterator itr=vlist.begin();

#ifdef DEBUG_VTREE
	PL_DBGOSH << "VTree::search_recursive end" << endl;
	for (; itr != vlist.end(); itr++) {
		PL_DBGOSH << "VTree::search:tid=" << (*itr)->get_triangle()->get_id() 
				  << endl;
	}
	itr=vlist.begin();
#endif

	vector<PrivateTriangle*> *tri_list = new vector<PrivateTriangle*>;
	for (; itr != vlist.end(); itr++) {
#ifdef MEMCOPY_TYPE
		PrivateTriangle *tri = new PrivateTriangle(
				(*itr)->getVertex(), (*itr)->get_normal(),
				(*itr)->getArea(), (*itr)->get_id()) ;

#ifdef DEBUG_VTREE
		PL_DBGOSH << "VTree::get_id:" << tri->get_id() << endl;
#endif
		tri_list->push_back(*tri);
		delete(tri);
#endif
		tri_list->push_back((*itr)->get_triangle());
	}
	return tri_list;
}

// public /////////////////////////////////////////////////////////////////////
POLYLIB_STAT VTree::search(
	BBox						*bbox, 
	bool						every, 
	vector<PrivateTriangle*>	*tri_list
) const {
#ifdef DEBUG_VTREE
	PL_DBGOSH << "VTree::search2:@@@------------------------@@@" << endl;
	Vec3f min = bbox->getPoint(0);
	Vec3f max = bbox->getPoint(7);
	PL_DBGOSH << "VTree::min(" << min << "),max(" << max << ")" << endl;
#endif

	if (m_root == 0) {
		cout << "Error" << endl;
		return PLSTAT_ROOT_NODE_NOT_EXIST;
	}
	vector<VElement*> vlist;
	search_recursive(m_root, *bbox, every, &vlist);

	vector<VElement*>::iterator itr=vlist.begin();

#ifdef DEBUG_VTREE
	for (; itr != vlist.end(); itr++) {
		PL_DBGOSH << "VTree::search:tid=" << (*itr)->get_triangle()->get_id() 
				  << endl;
	}
#endif

// 木を検索して返ってきた結果を返り値に変換する。

	itr = vlist.begin();
	for (; itr != vlist.end(); itr++) {
#ifdef MEMCOPY_TYPE
		PrivateTriangle *tri = new PrivateTriangle;

		tri->m_id =  (*itr)->m_id;
		tri->m_area =  (*itr)->m_area;
		for (int j=0; j<3; j++) {
			(tri->m_vertex[j])[AXIS_X] = (((*itr)->m_vertex)[j])[AXIS_X];
			(tri->m_vertex[j])[AXIS_Y] = (((*itr)->m_vertex)[j])[AXIS_Y];
			(tri->m_vertex[j])[AXIS_Z] = (((*itr)->m_vertex)[j])[AXIS_Z];
		}
		tri->m_normal[0] = ((*itr)->m_normal)[0];
		tri->m_normal[1] = ((*itr)->m_normal)[1];
		tri->m_normal[2] = ((*itr)->m_normal)[2];

		tri_list->push_back(tri);
#endif
		tri_list->push_back( (*itr)->get_triangle() );
	}
	return PLSTAT_OK;
}

// public /////////////////////////////////////////////////////////////////////
unsigned int VTree::memory_size() {
	VNode			*vnode;
	unsigned int	node_cnt = 1;		// ノード数
	unsigned int	poly_cnt = 0;		// ポリゴン数
	unsigned int	size;

	if ((vnode = m_root->get_left()) != NULL) {; 
		node_count(vnode, &node_cnt, &poly_cnt);
	}
#ifdef DEBUG_VTREE
	PL_DBGOSH << "VTree::memory_size1():node,poly=" << node_cnt << "," 
			  << poly_cnt << endl;
#endif
		
	if ((vnode = m_root->get_right()) != NULL) {; 
		node_count(vnode, &node_cnt, &poly_cnt);
	}
#ifdef DEBUG_VTREE
	PL_DBGOSH << "VTree::memory_size2():node,poly=" << node_cnt << "," 
			  << poly_cnt << endl;
#endif

	size  = sizeof(VTree);
	size += sizeof(VNode)	 * node_cnt;
	size += sizeof(VElement) * poly_cnt;
	return size;
}

// private ////////////////////////////////////////////////////////////////////
void VTree::traverse(VNode* vn, VElement* elm, VNode** vnode) const
{
// --- ims ---<
	// set bbox for search triangle
	vn->set_bbox_search(elm);
// --- ims --->

	if (vn->is_leaf()) {
		if (*vnode == 0) {
			*vnode = vn;
		}
		return;
	}

	Vec3f vtx = elm->get_pos();
#ifdef SQ_RADIUS
	float& sqdist = q->m_sqdist;
#endif
	AxisEnum axis = vn->get_axis();
	float x = vn->get_left()->get_bbox().max[axis];
	if (vtx[axis] < x) {
		traverse(vn->get_left(), elm, vnode);
#ifdef SQ_RADIUS
		float d = x - vtx[axis];
		if (d*d < sqdist) {
			traverse(vn->get_right(), elm, vnode);
		}
#endif
	}
	else {
		traverse(vn->get_right(), elm, vnode);
#ifdef SQ_RADIUS
		float d = vtx[axis] - x;
		if (d*d < sqdist) {
			traverse(vn->get_left(), elm, vnode);
		}
#endif
	}
}

// private ////////////////////////////////////////////////////////////////////
void VTree::search_recursive(
	VNode				*vn, 
	const BBox			&bbox, 
	bool				every, 
	vector<VElement*>	*vlist
) const {
#ifdef DEBUG_VTREE
try{
	PL_DBGOSH << "VTree::search_recursive:@@@----------------------@@@" << endl;
#endif
	if (vn->is_leaf()) {
		vector<VElement*>::const_iterator itr = vn->get_vlist().begin();
		for (; itr != vn->get_vlist().end(); itr++) {
			// determine between bbox and 3 vertices of each triangle.
			if (every == true) {
				bool iscontain = true;
				const Vec3f *temp = (*itr)->get_triangle()->get_vertex();
				for (int i = 0; i < 3; i++) {
					if (bbox.contain(temp[i]) == false)  {
						iscontain = false;
						break;
					}
				}
				if (iscontain == true) {
					vlist->push_back(*itr);
				}
			}
			else{
				// determine between bbox and bbox crossed
				BBox e_bbox = (*itr)->get_bbox();

				if (e_bbox.crossed(bbox) == true) {
						vlist->push_back(*itr);
				}
			}
		}
#ifdef USE_DEPTH
		PL_DBGOSH << "VTree::search_recursive:depth=" << vn->get_depth() 
				  << ",elem= " << vn->get_vlist().size() << endl;
#endif
#ifdef DEBUG_VTREE
	PL_DBGOSH << "VTree::search_recursive:@@@-return---------------@@@" << endl;
#endif
		return;
	}

#ifdef DEBUG_VTREE
	Vec3f min = bbox.getPoint(0);
	Vec3f max = bbox.getPoint(7);
	PL_DBGOSH << "VTree::min(" << min << "),max(" << max << ")" << endl;
#endif

	BBox lbox = vn->get_left()->get_bbox_search();
	BBox rbox = vn->get_right()->get_bbox_search();

	if (lbox.crossed(bbox) == true) {
#ifdef USE_DEPTH
		PL_DBGOSH << "VTree::search_recursive:left=" << vn->get_depth() << endl;
#endif
		search_recursive(vn->get_left(), bbox, every, vlist);
	}

	if (rbox.crossed(bbox) == true) {
#ifdef USE_DEPTH
	PL_DBGOSH << "VTree::search_recursive:right=" << vn->get_depth() << endl;
#endif
		search_recursive(vn->get_right(), bbox, every, vlist);
	}
#ifdef DEBUG_VTREE
}
catch(char *str) {
	cout << str;
}
#endif
}

// private ////////////////////////////////////////////////////////////////////
#ifdef SQ_RADIUS
POLYLIB_STAT VTree::create(float sqradius) {
#else
POLYLIB_STAT VTree::create(
	int							max_elem, 
	const BBox					bbox, 
	vector<PrivateTriangle*>	*tri_list
) {
#endif
	destroy();

	m_max_elements = max_elem;
	m_root = new VNode();
	m_root->set_bbox(bbox);
	m_root->set_axis(AXIS_X);

	vector<PrivateTriangle*>::iterator itr;
	for (itr = tri_list->begin(); itr != tri_list->end(); itr++) {
		// make a new triangle
		VElement* elm = NULL;
		elm = new VElement(*itr);

		VNode* vnode = NULL;
		traverse(m_root, elm, &vnode);

		// the vtx didn't find in the tree
		// add a new vertex
		if (vnode == NULL) {
			PL_ERROSH << "[ERROR]VTree::create():Can't find appropriate node" 
					  << endl;
			return PLSTAT_NODE_NOT_FIND;
		}

		// find node to add a new triangle
		vnode->set_element(elm);

		// set bbox for search triangle
		vnode->set_bbox_search(elm);
		if (vnode->get_elements_num() > m_max_elements) {
			vnode->split(m_max_elements);
#ifdef DEBUG_VTREE
			m_vnode.push_back(vnode->get_left());
			m_vnode.push_back(vnode->get_right());
#endif
		}
	}

#ifdef USE_DEPTH
	m_root->dump_depth(0);
#endif
	return PLSTAT_OK;
}

// private ////////////////////////////////////////////////////////////////////
void VTree::node_count(
	VNode			*parent, 
	unsigned int	*node_cnt, 
	unsigned int	*tri_cnt
) {
	VNode	*vnode;

#ifdef DEBUG_VTREE
	PL_DBGOSH << "VTree::node_count1():" << *node_cnt << endl;
#endif
	if ((vnode = parent->get_left()) != NULL) {; 
		(*node_cnt)++;
		node_count(vnode, node_cnt, tri_cnt);
	}
	else {
		(*tri_cnt) += parent->get_elements_num();
	}
#ifdef DEBUG_VTREE
	PL_DBGOSH << "VTree::node_count2():" << *node_cnt << endl;
#endif
	if ((vnode = parent->get_right()) != NULL) {; 
		(*node_cnt)++;
		node_count(vnode, node_cnt, tri_cnt);
	}
}

} //namespace PolylibNS
// --- ims --->
