/*
 * Polylib - Polygon Management Library.
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2010-
 *
 */
#include <vector>
#include <algorithm>
#include "Polylib.h"
#include "polygons/Polygons.h"
#include "polygons/Triangle.h"
#include "polygons/TriMesh.h"
#include "polygons/VTree.h"
#include "common/PolylibCommon.h"
#include "common/tt.h"
#include "common/Vec3.h"
#include "common/BBox.h"
#include "file_io/TriMeshIO.h"

using namespace std;

namespace PolylibNS {

#define M_MAX_ELEMENTS 15	/// VTreeのノードが持つ最大要素数

/************************************************************************
 *  
 * TriMeshクラス
 *  @attention 三角形ポリゴン集合を管理するクラス（KD木用に特化したクラス)
 *  
 ***********************************************************************/
// public /////////////////////////////////////////////////////////////////////
TriMesh::TriMesh()
{
	m_vtree = NULL;
	m_tri_list = NULL;
	m_max_elements = M_MAX_ELEMENTS;
}

// public /////////////////////////////////////////////////////////////////////
TriMesh::~TriMesh()
{
	delete m_vtree;
	if (m_tri_list != NULL) {
		vector<PrivateTriangle*>::iterator itr;
		for (itr = m_tri_list->begin(); itr != m_tri_list->end(); itr++) {
			delete *itr;
		}
		m_tri_list->clear();
	}
	delete m_tri_list;
}

// public /////////////////////////////////////////////////////////////////////
void TriMesh::init(const vector<PrivateTriangle*>* trias)
{
	init_tri_list();
	vector<PrivateTriangle*>::const_iterator itr;
	for (itr = trias->begin(); itr != trias->end(); itr++) {
		m_tri_list->push_back(
			new PrivateTriangle((*itr)->get_vertex(),	(*itr)->get_normal(), 
								(*itr)->get_area(),		(*itr)->get_id())
		);
	}
}

// public /////////////////////////////////////////////////////////////////////
// std::sort用ファンクタ
struct PrivTriaLess{
	bool operator()( const PrivateTriangle *l, const PrivateTriangle *r ) const
	{
		return l->get_id() < r->get_id();
	}
};
// std::equal用ファンクタ
struct PrivTriaEqual{
	bool operator()( const PrivateTriangle *l, const PrivateTriangle *r ) const
	{
		return l->get_id() == r->get_id();
	}
};
void
TriMesh::add(
	const vector<PrivateTriangle*> *trias
)
{
#ifdef DEBUG
	PL_DBGOSH << "TriMesh::add_triangles() in." << endl;
#endif
	unsigned int i;

	if (m_tri_list == NULL) {
		m_tri_list = new vector<PrivateTriangle*>;
	}

	// ひとまず全部追加
	for( i=0; i<trias->size(); i++ ) {
		m_tri_list->push_back( new PrivateTriangle(*(trias->at(i))) );
	}

	// 三角形リストをID順にソート
	std::sort( m_tri_list->begin(), m_tri_list->end(), PrivTriaLess() );

	// ID重複ぶんを削除
	m_tri_list->erase(
		std::unique(m_tri_list->begin(), m_tri_list->end(), PrivTriaEqual()),
		m_tri_list->end());
}

// public /////////////////////////////////////////////////////////////////////
POLYLIB_STAT TriMesh::import(const map<string, string> fmap)
{
	init_tri_list();
	return TriMeshIO::load(m_tri_list, fmap);
}

// public /////////////////////////////////////////////////////////////////////
POLYLIB_STAT TriMesh::build()
{
	BBox bbox;
	vector<PrivateTriangle*>::iterator itr;

	/// TriMeshクラスに含まれる全三角形ポリゴンを外包するBoundingBoxを計算
	bbox.init();
	for (itr = m_tri_list->begin(); itr != m_tri_list->end(); itr++) {
		const Vec3f* vtx_arr = (*itr)->get_vertex();
		for (int i = 0; i < 3; i++) {
			bbox.add(vtx_arr[i]);
		}
	}
	m_bbox = bbox;

#ifdef DEBUG
	Vec3f min = m_bbox.getPoint(0);
	Vec3f max = m_bbox.getPoint(7);
	PL_DBGOSH << "TriMesh::build:min=(" <<min<< "),max=(" <<max<< ")" << endl;
#endif

	// 木構造作成
	if (m_vtree != NULL) delete m_vtree;
	m_vtree = new VTree(m_max_elements, m_bbox, m_tri_list);
	return PLSTAT_OK;
}

// public /////////////////////////////////////////////////////////////////////
int TriMesh::triangles_num() {
	if (m_tri_list == NULL)		return 0;
	else						return m_tri_list->size();
}

// public /////////////////////////////////////////////////////////////////////
const vector<PrivateTriangle*> *TriMesh::search(
	BBox	*bbox, 
	bool	every
) const {
#ifdef DEBUG
	Vec3f min = bbox->getPoint(0);
	Vec3f max = bbox->getPoint(7);
	PL_DBGOSH << "TriMesh::search:min=(" <<min<< "),max=(" <<max<< ")" << endl;
#endif

	return m_vtree->search(bbox, every);
}

// public /////////////////////////////////////////////////////////////////////
POLYLIB_STAT TriMesh::search(
	BBox						*bbox, 
	bool						every, 
	vector<PrivateTriangle*>	*tri_list
) const {
	return m_vtree->search(bbox, every, tri_list);
}

// public /////////////////////////////////////////////////////////////////////
const vector<PrivateTriangle*>* TriMesh::linear_search(
	BBox	*q_bbox, 
	bool	every
) const {
	vector<PrivateTriangle*>		   *tri_list = new vector<PrivateTriangle*>;
	vector<PrivateTriangle*>::iterator itr;

	for (itr = m_tri_list->begin(); itr != m_tri_list->end(); itr++) {
		BBox bbox;
		bbox.init();
		const Vec3f* vtx_arr = (*itr)->get_vertex();
		for (int i = 0; i < 3; i++) {
			bbox.add(vtx_arr[i]);
		}
		if (every == true) {
			if (q_bbox->contain(vtx_arr[0]) == true && 
				q_bbox->contain(vtx_arr[1]) == true &&
				q_bbox->contain(vtx_arr[2]) == true)
			{
				tri_list->push_back(*itr);
			}
		}
		else {
#ifdef OLD_DEF
			if (bbox.crossed(*q_bbox) == true				||
				bbox.contain(q_bbox->getPoint(0)) == true	||
				q_bbox->crossed(bbox) == true				||
				q_bbox->contain(bbox.getPoint(0)) == true) {
#else
			if (bbox.crossed(*q_bbox) == true) {
#endif
				tri_list->push_back(*itr);
			}
		}
	}
	return tri_list;
}

// public /////////////////////////////////////////////////////////////////////
POLYLIB_STAT TriMesh::linear_search(
	BBox						*q_bbox, 
	bool						every, 
	vector<PrivateTriangle*>	*tri_list
) const {
	if (tri_list == NULL) return PLSTAT_ARGUMENT_NULL;

	vector<PrivateTriangle*>::iterator itr;

	for (itr = m_tri_list->begin(); itr != m_tri_list->end(); itr++) {
		BBox bbox;
		bbox.init();
		const Vec3f* vtx_arr = (*itr)->get_vertex();
		for (int i = 0; i < 3; i++) {
			bbox.add(vtx_arr[i]);
		}
		if (every == true) {
			if (q_bbox->contain(vtx_arr[0]) == true	&&
				q_bbox->contain(vtx_arr[1]) == true	&&
				q_bbox->contain(vtx_arr[2]) == true) {
				tri_list->push_back(*itr);
#ifdef DEBUG
			PL_DBGOSH << "TriMesh::linear_search:IN TRUE" << endl;
			PL_DBGOSH << "     vertex 0:" << vtx_arr[0] << endl;
			PL_DBGOSH << "     vertex 1:" << vtx_arr[1] << endl;
			PL_DBGOSH << "     vertex 2:" << vtx_arr[2] << endl;
#endif
			}
		}
		else {
#ifdef OLD_DEF
			if (bbox.crossed(*q_bbox) == true				||
				q_bbox->crossed(bbox) == true				||
				bbox.contain(q_bbox->getPoint(0)) == true	||
				q_bbox->contain(bbox.getPoint(0)) == true) {
#else
			if (bbox.crossed(*q_bbox) == true) {
#endif
				tri_list->push_back(*itr);
#ifdef DEBUG
				PL_DBGOSH << "TriMesh::linear_search:IN FALSE" << endl;
#endif
			}
		}
#ifdef DEBUG
		for (int i=0; i<8; i++) {
			PL_DBGOSH << "TriMesh::linear_search:q_box[" << i << "]:" 
					  << q_bbox->getPoint(i) << endl;
		}
	    PL_DBGOSH << "TriMesh::linear_searc:" << " id:" << (*itr)->get_id()
				  << ",v(" << vtx_arr << ")" << endl;
#endif
	}
	return PLSTAT_OK;
}

// private ////////////////////////////////////////////////////////////////////
void TriMesh::init_tri_list()
{
	if (m_tri_list == NULL) {
		m_tri_list = new vector<PrivateTriangle*>;
	}
	else {
		vector<PrivateTriangle*>::iterator itr;
		for (itr = m_tri_list->begin(); itr != m_tri_list->end(); itr++) {
			delete *itr;
		}
		m_tri_list->clear();
	}
}

} //namespace PolylibNS
