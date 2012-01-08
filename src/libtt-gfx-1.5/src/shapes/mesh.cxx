#include <vector>
#include <list>
#include <algorithm>
#include <tt/util/time.h>
#include <tt/util/filename.h>
#include <tt/util/message.h>
#include "mesh.h"
#include <tt/shapes/mesh/mesh_edit.h>
#include <tt/shapes/mesh/mesh_io.h>
#include <tt/shapes/mesh/vtree.h>
#include <tt/shapes/mesh/halfedge.h>

using namespace std;

//=========================================================================
// class TriMesh
//=========================================================================

TriMesh::TriMesh()
{
	m_axis = 0;
	m_size = 0;
	m_scale = 0;
	m_smooth_angle = 0;

	m_vtree = 0;
	m_halfedge = 0;
}

TriMesh::~TriMesh()
{
	delete m_vtree;
	delete m_halfedge;
}

BBox
TriMesh::getTriangleBBox(int tid) const
{
	BBox bbox;
	bbox.add(this->vertex(tid, 0));
	bbox.add(this->vertex(tid, 1));
	bbox.add(this->vertex(tid, 2));
	return bbox;
}

int
TriMesh::getGroupID(int tid) const
{
	int sz = m_groups.size();
	for (int i=0; i<sz; i++)
	{
		TriMeshGroup* grp = m_groups.getData(i);
		if (tid <= grp->m_last) return i;
	}
	return -1;
}

const std::string&
TriMesh::getMaterialName(int tid) const
{
	static std::string none = "";
	int gid = getGroupID(tid);
	TriMeshGroup* grp = m_groups.getData(gid);
	if (grp)
		return grp->m_name;
	else
		return none;
}

void
TriMesh::computeBBox()
{
	BBox bbox;

	bbox.min = m_vdata[0];
	bbox.max = m_vdata[0];

	int nvtx = this->nVertices();

	for (int i=1; i<nvtx; i++)
	{
		bbox.add(m_vdata[i]);
	}

	m_bbox = bbox;
}

void
TriMesh::computeVTree(float sqdist)
{
	delete m_vtree;
	m_vtree = new VTree(this);
	m_vtree->create(sqdist);
}

void
TriMesh::computeHalfedge()
{
	delete m_halfedge;
	m_halfedge = new HalfEdge(this);
	m_halfedge->create();
}

void
TriMesh::packVID()
{
	int ntri = nTriangles();
	int ngrp = nGroups();
	int seq = 0;

	TtArray<Vec3i> tmparray(ntri);
	std::copy(m_trivid.begin(), m_trivid.end(), tmparray.begin());

	for (int i=0; i<ngrp; i++)
	{
		TriMeshGroup* grp = getGroup(i);
		int ntid = grp->m_tidlist.size();

		// this is an option which preserves the original order of face indices
		std::sort(grp->m_tidlist.begin(), grp->m_tidlist.end());

		grp->m_first = seq;
		for (int j=0; j<ntid; j++)
		{
			int tid = grp->m_tidlist[j];
			m_trivid[seq] = tmparray[tid];
			seq++;
		}
		grp->m_last = seq - 1;
		grp->m_tidlist.clear();
	}
}

void
TriMesh::unpackVID()
{
	int ngrp = nGroups();

	for (int i=0; i<ngrp; i++)
	{
		TriMeshGroup* grp = getGroup(i);
		int first = grp->m_first;
		int last  = grp->m_last;
		for (int j=first; j<=last; j++)
		{
			grp->m_tidlist.push_back(j);
		}
		grp->m_first = 0;
		grp->m_last = 0;
	}
}

void
TriMesh::newGroup(const std::string& mtl_name, std::vector<int>& tid_list)
{
	unpackVID();

	std::sort(tid_list.begin(), tid_list.end());
	cout << tid_list.size() << endl;

	// subtract triangle ids from existing groups
	int ngroups = nGroups();
	for (int i=0; i<ngroups; i++)
	{
		TriMeshGroup* grp = getGroup(i);
		std::vector<int> tmp_list;
		std::set_difference(grp->m_tidlist.begin(), grp->m_tidlist.end(), tid_list.begin(), tid_list.end(), back_inserter(tmp_list));
		cout << i << ": " << grp->m_tidlist.size() << " " << tmp_list.size() << endl;
		grp->m_tidlist = tmp_list;
	}

	// add a new group
	TriMeshGroup* grp = new TriMeshGroup();
	grp->m_name = mtl_name;
	grp->m_tidlist = tid_list;
	m_groups.addData(grp);

	tid_list.clear();

	packVID();
	tm_computeFlatNormals(this);
}

void
TriMesh::renameGroup(int gid, const std::string& mtl_name)
{
	TriMeshGroup* grp = getGroup(gid);
	if (grp)
		grp->m_name = mtl_name;
}

int
TriMesh::load(const TtFileName& fname, std::string fmt)
{
	TriMeshIO io(this);
	return io.load(fname, fmt);
}

int
TriMesh::save(const TtFileName& fname, std::string fmt) const
{
	TriMeshIO io(this);
	return io.save(fname, fmt);
}

//=========================================================================

void
tm_computeFlatNormals(TriMesh* mesh)
{
	tt_info("[mesh] flat normals...\n");

	int i, j;

	int ntri = mesh->nTriangles();

	mesh->m_ndata.resize(3 * ntri);
	mesh->m_trinid.resize(ntri);

	// set vertex normal
	for (i=0; i<ntri; i++)
	{
		Vec3f nml = mesh->fnormal(i);
		for (j=0; j<3; j++)
		{
			mesh->m_ndata[i * 3 + j] = nml;
		}
	}

	// set vertex normal indices
	for (i=0; i<ntri; i++)
	{
		int x = i * 3;
		mesh->m_trinid[i] = Vec3i(x, x + 1, x + 2);
	}

	tt_info("[mesh] flat normals... done\n");
}

void
tm_finalizeMesh(TriMesh* mesh)
{
	tt_info("[mesh] finalize...\n");
	mesh->computeBBox();

	if (mesh->nNormals() == 0)
	{
		TriMeshEditor me(mesh);

		float angle = mesh->m_smooth_angle;
		if (angle == 0)
		{
			me.computeNormals(angle);
		}
		else
		{
			me.mergeVertices(0);
			me.computeNormals(angle);
			me.mergeNormals();
		}
	}

	if (mesh->nGroups() == 0)
	{
		tt_info("[mesh] add a default group...\n");
		TriMeshGroup* grp = new TriMeshGroup();
		grp->m_name = "default";
		grp->m_first = 0;
		grp->m_last = mesh->nTriangles() - 1;
		mesh->addGroup(grp);
	}

	tt_info("[mesh] finalize... done\n");
}

