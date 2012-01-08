#include <tt/gm.h>
#include <tt/util/message.h>
#include <tt/gfx/transform.h>
#include <tt/shapes/mesh.h>
#include "vtree.h"
#include "halfedge.h"
#include "mesh_edit.h"

inline int
tt_prev3(int x)
{
	return (x == 0) ? 2 : x - 1;
}

inline int
tt_next3(int x)
{
	return (x == 2) ? 0 : x + 1;
}

//=========================================================================

TriMeshEditor::TriMeshEditor(TriMesh* mesh)
{
	m_mesh = mesh;
}

TriMesh*
TriMeshEditor::getGroup(int gid)
{
	TriMesh* mesh = new TriMesh();

	TriMeshGroup* grp = m_mesh->getGroup(gid);

	int ntri = grp->nTriangles();

	mesh->m_vdata.reallocate(0, 3*ntri);
	mesh->m_trivid.reallocate(0, ntri);

	mesh->m_ndata.reallocate(0, 3*ntri);
	mesh->m_trinid.reallocate(0, ntri);

	for (int i=0; i<ntri; i++)
	{
		for (int j=0; j<3; j++)
		{
			mesh->m_vdata.push_back(m_mesh->vertex(i, j));
			mesh->m_ndata.push_back(m_mesh->vnormal(i, j));
		}
		Vec3i idx(i*3, i*3+1, i*3+2);
		mesh->m_trivid.push_back(idx);
		mesh->m_trinid.push_back(idx);
	}

	return mesh;
}

void
TriMeshEditor::addMesh(const TriMesh* mesh)
{
	int o_nvtx = m_mesh->nVertices();
	int o_nnml = m_mesh->nNormals();
	int o_ntri = m_mesh->nTriangles();

	int nvtx = mesh->nVertices();
	int nnml = mesh->nNormals();
	int ntri = mesh->nTriangles();
	int ngrp = mesh->nGroups();

	if (nvtx > 0)
	{
		m_mesh->m_vdata.resizeConservatively(o_nvtx + nvtx);
		m_mesh->m_trivid.resizeConservatively(o_ntri + ntri);
		for (int i=0; i<nvtx; i++)
		{
			m_mesh->m_vdata[o_nvtx + i] = mesh->m_vdata[i];
		}
		Vec3i offset(o_nvtx, o_nvtx, o_nvtx);
		for (int i=0; i<ntri; i++)
		{
			m_mesh->m_trivid[o_ntri + i] = mesh->m_trivid[i] + offset;
		}
	}
	if (nnml > 0)
	{
		m_mesh->m_ndata.resizeConservatively(o_nnml + nnml);
		m_mesh->m_trinid.resizeConservatively(o_ntri + ntri);
		for (int i=0; i<nnml; i++)
		{
			m_mesh->m_ndata[o_nnml + i] = mesh->m_ndata[i];
		}
		Vec3i offset(o_nnml, o_nnml, o_nnml);
		for (int i=0; i<ntri; i++)
		{
			m_mesh->m_trinid[o_ntri + i] = mesh->m_trinid[i] + offset;
		}
	}
	if (ngrp > 0)
	{
		for (int i=0; i<ngrp; i++)
		{
			TriMeshGroup* grp = mesh->getGroup(i);
			TriMeshGroup* new_grp = new TriMeshGroup();
			*new_grp = *grp;
			int& offset = o_ntri;
			new_grp->m_first += offset;
			new_grp->m_last += offset;
			m_mesh->addGroup(new_grp);
		}
	}
	tm_finalizeMesh(m_mesh);
}

float
TriMeshEditor::findMinEdgeLengthSquared()
{
	int ntri = m_mesh->nTriangles();
	float min_sqdist;

	// init min_sqdist
	const Vec3f& p0 = m_mesh->vertex(0, 0);
	const Vec3f& p1 = m_mesh->vertex(0, 1);
	min_sqdist = distanceSquared(p0, p1);

	for (int i=0; i<ntri; i++)
	{
		for (int j=0; j<3; j++)
		{
			int k = tt_next3(j);
			const Vec3f& p0 = m_mesh->vertex(i, j);
			const Vec3f& p1 = m_mesh->vertex(i, k);
			float sqdist = distanceSquared(p0, p1);
			min_sqdist = gmMin(sqdist, min_sqdist);
		}
	}

	return min_sqdist;
}

void
TriMeshEditor::mergeVertices(float _sqdist)
{
	int i, j;

	float min_sqdist = .5 * findMinEdgeLengthSquared();
	float sqdist = gmMin(_sqdist, min_sqdist);

	tt_info("[mesh] merge vertices... (dist: %g -> %g)\n", sqrt(_sqdist), sqrt(sqdist));

	int orig_size = m_mesh->m_vdata.size();

	VTree* vtree;

	vtree = m_mesh->getVTree();
	if (!vtree)
	{
		m_mesh->computeVTree(sqdist);
		vtree = m_mesh->getVTree();
	}

	int nvtx = vtree->m_vdata.size();
	int ntri = vtree->m_trivtx.size();

	m_mesh->m_vdata.reallocate(nvtx, nvtx);
	m_mesh->m_v_tidlist.resize(nvtx);
	m_mesh->m_trivid.reallocate(ntri, ntri);

	std::list<VElement*>::iterator itr = vtree->m_vdata.begin();
	for (i=0; i<nvtx; i++, itr++)
	{
		VElement* elm = *itr;
		m_mesh->m_vdata[i] = elm->m_pos;

		std::list<int>::iterator itr2 = elm->m_tidlist.begin();
		int sz = elm->m_tidlist.size();
		m_mesh->m_v_tidlist[i].resize(sz);
		for (j=0; j<sz; j++, itr2++)
		{
			m_mesh->m_v_tidlist[i][j] = *itr2;
		}
	}

	for (i=0; i<ntri; i++)
	{
		const std::vector<VElement*>& elms = vtree->m_trivtx[i];
		m_mesh->m_trivid[i] = Vec3i(elms[0]->m_vid, elms[1]->m_vid, elms[2]->m_vid);
	}

	int curr_size = m_mesh->m_vdata.size();

	tt_info("[mesh] merge vertices... done (%d -> %d vertices)\n", orig_size, curr_size);
}

void
TriMeshEditor::uniformFrontFaces()
{
	tt_info("[mesh] uniform front faces...\n");

	if (m_mesh->m_v_tidlist.size() == 0)
	{
		tt_info("[mesh] uniform front faces... fail\n");
	}

	HalfEdge* halfedge;

	halfedge = m_mesh->getHalfEdge();
	if (!halfedge)
	{
		m_mesh->computeHalfedge();
		halfedge = m_mesh->getHalfEdge();
	}

	int ntri = halfedge->m_triedge.size();

	m_mesh->m_trivid.reallocate(0, ntri);
	for (int i=0; i<ntri; i++)
	{
		const std::vector<HE_Edge*>& edges = halfedge->m_triedge[i];
		Vec3i idx = Vec3i(edges[0]->m_vtx->m_vid, edges[1]->m_vtx->m_vid, edges[2]->m_vtx->m_vid);
		m_mesh->m_trivid.push_back(idx);
	}

	tt_info("[mesh] uniform front faces... done\n");
}

void
TriMeshEditor::heuristicFrontFaces()
{
	float z[3];
	float zcenter;
	float zmax;
	int zmax_tid;
	int ngrp = m_mesh->nGroups();

	for (int i=0; i<ngrp; i++)
	{
		TriMeshGroup* grp = m_mesh->getGroup(i);
		int tid_begin = grp->m_first;
		int tid_last = grp->m_last;

		z[0] = m_mesh->vertex(tid_begin, 0)[2];
		z[1] = m_mesh->vertex(tid_begin, 1)[2];
		z[2] = m_mesh->vertex(tid_begin, 2)[2];
		zmax = (z[0] + z[1] + z[2]) / 3;
		zmax_tid = tid_begin;
		tid_begin += 1;

		for (int j=tid_begin; j<=tid_last; j++)
		{
			z[0] = m_mesh->vertex(j, 0)[2];
			z[1] = m_mesh->vertex(j, 1)[2];
			z[2] = m_mesh->vertex(j, 2)[2];
			zcenter = (z[0] + z[1] + z[2]) / 3;
			if (zcenter > zmax)
			{
				zmax = zcenter;
				zmax_tid = j;
			}
		}

		Vec3f nml = m_mesh->fnormal(zmax_tid);
		if (dot(nml, Vec3f(0, 0, 1)) < 0)
		{
			reverseNormals(i);
		}
	}
}

void
TriMeshEditor::computeNormals(float angle)
{
	m_mesh->m_smooth_angle = angle;

	if (angle == 0)
	{
		tm_computeFlatNormals(m_mesh);
	}
	else
	{
		computeSmoothNormals(angle);
	}
}

struct TriId
{
	TriId(int _idx, int _nth)
	{
		idx = _idx;
		nth = _nth;
	}

	int idx;
	int nth;
};

// compute normals for each vertex
// This function must be called after mergeVertices().
void
TriMeshEditor::computeSmoothNormals(float angle)
{
	tt_info("[mesh] smoothing normals... (angle: %g)\n", angle);

	float coplanar = cos(gmRadians(1));
	float min_dot = cos(gmRadians(angle));
	int nvtx = m_mesh->nVertices();
	int ntri = m_mesh->nTriangles();
	int nnml = 0;

	m_mesh->m_ndata.resize(3 * ntri);
	m_mesh->m_trinid.resize(ntri);

	// build an array of a triangle list shared a vertex
	TtArray< std::list<TriId> > tri_list;
	tri_list.resize(nvtx);
	for (int i=0; i<ntri; i++)
	{
		for (int j=0; j<3; j++)
		{
			int vid = m_mesh->m_trivid[i][j];
			tri_list[vid].push_back(TriId(i, j));
		}
	}

	for (int i=0; i<nvtx; i++)
	{
		// check triangles around a shared vertex
		std::list<TriId>& lst = tri_list[i];
		while (lst.size() > 0)
		{
			TriId tid = lst.front();
			lst.pop_front();

			std::list<TriId> lst_merge;
			lst_merge.push_back(tid);

			Vec3f nml0 = m_mesh->fnormal(tid.idx);
			Vec3f nml = nml0;

			// find to be merged triangles
			std::list<TriId>::iterator itr = lst.begin();
			for (; itr != lst.end();)
			{
				Vec3f test_nml = m_mesh->fnormal(itr->idx);
				// check smoothing angle
				if (dot(nml0, test_nml) >= min_dot)
				{
					// don't include the triangle if a coplanar triangle has been already included
					std::list<TriId>::const_iterator itr3 = lst_merge.begin();
					for (; itr3 != lst_merge.end(); itr3++)
					{
						if (dot(m_mesh->fnormal(itr3->idx), test_nml) >= coplanar)
						{
							break;
						}
					}

					if (itr3 == lst_merge.end())
					{
						nml += test_nml;
					}
					lst_merge.push_back(*itr);
					itr = lst.erase(itr);
				}
				else
				{
					itr++;
				}
			}

			nml.normalize();

			// assign a merged normal to shared vertex normals
			std::list<TriId>::const_iterator itr2 = lst_merge.begin();
			for (; itr2 != lst_merge.end(); itr2++)
			{
				m_mesh->vnormal(itr2->idx, itr2->nth) = nml;
			}

			nnml++;
		}
	}

	// set vertex normal indices
	for (int i=0; i<ntri; i++)
	{
		int x = i * 3;
		m_mesh->m_trinid[i] = Vec3i(x, x + 1, x + 2);
	}

	tt_info("[mesh] smoothing normals... done (%d -> %d normals)\n", ntri*3, nnml);
}

void
TriMeshEditor::mergeNormals()
{
	tt_info("[mesh] merge normals...\n");

	int orig_size = m_mesh->m_ndata.size();
	int curr_size = m_mesh->m_ndata.size();

	tt_info("[mesh] merge normals... done (%d -> %d normals)\n", orig_size, curr_size);
}

void
TriMeshEditor::reverseNormals()
{
	int i;

	int nnml = m_mesh->nNormals();
	int ntri = m_mesh->nTriangles();

	for (i=0; i<nnml; i++)
	{
		m_mesh->m_ndata[i] *= -1;
	}

	for (i=0; i<ntri; i++)
	{
		Vec3i x = m_mesh->m_trivid[i];
		m_mesh->m_trivid[i] = Vec3i(x[0], x[2], x[1]);
	}
}

void
TriMeshEditor::reverseNormals(int gid)
{
	TriMeshGroup* grp = m_mesh->getGroup(gid);
	if (grp->hasRandomTID())
	{
		int sz = grp->m_tidlist.size();
		for (int i=0; i<sz; i++)
		{
			int tid = grp->m_tidlist[i];

			m_mesh->vnormal(tid, 0) *= -1;
			m_mesh->vnormal(tid, 1) *= -1;
			m_mesh->vnormal(tid, 2) *= -1;

			Vec3i x = m_mesh->m_trivid[tid];
			m_mesh->m_trivid[tid] = Vec3i(x[0], x[2], x[1]);
		}
	}
	else
	{
		int first = grp->m_first;
		int last  = grp->m_last;
		for (int i=first; i<=last; i++)
		{
			int tid = i;

			m_mesh->vnormal(tid, 0) *= -1;
			m_mesh->vnormal(tid, 1) *= -1;
			m_mesh->vnormal(tid, 2) *= -1;

			Vec3i x = m_mesh->m_trivid[tid];
			m_mesh->m_trivid[tid] = Vec3i(x[0], x[2], x[1]);
		}
	}
}

void
TriMeshEditor::transform(const Mat4x4f& M)
{
	int nvtx = m_mesh->nVertices();
	int nnml = m_mesh->nNormals();

	// transform vertices
	for (int i=0; i<nvtx; i++)
	{
		m_mesh->m_vdata[i] = M.transform(m_mesh->m_vdata[i]);
	}

	// transform normals
	for (int i=0; i<nnml; i++)
	{
		m_mesh->m_ndata[i] = M.transformVector(m_mesh->m_ndata[i]);
	}

	Transform T(M);
	m_mesh->m_bbox = T.transform(m_mesh->m_bbox);
}

void
TriMeshEditor::grouping()
{
	HalfEdge* halfedge = m_mesh->getHalfEdge();

	tt_info("[mesh] grouping...\n");

	if (!halfedge)
	{
		tt_info("[mesh] grouping... fail\n");
		return;
	}

	int nclusters = halfedge->m_cluster.size();

	m_mesh->clearGroups();
	for (int i=0; i<nclusters; i++)
	{
		TriMeshGroup* grp = new TriMeshGroup();
		grp->m_name = "default";
		grp->m_tidlist = halfedge->m_cluster[i];
		m_mesh->addGroup(grp);
	}

	m_mesh->packVID();

	tt_info("[mesh] grouping... done (%d groups)\n", nclusters);
}

void
TriMeshEditor::selectCluster(std::vector<int>& tid_list, int tid, float angle)
{
	HalfEdge* halfedge = m_mesh->getHalfEdge();

	tt_info("[mesh] select cluster...\n");

	if (!halfedge)
	{
		tt_info("[mesh] select cluster... fail\n");
		return;
	}

	halfedge->selectCluster(tid_list, tid, angle);

	tt_info("[mesh] select cluster... done (%d triangles)\n", tid_list.size());
}

void
TriMeshEditor::fix(float sqdist, float angle)
{
	mergeVertices(sqdist);
	uniformFrontFaces();
	grouping();
	heuristicFrontFaces();
	computeNormals(angle);
	mergeNormals();
}

