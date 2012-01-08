#include <tt/util/message.h>
#include <tt/shapes/mesh.h>
#include "halfedge.h"

//=========================================================================
//
// class HE_Edge
//
//=========================================================================

HE_Edge::HE_Edge(int eid, HE_Vertex* vtx, HE_Face* face)
{
	m_eid = eid;
	m_vtx = vtx;
	m_face = face;
	m_pair = 0;
	m_next = 0;
	m_rdir = 0;
	m_nonmanifold = 0;

	vtx->m_edge = this;
	face->m_edge = this;
}

HE_Edge::~HE_Edge()
{
}

// check whether an edge 'a' and an edge 'b' are the pair
// RETURN
//    1: a pair, the front faces are the same directions
//   -1: a pair, the front faces are different directions
//    0: not a pair
int
isPair(HE_Edge* a, HE_Edge* b)
{
	if (a->m_face != b->m_face)
	{
		// the pair
		if (a->m_vtx == b->m_next->m_vtx && a->m_next->m_vtx == b->m_vtx)
		{
			return 1;
		}

		// the pair, but a different order
		if (a->m_vtx == b->m_vtx && a->m_next->m_vtx == b->m_next->m_vtx)
		{
			return -1;
		}
	}

	// not the pair
	return 0;
}

//=========================================================================

class UniformFrontFace
{
	public:
		UniformFrontFace(HalfEdge* o) : m_halfedge(o) {}

		int operator ()(HE_Edge* edge, HE_Edge* pair)
		{
			if (pair->m_face->m_flag == 0)
			{
				if (edge->m_rdir == -1)
				{
					m_halfedge->reverseFrontFace(pair->m_face->m_tid);
				}
				return 1;
			}
			return 0;
		}
	private:
		HalfEdge* m_halfedge;
};

//=========================================================================

class QueryCluster
{
	public:
		QueryCluster(HalfEdge* o, float angle) : m_halfedge(o)
		{
			m_min_dot = cos(gmRadians(angle));
		}

		int operator ()(HE_Edge* edge, HE_Edge* pair)
		{
			if (pair->m_face->m_flag == 0)
			{
				Vec3f nml = edge->m_face->m_nml;
				Vec3f pairnml = pair->m_face->m_nml;
				float d = dot(nml, pairnml);
				if (d >= m_min_dot)
				{
					return 1;
				}
			}
			return 0;
		}
	private:
		HalfEdge* m_halfedge;
		float m_min_dot;
};

//=========================================================================
//
// class HalfEdge
//
//=========================================================================

HalfEdge::HalfEdge(TriMesh* o)
{
	m_mesh = o;
}

HalfEdge::~HalfEdge()
{
	int sz;

	sz = m_vdata.size();
	for (int i=0; i<sz; i++)
	{
		delete m_vdata[i];
	}

	sz = m_fdata.size();
	for (int i=0; i<sz; i++)
	{
		delete m_fdata[i];
	}

	int ntri = m_triedge.size();
	for (int i=0; i<ntri; i++)
	{
		for (int j=0; j<3; j++)
		{
			delete m_triedge[i][j];
		}
	}
}

void
HalfEdge::clearFaceFlag()
{
	int ntri = m_triedge.size();
	for (int i=0; i<ntri; i++)
	{
		m_fdata[i]->m_flag = 0;
	}
}

void
HalfEdge::create()
{
	int i, j;
	int ntri = m_mesh->nTriangles();
	int nvtx = m_mesh->nVertices();

	m_triedge.resize(ntri);
	m_fdata.resize(ntri);
	m_vdata.resize(nvtx);

	// create vertices
	for (int i=0; i<nvtx; i++)
	{
		HE_Vertex* vtx = new HE_Vertex();
		vtx->m_vid = i;
		m_vdata[i] = vtx;
	}

	// create faces
	for (int i=0; i<ntri; i++)
	{
		HE_Face* face = new HE_Face();
		face->m_tid = i;
		face->m_flag = 0;
		m_fdata[i] = face;
	}

	// create edges
	int eid = 0;
	for (i=0; i<ntri; i++)
	{
		m_triedge[i].resize(3);
		for (j=0; j<3; j++)
		{
			int vid = m_mesh->m_trivid[i][j];
			int tid = i;
			HE_Vertex* vtx = m_vdata[vid];
			HE_Face* face = m_fdata[tid];
			HE_Edge* edge = new HE_Edge(eid, vtx, face);
			m_triedge[i][j] = edge;
			eid++;
		}
	}

	// set the next edge
	for (i=0; i<ntri; i++)
	{
		for (j=0; j<3; j++)
		{
			HE_Edge* edge = m_triedge[i][j];
			int x = (j == 2) ? 0 : j + 1;
			edge->m_next = m_triedge[i][x];
		}
	}

	// set the pair edge
	for (i=0; i<ntri; i++)
	{
		for (j=0; j<3; j++)
		{
			HE_Edge* edge = m_triedge[i][j];
			if (edge->m_pair == 0)
			{
				std::vector<int>& tid_list0 = m_mesh->m_v_tidlist[edge->m_next->m_vtx->m_vid];
				std::vector<int> tid_list;

				int npair = findPairTri(edge->m_vtx->m_vid, i, tid_list0, tid_list);
				if (npair == 0)
				{
					edge->m_nonmanifold = 2;
				}
				else if(npair > 1)
				{
					edge->m_nonmanifold = 1;
				}
				if (npair == 1)
				{
					int pairtid = tid_list[0];
					for (int k=0; k<3; k++)
					{
						HE_Edge* pair = m_triedge[pairtid][k];
						int r = isPair(edge, pair);
						if (r != 0)
						{
							edge->m_pair = pair;
							pair->m_pair = edge;
							edge->m_rdir = r;
							pair->m_rdir = r;
							break;
						}
					}
				}
			}
		}
	}

	// uniform front faces
	UniformFrontFace query(this);
	for (i=0; i<ntri; i++)
	{
		if (m_fdata[i]->m_flag == 0)
		{
			// create a new cluster
			std::vector<int> tid_list;
			traverse(tid_list, i, query);
			m_cluster.push_back(tid_list);
		}
	}

	// sort triedge order
	for (i=0; i<ntri; i++)
	{
		HE_Edge* edge = m_triedge[i][0];
		for (j=0; j<3; j++)
		{
			m_triedge[i][j] = edge;
			edge = edge->m_next;
		}
	}

	// compute face normals
	for (i=0; i<ntri; i++)
	{
		m_fdata[i]->m_nml = m_mesh->fnormal(i);
	}
}

// find the triangle 
int
HalfEdge::findPairTri(int vid0, int tid0, const std::vector<int>& tid_list0, std::vector<int>& tid_list)
{
	std::vector<int>::const_iterator itr = tid_list0.begin();

	for (; itr != tid_list0.end(); itr++)
	{
		int pairtid = *itr;
		if (pairtid != tid0)
		{
			for (int j=0; j<3; j++)
			{
				if (m_mesh->m_trivid[pairtid][j] == vid0)
				{
					tid_list.push_back(pairtid);
					break;
				}
			}
		}
	}

	return tid_list.size();
}

// reverse the front face
// change edge->m_vtx and edge->m_next
// keep triedge array (i.e. keep the pair)
void
HalfEdge::reverseFrontFace(int tid)
{
	std::vector<HE_Edge*>& triedge = m_triedge[tid];
	std::vector<HE_Vertex*> trivtx(3);
	for (int i=0; i<3; i++)
	{
		trivtx[i] = triedge[i]->m_vtx;
	}
	for (int i=0; i<3; i++)
	{
		int prev = (i == 0) ? 2 : i - 1;
		int next = (i == 2) ? 0 : i + 1;
		triedge[i]->m_next = triedge[prev];
		triedge[i]->m_vtx = trivtx[next];
		if (triedge[i]->m_pair)
		{
			triedge[i]->m_rdir *= -1;
			triedge[i]->m_pair->m_rdir *= -1;
		}
	}
}

void
HalfEdge::selectCluster(std::vector<int>& tid_list, int tid, float angle)
{
	clearFaceFlag();
	QueryCluster query(this, angle);
	traverse(tid_list, tid, query);
}

#if 0
void
HalfEdge::traverseCluster(std::vector<int>& tid_list, int tid, float min_dot)
{
	HE_Edge* edge;
	HE_Face* pairface;

	m_fdata[tid]->m_flag = 1;
	tid_list.push_back(tid);
	Vec3f nml = m_fdata[tid]->m_nml;

	for (int i=0; i<3; i++)
	{
		edge = m_triedge[tid][i];
		if (edge->m_pair == 0)
		{
			continue;
		}
		pairface = edge->m_pair->m_face;
		if (pairface->m_flag == 0)
		{
			Vec3f pairnml = pairface->m_nml;
			float d = dot(nml, pairnml);
			if (gmFuzEQ(d, 1))
			{
				traverseCluster(tid_list, pairface->m_tid, min_dot);
			}
			else
			{
				/*
				Vec3f v0 = m_vtree->m_vdata[edge->m_vtx->m_vid]->m_pos;
				Vec3f v1 = m_vtree->m_vdata[edge->m_next->m_vtx->m_vid]->m_pos;
				Vec3f e1 = v0 - v1;
				Vec3f x1 = cross(pairnml, e1).normalize();
				*/
				//if (dot(nml, x1) < 0 && d > min_dot)
				if (d > min_dot)
				{
					traverseCluster(tid_list, pairface->m_tid, min_dot);
				}
				/*
				Vec3f v0 = m_vtree->m_vdata[edge->m_vtx->m_vid]->m_pos;
				Vec3f v1 = m_vtree->m_vdata[edge->m_next->m_vtx->m_vid]->m_pos;
				Vec3f e0 = v1 - v0;
				Vec3f e1 = -e0;
				Vec3f x0 = cross(nml, e0).normalize();
				Vec3f x1 = cross(pairnml, e1).normalize();
				Vec3f y0 = (nml + pairnml).normalize();
				Vec3f y1 = x0 + x1;
				if (dot(y0, y1) < 0 && d > min_dot)
				{
					traverseCluster(tid_list, pairface->m_tid, min_dot);
				}
				*/
			}
		}
	}
}
#endif

