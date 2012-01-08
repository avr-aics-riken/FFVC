#ifndef halfedge_h
#define halfedge_h

#include <tt/gm.h>
#include <vector>
#include <stack>

class HE_Edge;
class TriMesh;

//=========================================================================

class HE_Vertex
{
public:
	HE_Edge* m_edge;
	int m_vid;
};

//=========================================================================

class HE_Face
{
public:
	HE_Edge* m_edge;
	int m_tid;
	Vec3f m_nml;
	int m_flag;
};

//=========================================================================

class HE_Edge
{
public:
	HE_Edge(int eid, HE_Vertex* vtx, HE_Face* face);
	~HE_Edge();

	int m_eid;
	HE_Vertex* m_vtx;
	HE_Face* m_face;
	HE_Edge* m_pair;
	HE_Edge* m_next;

	// relative direction
	int m_rdir;
	// flag for a non-maniforld edge
	int m_nonmanifold;
};

int isPair(HE_Edge* a, HE_Edge* b);

//=========================================================================

class HalfEdge
{
public:
	HalfEdge(TriMesh* o);
	~HalfEdge();

	void create();
	void selectCluster(std::vector<int>& tid_list, int tid, float angle);
	void reverseFrontFace(int tid);

	std::vector< std::vector<HE_Edge*> > m_triedge;
	std::vector< std::vector<int> > m_cluster;

private:
	void clearFaceFlag();
	int findPairTri(int vid0, int tid0, const std::vector<int>& tid_list0, std::vector<int>& tid_list);
	template<class F>
	void traverse(std::vector<int>& tid_list, int tid, F& query);

	TriMesh* m_mesh;
	std::vector<HE_Vertex*> m_vdata;
	std::vector<HE_Face*> m_fdata;
};

//=========================================================================

template<class F>
void
HalfEdge::traverse(std::vector<int>& tid_list, int tid0, F& query)
{
	HE_Edge* edge;
	HE_Edge* pair;
	int tid;
	std::stack<int> tid_stack;

	tid_stack.push(tid0);
	tid_list.push_back(tid0);
	m_fdata[tid0]->m_flag = 1;

	while (tid_stack.size() > 0)
	{
		tid = tid_stack.top();
		tid_stack.pop();

		for (int i=0; i<3; i++)
		{
			edge = m_triedge[tid][i];
			pair = edge->m_pair;
			if (pair == 0)
			{
				continue;
			}
			if (query(edge, pair))
			{
				int pairtid = pair->m_face->m_tid;
				tid_stack.push(pairtid);
				tid_list.push_back(pairtid);
				m_fdata[pairtid]->m_flag = 1;
			}
		}
	}
}

#endif  // halfedge_h

