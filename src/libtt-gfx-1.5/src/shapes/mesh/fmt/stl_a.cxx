#include <fstream>
#include <vector>
#include <iomanip>
#include <tt/shapes/mesh.h>
#include "stl.h"

#define SCIENTIFIC_OUT 0

using namespace std;

int
stl_a_scan(TriMeshInfo* info, const char* fname, const char* opt)
{
	ifstream is(fname);
	if (is.fail()) return 0;

	string token;

	while (is >> token && !is.eof())
	{
		if (token == "vertex")
		{
			info->vtx_size++;
		}
		else if (token == "facet")
		{
			info->tri_size++;
		}
	}
	info->nml_size = info->vtx_size;

	return 1;
}

int
stl_a_load(TriMesh* mesh, const char* fname, const char* opt)
{
	ifstream is(fname);
	if (is.fail()) return 0;

	TriMeshInfo info;
	stl_a_scan(&info, fname);

	mesh->m_vdata.reallocate(0, info.vtx_size);
	mesh->m_ndata.reallocate(0, info.nml_size);
	mesh->m_trivid.reallocate(0, info.tri_size);
	mesh->m_trinid.reallocate(0, info.tri_size);

	string token;

	while (is >> token && !is.eof())
	{
		if (token == "solid")
		{
			string s;
			is >> s;
		}
		else if (token == "facet")
		{
			string s;
			is >> s;
			Vec3f nml;
			is >> nml;
			nml.normalize();
			mesh->m_ndata.push_back(nml);
			mesh->m_ndata.push_back(nml);
			mesh->m_ndata.push_back(nml);
		}
		else if (token == "vertex")
		{
			Vec3f vtx;
			is >> vtx;
			mesh->m_vdata.push_back(vtx);
		}
		else if (token == "outer")
		{
			string s;
			is >> s;
		}
		else if (token == "endloop")
		{
		}
		else if (token == "endfacet")
		{
			int sz = mesh->m_trivid.size();
			Vec3i vtxidx(sz*3, sz*3+1, sz*3+2);
			Vec3i nmlidx(sz*3, sz*3+1, sz*3+2);
			mesh->m_trivid.push_back(vtxidx);
			mesh->m_trinid.push_back(nmlidx);
		}
		else if (token == "endsolid")
		{
			string s;
			is >> s;
		}
	}

	return 1;
}

int
stl_a_save(const TriMesh* mesh, const char* fname, const char* opt)
{
	ofstream os(fname);
	if (os.fail()) return 0;

	int n_triangles = mesh->nTriangles();

	os << "solid " << "model1" << endl;
	for (int i=0; i<n_triangles; i++)
	{
#if SCIENTIFIC_OUT
		os << "  facet " << "normal " << setprecision(6) << scientific << mesh->fnormal(i) << endl;
#else
		os << "  facet " << "normal " << setprecision(6) << mesh->fnormal(i) << endl;
#endif
		os << "    outer " << "loop" << endl;
		for (int j=0; j<3; j++)
		{
#if SCIENTIFIC_OUT
			os << "      vertex " << setprecision(6) << scientific << mesh->vertex(i, j) << endl;
#else
			os << "      vertex " << setprecision(6) << mesh->vertex(i, j) << endl;
#endif
		}
		os << "    endloop" << endl;
		os << "  endfacet" << endl;
	}
	os << "endsolid " << "model1" << endl;

	return 1;
}

