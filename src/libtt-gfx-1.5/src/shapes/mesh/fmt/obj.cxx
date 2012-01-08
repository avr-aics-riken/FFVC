#include <fstream>
#include <sstream>
#include <vector>
#include <tt/shapes/mesh.h>
#include "obj.h"

#define BUFFER_SIZE 1024

using namespace std;

///////////////////////////////////////////////////////////////////////////
//
// class Index3
//
///////////////////////////////////////////////////////////////////////////

class Index3
{
public:
	Index3() {}
	Index3(int v) : idx3(v, v, v) {}
	Index3(int v, int vt, int vn) : idx3(v, vt, vn) {}
	Index3(const Vec3i& x) : idx3(x) {}

	Index3& operator +=(const Index3& x)
	{
		idx3 += x.idx3;
		return *this;
	}
	Index3& operator -=(const Index3& x)
	{
		idx3 -= x.idx3;
		return *this;
	}
	Index3  operator + (const Index3& x) const
	{
		return Index3(idx3 + x.idx3);
	}
	Index3  operator - (const Index3& x) const
	{
		return Index3(idx3 - x.idx3);
	}
	int& operator [](int i)
	{
		return idx3[i];
	}
	const int& operator [](int i) const
	{
		return idx3[i];
	}

private:
	Vec3i idx3;
};

inline
std::istream& operator >>(std::istream& is, Index3& o)
{
	is >> o[0];
	if (is.get() == '/')
	{
		if (is.peek() != '/')
		{
			is >> o[1];
		}
		if (is.get() == '/')
		{
			is >> o[2];
		}
	}

	return is;
}

inline
std::ostream& operator <<(std::ostream& os, const Index3& o)
{
	os << o[0];
	if (o[1] != 0)
	{
		os << "/" << o[1];
	}
	else if (o[2] != 0)
	{
		os << "/";
	}
	if (o[2] != 0)
	{
		os << "/" << o[2];
	}
	return os;
}

///////////////////////////////////////////////////////////////////////////

int
obj_scan(TriMeshInfo* info, const char* fname, const char* opt)
{
	ifstream _is(fname);
	if (_is.fail()) return 0;

	char line[BUFFER_SIZE];

	while (_is.getline(line, BUFFER_SIZE) && !_is.eof())
	{
		istringstream is(line);
		string token;

		is >> token;

		if (token == "v")
		{
			info->vtx_size++;
		}
		else if (token == "vt")
		{
			info->tex_size++;
		}
		else if (token == "vn")
		{
			info->nml_size++;
		}
		else if (token == "f")
		{
			info->fce_size++;

			char buffer[BUFFER_SIZE];
			is.getline(buffer, BUFFER_SIZE);
			istringstream istr(buffer);

			int n = 0;
			while (1)
			{
				Index3 idx3; 
				istr >> idx3;
				if (idx3[0] == 0) break;
				n++;
			}
			info->tri_size += n - 2;
		}
		else if (token == "g")
		{
			info->grp_size++;
		}
		else if (token == "usemtl")
		{
			info->mtl_size++;
		}
	}

	return 1;
}

int
obj_load(TriMesh* mesh, const char* fname, const char* opt)
{
	ifstream _is(fname);
	if (_is.fail()) return 0;

	int     linenr = 0;
	char    line[BUFFER_SIZE];
	int     nfce = 0;
	int     clockwise = 0;

	string grp_name;
	string mtl_name;

	TriMeshInfo info;
	obj_scan(&info, fname);

	mesh->m_vdata.reallocate(0, info.vtx_size);
	mesh->m_texdata.reallocate(0, info.tex_size);
	mesh->m_ndata.reallocate(0, info.nml_size);
	if (info.vtx_size > 0)
	{
		mesh->m_trivid.reallocate(0, info.tri_size);
	}
	if (info.tex_size > 0)
	{
		mesh->m_tritexid.reallocate(0, info.tri_size);
	}
	if (info.nml_size > 0)
	{
		mesh->m_trinid.reallocate(0, info.tri_size);
	}

	while (_is.getline(line, BUFFER_SIZE) && !_is.eof())
	{
		linenr++;
		istringstream is(line);
		string token;

		is >> token;

		if (token == "v")
		{
			Vec3f vtx;
			is >> vtx;
			mesh->m_vdata.push_back(vtx);
		}
		else if (token == "vt")
		{
			Vec2f tex;
			is >> tex;
			mesh->m_texdata.push_back(tex);
		}
		else if (token == "vn")
		{
			Vec3f nml;
			is >> nml;
			mesh->m_ndata.push_back(nml);
		}
		else if (token == "f")
		{
			nfce++;

			char buffer[BUFFER_SIZE];
			is.getline(buffer, BUFFER_SIZE);
			istringstream istr(buffer);

			Index3 idx3[3];
			Index3 offset3(-1);
			istr >> idx3[0] >> idx3[1] >> idx3[2];
			idx3[0] += offset3;
			idx3[1] += offset3;
			idx3[2] += offset3;

			while (1)
			{
				Vec3i vid3, texid3, nid3;

				if (clockwise)
				{
					vid3.assign(idx3[0][0], idx3[2][0], idx3[1][0]);
					texid3.assign(idx3[0][1], idx3[2][1], idx3[1][1]);
					nid3.assign(idx3[0][2], idx3[2][2], idx3[1][2]);
				}
				else
				{
					vid3.assign(idx3[0][0], idx3[1][0], idx3[2][0]);
					texid3.assign(idx3[0][1], idx3[1][1], idx3[2][1]);
					nid3.assign(idx3[0][2], idx3[1][2], idx3[2][2]);
				}

				if (vid3[0] != -1)
				{
					mesh->m_trivid.push_back(vid3);
				}
				if (texid3[0] != -1)
				{
					mesh->m_tritexid.push_back(texid3);
				}
				if (nid3[0] != -1)
				{
					mesh->m_trinid.push_back(nid3);
				}

				idx3[1] = idx3[2];
				idx3[2] = Index3();
				istr >> idx3[2];
				if (idx3[2][0] == 0)
					break;
				idx3[2] += offset3;
			}
		}
		else if (token == "g")
		{
			is >> grp_name;
		}
		else if (token == "usemtl")
		{
			is >> mtl_name;

			TriMeshGroup* grp = new TriMeshGroup();
			grp->m_name = mtl_name;
			grp->m_first = mesh->m_trivid.size();
			mesh->addGroup(grp);
		}
	}

	// finalize material groups
	int ngrp = mesh->nGroups();
	for (int i=0; i<ngrp; i++)
	{
		int end;
		if (i == ngrp - 1)
		{
			end = mesh->m_trivid.size();
		}
		else
		{
			TriMeshGroup* grp = mesh->getGroup(i + 1);
			end = grp->m_first;
		}
		TriMeshGroup* grp = mesh->getGroup(i);
		grp->m_last = end - 1;
	}

	return 1;
}

int
obj_save(const TriMesh* mesh, const char* fname, const char* opt)
{
	ofstream os(fname);
	if (os.fail()) return 0;

	int i, j;
	int k = 0;

	int save_v  = 1;
	int save_vt = 0;
	int save_vn = 0;

	int nvtx = mesh->nVertices();
	int nnml = mesh->nNormals();
	int ntex = mesh->nTexCoords();
	int ntri = mesh->nTriangles();
	int ngrp = mesh->nGroups();

	os << "# Wavefront obj format" << endl;
	os << endl;

	// print vertex data
	if (save_v && nvtx > 0)
	{
		for (i=0; i<nvtx; i++)
		{
			os << "v " << mesh->m_vdata[i] << endl;
		}
		os << endl;
	}

	// print texture coordinates data
	if (save_vt && ntex > 0)
	{
		for (i=0; i<ntex; i++)
		{
			os << "vt " << mesh->m_texdata[i] << endl;
		}
		os << endl;
	}

	// print normal data
	if (save_vn && nnml > 0)
	{
		for (i=0; i<nnml; i++)
		{
			os << "vn " << mesh->m_ndata[i] << endl;
		}
		os << endl;
	}

	for (i=0; i<ntri; i++)
	{
		// print material name
		if (ngrp > 0 && k < ngrp)
		{
			TriMeshGroup* grp = mesh->getGroup(k);
			if (i == grp->m_first)
			{
				os << "usemtl " << grp->m_name << endl;
				k++;
			}
		}

		// print face
		os << "f";
		for (j=0; j<3; j++)
		{
			int vid = 0;
			int texid = 0;
			int nid = 0;

			if (save_v && nvtx > 0)
			{
				vid = mesh->m_trivid[i][j] + 1;
			}
			if (save_vt && ntex > 0)
			{
				texid = mesh->m_tritexid[i][j] + 1;
			}
			if (save_vn && nnml > 0)
			{
				nid = mesh->m_trinid[i][j] + 1;
			}
			Index3 idx3(vid, texid, nid);
			os << " " << idx3;
		}
		os << endl;
	}

	return 1;
}

