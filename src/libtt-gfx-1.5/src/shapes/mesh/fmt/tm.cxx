#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <tt/util/file.h>
#include <tt/shapes/mesh.h>
#include "tm.h"

#define BUFFER_SIZE 1024

using namespace std;

int
tm_scan(TriMeshInfo* info, const char* fname, const char* opt)
{
	ifstream is(fname);
	if (is.fail()) return 0;

	char line[BUFFER_SIZE];

	while (is.getline(line, BUFFER_SIZE) && !is.eof())
	{
		istringstream istr(line);
		string token;

		istr >> token;

		if (token == "nvtx")
		{
			istr >> info->vtx_size;
		}
		else if (token == "nnml")
		{
			istr >> info->nml_size;
		}
		else if (token == "ntex")
		{
			istr >> info->tex_size;
		}
		else if (token == "namb")
		{
			istr >> info->amb_size;
		}
		else if (token == "ntri")
		{
			istr >> info->tri_size;
		}
		else if (token == "nmtl")
		{
			istr >> info->mtl_size;
		}
		else if (token == "end_header")
		{
			break;
		}
	}

	return 1;
}

int
tm_load(TriMesh* mesh, const char* fname, const char* opt)
{
	ifstream is(fname);
	if (is.fail()) return 0;

	int inv = tt_check_machine_endian() == TT_LITTLE_ENDIAN ? 0 : 1;

	TriMeshInfo info;
	tm_scan(&info, fname);

	int nvtx = info.vtx_size;
	int nnml = info.nml_size;
	int ntex = info.tex_size;
	int namb = info.amb_size;
	int ntri = info.tri_size;

	char line[BUFFER_SIZE];

	while (is.getline(line, BUFFER_SIZE) && !is.eof())
	{
		istringstream istr(line);
		string token;

		istr >> token;

		if (token == "mtl")
		{
			TriMeshGroup* grp = new TriMeshGroup();
			istr >> grp->m_name >> grp->m_first >> grp->m_last;
			mesh->addGroup(grp);
		}
		else if (token == "smooth")
		{
			istr >> mesh->m_smooth_angle;
		}
		else if (token == "end_header")
		{
			break;
		}
	}

	if (nvtx > 0)
	{
		mesh->m_vdata.resize(nvtx);
		mesh->m_trivid.resize(ntri);
		tt_read(is, mesh->m_vdata.ptr(), sizeof(float), 3*nvtx, inv);
		tt_read(is, mesh->m_trivid.ptr(), sizeof(int), 3*ntri, inv);
	}

	if (nnml > 0)
	{
		mesh->m_ndata.resize(nnml);
		mesh->m_trinid.resize(ntri);
		tt_read(is, mesh->m_ndata.ptr(), sizeof(float), 3*nnml, inv);
		tt_read(is, mesh->m_trinid.ptr(), sizeof(int), 3*ntri, inv);
	}

	if (ntex > 0)
	{
		mesh->m_texdata.resize(ntex);
		mesh->m_tritexid.resize(ntri);
		tt_read(is, mesh->m_texdata.ptr(), sizeof(float), 2*ntex, inv);
		tt_read(is, mesh->m_tritexid.ptr(), sizeof(int), 3*ntri, inv);
	}

	if (namb > 0)
	{
		mesh->m_adata.resize(namb);
		tt_read(is, mesh->m_adata.ptr(), sizeof(float), 4*namb, inv);
	}

	return 1;
}

int
tm_save(const TriMesh* mesh, const char* fname, const char* opt)
{
	ofstream os(fname);
	if (os.fail()) return 0;

	int inv = tt_check_machine_endian() == TT_LITTLE_ENDIAN ? 0 : 1;

	int save_vtx = 1;
	int save_nml = 0;
	int save_tex = 0;
	int save_amb = 1;
	int save_mtl = 0;

	int nvtx = mesh->nVertices();
	int nnml = mesh->nNormals();
	int ntex = mesh->nTexCoords();
	int namb = mesh->m_adata.size();
	int ntri = mesh->nTriangles();
	int nmtl = mesh->nMaterials();

	if (save_vtx && nvtx > 0)
	{
		os << "nvtx " << nvtx << endl;
	}
	if (save_nml && nnml > 0)
	{
		os << "nnml " << nnml << endl;
	}
	if (save_tex && ntex > 0)
	{
		os << "ntex " << ntex << endl;
	}
	if (save_amb && namb > 0)
	{
		os << "namb " << namb << endl;
	}
	if (ntri > 0)
	{
		os << "ntri " << ntri << endl;
	}
	if (save_mtl && nmtl > 0)
	{
		os << "nmtl " << nmtl << endl;
		for (int i=0; i<nmtl; i++)
		{
			const TriMeshGroup* grp = mesh->getGroup(i);
			os << "mtl " << grp->m_name << " " << grp->m_first << " " << grp->m_last << endl;
		}
	}

	float smooth_angle = mesh->m_smooth_angle;
	/*
	if (smooth_angle > 0)
	{
		os << "smooth " << smooth_angle << endl;
	}
	*/

	os << "end_header" << endl;

	if (save_vtx && nvtx > 0)
	{
		tt_write(os, mesh->m_vdata.ptr(), sizeof(float), 3*nvtx, inv);
		tt_write(os, mesh->m_trivid.ptr(), sizeof(int), 3*ntri, inv);
	}

	if (save_nml && nnml > 0)
	{
		tt_write(os, mesh->m_ndata.ptr(), sizeof(float), 3*nnml, inv);
		tt_write(os, mesh->m_trinid.ptr(), sizeof(int), 3*ntri, inv);
	}

	if (save_tex && ntex > 0)
	{
		tt_write(os, mesh->m_texdata.ptr(), sizeof(float), 2*ntex, inv);
		tt_write(os, mesh->m_tritexid.ptr(), sizeof(int), 3*ntri, inv);
	}

	if (save_amb && namb > 0)
	{
		tt_write(os, mesh->m_adata.ptr(), sizeof(float), 4*namb, inv);
	}

	return 1;
}

