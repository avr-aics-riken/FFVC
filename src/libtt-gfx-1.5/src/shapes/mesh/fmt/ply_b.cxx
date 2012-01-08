#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <tt/util/file.h>
#include <tt/shapes/mesh.h>
#include "ply.h"

#define BUFFER_SIZE 1024

using namespace std;

int
ply_b_scan(TriMeshInfo* info, const char* fname, const char* opt)
{
	ifstream is(fname);
	if (is.fail()) return 0;

	char line[BUFFER_SIZE];

	while (is.getline(line, BUFFER_SIZE) && !is.eof())
	{
		istringstream istr(line);
		string token;

		istr >> token;
		if (token == "element")
		{
			istr >> token;
			if (token == "vertex")
			{
				istr >> info->vtx_size;
			}
			else if (token == "face")
			{
				istr >> info->fce_size;
			}
		}
		else if (token == "end_header")
		{
			break;
		}
	}

	return 1;
}

int
ply_b_load(TriMesh* mesh, const char* fname, const char* opt)
{
	ifstream is(fname);
	if (is.fail()) return 0;

	int inv = tt_check_machine_endian() == TT_BIG_ENDIAN ? 0 : 1;

	TriMeshInfo info;
	ply_b_scan(&info, fname);

	mesh->m_vdata.resize(info.vtx_size);

	char line[BUFFER_SIZE];

	while (is.getline(line, BUFFER_SIZE) && !is.eof())
	{
		if (string(line) == "end_header")
		{
			break;
		}
	}

	// read vertices
	for (int i=0; i<info.vtx_size; i++)
	{
		float vtx[3];
		tt_read(is, vtx, sizeof(float), 3, inv);
		mesh->m_vdata[i] = Vec3f(vtx[0], vtx[1], vtx[2]);
	}

	// read faces
	vector<Vec3i> idx_array;
	for (int i=0; i<info.fce_size; i++)
	{
		uchar m;
		tt_read(is, &m, sizeof(uchar), 1, inv);

		int* idx = new int[m];
		tt_read(is, idx, sizeof(int), m, inv);

		// ply format is counter clock-wise
		for (int j=2; j<m; j++)
		{
			idx_array.push_back(Vec3i(idx[0], idx[j-1], idx[j]));
		}
		delete [] idx;
	}

	copy_data(mesh->m_trivid, idx_array);

	return 1;
}

int
ply_b_save(const TriMesh* mesh, const char* fname, const char* opt)
{
	cerr << "ply_b_save: not implemented" << endl;
	return 0;
}

