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
ply_a_scan(TriMeshInfo* info, const char* fname, const char* opt)
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
ply_a_load(TriMesh* mesh, const char* fname, const char* opt)
{
	ifstream is(fname);
	if (is.fail()) return 0;

	TriMeshInfo info;
	ply_a_scan(&info, fname);

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
	Vec3f vtx;
	for (int i=0; i<info.vtx_size; i++)
	{
		is.getline(line, BUFFER_SIZE);
		istringstream istr(line);

		istr >> vtx;
		mesh->m_vdata[i] = vtx;
	}

	// read faces
	vector<Vec3i> idx_array;
	Vec3i idx;
	int m;
	for (int i=0; i<info.fce_size; i++)
	{
		is.getline(line, BUFFER_SIZE);
		istringstream istr(line);

		// ply format is counter clock-wise
		istr >> m >> idx[0] >> idx[1] >> idx[2];
		idx_array.push_back(idx);
		for (int j=0; j<m-3; j++)
		{
			idx[1] = idx[2];
			istr >> idx[2];
			idx_array.push_back(idx);
		}
	}

	copy_data(mesh->m_trivid, idx_array);

	return 1;
}

int
ply_a_save(const TriMesh* mesh, const char* fname, const char* opt)
{
	cerr << "ply_a_save: not implemented" << endl;
	return 0;
}

