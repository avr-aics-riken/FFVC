#include <fstream>
#include <vector>
#include <tt/shapes/mesh.h>
#include "off.h"

using namespace std;

int
off_scan(TriMeshInfo* info, const char* fname, const char* opt) {
	cerr << "off_scan: not implemented" << endl;
	return 0;
}

int
off_load(TriMesh* mesh, const char* fname, const char* opt) {
	ifstream is(fname);
	if (is.fail()) return 0;

	string token;
	int    n_vertices, n_faces, dummy;
	int    i;

	// read OFF
	is >> token;
	is >> n_vertices >> n_faces >> dummy;

	mesh->m_vdata.resize(n_vertices);

	// vertices
	Vec3f vtx;
	for (i=0; i<n_vertices; i++) {
		is >> vtx;
		mesh->m_vdata[i] = vtx;
	}

	// faces
	Vec3i         idx;
	vector<Vec3i> idx_array;
	int          m;
	for (i=0; i<n_faces; i++) {
		is >> m >> idx;
		idx_array.push_back(idx);
		for (int j=0; j<m-3; j++) {
			idx[1] = idx[2];
			is >> idx[2];
			idx_array.push_back(idx);
		}
	}

	copy_data(mesh->m_trivid, idx_array);

	return 1;
}

int
off_save(const TriMesh* mesh, const char* fname, const char* opt) {
	cerr << "off_save: not implemented" << endl;
	return 0;
}

