#include <fstream>
#include <tt/util/file.h>
#include <tt/shapes/mesh.h>
#include "stl.h"

#define STLHEAD 80	// header size for STL binary

using namespace std;

int
stl_b_scan(TriMeshInfo* info, const char* fname, const char* opt)
{
    ifstream ifs(fname, ios::in | ios::binary); 
	if (ifs.fail()) return 0;

	int inv = tt_check_machine_endian() == TT_LITTLE_ENDIAN ? 0 : 1;

	uint element = 0;

	char buf[STLHEAD];
    for (int i=0; i<STLHEAD; i++) buf[i] = 0;
	tt_read(ifs, buf, sizeof(char), STLHEAD, inv);
	tt_read(ifs, &element, sizeof(uint), 1, inv);

	info->vtx_size = element * 3;
	info->tri_size = element;

	return 1;
}

int
stl_b_load(TriMesh* mesh, const char* fname, const char* opt)
{
    ifstream ifs(fname, ios::in | ios::binary); 
	if (ifs.fail()) return 0;

	int inv = tt_check_machine_endian() == TT_LITTLE_ENDIAN ? 0 : 1;

	uint element = 0;

	char buf[STLHEAD];
    for (int i=0; i<STLHEAD; i++) buf[i] = 0;
	tt_read(ifs, buf, sizeof(char), STLHEAD, inv);
	tt_read(ifs, &element, sizeof(uint), 1, inv);

	mesh->m_vdata.resize(element*3);
	mesh->m_ndata.resize(element*3);
	mesh->m_trivid.resize(element);
	mesh->m_trinid.resize(element);

    float nml[3];
	float vtx[3];
    ushort padding=0;
    for (uint i=0; i<element; i++)
	{
		int t = i * 3;

		// one plane normal
		tt_read(ifs, nml, sizeof(float), 3, inv);
		mesh->m_ndata[t+0] = Vec3f(nml[0], nml[1], nml[2]);
		mesh->m_ndata[t+1] = Vec3f(nml[0], nml[1], nml[2]);
		mesh->m_ndata[t+2] = Vec3f(nml[0], nml[1], nml[2]);
		mesh->m_trinid[i] = Vec3i(t, t+1, t+2);

		// three vertices
		for (int j=0; j<3; j++)
		{
			tt_read(ifs, vtx, sizeof(float), 3, inv);
			mesh->m_vdata[t+j] = Vec3f(vtx[0], vtx[1], vtx[2]);
		}
		mesh->m_trivid[i] = Vec3i(t, t+1, t+2);

		// padding 2 bytes
		tt_read(ifs, &padding, sizeof(ushort), 1, inv);
    }

	return 1;
}

int
stl_b_save(const TriMesh* mesh, const char* fname, const char* opt)
{
    ofstream ofs(fname, ios::out | ios::binary); 
	if (ofs.fail()) return 0;
    
	int inv = tt_check_machine_endian() == TT_LITTLE_ENDIAN ? 0 : 1;

	uint element = mesh->nTriangles();

    char buf[STLHEAD];
    for (int i=0; i<STLHEAD; i++) buf[i] = 0;
    strcpy(buf, "default");
	tt_write(ofs, buf, 1, STLHEAD, inv);
	tt_write(ofs, &element, sizeof(uint), 1, inv);
    
    ushort padding=0;
	for (uint m=0; m<element; m++)
	{
		// one plane normal
		tt_write(ofs, mesh->fnormal(m).ptr(), sizeof(float), 3, inv);

		// three vertices
		tt_write(ofs, mesh->vertex(m, 0).ptr(), sizeof(float), 3, inv);
		tt_write(ofs, mesh->vertex(m, 1).ptr(), sizeof(float), 3, inv);
		tt_write(ofs, mesh->vertex(m, 2).ptr(), sizeof(float), 3, inv);

		// padding 2 bytes
		tt_write(ofs, &padding, sizeof(ushort), 1, inv);
	}

	return 1;
}

