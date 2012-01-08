#include <string>
#include <tt/util/message.h>
#include <tt/util/time.h>
#include <tt/shapes/mesh.h>
#include "mesh_io.h"
#include "fmt/obj.h"
#include "fmt/off.h"
#include "fmt/ply.h"
#include "fmt/stl.h"
#include "fmt/tm.h"

using namespace std;

typedef int (*f_mesh_scanner_type)(TriMeshInfo* info, const char* fname, const char* opt);
typedef int (*f_mesh_loader_type)(TriMesh* mesh, const char* fname, const char* opt);
typedef int (*f_mesh_saver_type)(const TriMesh* mesh, const char* fname, const char* opt);

struct MeshFunc
{
	const char* format;
	f_mesh_scanner_type scanner;
	f_mesh_loader_type loader;
	f_mesh_saver_type saver;
	const char* extensions;
	const char* comments;
};

MeshFunc mesh_func[] =
{
	{"obj_a", &obj_scan,    &obj_load,    &obj_save,   ".obj",  "wavefront ascii format"},
	{"vobj",  &obj_scan,    &obj_load,    &obj_save,   ".vobj", "vobj ascii format"},
	{"off",   &off_scan,    &off_load,    &off_save,   ".off",  "off format"},
	{"ply_a", &ply_a_scan,  &ply_a_load,  &ply_a_save, ".plya", "stanford ply ascii format"},
	{"ply_b", &ply_b_scan,  &ply_b_load,  &ply_b_save, ".ply",  "stanford ply binary format"},
	{"stl_a", &stl_a_scan,  &stl_a_load,  &stl_a_save, ".stla", "stl ascii format"},
	{"stl_b", &stl_b_scan,  &stl_b_load,  &stl_b_save, ".stl",  "stl binary format"},
	{"tm",    &tm_scan,     &tm_load,     &tm_save,    ".tm",   "tm binary format"},
	{0},
};

// guess the format from a given file name.
std::string
f_mesh_guess_format(const TtFileName& fname)
{
	string ext = fname.ext();

	for (int i=0; mesh_func[i].format; i++)
	{
		if (string(mesh_func[i].extensions) == ext)
		{
			return mesh_func[i].format;
		}
	}

	return "";
}

int
f_mesh_scan(TriMeshInfo* info, const TtFileName& fname, std::string fmt)
{
	if (fmt == "")
	{
		fmt = f_mesh_guess_format(fname);
		if (fmt == "")
		{
			cerr << "f_mesh_scan: unknown file extension" << endl;
			return 0;
		}
	}

	for (int i=0; mesh_func[i].format; i++)
	{
		if (string(mesh_func[i].format) == fmt)
		{
			return mesh_func[i].scanner(info, fname.c_str(), 0);
		}
	}

	cerr << "f_mesh_scan: unknown file format: " << fmt << endl;

	return 0;
}

int
f_mesh_load(TriMesh* mesh, const TtFileName& fname, std::string fmt)
{
	if (fmt == "")
	{
		fmt = f_mesh_guess_format(fname);
		if (fmt == "")
		{
			cerr << "f_mesh_load: unknown file extension" << endl;
			return 0;
		}
	}

	for (int i=0; mesh_func[i].format; i++)
	{
		if (string(mesh_func[i].format) == fmt)
		{
			return mesh_func[i].loader(mesh, fname.c_str(), 0);
		}
	}

	cerr << "f_mesh_load: unknown file format: " << fmt << endl;

	return 0;
}

int
f_mesh_save(const TriMesh* mesh, const TtFileName& fname, std::string fmt)
{
	if (fmt == "")
	{
		fmt = f_mesh_guess_format(fname);
		if (fmt == "")
		{
			cerr << "f_mesh_save: unknown file extension" << endl;
			return 0;
		}
	}

	for (int i=0; mesh_func[i].format; i++)
	{
		if (string(mesh_func[i].format) == fmt)
		{
			return mesh_func[i].saver(mesh, fname.c_str(), 0);
		}
	}

	cerr << "f_mesh_save: unknown file format: " << fmt << endl;

	return 0;
}

//=========================================================================

TriMeshIO::TriMeshIO(TriMesh* mesh)
{
	m_mesh = mesh;
}

TriMeshIO::TriMeshIO(const TriMesh* mesh)
{
	m_mesh = (TriMesh*)mesh;
}

int
TriMeshIO::load(const TtFileName& fname, std::string fmt)
{
	TtTime tm;
	int r;

	tt_info("[mesh] loading... (%s %s)\n", fname.c_str(), fmt.c_str());

	tm.start();

	r = f_mesh_load(m_mesh, fname, fmt);

	if (r)
	{
		tm_finalizeMesh(m_mesh);
	}

	tm.end();

	tt_info("[mesh] loading... done (%.1f sec)\n", tm.getElapsedSec());

	return r;
}

int
TriMeshIO::save(const TtFileName& fname, std::string fmt) const
{
	TtTime tm;
	int r;

	tt_info("[mesh] saving... (%s %s)\n", fname.c_str(), fmt.c_str());

	tm.start();
	r = f_mesh_save(m_mesh, fname, fmt);
	tm.end();

	tt_info("[mesh] saving... done (%.1f sec)\n", tm.getElapsedSec());

	return r;
}

