#include <sstream>
#include <tt/util/time.h>
#include <tt/util/message.h>
#include <tt/shapes/mesh.h>
#include <tt/shapes/mesh/mesh_edit.h>
#include <tt/vx/voxelizer_octree.h>
#include <tt/vx/fmt/svx.h>
#include <tt/vx/fmt/ovx.h>
#include <tt/vx/fmt/sph.h>
#include "app.h"

using namespace std;


Application::Application()
{
	m_mesh = 0;
	m_vox = 0;
}

Application::~Application()
{
	if (m_mesh) delete m_mesh;
	if (m_vox) delete m_vox;
}

TriMesh*
Application::getCurrMesh()
{
	return getMesh(0);
}

TriMesh*
Application::getMesh(int i)
{
	return m_mesh;
}

int
Application::loadObject(const TtFileName& fname, std::string opt)
{
	if (m_mesh) {
		delete m_mesh; m_mesh = 0;
	}

	m_mesh = new TriMesh();

	int r = m_mesh->load(fname, opt);
	if (!r) {
		delete m_mesh; m_mesh = 0;
		return 0;
	}

	return 1;
}

int
Application::saveObject(const TtFileName& fname, std::string fmt)
{
	TriMesh* mesh = getCurrMesh();
	if (mesh)
		return mesh->save(fname);

	return 0;
}

int
Application::saveVoxel(const TtFileName& fname, std::string opt)
{
	int r = 0;

	tt_info("[app] saving... voxel (%s %s)\n", fname.c_str(), opt.c_str());
	TtTime tm;
	tm.start();
	if (m_vox)
	{
		if (fname.ext() == ".svx")
		{
			r = svx_save(m_vox, fname.c_str(), opt.c_str());
		}
		else if (fname.ext() == ".ovx")
		{
			r = ovx_save(m_vox, fname.c_str(), opt.c_str());
		}
		else if (fname.ext() == ".sph")
		{
			r = sph_save(m_vox, fname.c_str(), opt.c_str());
		}
	}
	tm.end();
	tt_info("[app] saving... voxel (%.2f)\n", tm.getElapsedSec());

	return r;
}

void
Application::fixObject(int oid, float dist, float angle)
{
	TriMesh* mesh = getMesh(oid);
	if (!mesh) return;
	TriMeshEditor me(mesh);
	me.fix(dist * dist, angle);
}

void
Application::mergeObjectVertices(int oid, float dist)
{
	TriMesh* mesh = getMesh(oid);
	if (!mesh) return;
	TriMeshEditor me(mesh);
	me.mergeVertices(dist * dist);
}

void
Application::mergeObjectNormals(int oid)
{
	TriMesh* mesh = getMesh(oid);
	if (!mesh) return;
	TriMeshEditor me(mesh);
	me.mergeNormals();
}

void
Application::uniformObjectFrontFaces(int oid)
{
	TriMesh* mesh = getMesh(oid);
	if (!mesh) return;
	TriMeshEditor me(mesh);
	me.uniformFrontFaces();
}

void
Application::computeObjectNormals(int oid, float angle)
{
	TriMesh* mesh = getMesh(oid);
	if (!mesh) return;
	TriMeshEditor me(mesh);
	me.computeNormals(angle);
}

void
Application::reverseGroupNormals(int gid)
{
	TriMesh* mesh = getCurrMesh();
	if (!mesh) return;
	TriMeshEditor me(mesh);
	me.reverseNormals(gid);
}

std::string
Application::computeVoxelParam(int axis, int axis_nelm)
{
	TriMesh* mesh = getCurrMesh();
	if (!mesh) return "";

	BBox bbox = mesh->getBBox();

	float pitch = bbox.length(axis) / axis_nelm;
	Vec3f nelm3 = (bbox.size() + Vec3f(.5*pitch)) / pitch;
	int max_nelm = gmMax3(nelm3[0], nelm3[1], nelm3[2]);
	int max_depth = gmCeil(log2f(max_nelm));
	int nelm = int(powf(2, max_depth));

	ostringstream os;
	os << bbox.min << " " << Vec3f(pitch) << " " << Vec3f(nelm);
	return os.str();
}

void
Application::voxelizePolygon(const Vec3f& orig, const Vec3f& pitch, const Vec3i& nelm)
{
	TtTime tm;
	float t;
	int nx, ny, nz;
	nx = nelm[0];
	ny = nelm[1];
	nz = nelm[2];

	tt_info("[app] voxelizing... (%d x %d x %d = %d)\n", nx, ny, nz, nx*ny*nz);

	TriMesh* mesh = getCurrMesh();
	if (!mesh) return;

	Vec3f min = orig;
	Vec3f max = orig + multi(pitch, Vec3f(nx, ny, nz));
	BBox bbox = BBox(min, max);
	int max_depth = int(log2(nx));

	tm.start();

	int id = 1;
	delete m_vox;
	m_vox = new VxVoxelizerOctree(bbox, max_depth);
	m_vox->voxelize(mesh, id);
	id = m_vox->getNumId();

	tm.end();
	t = tm.getElapsedSec();

	tt_info("[app] voxelizing... done (%.2f sec, %d ids)\n", t, id);
}

void
Application::setRatioType(int type)
{
	m_vox->setRatioType(type);
}

