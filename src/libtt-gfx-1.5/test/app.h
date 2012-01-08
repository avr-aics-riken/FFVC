#ifndef app_h
#define app_h

#include <tt/tt.h>
#include <tt/gm.h>
#include <tt/util/filename.h>

class TriMesh;
class VxVoxelizerOctree;

class Application
{
public:
	Application();
	virtual ~Application();

	TriMesh* getCurrMesh();
	TriMesh* getMesh(int i);

	int  loadObject(const TtFileName& fname, std::string opt = "");
	int  saveObject(const TtFileName& fname, std::string opt = "");
	int  saveVoxel(const TtFileName& fname, std::string opt = "");

	void fixObject(int oid, float dist, float angle);
	void mergeObjectVertices(int oid, float dist);
	void mergeObjectNormals(int oid);
	void uniformObjectFrontFaces(int oid);
	void computeObjectNormals(int oid, float angle);
	void reverseGroupNormals(int gid);

	std::string computeVoxelParam(int axis, int axis_nelm);
	std::string computeVoxelParamByPitch(const Vec3f& pitch);
	void voxelizePolygon(const Vec3f& orig, const Vec3f& pitch, const Vec3i& elm);

	void setRatioType(int type);


	TriMesh* m_mesh;
	VxVoxelizerOctree* m_vox;

private:
};

#endif  // app_h

