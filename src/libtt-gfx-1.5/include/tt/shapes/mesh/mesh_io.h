#ifndef mesh_io_h
#define mesh_io_h

class TriMesh;
class TriMeshInfo;
class TtFileName;

int f_mesh_scan(TriMeshInfo* info, const TtFileName& fname, std::string fmt="");
int f_mesh_load(TriMesh* mesh, const TtFileName& fname, std::string fmt="");
int f_mesh_save(const TriMesh* mesh, const TtFileName& fname, std::string fmt="");

//=========================================================================

class TriMeshIO
{
public:
	TriMeshIO(TriMesh* mesh);
	TriMeshIO(const TriMesh* mesh);

	int load(const TtFileName& fname, std::string fmt="");
	int save(const TtFileName& fname, std::string fmt="") const;

private:
	TriMesh* m_mesh;
};

#endif  // mesh_io_h

