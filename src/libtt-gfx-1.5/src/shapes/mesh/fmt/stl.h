#ifndef stl_h
#define stl_h

class TriMesh;
class TriMeshInfo;

// STL ascii format
int stl_a_scan(TriMeshInfo* info, const char* fname, const char* opt="");
int stl_a_load(TriMesh* mesh, const char* fname, const char* opt="");
int stl_a_save(const TriMesh* mesh, const char* fname, const char* opt="");

// STL binary format
int stl_b_scan(TriMeshInfo* info, const char* fname, const char* opt="");
int stl_b_load(TriMesh* mesh, const char* fname, const char* opt="");
int stl_b_save(const TriMesh* mesh, const char* fname, const char* opt="");

#endif  // stl_h

