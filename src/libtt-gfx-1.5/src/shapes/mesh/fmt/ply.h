#ifndef ply_h
#define ply_h

class TriMesh;
class TriMeshInfo;

// stanford ply ascii format
int ply_a_scan(TriMeshInfo* info, const char* fname, const char* opt="");
int ply_a_load(TriMesh* mesh, const char* fname, const char* opt="");
int ply_a_save(const TriMesh* mesh, const char* fname, const char* opt="");

// stanford ply binary format
int ply_b_scan(TriMeshInfo* info, const char* fname, const char* opt="");
int ply_b_load(TriMesh* mesh, const char* fname, const char* opt="");
int ply_b_save(const TriMesh* mesh, const char* fname, const char* opt="");

#endif  // ply_h

