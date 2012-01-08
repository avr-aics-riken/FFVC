#ifndef obj_h
#define obj_h

class TriMesh;
class TriMeshInfo;

int obj_scan(TriMeshInfo* info, const char* fname, const char* opt="");
int obj_load(TriMesh* mesh, const char* fname, const char* opt="");
int obj_save(const TriMesh* mesh, const char* fname, const char* opt="");

#endif  // obj_h

