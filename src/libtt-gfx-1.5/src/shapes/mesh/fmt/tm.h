#ifndef tm_h
#define tm_h

class TriMesh;
class TriMeshInfo;

int tm_scan(TriMeshInfo* info, const char* fname, const char* opt="");
int tm_load(TriMesh* mesh, const char* fname, const char* opt="");
int tm_save(const TriMesh* mesh, const char* fname, const char* opt="");

#endif  // tm_h

