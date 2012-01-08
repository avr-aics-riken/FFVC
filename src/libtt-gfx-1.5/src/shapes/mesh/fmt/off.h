#ifndef off_h
#define off_h

class TriMesh;
class TriMeshInfo;

int off_scan(TriMeshInfo* info, const char* fname, const char* opt="");
int off_load(TriMesh* mesh, const char* fname, const char* opt="");
int off_save(const TriMesh* mesh, const char* fname, const char* opt="");

#endif  // off_h

