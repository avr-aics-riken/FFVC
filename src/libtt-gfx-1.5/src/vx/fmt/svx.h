#ifndef svx_h
#define svx_h

class VxVoxelizer;

int svx_load(VxVoxelizer* vox, const char* fname, const char* opt = "");
int svx_save(const VxVoxelizer* vox, const char* fname, const char* opt = "");

#endif  // svx_h

