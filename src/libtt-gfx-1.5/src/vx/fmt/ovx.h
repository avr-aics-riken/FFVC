#ifndef ovx_h
#define ovx_h

class VxVoxelizerOctree;

int ovx_load(VxVoxelizerOctree* vox, const char* fname, const char* opt = "");
int ovx_save(const VxVoxelizerOctree* vox, const char* fname, const char* opt = "");

#endif  // ovx_h

