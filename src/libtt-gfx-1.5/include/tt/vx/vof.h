#ifndef vof_h
#define vof_h

void f_compute_ratio(const TriMesh* mesh, const VxOctree* octree, const VxOctreeNode* node, float* volumeRatio, float areaRatio[6]);

#endif  // vof_h

