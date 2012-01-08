#include <math.h>
#include <tt/shapes/mesh.h>
#include "octree.h"
#include "itsect_axis.h"

void f_compute_itsect_on_edge(const TriMesh* mesh, const VxOctree* octree, const VxOctreeNode* node, int hit[12], float P[12][3], float N[12][3])
{
	if (!node->isBoundary()) return;

	BBox bbox = octree->getBBox(node);
	int edge_org[12] = {0, 2, 6, 4, 0, 4, 5, 1, 0, 1, 3, 2};
	for (int i=0; i<12; i++)
	{
		hit[i] = 0;
		vec3f_copy(P[i], bbox.getPoint(edge_org[i]));
	}
	Vec3f length = bbox.size();

	VxOctreeBoundaryData* data = node->getBoundaryData();
	int nelm = data->m_nelm;
	int* elms = data->m_elms;
	for (int i=0; i<nelm; i++)
	{
		int tid = elms[i];
		const float* p0 = mesh->vertex(tid, 0);
		const float* p1 = mesh->vertex(tid, 1);
		const float* p2 = mesh->vertex(tid, 2);
		int r;
		float* p;
		float t;

		for (int j=0; j<12; j++)
		{
			int axis = j / 4;
			p = P[j];
			r = f_itsect_tri_axis(axis, p0, p1, p2, &p[0], &p[1], &p[2]);
			if (r)
			{
				t = (p[axis] - bbox.min[axis])/length[axis];
				if (t >= 0 && t < 1)
				{
					hit[j] = 1;
					vec3f_copy(N[j], mesh->fnormal(tid));
				}
			}
		}
	}
}

void f_compute_tangent_plane(const int hit[12], const float P[12][3], const float N[12][3], float _p[3], float _n[3])
{
	Vec3f p(0, 0, 0);
	Vec3f n(0, 0, 0);
	int count = 0;
	for (int i=0; i<12; i++)
	{
		if (hit[i])
		{
			p += Vec3f(P[i]);
			n += Vec3f(N[i]);
			count++;
		}
	}
	p /= count;
	n.normalize();
	vec3f_copy(_p, p);
	vec3f_copy(_n, n);
}

// p[4][3]: position of corners
// c[3]: position of a tangent plane
// n[3]: normal of a tangent plane
double f_compute_area(const float p[4][3], const float c[3], const float n[3]){
    double area = 0;
    int i;
	int start = 0;
	float _f1[4];
	float _f2[4];
    for(i=0; i<4; i++){
      const float* p1 = p[i];
      const float* p2 = p[(i+1)%4];
      
      _f1[i] = n[0]*(p1[0]-c[0]) + n[1]*(p1[1]-c[1]) + n[2]*(p1[2]-c[2]);
      _f2[i] = n[0]*(p2[0]-c[0]) + n[1]*(p2[1]-c[1]) + n[2]*(p2[2]-c[2]);

      if(_f1[i] < 0 && _f2[i] >= 0){
		  start = i;
	  }
	}
    bool first = true;
    float q0[3];
    for(i=start; i<start+4; i++){
      const float* p1 = p[i%4];
      const float* p2 = p[(i+1)%4];
      
      float f1 = _f1[i%4];
      float f2 = _f2[i%4];
      
      float q1[3], q2[3];
      if(f1 >= 0 && f2 >= 0){
        q1[0] = p1[0];
        q1[1] = p1[1];
        q1[2] = p1[2];
        
        q2[0] = p2[0];
        q2[1] = p2[1];
        q2[2] = p2[2];
      }
      else if(f1 >= 0 && f2 < 0){
        q1[0] = p1[0];
        q1[1] = p1[1];
        q1[2] = p1[2];
        
        float w1 = (float)(fabs(f2)/(fabs(f1)+fabs(f2)));
        float w2 = (float)(fabs(f1)/(fabs(f1)+fabs(f2)));
        q2[0] = w1*p1[0] + w2*p2[0];
        q2[1] = w1*p1[1] + w2*p2[1];
        q2[2] = w1*p1[2] + w2*p2[2];
      }
      else if(f1 < 0 && f2 >= 0){
        float w1 = (float)(fabs(f2)/(fabs(f1)+fabs(f2)));
        float w2 = (float)(fabs(f1)/(fabs(f1)+fabs(f2)));
        q1[0] = w1*p1[0] + w2*p2[0];
        q1[1] = w1*p1[1] + w2*p2[1];
        q1[2] = w1*p1[2] + w2*p2[2];
        
        q2[0] = p2[0];
        q2[1] = p2[1];
        q2[2] = p2[2];
      }
      else
        continue;
      
      if(first){
        first = false;
        q0[0] = q1[0];
        q0[1] = q1[1];
        q0[2] = q1[2];
      }
      else{
        float v1[3] = {q1[0] - q0[0], q1[1] - q0[1], q1[2] - q0[2]};
        float v2[3] = {q2[0] - q0[0], q2[1] - q0[1], q2[2] - q0[2]};
        
        double nx = v1[1]*v2[2] - v1[2]*v2[1];
        double ny = v1[2]*v2[0] - v1[0]*v2[2];
        double nz = v1[0]*v2[1] - v1[1]*v2[0];
        
        area += sqrt(nx*nx + ny*ny + nz*nz);
      }
    }
    
    return 0.5*area;
}

void f_compute_ratio(const BBox& bbox, const float p[3], const float n[3], float* volumeRatio, float areaRatio[6])
{
	int i, j;
	double area[6];
	double volume = 0;
	float c[4][3];
	int face[6][4] = {
		{0, 2, 6, 4},
		{1, 3, 7, 5},
		{0, 4, 5, 1},
		{2, 6, 7, 3},
		{0, 1, 3, 2},
		{4, 5, 7, 6},
	};

	for (j=0; j<6; j++)
	{
		int axis = j / 2;
		for (i=0; i<4; i++)
		{
			vec3f_copy(c[i], bbox.getPoint(face[j][i]));
		}
		area[j] = f_compute_area(c, p, n);
		volume += area[j]*fabs(c[0][axis]-p[axis]);
	}

	volume /= 3.;

	Vec3f space = bbox.size();
	*volumeRatio = (float)(volume/(space[0]*space[1]*space[2]));
	areaRatio[0] = (float)(area[0]/(space[1]*space[2]));
	areaRatio[1] = (float)(area[1]/(space[1]*space[2]));
	areaRatio[2] = (float)(area[2]/(space[2]*space[0]));
	areaRatio[3] = (float)(area[3]/(space[2]*space[0]));
	areaRatio[4] = (float)(area[4]/(space[0]*space[1]));
	areaRatio[5] = (float)(area[5]/(space[0]*space[1]));
}

void f_compute_ratio(const TriMesh* mesh, const VxOctree* octree, const VxOctreeNode* node, float* volumeRatio, float areaRatio[6])
{
	int hit[12];
	float P[12][3];
	float N[12][3];
	f_compute_itsect_on_edge(mesh, octree, node, hit, P, N);

	float p[3];
	float n[3];
	f_compute_tangent_plane(hit, P, N, p, n);

	BBox bbox = octree->getBBox(node);
	f_compute_ratio(bbox, p, n, volumeRatio, areaRatio);
}

