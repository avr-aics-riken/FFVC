#ifndef mesh_edit_h
#define mesh_edit_h

#include <vector>

class TriMesh;
class HalfEdge;
class Mat4x4f;

class TriMeshEditor
{
public:
	TriMeshEditor(TriMesh* mesh);

	TriMesh* getGroup(int gid);
	void addMesh(const TriMesh* mesh);

	float findMinEdgeLengthSquared();

	void fix(float sqdist, float angle);

	void mergeVertices(float sqdist);
	void mergeNormals();
	void computeNormals(float angle);
	void reverseNormals();
	void reverseNormals(int gid);
	void uniformFrontFaces();
	void heuristicFrontFaces();
	void grouping();
	void selectCluster(std::vector<int>& tid_list, int tid, float angle);

	void transform(const Mat4x4f& M);

private:
	void computeSmoothNormals(float angle);

	TriMesh* m_mesh;
};

#endif  // mesh_edit_h

