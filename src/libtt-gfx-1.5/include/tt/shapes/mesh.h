#ifndef tt_mesh_h
#define tt_mesh_h

#include <vector>
#include <tt/array.h>
#include <tt/util/filename.h>
#include <tt/util/ptr_set.h>
#include <tt/gfx/shape.h>

class VTree;
class HalfEdge;

// counter clock-wise
inline Vec3f f_triangle_normal(const Vec3f& p0, const Vec3f& p1, const Vec3f& p2)
{
	return cross(p1-p0, p2-p0).normalize();
}

//=========================================================================
// class TriMeshInfo
//=========================================================================

struct TriMeshInfo
{
	int vtx_size;
	int tex_size;
	int nml_size;
	int amb_size;
	int fce_size;
	int tri_size;
	int grp_size;
	int mtl_size;

	TriMeshInfo()
	{
		vtx_size = 0;
		tex_size = 0;
		nml_size = 0;
		amb_size = 0;
		fce_size = 0;
		tri_size = 0;
		grp_size = 0;
		mtl_size = 0;
	}
};

//=========================================================================
// class TriMeshGroup
//=========================================================================

struct TriMeshGroup
{
	TriMeshGroup() : m_first(0), m_last(0) {}

	bool hasRandomTID() const {
		return m_tidlist.size() > 0 ? true : false;
	}
	int nTriangles() {
		if (hasRandomTID()) return m_tidlist.size();
		return m_last - m_first + 1;
	}

	std::string m_name;
	std::vector<int> m_tidlist;
	int m_first, m_last;
};

//=========================================================================
// class TriMesh
//=========================================================================

class TriMesh : public Shape
{
public:
	TriMesh();
	~TriMesh();

	ShapeType   getTypeId()   const { return TRIANGLE_MESH; }
	std::string getTypeName() const { return "TriMesh"; }

	BBox getBBox() const { return m_bbox; }

	bool isConvex() const { return false; }

	void draw() const;

	int load(const TtFileName& fname, std::string fmt="");
	int save(const TtFileName& fname, std::string fmt="") const;

	int nVertices()  const { return m_vdata.size(); }
	int nNormals()   const { return m_ndata.size(); }
	int nTexCoords() const { return m_texdata.size(); }
	int nTriangles() const { return m_trivid.size(); }
	int nGroups()    const { return m_groups.size(); }
	int nMaterials() const { return m_groups.size(); }
	int nAmbients()  const { return m_adata.size(); }

	// vertex of the triangle
	const Vec3f& vertex(int tid, int n) const
	{ return m_vdata[m_trivid[tid][n]]; }

	// vertex texture coordinate
	const Vec2f& vtexcoord(int tid, int n) const
	{ return m_texdata[m_tritexid[tid][n]]; }

	// vertex normal
	const Vec3f& vnormal(int tid, int n) const
	{ return m_ndata[m_trinid[tid][n]]; }

	// vertex normal
	Vec3f& vnormal(int tid, int n)
	{ return m_ndata[m_trinid[tid][n]]; }

	// face normal
	Vec3f fnormal(int tid) const
	{ return f_triangle_normal(vertex(tid, 0), vertex(tid, 1), vertex(tid, 2)); }

	// get a bounding box of a triangle
	BBox getTriangleBBox(int tid) const;

	TriMeshGroup* getGroup(int i) const
	{ return m_groups.getData(i); }
	void addGroup(TriMeshGroup* grp)
	{ m_groups.addData(grp); }
	void clearGroups()
	{ m_groups.clearData(); }
	void newGroup(const std::string& mtl_name, std::vector<int>& tid_list);
	void renameGroup(int gid, const std::string& mtl_name);
	void packVID();
	void unpackVID();

	int getGroupID(int tid) const;
	const std::string& getMaterialName(int tid) const;

	VTree*    getVTree() const { return m_vtree; }
	HalfEdge* getHalfEdge() const { return m_halfedge; }

	void computeBBox();
	void computeVTree(float sqdist);
	void computeHalfedge();

	// vertex and index
	TtArray<Vec3f> m_vdata;
	TtArray<Vec3i> m_trivid;

	// normal and index
	TtArray<Vec3f> m_ndata;
	TtArray<Vec3i> m_trinid;

	// texture coordinates and index
	TtArray<Vec2f> m_texdata;
	TtArray<Vec3i> m_tritexid;

	// ambient occlusion data
	TtArray<Vec4f> m_adata;

	// tid list shared a vertex
	std::vector< std::vector<int> > m_v_tidlist;

	float      m_smooth_angle;

private:
	PointerSet<TriMeshGroup> m_groups;
	BBox       m_bbox;
	VTree*     m_vtree;
	HalfEdge*  m_halfedge;

	TtFileName m_fname;
	int        m_axis;
	float      m_size;
	float      m_scale;

friend class TriMeshEditor;
};

void tm_computeFlatNormals(TriMesh* mesh);
void tm_finalizeMesh(TriMesh* mesh);

#endif  // tt_mesh_h

