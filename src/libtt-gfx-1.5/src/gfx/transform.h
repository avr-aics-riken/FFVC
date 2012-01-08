#ifndef tt_transform_h
#define tt_transform_h

#include <tt/gm.h>

//=========================================================================
// class Transform
//=========================================================================

class Transform {
public:
	Transform(const Mat4x4f& M = Mat4x4f::identity());

	void assign(const Mat4x4f& M);
	void assign(const Mat4x4f& M, const Mat4x4f& N);

	Vec3f transformP(const Vec3f& o) const;
	Vec3f transformV(const Vec3f& o) const;
	Vec3f transformN(const Vec3f& o) const;
	Ray   transform (const Ray& o, float& len) const;
	BBox  transform (const BBox& o) const;

	static Mat4x4f getCenteringMatrix(const BBox& bbox);
	static Mat4x4f getResizingMatrix(const BBox& bbox, int axis, float size);
	static Mat4x4f getAligningMatrix(const BBox& bbox, int axis, int sign);

	Mat4x4f M, N;
};

//=========================================================================
// class Coordinates
//=========================================================================

class Coordinates {
public:
	Coordinates(const Vec3f& x, const Vec3f& y, const Vec3f& z);

	Vec3f w2o(const Vec3f& v) const;
	Vec3f o2w(const Vec3f& v) const;

	Vec3f X, Y, Z;
};

//=========================================================================
// class HPoint
//=========================================================================

class HPoint {
public:
	HPoint() : t(0), u(0), v(0), tid(0) {}

	Vec3f w2o(const Vec3f& v) const;
	Vec3f o2w(const Vec3f& v) const;

	float t, u, v;
	uint  tid;
	Vec3f P;
	Vec3f N;
	Vec3f S, T;
	Vec2f texcoord;
};

#endif  // tt_transform_h

