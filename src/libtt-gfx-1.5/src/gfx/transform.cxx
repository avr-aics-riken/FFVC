#include "transform.h"

//=========================================================================
// class Transform
//=========================================================================

Transform::Transform(const Mat4x4f& _M) : M(_M) {
	N = M.inverse();
}

void
Transform::assign(const Mat4x4f& _M) {
	M = _M;
	N = M.inverse();
}

void
Transform::assign(const Mat4x4f& _M, const Mat4x4f& _n) {
	M = _M;
	N = _n;
}

Vec3f
Transform::transformP(const Vec3f& o) const {
	return M.transform(o);
}

Vec3f
Transform::transformV(const Vec3f& o) const {
	return M.transformVector(o);
}

Vec3f
Transform::transformN(const Vec3f& o) const {
	return Vec3f(
	N[0][0]*o[0] + N[1][0]*o[1] + N[2][0]*o[2],
	N[0][1]*o[0] + N[1][1]*o[1] + N[2][1]*o[2],
	N[0][2]*o[0] + N[1][2]*o[1] + N[2][2]*o[2]).normalize();
}

Ray
Transform::transform(const Ray& o, float& len) const {
	Vec3f org = this->transformP(o.org);
	Vec3f dir = this->transformV(o.dir);
	len = dir.length();
	assert(!gmIsZero(len));
	Ray ray(org, dir);
	ray.min_t = o.min_t * len;
	ray.max_t = o.max_t * len;
	return ray;
}

BBox
Transform::transform(const BBox& o) const {
	const Vec3f& p = o.min;
	const Vec3f& q = o.max;
	BBox bbox;
	bbox.add(this->transformP(Vec3f(p[0], p[1], p[2])));
	bbox.add(this->transformP(Vec3f(q[0], p[1], p[2])));
	bbox.add(this->transformP(Vec3f(q[0], q[1], p[2])));
	bbox.add(this->transformP(Vec3f(p[0], q[1], p[2])));

	bbox.add(this->transformP(Vec3f(p[0], p[1], q[2])));
	bbox.add(this->transformP(Vec3f(q[0], p[1], q[2])));
	bbox.add(this->transformP(Vec3f(q[0], q[1], q[2])));
	bbox.add(this->transformP(Vec3f(p[0], q[1], q[2])));
	return bbox;
}

Mat4x4f
Transform::getCenteringMatrix(const BBox& bbox)
{
	Vec3f center = bbox.center();
	return Mat4x4f::translate(-center);
}

Mat4x4f
Transform::getResizingMatrix(const BBox& bbox, int axis, float size)
{
	float scale  = size / (bbox.max[axis] - bbox.min[axis]);
	return Mat4x4f::scale(scale, scale, scale);
}

Mat4x4f
Transform::getAligningMatrix(const BBox& bbox, int axis, int sign)
{
	Vec3f vec(0, 0, 0);

	if (sign > 0)
		vec[axis] = bbox.max[axis];
	else
		vec[axis] = bbox.min[axis];

	return Mat4x4f::translate(-vec);
}

//=========================================================================
// class Coordinates
//=========================================================================

Coordinates::Coordinates(const Vec3f& x, const Vec3f& y, const Vec3f& z)
: X(x), Y(y), Z(z) {
}

Vec3f
Coordinates::w2o(const Vec3f& v) const {
	return Vec3f(dot(X, v), dot(Y, v), dot(Z, v));
}

Vec3f
Coordinates::o2w(const Vec3f& v) const {
	return Vec3f(
		X[0] * v[0] + Y[0] * v[1] + Z[0] * v[2],
		X[1] * v[0] + Y[1] * v[1] + Z[1] * v[2],
		X[2] * v[0] + Y[2] * v[1] + Z[2] * v[2]);
}

//=========================================================================
// class HPoint
//=========================================================================

Vec3f
HPoint::w2o(const Vec3f& v) const {
	return Coordinates(S, T, N).w2o(v);
}

Vec3f
HPoint::o2w(const Vec3f& v) const {
	return Coordinates(S, T, N).o2w(v);
}

