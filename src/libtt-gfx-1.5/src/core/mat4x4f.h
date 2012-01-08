#ifndef mat4x4f_h
#define mat4x4f_h

#include <assert.h>
#include <iostream>
#include "vec3.h"
#include "vec4.h"

class Mat4x4f {
public:
	Mat4x4f();
	Mat4x4f(float, float, float, float,
	float, float, float, float,
	float, float, float, float,
	float, float, float, float);
	// use default copy constructor
	// use default operator =()

	      float* operator [](int i)       { return &m_M[i][0]; }
	const float* operator [](int i) const { return &m_M[i][0]; }

	Mat4x4f& assign(float, float, float, float,
	float, float, float, float,
	float, float, float, float,
	float, float, float, float);

	Mat4x4f& operator *=(const Mat4x4f& o);
	Mat4x4f& operator *=(float c);
	Mat4x4f& operator /=(float c);

	Mat4x4f  operator * (const Mat4x4f& o) const;
	Mat4x4f  operator * (float c) const;
	Mat4x4f  operator / (float c) const;
	Mat4x4f  operator - () const;
	Vec4f    operator * (const Vec4f& v) const;

friend Mat4x4f operator * (float c, const Mat4x4f& o);
friend Vec4f   operator * (const Vec4f& v, const Mat4x4f& o);

	bool operator ==(const Mat4x4f& o) const;
	bool operator !=(const Mat4x4f& o) const;

	Mat4x4f inverse() const;
	Mat4x4f transpose() const;
	Mat4x4f adjoint() const;
	float   determinant() const;
	bool    isSingular() const;

	// transform point
	Vec3f transform(const Vec3f& v) const;
	// transform vector
	Vec3f transformVector(const Vec3f& v) const;

	// transformation matrices
	static Mat4x4f identity();
	static Mat4x4f translate(float x, float y, float z);
	static Mat4x4f rotate(float angle, float x, float y, float z);
	static Mat4x4f scale(float x, float y, float z);
	static Mat4x4f translate(const Vec3f& v);
	static Mat4x4f rotate(float angle, const Vec3f& axis);
	static Mat4x4f scale(const Vec3f& v);

	// cubic basis matrices
	static Mat4x4f hermiteBasis();
	static Mat4x4f bezierBasis();
	static Mat4x4f bsplineBasis();
	static Mat4x4f catmullromBasis();

private:
	float m_M[4][4];
};

std::istream& operator >>(std::istream& is, Mat4x4f& M);
std::ostream& operator <<(std::ostream& os, const Mat4x4f& M);
std::ostream& fancy_print(std::ostream& os, const Mat4x4f& M);

#endif  // mat4x4f_h

