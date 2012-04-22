/*
 * Polylib - Polygon Management Library.
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2010-
 *
 */
#ifndef vec3_h
#define vec3_h

#include <iostream>
#include <math.h>
#include "common/vec3_func.h"
#include "common/vec3f_func.h"
#include "common/axis.h"

namespace PolylibNS {

////////////////////////////////////////////////////////////////////////////
///
/// クラス:Vec3<T>
///
////////////////////////////////////////////////////////////////////////////
template<typename T>
class Vec3 {
public:
	Vec3(T v = 0)			{ t[0] = t[1] = t[2] = v; }
	Vec3(T _x, T _y, T _z)	{ t[0]=_x; t[1]=_y; t[2]=_z; }
	Vec3(const T v[3])		{ t[0] = v[0]; t[1] = v[1]; t[2] = v[2]; }
	// use default copt[1] constructor
	// use default operator =()

	Vec3<T>& assign(T _x, T _y, T _z) { 
		t[0] = _x; t[1] = _y; t[2] = _z; return *this; 
	}
 
	operator       T*()       { return &t[0]; }
	operator const T*() const { return &t[0]; }
	      T* ptr()       { return &t[0]; }
	const T* ptr() const { return &t[0]; }
	      T& operator [](const AxisEnum& axis) { 
              return t[axis];
          }
	const T& operator [](const AxisEnum& axis) const {
              return t[axis];
          }
//	      T& operator [](int i) { 
//              if (i == 0) return t[0];
//              else if (i == 1) return t[1];
//              return t[2];
//          }
//	const T& operator [](int i) const {
//              if (i == 0) return t[0];
//              else if (i == 1) return t[1];
//              return t[2];
//          }

	Vec3<T>& operator+=(const Vec3<T>& v) {
		t[0] += v.t[0]; t[1] += v.t[1]; t[2] += v.t[2]; 
		return *this;
	}
	Vec3<T>& operator-=(const Vec3<T>& v) {
		t[0] -= v.t[0]; t[1] -= v.t[1]; t[2] -= v.t[2]; 
		return *this;
	}
	Vec3<T>& operator*=(const Vec3<T>& v) {
		t[0] *= v.t[0]; t[1] *= v.t[1]; t[2] *= v.t[2]; 
		return *this;
	}
	Vec3<T>& operator/=(const Vec3<T>& v) {
		t[0] /= v.t[0]; t[1]/= v.t[1]; t[2] /= v.t[2]; 
		return *this;
	}
	Vec3<T>& operator*=(T s) {
		t[0] *= s; t[1] *= s; t[2] *= s; 
		return *this;
	}
	Vec3<T>& operator/=(T s) {
		T inv = 1./s;
		t[0] *= inv; t[1] *= inv; t[2] *= inv; 
		return *this;
	}

	Vec3<T> operator+(const Vec3<T>& v) const {
		return Vec3<T>(t[0] + v.t[0], t[1] + v.t[1], t[2] + v.t[2]);
	}
	Vec3<T> operator-(const Vec3<T>& v) const {
		return Vec3<T>(t[0] - v.t[0], t[1] - v.t[1], t[2] - v.t[2]);
	}
	Vec3<T> operator*(const Vec3<T>& v) const {
		return Vec3<T>(t[0] * v.t[0], t[1] * v.t[1], t[2] * v.t[2]);
	}
	Vec3<T> operator/(const Vec3<T>& v) const {
		return Vec3<T>(t[0] / v.t[0], t[1] / v.t[1], t[2] / v.t[2]);
	}
	Vec3<T> operator*(T s) const {
		return Vec3<T>(t[0] * s, t[1] * s, t[2] * s);
	}
	Vec3<T> operator/(T s) const {
		T inv = 1./s;
		return Vec3<T>(t[0] * inv, t[1] * inv, t[2] * inv);
	}
	Vec3<T> operator-() const {
		return Vec3<T>(-t[0], -t[1], -t[2]);
	}

	bool operator==(const Vec3<T>& v) const {
		return t[0] == v.t[0] && t[1] == v.t[1] && t[2] == v.t[2];
	}
	bool operator!=(const Vec3<T>& v) const {
		return !(*this == v);
	}

	static Vec3<T> xaxis() { return Vec3<T>(1, 0, 0); }
	static Vec3<T> yaxis() { return Vec3<T>(0, 1, 0); }
	static Vec3<T> zaxis() { return Vec3<T>(0, 0, 1); }

	float lengthSquared() const { return t[0] * t[0] + t[1] * t[1] + t[2] *t [2]; }
	float length() const { return sqrtf(lengthSquared()); }
	Vec3<T>& normalize() {
		float len = length();
		if (len != 0)
			return *this /= len;
		else
			return *this;
	}
	Vec3<T>& normalize(float* len) {
		*len = length();
		if (*len != 0)
			return *this /= *len;
		else
			return *this;
	}
	float average() const { return (t[0] + t[1] + t[2])/3.f; }

	T t[3];
};

//=========================================================================
// typedef
//=========================================================================
typedef Vec3<unsigned char>	Vec3uc;
typedef Vec3<int>			Vec3i;
typedef Vec3<float>			Vec3f;

//=========================================================================
// inline
//=========================================================================
inline Vec3f operator*(float s, const Vec3f& v) {
	return Vec3f(s*v.t[0], s*v.t[1], s*v.t[2]);
}
inline Vec3f multi(const Vec3f& a, const Vec3f& b) {
	return a * b;
}

inline float dot(const Vec3f& a, const Vec3f& b) {
	return a.t[0] * b.t[0] + a.t[1] * b.t[1] + a.t[2] * b.t[2];
}
inline Vec3f cross(const Vec3f& a, const Vec3f& b) {
	return Vec3f(a.t[1] * b.t[2] - a.t[2] * b.t[1],
	             a.t[2] * b.t[0] - a.t[0] * b.t[2],
	             a.t[0] * b.t[1] - a.t[1] * b.t[0]);
}
inline float distanceSquared(const Vec3f& a, const Vec3f& b) {
	return (a - b).lengthSquared();
}
inline float distance(const Vec3f& a, const Vec3f& b) {
	return (a - b).length();
}

template<typename T>
inline std::istream& operator>>(std::istream& is, Vec3<T>& v) {
	return is >> v[0] >> v[1] >> v[2];
}
template<typename T>
inline std::ostream& operator<<(std::ostream& os, const Vec3<T>& v) {
	return os << v[0] << " " << v[1] << " " << v[2];
}
inline std::istream& operator>>(std::istream& is, Vec3uc& v) {
	int x[3];
	is >> x[0] >> x[1] >> x[2];
	v[0] = x[0]; v[1] = x[1]; v[2] = x[2];
	return is;
}
inline std::ostream& operator<<(std::ostream& os, const Vec3uc& v) {
	int x[3];
	x[0] = v[0]; x[1] = v[1]; x[2] = v[2];
	os << x[0] << " " << x[1] << " " << x[2];
	return os;
}

inline bool lessVec3f(const Vec3f& a, const Vec3f& b) {
	return (a.lengthSquared() < b.lengthSquared()) ? true : false;
}

} //namespace PolylibNS

#endif  // vec3_h

