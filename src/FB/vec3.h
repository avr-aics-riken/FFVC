#ifndef _VEC3_H_
#define _VEC3_H_

//##################################################################################
//
// vec3 class
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2014 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

/**
 * @file   Vec3.h
 * @brief  Vec3<T> class Header
 * @author aics
 */

//##################################################################################
//
// ATTENTION : If you want modify this header, please consult with organizer.
//             Since this header class is being used for several libraries,
//             we need to take care fo consistency between them.
//
//##################################################################################

#include <iostream>
#include <math.h>


namespace Vec3class {
  
typedef enum {
	AXIS_X = 0,
	AXIS_Y,
	AXIS_Z,
	AXIS_ERROR
} AxisEnum;

//=========================================================================
// class Vec3<T>
//=========================================================================

template<typename T>
class Vec3 {
public:
	Vec3(T v = 0)			{ t[0] = t[1] = t[2] = v; }
	Vec3(T _x, T _y, T _z)	{ t[0]=_x; t[1]=_y; t[2]=_z; }
	Vec3(const T v[3])		{ t[0] = v[0]; t[1] = v[1]; t[2] = v[2]; }

	Vec3<T>& assign(T _x, T _y, T _z) { 
		t[0]=_x; t[1]=_y; t[2]=_z; 
		return *this; 
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

	T lengthSquared() const { 
		return t[0] * t[0] + t[1] * t[1] + t[2] *t [2]; 
	}

	T length() const { return sqrt(lengthSquared()); }

	Vec3<T>& normalize() {
		T len = length();
		if (len != 0)
			return *this /= len;
		else
			return *this;
	}

	Vec3<T>& normalize(T* len) {
		*len = length();
		if (*len != 0)
			return *this /= *len;
		else
			return *this;
	}

	T average() const { return (t[0] + t[1] + t[2])/3.f; }
  
	T t[3];
  T x, y, z;
};

//=========================================================================
// typedef
//=========================================================================

typedef Vec3<unsigned char> Vec3uc;
typedef Vec3<int>           Vec3i;
typedef Vec3<float>         Vec3f;
typedef Vec3<double>        Vec3d;

//=========================================================================
// inline
//=========================================================================  

template <typename T>
inline Vec3<T> operator*(T s, const Vec3<T>& v) {
	return Vec3<T>(s*v.t[0], s*v.t[1], s*v.t[2]);
}
  
template <typename T>
inline Vec3<T> multi(const Vec3<T>& a, const Vec3<T>& b) {
	return a * b;
}
  
template <typename T>
inline T dot(const Vec3<T>& a, const Vec3<T>& b) {
	return a.t[0] * b.t[0] + a.t[1] * b.t[1] + a.t[2] * b.t[2];
}
  
template <typename T>
inline Vec3<T> cross(const Vec3<T>& a, const Vec3<T>& b) {
	return Vec3<T>(a.t[1] * b.t[2] - a.t[2] * b.t[1],
		a.t[2] * b.t[0] - a.t[0] * b.t[2],
		a.t[0] * b.t[1] - a.t[1] * b.t[0]);
}
  
template <typename T>
inline T distanceSquared(const Vec3<T>& a, const Vec3<T>& b) {
	return (a - b).lengthSquared();
}
  
template <typename T>
inline T distance(const Vec3<T>& a, const Vec3<T>& b) {
	return (a - b).length();
}

// @brief compare length between a and b, if a<b return true
inline bool lessVec3f(const Vec3f& a, const Vec3f& b) 
{
	return (a.lengthSquared() < b.lengthSquared()) ? true : false;
}

//=============================


template<typename T>
inline std::istream& operator>>(std::istream& is, Vec3<T>& v) 
{
	return is >> v[0] >> v[1] >> v[2];
}

template<typename T>
inline std::ostream& operator<<(std::ostream& os, const Vec3<T>& v) 
{
	return os << v[0] << " " << v[1] << " " << v[2];
}



inline std::istream& operator>>(std::istream& is, Vec3uc& v) 
{
	int x[3];
	is >> x[0] >> x[1] >> x[2];
	v[0]=x[0]; v[1]=x[1]; v[2]=x[2];
	return is;
}

inline std::ostream& operator<<(std::ostream& os, const Vec3uc& v) 
{
	int x[3];
	x[0]=v[0]; x[1]=v[1]; x[2]=v[2];
	os << x[0] << " " << x[1] << " " << x[2];
	return os;
}

} // namespace Vec3

#endif  // _VEC3_H_
