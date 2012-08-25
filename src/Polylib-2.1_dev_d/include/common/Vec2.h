/*
 * Polylib - Polygon Management Library.
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2010-
 *
 */
#ifndef vec2_h
#define vec2_h

#include <iostream>
#include <math.h>

namespace PolylibNS {

////////////////////////////////////////////////////////////////////////////
///
/// クラス:Vec2<T>
///
////////////////////////////////////////////////////////////////////////////
template<typename T>
class Vec2 {
public:
	Vec2(T v = 0)		{ x = y = v; }
	Vec2(T _x, T _y)	{ x = _x; y = _y; }
	Vec2(const T v[2])	{ x = v[0]; y = v[1]; }
	// use default copy constructor
	// use default operator =()

	Vec2<T>& assign(T _x, T _y) { x = _x; y = _y; return *this; }

	operator       T*()       { return &x; }
	operator const T*() const { return &x; }
	      T* ptr()       { return &x; }
	const T* ptr() const { return &x; }
	      T& operator [](int i)       { return (&x)[i]; }
	const T& operator [](int i) const { return (&x)[i]; }

	Vec2<T>& operator+=(const Vec2<T>& v) {
		x += v.x; y += v.y;
		return *this;
	}
	Vec2<T>& operator-=(const Vec2<T>& v) {
		x -= v.x; y -= v.y;
		return *this;
	}
	Vec2<T>& operator*=(const Vec2<T>& v) {
		x *= v.x; y *= v.y;
		return *this;
	}
	Vec2<T>& operator/=(const Vec2<T>& v) {
		x /= v.x; y /= v.y;
		return *this;
	}
	Vec2<T>& operator*=(T s) {
		x *= s; y *= s;
		return *this;
	}
	Vec2<T>& operator/=(T s) {
		T inv = 1./s;
		x *= inv; y *= inv;
		return *this;
	}

	Vec2<T> operator+(const Vec2<T>& v) const {
		return Vec2<T>(x + v.x, y + v.y);
	}
	Vec2<T> operator-(const Vec2<T>& v) const {
		return Vec2<T>(x - v.x, y - v.y);
	}
	Vec2<T> operator*(const Vec2<T>& v) const {
		return Vec2<T>(x * v.x, y * v.y);
	}
	Vec2<T> operator/(const Vec2<T>& v) const {
		return Vec2<T>(x / v.x, y / v.y);
	}
	Vec2<T> operator*(T s) const {
		return Vec2<T>(x * s, y * s);
	}
	Vec2<T> operator/(T s) const {
		T inv = 1./s;
		return Vec2<T>(x * inv, y * inv);
	}
	Vec2<T> operator-() const {
		return Vec2<T>(-x, -y);
	}

	bool operator==(const Vec2<T>& v) const {
		return x == v.x && y == v.y;
	}
	bool operator!=(const Vec2<T>& v) const {
		return !(*this == v);
	}

	static Vec2<T> xaxis() { return Vec2<T>(1, 0); }
	static Vec2<T> yaxis() { return Vec2<T>(0, 1); }

	float lengthSquared() const { return x*x + y*y; }
	float length() const { return sqrtf(lengthSquared()); }
	Vec2<T>& normalize() {
		float len = length();
		if (len != 0) {
			return *this /= len;
		}
		else {
			return *this;
		}
	}
	Vec2<T>& normalize(float* len) {
		*len = length();
		if (*len != 0) {
			return *this /= *len;
		}
		else {
			return *this;
		}
	}
	float average() const { return (x + y)/2.f; }

	T x, y;
};

//=========================================================================
// typedef
//=========================================================================
typedef Vec2<unsigned char>	Vec2uc;
typedef Vec2<int>			Vec2i;
typedef Vec2<float>			Vec2f;

//=========================================================================
// inline
//=========================================================================
inline Vec2f operator*(float s, const Vec2f& v) {
	return Vec2f(s*v.x, s*v.y);
}

inline float distanceSquared(const Vec2f& a, const Vec2f& b) {
	return (a - b).lengthSquared();
}
inline float distance(const Vec2f& a, const Vec2f& b) {
	return (a - b).length();
}

template<typename T>
inline std::istream& operator>>(std::istream& is, Vec2<T>& v) {
	return is >> v[0] >> v[1];
}
template<typename T>
inline std::ostream& operator<<(std::ostream& os, const Vec2<T>& v) {
	return os << v[0] << " " << v[1];
}
inline std::istream& operator>>(std::istream& is, Vec2uc& v) {
	int x[2];
	is >> x[0] >> x[1];
	v[0] = x[0]; v[1] = x[1];
	return is;
}
inline std::ostream& operator<<(std::ostream& os, const Vec2uc& v) {
	int x[2];
	x[0] = v[0]; x[1] = v[1];
	os << x[0] << " " << x[1];
	return os;
}

inline bool lessVec2f(const Vec2f& a, const Vec2f& b) {
	return (a.lengthSquared() < b.lengthSquared()) ? true : false;
}

} //namespace PolylibNS

#endif  // vec2_h

