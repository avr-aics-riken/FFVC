#ifndef vec4_h
#define vec4_h

#include <iostream>
#include <math.h>
#include "vec4_func.h"

//=========================================================================
// class Vec4<T>
//=========================================================================

template<typename T>
class Vec4
{
public:
	Vec4(T v=0) { x=y=z=w=v; }
	Vec4(T _x, T _y, T _z, T _w) { x=_x; y=_y; z=_z; w=_w; }
	Vec4(const T v[4]) { x=v[0]; y=v[1]; z=v[2]; w=v[3]; }
	// use default copy constructor
	// use default operator =()

	Vec4<T>& assign(T _x, T _y, T _z, T _w) { x=_x; y=_y; z=_z; w=_w; return *this; }

	operator       T*()       { return &x; }
	operator const T*() const { return &x; }
	      T* ptr()       { return &x; }
	const T* ptr() const { return &x; }
	      T& operator [](int i)       { return (&x)[i]; }
	const T& operator [](int i) const { return (&x)[i]; }

	Vec4<T>& operator+=(const Vec4<T>& v) {
		x+=v.x; y+=v.y; z+=v.z; w+=v.w;
		return *this;
	}
	Vec4<T>& operator-=(const Vec4<T>& v) {
		x-=v.x; y-=v.y; z-=v.z; w-=v.w;
		return *this;
	}
	Vec4<T>& operator*=(const Vec4<T>& v) {
		x*=v.x; y*=v.y; z*=v.z; w*=v.w;
		return *this;
	}
	Vec4<T>& operator/=(const Vec4<T>& v) {
		x/=v.x; y/=v.y; z/=v.z; w/=v.w;
		return *this;
	}
	Vec4<T>& operator*=(T s) {
		x*=s; y*=s; z*=s; w*=s;
		return *this;
	}
	Vec4<T>& operator/=(T s) {
		T inv = 1./s;
		x*=inv; y*=inv; z*=inv; w*=inv;
		return *this;
	}

	Vec4<T> operator+(const Vec4<T>& v) const {
		return Vec4<T>(x+v.x, y+v.y, z+v.z, w+v.w);
	}
	Vec4<T> operator-(const Vec4<T>& v) const {
		return Vec4<T>(x-v.x, y-v.y, z-v.z, w-v.w);
	}
	Vec4<T> operator*(const Vec4<T>& v) const {
		return Vec4<T>(x*v.x, y*v.y, z*v.z, w*v.w);
	}
	Vec4<T> operator/(const Vec4<T>& v) const {
		return Vec4<T>(x/v.x, y/v.y, z/v.z, w/v.w);
	}
	Vec4<T> operator*(T s) const {
		return Vec4<T>(x*s, y*s, z*s, w*s);
	}
	Vec4<T> operator/(T s) const {
		T inv = 1./s;
		return Vec4<T>(x*inv, y*inv, z*inv, w*inv);
	}
	Vec4<T> operator-() const {
		return Vec4<T>(-x, -y, -z, -w);
	}

	bool operator==(const Vec4<T>& v) const {
		return x == v.x && y == v.y && z == v.z && w == v.w;
	}
	bool operator!=(const Vec4<T>& v) const {
		return !(*this == v);
	}

	static Vec4<T> xaxis() { return Vec4<T>(1, 0, 0, 0); }
	static Vec4<T> yaxis() { return Vec4<T>(0, 1, 0, 0); }
	static Vec4<T> zaxis() { return Vec4<T>(0, 0, 1, 0); }
	static Vec4<T> waxis() { return Vec4<T>(0, 0, 0, 1); }

	float average() const { return (x + y + z + w)/4.f; }

	T x, y, z, w;
};

//=========================================================================
// typedef
//=========================================================================

typedef Vec4<unsigned char> Vec4uc;
typedef Vec4<int> Vec4i;
typedef Vec4<float> Vec4f;

//=========================================================================

inline Vec4f operator*(float s, const Vec4f& v) {
	return Vec4f(s*v.x, s*v.y, s*v.z, s*v.w);
}

template<typename T>
inline std::istream& operator>>(std::istream& is, Vec4<T>& v) {
	return is >> v[0] >> v[1] >> v[2] >> v[3];
}
template<typename T>
inline std::ostream& operator<<(std::ostream& os, const Vec4<T>& v) {
	return os << v[0] << " " << v[1] << " " << v[2] << " " << v[3];
}
inline std::istream& operator>>(std::istream& is, Vec4uc& v) {
	int x[4];
	is >> x[0] >> x[1] >> x[2] >> x[3];
	v[0]=x[0]; v[1]=x[1]; v[2]=x[2]; v[3]=x[3];
	return is;
}
inline std::ostream& operator<<(std::ostream& os, const Vec4uc& v) {
	int x[4];
	x[0]=v[0]; x[1]=v[1]; x[2]=v[2]; x[3]=v[3];
	os << x[0] << " " << x[1] << " " << x[2] << " " << x[3];
	return os;
}

#endif  // vec4_h

