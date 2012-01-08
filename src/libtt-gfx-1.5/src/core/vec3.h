#ifndef vec3_h
#define vec3_h

#include <iostream>
#include <math.h>
#include "vec3_func.h"
#include "vec3f_func.h"

//=========================================================================
// class Vec3<T>
//=========================================================================

template<typename T>
class Vec3
{
public:
	Vec3(T v=0) { x=y=z=v; }
	Vec3(T _x, T _y, T _z) { x=_x; y=_y; z=_z; }
	Vec3(const T v[3]) { x=v[0]; y=v[1]; z=v[2]; }
	// use default copy constructor
	// use default operator =()

	Vec3<T>& assign(T _x, T _y, T _z) { x=_x; y=_y; z=_z; return *this; }

	operator       T*()       { return &x; }
	operator const T*() const { return &x; }
	      T* ptr()       { return &x; }
	const T* ptr() const { return &x; }
	      T& operator [](int i)       { return (&x)[i]; }
	const T& operator [](int i) const { return (&x)[i]; }

	Vec3<T>& operator+=(const Vec3<T>& v) {
		x+=v.x; y+=v.y; z+=v.z; 
		return *this;
	}
	Vec3<T>& operator-=(const Vec3<T>& v) {
		x-=v.x; y-=v.y; z-=v.z; 
		return *this;
	}
	Vec3<T>& operator*=(const Vec3<T>& v) {
		x*=v.x; y*=v.y; z*=v.z; 
		return *this;
	}
	Vec3<T>& operator/=(const Vec3<T>& v) {
		x/=v.x; y/=v.y; z/=v.z; 
		return *this;
	}
	Vec3<T>& operator*=(T s) {
		x*=s; y*=s; z*=s; 
		return *this;
	}
	Vec3<T>& operator/=(T s) {
		T inv = 1./s;
		x*=inv; y*=inv; z*=inv; 
		return *this;
	}

	Vec3<T> operator+(const Vec3<T>& v) const {
		return Vec3<T>(x+v.x, y+v.y, z+v.z);
	}
	Vec3<T> operator-(const Vec3<T>& v) const {
		return Vec3<T>(x-v.x, y-v.y, z-v.z);
	}
	Vec3<T> operator*(const Vec3<T>& v) const {
		return Vec3<T>(x*v.x, y*v.y, z*v.z);
	}
	Vec3<T> operator/(const Vec3<T>& v) const {
		return Vec3<T>(x/v.x, y/v.y, z/v.z);
	}
	Vec3<T> operator*(T s) const {
		return Vec3<T>(x*s, y*s, z*s);
	}
	Vec3<T> operator/(T s) const {
		T inv = 1./s;
		return Vec3<T>(x*inv, y*inv, z*inv);
	}
	Vec3<T> operator-() const {
		return Vec3<T>(-x, -y, -z);
	}

	bool operator==(const Vec3<T>& v) const {
		return x == v.x && y == v.y && z == v.z;
	}
	bool operator!=(const Vec3<T>& v) const {
		return !(*this == v);
	}

	static Vec3<T> xaxis() { return Vec3<T>(1, 0, 0); }
	static Vec3<T> yaxis() { return Vec3<T>(0, 1, 0); }
	static Vec3<T> zaxis() { return Vec3<T>(0, 0, 1); }

	float lengthSquared() const { return x*x + y*y + z*z; }
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
	float average() const { return (x + y + z)/3.f; }

	T x, y, z;
};

//=========================================================================
// typedef
//=========================================================================

typedef Vec3<unsigned char> Vec3uc;
typedef Vec3<int> Vec3i;
typedef Vec3<float> Vec3f;

//=========================================================================

inline Vec3f operator*(float s, const Vec3f& v) {
	return Vec3f(s*v.x, s*v.y, s*v.z);
}
inline Vec3f multi(const Vec3f& a, const Vec3f& b) {
	return a * b;
}

inline float dot(const Vec3f& a, const Vec3f& b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}
inline Vec3f cross(const Vec3f& a, const Vec3f& b) {
	return Vec3f(a.y * b.z - a.z * b.y,
	             a.z * b.x - a.x * b.z,
	             a.x * b.y - a.y * b.x);
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
	v[0]=x[0]; v[1]=x[1]; v[2]=x[2];
	return is;
}
inline std::ostream& operator<<(std::ostream& os, const Vec3uc& v) {
	int x[3];
	x[0]=v[0]; x[1]=v[1]; x[2]=v[2];
	os << x[0] << " " << x[1] << " " << x[2];
	return os;
}

inline bool lessVec3f(const Vec3f& a, const Vec3f& b) {
	return (a.lengthSquared() < b.lengthSquared()) ? true : false;
}

#endif  // vec3_h

