#ifndef _FB_vec3_h
#define _FB_vec3_h

// #################################################################
//
// CAERU Library
//
// Copyright (c) 2012 All right reserved.
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
//
// #################################################################

/**
 * @file   vec3.h
 * @brief  FlowBase Vec3<T> class Header
 * @author T. Tawara and kero
 */

#include <iostream>
#include <math.h>
#include "cpm_Define.h"
#include "vec3_func.h"
#include "vec3f_func.h"

namespace FB {
  
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

  /**
   * @brief add component: Vec3.x += Vec3.x
   */
	Vec3<T>& operator+=(const Vec3<T>& v) 
  {
		x+=v.x; y+=v.y; z+=v.z; 
		return *this;
	}
  
  /**
   * @brief subtract component: Vec3.x -= Vec3.x
   */
	Vec3<T>& operator-=(const Vec3<T>& v) 
  {
		x-=v.x; y-=v.y; z-=v.z; 
		return *this;
	}
  
  /**
   * @brief multiply component: Vec3.x *= Vec3.x
   */
	Vec3<T>& operator*=(const Vec3<T>& v) 
  {
		x*=v.x; y*=v.y; z*=v.z; 
		return *this;
	}
  
  /**
   * @brief divide by component: Vec3.x /= Vec3.x
   */
	Vec3<T>& operator/=(const Vec3<T>& v) 
  {
		x/=v.x; y/=v.y; z/=v.z; 
		return *this;
	}
  
  /**
   @brief multiply scalar: Vec3.x *= s
   */
	Vec3<T>& operator*=(T s) 
  {
		x*=s; y*=s; z*=s; 
		return *this;
	}
  
  /**
   * @brief divided by scalar: Vec3.x /= s
   */
	Vec3<T>& operator/=(T s) 
  {
		T inv = 1./s;
		x*=inv; y*=inv; z*=inv; 
		return *this;
	}

  /**
   * @fn Vec3<T> operator+(const Vec3<T>& v) const
   * @brief add: Vec3.x += v.x
   */
	Vec3<T> operator+(const Vec3<T>& v) const 
  {
		return Vec3<T>(x+v.x, y+v.y, z+v.z);
	}
  
  /**
   * @brief add: Vec3.x -= v.x
   */
	Vec3<T> operator-(const Vec3<T>& v) const 
  {
		return Vec3<T>(x-v.x, y-v.y, z-v.z);
	}
  
  /**
   * @brief multiply: Vec3.x *= v.x
   */
	Vec3<T> operator*(const Vec3<T>& v) const 
  {
		return Vec3<T>(x*v.x, y*v.y, z*v.z);
	}
  
  /**
   * @brief divide: Vec3.x /= v.x
   */
	Vec3<T> operator/(const Vec3<T>& v) const 
  {
		return Vec3<T>(x/v.x, y/v.y, z/v.z);
	}
  
  /**
   * @brief multiply scalar: Vec3.x *= s
   */
	Vec3<T> operator*(T s) const 
  {
		return Vec3<T>(x*s, y*s, z*s);
	}
  
  /**
   * @brief divided by scalar: Vec3.x /= s
   */
	Vec3<T> operator/(T s) const 
  {
		T inv = 1./s;
		return Vec3<T>(x*inv, y*inv, z*inv);
	}
  
  /**
   * @brief minus: -Vec3.x
   */
	Vec3<T> operator-() const 
  {
		return Vec3<T>(-x, -y, -z);
	}

  /**
   * @brief return true iff equivalent
   */
	bool operator==(const Vec3<T>& v) const 
  {
		return x == v.x && y == v.y && z == v.z;
	}
  
  /**
   * @brief return true iff not equivalent
   */
	bool operator!=(const Vec3<T>& v) const 
  {
		return !(*this == v);
	}

	static Vec3<T> xaxis() { return Vec3<T>(1, 0, 0); }
	static Vec3<T> yaxis() { return Vec3<T>(0, 1, 0); }
	static Vec3<T> zaxis() { return Vec3<T>(0, 0, 1); }

  /**
   * @brief return squared length of Vec
   */
	float lengthSquared() const 
  { 
    return x*x + y*y + z*z; 
  }
  
  /**
   * @brief return length of Vec
   */
	float length() const 
  { 
    return sqrtf(lengthSquared()); 
  }
  
  /**
   * @brief normalized by length
   */
	Vec3<T>& normalize() 
  {
		float len = length();
		if (len != 0)
			return *this /= len;
		else
			return *this;
	}
  
  /**
   * @brief normalized by len
   */
	Vec3<T>& normalize(float* len) 
  {
		*len = length();
		if (*len != 0)
			return *this /= *len;
		else
			return *this;
	}
  
  /**
   * @brief return averaged value
   */
	float average() const 
  { 
    return (x + y + z)/3.f; 
  }

	T x, y, z;
};

//=========================================================================
// typedef
//=========================================================================

typedef Vec3<unsigned char> Vec3uc;
typedef Vec3<int> Vec3i;
typedef Vec3<float> Vec3f;
typedef Vec3<REAL_TYPE> Vec3r;

//=========================================================================
  
/**
 * @brief multiply s
 */
inline Vec3r operator*(double s, const Vec3r& v)
{
  return Vec3r(s*v.x, s*v.y, s*v.z);
}
  
/**
 * @brief multiply s
 */
inline Vec3f operator*(float s, const Vec3f& v) 
{
	return Vec3f(s*v.x, s*v.y, s*v.z);
}

/**
 * @brief multiply
 */
inline Vec3f multi(const Vec3f& a, const Vec3f& b) 
{
	return a * b;
}

/**
 * @brief inner product
 */
inline float dot(const Vec3f& a, const Vec3f& b) 
{
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

/**
 * @brief cross product
 */
inline Vec3f cross(const Vec3f& a, const Vec3f& b) 
{
	return Vec3f(a.y * b.z - a.z * b.y,
	             a.z * b.x - a.x * b.z,
	             a.x * b.y - a.y * b.x);
}

/**
 * @brief squared distance
 */
inline float distanceSquared(const Vec3f& a, const Vec3f& b) 
{
	return (a - b).lengthSquared();
}

/**
 * @brief distance
 */
inline float distance(const Vec3f& a, const Vec3f& b) 
{
	return (a - b).length();
}


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

/**
 * @brief compare length between a and b, if a<b return true
 */
inline bool lessVec3f(const Vec3f& a, const Vec3f& b) 
{
	return (a.lengthSquared() < b.lengthSquared()) ? true : false;
}
  
} // namespace FB

#endif  // _FB_vec3_h
