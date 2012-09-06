/*
 * Polylib - Polygon Management Library.
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2010-
 *
 */
#ifndef vec3f_func_h
#define vec3f_func_h

namespace PolylibNS {

inline
void vec3f_copy(float to[3], const float from[3])
{
	to[0] = from[0];
	to[1] = from[1];
	to[2] = from[2];
}

inline
void vec3f_set(float v[3], float x, float y, float z)
{
	v[0] = x;
	v[1] = y;
	v[2] = z;
}

inline
void vec3f_min(float c[3], const float a[3], const float b[3])
{
	c[0] = (a[0] < b[0]) ? a[0] : b[0];
	c[1] = (a[1] < b[1]) ? a[1] : b[1];
	c[2] = (a[2] < b[2]) ? a[2] : b[2];
}

inline
void vec3f_max(float c[3], const float a[3], const float b[3])
{
	c[0] = (a[0] > b[0]) ? a[0] : b[0];
	c[1] = (a[1] > b[1]) ? a[1] : b[1];
	c[2] = (a[2] > b[2]) ? a[2] : b[2];
}

inline
void vec3f_plus(float c[3], const float a[3], const float b[3])
{
	c[0] = a[0] + b[0];
	c[1] = a[1] + b[1];
	c[2] = a[2] + b[2];
}

inline
void vec3f_minus(float c[3], const float a[3], const float b[3])
{
	c[0] = a[0] - b[0];
	c[1] = a[1] - b[1];
	c[2] = a[2] - b[2];
}

inline
void vec3f_multi(float c[3], const float a[3], const float b[3])
{
	c[0] = a[0] * b[0];
	c[1] = a[1] * b[1];
	c[2] = a[2] * b[2];
}

inline
void vec3f_multi(float c[3], const float a[3], const float b)
{
	c[0] = a[0] * b;
	c[1] = a[1] * b;
	c[2] = a[2] * b;
}

inline
void vec3f_div(float c[3], const float a[3], const float b[3])
{
	c[0] = a[0] / b[0];
	c[1] = a[1] / b[1];
	c[2] = a[2] / b[2];
}

inline
void vec3f_div(float c[3], const float a[3], const float b)
{
	float inv = 1.f / b;
	c[0] = a[0] * inv;
	c[1] = a[1] * inv;
	c[2] = a[2] * inv;
}

inline
float vec3f_dot(const float a[3], const float b[3])
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

inline
void vec3f_cross(float c[3], const float a[3], const float b[3])
{
	c[0] = a[1] * b[2] - a[2] * b[1];
	c[1] = a[2] * b[0] - a[0] * b[2];
	c[2] = a[0] * b[1] - a[1] * b[0];
}

inline
float vec3f_sqdist(const float a[3], const float b[3])
{
	float x = a[0] - b[0];
	float y = a[1] - b[1];
	float z = a[2] - b[2];
	return x*x + y*y + z*z;
}

inline
float vec3f_dist(const float a[3], const float b[3])
{
	return sqrtf(vec3f_sqdist(a, b));
}

} //namespace PolylibNS

#endif  // vec3f_func_h

