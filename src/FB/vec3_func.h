#ifndef _FB_vec3_func_h
#define _FB_vec3_func_h

// #################################################################
//
// CAERU Library
//
// Copyright (c) All right reserved. 2012
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
//
// #################################################################

//@file vec3_func.h
//@brief FlowBase class geometry functions
//@author T. Tawara and kero

namespace FB {
  
/**
 @fn template<typename T1, typename T2> inline void vec3_copy(T1 to[3], const T2 from[3])
 @brief copy data
 */
template<typename T1, typename T2> inline
void vec3_copy(T1 to[3], const T2 from[3])
{
	to[0] = from[0];
	to[1] = from[1];
	to[2] = from[2];
}

/**
 @fn template<typename T1, typename T2> inline void vec3_set(T1 v[3], T2 x, T2 y, T2 z)
 @brief set data
 */
template<typename T1, typename T2> inline
void vec3_set(T1 v[3], T2 x, T2 y, T2 z)
{
	v[0] = x;
	v[1] = y;
	v[2] = z;
}

/**
 @fn template<typename T> inline void vec3_min(T c[3], const T a[3], const T b[3])
 @brief get minimum
 */
template<typename T> inline
void vec3_min(T c[3], const T a[3], const T b[3])
{
	c[0] = (a[0] < b[0]) ? a[0] : b[0];
	c[1] = (a[1] < b[1]) ? a[1] : b[1];
	c[2] = (a[2] < b[2]) ? a[2] : b[2];
}

/**
 @fn template<typename T> inline void vec3_max(T c[3], const T a[3], const T b[3])
 @brief get maximum
 */
template<typename T> inline
void vec3_max(T c[3], const T a[3], const T b[3])
{
	c[0] = (a[0] > b[0]) ? a[0] : b[0];
	c[1] = (a[1] > b[1]) ? a[1] : b[1];
	c[2] = (a[2] > b[2]) ? a[2] : b[2];
}

/**
 @fn template<typename T1, typename T2, typename T3> inline void vec3_plus(T1 c[3], const T2 a[3], const T3 b[3])
 @brief c[] = a[] + b[]
 */
template<typename T1, typename T2, typename T3> inline
void vec3_plus(T1 c[3], const T2 a[3], const T3 b[3])
{
	c[0] = a[0] + b[0];
	c[1] = a[1] + b[1];
	c[2] = a[2] + b[2];
}

/**
 @fn template<typename T1, typename T2, typename T3> inline void vec3_minus(T1 c[3], const T2 a[3], const T3 b[3])
 @brief c[] = a[] - b[]
 */
template<typename T1, typename T2, typename T3> inline
void vec3_minus(T1 c[3], const T2 a[3], const T3 b[3])
{
	c[0] = a[0] - b[0];
	c[1] = a[1] - b[1];
	c[2] = a[2] - b[2];
}

/**
 @fn template<typename T1, typename T2, typename T3> inline void vec3_multi(T1 c[3], const T2 a[3], const T3 b[3])
 @brief c[] = a[] * b[]
 */
template<typename T1, typename T2, typename T3> inline
void vec3_multi(T1 c[3], const T2 a[3], const T3 b[3])
{
	c[0] = a[0] * b[0];
	c[1] = a[1] * b[1];
	c[2] = a[2] * b[2];
}

/**
 @fn template<typename T1, typename T2, typename T3> inline void vec3_multi(T1 c[3], const T2 a[3], const T3 b)
 @brief c[] = a[] * b
 */
template<typename T1, typename T2, typename T3> inline
void vec3_multi(T1 c[3], const T2 a[3], const T3 b)
{
	c[0] = a[0] * b;
	c[1] = a[1] * b;
	c[2] = a[2] * b;
}

/**
 @fn template<typename T1, typename T2, typename T3> inline void vec3_div(T1 c[3], const T2 a[3], const T3 b[3])
 @brief c[] = a[] / b[]
 */
template<typename T1, typename T2, typename T3> inline
void vec3_div(T1 c[3], const T2 a[3], const T3 b[3])
{
	c[0] = a[0] / b[0];
	c[1] = a[1] / b[1];
	c[2] = a[2] / b[2];
}

/**
 @fn template<typename T1, typename T2, typename T3> inline void vec3_div(T1 c[3], const T2 a[3], const T3 b)
 @brief c[] = a[] / b
 */
template<typename T1, typename T2, typename T3> inline
void vec3_div(T1 c[3], const T2 a[3], const T3 b)
{
	T1 inv = 1. / b;
	c[0] = a[0] * inv;
	c[1] = a[1] * inv;
	c[2] = a[2] * inv;
}

} // end of namspace FB

#endif  // _FB_vec3_func_h
