///////////////////////////////////////////////////////////////////////////
//
// graphics math constants and utility functions
//
// libgm++: Graphics Math Library
// Ferdi Scheepers and Stephen F May
// 15 June 1994
//
// Reference:
// - Graphics Gems V, Chapter 7 - 7
//
///////////////////////////////////////////////////////////////////////////

#ifndef GMUTILS_H
#define GMUTILS_H

#include <iostream>
#include <math.h>

//=========================================================================
// constants
//=========================================================================

float const gmEPSILON =   1.e-5;
float const gmPI =        3.14159265358979323846;
float const gm2PI =       6.28318530717958623200;
float const gmINVPI =     0.31830988618379069122;
float const gmPIDIV2 =    1.57079632679489655800;
float const gmPIDIV4 =    0.78539816339744827900;
float const gmINF =       1.e30;

float const gmDEGTORAD =  0.01745329251994329547;
float const gmRADTODEG = 57.29577951308232286465;

//=========================================================================

inline float gmAbs(float f)
{
  return (f >= 0) ? f : -f;
}

inline float gmSign(float f)
{
  return (f < 0) ? -1 : 1;
}

inline float gmZSign(float f)
{
  return (f > 0) ? 1 : (f < 0) ? -1 : 0;
}

//=========================================================================

inline float gmSqr(float f)
{
  return f * f;
}

inline float gmCube(float f)
{
  return f * f * f;
}

inline float gmInv(float f)
{
  return 1. / f;
}

//=========================================================================
// rounding
//=========================================================================

inline float gmTrunc(float f)
{
  return float(int(f));
}

inline float gmRound(float f)
{
  return (f >= 0) ? float(int(f + .5)) : float(- int(.5 - f));
}

inline float gmFloor(float f)
{
  return (f == int(f)) ? f : (f > 0) ? float(int(f)) : float(int(f) - 1);
}

inline float gmCeil(float f)
{
  return (f == int(f)) ? f : (f > 0) ? float(int(f) + 1) : float(int(f));
}

//=========================================================================
// degrees and radians
//=========================================================================

inline float gmDegrees(float f)
{
  return f * gmRADTODEG;
}

inline float gmRadians(float f)
{
  return f * gmDEGTORAD;
}

//=========================================================================
// comparison
//=========================================================================

inline bool gmIsZero(float f, float epsilon = gmEPSILON)
{
  return (gmAbs(f) < epsilon);
}

// f ~= g
inline bool gmFuzEQ(float f, float g, float epsilon = gmEPSILON)
{
  return (f <= g) ? (f >= g - epsilon) : (f <= g + epsilon);
}

// f ~>= g
inline bool gmFuzGEQ(float f, float g, float epsilon = gmEPSILON)
{
  return (f >= g - epsilon);
}

// f ~<= g
inline bool gmFuzLEQ(float f, float g, float epsilon = gmEPSILON)
{
  return (f <= g + epsilon);
}

//=========================================================================
// min & max
//=========================================================================

inline float gmMin(float f, float g)
{
  return (f < g) ? f : g;
}

inline float gmMax(float f, float g)
{
  return (f > g) ? f : g;
}

inline float gmMin3(float f, float g, float h)
{
  return (f < g) ? gmMin(f, h) : gmMin(g, h);
}

inline float gmMax3(float f, float g, float h)
{
  return (f > g) ? gmMax(f, h) : gmMax(g, h);
}

inline int gmMinAxis(float x, float y, float z)
{
	if (x < y)
	{
		if (x < z) return 0;
		else       return 2;
	}
	else
	{
		if (y < z) return 1;
		else       return 2;
	}
}

inline int gmMaxAxis(float x, float y, float z)
{
	if (x > y)
	{
		if (x > z) return 0;
		else       return 2;
	}
	else
	{
		if (y > z) return 1;
		else       return 2;
	}
}

//=========================================================================

template<typename T>
inline void gmSwap(T& f, T& g)
{
  T tmp = f;
  f = g;
  g = tmp;
}

inline float gmClamp(float val, float low, float high)
{
	if (val < low)
		return low;
	else if (val > high)
		return high;
	else
		return val;
}

inline float gmLerp(float f, float l, float h)
{
  return l + ((h - l) * f );
}

inline float gmSmooth(float f)
{
  return (3. - 2. * f) * f * f;
}



//=========================================================================
// functions for integer parameters
//=========================================================================

inline int gmAbs_i(int f)
{
  return (f >= 0) ? f : -f;
}

//=========================================================================
// min & max
//=========================================================================

inline int gmMin_i(int f, int g)
{
  return (f < g) ? f : g;
}

inline int gmMax_i(int f, int g)
{
  return (f > g) ? f : g;
}

inline int gmMin3_i(int f, int g, int h)
{
  return (f < g) ? gmMin_i(f, h) : gmMin_i(g, h);
}

inline int gmMax3_i(int f, int g, int h)
{
  return (f > g) ? gmMax_i(f, h) : gmMax_i(g, h);
}

//=========================================================================

inline int gmClamp_i(int val, int low, int high)
{
	if (val < low)
		return low;
	else if (val > high)
		return high;
	else
		return val;
}

#endif  // GMUTILS_H

