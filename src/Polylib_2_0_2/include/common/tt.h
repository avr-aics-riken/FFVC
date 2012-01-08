/*
 * Polylib - Polygon Management Library.
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2010-
 *
 */
#ifndef polylib_tt_h
#define polylib_tt_h

typedef unsigned char  uchar;
typedef unsigned short ushort;
typedef unsigned int   uint;
typedef unsigned long  ulong;

#ifndef M_PI
# define M_PI		3.14159265358979323846	/* pi */
#endif

#ifndef __GNUC__
# include <math.h>
# define sqrtf(x) (float)sqrt((double)(x))
#endif

#define GL_GLEXT_PROTOTYPES 1

#endif  // polylib_tt_h

