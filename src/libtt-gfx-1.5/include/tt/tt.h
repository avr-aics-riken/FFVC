#ifndef tt_h
#define tt_h

typedef unsigned char  uchar;
typedef unsigned short ushort;
typedef unsigned int   uint;
typedef unsigned long  ulong;

#ifndef M_PI
# define M_PI		3.14159265358979323846	/* pi */
#endif // !M_PI

#ifndef __GNUC__
#include <math.h>
#define sqrtf(x) (float)sqrt((double)(x))
#define log2f(x) (float)(log((double)(x))/log(2.0))

inline double cbrt(double x, double e =0.0001) {
	register double p = 0.0;
	register double q = x;
	while ( fabs(q - p) > e ) {
		p = q;
		q = p - (p*p*p - x)/(3.*p*p);
	}
	return q;
}
#endif // !__GNUC__

#define GL_GLEXT_PROTOTYPES 1

#endif  // tt_h

