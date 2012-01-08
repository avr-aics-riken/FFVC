#ifndef Ray_h
#define Ray_h

#include <tt/gm.h>

#define RAY_EPSILON  1e-3
#define RAY_INFINITY 1e09

class Ray
{
public:
	Ray() : min_t(0), max_t(RAY_INFINITY) {}
	Ray(const Vec3f& o, const Vec3f& d) : org(o), dir(d), min_t(0), max_t(RAY_INFINITY) {}

	Vec3f pointAt(float t) const { return org + t * dir; }

	Vec3f org, dir;
	mutable float min_t, max_t;
};


#endif  // Ray_h

