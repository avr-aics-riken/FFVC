#include "shape.h"

Shape::Shape()
: m_sides(2), m_smooth(1), m_revnml(0), m_minvt(0, 0), m_maxvt(1, 1)
{
}

Vec2f
Shape::uv2texcoord(float u, float v) const
{
	float s = (1 - u) * m_minvt[0] + u * m_maxvt[0];
	float t = (1 - v) * m_minvt[1] + v * m_maxvt[1];
	return Vec2f(s, t);
}

float
Shape::area() const
{
	return 0.;
}

Vec3f
Shape::sample(float u, float v) const
{
	return Vec3f(0, 0, 0);
}

