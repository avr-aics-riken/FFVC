#ifndef tt_shape_h
#define tt_shape_h

#include <tt/gm.h>
#include <tt/gfx/transform.h>

class Shape
{
public:
	enum ShapeType
	{
		BOX,
		SPHERE,
		CYLINDER,
		RECTANGLE,
		DISK,
		TRIANGLE_MESH,
		TRIANGLE,
	};

	Shape();
	virtual ~Shape() {}

	virtual ShapeType getTypeId() const =0;
	virtual std::string getTypeName() const =0;

	virtual bool isConvex() const =0;
	virtual BBox getBBox() const =0;

	virtual float area() const;
	virtual Vec3f sample(float u, float v) const;

	int m_sides;
	int m_smooth;
	int m_revnml;
	Vec2f m_minvt, m_maxvt;

protected:
	Vec2f uv2texcoord(float u, float v) const;
};

#endif  // tt_shape_h

