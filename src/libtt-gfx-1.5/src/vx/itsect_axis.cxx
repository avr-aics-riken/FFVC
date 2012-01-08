// compute an intersection of a triangle and the x-axis.
// return the intersection point in the x-coordinate.
bool f_itsect_tri_xaxis(const float p0[3], const float p1[3], const float p2[3], float* x, float y, float z)
{
	float a0 = (p1[1]-y)*(p2[2]-z) - (p1[2]-z)*(p2[1]-y);
	float a1 = (p2[1]-y)*(p0[2]-z) - (p2[2]-z)*(p0[1]-y);
	float a2 = (p0[1]-y)*(p1[2]-z) - (p0[2]-z)*(p1[1]-y);
	if(a0*a1 < 0 || a1*a2 < 0 || a2*a0 < 0)
		return false;
	*x = (a0*p0[0] + a1*p1[0] + a2*p2[0])/(a0 + a1 + a2);
	return true;
}

// compute an intersection of a triangle and the y-axis.
// return the intersection point in the y-coordinate.
bool f_itsect_tri_yaxis(const float p0[3], const float p1[3], const float p2[3], float x, float* y, float z)
{
	float a0 = (p1[2]-z)*(p2[0]-x) - (p1[0]-x)*(p2[2]-z);
	float a1 = (p2[2]-z)*(p0[0]-x) - (p2[0]-x)*(p0[2]-z);
	float a2 = (p0[2]-z)*(p1[0]-x) - (p0[0]-x)*(p1[2]-z);
	if(a0*a1 < 0 || a1*a2 < 0 || a2*a0 < 0)
		return false;
	*y = (a0*p0[1] + a1*p1[1] + a2*p2[1])/(a0 + a1 + a2);
	return true;
}

// compute an intersection of a triangle and the z-axis.
// return the intersection point in the z-coordinate.
bool f_itsect_tri_zaxis(const float p0[3], const float p1[3], const float p2[3], float x, float y, float* z)
{
	float a0 = (p1[0]-x)*(p2[1]-y) - (p1[1]-y)*(p2[0]-x);
	float a1 = (p2[0]-x)*(p0[1]-y) - (p2[1]-y)*(p0[0]-x);
	float a2 = (p0[0]-x)*(p1[1]-y) - (p0[1]-y)*(p1[0]-x);
	if(a0*a1 < 0 || a1*a2 < 0 || a2*a0 < 0)
		return false;
	*z = (a0*p0[2] + a1*p1[2] + a2*p2[2])/(a0 + a1 + a2);
	return true;
}

bool f_itsect_tri_axis(int axis, const float p0[3], const float p1[3], const float p2[3], float* x, float* y, float* z)
{
	switch (axis)
	{
		case 0: return f_itsect_tri_xaxis(p0, p1, p2, x, *y, *z);
		case 1: return f_itsect_tri_yaxis(p0, p1, p2, *x, y, *z);
		case 2: return f_itsect_tri_zaxis(p0, p1, p2, *x, *y, z);
	}
	return false;
}

