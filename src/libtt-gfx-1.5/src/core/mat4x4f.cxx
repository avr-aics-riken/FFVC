//=========================================================================
// SEE:
// Graphics Gems V
// libgm++: Graphics Math Library
// Ferdi Scheepers and Stephen F May
// 15 June 1994
//=========================================================================

#include "mat4x4f.h"
#include "gm_util.h"

// dot product of row i of matrix A and column j of matrix B
inline float
RCD(const Mat4x4f& A, const Mat4x4f& B, int i, int j)
{
	return A[i][0] * B[0][j] + A[i][1] * B[1][j] +
	       A[i][2] * B[2][j] + A[i][3] * B[3][j];
}

// MINOR of M[r][c]; r in {0,1,2,3}-{r0,r1,r2}; c in {0,1,2,3}-{c0,c1,c2}
inline float
MINOR(const Mat4x4f& M, int r0, int r1, int r2, int c0, int c1, int c2)
{
	return M[r0][c0] * (M[r1][c1] * M[r2][c2] - M[r2][c1] * M[r1][c2]) -
	       M[r0][c1] * (M[r1][c0] * M[r2][c2] - M[r2][c0] * M[r1][c2]) +
	       M[r0][c2] * (M[r1][c0] * M[r2][c1] - M[r2][c0] * M[r1][c1]);
}

//=========================================================================

Mat4x4f::Mat4x4f()
{
	assign(0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0);
}

Mat4x4f::Mat4x4f(float a00, float a01, float a02, float a03,
	float a10, float a11, float a12, float a13,
	float a20, float a21, float a22, float a23,
	float a30, float a31, float a32, float a33)
{
	assign(a00, a01, a02, a03,
	       a10, a11, a12, a13,
	       a20, a21, a22, a23,
	       a30, a31, a32, a33);
}

Mat4x4f&
Mat4x4f::assign(float a00, float a01, float a02, float a03,
	float a10, float a11, float a12, float a13,
	float a20, float a21, float a22, float a23,
	float a30, float a31, float a32, float a33)
{
	m_M[0][0] = a00; m_M[0][1] = a01; m_M[0][2] = a02; m_M[0][3] = a03;
	m_M[1][0] = a10; m_M[1][1] = a11; m_M[1][2] = a12; m_M[1][3] = a13;
	m_M[2][0] = a20; m_M[2][1] = a21; m_M[2][2] = a22; m_M[2][3] = a23;
	m_M[3][0] = a30; m_M[3][1] = a31; m_M[3][2] = a32; m_M[3][3] = a33;
	return *this;
}

Mat4x4f&
Mat4x4f::operator *=(const Mat4x4f& M)
{
	assign(
	RCD(*this, M, 0, 0), RCD(*this, M, 0, 1), 
	RCD(*this, M, 0, 2), RCD(*this, M, 0, 3),
	RCD(*this, M, 1, 0), RCD(*this, M, 1, 1),
	RCD(*this, M, 1, 2), RCD(*this, M, 1, 3),
	RCD(*this, M, 2, 0), RCD(*this, M, 2, 1),
	RCD(*this, M, 2, 2), RCD(*this, M, 2, 3),
	RCD(*this, M, 3, 0), RCD(*this, M, 3, 1),
	RCD(*this, M, 3, 2), RCD(*this, M, 3, 3));
	return *this;
}

Mat4x4f&
Mat4x4f::operator *=(float d)
{
	m_M[0][0] *= d; m_M[0][1] *= d; m_M[0][2] *= d; m_M[0][3] *= d;
	m_M[1][0] *= d; m_M[1][1] *= d; m_M[1][2] *= d; m_M[1][3] *= d;
	m_M[2][0] *= d; m_M[2][1] *= d; m_M[2][2] *= d; m_M[2][3] *= d;
	m_M[3][0] *= d; m_M[3][1] *= d; m_M[3][2] *= d; m_M[3][3] *= d;
	return *this;
}

Mat4x4f&
Mat4x4f::operator /=(float d)
{
	assert(!gmIsZero(d));
	float di = 1 / d;
	return *this *= di;
}

Mat4x4f
Mat4x4f::operator -() const
{
	Mat4x4f M(*this);
	return M *= -1;
}

Mat4x4f
Mat4x4f::operator *(const Mat4x4f& N) const
{
	Mat4x4f M(*this);
	return M *= N;
}

Mat4x4f
Mat4x4f::operator *(float d) const
{
	Mat4x4f M(*this);
	return M *= d;
}

Mat4x4f
Mat4x4f::operator /(float d) const
{
	Mat4x4f M(*this);
	return M /= d;
}

Mat4x4f
operator *(float d, const Mat4x4f& N)
{
	Mat4x4f M(N);
	return M *= d;
}

bool
Mat4x4f::operator ==(const Mat4x4f& M) const
{
	for (int i=0; i<4; i++)
	{
		for (int j=0; j<4; j++)
		{
			if (!gmFuzEQ(m_M[i][j], M[i][j])) return false;
		}
	}
	return true;
}

bool
Mat4x4f::operator !=(const Mat4x4f& M) const
{
	return !(*this == M);
}

Vec4f
Mat4x4f::operator *(const Vec4f& v) const
{
	return Vec4f(
	m_M[0][0] * v[0] + m_M[0][1] * v[1] + m_M[0][2] * v[2] + m_M[0][3] * v[3],
	m_M[1][0] * v[0] + m_M[1][1] * v[1] + m_M[1][2] * v[2] + m_M[1][3] * v[3],
	m_M[2][0] * v[0] + m_M[2][1] * v[1] + m_M[2][2] * v[2] + m_M[2][3] * v[3],
	m_M[3][0] * v[0] + m_M[3][1] * v[1] + m_M[3][2] * v[2] + m_M[3][3] * v[3]);
}

Vec4f
operator *(const Vec4f& v, const Mat4x4f& M)
{
	return Vec4f(
	v[0] * M[0][0] + v[1] * M[1][0] + v[2] * M[2][0] + v[3] * M[3][0],
	v[0] * M[0][1] + v[1] * M[1][1] + v[2] * M[2][1] + v[3] * M[3][1],
	v[0] * M[0][2] + v[1] * M[1][2] + v[2] * M[2][2] + v[3] * M[3][2],
	v[0] * M[0][3] + v[1] * M[1][3] + v[2] * M[2][3] + v[3] * M[3][3]);
}

Mat4x4f
Mat4x4f::transpose() const
{
	return Mat4x4f(
	m_M[0][0], m_M[1][0], m_M[2][0], m_M[3][0],
	m_M[0][1], m_M[1][1], m_M[2][1], m_M[3][1],
	m_M[0][2], m_M[1][2], m_M[2][2], m_M[3][2],
	m_M[0][3], m_M[1][3], m_M[2][3], m_M[3][3]);
}

Mat4x4f
Mat4x4f::inverse() const
{
	//assert(!isSingular());
	return adjoint() * gmInv(determinant());
}

Mat4x4f
Mat4x4f::adjoint() const
{
	Mat4x4f A;

	A[0][0] =  MINOR(*this, 1, 2, 3, 1, 2, 3);
	A[0][1] = -MINOR(*this, 0, 2, 3, 1, 2, 3);
	A[0][2] =  MINOR(*this, 0, 1, 3, 1, 2, 3);
	A[0][3] = -MINOR(*this, 0, 1, 2, 1, 2, 3);
	A[1][0] = -MINOR(*this, 1, 2, 3, 0, 2, 3);
	A[1][1] =  MINOR(*this, 0, 2, 3, 0, 2, 3);
	A[1][2] = -MINOR(*this, 0, 1, 3, 0, 2, 3);
	A[1][3] =  MINOR(*this, 0, 1, 2, 0, 2, 3);
	A[2][0] =  MINOR(*this, 1, 2, 3, 0, 1, 3);
	A[2][1] = -MINOR(*this, 0, 2, 3, 0, 1, 3);
	A[2][2] =  MINOR(*this, 0, 1, 3, 0, 1, 3);
	A[2][3] = -MINOR(*this, 0, 1, 2, 0, 1, 3);
	A[3][0] = -MINOR(*this, 1, 2, 3, 0, 1, 2);
	A[3][1] =  MINOR(*this, 0, 2, 3, 0, 1, 2);
	A[3][2] = -MINOR(*this, 0, 1, 3, 0, 1, 2);
	A[3][3] =  MINOR(*this, 0, 1, 2, 0, 1, 2);

	return A;
}

float
Mat4x4f::determinant() const
{
	return
	m_M[0][0] * MINOR(*this, 1, 2, 3, 1, 2, 3) -
	m_M[0][1] * MINOR(*this, 1, 2, 3, 0, 2, 3) +
	m_M[0][2] * MINOR(*this, 1, 2, 3, 0, 1, 3) -
	m_M[0][3] * MINOR(*this, 1, 2, 3, 0, 1, 2);
}

bool
Mat4x4f::isSingular() const
{
	return gmIsZero(determinant());
}

Vec3f
Mat4x4f::transform(const Vec3f& v) const
{
	return Vec3f(
	v[0] * m_M[0][0] + v[1] * m_M[0][1] + v[2] * m_M[0][2] + m_M[0][3],
	v[0] * m_M[1][0] + v[1] * m_M[1][1] + v[2] * m_M[1][2] + m_M[1][3],
	v[0] * m_M[2][0] + v[1] * m_M[2][1] + v[2] * m_M[2][2] + m_M[2][3]);
}

Vec3f
Mat4x4f::transformVector(const Vec3f& v) const
{
	return Vec3f(
	v[0] * m_M[0][0] + v[1] * m_M[0][1] + v[2] * m_M[0][2],
	v[0] * m_M[1][0] + v[1] * m_M[1][1] + v[2] * m_M[1][2],
	v[0] * m_M[2][0] + v[1] * m_M[2][1] + v[2] * m_M[2][2]);
}

//=========================================================================
// transformation matrices
//=========================================================================
Mat4x4f
Mat4x4f::identity()
{
	return Mat4x4f(
	1, 0, 0, 0,
	0, 1, 0, 0,
	0, 0, 1, 0,
	0, 0, 0, 1);
}

Mat4x4f
Mat4x4f::translate(float x, float y, float z)
{
	return Mat4x4f(
	1, 0, 0, x,
	0, 1, 0, y,
	0, 0, 1, z,
	0, 0, 0, 1);
}

Mat4x4f
Mat4x4f::rotate(float angle, float x, float y, float z)
{
	Mat4x4f M;

	float len = sqrtf(x*x + y*y + z*z);
	float a = x / len;
	float b = y / len;
	float c = z / len;
	float aa = a * a;
	float bb = b * b;
	float cc = c * c;
	float sine = sin(gmRadians(angle));
	float cosine = cos(gmRadians(angle));
	float omcos = 1 - cosine;

	M[0][0] = aa + (1 - aa) * cosine;
	M[1][1] = bb + (1 - bb) * cosine;
	M[2][2] = cc + (1 - cc) * cosine;
	M[1][0] = a * b * omcos + c * sine;
	M[2][0] = a * c * omcos - b * sine;
	M[0][1] = a * b * omcos - c * sine;
	M[2][1] = b * c * omcos + a * sine;
	M[0][2] = a * c * omcos + b * sine;
	M[1][2] = b * c * omcos - a * sine;
	M[0][3] = M[1][3] = M[2][3] = M[3][0] = M[3][1] = M[3][2] = 0;
	M[3][3] = 1;

	return M;
}

Mat4x4f
Mat4x4f::scale(float x, float y, float z)
{
	return Mat4x4f(
	x, 0, 0, 0,
	0, y, 0, 0,
	0, 0, z, 0,
	0, 0, 0, 1);
}

Mat4x4f
Mat4x4f::translate(const Vec3f& v)
{
	return translate(v[0], v[1], v[2]);
}

Mat4x4f
Mat4x4f::rotate(float angle, const Vec3f& axis)
{
	return rotate(angle, axis[0], axis[1], axis[2]);
}

Mat4x4f
Mat4x4f::scale(const Vec3f& v)
{
	return scale(v[0], v[1], v[2]);
}

//=========================================================================
// cubic basis matrices
//
// the same convention as the book
// "Computer Graphics PRINCIPLES AND PRACTICE"
//
// Q(t) = [x(t) y(t) z(t)] = T * M * G
//
//                |P1|
// Q(t) = T * M * |P4|
//                |R1|
//                |R4|
//=========================================================================
Mat4x4f
Mat4x4f::hermiteBasis()
{
	return Mat4x4f(
	 2, -2,  1,  1,
	-3,  3, -2, -1,
	 0,  0,  1,  0,
	 1,  0,  0,  0);
}

Mat4x4f
Mat4x4f::bezierBasis()
{
	return Mat4x4f(
	-1,  3, -3,  1,
	 3, -6,  3,  0,
	-3,  3,  0,  0,
	 1,  0,  0,  0);
}

Mat4x4f
Mat4x4f::bsplineBasis()
{
	return Mat4x4f(
	-1,  3, -3,  1,
	 3, -6,  3,  0,
	-3,  0,  3,  0,
	 1,  4,  1,  0) / 6;
}

Mat4x4f
Mat4x4f::catmullromBasis()
{
	return Mat4x4f(
	-1,  3, -3,  1,
	 2, -5,  4, -1,
	-1,  0,  1,  0,
	 0,  2,  0,  0) / 2;
}

//=========================================================================

std::istream&
operator >>(std::istream& is, Mat4x4f& M)
{
	is >> M[0][0] >> M[0][1] >> M[0][2] >> M[0][3];
	is >> M[1][0] >> M[1][1] >> M[1][2] >> M[1][3];
	is >> M[2][0] >> M[2][1] >> M[2][2] >> M[2][3];
	is >> M[3][0] >> M[3][1] >> M[3][2] >> M[3][3];
	return is;
}

std::ostream&
operator <<(std::ostream& os, const Mat4x4f& M)
{
	os << M[0][0] << " " << M[0][1] << " " << M[0][2] << " " << M[0][3] << " ";
	os << M[1][0] << " " << M[1][1] << " " << M[1][2] << " " << M[1][3] << " ";
	os << M[2][0] << " " << M[2][1] << " " << M[2][2] << " " << M[2][3] << " ";
	os << M[3][0] << " " << M[3][1] << " " << M[3][2] << " " << M[3][3];
	return os;
}

std::ostream&
fancy_print(std::ostream& os, const Mat4x4f& M)
{
	os << M[0][0] << " " << M[0][1] << " " << M[0][2] << " " << M[0][3] << std::endl;
	os << M[1][0] << " " << M[1][1] << " " << M[1][2] << " " << M[1][3] << std::endl;
	os << M[2][0] << " " << M[2][1] << " " << M[2][2] << " " << M[2][3] << std::endl;
	os << M[3][0] << " " << M[3][1] << " " << M[3][2] << " " << M[3][3] << std::endl;
	return os;
}

