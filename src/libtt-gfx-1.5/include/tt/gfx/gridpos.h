#ifndef gridpos_h
#define gridpos_h

#include <tt/gfx/bbox.h>

// nx, ny, nz: number of elements
// ix, iy, iz: indexed position: [0, nx-1], [0, ny-1], [0, nz-1]
// x, y, z: normalized position: [0, 1), [0, 1), [0, 1)
// idx: sequential index: [0, nx*ny*nz-1]

class GridIdx
{
public:
	GridIdx() {}
	GridIdx(int nx, int ny, int nz) {
		m_nelm[0] = nx;
		m_nelm[1] = ny;
		m_nelm[2] = nz;
		m_area = nx * ny;
	}

	int isInside(int ix, int iy, int iz) const {
		if (ix < 0 || ix >= m_nelm[0] ||
			iy < 0 || iy >= m_nelm[1] ||
			iz < 0 || iz >= m_nelm[2])
			return 0;
		return 1;
	}

	int xidx(int dx) const {
		return dx;
	}
	int yidx(int dy) const {
		return dy * m_nelm[0];
	}
	int zidx(int dz) const {
		return dz * m_area;
	}

	int xyz2idx(int ix, int iy, int iz) const {
		return ix + iy * m_nelm[0] + iz * m_area;
	}

	void idx2xyz(int idx, int* ix, int* iy, int* iz) const {
		*iz = idx / m_area;
		int m = idx % m_area;
		*iy = m / m_nelm[0];
		*ix = m % m_nelm[0];
	}

	// an indexed position to a normalized position
	void xyz2rpos(int ix, int iy, int iz, float* x, float* y, float* z) const {
		*x = float(ix) / m_nelm[0];
		*y = float(iy) / m_nelm[1];
		*z = float(iz) / m_nelm[2];
	}

	// a normalized position to an indexed position
	void rpos2xyz(float x, float y, float z, int* ix, int* iy, int* iz) const {
		*ix = int(x * m_nelm[0]);
		*iy = int(y * m_nelm[1]);
		*iz = int(z * m_nelm[2]);
	}

	void getNElements(int* nx, int* ny, int* nz) const {
		*nx = m_nelm[0];
		*ny = m_nelm[1];
		*nz = m_nelm[2];
	}
	void setNElements(int nx, int ny, int nz) {
		m_nelm[0] = nx;
		m_nelm[1] = ny;
		m_nelm[2] = nz;
		m_area = nx * ny;
	}

protected:
	int m_nelm[3];
	int m_area;
};

// x, y, z: world position: [minx, maxx), [miny, maxy), [minz, maxz)
class GridPos : public GridIdx
{
public:
	GridPos() {}
	GridPos(const BBox& bbox, int nx, int ny, int nz)
	: GridIdx(nx, ny, nz) {
		m_org[0] = bbox.min[0];
		m_org[1] = bbox.min[1];
		m_org[2] = bbox.min[2];
		m_len[0] = bbox.length(0);
		m_len[1] = bbox.length(1);
		m_len[2] = bbox.length(2);
	}

	// an indexed position to a world position
	void xyz2pos(int ix, int iy, int iz, float* x, float* y, float* z) const {
		*x = ix * m_len[0] / m_nelm[0] + m_org[0];
		*y = iy * m_len[1] / m_nelm[1] + m_org[1];
		*z = iz * m_len[2] / m_nelm[2] + m_org[2];
	}

	// a world position to an indexed position
	void pos2xyz(float x, float y, float z, int* ix, int* iy, int* iz) const {
		*ix = int(m_nelm[0] * (x - m_org[0]) / m_len[0]);
		*iy = int(m_nelm[1] * (y - m_org[1]) / m_len[1]);
		*iz = int(m_nelm[2] * (z - m_org[2]) / m_len[2]);
	}

	void getOrigin(float* ox, float* oy, float* oz) const {
		*ox = m_org[0];
		*oy = m_org[1];
		*oz = m_org[2];
	}

	void getLength(float* lx, float* ly, float* lz) const {
		*lx = m_len[0];
		*ly = m_len[1];
		*lz = m_len[2];
	}

	void getPitch(float* px, float* py, float* pz) const {
		*px = m_len[0] / m_nelm[0];
		*py = m_len[1] / m_nelm[1];
		*pz = m_len[2] / m_nelm[2];
	}

	BBox getBBox() const {
		return BBox(m_org[0], m_org[1], m_org[2], 
				m_org[0] + m_len[0], m_org[1] + m_len[1], m_org[2] + m_len[2]); 
	}
	void setBBox(const BBox& bbox) {
		m_org[0] = bbox.min[0];
		m_org[1] = bbox.min[1];
		m_org[2] = bbox.min[2];
		m_len[0] = bbox.length(0);
		m_len[1] = bbox.length(1);
		m_len[2] = bbox.length(2);
	}

private:
	float m_org[3];
	float m_len[3];
};

#endif  // gridpos_h

