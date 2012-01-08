#ifndef vxgrid_h
#define vxgrid_h

#include <queue>
#include <tt/gfx/gridpos.h>

class VxGrid
{
public:
	VxGrid(int nx, int ny, int nz);
	~VxGrid();

	void getSize(int* nx, int* ny, int* nz) const {
		*nx = m_nelm[0]; *ny = m_nelm[1]; *nz = m_nelm[2];
	}

	void clear();

	int countVoxels(int x0, int y0, int z0, int x1, int y1, int z1, int id) const;
	int fill(int ix, int iy, int iz, int id);

	int getMark(int idx) const {
		return m_data[idx] & 1;
	}
	void setMark(int idx, int b) {
		m_data[idx] &= ~1;
		m_data[idx] |= b;
	}

	int isBoundary(int idx) const {
		return m_data[idx] & 2;
	}
	void setBoundary(int idx, int b) {
		m_data[idx] &= ~2;
		m_data[idx] |= b << 1;
	}

	int getId(int idx) const {
		return m_data[idx] >> 2;
	}
	void setId(int idx, int id) {
		m_data[idx] &= 3;
		m_data[idx] |= id << 2;
	}

private:
	void pushVoxel(std::queue<int>& idx_queue, int idx, int id)
	{
		idx_queue.push(idx);
		setMark(idx, 1);
		setId(idx, id);
	}

	int m_nelm[3];
	char* m_data;
	GridIdx m_gp;
};

#endif  // vxgrid_h

