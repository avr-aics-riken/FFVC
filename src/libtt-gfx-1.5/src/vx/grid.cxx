#include "grid.h"

VxGrid::VxGrid(int nx, int ny, int nz)
: m_gp(nx, ny, nz)
{
	m_nelm[0] = nx;
	m_nelm[1] = ny;
	m_nelm[2] = nz;
	m_data = new char[nx * ny * nz];
	clear();
}

VxGrid::~VxGrid()
{
	delete [] m_data;
}

void
VxGrid::clear()
{
	int sz = m_nelm[0] * m_nelm[1] * m_nelm[2];
	for (int i=0; i<sz; i++) m_data[i] = 0;
}

int
VxGrid::countVoxels(int x0, int y0, int z0, int x1, int y1, int z1, int id) const
{
	int x, y, z;
	int count = 0;

	for (z=z0; z<z1; z++)
	{
		for (y=y0; y<y1; y++)
		{
			for (x=x0; x<x1; x++)
			{
				int idx = m_gp.xyz2idx(x, y, z);
				if (getId(idx) == id) count++;
			}
		}
	}
	return count;
}

int
VxGrid::fill(int _ix, int _iy, int _iz, int id)
{
	int count = 0;
	std::queue<int> idx_queue;
	int idx;
	int ix, iy, iz;
	int next;
	int min[3];
	int max[3];
	int dx = 1;
	int dy = m_nelm[0];
	int dz = m_nelm[0] * m_nelm[1];

	min[0] = min[1] = min[2] = 0;
	max[0] = m_nelm[0] - 1;
	max[1] = m_nelm[1] - 1;
	max[2] = m_nelm[2] - 1;

	// push the root voxel
	if (m_gp.isInside(_ix, _iy, _iz))
	{
		idx = m_gp.xyz2idx(_ix, _iy, _iz);
		if (!getMark(idx))
		{
			pushVoxel(idx_queue, idx, id);
			count++;
		}
	}

	while (idx_queue.size() > 0)
	{
		idx = idx_queue.front();
		idx_queue.pop();
		m_gp.idx2xyz(idx, &ix, &iy, &iz);

		if (ix > min[0])
		{
			next = idx - dx;
			if (!getMark(next))
			{
				pushVoxel(idx_queue, next, id);
				count++;
			}
		}
		if (ix < max[0])
		{
			next = idx + dx;
			if (!getMark(next))
			{
				pushVoxel(idx_queue, next, id);
				count++;
			}
		}
		if (iy > min[1])
		{
			next = idx - dy;
			if (!getMark(next))
			{
				pushVoxel(idx_queue, next, id);
				count++;
			}
		}
		if (iy < max[1])
		{
			next = idx + dy;
			if (!getMark(next))
			{
				pushVoxel(idx_queue, next, id);
				count++;
			}
		}
		if (iz > min[2])
		{
			next = idx - dz;
			if (!getMark(next))
			{
				pushVoxel(idx_queue, next, id);
				count++;
			}
		}
		if (iz < max[2])
		{
			next = idx + dz;
			if (!getMark(next))
			{
				pushVoxel(idx_queue, next, id);
				count++;
			}
		}
	}

	return count;
}

