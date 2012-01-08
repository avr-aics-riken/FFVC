#include <stdio.h>
#include <iostream>
#include <sstream>
#include <tt/tt.h>
#include <tt/util/message.h>
#include <tt/util/file.h>
#include <tt/vx/voxelizer_octree.h>
#include "svx.h"

using namespace std;

int
svx_load(VxVoxelizer* _vox, const char* fname, const char* opt)
{
	tt_info("[svx] loading... (%s)\n", fname);

	const VxVoxelizerOctree* vox = dynamic_cast<const VxVoxelizerOctree*>(_vox);

	TtInFile ifile(fname);

	int   x, y, z;
	uint  padding;
	float origin[3];
	float pitch[3];
	int   nelm[3];
	int   type = 0;

	ifile.read(&padding, sizeof(uint), 1);
	if (padding != 12)
	{
		ifile.setInverseFlag(1);
	}
	ifile.read(nelm, sizeof(uint), 3);
	ifile.read(&padding, sizeof(uint), 1);

	ifile.read(&padding, sizeof(uint), 1);
	ifile.read(origin, sizeof(float), 3);
	ifile.read(&padding, sizeof(uint), 1);

	ifile.read(&padding, sizeof(uint), 1);
	ifile.read(pitch, sizeof(float), 3);
	ifile.read(&padding, sizeof(uint), 1);

	ifile.read(&padding, sizeof(uint), 1);
	ifile.read(&type, sizeof(int), 1);
	ifile.read(&padding, sizeof(uint), 1);

	Vec3f min = Vec3f(origin[0], origin[1], origin[2]);
	Vec3f max = Vec3f(pitch[0] * nelm[0], pitch[1] * nelm[1], pitch[2] * nelm[2]) + min;
	BBox bbox(min, max);

	if (type == 0 || type & (0x1<<0))
	{
		// read volume ratio
		ifile.read(&padding, sizeof(uint), 1);
		for (z=0; z<nelm[2]; z++)
		{
			for (y=0; y<nelm[1]; y++)
			{
				for (x=0; x<nelm[0]; x++)
				{
					float data;
					ifile.read(&data, sizeof(float), 1);
				}
			}
		}
		ifile.read(&padding, sizeof(uint), 1);
	}

	if (type == 0 || type & (0x1<<1))
	{
		// read area ratio
	}

	if (type == 0 || type & (0x1<<2))
	{
		// read medium ID
		ifile.read(&padding, sizeof(uint), 1);
		for (z=0; z<nelm[2]; z++)
		{
			for (y=0; y<nelm[1]; y++)
			{
				for (x=0; x<nelm[0]; x++)
				{
					int data;
					ifile.read(&data, sizeof(int), 1);
				}
			}
		}
		ifile.read(&padding, sizeof(uint), 1);
	}

	if (type == 0 || type & (0x1<<3))
	{
		// read cell boundary condition ID
	}

	if (type == 0 || type & (0x1<<4))
	{
		// read face boundary condition ID
	}

	return 1;
}

int
svx_save(const VxVoxelizer* _vox, const char* fname, const char* opt)
{
	tt_info("[svx] saving... (%s)\n", fname);

	const VxVoxelizerOctree* vox = dynamic_cast<const VxVoxelizerOctree*>(_vox);
	std::istringstream is_opt(opt);

	TtOutFile ofile(fname);

	int   x, y, z;
	uint  padding;
	float origin[3];
	float pitch[3];
	int   nelm[3];
	int   type = 5;

	vox->getOrigin(&origin[0], &origin[1], &origin[2]);
	vox->getPitch(&pitch[0], &pitch[1], &pitch[2]);
	vox->getSize(&nelm[0], &nelm[1], &nelm[2]);

	string token;
	while (is_opt >> token)
	{
		if (token == "type")
		{
			is_opt >> type;
		}
		else if (token == "nelm")
		{
			is_opt >> nelm[0] >> nelm[1] >> nelm[2];
		}
	}

	padding = sizeof(uint) * 3;
	ofile.write(&padding, sizeof(uint), 1);
	ofile.write(nelm, sizeof(uint), 3);
	ofile.write(&padding, sizeof(uint), 1);

	padding = sizeof(float) * 3;
	ofile.write(&padding, sizeof(uint), 1);
	ofile.write(origin, sizeof(float), 3);
	ofile.write(&padding, sizeof(uint), 1);

	padding = sizeof(float) * 3;
	ofile.write(&padding, sizeof(uint), 1);
	ofile.write(pitch, sizeof(float), 3);
	ofile.write(&padding, sizeof(uint), 1);

	padding = sizeof(int);
	ofile.write(&padding, sizeof(uint), 1);
	ofile.write(&type, sizeof(int), 1);
	ofile.write(&padding, sizeof(uint), 1);

	if (type == 0 || type & (0x1<<0))
	{
		// write volume ratio
		padding = sizeof(float) * nelm[0] * nelm[1] * nelm[2];
		ofile.write(&padding, sizeof(uint), 1);
		for (z=0; z<nelm[2]; z++)
		{
			for (y=0; y<nelm[1]; y++)
			{
				for (x=0; x<nelm[0]; x++)
				{
					float volume;
					float area[6];
					vox->computeRatio(x, y, z, &volume, area);
					ofile.write(&volume, sizeof(float), 1);
				}
			}
		}
		ofile.write(&padding, sizeof(uint), 1);
	}

	if (type == 0 || type & (0x1<<1))
	{
		// write area ratio for the x-axis
		padding = sizeof(float) * (nelm[0]+1) * nelm[1] * nelm[2];
		ofile.write(&padding, sizeof(uint), 1);
		for (z=0; z<nelm[2]; z++)
		{
			for (y=0; y<nelm[1]; y++)
			{
				for (x=0; x<nelm[0]; x++)
				{
					float volume;
					float area[6];
					vox->computeRatio(x, y, z, &volume, area);
					ofile.write(&area[0], sizeof(float), 1);
					if (x == nelm[0]-1)
					{
						ofile.write(&area[1], sizeof(float), 1);
					}
				}
			}
		}
		ofile.write(&padding, sizeof(uint), 1);

		// write area ratio for the y-axis
		padding = sizeof(float) * nelm[0] * (nelm[1]+1) * nelm[2];
		ofile.write(&padding, sizeof(uint), 1);
		for (z=0; z<nelm[2]; z++)
		{
			for (y=0; y<nelm[1]; y++)
			{
				for (x=0; x<nelm[0]; x++)
				{
					float volume;
					float area[6];
					vox->computeRatio(x, y, z, &volume, area);
					ofile.write(&area[2], sizeof(float), 1);
					if (y == nelm[1]-1)
					{
						ofile.write(&area[3], sizeof(float), 1);
					}
				}
			}
		}
		ofile.write(&padding, sizeof(uint), 1);

		// write area ratio for the z-axis
		padding = sizeof(float) * nelm[0] * nelm[1] * (nelm[2]+1);
		ofile.write(&padding, sizeof(uint), 1);
		for (z=0; z<nelm[2]; z++)
		{
			for (y=0; y<nelm[1]; y++)
			{
				for (x=0; x<nelm[0]; x++)
				{
					float volume;
					float area[6];
					vox->computeRatio(x, y, z, &volume, area);
					ofile.write(&area[4], sizeof(float), 1);
					if (z == nelm[2]-1)
					{
						ofile.write(&area[5], sizeof(float), 1);
					}
				}
			}
		}
		ofile.write(&padding, sizeof(uint), 1);
	}

	if (type == 0 || type & (0x1<<2))
	{
		// write medium ID
		padding = sizeof(int) * nelm[0] * nelm[1] * nelm[2];
		ofile.write(&padding, sizeof(uint), 1);
		for (z=0; z<nelm[2]; z++)
		{
			for (y=0; y<nelm[1]; y++)
			{
				for (x=0; x<nelm[0]; x++)
				{
					int id;
					id = vox->getId(x, y, z);
					ofile.write(&id, sizeof(int), 1);
				}
			}
		}
		ofile.write(&padding, sizeof(uint), 1);
	}

	if (type == 0 || type & (0x1<<3))
	{
		// write cell boundary condition ID
	}

	if (type == 0 || type & (0x1<<4))
	{
		// write face boundary condition ID
	}

	return 1;
}
