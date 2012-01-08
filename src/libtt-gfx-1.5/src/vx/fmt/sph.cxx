#include <string>
#include <sstream>
#include <tt/util/file.h>
#include <tt/util/message.h>
#include <tt/vx/voxelizer_octree.h>
#include "sph.h"

int sph_save(const VxVoxelizer* _vox, const char* fname, const char* opt)
{
	tt_info("[sph] saving... (%s)\n", fname);

	const VxVoxelizerOctree* vox = dynamic_cast<const VxVoxelizerOctree*>(_vox);
	std::istringstream is_opt(opt);

	TtOutFile ofile(fname);

	int len;
	int svType;
	int dType;
	int imax, jmax, kmax;
	float xorg, yorg, zorg;
	float xpitch, ypitch, zpitch;
	int step;
	float time;

	svType = 1;
	dType = 1;
	vox->getSize(&imax, &jmax, &kmax);
	vox->getOrigin(&xorg, &yorg, &zorg);
	vox->getPitch(&xpitch, &ypitch, &zpitch);
	step = 0;
	time = 0;

	std::string token;
	while (is_opt >> token)
	{
		if (token == "nelm")
		{
			is_opt >> imax >> jmax >> kmax;
		}
	}

	len = 4 * 2;
	ofile.write(&len, 4, 1);
	ofile.write(&svType, 4, 1);
	ofile.write(&dType, 4, 1);
	ofile.write(&len, 4, 1);

	len = 4 * 3;
	ofile.write(&len, 4, 1);
	ofile.write(&imax, 4, 1);
	ofile.write(&jmax, 4, 1);
	ofile.write(&kmax, 4, 1);
	ofile.write(&len, 4, 1);

	len = 4 * 3;
	ofile.write(&len, 4, 1);
	ofile.write(&xorg, 4, 1);
	ofile.write(&yorg, 4, 1);
	ofile.write(&zorg, 4, 1);
	ofile.write(&len, 4, 1);

	len = 4 * 3;
	ofile.write(&len, 4, 1);
	ofile.write(&xpitch, 4, 1);
	ofile.write(&ypitch, 4, 1);
	ofile.write(&zpitch, 4, 1);
	ofile.write(&len, 4, 1);

	len = 4 * 2;
	ofile.write(&len, 4, 1);
	ofile.write(&step, 4, 1);
	ofile.write(&time, 4, 1);
	ofile.write(&len, 4, 1);

	int x, y, z;
	float volume;
	float area[6];
	len = 4 * imax * jmax * kmax;
	ofile.write(&len, 4, 1);
	for (z=0; z<kmax; z++)
	{
		for (y=0; y<jmax; y++)
		{
			for (x=0; x<imax; x++)
			{
				vox->computeRatio(x, y, z, &volume, area);
				ofile.write(&volume, 4, 1);
			}
		}
	}
	ofile.write(&len, 4, 1);

	return 1;
}

