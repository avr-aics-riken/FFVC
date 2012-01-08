#include <stdio.h>
#include <iostream>
#include <sstream>
#include <tt/tt.h>
#include <tt/util/message.h>
#include <tt/util/file.h>
#include <tt/vx/voxelizer_octree.h>
#include <tt/vx/octree.h>
#include "ovx.h"

using namespace std;

int
ovx_load(VxVoxelizerOctree* vox, const char* fname, const char* opt)
{
	return 0;
}

class VxOctreeWriteOVX : public VxOctreeTravFunc
{
public:
	VxOctreeWriteOVX(TtOutFile* ofile, const VxVoxelizerOctree* vox, uint datatype) : m_ofile(ofile), m_vox(vox), m_datatype(datatype)
	{
	}
	void doCommon(VxOctreeNode* node)
	{
		uint padding;
		uint rpos;
		uint hasChild;

		rpos = pdg_top(node->getPosi());
		hasChild = node->isLeaf() ? 0 : 1;

		padding = sizeof(uint) * 2;
		m_ofile->write(&padding, sizeof(uint), 1);
		m_ofile->write(&rpos, sizeof(uint), 1);
		m_ofile->write(&hasChild, sizeof(uint), 1);
		m_ofile->write(&padding, sizeof(uint), 1);

		float volume;
		float area[6];
		int id;

		if (node->isLeaf())
		{
			m_vox->computeRatio(node, &volume, area);
			id = node->getId();
		}
		else
		{
			volume = 0;
			for (int i=0; i<6; i++) area[i] = 0;
			id = 0;
		}

		if (m_datatype == 0 || m_datatype & (0x1<<0))
		{
			// write volume ratio
			padding = sizeof(float);
			m_ofile->write(&padding, sizeof(uint), 1);
			m_ofile->write(&volume, sizeof(float), 1);
			m_ofile->write(&padding, sizeof(uint), 1);
		}

		if (m_datatype == 0 || m_datatype & (0x1<<1))
		{
			// write area ratio
			padding = sizeof(float) * 6;
			m_ofile->write(&padding, sizeof(uint), 1);
			m_ofile->write(area, sizeof(float), 6);
			m_ofile->write(&padding, sizeof(uint), 1);
		}

		if (m_datatype == 0 || m_datatype & (0x1<<2))
		{
			// write medium ID
			padding = sizeof(int);
			m_ofile->write(&padding, sizeof(uint), 1);
			m_ofile->write(&id, sizeof(int), 1);
			m_ofile->write(&padding, sizeof(uint), 1);
		}

		if (m_datatype == 0 || m_datatype & (0x1<<3))
		{
			// write cell boundary condition ID
		}

		if (m_datatype == 0 || m_datatype & (0x1<<4))
		{
			// write face boundary condition ID
		}
	}
	void doInterior(VxOctreeNode* node)
	{
		doCommon(node);
	}
	void doLeaf(VxOctreeNode* node)
	{
		doCommon(node);
	}

	TtOutFile* m_ofile;
	const VxVoxelizerOctree* m_vox;
	uint m_datatype;
};

int
ovx_save(const VxVoxelizerOctree* vox, const char* fname, const char* opt)
{
	tt_info("[ovx] saving... (%s)\n", fname);

	std::istringstream is_opt(opt);

	TtOutFile ofile(fname);

	uint  padding;
	uint  relm[3];
	uint  datatype;
	uint  null;
	uint  areatype;
	float bbox_min[3];
	float bbox_max[3];

	relm[0] = relm[1] = relm[2] = 1;
	datatype = 5;
	null = 0;
	areatype = 2;

	string token;
	while (is_opt >> token)
	{
		if (token == "type")
		{
			is_opt >> datatype;
		}
	}

	BBox bbox = vox->getBBox();
	vec3f_copy(bbox_min, bbox.min);
	vec3f_copy(bbox_max, bbox.max);

	// the number of root nodes
	padding = sizeof(uint) * 3;
	ofile.write(&padding, sizeof(uint), 1);
	ofile.write(relm, sizeof(uint), 3);
	ofile.write(&padding, sizeof(uint), 1);

	// data length
	padding = sizeof(uint) * 2;
	ofile.write(&padding, sizeof(uint), 1);
	ofile.write(&datatype, sizeof(uint), 1);
	ofile.write(&null, sizeof(uint), 1);
	ofile.write(&padding, sizeof(uint), 1);

	padding = sizeof(uint) * 7;
	ofile.write(&padding, sizeof(uint), 1);
	ofile.write(&areatype, sizeof(uint), 1);
	ofile.write(bbox_min, sizeof(float), 3);
	ofile.write(bbox_max, sizeof(float), 3);
	ofile.write(&padding, sizeof(uint), 1);

	VxOctree* octree = vox->getOctree();
	VxOctreeWriteOVX func(&ofile, vox, datatype);
	octree->traverse(octree->root(), &func);

	return 1;
}
