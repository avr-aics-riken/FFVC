#ifndef pedig2_h
#define pedig2_h

#define PDG_OTHER_MASK 0xfff0000000000000LL
#define PDG_LEVEL_MASK 0x000f000000000000LL
#define PDG_PEDIG_MASK 0x0000ffffffffffffLL
#define PDG_LEVEL_SHIFT 48

typedef unsigned long long VxOctreePedigree;

struct VxOctreePosi
{
	VxOctreePosi() : x(0), y(0), z(0), w(0) {}
	unsigned short x, y, z, w;
};

//=========================================================================
// VxOctreePedigree
//=========================================================================

inline void pdg_push(VxOctreePedigree& pedig, int oct)
{
	VxOctreePedigree other = pedig & PDG_OTHER_MASK;
	VxOctreePedigree lv = (pedig & PDG_LEVEL_MASK) >> PDG_LEVEL_SHIFT;
	VxOctreePedigree pd = pedig & PDG_PEDIG_MASK;
	lv++;
	pd = (pd << 3) + oct;
	pedig = other | lv<<PDG_LEVEL_SHIFT | pd;
}

inline void pdg_pop(VxOctreePedigree& pedig)
{
	VxOctreePedigree other = pedig & PDG_OTHER_MASK;
	VxOctreePedigree lv = (pedig & PDG_LEVEL_MASK) >> PDG_LEVEL_SHIFT;
	VxOctreePedigree pd = pedig & PDG_PEDIG_MASK;
	lv--;
	pd = pd >> 3;
	pedig = other | lv<<PDG_LEVEL_SHIFT | pd;
}

inline int pdg_top(const VxOctreePedigree& pedig)
{
	return pedig & 0x7;
}

inline int pdg_get_depth(const VxOctreePedigree& pedig)
{
	return (pedig & PDG_LEVEL_MASK) >> PDG_LEVEL_SHIFT;
}

inline void pdg_pedigree2xyz(VxOctreePedigree pedig, int* x, int* y, int* z, int* depth)
{
	int i, lv, oct;

	lv = pdg_get_depth(pedig);
	*x = *y = *z = 0;
	for (i=0; i<lv; i++)
	{
		oct = pdg_top(pedig);
		if (oct & 0x1) *x += (0x1<<i);
		if (oct & 0x2) *y += (0x1<<i);
		if (oct & 0x4) *z += (0x1<<i);
		pdg_pop(pedig);
	}
	*depth = lv;
}

inline void pdg_xyz2pedigree(VxOctreePedigree* _pedig, int x, int y, int z, int level)
{
	VxOctreePedigree pedig = 0;
	int i, oct;

	for (i=level-1; i>=0; i--)
	{
		oct = 0;
		oct |=  (x>>i) & 0x1;
		oct |= ((y>>i) & 0x1) << 1;
		oct |= ((z>>i) & 0x1) << 2;
		pdg_push(pedig, oct);
	}
	*_pedig = pedig;
}

//=========================================================================
// VxOctreePosi
//=========================================================================

inline void pdg_push(VxOctreePosi& posi, int oct)
{
	posi.x = (posi.x << 1) | (oct    & 0x1);
	posi.y = (posi.y << 1) | (oct>>1 & 0x1);
	posi.z = (posi.z << 1) | (oct>>2 & 0x1);
	unsigned short other = posi.w & 0xfff0;
	unsigned short lv = posi.w & 0xf;
	lv++;
	posi.w = other | lv;
}

inline int pdg_top(const VxOctreePosi& posi)
{
	int oct = 0;
	oct |=  posi.x & 0x1;
	oct |= (posi.y & 0x1) << 1;
	oct |= (posi.z & 0x1) << 2;
	return oct;
}

inline int pdg_get_depth(const VxOctreePosi& posi)
{
	return posi.w & 0xf;
}

inline void pdg_posi2xyz(const VxOctreePosi& posi, int* x, int* y, int* z, int* depth)
{
	*x = posi.x;
	*y = posi.y;
	*z = posi.z;
	*depth = pdg_get_depth(posi);
}

#endif  // pedig2_h

