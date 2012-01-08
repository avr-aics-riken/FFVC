#include <math.h>
#include <stdlib.h>
#include "table.h"

namespace gm {
	// [0, pi)
	float cos_pi [256];
	float sin_pi [256];
	// [0, 2pi)
	float cos_2pi[256];
	float sin_2pi[256];
}

namespace tt
{
	float g_colortab[256][3];
}

void
init_trifunc_table() {
	for (int i=0; i<256; i++) {
		float radian = M_PI * i / 256;
		gm::cos_pi [i] = cosf(radian);
		gm::sin_pi [i] = sinf(radian);
		gm::cos_2pi[i] = cosf(radian * 2);
		gm::sin_2pi[i] = sinf(radian * 2);
	}
}

void
init_color_table() {
	srand(1);
	for (int i=0; i<256; i++) {
		tt::g_colortab[i][0] = (float)rand() / (RAND_MAX + 1.);
		tt::g_colortab[i][1] = (float)rand() / (RAND_MAX + 1.);
		tt::g_colortab[i][2] = (float)rand() / (RAND_MAX + 1.);
	}
}

class MakeLookupTable {
public:
	MakeLookupTable() {
		init_trifunc_table();
		init_color_table();
	}
};

MakeLookupTable lookup_table;

