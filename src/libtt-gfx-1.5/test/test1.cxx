#include <stdlib.h>
#include <string>
#include <sstream>
#include "app.h"

using namespace std;

int
main(int argc, char *argv[])
{
	string ifname_obj = "flamingo.obj";
	string ofname_stl = "out.stl";
	string ofname_svx = "out.svx";
	string ofname_ovx = "out.ovx";
	string ofname_sph = "out.sph";
	Vec3f orig, pitch;
	Vec3i nelm;

	Application app;

	// load mesh data
	app.loadObject(ifname_obj);

	// clean mesh data
	app.fixObject(0, 0, 30);

	// save mesh data
	app.saveObject(ofname_stl);

	// compute voxel parameters
	// ex.) 100 voxels along the y axis
	string s = app.computeVoxelParam(1, 100);
	cout << "voxel param: " << s << endl;

	istringstream is(s);
	is >> orig >> pitch >> nelm;

	// voxelize mesh data
	app.voxelizePolygon(orig, pitch, nelm);

	// save voxel data in svx format
	// it is possible to set options
	app.saveVoxel(ofname_svx, "type 7 nelm 50 100 100");

	// save voxel data in ovx format
	app.saveVoxel(ofname_ovx, "type 7");

	// save voxel data in sph format
	app.saveVoxel(ofname_sph, "nelm 50 100 100");

	return 0;
}

