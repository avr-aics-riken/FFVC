#include <stdlib.h>
#include <string>
#include <sstream>
#include "app.h"

using namespace std;

int
main(int argc, char *argv[])
{
	string ifname_obj = "box.stla";
	string ofname_stl = "out2.stl";
	string ofname_svx = "out2.svx";
	string ofname_ovx = "out2.ovx";
	string ofname_sph = "out2.sph";
	Vec3f orig, pitch;
	Vec3i nelm;

	Application app;

	// load mesh data
	app.loadObject(ifname_obj);

	// clean mesh data
	app.fixObject(0, 0, 30);

	// save mesh data
	app.saveObject(ofname_stl);

	orig = Vec3f(1, -.5, -.5);
	pitch = Vec3f(0.01, 0.01, 0.01);
	nelm = Vec3i(200, 200, 200);

	// voxelize mesh data
	app.voxelizePolygon(orig, pitch, nelm);

	// set a ratio type to a binary mode
	app.setRatioType(0);

	// save voxel data in svx format
	// it is possible to set options
	app.saveVoxel(ofname_svx, "type 7");

	// save voxel data in ovx format
	app.saveVoxel(ofname_ovx, "type 7");

	// save voxel data in sph format
	app.saveVoxel(ofname_sph);

	return 0;
}

