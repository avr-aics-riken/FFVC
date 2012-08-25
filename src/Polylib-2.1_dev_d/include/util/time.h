/*
 * Polylib - Polygon Management Library.
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2010-
 *
 */
#ifndef polylib_time_h
#define polylib_time_h

namespace PolylibNS {
bool getrusage_sec(
	double	*usr_time, 
	double	*sys_time, 
	double	*total
);
}

#endif //polylib_time_h
