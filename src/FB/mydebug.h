#ifndef _SKL_SOLVER_DEBUG_H_
#define _SKL_SOLVER_DEBUG_H_

/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file mydebug.h
//@brief FlowBase mydebug Header
//@author keno, FSI Team, VCAD, RIKEN

#include <stdlib.h>

#ifdef DEBUG
#define debug_message() printf("%s (%d):\n",__FILE__, __LINE__)
#else
#define debug_message()
#endif

#define Exit(x) \
((void)printf("exit at %s:%u\n", __FILE__, __LINE__), exit((x)))

#define message() printf("\t%s (%d):\n",__FILE__, __LINE__)
#define mark() printf("%s (%d):\n",__FILE__, __LINE__)

#define stamped_printf printf("%s (%d):  ",__FILE__, __LINE__), printf
#define stamped_fprintf fprintf(fp, "%s (%d):  ",__FILE__, __LINE__), fprintf
#define stamped_mprintf fprintf(mp, "%s (%d):  ",__FILE__, __LINE__), fprintf

#define Hostonly_ if(pn.ID==0)

#define TIMING__ if ( ModeTiming == ON )

#endif // _SKL_SOLVER_DEBUG_H_
