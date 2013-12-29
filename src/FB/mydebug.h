#ifndef _FB_DEBUG_H_
#define _FB_DEBUG_H_

//##################################################################################
//
// Flow Base class
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2013 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2013 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

//@file   mydebug.h
//@brief  FlowBase mydebug Header
//@author kero

#include <stdlib.h>

#ifdef DEBUG
#define debug_message() printf("%s (%d):\n",__FILE__, __LINE__)
#else
#define debug_message()
#endif

#define Exit(x) \
((void)printf("exit at %s:%u\n", __FILE__, __LINE__), exit((x)))

#define message() printf("\t%s (%d):\n",__FILE__, __LINE__)
#define mark() printf("%s (%d) [%d]:\n",__FILE__, __LINE__, paraMngr->GetMyRankID())

#define stamped_printf printf("%s (%d):  ",__FILE__, __LINE__), printf
#define stamped_fprintf fprintf(fp, "%s (%d):  ",__FILE__, __LINE__), fprintf

#define Hostonly_ if(paraMngr->GetMyRankID()==0) 

#define TIMING__ if ( ModeTiming == ON )

#endif // _FB_DEBUG_H_
