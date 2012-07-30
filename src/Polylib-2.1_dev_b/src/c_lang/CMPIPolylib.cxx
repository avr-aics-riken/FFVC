/*
 * Polylib - Polygon Management Library.
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2010-
 *
 */
#include <iostream>
#include <vector>
#include "common/PolylibCommon.h"
#include "common/PolylibStat.h"
#include "MPIPolylib.h"
#include "c_lang/CPolylib.h"
#include "c_lang/CMPIPolylib.h"

///
/// C言語用MPIPolylib-API（MPI版）
///

using namespace std;
using namespace PolylibNS;

// init_parallel_info
POLYLIB_STAT
mpipolylib_init_parallel_info(
	MPI_Comm comm,
	float bpos[3],
	unsigned int bbsize[3],
	unsigned int gcsize[3],
	float dx[3]
)
{
	return (MPIPolylib::get_instance())->init_parallel_info( comm, bpos, bbsize, gcsize, dx );
}



// load_rank0
POLYLIB_STAT
mpipolylib_load_rank0(char* config_name)
{
	if( config_name == NULL ) {
		return (MPIPolylib::get_instance())->load_rank0();
	}
	string fname = config_name;
	return (MPIPolylib::get_instance())->load_rank0( fname );
}



// load_parallel
POLYLIB_STAT
mpipolylib_load_parallel(char* config_name)
{
	if( config_name == NULL ) {
		return (MPIPolylib::get_instance())->load_parallel();
	}
	string fname = config_name;
	return (MPIPolylib::get_instance())->load_parallel( fname );
}



// save_rank0
POLYLIB_STAT
mpipolylib_save_rank0(
	char	**p_fname,
	char	*format,
	char	*extend
)
{
	string s_fname;
	string s_format = format;
	string s_extend;
	POLYLIB_STAT stat;

	if( extend ) s_extend = extend;

	if( extend==NULL ) {
		stat = (MPIPolylib::get_instance())->save_rank0( &s_fname, s_format );
	}
	else {
		stat = (MPIPolylib::get_instance())->save_rank0( &s_fname, s_format, s_extend );
	}
	*p_fname = (char*)malloc( s_fname.size()+1 );
	if(p_fname == NULL){
		fprintf(stderr,"mpipolylib_save_rank0: Can not allocate memory.\n");
		return PLSTAT_MEMORY_NOT_ALLOC;
	}
	strcpy( *p_fname, s_fname.c_str() );
	return stat;
}



// save_parallel
POLYLIB_STAT
mpipolylib_save_parallel(
	char	**p_fname,
	char	*format,
	char	*extend
)
{
	string s_fname;
	string s_format = format;
	string s_extend;
	POLYLIB_STAT stat;

	if( extend ) s_extend = extend;

	if( extend==NULL ) {
		stat = (MPIPolylib::get_instance())->save_parallel( &s_fname, s_format );
	}
	else {
		stat = (MPIPolylib::get_instance())->save_parallel( &s_fname, s_format, s_extend );
	}
	*p_fname = (char*)malloc( s_fname.size()+1 );
	if(p_fname == NULL){
		fprintf(stderr,"mpipolylib_save_parallel: Can not allocate memory.\n");
		return PLSTAT_MEMORY_NOT_ALLOC;
	}
	strcpy( *p_fname, s_fname.c_str() );
	return stat;
}



// search_polygons
TriangleStruct** mpipolylib_search_polygons(
	char* group_name,
	float min_pos[3],
	float max_pos[3],
	int every, 
	int *num_tri,
	POLYLIB_STAT *err
)
{
	Vec3f c_min_pos,c_max_pos;
	string c_group_name(group_name);

	for(int i=0; i<3; i++){
		c_min_pos[i] = min_pos[i];
		c_max_pos[i] = max_pos[i];
	}

	bool b_every;
	if(every == POLYLIB_TRUE) b_every = true;
	else b_every = false;

	//Polylibから三角形リストを抽出
	std::vector<Triangle*>*  tri_list =
			(MPIPolylib::get_instance())->search_polygons(
						c_group_name, c_min_pos, c_max_pos, b_every);
	*num_tri  = tri_list->size();

	//三角形リストのポインタ配列の確保
	TriangleStruct **p_tri =
			(TriangleStruct**)malloc( sizeof(TriangleStruct*)*(*num_tri) );
	if(p_tri == NULL){
		fprintf(stderr,"mpipolylib_serch_polygons: Can not allocate memory.\n");
		*err = PLSTAT_MEMORY_NOT_ALLOC;
		return NULL;
	}

	std::vector<Triangle*>::iterator itr;
	int num = 0;
	for(itr=tri_list->begin(); itr!=tri_list->end(); itr++){
		//TriangleクラスインスタンスをTriangleStruct*にキャストする
		p_tri[num] = (TriangleStruct*)( *itr );
		num++;
	}
	delete(tri_list);
	*err = PLSTAT_OK;
	return p_tri;
}



// show_group_hierarchy
void mpipolylib_show_group_hierarchy()
{
	(MPIPolylib::get_instance())->show_group_hierarchy();
}



// show_group_info
POLYLIB_STAT mpipolylib_show_group_info(char* group_name)
{
	string c_group_name(group_name);
	return (MPIPolylib::get_instance())->show_group_info(c_group_name);
}



// used_memory_size
unsigned int mpipolylib_used_memory_size()
{
	return (MPIPolylib::get_instance())->used_memory_size();
}

// eof
