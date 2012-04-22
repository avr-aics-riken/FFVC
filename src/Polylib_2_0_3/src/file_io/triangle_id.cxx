/*
 * Polylib - Polygon Management Library.
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2010-
 *
 */
#include <fstream>
#include <vector>
#include <iomanip>
#include "polygons/Triangle.h"
#include "file_io/triangle_id.h"

#define AS_BINARY	0
#define AS_ASCII	1

namespace PolylibNS {

using namespace std;

//////////////////////////////////////////////////////////////////////////////
// 変更:ポリゴンIDのバイナリ入力対応 2010.10.19
POLYLIB_STAT load_id(
	vector<PrivateTriangle*>	*tri_list, 
	string 						fname,
	ID_FORMAT					id_format
) {
	ifstream	is;

	if (id_format == ID_BIN)	is.open(fname.c_str(), ios::binary);
	else						is.open(fname.c_str());

	if (is.fail()) {
		PL_ERROSH << "[ERROR]triangle_id:load_id():Can't open " << fname << endl;
		return PLSTAT_STL_IO_ERROR;
	}

	vector<PrivateTriangle*>::iterator itr = tri_list->begin();

	if (id_format == ID_BIN) {
		int		id;
		while (is.read((char*)&id, sizeof(int)) && !is.eof()) {
			// ポリゴン数とIDファイルの行数が一致しているか?
			if (itr == tri_list->end()) {
				PL_ERROSH << "[ERROR]triangle_id::load_id():Triangle number "
						  << "is short:" << fname << endl;
				return PLSTAT_STL_IO_ERROR;
			}
			(*itr)->set_id(id);
			itr++;
		}
	}
	else {
		string	id;
		while (is >> id && !is.eof()) {
			// ポリゴン数とIDファイルの行数が一致しているか?
			if (itr == tri_list->end()) {
				PL_ERROSH << "[ERROR]triangle_id::load_id():Triangle number "
						  << "is short:" << fname << endl;
				return PLSTAT_STL_IO_ERROR;
			}
			(*itr)->set_id(atoi(id.c_str()));
			itr++;
		}
	}

	// ポリゴン数とIDファイルの行数が一致しているか?
	if (itr != tri_list->end()) {
		PL_ERROSH << "[ERROR]triangle_id::load_id():ID number is short:" 
				  << fname << endl;
		return PLSTAT_STL_IO_ERROR;
	}

	if (!is.eof() && is.fail()) {
		PL_ERROSH << "[ERROR]triangle_id:load_id():Error in loading:" 
				  << fname << endl;
		return PLSTAT_STL_IO_ERROR;
	}

	return PLSTAT_OK;
}

//////////////////////////////////////////////////////////////////////////////
// 変更:ポリゴンIDのバイナリ出力対応 2010.10.19
POLYLIB_STAT save_id(
	vector<PrivateTriangle*>	*tri_list, 
	string 						fname,
	ID_FORMAT					id_format
) {
	ofstream os;

	if (id_format == ID_BIN)	os.open(fname.c_str(), ios::binary);
	else						os.open(fname.c_str());

	if (os.fail()) {
		PL_ERROSH << "[ERROR]triangle_id:save_id():Can't open " << fname << endl;
		return PLSTAT_STL_IO_ERROR;
	}

	vector<PrivateTriangle*>::iterator itr;
	if (id_format == ID_BIN) {
		for (itr = tri_list->begin(); itr != tri_list->end(); itr++) {
			int		id = (*itr)->get_id();
			os.write((char *)&id, sizeof(int));
		}
	}
	else {
		for (itr = tri_list->begin(); itr != tri_list->end(); itr++) {
			os	<< (*itr)->get_id() << endl;
		}
	}

	if (!os.eof() && os.fail()) {
		PL_ERROSH << "[ERROR]triangle_id:save_id():Error in saving:" 
				  << fname << endl;
		return PLSTAT_STL_IO_ERROR;
	}

	return PLSTAT_OK;
}

} //namespace PolylibNS
