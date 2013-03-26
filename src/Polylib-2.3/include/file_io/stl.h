/*
 * Polylib - Polygon Management Library.
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2010-
 *
 */
#ifndef stl_h
#define stl_h

#include <vector>
#include "common/PolylibCommon.h"

namespace PolylibNS {

///
/// ASCIIモードのSTLファイルを読み込み、tri_listに三角形ポリゴン情報を設定する。
///
///  @param[in,out] tri_list	三角形ポリゴンリストの領域。
///  @param[in]		fname		STLファイル名。
///  @param[in,out] total		ポリゴンIDの通番。
///  @return	POLYLIB_STATで定義される値が返る。
///
POLYLIB_STAT stl_a_load(
	std::vector<PrivateTriangle*>	*tri_list, 
	std::string 					fname,
	int								*total,
	float							scale=1.0
);

///
/// 三角形ポリゴン情報をASCIIモードでSTLファイルに書き出す。
///
///  @param[in] tri_list	三角形ポリゴン情報。
///  @param[in] fname		STLファイル名。
///  @return	POLYLIB_STATで定義される値が返る。
///
POLYLIB_STAT stl_a_save(
	std::vector<PrivateTriangle*>	*tri_list, 
	std::string 					fname
);

///
/// バイナリモードのSTLファイルを読み込み、tri_listに三角形ポリゴン情報を設定
/// する。
///
///  @param[in,out] tri_list	三角形ポリゴンリストの領域。
///  @param[in]		fname		ファイル名。
///  @param[in,out] total		ポリゴンIDの通番。
///  @return	POLYLIB_STATで定義される値が返る。
///
POLYLIB_STAT stl_b_load(
	std::vector<PrivateTriangle*>	*tri_list, 
	std::string						fname,
	int								*total,
	float							scale=1.0
);

///
/// 三角形ポリゴン情報をバイナリモードでSTLファイルに書き出す。
///
///  @param[in] tri_list	三角形ポリゴン情報。
///  @param[in] fname		STLファイル名。
///  @return	POLYLIB_STATで定義される値が返る。
///
POLYLIB_STAT stl_b_save(
	std::vector<PrivateTriangle*>	*tri_list, 
	std::string						fname
);

///
/// STLファイルを読み込みバイナリかアスキーかを判定する。
///
///  @param[in] STLファイルのフルパス名。
///  @return	true:アスキー形式 / false:バイナリ形式。
/// 
bool is_stl_a(
	std::string		path
);

///
/// STLファイル名から名称(拡張子を除いた部分)を取得する。
///
///  @param[in] STLファイルのフルパス名。
///  @return	拡張子を除いた名称。
///  @attention	戻り値のchar *は解放不要。
///
char *stl_get_fname(
	std::string		path
);

///
/// STLファイル名から拡張子のみを取得する。
///
///  @param[in] STLファイルのフルパス名。
///  @return	拡張子。
///  @attention	戻り値のchar *は解放不要。
///
char *stl_get_ext(
	std::string		path
);

} //namespace PolylibNS

#endif  // stl_h

