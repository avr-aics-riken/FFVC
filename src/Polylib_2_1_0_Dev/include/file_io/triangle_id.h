/*
 * Polylib - Polygon Management Library.
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2010-
 *
 */
#ifndef stl_h
#define stl_h

#include <vector>
#include "polygons/Triangle.h"
#include "common/PolylibStat.h"
#include "common/PolylibCommon.h"

namespace PolylibNS {

///
/// 三角形ポリゴンIDをIDファイルから読み込み、m_idにセットする。
///
///  @param[in,out] tri_list	三角形ポリゴン情報。
///  @param[in]		fname		三角形ポリゴンIDファイル名。
///  @param[in]		id_format	三角形ポリゴンIDファイルの入力形式。
///  @return	POLYLIB_STATで定義される値が返る
///
POLYLIB_STAT load_id(
	std::vector<PrivateTriangle*>	*tri_list, 
	std::string 					fname,
	ID_FORMAT						id_format
);

///
/// 三角形ポリゴンIDをファイルに出力する。
///
///  @param[in] tri_list	三角形ポリゴン情報。
///  @param[in] fname		三角形ポリゴンIDファイル名。
///  @param[in]	id_format	三角形ポリゴンIDファイルの出力形式。
///  @return	POLYLIB_STATで定義される値が返る
///
POLYLIB_STAT save_id(
	std::vector<PrivateTriangle*>	*tri_list, 
	std::string 					fname,
	ID_FORMAT						id_format
);

} //namespace PolylibNS

#endif  // stl_h

