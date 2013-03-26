/*
 * Polylib - Polygon Management Library.
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2010-
 *
 */
#ifndef mesh_io_h
#define mesh_io_h

#include <vector>
#include <map>
#include "common/PolylibStat.h"
#include "common/PolylibCommon.h"

namespace PolylibNS {

////////////////////////////////////////////////////////////////////////////
///
/// クラス:TriMeshIO
/// 三角形ポリゴン入出力管理。
///
////////////////////////////////////////////////////////////////////////////
class TriMeshIO {
public:
	///
	/// STLファイルを読み込み、tri_listにセットする。
	///
	///  @param[in,out] tri_list	三角形ポリゴンリストの領域。
	///  @param[in]		fmap		ファイル名、ファイルフォーマットのセット。
	///  @return	POLYLIB_STATで定義される値が返る。
	///
	static POLYLIB_STAT load(
		std::vector<PrivateTriangle*>				*tri_list,
		const std::map<std::string, std::string>	&fmap,
		float scale = 1.0
	);

	///
	/// tri_listの内容をSTL形式でファイルへ保存。
	///
	///  @param[in] tri_list	三角形ポリゴンのリスト(出力内容)。
	///  @param[in] fname		ファイル名。
	///  @param[in] fmt			ファイルフォーマット。
	///  @return	POLYLIB_STATで定義される値が返る。
	///
	static POLYLIB_STAT save(
		std::vector<PrivateTriangle*>	*tri_list,
		std::string						fname, 
		std::string 					fmt = ""
	);

	///
	/// ファイル名を元に入力ファイルのフォーマットを取得する。
	///
	///  @param[in] filename		入力ファイル名。
	///  @return	判定したファイルフォーマット。
	///  @attention	ファイル拡張子が"stl"の場合、ファイルを読み込んで判定する。
	///
	static std::string input_file_format(
		const std::string &filename
	);

	/// STLファイルのフォーマット種別
	///
	///  @attention STLファイルの拡張子とは異なるので注意すること。
	///
	static const std::string FMT_STL_A;		///< アスキーファイル
	static const std::string FMT_STL_AA;	///< アスキーファイル
	static const std::string FMT_STL_B;		///< バイナリファイル
	static const std::string FMT_STL_BB;	///< バイナリファイル
	static const std::string DEFAULT_FMT;	///< TrimeshIO.cxxで定義している値
};

} //namespace PolylibNS

#endif  // mesh_io_h
