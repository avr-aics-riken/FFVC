/*
 * Polylib - Polygon Management Library.
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2010-
 *
 */
#ifndef c_polylib_h
#define c_polylib_h

#ifdef __cplusplus
extern "C" {  // for C++
#else
#endif

#include "common/PolylibStat.h"

///
/// C言語用Polylib（単一プロセス版）
///
/// 注意：C言語用Polylibでは、以下のC++版APIの代替機能はありません。
///            Polylib::set_factory()
///            Polylib::move()
///            Polylib::get_root_groups()
///

#define POLYLIB_FALSE 0
#define POLYLIB_TRUE  1

///
/// 三角形ポリゴン情報構造体
///
typedef struct {
	float m_vertex[9];	///< ３頂点座標
	float m_normal[3];	///< 法線ベクトル
	float m_area;		///< 面積
} TriangleStruct;

///
/// Polylib::loadメソッドのラッパー関数。
/// 引数で指定された設定ファイルを読み込み、グループツリーを作成する。
/// 続いて設定ファイルで指定されたSTLファイルを読み込み、KD木を作成する。
/// @param[in] fname 設定ファイル名。デフォルト値は、polylib_config.xml。
///  @return	POLYLIB_STATで定義される値が返る。
///
POLYLIB_STAT polylib_load(char* config_name);

///
/// Polylib::saveメソッドのラッパー関数。
/// PolygoGroupツリー、三角形ポリゴン情報の保存。
/// グループツリーの情報を設定ファイルへ出力。三角形ポリゴン情報をSTL
///	ファイルへ出力。三角形ポリゴンIDをIDファイルへ出力。
///
///  @param[out]    p_fname	設定ファイル名。
///  @param[in]     format	TriMeshIOクラスで定義されているSTLファイルの
///							フォーマット。
///  @param[in]     extend	ファイル名に付加する文字列。NULLを指定した
///							場合は、付加文字列として本メソッド呼び出し時の
///							年月日時分秒(YYYYMMDD24hhmmss)を用いる。
///  @return	POLYLIB_STATで定義される値が返る。
///  @attention	ファイル名命名規約は次の通り。
///			定義ファイル : polylib_config_付加文字.xml。
///			STLファイル  : ポリゴングループ名_付加文字.拡張子。
///			IDファイル   : ポリゴングループ名_付加文字.ID。
///
POLYLIB_STAT polylib_save(
	char	**p_fname,
	char	*format,
	char	*extend
);

///
/// Polylib::search_polygonsメソッドのラッパー関数。
/// 位置ベクトルmin_posとmax_posにより特定される矩形領域に含まれる、
/// 特定のグループとその子孫グループに属する三角形ポリゴンをKD探索に
/// より抽出する。Polylib内でメモリ領域が確保される。
///  @param[in]		group_name	抽出グループ名。
///  @param[in]		min_pos		抽出する矩形領域の最小値。(x,y,z順の配列)
///  @param[in]		max_pos		抽出する矩形領域の最大値。(x,y,z順の配列)
///  @param[in]		every		抽出オプション。
///   1：3頂点が全て検索領域に含まれるポリゴンを抽出する。
///   0：三角形のBBoxが一部でも検索領域と交差するものを抽出する。
///  @param[out]	num_tri		抽出された三角形ポリゴン数
///	 @param[out]	err			POLYLIB_STATで定義される値。通常はPLSTAT_OK。
///  @return 抽出された三角形ポリゴン構造体配列へのポインタを返す。
///  @attention 検索結果の三角形ポリゴン構造体はfreeしないでください。構造体を格納する
///             配列は呼び出し側でfreeして下さい。
///
TriangleStruct** polylib_search_polygons(
	char* group_name,
	float min_pos[3],
	float max_pos[3],
	int every, 
	int *num_tri,
	POLYLIB_STAT *err);

///
/// Polylib::show_group_hierarchyメソッドのラッパー関数。
/// グループ階層構造リストを標準出力に出力する。
///
void polylib_show_group_hierarchy();

///
/// Polylib::show_group_infoメソッドのラッパー関数。
/// グループの情報を出力する。(親グループ名、自信の名前、ファイル名、
///   登録三角形数、3頂点ベクトルの座標、法線ベクトルの座標、面積)
///  @param[in] group_name グループ名
///  @return POLYLIB_STATで定義される値が返る。
///
POLYLIB_STAT polylib_show_group_info(char* group_name);

///
///
///
unsigned int polylib_used_memory_size();

#ifdef __cplusplus
} // extern "C" or extern
#else
#endif

#endif //c_polylib_h
