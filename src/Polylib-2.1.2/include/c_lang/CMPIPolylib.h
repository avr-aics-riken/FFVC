/*
 * Polylib - Polygon Management Library.
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2010-
 *
 */
#ifndef c_mpi_polylib_h
#define c_mpi_polylib_h

#ifdef __cplusplus
extern "C" {  // for C++
#else
#endif

#include "mpi.h"
#include "common/PolylibStat.h"
#include "c_lang/CPolylib.h"

///
/// C言語用MPIPolylib（MPI版）
///
/// 注意：C言語用MPIPolylibでは、以下のC++版APIの代替機能はありません。
///            MPIPolylib::set_factory()
///            MPIPolylib::move()
///            MPIPolylib::migrate()
///            MPIPolylib::get_root_groups()
///

///
/// MPIPolylib::init_parallel_infoメソッドのラッパー関数。
/// 並列計算関連情報の設定と初期化を行う。
/// 全rankで各々設定を行い、その領域情報を全rankへ配信する。
///
///  @param[in] comm	MPIコミュニケーター
///  @param[in] bpos	自PE担当領域の基点座標
///  @param[in] bbsize	同、計算領域のボクセル数
///  @param[in] gcsize	同、ガイドセルのボクセル数
///  @param[in] dx		同、ボクセル１辺の長さ
///  @return POLYLIB_STATで定義される値が返る。
///
POLYLIB_STAT
mpipolylib_init_parallel_info(
	MPI_Comm comm,
	float bpos[3],
	unsigned int bbsize[3],
	unsigned int gcsize[3],
	float dx[3]
);

///
/// MPIPolylib::load_rank0メソッドのラッパー関数。
/// rank0によるデータ構築。
/// 指定された設定ファイルをrank0にて読み込み、グループ階層構造の構築
/// およびポリゴンデータの構築を行う。
/// グループ階層構造は全rankにb_castされ、情報を共有する。
/// ポリゴンデータは各rank領域毎のデータが分配される。
/// @param[in] fname 設定ファイル名。NULLならば、polylib_config.xml。
///  @return	POLYLIB_STATで定義される値が返る。
///
POLYLIB_STAT
mpipolylib_load_rank0(char* config_name);

///
/// MPIPolylib::load_parallelメソッドのラッパー関数。
/// 全rank並列でのデータ構築。
/// 指定された設定ファイルを各rankにて読み込み、グループ階層構造の構築、
/// およびポリゴンデータの構築を行う。
/// @attension 各rankが読み込むファイルに記述されたグループ階層構造が一致している必要がある。
///
/// @param[in] config_filename	初期化ファイル名。NULLならば、polylib_config.xmlを読む。
/// @return	POLYLIB_STATで定義される値が返る。
///
POLYLIB_STAT
mpipolylib_load_parallel(char* config_name);

///
/// MPIPolylib::save_rank0メソッドのラッパー関数。
/// rank0によるデータ保存。
/// rank0のMPIPolylibインスタンスが保持するグループ階層構造を設定ファイルに書き出す。
/// 同時に各rankに分散するポリゴンデータもrank0に集められ、指定されたフォーマットの
/// STLファイルにrank0で書き出す。
///
///  @param[out]    p_fname	設定ファイル名返却用
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
POLYLIB_STAT
mpipolylib_save_rank0(
	char	**p_fname,
	char	*format,
	char	*extend
);

///
/// MPIPolylib::save_parallelメソッドのラッパー関数。
/// 全rank並列でのデータ保存。
/// 各rankのMPIPolylibインスタンスが保持するグループ階層構造を設定ファイルに各rank毎に書き出す。
/// 同時にポリゴンデータも指定されたフォーマットのSTLファイルに各rank毎に書き出す。
///
///  @param[out]    p_fname	設定ファイル名返却用
///  @param[in]     format	TriMeshIOクラスで定義されているSTLファイルの
///							フォーマット。
///  @param[in]     extend	ファイル名に付加する文字列。NULLを指定した
///							場合は、付加文字列として本メソッド呼び出し時の
///							年月日時分秒(YYYYMMDD24hhmmss)を用いる。
///  @return	POLYLIB_STATで定義される値が返る。
///  @attention	ファイル名命名規約は次の通り。
///			定義ファイル : polylib_config_ランク番号_付加文字.xml。
///			STLファイル  : ポリゴングループ名_ランク番号_付加文字.拡張子。
///			IDファイル   : ポリゴングループ名_ランク番号_付加文字.ID。
///
POLYLIB_STAT
mpipolylib_save_parallel(
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
TriangleStruct** mpipolylib_search_polygons(
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
void mpipolylib_show_group_hierarchy();

///
/// Polylib::show_group_infoメソッドのラッパー関数。
/// グループの情報を出力する。(親グループ名、自信の名前、ファイル名、
///   登録三角形数、3頂点ベクトルの座標、法線ベクトルの座標、面積)
///  @param[in] group_name グループ名
///  @return POLYLIB_STATで定義される値が返る。
///
POLYLIB_STAT mpipolylib_show_group_info(char* group_name);

///
/// MPIPolylibが利用中の概算メモリ量を返す
///
/// @return 利用中のメモリ量(byte)
///
unsigned int mpipolylib_used_memory_size();




#ifdef __cplusplus
} // extern "C" or extern
#else
#endif

#endif //cmpi_polylib_h
