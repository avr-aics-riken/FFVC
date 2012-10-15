/*
 * Polylib - Polygon Management Library.
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2010-
 *
 */
#ifndef polylib_stat_h
#define polylib_stat_h

////////////////////////////////////////////////////////////////////////////
///
/// Polylibで利用するEnumの定義
///
////////////////////////////////////////////////////////////////////////////
typedef enum {
	PLSTAT_OK,					///< 処理が成功した。
	PLSTAT_NG,					///< 一般的なエラー。
	PLSTAT_INSTANCE_EXISTED,	///< Polylibインスタンスがすでに存在している。
	PLSTAT_INSTANCE_NOT_EXIST,	///< Polylibインスタンスが存在しない。
	PLSTAT_STL_IO_ERROR,		///< STLファイルIOエラー
	PLSTAT_UNKNOWN_STL_FORMAT,	///< ファイル拡張子が.stla、.stlb、.stl以外。
	PLSTAT_FILE_NOT_SET,		///< リーフグループにファイル名が未設定。
	PLSTAT_CONFIG_ERROR,		///< 定義ファイルでエラー発生
	PLSTAT_GROUP_NOT_FOUND,		///< グループ名がPolylibに未登録。
	PLSTAT_GROUP_NAME_EMPTY,	///< グループ名が空である。
	PLSTAT_GROUP_NAME_DUP,		///< グループ名が重複している。
	PLSTAT_MEMORY_NOT_ALLOC,	///< メモリ確保に失敗した。
	PLSTAT_POLYGON_NOT_EXIST,	///< PolygonGroupにPolygonsが未設定。
	PLSTAT_TRIANGLE_NOT_EXIST,	///< Polygonsにポリゴンリストが未設定。
	PLSTAT_NODE_NOT_FIND,		///< KD木生成時に検索点が見つからなかった。
	PLSTAT_ROOT_NODE_NOT_EXIST,	///< KD木のルートノードが存在しない。
	PLSTAT_ARGUMENT_NULL,		///< 引数のメモリ確保が行われていない。
	PLSTAT_MPI_ERROR,			///< MPI関数がエラーを戻した。
// 以下は未使用
//	PLSTAT_GROUP_UNMATCH,		///< グループ並びがランク0と一致しなかった。
//	PLSTAT_UNkNOWN_ERROR,		///< 予期せぬエラー。
//	PLSTAT_FILE_NOT_FOUND,		///< 読み込みファイルが存在しない。
//	PLSTAT_FILE_NOT_OPEN,		///< ファイルが開けなかった。
//	PLSTAT_UNKNOWN_FILE_FORMAT,	///< stl_a、stl_b以外が指定された。
//	PLSTAT_FILE_NOT_SCAN,		///< ファイル読み込みエラー。
//	PLSTAT_PARAMETER_NOT_FOUND,	///< ConfigファイルにParamタグがない。
} POLYLIB_STAT;

#ifndef C_LANG // C言語版でも本ヘッダを使用しているため#ifndefを追加 2010.11.04
namespace PolylibNS {
////////////////////////////////////////////////////////////////////////////
///
/// PolylibStat文字列出力用クラス
///
////////////////////////////////////////////////////////////////////////////
class PolylibStat2 {
public:
///
/// PolylibStat文字列出力。
///
///  @param[in] stat	PolylibStat値。
///  @return    PolylibStat値を文字列化したもの。
///
static std::string String(
	POLYLIB_STAT	stat
) {
	if (stat == PLSTAT_OK)							return "PLSTAT_OK";
	else if (stat == PLSTAT_NG)						return "PLSTAT_NG";
	else if (stat == PLSTAT_INSTANCE_EXISTED) 		return "PLSTAT_INSTANCE_EXISTED";
	else if (stat == PLSTAT_INSTANCE_NOT_EXIST) 	return "PLSTAT_INSTANCE_NOT_EXIST";
	else if (stat == PLSTAT_STL_IO_ERROR) 			return "PLSTAT_STL_IO_ERROR";
	else if (stat == PLSTAT_UNKNOWN_STL_FORMAT) 	return "PLSTAT_UNKNOWN_STL_FORMAT";
	else if (stat == PLSTAT_FILE_NOT_SET) 			return "PLSTAT_FILE_NOT_SET";
	else if (stat == PLSTAT_CONFIG_ERROR) 			return "PLSTAT_CONFIG_ERROR";
	else if (stat == PLSTAT_GROUP_NOT_FOUND) 		return "PLSTAT_GROUP_NOT_FOUND";
	else if (stat == PLSTAT_GROUP_NAME_EMPTY) 		return "PLSTAT_GROUP_NAME_EMPTY";
	else if (stat == PLSTAT_GROUP_NAME_DUP) 		return "PLSTAT_GROUP_NAME_DUP";
	else if (stat == PLSTAT_MEMORY_NOT_ALLOC) 		return "PLSTAT_MEMORY_NOT_ALLOC";
	else if (stat == PLSTAT_POLYGON_NOT_EXIST) 		return "PLSTAT_POLYGON_NOT_EXIST";
	else if (stat == PLSTAT_TRIANGLE_NOT_EXIST) 	return "PLSTAT_TRIANGLE_NOT_EXIST";
	else if (stat == PLSTAT_NODE_NOT_FIND) 			return "PLSTAT_NODE_NOT_FIND";
	else if (stat == PLSTAT_ROOT_NODE_NOT_EXIST) 	return "PLSTAT_ROOT_NODE_NOT_EXIST";
	else if (stat == PLSTAT_ARGUMENT_NULL) 			return "PLSTAT_ARGUMENT_NULL";
	else if (stat == PLSTAT_MPI_ERROR) 				return "PLSTAT_MPI_ERROR";
	else											return "UNKNOW_STATUS";
}
};
}
#endif // C_LANG

#endif // polylib_stat_h
