/*
 * Polylib - Polygon Management Library.
 * Version     : 2.0.3
 * Release date: Nov.29.2010
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2010-
 *
 */
#ifndef polylib_h
#define polylib_h

#include <vector>
#include <iostream>
#include "polygons/Polygons.h"
#include "polygons/TriMesh.h"
#include "polygons/Triangle.h"
#include "groups/PolygonGroup.h"
#include "groups/PolygonGroupFactory.h"
#include "file_io/PolylibConfig.h"
#include "common/PolylibStat.h"
#include "common/PolylibCommon.h"
#include "common/BBox.h"
#include "common/Vec3.h"

namespace PolylibNS {

////////////////////////////////////////////////////////////////////////////
///
/// クラス:CalcAreaInfo
/// 計算領域情報。
///
////////////////////////////////////////////////////////////////////////////
struct CalcAreaInfo {
	/// 基点座標
	Vec3f m_bpos;

	/// 計算領域のボクセル数
	Vec3f m_bbsize;

	/// ガイドセルのボクセル数
	Vec3f m_gcsize;

	/// ボクセル１辺の長さ
	Vec3f m_dx;

	/// ガイドセルを含めた担当領域の最小位置
	Vec3f m_gcell_min;

	/// ガイドセルを含めた担当領域の最大位置
	Vec3f m_gcell_max;

	/// ガイドセルを含めたBounding Box
	BBox m_gcell_bbox;
};

////////////////////////////////////////////////////////////////////////////
///
/// クラス:PolylibMoveParams
/// Polylib::move()の引数として利用するパラメタセットクラスです。
/// 本クラスメンバ変数ではパラメタが不足する場合は、継承クラスをユーザ定義
/// してください。
///
////////////////////////////////////////////////////////////////////////////
class PolylibMoveParams {
public:
	/// 現在の計算ステップ番号
	int	m_current_step;

	/// 移動後の計算ステップ番号
	int m_next_step;

	/// １計算ステップあたりの時間変異
	double m_delta_t;
};

////////////////////////////////////////////////////////////////////////////
///
/// クラス:Polylib
///	ポリゴンを管理する為のクラスライブラリです。
///
////////////////////////////////////////////////////////////////////////////
class Polylib
{
public:
	///
	/// singletonのPolylibインスタンス取得。
	/// デフォルトのFactoryクラスであるPolygonGroupFactoryを使用してインスタンス
	/// を生成する。
	///
	///  @return	Polylibクラスのインスタンス。
	///  @attention	呼び出し側でdeleteはできません。
	///
	static Polylib* get_instance();

	///
	/// PolygonGroupクラスを生成するためのFactoryクラスを登録。
	/// 本メソッドは、独自のFactoryクラスを登録しない限り、呼び出し不要である。
	/// コンストラクタで生成したFactoryクラスを破棄し、代わりに引数で指定された
	/// Factoryクラスを登録する。
	///
	///  @param[in] factory	Factoryクラス。
	///  @attention	PolygonGroupを拡張した場合、拡張後のPolygonGroupのFactory
	///				クラスを登録する。
	///
	void set_factory(
		PolygonGroupFactory		*factory = NULL
	);

	///
	/// PolygoGroup、三角形ポリゴン情報の読み込み。
	/// 引数で指定された設定ファイルを読み込み、グループツリーを作成する。
	/// 続いて設定ファイルで指定されたSTLファイルを読み込み、KD木を作成する。
	///
	///  @param[in] config_name 設定ファイル名。
	///  @return	POLYLIB_STATで定義される値が返る。
	///
	POLYLIB_STAT load(
		std::string			config_name = "polylib_config.xml"
	);

	///
	/// PolygoGroupツリー、三角形ポリゴン情報の保存。
	/// グループツリーの情報を設定ファイルへ出力。三角形ポリゴン情報をSTL
	///	ファイルへ出力。
	///
	///	 @param[out] p_config_name	保存した設定ファイル名の返却用。
	///  @param[in]	 stl_format		TriMeshIOクラスで定義されているSTLファイルの
	///								フォーマット。
	///  @param[in]	 extend			ファイル名に付加する文字列。省略可。省略した
	///								場合は、付加文字列として本メソッド呼び出し時
	///								の年月日時分秒(YYYYMMDD24hhmmss)を用いる。
	///  @return	POLYLIB_STATで定義される値が返る。
	///  @attention	ファイル名命名規約は次の通り。
	///			定義ファイル : polylib_config_ランク番号_付加文字.xml。
	///			STLファイル  : ポリゴングループ名_ランク番号_付加文字.拡張子。
	///
	POLYLIB_STAT save(
		std::string			*p_config_name,
		std::string			stl_format,
		std::string			extend = ""
	);

	///
	/// 三角形ポリゴン座標の移動。
	/// 本クラスインスタンス配下の全PolygonGroupのmoveメソッドが呼び出される。
	/// moveメソッドは、PolygonGroupクラスを拡張したクラスに利用者が記述する。
	///
	///  @param[in] params	Polylib.hで宣言された移動計算パラメータセット。
	///  @return	POLYLIB_STATで定義される値が返る。
	///
	POLYLIB_STAT move(
		PolylibMoveParams	&params
	);

	///
	/// PolygoGroupツリーの最上位ノードの取得。
	///
	///  @return	最上位ノードのvector。
	///  @attention 返却したPolygonGroupは、削除不可。vectorは要削除。
	///
	std::vector<PolygonGroup *> *get_root_groups() const;

	///
	/// 三角形ポリゴンの検索。
	/// 位置ベクトルmin_posとmax_posにより特定される矩形領域に含まれる、
	/// 三角形ポリゴンをgroup_nameで指定されたグループの下から探索する。
	///
	///  @param[in] group_name	抽出グループ名。
	///  @param[in] min_pos		抽出する矩形領域の最小値。
	///  @param[in] max_pos		抽出する矩形領域の最大値。
	///  @param[in] every		true:3頂点が全て検索領域に含まれるものを抽出。
	///   						false:3頂点の一部でも検索領域と重なるものを抽出。
	///  @return	抽出した三角形ポリゴンのvector。
	///  @attention 返却した三角形ポリゴンは、削除不可。vectorは要削除。
	///
	std::vector<Triangle*>* search_polygons(
		std::string		group_name, 
		Vec3f			min_pos, 
		Vec3f			max_pos, 
		bool			every
	) const;

	///
	/// 引数のグループ名が既存グループと重複しないかチェック。
	///
	///  @param[in] pg_name		グループ名
	///  @param[in] parent_path	親グループまでのフルパス
	///  @return	POLYLIB_STATで定義される値が返る。
	///  @attention	Polylib内部で使用する関数であり、通常は利用者が用いるもの
	///				ではない。
	///
	POLYLIB_STAT check_group_name(
		const std::string	&pg_name, 
		const std::string	&parent_path
	);

	///
	/// PolygonGroupのインスタンスの生成。
	/// 本クラスが管理しているFactoryクラスを利用して、引数で渡されたクラス名
	/// に応じたPolygonGroupのインスタンスを生成する。
	///
	///  @param[in] class_name		クラス名
	///  @return	生成したPolygonGroup
	///  @attention	Polylib内部で使用する関数であり、通常は利用者が用いるもの
	///				ではない。
	///
	PolygonGroup *create_polygon_group(
		std::string		class_name
	);

	///
	/// PolygonGroupの追加。
	/// 本クラスが管理しているPolygonGroupのリストにPolygonGroupを追加する。
	///
	///  @param[in] pg		PolygonGroup
	///  @attention	Polylib内部で使用する関数であり、通常は利用者が用いるもの
	///				ではない。
	///
	void add_pg_list(
		PolygonGroup	*pg
	);

	///
	/// グループ階層構造を標準出力に出力。
	/// 2010.10.20 引数FILE *追加。
	///  @param[in] fp	出力先ファイル。指定されて行ければ、標準出力へ出力する。
	///
	///  @attention	テスト用の関数であり、通常は利用者が用いるものではない。
	///
	void show_group_hierarchy(
		FILE	*fp	= NULL
	);

	///
	/// グループの情報と配下の三角形ポリゴン情報を標準出力に出力。
	/// 	親グループ名、自身の名前、STLファイル名、登録三角形数、3頂点ベクト
	///		ルの座標、法線ベクトルの座標、面積。
	///
	///  @param[in] group_name グループ名。
	///  @return	POLYLIB_STATで定義される値が返る。
	///  @attention	テスト用の関数であり、通常は利用者が用いるものではない。
	///
	POLYLIB_STAT show_group_info(
		std::string		group_name
	);

	///
	/// Polylibが利用中の概算メモリ量を返す
	///
	/// @return 利用中のメモリ量(byte)
	///
	unsigned int used_memory_size();

	///
	/// グループの取得。
	/// nameで与えられた名前のPolygonGroupを返す。
	///  @param[in] name グループ名
	///  @return	ポリゴングループクラスのポインタ。エラー時はNULLが返る。
	///  @attention オーバーロードメソッドあり。
	///
	PolygonGroup* get_group(
		std::string		name
	) const;

protected:
	///
	/// コンストラクタ
	///
	/// @attention
	///   singletonのため、子クラス以外からの呼び出し不可とする
	///
	Polylib();//PolygonGroupFactory* factory = 0);

	///
	/// デストラクタ
	///
	~Polylib();

	///
	/// グループツリー作成。
	/// 設定ファイルを管理するPolylibConfigクラスからXMLタグを得て、適切な
	/// PolygonGroupを作成し、グループツリーに登録する。
	///
	///  @param[in] config	設定ファイル管理クラス
	///  @return	POLYLIB_STATで定義される値が返る。
	///  @attention	オーバーロードメソッドあり。
	///
	/// 以下のコメントはDoxygenには出力したくないのだが...
	/// configに実体を渡すとPolylibConfigのデストラクタが2回(1回目は本関数を
	/// 抜けるとき、2回目は本関数を呼び出した関数(load_config、make_group_tree)
	/// から抜けるとき)呼ばれてしまい、結果的にSegmentation Faultで落ちてしま
	/// う。
	///
    POLYLIB_STAT make_group_tree(
        PolylibConfig	*config
    );

	///
	/// 引数の内容でグループ階層構造を構築。
	///
	///  @param[in] config_contents	設定ファイルの内容(XML形式)。
	///  @return	POLYLIB_STATで定義される値が返る。
	///  @attention	MPIPolylibクラスがMPI環境で利用することを想定している。
	///  @attention	オーバーロードメソッドあり。
	///
	POLYLIB_STAT make_group_tree(
		std::string		config_contents
	);

	///
	/// 設定ファイルを読み込み、内容をcontentsに設定。
	///
	///  @param[out] contents	設定ファイルの内容(XML形式)。
	///  @param[in]  fname		設定ファイル名。
	///  @return	POLYLIB_STATで定義される値が返る。
	///  @attention	MPIPolylibクラスがMPI環境で利用することを想定している。
	///
	POLYLIB_STAT load_config_file(
		std::string			*contents,
		std::string			fname = "" 
	);

	///
	/// 三角形IDファイルの存在が必須なload関数。
	/// loadと同様の動作を行う。但し読み込み時には、三角形IDファイルが必要で
	/// あり、このファイルに記述されているIDを用いてm_idを設定する。
	///
	///  @param[in] config_name	設定ファイル名。
	///  @param[in]	id_format	三角形IDファイルの入力形式。
	///  @return	POLYLIB_STATで定義される値が返る。
	///  @attention	MPIPolylibクラスがMPI環境で利用することを想定している。
	///
	POLYLIB_STAT load_with_idfile(
		std::string		config_name,
		ID_FORMAT		id_format
	);

	///
	/// STLファイルの読み込み。
	/// グループツリーの全リーフについて、設定されているSTLファイルから
	/// ポリゴン情報を読み込む。読み込んだ後、KD木の生成、法線の計算、面積の
	/// 計算を行う。
	///
	///  @param[in] with_id_file	trueならば、三角形ポリゴンIDファイルを読み
	///								込んでm_idを設定する。
	///								falseならば、STL読み込み時にm_idを自動生成。
	///  @param[in]	id_format		三角形IDファイルの入力形式。
	///  @return	POLYLIB_STATで定義される値が返る。
	///
	POLYLIB_STAT load_polygons(
		bool		with_id_file,
		ID_FORMAT	id_format
	);

	///
	/// 設定ファイルの保存。
	/// メモリに展開しているグループツリー情報から設定ファイルを生成する。
	///
	///  @param[in] rank_no	ランク番号
	///  @param[in] extend	ファイル名に付加する文字列
	///  @param[in] format	TriMeshIOクラスで定義されているSTLファイルのフォー
	///						マット。
	///  @return	作成した設定ファイルの名称。エラー時はNULLが返る。
	///
	char *save_config_file(
		std::string	rank_no,
		std::string	extend,
		std::string	format
	);

	/// PolygoGroupツリー、三角形ポリゴン情報の保存。
	/// グループツリー情報を設定ファイルへ出力。三角形ポリゴン情報をSTLファイル
	/// へ出力。ID情報をIDファイルへ出力。ファイル名にランク番号を付加する。
	///
	///	 @param[out] p_config_name	保存した設定ファイル名の返却用。
	///	 @param[in]  myrank			自ランク番号。
	///	 @param[in]	 maxrank		最大ランク番号。
	///	 @param[in]	 extend			ファイ名に付加される文字列。
	///	 @param[in]	 stl_format		STLファイルフォーマット指定。
	///  @param[in]	 id_format		三角形IDファイルの出力形式。
	///  @return	POLYLIB_STATで定義される値が返る。
	///  @attention	ファイル名命名規約は次の通り。
	///			定義ファイル : polylib_config_ランク番号_付加文字.xml。
	///			STLファイル  : ポリゴングループ名_ランク番号_付加文字.拡張子。
	///			IDファイル   : ポリゴングループ名_ランク番号_付加文字.ID。
	///  @attention	MPIPolylibクラスがMPI環境で利用することを想定している。
	POLYLIB_STAT save_with_rankno(
		std::string		*p_config_name,
		int				myrank,
		int				maxrank,
		std::string		extend,
		std::string		stl_format,
		ID_FORMAT		id_format
	);

	///
	/// グループ名の表示。
	/// 指定されたグループ以下の階層構造をツリー形式で標準出力に出力する。
	/// 2010.10.20 引数FILE *追加。
	///
	///  @param[in] p	検索の基点となるPolygonGroupのポインタ
	///  @param[in] tab 階層の深さを示すスペース
	///  @param[in] fp	出力先ファイル。指定されて行ければ、標準出力へ出力する。
	///
	void show_group_name(
		PolygonGroup	*p, 
		std::string		tab,
		FILE			*fp
	);

	///
	/// グループの取得。
	/// internal_idで与えられたm_internal_idを持つPolygonGroupを返す。
	///  @param[in] internal_id ポリゴングループID
	///  @return	ポリゴングループクラスのポインタ。エラー時はNULLが返る。
	///  @attention オーバーロードメソッドあり。
	///
	PolygonGroup* get_group(
		int	internal_id
	) const;

private:
	///
	/// 三角形ポリゴンの検索。
	/// 位置ベクトルmin_posとmax_posにより特定される矩形領域に含まれる、
	/// 三角形ポリゴンをgroup_nameで指定されたグループの下から探索する。
	/// 
	///  @param[in]  group_name	抽出グループ名
	///  @param[in]  min_pos	抽出する矩形領域の最小値
	///  @param[in]  max_pos	抽出する矩形領域の最大値
	///  @param[in]  every		true:3頂点が全て検索領域に含まれるものを抽出。
	///   						false:3頂点のBBoxが一部でも検索領域と重なるものを抽出。
	///  @param[in]	 linear		true:線形探索を行う。
	///							false:KD木探索を行う。
	///  @param[out] ret		POLYLIB_STATで定義される値(デバック用)。
	///  @return	抽出した三角形ポリゴンリスト
	///  @attention	publicなsearch_polygons()は内部で本関数を利用している。
	///
	std::vector<PrivateTriangle*> *search_polygons(
		std::string		group_name, 
		Vec3f			min_pos, 
		Vec3f			max_pos,
		bool			every,
		bool			linear, 
		POLYLIB_STAT	*ret
	) const;

	///
	/// グループの検索。
	/// 基点となるポリゴングループに連なる子孫ポリゴングループを全て抽出する。
	///  @param[in]  p	探索の基点となるポリゴングループへのポインタ
	///  @param[out] pg	抽出した子孫ポリゴングループのリスト
	///
	void search_group(
		PolygonGroup				*p, 
		std::vector<PolygonGroup*>	*pg
	) const;

	///
	/// XMLタグの作成。
	/// グループツリーから設定ファイルに出力するXMLタグを作成する。実際のタグ
	/// 作成はPolylibConfigクラスが受け持つ。
	/// 
	///  @param[in,out]	elem		XMLノード
	///  @param[in]		pg			出力するグループ
	///  @param[in]		rank_no		ランク番号
	///  @param[in]		extend		ファイル名に付加する文字列
	///  @param[in]		format		TriMeshIOクラスで定義されているSTLファイルの
	///								フォーマット。
	///  @return	POLYLIB_STATで定義される値が返る。
	///
	POLYLIB_STAT mk_xml(
		xmlNodePtr		elem,
		PolygonGroup	*pg,
		std::string		rank_no,
		std::string		extend,
		std::string		format
	);

protected:
	//=======================================================================
	// クラス変数
	//=======================================================================
	/// 自クラスのインスタンス(singleton)
	static Polylib				*m_instance;

	/// PolygonGroupのファクトリークラス
	PolygonGroupFactory			*m_factory;

	/// ポリゴングループリスト
	std::vector<PolygonGroup*>	m_pg_list;
};

} //namespace PolylibNS

#endif // polylib_h
