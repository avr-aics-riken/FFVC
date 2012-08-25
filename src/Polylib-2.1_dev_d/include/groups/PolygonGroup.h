/*
 * Polylib - Polygon Management Library.
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2010-
 *
 */
#ifndef polylib_polygongroup_h
#define polylib_polygongroup_h

#include <map>
//#include <libxml/tree.h>
#include "polygons/Triangle.h"
#include "common/PolylibStat.h"
#include "common/PolylibCommon.h"
#include "TextParser.h"

namespace PolylibNS {

class Polylib;
class Polygons;
//class PolylibCfgElem;
class PolylibMoveParams;

////////////////////////////////////////////////////////////////////////////
///
/// クラス:PolygonGroup
/// ポリゴングループを管理するクラスです。
///
////////////////////////////////////////////////////////////////////////////
class PolygonGroup {
public:
	///
	/// コンストラクタ
	///
	PolygonGroup();

	///
 	/// デストラクタ
	///
	virtual ~PolygonGroup();

	///
	/// 引数で与えられる三角形ポリゴンリストを複製し、KD木の生成を行う。
	///
	///  @param[in] tri_list	設定する三角形ポリゴンリスト。
	///  @param[in] clear		true:ポリゴン複製、面積計算、KD木生成を行う。
	///  						false:面積計算、KD木生成だけを行う。
	///  @return	POLYLIB_STATで定義される値が返る。
	///  @attention TriMeshクラスのinit()参照。オーバーロードメソッドあり。
	///
	POLYLIB_STAT init(
		const std::vector<PrivateTriangle*>		*tri_list, 
		bool									clear = true
	);

	///
	/// PolygonGroupツリーの作成。
	/// 設定ファイルの内容を再帰的に呼び出し、PolygonGroupツリーを作成する。
	///
	///  @param[in] polylib		Polygonクラスのインスタンス
	///  @param[in] parent		親グループ
	///  @param[in] tp 　　　　　　　　　　　　TextParser のインスタンス
	///  @return	POLYLIB_STATで定義される値が返る。
	///
	virtual POLYLIB_STAT build_group_tree(
		Polylib					*polylib,
		PolygonGroup			*parent,
		TextParser* tp
	);

	///
	/// 三角形ポリゴンの法線ベクトルの計算、面積の計算、KD木の生成を行う。
	/// 三角形ポリゴンはTriMeshクラスが管理している。
	///
	///  @return	POLYLIB_STATで定義される値が返る。
	///  @attention	TriMeshクラスのbuild()参照。
	///
	POLYLIB_STAT build_polygon_tree();

	///
	/// STLファイルからポリゴン情報を読み込み、TriMeshクラスに登録する。
	///
	///  @return POLYLIB_STATで定義される値が返る。
	///  @attention TriMeshクラスのimport()参照。
	///
	POLYLIB_STAT load_stl_file();

	///
	/// 三角形ポリゴンIDファイルからポリゴンIDを読み込み、m_internal_idに登録する。
	///
	///  @param[in] id_format	三角形IDファイルの入力形式。
	///  @return	POLYLIB_STATで定義される値が返る。
	///
	POLYLIB_STAT load_id_file(
		ID_FORMAT		id_format
	);

	///
	/// TriMeshクラスが管理しているポリゴン情報をSTLファイルに出力する。
	///
	///  @param[in] rank_no	ファイル名に付加するランク番号。
	///  @param[in] extend	ファイル名に付加する自由文字列。
	///  @param[in] format	STLファイルフォーマット。
	///  @return	POLYLIB_STATで定義される値が返る。
	///  @attention TriMeshIOクラスのsave()参照。オーバーロードメソッドあり。
	///
	POLYLIB_STAT save_stl_file(
		std::string		rank_no,
		std::string		extend,
		std::string		format
	);

	///
	/// TriMeshクラスが管理しているポリゴン情報をSTLファイルに出力する。
	/// TextParser 対応版
	///  @param[in] rank_no	ファイル名に付加するランク番号。
	///  @param[in] extend	ファイル名に付加する自由文字列。
	///  @param[in] format	STLファイルフォーマット。
	///  @param[in,out] stl_fname_map stl ファイル名とポリゴングループのパス
	///  @return	POLYLIB_STATで定義される値が返る。
	///  @attention TriMeshIOクラスのsave()参照。オーバーロードメソッドあり。
	///
	POLYLIB_STAT save_stl_file(
		std::string		rank_no,
		std::string		extend,
		std::string		format,
		std::map<std::string,std::string>& stl_fname_map
	);

	///
	/// 三角形ポリゴンIDファイルにポリゴンIDを出力する。IDファイル名は、
	/// 階層化されたグループ名_ランク番号_自由文字列.id。
	///
	///  @param[in] rank_no		ファイル名に付加するランク番号。
	///  @param[in] extend		ファイル名に付加する自由文字列。
	///  @param[in] id_format	三角形IDファイルの出力形式。
	///  @return	POLYLIB_STATで定義される値が返る。
	/// 
	POLYLIB_STAT save_id_file(
		std::string 	rank_no,
		std::string		extend,
		ID_FORMAT		id_format
	);


	///
	/// 設定ファイルに出力するTextParserのリーフを編集する.
	/// デフォルトでは何もしない。
	/// CarGroup.cxx の例を参照.
	///
	///  @param[in] pointer to TextParser 
	///  @param[in] rank_no	ファイル名に付加するランク番号。
	///  @param[in] extend	ファイル名に付加する自由文字列。
	///  @param[in] format	STLファイルフォーマット。
	///  @return	POLYLIB_STATで定義される値が返る。
	///  @attention  do nothing by default
	///

	virtual POLYLIB_STAT mk_param_tag(
					  TextParser* pt,
		std::string		rank_no,
		std::string		extend,
		std::string		format
	);


	///
	/// 三角形ポリゴン移動メソッド。virtual用の関数なので処理はない。
	///
	///  @param[in] params	Polylib.hで宣言しているパラメタセットクラス。
	///  @return	POLYLIB_STATで定義される値が返る。
	///
	virtual POLYLIB_STAT move(
		PolylibMoveParams	&params
	);

	///
	/// KD木探索により、指定矩形領域に含まれるポリゴンを抽出する。
	///
	///  @param[in] bbox	矩形領域。
	///  @param[in]	every	true:3頂点が全て検索領域に含まれるものを抽出。
	///  					false:1頂点でも検索領域に含まれるものを抽出。
	///  @return	抽出したポリゴンリストのポインタ。
	///  @attention オーバーロードメソッドあり。
	///
	const std::vector<PrivateTriangle*>* search(
		BBox	*bbox, 
		bool	every
	) const;

	///
	/// KD木探索により、指定矩形領域に含まれるポリゴンを抽出する。
	///
	///  @param[in]		bbox		矩形領域。
	///  @param[in]		every		true:3頂点が全て検索領域に含まれるものを抽出。
	///  							false:1頂点でも検索領域に含まれるものを抽出。
	///  @param[in,out] tri_list	抽出した三角形ポリゴンリストのポインタ。
	///  @return	POLYLIB_STATで定義される値が返る。
	///  @attention オーバーロードメソッドあり。
	///
	POLYLIB_STAT search(
		BBox							*bbox, 
		bool							every, 
		std::vector<PrivateTriangle*>	*tri_list
	) const;

	///
	/// 線形探索により、指定矩形領域に含まれるポリゴンを抽出する。
	///
	///  @param[in] bbox	矩形領域。
	///  @param[in]	every	true:3頂点が全て検索領域に含まれるものを抽出。
	///  					false:1頂点でも検索領域に含まれるものを抽出。
	///  @return	抽出したポリゴンリストのポインタ。
	///  @attention	オーバーロードメソッドあり。
	///
	const std::vector<PrivateTriangle*>* linear_search(
		BBox	*bbox, 
		bool	every
	) const;

	///
	/// 線形探索により、指定矩形領域に含まれるポリゴンを抽出する。
	///
	///  @param[in]		bbox		 矩形領域。
	///  @param[in]		every	 	true:3頂点が全て検索領域に含まれるものを抽出。
	///  						 	false:1頂点でも検索領域に含まれるものを抽出。
	///  @param[in,out] tri_list	抽出した三角形ポリゴンリストのポインタ。
	///  @return	POLYLIB_STATで定義される値が返る。
	///  @attention	オーバーロードメソッドあり。
	///
	POLYLIB_STAT linear_search(
		BBox							*bbox, 
		bool							every, 
		std::vector<PrivateTriangle*>	*tri_list
	) const;

	///
	/// PolygonGroupのフルパス名を取得する。
	///
	///  @return フルパス名。
	///
	std::string acq_fullpath();

	///
	/// カンマ区切りでSTLファイル名リストを取得。
	///
	///  @return ファイル名リスト。
	///
	std::string acq_file_name();

	///
	/// PE領域間移動する三角形ポリゴンリストの取得。
	///
	///  @param[in]	neibour_bbox		隣接PE領域バウンディングボックス。
	///  @param[in]	exclude_tria_ids	領域移動対象外三角形IDリスト。
	///  @return	検索結果三角形リスト。
	///
	const std::vector<PrivateTriangle*>* search_outbounded(
		BBox				neibour_bbox,
		std::vector<int>	*exclude_tria_ids
	);

	///
	/// 三角形リストの追加。
	///
	///  @param[in]	tri_list	三角形ポリゴンリストのポインタ。
	///  @return	POLYLIB_STATで定義される値が返る。
	///  @attention	三角形IDが重複した三角形は追加しない。KD木の再構築はしない。
	///
	POLYLIB_STAT add_triangles(
		std::vector<PrivateTriangle*>	*tri_list
	);

	///
	/// ポリゴン情報を再構築する。（KD木の再構築をおこなう）
	///
	///  @return	POLYLIB_STATで定義される値が返る。
	///
	POLYLIB_STAT rebuild_polygons();

	///
	/// グループ情報（ランク番号、親グループ名、自分のグループ名、ファイル名、
	/// 頂点数、各頂点のXYZ座標値、法線ベクトルのXYZ座標値、面積）を出力する。
	///
	///  @param[in] irank ランク数。
	///  @return	POLYLIB_STATで定義される値が返る。
	///
	POLYLIB_STAT show_group_info(
		int irank = -1
	);

  // add keno 20120331
  /// ポリゴングループの要素数を返す
  int get_group_num_tria( void );
  
  /// ポリゴンの面積を積算して返す
  float get_group_area( void );
  
	//=======================================================================
	// Setter/Getter
	//=======================================================================
	///
	/// クラス名を取得。
	///
	/// @return		クラス名。
	/// @attention	本クラスを継承する場合、継承後のクラス名を返すように変更す
	/// 			ることる。
	///
	static std::string get_class_name() {return "PolygonGroup";};

	///
	/// クラス名を取得。
	///
	/// @return		クラス名。
	/// @attention	継承するクラスのクラス名取得関数get_class_name()を呼び出す。
	///
	virtual std::string whoami()		{return get_class_name();};

	///
	/// STLファイル名とファイルフォーマットを設定。
	///
	///  @param[in] fname STLファイル名とファイルフォーマットの対応マップ。
	///
	void set_file_name(std::map<std::string, std::string> fname) {
		m_file_name = fname;
	}

	///
	/// STLファイル名とファイルフォーマットの対応マップ取得。
	///
	///  @return STLファイル名とファイルフォーマットの対応マップ。
	///
	std::map<std::string, std::string> get_file_name() const {
		return m_file_name;
	}

	///
	/// グループ名を設定。
	///
	/// @param[in] name グループ名。
	///
	void set_name(std::string name){
		m_name = name;
	}

	///
	/// グループ名を取得。
	///
	/// @return グループ名。
	///
	std::string get_name(void){
		return m_name;
	}

	///
	/// 親グループのフルパス名を設定。
	///
	/// @param[in] ppath 親グループのフルパス名。
	///
	void set_parent_path(std::string ppath){
		m_parent_path = ppath;
	}

	///
	/// 親グループのフルパス名を取得。
	///
	/// @return 親グループのフルパス名。
	///
	std::string get_parent_path(void){
		return m_parent_path;
	}

	///
	/// 親グループを設定。
	///
	/// @param[in] p 親グループのポインタ。
	///
	void set_parent(PolygonGroup* p) {
		m_parent = p;
	}

	///
	/// 親グループを取得。
	///
	/// @return 親グループのポインタ。
	///
	PolygonGroup* get_parent(void) {
		return m_parent;
	}

	///
	/// 子グループを設定。
	///
	/// @param[in] p	子グループのリスト。
	///
	void set_children(std::vector<PolygonGroup*>& p) {
		m_children = p;
	}

	///
	/// 子グループを取得。
	///
	/// @return 子グループのリスト。
	///
	std::vector<PolygonGroup*>& get_children(void) {
		return m_children;
	}

	///
	/// 子グループを追加。
	///
	/// @param[in] p	子グループ。
	///
	void add_children(PolygonGroup* p) {
		m_children.push_back(p);
	}

	///
	/// Polygonクラスが管理する三角形ポリゴンリストを取得。
	///
	/// @return 三角形ポリゴンリスト。
	///
	std::vector<PrivateTriangle*>* get_triangles() {
		return m_polygons->get_tri_list();
	}

	///
	/// Polygonクラスが管理するKD木クラスを取得。
	///
	/// @return KD木ポリゴンリスト。
	///
	VTree *get_vtree() {
		return m_polygons->get_vtree();
	}

	///
	/// ポリゴングループIDを取得。
	/// メンバー名修正( m_id -> m_internal_id) 2010.10.20
	///
	///  @return ポリゴングループID。
	///
	int get_internal_id() {
		return m_internal_id;
	}

	///
	/// ユーザ定義IDを取得。
	/// 処理追加 2010.10.20
	///
	///  @return ユーザ定義ID。
	///
	int get_id() {
		return m_id;
	}

	///
	/// 移動対象フラグを取得。
	///
	///  @return 移動対象フラグ。
	///
	int get_movable() {
		return m_movable;
	}

	///
	/// move()による移動前三角形一時保存リストの個数を取得。
	///
	///  @return 一時保存リストサイズ。
	///
	size_t get_num_of_trias_before_move() {
		if (m_trias_before_move == NULL)	return 0;
		else								return m_trias_before_move->size();
	}

	///
	/// configファイルに記述するParamタグのクラス名(value="...")。
	///
	static const char *ATT_NAME_CLASS;

protected:

	///
	/// 設定ファイルから取得したPolygonGroupの情報をインスタンスにセットする。
	///
	///  @param[in] polylib		Polygonクラスのインスタンス。
	///  @param[in] parent		親グループ。
	///  @param[in] tp              TextParserクラスのインスタンス
	///  @return	POLYLIB_STATで定義される値が返る。
	///
	POLYLIB_STAT setup_attribute (
		Polylib					*polylib,
		PolygonGroup			*parent,
		TextParser *tp
	);

	///
	/// move()メソッド実行により、頂点が隣接セルよりも遠くへ移動した三角形情報
	/// を報告（前処理）。
	///
	///  @return	POLYLIB_STATで定義される値が返る。
	///  @attention 本メソッドはデバッグ用です。
	///  			派生クラスでオーバーライドしたmove()メソッド内で、座標移動
	///				処理前に呼ぶこと。
	///
	POLYLIB_STAT init_check_leaped();

	///
	/// move()メソッド実行により、頂点が隣接セルよりも遠くへ移動した三角形情報
	/// を報告（後処理）。該当する三角形について、以下の情報をcerrへ出力する。
	///		・ポリゴングループID
	///		・三角形ID
	///		・移動前/後の頂点座標
	///
	///  @param[in] origin		計算領域起点座標
	///  @param[in] cell_size	ボクセルサイズ
	///  @return	POLYLIB_STATで定義される値が返る。
	///  @attention	本メソッドはデバッグ用です。
	///  			派生クラスでオーバーライドしたmove()メソッド内で、座標移動
	///				処理後に呼ぶこと。
	///
	POLYLIB_STAT check_leaped(
		Vec3f origin,
		Vec3f cell_size
	);

	///
	/// 2点が隣接ボクセルよりも離れているか？
	///
	///  @param[in] origin		計算領域起点座標。
	///  @param[in] cell_size	ボクセルサイズ。
	///  @param[in] pos1			点(1)。
	///  @param[in] pos2			点(2)。
	///  @return	true:2点が隣接ボクセルよりも離れている。
	///
	bool is_far(
		Vec3f origin,
		Vec3f cell_size,
		Vec3f pos1,
		Vec3f pos2
	);

private:
	///
	/// STLファイル名を作成。ファイル名は、以下の通り。
	/// グループ名のフルパス_ランク番号_自由文字列.フォーマット文字列。
	///
	///  @param[in] rank_no	ファイル名に付加するランク番号。
	///  @param[in] extend	ファイル名に付加する自由文字列。
	///  @param[in] format	STLファイルフォーマット。
	///  @return	POLYLIB_STATで定義される値が返る。
	///
	char *mk_stl_fname(
		std::string		rank_no,
		std::string		extend,
		std::string		format
	);


	///
	/// STLファイル名を作成。ファイル名は、以下の通り。
	/// グループ名のフルパス_ランク番号_自由文字列.フォーマット文字列。
	/// TextParser 対応版
	///
	///  @param[in] rank_no	ファイル名に付加するランク番号。
	///  @param[in] extend	ファイル名に付加する自由文字列。
	///  @param[in] format	STLファイルフォーマット。
	///  @param[in,out] stl_fname_map stl ファイル名とポリゴングループのパス
	///  @return	POLYLIB_STATで定義される値が返る。
	///
	char *mk_stl_fname(
		std::string		rank_no,
		std::string		extend,
		std::string		format,
		std::map<std::string,std::string>& stl_fname_map
	);

	///
	/// 三角形ポリゴンIDファイル名を作成。ファイル名は、以下の通り。
	/// グループ名のフルパス_ランク番号_自由文字列.id。
	///
	///  @param[in] rank_no	ファイル名に付加するランク番号。
	///  @param[in] extend	ファイル名に付加する自由文字列。
	///  @return	POLYLIB_STATで定義される値が返る。
	///
	char *mk_id_fname(
		std::string		extend,
		std::string		rank_no
	);

	///
	/// 全PolygonGroupに一意のグループIDを作成する。
	///
	///  @return	グループID。
	///
	int create_global_id();

protected:
	//=======================================================================
	// クラス変数
	//=======================================================================
	/// グループID。
	// メンバー名修正( m_id -> m_internal_id) 2010.10.20
	int								m_internal_id;

	/// 自グループ名。
	std::string						m_name;

	/// 親グループのパス名。
	std::string							m_parent_path;

	/// 親グループへのポインタ。
	PolygonGroup						*m_parent;

	/// 子グループへのポインタリスト。
	std::vector<PolygonGroup*>			m_children;

	/// STLファイル名とファイル形式。
	std::map<std::string, std::string>	m_file_name;

	/// 三角形Polygonsクラス。
	Polygons							*m_polygons;

	/// moveメソッドにより移動するグループか？
	bool								m_movable;

	/// KD木の再構築が必要か？
	bool								m_need_rebuild;

	/// move()による移動前三角形一時保存リスト。
	std::vector<PrivateTriangle*>		*m_trias_before_move;
private:
	/// ユーザ定義id : 追加 2010.10.20
	int									m_id;
};

} //namespace PolylibNS

#endif //polylib_polygongroup_h
