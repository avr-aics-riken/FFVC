/*
 * Polylib - Polygon Management Library.
 * Version     : 2.0.3
 * Release date: Nov.29.2010
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2010-
 *
 */
#ifndef MPIPolylib_h
#define MPIPolylib_h

#include <vector>
#include <map>
#include "mpi.h"
#include "Polylib.h"
#include "groups/PolygonGroup.h"

namespace PolylibNS {
////////////////////////////////////////////////////////////////////////////
///
/// クラス:ParallelInfo
/// 並列プロセス情報。
///
////////////////////////////////////////////////////////////////////////////
struct ParallelInfo {
	/// MPIコミュニケータ
	MPI_Comm m_comm;

	/// ランク数
	int m_rank;

	/// 計算領域情報
	CalcAreaInfo m_area;

	/// migrate除外三角形IDマップ(k:グループID, v:三角形IDリスト)
	std::map< int, std::vector<int> > m_exclusion_map;
};

////////////////////////////////////////////////////////////////////////////
///
/// クラス:MPIPolylib
/// ポリゴンを管理する為の並列版クラスライブラリです。
///
////////////////////////////////////////////////////////////////////////////
class MPIPolylib : public Polylib {
public:
	///
	/// インスタンス取得。本クラスはsingltonクラスです。
	///
	/// @return MPIPolylibクラスのインスタンス
	///
	static MPIPolylib* get_instance();

	///
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
	init_parallel_info(
		MPI_Comm comm,
		float bpos[3],
		unsigned int bbsize[3],
		unsigned int gcsize[3],
		float dx[3]
	);

	///
	/// Polylib::load()のオーバライドメソッド。
	/// @attention 並列環境では利用できません。
	///
	/// @param[in] config_filename	初期化ファイル名。
	/// @return 常に PLSTAT_NG が返ります。
	///
	POLYLIB_STAT
	load(
		std::string config_filename
	){ return PLSTAT_NG; };

	///
	/// rank0によるデータ構築。
	/// 指定された設定ファイルをrank0にて読み込み、グループ階層構造の構築
	/// およびポリゴンデータの構築を行う。
	/// グループ階層構造は全rankにb_castされ、情報を共有する。
	/// ポリゴンデータは各rank領域毎のデータが分配される。
	///
	/// @param[in] config_filename	初期化ファイル名。未指定時はデフォルトファイルを読む。
	/// @return	POLYLIB_STATで定義される値が返る。
	///
	POLYLIB_STAT
	load_rank0(
		std::string config_filename = "",
		float scale = 1.0
	);

	///
	/// 全rank並列でのデータ構築。
	/// 指定された設定ファイルを各rankにて読み込み、グループ階層構造の構築、
	/// およびポリゴンデータの構築を行う。
	/// @attention 各rankが読み込むファイルに記述されたグループ階層構造が一致している必要がある。
	///
	/// @param[in] config_filename	初期化ファイル名。未指定時はデフォルトファイルを読む。
	/// @param[in] id_format		三角形IDファイルの入力形式。
	/// @return	POLYLIB_STATで定義される値が返る。
	///
	POLYLIB_STAT
	load_parallel( 
		std::string config_filename = "",
		ID_FORMAT	id_format = ID_BIN
	);

	///
	/// Polylib::save()のオーバライドメソッド。
	/// @attention 並列環境では利用できません。
	///
	/// @param[out] p_config_filename	初期化ファイル名。
	/// @return 常に PLSTAT_NG が返ります。
	///
	POLYLIB_STAT
	save(
		std::string *p_config_filename
	){ return PLSTAT_NG; };

	///
	/// rank0によるデータ保存。
	/// rank0の本クラスインスタンスが保持するグループ階層構造を設定ファイルに書き出す。
	/// 同時に各rankに分散するポリゴンデータもrank0に集められ、指定されたフォーマットの
	/// STLファイルにrank0で書き出す。
	/// 設定ファイル命名規則は以下の通り
	///   polylib_config_付加文字列.xml
	///   polylib_config_付加文字列.tpp
	/// STLファイル命名規則は以下の通り
	///   ポリゴングループ名称_付加文字列.拡張子
	///
	/// @param[out] p_config_filename	設定ファイル名返却用stringインスタンスへのポインタ
	/// @param[in]  stl_format			STLファイルフォーマット。 "stl_a":アスキー形式　"stl_b":バイナリ形式
    /// @param[in]  extend				ファイル名に付加する文字列。省略可。省略
	///									した場合は、付加文字列として本メソッド呼
	///									び出し時の年月日時分秒(YYYYMMDD24hhmmss)
	///									を用いる。
	/// @return	POLYLIB_STATで定義される値が返る。
	/// @attention 出力引数p_config_filenameの返却値はrank0でのみ有効
	///
	POLYLIB_STAT
	save_rank0(
		std::string *p_config_filename,
		std::string stl_format,
		std::string extend = ""
	);

	///
	/// 全rank並列でのデータ保存。
	/// 各rankの本クラスインスタンスが保持するグループ階層構造を設定ファイルに各rank毎に書き出す。
	/// 同時にポリゴンデータも指定されたフォーマットのSTLファイルに各rank毎に書き出す。
	/// 設定ファイル命名規則は以下の通り
	///   polylib_config_ランク番号_付加文字列.xml
	///   polylib_config_ランク番号_付加文字列.tpp
	/// STLファイル命名規則は以下の通り
	///   ポリゴングループ名称_ランク番号_付加文字列.拡張子
	///
	/// @param[out] p_config_filename	設定ファイル名返却用stringインスタンスへのポインタ
	/// @param[in] stl_format	STLファイルフォーマット。 "stl_a":アスキー形式　"stl_b":バイナリ形式
    /// @param[in]  extend				ファイル名に付加する文字列。省略可。省略
	///									した場合は、付加文字列として本メソッド呼
	///									び出し時の年月日時分秒(YYYYMMDD24hhmmss)
	///									を用いる。
	/// @param[in] id_format	三角形IDファイルの出力形式。
	/// @return	POLYLIB_STATで定義される値が返る。
	///
	POLYLIB_STAT
	save_parallel(
		std::string *p_config_filename,
		std::string stl_format,
		std::string extend = "",
		ID_FORMAT	id_format = ID_BIN
	);

	///
	/// ポリゴン座標の移動。
	/// 本クラスインスタンス配下の全PolygonGroupのmoveメソッドが呼び出される。
	///
	/// @param[in] params	移動計算要パラメタセット。
	/// @return	POLYLIB_STATで定義される値が返る。
	///
	POLYLIB_STAT
	move(
		PolylibMoveParams &params
	);

	///
	/// ポリゴンデータのPE間移動。
	/// 本クラスインスタンス配下の全PolygonGroupのポリゴンデータについて、
	/// moveメソッドにより移動した三角形ポリゴン情報を隣接PE間でやり取りする。
	///
	/// @return	POLYLIB_STATで定義される値が返る。
	///
	POLYLIB_STAT
	migrate();

	///
	/// m_myprocの内容をget
	/// @return 自PE領域情報
	///
	ParallelInfo get_myproc(){ return m_myproc; };

	///
	/// MPIPolylibが利用中の概算メモリ量を返す
	///
	/// @return 利用中のメモリ量(byte)
	///
	unsigned int used_memory_size();

protected:
	///
	/// コンストラクタ。
	/// singletonのため非公開。本クラスインスタンス取得にはget_instance()を利用する。
	///
	MPIPolylib();

	///
	/// デストラクタ。
	///
	~MPIPolylib();

	///
	/// 指定されたグループ以下の階層構造をツリー形式で標準出力に出力する。
	///  @param p	表示対象となるグループのポインタ。
	///  @param tab	階層の深さを示すスペース。
	///  @attention プロセス毎に動作する。
    ///   出力にランク数が加わる以外は非並列版と同じ。
	///
	void show_group_name(PolygonGroup* p, std::string tab);

	///
	/// 設定ファイル内容を他rankへbroadcastする。
	///
	/// @param[in] config_contents 初期化ファイル内容。
	/// @return	POLYLIB_STATで定義される値が返る。
	///
	POLYLIB_STAT
	broadcast_config(
		std::string config_contents
	);

	///
	/// 各PE領域内ポリゴン情報を全rankに送信
	///
	/// @return	POLYLIB_STATで定義される値が返る。
	///
	POLYLIB_STAT
	send_polygons_to_all();

	///
	/// グループID＆グループ内三角形数の送信情報を作成。
	/// 
	/// @param[in,out] p_vec 情報追加先ベクタ
	/// @param[in] group_id グループID
	/// @param[in] p_trias グループ内三角形リスト
	/// @return	POLYLIB_STATで定義される値が返る。
	///
	POLYLIB_STAT
	pack_num_trias(
		std::vector<int>* p_vec,
		int group_id,
		const std::vector<PrivateTriangle*>* p_trias
	);

	///
	/// 三角形の送信情報を作成。
	/// 
	/// @param[in,out] p_vec 情報追加先ベクタ
	/// @param[in] p_trias グループ内三角形リスト
	/// @return	POLYLIB_STATで定義される値が返る。
	///
	POLYLIB_STAT
	pack_trias(
		std::vector<float>* p_vec,
		const std::vector<PrivateTriangle*>* p_trias
	);

	///
	/// 三角形IDの送信情報を作成。
	/// 
	/// @param[in,out] p_vec 情報追加先ベクタ
	/// @param[in] p_trias グループ内三角形リスト
	/// @return	POLYLIB_STATで定義される値が返る。
	///
	POLYLIB_STAT
	pack_tria_ids(
		std::vector<int>* p_vec,
		const std::vector<PrivateTriangle*>* p_trias
	);

	///
	/// 自領域内ポリゴンのみ抽出してポリゴン情報を再構築。
	/// migrate実行後に行う。
	/// 
	/// @return	POLYLIB_STATで定義される値が返る。
	///
	POLYLIB_STAT
	erase_outbounded_polygons();

	///
	/// ポリゴングループ定義情報をrank0から受信し、グループ階層構造を構築。
	///
	/// @return	POLYLIB_STATで定義される値が返る。
	///
	POLYLIB_STAT
	broadcast_config_from_rank0();

	///
	/// 自領域に必要なポリゴン情報をrank0から受信
	/// 
	/// @return	POLYLIB_STATで定義される値が返る。
	///
	POLYLIB_STAT
	receive_polygons_from_rank0();

	///
	/// 他rankからポリゴン情報をrank0で受信
	///
	POLYLIB_STAT
	gather_polygons();

	///
	/// rank0へポリゴン情報を送信
	///
	POLYLIB_STAT
	send_polygons_to_rank0();

	///
	/// 移動除外三角形IDリストの作成
	///
	POLYLIB_STAT
	select_excluded_trias( PolygonGroup *p_pg );

protected:
	///
	/// プロセス担当領域クラスのポインタを返す
	///  @param[in] rank ランク数
	///  @return プロセス担当領域クラスのポインタ
	///
	ParallelInfo* get_proc(int rank);

	/// 自PE担当領域情報
	ParallelInfo m_myproc;

	/// 自PEを除く全PE担当領域情報リスト
	std::vector<ParallelInfo*> m_other_procs;

	/// 隣接PE担当領域情報リスト
	std::vector<ParallelInfo*> m_neibour_procs;

	/// 自プロセスのランク数
	int m_myrank;

	/// 全プロセス数
	int m_numproc;

	/// 自プロセスが利用するコミュニケーター
	MPI_Comm m_mycomm;
};

} //namespace PolylibNS

#endif // MPIPolylib_h
