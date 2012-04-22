/*
 * Polylib - Polygon Management Library.
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2010-
 *
 */
#ifndef polylib_config_h
#define polylib_config_h

#include <vector>
#include <libxml/tree.h>
#include <common/PolylibStat.h>
#include <common/PolylibCommon.h>

namespace PolylibNS {

///
/// Paramで用いるデータ型
///
enum PolylibCfgParamType
{
	INT,
	REAL,
	STRING
};

////////////////////////////////////////////////////////////////////////////
///
/// クラス:PolylibCfgParam
///  @attention	定義ファイルのParam要素の管理。 Param要素の形式は
///				"<Param name="" dtype="" value="">"である。
///
////////////////////////////////////////////////////////////////////////////
class PolylibCfgParam {
public:
	///
	/// コンストラクタ
	///
	///  @param[in] name	Param名。
	///  @param[in] type	Paramデータ型（STRING/INT/REAL のいずれか）。
	///  @param[in] value	Paramのvalue値。
	///
	PolylibCfgParam(
		const std::string	name, 
		const std::string	type, 
		const std::string	value
	);

	//=======================================================================
	// Setter/Getter
	//=======================================================================
	///
	/// Paramのデータ型。
	///  @return データ型。
	///
	const PolylibCfgParamType get_data_type() const {
		return m_type;
	}

	///
	/// 文字列データ取得。
	///  @return 文字列型データ。
	///
	const std::string get_string_data() const {
		return m_str;
	}

	///
	/// INTEGER型データ取得。
	///  @return INTEGER型データ。
	///
	const int get_int_data() const {
		return m_int;
	}

	///
	/// FLOAT型データ取得。
	///  @return float型のデータ。
	///
	const float get_real_data() const {
		return m_real;
	}

	///
	/// Param名 。
	///  @return パラメータ名。
	///
	const std::string get_name() const {
		return m_name;
	}

private:
	//=======================================================================
	// クラス数
	//=======================================================================
	/// Param名 。
	std::string			m_name;

	/// Paramのデータ型。
	PolylibCfgParamType	m_type;

	/// Paramのvalue値(INTEGER型)。
	int					m_int;

	/// Paramのvalue値(FLOAT型)。
	float				m_real;

	/// Paramのvalue値(STRING型)。
	std::string			m_str;
};

////////////////////////////////////////////////////////////////////////////
///
/// クラス:PolylibCfgElem
///  @attention	configファイルのElem要素の管理。 Elem要素の形式は
///				"<Elem name="">"である。
///
////////////////////////////////////////////////////////////////////////////
class PolylibCfgElem {
public:
	///
	/// コンストラクタ。
	///
	PolylibCfgElem(
		const std::string name
	);

	///
 	/// デストラクタ。
	///
	~PolylibCfgElem();

	///
	/// nameで指定される最初のElem要素を返す。
	///
	///  @param[in] name	要素名(指定しない場合最初の要素)。
	///  @return Elem要素（存在しない場合NULLが返る）。
	///
	const PolylibCfgElem* first_element(
		const std::string name = ""
	) const;

	///
	/// nameで指定される次のElem要素を返す。
	///
	///  @param[in] param	前の要素。
	///  @param[in] name	要素名(指定しない場合すぐ次の要素)。
	///  @return Elem要素（存在しない場合NULLが返る）。
	///
	const PolylibCfgElem* next_element(
		const PolylibCfgElem	*param, 
		const std::string		name = ""
	) const;

	///
	/// nameで指定される次のElem要素を返す。
	///
	///  @param[in] name	要素名(指定しない場合最初の要素)。
	///  @return Param要素（存在しない場合NULLが返る）。
	///
	const PolylibCfgParam* first_param(
		const std::string	name = ""
	) const;

	///
	/// nameで指定される次のParam要素を返す。
	///
	///  @param[in] param	前の要素。
	///  @param[in] name	要素名(指定しない場合すぐ次の要素)。
	///  @return Param要素（存在しない場合NULLが返る）
	///
	const PolylibCfgParam* next_param(
		const PolylibCfgParam	*param, 
		const std::string		name = ""
	) const;

	//=======================================================================
	// Setter/Getter
	//=======================================================================
	///
	/// Elem名を返す。
	///  @return Elem名。
	///
	const std::string get_name() const {
		return m_name;
	}

	///
	/// 子要素のElemを追加。
	///  @param[in] elem 追加するElem要素。
	///
	void set_elem(PolylibCfgElem* elem) {
		m_elem->push_back(elem);
	}

	///
	/// 子要素のParamを追加。
	///  @param[in] param 追加するParam要素。
	///
	void set_param(PolylibCfgParam* param) {
		m_param->push_back(param);
	}

private:
	//=======================================================================
	// クラス変数
	//=======================================================================
	/// Elem名。
	std::string						m_name;

	/// 子要素のElemリスト。
	std::vector<PolylibCfgElem*>	*m_elem;

	/// 子要素のParamリスト。
	std::vector<PolylibCfgParam*>	*m_param;
};

////////////////////////////////////////////////////////////////////////////
///
/// クラス:PolylibConfig
/// configファイルの管理
///  @attention XMLの書式はV-Sphereに準拠する。
///
////////////////////////////////////////////////////////////////////////////
class PolylibConfig {
public:
	///
	/// コンストラクタ。
	///
	///  @attention	オーバーロードメソッドあり。
	///
	PolylibConfig(std::string fname);

	///
	/// コンストラクタ。
	///
	///  @Exception 設定ファイルに不備があった場合は例外PLSTAT_NGを投げる。
	///  @attention	オーバーロードメソッドあり。
	///
	PolylibConfig();

	///
 	/// デストラクタ。
	///
	~PolylibConfig();

	///
	/// 引数で渡されたXML形式のデータをメモリ展開する。
	///
	///  @param[in] contents XML形式の文字列。
	///  @return	POLYLIB_STATで定義される値が返る。
	///
	POLYLIB_STAT parse_xml_on_memory(
		std::string		contents
	);

	///
	/// 設定ファイル読み込み。
	/// 設定ファイルを読み込み、libxml2ライブラリを用いて構文解析した後、
	/// 引数contentsに代入して上位に戻す。
	/// 
	///  @param[out] contents	XML形式の文字列。
	///  @param[in]  fname		設定ファイル名。
	///  @return	POLYLIB_STATで定義される値が返る。
	///
	static POLYLIB_STAT load_config_file(
		std::string		*contents,
		std::string		fname
	);

	///
	/// 設定ファイルに出力するParameterタグを作成する。
	///
	///  @param[in] doc	libxml2ライブラリで定義しているXML文書構造体。
	///  @return	作成したParameterタグ。エラー時にはNULLが返る。
	///
	static xmlNodePtr mk_parameter_tag(
		xmlDocPtr	doc
	);

	///
	/// 設定ファイルに出力するElemタグを作成する。
	///
	///  @param[in] elem	Parameterタグ構造体、またはElemタグ構造体。
	///  @return 作成しElemタグ。エラー時にはNULLが返る。
	///
	static xmlNodePtr mk_elem_tag(
		xmlNodePtr	elem
	);

	///
	/// 設定ファイルに出力するParamタグを作成する。出力する属性値は文字列型。
	///
	///  @param[in] elem	親のElemタグ。
	///  @param[in] name	Paramタグの属性名。
	///  @param[in] value	Paramタグの属性値。
	///  @return	作成しParamタグ。エラー時にはNULLが返る。
	///  @attention	オーバーロードメソッドあり。
	///
	static POLYLIB_STAT mk_param_tag(
		xmlNodePtr		elem,
		std::string 	name,
		std::string 	value
	);

	///
	/// 設定ファイルに出力するParamタグを作成する。出力する属性値は整数型。
	///
	///  @param[in] elem	親のElemタグ。
	///  @param[in] name	Paramタグの属性名。
	///  @param[in] value	Paramタグの属性値。
	///  @return	作成しParamタグ。エラー時にはNULLが返る。
	///  @attention	オーバーロードメソッドあり。
	///
	static POLYLIB_STAT mk_param_tag(
		xmlNodePtr		elem,
		std::string		name,
		int				value
	);

	///
	/// 設定ファイルに出力するParamタグを作成する。出力する属性値は実数型。
	///
	///  @param[in] elem	親のElemタグ。
	///  @param[in] name	Paramタグの属性名。
	///  @param[in] value	Paramタグの属性値。
	///  @return	作成しParamタグ。エラー時にはNULLが返る。
	///  @attention	オーバーロードメソッドあり
	///
	static POLYLIB_STAT mk_param_tag(
		xmlNodePtr		elem,
		std::string		name,
		double			value
	);

	///
	/// 設定情報をXML形式でファイルに出力する。
	/// 設定ファイルのファイル名は、polylib_config_ランク番号_自由文字列.xml。
	///
	///  @param[in] doc		libxml2ライブラリで定義しているXML構造体。
	///  @param[in] rank_no	ファイル名に付加するランク番号。
	///  @param[in] extend	ファイル名に付加する自由文字列。
	///  @return	出力ファイル名。エラー時にはNULLが返る。
	///  @attention	戻り値のchar *はフリー不要。
	static char *save_file(
		xmlDocPtr		doc,
		std::string		rank_no,
		std::string		extend
	);

	//=======================================================================
	// Setter/Getter
	//=======================================================================
	///
 	/// ルートノード取得。
	///
	///  @return	Elemタグ構造体。
	///
	const PolylibCfgElem* get_root_elem() const {
		return m_elem;
	}

private:
	///
	/// libxml2を使用して、XML形式の設定情報を構文解析し、PolylibCfgElemまたは
	/// PolylibCfgParamのインスタンスを生成する。
	///
	///  @param[in] parent	親のタグ情報。
	///  @param[in] node	子ノード。
	///  @return	POLYLIB_STATで定義される値が返る。
	///
	POLYLIB_STAT xml_parse(
		PolylibCfgElem	*parent, 
		xmlNodePtr		node
	);

	//=======================================================================
	// クラス変数
	//=======================================================================
	/// ルートノード。
	PolylibCfgElem	*m_elem;
};

} //namespace PolylibNS

#endif //polylib_config_h
