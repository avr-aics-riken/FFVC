/*
 * Polylib - Polygon Management Library.
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2010-
 *
 */
#include <fstream>
#include "Polylib.h"

using namespace std;
using namespace PolylibNS;

/// デバッグ用ランク番号グローバル文字列
std::string PolylibNS::gs_rankno = "";

/************************************************************************
 *
 * Polylibクラス
 *
 ***********************************************************************/
// public /////////////////////////////////////////////////////////////////////
Polylib* Polylib::get_instance() {
	// この方法ならば、デストラクタを呼び出さなくてもクラスインスタンスの領域
	// が解放できるし、もし本関数が複数回呼び出されても、クラスインスタンスが
	// 複数作成されることはない(=singletonになる)
	static Polylib m_instance;
	return &m_instance;
}

// public /////////////////////////////////////////////////////////////////////
void Polylib::set_factory(
	PolygonGroupFactory		*factory
) {
	if (factory == NULL) return;

	// 明示的に指定された場合は、そのFactoryクラスを使用する
	delete m_factory;
	m_factory = factory;
}

// public /////////////////////////////////////////////////////////////////////
POLYLIB_STAT Polylib::load(
	string	config_name
) {
#ifdef DEBUG
	PL_DBGOSH << "Polylib::load() in." << endl;
#endif
	// 設定ファイル読み込み
	try {
		PolylibConfig base(config_name);

		// グループツリー作成
		POLYLIB_STAT stat = make_group_tree(&base);
		if (stat != PLSTAT_OK)	return stat;

		// STLファイル読み込み (三角形IDファイルは不要なので、第二引数はダミー)
		return load_polygons(false, ID_BIN);
	}
	catch (POLYLIB_STAT e) {
		return e;
	}
}

// public /////////////////////////////////////////////////////////////////////
POLYLIB_STAT Polylib::save(
	string	*p_config_name,
	string	stl_format,
	string	extend
) {
#ifdef DEBUG
	PL_DBGOSH << "Polylib::save() in." << endl;
#endif
	char	my_extend[128];

	// 拡張文字列がカラであれば、現在時刻から作成
	if (extend == "") {
		time_t		timer = time(NULL);
		struct tm	*date = localtime(&timer);
		sprintf(my_extend, "%04d%02d%02d%02d%02d%02d",
			date->tm_year+1900,	date->tm_mon+1,	date->tm_mday,
			date->tm_hour,		date->tm_min,	date->tm_sec);
	}
	else {
		sprintf(my_extend, "%s", extend.c_str());
	}

	char	*config_name = save_config_file("", my_extend, stl_format);
	if (config_name == NULL)	return PLSTAT_NG;
	else						*p_config_name = string(config_name);

	vector<PolygonGroup*>::iterator	it;
	for (it = m_pg_list.begin(); it != m_pg_list.end(); it++) {
		//リーフのみがポリゴン情報を持っている
		if ((*it)->get_children().empty() == false)	continue;

		// ポリゴン数が0ならばファイル出力不要 2010.10.19
		if ((*it)->get_triangles()->size() == 0)	continue;

		// STLファイル保存 (第一引数のランク番号は不要)
		POLYLIB_STAT stat = (*it)->save_stl_file("", my_extend, stl_format);
		if (stat != PLSTAT_OK) return stat;
	}
	return PLSTAT_OK;
}

// public /////////////////////////////////////////////////////////////////////
POLYLIB_STAT Polylib::move(
	PolylibMoveParams	&params
) {
#ifdef DEBUG
	PL_DBGOSH << "Polylib::move() in." << endl;
#endif
	POLYLIB_STAT ret;
	vector<PolygonGroup*>::iterator it;

	for (it = m_pg_list.begin(); it != m_pg_list.end(); it++) {

		// リーフグループで、movableフラグONのポリゴンを移動
		if ((*it)->get_children().empty() == true && (*it)->get_movable() ) {
			ret = (*it)->move(params);
			if (ret != PLSTAT_OK)	return ret;

			// 座標移動したのでKD木の再構築
			ret = (*it)->rebuild_polygons();
			if( ret != PLSTAT_OK ) return ret;
		}

	}
	return PLSTAT_OK;
}

// public /////////////////////////////////////////////////////////////////////
vector<PolygonGroup *> *Polylib::get_root_groups() const {
#ifdef DEBUG
	PL_DBGOSH << "Polylib::get_root_groups() in." << endl;
#endif
	vector<PolygonGroup*>					*root = new vector<PolygonGroup*>;
	vector<PolygonGroup*>::const_iterator	it;

	for (it = m_pg_list.begin(); it != m_pg_list.end(); it++) {
		if (((*it)->get_parent()) == NULL) {
			root->push_back(*it);
		}
	}
	return root;
}

// public /////////////////////////////////////////////////////////////////////
vector<Triangle*>* Polylib::search_polygons(
	string		group_name, 
	Vec3f		min_pos, 
	Vec3f		max_pos, 
	bool		every
) const {
#ifdef DEBUG
	PL_DBGOSH << "Polylib::search_polygons() in." << endl;
#endif
	POLYLIB_STAT ret;
	return (vector<Triangle*>*)
		search_polygons(group_name, min_pos, max_pos,every, false, &ret);
}

// public /////////////////////////////////////////////////////////////////////
POLYLIB_STAT Polylib::check_group_name(
	const string	&name, 
	const string	&path
) {
#ifdef DEBUG
	PL_DBGOSH << "Polylib::check_group_name() in." << endl;
#endif
	if (name.empty() == true) {
		PL_ERROSH << "[ERROR]Polylib::check_group_name():Group name is empty." 
				  << endl;
		return PLSTAT_GROUP_NAME_EMPTY;
	}

	vector<PolygonGroup*>::iterator it;
	for (it = m_pg_list.begin(); it != m_pg_list.end(); it++) {
		if ((*it)->get_name() == name && (*it)->get_parent_path() == path) {
			PL_ERROSH << "[ERROR]Polylib::check_group_name():Group name is "
				<< "duplicate:name:" << name << "," << "path:" << path 
				<< endl;
			return PLSTAT_GROUP_NAME_DUP;
		}
	}
	return PLSTAT_OK;
}

// public /////////////////////////////////////////////////////////////////////
PolygonGroup *Polylib::create_polygon_group(string class_name)
{
#ifdef DEBUG
	PL_DBGOSH << "Polylib::create_polygon_group() in." << endl;
#endif
	return m_factory->create_instance(class_name);
}

// public /////////////////////////////////////////////////////////////////////
void Polylib::add_pg_list(PolygonGroup *pg)
{
#ifdef DEBUG
	PL_DBGOSH << "Polylib::add_pg_list() in." << endl;
#endif
	m_pg_list.push_back(pg);
}

// public /////////////////////////////////////////////////////////////////////
/// 2010.10.20 引数FILE *追加。
void Polylib::show_group_hierarchy(FILE	*fp)
{
#ifdef DEBUG
	PL_DBGOSH << "Polylib::show_group_hierarchy() in." << endl;
#endif
	string		tab;
	vector<PolygonGroup*>::iterator	it;
	for (it = m_pg_list.begin(); it != m_pg_list.end(); it++) {
		if ((*it)->get_parent() != NULL) {
			// Not Use
		}
		else{
			show_group_name(*it, tab, fp);
		}
	}
}

// public /////////////////////////////////////////////////////////////////////
POLYLIB_STAT Polylib::show_group_info(string group_name)
{
#ifdef DEBUG
	PL_DBGOSH << "Polylib::show_group_info() in." << endl;
#endif
	PolygonGroup *p = get_group(group_name);
	if (p == NULL) {
		PL_ERROSH << "[ERROR]Polylib::show_group_info():Group not found:"
			 << group_name << endl;
		return PLSTAT_GROUP_NOT_FOUND;
	}

	return p->show_group_info();
}

// public /////////////////////////////////////////////////////////////////////
unsigned int Polylib::used_memory_size()
{
	unsigned int						size;
	vector<PolygonGroup*>::iterator		pg;
	vector<PrivateTriangle*>::iterator	pt;

	// 自クラスとFactoryクラス
	size = sizeof(Polylib) + sizeof(PolygonGroupFactory);

	// ポリゴングループ
#ifdef DEBUG
PL_DBGOSH << "Polylib::used_memory_size:PolygonGroup num=" << m_pg_list.size() << endl;
#endif
	for (pg = m_pg_list.begin(); pg != m_pg_list.end(); pg++) {

		// PolygonGroupクラス
		size += sizeof(PolygonGroup);

		// TriMeshクラス
		size += sizeof(TriMesh);

		// 三角形移動前一時リスト
		size += (*pg)->get_num_of_trias_before_move() * sizeof(vector<void> *);

		// リーフにはポリゴンがある
		if ((*pg)->get_children().empty()) {

			// 三角形ポリゴン
			vector<PrivateTriangle*>	*tri_list = (*pg)->get_triangles();
#ifdef DEBUG
PL_DBGOSH << "Polylib::used_memory_size:PrivateTriangle num=" << tri_list->size() << endl;
#endif
			size += tri_list->size() * sizeof(PrivateTriangle);

			// KD木
			VTree	*vtree = (*pg)->get_vtree();
			size += vtree->memory_size();
		}

	}

	return size;
}

// public /////////////////////////////////////////////////////////////////////
PolygonGroup* Polylib::get_group(string name) const
{
#ifdef DEBUG
	PL_DBGOSH << "Polylib::get_group(" << name << ") in." << endl;
#endif
	vector<PolygonGroup*>::const_iterator it;
	for (it = m_pg_list.begin(); it != m_pg_list.end(); it++) {
		if (name == (*it)->acq_fullpath()) {
#ifdef DEBUG
			PL_DBGOS	<< "get_group: " << (*it)->get_parent_path()
						<< " name: " << (*it)->get_name()
						<< " size: " << (*it)->get_children().size() << endl;
#endif
			return *it;
		}
	}
	for (it = m_pg_list.begin(); it != m_pg_list.end(); it++) {
		if (name == (*it)->get_name()) {
			return *it;
		}
	}
#ifdef DEBUG
	PL_DBGOS << "Polylib::get_group(" << name << ") returns NULL" << endl;
#endif
	return NULL;
}

// protected //////////////////////////////////////////////////////////////////
Polylib::Polylib()
{
	gs_rankno = "";
	// デフォルトのファクトリークラスを登録する 2010.08.16
	m_factory = new PolygonGroupFactory();
}

// protected //////////////////////////////////////////////////////////////////
Polylib::~Polylib()
{
	vector<PolygonGroup*>::iterator it;
	for (it = m_pg_list.begin(); it != m_pg_list.end();) {
#ifdef DEBUG
PL_DBGOSH << "~Polylib():" << (*it)->get_name() << endl;
#endif
		(*it)->get_children().clear();
		delete *it;
		it = m_pg_list.erase(it);
	}
}

// protected //////////////////////////////////////////////////////////////////
POLYLIB_STAT Polylib::make_group_tree(
	PolylibConfig	*base
) {
#ifdef DEBUG
	PL_DBGOSH << "Polylib::make_group_tree() in." << endl;
#endif
	const PolylibCfgElem* root = base->get_root_elem();
	if (root == NULL) {
		PL_ERROSH << "[ERROR]Polylib::make_group_tree():Root tag not found." 
				  << endl;
		return PLSTAT_CONFIG_ERROR;
	}

	// ルート直下のノード分ループを回す
	const PolylibCfgElem* elem = root->first_element();
	if (elem == NULL) {
		PL_ERROSH << "[ERROR]Polylib::make_group_tree():Elem tag not found." 
				  << endl;
		return PLSTAT_CONFIG_ERROR;
	}

	while (elem != NULL) {
		const PolylibCfgParam	*param_class;
		PolygonGroup			*pg;
		string					class_name;

		// paramタグを取得し、PolygonGroupのインスタンスを作成
		param_class	= elem->first_param(PolygonGroup::ATT_NAME_CLASS);
		if (param_class == NULL) {
			PL_ERROSH << "[ERROR]Polylib::make_group_tree():Class name not found."
					  << endl;
			return PLSTAT_CONFIG_ERROR;
		}
		class_name  = param_class->get_string_data();
		pg			= m_factory->create_instance(class_name);
		add_pg_list(pg);	// グループリストに追加

		if (pg == NULL) {
			PL_ERROSH << "[ERROR]Polylib::make_group_tree():Class name not found."
					  << endl;
			return PLSTAT_CONFIG_ERROR;
		}

		// 配下のタグを取得して、PolygonGroupツリーを作成
		POLYLIB_STAT res = pg->build_group_tree(this, NULL, elem);
		if (res != PLSTAT_OK) return res;

		// 次ルートノードの取得
		elem = root->next_element(elem);
	}

	return PLSTAT_OK;
}

// protected //////////////////////////////////////////////////////////////////
POLYLIB_STAT Polylib::make_group_tree(
	string		config_contents
) {
#ifdef DEBUG
	PL_DBGOSH << "Polylib::make_group_tree() in." << endl;
#endif
	POLYLIB_STAT	stat;

	try {
		PolylibConfig base;
		stat = base.parse_xml_on_memory(config_contents);
		if (stat != PLSTAT_OK)	return PLSTAT_NG;

		return make_group_tree(&base);;
	}
	catch (POLYLIB_STAT e) {
		return e;
	}
}

// protected //////////////////////////////////////////////////////////////////
POLYLIB_STAT Polylib::load_config_file(
	string		*contents,
	string		fname
) {
#ifdef DEBUG
	PL_DBGOSH << "Polylib::load_config_file() in." << endl;
#endif
	// 設定ファイルを読み込みstring型で返す
	return PolylibConfig::load_config_file(contents, fname);
}

// protected //////////////////////////////////////////////////////////////////
POLYLIB_STAT Polylib::load_with_idfile(
	string		config_name,
	ID_FORMAT	id_format
) {
#ifdef DEBUG
	PL_DBGOSH << "Polylib::load_with_idfile() in." << endl;
#endif
	// 設定ファイル読み込み
	PolylibConfig base(config_name);

	// グループツリー作成
	POLYLIB_STAT stat = make_group_tree(&base);
	if (stat != PLSTAT_OK)	return stat;

	// STLファイルとIDファイル読み込み
	return load_polygons(true, id_format);
}

// protected //////////////////////////////////////////////////////////////////
POLYLIB_STAT Polylib::load_polygons(
	bool		with_id_file,
	ID_FORMAT	id_format
)
{
#ifdef DEBUG
	PL_DBGOSH << "Polylib::load_polygons() in." << endl;
#endif
	vector<PolygonGroup*>::iterator it;
	for (it = m_pg_list.begin(); it != m_pg_list.end(); it++) {
		// リーフの場合
		if ((*it)->get_children().empty() == true) {

			//STLファイルを読み込む
			POLYLIB_STAT ret = (*it)->load_stl_file();
			if (ret != PLSTAT_OK)		return ret;

			// 必要であればIDファイルを読み込んでm_idを設定
			if (with_id_file == true) {
				POLYLIB_STAT ret = (*it)->load_id_file(id_format);
				if (ret != PLSTAT_OK)		return ret;
			}
		}
	}
	return PLSTAT_OK;
}

// protected //////////////////////////////////////////////////////////////////
char *Polylib::save_config_file(
	string	rank_no,
	string	extend,
	string	format
){
#ifdef DEBUG
	PL_DBGOSH << "Polylib::save_config_file() in." << endl;
#endif
	vector<PolygonGroup*>::iterator	it;
	POLYLIB_STAT					stat;

	// DOM作成
	xmlDocPtr doc = xmlNewDoc((xmlChar *)"1.0");
	doc->children = PolylibConfig::mk_parameter_tag(doc);
	if (doc->children == NULL) return NULL;

	for (it = m_pg_list.begin(); it != m_pg_list.end(); it++) {
		if ((*it)->get_parent() == NULL) {
			// Elemタグ作成
			xmlNodePtr elem = PolylibConfig::mk_elem_tag(doc->children);
			if (elem == NULL)		return NULL;
			stat = mk_xml(elem, *it, rank_no, extend, format);
			if (stat != PLSTAT_OK)	return NULL;
		}
	}

	// ファイル出力
	char	*config_name;
	if (rank_no == "") {
		// ランク番号不要
		config_name = PolylibConfig::save_file(doc, "", extend);
	}
	else {
		// 要ランク番号(MPI版)
		config_name = PolylibConfig::save_file(doc, rank_no, extend);
	}
#ifdef DEBUG
PL_DBGOSH << "save_config_file():" << config_name << endl;
#endif
	xmlFreeDoc(doc);

	return config_name;
}

// protected //////////////////////////////////////////////////////////////////
POLYLIB_STAT Polylib::save_with_rankno(
	string		*p_config_name,
	int			myrank,
	int 		maxrank,
	string		extend,
	string		stl_format,
	ID_FORMAT	id_format
){
#ifdef DEBUG
	PL_DBGOSH << "Polylib::save_with_rankno() in." << endl;
#endif
	char	my_extend[128];

	// 拡張文字列がカラであれば、現在時刻から作成
	if (extend == "") {
		time_t		timer = time(NULL);
		struct tm	*date = localtime(&timer);
		sprintf(my_extend, "%04d%02d%02d%02d%02d%02d",
			date->tm_year+1900, date->tm_mon+1, date->tm_mday,
			date->tm_hour,      date->tm_min,   date->tm_sec);
	}
	else {
		sprintf(my_extend, "%s", extend.c_str());
	}

	// ランク番号の整形
	char	rank_no[16];
	int		fig = (int)log10((double)maxrank) + 1;
	sprintf(rank_no, "%0*d", fig, myrank);

	// 定義ファイルの保存
	char	*config_name = save_config_file(rank_no, my_extend, stl_format);
	if (config_name == NULL)	return PLSTAT_NG;
	else						*p_config_name = string(config_name);

	// STLファイルとIDファイルの保存
	vector<PolygonGroup*>::iterator	it;
	for (it = m_pg_list.begin(); it != m_pg_list.end(); it++) {
		POLYLIB_STAT	stat;

		//リーフのみがポリゴン情報を持っている
		if ((*it)->get_children().empty() == false) continue;

		// ポリゴン数が0ならばファイル出力不要 2010.10.19
		if ((*it)->get_triangles()->size() == 0) continue;

		stat = (*it)->save_stl_file(rank_no, my_extend, stl_format);
		if (stat != PLSTAT_OK)	return stat;
		stat = (*it)->save_id_file(rank_no, my_extend, id_format);
		if (stat != PLSTAT_OK)	return stat;
	}
	return PLSTAT_OK;
}

// protected //////////////////////////////////////////////////////////////////
// 2010.10.20 引数FILE *追加。
void Polylib::show_group_name(
	PolygonGroup	*pg, 
	string			tab,
	FILE			*fp
){
	vector<PolygonGroup*>::iterator it;

	// ユーザ定義id出力 2010.10.20
	if (fp == NULL) {
		PL_DBGOSH << "Polylib::show_group_name: " << tab ;
		if (pg->get_parent_path().empty() == true)  PL_DBGOS << "+";
		PL_DBGOS << pg->get_name() << ":" << pg->acq_file_name() << ":"
				 << pg->get_id() << endl;

	}
	else {
		// 出力先が指定されていれば、そちらに出力
		fprintf(fp, "%sPolylib::show_group_name:%s%s%s:%s:%d\n", 
			gs_rankno.c_str(),		tab.c_str(),
			(pg->get_parent_path().empty() == true) ? "+" : "",
			pg->get_name().c_str(),	pg->acq_file_name().c_str(),
			pg->get_id());
	}

	tab = tab + "    ";
	for (it = pg->get_children().begin(); it != pg->get_children().end(); it++){
		show_group_name(*it, tab, fp);
	}
}

// protected //////////////////////////////////////////////////////////////////
PolygonGroup* Polylib::get_group(int internal_id) const
{
#ifdef DEBUG
	PL_DBGOSH << "Polylib::get_group(" << internal_id << ") in." << endl;
#endif
	vector<PolygonGroup*>::const_iterator it;
	for (it = m_pg_list.begin(); it != m_pg_list.end(); it++) {
		if (internal_id == (*it)->get_internal_id()) {
			return *it;
		}
	}
#ifdef DEBUG
	PL_DBGOS << "Polylib::get_group(" << internal_id << ") returns NULL" << endl;
#endif
	return NULL;
}

// private ////////////////////////////////////////////////////////////////////
vector<PrivateTriangle*>* Polylib::search_polygons(
	string			group_name, 
	Vec3f			min_pos, 
	Vec3f			max_pos,
	bool			every,
	bool			linear, 
	POLYLIB_STAT	*ret
) const {
#ifdef DEBUG
	PL_DBGOSH << "Polylib::search_polygons() in." << endl;
#endif
	vector<PrivateTriangle*>* tri_list = new vector<PrivateTriangle*>;
	PolygonGroup* pg = get_group(group_name);

	if (pg == 0) {
		PL_ERROSH << "[ERROR]Polylib::search_polygons():Group not found: " 
				  << group_name << endl;
		*ret = PLSTAT_GROUP_NOT_FOUND;
		return tri_list;
	}
	vector<PolygonGroup*>* pg_list2 = new vector<PolygonGroup*>;

#ifdef BENCHMARK
	double st1, st2, ut1, ut2, tt1, tt2;
	bool ret1, ret2;
	ret1 = getrusage_sec(&ut1, &st1, &tt1);
#endif

	//子孫を検索
	search_group(pg, pg_list2);

	//自身を追加
	pg_list2->push_back(pg);

	// 検索範囲
	BBox bbox;
	bbox.init();
	bbox.add(min_pos);
	bbox.add(max_pos);

	//全ポリゴングループを検索
	vector<PolygonGroup*>::iterator it;
	for (it = pg_list2->begin(); it != pg_list2->end(); it++) {

		//リーフ構造からのみ検索を行う
		if ((*it)->get_children().size()==0) {
			POLYLIB_STAT ret2;
#ifdef DEBUG
			if (linear == false) {
#endif
				ret2 = (*it)->search (&bbox,every,tri_list);
				if (ret2 != PLSTAT_OK) {
					*ret = ret2;
					return tri_list;
				}
#ifdef DEBUG
			}
			else{
				ret2 = (*it)->linear_search (&bbox,every,tri_list);
				if (ret2 != PLSTAT_OK) {
					*ret = ret2;
					return tri_list;
				}
			}
#endif
		}
	}
#ifdef BENCHMARK
	ret2 = getrusage_sec(&ut2,&st2,&tt2);
	if (ret1 == false || ret2 == false) {
		PL_DBGOSH << "Search SYS   Time Error" << endl;
		PL_DBGOSH << "Search CPU   Time Error" << endl;
		PL_DBGOSH << "Search Total Time Error" << endl;
	}
	else{
		cout.setf(ios::scientific, ios::floatfield);
		PL_DBGOSH << "Search SYS   Time:" << st2 - st1 << endl;
		PL_DBGOSH << "Search CPU   Time:" << ut2 - ut1 << endl;
		PL_DBGOSH << "Search Total Time:" << tt2 - tt1 << endl;
		std::cout.unsetf(ios::scientific);
	}
#endif

	delete pg_list2;
	*ret = PLSTAT_OK;
	return tri_list;
}

// private ////////////////////////////////////////////////////////////////////
void Polylib::search_group(
	PolygonGroup			*p, 
	vector<PolygonGroup*>	*pg
) const {
#ifdef DEBUG
	PL_DBGOSH << "Polylib::search_group() in." << endl;
#endif
	vector<PolygonGroup*>::iterator it;
	for (it = p->get_children().begin(); it != p->get_children().end(); it++) {
		pg->push_back(*it);
		search_group(*it, pg);
	}
}

// private //////////////////////////////////////////////////////////////////
POLYLIB_STAT Polylib::mk_xml(
	xmlNodePtr			elem,
	PolygonGroup		*pg,
	string				rank_no,
	string				extend,
	string				format
)
{
#ifdef DEBUG
	PL_DBGOSH << "Polylib::mk_xml() in." << endl;
#endif
	// Paramタグ作成
	pg->mk_param_tag(elem, rank_no, extend, format);

	vector<PolygonGroup*>::iterator it;
	POLYLIB_STAT					stat;
	for (it = pg->get_children().begin(); it != pg->get_children().end(); it++){
		// Elemタグ作成
		xmlNodePtr elem2 = PolylibConfig::mk_elem_tag(elem);
		if (elem2 == NULL)			return PLSTAT_NG;
		// Elemタグ配下のParamタグを作成する
		stat = mk_xml(elem2, *it, rank_no, extend, format);
		if (stat != PLSTAT_OK)		return stat;
			
	}
	return PLSTAT_OK;
}

