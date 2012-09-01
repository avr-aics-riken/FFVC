/*
 * Polylib - Polygon Management Library.
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2010-
 *
 */
#include <fstream>
#include <map>
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

/// TextParser
// public /////////////////////////////////////////////////////////////////////
// 
POLYLIB_STAT Polylib::load(
	string	config_name
) {
#ifdef DEBUG
	PL_DBGOSH << "Polylib::load_test() in." << endl;
#endif

	// 設定ファイル読み込み
	try {
	  //PolylibConfig base(config_name);

	  tp->read(config_name);
	  //  tp->write("tmp.tpp");
	  // グループツリー作成
	  POLYLIB_STAT stat = make_group_tree(tp);
	  if (stat != PLSTAT_OK)	return stat;

	  // STLファイル読み込み (三角形IDファイルは不要なので、第二引数はダミー)
	  return load_polygons(false, ID_BIN);
	}
	catch (POLYLIB_STAT e) {
		return e;
	}
}

// public /////////////////////////////////////////////////////////////////////
//TextParser 版

POLYLIB_STAT Polylib::save(
	string	*p_config_name,
	string	stl_format,
	string	extend
) {
  //#ifdef DEBUG
	PL_DBGOSH << "Polylib::save() in." << endl;
	//#endif
	char	my_extend[128];
	POLYLIB_STAT stat=PLSTAT_OK;


	// 拡張文字列がカラであれば、現在時刻から作成
	if (extend == "") {
	  time_t timer = time(NULL);
	  struct tm	*date = localtime(&timer);
	  sprintf(my_extend, "%04d%02d%02d%02d%02d%02d",
		  date->tm_year+1900,date->tm_mon+1,date->tm_mday,
		  date->tm_hour,date->tm_min,date->tm_sec);
	}
	else {
	  sprintf(my_extend, "%s", extend.c_str());
	}
	
	//PL_DBGOSH << __FUNCTION__ << " extend "<< my_extend << endl;

	map<string,string> stl_fname_map;

	vector<PolygonGroup*>::iterator	it;
	for (it = m_pg_list.begin(); it != m_pg_list.end(); it++) {
	  //リーフのみがポリゴン情報を持っている
	  if ((*it)->get_children().empty() == false)	continue;

	  // ポリゴン数が0ならばファイル出力不要 2010.10.19
	  if ((*it)->get_triangles()->size() == 0)	continue;

	  // STLファイル保存 (第一引数のランク番号は不要)
	  stat = (*it)->save_stl_file("", my_extend, stl_format,
						   stl_fname_map);
	  if (stat != PLSTAT_OK) return stat;

	  stat = (*it)->mk_param_tag(tp, "", "", "");
	  if (stat != PLSTAT_OK) return stat;
	}


	// update stl filepath 
	// clear file path first
	stat=clearfilepath(tp);
	//tp->write("tmp2.tpp");
	// set filepath
	stat=setfilepath(stl_fname_map);

	//	string tmp_extend = my_extend;
	//	char	*config_name = save_config_file("", tmp_extend, stl_format);
	char	*config_name = save_config_file("", my_extend, stl_format);
	//	PL_DBGOSH << __FUNCTION__ << " config_name "<< config_name << endl;

	if (config_name == NULL)	return PLSTAT_NG;
	else	*p_config_name = string(config_name);
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

	
	//Polylib にTextParser クラスを持たせる。
	tp = new TextParser;
	
	
	//PL_DBGOS<< __FUNCTION__ <<" m_factory "<< m_factory << " tp " << tp<<std::endl;

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
	if(tp !=0) delete tp;

}

// protected //////////////////////////////////////////////////////////////////
// TextParser 版
POLYLIB_STAT Polylib::make_group_tree(
    TextParser* tp
) {
#ifdef DEBUG
	PL_DBGOSH << "Polylib::make_group_tree(TextParser) in." << endl;
#endif
	TextParserError status = TP_NO_ERROR;

	// 念のため階層構造の最上位へ
	string cur;	
	tp->currentNode(cur);
	if(cur != "/Polylib") {
	  status=tp->changeNode("/Polylib");
	  tp->currentNode(cur);

	  if(status!=TP_NO_ERROR){
	    PL_ERROSH << 
	      "[ERROR]Polylib::make_group_tree(TextParser):Root node not found."
		      << endl;
	    return PLSTAT_CONFIG_ERROR;
	  }

	}
#if 0 
	if(cur != "/") {
	  status=tp->changeNode("/");
	  tp->currentNode(cur);

	  if(status!=TP_NO_ERROR){
	    PL_ERROSH << 
	      "[ERROR]Polylib::make_group_tree(TextParser):Root node not found."
		      << endl;
	    return PLSTAT_CONFIG_ERROR;
	  }

	}
#endif

	// ノードとリーフのリストを取る

	//ノードを取得
	vector<string> nodes;
	tp->getNodes(nodes);
	string current_node;


	//	tp->currentNode(current_node);
	// vector<string> leaves;
	// tp->getLabels(leaves);
	// if(nodes.size()==0 && leaves.size()==0){ return PLSTAT_CONFIG_ERROR;}

	// string class_name = "PolygonGroup";
	// if(leaves.size()!=0){
	//   vector<string>::iterator leaf_iter=find(leaves.begin(),
	// 					  leaves.end(),
	// 					  PolygonGroup::ATT_NAME_CLASS);
	//   if(leaf_iter != leaves.end()){
	//     //	    class_name = *leaf_iter;
	//     string value;
	//     status=tp->getValue((*leaf_iter),value);
	//     class_name=value;

	//   }
	// }


	// loop over node recurcively.
	PolygonGroup *pg;
	for(vector<string>::iterator nodes_iter = nodes.begin(); 
	    nodes_iter != nodes.end();
	    nodes_iter++){

	  status = tp->changeNode(*nodes_iter);
	  tp->currentNode(current_node);
#ifdef DEBUG	  
	  PL_DBGOS<<"current_node "<< current_node <<  endl;
#endif // DEBUG	  
	  vector<string> nodes;
	  vector<string> leaves;

	  tp->getNodes(nodes);
	  tp->getLabels(leaves);
	  if(nodes.size()==0 && leaves.size()==0){ return PLSTAT_CONFIG_ERROR;}

	  string class_name = "PolygonGroup"; //default

	  if(leaves.size()!=0){
	    vector<string>::iterator leaf_iter=find(leaves.begin(),
						    leaves.end(),
						    PolygonGroup::ATT_NAME_CLASS);
	    if(leaf_iter != leaves.end()){
	      //     class_name = *leaf_iter;
	      string value;
	      status=tp->getValue((*leaf_iter),value);
	      class_name=value;
	    }
	  }
	  pg = m_factory->create_instance(class_name);
	  add_pg_list(pg);
	  if (pg == NULL) {
	    PL_ERROSH << "[ERROR]Polylib::make_group_tree():Class name not found."
		      << class_name
		      << endl;
	    return PLSTAT_CONFIG_ERROR;
	  }


	  // 配下のタグを取得して、PolygonGroupツリーを作成
	  //	  POLYLIB_STAT res = pg->build_group_tree(this, NULL, elem);
	  //  
	  POLYLIB_STAT res = pg->build_group_tree(this, NULL, tp);
	  if (res != PLSTAT_OK) return res;
	  status = tp->changeNode("..");

	} 
	

	return PLSTAT_OK;
}

//TextParser Version
// protected //////////////////////////////////////////////////////////////////

POLYLIB_STAT Polylib::make_group_tree(
	string		config_contents
) {
#ifdef DEBUG
	PL_DBGOSH << "Polylib::make_group_tree() in." << endl;
#endif
	try {
	 tp->read(config_contents);
	  return make_group_tree(tp);
	}
	catch (POLYLIB_STAT e) {
		return e;
	}
	return PLSTAT_NG;
}

// TextParser 版
// protected //////////////////////////////////////////////////////////////////
POLYLIB_STAT Polylib::load_config_file(
	string		*contents,
	string		fname
) {
#ifdef DEBUG
 	PL_DBGOSH << "Polylib::load_config_file() in." << endl;
#endif
	// 設定ファイルを読み込みstring型で返す
	//	return PolylibConfig::load_config_file(contents, fname);
	//	return PolylibConfig::load_config_file(contents, fname);

	// right now just return no error
	return PLSTAT_OK;

}
/// textparser 版
// protected //////////////////////////////////////////////////////////////////
POLYLIB_STAT Polylib::load_with_idfile(
	string		config_name,
	ID_FORMAT	id_format
) {
#ifdef DEBUG
	PL_DBGOSH << "Polylib::load_with_idfile() in." << endl;
#endif
	// 設定ファイル読み込み
	tp->read_local(config_name);

	// グループツリー作成
	
	POLYLIB_STAT stat = make_group_tree(tp);
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
//TextPArser 版
char *Polylib::save_config_file(
	string	rank_no,
	string	extend,
	string	format
){
#ifdef DEBUG
  //  PL_DBGOSH << "Polylib::save_config_file() in. " << endl;
#endif

	vector<PolygonGroup*>::iterator	it;
	//POLYLIB_STAT					stat;

	// ファイル出力
	char	*config_name;

	if (rank_no == "") {
		// ランク番号不要
	  config_name = polylib_config_save_file("", extend);
	  //config_name = PolylibConfig::save_file(doc, "", extend);
	}
	else {
	  // 要ランク番号(MPI版)
	  //	config_name = PolylibConfig::save_file(doc, rank_no, extend);
	  config_name = polylib_config_save_file(rank_no, extend);
	}

#ifdef DEBUG
	//PL_DBGOSH << "save_config_file(): " << config_name << endl;
#endif
//	xmlFreeDoc(doc);

	return config_name;
}

/////////////////////////////////////////　　
// 追加　２０１２ー０８
//////////////////////////////////

char* Polylib::polylib_config_save_file(
					string rank_no,
					string extend){

#define POLYLIB_CONFIG_NAME		"polylib_config"  
#define POLYLIB_CONFIG_EXT		"tpp"

  //  cout <<__FUNCTION__<<" in "<< rank_no <<" "<< extend << endl;

  static char fname[1024];
  if (rank_no == ""){
    sprintf(fname,"%s_%s.%s",POLYLIB_CONFIG_NAME,extend.c_str(),POLYLIB_CONFIG_EXT);
  } else {
    sprintf(fname, "%s_%s_%s.%s", POLYLIB_CONFIG_NAME, rank_no.c_str(),
	    extend.c_str(), POLYLIB_CONFIG_EXT);
  }
  string fname_string = fname;

  //config_file の書き込み

  tp->write(fname_string,1);

  return fname;

}


//TextParser version 
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
	PL_DBGOSH << "Polylib::save_with_rankno() in. " << endl;
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

#ifdef DEBUG
	PL_DBGOSH << "Polylib::save_with_rankno() rank_no" << rank_no<< endl;
#endif

	map<string,string> stl_fname_map;
	// STLファイルとIDファイルの保存
	vector<PolygonGroup*>::iterator	it;
	for (it = m_pg_list.begin(); it != m_pg_list.end(); it++) {
		POLYLIB_STAT	stat;
		//リーフのみがポリゴン情報を持っている
		if ((*it)->get_children().empty() == false) continue;

		// ポリゴン数が0ならばファイル出力不要 2010.10.19
		if ((*it)->get_triangles()->size() == 0) continue;

		//stat = (*it)->save_stl_file(rank_no, my_extend, stl_format);
		stat = (*it)->save_stl_file(rank_no, my_extend, stl_format,stl_fname_map);
		if (stat != PLSTAT_OK)	return stat;
		stat = (*it)->save_id_file(rank_no, my_extend, id_format);
		if (stat != PLSTAT_OK)	return stat;
		string rank_string,my_extend_string;
		rank_string=rank_no;
		my_extend_string = my_extend;
		stat = (*it)->mk_param_tag(tp, "", "", "");
		if (stat != PLSTAT_OK) return stat;
	  
	}
	tp->changeNode("/"); //
	//	cout << "before cleanfilepath" <<endl;
	clearfilepath(tp); //clear whole file path
	//	cout << "after cleanfilepath" <<endl;

	//	cout << "before setfilepath" <<endl;
	setfilepath(stl_fname_map);
	//cout << "after setfilepath" <<endl;

	// 定義ファイルの保存	


	char	*config_name = save_config_file(rank_no, my_extend, stl_format);
#ifdef DEBUG	
	PL_DBGOSH << __FUNCTION__ << " config_name "<< config_name << endl;
#endif 
#ifdef DEBUG
	PL_DBGOSH << "Polylib::save_with_rankno() config_name " << config_name << endl;
#endif
	if (config_name == NULL)	return PLSTAT_NG;
	else	*p_config_name = string(config_name);
	return PLSTAT_OK;

}

/////////////////////////////////////////　　
// 追加　２０１２ー０８
//////////////////////////////////

POLYLIB_STAT Polylib::setfilepath(map<string,string>& stl_fname_map){

#ifdef DEBUG
  PL_DBGOS << "stl_map size " <<  stl_fname_map.size()<<endl;
#endif //DEBUG
  //  tp->changeNode("/");
  tp->changeNode("/Polylib");
  for(map<string,string>::iterator map_iter=stl_fname_map.begin();
      map_iter != stl_fname_map.end();
      map_iter++){
#ifdef DEBUG
    PL_DBGOS << "stl_map " <<  map_iter->first <<" "<< map_iter->second<<endl;
#endif // DEBUG
    tp->changeNode(map_iter->first); //
    string cur;
    tp->currentNode(cur);
    
    //    cout << __FUNCTION__ << " " <<  cur << endl;
    string value = "\"" + map_iter->second + "\"";
    //    cout << "before leaf createion" <<endl;
    tp->createLeaf("filepath",value);
    //    cout << "after leaf createion" <<endl;
    //check
    //    vector<string> leaves;
    //   tp->getLabels(leaves);  //label　の取り直し
    //    for(vector<string>::iterator leaf_iter=leaves.begin();
    //	leaf_iter!=leaves.end();
    //	leaf_iter++){
    //      cout << __FUNCTION__ << " "  << *leaf_iter << endl;
    //    }
    //  	  string value;
    //    tp->getValue("filepath",value);
    //    cout << __FUNCTION__ << " " << cur << " "<< value <<endl;

    //    tp->changeNode("/");

    tp->changeNode("/Polylib");
  }
  return PLSTAT_OK;
}
/////////////////////////////////////////　　
// 追加　２０１２ー０８
//////////////////////////////////
POLYLIB_STAT Polylib::clearfilepath(TextParser* tp_ptr){

  // recursive にするため、TextParserのポインタを引数に取る。

  vector<string> leaves;
  tp_ptr->getLabels(leaves,1);

  vector<string>::iterator leaf_iter=find(leaves.begin(),leaves.end(),"filepath");
  if(leaf_iter!=leaves.end()){ // 見つかったら
    tp_ptr->deleteLeaf("filepath");
  } 
  leaves.clear();
  tp_ptr->getLabels(leaves);


  leaf_iter=find(leaves.begin(),leaves.end(),"filepath[0]");
  if(leaf_iter!=leaves.end()){ // 見つかったら

    int index=0;
    leaf_iter = leaves.begin();
    while(leaf_iter != leaves.end()){
      stringstream ss;
      string tmpstring="filepath";
      ss << tmpstring <<"["<<index<<"]";
      ss >> tmpstring;
      leaf_iter = find(leaf_iter,leaves.end(),tmpstring);
      if(leaf_iter == leaves.end()) break;
      TextParserError tp_error=tp_ptr -> deleteLeaf(tmpstring);
      if (tp_error!=TP_NO_ERROR) {
	PL_ERROSH << "[ERROR]Polylib::save() "
		  << "can not remove leaf = " << tmpstring << endl;
	return PLSTAT_NG;
      }
      index++;
      leaf_iter++;
    }
  }

  vector<string> nodes;
  tp_ptr->getNodes(nodes);
  for( vector<string>::iterator node_iter=nodes.begin();
       node_iter!=nodes.end();
       node_iter++){
    tp_ptr->changeNode(*node_iter);
    clearfilepath(tp_ptr);
    tp_ptr->changeNode("..");
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
	// ユーザ定義ラベル出力 2012.08.31
	if (fp == NULL) {
		PL_DBGOSH << "Polylib::show_group_name: " << tab ;
		if (pg->get_parent_path().empty() == true)  PL_DBGOS << "+";
		PL_DBGOS << pg->get_name() << ":" << pg->acq_file_name() << ":"
				 << pg->get_id() << ":" << pg->get_label() << endl;

	}
	else {
		// 出力先が指定されていれば、そちらに出力
		fprintf(fp, "%sPolylib::show_group_name:%s%s%s:%s:%d:%s\n", 
			gs_rankno.c_str(),		tab.c_str(),
			(pg->get_parent_path().empty() == true) ? "+" : "",
			pg->get_name().c_str(),	pg->acq_file_name().c_str(),
			pg->get_id(), pg->get_label().c_str());
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

