/*
 * Polylib - Polygon Management Library.
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2010-
 *
 */
#include <iostream>
#include <sstream>
#include <vector>
#include "file_io/PolylibConfig.h"

namespace PolylibNS {

using namespace std;

#define POLYLIB_CONFIG_NAME		"polylib_config"
#define POLYLIB_CONFIG_EXT		"xml"
#define PARAMETER				"Parameter"
#define ELEM					"Elem"
#define PARAM					"Param"
#define ELEM_NAME_POLYGONGROUP	"PolygonGroup"
#define ATT_NAME				"name"
#define ATT_DTYPE				"dtype"
#define ATT_VALUE				"value"
#define DTYPE_STR				"STRING"
#define DTYPE_INT				"INT"
#define DTYPE_REAL				"REAL"

/************************************************************************
 *
 * PolylibConfigクラス
 *
 ***********************************************************************/
// public /////////////////////////////////////////////////////////////////////
PolylibConfig::PolylibConfig() {
}

PolylibConfig::PolylibConfig(
	string		fname
) {
	xmlDocPtr doc;

	// 設定ファイル読み込み
#ifdef DEBUG
PL_DBGOSH << "PolylibConfig::PolylibConfig:" << fname << endl;
#endif
	doc = xmlParseFile(fname.c_str());
	if (doc == NULL) {
		PL_ERROSH << "[ERROR]PolylibConfig::PolylibConfig():Can't read config"
				  << " file:" << fname << endl;
		throw PLSTAT_NG;
	}

	// ルートノードを取得 & XMLを解析
	xmlNodePtr cur = xmlDocGetRootElement(doc);
	if (cur == NULL) {
		PL_ERROSH << "[ERROR]PolylibConfig::PolylibConfig():Can't get root"
				  << " node:" << fname << endl;
		xmlFreeDoc(doc);
		throw PLSTAT_NG;
	}
	m_elem = new PolylibCfgElem("");
	if (xml_parse(m_elem, cur) != PLSTAT_OK) {
		throw PLSTAT_NG;
	}

	xmlFreeDoc(doc);
}

// public /////////////////////////////////////////////////////////////////////
PolylibConfig::~PolylibConfig()
{
	delete m_elem;
}

// public /////////////////////////////////////////////////////////////////////
POLYLIB_STAT PolylibConfig::parse_xml_on_memory(
	string		contents
) {
	// DOM生成
	xmlDocPtr doc = xmlParseMemory((char *)(contents.c_str()), 
													strlen(contents.c_str()));
	if (doc == NULL) {
		PL_ERROSH << "[ERROR]PolylibConfig::parse_xml_on_memory():Can't create "
				  << "xml doc" << endl;
		return PLSTAT_NG;
	}

	// ルートノードを取得 & XMLを解析
	xmlNodePtr cur = xmlDocGetRootElement(doc);
	if (cur == NULL) {
		PL_ERROSH << "[ERROR]PolylibConfig::parse_xml_on_memory():Can't get "
				  << "root node." << endl;
		xmlFreeDoc(doc);
		return PLSTAT_NG;
	}
	m_elem = new PolylibCfgElem("");
	if (xml_parse(m_elem, cur) != PLSTAT_OK) {
		xmlFreeDoc(doc);
		return PLSTAT_NG;
	}
	else {
		xmlFreeDoc(doc);
		return PLSTAT_OK;
	}
}

// public /////////////////////////////////////////////////////////////////////
POLYLIB_STAT PolylibConfig::load_config_file(
	string		*contents,
	string		fname
) {
	xmlDocPtr	doc;
	xmlChar		*buff;
	int			size;

	// 設定ファイル読み込み
	doc = xmlParseFile(fname.c_str());
	if (doc == NULL) {
		PL_ERROSH << "[ERROR]PolylibConfig::load_config_file():Can't read "
				  << fname << endl;
		return PLSTAT_NG;
	}
	xmlDocDumpMemory(doc, &buff, &size);
	*contents = string((char *)buff);
	free(buff);
	xmlFreeDoc(doc);

	return PLSTAT_OK;
}

// public /////////////////////////////////////////////////////////////////////
xmlNodePtr PolylibConfig::mk_parameter_tag(
	xmlDocPtr	doc
){
	xmlNodePtr node;
	node = xmlNewDocNode(doc, NULL, (xmlChar *)PARAMETER, NULL);
	if (node == NULL) {
		PL_ERROSH << "[ERROR]PolylibConfig::mk_parameter_tag():" 
				  << "xmlNewDocNode return NULL." << endl;
	}
	return node;
}

// public /////////////////////////////////////////////////////////////////////
xmlNodePtr PolylibConfig::mk_elem_tag(
	xmlNodePtr	elem
) {
	xmlNodePtr	elem2 = xmlNewChild(elem, NULL, (xmlChar *)ELEM, NULL);
	if (elem2 == NULL) {
		PL_ERROSH << "[ERROR]PolylibConfig::mk_elem_tag():" 
				  << "xmlNewChild return NULL." << endl;
		return elem2;
	}
	xmlSetProp(elem2, (xmlChar *)ATT_NAME, (xmlChar *)ELEM_NAME_POLYGONGROUP);
	return elem2;
}

// public /////////////////////////////////////////////////////////////////////
POLYLIB_STAT PolylibConfig::mk_param_tag(
	xmlNodePtr		elem,
	string 			name,
	string 			value
) {
	xmlNodePtr  param = xmlNewChild(elem, NULL, (xmlChar *)PARAM, NULL);
	if (param == NULL) {
		PL_ERROSH << "[ERROR]PolylibConfig::mk_param_tag():xmlNewChild return "
				  << "NULL." << endl;
		return PLSTAT_NG;
	}
	xmlSetProp(param, (xmlChar *)ATT_NAME,  (xmlChar *)(name.c_str()));
	xmlSetProp(param, (xmlChar *)ATT_DTYPE, (xmlChar *)DTYPE_STR);
	xmlSetProp(param, (xmlChar *)ATT_VALUE, (xmlChar *)(value.c_str()));
	return PLSTAT_OK;
}

// public /////////////////////////////////////////////////////////////////////
POLYLIB_STAT PolylibConfig::mk_param_tag(
	xmlNodePtr		elem,
	string			name,
	int				value
) {
	char	ival[32];
	sprintf(ival, "%d", value);
	xmlNodePtr  param = xmlNewChild(elem, NULL, (xmlChar *)PARAM, NULL);
	if (param == NULL) {
		PL_ERROSH << "[ERROR]PolylibConfig::mk_param_tag():xmlNewChild return "
				  << "NULL." << endl;
		return PLSTAT_NG;
	}
	xmlSetProp(param, (xmlChar *)ATT_NAME,  (xmlChar *)(name.c_str()));
	xmlSetProp(param, (xmlChar *)ATT_DTYPE, (xmlChar *)DTYPE_INT);
	xmlSetProp(param, (xmlChar *)ATT_VALUE, (xmlChar *)ival);
	return PLSTAT_OK;
}

// public /////////////////////////////////////////////////////////////////////
POLYLIB_STAT PolylibConfig::mk_param_tag(
	xmlNodePtr	elem,
	string		name,
	double		value
) {
	char	fval[32];
	sprintf(fval, "%f", value);
	xmlNodePtr  param = xmlNewChild(elem, NULL, (xmlChar *)PARAM, NULL);
	if (param == NULL) {
		PL_ERROSH << "[ERROR]PolylibConfig::mk_param_tag():xmlNewChild return "
				  << "NULL." << endl;
		return PLSTAT_NG;
	}
	xmlSetProp(param, (xmlChar *)ATT_NAME,  (xmlChar *)(name.c_str()));
	xmlSetProp(param, (xmlChar *)ATT_DTYPE, (xmlChar *)DTYPE_REAL);
	xmlSetProp(param, (xmlChar *)ATT_VALUE, (xmlChar *)fval);
	return PLSTAT_OK;
}

// public /////////////////////////////////////////////////////////////////////
char *PolylibConfig::save_file(
	xmlDocPtr	doc,
	string 		rank_no,
	string		extend
) {
	static char	fname[1024];
	if (rank_no == "") {
		sprintf(fname, "%s_%s.%s", POLYLIB_CONFIG_NAME, extend.c_str(), 
														POLYLIB_CONFIG_EXT);
	}
	else {
		sprintf(fname, "%s_%s_%s.%s", POLYLIB_CONFIG_NAME, rank_no.c_str(), 
											extend.c_str(), POLYLIB_CONFIG_EXT);
	}
	// 整形して出力しよう
	//xmlSaveFile("test.xml", doc);
	if (xmlSaveFormatFile(fname, doc, xmlIndentTreeOutput) == -1) {
		PL_ERROSH << "[ERROR]PolylibConfig::save_file():xmlSaveFormatFile:"
				  << fname << endl;
		return NULL;
	}
	return fname;
}

// private ////////////////////////////////////////////////////////////////////
POLYLIB_STAT PolylibConfig::xml_parse(
	PolylibCfgElem	*parent,
	xmlNodePtr		node
) {
	xmlNodePtr cur = NULL;

	for (cur = node; cur; cur = cur->next) {
		ostringstream oss;
		oss << cur->name;
		if (oss.str() == ELEM) {
			xmlChar* name = xmlGetProp(cur, (const xmlChar *)ATT_NAME);
			ostringstream oss_name;
			oss_name << name;
#ifdef DEBUG
PL_DBGOSH << "PolylibConfig::xml_parse():Elem:name=" << name << endl;
#endif
			if (oss_name.str() != ELEM_NAME_POLYGONGROUP) {
				continue;
			}
			PolylibCfgElem* elem = new PolylibCfgElem(oss_name.str());
			if (xml_parse(elem, cur->children) != PLSTAT_OK) {
				return PLSTAT_NG;
			}
			parent->set_elem(elem);
		}
		else if (oss.str() == PARAM) {
			ostringstream oss_name;
			ostringstream oss_type;
			ostringstream oss_value;
			xmlChar* name	= xmlGetProp(cur, (const xmlChar *)ATT_NAME);
			xmlChar* dtype	= xmlGetProp(cur, (const xmlChar *)ATT_DTYPE);
			xmlChar* value	= xmlGetProp(cur, (const xmlChar *)ATT_VALUE);
			oss_name	<< name;
			oss_type	<< dtype;
			oss_value	<< value;
#ifdef DEBUG
PL_DBGOSH << "PolylibConfig::xml_parse():Param:name,dtype,value=" << name << "," << dtype << "," << value << endl;
#endif
			PolylibCfgParam* param = new PolylibCfgParam(
				oss_name.str(), oss_type.str(), oss_value.str());
			PolylibCfgParamType	data_type = param->get_data_type();
			if (data_type != INT && data_type != REAL && data_type != STRING) {
				PL_ERROSH << "[ERROR]PolylibConfig::xml_parse():Unexpected dtype"
						  << ":name,dtype,value=" << name << "," << dtype << ","
						  << value << endl;
				return PLSTAT_NG;
			}
			parent->set_param(param);
		}
		else if (oss.str() == PARAMETER) {
			if (xml_parse(parent, cur->children) != PLSTAT_OK) {
				return PLSTAT_NG;
			}
		}
// 不要 2010.07.15
//		else if (oss.str() == "SphereConfig") {
//			xml_parse(parent, cur->children);
//		}
	}
	return PLSTAT_OK;
}

/************************************************************************
 *
 * PolylibCfgElemクラス
 *
 ***********************************************************************/
// public /////////////////////////////////////////////////////////////////////
PolylibCfgElem::PolylibCfgElem(
	const string name
) {
	m_elem  = new vector<PolylibCfgElem*>;
	m_param = new vector<PolylibCfgParam*>;
	m_name = name;
}

// public /////////////////////////////////////////////////////////////////////
PolylibCfgElem::~PolylibCfgElem()
{
	vector<PolylibCfgParam*>::const_iterator pitr;
	for (pitr = m_param->begin(); pitr != m_param->end(); pitr++) {
		delete *pitr;
	}
	delete m_param;

	vector<PolylibCfgElem*>::const_iterator eitr;
	for (eitr = m_elem->begin(); eitr != m_elem->end(); eitr++) {
		delete *eitr;
	}
	delete m_elem;
}

// public /////////////////////////////////////////////////////////////////////
const PolylibCfgElem* PolylibCfgElem::first_element (
	const string name
) const {
	if (m_elem->size() == 0) {
		return NULL;
	}
	vector<PolylibCfgElem*>::const_iterator itr = m_elem->begin();
	if (name.empty() == true) {
		return *itr;
	}
	for (itr = m_elem->begin(); itr != m_elem->end(); itr++) {
		if ((*itr)->get_name() == name) {
			return *itr;
		}
	}
	return NULL;
}

// public /////////////////////////////////////////////////////////////////////
const PolylibCfgElem* PolylibCfgElem::next_element (
	const PolylibCfgElem	*elem, 
	const string			name
) const {
	if (elem == NULL) return NULL;

	vector<PolylibCfgElem*>::const_iterator itr = m_elem->begin();
	bool flg = false;
	for (itr = m_elem->begin(); itr != m_elem->end(); itr++) {
		if (*itr == elem) {
			flg = true;
			continue;
		}
		if (flg == false) {
			continue;
		}
		if (name.empty() == true && (*itr)->get_name() == elem->get_name()) {
			return *itr;
		}
		else if (name.empty() == false && (*itr)->get_name() == name) {
			return *itr;
		}
	}
	return NULL;
}

// public /////////////////////////////////////////////////////////////////////
const PolylibCfgParam* PolylibCfgElem::first_param (
	const string	name
) const {
	if (m_param->size() == 0) {
		return NULL;
	}
	vector<PolylibCfgParam*>::const_iterator itr = m_param->begin();
	if (name.empty() == true) {
		return *itr;
	}
	for (itr = m_param->begin(); itr != m_param->end(); itr++) {
		if ((*itr)->get_name() == name) {
			return *itr;
		}
	}
	return NULL;
}

// public /////////////////////////////////////////////////////////////////////
const PolylibCfgParam* PolylibCfgElem::next_param (
	const PolylibCfgParam	*param,
	const string			name
) const {
	if (param == NULL) return NULL;

	vector<PolylibCfgParam*>::const_iterator itr = m_param->begin();
	bool flg = false;
	for (itr = m_param->begin(); itr != m_param->end(); itr++) {
		if (*itr == param) {
			flg = true;
			continue;
		}
		if (flg == false) {
			continue;
		}
		if (name.empty() == true && (*itr)->get_name() == param->get_name()) {
			return *itr;
		}
		else if (name.empty() == false && (*itr)->get_name() == name) {
			return *itr;
		}
	}
	return NULL;
}

/************************************************************************
 *
 * PolylibCfgParamクラス
 *
 ***********************************************************************/
// public /////////////////////////////////////////////////////////////////////
PolylibCfgParam::PolylibCfgParam(
	const string	name, 
	const string	type, 
	const string	value
) {
	char		capital[256];
	const char	*c_str = type.c_str();

	// 大文字でも小文字でもOK 2010.10.15
	memset(capital, '\0', sizeof(capital));
	for (int i = 0; c_str[i]; i++) {
		capital[i] = toupper(c_str[i]);
	}
	m_name = name;
	if (!strcmp(capital, DTYPE_STR)) {
		m_type	= STRING;
		m_str	= value;
	}
	else if (!strcmp(capital, DTYPE_INT)) {
		m_type	= INT;
		m_int	= atoi(value.c_str());
	}
	else if (!strcmp(capital, DTYPE_REAL)) {
		m_type	= REAL;
		m_real	= atof(value.c_str());
	}
}
}
