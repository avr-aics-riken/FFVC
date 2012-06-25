/****************************************************************************
**
** Copyright (C) 2012 Tokyo University.
**
****************************************************************************/

/** @file TextParser.cpp
 * ここには TextParser クラスと
 * C言語用APIが実装されています。
 *
 */

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif //HAVE_CONFIG_H

#include <string.h>

#include "TextParser.h"

/** 
 * プロセス内唯一の TextParser インスタンスを返します。
 *
 * @return TextParser インスタンスのアドレス
 * 
 */

TextParser* TextParser::get_instance(){
#ifdef MYDEBUG
  std::cout<< "TextParser::get_instance() called."<<std::endl;
#endif
  static TextParser instance;
  return &instance;
}

/** 
 * デフォルトコンストラクタで、
 * パラメータのデータ構造へのインスタンスを取得します。
 *
 *
 */
TextParser::TextParser(){
  if( dataTree()==0) _data_tree=TextParserTree::get_instance();
}

/** 指定されたパラメータファイルを読み込み、字句解析を行って
 *  その結果をtree構造データとして格納します。
 *
 * @param[in] file 入力するパラメータファイル名
 * 
 * @return エラーコード 
 */

TextParserError TextParser::read(const std::string& file){
  //std::cout << file << std::endl;
  
  TextParserError ret = TP_NO_ERROR;

  try {
    ret=dataTree()->readParameters(file);
  } catch (std::exception ex) {
    ret = TextParserErrorHandler(TP_FILEINPUT_ERROR, file);
  }
   return ret;
}



/** 指定されたファイルに現在保持しているtree構造を、
 *  パラメータファイルの書式で書き出します。
 *  依存関係付き値や、配列型ラベルの引数は、
 *  展開し確定した値を書き出します。
 *
 * @param[in] filename 出力するパラメータファイル名
 * @return エラーコード 
 *
 */

TextParserError TextParser::write(const std::string &file,int order)
{
    TextParserError ret = TP_NO_ERROR;

    try {
      dataTree()->writeParameters(file,order);
      
    } catch (std::exception ex) {
      ret = TextParserErrorHandler(TP_FILEOUTPUT_ERROR, file);
    }
    return ret;
    //    *error = ret;
}



/** 全てのパラメータのラベルパスを取得
 *
 * @param[out] labels 全てのパラメータへのラベルパス
 * @return エラーコード
 */
TextParserError TextParser::getAllLabels(std::vector<std::string>& labels)
{
  //    std::vector<std::string> labels;
  //    MgppError ret = MGPP_NO_ERROR;
  TextParserError ret = TP_NO_ERROR;
  
    try {
      if (!dataTree()->isReady()) {
	ret = TextParserErrorHandler(TP_DATABASE_NOT_READY_ERROR, "");
      } else {
	std::map<unsigned int, std::string>::iterator 
	  il = dataTree()->_leaf_paths.begin();
	while (il != dataTree()->_leaf_paths.end()) {
	  labels.push_back(il->second);
	  il++;
	}
      }
    } catch (std::exception ex) {
        ret = TextParserErrorHandler(TP_GET_PARAMETER_ERROR, "");
    }


    //    *error = ret;
     //    return labels;
    return ret;
}



/** リーフのラベルが示す値のタイプを返す関数
 * 
 *　@param[in] label リーフのラベル
 *  @param[out] error エラーコード
 *  @return labelが指す値のタイプ 型はTextParserValueType
 */
TextParserValueType TextParser::getType(const std::string& label, int *error)
{
    TextParserValueType returner;
    TextParserError ret = TP_NO_ERROR;

    try {
      if (!dataTree()->isReady()) {
	ret = TextParserErrorHandler(TP_DATABASE_NOT_READY_ERROR, "");
        } else {
            TextParserValue *value;
            ret = dataTree()->getLeafValue(label, &value);
            if (ret != TP_NO_ERROR) {
                if (ret == TP_MISSING_PATH_ELEMENT_ERROR)
		  TextParserErrorHandler(ret, label);
            } else {
                    returner = value->_value_type;
                if (value->_value_type == TP_UNDEFFINED_VALUE) {
                    //*type = 0;
                } else if (value->_value_type == TP_NUMERIC_VALUE) {
                    //*type = 1;
                } else if (value->_value_type == TP_STRING_VALUE) {
                    //*type = 2;
                } else if (value->_value_type == TP_VECTOR_UNDEFFINED) {
                        //*type = 4;
                } else if (value->_value_type == TP_VECTOR_NUMERIC) {
                        //*type = 5;
                } else if (value->_value_type == TP_VECTOR_STRING) {
                        //*type = 6;
                } else {
                    ret = TextParserErrorHandler(TP_ILLEGAL_VALUE_TYPE_ERROR, "");
                }
            }
        }


    } catch (std::exception ex) {
      ret = TextParserErrorHandler(TP_GET_PARAMETER_ERROR, "");
    }
    
    *error = ret;
    
    return returner;
}




/** パラメータの値を取得
 *
 * @param[in] label リーフのラベルパス
 * @param[out] value_string リーフの値
 * @param[out] ierror エラーコード
 * 
 *
 */

void TextParser::getValue(const std::string& label,std::string& value_string,int* ierror)
{
  //  std::string value_string;
  TextParserError ret = TP_NO_ERROR;
  
  try {
    if (!dataTree()->isReady()) {
      ret = TextParserErrorHandler(TP_DATABASE_NOT_READY_ERROR, "");
    } else {
      TextParserValue *value;
      ret = dataTree()->getLeafValue(label, &value);
      value_string = value->value();
      if (ret == TP_MISSING_PATH_ELEMENT_ERROR)
	TextParserErrorHandler(ret, label);
    }
    
  } catch (std::exception ex) {
    ret = TextParserErrorHandler(TP_GET_PARAMETER_ERROR, "");
  }

  //  return ret;
    *ierror = ret;
    
  //  return value_string;
}

/** パラメータの値を取得
 *
 * @param[in] label パラメータのラベルパス
 * @param[out] value_string パラメータの値
 * @return エラーコード
 *
 */

TextParserError TextParser::getValue(const std::string& label,std::string& value_string)
{
  //  std::string value_string;
  TextParserError ret = TP_NO_ERROR;
  
  try {
    if (!dataTree()->isReady()) {
      ret = TextParserErrorHandler(TP_DATABASE_NOT_READY_ERROR, "");
    } else {
      TextParserValue *tmpvalue;
      ret = dataTree()->getLeafValue(label, &tmpvalue);
      value_string = tmpvalue->value();
      if (ret == TP_MISSING_PATH_ELEMENT_ERROR)
	TextParserErrorHandler(ret, label);
    }
    
  } catch (std::exception ex) {
    ret = TextParserErrorHandler(TP_GET_PARAMETER_ERROR, "");
  }

  return ret;
    

}

/** 現在のノードの取得
 *
 * @param[out] returner ノードのパスラベル(std::string).
 * @return エラーコード.
 *
 */
TextParserError TextParser::currentNode(std::string& returner)
{
    TextParserError ret = TP_NO_ERROR;
    //    std::string returner = "";

    try {
        //! インスタンスの取得

        if (!dataTree()->isReady()) {
                ret = TextParserErrorHandler(TP_DATABASE_NOT_READY_ERROR, "");
        } else {
            if (dataTree()->_current_element == 0) {
                    // do nothing
                    //            *label = 0;
                    // is this right. Check this later.
                    returner = dataTree()->_label;
            } else {
                    //*label = (char *)dataTree()->_current_element->_label.c_str();
                    returner = dataTree()->_current_element->_label.c_str();
            }
        }
    } catch (std::exception ex) {
        ret = TextParserErrorHandler(TP_GET_PARAMETER_ERROR, "");
    }

    //    *error = ret;
    return ret;
    //    return returner;
}


/** 現在のノード内の全ノードを取得
 *
 * @param[out] node_list 現在のノード内の全ノード std::vector<std::sring>型
 * @return error エラーコード
 */

TextParserError TextParser::getNodes(std::vector<std::string>& node_list,int oswitch){
  
  TextParserError ret = TP_NO_ERROR;
  
  try {
    
    if (!dataTree()->isReady()) {
      ret = TextParserErrorHandler(TP_DATABASE_NOT_READY_ERROR, "");
    } else {
      ret =dataTree()->getCurrentNodeLabels(node_list);
    }
  } catch (std::exception ex) {
    ret = TextParserErrorHandler(TP_GET_PARAMETER_ERROR, "");
  }
  
  if(oswitch!=0){
    std::vector< std::string> output;
  dataTree()->nodeSort(node_list,output,oswitch);
  node_list=output;
}
//  *error = ret;
//    return labels;
return ret;

}


/** リーフのストリングをcharに変換
 * @param[in] value リーフの値（文字列).
 * @param[out] error エラーコード.
 * @return char に変換された値.
 */
char TextParser::convertChar(const std::string& value, int *error)
{
    *error = 0;

    char returner;
    int int_recieve;
    try {
      std::string val = CorrectValueString(value);
        std::stringstream ss;
        ss << val;
	//        ss >> (int_recieve);
	//	returner=int_recieve;
	ss >> (returner);
    } catch (std::exception ex) {
        *error = TextParserErrorHandler(TP_MEMORY_ALLOCATION_ERROR, value);
    }
  
    return returner;
}


/** リーフのストリングをshortに変換
 * @param[in] value リーフの値（文字列).
 * @param[out] error エラーコード.
 * @return short に変換された値.
 */
short TextParser::convertShort(const std::string& value, int *error)
{
    *error = 0;

    short returner;
    try {
        std::string val = CorrectValueString(value);
        std::stringstream ss;
        ss << val;
        ss >> (returner);
    } catch (std::exception ex) {
        *error = TextParserErrorHandler(TP_MEMORY_ALLOCATION_ERROR, value);
    }
  
    return returner;
}



/** リーフのストリングをintに変換
 * @param[in] value リーフの値（文字列).
 * @param[out] error エラーコード.
 * @return int に変換された値.
 */
int TextParser::convertInt(const std::string& value, int *error)
{
    *error = 0;

    int returner;
    try {
        std::string val = CorrectValueString(value);
        std::stringstream ss;
        ss << val;
        ss >> (returner);
    } catch (std::exception ex) {
        *error = TextParserErrorHandler(TP_MEMORY_ALLOCATION_ERROR, value);
    }
  
    return returner;
}

/** リーフのストリングをlongに変換
 * @param[in] value リーフの値（文字列).
 * @param[out] error エラーコード.
 * @return long に変換された値.
 */
long TextParser::convertLong(const std::string& value, int *error)
{
    *error = 0;

    long returner;
    try {
        std::string val = CorrectValueString(value);
        std::stringstream ss;
        ss << val;
        ss >> (returner);
    } catch (std::exception ex) {
        *error = TextParserErrorHandler(TP_MEMORY_ALLOCATION_ERROR, value);
    }
  
    return returner;
}


/** リーフのストリングをlong longに変換
 * @param[in] value リーフの値（文字列).
 * @param[out] error エラーコード.
 * @return long longに変換された値.
 */

long long TextParser::convertLongLong(const std::string& value, int *error)
{
    *error = 0;

    long long returner;
    try {
        std::string val = CorrectValueString(value);
        std::stringstream ss;
        ss << val;
        ss >> (returner);
    } catch (std::exception ex) {
        *error = TextParserErrorHandler(TP_MEMORY_ALLOCATION_ERROR, value);
    }
  
    return returner;
}



/** リーフのストリングをfloatに変換
 * @param[in] value リーフの値（文字列).
 * @param[out] error エラーコード.
 * @return float に変換された値.
 */

float TextParser::convertFloat(const std::string& value, int *error)
{
    *error=0;

    float returner;
    try {
        std::string val = CorrectValueString(value);
        std::stringstream ss;
        ss << val;
        ss >> (returner);
    } catch (std::exception ex) {
        *error = TextParserErrorHandler(TP_MEMORY_ALLOCATION_ERROR, value);
    }
  
    return returner;
}


/** リーフのストリングをdoubleに変換
 * @param[in] value リーフの値（文字列).
 * @param[out] error エラーコード.
 * @return double に変換された値.
 */
double TextParser::convertDouble(const std::string& value, int *error)
{
    *error=0;

    double returner;
    try {
        std::string val = CorrectValueString(value);
        std::stringstream ss;
        ss << val;
        ss >> (returner);
    } catch (std::exception ex) {
        *error = TextParserErrorHandler(TP_MEMORY_ALLOCATION_ERROR, value);
    }
  
    return returner;
}
/** リーフのストリングをboolに変換
 *
 * true : "true"(大文字小文字は関係無し)  
 * false : "false"（大文字小文字は関係無し）  
 * 上記を満たさない文字列の場合は、falseを返す。
 *
 * @param[in] value リーフの値（文字列).
 * @param[out] error エラーコード.
 * @return bool に変換された値.
 */
bool TextParser::convertBool(const std::string& value, int *error)
{
   *error=0;
    bool returner=false;
#ifdef MYDEBUG3
    std::cout << "ConvertBool" << " value |" << value << "|"<<std::endl;
#endif 
    std::string tstring="true";
    if(TextParserStringCompare(value,tstring)){
      returner=true;
      return returner;
    }
    tstring="1";
    if(TextParserStringCompare(value,tstring)){
      returner=true;
      return returner;
    }
    tstring="false";
    if(TextParserStringCompare(value,tstring)){
      return returner;
    }
    tstring="0";
    if(TextParserStringCompare(value,tstring)){
      return returner;
    }
    return returner;
}


/** ベクトル値を要素ごとに文字列に分割する。
 * @param[in] vector_value ベクトル値.
 * @param[out] velem 分離された値の文字列のベクトル
 * @return エラーコード.
 */
TextParserError TextParser::splitVector(const std::string& vector_value,
					std::vector<std::string>& velem )
{
    TextParserError ret = TP_NO_ERROR;

    try {
        if (!dataTree()->isReady()) {
                ret = TextParserErrorHandler(TP_DATABASE_NOT_READY_ERROR, "");
        } else {
            ret = dataTree()->splitVectorValue(vector_value, velem);
        }
    } catch (std::exception ex) {
        ret = TextParserErrorHandler(TP_GET_PARAMETER_ERROR, "");
    }
    //    std::cout << velem.size() <<std::endl;
    return ret;
}




/** パラメータの削除
 *
 * @return エラーコード
 *
 */
TextParserError TextParser::remove()
{
    TextParserError ret = TP_NO_ERROR;
    try {
        ret = dataTree()->removeElement();
    } catch (std::exception ex) {
        ret = TextParserErrorHandler(TP_REMOVE_ELEMENT_ERROR, "");
    }
    return ret;
}
/** ノードの移動
 *
 * @param[in] label ノードのパスラベル.
 * @return エラーコード.
 *
 */
TextParserError TextParser::changeNode(const std::string& label)
{
    TextParserError ret = TP_NO_ERROR;

    try {
        //! インスタンスの取得

        if (!dataTree()->isReady()) {
                ret = TextParserErrorHandler(TP_DATABASE_NOT_READY_ERROR, "");
        } else {
            ret = dataTree()->changeNode(label);
        }
    } catch (std::exception ex) {
        ret = TextParserErrorHandler(TP_GET_PARAMETER_ERROR, "");
    }
    return ret;

}


/** 現在のノード内の全ラベルを取得
 *
 * @param[out] labels 現在のディレクトリ内の全ラベル std::vector<std::sring>型
 * @return エラーコード
 */
TextParserError TextParser::getLabels(std::vector<std::string>& labels,int oswitch)
{
  //    std::vector<std::string> labels;
    TextParserError ret = TP_NO_ERROR;
  
    try {
        //! インスタンスの取得
        if (!dataTree()->isReady()) {
            ret = TextParserErrorHandler(TP_DATABASE_NOT_READY_ERROR, "");
        } else {
            ret = dataTree()->getCurrentLeafLabels(labels);
        }
    } catch (std::exception ex) {
        ret = TextParserErrorHandler(TP_GET_PARAMETER_ERROR, "");
    }

    if(oswitch==0){
      std::vector<std::string> output;
      //    dataTree()->labelSort(labels,output,0);
      dataTree()->labelSort(labels,output,oswitch);
      labels=output;
    }
    return ret;
    //    *error = ret;
    
    //    return labels;
}


 // global functions for C API 

 /** 入力ファイルを読み込み、各パラメータをtree構造のデータとして格納する。C用API
  * 
  * @param[in] cfile 入力ファイル名
  * @return エラーコード、TextParserErrorによる。intで取得。
  */
int tp_read(char* cfile){
  std::string filestring(cfile);
  TextParser* tp_ptr=TextParser::get_instance();
  int error = tp_ptr->read(filestring);
  return  error;
}

 /** ファイルにデータ構造を格納する。C用API
  * 
  * @param[in] cfile 出力ファイル名
  * @return エラーコード、TextParserErrorによる。intで取得。
  */
int tp_write(char* cfile){
  std::string filestring(cfile);
  TextParser* tp_ptr=TextParser::get_instance();
  int error = tp_ptr->write(filestring);
  return  error;
}

 /** ファイルにデータ構造を破棄する。C用API
  * 
  * @return エラーコード、TextParserErrorによる。intで取得。
  */

int tp_remove(){
  TextParser* tp_ptr=TextParser::get_instance();
  int error = tp_ptr->remove();
  return  error;
}


 /** リーフの総数を取得する。C用API
  *
  * @param[out] Nleaves リーフの総数
  * @return エラーコード、TextParserErrorによる。intで取得。
  */

int tp_getNumberOfLeaves(unsigned int* Nleaves){
  TextParser* tp_ptr=TextParser::get_instance();
  int error = 0;
  if( !tp_ptr->dataTree()->isReady() ){
    error=TextParserErrorHandler(TP_DATABASE_NOT_READY_ERROR, "");
  } else {
    int size= tp_ptr->dataTree()->_leaf_paths.size();
    *Nleaves = size;
  }
  return error;
}


 /** ラベルの取得。C用API
  *
  * @param[in] ilabel 取得するラベルの番号
  * @param[out] label ラベル
  * @return エラーコード、TextParserErrorによる。intで取得。
  */

int tp_getLabel(int ilabel,char* label){

  TextParser* tp_ptr=TextParser::get_instance();
  int error = 0;
  if( !tp_ptr->dataTree()->isReady() ){
    error=TextParserErrorHandler(TP_DATABASE_NOT_READY_ERROR, "");
  } else {
    if (ilabel >= tp_ptr->dataTree()->_leaf_paths.size() ) {
      error = TextParserErrorHandler(TP_ID_OVER_ELEMENT_NUMBER_ERROR, "");
    } else {
      std::map<unsigned int, std::string>::iterator
	il = tp_ptr->dataTree()->_leaf_paths.find(ilabel);
      strcpy(label,il->second.c_str());
    }
  }
  return error;
}
 /** 値の取得。C用API
  *
  * @param[in] label typeを取得するラベル
  * @param[out] value 値 
  * @return エラーコード、TextParserErrorによる。intで取得。
  */

int tp_getValue(char* label,char* value){

  TextParser* tp_ptr=TextParser::get_instance();
  int error = 0;
  std::string slabel(label);
  std::string svalue;
  error=tp_ptr->getValue(slabel,svalue);
  strcpy(value,svalue.c_str());
  return error;
}

 /** タイプの取得。C用API
  *
  * @param[in] label typeを取得するラベル
  * @param[out] type 値のタイプ。TextParserValuetype による。int で取得。
  * @return エラーコード、TextParserErrorによる。intで取得。
  */

int tp_getType(char* label,int* type){
  *type=-10;
  TextParser* tp_ptr=TextParser::get_instance();
  int error = 0;
  std::string slabel(label);
  *type=tp_ptr->getType(slabel,&error);
  return error;
}

int typeToInt(TextParserValueType type){
  int ret=-1;
  
  switch(type){
    
  TP_UNDEFFINED_VALUE:
    ret = 0; break;  //!< 不定
  TP_NUMERIC_VALUE:
    ret = 1; break;     //!< 数値
  TP_STRING_VALUE:
    ret = 2; break;      //!< 文字列
  TP_DEPENDENCE_VALUE:
    ret = 3; break;  //!< 依存関係付き値
  TP_VECTOR_UNDEFFINED:
    ret = 4; break; //!< ベクトル型不定
  TP_VECTOR_NUMERIC:
    ret = 5; break;    //!< ベクトル型数値
  TP_VECTOR_STRING:
    ret = 6; break;     //!< ベクトル型文字列
  default:
    break;
  }
}

/** パラメータの値を文字列から char 型へ変換します。C用API.
 *
 * @param[in] value 値の文字列
 * @param[out] error エラーコード
 * @return char型への変換値
 */

char tp_convertChar(char* value,int *error){
  TextParser* tp_ptr=TextParser::get_instance();
  std::string svalue(value);

  int int_recieve=tp_ptr->convertInt(svalue,error);
  char returner = int_recieve;
  std::cout << __FUNCTION__ << svalue  << " int "<< int_recieve<< " char "<< (int)returner <<std::endl; 
  return returner;

}
/** パラメータの値を文字列から short 型へ変換します。C用API.
 *
 * @param[in] value 値の文字列
 * @param[out] error エラーコード
 * @return short型への変換値
 */

short tp_convertShort(char* value,int *error){
  TextParser* tp_ptr=TextParser::get_instance();
  std::string svalue(value);
  return tp_ptr->convertShort(svalue,error);
}
/** パラメータの値を文字列から int 型へ変換します。C用API.
 *
 * @param[in] value 値の文字列
 * @param[out] error エラーコード
 * @return short型への変換値
 */

int tp_convertInt(char* value,int *error){
  TextParser* tp_ptr=TextParser::get_instance();
  std::string svalue(value);
  return tp_ptr->convertInt(svalue,error);
}

/** パラメータの値を文字列から long 型へ変換します。C用API.
 *
 * @param[in] value 値の文字列
 * @param[out] error エラーコード
 * @return short型への変換値
 */

long tp_convertLong(char* value,int *error){
  TextParser* tp_ptr=TextParser::get_instance();
  std::string svalue(value);
  return tp_ptr->convertLong(svalue,error);
}

/** パラメータの値を文字列から long long 型へ変換します。C用API.
 *
 * @param[in] value 値の文字列
 * @param[out] error エラーコード
 * @return short型への変換値
 */

long long tp_convertLongLong(char* value,int *error){
  TextParser* tp_ptr=TextParser::get_instance();
  std::string svalue(value);
  return tp_ptr->convertLongLong(svalue,error);
}

/** パラメータの値を文字列から float 型へ変換します。C用API.
 *
 * @param[in] value 値の文字列
 * @param[out] error エラーコード
 * @return short型への変換値
 */
float tp_convertFloat(char* value,int *error){
  TextParser* tp_ptr=TextParser::get_instance();
  std::string svalue(value);
  return tp_ptr->convertFloat(svalue,error);
}

/** パラメータの値を文字列から double 型へ変換します。C用API.
 *
 * @param[in] value 値の文字列
 * @param[out] error エラーコード
 * @return short型への変換値
 */
double tp_convertDouble(char* value,int *error){
  TextParser* tp_ptr=TextParser::get_instance();
  std::string svalue(value);
  return tp_ptr->convertDouble(svalue,error);
}

/** パラメータの値を文字列から bool型をint 型の形式で変換します。C用API.
 *
 * true --> 1 
 * false --> 0
 *
 * @param[in] value 値の文字列
 * @param[out] error エラーコード
 * @return int型への変換値
 */
int tp_convertBool(char* value,int *error){
  TextParser* tp_ptr=TextParser::get_instance();
  std::string svalue(value);
  int ret =0;
  if(tp_ptr->convertBool(svalue,error)) ret=1;
  return ret ;
}

/** ベクトル型の値の要素数を取得する。C用API.
 *
 * @param[in] vector_value ベクトル型の値（文字列）
 * @param[out] nvec 要素数 
 * @return エラーコード、TextParserErrorによる。intで取得。
 *
 */
int tp_getNumberOfElements(char* vector_value,unsigned int* nvec ){

  TextParser* tp_ptr=TextParser::get_instance();

  std::string svalue(vector_value);
  std::vector<std::string> tmp_vector_string;
  int error=tp_ptr->splitVector(svalue,tmp_vector_string);
  int size=tmp_vector_string.size();
  *nvec =  size;
  return error;
}

/** 指定したindexのベクトル型の値の要素を取得する. C用API.
 *
 * @param[in] vector_value ベクトル型の値（文字列）
 * @param[in] ivec 取得する要素のインデックス
 * @param[out] velem ベクトル値の要素（文字列）
 * @return エラーコード、TextParserErrorによる。intで取得。
 *
 */
int tp_getIthElement(char* vector_value,unsigned int ivec,char* velem){
  TextParser* tp_ptr=TextParser::get_instance();
  std::string svalue(vector_value);
  std::vector<std::string> tmp_vector_string;
  int error=tp_ptr->splitVector(svalue,tmp_vector_string);
  strcpy(velem,tmp_vector_string[ivec].c_str());
  return error;
}

 /** カレントノードを取得します。C用API
  *
  * @param[out] label ノードのラベルが返ります。
  * @return エラーコード、TextParserErrorによる。intで取得。
  */
int tp_currentNode(char* label){

  TextParser* tp_ptr=TextParser::get_instance();
  int error = 0;
  std::string slabel(label);
  error=tp_ptr->currentNode(slabel);
  strcpy(label,slabel.c_str());
  return error;
}
 /** labelで指定されたノードをカレントノードに設定します。C用API
  *
  * @param[in] label 移動するノード
  * @return エラーコード、TextParserErrorによる。intで取得。
  */
int tp_changeNode(char* label){
  TextParser* tp_ptr=TextParser::get_instance();
  int error = 0;
  std::string slabel(label);
  error=tp_ptr->changeNode(slabel);
  return error;
}

 /** カレントノードの子ノードの数を取得します。C用API
  *
  * @param[out] nnodes 子ノードの数
  * @return エラーコード、TextParserErrorによる。intで取得。
  */
int tp_getNumberOfCNodes(int* nnodes){
  TextParser* tp_ptr=TextParser::get_instance();
  int error = 0;
  std::vector<std::string> node_labels;
  error=tp_ptr->getNodes(node_labels);
  *nnodes=node_labels.size();
  return error;
}

 /** カレントノードのリーフの数を取得します。C用API
  *
  * @param[out] nleaves リーフの数
  * @return エラーコード、TextParserErrorによる。intで取得。
  */
int tp_getNumberOfCLeaves(int* nleaves){
  TextParser* tp_ptr=TextParser::get_instance();
  int error = 0;
  std::vector<std::string> leaf_labels;
  error=tp_ptr->getLabels(leaf_labels);
  *nleaves=leaf_labels.size();
  return error;
}

 /** カレントノードにあるインデックスで指定したノードのラベルを取得します。C用API
  * 
  * @param[in] inode ノードのインデックス
  * @param[out] node ノードのラベル
  * @return エラーコード、TextParserErrorによる。intで取得。
  */
int tp_getIthNode(int inode,char* node){

  TextParser* tp_ptr=TextParser::get_instance();
  int error = 0;
  std::vector<std::string> node_labels;
  error=tp_ptr->getNodes(node_labels);
  strcpy(node,node_labels[inode].c_str());
  return error;
}
 /** カレントノードにあるインデックスで指定したノードのラベルを取得します。C用API
  * 
  * @param[in] inode ノードのインデックス
  * @param[out] node ノードのラベル
  * @param[in] order ラベルの出力順　0:getIthNode 同様　1:配列型ラベルのインデックス順　2:ラベルの出現順.
  * @return エラーコード、TextParserErrorによる。intで取得。
  */
int tp_getIthNodeOrder(int inode,char* node,int order){


  TextParser* tp_ptr=TextParser::get_instance();
  int error = 0;
  std::vector<std::string> node_labels;
  error=tp_ptr->getNodes(node_labels,order);
  strcpy(node,node_labels[inode].c_str());
  return error;
}

 /** カレントノードにあるインデックスで指定したリーフを取得します。C用API
  * 
  * @param[in] ileaf リーフのインデックス
  * @param[out] leaf リーフのラベル
  * @return エラーコード、TextParserErrorによる。intで取得。
  */
int tp_getIthLeaf(int ileaf,char* leaf){
  TextParser* tp_ptr=TextParser::get_instance();
  int error = 0;
  std::vector<std::string> leaf_labels;
  error=tp_ptr->getLabels(leaf_labels);
  strcpy(leaf,leaf_labels[ileaf].c_str());
  return error;
}

 /** カレントノードにあるインデックスで指定したリーフを取得します。C用API
  * 
  * @param[in] ileaf リーフのインデックス
  * @param[out] leaf リーフのラベル
  * @param[in] order ラベルの出力順　0:getIthNode 同様　1:配列型ラベルのインデックス順　2:ラベルの出現順.
  * @return エラーコード、TextParserErrorによる。intで取得。
  */
int tp_getIthLeafOrder(int ileaf,char* leaf,int order){
  TextParser* tp_ptr=TextParser::get_instance();
  int error = 0;
  std::vector<std::string> leaf_labels;
  error=tp_ptr->getLabels(leaf_labels,order);
  strcpy(leaf,leaf_labels[ileaf].c_str());
  return error;
}

/** ファイルからパラメータを読み込む.　Fortran用API
 *
 * @param[in] file ファイル名
 * @param[in] length ファイル名文字列の長さ
 * @return エラーコード、TextParserErrorによる。intで取得。
 */

int tp_read_fort_(char* file,int* length){
  return tp_read(file);
}

/** パラメータをファイルに書き出す.　Fortran用API
 *
 * @param[in] file ファイル名
 * @param[in] length ファイル名の長さ 
 * @return エラーコード、TextParserErrorによる。intで取得。
 */
int tp_write_fort_(char* file,int* length){
  return tp_write(file);
}
/** パラメータデータをメモリから消去　Fortran用API
 *
 * @return エラーコード、TextParserErrorによる。intで取得。
 */
int tp_remove_fort_(){
  return tp_remove();
}

/** リーフの総数を取得する.　Fortran用API
 *
 * @param[out] nleaves リーフの総数
 * @return エラーコード、TextParserErrorによる。intで取得。
 */
int tp_get_number_of_leaves_fort_(int* nleaves ){
  unsigned int i;
  //  std::cout << "tp_get_number_of_leaves_fort_"<< std::endl;
  int error=tp_getNumberOfLeaves(&i);
  //  std::cout << "tp_get_number_of_leaves_fort_"<< std::endl;
  *nleaves=i;

  return error;
}

/** リーフのラベルを取得する.　Fortran用API
 *
 * @param[in] ileaf インデックス
 * @param[out] label リーフのラベル
 * @param[out] length 文字列長
 * @return エラーコード、TextParserErrorによる。intで取得。
 */
int tp_get_label_fort_(int* ileaf ,char* label,int* length){

  unsigned int i;
  i=*ileaf;
  int j;
  char tmp_label[TP_FORTRAN_BUFFER_SIZE];
  int error = tp_getLabel(i,tmp_label);
  int llength = strlen(tmp_label);
  for (j=0;j<llength;j++){
    label[j]=tmp_label[j];
  }

    //return error;

  return error;
}

/** リーフの値を取得する.　Fortran用API
 *
 * @param[in] label リーフのラベル
 * @param[in] label_length 文字列長
 * @param[out] value 値
 * @param[out] value_length 文字列長
 * @return エラーコード、TextParserErrorによる。intで取得。
 */
int tp_get_value_fort_(char* label,char* value,int* label_length,int* value_length){

  char tmp_value[TP_FORTRAN_BUFFER_SIZE];
  int error=tp_getValue(label,tmp_value);
  int vlen=strlen(tmp_value);
  int j;

  //  std::cout <<error<<"length "<<vlen<< "|"<<tmp_value<<"|"<<std::endl;

  for (j=0;j<vlen;j++){
    value[j]=tmp_value[j];
  }

  return error;
}
/** 値のタイプを取得する.　Fortran用API
 *
 * @param[in] label リーフのラベル
 * @param[in] label_length 文字列長
 * @param[out] type
 * @return エラーコード、TextParserErrorによる。intで取得。
 */
int tp_get_type_fort_(char* label,int* type,int* label_length){

  return tp_getType(label,type);
}



/** パラメータの値を文字列から char 型へ変換します。fortran 用API.
 *
 * @param[in] value 値の文字列
 * @param[out] error エラーコード
 * @param[out] value_length 文字列長
 * @return char型への変換値
 */

char tp_convert_char_fort_(char* value,int *error,int* value_length){
  
  return tp_convertChar(value,error);
  
}

/** パラメータの値を文字列から short 型へ変換します。fortran 用API.
 *
 * @param[in] value 値の文字列
 * @param[out] error エラーコード
 * @param[out] value_length 文字列長
 * @return short型への変換値
 */

short tp_convert_short_fort_(char* value,
			     int *error,int* value_length){
  return tp_convertShort(value,error);
}

/** パラメータの値を文字列から int 型へ変換します。fortran 用API.
 *
 * @param[in] value 値の文字列
 * @param[out] error エラーコード
 * @param[out] value_length 文字列長
 * @return int型への変換値
 */

int tp_convert_int_fort_(char* value,
			 int *error,int* value_length){
  return tp_convertInt(value,error);
  
} 

/** パラメータの値を文字列から long 型へ変換します。fortran 用API.
 *
 * @param[in] value 値の文字列
 * @param[out] error エラーコード
 * @param[out] value_length 文字列長
 * @return long型への変換値
 */

long tp_convert_long_fort_(char* value,
			  int *error,int* value_length){
  return tp_convertLong(value,error);
} 

/** パラメータの値を文字列から long long 型へ変換します。fortran 用API.
 *
 * @param[in] value 値の文字列
 * @param[out] error エラーコード
 * @param[out] value_length 文字列長
 * @return long long型への変換値
 */

long long tp_convert_long_long_fort_(char* value,
				     int *error,int* value_length){
  return tp_convertLongLong(value,error);
}
/** パラメータの値を文字列から float 型へ変換します。fortran 用API.
 *
 * @param[in] value 値の文字列
 * @param[out] error エラーコード
 * @param[out] value_length 文字列長
 * @return float型への変換値
 */

float tp_convert_float_fort_(char* value,
			     int *error,int* value_length){
  return tp_convertFloat(value,error);
}
/** パラメータの値を文字列から double 型へ変換します。fortran 用API.
 *
 * @param[in] value 値の文字列
 * @param[out] error エラーコード
 * @param[out] value_length 文字列長
 * @return double型への変換値
 */

double tp_convert_double_fort_(char* value,
			       int *error,int* value_length){
  return tp_convertDouble(value,error);
}

/** パラメータの値を文字列から int(boolian) 型へ変換します。fortran 用API.
 *
 * @param[in] value 値の文字列
 * @param[out] error エラーコード
 * @param[out] value_length 文字列長
 * @return int 型への変換値
 */
int tp_convert_logical_fort_(char* value,
			     int *error,int* value_length){
  return tp_convertBool(value,error);
}


/** ベクトル型の値の要素数を取得する。Fortrean用API.
 *
 * @param[in] vector_value ベクトル型の値（文字列）
 * @param[out] nvec 要素数 
 * @param[in] vector_length 文字列長
 * @return エラーコード、TextParserErrorによる。intで取得。
 *
 */
int tp_get_number_of_elements_fort_(char* vector_value,unsigned int* nvec,int* vector_length){
  return tp_getNumberOfElements(vector_value,nvec);
} 


/** 指定したindexのベクトル型の値の要素を取得する. Fortran用API.
 *
 * @param[in] vector_value ベクトル型の値（文字列）
 * @param[in] ivec 取得する要素のインデックス
 * @param[out] velem ベクトル値の要素（文字列）
 * @param[in] vector_value_length ベクトル値の文字列長
 * @param[out] velem_length ベクトル値の要素の文字列長
 * @return エラーコード、TextParserErrorによる。intで取得。
 *
 */
int tp_get_ith_element_fort_(char* vector_value,
			     int* ivec,
			     char* velem,
			     int* vector_value_length,
			     int* velem_length){

  unsigned int id_vec;
  char tmp_velem[TP_FORTRAN_BUFFER_SIZE];
  int i;
  int error;
  id_vec=*ivec;
  error = tp_getIthElement(vector_value,id_vec,tmp_velem);
  for(i=0;i<strlen(tmp_velem);++i){
    velem[i]=tmp_velem[i];
  }
  

  return error;

}


  // 
  //カレントノードの取得、ノードの移動 fortran　用 API


int tp_current_node_fort_(char* label,int* label_length){
  return tp_currentNode(label);
}
int tp_change_node_fort_(char* label,int* label_length){
  return tp_changeNode(label);
}
int tp_get_number_of_cnodes_fort_(int* nnodes){
  return tp_getNumberOfCNodes(nnodes);
}
int tp_get_number_of_cleaves_fort_(int* nleaves){
  return tp_getNumberOfCLeaves(nleaves);
}
int tp_get_ith_node_fort_(int* inode,char* node,int* node_length){
  //  std::cout << "tp_get_ith_node_fort_" << *inode << std::endl;
  char tmp_node[TP_FORTRAN_BUFFER_SIZE];
  int i;
  int error=tp_getIthNode(*inode,tmp_node);
  for(i=0;i<strlen(tmp_node);++i){
    node[i]=tmp_node[i];
  }

  return error;
}

int tp_get_ith_node_order_fort_(int* inode,char* node,int* order,int* node_length){
  //  std::cout << "tp_get_ith_node_fort_" << *inode << std::endl;
  char tmp_node[TP_FORTRAN_BUFFER_SIZE];
  int i;
  int error=tp_getIthNodeOrder(*inode,tmp_node,*order);
  for(i=0;i<strlen(tmp_node);++i){
    node[i]=tmp_node[i];
  }

  return error;
}

int tp_get_ith_leaf_fort_(int* ileaf,char* leaf,int* leaf_length){

  char tmp_leaf[TP_FORTRAN_BUFFER_SIZE];
  int i;
  int error=tp_getIthLeaf(*ileaf,tmp_leaf);
  for(i=0;i<strlen(tmp_leaf);++i){
    leaf[i]=tmp_leaf[i];
  }

  return error;
}
int tp_get_ith_leaf_order_fort_(int* ileaf,char* leaf,int* order,int* leaf_length){

  char tmp_leaf[TP_FORTRAN_BUFFER_SIZE];
  int i;
  int error=tp_getIthLeafOrder(*ileaf,tmp_leaf,*order);
  for(i=0;i<strlen(tmp_leaf);++i){
    leaf[i]=tmp_leaf[i];
  }

  return error;
}
