// -*- mode: c++ -*-
/****************************************************************************
**
** Copyright (C) 2012 Tokyo University.
**
****************************************************************************/

/** @file TextParser.h
 * ここには TextParser クラス　と C用 Fortran 用 API 関数が定義されています。
 *
 */

#ifndef __TEXTPARSER_H__
#define __TEXTPARSER_H__


//# define TP_FORTRAN_BUFFER_SIZE 128
# define TP_FORTRAN_BUFFER_SIZE 1024
//#include <varargs.h>
//#include <stdarg.h>

#include "TextParserCommon.h"
#ifdef __cplusplus
#include <iostream>
#include <string>
#include <sstream>
#include "TextParserElement.h"
#include "TextParserTree.h"
#endif // __cplusplus

#ifdef __cplusplus


/** 
 * @class TextParser TextParser.h
 * @brief TextParserライブラリのAPIクラス TextParser が定義されています.
 * パラメータのファイルからの読み込み、
 * パラメータのデータ構造のファイルへの書き出し、
 * パラメータのデータ構造へのアクセス、
 * パラメータ（文字列）の変換、
 * 及びそのヘルパー関数が定義されています。
 * 現時点ではSingletonパターンで実装されています。
 */


class TextParser{

 public:
  static TextParser* get_instance(); //!< TextParser のインスタンスを取得します。
  TextParserError read(const std::string& file); //!< file からパラメータを読み込みます。
  TextParserError write(const std::string& file,int order=0); //!< file にパラメータのデータ構造全体を書き込みます。
  TextParserError remove();//!< 格納しているパラメータのデータを破棄します。

  TextParserError getAllLabels(std::vector<std::string>& labels);//!< 格納されているパラメータへの(リーフへの)全てのラベルを取得します。
  TextParserValueType getType(const std::string& label, int *error);//!<ラベルで示されるパラメータの値のタイプを返します。（TextParserError参照） 
  TextParserError getValue(const std::string& label,std::string& value);//!< ラベルで示されるパラメータの値を取得します。
  void getValue(const std::string& label,std::string& value, int*ierror);//!< ラベルで示されるパラメータの値を取得します。
  TextParserError currentNode(std::string& returner); //!<現在のカレントノードを取得します。
  //  std::string CurrentNode(int *error);
  TextParserError getNodes(std::vector<std::string>& node_list,int oswitch=0);//!<カレントノード内の子ノードのラベル（相対パス）のリストを取得します。
  TextParserError getLabels(std::vector<std::string>& labels,int oswitch=0);//!<カレントノード内のリーフのラベル（相対パス）のリストを取得します。

  TextParserError changeNode(const std::string& label);//!< ラベルで指定したノードをカレントノードにします。

  char convertChar(const std::string& value, int *error); //!< パラメータの値を文字列から char 型へ変換します。 
  short convertShort(const std::string& value, int *error);//!< パラメータの値を文字列から short 型へ変換します。 
  int convertInt(const std::string& value, int *error); //!< パラメータの値を文字列から int 型へ変換します。 
  long convertLong(const std::string& value, int *error); //!< パラメータの値を文字列から long 型へ変換します。 
  long long convertLongLong(const std::string& value, int *error);//!< パラメータの値を文字列から long long 型へ変換します。 
  float convertFloat(const std::string& value, int *error);//!< パラメータの値を文字列から float 型へ変換します。 
  double convertDouble(const std::string& value, int *error);//!< パラメータの値を文字列から double 型へ変換します。 
  bool convertBool(const std::string& value, int *error);//!< パラメータの値を文字列から bool 型へ変換します。 

  TextParserError splitVector(const std::string& vector_value,
			      std::vector<std::string>& velem );//!<ベクトル値を分割する。

 TextParserTree* dataTree() const {return _data_tree;} //!< パラメータデータ構造へのポインタ。
 
protected:
  
private:
  //  TextParser(); //!< デフォルトコンストラクタ 利用禁止
  TextParser(const TextParser& rhs){} //!< コピーコンストラクタ 利用禁止
  TextParser& operator=(const TextParser& rhs){} //!< 代入演算子 利用禁止
  TextParser(); //!< デフォルトコンストラクタ 利用禁止
  ~TextParser(){  remove();  } //!< デストラクタ。パラメータデータ構造を消去します。

  TextParserTree* _data_tree; //!< パラメータデータ構造へのポインタ。
};


#endif /* __cplusplus */


// C / Fortran 用　API 関数 

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
  // for c_interfaces
  int tp_read(char* file); //!< file からパラメータを読み込みます.  C用API.
  int tp_write(char* file); //!< file にパラメータのデータ構造全体を書き込みます. C用API.
  int tp_remove();//!< 格納しているパラメータのデータを破棄します。 C用API.
  // full path access
  int tp_getNumberOfLeaves(unsigned int* Nleaves);//!< リーフの総数を取得します。C用API.
  int tp_getLabel(int ilabel,char* label);//!< インデックスで指定したリーフのラベルを取得します。C用API.

  // value and value type
  int tp_getType(char* label,int* type);//!< ラベルで指定したリーフのタイプを取得します。C用API.
  int tp_getValue(char* label,char* value);//!< ラベルで指定したリーフの値を取得します。C用API.

  //値の変換
  char tp_convertChar(char* value,
		      int *error); //!< パラメータの値を文字列から char 型へ変換します。C用API.
  short tp_convertShort(char* value,
			int *error);//!< パラメータの値を文字列から short 型へ変換します。 C用API.
  int tp_convertInt(char* value,
		    int *error); //!< パラメータの値を文字列から int 型へ変換します。C用API. 
  long tp_convertLong(char* value,
		      int *error); //!< パラメータの値を文字列から long 型へ変換します。 C用API.
  long long tp_convertLongLong(char* value,
			       int *error);//!< パラメータの値を文字列から long long 型へ変換します。C用API. 
  float tp_convertFloat(char* value,
			int *error);//!< パラメータの値を文字列から float 型へ変換します。 C用API.
  double tp_convertDouble(char* value,
			  int *error);//!< パラメータの値を文字列から double 型へ変換します。 C用API.
  int tp_convertBool(char* value,
		     int *error);//!< パラメータの値を文字列から int 型 でbool型を表現したものへ変換します。 C用API.

  //ベクトル値の操作
  int tp_getNumberOfElements(char* vector_value,unsigned int* nvec ); //!<ベクトル値の要素数　C用API.
  int tp_getIthElement(char* vector_value,unsigned int ivec,char* velem );//!<ベクトル値のivec番目の要素を取得する。C用API.

  //カレントノードの取得、ノードの移動
  int tp_currentNode(char*); //!<現在のカレントノードを取得します。C用API.
  int tp_changeNode(char*);//!< ラベルで指定したノードをカレントノードにします。C用API.
  int tp_getNumberOfCNodes(int* nnodes);//!< カレントノードにある子ノードの数を取得します。C用API.
  int tp_getNumberOfCLeaves(int* nleaves);//!< カレントノードにあるリーフの数を取得します。C用API.
  int tp_getIthNode(int inode,char* node);//!< カレントノードにあるinode番目のノードのラベルを取得します。C用API.
  int tp_getIthLeaf(int ileaf,char* leaf);//!< カレントノードにあるileaf番目のリーフのラベルを取得します。C用API.
  int tp_getIthNodeOrder(int inode,char* node,int order);//!< カレントノードにあるinode番目のノードのラベルを取得します。C用API.
  int tp_getIthLeafOrder(int ileaf,char* leaf,int order);//!< カレントノードにあるileaf番目のリーフのラベルを取得します。C用API.


  //fortran 用 API
  int tp_read_fort_(char* file,int* length);//!< file からパラメータを読み込みます. Fortran用API.
  int tp_write_fort_(char* file,int* length);//!< file にパラメータのデータ構造全体を書き込みます. Fortran用API. 
  int tp_remove_fort_();//!< 格納しているパラメータのデータを破棄します。 Fortran用API. 
  int tp_get_number_of_leaves_fort_(int* nleaves );//!< リーフの総数を取得します。 Fortran用API. 
  int tp_get_label_fort_(int* ileaf,char* label,int* length);//!< インデックスで指定したリーフのラベルを取得します。 Fortran用API. 

  int tp_get_type_fort_(char* label,int* type,int* label_length);//!< ラベルで指定したリーフのタイプを取得します。 Fortran用API. 
  int tp_get_value_fort_(char* label,char* value,int* label_length,int* value_length);//!< ラベルで指定したリーフの値を取得します。 Fortran用API. 

  // 型変換用関数
  char tp_convert_char_fort_(char* value,
			     int *error,int* value_length); //!< パラメータの値を文字列から char 型へ変換します。Fortran用API.
  short tp_convert_short_fort_(char* value,
			int *error,int* value_length);//!< パラメータの値を文字列から short 型へ変換します。 Fortran用API.
  int tp_convert_int_fort_(char* value,
		     int *error,int* value_length); //!< パラメータの値を文字列から int 型へ変換します。Fortran用API. 
  long tp_convert_long_fort_(char* value,
		      int *error,int* value_length); //!< パラメータの値を文字列から long 型へ変換します。 Fortran用API.
   long long tp_convert_long_long_fort_(char* value,
  			       int *error,int* value_length);//!< パラメータの値を文字列から long long 型へ変換します。Fortran用API. 
  float tp_convert_float_fort_(char* value,
			       int *error,int* value_length);//!< パラメータの値を文字列から float 型へ変換します。 Fortran用API.
  double tp_convert_double_fort_(char* value,
			  int *error,int* value_length);//!< パラメータの値を文字列から double 型へ変換します。 Fortran用API.
  int tp_convert_logical_fort_(char* value,
		     int *error,int* value_length);//!< パラメータの値を文字列から int 型 でbool型を表現したものへ変換します。 Fortran用API.
  //ベクトル値の操作
  int tp_get_number_of_elements_fort_(char* vector_value,unsigned int* nvec,int* vector_length); //!<ベクトル値の要素数の取得  Fortran用API. 
  int tp_get_ith_element_fort_(char* vector_value, int* ivec,char* velem ,int* vector_value_length,int* velem_length);//!<ベクトル値のivec番目の要素を取得する。 Fortran用API. 
  


  //カレントノードの取得、ノードの移動 fortran　用 API
  int tp_current_node_fort_(char* label,int* label_length); //!<現在のカレントノードを取得します。 Fortran用API. 
  int tp_change_node_fort_(char* lable,int* label_length);//!< ラベルで指定したノードをカレントノードにします。 Fortran用API. 
  int tp_get_number_of_cnodes_fort_(int* nnodes);//!< カレントノードにある子ノードの数を取得します。 Fortran用API. 
  int tp_get_number_of_cleaves_fort_(int* nleaves);//!< カレントノードにあるリーフの数を取得します。 Fortran用API. 
  int tp_get_ith_node_fort_(int* inode,char* node,int* node_length);//!< カレントノードにあるinode番目のノードのラベルを取得します。 Fortran用API. 
  int tp_get_ith_leaf_fort_(int* ileaf,char* leaf,int* leaf_length);//!< カレントノードにあるileaf番目のリーフのラベルを取得します。 Fortran用API. 

  int tp_get_ith_node_order_fort_(int* inode,char* node,int* order,int* node_length);//!< カレントノードにあるinode番目のノードのラベルを取得します。 Fortran用API. 
  int tp_get_ith_leaf_order_fort_(int* ileaf,char* leaf,int* order,int* leaf_length);//!< カレントノードにあるileaf番目のリーフのラベルを取得します。 Fortran用API. 



#ifdef __cplusplus
}
#endif /*__cplusplus*/

#endif // __TEXTPARSER_H__


