// -*- mode: c++ -*-
/****************************************************************************
 **
 ** Copyright (C) 2012 Tokyo University.
 **
 ****************************************************************************/
/** @file TextParserTree.h
 * ここには TextParserTree クラスが定義されています。
 *
 */

#ifndef __TEXTPARSER_TREE_H__
#define __TEXTPARSER_TREE_H__


#include <iostream>
#include <string>
#include <string.h>
#include <sstream>
#include <map>
#include <algorithm>
#include <vector>

#include "TextParserCommon.h"
#include "TextParserElement.h"

/** パラメータファイルのデータ構造を保持するクラス
 *
 *
 *
 */

class TextParserTree{

 private:
  // constructors. not for using.
  TextParserTree(){} //!< デフォルトコンストラクタ 利用禁止
  TextParserTree(const TextParserTree& rhs){} //!< コピーコンストラクタ 利用禁止
  TextParserTree& operator=(const TextParserTree& rhs); //!< 代入演算子 利用禁止

  // private valiables

  static  bool _initialized;    //!< 初期化判定フラッグ
  bool _is_ready;    //!< アクセス可否判定フラッグ
  std::string _input_filename;   //!< 入力ファイル名
  TextParserParseMode _parse_mode;   //!< 現在の解析モード
  TextParserParseMode _return_parse_mode; //!< 戻り先のパースモード
  std::map<std::string, unsigned int> _array_label_number; //!< 配列形式ラベルの使用数
  TextParserBool _result_bool;   //!< 条件式の評価結果
  TextParserBool _current_bool; //!< 依存関係の : の左又は右を表す（左はTextParserTRUE、右はTestParserFALSE）
  std::map<unsigned int, TextParserLeaf *> _unresolved_leaves;  //!< 未解決の依存関係付き値
  bool _node_open;    //!< ノードの開閉状態
  

 public:
  static TextParserTree* get_instance(); //!< インスタンスの取得
  bool isReady();//!< データ構造にアクセス可能かどうか
  void read(const std::string* file); //!< 指定したファイルからのパラメータ読み込み


  // スタティックな関数
  static TextParserError SetArrayLabelIndex(std::string& label, std::map<std::string, unsigned int>& array_label_number);
  static TextParserError GetElementAbsolutePath(TextParserElement* element, std::string& path);//!< エレメント（相対パス）の絶対パスを取得
  static unsigned int GetLeafID();//!< リーフ固有のIDを取得
 
  // these should be private...
  std::string _label;  //!< ルートディレクトリのラベル
  std::map<std::string, TextParserNode *>_nodes; //!< TextParserNodeクラスのキー付き配列
  std::map<unsigned int, std::string> _leaf_paths; //!< リーフの絶対パスリスト
  TextParserElement* _current_element;  //!< 現在のエレメント
  bool _debug_write;  //!< デバッグ文の表示スイッチ

  //  unsigned int _current_line;   //!< 入力中のファイルの現在の行
  //  unsigned int _current_leaf_id;  //!< 現在のリーフのID
  // static化
  static  unsigned int _current_line; //!< 入力中のファイルの現在の行
  static  unsigned int _current_leaf_id; //!< 現在のリーフのID

 
   // エレメントの追加、削除と取得
  
  TextParserError addElement(TextParserNode *directory);//!< Elementの追加
  TextParserError removeElement();//!< Elementの削除
  TextParserNode *getNode(const std::string& label);//!< ノードの取得
  TextParserError getNode(std::string& label,
			  TextParserElement *parent_element,
			  TextParserNode **dirtectory);//!< ノードの取得
  
  
 // パラメータの入力
  void initialize();//!< パラメータデータ構造の初期化
  TextParserError readParameters(const std::string& filename);//<! 指定したファイルからパラメータを読み込む。


  // パラメータのパース
private:
  TextParserError parseLine(std::stringstream& ss, std::string& buffer); //!< 一行あるいは、バッファー内の文字列が解析可能になるまで読み込んで、バッファーを処理する。
  TextParserError removeCommentsAndSpaces(std::stringstream& ss, std::string& buffer);//!< コメントと空白を削除する。
  TextParserError removeComment(std::stringstream& ss, std::string& buffer, bool& answer);//!< 空白を削除する。
  TextParserError parseLabel(std::stringstream& ss, std::string& buffer, bool& answer, std::string& label); //!< バッファーの内容をラベルかどうか判定する。
  TextParserError parseLabel(std::string& buffer, bool& answer, std::string& label);//!< バッファーの内容をラベルかどうか判定する。
  bool isCorrectLabel(std::string& label, bool pathlabel);//!< ラベルの正しさを検証する。
  bool isPathLabel(const std::string& label); //!< ラベルがパスであるかどうかを判定する。
  bool isIncludedInPath(const std::string& label, TextParserNode *node); //!< ノードが上位のパスに含まれていないか検証する。
  bool isArrayLabelExist(const std::string& label, TextParserElement *parent_element, const TextParserElementType type);//!< 配列ラベルを示す同じエレメントが存在するか判定
  TextParserError setArrayLabelIndex(std::string& label); //!< 配列型ラベルならばインデックスを割り当てる。
  TextParserError parseNode(std::stringstream& ss, std::string& buffer, bool& answer); //!< ノードかどうかパースする。
  TextParserError parseEndOfNode(std::stringstream& ss, std::string& buffer, bool& answer); //!< ノードの終わりかどうかパースする。
  TextParserError openNode(const std::string& label); //!<ノードの開始
  TextParserError closeNode();//!< ノードの終了
  TextParserError getOrAddNode(std::string& label, TextParserElement *parent_element, TextParserNode **node);//!< ディレクトリの取得又は追加
  TextParserError parseLeaf(std::stringstream& ss, std::string& buffer, bool& answer);//!<リーフの判定
  TextParserError openLeaf(const std::string& label);
  TextParserError closeLeaf();
  TextParserError parseDependValue(std::stringstream& ss, std::string& buffer, bool& answer);
  TextParserError parseVectorValue(std::stringstream& ss, std::string& buffer, bool& answer);
  TextParserError parseVectorValue(std::string& buffer, bool& answer);
  TextParserError parseUndefinedValue(std::stringstream& ss, std::string& buffer, bool& answer);
  TextParserError parseUndefinedValue(std::string& buffer, bool& answer);
  TextParserError parseValue(std::stringstream& ss, std::string& buffer, bool& answer, std::string& value, TextParserValueType& value_type);
  TextParserError parseValue(std::string& buffer, bool& answer, std::string& value, TextParserValueType& value_type);
  TextParserError parseOpenBrancket(std::stringstream& ss, std::string& buffer, bool& answer);
  TextParserError parseOpenBrancket(std::string& buffer, bool& answer);
  TextParserError parseClosedBrancket(std::stringstream& ss, std::string& buffer, bool& answer);
  TextParserError parseClosedBrancket(std::string& buffer, bool& answer);
  TextParserError parseEqual(std::stringstream& ss, std::string& buffer, bool& answer);
  TextParserError parseEqual(std::string& buffer, bool& answer);
  TextParserError parseNotEqual(std::stringstream& ss, std::string& buffer, bool& answer);
  TextParserError parseNotEqual(std::string& buffer, bool& answer);
  TextParserError parseAnd(std::stringstream& ss, std::string& buffer, bool& answer);
  TextParserError parseAnd(std::string& buffer, bool& answer);
  TextParserError parseOr(std::stringstream& ss, std::string& buffer, bool& answer);
  TextParserError parseOr(std::string& buffer, bool& answer);
  TextParserError parseQuestion(std::stringstream& ss, std::string& buffer, bool& answer);
  TextParserError parseQuestion(std::string& buffer, bool& answer);
  TextParserError parseColon(std::stringstream& ss, std::string& buffer, bool& answer);
  TextParserError parseColon(std::string& buffer, bool& answer);
  TextParserError parseDelimiter(std::stringstream& ss, std::string& buffer, bool& answer);
  TextParserError parseDelimiter(std::string& buffer, bool& answer);
  TextParserError parseEndOfVector(std::stringstream& ss, std::string& buffer, bool& answer);
  TextParserError parseEndOfVector(std::string& buffer, bool& answer);
  bool isCorrectValue(std::string& value, TextParserValueType& value_type);
  TextParserError setValue(std::string& value, TextParserValueType& value_type);
  TextParserError setUndefinedValue();
  TextParserError setVectorValue(std::stringstream& ss, std::string& buffer);
  TextParserError setDependenceExpression(std::stringstream& ss, std::string& buffer);
  TextParserError setConditionalExpression(std::stringstream& ss, std::string& buffer);
  TextParserError setDependenceValue(std::stringstream& ss, std::string& buffer);
  TextParserError parseDependenceExpression(TextParserLeaf *leaf);
  TextParserError parseConditionalExpression(std::string& buffer, TextParserBool& result);
  TextParserError parseDependenceValue(std::string& buffer, TextParserBool set);
  TextParserError resolveConditionalExpression(std::string& label, std::string& value, TextParserValueType& value_type, const TextParserBool is_equal, TextParserBool& result);
  TextParserBool resolveAnd(TextParserBool& left, TextParserBool& right);
  TextParserBool resolveOr(TextParserBool& left, TextParserBool& right);
  
  int array_label_test(const std::string &label,std::string &key);


public:
  TextParserError getElementRelativePath (std::string& path, bool add, TextParserElement **parent_element);
  TextParserError getElement(const std::string& path_label, const TextParserElementType type, TextParserElement** element);
  
  // original...
  //  TextParserError getLeafValue(std::string& path, TextParserValue **value); 
  
  TextParserError getLeafValue(const std::string &path, TextParserValue **value);
  
  // パラメータの出力
  TextParserError writeParameters(const std::string& filename,int order=0);
  // パラメータの取得
  TextParserError getLeafLabels(std::map<unsigned int, std::string>& labels);
  TextParserError changeNode(const std::string& label);
  TextParserError getCurrentNode(std::string& label);
  TextParserError getCurrentNodeNumber(unsigned int *number);
  TextParserError getCurrentNodeLabel(unsigned int id, char **label);
  TextParserError getCurrentNodeLabels(std::vector<std::string>& labels);
  TextParserError getCurrentLeafNumber(unsigned int *number);
  TextParserError getCurrentLeafLabel(unsigned int id, char **label);
  TextParserError getCurrentLeafLabels(std::vector<std::string>& labels);
  TextParserError splitVectorValue(const std::string &vector_value,
				   std::vector<std::string>& values);
  
  // デバッグ表示
    void debugWrite(bool swt);
    void debugWriteValue(const std::string& value, const TextParserValueType& value_type);
    void debugWriteTextParserBool(const TextParserBool& result);

  // 並べ替え
  TextParserError labelSort(const std::vector<std::string>& input,
			    std::vector<std::string>& output,
			    int iswitch);
  TextParserError labelSort_1(const std::vector<std::string>& input,
			      std::vector<std::string>& output);
  TextParserError labelSort_2(const std::vector<std::string>& input,
			      std::vector<std::string>& output);
  TextParserError nodeSort(const std::vector<std::string>& input,
			   std::vector<std::string>& output,
			   int iswitch);
  TextParserError nodeSort_2(const std::vector<std::string>& input,
			     std::vector<std::string>& output);
  
};

//文字列処理
std::string CorrectValueString(std::string buffer);
std::string TextParserRemoveHeadSpaces(std::string buffer);
std::string TextParserRemoveTailSpaces(std::string buffer);
std::string TextParserStringToLower(const std::string& str);
bool TextParserStringCompare(const std::string& str0, const std::string& str1);

// エラー処理
TextParserError TextParserErrorHandler(const TextParserError error_code, const std::string& sub_message);



#endif // __TEXTPARSER_TREE_H__
