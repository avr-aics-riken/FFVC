/****************************************************************************
 **
 ** Copyright (C) 2012 Tokyo University.
 **
 ****************************************************************************/
/** @file TextParserTree.cpp
 * ここには TextParserTree クラスが実装されています。
 *
 */
#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif //HAVE_CONFIG_H


#ifdef ENABLE_MPI 
#ifdef BUILD_MPI
#include "mpi.h"
// nothig to do ...
#endif
#else
#endif



#include "TextParserTree.h"

unsigned int TextParserTree::_current_line=0;
unsigned int TextParserTree::_current_leaf_id=0;
bool TextParserTree::_initialized=false;

// 値の次に来る可能性のある文字
// 　通常だと'\t' '\n' '\r'
// 　ベクトルの次の値だと','ベクトルの終わりだと ')'
// 　次のディレクトリやリーフだと '\"'
// 　依存関係付き値だと ':"'
//   コメントだと'/'
// std::string MgppValueNext[] = {"\t", "\n", "\r", ",", ")", "\"", ":", "/"};
// unsigned int MgppNumValueNext = 8;
std::string TextParserValueNext[] = {"\t", "\n", "\r", ",", ")", "\"", ":", "/"};
unsigned int TextParserNumValueNext = 8;
char charTextParserValueNext[] = {'\t', '\n', '\r', ',', ')', '\"', ':', '/'};




/** プロセス内唯一の TextParserTree インスタンスを返します。
 *
 * @return TextParserTree インスタンスのアドレス
 * 
 */
TextParserTree* TextParserTree::get_instance(){

#ifdef MYDEBUG
  std::cout<< "TextParserTree::get_instance() called."<<std::endl;
#endif //MYDEBUG

  static TextParserTree instance;

#ifdef MYDEBUG
  std::cout<< "address to TextParserTree :"<< &instance <<std::endl;
#endif //MYDEBUG

  if(!instance._initialized ){
    instance._initialized=true;
    instance._is_ready = false;
    instance._label = "/";
    instance._parse_mode = TP_NO_PARSE;
    instance._return_parse_mode = TP_NO_PARSE;
    instance._current_element = 0;
    instance._current_line = 0;
    instance._current_leaf_id = 0;
    instance._node_open = false;
    instance._result_bool = TP_UNDEFINED_BOOL;
    instance._current_bool = TP_UNDEFINED_BOOL;
    instance._debug_write = false;
#ifdef MYDEBUG
    std::cout<< "TextParserTree is initialized."<< std::endl;
#endif //MYDEBUG
    

  }
  return &instance;
}



/** 配列ラベルのインデックス設定
 *
 * @param[in] label ラベル
 * @param[in] array_label_number 配列形式ラベルの使用数格納データ
 * @return リターンコード
 *
 * 配列ラベル（末尾が[@]）だったらインデックスを割り当てる
 *
 */
TextParserError TextParserTree::SetArrayLabelIndex(std::string& label, std::map<std::string, unsigned int>& array_label_number)
{
  int a = label.find("[@]");
  if (a >= 0) { // 配列タイプのラベル
    if (a == label.size() - 3) {
      std::string label0 = TextParserStringToLower(label.substr(0, a)); // ラベルの文字列部分を取り出し小文字に変換
      unsigned int number = 0;
      std::map<std::string, unsigned int>::iterator li = array_label_number.find(label0);
      if (li != array_label_number.end()) {  // ラベルが既に存在する
	number = li->second;            // 配列添え字を更新
	number++;
	li->second = number;
      } else {                                // ラベルを登録
	array_label_number.insert(str_int(label0, number));
      }
#if 0
      std::ostringstream os;
      os << number;
      label = label0 + "[" + os.str() + "]";
#else
      char os[20];
      sprintf(os, "%d", number);
      label = label0 + "[" + os + "]";
#endif
    } else {
      return TP_ILLEGAL_ARRAY_LABEL_ERROR;
    }
  }

  return TP_NO_ERROR;
}




/** パラメータファイルの読み込み
 *
 * @param[in] filename ファイル名
 * @return エラーコード
 *
 */
TextParserError TextParserTree::readParameters(const std::string& filename)
{

  TextParserError ret = TP_NO_ERROR;  // 戻り値
  if (_is_ready) {
    return TextParserErrorHandler(TP_DATABASE_ALREADY_SET_ERROR, "");
  }
  try{
    std::string line;
    std::stringstream ss;
    
#ifdef BUILD_MPI
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    if(myrank==0){
#endif
      
      std::ifstream ifs(filename.c_str());
      while (getline(ifs,line)){
	ss << line << '\n';
      }
      ifs.close();
      
      
#ifdef BUILD_MPI
    } // myrank==0 end 
#endif
    
#ifdef BUILD_MPI
    int ssize=ss.str().size();
    MPI_Bcast(&ssize,1,MPI_INT,0,MPI_COMM_WORLD);
    
    char mpi_buffer[ssize+1];
    if(myrank==0){
      strcpy(mpi_buffer,ss.str().c_str());
    }
    
    MPI_Bcast(mpi_buffer,ssize+1,MPI_CHAR,0,MPI_COMM_WORLD);
    
    if(myrank!=0){
      ss<<mpi_buffer;
    }
#endif
    
    // now ss has an entire file contents.
    
    if (ss)
      {
	_unresolved_leaves.clear();
	_leaf_paths.clear();
	_parse_mode = TP_NODE_PARSE; // ノード解析モード
	_current_element = 0;   // カレントエレメントをルートディレクトリに設定
	_current_line = 0;       // ファイルのカレント行を初期化
	_current_leaf_id = 0;   // 現在のリーフのIDを初期化
	_node_open = false; // ノードの開閉状態を閉に設定
	
	
	std::string buffer;
	while (getline(ss, buffer)) {
	  _current_line++;

//	  TextParserError ret = parseLine(ifs, buffer);     // 一行のパラメータ解析
	  TextParserError ret = parseLine(ss, buffer);     // 一行のパラメータ解析

	  if (ret != TP_NO_ERROR) {
	    initialize();
	    return ret;
	  }
	} //while(getline(loop))

	//! 未解決の依存関係付き値の解決
	if (_unresolved_leaves.size() > 0) {
#ifdef TP_DEBUG
	  if (_debug_write) {
	    std::cout << "☆ 未解決の依存関係付き値の総数 = " << _unresolved_leaves.size() << std::endl;
	  }
#endif
	  while (_unresolved_leaves.size() > 0) {
	    unsigned int unresolved_number = _unresolved_leaves.size();
	    std::map<unsigned int, TextParserLeaf *>::iterator li = _unresolved_leaves.begin();
	    while (li != _unresolved_leaves.end()) {
	      _current_line = li->first;
	      _current_element = li->second;
	      ret = parseDependenceExpression(li->second);
	      if (ret == TP_NO_ERROR) {
#ifdef TP_DEBUG
		if (_debug_write) {
		  std::cout << "★ 解決した依存関係付き値 = " << li->second->_label << std::endl;
		}
#endif
		std::map<unsigned int, TextParserLeaf *>::iterator lt = li;
		li++;
		_unresolved_leaves.erase(lt);   // 削除
	      } else {
		if (ret == TP_MISSING_PATH_ELEMENT_ERROR
		    || ret == TP_UNRESOLVED_LABEL_USED_WARNING) {
		  li++;
		} else {
		  initialize();
		  return ret;
		}
	      }
	    }
	    // 解決されなかったら終了
	    if (unresolved_number == _unresolved_leaves.size()) break;
	  }
	  if (_unresolved_leaves.size() > 0) {
	    // 未解決の依存関係付き値が残った場合未定義値を設定
	    std::map<unsigned int, TextParserLeaf *>::iterator li = _unresolved_leaves.begin();
	    while (li != _unresolved_leaves.end()) {
	      _current_line = li->first;
	      _current_element = li->second;
	      ret = parseDependenceExpression(li->second);
	      std::string expression = "";
	      TextParserValue *value = li->second->_value;
	      expression = li->second->_label + " = " + value->_value;
	      if (ret != TP_NO_ERROR) {
		if (ret == TP_UNRESOLVED_LABEL_USED_WARNING) {
		  TextParserErrorHandler(TP_UNRESOLVED_LABEL_USED_WARNING, expression);
		  ret = setUndefinedValue();                  // 未定義値の設定
		  if (ret != TP_NO_ERROR) {
		    initialize();
		    return ret;
		  }
		} else {
		  if (ret == TP_MISSING_PATH_ELEMENT_ERROR) {
		    TextParserErrorHandler(TP_MISSING_PATH_ELEMENT_ERROR, expression);
		  }
		  initialize();
		  return ret;
		}
	      }
	      li++;   // 最初のリーフでエラーになるはずだから本来不要
	    }
	  }
	}
      } else {
      ret = TextParserErrorHandler(TP_FILEOPEN_ERROR, filename);
      initialize();
      return ret;
    }
    if (_node_open) {
      ret = TextParserErrorHandler(TP_NODE_END_MISSING_ERROR, filename);
      initialize();
      return ret;
    }
    ret = getLeafLabels(_leaf_paths);
    if (ret != TP_NO_ERROR) {
      initialize();
      return ret;
    }
    _input_filename = filename;
    _is_ready = true;  // データへのアクセスを許可
  } catch (std::ios_base::failure) {
    ret = TextParserErrorHandler(TP_FILEINPUT_ERROR, filename);
  } catch (std::exception ex) {
    ret = TextParserErrorHandler(TP_FILEINPUT_ERROR, filename);
  }
  initialize();

#ifdef MYDEBUG
  std::cout<< "readParameters(const std::string& filename) end." <<std::endl;
#endif
  return ret;
}


/** エレメントのパスを取得
 *
 * @param[in] element エレメント
 * @param[out] path パス
 * @return リターンコード
 *
 */
TextParserError TextParserTree::GetElementAbsolutePath(TextParserElement* element, std::string& path)
{
  std::vector<std::string> nodes;
  nodes.push_back(element->_label);
  TextParserElement *parent = element->_parent;
  while (parent != 0) {
    nodes.push_back(parent->_label);
    parent = parent->_parent;
  }
  path = "";
  for (int id = nodes.size() - 1; id >= 0; id--) {
    path += "/" + nodes[id];
  }


  return TP_NO_ERROR;
}


/** リーフのIDを取得
 *
 * @return リーフID
 *
 */
unsigned int TextParserTree::GetLeafID()
{
#ifdef MYDEBUG
  std::cout <<"TextParserTree::GetLeafID() started"<<std::endl;
#endif //MYDEBUG
  //! インスタンスの取得

  // TextParserTree *db = TextParserTree::get_instance();
  // unsigned int id = db->_current_leaf_id;
  // db->_current_leaf_id++;

  unsigned int id = _current_leaf_id;
  _current_leaf_id++;

#ifdef MYDEBUG
  std::cout <<"extparserTree::GetLeafID() end"<<std::endl;
#endif //MYDEBUG

  return id;

}




// 文字列処理ルーチン

/** 数値の文字列を修正する
 *
 * @param[in] buffer 入力文字列
 * @return 出力文字列
 *
 * 空白文字を削除し、指数部のdやDはeに変える
 *
 */
std::string CorrectValueString(std::string buffer)
{
  unsigned int c = 0;
  while (buffer.size() > 0) {
    if (buffer[c] == ' ')
      {
	buffer = buffer.erase(c, 1);
      } else {
      if (buffer[c] == 'd') buffer[c] = 'e';  // 指数をeに変換
      if (buffer[c] == 'D') buffer[c] = 'e';  // 指数をeに変換
      c++;
      if (c >= buffer.size()) break;
    }
  }

  return buffer;
}

/** 頭の空白文字を削除する
 *
 * @param[in] buffer 入力文字列
 * @return 出力文字列
 *
 */
std::string TextParserRemoveHeadSpaces(std::string buffer)
{
  if (buffer.size() == 0) return buffer;

  for (unsigned int c = 0; c < buffer.size(); c++) {
    if (buffer[c] != ' ' && buffer[c] != '\t' && buffer[c] != '\r' && buffer[c] != '\n')
      {
	buffer = buffer.substr(c);
	return buffer;
      }
  }
  return "";
}

/** 末尾の空白文字を削除する
 *
 * @param[in] buffer 入力文字列
 * @return 出力文字列
 *
 */
std::string TextParserRemoveTailSpaces(std::string buffer)
{
  if (buffer.size() == 0) return buffer;

  for (unsigned int c = buffer.size() - 1; c >= 0; c--) {
    if (buffer[c] != ' ' && buffer[c] != '\t' && buffer[c] != '\r' && buffer[c] != '\n')
      {
	if (c < buffer.size() - 1) {
	  buffer = buffer.substr(0, c + 1);
	}
	return buffer;
      }
  }
  return "";
}

/** 文字列を小文字に変換
 *
 * @param[in] str 文字列
 * @return 文字列
 *
 */
std::string TextParserStringToLower(const std::string& str)
{
  std::string str_cpy(str);
  // 小文字に変換
  transform(str_cpy.begin(), str_cpy.end(), str_cpy.begin(), ::tolower);
    
  return str_cpy;
}

/** 大文字小文字を無視して文字列を比較
 *
 * @param[in] str0 文字列
 * @param[in] str1 文字列
 * @return 比較結果　（true : 等しい、false : 異なる）
 *
 */
bool TextParserStringCompare(const std::string& str0, const std::string& str1)
{
  if (str0.size() != str1.size()) return false;   // 長さが違う

  std::string str0_cpy = TextParserStringToLower(str0);
  std::string str1_cpy = TextParserStringToLower(str1);
    
  return (str0_cpy == str1_cpy);
}

/** エラーメッセージを表示する
 *
 * @param[in] error_code エラーコード
 * @param[in] sub_message サブメッセージ
 * @return エラーコード
 *
 */
TextParserError TextParserErrorHandler(const TextParserError error_code, const std::string& sub_message)
{
#ifdef MYDEBUG    
  std::cout<< "TextParserErrorHandler() start"<<std::endl;
  std::cout<< "error_code " <<error_code <<" submessage " <<sub_message<< std::endl;
#endif // MYDEBUG    


  if (error_code > TP_NO_ERROR) {
    if (error_code < TP_WARNING ) {
      std::cerr << "*Error #" << error_code << ": ";
    } else {
      std::cerr << "*Warning #" << error_code << ": ";
    }
    switch (error_code) {
    case TP_DATABASE_NOT_READY_ERROR:
      std::cerr << "Database is not ready";
      break;
    case TP_DATABASE_ALREADY_SET_ERROR:
      std::cerr << "Database has been already set";
      break;
    case TP_FILEOPEN_ERROR:
      std::cerr << "File open failed";
      break;
    case TP_FILEINPUT_ERROR:
      std::cerr << "File input failed";
      break;
    case TP_FILEOUTPUT_ERROR:
      std::cerr << "File outnput failed";
      break;
    case TP_ENDOF_FILE_ERROR:
      std::cerr << "End of file";
      break;
    case TP_ILLEGAL_TOKEN_ERROR:
      std::cerr << "Illegal token";
      break;
    case TP_MISSING_LABEL_ERROR:
      std::cerr << "Missing label";
      break;
    case TP_ILLEGAL_LABEL_ERROR:
      std::cerr << "Illegal label";
      break;
    case TP_ILLEGAL_ARRAY_LABEL_ERROR:
      std::cerr << "Illegal array type label";
      break;
    case TP_MISSING_ELEMENT_ERROR:
      std::cerr << "Missing element";
      break;
    case TP_ILLEGAL_ELEMENT_ERROR:
      std::cerr << "Illegal element";
      break;
    case TP_NODE_END_ERROR:
      std::cerr << "too much Node end";
      break;
    case TP_NODE_END_MISSING_ERROR:
      std::cerr << "Node termination is Missing";
      break;
    case TP_NODE_NOT_FOUND_ERROR:
      std::cerr << "The Node is not found";
      break;
    case TP_LABEL_ALREADY_USED_ERROR:
      std::cerr << "Label is already used";
      break;
    case TP_LABEL_ALREADY_USED_PATH_ERROR:
      std::cerr << "Label is already used in path";
      break;
    case TP_ILLEGAL_CURRENT_ELEMENT_ERROR:
      std::cerr << "Illegal current element ";
      break;
    case TP_ILLEGAL_PATH_ELEMENT_ERROR:
      std::cerr << "Illegal path element ";
      break;
    case TP_MISSING_PATH_ELEMENT_ERROR:
      std::cerr << "Missing path element";
      break;
    case TP_ILLEGAL_LABEL_PATH_ERROR:
      std::cerr << "Illegal label path";
      break;
    case TP_UNKNOWN_ELEMENT_ERROR:
      std::cerr << "Unknown element";
      break;
    case TP_MISSING_EQUAL_NOT_EQUAL_ERROR:
      std::cerr << "Missing both == and !=";
      break;
    case TP_MISSING_AND_OR_ERROR:
      std::cerr << "Missing both && and ||";
      break;
    case TP_MISSING_CONDITION_EXPRESSION_ERROR:
      std::cerr << "Missing condition expression";
      break;
    case TP_ILLEGAL_DEPENDENCE_EXPRESSION_ERROR:
      std::cerr << "Illegal dependence expression";
      break;
    case TP_MISSING_CLOSED_BRANCKET_ERROR:
      std::cerr << "Missing closed brancket";
      break;
    case TP_ILLEGAL_CONDITION_EXPRESSION_ERROR:
      std::cerr << "Illegal condition expression";
      break;
    case TP_MISSING_VALUE_ERROR:
      std::cerr << "Missing value";
      break;
    case TP_ILLEGAL_VALUE_ERROR:
      std::cerr << "Illegal value";
      break;
    case TP_ILLEGAL_NUMERIC_VALUE_ERROR:
      std::cerr << "Illegal numeric value";
      break;
    case TP_ILLEGAL_VALUE_TYPE_ERROR:
      std::cerr << "Illegal value type";
      break;
    case TP_MISSING_VECTOR_END_ERROR:
      std::cerr << "Missing vector end";
      break;
    case TP_VALUE_CONVERSION_ERROR:
      std::cerr << "Value conversion failed";
      break;
    case TP_MEMORY_ALLOCATION_ERROR:
      std::cerr << "Memory allocation failed";
      break;
    case TP_MISSING_COMMENT_END_ERROR:
      std::cerr << "Missing comment end";
      break;
    case TP_ID_OVER_ELEMENT_NUMBER_ERROR:
      std::cerr << "ID is over the element number";
      break;
    case TP_GET_PARAMETER_ERROR:
      std::cerr << "Get parameter failed";
      break;
    case TP_UNSUPPORTED_ERROR:
      std::cerr << "Unsupported function";
      break;
    case TP_UNDEFINED_VALUE_USED_WARNING:
      std::cerr << "Undefined value used";
      break;
    case TP_UNRESOLVED_LABEL_USED_WARNING:
      std::cerr << "Unresolved label used";
      break;
    default:
      std::cerr << "Undefined error code";
      break;
    }
    unsigned int current_line = (TextParserTree::get_instance())->_current_line;
    //unsigned int current_line = TextParserTree::getInstance()->_current_line;
    if (sub_message.size() > 0) std::cerr << " : " + sub_message;
    if (current_line > 0) std::cerr << " : line " << current_line;
    std::cerr << std::endl;
  }


#ifdef MYDEBUG    
  std::cout<< "TextParserErrorHandler() end"<<std::endl;
#endif // MYDEBUG    

  return error_code;
}


/** 一行のパラメータ解析
 *
 * @param[in] ss 入力ファイルポインタ
 * @param[in,out] buffer 入力文字列
 * @return エラーコード
 *
 */
TextParserError TextParserTree::parseLine(std::stringstream& ss, std::string& buffer)
{
#ifdef MYDEBUG
  std::cout << "parseLine(std::stringstream& ss, std::string& buffer) start" 
	    <<std::endl;
  std::cout <<"buffer "<<buffer << std::endl;
#endif //MYDEBUG
  TextParserError ret;                                          // 戻り値
  bool answer;
  std::string comment;                                    // コメント
  std::string label;                                      // ラベル
  std::string value;                                      // 値
  TextParserValueType value_type;                               // 値のタイプ
  unsigned int line;

  while (buffer.size() > 0) {
    ret = removeCommentsAndSpaces(ss, buffer);
    if (ret != TP_NO_ERROR) {
      if (ret == TP_ENDOF_FILE_ERROR) ret = TP_NO_ERROR;
      return ret;
    }
    switch (_parse_mode) {
    case TP_NODE_PARSE:
      ret = parseLabel(ss, buffer, answer, label);
      if (ret != TP_NO_ERROR) return ret;
      if (answer) {                   // ラベル?
	// ラベルの異常検出
	if (!isCorrectLabel(label, true)) {
	  return TextParserErrorHandler(TP_ILLEGAL_LABEL_ERROR, label);
	}
	// エレメントの解析(ディレクトリ又はリーフが必要)
	ret = parseNode(ss, buffer, answer);
	if (ret != TP_NO_ERROR) return ret;
	if (answer) {                  // ディレクトリ?
	  ret = openNode(label);             // ディレクトリの開始
	  if (ret != TP_NO_ERROR) return ret;
	  _parse_mode = TP_NODE_PARSE;     // ディレクトリ解析モード
	  break;
	}
	ret = parseLeaf(ss, buffer, answer);
	if (ret != TP_NO_ERROR) {
	  if (ret == TP_ENDOF_FILE_ERROR
	      && _current_element == 0) {
	  } else {
	    return ret;
	  }
	}
	if (answer) { // リーフ?
	  ret = openLeaf(label); // リーフの開始
	  if (ret != TP_NO_ERROR) return ret;
	  _parse_mode = TP_LEAF_PARSE; // ディレクトリ解析モード
	  break;
	}
	return TextParserErrorHandler(TP_MISSING_ELEMENT_ERROR, buffer);
      }
      ret = parseEndOfNode(ss, buffer, answer);
      if (ret != TP_NO_ERROR) {
	if (ret == TP_ENDOF_FILE_ERROR
	    && _current_element != 0 && _current_element->_parent == 0) {
	} else {
	  return ret;
	}
      }
      if (answer) {  // ディレクトリの終わり?
	ret = closeNode();  // ディレクトリの終了
	if (ret != TP_NO_ERROR) return ret;
	break;
      }
      ret = parseDelimiter(ss, buffer, answer);
      if (answer) {                                       // デリミタ?
	// ディレクトリ間、リーフ間のデリミタを無視
	break;
      }
      if (buffer.size() > 0) {
	return TextParserErrorHandler(TP_ILLEGAL_TOKEN_ERROR, buffer);
      }
      break;
    case TP_LEAF_PARSE:
      line = _current_line; // 現在の行数の保存（未解決の依存関係付き値の保存時に使用するため）
      ret = parseDependValue(ss, buffer, answer);
      if (ret != TP_NO_ERROR) return ret;
      if (answer) {   // 依存関係付き値?
	ret = setDependenceExpression(ss, buffer);
	if (ret != TP_NO_ERROR) return ret;
	ret = parseDependenceExpression((TextParserLeaf *)_current_element);
	if (ret != TP_NO_ERROR) {
	  if (ret == TP_MISSING_PATH_ELEMENT_ERROR
	      || ret == TP_UNRESOLVED_LABEL_USED_WARNING) {
	    _unresolved_leaves.insert(uint_leaf(line, (TextParserLeaf *)_current_element));
	  } else {
	    return ret;
	  }
	} 
	ret = closeLeaf();                          // リーフの終了
	if (ret != TP_NO_ERROR) return ret;
	break;
      } // parseDependValue answer end
      ret = parseVectorValue(ss, buffer, answer);
      if (ret != TP_NO_ERROR) return ret;
      if (answer) {             // ベクトル値?
	ret = setVectorValue(ss, buffer);
	if (ret != TP_NO_ERROR) return ret;
	ret = closeLeaf();                          // リーフの終了
	if (ret != TP_NO_ERROR) return ret;
	break;
      }
      ret = parseUndefinedValue(ss, buffer, answer);
      if (ret != TP_NO_ERROR) return ret;
      if (answer) {          // 未定義の値?
	ret = setUndefinedValue();                  // 未定義値の設定
	if (ret != TP_NO_ERROR) return ret;
	ret = closeLeaf();                          // リーフの終了
	if (ret != TP_NO_ERROR) return ret;
	break;
      }
      ret = parseValue(ss, buffer, answer, value, value_type);
      if (ret != TP_NO_ERROR) return ret;
      if (answer) {// スカラー値?
	if (!isCorrectValue(value, value_type)) {    // 値のチェック
	  return TextParserErrorHandler(TP_ILLEGAL_VALUE_ERROR, value);
	}
	ret = setValue(value, value_type);          // 値の設定
	if (ret != TP_NO_ERROR) return ret;
	ret = closeLeaf();                          // リーフの終了
	if (ret != TP_NO_ERROR) return ret;
	break;
      }
      return TextParserErrorHandler(TP_ILLEGAL_VALUE_ERROR, buffer);
      break;
    default:
      break;
    }
  }
#ifdef MYDEBUG
  std::cout << "parseLine(std::stringstream& ss, std::string& buffer) start" 
	    <<std::endl;
#endif //MYDEBUG

  return TP_NO_ERROR;
}



/** データベースの状態
 *
 * @return データベースの状態
 *
 */
bool TextParserTree::isReady()
{
  return _is_ready;
}


/** コメントと空白の削除
 * @param[in] ss 入力ファイルストリーム
 * @param[in,out] buffer 入力文字列
 * @return エラーコード
 *
 */
TextParserError TextParserTree::removeCommentsAndSpaces(std::stringstream& ss, std::string& buffer)
{
#ifdef MYDEBUG
  std::cout << 
    "removeCommentsAndSpaces(std::stringstream& ss, std::string& buffer) start"
	    <<  std::endl;
  std::cout << "buffer "<<buffer <<std::endl;
#endif //MYDEBUG
  bool answer = false;
  while (true) {
    buffer = TextParserRemoveHeadSpaces (buffer);             // 先頭の空白削除
    TextParserError ret = removeComment(ss, buffer, answer); // コメントの削除
    if (ret != TP_NO_ERROR) return ret;
    while (buffer.size() == 0) {
      if (getline(ss, buffer)) {
	_current_line++;
      } else {
	return TP_NO_ERROR;
      }
      buffer = TextParserRemoveHeadSpaces (buffer);             // 先頭の空白削除
    }
    if (!answer) break; // コメントを削除しなければ新たな空白は生まれないので終了する
  }

#ifdef MYDEBUG
  std::cout << 
    "removeCommentsAndSpaces(std::stringstream& ss, std::string& buffer) end"
	    <<  std::endl;
  std::cout << "buffer "<<buffer <<std::endl;
#endif //MYDEBUG

  return TP_NO_ERROR;
}




//:::::::::::::::::::::::::::::::
/*
// 値の次に来る可能性のある文字
// 　通常だと'\t' '\n' '\r'
// 　ベクトルの次の値だと','ベクトルの終わりだと ')'
// 　次のディレクトリやリーフだと '\"'
// 　依存関係付き値だと ':"'
//   コメントだと'/'
std::string MgppValueNext[] = {"\t", "\n", "\r", ",", ")", "\"", ":", "/"};
unsigned int MgppNumValueNext = 8;
*/


#if 1
/** エレメントの追加
 *
 * @param[in] node TextParserNodeのポインタ
 * @return エラーコード
 *
 * ディレクトリエレメントを追加します
 * キーとなるラベルは小文字に変換します
 *
 */
TextParserError TextParserTree::addElement(TextParserNode *node)
{
  std::string label = TextParserStringToLower(node->_label);
  _nodes.insert(str_node(label, node));

  return TP_NO_ERROR;
}
#endif
/** エレメントの削除
 *
 * @return エラーコード
 *
 */
TextParserError TextParserTree::removeElement()
{
  // データベースは読み込みが失敗しても削除できるようにする
  //if (!_is_ready) {
  //    return TextParserErrorHandler(TP_DATABASE_NOT_READY_ERROR, "");
  //}
  TextParserError ret = TP_NO_ERROR;
  _is_ready = false;
  try {
    std::map<std::string, TextParserNode *>::iterator di = _nodes.begin();
    while (di != _nodes.end()) {
      TextParserNode *dir = (di)->second;
      ret = dir->removeElement();
      if (ret != TP_NO_ERROR) {
	TextParserErrorHandler(ret, dir->_label);
	return ret;
      }
      delete dir;
      di++;
    }
    _nodes.clear();

    _array_label_number.clear();
    _unresolved_leaves.clear();
    _leaf_paths.clear();
  } catch (std::exception ex) {
    ret = TextParserErrorHandler(TP_REMOVE_ELEMENT_ERROR, "");
  }

  return ret;
}

/** ディレクトリの取得
 *
 * @param[in] label ラベル
 * @return ディレクトリ
 *
 * ラベルの一致するディレクトリを取得
 *
 */
TextParserNode *TextParserTree::getNode(const std::string& label)
{
  std::string label_cpy = TextParserStringToLower(label);
  std::map<std::string, TextParserNode *>::iterator di = _nodes.find(label_cpy);
  if (di != _nodes.end()) {
    return di->second;
  } else {
    return 0;
  }
}

/** ディレクトリの取得
 *
 * @param[in] label ラベル
 * @param[in] parent_element 親エレメント
 * @param[out] node ディレクトリ
 * @return リターンコード
 *
 * ラベルに一致するディレクトリを取得する
 *
 */
TextParserError TextParserTree::getNode(std::string& label, TextParserElement *parent_element, TextParserNode **node)
{
  *node = 0;
  if (parent_element == 0) {							// 親エレメントがルートディレクトリ
    *node = getNode(label);               // 既存のディレクトリを取得
    if (*node == 0) {   // 存在しない
      return TextParserErrorHandler(TP_NODE_NOT_FOUND_ERROR, label);
    }
  } else if (parent_element->_type == TP_NODE_ELEMENT) {
    TextParserNode *parent_dir = (TextParserNode *)parent_element;
    *node = parent_dir->getNode(label);	// 既存のディレクトリを取得
    if (*node == 0) {   // 存在しない
      return TextParserErrorHandler(TP_NODE_NOT_FOUND_ERROR, label);
    }
  } else {
    return TextParserErrorHandler(TP_ILLEGAL_CURRENT_ELEMENT_ERROR, "");
  }

  return TP_NO_ERROR;
}

/** パラメータファイルの読み込み後の変数の初期化
 *
 */
void TextParserTree::initialize()
{
  _parse_mode = TP_NO_PARSE;
  _return_parse_mode = TP_NO_PARSE;
  _current_element = 0;
  _current_line = 0;
  _current_leaf_id = 0;
  _node_open = false;
}


/** コメントの削除
 *
 * @param [in] ss ファイルストリーム。
 * @param [in,out] buffer 入力文字列。
 * @param [out] answer コメントを削除したかどうか。
 * @return エラーコード
 *
 */
TextParserError TextParserTree::removeComment(std::stringstream& ss, std::string& buffer, bool& answer)
{
  answer = false;

  // 最初が "//" ?
  if (!buffer.compare(0, 2, "//")) {
    answer = true;
    buffer = "";
  } else if (!buffer.compare(0, 2, "/*")) {
    answer = true;
    while (true) {
      int e = buffer.find("*/");
      if (e >= 0) {
	buffer = buffer.substr(e + 2);
	break;
      } else {
	if (getline(ss, buffer)) {
	  _current_line++;
	} else {
	  return TextParserErrorHandler(TP_MISSING_COMMENT_END_ERROR, buffer);
	}
      }
    }
    buffer = TextParserRemoveHeadSpaces (buffer);                         // 先頭の空白削除
  }

  return TP_NO_ERROR;
}

/** ラベルの判定
 *
 * @param[in] ss 入力ファイルポインタ
 * @param[in,out] buffer 入力文字列
 * @param[out] answer ラベルか否か
 * @param[out] label ラベル
 * @return リターンコード
 *
 */
TextParserError TextParserTree::parseLabel(std::stringstream& ss, std::string& buffer, bool& answer, std::string& label)
{
#ifdef MYDEBUG
  std::cout << "parseLabel(std::stringstream& ss, std::string& buffer, bool& answer, std::string& label) start"<<std::endl;
  std::cout << "label" << label <<std::endl;
#endif

  TextParserError ret;
  answer=false;
 loop:
#ifdef MYDEBUG
    std::cout << "parseLabel 1\""<< buffer <<"\"" <<std::endl;
#endif //MYDEBUG
    ret = removeCommentsAndSpaces(ss, buffer);
#ifdef MYDEBUG
  std::cout << "parseLabel 2\""<< buffer <<"\"" <<std::endl;
#endif //MYDEBUG
  if (ret != TP_NO_ERROR) return ret;


  int pos_amp = buffer.find('&');
  int pos_pipe = buffer.find('|');
  int pos_lpar = buffer.find('(');

  if( pos_amp==0  || pos_pipe==0 || pos_lpar==0) return ret;
 

  // new version
  int pos_eq=buffer.find('=');
  // left curly brace
  int pos_lcb=buffer.find('{');
  // right curly brace
  int pos_rcb=buffer.find('}');
#ifdef MYDEBUG2
  std::cout << "parseLabel pos_eq "<< pos_eq 
	    << " pos_lcb"<< pos_lcb 
	    << " pos_rcb"<< pos_rcb
	    <<std::endl;
#endif // MYDEBUG2

  if ( ( pos_eq <0 ) && (pos_lcb < 0) && (pos_rcb < 0) ) {
    std::string tmp_buffer;
    if(!getline(ss,tmp_buffer)) goto end;
    _current_line++;
    buffer+=tmp_buffer;
#ifdef MYDEBUG2
    std::cout << "parseLabel tmp_buffer\""<< tmp_buffer <<"\"" 
	      <<" buffer \""<< buffer 
	      <<"\" _current_line "<< _current_line
	      <<std::endl;
    //      getline(std::cin,tmp_buffer);
    std::cout << "buffer in parseLabel \""<< buffer <<"\"" <<std::endl;
#endif // MYDEBUG2
    goto loop;
  } else {
    goto end;
  }
 end:     


  ret = parseLabel(buffer, answer, label);
  if (ret != TP_NO_ERROR) return ret;

#ifdef MYDEBUG
  std::cout << "parseLabel(std::stringstream& ss, std::string& buffer, bool& answer, std::string& label) end"<<std::endl;
#endif


  return TP_NO_ERROR;
}

/** ラベルの判定
 *
 * @param[in,out] buffer 入力文字列
 * @param[out] answer ラベルか否か
 * @param[out] label ラベル
 * @return リターンコード
 *
 */
TextParserError TextParserTree::parseLabel(std::string& buffer, bool& answer, std::string& label)
{

#ifdef MYDEBUG2
  std::cout<< "parseLabel(std::string& buffer, bool& answer, std::string& label) start" << std::endl;
  std::cout<< "buffer :" << buffer <<    std::endl;
  std::cout<< "label :" << label <<    std::endl;
  std::cout<< "answer :" << answer <<    std::endl;

#endif
  answer = false;

  buffer = TextParserRemoveHeadSpaces(buffer);

  int pos_amp = buffer.find('&');
  int pos_pipe = buffer.find('|');
  int pos_lpar = buffer.find('(');

  if( pos_amp==0 || pos_pipe==0 || pos_lpar==0) return TP_NO_ERROR;



#if 1 // new version
  // new version
  int pos_eq=buffer.find('=');
  // left curly brace
  int pos_lcb=buffer.find('{');
  int pos_exc=buffer.find('!');
  int pos=-100;
#ifdef MYDEBUG2    
  std::cout << "pos_eq " << pos_eq 
	    << " pos_lcb "<<pos_lcb 
	    << " pos_exc "<<pos_exc
	    <<std::endl;
#endif // MYDEBUG2    
    
  std::vector<int> tmp_pos_list;
  if(pos_eq >=0) {
    tmp_pos_list.push_back(pos_eq);
#ifdef MYDEBUG2    
    std::cout << "pos_eq label? " << buffer.substr(0,pos_eq) 
	      << " buffer? " << buffer.substr(pos_eq) << std::endl;
#endif // MYDEBUG2    
  }
  if(pos_lcb >=0) {
    tmp_pos_list.push_back(pos_lcb);
#ifdef MYDEBUG2    
    std::cout << "pos_lcb label? " << buffer.substr(0,pos_lcb) 
	      << " buffer? " << buffer.substr(pos_lcb) << std::endl;
#endif // MYDEBUG2    
  }
  if(pos_exc >=0) {
    tmp_pos_list.push_back(pos_exc);
#ifdef MYDEBUG2    
    std::cout << "pos_exc label? " << buffer.substr(0,pos_exc) 
	      << " buffer? " << buffer.substr(pos_exc) << std::endl;
#endif // MYDEBUG2    
  }
    
  if(tmp_pos_list.size()==0){
#ifdef MYDEBUG2    
    std::cout << "It may not be a label." << std::endl;
#endif // MYDEBUG2    
    return TP_NO_ERROR;
  } else {
    answer = true;
#ifdef MYDEBUG2    
    std::cout << "It may be a label." << std::endl;
#endif //MYDEBUG2
    pos=*( std::min_element( tmp_pos_list.begin(), tmp_pos_list.end() ) );
  }

#if 0      
  if (pos_eq < 0){
    pos = pos_lcb;
  } else if(pos_lcb < 0){
    pos = pos_eq;
  } else { // both pos_lcb and pos_eq are positive. take lower. 
    pos=pos_eq;
    if(pos_eq>pos_lcb)pos=pos_lcb;
  }
#endif
#ifdef MYDEBUG
  std::cout << " pos " << pos <<std::endl;	    
#endif // MYDEBUG2
  //std::string tmp_label=buffer.substr(0,pos-1);
  //TextParserRemoveHeadSpaces(tmp_label);

  label=buffer.substr(0,pos);
  buffer=buffer.substr(pos);
#ifdef MYDEBUG
  std::cout<< "label :\"" << label << "\""<<   std::endl;
  std::cout<< "buffer :\"" << buffer <<"\""   << std::endl;
#endif // MYDEBUG
  //    TextParserRemoveHeadSpaces(label);

  // remove double quote
  std::string tmp_label=label;
  int fdq=tmp_label.find('\"');

  //  std::cout<< "parseLabel fdq: "<<fdq<<std::endl;
  if(fdq>=0){
    tmp_label=tmp_label.substr(fdq+1);
    //    std::cout<< "parseLabel tmp_label 1:\""<<tmp_label<<std::endl;
    int sdq=tmp_label.find('\"');
    //    std::cout<< "parseLabel sdq: "<<sdq<<std::endl;
    if(sdq>=0){
      tmp_label=tmp_label.substr(0,sdq);
      //      std::cout<< "parseLabel tmp_label 2:\""<<tmp_label<<std::endl;
    } else {	
      TextParserErrorHandler(TP_ILLEGAL_TOKEN_ERROR,label);
      return TP_ILLEGAL_TOKEN_ERROR;
    }
  }
  //  std::cout<<"label :|"<<label
  //	   <<"| tmp_label |"<<tmp_label<<"|"<<std::endl;
  label=TextParserRemoveTailSpaces(tmp_label);


#if 0
  if(label[0]!='\"'){
    std::cout << "no  double quote" <<std::endl;
    label=TextParserRemoveTailSpaces(label);
  }else{
    std::cout << "has a pare of double quotes " <<std::endl;
    label=label.substr(1);
    int e = label.find("\""); // second double quote.
    label=label.substr(0,e);
  } 
#endif
  //  std::cout <<  "final label |" <<label<<"|"<<std::endl;

#endif // new version


#if 0 // old version
  // 最初が " ?
  if (buffer[0] != '\"') return TP_NO_ERROR;    // "から始まらないのでラベルではない

  buffer = buffer.substr(1);
  // ２つ目の"を探す
  int e = buffer.find("\"");

  if (e >= 0) {
    answer = true;
    label = buffer.substr(0, e);
    buffer = buffer.substr(e + 1);
  } else {
    return TP_NO_ERROR;                       // "で終わらないのでラベルではない
  }
#endif
    

#ifdef TP_DEBUG
  if (_debug_write) {
    std::cout << "ラベル : " + label << std::endl;
  }
#endif

#ifdef MYDEBUG2
  std::cout<< "parseLabel(std::string& buffer, bool& answer, std::string& label) end" << std::endl;
  std::cout<< "buffer :\"" << buffer <<"\""   << std::endl;
  std::cout<< "label :\"" << label << "\""<<   std::endl;
  std::cout<< "answer :\"" << answer << "\"" <<  std::endl;
#endif

  return TP_NO_ERROR;
}

/** ラベルが正常か判定
 *
 * @param[in] label ラベル
 * @param[in] pathlabel パスラベルか否か
 * @return ラベルが正常か否か
 *
 */
bool TextParserTree::isCorrectLabel(std::string& label, bool pathlabel)
{

#ifdef MYDDEEBUG2
  std::cout <<" isCorrectLabel(std::string& label, bool pathlabel) !" << label<<"| start"<<std::endl
#endif // MYDDEEBUG2

    std::string label0 = label;
  int a = label0.find("[@]");
  if (a >= 0) { // 配列タイプのラベル
    if (a == label0.size() - 3) {    // ラベルの最後以外に[@]が有る
      label0 = label.substr(0, a);
    } else {
#ifdef MYDDEEBUG2
      std::cout <<" isCorrectLabel(std::string& label, bool pathlabel) !" << label<<"| end 2 "<<std::endl
#endif // MYDDEEBUG2
	return true;
    }
  }
  for (unsigned int i = 0; i < label0.size(); i++) {
    if (label0[i] >= 'a' && label0[i] <= 'z') continue;
    if (label0[i] >= 'A' && label0[i] <= 'Z') continue;
    if (label0[i] >= '0' && label0[i] <= '9') continue;
    if (label0[i] == '_') continue;
    if (label0[i] == '-') continue;
    if (label0[i] == '[') continue;
    if (label0[i] == ']') continue;
    if (pathlabel) {    // パスラベル
      if (label0[i] == '.') continue;
      if (label0[i] == '/') continue;
    }
    // 使用できない文字
#ifdef MYDDEEBUG2
    std::cout <<" isCorrectLabel(std::string& label, bool pathlabel) !" << label<<"| end 1 "<<std::endl
#endif // MYDDEEBUG2
      return false;
  }

#ifdef MYDDEEBUG2
  std::cout <<" isCorrectLabel(std::string& label, bool pathlabel) !" << label<<"| end 0 "<<std::endl
#endif // MYDDEEBUG2
    return true;
}

/** ラベルがパスか判定
 *
 * @param[in] label ラベル
 * @return ラベルがパスか否か
 *
 * 一つでも'/'を含めばパスと認識
 * ".."から始まればパスと認識
 *
 */
bool TextParserTree::isPathLabel(const std::string& label)
{
  if (!label.compare(0, 2, "..")) {
    return true;
  }
  for (unsigned int i = 0; i < label.size(); i++) {
    if (label[i] == '/') {
      return true;
    }
  }

  return false;
}

/** ラベルがパス内に存在するか判定
 *
 * @param[in] label ラベル
 * @param[in] node ディレクトリ
 * @return ラベルがパス内に存在するか否か
 *
 * パス内に同じラベルのディレクトリが存在するかチェックする
 * 配列型のラベルの場合は[@]の前の部分が同一のラベルの存在もチェックする
 *
 */
bool TextParserTree::isIncludedInPath(const std::string& label, TextParserNode *node)
{
  int a = label.find("[@]");
  if (a < 0) return false;                    // 配列ラベルではない
  if (a != label.size() - 3) return false;    // [@]の後に余計な文字がある
  std::string label0 = label.substr(0, a);
  TextParserElement *element = node;
  while (element != 0) {
    std::string label1 = element->_label;
    int a = label1.find("[");
    if (a >= 0) { // 配列タイプのラベル
      label1 = label1.substr(0, a);
    }
    if (TextParserStringCompare(label0, label1)) return true;
    element = element->_parent;
  }

  return false;
}

/** 配列ラベルが同じエレメントが存在するか判定
 *
 * @param[in] label ラベル
 * @param[in] parent_element エレメント
 * @param[in] type エレメントのタイプ
 * @return 配列ラベルが同じエレメントが存在するか否か
 *
 */
bool TextParserTree::isArrayLabelExist(const std::string& label, TextParserElement *parent_element, const TextParserElementType type)
{
  int a = label.find("[@]");
  if (a < 0) return false;    // 配列ラベルではない
  if (a != label.size() - 3) return false;    // [@]の後に余計な文字がある
  std::string label0 = label.substr(0, a);
  if (parent_element == 0) {  // ルートディレクトリ
    if (type == TP_NODE_ELEMENT) {
      std::map<std::string, TextParserNode *>::iterator di = _nodes.begin();
      while (di != _nodes.end()) {
	TextParserNode *node = di->second;
	std::string label1 = node->_label;
	int a = label1.find("[");
	if (a >= 0) { // 配列タイプのラベル
	  label1 = label1.substr(0, a);
	}
	if (TextParserStringCompare(label0, label1)) return true;
	di++;
      }
    }
  } else {
    if (type == TP_NODE_ELEMENT) {
      TextParserNode *parent_dir = (TextParserNode *)parent_element;
      std::map<std::string, TextParserNode *>::iterator di = parent_dir->_nodes.begin();
      while (di != parent_dir->_nodes.end()) {
	TextParserNode *node = di->second;
	std::string label1 = node->_label;
	int a = label1.find("[");
	if (a >= 0) { // 配列タイプのラベル
	  label1 = label1.substr(0, a);
	}
	if (TextParserStringCompare(label0, label1)) return true;
	di++;
      }
    } else if (type == TP_LEAF_ELEMENT) {
      TextParserNode *parent_dir = (TextParserNode *)parent_element;
      std::map<std::string, TextParserLeaf *>::iterator li = parent_dir->_leaves.begin();
      while (li != parent_dir->_leaves.end()) {
	TextParserLeaf *leaf = li->second;
	std::string label1 = leaf->_label;
	int a = label1.find("[");
	if (a >= 0) { // 配列タイプのラベル
	  label1 = label1.substr(0, a);
	}
	if (TextParserStringCompare(label0, label1)) return true;
	li++;
      }
    }
  }

  return false;
}

/** 配列ラベルのインデックス設定
 *
 * @param[in,out] label ラベル
 * @return リターンコード
 *
 * 配列ラベル（末尾が[@]）だったらインデックスを割り当てる
 *
 */
TextParserError TextParserTree::setArrayLabelIndex(std::string& label)
{
  return SetArrayLabelIndex(label, _array_label_number);
}

/** ディレクトリの判定
 *
 * @param[in] ss 入力ファイルポインタ
 * @param[in,out] buffer 入力文字列
 * @param[out] answer ディレクトリか否か
 * @return リターンコード
 *
 */
TextParserError TextParserTree::parseNode(std::stringstream& ss, std::string& buffer, bool& answer)
{
  TextParserError ret;
  answer = false;

  ret = removeCommentsAndSpaces(ss, buffer);
  if (ret != TP_NO_ERROR) return ret;

  // 最初が { ?
  if (buffer[0] != '{') return TP_NO_ERROR;     // {から始まらないのでディレクトリではない

  answer = true;
  buffer = buffer.substr(1);

#ifdef TP_DEBUG
  /*
    if (_debug_write) {
    std::cout << "<ディレクトリ>" << std::endl;
    }
  */
#endif
  return TP_NO_ERROR;
}

/** ディレクトリの終わり判定
 *
 * @param[in] ss 入力ファイルポインタ
 * @param[in,out] buffer 入力文字列
 * @param[out] answer ディレクトリの終わりか否か
 * @return リターンコード
 *
 */
TextParserError TextParserTree::parseEndOfNode(std::stringstream& ss, std::string& buffer, bool& answer)
{
  TextParserError ret;
  answer = false;

  ret = removeCommentsAndSpaces(ss, buffer);
  if (ret != TP_NO_ERROR) return ret;

  // 最初が } ?
  if (buffer[0] != '}') return TP_NO_ERROR;     // }から始まらないのでディレクトリの終わりではない

  answer = true;
  buffer = buffer.substr(1);

#ifdef TP_DEBUG
  /*
    if (_debug_write) {
    std::cout << "<ディレクトリの終わり>" << std::endl;
    }
  */
#endif
  return TP_NO_ERROR;
}

/** ディレクトリの開始
 *
 * @param[in] label ラベルパス
 * @return リターンコード
 *
 * カレントエレメントに存在しなければディレクトリを作成し、それを
 * カレントエレメントにする
 *
 */
TextParserError TextParserTree::openNode(const std::string& label)
{
#ifdef MYDEBUG
  std::cout<<"openNode(const std::string& label) start" <<std::endl;
  std::cout<<"label :"<<label <<std::endl;
#endif
  TextParserError ret;
  TextParserNode *node = 0;
  TextParserElement *parent_element = 0;
  std::string path = label;

  if(label.size()==0){
    // そもそもラベルのサイズが0なのはエラー
    return TextParserErrorHandler(TP_ILLEGAL_LABEL_ERROR,label);
  }

  ret = getElementRelativePath(path, true, &parent_element);

  if (ret != TP_NO_ERROR) return ret;
  if (path.size() == 0) {
    node->_return_elements.push_back(_current_element);    // カレントエレメントを戻り先にスタック
    _current_element = parent_element;                               // カレントエレメントに設定
  } else {
    // ラベルの異常検出
    if (!isCorrectLabel(path, false)) {
      std::cout << "test "<< std::endl;
      return TextParserErrorHandler(TP_ILLEGAL_LABEL_ERROR, label);
    }
    ret = getOrAddNode(path, parent_element, &node);   // ディレクトリの取得又は追加
    if (ret != TP_NO_ERROR) return ret;
    node->_return_elements.push_back(_current_element);    // カレントエレメントを戻り先にスタック
    _current_element = node;                               // カレントエレメントに設定
    node->setLineN(_current_line);
  }
  _node_open = true;
    
#ifdef TP_DEBUG
  if (_debug_write) {
    if (_current_element == 0) {
      std::cout << "ディレクトリの開始 : /" << std::endl;
    } else {
      std::cout << "ディレクトリの開始 : " + _current_element->_label << std::endl;
    }
  }
#endif
    
#ifdef MYDEBUG
  std::cout<<"openNode(const std::string& label) end"
	   <<std::endl;
#endif
  return TP_NO_ERROR;
}

/** ノードの終了
 *
 * @return リターンコード
 *
 * パースモードを戻し、ディレクトリの親をカレントエレメントにする
 *
 */
TextParserError TextParserTree::closeNode()
{
  TextParserNode *current_dir = (TextParserNode *)_current_element;
  if (current_dir != 0) {
    _current_element = current_dir->_return_elements.back();    // 最新の戻り先エレメントをカレントエレメントに設定
    current_dir->_return_elements.pop_back();                   // 最新の戻り先エレメントを削除
    if (_current_element != 0 && _current_element->_type != TP_NODE_ELEMENT) {
      return TextParserErrorHandler(TP_ILLEGAL_CURRENT_ELEMENT_ERROR, "");
    }
    _parse_mode = TP_NODE_PARSE;                         // ディレクトリ解析モード
  } else {
    return TextParserErrorHandler(TP_NODE_END_ERROR, "");
  }
  _node_open = false;

#ifdef TP_DEBUG
  if (_debug_write) {
    std::cout << "ディレクトリの終了 : " + current_dir->_label << std::endl;
  }
#endif
  return TP_NO_ERROR;
}

/** ディレクトリの取得又は追加
 *
 * @param[in] label ラベル
 * @param[in] parent_element 親エレメント
 * @param[out] node ディレクトリ
 * @return リターンコード
 *
 * ラベルに一致するディレクトリが有れば取得し、無ければ追加する
 *
 */
TextParserError TextParserTree::getOrAddNode(std::string& label, TextParserElement *parent_element, TextParserNode **node)
{
  TextParserError ret;
  *node = 0;
  if (parent_element == 0) {							// 親エレメントがルートディレクトリ
    if (isArrayLabelExist(label, parent_element, TP_LEAF_ELEMENT)) {
      return TextParserErrorHandler(TP_LABEL_ALREADY_USED_ERROR, label);
    }
    ret = setArrayLabelIndex(label);                // 配列型ラベルのインデックス設定
    if (ret != TP_NO_ERROR) {
      return TextParserErrorHandler(ret, label);
    }
    *node = getNode(label);               // 既存のディレクトリを取得
    if (*node == 0) {   // 存在しない
      try {
	*node = new TextParserNode(label);      // 新規にディレクトリを作成
	if (*node == 0) {
	  return TextParserErrorHandler(TP_MEMORY_ALLOCATION_ERROR, label);
	}
      } catch (std::exception ex) {
	return TextParserErrorHandler(TP_MEMORY_ALLOCATION_ERROR, label);
      }
      addElement(*node);                     // ルートディレクトリに追加
      (*node)->_parent = parent_element;		// 親エレメントを設定
#ifdef TP_DEBUG
      if (_debug_write) {
	std::cout << "★ ディレクトリの作成 : " + label << std::endl;
      }
#endif
    }
  } else if (parent_element->_type == TP_NODE_ELEMENT) {
    TextParserNode *parent_dir = (TextParserNode *)parent_element;
    if (isIncludedInPath(label, parent_dir)) {
      return TextParserErrorHandler(TP_LABEL_ALREADY_USED_PATH_ERROR, label);
    }
    if (isArrayLabelExist(label, parent_element, TP_LEAF_ELEMENT)) {
      return TextParserErrorHandler(TP_LABEL_ALREADY_USED_ERROR, label);
    }
    ret = parent_dir->setArrayLabelIndex(label);    // 配列型ラベルのインデックス設定
    if (ret != TP_NO_ERROR) {
      return TextParserErrorHandler(ret, label);
    }
    *node = parent_dir->getNode(label);	// 既存のディレクトリを取得
    if (*node == 0) {   // 存在しない
      // 同一ラベルのリーフを探す
      TextParserLeaf *leaf = parent_dir->getLeaf(label);	// 既存のリーフを取得
      if (leaf != 0) {                            // 同一名のリーフがあればエラー
	return TextParserErrorHandler(TP_LABEL_ALREADY_USED_ERROR, label);
      }
      try {
	*node = new TextParserNode(label);      // 新規にディレクトリを作成
	if (*node == 0) {
	  return TextParserErrorHandler(TP_MEMORY_ALLOCATION_ERROR, label);
	}
      } catch (std::exception ex) {
	return TextParserErrorHandler(TP_MEMORY_ALLOCATION_ERROR, label);
      }
      parent_dir->addElement(*node);			// 親ディレクトリに追加
      (*node)->_parent = parent_element;		// 親エレメントを設定
#ifdef TP_DEBUG
      if (_debug_write) {
	std::cout << "★ ディレクトリの作成 : " + label << std::endl;
      }
#endif
    }
  } else {
    return TextParserErrorHandler(TP_ILLEGAL_CURRENT_ELEMENT_ERROR, "");
  }

  return TP_NO_ERROR;
}

/** リーフの判定
 *
 * @param[in] ss 入力ファイルポインタ
 * @param[in,out] buffer 入力文字列
 * @param[out] answer リーフか否か
 * @return リターンコード
 *
 */
TextParserError TextParserTree::parseLeaf(std::stringstream& ss, std::string& buffer, bool& answer)
{
  TextParserError ret;
  answer = false;

  ret = removeCommentsAndSpaces(ss, buffer);
  if (ret != TP_NO_ERROR) return ret;

  // 最初が = ?
  if (buffer[0] != '=') return TP_NO_ERROR; // =から始まらないのでリーフではない
  // =から始まらないのでディレクトリではない

  answer = true;
  buffer = buffer.substr(1);

#ifdef TP_DEBUG
  /*
    if (_debug_write) {
    std::cout << "<リーフ>" << std::endl;
    }
  */
#endif
  return TP_NO_ERROR;
}

/** リーフの開始
 *
 * @param[in] label ラベル
 * @return リターンコード
 *
 * カレントエレメントに存在しなければリーフを作成し、それを
 * カレントエレメントにする
 *
 */
TextParserError TextParserTree::openLeaf(const std::string& label)
{
#ifdef MYDEBUG
  std::cout << "openLeaf(const std::string& label) start"<<  std::endl;
  std::cout << "label :"<<label <<std::endl;
#endif
  TextParserError ret = TP_NO_ERROR;
  TextParserElement *parent_element = 0;
  std::string path = label;
  ret = getElementRelativePath(path, true, &parent_element);
  if (ret != TP_NO_ERROR) return ret;
  if (path.size() == 0) {
    return TextParserErrorHandler(TP_ILLEGAL_LABEL_PATH_ERROR, label);
  }
  TextParserLeaf *leaf = 0;
  if (parent_element != 0 && parent_element->_type == TP_NODE_ELEMENT) {
    // ラベルの異常検出
    if (!isCorrectLabel(path, false)) {
      return TextParserErrorHandler(TP_ILLEGAL_LABEL_ERROR, label);
    }
    TextParserNode *parent_dir = (TextParserNode *)parent_element;
    if (isIncludedInPath(path, parent_dir)) {
      return TextParserErrorHandler(TP_LABEL_ALREADY_USED_PATH_ERROR, path);
    }
    if (isArrayLabelExist(path, parent_element, TP_NODE_ELEMENT)) {
      return TextParserErrorHandler(TP_LABEL_ALREADY_USED_ERROR, path);
    }
    ret = parent_dir->setArrayLabelIndex(path);// 配列型ラベルのインデックス設定
    if (ret != TP_NO_ERROR) {
      return TextParserErrorHandler(ret, path);
    }
    leaf = parent_dir->getLeaf(path);          // 既存のリーフを取得
    if (leaf == 0) {    // 存在しない
      // 同一ラベルのディレクトリを探す
      TextParserNode *node = parent_dir->getNode(path);	// 既存のディレクトリを取得
      if (node != 0) {                            // 同一名のリーフがあればエラー
	return TextParserErrorHandler(TP_LABEL_ALREADY_USED_ERROR, path);
      }
      try {
	leaf = new TextParserLeaf(path);             // 新規にリーフを作成
	if (leaf == 0) {
	  return TextParserErrorHandler(TP_MEMORY_ALLOCATION_ERROR, "leaf");
	}
      } catch (std::exception ex) {
	return TextParserErrorHandler(TP_MEMORY_ALLOCATION_ERROR, "leaf");
      }
      parent_dir->addElement(leaf);               // 親ディレクトリに追加
      leaf->_parent = parent_element;             // 親エレメントを設定
      leaf->setLineN(_current_line);
    } else {
      return TextParserErrorHandler(TP_LABEL_ALREADY_USED_ERROR, path);
    }
  } else {
    return TextParserErrorHandler(TP_ILLEGAL_CURRENT_ELEMENT_ERROR, "");
  }
  leaf->_return_elements.push_back(_current_element);         // カレントエレメントを戻り先にスタック
  _current_element = leaf;                                    // カレントエレメントに設定

#ifdef TP_DEBUG
  if (_debug_write) {
    std::cout << "★ リーフの作成 : " + leaf->_label << std::endl;
  }
#endif //TP_DEBUG

#ifdef MYDEBUG
  std::cout << "openLeaf(const std::string& label) end"<<  std::endl;
  std::cout << "label :"<<label <<std::endl;
#endif //MYDEBUG


  return ret;
}

/** リーフの終了
 *
 * @return リターンコード
 *
 * パースモードを戻し、リーフの親をカレントエレメントにする
 *
 */
TextParserError TextParserTree::closeLeaf()
{
#ifdef MYDEBUG
  std::cout<<  "closeLeaf() start"<<std::endl;
#endif //MYDEBUG
  TextParserLeaf *current_leaf = (TextParserLeaf *)_current_element;
  if (current_leaf != 0) {
    _current_element = current_leaf->_return_elements.back();  // 最新の戻り先エレメントをカレントエレメントに設定
    current_leaf->_return_elements.pop_back();                  // 最新の戻り先エレメントを削除
    if (_current_element == 0) {                                // ルートディレクトリ    
      _parse_mode = TP_NODE_PARSE;                     // ディレクトリ解析モード
    } else if (_current_element->_type == TP_NODE_ELEMENT) {  // ディレクトリ    
      _parse_mode = TP_NODE_PARSE;                     // ディレクトリ解析モード
    } else {
      return TextParserErrorHandler(TP_ILLEGAL_CURRENT_ELEMENT_ERROR, "");
    }
  } else {
    return TextParserErrorHandler(TP_ILLEGAL_CURRENT_ELEMENT_ERROR, "");
  }

#ifdef TP_DEBUG
  if (_debug_write) {
    std::cout << "リーフの終了 : " + current_leaf->_label << std::endl;
  }
#endif //TP_DEBUG
#ifdef MYDEBUG
  std::cout<<  "closeLeaf() end"<<std::endl;
#endif //MYDEBUG
  return TP_NO_ERROR;
}

/** 依存関係付き値の判定
 *
 * @param[in] ss 入力ファイルポインタ
 * @param[in,out] buffer 入力文字列
 * @param[out] answer 依存関係付き否か
 * @return リターンコード
 *
 */
TextParserError TextParserTree::parseDependValue(std::stringstream& ss, std::string& buffer, bool& answer)
{
  TextParserError ret;
  answer = false;

  ret = removeCommentsAndSpaces(ss, buffer);
  if (ret != TP_NO_ERROR) return ret;

  // 最初が @dep ?
  if (!TextParserStringCompare(buffer.substr(0, 4), "@dep")) return TP_NO_ERROR;  // @depから始まらないので依存関係付き値ではない 

  answer = true;
  buffer = buffer.substr(4);

#ifdef TP_DEBUG
  if (_debug_write) {
    std::cout << "<依存関係付き値>" << std::endl;
  }
#endif //TP_DEBUG
  return TP_NO_ERROR;
}

/** ベクトル値の判定
 *
 * @param[in] ss 入力ファイルポインタ
 * @param[in,out] buffer 入力文字列
 * @param[out] answer ベクトル値か否か
 * @return リターンコード
 *
 */
TextParserError TextParserTree::parseVectorValue(std::stringstream& ss, std::string& buffer, bool& answer)
{
  TextParserError ret;

  ret = removeCommentsAndSpaces(ss, buffer);
  if (ret != TP_NO_ERROR) return ret;

  ret = parseVectorValue(buffer, answer);
  if (ret != TP_NO_ERROR) return ret;

  return TP_NO_ERROR;
}

/** ベクトル値の判定
 *
 * @param[in,out] buffer 入力文字列
 * @param[out] answer ベクトル値か否か
 * @return リターンコード
 *
 */
TextParserError TextParserTree::parseVectorValue(std::string& buffer, bool& answer)
{
  answer = false;

  // 最初が ( ?
  if (buffer[0] != '(') return TP_NO_ERROR;     // (から始まらないのでベクトル値ではない

  answer = true;
  buffer = buffer.substr(1);

#ifdef TP_DEBUG
  if (_debug_write) {
    std::cout << "<ベクトル値>" << std::endl;
  }
#endif
  return TP_NO_ERROR;
}

/** 未定義値の判定
 *
 * @param[in] ss 入力ファイルポインタ
 * @param[in,out] buffer 入力文字列
 * @param[out] answer 未定義値か否か
 * @return リターンコード
 *
 */
TextParserError TextParserTree::parseUndefinedValue(std::stringstream& ss, std::string& buffer, bool& answer)
{
  TextParserError ret;

  ret = removeCommentsAndSpaces(ss, buffer);
  if (ret != TP_NO_ERROR) return ret;

  ret = parseUndefinedValue(buffer, answer);
  if (ret != TP_NO_ERROR) return ret;

  return TP_NO_ERROR;
}

/** 未定義値の判定
 *
 * @param[in,out] buffer 入力文字列
 * @param[out] answer 未定義値か否か
 * @return リターンコード
 *
 */
TextParserError TextParserTree::parseUndefinedValue(std::string& buffer, bool& answer)
{
  answer = false;

  // 最初が UNDEF ?
  if (buffer.compare(0, 5, "UNDEF")) return TP_NO_ERROR;  // UNDEFから始まらないので未定義値ではない

  if (buffer.size() > 5) {
    // 次に続く文字が値の区切りか？ UNDEFの後はスペースも可
    bool found = false;
    for (unsigned int i = 0; i < TextParserNumValueNext; i++) {
      if (!buffer.compare(5, 1 , TextParserValueNext[i]) || !buffer.compare(5, 1 , " ")) {
	found = true;
	break;
      }
    }
    if (!found) return TP_NO_ERROR;     // UNDEFの後に区切り以外の文字があるので未定義値ではない
  }
  answer = true;
  buffer = buffer.substr(5);

#ifdef TP_DEBUG
  if (_debug_write) {
    std::cout << "<未定義値>" << std::endl;
  }
#endif
  return TP_NO_ERROR;
}

/** 値の判定
 *
 * @param[in] ss 入力ファイルポインタ
 * @param[in,out] buffer 入力文字列
 * @param[out] answer 値か否か
 * @param[out] value 値
 * @param[out] value_type 値のタイプ
 * @return リターンコード
 *
 */
TextParserError TextParserTree::parseValue(std::stringstream& ss,
					   std::string& buffer,
					   bool& answer,
					   std::string& value,
					   TextParserValueType& value_type)
{
  TextParserError ret;

  ret = removeCommentsAndSpaces(ss, buffer);
  if (ret != TP_NO_ERROR) return ret;

  ret = parseValue(buffer, answer, value, value_type);
  if (ret != TP_NO_ERROR) return ret;

  return TP_NO_ERROR;
}

/** 値の判定
 *
 * @param[in,out] buffer 入力文字列
 * @param[out] answer 値か否か
 * @param[out] value 値
 * @param[out] value_type 値のタイプ
 * @return リターンコード
 *
 */
TextParserError TextParserTree::parseValue(std::string& buffer,
					   bool& answer,
					   std::string& value,
					   TextParserValueType& value_type)
{
#ifdef MYDEBUG3
  std::cout << "parseValue start buffer|"<<buffer<<"|"<<std::endl;
#endif //MYDEBUG3

  answer = false;

  value_type = TP_UNDEFFINED_VALUE;
  // 最初が " ?
  if (buffer[0] == '\"') {    // 値が文字列
    buffer = buffer.substr(1);
    int e = buffer.find("\"");
    if (e >= 0) {
      value = buffer.substr(0, e);
      buffer = buffer.substr(e + 1);
      value_type = TP_STRING_VALUE;
    } else {
      return TP_NO_ERROR;                               // "で終わらないので文字列ではない
    }
  } else {                    // 値が数値
    // 数値の区切りで切り取る
    // 最初は'+' '-' '.' '0'-'9'

    if (buffer[0] != '+' && buffer[0] != '-' && buffer[0] != '.'
        && (buffer[0] < '0' || buffer[0] > '9')) {
      return TP_NO_ERROR;
    }
		
    // 値の区切りの検出
    unsigned int end = buffer.size();

    for (unsigned int i = 0; i < TextParserNumValueNext; i++) {
      unsigned int e = buffer.find(charTextParserValueNext[i]);
      //      unsigned int e = buffer.find(TextParserValueNext[i]);
      if (e > 0 && e < end) end = e;
    }

    // special case 
    // ダブルクオート無しのラベルにも対応する。カンマは直前でチェックされているので、
    // '=' の前の スペース区切りを検出する。
    int pos_eq = buffer.find('=');
    if (pos_eq != std::string::npos){ // 次のラベルが後ろに続いている。
      std::string candidate=buffer.substr(0,pos_eq);
      candidate = TextParserRemoveTailSpaces(candidate); // 末尾の空白を削除
      int lspace=candidate.rfind(' ');
#ifdef MYDEBUG3
      std::cout << "parseValue: have next label on back in the line."
		<< std::endl;
      std::cout << "pos_eq "<<pos_eq << std::endl;
      std::cout << "end "<<end << std::endl;
      std::cout << "candidate |"<<candidate<<"|"<<std::endl;
      std::cout << "lspace "<<lspace << std::endl;
#endif //MYDEBUG3
      if(lspace!=std::string::npos){
	if(end>lspace) end=lspace;
       }
#ifdef MYDEBUG3      
      std::cout<< "test value  |"<< buffer.substr(0,end)<<"|"<<std::endl;
      std::cout<< "test buffer |"<< buffer.substr(end) <<"|"<<std::endl;
#endif //MYDEBUG3      
    }

    if (end < buffer.size()) {
      value = buffer.substr(0, end);
      buffer = buffer.substr(end);                        // ベクトルの','や')'、ラベルの'\"'を先頭に残すため検出した文字を削除しない
    } else {
      value = buffer;
      buffer = "";
    }
    value = TextParserRemoveTailSpaces(value);                    // 末尾の空白を削除
		
    // 数値のチェック
    // 許される表現 (+/-) 0-9(.(0-9)) (e/d (+/-) 0-9)
    std::string nvalue = value;
    if (nvalue[0] == '+' || nvalue[0] == '-') {             // 先頭が+-
      nvalue = nvalue.substr(1);                          // +-を削除
      nvalue = TextParserRemoveHeadSpaces(nvalue);              // 先頭の空白を削除
    }
    int e = nvalue.find("e");
    if (e < 0) e = nvalue.find("d");
    if (e < 0) e = nvalue.find("E");
    if (e < 0) e = nvalue.find("D");
    if (e >= 0) {   // 指数？
      std::string exponent = nvalue.substr(e + 1);        // 指数部分
      exponent = TextParserRemoveHeadSpaces(exponent);          // 先頭の空白を削除
      if (exponent[0] == '+' || exponent[0] == '-') {     // 先頭が+-
	exponent = exponent.substr(1);                  // +-を削除
	exponent = TextParserRemoveHeadSpaces(exponent);      // 先頭の空白を削除
      }
      for (unsigned int i = 0; i < exponent.size(); i++) {
	if (exponent[i] < '0' || exponent[i] > '9') {   // 指数部に数値以外がある
	  return TextParserErrorHandler(TP_ILLEGAL_NUMERIC_VALUE_ERROR, value);
	}
      }
      nvalue = nvalue.substr(0, e);                       // 仮数部分
      nvalue = TextParserRemoveTailSpaces(nvalue);              // 末尾の空白を削除
    }

    for (unsigned int i = 0; i < nvalue.size(); i++) {
      if ((nvalue[i] < '0' || nvalue[i] > '9') && nvalue[i] != '.') {   // 仮数部に数値と.以外がある
#ifdef MYDEBUG3
  std::cout << "parseValue 0 buffer=|"<<buffer
	    << "| value=|"<<value<<"|"<<std::endl;
#endif //MYDEBUG3

	return TextParserErrorHandler(TP_ILLEGAL_NUMERIC_VALUE_ERROR, value);
      }
    }
    // 数値変換チェック
    double dval;
#if 0
    try {
      std::stringstream ss;
      ss << value;
      ss >> dval;
    } catch (std::exception ex) {
      return TextParserErrorHandler(TP_ILLEGAL_NUMERIC_VALUE_ERROR, value);
    }
#else
    try {
      sscanf(value.c_str(), "lf", &dval);
    } catch (std::exception ex) {

#ifdef MYDEBUG3
  std::cout << "parseValue 1 buffer=|"<<buffer
	    << "| value=|"<<value<<"|"<<std::endl;
#endif //MYDEBUG3

      return TextParserErrorHandler(TP_ILLEGAL_NUMERIC_VALUE_ERROR, value);
    }
#endif
		
    value_type = TP_NUMERIC_VALUE;
  }


#ifdef MYDEBUG3
  std::cout << "parseValue end buffer=|"<<buffer
	    << "| value=|"<<value<<"|"<<std::endl;
#endif //MYDEBUG3



  answer = true;

#ifdef TP_DEBUG
  if (_debug_write) {
    if (value_type == TP_STRING_VALUE) {
      std::cout << "値（文字列） : " + value << std::endl;
    } else if (value_type == TP_NUMERIC_VALUE) {
      std::cout << "値（数値） : " + value << std::endl;
    }
  }
#endif //TP_DEBUG
  return TP_NO_ERROR;
}

/** 右括弧の判定
 *
 * @param[in] ss 入力ファイルポインタ
 * @param[in,out] buffer 入力文字列
 * @param[out] answer 右括弧か否か
 * @return リターンコード
 *
 */
TextParserError TextParserTree::parseOpenBrancket(std::stringstream& ss, std::string& buffer, bool& answer)
{
  TextParserError ret;

  ret = removeCommentsAndSpaces(ss, buffer);
  if (ret != TP_NO_ERROR) return ret;

  ret = parseOpenBrancket(buffer, answer);
  if (ret != TP_NO_ERROR) return ret;

  return TP_NO_ERROR;
}
    
/** 右括弧の判定
 *
 * @param[in,out] buffer 入力文字列
 * @param[out] answer 右括弧か否か
 * @return リターンコード
 *
 */
TextParserError TextParserTree::parseOpenBrancket(std::string& buffer, bool& answer)
{
  answer = false;

  // 最初が ( ?
  if (buffer[0] != '(') return TP_NO_ERROR; 

  answer = true;
  buffer = buffer.substr(1);

#ifdef TP_DEBUG
  if (_debug_write) {
    std::cout << "(" << std::endl;
  }
#endif //TP_DEBUG
  return TP_NO_ERROR;
}
    
/** 左括弧の判定
 *
 * @param[in] ss 入力ファイルポインタ
 * @param[in,out] buffer 入力文字列
 * @param[out] answer 左括弧か否か
 * @return リターンコード
 *
 */
TextParserError TextParserTree::parseClosedBrancket(std::stringstream& ss, std::string& buffer, bool& answer)
{
  TextParserError ret;

  ret = removeCommentsAndSpaces(ss, buffer);
  if (ret != TP_NO_ERROR) return ret;

  ret = parseClosedBrancket(buffer, answer);
  if (ret != TP_NO_ERROR) return ret;

  return TP_NO_ERROR;
}
    
/** 左括弧の判定
 *
 * @param[in,out] buffer 入力文字列
 * @param[out] answer 左括弧か否か
 * @return リターンコード
 *
 */
TextParserError TextParserTree::parseClosedBrancket(std::string& buffer, bool& answer)
{
  answer = false;

  // 最初が ) ?
  if (buffer[0] != ')') return TP_NO_ERROR; 

  answer = true;
  buffer = buffer.substr(1);

#ifdef TP_DEBUG
  if (_debug_write) {
    std::cout << ")" << std::endl;
  }
#endif //TP_DEBUG
  return TP_NO_ERROR;
}
    
/** ==の判定
 *
 * @param[in] ss 入力ファイルポインタ
 * @param[in,out] buffer 入力文字列
 * @param[out] answer ==か否か
 * @return リターンコード
 *
 */
TextParserError TextParserTree::parseEqual(std::stringstream& ss, std::string& buffer, bool& answer)
{
  TextParserError ret;

  ret = removeCommentsAndSpaces(ss, buffer);
  if (ret != TP_NO_ERROR) return ret;

  ret = parseEqual(buffer, answer);
  if (ret != TP_NO_ERROR) return ret;

  return TP_NO_ERROR;
}
    
/** ==の判定
 *
 * @param[in,out] buffer 入力文字列
 * @param[out] answer ==か否か
 * @return リターンコード
 *
 */
TextParserError TextParserTree::parseEqual(std::string& buffer, bool& answer)
{
  answer = false;

  // 最初が == ?
  if (buffer.compare(0, 2, "==")) return TP_NO_ERROR; 

  answer = true;
  buffer = buffer.substr(2);

#ifdef TP_DEBUG
    if (_debug_write) {
    std::cout << "==" << std::endl;
    }
#endif //TP_DEBUG
  return TP_NO_ERROR;
}
    
/** !=の判定
 *
 * @param[in] ss 入力ファイルポインタ
 * @param[in,out] buffer 入力文字列
 * @param[out] answer !=か否か
 * @return リターンコード
 *
 */
TextParserError TextParserTree::parseNotEqual(std::stringstream& ss, std::string& buffer, bool& answer)
{
  TextParserError ret;

  ret = removeCommentsAndSpaces(ss, buffer);
  if (ret != TP_NO_ERROR) return ret;

  ret = parseNotEqual(buffer, answer);
  if (ret != TP_NO_ERROR) return ret;

  return TP_NO_ERROR;
}
    
/** !=の判定
 *
 * @param[in,out] buffer 入力文字列
 * @param[out] answer !=か否か
 * @return リターンコード
 *
 */
TextParserError TextParserTree::parseNotEqual(std::string& buffer, bool& answer)
{
  answer = false;

  // 最初が != ?
  if (buffer.compare(0, 2, "!=")) return TP_NO_ERROR; 

  answer = true;
  buffer = buffer.substr(2);

#ifdef TP_DEBUG
  /*
    if (_debug_write) {
    std::cout << "!=" << std::endl;
    }
  */
#endif //TP_DEBUG
  return TP_NO_ERROR;
}
    
/** &&の判定
 *
 * @param[in] ss 入力ファイルポインタ
 * @param[in,out] buffer 入力文字列
 * @param[out] answer &&か否か
 * @return リターンコード
 *
 */
TextParserError TextParserTree::parseAnd(std::stringstream& ss, std::string& buffer, bool& answer)
{
  TextParserError ret;

  ret = removeCommentsAndSpaces(ss, buffer);
  if (ret != TP_NO_ERROR) return ret;

  ret = parseAnd(buffer, answer);
  if (ret != TP_NO_ERROR) return ret;

  return TP_NO_ERROR;
}
    
/** &&の判定
 *
 * @param[in,out] buffer 入力文字列
 * @param[out] answer &&か否か
 * @return リターンコード
 *
 */
TextParserError TextParserTree::parseAnd(std::string& buffer, bool& answer)
{
  answer = false;
#ifdef MYDEBUG3
  std::cout<< "parseAnd |"<<buffer<<std::endl;
#endif // MYDEBUG3

  // 最初が && ?
  if (buffer.compare(0, 2, "&&")) return TP_NO_ERROR; 

  answer = true;
  buffer = buffer.substr(2);

#ifdef TP_DEBUG
  if (_debug_write) {
    std::cout << "&&" << std::endl;
  }
#endif //TP_DEBUG
  return TP_NO_ERROR;
}
    
/** ||の判定
 *
 * @param[in] ss 入力ファイルポインタ
 * @param[in,out] buffer 入力文字列
 * @param[out] answer ||か否か
 * @return リターンコード
 *
 */
TextParserError TextParserTree::parseOr(std::stringstream& ss, std::string& buffer, bool& answer)
{
  TextParserError ret;

  ret = removeCommentsAndSpaces(ss, buffer);
  if (ret != TP_NO_ERROR) return ret;

  ret = parseOr(buffer, answer);
  if (ret != TP_NO_ERROR) return ret;

  return TP_NO_ERROR;
}

/** ||の判定
 *
 * @param[in,out] buffer 入力文字列
 * @param[out] answer ||か否か
 * @return リターンコード
 *
 */
TextParserError TextParserTree::parseOr(std::string& buffer, bool& answer)
{
  answer = false;

  // 最初が || ?
  if (buffer.compare(0, 2, "||")) return TP_NO_ERROR; 

  answer = true;
  buffer = buffer.substr(2);

#ifdef TP_DEBUG
  if (_debug_write) {
    std::cout << "||" << std::endl;
  }
#endif //TP_DEBUG
  return TP_NO_ERROR;
}

/** ?の判定
 *
 * @param[in] ss 入力ファイルポインタ
 * @param[in,out] buffer 入力文字列
 * @param[out] answer ?か否か
 * @return リターンコード
 *
 */
TextParserError TextParserTree::parseQuestion(std::stringstream& ss, std::string& buffer, bool& answer)
{
  TextParserError ret;

  ret = removeCommentsAndSpaces(ss, buffer);
  if (ret != TP_NO_ERROR) return ret;

  ret = parseQuestion(buffer, answer);
  if (ret != TP_NO_ERROR) return ret;

  return TP_NO_ERROR;
}
    
/** ?の判定
 *
 * @param[in,out] buffer 入力文字列
 * @param[out] answer ?か否か
 * @return リターンコード
 *
 */
TextParserError TextParserTree::parseQuestion(std::string& buffer, bool& answer)
{
  answer = false;

#ifdef MYDEBUG3    
  std::cout << "@@@ parseQuestion buffer="<<buffer<<std::endl;
#endif //MYDEBUG3    

  // 最初が ? ?
  if (buffer[0] != '?') return TP_NO_ERROR; 

  answer = true;
  buffer = buffer.substr(1);

#ifdef TP_DEBUG
  if (_debug_write) {
    std::cout << "?" << std::endl;
  }
#endif //TP_DEBUG
  return TP_NO_ERROR;
}
    
/** :の判定
 *
 * @param[in] ss 入力ファイルポインタ
 * @param[in,out] buffer 入力文字列
 * @param[out] answer :か否か
 * @return リターンコード
 *
 */
TextParserError TextParserTree::parseColon(std::stringstream& ss, std::string& buffer, bool& answer)
{
  TextParserError ret;

  ret = removeCommentsAndSpaces(ss, buffer);
  if (ret != TP_NO_ERROR) return ret;

  ret = parseColon(buffer, answer);
  if (ret != TP_NO_ERROR) return ret;

  return TP_NO_ERROR;
}
    
/** :の判定
 *
 * @param[in,out] buffer 入力文字列
 * @param[out] answer :か否か
 * @return リターンコード
 *
 */
TextParserError TextParserTree::parseColon(std::string& buffer, bool& answer)
{
  answer = false;

  // 最初が : ?
  if (buffer[0] != ':') return TP_NO_ERROR; 

  answer = true;
  buffer = buffer.substr(1);

#ifdef TP_DEBUG
  if (_debug_write) {
    std::cout << ":" << std::endl;
  }
#endif //TP_DEBUG
  return TP_NO_ERROR;
}
    
/** デリミタの判定
 *
 * @param[in] ss 入力ファイルポインタ
 * @param[in,out] buffer 入力文字列
 * @param[out] answer デリミタか否か
 * @return リターンコード
 *
 */
TextParserError TextParserTree::parseDelimiter(std::stringstream& ss, std::string& buffer, bool& answer)
{
  TextParserError ret;

  ret = removeCommentsAndSpaces(ss, buffer);
  if (ret != TP_NO_ERROR) return ret;

  ret = parseDelimiter(buffer, answer);
  if (ret != TP_NO_ERROR) return ret;

  return TP_NO_ERROR;
}

/** デリミタの判定
 *
 * @param[in,out] buffer 入力文字列
 * @param[out] answer デリミタか否か
 * @return リターンコード
 *
 */
TextParserError TextParserTree::parseDelimiter(std::string& buffer, bool& answer)
{
  answer = false;

  // 最初が , ?
  if (buffer[0] != ',') return TP_NO_ERROR;

  answer = true;
  buffer = buffer.substr(1);

#ifdef TP_DEBUG
  if (_debug_write) {
    std::cout << "デリミタ" << std::endl;
  }
#endif
  return TP_NO_ERROR;
}

/** ベクトルの終わりの判定
 *
 * @param[in] ss 入力ファイルポインタ
 * @param[in,out] buffer 入力文字列
 * @param[out] answer ベクトルの終わりか否か
 * @return リターンコード
 *
 */
TextParserError TextParserTree::parseEndOfVector(std::stringstream& ss, std::string& buffer, bool& answer)
{
  TextParserError ret;

  ret = removeCommentsAndSpaces(ss, buffer);
  if (ret != TP_NO_ERROR) return ret;

  ret = parseEndOfVector(buffer, answer);
  if (ret != TP_NO_ERROR) return ret;

  return TP_NO_ERROR;
}

/** ベクトルの終わりの判定
 *
 * @param[in,out] buffer 入力文字列
 * @param[out] answer ベクトルの終わりか否か
 * @return リターンコード
 *
 */
TextParserError TextParserTree::parseEndOfVector(std::string& buffer, bool& answer)
{
  answer = false;

  // 最初が ) ?
  if (buffer[0] != ')') return TP_NO_ERROR; 

  answer = true;
  buffer = buffer.substr(1);

#ifdef TP_DEBUG
  if (_debug_write) {
    std::cout << "ベクトルの終わり" << std::endl;
  }
#endif //TP_DEBUG
  return TP_NO_ERROR;
}

/** 値が正常か判定
 *
 * @param[in] value 値
 * @param[in] value_type 値のタイプ
 * @return 値が正常か否か
 *
 */
bool TextParserTree::isCorrectValue(std::string& value, TextParserValueType& value_type)
{
  if (value_type == TP_STRING_VALUE) {          // 文字列
    for (unsigned int i = 0; i < value.size(); i++) {
      if (value[i] >= 'a' && value[i] <= 'z') continue;
      if (value[i] >= 'A' && value[i] <= 'Z') continue;
      if (value[i] >= '0' && value[i] <= '9') continue;
      if (value[i] == '_') continue;
      if (value[i] == '-') continue;
      if (value[i] == '.') continue;
      if (value[i] == '/') continue;
      // 使用できない文字
      return false;
    }
  } else if (value_type == TP_NUMERIC_VALUE) {  // 数値
    for (unsigned int i = 0; i < value.size(); i++) {
      if (value[i] >= '0' && value[i] <= '9') continue;
      if (value[i] == '+') continue;
      if (value[i] == '-') continue;
      if (value[i] == '.') continue;
      if (value[i] == ' ') continue;
      if (value[i] == 'e') continue;
      if (value[i] == 'd') continue;
      if (value[i] == 'E') continue;
      if (value[i] == 'D') continue;
      // 使用できない文字
      return false;
    }
  } else {
    return false;
  }

  return true;
}

/** 値の設定
 *
 * @param[in] value 値
 * @param[in] value_type 値のタイプ
 * @return リターンコード
 *
 */
TextParserError TextParserTree::setValue(std::string& value, TextParserValueType& value_type)
{
  TextParserError ret;
  if (_current_element != 0 && _current_element->_type == TP_LEAF_ELEMENT) {
    TextParserLeaf *current_leaf = (TextParserLeaf *)_current_element;
    if (current_leaf->_value_type == TP_DEPENDENCE_VALUE) {   // 依存関係値
      TextParserValue *value_element = current_leaf->_value;    // 値のエレメントを取得
      if (value_element == 0) {
	return TextParserErrorHandler(TP_ILLEGAL_VALUE_ERROR, current_leaf->_label);
      }
      value_element->_value = value;
      value_element->_value_type = value_type;
      current_leaf->_value_type = value_type;
    } else {
      TextParserValue *value_element;
      try {
	value_element = new TextParserValue(value, value_type);    // 値のエレメントを新規に作成
	if (value_element == 0) {
	  return TextParserErrorHandler(TP_MEMORY_ALLOCATION_ERROR, "value_element");
	}
      } catch (std::exception ex) {
	return TextParserErrorHandler(TP_MEMORY_ALLOCATION_ERROR, "value_element");
      }
      ret = current_leaf->setElement(value_element);
      if (ret != TP_NO_ERROR) { // カレントリーフに追加
	return TextParserErrorHandler(ret, value);
      }
    }
  } else {
    return TextParserErrorHandler(TP_ILLEGAL_CURRENT_ELEMENT_ERROR, "");
  }

#ifdef TP_DEBUG
  if (_debug_write) {
    std::cout << "★ 値の設定 : ";
    debugWriteValue(value, value_type);
    std::cout << std::endl;
  }
#endif
  return TP_NO_ERROR;
}

/** 未定義値の設定
 *
 * @return リターンコード
 *
 */
TextParserError TextParserTree::setUndefinedValue()
{
  TextParserError ret;
  if (_current_element != 0 && _current_element->_type == TP_LEAF_ELEMENT) {
    TextParserLeaf *current_leaf = (TextParserLeaf *)_current_element;
    if (current_leaf->_value_type == TP_DEPENDENCE_VALUE) {   // 依存関係値
      TextParserValue *value_element = current_leaf->_value;    // 値のエレメントを取得
      if (value_element == 0) {
	return TextParserErrorHandler(TP_ILLEGAL_VALUE_ERROR, current_leaf->_label);
      }
      value_element->_value = "UNDEF";
      value_element->_value_type = TP_UNDEFFINED_VALUE;
      current_leaf->_value_type = TP_UNDEFFINED_VALUE;
    } else {
      TextParserValue *value_element;
      try {
	value_element = new TextParserValue("UNDEF", TP_UNDEFFINED_VALUE);    // 値のエレメントを新規に作成
	if (value_element == 0) {
	  return TextParserErrorHandler(TP_MEMORY_ALLOCATION_ERROR, "value_element");
	}
      } catch (std::exception ex) {
	return TextParserErrorHandler(TP_MEMORY_ALLOCATION_ERROR, "value_element");
      }
      ret = current_leaf->setElement(value_element);
      if (ret != TP_NO_ERROR) { // カレントリーフに追加
	return TextParserErrorHandler(ret, "UNDEF");
      }
    }
  } else {
    return TextParserErrorHandler(TP_ILLEGAL_CURRENT_ELEMENT_ERROR, "");
  }

#ifdef TP_DEBUG
  if (_debug_write) {
    std::cout << "★ 未定義値の設定" << std::endl;
  }
#endif
  return TP_NO_ERROR;
}

/** ベクトル値の設定
 *
 * @param[in] ss 入力ファイルポインタ
 * @param[in,out] buffer 入力文字列
 * @return リターンコード
 *
 * 値としてベクトル値の文字列を設定する
 *
 */
TextParserError TextParserTree::setVectorValue(std::stringstream& ss, std::string& buffer)
{
  TextParserError ret;
  bool answer = false;
  TextParserValue *value_element = 0;

  // 値の追加
  if (_current_element != 0 && _current_element->_type == TP_LEAF_ELEMENT) {
    TextParserLeaf *current_leaf = (TextParserLeaf *)_current_element;
    try {
      value_element = new TextParserValue("", TP_DEPENDENCE_VALUE);    // 値のエレメントを新規に作成
      if (value_element == 0) {
	return TextParserErrorHandler(TP_MEMORY_ALLOCATION_ERROR, "value_element");
      }
    } catch (std::exception ex) {
      return TextParserErrorHandler(TP_MEMORY_ALLOCATION_ERROR, "value_element");
    }
    ret = current_leaf->setElement(value_element);
    if (ret != TP_NO_ERROR) return TextParserErrorHandler(ret, "");
  } else {
    return TextParserErrorHandler(TP_ILLEGAL_CURRENT_ELEMENT_ERROR, "");
  }
  value_element->_value = "(";
  TextParserValueType value_type;
  std::string value;
  bool first = true;              // 一つ目？
  TextParserValueType first_value_type; // 一つ目の値のタイプ
  while(true) {
    // 値の取得
    ret = parseUndefinedValue(ss, buffer, answer);
    if (ret != TP_NO_ERROR) return ret;
    if (answer) {          // 未定義の値?
      value_element->_value += "UNDEF";
      value_element->_value_type = TP_VECTOR_UNDEFFINED;
      value_type = TP_UNDEFFINED_VALUE;
    } else {
      ret = parseValue(ss, buffer, answer, value, value_type);
      if (ret != TP_NO_ERROR) return ret;
      if (answer) {       // スカラー値?
	if (value_type == TP_NUMERIC_VALUE) {
	  value_element->_value += value;
	  value_element->_value_type = TP_VECTOR_NUMERIC;
	} else if (value_type == TP_STRING_VALUE) {
	  value_element->_value += "\"" + value + "\"";
	  value_element->_value_type = TP_VECTOR_STRING;
	} else {
	  return TextParserErrorHandler(TP_ILLEGAL_VALUE_TYPE_ERROR, buffer);
	}
      } else {
	return TextParserErrorHandler(TP_ILLEGAL_VALUE_ERROR, buffer);
      }
    }
    // 値の種類のチェックのため最初の値のタイプを保存する
    if (first) {
      first_value_type = value_type;
      first = false;
    }
    ret = parseDelimiter(ss, buffer, answer);
    if (answer) {   // デリミタ?
      value_element->_value += ",";
    } else {
      ret = parseEndOfVector(ss, buffer, answer);
      if (ret != TP_NO_ERROR) return ret;
      if (answer) {   // ベクトルの終わり
	value_element->_value += ")";
	break;
      } else {
	return TextParserErrorHandler(TP_MISSING_VECTOR_END_ERROR, buffer);
      }
    }
  }
#ifdef TP_DEBUG
  if (_debug_write) {
    std::cout << "<ベクトル値>：" + value_element->_value << std::endl;
  }
#endif
  return TP_NO_ERROR;
}

/** 依存関係式の設定
 *
 * @param[in] ss 入力ファイルポインタ
 * @param[in,out] buffer 入力文字列
 * @return リターンコード
 *
 * 値として依存関係式の文字列を設定する
 *
 */
TextParserError TextParserTree::setDependenceExpression(std::stringstream& ss, std::string& buffer)
{
  TextParserError ret;
  bool answer = false;
  TextParserValue *value_element = 0;

#ifdef MYDEBUG3
  std::cout <<"aaa " << buffer<< std::endl;
#endif //MYDEBUG3

  // 値の追加
  if (_current_element != 0 && _current_element->_type == TP_LEAF_ELEMENT) {
    TextParserLeaf *current_leaf = (TextParserLeaf *)_current_element;
    try {
      value_element = new TextParserValue("", TP_DEPENDENCE_VALUE);    // 値のエレメントを新規に作成
      if (value_element == 0) {
	return TextParserErrorHandler(TP_MEMORY_ALLOCATION_ERROR, "value_element");
      }
    } catch (std::exception ex) {
      return TextParserErrorHandler(TP_MEMORY_ALLOCATION_ERROR, "value_element");
    }
    ret = current_leaf->setElement(value_element);
    if (ret != TP_NO_ERROR) return TextParserErrorHandler(ret, "");
  } else {
    return TextParserErrorHandler(TP_ILLEGAL_CURRENT_ELEMENT_ERROR, "");
  }
  ret = parseOpenBrancket(ss, buffer, answer);
#ifdef MYDEBUG3
  std::cout <<"bbb " << buffer<< std::endl;
#endif //MYDEBUG3

  if (ret != TP_NO_ERROR) return ret;
  if (answer) {                   // 条件式の開始？
    value_element->_value = "(";
    ret = setConditionalExpression(ss, buffer);  // 条件式の登録
#ifdef MYDEBUG3
    std::cout <<"ccc " << buffer<< std::endl;
#endif //MYDEBUG3
    if (ret != TP_NO_ERROR) return ret;
    ret = parseQuestion(ss, buffer, answer);
    if (ret != TP_NO_ERROR) return ret;
    if (answer) {
      value_element->_value += "?";
      ret = setDependenceValue(ss, buffer);    // 値の登録
      if (ret != TP_NO_ERROR) return ret;
      ret = parseColon(ss, buffer, answer);
      if (ret != TP_NO_ERROR) return ret;
      if (answer) {
	value_element->_value += ":";
	ret = setDependenceValue(ss, buffer);    // 値の登録
	if (ret != TP_NO_ERROR) return ret;
      } else {
#ifdef MYDEBUG3	      
	std::cout << "setDependenceExpression1 buffer=" 
		  <<buffer<<std::endl;
#endif
	return TextParserErrorHandler(TP_ILLEGAL_DEPENDENCE_EXPRESSION_ERROR, buffer);
      }
    } else {
#ifdef MYDEBUG3	      
      std::cout << "setDependenceExpression2 buffer|" 
		<<buffer <<" answer "<< answer <<std::endl;
#endif
      return TextParserErrorHandler(TP_ILLEGAL_DEPENDENCE_EXPRESSION_ERROR, buffer);
    }
  } else {
#ifdef MYDEBUG3	      
    std::cout << "setDependenceExpression3 buffer=" 
	      <<buffer <<" answer "<< answer <<std::endl;
#endif
    return TextParserErrorHandler(TP_ILLEGAL_DEPENDENCE_EXPRESSION_ERROR, buffer);
  }

#ifdef TP_DEBUG
  if (_debug_write) {
    std::cout << "<依存関係式>：" + value_element->_value << std::endl;
  }
#endif
  return TP_NO_ERROR;
}

/** 条件式の登録
 *
 * @param[in] ss 入力ファイルポインタ
 * @param[in,out] buffer 入力文字列
 * @return リターンコード
 *
 */
TextParserError TextParserTree::setConditionalExpression(std::stringstream& ss, std::string& buffer)
{
  TextParserError ret;
  bool answer = false;
  TextParserValue *value_element = 0;
  std::string value;
  TextParserValueType value_type;


#ifdef MYDEBUG3
  std::cout << "setConditionalExpression 1 |"<<buffer <<std::endl;
#endif //MYDEBUG3

  // 値の取得
  if (_current_element != 0 && _current_element->_type == TP_LEAF_ELEMENT) {
    TextParserLeaf *current_leaf = (TextParserLeaf *)_current_element;
    value_element = current_leaf->_value;
    if (value_element == 0) {
      return TextParserErrorHandler(TP_ILLEGAL_VALUE_ERROR,
				    current_leaf->_label);
    }
  } else {
    return TextParserErrorHandler(TP_ILLEGAL_CURRENT_ELEMENT_ERROR, "");
  }

 loop:

  std::string label;
  ret = parseLabel(ss, buffer, answer, label);
#ifdef MYDEBUG3
  std::cout << "setConditionalExpression 2 |"<<buffer <<std::endl;
#endif //MYDEBUG3

  if (ret != TP_NO_ERROR) return ret;
  if (answer) { // parseLabel yes
    value_element->_value += "\"" + label + "\"";
    ret = parseEqual(ss, buffer, answer);
#ifdef MYDEBUG3
    std::cout << "setConditionalExpression 3 |"<<buffer <<std::endl;
#endif //MYDEBUG3

    if (ret != TP_NO_ERROR) return ret;
    if (answer) {  // "==" ?
      value_element->_value += "==";
      ret = parseUndefinedValue(ss, buffer, answer);
#ifdef MYDEBUG3
      std::cout << "setConditionalExpression 4 |"<<buffer <<std::endl;
#endif //MYDEBUG3

      if (ret != TP_NO_ERROR) return ret;
      if (answer) {
	if (_result_bool == TP_UNDEFINED_BOOL) {
	  TextParserErrorHandler(TP_UNDEFINED_VALUE_USED_WARNING, label);
	}
	value_element->_value += "UNDEF";
      } else {
	ret = parseValue(ss, buffer, answer, value, value_type);
#ifdef MYDEBUG3
	std::cout << "setConditionalExpression 5 |"<<buffer <<std::endl;
#endif //MYDEBUG3

	if (ret != TP_NO_ERROR) return ret;
	if (answer) {
	  if (!isCorrectValue(value, value_type)) {    // 値のチェック
	    return TextParserErrorHandler(TP_ILLEGAL_VALUE_ERROR, value);
	  }
	  if (value_type == TP_STRING_VALUE) {
	    value_element->_value += "\"" + value + "\"";
	  } else if (value_type == TP_NUMERIC_VALUE) {
	    value_element->_value += value;
	  }
	} else {
	  return TextParserErrorHandler(TP_MISSING_VALUE_ERROR, buffer);
	}
      }
    } else {
      ret = parseNotEqual(ss, buffer, answer);
#ifdef MYDEBUG3
      std::cout << "setConditionalExpression 6 |"<<buffer <<std::endl;
#endif //MYDEBUG3


      if (ret != TP_NO_ERROR) return ret;
      if (answer) {    // "!=" ?
	value_element->_value += "!=";
	ret = parseUndefinedValue(ss, buffer, answer);
#ifdef MYDEBUG3
	std::cout << "setConditionalExpression 7 |"<<buffer <<std::endl;
#endif //MYDEBUG3

	if (ret != TP_NO_ERROR) return ret;
	if (answer) { //undefined value
	  if (_result_bool == TP_UNDEFINED_BOOL) {
	    TextParserErrorHandler(TP_UNDEFINED_VALUE_USED_WARNING, label);
	  }
	  value_element->_value += "UNDEF";
	} else { // not undefined value
	  ret = parseValue(ss, buffer, answer, value, value_type);
#ifdef MYDEBUG3
	  std::cout << "setConditionalExpression 8 |"<<buffer <<std::endl;
#endif //MYDEBUG3
	  if (ret != TP_NO_ERROR) return ret;
	  if (answer) { // parseValue answer 
	    if (!isCorrectValue(value, value_type)) {    // 値のチェック
	      return TextParserErrorHandler(TP_ILLEGAL_VALUE_ERROR, value);
	    }
	    if (value_type == TP_STRING_VALUE) { //string
	      value_element->_value += "\"" + value + "\"";
	    } else if (value_type == TP_NUMERIC_VALUE) { //numerical
	      value_element->_value += value;
	    } else {
	      return TextParserErrorHandler(TP_ILLEGAL_VALUE_TYPE_ERROR, 
					    value);
	    }
	  } else { //parseValue answer false
	    return TextParserErrorHandler(TP_MISSING_VALUE_ERROR, buffer);
	  } //parseValue answer 
	} // undefined Value answer
      } else { // "!=" answer false 
	return TextParserErrorHandler(TP_MISSING_EQUAL_NOT_EQUAL_ERROR, buffer);
      } // "!=" answer end
    }
    ret = parseClosedBrancket(ss, buffer, answer);

#ifdef MYDEBUG3
    std::cout << "setConditionalExpression 9 |"
	      <<buffer << "| answer "<< answer <<std::endl;
#endif //MYDEBUG3

    if (ret != TP_NO_ERROR) return ret;
    if (answer) {   // 条件式の終了？
      value_element->_value += ")";
    } else {
      return TextParserErrorHandler(TP_MISSING_CLOSED_BRANCKET_ERROR, buffer);
    }
#ifdef MYDEBUG3
    std::cout << value_element->_value << std::endl;
#endif //MYDEBUG3

    buffer=TextParserRemoveHeadSpaces(buffer);
    int qpos=buffer.find("?");
#ifdef MYDEBUG3
    std::cout << "qpos "<< qpos << " buffer "<< buffer<< std::endl;
#endif //MYDEBUG3
    if (qpos<0) 
      return TextParserErrorHandler(
				    TP_ILLEGAL_CONDITION_EXPRESSION_ERROR,
				    buffer);

    //    if (qpos!=0) goto loop;
    //    if(qpos!=0){
    //      ret = setConditionalExpression(ss, buffer);
    //    }

  } else {   // parseLabel no ???
    ret = parseOpenBrancket(ss, buffer, answer);

#ifdef MYDEBUG3
    std::cout << "setConditionalExpression 10 |"<< buffer <<std::endl;
#endif //MYDEBUG3

    if (ret != TP_NO_ERROR) return ret;
    if (answer) {  // 条件式の開始？
      value_element->_value += "(";
      ret = setConditionalExpression(ss, buffer);

#ifdef MYDEBUG3
      std::cout << "setConditionalExpression 11 |"<<buffer <<std::endl;
#endif //MYDEBUG3

      if (ret != TP_NO_ERROR) return ret;
      ret = parseAnd(ss, buffer, answer);

#ifdef MYDEBUG3
      std::cout << "setConditionalExpression 12 |"<<buffer <<std::endl;
#endif //MYDEBUG3

      if (ret != TP_NO_ERROR) return ret;
      if (answer) {   // &&
	value_element->_value += "&&";
	ret = parseOpenBrancket(ss, buffer, answer);
	if (ret != TP_NO_ERROR) return ret;
	if (answer) {  // 条件式の開始？
	  value_element->_value += "(";
	  ret = setConditionalExpression(ss, buffer);
	  if (ret != TP_NO_ERROR) return ret;
	  ret = parseClosedBrancket(ss, buffer, answer);
	  if (ret != TP_NO_ERROR) return ret;
	  if (answer) {   // ')' ?
	    value_element->_value += ")";
	  } else { // ')' ?
	    return TextParserErrorHandler(TP_MISSING_CLOSED_BRANCKET_ERROR, buffer);
	  } // ')' ?
	} else { // 条件式の開始？
	  return TextParserErrorHandler(TP_MISSING_CONDITION_EXPRESSION_ERROR, buffer);
	} // 条件式の開始？
      } else { // not &&
	ret = parseOr(ss, buffer, answer);
	if (ret != TP_NO_ERROR) return ret;
	if (answer) { // parseOr true
	  value_element->_value += "||";
	  ret = parseOpenBrancket(ss, buffer, answer);
	  if (ret != TP_NO_ERROR) return ret;
	  if (answer) {  // 条件式の開始？ true
	    value_element->_value += "(";
	    ret = setConditionalExpression(ss, buffer);
	    if (ret != TP_NO_ERROR) return ret;
	    ret = parseClosedBrancket(ss, buffer, answer);
	    if (ret != TP_NO_ERROR) return ret;
	    if (answer) {   // ')' true
	      value_element->_value += ")";
	    } else {   // ')' false 
	      return
		TextParserErrorHandler(TP_MISSING_CLOSED_BRANCKET_ERROR,
				       buffer);
	    }   // ')' end
	  } else { // 条件式の開始？ false 
	    return 
	      TextParserErrorHandler(TP_MISSING_CONDITION_EXPRESSION_ERROR,
				     buffer);
	  } // 条件式の開始？ end
	} else {    // 無駄な括弧を許す ((条件式)) false
	  ret = parseClosedBrancket(ss, buffer, answer);
	  if (ret != TP_NO_ERROR) return ret;
	  if (answer) {   // ')'  true
	    value_element->_value += ")";
	  } else {  // ')'  false
	    return TextParserErrorHandler(TP_MISSING_AND_OR_ERROR, buffer);
	  }  // ')'  end
	} // 無駄な括弧を許す ((条件式)) end
      } 
    } else { //条件式の開始？　
#ifdef MYDEBUG3
      std::cout << "no no no ( " << std::endl;
#endif
     
      return
	TextParserErrorHandler(TP_ILLEGAL_CONDITION_EXPRESSION_ERROR,
			       buffer);
    }
  }

#ifdef MYDEBUG3
  std::cout << "no error:set conditional " << std::endl;
#endif

  return TP_NO_ERROR;
}

/** 依存関係式の値の登録
 *
 * @param[in] ss 入力ファイルポインタ
 * @param[in,out] buffer 入力文字列
 * @return リターンコード
 *
 */
TextParserError TextParserTree::setDependenceValue(std::stringstream& ss, std::string& buffer)
{
  TextParserError ret;
  bool answer = false;
  TextParserValue *value_element = 0;
  std::string value;
  TextParserValueType value_type;

  // 値の取得
  if (_current_element != 0 && _current_element->_type == TP_LEAF_ELEMENT) {
    TextParserLeaf *current_leaf = (TextParserLeaf *)_current_element;
    value_element = current_leaf->_value;
    if (value_element == 0) {
      return TextParserErrorHandler(TP_ILLEGAL_VALUE_ERROR, current_leaf->_label);
    }
  } else {
    return TextParserErrorHandler(TP_ILLEGAL_CURRENT_ELEMENT_ERROR, "");
  }
  while (true) {
    ret = parseVectorValue(ss, buffer, answer);
    if (ret != TP_NO_ERROR) return ret;
    if (answer) {                // ベクトル値?
      value_element->_value += "(";
      bool vector_end = false;
      while (true) {
	ret = parseValue(ss, buffer, answer, value, value_type);
	if (ret != TP_NO_ERROR) return ret;
	if (answer) {       // スカラー値?
	  if (!isCorrectValue(value, value_type)) {   // 値のチェック
	    return TextParserErrorHandler(TP_ILLEGAL_VALUE_ERROR, value);
	  }
	  if (value_type == TP_STRING_VALUE) {
	    value_element->_value += "\"" + value + "\"";
	  } else if (value_type == TP_NUMERIC_VALUE) {
	    value_element->_value += value;
	  }
	  while (true) {
	    ret = parseDelimiter(ss, buffer, answer);
	    if (ret != TP_NO_ERROR) return ret;
	    if (answer) {                                       // デリミタ?
	      value_element->_value += ",";
	      break;
	    }
	    ret = parseEndOfVector(ss, buffer, answer);
	    if (ret != TP_NO_ERROR) return ret;
	    if (answer) {             // ベクトルの終わり
	      value_element->_value += ")";
	      vector_end = true;
	      break;
	    }
#ifdef MYDEBUG3	      
	    std::cout << "setDependenceValue1 buffer=" 
		      <<buffer<<std::endl;
#endif
	    return TextParserErrorHandler(TP_ILLEGAL_DEPENDENCE_EXPRESSION_ERROR, buffer);
	  }
	} else {
#ifdef MYDEBUG3	      
	  std::cout << "setDependenceValue2 buffer=" 
		    <<buffer<<std::endl;
#endif
	  return TextParserErrorHandler(TP_ILLEGAL_DEPENDENCE_EXPRESSION_ERROR, buffer);
	}
	if (vector_end) break;
      }
      break;
    }
    ret = parseUndefinedValue(ss, buffer, answer);
    if (ret != TP_NO_ERROR) return ret;
    if (answer) {      // 未定義の値?
      value_element->_value += "UNDEF";
      break;
    }
    ret = parseValue(ss, buffer, answer, value, value_type);
    if (ret != TP_NO_ERROR) return ret;
    if (answer) {    // スカラー値?
      if (!isCorrectValue(value, value_type)) {       // 値のチェック
	return TextParserErrorHandler(TP_ILLEGAL_VALUE_ERROR, value);
      }
      if (value_type == TP_STRING_VALUE) {
	value_element->_value += "\"" + value + "\"";
      } else if (value_type == TP_NUMERIC_VALUE) {
	value_element->_value += value;
      }
      break;
    }
#ifdef MYDEBUG3	      
    std::cout << "setDependenceValue3 buffer=" 
	      <<buffer<<std::endl;
#endif
    return TextParserErrorHandler(TP_ILLEGAL_DEPENDENCE_EXPRESSION_ERROR, buffer);
  }

  return TP_NO_ERROR;
}

/** 依存関係式の解析
 *
 * @param[in] leaf リーフ
 * @return リターンコード
 *
 */
TextParserError TextParserTree::parseDependenceExpression(TextParserLeaf *leaf)
{
  TextParserError ret;
  bool answer;
  TextParserValue *value_element = 0;
  TextParserBool result;
#ifdef MYDEBUG3
  std::cout << "parseDependenceExpression start" <<std::endl;
#endif //MYDEBUG3

  // 値の取得
  if (_current_element != 0 && _current_element->_type == TP_LEAF_ELEMENT) {
    TextParserLeaf *current_leaf = (TextParserLeaf *)_current_element;
    value_element = current_leaf->_value;
    if (value_element == 0) {
      return TextParserErrorHandler(TP_ILLEGAL_VALUE_ERROR, current_leaf->_label);
    }
  } else {
    return TextParserErrorHandler(TP_ILLEGAL_CURRENT_ELEMENT_ERROR, "");
  }

  std::string buffer = value_element->_value;

#ifdef MYDEBUG3
  std::cout << "buffer0 |"<<buffer<<"|" <<std::endl;
#endif //MYDEBUG3


#ifdef TP_DEBUG
  if (_debug_write) {
    std::cout << "依存関係式：" + buffer << std::endl;
  }
#endif

  ret = parseOpenBrancket(buffer, answer);
#ifdef MYDEBUG3
  std::cout << "buffer1 |"<<buffer<<"|" <<std::endl;
#endif //MYDEBUG3

  if (ret != TP_NO_ERROR) return ret;
  if (answer) {   //parseOpenBrancket true                 // 条件式の開始？
    ret = parseConditionalExpression(buffer, result);
#ifdef MYDEBUG3
    std::cout << "buffer2 |"<<buffer<<"|" <<std::endl;
#endif //MYDEBUG3

    if (ret != TP_NO_ERROR) {
      return ret;
    }
    ret = parseQuestion(buffer, answer);    // ?
#ifdef MYDEBUG3
    std::cout << "buffer3 |"<<buffer<<"|" <<std::endl;
#endif //MYDEBUG3

    if (ret != TP_NO_ERROR) return ret;
    if (answer) { //parseQuetion true
      ret = parseDependenceValue(buffer, result);
      // resultがtrueなら値を設定
#ifdef MYDEBUG3
      std::cout << "buffer4 |"<<buffer<<"|" <<std::endl;
#endif //MYDEBUG3
                         


      if (ret != TP_NO_ERROR) return ret;
      ret = parseColon(buffer, answer);
#ifdef MYDEBUG3
      std::cout << "buffer5 |"<<buffer<<"|" <<std::endl;
#endif //MYDEBUG3

      if (ret != TP_NO_ERROR) return ret;
      if (answer) { //parseColon true
	ret = parseDependenceValue(buffer, result == TP_TRUE ? TP_FALSE : TP_TRUE);   // resultがfalseなら値を設定
	if (ret != TP_NO_ERROR) return ret;
      } else { //parseColon false 
#ifdef MYDEBUG3	      
	std::cout << "parseDependenceExpression1a buffer=" 
		  <<buffer<<std::endl;
#endif
	return TextParserErrorHandler(TP_ILLEGAL_DEPENDENCE_EXPRESSION_ERROR, buffer);
      } //parseColon end
    } else { //parseQuetion false
#ifdef MYDEBUG3	      
      std::cout << "parseDependenceExpression2a buffer=" 
		<<buffer<<std::endl;
#endif
      return TextParserErrorHandler(TP_ILLEGAL_DEPENDENCE_EXPRESSION_ERROR, buffer);
    } //parseQuetion end
  } else {  //parseOpenBrancket false     // 条件式の開始？
#ifdef MYDEBUG3	      
    std::cout << "parseDependenceExpression3a buffer=" 
	      <<buffer<<std::endl;
#endif
    return TextParserErrorHandler(TP_ILLEGAL_DEPENDENCE_EXPRESSION_ERROR, buffer);
  }

  return TP_NO_ERROR;
}

/** 条件式の解析
 *
 * @param[in,out] buffer 入力文字列
 * @param[out] result 結果の論理値
 * @return リターンコード
 *
 */
TextParserError TextParserTree::parseConditionalExpression(std::string& buffer, TextParserBool& result)
{
#ifdef TP_DEBUG
  /*
    if (_debug_write) {
    std::cout << "<条件式の解析>" << std::endl;
    }
  */
#endif
  TextParserError ret;
  bool answer = false;
  std::string value;
  TextParserValueType value_type;
  TextParserBool left_result = TP_UNDEFINED_BOOL;
  TextParserBool right_result = TP_UNDEFINED_BOOL;

  std::string label;
  ret = parseLabel(buffer, answer, label);
  if (ret != TP_NO_ERROR) return ret;
  if (answer) {
    ret = parseEqual(buffer, answer);
    if (ret != TP_NO_ERROR) return ret;
    if (answer) {  // "==" ?
      ret = parseUndefinedValue(buffer, answer);
      if (ret != TP_NO_ERROR) return ret;
      if (answer) {
	if (_result_bool == TP_UNDEFINED_BOOL) {
	  TextParserErrorHandler(TP_UNDEFINED_VALUE_USED_WARNING, label);
	}
	result = TP_UNDEFINED_BOOL;
      } else {
	ret = parseValue(buffer, answer, value, value_type);
	if (ret != TP_NO_ERROR) return ret;
	if (answer) {
	  if (!isCorrectValue(value, value_type)) {    // 値のチェック
	    return TextParserErrorHandler(TP_ILLEGAL_VALUE_ERROR, value);
	  }
	  ret = resolveConditionalExpression(label, value, value_type, TP_TRUE, result);
	  if (ret != TP_NO_ERROR) return ret;
	  if (result == TP_UNDEFINED_BOOL) {
	    TextParserErrorHandler(TP_UNDEFINED_VALUE_USED_WARNING, label);
	  }
	} else {
	  return TextParserErrorHandler(TP_MISSING_VALUE_ERROR, buffer);
	}
      }
    } else {
      ret = parseNotEqual(buffer, answer);
      if (ret != TP_NO_ERROR) return ret;
      if (answer) {    // "!=" ?
	ret = parseUndefinedValue(buffer, answer);
	if (ret != TP_NO_ERROR) return ret;
	if (answer) {
	  if (_result_bool == TP_UNDEFINED_BOOL) {
	    TextParserErrorHandler(TP_UNDEFINED_VALUE_USED_WARNING, label);
	  }
	  result = TP_UNDEFINED_BOOL;
	} else {
	  ret = parseValue(buffer, answer, value, value_type);
	  if (ret != TP_NO_ERROR) return ret;
	  if (answer) {
	    if (!isCorrectValue(value, value_type)) {    // 値のチェック
	      return TextParserErrorHandler(TP_ILLEGAL_VALUE_ERROR, value);
	    }
	    ret = resolveConditionalExpression(label, value, value_type, TP_FALSE, result);
	    if (ret != TP_NO_ERROR) return ret;
	    if (result == TP_UNDEFINED_BOOL) {
	      TextParserErrorHandler(TP_UNDEFINED_VALUE_USED_WARNING, label);
	    }
	  } else {
	    return TextParserErrorHandler(TP_MISSING_VALUE_ERROR, buffer);
	  }
	}
      } else {
	return TextParserErrorHandler(TP_MISSING_EQUAL_NOT_EQUAL_ERROR, buffer);
      }
    }
    ret = parseClosedBrancket(buffer, answer);
    if (ret != TP_NO_ERROR) return ret;
    if (answer) {   // 条件式の終了？
    } else {
      return TextParserErrorHandler(TP_MISSING_CLOSED_BRANCKET_ERROR, buffer);
    }
  } else {
    ret = parseOpenBrancket(buffer, answer);
    if (ret != TP_NO_ERROR) return ret;
    if (answer) {  // 条件式の開始？
      ret = parseConditionalExpression(buffer, left_result);
      if (ret != TP_NO_ERROR) return ret;
      ret = parseAnd(buffer, answer);
      if (ret != TP_NO_ERROR) return ret;
      if (answer) {   // &&
	ret = parseOpenBrancket(buffer, answer);
	if (ret != TP_NO_ERROR) return ret;
	if (answer) {  // 条件式の開始？
	  ret = parseConditionalExpression(buffer, right_result);
	  if (ret != TP_NO_ERROR) return ret;
	  result = resolveAnd(left_result, right_result);
	  ret = parseClosedBrancket(buffer, answer);
	  if (ret != TP_NO_ERROR) return ret;
	  if (answer) {   // ')' ?
	  } else {
	    return TextParserErrorHandler(TP_MISSING_CLOSED_BRANCKET_ERROR, buffer);
	  }
	} else {
	  return TextParserErrorHandler(TP_MISSING_CONDITION_EXPRESSION_ERROR, buffer);
	}
      } else {
	ret = parseOr(buffer, answer);
	if (ret != TP_NO_ERROR) return ret;
	if (answer) {
	  ret = parseOpenBrancket(buffer, answer);
	  if (ret != TP_NO_ERROR) return ret;
	  if (answer) {  // 条件式の開始？
	    ret = parseConditionalExpression(buffer, right_result);
	    if (ret != TP_NO_ERROR) return ret;
	    result = resolveOr(left_result, right_result);
	    ret = parseClosedBrancket(buffer, answer);
	    if (ret != TP_NO_ERROR) return ret;
	    if (answer) {   // ')' ?
	    } else {
	      return TextParserErrorHandler(TP_MISSING_CLOSED_BRANCKET_ERROR, buffer);
	    }
	  } else {
	    return TextParserErrorHandler(TP_MISSING_CONDITION_EXPRESSION_ERROR, buffer);
	  }
	} else {    // 無駄な括弧を許す ((条件式))
	  ret = parseClosedBrancket(buffer, answer);
	  if (ret != TP_NO_ERROR) return ret;
	  if (answer) {   // ')' ?
	    result = left_result;   // 左の式の結果を設定
	  } else {
	    return TextParserErrorHandler(TP_MISSING_AND_OR_ERROR, buffer);
	  }
	}
      }
    } else {
      return TextParserErrorHandler(TP_ILLEGAL_CONDITION_EXPRESSION_ERROR, buffer);
    }
  }

  return TP_NO_ERROR;
}

/** 依存関係の値の解析
 *
 * @param[in,out] buffer 文字列
 * @param[in] set 文字列
 * @return リターンコード
 *
 */
TextParserError TextParserTree::parseDependenceValue(std::string& buffer, TextParserBool set)
{
  TextParserError ret;
  bool answer = false;
  std::string value;
  TextParserValueType value_type;

  while (true) {
    ret = parseVectorValue(buffer, answer);
    if (ret != TP_NO_ERROR) return ret;
    if (answer) {                // ベクトル値?
      bool vector_end = false;
      while (true) {
	ret = parseValue(buffer, answer, value, value_type);
	if (ret != TP_NO_ERROR) return ret;
	if (answer) {       // スカラー値?
	  if (!isCorrectValue(value, value_type)) {   // 値のチェック
	    return TextParserErrorHandler(TP_ILLEGAL_VALUE_ERROR, value);
	  }
	  if (set == TP_TRUE) {
	    ret = setValue(value, value_type);          // 値の設定
	    if (ret != TP_NO_ERROR) return ret;
	  }
	  while (true) {
	    ret = parseDelimiter(buffer, answer);
	    if (ret != TP_NO_ERROR) return ret;
	    if (answer) {                                       // デリミタ?
	      break;
	    }
	    ret = parseEndOfVector(buffer, answer);
	    if (ret != TP_NO_ERROR) return ret;
	    if (answer) {             // ベクトルの終わり
	      vector_end = true;
	      break;
	    }
#ifdef MYDEBUG3	      
	    std::cout << "parseDependenceValue1 buffer=" 
		      <<buffer<<std::endl;
#endif
	    return TextParserErrorHandler(TP_ILLEGAL_DEPENDENCE_EXPRESSION_ERROR, buffer);
	  }
	} else {
#ifdef MYDEBUG3	      
	  std::cout << "parseDependenceValue2 buffer=" 
		    <<buffer<<std::endl;
#endif
	  return TextParserErrorHandler(TP_ILLEGAL_DEPENDENCE_EXPRESSION_ERROR, buffer);
	}
	if (vector_end) break;
      }
      break;
    }
    ret = parseUndefinedValue(buffer, answer);
    if (ret != TP_NO_ERROR) return ret;
    if (answer) {      // 未定義の値?
      if (set == TP_TRUE) {
	ret = setUndefinedValue();          // 未定義値の設定
	if (ret != TP_NO_ERROR) return ret;
      }
      break;
    }
    ret = parseValue(buffer, answer, value, value_type);
    if (ret != TP_NO_ERROR) return ret;
    if (answer) {    // スカラー値?
      if (!isCorrectValue(value, value_type)) {       // 値のチェック
	return TextParserErrorHandler(TP_ILLEGAL_VALUE_ERROR, value);
      }
      if (set == TP_TRUE) {
	ret = setValue(value, value_type);          // 値の設定
	if (ret != TP_NO_ERROR) return ret;
      }
      break;
    }
#ifdef MYDEBUG3	      
    std::cout << "parseDependenceValue3 buffer=" 
	      <<buffer<<std::endl;
#endif
    return TextParserErrorHandler(TP_ILLEGAL_DEPENDENCE_EXPRESSION_ERROR, buffer);
  }

  return TP_NO_ERROR;
}

/** 条件式の論理値計算
 *
 * @param[in] label ラベル
 * @param[in] value 値
 * @param[in] value_type 値のタイプ
 * @param[in] is_equal 等号か？
 * @param[out] result 論理値
 * @return リターンコード
 *
 */
TextParserError TextParserTree::resolveConditionalExpression(std::string& label, std::string& value, TextParserValueType& value_type, const TextParserBool is_equal, TextParserBool& result)
{
#ifdef TP_DEBUG
  if (_debug_write) {
    std::cout << "条件式 (" + label;
    if (is_equal == TP_TRUE) {
      std::cout << " == ";
    } else {
      std::cout << " != ";
    }
    debugWriteValue(value, value_type);
    std::cout << ")" << std::endl;
  }
#endif
  TextParserError ret;
  std::string path = label;
  if (path[0] != '/') {  // ラベルが相対パスなら絶対パスに変換
    if (_current_element->_type == TP_LEAF_ELEMENT) {
      TextParserLeaf *current_leaf = (TextParserLeaf *)_current_element;
      std::string parent_path;
      if (!path.compare(0, 3, "../")) {        // ラベルが”../”から始まる場合
	path = path.substr(3);
	if (current_leaf->_parent != 0) {
	  TextParserElement *parent = current_leaf->_parent;
	  ret = GetElementAbsolutePath(parent->_parent, parent_path);
	  if (ret != TP_NO_ERROR) return ret;
	  path = parent_path + "/" + path;
	}
      } else {
	if (!path.compare(0, 2, "./")) {   // ラベルが”./”から始まる場合
	  path = path.substr(2);
	}
	current_leaf = (TextParserLeaf *)_current_element;
	ret = GetElementAbsolutePath(current_leaf->_parent, parent_path);
	if (ret != TP_NO_ERROR) return ret;
	path = parent_path + "/" + path;
      }
    } else {
      return TextParserErrorHandler(TP_ILLEGAL_CURRENT_ELEMENT_ERROR, label);
    }
  }
  TextParserValue *leaf_value = 0;
  ret = getLeafValue(path, &leaf_value);
  if (ret != TP_NO_ERROR) return ret;
  if (leaf_value == 0) {
    return TextParserErrorHandler(TP_ILLEGAL_VALUE_ERROR, label);
  }

  if (leaf_value->_value_type == TP_DEPENDENCE_VALUE) {
#ifdef TP_DEBUG
    if (_debug_write) {
      std::cout << "未解決の依存関係付き値：" + label << std::endl;
    }
#endif
    return TP_UNRESOLVED_LABEL_USED_WARNING;
  } else if (leaf_value->_value_type == TP_STRING_VALUE) {
    if (TextParserStringCompare(value, leaf_value->_value)) {
      result = TP_TRUE;
    } else {
      result = TP_FALSE;
    }
  } else if (leaf_value->_value_type == TP_NUMERIC_VALUE) {
    // 数値はdoubleに変換して比較
    std::string val0 = CorrectValueString(value);
    std::string val1 = CorrectValueString(leaf_value->_value);
    double dval0, dval1;
#if 0
    try {
      std::stringstream ss;
      ss << val0;
      ss >> dval0;
      ss.clear();
      ss << val1;
      ss >> dval1;
      if (dval0 == dval1) {
	result = TP_TRUE;
      } else {
	result = TP_FALSE;
      }
    } catch (std::exception ex) {
      return TextParserErrorHandler(TP_ILLEGAL_NUMERIC_VALUE_ERROR, value);
    }
#else
    try {
      sscanf(val0.c_str(), "%lf", &dval0);
      sscanf(val1.c_str(), "%lf", &dval1);
      if (dval0 == dval1) {
	result = TP_TRUE;
      } else {
	result = TP_FALSE;
      }
    } catch (std::exception ex) {
      return TextParserErrorHandler(TP_ILLEGAL_NUMERIC_VALUE_ERROR, value);
    }
#endif
  } else if (leaf_value->_value_type == TP_UNDEFFINED_VALUE) {
    result = TP_UNDEFINED_BOOL;
  } else {
    return TextParserErrorHandler(TP_ILLEGAL_VALUE_TYPE_ERROR, value);
  }
  if (result != TP_UNDEFINED_BOOL) {
    if (is_equal == TP_FALSE) {
      result = (result == TP_TRUE) ? TP_FALSE : TP_TRUE;   // 不等号の場合は反転
    }
  }

#ifdef TP_DEBUG
  if (_debug_write) {
    std::cout << "☆ 結果 == ";
    debugWriteTextParserBool(result);
    std::cout << std::endl;
  }
#endif
  return TP_NO_ERROR;
}

/** 論理積の解析
 *
 * @param[in] left 左辺値
 * @param[in] right 右辺値
 * @return 論理積
 *
 */
TextParserBool TextParserTree::resolveAnd(TextParserBool& left, TextParserBool& right)
{
  if (left == TP_UNDEFINED_BOOL || right == TP_UNDEFINED_BOOL) {
#ifdef TP_DEBUG
    if (_debug_write) {
      std::cout << "☆ 論理積 = 不定" << std::endl;
    }
#endif
    return TP_UNDEFINED_BOOL;
  } else if (left == TP_TRUE && right == TP_TRUE) {
#ifdef TP_DEBUG
    if (_debug_write) {
      std::cout << "☆ 論理積 = True" << std::endl;
    }
#endif
    return TP_TRUE;
  } else {
#ifdef TP_DEBUG
    if (_debug_write) {
      std::cout << "☆ 論理積 = False" << std::endl;
    }
#endif
    return TP_FALSE;
  }
}

/** 論理和の解析
 *
 * @param[in] left 左辺値
 * @param[in] right 右辺値
 * @return 論理和
 *
 */
TextParserBool TextParserTree::resolveOr(TextParserBool& left, TextParserBool& right)
{
  if (left == TP_UNDEFINED_BOOL || right == TP_UNDEFINED_BOOL) {
#ifdef TP_DEBUG
    if (_debug_write) {
      std::cout << "☆ 論理和 = 不定" << std::endl;
    }
#endif
    return TP_UNDEFINED_BOOL;
  } else if (left == TP_TRUE || right == TP_TRUE) {
#ifdef TP_DEBUG
    if (_debug_write) {
      std::cout << "☆ 論理和 = True" << std::endl;
    }
#endif
    return TP_TRUE;
  } else {
#ifdef TP_DEBUG
    if (_debug_write) {
      std::cout << "☆ 論理和 = False" << std::endl;
    }
#endif
    return TP_FALSE;
  }
}

/** エレメントの相対パス取得
 *
 * @param[in,out] path パスラベル
 * @param[in] add パスラベル
 * @param[out] parent_element 親エレメント
 * @return リターンコード
 *
 */
TextParserError TextParserTree::getElementRelativePath (std::string& path, bool add, TextParserElement **parent_element)
{
#ifdef MYDEBUG
  std::cout << "getElementRelativePath (std::string& path, bool add, TextParserElement **parent_element) start" <<std::endl;
#endif //MYDEBUG

  TextParserError ret = TP_NO_ERROR;
  TextParserNode *node;
  if (isPathLabel(path)) {                    // ラベルがパス？
    if (path[0] == '/') {                   // 絶対パス
      path = path.substr(1);
      *parent_element = 0;                // 親エレメントをルートエレメントに設定
    } else {                                // 相対パス
      if (!path.compare(0, 2, "..")) {
	path = path.substr(2);
	if (_current_element != 0 && _current_element->_type == TP_NODE_ELEMENT) {   // ルート以外のディレクトリ
	  *parent_element = _current_element->_parent;                                    // 親ディレクトリに移動
	} else {
	  return TextParserErrorHandler(TP_ILLEGAL_LABEL_PATH_ERROR, path);
	}
	if (path.size() > 0) {
	  if (path[0] == '/') {
	    path = path.substr(1);
	  } else {    // ".."の次に"/"以外が来ることは無い
	    return TextParserErrorHandler(TP_ILLEGAL_LABEL_PATH_ERROR, path);
	  }
	}
      } else {
	if (!path.compare(0, 2, "./")) {        // "./"は無くて良い
	  path = path.substr(2);
	}
	*parent_element = _current_element;     // 親エレメントをカレントエレメントに設定
      }
    }
    while (path.size() > 0) {        // パスをラベルに分解して対応するディレクトリが無ければ作成する
      int s = path.find("/");
      if (s < 0) break;
      std::string label = path.substr(0, s);
      if (label.size() == 0 || !isCorrectLabel(label, false)) {
	return TextParserErrorHandler(TP_ILLEGAL_LABEL_PATH_ERROR, label);
      }
      if (add) {
	ret = getOrAddNode(label, *parent_element, &node);// ディレクトリの取得又は追加
      } else {
	ret = getNode(label, *parent_element, &node);     // ディレクトリの取得又は追加
      }
      if (ret != TP_NO_ERROR) return ret;
      *parent_element = node;
      path = path.substr(s + 1);
    }
  } else {
    *parent_element = _current_element;				        // 親エレメントをカレントエレメントに設定
  }

#ifdef MYDEBUG
  std::cout << "getElementRelativePath (std::string& path, bool add, TextParserElement **parent_element) end" <<std::endl;
#endif //MYDEBUG


  return ret;
}

/** エレメントの取得
 *
 * @param[in] label パスラベル
 * @param[in] type エレメントタイプ
 * @param[out] element エレメント
 * @return リターンコード
 *
 */
TextParserError TextParserTree::getElement(const std::string& label, const TextParserElementType type, TextParserElement** element)
{
  TextParserError ret = TP_NO_ERROR;
  *element = 0;
  //  std::cout << "getElement 0" << label <<std::endl;
  if (type != TP_NODE_ELEMENT && type != TP_LEAF_ELEMENT) {
    return TextParserErrorHandler(TP_ILLEGAL_ELEMENT_ERROR, label);
  }
  //  std::cout << "getElement 1" << label <<std::endl;
  TextParserElement *parent_element = 0;
  std::string path = label;
  ret = getElementRelativePath (path, false, &parent_element);
  //  std::cout << "getElement 2" << label <<std::endl;
  if (ret != TP_NO_ERROR) return ret;
  //  std::cout << "getElement 3" << label <<std::endl;
  if (path.size() == 0) { // "/"や".."の場合
    if (type == TP_NODE_ELEMENT) {
      *element = 0;
      return TP_NO_ERROR;
    } else {
      return TP_MISSING_PATH_ELEMENT_ERROR;
    }
  }
  //  std::cout << "getElement 4" << label <<std::endl;
  TextParserNode *parent_dir = 0;
  if (parent_element == 0) {  // ルートディレクトリ
    *element = getNode(path);
    //  std::cout << "getElement 5" << label <<std::endl;
    if (*element == 0) {
      //  std::cout << "getElement 6" << label <<std::endl;
      return TP_MISSING_PATH_ELEMENT_ERROR;
    }
  } else if (parent_element->_type == TP_NODE_ELEMENT) {
    parent_dir = (TextParserNode *)parent_element;
    if (type == TP_NODE_ELEMENT) {
      //  std::cout << "getElement 7" << label <<std::endl;
      *element = parent_dir->getNode(path);
    } else if (type == TP_LEAF_ELEMENT) {
      //      std::cout << "getElement 8" << label <<std::endl;
      *element = parent_dir->getLeaf(path);
      //      std::cout << "getElement 8" << path << *element <<std::endl;
    }
    if (*element == 0) {
      //      std::cout << "getElement 9" << label <<std::endl;
      return TP_MISSING_PATH_ELEMENT_ERROR;
    }
  } else {
    //  std::cout << "getElement 10" << label <<std::endl;
    return TextParserErrorHandler(TP_ILLEGAL_PATH_ELEMENT_ERROR, "");
  }
  //  std::cout << "getElement 11" << label <<std::endl;
  //  std::cout << "getElement " << path << " " <<(*element) << std::endl;
  return ret;
}





/** リーフの値を取得
 *
 * @param[in] path リーフのパス
 * @param[out] value 値
 * @return リターンコード
 *
 */
TextParserError TextParserTree::getLeafValue(const  std::string& path, TextParserValue **value)
//TextParserError TextParserTree::getLeafValue(cosnt std::string* path, TextParserValue **value)
{
  TextParserError ret;
  TextParserElement* element = 0;
  ret = getElement(path, TP_LEAF_ELEMENT, &element); // パスからエレメントを取得
  if (ret != TP_NO_ERROR) return ret;
  if (element->_type == TP_LEAF_ELEMENT) {
    TextParserLeaf *leaf = (TextParserLeaf *)element;
    if (leaf->_value == 0) {
      return TextParserErrorHandler(TP_ILLEGAL_VALUE_ERROR, "");
    } else {
      *value = leaf->_value;
    }
  } else {
    return TextParserErrorHandler(TP_ILLEGAL_PATH_ELEMENT_ERROR, "");
  }

#ifdef TP_DEBUG
  if (_debug_write) {
    std::cout << "リーフの値（" + path + " = " + (*value)->_value + ")" << std::endl;
  }
#endif
  return TP_NO_ERROR;
}

/** パラメータファイルの書き出し
 *
 * @param[in] filename ファイル名
 * @return エラーコード
 *
 */
TextParserError TextParserTree::writeParameters(const std::string& filename,int order)
{
#ifdef MYDEBUG
  std::cout<<"writeParameters(const std::string& filename) start"<<std::endl;
#endif //MYDEBUG

  TextParserError ret = TP_NO_ERROR;
  if (!_is_ready) {
    return TextParserErrorHandler(TP_DATABASE_NOT_READY_ERROR, "");
  }

  // 出力ファイルのオープン
  std::ostream *ofs;
  std::ofstream *fp = 0;
  if (!filename.compare("stdout")) {
    ofs = &std::cout;
  } else {
    try {
      fp = new std::ofstream(filename.c_str());
      if (fp == 0) {
	return TextParserErrorHandler(TP_FILEOPEN_ERROR, filename);
      }
      ofs = fp;
    } catch (std::exception ex) {
      ret = TextParserErrorHandler(TP_FILEOPEN_ERROR, filename);
    }
  }
  if (!(*ofs)) {
    ret = TextParserErrorHandler(TP_FILEOPEN_ERROR, filename);
  }
  if (ret != TP_NO_ERROR) return ret;

  // パラメータの出力
  try {
    std::map<std::string, TextParserNode *>::iterator di = _nodes.begin();
    while (di != _nodes.end()) {
      TextParserNode *dir = (di)->second;
      ret = dir->writeNode(*ofs, 0,order);
      if (ret != TP_NO_ERROR) return ret;
      di++;
    }
  } catch (std::ios_base::failure) {
    ret = TextParserErrorHandler(TP_FILEOUTPUT_ERROR, filename);
  } catch (std::exception ex) {
    ret = TextParserErrorHandler(TP_FILEOUTPUT_ERROR, filename);
  }

  if (fp != 0) delete fp;


#ifdef MYDEBUG
  std::cout<<"writeParameters(const std::string& filename) endl"<<std::endl;
#endif //MYDEBUG
  return ret;
}

/** リーフのラベル取得
 *
 * @param[in,out] labels リーフのラベル
 * @return エラーコード
 *
 */
TextParserError TextParserTree::getLeafLabels(std::map<unsigned int, std::string>& labels)
{
  TextParserError ret;

  std::map<std::string, TextParserNode *>::iterator di = _nodes.begin();
  while (di != _nodes.end()) {
    TextParserNode *dir = (di)->second;
    ret = dir->getLeafLabels(labels);
    if (ret != TP_NO_ERROR) return ret;
    di++;
  }

  return TP_NO_ERROR;
}

/** ディレクトリの移動
 *
 * @param[in] label ディレクトリのラベルパス
 * @return エラーコード
 *
 */
TextParserError TextParserTree::changeNode(const std::string& label)
{
  TextParserError ret;
  TextParserNode *node = 0;
  TextParserElement *parent_element = 0;
  std::string path = label;
  ret = getElementRelativePath (path, false, &parent_element);           // パスで指定したエレメントの相対パスを取得
  if (ret != TP_NO_ERROR) return ret;
  if (path.size() == 0) { // "/"や".."の場合
    _current_element = parent_element;                          // パスが""の時は親のエレメントをカレントに設定（".."のケース）
  } else {
    ret = getNode(path, parent_element, &node);       // ディレクトリの取得又は追加
    if (ret != TP_NO_ERROR) return ret;
    _current_element = node;                               // カレントエレメントに設定
  }

#ifdef TP_DEBUG
  if (_debug_write) {
    if (_current_element == 0) {
      std::cout << "ディレクトリの移動 : /" << std::endl;
    } else {
      std::cout << "ディレクトリの移動 : " + _current_element->_label << std::endl;
    }
  }
#endif
  return TP_NO_ERROR;
}

/** 現在のディレクトリの取得
 *
 * @param[out] label ディレクトリの絶対パス
 * @return エラーコード
 *
 */
TextParserError TextParserTree::getCurrentNode(std::string& label)
{
  TextParserError ret;
  if (_current_element == 0) {
    label = "/";
  } else {
    if (_current_element->_type != TP_NODE_ELEMENT) {
      return TextParserErrorHandler(TP_ILLEGAL_CURRENT_ELEMENT_ERROR, "");
    }
    ret = GetElementAbsolutePath(_current_element, label);
  }

  return TP_NO_ERROR;
}

/** 現在のディレクトリ内のディレクトリ数を取得
 *
 * @param[out] number ディレクトリ数
 * @return エラーコード
 *
 */
TextParserError TextParserTree::getCurrentNodeNumber(unsigned int *number)
{
  if (_current_element == 0) {
    *number = _nodes.size();
  } else {
    if (_current_element->_type != TP_NODE_ELEMENT) {
      return TextParserErrorHandler(TP_ILLEGAL_CURRENT_ELEMENT_ERROR, "");
    }
    TextParserNode *current_dir = (TextParserNode *)_current_element;
    *number = current_dir->_nodes.size();
  }

  return TP_NO_ERROR;
}

/** 現在のディレクトリ内のディレクトリのラベルを取得
 *
 * @param[in] id ディレクトリのID
 * @param[out] label ディレクトリのラベル（ポインタ）
 * @return エラーコード
 *
 */
TextParserError TextParserTree::getCurrentNodeLabel(unsigned int id, char **label)
{
  if (_current_element == 0) {
    std::map<std::string, TextParserNode *>::iterator di = _nodes.begin();
    for (unsigned int i = 0; i < id; i++) di++;
    if (di != _nodes.end()) {
      TextParserNode *node = di->second;
      *label = (char *)node->_label.c_str();
    } else {
      return TextParserErrorHandler(TP_ID_OVER_ELEMENT_NUMBER_ERROR, "");
    }
  } else {
    if (_current_element->_type != TP_NODE_ELEMENT) {
      return TextParserErrorHandler(TP_ILLEGAL_CURRENT_ELEMENT_ERROR, "");
    }
    TextParserNode *current_dir = (TextParserNode *)_current_element;
    std::map<std::string, TextParserNode *>::iterator di = current_dir->_nodes.begin();
    for (unsigned int i = 0; i < id; i++) di++;
    if (di != current_dir->_nodes.end()) {
      TextParserNode *node = di->second;
      *label = (char *)node->_label.c_str();
    } else {
      return TextParserErrorHandler(TP_ID_OVER_ELEMENT_NUMBER_ERROR, "");
    }
  }

  return TP_NO_ERROR;
}

/** 現在のディレクトリ内の全ディレクトリのラベルを取得
 *
 * @param[out] labels ディレクトリのラベル
 * @return エラーコード
 *
 */
TextParserError TextParserTree::getCurrentNodeLabels(std::vector<std::string>& labels)
{
  if (_current_element == 0) {
    std::map<std::string, TextParserNode *>::iterator di = _nodes.begin();
    while (di != _nodes.end()) {
      TextParserNode *node = di->second;
      labels.push_back(node->_label);
      di++;
    }
  } else {
    if (_current_element->_type != TP_NODE_ELEMENT) {
      return TextParserErrorHandler(TP_ILLEGAL_CURRENT_ELEMENT_ERROR, "");
    }
    TextParserNode *current_dir = (TextParserNode *)_current_element;
    std::map<std::string, TextParserNode *>::iterator di = current_dir->_nodes.begin();
    while (di != current_dir->_nodes.end()) {
      TextParserNode *node = di->second;
      labels.push_back(node->_label);
      di++;
    }
  }

  return TP_NO_ERROR;
}

/** 現在のディレクトリ内のリーフ数を取得
 *
 * @param[out] number リーフ数
 * @return エラーコード
 *
 */
TextParserError TextParserTree::getCurrentLeafNumber(unsigned int *number)
{
  if (_current_element == 0) {
    *number = 0;
  } else {
    if (_current_element->_type != TP_NODE_ELEMENT) {
      return TextParserErrorHandler(TP_ILLEGAL_CURRENT_ELEMENT_ERROR, "");
    }
    TextParserNode *current_dir = (TextParserNode *)_current_element;
    *number = current_dir->_leaves.size();
  }

  return TP_NO_ERROR;
}

/** 現在のディレクトリ内のリーフのラベルを取得
 *
 * @param[in] id リーフのID
 * @param[out] label リーフのラベル（ポインタ）
 * @return エラーコード
 *
 */
TextParserError TextParserTree::getCurrentLeafLabel(unsigned int id, char **label)
{
  if (_current_element == 0) {
    return TextParserErrorHandler(TP_ILLEGAL_CURRENT_ELEMENT_ERROR, "");
  } else {
    if (_current_element->_type != TP_NODE_ELEMENT) {
      return TextParserErrorHandler(TP_ILLEGAL_CURRENT_ELEMENT_ERROR, "");
    }
    TextParserNode *current_dir = (TextParserNode *)_current_element;
    std::map<std::string, TextParserLeaf *>::iterator li = current_dir->_leaves.begin();
    for (unsigned int i = 0; i < id; i++) li++;
    if (li != current_dir->_leaves.end()) {
      TextParserLeaf *leaf = li->second;
      *label = (char *)leaf->_label.c_str();
    } else {
      return TextParserErrorHandler(TP_ID_OVER_ELEMENT_NUMBER_ERROR, "");
    }
  }

  return TP_NO_ERROR;
}

/** 現在のディレクトリ内の全リーフのラベルを取得
 *
 * @param[out] labels リーフのラベル
 * @return エラーコード
 *
 */
TextParserError TextParserTree::getCurrentLeafLabels(std::vector<std::string>& labels)
{
  if (_current_element == 0) {
    labels.clear();
  } else {
    if (_current_element->_type != TP_NODE_ELEMENT) {
      return TextParserErrorHandler(TP_ILLEGAL_CURRENT_ELEMENT_ERROR, "");
    }
    TextParserNode *current_dir = (TextParserNode *)_current_element;
    std::map<std::string, TextParserLeaf *>::iterator li = current_dir->_leaves.begin();
    while (li != current_dir->_leaves.end()) {
      TextParserLeaf *leaf = li->second;
      labels.push_back(leaf->_label);
      li++;
    }
  }

  return TP_NO_ERROR;
}

/** ベクトル値を分離
 * @param[in] vector_value ベクトル値.
 * @param[out] values 分離された値の文字列の配列.
 * @return エラーコード.
 */
TextParserError TextParserTree::splitVectorValue(const std::string &vector_value, std::vector<std::string>& values)
{
  TextParserError ret;
  bool answer = false;
  
  std::string buffer = vector_value;
  ret = parseVectorValue(buffer, answer);
  if (ret != TP_NO_ERROR) return ret;
  TextParserValueType value_type;
  std::string value;
  if (answer) {             // ベクトル値?
    while(true) {
      // 値の取得
      ret = parseUndefinedValue(buffer, answer);
      if (ret != TP_NO_ERROR) return ret;
      if (answer) {          // 未定義の値?
	values.push_back("UNDEF");
      } else {
	ret = parseValue(buffer, answer, value, value_type);
	if (ret != TP_NO_ERROR) return ret;
	if (answer) {       // スカラー値?
	  values.push_back(value);
	} else {
	  return TextParserErrorHandler(TP_ILLEGAL_VALUE_ERROR, buffer);
	}
      }
      ret = parseDelimiter(buffer, answer);
      if (answer) {   // デリミタ?
      } else {
	ret = parseEndOfVector(buffer, answer);
	if (ret != TP_NO_ERROR) return ret;
	if (answer) {   // ベクトルの終わり
	  break;
	} else {
	  return TextParserErrorHandler(TP_MISSING_VECTOR_END_ERROR, buffer);
	}
      }
    }
  } else {
    ret = TextParserErrorHandler(TP_ILLEGAL_VALUE_ERROR, buffer);
  }

  return TP_NO_ERROR;
}




/** デバッグ文の表示スイッチ設定
 *
 * @param[in] swt 表示スイッチ
 *
 */
void TextParserTree::debugWrite(bool swt)
{
  _debug_write = swt;
}

/** 値のデバッグ表示
 *
 * @param[in] value 値
 * @param[in] value_type 値のタイプ
 *
 */
void TextParserTree::debugWriteValue(const std::string& value, const TextParserValueType& value_type)
{
#ifdef TP_DEBUG
  if (_debug_write) {
    if (value_type == TP_STRING_VALUE) std::cout << "\"";
    std::cout << value;
    if (value_type == TP_STRING_VALUE) std::cout << "\"";
  }
#endif
}

/** 論理値のデバッグ表示
 *
 * @param[in] result 論理値の
 *
 */
void TextParserTree::debugWriteTextParserBool(const TextParserBool& result)
{
#ifdef TP_DEBUG
  if (_debug_write) {
    if (result == TP_UNDEFINED_BOOL) {
      std::cout << "不定";
    } else if (result == TP_TRUE) {
      std::cout << "True";
    } else if (result == TP_FALSE) {
      std::cout << "False";
    }
  }
#endif
}




/** リーフラベルのソート
 * @param[in] input 入力ラベル
 * @param[out] output ソートしたラベル
 * @param[in] iswitch ソートスイッチ 0:何もしない 1:配列ラベルのインデックス順 2:パラメータファイルの記述順
 * @return エラーコード
 */

TextParserError TextParserTree::labelSort(const std::vector<std::string>& input,
					  std::vector<std::string>& output,
					  int iswitch)
{
  TextParserError ret =TP_NO_ERROR;
  //  std::cout << "labelSort "<< iswitch<< std::endl;  
  if(iswitch==0){
    output=input;
    return ret;
  } else if(iswitch==1) {
    ret=labelSort_1(input,output);
    return ret;
  } else if(iswitch==2) {
    ret=labelSort_2(input,output);

    return ret;
  }


return ret ;
}


/** リーフラベルのソート. 配列ラベルをインデックス順にソートする.
 * @param[in] input 入力ラベルリスト
 * @param[out] output ソートしたラベルリスト
 * @return エラーコード
 */


TextParserError
TextParserTree::labelSort_1(const std::vector<std::string>& input,
			    std::vector<std::string>& output)
{
  TextParserError ret=TP_NO_ERROR;
  
  std::vector<TextParserElement*> velement; 
  std::string key;
  std::vector<int> order_switch;
  std::vector<std::string> key_list;
  std::vector<int> key_number;
  std::vector< std::map<int,std::string> > map_buffer;


  std::vector< std::string >::const_iterator iter=input.begin();
  while ( iter !=input.end() ){
    std::string label=(*iter);
    // std::cout <<label <<std::endl;
    int number=array_label_test(label,key);
    
    if( number==-10000 || number == -1000){
      return TextParserErrorHandler(TP_ILLEGAL_LABEL_PATH_ERROR,label);
    }
    
    if (number < 0){
      order_switch.push_back(0);
    } else {
      std::vector<std::string>::iterator
	kl_iter=find(key_list.begin(),key_list.end(),key);
      if(kl_iter != key_list.end()){
	int index=kl_iter-key_list.begin();
	order_switch.push_back(key_number[index]);
	map_buffer[order_switch.back()-1].
	  insert(std::map<int,std::string>::value_type(number,label) );
      } else {

	key_list.push_back(key);
	key_number.push_back(key_list.size());
	order_switch.push_back(key_list.size());
	//std::cout << "new key "<< key <<" "<<key_list.size()<<std::endl;
	
	std::map<int,std::string> dummy;
	dummy.insert(std::map<int,std::string>::value_type(number,label) );

	map_buffer.push_back(dummy);
      }
    }
    
    iter++;
  }

   std::vector< std::map<int,std::string>::iterator > map_buffer_iter;
   std::vector< std::map<int,std::string>::iterator > map_buffer_iter_end;

   for (int k=0;k<key_list.size();k++){
     map_buffer_iter.push_back(map_buffer[k].begin());
     map_buffer_iter_end.push_back(map_buffer[k].end());
   }
   
   iter=input.begin();
   int i=0;
   
   while ( iter !=input.end() ){
     
     std::string label=(*iter);
     if(order_switch[i]==0) {
       output.push_back(label);
     } else {
       //std::cout<<"order_switch " << order_switch[i] << std::endl;
       if(map_buffer_iter[order_switch[i]-1]
	  !=map_buffer_iter_end[order_switch[i]-1]){
	 //	 std::cout << "aaa " <<  std::endl; 
	 output.push_back(map_buffer_iter[order_switch[i]-1]->second);
	 //	 std::cout << "bbb " <<  std::endl; 
	 map_buffer_iter[order_switch[i]-1]++;
       } else {
	 return TextParserErrorHandler(TP_ERROR,label);
       }
     }
     ++iter;
     ++i;
   }

   //std::cout << "labelSort_1 end" <<std::endl;

 return ret;

}


/** リーフラベルのソート. ラベルをパラメータファイルでの出現順に並べる.
 * @param[in] input 入力ラベルリスト
 * @param[out] output ソートしたラベルリスト
 * @return エラーコード
 */

TextParserError
TextParserTree::labelSort_2(const std::vector<std::string>& input,
			    std::vector<std::string>& output)
{
  TextParserError ret=TP_NO_ERROR;
  
  std::map<int,std::string>  map_buffer;


  std::vector< std::string >::const_iterator iter=input.begin();
  while ( iter !=input.end() ){
    std::string label=(*iter);
    TextParserElement* element=0; 
    //    std::cout << label  <<std::endl; 
    ret=getElement(label,TP_LEAF_ELEMENT,&element);
    TextParserLeaf* leaf = (TextParserLeaf*) element; 
    if (ret!=TP_NO_ERROR) {
      return TextParserErrorHandler(ret,label);
    }    
    //    std::cout << label  << " "<< leaf->_id<<std::endl; 
    map_buffer.insert(std::map<int,std::string>::value_type(leaf->_id,label) );
    
    iter++;
  }
  std::map<int,std::string>::iterator map_iter=map_buffer.begin();
  while( map_iter != map_buffer.end()){
    output.push_back(map_iter->second);
    ++map_iter;
  }



return ret;
}



/** check array-ed label number.
 * @param[in] label label.
 * @param[out] key label string without number.
 * @return label id.
 */

int TextParserTree::array_label_test(const std::string& label,std::string& key){

  int number=-10000 ;
  key="";

  int lsb=label.find("[");
  if (lsb<0){
    return lsb;
    //  output.push_back(label);
  } else if (lsb == 0) { 
    TextParserErrorHandler(TP_ILLEGAL_LABEL_PATH_ERROR,label);
    return lsb;
  } else {
    int rsb=label.find("]");
    if(rsb<0 || lsb >= rsb ){
      TextParserErrorHandler(TP_ILLEGAL_LABEL_PATH_ERROR,label); 
      return -1000;
    } else {
      key=label.substr(0,lsb);
      std::string num_string=label.substr(lsb+1,rsb-lsb-1);
      std::stringstream ts;
      ts << num_string;
      
      ts >>number;
       // std::cout << "array_label_test lsb rsb "<< lsb << " "<< rsb 
       // 		 << " key " <<key 
       // 		<< " num_string "<< num_string 
       // 		<<" number "<< number<<std::endl;
    }
  }
  
  return number;
}


/** ノードラベルのソート
 * @param[in] input 入力ラベル
 * @param[out] output ソートしたラベル
 * @param[in] iswitch ソートスイッチ 0:何もしない 1:配列ラベルのインデックス順 2:パラメータファイルの記述順
 * @return エラーコード
 */

TextParserError TextParserTree::nodeSort(const std::vector<std::string>& input,
					  std::vector<std::string>& output,
					  int iswitch)
{
  TextParserError ret =TP_NO_ERROR;
  //  std::cout << "labelSort "<< iswitch<< std::endl;  
  if(iswitch==0){
    output=input;
    return ret;
  } else if(iswitch==1) {
    // same as leaf label
    ret=labelSort_1(input,output);
    return ret;
  } else if(iswitch==2) {
    ret=nodeSort_2(input,output);

    return ret;
  }


return ret ;
}



/** ノードラベルのソート. ラベルをパラメータファイルでの出現順に並べる.
 * @param[in] input 入力ラベルリスト
 * @param[out] output ソートしたラベルリスト
 * @return エラーコード
 */

TextParserError
TextParserTree::nodeSort_2(const std::vector<std::string>& input,
			    std::vector<std::string>& output)
{
  TextParserError ret=TP_NO_ERROR;
  
  std::multimap<int,std::string>  map_buffer;


  std::vector< std::string >::const_iterator iter=input.begin();
  while ( iter !=input.end() ){
    std::string label=(*iter);
    TextParserElement* element=0; 
    //    std::cout << label  <<std::endl; 
    ret=getElement(label,TP_NODE_ELEMENT,&element);
    TextParserNode* node = (TextParserNode*) element; 
    if (ret!=TP_NO_ERROR) {
      return TextParserErrorHandler(ret,label);
    }    
    //    std::cout << label  << " "<< leaf->_id<<std::endl; 
    map_buffer.
      insert(std::multimap<int,std::string>::value_type(node->line(),label) );
    
    iter++;
  }
  std::multimap<int,std::string>::iterator map_iter=map_buffer.begin();
  while( map_iter != map_buffer.end()){
    output.push_back(map_iter->second);
    ++map_iter;
  }

  return ret;
}
