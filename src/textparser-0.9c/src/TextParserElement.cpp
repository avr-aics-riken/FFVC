/****************************************************************************
**
** Copyright (C) 2012 Tokyo University.
**
****************************************************************************/
/** @file TextParserElement.cpp
 * ここには TextParserElement クラス及びその派生クラス TextParserNode ,
 * TextParserLeaf が実装されています。
 *
 */

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif //HAVE_CONFIG_H

#include <exception>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <utility>
#include "TextParserCommon.h"
#include "TextParserElement.h"
#include "TextParserTree.h"

#include "../config.h"

/** TextParserElement のコンストラクタ
 *
 */
TextParserElement::TextParserElement()
{
    _label = "";
    _parent = 0;
    _type = TP_UNKNOWN_ELEMENT;
}

/** TextParserNodeのコンストラクタ
 *
 * @param[in] label ディレクトリのラベル
 *
 * ラベルは小文字に変換して保存
 *
 */
TextParserNode::TextParserNode(const std::string label)
{
    _label = label;
    _parent = 0;
    _type = TP_NODE_ELEMENT;
}

/** エレメントの追加
 *
 * @param[in] node TextParserNodeのポインタ
 * @return エラーコード
 *
 * ディレクトリエレメントを追加します
 *
 */
TextParserError TextParserNode::addElement(TextParserNode *node)
{
    std::string label = TextParserStringToLower(node->_label);
    _nodes.insert(str_node(label, node));

    return TP_NO_ERROR;
}

/** エレメントの追加
 *
 * @param[in] leaf TextParserLeafのポインタ
 * @return エラーコード
 *
 * リーフエレメントを追加します
 *
 */
TextParserError TextParserNode::addElement(TextParserLeaf *leaf)
{
    std::string label = TextParserStringToLower(leaf->_label);
    _leaves.insert(str_leaf(label, leaf));
    leaf->_id = TextParserTree::GetLeafID();

    return TP_NO_ERROR;
}

/** エレメントの削除
 *
 * @return エラーコード
 *
 */
TextParserError TextParserNode::removeElement()
{
    TextParserError ret;
    try {
        std::map<std::string, TextParserNode *>::iterator ni = _nodes.begin();
        while (ni != _nodes.end()) {
            TextParserNode *node = (ni)->second;
            ret = node->removeElement();
            if (ret != TP_NO_ERROR) return ret;
            delete node;
            ni++;
        }
        _nodes.clear();

        std::map<std::string, TextParserLeaf *>::iterator li = _leaves.begin();
        while (li != _leaves.end()) {
            TextParserLeaf *leaf = (li)->second;
            ret = leaf->removeElement();
            if (ret != TP_NO_ERROR) return ret;
            delete leaf;
            li++;
        }
        _leaves.clear();

        _array_label_number.clear();
    } catch (std::exception ex) {
        return TP_REMOVE_ELEMENT_ERROR;
    }

    return TP_NO_ERROR;
}

/** ディレクトリの取得
 *
 * @param[in] label ラベル
 * @return ディレクトリ
 *
 * ラベルの一致するディレクトリを取得
 *
 */
TextParserNode *TextParserNode::getNode(const std::string& label)
{
    std::string label_cpy = TextParserStringToLower(label);
    std::map<std::string, TextParserNode *>::iterator di = _nodes.find(label_cpy);
    if (di != _nodes.end()) {
        return di->second;
    } else {
        if (!label.compare("..")) {
            return (TextParserNode *)_parent;
        } else {
            return 0;
        }
    }
}

/** リーフの取得
 *
 * @param[in] label ラベル
 * @return リーフ
 *
 * ラベルの一致するリーフを取得
 *
 */
TextParserLeaf *TextParserNode::getLeaf(const std::string& label)
{
    std::string label_cpy = TextParserStringToLower(label);
    //    std::cout << "getLeaf " << label << " "<< label_cpy << std::endl; 

    std::map<std::string, TextParserLeaf *>::iterator li = _leaves.find(label_cpy);
    if (li != _leaves.end()) {
        return li->second;
    } else {
        return 0;
    }
}

/** 配列ラベルのインデックス設定
 *
 * @param[in] label ラベル
 * @return リターンコード
 *
 * 配列ラベル（末尾が[@]）だったらインデックスを割り当てる
 *
 */
TextParserError TextParserNode::setArrayLabelIndex(std::string& label)
{
    return TextParserTree::SetArrayLabelIndex(label, _array_label_number);
}

/** ディレクトリの書き出し
 *
 * @param[in] ofs 出力ファイルポインタ
 * @param[in] level 階層のレベル
 * @return エラーコード
 *
 */
TextParserError TextParserNode::writeNode(std::ostream& ofs, unsigned int level,int order)
{
    for (unsigned int l = 0; l < level; l++) ofs << TP_INDENT_STRING;
#ifdef TP_NODQOUT    
        ofs <<  _label + " {\n";
#else 
        ofs << "\"" + _label + "\" {\n";
#endif //TP_NODQOUT    
	if(order!=0){
	  // TextParserTree* tree_instance=TextParserTree::get_instance();
	{
	  std::vector<std::string> tmplabel1;
	  std::vector<std::string> tmplabel2;
	  std::map<std::string, TextParserNode *>::iterator di = _nodes.begin();
	  
	  while (di != _nodes.end()) {
	    TextParserNode *dir = di->second;
	    tmplabel1.push_back(dir->_label);
	    di++;
	  }
	  //	  tree_instance->nodeSort(tmplabel1,tmplabel2,order);
	  element_node_sort(tmplabel1,tmplabel2,order);
	  std::vector<std::string>::iterator tmp_iter=tmplabel2.begin();
	  while(tmp_iter !=tmplabel2.end()){
	    TextParserNode *dir = _nodes[TextParserStringToLower(*tmp_iter)];
	    //std::cout << __func__ << " 1 "<<*tmp_iter << std::endl;
	    //std::cout << __func__ << " 2 "<<dir->_label << std::endl;
	    dir->writeNode(ofs, level + 1,order);
	    tmp_iter++;
	  }
	}
	
	{
	  std::vector<std::string> tmplabel1;
	  std::vector<std::string> tmplabel2;
	  
	  std::map<std::string, TextParserLeaf *>::iterator
	    li = _leaves.begin();
	  while (li != _leaves.end()) {
	    TextParserLeaf *leaf = li->second;
	    // std::cout << __func__<< __LINE__
	    // 	      <<" "<<leaf 
	    // 	      << " label " << leaf->_label <<std::endl;
	    tmplabel1.push_back(leaf->_label);
	    li++;
	  }
	  // tree_instance->labelSort(tmplabel1,tmplabel2,order);
	  element_label_sort(tmplabel1,tmplabel2,order);

	  std::vector<std::string>::iterator tmp_iter=tmplabel2.begin();
	  while(tmp_iter !=tmplabel2.end()){
	    
	    //	    std::string tmpstring = *tmp_iter;

	    TextParserLeaf 
	      *tmpleaf2 =_leaves[TextParserStringToLower(*tmp_iter)];
	    if (tmpleaf2==0){	    
	      return TextParserErrorHandler(TP_ERROR,*tmp_iter);
	    }
	    tmpleaf2->writeLeaf(ofs, level + 1);

	    tmp_iter++;
	  }
	  
	}

	} else {
	  
	  // order ==0 original code.
	  std::map<std::string, TextParserNode *>::iterator di = _nodes.begin();
	  while (di != _nodes.end()) {
	    TextParserNode *dir = di->second;
	    dir->writeNode(ofs, level + 1);
	    di++;
	  }
	  std::map<std::string, TextParserLeaf *>::iterator li = _leaves.begin();
	  while (li != _leaves.end()) {
	    TextParserLeaf *leaf = li->second;
	    std::cout << __func__ << " 5 leaf "<<leaf->_label << std::endl;
	    leaf->writeLeaf(ofs, level + 1);
	    li++;
	  }
	}
	
	for (unsigned int l = 0; l < level; l++) ofs << TP_INDENT_STRING;
	ofs << "}\n";
	
	return TP_NO_ERROR;
}

/** リーフのラベル取得
 *
 * @param[in,out] labels リーフのラベル
 * @return エラーコード
 *
 */
TextParserError TextParserNode::getLeafLabels(std::map<unsigned int, std::string>& labels)
{
    TextParserError ret;
 
	std::map<std::string, TextParserNode *>::iterator di = _nodes.begin();
	while (di != _nodes.end()) {
        TextParserNode *dir = di->second;
        ret = dir->getLeafLabels(labels);
        if (ret != TP_NO_ERROR) return ret;
        di++;
	}
	std::map<std::string, TextParserLeaf *>::iterator li = _leaves.begin();
	while (li != _leaves.end()) {
        TextParserLeaf *leaf = li->second;
        std::string path;
        TextParserTree::GetElementAbsolutePath(leaf, path);
        labels.insert(int_str(leaf->_id, path));
        li++;
    }
    

    return TP_NO_ERROR;
}

/** TextParserLeafのコンストラクタ
 *
 * @param[in] label リーフのラベル
 *
 */
TextParserLeaf::TextParserLeaf(const std::string label)
{
    _label = label;
    _parent = 0;
    _type = TP_LEAF_ELEMENT;
    _value_type = TP_UNDEFFINED_VALUE;
    _value = 0;
}

/** エレメントの設定
 *
 * @param[in] value TextParserValueのポインタ
 * @return エラーコード
 *
 * 値エレメントを設定します
 *
 */
TextParserError TextParserLeaf::setElement(TextParserValue *value)
{
    if (_value_type == TP_UNDEFFINED_VALUE) { // 値のタイプが未定義
        _value_type = value->_value_type;       // 値のタイプをコピー
    }
    _value = value;

    return TP_NO_ERROR;
}

/** エレメントの削除
 *
 * @return エラーコード
 *
 */
TextParserError TextParserLeaf::removeElement()
{
    try {
        if (_value != 0) {
            delete _value;
        }
    } catch (std::exception ex) {
        return TP_REMOVE_ELEMENT_ERROR;
    }

    return TP_NO_ERROR;
}

/** リーフの書き出し
 *
 * @param[in] ofs 出力ファイルポインタ
 * @param[in] level 階層のレベル
 * @return エラーコード
 *
 */
TextParserError TextParserLeaf::writeLeaf(std::ostream& ofs, unsigned int level)
{
    for (unsigned int l = 0; l < level; l++) ofs << TP_INDENT_STRING;

    if (_value == 0) {
    } else {
        std::string val = _value->_value;
        if (_value->_value_type == TP_STRING_VALUE) {
            val = "\"" + val + "\"";                    // 文字列の場合は""で囲む
        }
#ifdef TP_NODQOUT
        ofs <<  _label + " = " + val + "\n";
#else   //TP_NODQOUT
        ofs << "\"" + _label + "\" = " + val + "\n";
#endif  //TP_NODQOUT

    }

    return TP_NO_ERROR;
}

/** TextParserValueのコンストラクタ
 *
 * @param[in] value 値の文字列
 * @param[in] type 値のタイプ
 *
 */
TextParserValue::TextParserValue(const std::string& value, const TextParserValueType& type)
{
    _value = value;
    _value_type = type;
    _type = TP_VALUE_ELEMENT;
}

TextParserError element_node_sort(const std::vector<std::string>& input,
				  std::vector<std::string>& output,
		       int order)
{
 TextParserError ret =TP_NO_ERROR;
  if(order==0){
    output=input;
    return ret;
  } else if(order==1) {
    ret=element_labelSort_1(input,output);
    return ret;
  } else if(order==2) {
    //    ret=element_nodeSort_2(input,output);
    return ret;
  }
  return TP_ERROR ;
}

TextParserError element_label_sort(const std::vector<std::string>& input,
				  std::vector<std::string>& output,
		       int order)
{
 TextParserError ret =TP_NO_ERROR;
  if(order==0){
    output=input;
    return ret;
  } else if(order==1) {
    ret=element_labelSort_1(input,output);
    return ret;
  } else if(order==2) {
    //    ret=element_labelSort_2(input,output);
    return ret;
  }
  return TP_ERROR ;
}

TextParserError element_labelSort_1(const std::vector<std::string>& input,
				    std::vector<std::string>& output){

  std::cout << "element_labelSort_1" <<std::endl;
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
    int number=element_array_label_test(label,key);
    
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
   std::cout << "element_labelSort_1 end" <<std::endl;
 return ret;

}


/** check array-ed label number.
 * @param[in] label label.
 * @param[out] key label string without number.
 * @return label id.
 */

int element_array_label_test(const std::string& label,std::string& key){

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

