 /*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2012
 *
 */

//@file TPControl.C
//@brief FlowBase Text Parser Control class
//@author keno, FSI Team, VCAD, RIKEN

#include "TPControl.h"


/**
 @fn int TPControl::getTPinstance()
 @brief TextParserLibraryのインスタンス生成
 */
int TPControl::getTPinstance()
{
  int ierr = TP_NO_ERROR;

  tp=TextParser::get_instance();
  if( !tp )
  {
    cerr << "ERROR : instance of TextParser" << endl;
    return TP_ERROR;
  }
  return ierr;
}

/**
 @fn int TPControl::readTPfile( const string filename )
 @brief TextParserオブジェクトに入力ファイルをセットする
 @param filename 入力ファイル名
 @retval エラーコード
 */
int TPControl::readTPfile( const string filename )
{
  int ierr = TP_NO_ERROR;
  if( !tp ) return TP_ERROR;

  // read
  if( (ierr = tp->read(filename)) != TP_NO_ERROR )
  {
    cout << "ERROR : in input file: " << filename << endl
         << "  ERROR CODE = "<< ierr << endl;
    return ierr;
  }
  return ierr;
}


/**
 @fn bool TPControl::GetVector(string label, int *vec, const int nvec)
 @brief TextParser入力ファイルからベクトル値を取得する（整数型）
 @param label 取得するベクトルのラベル（絶対パス）
 @param *vec ベクトル格納配列ポインタ
 @param nvec ベクトルサイズ
 */
bool TPControl::GetVector(string label, int *vec, const int nvec)
{
  int i;
  int ierr = TP_NO_ERROR;
  string value;

  if( !tp ) return false;

  // ラベルがあるかチェック
  if( !chkLabel(label)){
	  //std::cout <<  "WARNING no label " << label << std::endl;
	  return false;
  }

  // get value
  if( (ierr = tp->getValue(label, value)) != TP_NO_ERROR ) return false;

  // get type
  TextParserValueType type = tp->getType(label, &ierr);
  if( ierr != TP_NO_ERROR ) return false;
  if( type != TP_VECTOR_NUMERIC ) return false;

  // split
  vector<string> vec_value;
  if( (ierr = tp->splitVector(value, vec_value)) != TP_NO_ERROR ) return false;

  // check number of vector element
  if( vec_value.size() != nvec ) return false;

  // string to real
  for( i=0;i<vec_value.size();i++ )
  {
    vec[i] = tp->convertInt(vec_value[i], &ierr);
    if( ierr != TP_NO_ERROR ) return false;
  }

  return true;
}

/**
 @fn bool TPControl::GetVector(string label, REAL_TYPE *vec, const int nvec)
 @brief TextParser入力ファイルからベクトル値を取得する（実数型）
 @param label 取得するベクトルのラベル（絶対パス）
 @param *vec ベクトル格納配列ポインタ
 @param nvec ベクトルサイズ
 */
bool TPControl::GetVector(string label, REAL_TYPE *vec, const int nvec)
{
  int i;
  int ierr = TP_NO_ERROR;
  string value;

  if( !tp ) return false;

  // ラベルがあるかチェック
  if( !chkLabel(label)){
	  //std::cout <<  "WARNING no label " << label << std::endl;
	  return false;
  }

  // get value
  if( (ierr = tp->getValue(label, value)) != TP_NO_ERROR ){
	  cout << " GetVector debug 333" << endl;
	  return false;
  }

  // get type
  TextParserValueType type = tp->getType(label, &ierr);
  if( ierr != TP_NO_ERROR ) return false;
  if( type != TP_VECTOR_NUMERIC ) return false;

  // split
  vector<string> vec_value;
  if( (ierr = tp->splitVector(value, vec_value)) != TP_NO_ERROR ) return false;

  // check number of vector element
  if( vec_value.size() != nvec ) return false;

  // string to real
  for( i=0;i<vec_value.size();i++ )
  {
    vec[i] = tp->convertDouble(vec_value[i], &ierr);
    if( ierr != TP_NO_ERROR ) return false;
  }

  return true;
}

/**
 @fn bool TPControl::GetVector(string label, string *vec, const int nvec)
 @brief TextParser入力ファイルからベクトル値を取得する（文字列型）
 @param label 取得するベクトルのラベル（絶対パス）
 @param *vec ベクトル格納配列ポインタ
 @param nvec ベクトルサイズ
 */
bool TPControl::GetVector(string label, string *vec, const int nvec)
{
  int i;
  int ierr = TP_NO_ERROR;
  string value;

  if( !tp ) return false;

  // ラベルがあるかチェック
  if( !chkLabel(label)){
	  //std::cout <<  "WARNING no label " << label << std::endl;
	  return false;
  }

  // get value
  if( (ierr = tp->getValue(label, value)) != TP_NO_ERROR ) return false;

  // get type
  TextParserValueType type = tp->getType(label, &ierr);
  if( ierr != TP_NO_ERROR ) return false;
  if( type != TP_VECTOR_NUMERIC ) return false;

  // split
  vector<string> vec_value;
  if( (ierr = tp->splitVector(value, vec_value)) != TP_NO_ERROR ) return false;

  // check number of vector element
  if( vec_value.size() != nvec ) return false;

  // string to string
  for( i=0;i<vec_value.size();i++ )
  {
    vec[i] = vec_value[i];
  }

  return true;
}


/**
 @fn bool TPControl::GetValue(const string label, int *ct)
 @brief TextParser入力ファイルから変数を取得する（整数型）
 @param label 取得する変数のラベル（絶対パス）
 @param *ct 変数格納ポインタ
 */
bool TPControl::GetValue(const string label, int *ct)
{
  int ierror;
  std::string value;
  
  if( !tp ) return false;

  // ラベルがあるかチェック
  if( !chkLabel(label)){
	  //std::cout <<  "WARNING no label " << label << std::endl;
	  return false;
  }

  //std::cout << " label: " << label << std::endl;

  //値の取得
  ierror=tp->getValue(label,value);//labelは絶対パスを想定
  if (ierror != TP_NO_ERROR){
	  std::cout << " label: " << label << std::endl;
	  std::cout <<  "ERROR no label " << label << ierror << std::endl;
	  return false;
  }

  //std::cout << " value: " << value << std::endl;

  //型の取得
  TextParserValueType type = tp->getType(label, &ierror);
  if (ierror != TP_NO_ERROR){
	  std::cout << " label: " << label << std::endl;
	  std::cout <<  "ERROR in TextParser::getType file: " << ierror << std::endl;
	  return false;
  }
  if( type != TP_NUMERIC_VALUE ){
	  std::cout << " label: " << label << std::endl;
	  std::cout <<  "ERROR in TextParser::Type error: " << ierror << std::endl;
	  return false;
  }

  //std::cout << " value type: " << type << std::endl;

  // string to real
  int val = tp->convertInt(value, &ierror);
  if (ierror != TP_NO_ERROR){
	  std::cout << " label: " << label << std::endl;
	  std::cout <<  "ERROR convertInt " << ierror << std::endl;
	  return false;
  }
  //std::cout << " val: " << val << std::endl;

  *ct=val;

  return true;
}

/**
 @fn bool TPControl::GetValue(const string label, REAL_TYPE *ct)
 @brief TextParser入力ファイルから変数を取得する（実数型）
 @param label 取得する変数のラベル（絶対パス）
 @param *ct 変数格納ポインタ
 */
bool TPControl::GetValue(const string label, REAL_TYPE *ct)
{
  int ierror;
  std::string value;
  std::string node;

  if( !tp ) return false;
  
  // ラベルがあるかチェック
  if( !chkLabel(label)){
	  //std::cout <<  "WARNING no label " << label << std::endl;
	  return false;
  }

  //値の取得
  ierror=tp->getValue(label,value);//labelは絶対パスを想定
  if (ierror != TP_NO_ERROR){
	  std::cout << " label: " << label << std::endl;
	  std::cout <<  "ERROR no label " << ierror << std::endl;
	  return false;
  }

  //std::cout << " value: " << value << std::endl;

  //型の取得
  TextParserValueType type = tp->getType(label, &ierror);
  if (ierror != TP_NO_ERROR){
	  std::cout << " label: " << label << std::endl;
	  std::cout <<  "ERROR in TextParser::getType file: " << ierror << std::endl;
	  return false;
  }
  if( type != TP_NUMERIC_VALUE ){
	  std::cout << " label: " << label << std::endl;
	  std::cout <<  "ERROR in TextParser::Type error: " << ierror << std::endl;
	  return false;
  }

  //std::cout << " value type: " << type << std::endl;

  // string to real
  REAL_TYPE val = tp->convertFloat(value, &ierror);
  if (ierror != TP_NO_ERROR){
	  std::cout << " label: " << label << std::endl;
	  std::cout <<  "ERROR convertInt " << ierror << std::endl;
	  return false;
  }
  //std::cout << " val: " << val << std::endl;

  *ct=val;

  return true;
}


/**
 @fn bool TPControl::GetValue(const string label, string *ct)
 @brief TextParser入力ファイルから変数を取得する（文字列型）
 @param label 取得する変数のラベル（絶対パス）
 @param *ct 変数格納ポインタ
 */
bool TPControl::GetValue(const string label, string *ct)
{
  int ierror;
  std::string value;
  
  if( !tp ) return false;

  // ラベルがあるかチェック
  if( !chkLabel(label)){
	  //std::cout <<  "WARNING no label " << label << std::endl;
	  return false;
  }

  //std::cout << " label: " << label << std::endl;

  //値の取得
  ierror=tp->getValue(label,value);//labelは絶対パスを想定
  if (ierror != TP_NO_ERROR){
	  std::cout << " label: " << label << std::endl;
	  std::cout <<  "ERROR no label " << label << std::endl;
	  return false;
  }

  //std::cout << " value: " << value << std::endl;

  //型の取得
  TextParserValueType type = tp->getType(label, &ierror);
  if (ierror != TP_NO_ERROR){
	  std::cout << " label: " << label << std::endl;
	  std::cout <<  "ERROR in TextParser::getType file: " << ierror << std::endl;
	  return false;
  }
  if( type != TP_STRING_VALUE ){
	  std::cout << " label: " << label << std::endl;
	  std::cout <<  "ERROR in TextParser::Type error: " << ierror << std::endl;
	  return false;
  }

  //std::cout << " value type: " << type << std::endl;

  *ct=value;

  return true;
}


/**
 @fn bool TPControl::chkLabel(const string label)
 @brief ラベルの有無をチェック
 @param label チェックするラベル（絶対パス）
 */
bool TPControl::chkLabel(const string label)
{
  int ierror;
  std::string value;
  
  if( !tp ) return false;
  
  //std::cout << " label: " << label << std::endl;

  // ラベルがあるかチェック
  std::vector<std::string> labels;
  ierror=tp->getAllLabels(labels);
  if (ierror != 0){
    std::cout <<  "ERROR in TextParser::getAllLabels file: "  
        << " ERROR CODE "<< ierror << std::endl;
    return false;
  }
  int flag=0;
  for (int i = 0; i < labels.size(); i++) {
	  if( !strcasecmp(label.c_str(), labels[i].c_str()) ){
		  //std::cout << i << " " << labels[i] << std::endl;
		  flag=1;
		  break;
	  }
  }
  if(flag==0){
	  //std::cout << " label: " << label << std::endl;
	  //std::cout <<  "ERROR no label " << label << std::endl;
	  return false;
  }
  return true;
}


/**
 @fn bool TPControl::chkNode(const string label)
 @brief ノードの有無をチェック
 @param label チェックするノード（絶対パス）
 */
bool TPControl::chkNode(const string label)
{
  int ierror;
  std::string node;
  std::vector<std::string> labels;
  int len=label.length();

  if( !tp ) return false;
  
  //std::cout << " label: " << label << std::endl;

  // Nodeがあるかチェック
  ierror=tp->getAllLabels(labels);
  if (ierror != 0){
    std::cout <<  "ERROR in TextParser::getAllLabels file: "  
        << " ERROR CODE "<< ierror << std::endl;
    return false;
  }
  int flag=0;
  for (int i = 0; i < labels.size(); i++) {
	  node=labels[i].substr(0,len);
	  //std::cout << i << " labels = " << labels[i] << std::endl;
	  //std::cout << i << " node   = " << node  << std::endl;
	  if( node==label ){
		  //std::cout << i << " " << labels[i] << std::endl;
		  flag=1;
		  break;
	  }
  }
  if(flag==0){
	  //std::cout << " label: " << label << std::endl;
	  //std::cout <<  "ERROR no label " << label << std::endl;
	  return false;
  }
  return true;
}

/**
 @fn bool TPControl::countLabels(const string label)
 @brief ノード以下のラベルの数を数える
 @param label ラベルを数えるノードの絶対パス
 @retval ラベルの数（エラー、もしくはない場合は-1を返す）
 */
int TPControl::countLabels(const string label)
{
  int ierror;
  std::string node,str,chkstr="";
  std::vector<std::string> labels;
  int len=label.length();
  int flag=0;
  int inode=0;
  int next=0;

  if( !tp ) return -1;
  
  // Nodeがあるかチェック
  ierror=tp->getAllLabels(labels);
  if (ierror != 0){
    std::cout <<  "ERROR in TextParser::getAllLabels file: "  
        << " ERROR CODE "<< ierror << std::endl;
    return -1;
  }
  for (int i = 0; i < labels.size(); i++) {
	  node=labels[i].substr(0,len);
	  //std::cout << i << " labels = " << labels[i] << std::endl;
	  //std::cout << i << " node   = " << node  << std::endl;
	  if( !strcasecmp(node.c_str(), label.c_str()) ){
		  str=labels[i].substr(len+1);
		  //std::cout << " str    = " << str  << std::endl;
		  //std::cout << " chkstr = " << chkstr << std::endl;
		  next=str.find("/");
		  //std::cout << " init inode = " << inode  << std::endl;
		  if(next==0) inode++;
		  else{
			  if(chkstr!=str.substr(0,next)){
				  chkstr=str.substr(0,next);
			  //if(chkstr!=str.substr(0,next-1)){
				 // chkstr=str.substr(0,next-1);
				  inode++;
			  }
		  }
		  //std::cout << " end  inode = " << inode  << std::endl;
	  }
  }

  return inode;
}


/**
 @fn bool TPControl::GetNodeStr(const string label, const int nnode, string *ct)
 @brief ノード以下のnnode番目の文字列を取得する
 @param label ノードの絶対パス
 @param nnode 取得する文字列が現れる順番
 @param *ct 取得した文字列
 */
bool TPControl::GetNodeStr(const string label, const int nnode, string *ct)
{
  int ierror;
  std::string node,str,chkstr="";
  std::vector<std::string> labels;
  int len=label.length();
  int flag=0;
  int inode=0;
  int next=0;

  if( !tp ) return -1;
  
  //std::cout << " label = " << label << std::endl;
  //std::cout << " nnode = " << nnode  << std::endl;

  // Nodeがあるかチェック
  ierror=tp->getAllLabels(labels);
  if (ierror != 0){
    std::cout <<  "ERROR in TextParser::getAllLabels file: "  
        << " ERROR CODE "<< ierror << std::endl;
    return false;
  }
  for (int i = 0; i < labels.size(); i++) {
	  node=labels[i].substr(0,len);
	  //std::cout << i << " labels = " << labels[i] << std::endl;
	  //std::cout << i << " node   = " << node  << std::endl;

	  if( !strcasecmp(node.c_str(), label.c_str()) ){
		  str=labels[i].substr(len+1);
		  next=str.find("/");
		  //std::cout << " str = " << str  << std::endl;
		  //std::cout << " next = " << next << std::endl;
		  if(next==0){
			  inode++;
		  }
		  else{
			  if(chkstr!=str.substr(0,next)){
				  chkstr=str.substr(0,next);
				  inode++;
			  }
		  }
		  //std::cout << " inode = " << inode  << std::endl;
		  //std::cout << " chkstr = " << chkstr << std::endl;
		  if(inode==nnode){
			  *ct=chkstr;
			  return true;
		  }
	  }
  }
  return false;
}

