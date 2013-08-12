//##################################################################################
//
// Flow Base class
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2013 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2013 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

/**
 * @file   TPControl.C
 * @brief  TextParser Control class
 * @author kero
 */

#include "TPControl.h"



// #################################################################
// ラベルの有無をチェック
bool TPControl::chkLabel(const string label)
{
  int ierror;
  string value;
  
  if( !tp ) return false;
  
  // ラベルがあるかチェック
  vector<string> labels;
  
  ierror=tp->getAllLabels(labels);
  
  if (ierror != 0)
  {
    cout <<  "ERROR in TextParser::getAllLabels file: "
    << " ERROR CODE "<< ierror << endl;
    Exit(0);
  }
  
  int flag=0;
  for (int i = 0; i < labels.size(); i++)
  {
	  if( !strcasecmp(label.c_str(), labels[i].c_str()) )
    {
		  flag=1;
		  break;
	  }
  }
  
  if (flag==0)
  {
	  return false;
  }
  
  return true;
}



// #################################################################
// ノードの有無をチェック
bool TPControl::chkNode(const string label)
{
  int ierror;
  string node;
  vector<string> labels;
  int len=label.length();
  
  if( !tp ) return false;
  
  // Nodeがあるかチェック
  ierror = tp->getAllLabels(labels);
  
  if (ierror != 0)
  {
    cout <<  "ERROR in TextParser::getAllLabels file: "
    << " ERROR CODE "<< ierror << endl;
    Exit(0);
  }
  
  int flag=0;
  for (int i = 0; i < labels.size(); i++) {
	  node = labels[i].substr(0,len);
    
    if ( !strcasecmp(node.c_str(), label.c_str()) )
    {
		  flag=1;
		  break;
	  }
  }
  
  if (flag==0)
  {
	  return false;
  }
  
  return true;
}



// #################################################################
// ノード以下のラベルの数を数える
int TPControl::countLabels(const string label)
{
  int ierror;
  string node,str,chkstr="";
  vector<string> labels;
  int len=label.length();
  int flag=0;
  int inode=0;
  int next=0;
  
  if( !tp ) return -1;
  
  // Nodeがあるかチェック
  ierror=tp->getAllLabels(labels);
  
  if (ierror != 0){
    cout <<  "ERROR in TextParser::getAllLabels file: "
    << " ERROR CODE "<< ierror << endl;
    return -1;
  }
  
  for (int i = 0; i < labels.size(); i++) {
	  node=labels[i].substr(0,len);
    
	  if( !strcasecmp(node.c_str(), label.c_str()) ){
		  str=labels[i].substr(len+1);
		  next=str.find("/");
      
		  if(next==0) inode++;
		  else{
			  if(chkstr!=str.substr(0,next)){
				  chkstr=str.substr(0,next);
				  inode++;
			  }
		  }
      
	  }
  }
  
  return inode;
}



// #################################################################
// TextParserLibraryのインスタンス生成
void TPControl::getTPinstance()
{
  tp = new TextParser;
}



// #################################################################
// ノード以下のnnode番目の文字列を取得する
bool TPControl::GetNodeStr(const string label, const int nnode, string *ct)
{
  if ( !tp ) return -1;
  
  int ierror;
  int len=label.length();
  int flag=0;
  int inode=0;
  int next=0;
  
  string node;
  string str;
  string chkstr="";
  vector<string> labels;

  
  // Nodeがあるかチェック
  ierror = tp->getAllLabels(labels);
  
  if (ierror != 0)
  {
    cout <<  "ERROR in TextParser::getAllLabels file: " << " ERROR CODE "<< ierror << endl;
    Exit(0);
  }
  
  for (int i = 0; i < labels.size(); i++) {
	  node = labels[i].substr(0, len);
    
	  if ( !strcasecmp(node.c_str(), label.c_str()) )
    {
		  str = labels[i].substr(len+1);
		  next = str.find("/");
      
		  if ( next == 0 )
      {
			  inode++;
		  }
		  else
      {
			  if ( chkstr != str.substr(0, next) )
        {
				  chkstr = str.substr(0, next);
				  inode++;
			  }
		  }
      
		  if ( inode == nnode )
      {
			  *ct = chkstr;
			  return true;
		  }
	  }
  }
  return false;
}



// #################################################################
// TextParser入力ファイルからベクトル値を取得する（整数型）
bool TPControl::GetVector(const string label, int *vec, const int nvec)
{
  int ierr = TP_NO_ERROR;
  string value;

  if( !tp ) return false;

  // ラベルがあるかチェック
  if ( !chkLabel(label)) return false;

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
  for(int i=0; i<vec_value.size(); i++ )
  {
    vec[i] = tp->convertInt(vec_value[i], &ierr);
    if( ierr != TP_NO_ERROR ) return false;
  }

  return true;
}


// #################################################################
// TextParser入力ファイルからベクトル値を取得する（実数型）
bool TPControl::GetVector(const string label, REAL_TYPE *vec, const int nvec)
{
  int ierr = TP_NO_ERROR;
  string value;

  if( !tp ) return false;

  // ラベルがあるかチェック
  if ( !chkLabel(label)) return false;

  // get value
  if ( (ierr = tp->getValue(label, value)) != TP_NO_ERROR )
  {
	  Exit(0);
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
  for(int i=0; i<vec_value.size(); i++ )
  {
    vec[i] = tp->convertDouble(vec_value[i], &ierr);
    if( ierr != TP_NO_ERROR ) return false;
  }

  return true;
}


// #################################################################
// TextParser入力ファイルからベクトル値を取得する（文字列型）
bool TPControl::GetVector(string label, string *vec, const int nvec)
{
  int ierr = TP_NO_ERROR;
  string value;

  if( !tp ) return false;

  // ラベルがあるかチェック
  if ( !chkLabel(label)) return false;

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
  for(int i=0;i<vec_value.size();i++ )
  {
    vec[i] = vec_value[i];
  }

  return true;
}


// #################################################################
// TextParser入力ファイルから変数を取得する（整数型）
bool TPControl::GetValue(const string label, int *ct)
{
  int ierror;
  string value;
  
  if( !tp ) return false;

  // ラベルがあるかチェック
  if ( !chkLabel(label)) return false;

  //値の取得
  ierror=tp->getValue(label,value);//labelは絶対パスを想定
  if (ierror != TP_NO_ERROR){
	  cout << " label: " << label << endl;
	  cout <<  "ERROR no label " << label << ierror << endl;
	  Exit(0);
  }

  //型の取得
  TextParserValueType type = tp->getType(label, &ierror);
  if (ierror != TP_NO_ERROR){
	  cout << " label: " << label << endl;
	  cout <<  "ERROR in TextParser::getType file: " << ierror << endl;
	  Exit(0);
  }
  if( type != TP_NUMERIC_VALUE ){
	  cout << " label: " << label << endl;
	  cout <<  "ERROR in TextParser::Type error: " << ierror << endl;
	  Exit(0);
  }

  // string to real
  int val = tp->convertInt(value, &ierror);
  if (ierror != TP_NO_ERROR){
	  cout << " label: " << label << endl;
	  cout <<  "ERROR convertInt " << ierror << endl;
	  Exit(0);
  }

  *ct=val;

  return true;
}


// #################################################################
// TextParser入力ファイルから変数を取得する（実数型）
bool TPControl::GetValue(const string label, float *ct)
{
  int ierror;
  string value;
  string node;

  if( !tp ) return false;
  
  // ラベルがあるかチェック
  if ( !chkLabel(label)) return false;

  //値の取得
  ierror=tp->getValue(label,value);//labelは絶対パスを想定
  if (ierror != TP_NO_ERROR){
	  cout << " label: " << label << endl;
	  cout <<  "ERROR no label " << ierror << endl;
	  Exit(0);
  }

  //型の取得
  TextParserValueType type = tp->getType(label, &ierror);
  if (ierror != TP_NO_ERROR){
	  cout << " label: " << label << endl;
	  cout <<  "ERROR in TextParser::getType file: " << ierror << endl;
	  Exit(0);
  }
  if( type != TP_NUMERIC_VALUE ){
	  cout << " label: " << label << endl;
	  cout <<  "ERROR in TextParser::Type error: " << ierror << endl;
	  Exit(0);
  }

  // string to real
  float val = tp->convertFloat(value, &ierror);
  if (ierror != TP_NO_ERROR){
	  cout << " label: " << label << endl;
	  cout <<  "ERROR convertInt " << ierror << endl;
	  Exit(0);
  }

  *ct=val;

  return true;
}


// #################################################################
// TextParser入力ファイルから変数を取得する（実数型）
bool TPControl::GetValue(const string label, double *ct)
{
  int ierror;
  string value;
  string node;
  
  if( !tp ) return false;
  
  // ラベルがあるかチェック
  if ( !chkLabel(label)) return false;
  
  //値の取得 labelは絶対パスを想定
  ierror=tp->getValue(label,value);
  
  if (ierror != TP_NO_ERROR){
	  cout << " label: " << label << endl;
	  cout <<  "ERROR no label " << ierror << endl;
	  Exit(0);
  }
  
  //型の取得
  TextParserValueType type = tp->getType(label, &ierror);
  if (ierror != TP_NO_ERROR){
	  cout << " label: " << label << endl;
	  cout <<  "ERROR in TextParser::getType file: " << ierror << endl;
	  Exit(0);
  }
  if( type != TP_NUMERIC_VALUE ){
	  cout << " label: " << label << endl;
	  cout <<  "ERROR in TextParser::Type error: " << ierror << endl;
	  Exit(0);
  }
  
  // string to real
  double val = tp->convertDouble(value, &ierror);
  if (ierror != TP_NO_ERROR){
	  cout << " label: " << label << endl;
	  cout <<  "ERROR convertInt " << ierror << endl;
	  Exit(0);
  }
  
  *ct=val;
  
  return true;
}


// #################################################################
// TextParser入力ファイルから変数を取得する（文字列型）
bool TPControl::GetValue(const string label, string *ct)
{
  int ierror;
  string value;
  
  if ( !tp ) return false;
  
  // ラベルがあるかチェック
  if ( !chkLabel(label)) return false;
  
  //値の取得
  ierror = tp->getValue(label, value); //labelは絶対パスを想定
  
  if (ierror != TP_NO_ERROR)
  {
	  cout << " label: " << label << endl;
	  cout <<  "ERROR no label " << label << endl;
	  Exit(0);
  }
  
  //型の取得
  TextParserValueType type = tp->getType(label, &ierror);
  if (ierror != TP_NO_ERROR)
  {
	  cout << " label: " << label << endl;
	  cout <<  "ERROR in TextParser::getType file: " << ierror << endl;
	  Exit(0);
  }
  
  if( type != TP_STRING_VALUE )
  {
	  cout << " label: " << label << endl;
	  cout <<  "ERROR in TextParser::Type error: " << ierror << endl;
	  Exit(0);
  }
  
  *ct=value;
  
  return true;
}


// #################################################################
// ラベルリストを作成する
// @note ラベルはパラメータファイルの出現順(2)
bool TPControl::getLabelVector(const string root, vector<string>& nodes)
{
  if ( !tp ) return false;
  
  if ( tp->changeNode(root) ) Exit(0); // 0 - no error
  
  tp->getNodes(nodes, 2);
  
  return true;
}

// #################################################################
// TextParserオブジェクトに入力ファイルをセットする
bool TPControl::readTPfile(const string filename)
{
  if( !tp )
  {
    cout << "No instanec of TextParser" << endl;
    return false;
  }
  
  TextParserError ierr;

  // read
  if ( (ierr=tp->read(filename)) != TP_NO_ERROR )
  {
    cout << "ERROR : in input file: " << filename << endl
    << "  ERROR CODE = "<< ierr << endl;
    return false;
  }
  return true;
}
