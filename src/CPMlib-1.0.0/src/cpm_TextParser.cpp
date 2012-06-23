/*
 * CPM - Cartesian Partition Manager
 *
 * Copyright (C) 2012 University of Tokyo.
 *
 */

/**
 * @file   cpm_TextParser.cpp
 * TextParserクラスのソースファイル
 * @author University of Tokyo
 * @date   2012/05/31
 */
#include "cpm_TextParser.h"

////////////////////////////////////////////////////////////////////////////////
// コンストラクタ
cpm_TextParser::cpm_TextParser()
  : cpm_Base()
{
  //インスタンス
  m_tp = TextParser::get_instance();
}

////////////////////////////////////////////////////////////////////////////////
// デストラクタ
cpm_TextParser::~cpm_TextParser()
{
  if( m_tp ) m_tp->remove();
}
////////////////////////////////////////////////////////////////////////////////
// 読み込み
int
cpm_TextParser::Read( std::string filename )
{
  int ret = 0;

  // ポインタのチェック
  if( !m_tp ) return TP_ERROR;

  // remove
  if( (ret = m_tp->remove()) != TP_NO_ERROR ) return ret;

  // read
  if( (ret = m_tp->read(filename)) != TP_NO_ERROR ) return ret;

  return ret;
}

////////////////////////////////////////////////////////////////////////////////
// ベクトル値の読み込み(単精度実数版)
int
cpm_TextParser::readVector(std::string label, float *vec, const int nvec)
{
  int ret = TP_NO_ERROR;
  std::string value;

  if( !m_tp ) return TP_ERROR;

  // get value
  if( (ret = m_tp->getValue(label, value)) != TP_NO_ERROR ) return ret;

  // get type
  TextParserValueType type = m_tp->getType(label, &ret);
  if( ret  != TP_NO_ERROR ) return ret;
  if( type != TP_VECTOR_NUMERIC ) return -1;

  // split
  std::vector<std::string> vec_value;
  if( (ret = m_tp->splitVector(value, vec_value)) != TP_NO_ERROR ) return ret;

  // check number of vector element
  if( vec_value.size() != nvec ) return -2;

  // string to real
  for( size_t i=0;i<vec_value.size();i++ )
  {
    vec[i] = m_tp->convertFloat(vec_value[i], &ret);
    if( ret != TP_NO_ERROR ) return ret;
  }

  return ret;
}

////////////////////////////////////////////////////////////////////////////////
// ベクトル値の読み込み(倍精度実数版)
int
cpm_TextParser::readVector(std::string label, double *vec, const int nvec)
{
  int ret = TP_NO_ERROR;
  std::string value;

  if( !m_tp ) return TP_ERROR;

  // get value
  if( (ret = m_tp->getValue(label, value)) != TP_NO_ERROR ) return ret;

  // get type
  TextParserValueType type = m_tp->getType(label, &ret);
  if( ret  != TP_NO_ERROR ) return ret;
  if( type != TP_VECTOR_NUMERIC ) return CPM_ERROR_TP_NOVECTOR;

  // split
  std::vector<std::string> vec_value;
  if( (ret = m_tp->splitVector(value, vec_value)) != TP_NO_ERROR ) return ret;

  // check number of vector element
  if( vec_value.size() != nvec ) return CPM_ERROR_TP_VECTOR_SIZE;

  // string to real
  for( size_t i=0;i<vec_value.size();i++ )
  {
    vec[i] = m_tp->convertDouble(vec_value[i], &ret);
    if( ret != TP_NO_ERROR ) return ret;
  }

  return ret;
}

////////////////////////////////////////////////////////////////////////////////
// ベクトル値の読み込み(整数版)
int
cpm_TextParser::readVector(std::string label, int *vec, const int nvec)
{
  int ret = TP_NO_ERROR;
  std::string value;

  if( !m_tp ) return TP_ERROR;

  // get value
  if( (ret = m_tp->getValue(label, value)) != TP_NO_ERROR ) return ret;

  // get type
  TextParserValueType type = m_tp->getType(label, &ret);
  if( ret  != TP_NO_ERROR ) return ret;
  if( type != TP_VECTOR_NUMERIC ) return -1;

  // split
  std::vector<std::string> vec_value;
  if( (ret = m_tp->splitVector(value, vec_value)) != TP_NO_ERROR ) return ret;

  // check number of vector element
  if( vec_value.size() != nvec ) return -2;

  // string to real
  for( size_t i=0;i<vec_value.size();i++ )
  {
    vec[i] = m_tp->convertInt(vec_value[i], &ret);
    if( ret != TP_NO_ERROR ) return ret;
  }

  return ret;
}

