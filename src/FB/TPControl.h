#ifndef _FB_TPCONTROL_H_
#define _FB_TPCONTROL_H_

// #################################################################
//
// CAERU Library
//
// Copyright (c) All right reserved. 2012
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
//
// #################################################################

/**
 @file TPControl.h
 @brief TextParser Control class Header
 @author kero
 */

#include <math.h>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "string.h"

#include "cpm_Define.h"
#include "cpm_ParaManager.h"

#include "FB_Define.h"
#include "TextParser.h"

using namespace std;

class SubDomain {
public:
  /** コンストラクタ */
  SubDomain()
  {
    m_pos[0] = m_pos[1] = m_pos[2] = 0;
    m_bcid[0] = m_bcid[1] = m_bcid[2] = m_bcid[3] = m_bcid[4] = m_bcid[5] = 0;
  }
  
  /**　デストラクタ */
  ~SubDomain()
  {
  }

  int m_pos[3];
  int m_bcid[6];
};



class DomainInfo {
public:
  /** コンストラクタ */
  DomainInfo()
  {
    m_globalOrigin[0] = m_globalOrigin[1] = m_globalOrigin[2] = REAL_TYPE(0.0);
    m_globalRegion[0] = m_globalRegion[1] = m_globalRegion[2] = REAL_TYPE(0.0);
    m_globalPitch[0]  = m_globalPitch[1]  = m_globalPitch[2]  = REAL_TYPE(0.0);
    m_domainDiv[0]    = m_domainDiv[1]    = m_domainDiv[2]    = 0;
  }
  
  /**　デストラクタ */
  ~DomainInfo()
  {
  }

  REAL_TYPE m_VoxelOrigin[3];
  REAL_TYPE m_VoxelSize[3];
  REAL_TYPE m_VoxelPitch[3];

  REAL_TYPE m_globalOrigin[3];
  REAL_TYPE m_globalRegion[3];
  REAL_TYPE m_globalPitch[3];
  int       m_domainDiv[3];
  vector<SubDomain> m_subDomain;
};



class MediumTableInfo {
public:
  int type;
  string label;
  map<string, REAL_TYPE> m_fval;

public:
  /** コンストラクタ */
  MediumTableInfo() 
  {
    type = -1;
  }
  
  /**　デストラクタ */
  ~MediumTableInfo() 
  {
  }
  
};

class TPControl {

private:
	TextParser* tp;  ///< テキストパーサ

public:
  /** コンストラクタ */
	TPControl(){};

  /**　デストラクタ */
	~TPControl(){};
	

  /**
   @brief TextParser入力ファイルからベクトル値を取得する（整数型）
   @param[int] label 取得するベクトルのラベル（絶対パス）
   @param[out] vec   ベクトル格納配列ポインタ
   @param[in]  nvec  ベクトルサイズ
   */
	bool GetVector(const string label, int *vec, const int nvec);
  
  
  /**
   @brief TextParser入力ファイルからベクトル値を取得する（実数型）
   @param[in]  label  取得するベクトルのラベル（絶対パス）
   @param[out] vec    ベクトル格納配列ポインタ
   @param[in]  nvec   ベクトルサイズ
   */
	bool GetVector(const string label, REAL_TYPE *vec, const int nvec);

  
  /**
   @brief TextParser入力ファイルからベクトル値を取得する（文字列型）
   @param[in]  label  取得するベクトルのラベル（絶対パス）
   @param[out] vec    ベクトル格納配列ポインタ
   @param[in]  nvec   ベクトルサイズ
   */
	bool GetVector(const string label, string *vec, const int nvec);
  
  
  /**
   @brief TextParser入力ファイルから変数を取得する（整数型）
   @param[in]  label 取得する変数のラベル（絶対パス）
   @param[out] ct    変数格納ポインタ
   */
	bool GetValue(const string label, int *ct);
  
  
  /**
   @brief TextParser入力ファイルから変数を取得する（実数型）
   @param[in]  label 取得する変数のラベル（絶対パス）
   @param[out] ct    変数格納ポインタ
   */
	bool GetValue(const string label, REAL_TYPE *ct);
  
  
  /**
   @brief TextParser入力ファイルから変数を取得する（文字列型）
   @param[in]  label 取得する変数のラベル（絶対パス）
   @param[out] ct    変数格納ポインタ
   */
	bool GetValue(const string label, string *ct);

	
  /**
   @brief ラベルの有無をチェック
   @param[in] label チェックするラベル（絶対パス）
   */
	bool chkLabel(const string label);
  
  
  /**
   @brief ノードの有無をチェック
   @param[in] label チェックするノード（絶対パス）
   */
	bool chkNode(const string label);
  
  
  /**
   @brief ノード以下のnnode番目の文字列を取得する
   @param[in]  label ノードの絶対パス
   @param[in]  nnode 取得する文字列が現れる順番
   @param[out] ct    取得した文字列
   */
	bool GetNodeStr(const string label, const int nnode, string *ct);
  
  
  /**
   @brief ノード以下のラベルの数を数える
   @param[in] label ラベルを数えるノードの絶対パス
   @retval ラベルの数（エラー、もしくはない場合は-1を返す）
   */
  int countLabels(const string label);
  
  /**
   @brief TextParserLibraryのインスタンス生成
   @retrun エラーコード
   */
  int getTPinstance();
  
  
  /**
   @brief TextParserオブジェクトに入力ファイルをセットする
   @param[in] filename 入力ファイル名
   @retval エラーコード
   */
	int readTPfile(const string filename);

};

#endif // _FB_TPCONTROL_H_
