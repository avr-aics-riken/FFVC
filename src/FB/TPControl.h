#ifndef _FB_TPCONTROL_H_
#define _FB_TPCONTROL_H_

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
 * @file   TPControl.h
 * @brief  TextParser Control class Header
 * @author kero
 */

#include <math.h>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "string.h"

#include "cpm_ParaManager.h"

#include "FB_Define.h"
#include "TextParser.h"

using namespace std;


class TPControl {

private:
	TextParser* tp;  ///< テキストパーサ

public:
  /** コンストラクタ */
	TPControl(){};

  /**　デストラクタ */
	~TPControl(){};
	

  /**
   * @brief ラベルリストを作成し，重複をチェックする
   * @param [in] root  テストするパラメータディレクトリ
   * @param [in] nodes ラベルvector
   */
  bool getLabelVector(const string root, vector<string>& nodes);
  
  
  /**
   * @brief TextParser入力ファイルからベクトル値を取得する（整数型）
   * @param [int] label 取得するベクトルのラベル（絶対パス）
   * @param [out] vec   ベクトル格納配列ポインタ
   * @param [in]  nvec  ベクトルサイズ
   */
	bool GetVector(const string label, int *vec, const int nvec);
  
  
  /**
   * @brief TextParser入力ファイルからベクトル値を取得する（実数型）
   * @param [in]  label  取得するベクトルのラベル（絶対パス）
   * @param [out] vec    ベクトル格納配列ポインタ
   * @param [in]  nvec   ベクトルサイズ
   */
	bool GetVector(const string label, REAL_TYPE *vec, const int nvec);

  
  /**
   * @brief TextParser入力ファイルからベクトル値を取得する（文字列型）
   * @param [in]  label  取得するベクトルのラベル（絶対パス）
   * @param [out] vec    ベクトル格納配列ポインタ
   * @param [in]  nvec   ベクトルサイズ
   */
	bool GetVector(const string label, string *vec, const int nvec);
  
  
  /**
   * @brief TextParser入力ファイルから変数を取得する（整数型）
   * @param [in]  label 取得する変数のラベル（絶対パス）
   * @param [out] ct    変数格納ポインタ
   */
	bool GetValue(const string label, int *ct);
  
  
  /**
   * @brief TextParser入力ファイルから変数を取得する（実数型）
   * @param [in]  label 取得する変数のラベル（絶対パス）
   * @param [out] ct    変数格納ポインタ
   */
	bool GetValue(const string label, float *ct);
  
  
  /**
   * @brief TextParser入力ファイルから変数を取得する（実数型）
   * @param [in]  label 取得する変数のラベル（絶対パス）
   * @param [out] ct    変数格納ポインタ
   */
	bool GetValue(const string label, double *ct);
  
  
  /**
   * @brief TextParser入力ファイルから変数を取得する（文字列型）
   * @param [in]  label 取得する変数のラベル（絶対パス）
   * @param [out] ct    変数格納ポインタ
   */
	bool GetValue(const string label, string *ct);

	
  /**
   * @brief ラベルの有無をチェック
   * @param [in] label チェックするラベル（絶対パス）
   */
	bool chkLabel(const string label);
  
  
  /**
   * @brief ノードの有無をチェック
   * @param [in] label チェックするノード（絶対パス）
   */
	bool chkNode(const string label);
  
  
  /**
   * @brief ノード以下のnnode番目の文字列を取得する
   * @param [in]  label ノードの絶対パス
   * @param [in]  nnode 取得する文字列が現れる順番
   * @param [out] ct    取得した文字列
   */
	bool GetNodeStr(const string label, const int nnode, string *ct);
  
  
  /**
   * @brief ノード以下のラベルの数を数える
   * @param [in] label ラベルを数えるノードの絶対パス
   * @retval ラベルの数（エラー、もしくはない場合は-1を返す）
   */
  int countLabels(const string label);
  
  
  /**
   * @brief TextParserLibraryのインスタンス生成
   * @retrun エラーコード
   */
  void getTPinstance();
  
  
  /**
   * @brief TextParserオブジェクトに入力ファイルをセットする
   * @param [in] filename 入力ファイル名
   * @retval エラーコード
   */
	bool readTPfile(const string filename);
  
  
  /** テキストパーサーの内容を破棄 */
  int remove()
  {
    return tp->remove();
  }
  
  string getVersionInfo()
  {
    return tp->getVersionInfo();
  }

};

#endif // _FB_TPCONTROL_H_
