/*
 * CPM - Cartesian Partition Manager
 *
 * Copyright (C) 2012 University of Tokyo.
 *
 */

/**
 * @file   cpm_TextParser.h
 * テキストパーサークラスのヘッダーファイル
 * @author University of Tokyo
 * @date   2012/05/31
 */

#ifndef _CPM_TEXTPARSER_H_
#define _CPM_TEXTPARSER_H_

#include "cpm_Base.h"
#include "TextParser.h"

/** CPMのテキストパーサークラス */
class cpm_TextParser : public cpm_Base
{
////////////////////////////////////////////////////////////////////////////////
// メンバー関数
////////////////////////////////////////////////////////////////////////////////
public:


private:


protected:
  /** コンストラクタ */
  cpm_TextParser();

  /** デストラクタ */
  virtual ~cpm_TextParser();

  /** 読み込み処理
   *  - ユーザは直接コールできない
   *  @param[in] filename 読み込むファイル名
   *  @return    TextParserクラスの終了コード
   */
  int Read( std::string filename ); 

  /** ベクトルデータの読み込み(単精度実数版)
   *  @param[in]  label ベクトルデータのテキストラベル
   *  @param[out] vec   読み込んだベクトルデータ(サイズはnvec確保されている必要がある)
   *  @param[in]  nvec  読み込んだベクトルデータの数
   *  @retval     1000未満 テキストパーサのエラーコード
   *  @retval     CPM_ERROR_TP_NOVECTOR(2001) 指定ラベルがベクトルデータではない
   *  @retval     CPM_ERROR_TP_VECTOR_SIZE(2002) ベクトルデータのサイズがnvecと一致しない
   */
  int readVector(std::string label, float *vec, const int nvec);

  /** ベクトルデータの読み込み(倍精度実数版)
   *  @param[in]  label ベクトルデータのテキストラベル
   *  @param[out] vec   読み込んだベクトルデータ(サイズはnvec確保されている必要がある)
   *  @param[in]  nvec  読み込んだベクトルデータの数
   *  @retval     1000未満 テキストパーサのエラーコード
   *  @retval     CPM_ERROR_TP_NOVECTOR(2001) 指定ラベルがベクトルデータではない
   *  @retval     CPM_ERROR_TP_VECTOR_SIZE(2002) ベクトルデータのサイズがnvecと一致しない
   */
  int readVector(std::string label, double *vec, const int nvec);

  /** ベクトルデータの読み込み(整数版)
   *  @param[in]  label ベクトルデータのテキストラベル
   *  @param[out] vec   読み込んだベクトルデータ(サイズはnvec確保されている必要がある)
   *  @param[in]  nvec  読み込んだベクトルデータの数
   *  @retval     1000未満 テキストパーサのエラーコード
   *  @retval     CPM_ERROR_TP_NOVECTOR(2001) 指定ラベルがベクトルデータではない
   *  @retval     CPM_ERROR_TP_VECTOR_SIZE(2002) ベクトルデータのサイズがnvecと一致しない
   */
  int readVector(std::string label, int *vec, const int nvec);

  


////////////////////////////////////////////////////////////////////////////////
// メンバー変数
////////////////////////////////////////////////////////////////////////////////
public:


private:


protected:
  /** テキストパーサークラスのインスタンス */
  TextParser *m_tp;


};

#endif /* _CPM_TEXTPARSER_H_ */

