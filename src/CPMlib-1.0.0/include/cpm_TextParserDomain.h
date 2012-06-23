/*
 * CPM - Cartesian Partition Manager
 *
 * Copyright (C) 2012 University of Tokyo.
 *
 */

/**
 * @file   cpm_TextParserDomain.h
 * 領域情報のテキストパーサークラスのヘッダーファイル
 * @author University of Tokyo
 * @date   2012/05/31
 */

#ifndef _CPM_TEXTPARSER_DOMAIN_H_
#define _CPM_TEXTPARSER_DOMAIN_H_

#include "cpm_TextParser.h"
#include "cpm_DomainInfo.h"

/** CPMの領域情報テキストパーサークラス */
class cpm_TextParserDomain : public cpm_TextParser
{
////////////////////////////////////////////////////////////////////////////////
// メンバー関数
////////////////////////////////////////////////////////////////////////////////
public:
  /** コンストラクタ */
  cpm_TextParserDomain();

  /** デストラクタ */
  virtual ~cpm_TextParserDomain();

  /** 読み込み処理
   *  - TextParserクラスを用いて領域分割情報ファイルを読み込む
   *  - TextParserクラスのインスタンスはクリア(remove)される
   *  @param[in]  filename 読み込むファイル名
   *  @param[out] errorcode CPMエラーコード
   *  @return 領域情報ポインタ
   */
  static cpm_GlobalDomainInfo* Read( std::string filename, int &errorcode ); 


private:
  /** 読み込み処理のメイン
   *  - TextParserクラスを用いて領域分割情報ファイルを読み込む
   *  - TextParserクラスのインスタンスはクリア(remove)される
   *  @param[in]  filename 読み込むファイル名
   *  @param[out] errorcode CPMエラーコード
   *  @return 領域情報ポインタ
   */
  cpm_GlobalDomainInfo* ReadMain( std::string filename, int &errorcode ); 

  /** DomainInfoの読み込み
   *  @param[inout] dInfo 領域情報
   *  @return CPMエラーコード
   */
  int ReadDomainInfo( cpm_GlobalDomainInfo* dInfo );

  /** ActiveSubDomainsの読み込み
   *  @param[inout] dInfo 領域情報
   *  @return CPMエラーコード
   */
  int ReadSubdomainInfo( cpm_GlobalDomainInfo* dInfo );


////////////////////////////////////////////////////////////////////////////////
// メンバー変数
////////////////////////////////////////////////////////////////////////////////
public:


private:


};

#endif /* _CPM_TEXTPARSER_DOMAIN_H_ */

