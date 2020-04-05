#ifndef _FB_DATA_HOLDER_H_
#define _FB_DATA_HOLDER_H_

//##################################################################################
//
// Flow Base class
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2015 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2016 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
// Copyright (c) 2016-2020 Research Institute for Information Research, Kyushu University.
// All rights reserved.
//
//##################################################################################

/*
 * @file   DataHolder.h
 * @brief  FlowBase DataHolder class Header
 * @author aics
 */

#include "DomainInfo.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <iostream>
#include <fstream>
#include "FB_Define.h"
#include "TextParser.h"
#include <stdarg.h> // for va_start()
#include "mydebug.h"

/* ------------- DataHolder -----------------------------------------*/


/// データ格納用コンテナ. 時刻値をキーとしたマップ.
typedef std::map<REAL_TYPE,REAL_TYPE> DATA_MAP;


/**
 * DataHolder データ格納クラス.
 */
class DataHolder {
protected:
  std::string m_label;  ///< ラベル
  std::string m_file;   ///< データファイル
  int m_id;        ///< MPIランク番号

  DATA_MAP m_data;  ///< データコンテナ

  int m_dataSize; ///< データ数

  REAL_TYPE m_minKey;   ///< 時刻最小値
  REAL_TYPE m_maxKey;   ///< 時刻最大値

public:

  /// デフォルトコンストラクタ.
  DataHolder() {}

  /// コンストラクタ.
  ///
  ///   @param[in] label ラベル
  ///   @param[in] file データファイルのパス
  ///   @param[in] id MPIランク番号
  ///
  DataHolder(const char* label, const char* file, int id);

  /// デストラクタ.
  ~DataHolder() {}

  /// 時刻を指定して値を取得する.
  ///
  ///   @param key 時刻値
  ///   @return 指定した時刻に対する値
  ///
  REAL_TYPE getValue(REAL_TYPE key) const;

  /// 時刻を指定して値を取得する.
  ///   getValueメソッドと同じ.
  ///   @param key 時刻値
  ///   @return 指定した時刻に対する値
  ///
  REAL_TYPE get_value(REAL_TYPE key) const {
    return getValue(key);
  }

  /// データ数取得.
  ///
  ///  @return データ数
  ///
  int getDataSize() const {
    return m_dataSize;
  }

  /// 値をスケーリングする．
  ///
  ///   指定された倍率を全てのデータに掛ける．
  ///   @param ratio 倍率
  ///
  void scale(REAL_TYPE ratio);

  /// 時刻の最小値を取得.
  ///
  ///   @return 時刻の最小値
  ///
  REAL_TYPE getMin() const {
    return m_minKey;
  }

  /// 時刻の最大値を取得.
  ///
  ///   @return 時刻の最大値
  ///
  REAL_TYPE getMax() const {
    return m_maxKey;
  }

  /// 格納データ情報を出力.
  ///
  ///   @param[in] fp 出力ファイルポインタ
  ///
  void printInfo(FILE* fp) const;

  /// 全格納データを表示(デバッグ用).
  ///
  ///   @param[in] fp 出力ファイルポインタ
  ///
  void printInfoDebug(FILE* fp) const;

protected:

  /// データファイル読み込み.
  ///
  ///   @param[in] fp 入力データファイルポインタ
  ///
  void readData(FILE* fp);

  /// コメント行をチェック.
  ///
  ///   @param[in] p 行バッファ先頭ポインタ
  ///   @return true:コメント行または空行/false:その他
  ///
  bool isCommentLine(const char* p);

  /// エラーメッセージ出力.
  ///
  ///   @param[in] fmt  出力フォーマット文字列
  ///
  void printError(const char* fmt, ...);

};


/* ------------- DataHolderManager ----------------------------------*/

/// DetataHolderコンテナ. ラベルをキーとしたマップ.
typedef std::map<std::string,DataHolder*> DATA_HOLDER_MAP;


/**
 * DataHolderManager データ格納管理クラス.
 */
class DataHolderManager : public DomainInfo {

protected:

  DATA_HOLDER_MAP m_dataHolders;  ///< DataHolderコンテナ(ラベルをキーとしたマップ)
  
  TextParser* tpCntl;
  
  /// エラーメッセージ出力.
  ///
  ///   @param[in] fmt  出力フォーマット文字列
  ///
  void printError(const char* fmt, ...);

public:

  /// コンストラクタ.
  DataHolderManager() {}
    
  /// デストラクタ.
  ~DataHolderManager() {
    DATA_HOLDER_MAP::iterator it;
    for (it = m_dataHolders.begin(); it != m_dataHolders.end(); it++) {
      delete it->second;
    }
  }

  /// ラベルを指定してDataHolderを取得.
  ///
  ///   @param [in] label ラベル
  ///   @return 対応するDataHolderへのポインタ
  ///
  ///   @note 対応するDataHolderがない場合は0を返す.
  ///
  DataHolder* getDataHolder(const std::string& label) {
    DATA_HOLDER_MAP::iterator it = m_dataHolders.find(label);
    if (it != m_dataHolders.end()) {
      return it->second;
    }
    else {
      return 0;
    }
  }


  /// 管理しているDataHolder数を取得.
  ///
  ///   @return DataHolder数
  ///
  int getNumDataHolder() {
    return m_dataHolders.size();
  }

  
  /// 管理しているDataHolderの情報を出力.
  ///
  ///   @param fp ファイルポインタ
  ///
  void printInfo(FILE* fp) const;

  
  /// 管理しているDataHolderの全データを出力(デバッグ用).
  ///
  ///   @param fp ファイルポインタ
  ///
  void printInfoDebug(FILE* fp) const;


  /// データファイル読み込み.
  ///
  ///   設定XMLファイルをパースし, 入力データファイルを読み込み,
  ///   DataHolderを生成し登録する.
  ///
  void readData();
  
  /**
   * @brief TPのポインタを受け取る
   * @param [in] tp  TextParserのポインタ
   * @return エラーコード
   */
  bool importTP(TextParser* tp);

};


#endif // _FB_DATA_HOLDER_H_
