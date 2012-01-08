/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2010
 *
 */

//@file DataHolder.C
//@brief DataHolder class
//@author keno, FSI Team, VCAD, RIKEN

#include "DataHolder.h"
#include <cassert>

#define BUF_SIZE 256

/* ------------- DataHolder -----------------------------------------*/


/// コンストラクタ.
///
///   @param[in] label ラベル
///   @param[in] file データファイルのパス
///   @param[in] id MPIランク番号
///
DataHolder::DataHolder(const char* label, const char* file, int id) :
  m_label(label), m_file(file), m_id(id)
{
  FILE* fp;
  if (!(fp = fopen(file, "r"))) {
    printError("cannot open data file.\n");
    assert(0);
  }
  readData(fp);

  fclose(fp);
}

/// データファイル読み込み.
///
///   @param[in] fp 入力データファイルポインタ
///
void DataHolder::readData(FILE* fp)
{
  char buf[BUF_SIZE];

  while (fgets(buf, BUF_SIZE, fp)) {
    if (isCommentLine(buf)) continue;

    // parse line
    double k, v;
    if (sscanf(buf, "%lf %lf", &k, &v) != 2 &&     // 空白区切り
        sscanf(buf, "%lf , %lf", &k, &v) != 2) {   // カンマ区切り
      printError("cnnnot parse line: %s\n", buf);
      assert(0);
    }
    SKL_REAL key = (SKL_REAL)k;
    SKL_REAL val = (SKL_REAL)v;

    // 重複を確認
    if (m_data.find(key) != m_data.end()) {
      printError("duplicated key value in line: %s\n", buf);
      assert(0);
    }

    m_data.insert(DATA_MAP::value_type(key, val));
  }

  // データサイズ設定
  m_dataSize = m_data.size();

  if (m_dataSize == 0) {
    printError("no data in file.\n");
    assert(0);
  }

  DATA_MAP::iterator it;

  // 最小キーの設定
  it = m_data.begin();
  m_minKey = (*it).first;

  // 最大キーの設定
  if (m_dataSize == 1) {
    m_maxKey = m_minKey;
  } else {
    it = m_data.end();
    it--;
    m_maxKey = (*it).first;
  }
}


/// コメント行をチェック
///
///   @param[in] p 行バッファ先頭ポインタ
///   @return true:コメント行または空行/false:その他
///
bool DataHolder::isCommentLine(const char* p)
{
  while (*p) {
    if (!isspace(*p)) {
      if (*p == '#') {
        return true;
      }
      else {
        return false;
      }
    }
    p++;
  }
  return true;   // 空行
}


/// 時刻を指定して値を取得する.
///
///   @param key 時刻値
///   @return 指定した時刻に対する値
///
SKL_REAL DataHolder::getValue(SKL_REAL key) const
{
  // キーが最小値以下の場合
  if (key <= m_minKey) {
    return m_data.find(m_minKey)->second;
  }

  // キーが最大値以上の場合
  if (key >= m_maxKey) {
    return m_data.find(m_maxKey)->second;
  }

  // 上記以外は線形補間
  DATA_MAP::const_iterator it = m_data.lower_bound(key);
  SKL_REAL key1 = it->first;
  SKL_REAL val1 = it->second;
  if (key == key1) return val1;
  it--;
  SKL_REAL key0 = it->first;
  SKL_REAL val0 = it->second;
  return (val1 - val0) * (key - key0) / (key1 - key0) + val0;
}


/// 値をスケーリングする．
///
///   指定された倍率を全てのデータに掛ける．
///   @param ratio 倍率
///
void DataHolder::scale(SKL_REAL ratio)
{
  DATA_MAP::iterator it;
  for (it = m_data.begin(); it != m_data.end(); it++) {
      it->second *= ratio;
  }
}


/// 格納データ情報を出力.
///
///   @param[in] fp 出力ファイルポインタ
///
void DataHolder::printInfo(FILE* fp) const
{
  fprintf(fp, "\nDataHolder: %s\n", m_label.c_str());
  fprintf(fp, "\tdata file:  %s\n", m_file.c_str());
  fprintf(fp, "\t# of data:  %d\n", m_dataSize);
  fprintf(fp, "\ttime range: [%g, %g]\n",  m_minKey, m_maxKey);
}


/// 全格納データを表示(デバッグ用).
///
///   @param[in] fp 出力ファイルポインタ
///
void DataHolder::printInfoDebug(FILE* fp) const
{
  printInfo(fp);
  DATA_MAP::const_iterator it;
  fprintf(fp, "\t(time, value):\n");
  for (it = m_data.begin(); it != m_data.end(); it++) {
    fprintf(fp, "\t\t(%g, %g)\n", (*it).first, (*it).second);
  }
}


/// エラーメッセージ出力.
///
///   @param[in] fmt  出力フォーマット文字列
///
void DataHolder::printError(const char* fmt, ...)
{
  if (m_id == 0) {
    fprintf(stderr, "DataHolder error: label=%s, file=%s\n", 
            m_label.c_str(), m_file.c_str());
    va_list ap;
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
  }
}


/* ------------- DataHolderManager ----------------------------------*/


/// 設定XMLファイルポインタを受け取る.
///
///   @param[in] cfg 設定XMLファイルポインタ
///
bool DataHolderManager::receiveCfgPtr(SklSolverConfig* cfg)
{
  if (!cfg) return false;
  m_cfg = cfg;
  return true;
}


/// データファイル読み込み.
///
///   設定XMLファイルをパースし, 入力データファイルを読み込み,
///   DataHolderを生成し登録する.
///
void DataHolderManager::readData()
{
  const CfgElem* elemTop = m_cfg->GetTop(PARAMETER);
  if (!elemTop) assert(0);
  const CfgElem* elemDH = elemTop->GetElemFirst("Data_Holder");
  while (elemDH) {
    const CfgParam* param = elemDH->GetParamFirst();
    while (param) {
      const char* label = param->GetName();
      const char* file;
      if (!param->GetData(&file)) {
        printError("name=\"%s\": cannnot parse 'value'\n", label);
        assert(0);
      }
//    printf("label=%s, file=%s\n", label, file);
      DataHolder* dh = new DataHolder(label, file, pn.ID);
      m_dataHolders.insert(DATA_HOLDER_MAP::value_type(label, dh));

      param = elemDH->GetParamNext(param);
    }
    elemDH = elemTop->GetElemNext(elemDH, "Data_Holder");
  }
}


/// 管理しているDataHolderの情報を出力.
///
///   @param fp ファイルポインタ
///
void DataHolderManager::printInfo(FILE* fp) const
{
  fprintf(fp, "\n-----------------------------------------\n");
  fprintf(fp, "DataHolder infomation\n");
  fprintf(fp, "\t# of DataHolders = %d\n", m_dataHolders.size());

  DATA_HOLDER_MAP::const_iterator it;
  for (it = m_dataHolders.begin(); it != m_dataHolders.end(); it++) {
    it->second->printInfo(fp);
  }
}


/// 管理しているDataHolderの全データを出力(デバッグ用).
///
///   @param fp ファイルポインタ
///
void DataHolderManager::printInfoDebug(FILE* fp) const
{
  fprintf(fp, "\n-----------------------------------------\n");
  fprintf(fp, "DataHolder infomation\n");
  fprintf(fp, "\t# of DataHolders = %d\n", m_dataHolders.size());

  DATA_HOLDER_MAP::const_iterator it;
  for (it = m_dataHolders.begin(); it != m_dataHolders.end(); it++) {
    it->second->printInfoDebug(fp);
  }
}


/// エラーメッセージ出力.
///
///   @param[in] fmt  出力フォーマット文字列
///
void DataHolderManager::printError(const char* fmt, ...)
{
  if (pn.ID == 0) {
    fprintf(stderr, "DataHolderManager error:\n");
    va_list ap;
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
  }
}
