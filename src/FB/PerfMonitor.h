#ifndef _SKL_FB_PERFMONITOR_H_
#define _SKL_FB_PERFMONITOR_H_

/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file PerfMonitor.h
//@brief FlowBase PerfMonitor class Header
//@author keno, FSI Team, VCAD, RIKEN

#include <cstdio>
using namespace std;

#include "FBDefine.h"
#include "Parallel_node.h"
#include "PerfWatch.h"


/// 排他測定用マクロ
#define PM_TIMING__         if (PerfMonitor::TimingLevel > 0)


/// 排他測定＋非排他測定用マクロ
#define PM_TIMING_DETAIL__  if (PerfMonitor::TimingLevel > 1)


/**
 * 計算性能測定管理クラス.
 */
class PerfMonitor : public Parallel_Node {
public:

  /// 測定対象タイプ
  enum Type {
    COMM,  ///< 通信
    CALC,  ///< 計算
  };

  /// 測定レベル制御変数.
  /// =0:測定なし/=1:排他測定のみ/=2:非排他測定も(ディフォルト)
  static unsigned TimingLevel;

private:
  unsigned m_nWatch;        ///< 測定区間数

  PerfWatch* m_watchArray;  ///< 測定時計配列

  PerfWatch  m_total;       ///< 全計算時間用測定時計

  bool m_gathered;          ///< 想定結果集計済みフラグ

public:
  /// コンストラクタ.
  PerfMonitor() : m_watchArray(0) {}

  /// デストラクタ.
  ~PerfMonitor() { if (m_watchArray) delete[] m_watchArray; }

  /// 初期化.
  ///
  /// 測定区間数分の測定時計を準備.
  /// 全計算時間用測定時計をスタート.
  /// @param[in] nWatch 測定区間数
  ///
  void initialize(unsigned nWatch) {
    m_nWatch = nWatch;
    m_watchArray = new PerfWatch[m_nWatch];
    m_gathered = false;

    m_total.setProperties("Total excution time", CALC, false);
    m_total.start();
  }

  /// 測定時計にプロパティを設定.
  ///
  ///   @param[in] key キー番号
  ///   @param[in] label ラベル
  ///   @param[in] type  測定対象タイプ(COMM or CALC)
  ///   @param[in] exclusive 排他測定フラグ(ディフォルトtrue)
  ///
  void setProperties(unsigned key, const string& label, Type type, bool exclusive=true) {
    if (key >= m_nWatch) {
      fprintf(stderr, "\tPerfMonitor::setProperties() error, out of range key\n");
      Exit(0);
    }
    m_watchArray[key].setProperties(label, type, exclusive);
    m_watchArray[key].setParallelInfo(pn);
  }

  /// 測定スタート.
  ///
  ///   @param[in] key キー番号
  ///
  void start(unsigned key) {
    if (key >= m_nWatch) {
      fprintf(stderr, "\tPerfMonitor::start() error, out of range key\n");
      Exit(0);
    }
    m_watchArray[key].start();
  }

  /// 測定ストップ.
  ///
  ///   @param[in] key キー番号
  ///   @param[in] flopPerTask 「タスク」あたりの計算量/通信量(バイト) (ディフォルト0)
  ///   @param[in] iterationCount  実行「タスク」数 (ディフォルト1)
  ///
  void stop(unsigned key, REAL_TYPE flopPerTask=0.0, unsigned iterationCount=1) {
    if (key >= m_nWatch) {
      fprintf(stderr, "\tPerfMonitor::stop() error, out of range key\n");
      Exit(0);
    }
    m_watchArray[key].stop(flopPerTask, iterationCount);
  }

  /// 測定結果情報をノード０に集約.
  ///
  ///   全計算時間用測定時計をストップ.
  ///
  void gather() {
    if (m_gathered) {
      fprintf(stderr, "\tPerfMonitor::gather() error, already gathered\n");
      Exit(0);
    }
    m_total.stop(0.0, 1);
    for (int i = 0; i < m_nWatch; i++) {
      m_watchArray[i].gather();
    }
    m_gathered = true;
  }

  /// 測定結果を出力.
  ///
  ///   排他測定区間のみ
  ///   @param[in] fp 出力ファイルポインタ
  ///
  ///   @note ノード0以外は, 呼び出されてもなにもしない
  ///
  void print(FILE* fp);

  /// 詳細な測定結果を出力.
  ///
  ///   ノード毎に非排他測定区間も出力
  ///   @param[in] fp 出力ファイルポインタ
  ///
  ///   @note ノード0以外は, 呼び出されてもなにもしない
  ///
  void printDetail(FILE* fp);

};

#endif // _SKL_FB_PERFMONITOR_H_
