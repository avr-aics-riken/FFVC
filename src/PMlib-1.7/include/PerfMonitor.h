#ifndef _CAERU_PERFMONITOR_H_
#define _CAERU_PERFMONITOR_H_

/* ##################################################################
 *
 * PMlib - Performance Monitor library
 *
 * Copyright (c) 2010-2013 Advanced Institute for Computational Science, RIKEN.
 * All right reserved.
 *
 * ###################################################################
 */

//@file   PerfMonitor.h
//@brief  PerfMonitor class Header
//@note MPI_Init(), MPI_Finalize()はライブラリ外で行う。

#include <cstdio>
#include <cstdlib>
#include "PerfWatch.h"

namespace pm_lib {
  
  /// 排他測定用マクロ
#define PM_TIMING__         if (PerfMonitor::TimingLevel > 0)
  
  
  /// 排他測定＋非排他測定用マクロ
#define PM_TIMING_DETAIL__  if (PerfMonitor::TimingLevel > 1)

  /// バージョン情報
#define PMLIB_VERS 17
  
  /**
   * 計算性能測定管理クラス.
   */
  class PerfMonitor {
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
    unsigned m_nWatch;         ///< 測定区間数
    bool m_gathered;           ///< 想定結果集計済みフラグ
    int num_threads;           ///< 並列スレッド数
    int num_process;           ///< 並列プロセス数
    int my_rank;               ///< 自ランク番号
    std::string parallel_mode; ///< 並列動作モード（"Serial", "OpenMP", "Flat MPI", "Hybrid"）
    PerfWatch* m_watchArray;   ///< 測定時計配列
    PerfWatch  m_total;        ///< 全計算時間用測定時計
    
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
      my_rank = 0;
      
      m_total.setProperties("Total excution time", CALC, my_rank, false);
      m_total.start();
    }
    
    /// ランク番号の通知
    void setRankInfo(const int myID) {
      my_rank = myID;
    }
    
    /// 測定時計にプロパティを設定.
    ///
    ///   @param[in] key キー番号
    ///   @param[in] label ラベル
    ///   @param[in] type  測定対象タイプ(COMM or CALC)
    ///   @param[in] exclusive 排他測定フラグ(ディフォルトtrue)
    ///
    void setProperties(unsigned key, const std::string& label, Type type, bool exclusive=true) {
      if (key >= m_nWatch) {
        fprintf(stderr, "\tPerfMonitor::setProperties() error, out of range key\n");
        PM_Exit(0);
      }
      m_watchArray[key].setProperties(label, type, my_rank, exclusive);
    }
    
    /// 並列モードを設定
    ///
    /// @param[in] p_mode 並列モード
    /// @param[in] n_thread
    /// @param[in] n_proc
    ///
    void setParallelMode(const std::string& p_mode, const int n_thread, const int n_proc) {
      parallel_mode = p_mode;
      num_threads   = n_thread;
      num_process   = n_proc;
    }
    
    /// 測定スタート.
    ///
    ///   @param[in] key キー番号
    ///
    void start(unsigned key) {
      if (key >= m_nWatch) {
        fprintf(stderr, "\tPerfMonitor::start() error, out of range key\n");
        PM_Exit(0);
      }
      m_watchArray[key].start();
    }
    
    /// 測定ストップ.
    ///
    ///   @param[in] key キー番号
    ///   @param[in] flopPerTask 「タスク」あたりの計算量/通信量(バイト) (ディフォルト0)
    ///   @param[in] iterationCount  実行「タスク」数 (ディフォルト1)
    ///
    void stop(unsigned key, double flopPerTask=0.0, unsigned iterationCount=1) {
      if (key >= m_nWatch) {
        fprintf(stderr, "\tPerfMonitor::stop() error, out of range key\n");
        PM_Exit(0);
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
        PM_Exit(0);
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
    ///   @param[in] fp           出力ファイルポインタ
    ///   @param[in] hostname     ホスト名
    ///   @param[in] operatorname 作業者名
    ///
    ///   @note ノード0以外は, 呼び出されてもなにもしない
    ///
    void print(FILE* fp, const std::string hostname, const std::string operatorname);
    
    /// 詳細な測定結果を出力.
    ///
    ///   ノード毎に非排他測定区間も出力
    ///   @param[in] fp 出力ファイルポインタ
    ///
    ///   @note ノード0以外は, 呼び出されてもなにもしない
    ///
    void printDetail(FILE* fp);
    
  };

} /* namespace pm_lib */

#endif // _CAERU_PERFMONITOR_H_
