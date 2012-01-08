/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2010
 *
 */

//@file PerfWatch.C
//@brief FlowBase PerfWatch class
//@author keno, FSI Team, VCAD, RIKEN

#include <cmath>
#include "PerfWatch.h"

bool PerfWatch::ExclusiveStarted = false;


/// 測定スタート.
void PerfWatch::start()
{
  if (m_label.empty()) {
    printError("PerfWatch::start()",  "properties not set\n");
    assert(0);
  }
  if (m_started) {
    printError("PerfWatch::start()",  "already started\n");
    assert(0);
  }
  if (m_exclusive && ExclusiveStarted) {
    printError("PerfWatch::start()",  "exclusive sections overlapped\n");
    assert(0);
  }
  m_started = true;
  if (m_exclusive) ExclusiveStarted = true;
  m_startTime = getTime();
}


/// 測定ストップ.
///
///   @param[in] flopPerTask 「タスク」あたりの計算量/通信量(バイト)
///   @param[in] iterationCount  実行「タスク」数
///
///   @note m_countには, iterationCountではなく, 「測定回数」を積算
///
void PerfWatch::stop(SKL_REAL flopPerTask, unsigned iterationCount)
{
  if (!m_started) {
    printError("PerfWatch::stop()",  "not started\n");
    assert(0);
  }
  m_time += getTime() - m_startTime;
  m_flop += (double)flopPerTask * iterationCount;
  m_count++;
  m_started = false;
  if (m_exclusive) ExclusiveStarted = false;
}


/// 測定結果情報をノード０に集約.
void PerfWatch::gather()
{
  if (m_gathered) {
    printError("PerfWatch::gather()",  "already gathered\n");
    assert(0);
  }
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  m_np = para_mng->GetNodeNum(pn.procGrp);
  
  if (!(m_timeArray = new double[m_np])) assert(0);
  if (!(m_flopArray = new double[m_np])) assert(0);
  if (!(m_countArray = new unsigned long[m_np])) assert(0);
  if (para_mng->IsParallel()) {
    if (!para_mng->Gather(&m_time, 1, SKL_ARRAY_DTYPE_DOUBLE,
                      m_timeArray, 1, SKL_ARRAY_DTYPE_DOUBLE, 0, pn.procGrp)) assert(0);
    if (!para_mng->Gather(&m_flop, 1, SKL_ARRAY_DTYPE_DOUBLE,
                      m_flopArray, 1, SKL_ARRAY_DTYPE_DOUBLE, 0, pn.procGrp)) assert(0);
    if (!para_mng->Gather(&m_count, 1, SKL_ARRAY_DTYPE_ULONG,
                      m_countArray, 1, SKL_ARRAY_DTYPE_ULONG, 0, pn.procGrp)) assert(0);
  }
  else {
    m_timeArray[0] = m_time;
    m_flopArray[0] = m_flop;
    m_countArray[0] = m_count;
  }

  m_gathered = true;

  if (pn.ID == 0) {
   
    if (m_exclusive) {

      // check call count
      int n = m_countArray[0];
      for (int i = 1; i < m_np; i++) {
        if (m_countArray[i] != n) m_valid = false;
      }
      // 排他測定(m_exclusive==true)では、全ノードでの測定回数が等しいことを仮定.
      // 等しくない場合(m_valid==fals)には、統計計算をスキップする.
      if (!m_valid) return;
  
      // 平均値
      m_time_av = 0.0;
      m_flop_av = 0.0;
      for (int i = 0; i < m_np; i++) {
        m_time_av += m_timeArray[i];
        m_flop_av += m_flopArray[i];
      }
      m_time_av /= m_np;
      m_flop_av /= m_np;
      
      // 標準偏差
      m_time_sd = 0.0;
      m_flop_sd = 0.0;
      if (m_np > 1) {
        for (int i = 0; i < m_np; i++) {
          double d_time = m_timeArray[i] - m_time_av;
          double d_flop = m_flopArray[i] - m_flop_av;
          m_time_sd += d_time*d_time;
          m_flop_sd += d_flop*d_flop;
        }
        m_time_sd = sqrt(m_time_sd / (m_np-1));
        m_flop_sd = sqrt(m_flop_sd / (m_np-1));
      }
      
      // 通信の場合，各ノードの通信時間の最大値
      m_time_comm = 0.0;
      if (!m_typeCalc) {
        double comm_max = 0.0;
        for (int i = 0; i < m_np; i++) {
          if (m_timeArray[i] > comm_max) comm_max = m_timeArray[i];
        }
        m_time_comm = comm_max; 
      }
    }

  } // pn.ID == 0

}


/// 詳細な測定結果を出力.
///
///   ノード毎に非排他測定区間も出力
///
///   @param[in] fp 出力ファイルポインタ
///   @param[in] totalTime 全排他測定区間での計算時間(平均値)の合計
///
///   @note ノード0からのみ, 呼び出し可能
///
void PerfWatch::printDatail(FILE* fp, double totalTime)
{
  if (pn.ID != 0) {
    printError("PerfWatch::printDetail()",  "call from non-root node\n");
    assert(0);
  }
  if (!m_gathered) {
    printError("PerfWatch::printDetail()",  "nt gathered information\n");
    assert(0);
  }

  double tMax = 0.0;
  for (int i = 0; i < m_np; i++) {
    tMax = (m_timeArray[i] > tMax) ? m_timeArray[i] : tMax;
  }

  // skip zero count
  unsigned long total_count = 0;
  for (int i = 0; i < m_np; i++) total_count += m_countArray[i];
  
  if ( total_count > 0 ) {
    fprintf(fp, "\n%s%s\n", m_exclusive ? "" : "*", m_label.c_str());
    fprintf(fp, "                call       accm[s] accm[%%]    waiting[s]  accm/call[s]    flop|msg     speed\n");
    for (int i = 0; i < m_np; i++) {
      string unit;
      fprintf(fp, "#%-5d: %12ld  %12.6e  %6.2f  %12.6e  %12.6e   %8.3e  %7.2f %s\n",
              i, 
              m_countArray[i], // コール回数
              m_timeArray[i], // ノードあたりの時間
              100*m_timeArray[i]/totalTime, // 非排他測定区間に対する割合
              tMax-m_timeArray[i], // ノード間の最大値を基準にした待ち時間
              (m_countArray[i]==0) ? 0.0: m_timeArray[i]/m_countArray[i], // 1回あたりの時間コスト
              m_flopArray[i], // ノードあたりの演算数
              flops((m_countArray[i]==0) ? 0.0 : m_flopArray[i]/m_timeArray[i], unit, m_typeCalc), unit.c_str() // スピード
              );
    }
  }
}


/// 単位変換.
///
///   @param[in] fops 浮動小数演算数/通信量(バイト)
///   @param[out] unit 単位の文字列
///   @param[in] mode 測定モード (通信か計算)
///   @return  単位変換後の数値
///
double PerfWatch::flops(double fops, string &unit, bool mode)
{
  double P, T, G, M, K, ret=0.0;
  K = 1024.0;
  M = 1024.0*K;
  G = 1024.0*M;
  T = 1024.0*G;
  P = 1024.0*T;
    
  if (mode) { // arithmetic - ture
    if      ( fops > P ) {
      ret = fops / P;
      unit = "Pflops";
    }
    else if ( fops > T ) {
      ret = fops / T;
      unit = "Tflops";
    }
    else if ( fops > G ) {
      ret = fops / G;
      unit = "Gflops";
    }
    else {
      ret = fops / M;
      unit = "Mflops";
    }
  }
  else { // communication - false
    if      ( fops > P ) {
      ret = fops / P;
      unit = "PB/sec";
    }
    else if ( fops > T ) {
      ret = fops / T;
      unit = "TB/sec";
    }
    else if ( fops > G ) {
      ret = fops / G;
      unit = "GB/sec";
    }
    else if ( fops > M ) {
      ret = fops / M;
      unit = "MB/sec";
    }
    else {
      ret = fops / K;
      unit = "KB/sec";
    }
  }
  return ret;
}


/// 時刻を取得.
///
///   Unix/Linux: gettimeofdayシステムコールを使用.
///   Windows: GetSystemTimeAsFileTime API(sph_win32_util.h)を使用.
///
///   @return 時刻値(秒)
///
double PerfWatch::getTime()
{
  struct timeval tv;
  gettimeofday(&tv, 0);
  return (double)tv.tv_sec + (double)tv.tv_usec * 1.0e-6;
}


/// エラーメッセージ出力.
///
///   @param[in] func メソッド名
///   @param[in] fmt  出力フォーマット文字列
///
void PerfWatch::printError(const char* func, const char* fmt, ...)
{
  if (pn.ID == 0) {
    fprintf(stderr, "%s error: \"%s\" ", func, m_label.c_str());
    va_list ap;
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
  }
}
