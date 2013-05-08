/* ##################################################################
 *
 * PMlib - Performance Monitor library
 *
 * Copyright (c) 2010-2013 Advanced Institute for Computational Science, RIKEN.
 * All right reserved.
 *
 * ###################################################################
 */

//@file   PerfMonitor.cpp
//@brief  PerfMonitor class

#include "PerfMonitor.h"
#include <time.h>
#include <mpi.h>

namespace pm_lib {
  
  /// 測定レベル制御変数.
  /// =0:測定なし
  /// =1:排他測定のみ
  /// =2:非排他測定も(ディフォルト)
  unsigned PerfMonitor::TimingLevel = 2;
  
  
  /// 測定結果を出力.
  ///
  ///   排他測定区間のみ
  ///   @param[in] fp           出力ファイルポインタ
  ///   @param[in] hostname     ホスト名
  ///   @param[in] operatorname 作業者名
  ///
  ///   @note ノード0以外は, 呼び出されてもなにもしない
  ///
  void PerfMonitor::print(FILE* fp, const std::string hostname, const std::string operatorname)
  {
    if (my_rank != 0) return;
    
    if (!m_gathered) {
      fprintf(stderr, "\tPerfMonitor::print() error, not gathered infomation\n");
      PM_Exit(0);
    }
    
    // タイムスタンプの取得
    struct tm *date;
    time_t now;
    int year, month, day;
    int hour, minute, second;
    
    time(&now);
    date = localtime(&now);
    
    year   = date->tm_year + 1900;
    month  = date->tm_mon + 1;
    day    = date->tm_mday;
    hour   = date->tm_hour;
    minute = date->tm_min;
    second = date->tm_sec;
    
    
    
    // 合計 > tot
    double tot = 0.0;
    int maxLabelLen = 0;
    for (int i = 0; i < m_nWatch; i++) {
      if (m_watchArray[i].m_exclusive) {
        tot +=  m_watchArray[i].m_time_av;
        int labelLen = m_watchArray[i].m_label.size();
        maxLabelLen = (labelLen > maxLabelLen) ? labelLen : maxLabelLen;
      }
    }
    maxLabelLen++;
    
    fprintf(fp, "\n\t-----------------------------------------------------------------\n");
    fprintf(fp, "\tReport of Timing Statistics PMlib version %3.1f\n", (float)PMLIB_VERS/10.0);
    fprintf(fp, "\n");
    fprintf(fp, "\tOperator  : %s\n", operatorname.c_str());
    fprintf(fp, "\tHost name : %s\n", hostname.c_str());
    fprintf(fp, "\tDate      : %04d/%02d/%02d : %02d:%02d:%02d\n", year, month, day, hour, minute, second);
    fprintf(fp, "\n");
    
    fprintf(fp,"\tParallel Mode                    :   %s ", parallel_mode.c_str());
    if (parallel_mode == "Serial") {
      fprintf(fp, "\n");
    } else if (parallel_mode == "FlatMPI") {
      fprintf(fp, "(%d processes)\n", num_process);
    } else if (parallel_mode == "OpenMP") {
      fprintf(fp, "(%d threads)\n", num_threads);
    } else if (parallel_mode == "Hybrid") {
      fprintf(fp, "(%d processes x %d threads)\n", num_process, num_threads);
    } else {
      fprintf(fp, "\n\tError : invalid Parallel mode \n");
      PM_Exit(0);
    }
    fprintf(fp, "\n");
    
    fprintf(fp, "\tTotal execution time            = %12.6e [sec]\n", m_total.m_time);
    fprintf(fp, "\tTotal time of measured sections = %12.6e [sec]\n", tot);
    fprintf(fp, "\n");
    
    fprintf(fp, "\tStatistics per MPI process [Node Average]\n");
    
    fprintf(fp, "\t%-*s|     call     |              accumulated time                    |          flop | messages[Bytes]\n", maxLabelLen, "Label");
    fprintf(fp, "\t%-*s|              |   avr[sec]    avr[%%]     sdv[sec]  avr/call[sec] |   avr         sdv         speed\n", maxLabelLen, "");
    fputc('\t', fp); for (int i = 0; i < maxLabelLen; i++) fputc('-', fp);
    fprintf(fp, "+--------------+--------------------------------------------------+--------------------------------------\n");
    
    // 反復あたりのタイムコスト
    double *m_tcost = NULL;
    if ( !(m_tcost = new double[m_nWatch]) ) PM_Exit(0);
    
    for (int i = 0; i < m_nWatch; i++) {
      PerfWatch& w = m_watchArray[i];
      if (!w.m_exclusive) continue;   // 
      if (w.m_label.empty()) continue;  //
      if (w.m_valid) {
        if ( w.m_count > 0 ) {
          m_tcost[i] = 100*w.m_time_av/tot;  // ノード平均/各排他測定セクションの実行時間のノード間平均値の合計値
        }
      }
    }
    
    // 降順ソート O(n^2) Brute force　
    unsigned *m_order = NULL;
    if ( !(m_order = new unsigned[m_nWatch]) ) PM_Exit(0);
    
    for (int i=0; i<m_nWatch; i++) {
      m_order[i] = i;
    }
    
    double tmp_d;
    unsigned tmp_u;
    for (int i=0; i<m_nWatch-1; i++) {
      PerfWatch& w = m_watchArray[i];
      if (!w.m_exclusive) continue;   // 
      if (w.m_label.empty()) continue;  //
      if (w.m_valid) {
        for (int j=i+1; j<m_nWatch; j++) {
          PerfWatch& q = m_watchArray[j];
          if (!q.m_exclusive) continue;   // 
          if (q.m_label.empty()) continue;  //
          if (q.m_valid) {
            if ( m_tcost[i] < m_tcost[j] ) {
              tmp_d=m_tcost[i]; m_tcost[i]=m_tcost[j]; m_tcost[j]=tmp_d;
              tmp_u=m_order[i]; m_order[i]=m_order[j]; m_order[j]=tmp_u;
            }
          }
        }
      }
    }
    
    
    // 表示
    double sum_time_av = 0.0;
    double sum_flop_av = 0.0;
    double sum_time_comm = 0.0;
    std::string unit;
    double flop_serial=0.0;
    
    /* 登録順表示
     for (int i = 0; i < m_nWatch; i++) {
     PerfWatch& w = m_watchArray[i];
     if (!w.m_exclusive) continue;   // 
     if (w.m_label.empty()) continue;  //
     if (w.m_valid) {
     if ( w.m_count > 0 ) {
     double flops_av = (w.m_count==0) ? 0.0 : w.m_flop_av/w.m_time_av;
     fprintf(fp, "\t%-*s: %12ld   %12.6e  %6.2f %5d %11.4e  %12.6e    %8.3e   %8.3e  %7.2f %s\n",
     maxLabelLen,
     w.m_label.c_str(), 
     w.m_count,            // コール回数
     w.m_time_av,          // ノード平均
     100*w.m_time_av/tot,  // ノード平均/各排他測定セクションの実行時間のノード間平均値の合計値
     m_order[i],           // 順位
     w.m_time_sd,          // 標準偏差
     (w.m_count==0) ? 0.0 : w.m_time_av/w.m_count, // 1回あたりの時間コスト
     w.m_flop_av,          // ノード平均
     w.m_flop_sd,          // 標準偏差
     w.flops(flops_av, unit, w.get_typeCalc()), unit.c_str());
     sum_time_av += w.m_time_av;
     // 計算セクションのみ
     if ( w.m_typeCalc ) sum_flop_av += w.m_flop_av;
     
     // 非計算セクションのみ
     if ( !w.m_typeCalc ) sum_time_comm += w.m_time_av; //w.m_time_comm;
     }
     }
     else {
     fprintf(fp, "\t%-*s: *** NA ***\n", maxLabelLen, w.m_label.c_str());
     }
     }
     */
    
    // タイムコスト順表示
    for (int j = 0; j < m_nWatch; j++) {
      int i = m_order[j];
      PerfWatch& w = m_watchArray[i];
      if (!w.m_exclusive) continue;   // 
      if (w.m_label.empty()) continue;  //
      if (w.m_valid) {
        if ( w.m_count > 0 ) {
          double flops_av;
          if (w.m_time_av == 0.0) {
            flops_av = 0.0;
          }
          else {
            flops_av = (w.m_count==0) ? 0.0 : w.m_flop_av/w.m_time_av;
          }
          fprintf(fp, "\t%-*s: %12ld   %12.6e  %6.2f  %11.4e  %12.6e    %8.3e   %8.3e  %7.2f %s\n",
                  maxLabelLen,
                  w.m_label.c_str(), 
                  w.m_count,            // コール回数
                  w.m_time_av,          // ノード平均
                  100*w.m_time_av/tot,  // ノード平均/各排他測定セクションの実行時間のノード間平均値の合計値
                  w.m_time_sd,          // 標準偏差
                  (w.m_count==0) ? 0.0 : w.m_time_av/w.m_count, // 1回あたりの時間コスト
                  w.m_flop_av,          // ノード平均
                  w.m_flop_sd,          // 標準偏差
                  w.flops(flops_av, unit, w.get_typeCalc()), unit.c_str());
          sum_time_av += w.m_time_av;
          
          // 計算セクションのみ
          if ( w.m_typeCalc ) sum_flop_av += w.m_flop_av;
          
          // 非計算セクションのみ
          if ( !w.m_typeCalc ) sum_time_comm += w.m_time_av; //w.m_time_comm;
        }
      }
      else {
        fprintf(fp, "\t%-*s: *** NA ***\n", maxLabelLen, w.m_label.c_str());
      }
    }
    
    fputc('\t', fp); for (int i = 0; i < maxLabelLen; i++) fputc('-', fp);
    fprintf(fp, "+--------------+--------------------------------------------------+--------------------------------------\n");
    flop_serial = PerfWatch::flops(sum_flop_av/sum_time_av, unit, true);
    fprintf(fp, "\t%-*s|%15s %12.6e                                       %8.3e              %7.2f %s\n", maxLabelLen, "Total", "",
            sum_time_av, sum_flop_av, flop_serial, unit.c_str());
    
    // 並列時
    int np;
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    
    if ( np > 1 ) {
      fprintf(fp,"\t%-*s| \t\t\t                                                                   %7.2f %s\n", 
              maxLabelLen, "   Performance", flop_serial*(double)np, unit.c_str());
    }
  }
  
  
  /// 詳細な測定結果を出力.
  ///
  ///   ノード毎に非排他測定区間も出力
  ///   @param[in] fp 出力ファイルポインタ
  ///
  ///   @note ノード0以外は, 呼び出されてもなにもしない
  ///
  void PerfMonitor::printDetail(FILE* fp)
  {
    if (my_rank != 0) return;
    
    if (!m_gathered) {
      fprintf(stderr, "\tPerfMonitor::printDetail() error, not gathered infomation\n");
      PM_Exit(0);
    }
    
    fprintf(fp, "\n-----------------------------------------------------------------\n");
    fprintf(fp, "Detail of Timing Statistics\n");
    
    double tot = 0.0;
    for (int i = 0; i < m_nWatch; i++) {
      if (m_watchArray[i].m_exclusive) {
        tot +=  m_watchArray[i].m_time_av;
      }
    }
    for (int i = 0; i < m_nWatch; i++) {
      m_watchArray[i].printDatail(fp, tot);
    }
    
  }
  
} /* namespace pm_lib */
