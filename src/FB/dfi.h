#ifndef _CBC_DFI_H_
#define _CBC_DFI_H_

/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file DFI.h
//@brief DFI class Header
//@author keno, Advanced Vis Team, AICS, RIKEN
//@note このクラスは、並列時のみにコールすること。MPI_Initialize(), Finalize()はクラス外で実行。
//      dfiファイルの書き出しは、マスターランクでのみ行う。

#include <string>
#include <stdio.h>
#include <mpi.h>
#include "FB_Define.h"

class DFI {
protected:
  int Num_Node;       /// MPI並列数
  int WriteCount;     /// FileInfoのファイル名を記述した回数（==時系列でコールされた回数）
  int my_id;          /// 自ノードのランク番号（dfiファイルの出力はランク0のみなので自明だが）
  int Gsize[3];       /// 計算領域全体の分割数
  int div_domain[3];  /// 全計算領域のノード分割数
  int procGrp;        /// プロセスグループ（デフォルト0）
  int guide;          /// ガイドセル数
  int* head;          /// bboxの開始インデクス(C index) [3*Num_Node]
  int* tail;          /// bboxの終端インデクス(C index) [3*Num_Node]
  char** hostname;    /// hostname char[m_np][LABEL]
  
public:
  DFI() {
    Num_Node   = 0;
    WriteCount = 0;
    my_id      = 0;
    procGrp    = 0;
    guide      = 0;
    head       = NULL;
    tail       = NULL;
    hostname   = NULL;
    
    for (int i=0; i<3; i++) {
      Gsize[i]      = 0;
      div_domain[i] = 0;
    }

  }
  ~DFI() {
    delete [] head;
    delete [] tail;
    for (int i=0; i<Num_Node; i++) {
      delete [] hostname[i];
    }
    delete [] hostname;
  }
  
protected:
  std::string Generate_DFI_Name(const std::string prefix, const int m_id);
  
  bool Write_File        (const std::string dfi_name, const std::string prefix, const int step, const bool mio);
  bool Write_Header      (FILE* fp, const unsigned tab, const std::string prefix);
  bool Write_Node        (FILE* fp, const unsigned tab, const int id, const std::string prefix);
  bool Write_NodeInfo    (FILE* fp, const unsigned tab, const std::string prefix);
  bool Write_OutFileInfo (FILE* fp, const unsigned tab, const std::string prefix, const int step, const bool mio);
  bool Write_OutFileName (FILE* fp, const unsigned tab, const std::string prefix, const int step, const int id, const bool mio);
  
  void Write_BaseName    (FILE* fp, const unsigned tab, const std::string prefix);
  void Write_FileFormat  (FILE* fp, const unsigned tab);
  void Write_GuideCell   (FILE* fp, const unsigned tab);
  void Write_MyID        (FILE* fp, const unsigned tab);
  void Write_NodeNum     (FILE* fp, const unsigned tab);
  void Write_Tab         (FILE* fp, const unsigned tab);
  void Write_NumDivDomain(FILE* fp, const unsigned tab);
  void Write_WholeSize   (FILE* fp, const unsigned tab);
  
public:
  bool init              (const int* g_size, const int* m_div, const int gc, const int* hidx, const int* tidx);
  bool Write_DFI_File    (const std::string prefix, const int step, const bool mio);
  
  std::string Generate_FileName(const std::string prefix, const int m_step, const int m_id, const bool mio=false);
  
  // @brief ホスト名のコピー
  // @param host
  void copy_hostname(const char* m_host, const int i) {
    strcpy(hostname[i], m_host);
  }

};

#endif // _CBC_DFI_H_
