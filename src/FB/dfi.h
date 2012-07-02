#ifndef _CBC_DFI_H_
#define _CBC_DFI_H_

// #################################################################
//
// CAERU Library
//
// Copyright (c) All right reserved. 2012
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
//
// #################################################################

//@file DFI.h
//@brief DFI class Header
//@author kero
//@note このクラスは、並列時のみにコールすること。MPI_Initialize(), Finalize()はクラス外で実行。
//      dfiファイルの書き出しは、マスターランクでのみ行う。

#include <string>
#include <stdio.h>
#include <mpi.h>
#include "FB_Define.h"

class DFI {
protected:
  int Num_Node;          ///< MPI並列数
  int my_id;             ///< 自ノードのランク番号（dfiファイルの出力はランク0のみなので自明だが）
  int Gsize[3];          ///< 計算領域全体の分割数
  int div_domain[3];     ///< 全計算領域のノード分割数
  int procGrp;           ///< プロセスグループ（デフォルト0）
  int guide;             ///< ガイドセル数
  int start_type;        ///< セッションのスタートモード（initial_start=0, restart=1, coarse_restart=2）
  int* head;             ///< bboxの開始インデクス(C index) [3*Num_Node]
  int* tail;             ///< bboxの終端インデクス(C index) [3*Num_Node]
  std::string  hostname; ///< ラベル
  
public:
  DFI() {
    Num_Node   = 0;
    my_id      = 0;
    procGrp    = 0;
    guide      = 0;
    start_type = -1;
    head       = NULL;
    tail       = NULL;
    
    for (int i=0; i<3; i++) {
      Gsize[i]      = 0;
      div_domain[i] = 0;
    }

  }
  ~DFI() {
    if ( head ) delete [] head;
    if ( tail ) delete [] tail;
  }
  
protected:
  std::string Generate_DFI_Name(const std::string prefix, const int m_id);
  
  bool Write_File        (const std::string dfi_name, const std::string prefix, const int step, int& dfi_mng, const bool mio);
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
  bool init              (const int* g_size, const int* m_div, const int gc, const int stype, const int* hidx, const int* tidx);
  bool Write_DFI_File    (const std::string prefix, const int step, int& dfi_mng, const bool mio);
  
  std::string Generate_FileName(const std::string prefix, const int m_step, const int m_id, const bool mio=false);

};

#endif // _CBC_DFI_H_
