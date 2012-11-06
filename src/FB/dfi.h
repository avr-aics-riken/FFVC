#ifndef _FB_DFI_H_
#define _FB_DFI_H_

// #################################################################
//
// CAERU Library
//
// Copyright (c) 2012 All right reserved.
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
//
// #################################################################

//@file   dfi.h
//@brief  DFI class Header
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
  int start_type;        ///< セッションのスタートモード（initial_start=0, restart=1, coarse_restart=2, restart_different_nproc=3）
  int* head;             ///< bboxの開始インデクス(C index) [3*Num_Node]
  int* tail;             ///< bboxの終端インデクス(C index) [3*Num_Node]
  std::string hostname;  ///< ラベル
  
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
  
  /**
   * @brief 出力DFIファイル名を作成する
   * @param [in] prefix ファイル接頭文字
   */
  std::string Generate_DFI_Name(const std::string prefix);
  
  
  /**
   * @brief DFIファイル:BaseName要素を出力する
   * @param [in] fp      ファイルポインタ
   * @param [in] tab     インデント
   * @param [in] prefix  ファイル接頭文字
   */
  void Write_BaseName(FILE* fp, const unsigned tab, const std::string prefix);

  
  /**
   * @brief DFIファイルを出力する
   * @param [in] dfi_name  DFIファイル名
   * @param [in] prefix    ファイル接頭文字
   * @param [in] step      ステップ数
   * @param [in] dfi_mng   出力管理カウンタ
   * @param [in] mio       出力時の分割指定　 true = local / false = gather
   */
  bool Write_File(const std::string dfi_name, const std::string prefix, const unsigned step, int& dfi_mng, const bool mio);
  
  /**
   * @brief DFIファイルを出力する
   * @param [in] dfi_name  DFIファイル名
   * @param [in] prefix    ファイル接頭文字
   * @param [in] step      ステップ数
   * @param [in] time      時間
   * @param [in] dfi_mng   出力管理カウンタ
   * @param [in] mio       出力時の分割指定　 true = local / false = gather
   */
  bool Write_File(const std::string dfi_name, const std::string prefix, const unsigned step, const double time, int& dfi_mng, const bool mio);
  
  
  /**
   * @brief DFIファイル:ファイルフォーマット要素を出力する
   * @param [in] fp      ファイルポインタ
   * @param [in] tab     インデント
   */
  void Write_FileFormat(FILE* fp, const unsigned tab);
  
  
  /**
   * @brief DFIファイル:ヘッダー要素を出力する
   * @param [in] fp      ファイルポインタ
   * @param [in] tab     インデント
   * @param [in] prefix  ファイル接頭文字
   */
  bool Write_Header(FILE* fp, const unsigned tab, const std::string prefix);
  
  
  /**
   * DFIファイル:ボクセル情報要素を出力する
   * @param [in] fp      ファイルポインタ
   * @param [in] tab     インデント
   * @param [in] n       対象ノードID
   * @param [in] prefix  ファイル接頭文字
   */
  bool Write_Node(FILE* fp, const unsigned tab, const int id, const std::string prefix);
  
  
  /**
   * @brief DFIファイル:ノード情報要素を出力する
   * @param [in] fp      ファイルポインタ
   * @param [in] tab     インデント
   * @param [in] prefix  ファイル接頭文字
   */
  bool Write_NodeInfo(FILE* fp, const unsigned tab, const std::string prefix);
  
  
  /**
   * @brief DFIファイル:出力ファイル情報要素を出力する
   * @param [in] fp      ファイルポインタ
   * @param [in] prefix  ファイル接頭文字
   * @param [in] tab     インデント
   * @param [in] step    ステップ数
   * @param [in] mio     出力時の分割指定　 true = local / false = gather
   */
  bool Write_OutFileInfo(FILE* fp, const unsigned tab, const std::string prefix, const unsigned step, const bool mio);
  
  
  /**
   * @brief DFIファイル:出力ファイル情報要素を出力する
   * @param [in] fp      ファイルポインタ
   * @param [in] prefix  ファイル接頭文字
   * @param [in] tab     インデント
   * @param [in] step    ステップ数
   * @param [in] time    時間
   * @param [in] mio     出力時の分割指定　 true = local / false = gather
   */
  bool Write_OutFileInfo(FILE* fp, const unsigned tab, const std::string prefix, const unsigned step, const double time, const bool mio);
  
  
  /**
   * @brief DFIファイル:ファイル名要素を出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   * @param [in] prefix ファイル接頭文字
   * @param [in] step   ステップ数
   * @param [in] id     対象ノードID
   * @param [in] mio    出力時の分割指定　 true = local / false = gather
   */
  bool Write_OutFileName(FILE* fp, const unsigned tab, const std::string prefix, const unsigned step, const int id, const bool mio);
  

  /**
   * @brief DFIファイル:ガイドセル要素を出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   * @param [in] gc     ガイドセル
   */
  void Write_GuideCell(FILE* fp, const unsigned tab);
  
  
  /**
   * @brief DFIファイル:ノード番号要素を出力する
   * @param [in] fp      ファイルポインタ
   * @param [in] tab     インデント
   */
  void Write_MyID(FILE* fp, const unsigned tab);
  
  
  /**
   * @brief DFIファイル:ノード数要素を出力する
   * @param [in] fp  ファイルポインタ
   * @param [in] tab インデント
   */
  void Write_NodeNum(FILE* fp, const unsigned tab);
  
  
  /**
   * @brief DFIファイル:I,J,K分割数要素を出力する
   * @param [in] fp      ファイルポインタ
   * @param [in] tab     インデント
   */
  void Write_NumDivDomain(FILE* fp, const unsigned tab);
  
  
  /**
   * @brief ステップ数を出力する
   * @param [in] fp      ファイルポインタ
   * @param [in] tab     インデント数
   * @param [in] step    ステップ数
   */
  void Write_Step(FILE* fp, const unsigned tab, const unsigned step);
  
  
  /**
   * @brief Tab(space２つ)を出力する
   * @param [in] fp      ファイルポインタ
   * @param [in] tab     インデント数
   */
  void Write_Tab(FILE* fp, const unsigned tab);
  
  
  /**
   * @brief DFIファイル:全体ボクセルサイズ要素を出力する
   * @param [in] fp      ファイルポインタ
   * @param [in] tab     インデント
   */
  void Write_WholeSize(FILE* fp, const unsigned tab);
  
  
  
public:
  
  /**
   * @brief ディレクトリ名を作成する
   * @param [in] prefix ファイル接頭文字
   * @param [in] m_step ステップ数
   */
  std::string Generate_DirName(const std::string prefix, const unsigned m_step);
  
  
  /**
   * @brief ファイル名を作成する
   * @param [in] prefix ファイル接頭文字
   * @param [in] m_step ステップ数
   * @param [in] m_id   ランク番号
   * @param [in] mio    出力時の分割指定　 true = local / false = gather(default)
   */
  std::string Generate_FileName(const std::string prefix, const unsigned m_step, const int m_id, const bool mio=false);
  

  /**
   * @brief ファイル名を作成する。（拡張子自由）
   * @param [in] prefix ファイル接頭文字
   * @param [in] xxx    拡張子
   * @param [in] m_step ステップ数
   * @param [in] m_id   ランク番号
   * @param [in] mio    出力時の分割指定　 true = local / false = gather(default)
   */
  std::string Generate_FileName_Free(const std::string prefix, const std::string xxx, const unsigned m_step, const int m_id, const bool mio=false);

  
  /**
   * @brief 初期化
   * @param [in] g_size グローバルサイズ
   * @param [in] m_div  ノード分割数
   * @param [in] gc     ガイドセル
   * @param [in] stype  スタートタイプ
   * @param [in] hidx   開始インデクス
   * @param [in] tidx   終端インデクス
   * @param [in] m_host ホスト名
   */
  bool init(const int* g_size, const int* m_div, const int gc, const int stype, const int* hidx, const int* tidx, const std::string m_host);
  
  
  /**
   * @brief データをファイルに書き込む
   * @param [in] prefix  ファイル接頭文字
   * @param [in] step    ステップ
   * @param [in] dfi_mng 出力管理カウンタ
   * @param [in] mio     出力時の分割指定　 true = local / false = gather
   */
  bool Write_DFI_File(const std::string prefix, const unsigned step, int& dfi_mng, const bool mio);
  
  
  /**
   * @brief データをファイルに書き込む
   * @param [in] prefix  ファイル接頭文字
   * @param [in] step    ステップ
   * @param [in] time    時間
   * @param [in] dfi_mng 出力管理カウンタ
   * @param [in] mio     出力時の分割指定　 true = local / false = gather
   */
  bool Write_DFI_File(const std::string prefix, const unsigned step, const double time, int& dfi_mng, const bool mio);
};

#endif // _FB_DFI_H_
