#ifndef _LAYOUT_H_
#define _LAYOUT_H_

// #################################################################
//
// output layoutfiles
//
// Copyright (c) All right reserved. 2012
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
//
// #################################################################

/** 
 * @file   layout.h
 * @brief  LAYOUT Class Header
 * @author kero
 */

#include <stdio.h>
//#include <stdlib.h>
#include <string.h>
#include <string>
#include <iostream>
#include <fstream>
#include <errno.h>
#include <sstream>

#include "mpi.h"

#include <sys/types.h>
#include <sys/stat.h>

#include <deque>

// for GetFullPathName
#if defined(_WIN32)
  #include <stdlib.h>
#else
  #include <sys/param.h>
  #include <stdlib.h>
#endif

#ifndef _WIN32
#include <unistd.h>           // for linux
#else
#include "sph_win32_util.h"   // for win32
#endif

#ifndef _WIN32
#include <dirent.h>
#else
#include "sph_win32_util.h"   // for win32
#endif

#include "PerfMonitor.h"
#include "cpm_ParaManager.h"

#include "FB_Define.h"
#include "mydebug.h"
#include "TPControl.h"
//#include "omp.h"

#include "dfiinfo.h"

#include "limits.h" // for UBUNTU

//#include "endianUtil.h"

//// FX10 profiler
//#if defined __K_FPCOLL
//#include "fjcoll.h"
//#elif defined __FX_FAPP
//#include "fj_tool/fjcoll.h"
//#endif

#include "LAYOUT_Define.h"

using namespace std;


class LAYOUT {

private:
  string fname;
  DfiInfo* DI;

  bool skip0;

  int ndfi;//number of dfi file list
  vector<string> dfi_name;
  vector<string> mname;
  vector<string> dname;
  vector<int> rankis;
  vector<int> rankie;

  string dirname;//出力ディレクトリ指定

public:
  /** コンストラクタ */
  LAYOUT();

  /**　デストラクタ */
  ~LAYOUT();

  /**  
   * @brief 引数のキープ
   * @param [in] m_fname  ファイル名
   */
  void SetInput(bool m_skip0, string m_fname);

  /**
   * @brief 入力ファイルの読み込み
   */
  void ReadInit();

  /**
   * @brief dfiファイルの読み込みとDfiInfoクラスデータの作成
   */
  void ReadDfiFiles();

  /**
   * @brief オプション入力からのディレクトリ名の反映
   */
  void SetDirName(string m_dname);

  /**
   * @brief オプション入力からのマシン名の反映
   */
  void SetMachineName(string m_mname);

  /**
   * @brief layoutファイルの出力
   */
  void OutputLayout();

  /**
   * @brief 出力指定ディレクトリのチェック
   */
  void CheckDir(string dirstr);

private:

  /**
   * @brief 入力ファイルの読み込み
   */
  void ReadInputFile(TPControl* tpCntl);

  /**
   * @brief layoutファイルの出力
   */
  void OutLay(DfiInfo *D);

  /**
   * @brief layoutファイル名を作成する
   * @param [in] prefix ファイル接頭文字
   * @param [in] str    xyz,func,,,
   * @param [in] m_step ステップ数
   */
  std::string Generate_LayoutFileName(const std::string prefix, const std::string str, const unsigned m_step);

  /**
   * @brief ファイル名を作成する。（拡張子自由）
   * @param [in] prefix ファイル接頭文字
   * @param [in] xxx    拡張子
   * @param [in] m_step ステップ数
   * @param [in] m_id   ランク番号
   * @param [in] mio    出力時の分割指定　 true = local / false = gather(default)
   */
  std::string Generate_FileName_Free(const std::string prefix, const std::string xxx, const unsigned m_step, const int m_id, const bool mio=false);

};

#endif // _LAYOUT_H_
