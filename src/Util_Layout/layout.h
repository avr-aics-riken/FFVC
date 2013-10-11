#ifndef _LAYOUT_H_
#define _LAYOUT_H_

//##################################################################################
//
// FFV-C : Frontflow / violet Cartesian
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2013 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2013 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

/** 
 * @file   layout.h
 * @brief  LAYOUT Class Header
 * @author kero
 */

#include "mpi.h"
#include <stdio.h>
//#include <stdlib.h>
#include <string.h>
#include <string>
#include <iostream>
#include <fstream>
#include <errno.h>
#include <sstream>

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
//#include "omp.h"

#include "../PLOT3D/dfiinfo.h"

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
  
public:
  cpm_ParaManager* paraMngr; ///< Cartesian Partition Manager
  
private:
  string fname;
  DfiInfo* DI;
  
  bool skip0;
  
  int ndfi;  //number of dfi file list
  vector<string> dfi_name;
  vector<string> mname;
  vector<string> dname;
  vector<int> rankis;
  vector<int> rankie;
  
  string dirname; //出力ディレクトリ指定
  
  int procGrp;           ///< プロセスグループ番号
  int myRank;            ///< 自ノードのランク番号
  int numProc;           ///< 全ランク数
  std::string HostName;  ///< ホスト名
  
  string basename_g;
  string basename_f;
  int IS_DivideFunc;
  
public:
  /** コンストラクタ */
  LAYOUT();
  
  /**　デストラクタ */
  ~LAYOUT();
  
  /**
   * @brief CPMのポインタをコピーし、ランク情報を設定
   * @param [in] m_paraMngr  cpm_ParaManagerクラス
   * @return  エラーコード
   */
  bool importCPM(cpm_ParaManager* m_paraMngr)
  {
    if ( !m_paraMngr ) return false;
    paraMngr = m_paraMngr;
    setRankInfo();
    return true;
  }
  
  /**
   * @brief ランク情報をセットする
   * @param [in] m_paraMngr  CPMlibポインタ
   * @param [in] m_proGrp    プロセスグループ番号
   */
  void setRankInfo()
  {
    procGrp = 0;
    myRank  = paraMngr->GetMyRankID();
    numProc = paraMngr->GetNumRank();
    HostName= paraMngr->GetHostName();
  }
  
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
  void ReadInputFile(TextParser* tpCntl);
  
  /**
   * @brief layout grid ファイルの出力
   */
  void OutLayGrid(DfiInfo *D);
  
  /**
   * @brief layout func ファイルの出力
   */
  void OutLayFunc(DfiInfo *D);
  
  
  /**
   * @brief layoutファイル名を作成する
   * @param [in] prefix ファイル接頭文字
   * @param [in] m_step ステップ数
   */
  std::string Generate_LayoutFileName(const std::string prefix, const unsigned m_step);
  
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
