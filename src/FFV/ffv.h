#ifndef _FFV_H_
#define _FFV_H_

// #################################################################
//
// FFV : Frontflow / violet
//
// Copyright (c) All right reserved. 2012
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
//
// #################################################################
//
// 以下のマクロはcpm_Define.hで定義されている
//   REAL_TYPE
//   X_MINUS, Y_MINUS, Z_MINUS, X_PLUS, Y_PLUS, Z_PLUS
//   X_DIR, Y_DIR, Z_DIR
//   PLUS2MINUS, MINUS2PLUS, BOTH

/** 
 * @file ffv.h
 * @brief FFV Class Header
 * @author kero
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <iostream>
#include <fstream>

#include "cpm_ParaManager.h"
#include "cpm_TextParserDomain.h"

#include "../FB/FB_Define.h"
#include "ffv_Define.h"
#include "mydebug.h"
#include "FBUtility.h"

using namespace std;

class FFV {
protected:
  
  
  int session_maxStep;     ///< セッションのステップ数
  int session_currentStep; ///< セッションの現在のステップ
  
  FILE *mp;  ///< 標準出力
  
public:
  cpm_ParaManager *paraMngr; ///< Cartesian Partition Maneger
  
  
public:
  /** コンストラクタ */
  FFV();
  
  /**　デストラクタ */
  virtual ~FFV();
  
  
public:
  /** 毎ステップ後に行う処理 */
  virtual bool stepPost();
  
  
  /** シミュレーションの終了時の処理
   * プロファイルの統計処理ほか
   */
  virtual bool Post();
  
  
  /** 初期化 
   * 格子生成、ビットフラグ処理ほか
   * @param[in] argc  main関数の引数の個数
   * @param[in] argv  main関数の引数リスト
   */
  virtual int Initialize(int argc, char **argv);
  
  
  /** 1ステップのコアの処理
   * @param[in] m_step   現在のステップ数
   */
  virtual int Loop(int m_step);
  
  
  /** シミュレーションの1ステップの処理
   *  Loop() + stepPost()
   */
  virtual int MainLoop();
  
  
  /** コマンドラインヘルプ */
  virtual void Usage();
  
  
  /** 
   * @brief マスターノードのみ trueを返す
   * @ret true(Rank==0) / false(Rank!=0)
   */
  bool IsMaster() 
  {
    return ( !paraMngr->IsParallel() || (paraMngr->GetMyRankID() == 0) ) ? true : false;
  }
  
};

#endif // _FFV_H_