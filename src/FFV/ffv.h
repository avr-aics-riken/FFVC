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

#include "FB_Define.h"
#include "ffv_Define.h"
#include "mydebug.h"
#include "FBUtility.h"
#include "Control.h"
#include "Alloc.h"
#include "FileIO.h"
#include "ParseBC.h"
#include "ParseMat.h"
//#include "ffv_SetBC.h"
//#include "VoxInfo.h"

#include "TPControl.h"

#include "omp.h"


#include "IP_Duct.h"
#include "IP_PPLT2D.h"
#include "IP_SHC1D.h"
#include "IP_PMT.h"
#include "IP_Rect.h"
#include "IP_Step.h"
#include "IP_Cylinder.h"
#include "IP_Polygon.h"
#include "IP_Sphere.h"

using namespace std;

class FFV {
private:
  
  int session_maxStep;     ///< セッションのステップ数
  int session_currentStep; ///< セッションの現在のステップ
  
  int G_size[3];           ///< 全ドメインの分割数
  REAL_TYPE G_org[3];      ///< 全ドメインの基点
  REAL_TYPE G_reg[3];      ///< 全ドメインのサイズ
  
  FILE *mp;    ///< 標準出力
  FILE *fp_b;  ///< 基本情報
  FILE *fp_w;  ///< 壁面情報
  FILE *fp_c;  ///< コンポーネント情報
  FILE *fp_d;  ///< 流量収支情報
  FILE *fp_i;  ///< 反復履歴情報
  FILE *fp_f;  ///< 力の履歴情報
  
  Control C;                 ///< 制御パラメータクラス
  FileIO F;                  ///< ファイル入出力クラス
  DTcntl DT;                 ///< 時間制御クラス
  Intrinsic* Ex;             ///< pointer to a base class
  ItrCtl IC[ItrCtl::ic_END]; ///< 反復情報管理クラス
  ReferenceFrame RF;         ///< 参照座標系クラス
  
//  SetBC3D BC;                ///< BCクラス
  
public:
  cpm_ParaManager *paraMngr; ///< Cartesian Partition Maneger
  
public:
  /** コンストラクタ */
  FFV();
  
  /**　デストラクタ */
  ~FFV();
  
  
public:
  
  /**
   @brief 組み込み例題のインスタンス
   @param [in] Cref Controlクラスのポインタ
   */
  void connectExample(Control* Cref);
  
  
  /** 計算領域情報を設定する
   * @param [in] dom_file  ドメインファイル名
   */
  void DomainInitialize(const string dom_file);
  
  
  /**
   * @brief 固定パラメータの設定
   */
  void fixed_parameters();
  
  /**
   * @brief 組み込み例題の設定
   * @param [in] Cref    コントロールクラス
   * @param [in] tpCntl  テキストパーサーのラッパー
   */
  void getExample(Control* Cref, TPControl* tpCntl);
  
  /** 初期化 
   * 格子生成、ビットフラグ処理ほか
   * @param[in] argc  main関数の引数の個数
   * @param[in] argv  main関数の引数リスト
   */
  int Initialize(int argc, char **argv);
  
  
  /** 
   * @brief マスターノードのみ trueを返す
   * @ret true(Rank==0) / false(Rank!=0)
   */
  bool IsMaster() 
  {
    return ( paraMngr->GetMyRankID() == 0 ) ? true : false;
  }
  
  
  /** 1ステップのコアの処理
   * @param[in] m_step   現在のステップ数
   */
  int Loop(int m_step);
  
  
  /** シミュレーションの1ステップの処理
   *  Loop() + stepPost()
   */
  int MainLoop();
  
  
  /** シミュレーションの終了時の処理
   * プロファイルの統計処理ほか
   */
  bool Post();
  
  
  /**
   * @brief 読み込んだ領域情報のデバッグライト
   * @param [in] dInfo  領域情報クラス
   */
  void printDomainInfo( cpm_GlobalDomainInfo* dInfo );
  
  
  /**
   * @brief 並列化と分割の方法を保持
   */
  void setParallelism();
  
  
  /** 毎ステップ後に行う処理 */
  bool stepPost();
  
  
  /** コマンドラインヘルプ */
  virtual void Usage();
  
  

};

#endif // _FFV_H_