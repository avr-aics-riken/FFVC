//##################################################################################
//
// FFV-C : Frontflow / violet Cartesian
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2014 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

/**
 * @file   ffv_Filter.C
 * @brief  FFV filter
 * @author aics
 */

#include "ffv.h"

int FFV::FilterLoop()
{
  printf("\n\n\tData sampling....\n\n");
}



// #################################################################
/* @brief フィルタ処理の初期化
 * @param [in] argc  main関数の引数の個数
 * @param [in] argv  main関数の引数リスト
 */
int FFV::FilterInitialize(int argc, char **argv)
{
  double TotalMemory   = 0.0;  ///< 計算に必要なメモリ量（ローカル）
  double PrepMemory    = 0.0;  ///< 初期化に必要なメモリ量（ローカル）
  double G_TotalMemory = 0.0;  ///< 計算に必要なメモリ量（グローバル）
  double G_PrepMemory  = 0.0;  ///< 初期化に必要なメモリ量（グローバル）
  double tmp_memory    = 0.0;  ///< 計算に必要なメモリ量（グローバル）
  double flop_task     = 0.0;  ///< flops計算用
  
  
  
  // cpm_ParaManagerのポインタをセット
  C.importCPM(paraMngr);
  V.importCPM(paraMngr);
  B.importCPM(paraMngr);
  BC.importCPM(paraMngr);
  MO.importCPM(paraMngr);
  
  // 入力ファイルの指定
  std::string input_file = argv[2];
  
  
  // ffvのパラメータローダのインスタンス生成
  TextParser tp_ffv;
  
  
  // パラメータのロードと保持
  {
    int ierror=0;
    
    if ( (ierror = tp_ffv.read(input_file)) != TP_NO_ERROR )
    {
      Hostonly_ stamped_printf("\tError at reading '%s' file : %d\n", input_file.c_str(), ierror);
      Exit(0);
    }
  }
  
  
  // TextParserクラスのポインタを各クラスに渡す
  C.importTP(&tp_ffv);
  B.importTP(&tp_ffv);
  M.importTP(&tp_ffv);
  MO.importTP(&tp_ffv);
  
  
  // 固定パラメータ
  fixedParameters();
  
  
  
  
  // ------------------------------------
  FILE* fp = NULL;
  
  // condition fileのオープン
  Hostonly_
  {
    if ( !(fp=fopen("filter-condition.txt", "w")) )
    {
      stamped_printf("\tSorry, can't open 'filter-condition.txt' file. Write failed.\n");
      return -1;
    }
  }
  
  // メッセージ表示
  Hostonly_
  {
    FBUtility::printVersion(fp,     "FFV Filter ", FFVC_VERSION_NO);
    FBUtility::printVersion(stdout, "FFV Filter ", FFVC_VERSION_NO);
  }
  
  
  // 反復制御クラスのインスタンス
  //C.getIteration();
  
  
  // 流体の解法アルゴリズムを取得
  C.getSolvingMethod4Flow();
  
  
  // 線形ソルバーの特定
  identifyLinearSolver(&tp_ffv);
  
  
  // 計算モデルの入力ソース情報を取得
  C.getGeometryModel();
  
  
  // Intrinsic classの同定
  identifyExample(fp);
  
  
  
  // パラメータの取得と計算領域の初期化，並列モードを返す
  std::string str_para = setupDomain(&tp_ffv);
  
  
  // mat[], cmp[]の作成
  createTable(fp);
  
  
  // 媒質情報をパラメータファイルから読み込み，媒質リストを作成する
  setMediumList(fp);
  
  
  V.setControlVars(C.NoCompo, Ex);
  
  
  // CompoListの設定，外部境界条件の読み込み保持
  setBCinfo();
  
  
  
  // タイミング測定の初期化
  if ( C.Mode.Profiling != OFF )
  {
    ModeTiming = ON;
    TIMING__ PM.initialize( tm_END );
    TIMING__ PM.setRankInfo( paraMngr->GetMyRankID() );
    TIMING__ PM.setParallelMode(str_para, C.num_thread, C.num_process);
    set_timing_label();
  }
  
  
  // タイミング測定開始
  TIMING_start(tm_init_sct);
  
  
  
  // 前処理に用いるデータクラスのアロケート -----------------------------------------------------

  allocArray_Prep(PrepMemory, TotalMemory);
  
  
  
  TIMING_start(tm_voxel_prep_sct);
  
  
  
  // 各問題に応じてモデルを設定 >> Polylib + Cutlib
  // 外部境界面およびガイドセルのカットとIDの処理
  setModel(PrepMemory, TotalMemory, fp);
  
  
  
  // 定義点上に交点がある場合の処理 >> カットするポリゴンのエントリ番号でフィルする
  unsigned long fill_cut = V.modifyCutOnCellCenter(d_bid, d_cut, C.FillID);
  
  Hostonly_
  {
    printf(    "\tCut on cell center = %16ld\n\n", fill_cut);
    fprintf(fp,"\tCut on cell center = %16ld\n\n", fill_cut);
  }
  
  
  
  // 領域情報の表示
  Hostonly_
  {
    printf("\n---------------------------------------------------------------------------\n");
    printf("\n\t>> Global Domain Information\n\n");
    C.printGlobalDomain(stdout, G_size, G_origin, G_region, pitch);
    
    fprintf(fp,"\n---------------------------------------------------------------------------\n");
    fprintf(fp,"\n\t>> Global Domain Information\n\n");
    C.printGlobalDomain(fp, G_size, G_origin, G_region, pitch);
  }
  
  
  // メモリ消費量の情報を表示
  Hostonly_
  {
    printf(    "\n---------------------------------------------------------------------------\n\n");
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
  }
  G_PrepMemory = PrepMemory;
  
  displayMemoryInfo(fp, G_PrepMemory, PrepMemory, "Preprocessor");
  
  
  
  
  // サンプリング準備
  setMonitorList();
  
  
  
  // Fill
  Hostonly_
  {
    printf(    "\n---------------------------------------------------------------------------\n\n");
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
    printf(    "\t>> Fill\n\n");
    fprintf(fp,"\t>> Fill\n\n");
  }
  
  fill(fp);
  
  
  // 全周カットのあるセルを固体セルIDで埋める
  V.replaceIsolatedFcell(d_bcd, C.FillID, d_bid);
  
  
  
  // ∆tの決め方とKindOfSolverの組み合わせで無効なものをはねる
  if ( !DT.chkDtSelect() )
  {
    Hostonly_
    {
      printf(    "\tCombination of specified 'TimeIncrement' and 'KindOfSolver' is not permitted.\n");
      fprintf(fp,"\tCombination of specified 'TimeIncrement' and 'KindOfSolver' is not permitted.\n");
    }
    return -1;
  }
  
  
  // パラメータファイルから得られた内部BCコンポーネント数を表示
  Hostonly_
  {
    printf("\n---------------------------------------------------------------------------\n\n");
    printf("\t>> Components\n\n");
    C.printNoCompo(stdout);
    printf("\n"); fflush(stdout);
    
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
    fprintf(fp,"\t>> Components\n\n");
    C.printNoCompo(fp);
    fprintf(fp,"\n"); fflush(fp);
  }
  
  
  
  
  // HEX/FANコンポーネントの形状情報からBboxと体積率を計算
  if ( C.EnsCompo.fraction )
  {
    TIMING_start(tm_init_alloc);
    allocArray_CompoVF(PrepMemory, TotalMemory);
    TIMING_stop(tm_init_alloc);
    
    setComponentVF();
  }
  
  
  
  // 内部周期境界の場合のガイドセルのコピー処理
  V.adjMediumPrdcInner(d_bcd, cmp);
  
  
  // 媒質数とKindOfSolverの整合性をチェックする
  if ( !chkMediumConsistency() )
  {
    Hostonly_
    {
      stamped_printf(    "\tchkMediumConsistency()\n");
      stamped_fprintf(fp,"\tchkMediumConsistency()\n");
    }
    return -1;
  }
  
  
  
  
  // BCIndexにビット情報をエンコードとコンポーネントインデクスの再構築
  encodeBCindex(fp);
  
  
  
  // Polygonモニタの数をcmp[]にセット
  if ( C.SamplingMode )
  {
    MO.setMonitorNpoint(cmp, C.NoCompo);
  }
  
  
  
  // 体積力を使う場合のコンポーネント配列の確保
  TIMING_start(tm_init_alloc);
  allocArray_Forcing(PrepMemory, TotalMemory, fp);
  TIMING_stop(tm_init_alloc);
  
  
  // コンポーネントの体積率を8bitで量子化し，圧力損失コンポの場合にはFORCING_BITをON > bcdにエンコード
  V.setCmpFraction(cmp, d_bcd, d_cvf);
  
  
  
  
  // 周期境界条件が設定されている場合のBCIndexの周期条件の強制同期
  BC.setBCIperiodic(d_bcd);
  BC.setBCIperiodic(d_bcp);
  BC.setBCIperiodic(d_cdf);
  
  
  // bcd/bcp/cdfの同期
  if ( numProc > 1 )
  {
    if ( paraMngr->BndCommS3D(d_bcd, size[0], size[1], size[2], guide, 1) != CPM_SUCCESS ) Exit(0);
    if ( paraMngr->BndCommS3D(d_bcp, size[0], size[1], size[2], guide, 1) != CPM_SUCCESS ) Exit(0);
    if ( paraMngr->BndCommS3D(d_cdf, size[0], size[1], size[2], guide, 1) != CPM_SUCCESS ) Exit(0);
  }
  
  
  // サンプリング点の整合性をチェック
  if ( C.SamplingMode == ON ) MO.checkStatus();
  
  
  // 時間積分幅 deltaT や物理パラメータの設定
  setParameters();
  
  
  // 必要なパラメータをSetBC3Dクラスオブジェクトにコピーする >> setParameters()の後
  BC.setControlVars(&C, mat, &RF, Ex);
  
  
  
  
  // コンポーネントのグローバルインデクス情報を取得し，CompoListの内容とセル数の情報を表示する
  setGlobalCompoIdx_displayInfo(fp);
  
  
  
  
  // 外部境界面の開口率を計算する
  V.countOpenAreaOfDomain(d_bcd, C.OpenDomain);
  
  
  Hostonly_
  {
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n\n");
    printf(    "\n---------------------------------------------------------------------------\n\n\n");
  }
  
  if (C.FIO.IO_Voxel == ON)
  {
    Ex->writeSVX(d_bcd, &C);
    Hostonly_
    {
      fprintf(fp,"\tVoxel file 'example.svx' was written.\n");
      printf(    "\tVoxel file 'example.svx' was written.\n");
      fprintf(fp,"\n---------------------------------------------------------------------------\n\n\n");
      printf(    "\n---------------------------------------------------------------------------\n\n\n");
    }
  }
  
  // mid[]を解放する  ---------------------------
  if ( d_mid ) delete [] d_mid;
  
  
  
  // 各ノードの領域情報をファイル出力
  gatherDomainInfo();
  
  
  
  
  TIMING_stop(tm_voxel_prep_sct);
  // ここまでが準備の時間セクション
  
  
  
  // 計算に用いる配列のアロケート ----------------------------------------------------------------------------------
  allocate_Main(TotalMemory);
  
  
  
  // 初期値とリスタート処理 瞬時値と平均値に分けて処理　------------------
  Hostonly_
  {
    printf(    "\n---------------------------------------------------------------------------\n\n");
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
  }
  
  // リスタートモードの選択
  if ( C.Start != initial_start)
  {
    selectRestartMode();
  }
  
  
  // 瞬時値のリスタート
  TIMING_start(tm_restart);
  Restart(fp);
  TIMING_stop(tm_restart);
  
  
  // 制御インターバルの初期化
  initInterval();
  
  
  // 平均値のリスタート
  if ( C.Mode.AverageRestart == ON )
  {
    TIMING_start(tm_restart);
    RestartAvrerage(fp, flop_task);
    TIMING_stop(tm_restart);
  }
  
  
  // リスタートの最大値と最小値の表示
  if ( C.Start != initial_start )
  {
    RestartDisplayMinmax(fp, flop_task);
  }
  
  
  
  // 利用ライブラリのバージョン番号取得
  C.ver_CPM = cpm_Base::getVersionInfo();
  C.ver_CIO = cio_DFI::getVersionInfo();
  C.ver_Poly= PL->getVersionInfo();
  C.ver_PM  = PM.getVersionInfo();
  C.ver_CUT = cutlib_VersionInfo();
  C.ver_TP  = tp_ffv.getVersionInfo();
  
  
  // 制御パラメータ，物理パラメータの表示
  Hostonly_
  {
    displayParameters(fp);
  }
  
  
  
  // ドライバ条件のチェック
  BC.checkDriver(fp);
  
  
  
  // 初期条件の条件設定
  setInitialCondition();
  
  
  // サンプリング元となるデータ配列の登録
  if ( C.SamplingMode == ON )
  {
    MO.setDataPtrs(d_v, d_p, d_ie, d_vrt);
  }
  
  
  // 出力ファイルの初期化
  initFileOut();
  
  
  
  // メモリ使用量の表示
  Hostonly_
  {
    printf(    "\n---------------------------------------------------------------------------\n\n");
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
  }
  
  G_TotalMemory = TotalMemory;
  
  displayMemoryInfo(fp, G_TotalMemory, TotalMemory, "Solver");
  
  
  
  // 履歴出力準備
  //prepHistoryOutput();
  
  
  Hostonly_ if ( fp ) fclose(fp);
  
  
  TIMING_stop(tm_init_sct);
  
  return 1;
}
