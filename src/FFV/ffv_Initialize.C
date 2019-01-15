//##################################################################################
//
// FFV-C : Frontflow / violet Cartesian
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2015 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2016 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

/**
 * @file   ffv_Initialize.C
 * @brief  FFV Class
 * @author aics
 */

#include "ffv.h"
#include "limits.h"
#include <algorithm>


// #################################################################
/* @brief 初期化格子生成、ビットフラグ処理ほか
 * @param [in] argc  main関数の引数の個数
 * @param [in] argv  main関数の引数リスト
 */
int FFV::Initialize(int argc, char **argv)
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
  GM.importCPM(paraMngr);

  for (int i=0; i<ic_END; i++)
  {
    LS[i].importCPM(paraMngr);
  }


  // 入力ファイルの指定
  std::string input_file = argv[1];


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

  //Hostonly_ printf("tp object= %d\n", sizeof(tp_ffv));


  // TextParserクラスのポインタを各クラスに渡す
  C.importTP(&tp_ffv);
  B.importTP(&tp_ffv);
  M.importTP(&tp_ffv);
  MO.importTP(&tp_ffv);


  // 固定パラメータ
  fixedParameters();



  // File IO classのインスタンス
  identifyFIO(&tp_ffv);


  F->importCPM(paraMngr);
  F->importExtClass(&tp_ffv, &RF, &C);




  // ------------------------------------
  FILE* fp = NULL;

  // condition fileのオープン
  Hostonly_
  {
    if ( !(fp=fopen("condition.txt", "w")) )
    {
      stamped_printf("\tSorry, can't open 'condition.txt' file. Write failed.\n");
      return -1;
    }
  }

  // メッセージ表示
  Hostonly_
  {
    FBUtility::printVersion(fp,     "Welcome to FFV  ", FFVC_VERSION_NO);
    FBUtility::printVersion(stdout, "Welcome to FFV  ", FFVC_VERSION_NO);

    FBUtility::printVersion(fp,     "FlowBase        ", FB_VERS);
    FBUtility::printVersion(stdout, "FlowBase        ", FB_VERS);
  }


  // 反復制御クラスのインスタンス
  C.getIteration();


  // 流体の解法アルゴリズムを取得
  C.getSolvingMethod4Flow();


  // 線形ソルバーの特定
  identifyLinearSolver(&tp_ffv);


  // 計算モデルの入力ソース情報を取得
  C.getGeometryModel();


  // Intrinsic classの同定
  identifyExample(fp);


  // パラメータの取得と計算領域の初期化，並列モードを返す
  // Polylibの基準値も設定
  std::string str_para = setDomain(&tp_ffv);

  // mat[], cmp[]の作成
  createTable(fp);


  // 媒質情報をパラメータファイルから読み込み，媒質リストを作成する
  setMediumList(fp);


  // フィルパラメータ
  GM.getFillParam(&tp_ffv, C.Unit.Param, C.RefLength, C.NoMedium, mat, fp);


  V.setControlVars(Ex);


  // CompoListの設定、境界条件の読み込み保持、パラメータの無次元化
  setBCinfo();


  // コンポーネントポインタのコピー
  GM.setCompoPtr(C.NoCompo, cmp);


  // タイミング測定の初期化
  if ( C.Mode.Profiling != OFF )
  {
    ModeTiming = ON;
    TIMING__ PM.initialize( PM_NUM_MAX );
    TIMING__ PM.setRankInfo( paraMngr->GetMyRankID(procGrp) );
    TIMING__ PM.setParallelMode(str_para, C.num_thread, C.num_process);
    set_timing_label();
  }


  // タイミング測定開始
  TIMING_start("Initialization_Section");


  // 前処理に用いるデータクラスのアロケート -----------------------------------------------------
  TIMING_start("Allocate_Arrays");


  // 配列アロケート前に一度コール
  setArraySize();
  allocArray_Prep(PrepMemory, TotalMemory);

  // カット情報保持領域
  allocArray_Cut(PrepMemory, TotalMemory);
  initCutInfo();

  TIMING_stop("Allocate_Arrays");



  TIMING_start("Voxel_Prep_Section");


  // 各問題に応じてモデルを設定 >> Polylib
  // 外部境界面およびガイドセルのカットとIDの処理
  setModel(PrepMemory, TotalMemory, fp);


  // 回転体
  setComponentSR();



  // 領域情報の表示
  Hostonly_
  {
    printf("\n----------\n");
    printf("\n\t>> Global Domain Information\n\n");
    printGlobalDomain(stdout);

    fprintf(fp,"\n----------\n");
    fprintf(fp,"\n\t>> Global Domain Information\n\n");
    printGlobalDomain(fp);
  }


  // メモリ消費量の情報を表示
  Hostonly_
  {
    printf(    "\n----------\n\n");
    fprintf(fp,"\n----------\n\n");
  }
  G_PrepMemory = PrepMemory;

  displayMemoryInfo(fp, G_PrepMemory, PrepMemory, "Preprocessor");
  displayMemoryInfo(stdout, G_PrepMemory, PrepMemory, "Preprocessor");



  // サンプリング準備
  setMonitorList();


  /* 再考
  if ( C.Mode.Example == id_Polygon )
  {

    // ガイドセルのIDをd_midに転写
    for (int face=0; face<NOFACE; face++)
    {
      if( nID[face] >= 0 ) continue;

      GM.copyIDonGuide(face, d_bcd, d_mid);
    }


    // Fill by Seed >> d_mid={-1, SeedID}
    Hostonly_
    {
      printf(    "\n----------\n\n");
      fprintf(fp,"\n----------\n\n");
      printf(    "\t>> Seed Filling Process\n\n");
      fprintf(fp,"\t>> Seed Filling Process\n\n");
    }

    TIMING_start("SeedFilling");
    GM.SeedFilling(fp, cmp, mat, d_mid, PL, PG, C.NoCompo);
    TIMING_stop("SeedFilling");

    //F->writeSVX(d_mid);

    // Sub-sampling
    Hostonly_
    {
      printf(    "\n----------\n\n");
      fprintf(fp,"\n----------\n\n");
      printf(    "\t>> Sub-sampling Process\n\n");
      fprintf(fp,"\t>> Sub-sampling Process\n\n");
    }

    TIMING_start("SubSampling");
    GM.SubSampling(fp, mat, d_mid, d_pvf, PL, C.NoCompo);
    TIMING_stop("SubSampling");
  }
  */


  // Fill
  Hostonly_
  {
    fprintf(fp,"\n----------\n\n");
    fprintf(fp,"\t>> Fill\n\n");
  }

  TIMING_start("Fill");
  if ( !GM.fill(d_bcd, d_bid, d_mid, d_cut) )
  {
    // debug routine
    //F->writeSVX(d_bcd);
    //F->writeRawSPH(d_mid);
    F->writeSVX(d_mid);
    Exit(0);
  }
  TIMING_stop("Fill");

  F->writeSVX(d_bcd);

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
    fprintf(fp,"\n----------\n\n");
    fprintf(fp,"\t>> Components\n\n");
    C.printNoCompo(fp);
    fprintf(fp,"\n"); fflush(fp);
  }




  // 体積率コンポーネントの配列確保
  if ( C.EnsCompo.fraction )
  {
    TIMING_start("Allocate_Arrays");
    allocArray_CompoVF(PrepMemory, TotalMemory);
    TIMING_stop("Allocate_Arrays");
  }


  // 形状情報からBboxと体積率を計算
  setComponentVF();



  // 内部周期境界の場合のガイドセルのコピー処理
  V.adjMediumPrdcInner(d_bcd, cmp, C.NoCompo);


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
  TIMING_start("Encode_BCindex");
  encodeBCindex(fp);
  TIMING_stop("Encode_BCindex");

  V.chkFlag(d_bcp, d_bid, d_cut, d_bcd);


  // Polygonモニタの数をcmp[]にセット
  if ( C.SamplingMode )
  {
    MO.setMonitorNpoint(cmp, C.NoCompo);
  }



  // 体積力を使う場合のコンポーネント配列の確保
  TIMING_start("Allocate_Arrays");
  allocArray_Forcing(PrepMemory, TotalMemory, fp, &C, cmp);
  TIMING_stop("Allocate_Arrays");


  // コンポーネントの体積率を8bitで量子化し，圧力損失コンポの場合にはFORCING_BITをON > bcdにエンコード
  TIMING_start("Compo_Fraction");
  V.setCmpFraction(cmp, d_bcd, d_cvf, C.NoCompo);
  TIMING_stop("Compo_Fraction");




  // 周期境界条件が設定されている場合のBCIndexの周期条件の強制同期
  BC.setBCIperiodic(d_bcd, ensPeriodic);
  BC.setBCIperiodic(d_bcp, ensPeriodic);
  BC.setBCIperiodic(d_cdf, ensPeriodic);


  // bcd/bcp/cdfの同期
  if ( numProc > 1 )
  {
    if ( paraMngr->BndCommS3D(d_bcd, size[0], size[1], size[2], guide, 1, procGrp) != CPM_SUCCESS ) Exit(0);
    if ( paraMngr->BndCommS3D(d_bcp, size[0], size[1], size[2], guide, 1, procGrp) != CPM_SUCCESS ) Exit(0);
    if ( paraMngr->BndCommS3D(d_cdf, size[0], size[1], size[2], guide, 1, procGrp) != CPM_SUCCESS ) Exit(0);
  }


  // サンプリング点の整合性をチェック
  if ( C.SamplingMode == ON ) MO.checkStatus();


  // 時間積分幅 deltaT や物理パラメータの設定
  setParameters();


  // 必要なパラメータをSetBC3Dクラスオブジェクトにコピーする >> setParameters()の後
  BC.setControlVars(&C, mat, &RF, Ex);


  // set phase
  if ( C.BasicEqs == INCMP_2PHASE )
  {
    B.get_Phase(cmp);
  }


  // コンポーネントのグローバルインデクス情報を取得し，CompoListの内容とセル数の情報を表示する
  dispGlobalCompoInfo(fp);




  // 外部境界面の開口率を計算する
  V.countOpenAreaOfDomain(d_bcd, C.OpenDomain);


  Hostonly_
  {
    fprintf(fp,"\n----------\n\n\n");
  }



  // 各ノードの領域情報をファイル出力
  TIMING_start("Gather_DomainInfo");
  gatherDomainInfo();
  TIMING_stop("Gather_DomainInfo");


  // 交点のグリフを出力
  if ( C.Hide.GlyphOutput != OFF )
  {
    TIMING_start("Generate_Glyph");
    /*
    int st[3], ed[3];
    st[0] = 145;
    ed[0] = 147;
    st[1] = 520;
    ed[1] = 545;
    st[2] = 230;
    ed[2] = 253;
    generateGlyph(d_cut, d_bid, fp, st, ed);
     */
    generateGlyph(d_cut, d_bid, fp);
    TIMING_stop("Generate_Glyph");
  }


  TIMING_stop("Voxel_Prep_Section");
  // ここまでが準備の時間セクション



  // IP_Cylinderの場合、物体の個数をチェック
  if ( C.Mode.Example == id_Cylinder )
  {
    if ( num_obstacle != ((IP_Cylinder*)Ex)->get_num_cyls() )
    {
      Hostonly_ printf("\tNumber of cylinders does not agree with the number of local boundaries %d %d\n",
                       num_obstacle, ((IP_Cylinder*)Ex)->get_num_cyls());
      Exit(0);
    }
  }



// ##########
#if 0
  write_distance(cut);
#endif
// ##########



  // 計算に用いる配列のアロケート ----------------------------------------------------------------------------------
  // SamplingでHelicityあるいはVorticityが指定されている場合には，Vorticityの配列を使う
  if ( MO.getStateVorticity() ) C.varState[var_Vorticity] = ON;

  TIMING_start("Allocate_Arrays");
  allocate_Main(TotalMemory);
  TIMING_stop("Allocate_Arrays");



  // File IO class への配列ポインタ
  F->setVarPointers(d_p,
                    d_v,
                    d_vf,
                    d_ie,
                    d_ws,
                    d_wv,
                    d_ap,
                    d_av,
                    d_ae,
                    d_dv,
                    d_rms_v,
                    d_rms_mean_v,
                    d_rms_p,
                    d_rms_mean_p,
                    d_rms_t,
                    d_rms_mean_t,
                    d_bcd,
                    d_cdf,
                    mat_tbl,
                    d_mid,
                    d_io_buffer,
                    d_av_mean,
                    d_arms_mean,
                    d_aR_mean,
                    d_aP_mean,
                    d_aE_mean,
                    d_aT_mean,
                    d_aPI_mean);

  F->getStartCondition();

  F->getStagingOption();

  F->getRestartDFI();



  // 初期値とリスタート処理 瞬時値と統計値に分けて処理　------------------
  Hostonly_
  {
    printf(    "\n----------\n\n");
    fprintf(fp,"\n----------\n\n");
  }


  // リスタートモードの選択と瞬時値のリスタート
  TIMING_start("Restart_Process");
  F->Restart(fp, CurrentStep, CurrentTime);
  TIMING_stop("Restart_Process");


  // 制御インターバルの初期化
  initInterval();


  // 統計値のリスタート
  if ( C.Mode.StatisticRestart == ON )
  {
    TIMING_start("Restart_Process");
    F->RestartStatistic(fp, CurrentStep, CurrentTime, CurrentStepStat, CurrentTimeStat, flop_task);
    TIMING_stop("Restart_Process");
  }


  // リスタートの最大値と最小値の表示
  if ( C.Start != initial_start )
  {
    F->RestartDisplayMinmax(fp, flop_task);
  }


  /* コンポーネントの統計値のリスタート
  if ( C.Mode.StatisticRestart == ON && num_obstacle > 0 )
  {
    if (C.Mode.Statistic == ON)
    {
      Hostonly_
      {
        if ( !(fp=fopen("component_statistic.txt", "r")) )
        {
          stamped_printf("\tSorry, can't open 'component_statistic.txt' file.\n");
          Exit(0);
        }

        H->loadCompoStatistics(fp, cmp, &C, cmp_force_avr);

        H->printCompoStatistics(cmp, cmp_force_avr);
      }
    }
  }
  */



  // 利用ライブラリのバージョン番号取得
  C.ver_CPM = cpm_Base::getVersionInfo();
  C.ver_CDM = cdm_DFI::getVersionInfo();
  C.ver_Poly= PL->getVersionInfo();
  C.ver_PM  = PM.getVersionInfo();
  C.ver_TP  = tp_ffv.getVersionInfo();



  // ドライバ条件のチェック
  BC.checkDriver(fp);



  // 初期条件の条件設定
  setInitialCondition();


  // サンプリング元となるデータ配列の登録
  if ( C.SamplingMode == ON )
  {
    MO.setDataPtrs(d_v, d_p, d_ie, d_vrt);
  }



  // CellIDとBCflagの出力 (guide cell=0)
  // 戻り値：全セルが同じ値cの場合にはcの値が戻り、異なるセルの値が存在する場合には-1
  int id_cell = F->writeCellID(0);
  int id_bcf  = F->writeBCflag(0);

  // 出力ファイルの初期化
  F->initFileOut(id_cell, id_bcf);


  // debug: obsolete format
  //F->writeSVX(d_bcd);


  // IBLANK 出力後に mid[]を解放する  ---------------------------
  if ( d_mid ) delete [] d_mid;





  // セッションを開始したときに、初期値をファイル出力  性能測定モードのときには出力しない
  if ( (C.Hide.PM_Test == OFF) && (0 == CurrentStep) )
  {
    flop_task = 0.0;
    F->OutputBasicVariables(CurrentStep, CurrentTime, flop_task);

    if ( (C.Mode.Statistic == ON) && (C.Start != initial_start) )
    {
      double flop_count=0.0;
      F->OutputStatisticalVarables(CurrentStep, CurrentTime, CurrentStepStat, CurrentTimeStat, flop_count);
    }
  }


  // 粗い格子を用いたリスタート時には出力
  if ( (C.Start == restart_sameDiv_refinement) || (C.Start == restart_diffDiv_refinement) )
  {
    flop_task = 0.0;
    F->OutputBasicVariables(CurrentStep, CurrentTime, flop_task);
  }


  // セルフェイス速度から領域境界平均速度を求める
  for (int face=0; face<NOFACE; face++)
  {
    REAL_TYPE vsum = 0.0;
    BoundaryOuter* T = BC.exportOBC(face);

    vobc_face_massflow_(&vsum, size, &guide, &face, d_vf, d_cdf, nID);
    T->setDomainMF(vsum);
  }


  // 反復法クラスの初期化
  LS_initialize(TotalMemory, &tp_ffv);


  // 制御パラメータ，物理パラメータの表示
  Hostonly_
  {
    displayParameters(fp);
  }


  // メモリ使用量の表示
  Hostonly_
  {
    printf(    "\n----------\n\n");
    fprintf(fp,"\n----------\n\n");
  }

  G_TotalMemory = TotalMemory;

  displayMemoryInfo(fp, G_TotalMemory, TotalMemory, "Solver");
  displayMemoryInfo(stdout, G_TotalMemory, TotalMemory, "Solver");


  // 粒子追跡 ------------------
  Hostonly_ {
    fprintf(fp,"\n----------\n\n");
    fprintf(stdout,"\n----------\n\n");
  }
  fflush(stdout);

  if (C.Mode.ParticleTracking == ON) {
    TR = new Cloud(d_bcd,
                   d_v,
                   deltaT,
                   &tp_ffv,
                   &PM);
    Hostonly_ {
      fprintf(fp,"\n\tParticle Tracking  ON\n\n");
      fprintf(stdout,"\n\tParticle Tracking  ON\n\n");
    }
      
    TR->importCPM(paraMngr);
    TR->setRankInfo(paraMngr, procGrp);
    TR->setDomainInfo(C.guide, C.RefLength);
    TR->initCloud(fp);

  }
  else
  {
    Hostonly_ {
      fprintf(fp,"\n\tParticle Tracking  OFF\n");
      fprintf(stdout,"\n\tParticle Tracking  OFF\n");
    }
  }


  TIMING_stop("Initialization_Section");

  // チェックモードの場合のコメント表示，前処理のみで中止---------------------------------------------------------
  if ( C.CheckParam == ON)
  {
		Hostonly_
    {
      printf(     "\n\tCheck mode --- Only pre-process\n\n");
      fprintf(fp, "\n\tCheck mode --- Only pre-process\n\n");
    }
    return 0;
	}

  if ( C.Hide.DryRun == ON )
  {
    Hostonly_
    {
      printf(     "\n\t#############  DRY RUN for DEBUGGING BC  #############\n\n");
    }
  }

  fflush(stdout);
  Hostonly_ fflush(fp);

  // 履歴出力準備
  prepHistoryOutput();


  Hostonly_ if ( fp ) fclose(fp);

  return 1;
}


// #################################################################
/* @brief 主計算部分に用いる配列のアロケーション
 * @param [in,out] total ソルバーに使用するメモリ量
 */
void FFV::allocate_Main(double &total)
{
  // 基本変数と必須領域
  allocArray_Main(total, &C);


  if ( C.LES.Calc == ON )
  {
    allocArray_LES(total);
  }


  if ( (C.AlgorithmF == Flow_FS_AB2) || (C.AlgorithmF == Flow_FS_AB_CN) )
  {
    allocArray_AB2(total);
  }

  if ( C.BasicEqs == INCMP_2PHASE )
  {
    allocArray_Interface(total);
  }


  // 時間平均・統計処理用の配列
  if ( C.Mode.Statistic == ON )
  {
    allocArray_Average(total, &C);

    allocArray_Statistic(total, &C);
  }
}



// #################################################################
/* @brief 全Voxelモデルの媒質数とKOSの整合性をチェック
 * @retval エラーコード
 */
bool FFV::chkMediumConsistency()
{
  int nmSolid = C.NoMediumSolid;
  int nmFluid = C.NoMediumFluid;

  if ( numProc > 1 )
  {
    int nms = nmSolid;
    int nmf = nmFluid;
    paraMngr->Allreduce(&nms, &nmSolid, 1, MPI_SUM, procGrp);
    paraMngr->Allreduce(&nmf, &nmFluid, 1, MPI_SUM, procGrp);
  }

  if ( (nmFluid == 0) && (nmSolid == 0) )
  {
    Hostonly_ printf("\tError : No medium\n");
    return false;
  }

  switch (C.KindOfSolver)
  {
    case FLOW_ONLY:
    case THERMAL_FLOW:
    case THERMAL_FLOW_NATURAL:

      if ( nmFluid == 0 )
      {
        Hostonly_ printf("\tError : No FLUID medium\n");
        return false;
      }
      break;

    case CONJUGATE_HT:
    case CONJUGATE_HT_NATURAL:
      if ( ( nmFluid == 0 ) || ( nmSolid == 0 ) )
      {
        Hostonly_ printf("\tError : Fluid/Solid should have at least one medium.\n");
        return false;
      }
      break;

    case SOLID_CONDUCTION:
      if ( nmSolid == 0 )
      {
        Hostonly_ printf("\tError : No Solid medium\n");
        return false;
      }
      break;
  };

  return true;
}



// #################################################################
/* @brief mat[], cmp[]のテーブルを作成
 * @param [in] fp  file pointer
 */
void FFV::createTable(FILE* fp)
{
  // コンポーネント数，境界条件数，媒質数を取得し，配列をアロケートする
  C.getNoOfComponent();


  // ポリゴンモデルの場合の局所境界条件数のチェック
  if ( C.Mode.Example == id_Polygon  &&  C.NoBC == 0 )
  {
    Hostonly_{
      printf("Error : In case of using polygon model, at least one local boundary condition is required.\n");
    }
    Exit(0);
  }


  // 媒質リストをインスタンス [0]はダミーとして利用しないので，配列の大きさはプラス１する
  if ( !(mat = new MediumList[C.NoCompo+1]) ) Exit(0);


  // CompoListクラスをインスタンス [0]はダミーとして利用しないので，配列の大きさはプラス１する
  if ( !(cmp = new CompoList[C.NoCompo+1]) ) Exit(0);


  // OBSTACLEコンポーネントの積算用 cmp_force_global[C.NoCompo+1][3]のイメージ
  if ( !(cmp_force_global = new REAL_TYPE[3*(C.NoCompo+1)]) ) Exit(0);
  if ( !(cmp_force_avr    = new REAL_TYPE[3*(C.NoCompo+1)]) ) Exit(0);
  if ( !(cmp_force_local  = new REAL_TYPE[3*(C.NoCompo+1)]) ) Exit(0);

  for (int i=0; i<3*(C.NoCompo+1); i++)
  {
    cmp_force_global[i] = 0.0;
    cmp_force_avr[i]    = 0.0;
    cmp_force_local[i]  = 0.0;
  }


  // 積算用バッファ
  if ( !(buffer_force = new REAL_TYPE[3*numProc]) ) Exit(0);
  for (int i=0; i<3*numProc; i++) buffer_force[i] = 0.0;


  // 各コンポーネントのOBSTACLEの有無
  if ( !(global_obstacle = new int[C.NoCompo+1]) ) Exit(0);
  for (int i=0; i<C.NoCompo+1; i++) global_obstacle[i] = 0;


  // CompoList, MediumListのポインタをセット
  BC.importCMP_MAT(cmp, mat);


  // Fortran用のデータ保持配列 >> mat_tbl[C.NoCompo+1][3]のイメージ
  if ( !(mat_tbl = new double[3*(C.NoCompo+1)]) ) Exit(0);
  for (int i=0; i<3*(C.NoCompo+1); i++) mat_tbl[i] = 1.0; // ゼロ割防止のため，1.0をいれておく

  // vec_tbl[C.NoCompo+1][7]のイメージ
  if ( !(vec_tbl = new REAL_TYPE[7*(C.NoCompo+1)]) ) Exit(0);
  for (int i=0; i<7*(C.NoCompo+1); i++) vec_tbl[i] = 0.0;


  Hostonly_
  {
    fprintf(fp,"\n----------\n\n");
    fprintf(fp,"\n\t>> Tables\n\n");
    fprintf(fp,"\tNo. of Medium    = %3d\n", C.NoMedium);
    fprintf(fp,"\tNo. of LocalBC   = %3d\n", C.NoBC);
    fprintf(fp,"\t----------------------\n");
    fprintf(fp,"\tNo. of Component = %3d\n", C.NoCompo);
    fprintf(fp,"\n");
  }

}



// #################################################################
/* @brief コンポーネントのグローバルなBboxを求め，CompoListの情報を表示
 * @param [in] fp file pointer
 */
void FFV::dispGlobalCompoInfo(FILE* fp)
{
  int st_i, st_j, st_k, ed_i, ed_j, ed_k;
  int node_st_i, node_st_j, node_st_k;
  int st[3], ed[3];


  // グローバルインデクスの配列インスタンス
  int* cgb = new int[6*(C.NoCompo+1)];

  for (int i=0; i<6*(C.NoCompo+1); i++) cgb[i] = 0;



  // ローカルインデクスからグローバルインデクスに変換
  for (int m=1; m<=C.NoCompo; m++)
  {

    if ( !cmp[m].existLocal() ) // コンポーネントが存在しないノードはゼロを代入
    {
      cgb[6*m+0] = 0;
      cgb[6*m+1] = 0;
      cgb[6*m+2] = 0;
      cgb[6*m+3] = 0;
      cgb[6*m+4] = 0;
      cgb[6*m+5] = 0;
    }
    else // コンポーネントが存在する場合
    {
      cmp[m].getBbox(st, ed);

      st_i = st[0];
      st_j = st[1];
      st_k = st[2];
      ed_i = ed[0];
      ed_j = ed[1];
      ed_k = ed[2];

      if ( numProc > 1 )
      {
        node_st_i = head[0];
        node_st_j = head[1];
        node_st_k = head[2];

        cgb[6*m+0] = node_st_i + st_i - 1;
        cgb[6*m+1] = node_st_j + st_j - 1;
        cgb[6*m+2] = node_st_k + st_k - 1;
        cgb[6*m+3] = node_st_i + ed_i - 1;
        cgb[6*m+4] = node_st_j + ed_j - 1;
        cgb[6*m+5] = node_st_k + ed_k - 1;
      }
      else
      {
        cgb[6*m+0] = st_i;
        cgb[6*m+1] = st_j;
        cgb[6*m+2] = st_k;
        cgb[6*m+3] = ed_i;
        cgb[6*m+4] = ed_j;
        cgb[6*m+5] = ed_k;
      }

    }
  }

  // 領域全体のbboxを求める
  if ( numProc > 1 )
  {
    int* m_gArray = NULL;
    int* m_eArray = NULL;
    int array_size = 6*(C.NoCompo+1);
    int st_x, st_y, st_z, ed_x, ed_y, ed_z, es;

    if ( !(m_gArray = new int[numProc*6]) ) Exit(0);
    if ( !(m_eArray = new int[numProc]  ) ) Exit(0);

    for (int i=0; i<numProc*6; i++) m_gArray[i] = 0;
    for (int i=0; i<numProc; i++) m_eArray[i] = 0;

    for (int n=1; n<=C.NoCompo; n++)
    {
      if ( numProc > 1 )
      {
        es = ( cmp[n].existLocal() ) ? 1 : 0;
        if ( paraMngr->Gather(&es, 1, m_eArray, 1, 0, procGrp) != CPM_SUCCESS ) Exit(0);
        if ( paraMngr->Gather(&cgb[6*n], 6, m_gArray, 6, 0, procGrp) != CPM_SUCCESS ) Exit(0);


        if (myRank == 0) // マスターノードのみ
        {
          // 初期値
          cgb[6*n+0] = 100000000;
          cgb[6*n+1] = 100000000;
          cgb[6*n+2] = 100000000;
          cgb[6*n+3] = 0;
          cgb[6*n+4] = 0;
          cgb[6*n+5] = 0;

          for (int m=0; m<numProc; m++)
          {
            if ( m_eArray[m]==1 ) // コンポーネントの存在ランクのみを対象とする
            {
              st_x = m_gArray[6*m+0]; // 各ランクのコンポーネント存在範囲のインデクス
              st_y = m_gArray[6*m+1];
              st_z = m_gArray[6*m+2];
              ed_x = m_gArray[6*m+3];
              ed_y = m_gArray[6*m+4];
              ed_z = m_gArray[6*m+5];

              if( st_x < cgb[6*n+0] ) { cgb[6*n+0] = st_x; } // 最大値と最小値を求める
              if( st_y < cgb[6*n+1] ) { cgb[6*n+1] = st_y; }
              if( st_z < cgb[6*n+2] ) { cgb[6*n+2] = st_z; }
              if( ed_x > cgb[6*n+3] ) { cgb[6*n+3] = ed_x; }
              if( ed_y > cgb[6*n+4] ) { cgb[6*n+4] = ed_y; }
              if( ed_z > cgb[6*n+5] ) { cgb[6*n+5] = ed_z; }
            }
          }
        }
      }

    }

    // destroy
    if (m_gArray)
    {
      delete[] m_gArray;
      m_gArray = NULL;
    }

    if (m_eArray)
    {
      delete[] m_eArray;
      m_eArray = NULL;
    }

  }


  // CompoListの内容とセル数の情報を表示する
  Hostonly_
  {
    fprintf(fp,"\n----------\n\n");
    fprintf(fp,"\t>> Component Information\n\n");
  }

  displayCompoInfo(cgb, fp);


  if ( cgb ) { delete [] cgb; cgb=NULL; }

}




// #################################################################
/* @brief CompoListの情報を表示する
 * @param [in]  fp   ファイルポインタ
 */
void FFV::displayCompoInfo(const int* cgb, FILE* fp)
{
  // サマリー
  Hostonly_
  {
    B.printCompoSummary(fp, cmp, C.BasicEqs);
  }

  double cr = (double)G_Wcell/ ( (double)G_size[0] * (double)G_size[1] * (double)G_size[2]) *100.0;

  // セル数の情報を表示する
  Hostonly_
  {
    fprintf(fp, "\tThis model includes %4d solid %s  [Solid cell ratio inside computational domain : %9.5f %%]\n\n",
            C.NoMediumSolid, (C.NoMediumSolid>1) ? "IDs" : "ID", cr);
  }



  // 詳細情報
  if ( C.NoBC >0 )
  {
    Hostonly_
    {
      B.printCompo( fp,     cgb, mat, cmp, BC.exportOBC() );
    }
  }


  // Check consistency of boundary condition
  for (int n=1; n<=C.NoCompo; n++)
  {
    if ( cmp[n].getType() == HT_SN )
    {
      if ( (C.KindOfSolver == FLOW_ONLY) || (C.KindOfSolver == THERMAL_FLOW) || (C.KindOfSolver == SOLID_CONDUCTION) )
      {
        Hostonly_ printf("\tInconsistent parameters of combination between Kind of Solver and Heat Transfer type SN. Check QBCF\n");
        fflush(stdout);
        Exit(0);
      }
    }
    if ( cmp[n].getType() == HT_SF )
    {
      if ( (C.KindOfSolver == FLOW_ONLY) || (C.KindOfSolver == THERMAL_FLOW_NATURAL) || (C.KindOfSolver == SOLID_CONDUCTION) )
      {
        Hostonly_ printf("\tInconsistent parameters of combination between Kind of Solver and Heat Transfer type SF. Check QBCF\n");
        fflush(stdout);
        Exit(0);
      }
    }

  }

}


// #################################################################
/* @brief 交点情報の表示（デバッグ）
 * @param [in,out] cut カット情報の配列
 * @param [in]     bid カット情報の配列
 */
void FFV::displayCutInfo(const long long* cut, const int* bid)
{

  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  FILE *fp=NULL;

  if ( !(fp=fopen("cutinfo.txt","w")) )
  {
    Hostonly_ printf("\tSorry, can't open 'cutinfo.txt', write failed.\n");
    Exit(0);
  }

  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {

        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        int bd = bid[m];

        const int b0 = getBit5(bd, X_minus); // (bd >> 0)  & MASK_5;
        const int b1 = getBit5(bd, X_plus);  // (bd >> 5)  & MASK_5;
        const int b2 = getBit5(bd, Y_minus); // (bd >> 10) & MASK_5;
        const int b3 = getBit5(bd, Y_plus);  // (bd >> 15) & MASK_5;
        const int b4 = getBit5(bd, Z_minus); // (bd >> 20) & MASK_5;
        const int b5 = getBit5(bd, Z_plus);  // (bd >> 25) & MASK_5;

        long long pos = cut[m];

        fprintf(fp, "%5d %5d %5d : %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %d %d %d %d %d %d\n",
                i,j,k,
                getCut9(pos, X_minus),
                getCut9(pos, X_plus),
                getCut9(pos, Y_minus),
                getCut9(pos, Y_plus),
                getCut9(pos, Z_minus),
                getCut9(pos, Z_plus),
                b0, b1, b2, b3, b4, b5);

      }
    }
  }

  fflush(fp);
  fclose(fp);
}




// #################################################################
/* @brief メモリ消費情報を表示
 * @param [in]     fp    ファイルポインタ
 * @param [in,out] G_mem グローバルメモリサイズ
 * @param [in]     L_mem ローカルメモリサイズ
 * @param [in]     str   表示用文字列
 */
void FFV::displayMemoryInfo(FILE* fp, double G_mem, double L_mem, const char* str)
{
  if ( numProc > 1 )
  {
    double tmp_memory = G_mem;
    if ( paraMngr->Allreduce(&tmp_memory, &G_mem, 1, MPI_SUM, procGrp) != CPM_SUCCESS ) Exit(0);
  }

  Hostonly_
  {
    FBUtility::MemoryRequirement(str, G_mem, L_mem, fp);
  }

  Hostonly_
  {
    printf("\n\n");
    fprintf(fp, "\n\n");
  }
}


// #################################################################
/* @brief 制御パラメータ，物理パラメータの表示
 * @param [in]  fp   ファイルポインタ
 * @note この関数はマスターノードのみ
 */
void FFV::displayParameters(FILE* fp)
{

  C.printSteerConditions(stdout, &DT, &RF);
  C.printSteerConditions(fp,     &DT, &RF);

  printCriteria(stdout);
  printCriteria(fp);

  F->printSteerConditions(stdout);
  F->printSteerConditions(fp);

  C.printParaConditions(fp,     mat);

  C.printInitValues(fp,     cmp);

  Ex->printPara(stdout, &C);
  Ex->printPara(fp, &C);


  // 外部境界面の開口率を表示
  C.printOuterArea(fp, G_Fcell, G_Acell, G_size);


  // 境界条件のリストと外部境界面のBC設定を表示

  fprintf(fp,"\n----------\n\n");
  fprintf(fp,"\t>> Outer Boundary Conditions\n\n");

  B.printFaceOBC(fp, G_region, BC.exportOBC(), mat);


  // モニタ情報の表示
  if ( C.SamplingMode == ON )
  {
    MO.printMonitorInfo(stdout, "sampling_info.txt", false); // ヘッダのみ

    FILE *fp_mon=NULL;

    if ( !(fp_mon=fopen("sampling_info.txt", "w")) )
    {
      stamped_printf("\tSorry, can't open 'sampling_info.txt' file. Write failed.\n");
      Exit(0);
    }

    MO.printMonitorInfo(fp_mon, "sampling_info.txt", true);  // 詳細モード
    if ( fp_mon ) fclose(fp_mon);
  }
}


// #################################################################
/* @brief BCIndexにビット情報をエンコードする
 */
void FFV::encodeBCindex(FILE* fp)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  // 基本ビット情報（Active, State, コンポ，媒質情報）を全領域についてエンコードする
  V.setBCIndexBase(d_bcd, d_bid, mat, cmp, L_Acell, G_Acell, C.KindOfSolver, C.NoMedium, C.NoCompo);



  // @attention bx[]の同期が必要 >> 以下の処理で隣接セルを参照するため
  if ( numProc > 1 )
  {
    if ( paraMngr->BndCommS3D(d_bcd, ix, jx, kx, gd, 1, procGrp) != CPM_SUCCESS ) Exit(0);
  }

  BC.setBCIperiodic(d_bcd, ensPeriodic);


  // STATEとACTIVEビットのコピー
  V.copyBCIbase(d_bcp, d_bcd);
  V.copyBCIbase(d_cdf, d_bcd);
  V.copyBCIbase(d_bid, d_bcd);



  // 圧力計算のビット情報をエンコードする -----
  V.setBCIndexP(d_bcd, d_bcp, &BC, cmp, C.Mode.Example, d_cut, d_bid, C.NoCompo);




  // 速度計算のビット情報をエンコードする -----
  V.setBCIndexV(d_cdf, &BC, cmp, C.Mode.Example, d_cut, d_bid, d_bcd, C.NoCompo, C.NoMedium, mat);



  // 温度計算のビット情報をエンコードする -----
  if ( C.isHeatProblem() )
  {
    V.setBCIndexH(d_cdf, d_bcd, &BC, C.KindOfSolver, cmp, d_cut, d_bid, C.NoCompo);
  }

  // 内部領域のFluid, Solidのセル数を数える C.Wcell(Local), G_Wcell(global)
  V.countCellState(L_Wcell, G_Wcell, d_bcd, SOLID);
  V.countCellState(L_Fcell, G_Fcell, d_bcd, FLUID);


  // エンコードされているエントリ番号のチェック 1<=order<=31
  V.chkOrder(d_bcd);

}



// #################################################################
/* @brief 固定パラメータの設定
 */
void FFV::fixedParameters()
{
  // 精度
  if ( sizeof(REAL_TYPE) == sizeof(double) )
  {
    C.Mode.Precision = sizeof(double);
  }
  else
  {
    C.Mode.Precision = sizeof(float);
  }

  // 定数
  C.Gravity =9.8; // gravity acceleration
  GM.setSubDivision(20);

}



// #################################################################
/* @brief 並列処理時の各ノードの分割数を集めてファイルに保存する
 */
void FFV::gatherDomainInfo()
{
  // 統計処理の母数
  double d = 1.0 /(double)numProc;
  double r;

  if ( numProc > 1 )
  {
    r = 1.0 /(double)(numProc-1);
  }
  else
  {
    r = 1.0;
  }

  int* m_size = NULL;           ///< 領域分割数
  unsigned long* bf_fcl = NULL; ///< Fluid cell
  unsigned long* bf_wcl = NULL; ///< Solid cell
  unsigned long* bf_acl = NULL; ///< Active cell

  REAL_TYPE* m_org = NULL;      ///< 基点
  REAL_TYPE* m_reg = NULL;      ///< 領域サイズ
  double* bf_srf = NULL;        ///< 表面数

  int* st_buf = NULL; ///<
  int* ed_buf = NULL; ///<


  if( !(m_size = new int[numProc*3]) ) Exit(0);
  if( !(bf_fcl = new unsigned long[numProc]) )   Exit(0);
  if( !(bf_wcl = new unsigned long[numProc]) )   Exit(0);
  if( !(bf_acl = new unsigned long[numProc]) )   Exit(0);

  if( !(m_org  = new REAL_TYPE[numProc*3]) ) Exit(0);
  if( !(m_reg  = new REAL_TYPE[numProc*3]) ) Exit(0);
  if( !(bf_srf = new double[numProc]) )   Exit(0);

  if( !(st_buf = new int [numProc*3]) ) Exit(0);
  if( !(ed_buf = new int [numProc*3]) ) Exit(0);

  // 領域情報の収集
  if ( numProc > 1 )
  {
    if ( paraMngr->Gather(size, 3, m_size, 3, 0, procGrp) != CPM_SUCCESS ) Exit(0);
    if ( paraMngr->Gather(origin, 3, m_org, 3, 0, procGrp) != CPM_SUCCESS ) Exit(0);
    if ( paraMngr->Gather(region, 3, m_reg, 3, 0, procGrp) != CPM_SUCCESS ) Exit(0);
    if ( paraMngr->Gather(&L_Fcell, 1, bf_fcl, 1, 0, procGrp) != CPM_SUCCESS ) Exit(0);
    if ( paraMngr->Gather(&L_Wcell, 1, bf_wcl, 1, 0, procGrp) != CPM_SUCCESS ) Exit(0);
    if ( paraMngr->Gather(&L_Acell, 1, bf_acl, 1, 0, procGrp) != CPM_SUCCESS ) Exit(0);
  }
  else // serial
  {
    memcpy(m_size, size, 3*sizeof(int));
    bf_fcl[0] = L_Fcell;
    bf_wcl[0] = L_Wcell;
    bf_acl[0] = L_Acell;
    memcpy(m_org, origin, 3*sizeof(REAL_TYPE));
    memcpy(m_reg, region, 3*sizeof(REAL_TYPE));
  }


  // Info. of computational domain
  double vol = (double)( (double)G_size[0] * (double)G_size[1] * (double)G_size[2]);
  double srf = 2.0   * ( (double)G_size[0] * (double)G_size[1]
                       + (double)G_size[1] * (double)G_size[2]
                       + (double)G_size[2] * (double)G_size[0]);

  // ローカルノード
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  double m_srf = (double)((ix*jx + jx*kx + kx*ix)) * 2.0;

  if ( nID[X_minus] < 0 ) m_srf -= (double)(jx*kx);  // remove face which does not join communication
  if ( nID[Y_minus] < 0 ) m_srf -= (double)(ix*kx);
  if ( nID[Z_minus] < 0 ) m_srf -= (double)(ix*jx);
  if ( nID[X_plus]  < 0 ) m_srf -= (double)(jx*kx);
  if ( nID[Y_plus]  < 0 ) m_srf -= (double)(ix*kx);
  if ( nID[Z_plus]  < 0 ) m_srf -= (double)(ix*jx);

  if ( numProc > 1 )
  {
    if ( paraMngr->Gather(&m_srf, 1, bf_srf, 1, 0, procGrp) != CPM_SUCCESS ) Exit(0);
  }
  else
  {
    bf_srf[0] = m_srf;
  }

  // mean of domain
  double m_vol = 0;
  double m_efv = 0;
  m_srf = 0;

  for (int i=0; i<numProc; i++) {
    ix = m_size[3*i];
    jx = m_size[3*i+1];
    kx = m_size[3*i+2];
    m_vol += (double)(ix*jx*kx);
    m_srf += bf_srf[i];
    m_efv += (double)bf_acl[i];
  }

  double d_vol = m_vol * d;
  double d_srf = m_srf * d;
  double d_efv = m_efv * d;

  // std. deviation of domain
  double vol_dv = 0.0;
  double srf_dv = 0.0;
  double efv_dv = 0.0;

  double d1, d2, d3;

  for (int i=0; i<numProc; i++)
  {
    ix = m_size[3*i];
    jx = m_size[3*i+1];
    kx = m_size[3*i+2];
    d1 = (double)(ix*jx*kx) - d_vol;
    d2 = bf_srf[i] - d_srf;
    d3 = (double)bf_acl[i] - d_efv;
    vol_dv += d1 * d1;
    srf_dv += d2 * d2;
    efv_dv += d3 * d3;
  }
  vol_dv = sqrt(vol_dv*r);
  srf_dv = sqrt(srf_dv*r);
  efv_dv = sqrt(efv_dv*r);


  FILE *fp=NULL;
  if ( !(fp=fopen("DomainInfo.txt", "w")) )
  {
    stamped_printf("\tSorry, can't open 'DomainInfo.txt' file. Write failed.\n");
    Exit(0);
  }

  // 全体情報の表示
  printGlobalDomain(fp);


  // ローカルノードの情報を表示
  for (int i=0; i<numProc; i++)
  {
    Hostonly_
    {
      fprintf(fp,"\nDomain %4d\n", i);
      fprintf(fp,"\t ix, jx,  kx        [-] =  %13ld %13ld %13ld\n",  m_size[i*3], m_size[i*3+1], m_size[i*3+2]);
      fprintf(fp,"\t(ox, oy, oz)  [m] / [-] = (%13.6e %13.6e %13.6e)  /  (%13.6e %13.6e %13.6e)\n",
              m_org[i*3]*C.RefLength,  m_org[i*3+1]*C.RefLength,  m_org[i*3+2]*C.RefLength, m_org[i*3],  m_org[i*3+1],  m_org[i*3+2]);
      fprintf(fp,"\t(Lx, Ly, Lz)  [m] / [-] = (%13.6e %13.6e %13.6e)  /  (%13.6e %13.6e %13.6e)\n",
              m_reg[i*3]*C.RefLength,  m_reg[i*3+1]*C.RefLength,  m_reg[i*3+2]*C.RefLength, m_reg[i*3],  m_reg[i*3+1],  m_reg[i*3+2]);

      if (C.NoCompo != 0) fprintf(fp, "\n\t no            Label    i_st    i_ed    j_st    j_ed    k_st    k_ed\n");
    }

    if ( numProc > 1 )
    {
      for (int n=1; n<=C.NoCompo; n++)
      {
        if ( !cmp[n].isKindMedium() )
        {
          if( paraMngr->Gather(cmp[n].getBbox_st(), 3, st_buf, 3, 0, procGrp) != CPM_SUCCESS ) Exit(0);
          if( paraMngr->Gather(cmp[n].getBbox_ed(), 3, ed_buf, 3, 0, procGrp) != CPM_SUCCESS ) Exit(0);
        }
        Hostonly_
        {
          fprintf(fp,"\t%3d %16s %7d %7d %7d %7d %7d %7d\n",
                  n, cmp[n].alias.c_str(), st_buf[i*3], ed_buf[i*3], st_buf[i*3+1], ed_buf[i*3+1], st_buf[i*3+2], ed_buf[i*3+2]);
        }
      }
    }
    else // serial
    {
      int *st, *ed;
      for (int n=1; n<=C.NoCompo; n++)
      {
        if ( !cmp[n].isKindMedium() )
        {
          st = cmp[n].getBbox_st();
          ed = cmp[n].getBbox_ed();

          Hostonly_
          {
            fprintf(fp,"\t%3d %16s %7d %7d %7d %7d %7d %7d\n",
                  n, cmp[n].alias.c_str(), st[0], ed[0], st[1], ed[1], st[2], ed[2]);
          }
        }
      }
    }
  }

  Hostonly_
  {
    fprintf(fp, "\n");
    fprintf(fp,"\n\t--------------------------------------------------\n");
    fprintf(fp,"\tReport of Whole Domain Statistics\n");
    fprintf(fp,"\tDomain size         = %7d %7d %7d\n", G_size[0], G_size[1], G_size[2]);
    fprintf(fp,"\tNumber of voxels    = %12.6e\n", vol);
    fprintf(fp,"\tNumber of surface   = %12.6e\n", srf);
    fprintf(fp,"\tEffective voxels    = %12.6e (%6.2f%%)\n", (REAL_TYPE)G_Acell, 100.0*(REAL_TYPE)G_Acell/vol);
    fprintf(fp,"\tFluid voxels        = %12.6e (%6.2f%%)\n", (REAL_TYPE)G_Fcell, 100.0*(REAL_TYPE)G_Fcell/vol);
    fprintf(fp,"\tWall  voxels        = %12.6e (%6.2f%%)\n", (REAL_TYPE)G_Wcell, 100.0*(REAL_TYPE)G_Wcell/vol);

    if ( numProc == 1 )
    {
      fprintf(fp,"\tDivision :          = %d : %s\n", numProc, "Serial");
    }
    else
    {
      fprintf(fp,"\tDivision :          = %d : %s\n", numProc, "Equal segregation");
    }

    fprintf(fp,"\n\t--------------------------------------------------\n");
    fprintf(fp,"\tDomain Statistics per MPI process\n");
    fprintf(fp,"\tMean volume in each domain           = %12.6e\n", d_vol);
    fprintf(fp,"\tStd. deviation of domain             = %12.6e\n", vol_dv);
    fprintf(fp,"\tMean comm. in each domain            = %12.6e\n", d_srf);
    fprintf(fp,"\tStd. deviation of surface            = %12.6e\n", srf_dv);
    fprintf(fp,"\tMean effective volume in each domain = %12.6e\n", d_efv);
    fprintf(fp,"\tStd. deviation of effective volume   = %12.6e\n", efv_dv);
    fprintf(fp,"\n");

    fprintf(fp,"\tDomain :     ix     jx     kx       Volume Vol_dv[%%]      Surface Srf_dv[%%] Fluid[%%] Solid[%%]      Eff_Vol Eff_Vol_dv[%%]      Eff_Srf Eff_srf_dv[%%]  Itr_scheme\n");
    fprintf(fp,"\t----------------------------------------------------------------------------------------------------\n");

    double tmp_vol, tmp_acl, tmp_fcl, tmp_wcl;

    for (int i=0; i<numProc; i++)
    {
      ix = m_size[3*i];
      jx = m_size[3*i+1];
      kx = m_size[3*i+2];
      tmp_vol = (double)(ix*jx*kx);
      tmp_acl = (double)bf_acl[i];
      tmp_fcl = (double)bf_fcl[i];
      tmp_wcl = (double)bf_wcl[i];
      fprintf(fp,"\t%6d : %6d %6d %6d ", i, ix, jx, kx);
      fprintf(fp,"%12.4e  %8.3f ", tmp_vol, 100.0*(tmp_vol-d_vol)/d_vol);
      fprintf(fp,"%12.4e  %8.3f ", bf_srf[i], (d_srf == 0.0) ? 0.0 : 100.0*(bf_srf[i]-d_srf)/d_srf);
      fprintf(fp,"%8.3f %8.3f ", 100.0*tmp_fcl/tmp_vol, 100.0*tmp_wcl/tmp_vol);
      fprintf(fp,"%12.4e      %8.3f \n", tmp_acl, 100.0*(tmp_acl-d_efv)/d_efv);
    }
    fprintf(fp,"\t---------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
  }

  if (fp) fclose(fp);


  if( m_size ) { delete [] m_size; m_size=NULL; }
  if( m_org  ) { delete [] m_org;  m_org =NULL; }
  if( m_reg  ) { delete [] m_reg;  m_reg =NULL; }
  if( bf_srf ) { delete [] bf_srf; bf_srf=NULL; }
  if( bf_fcl ) { delete [] bf_fcl; bf_fcl=NULL; }
  if( bf_wcl ) { delete [] bf_wcl; bf_wcl=NULL; }
  if( bf_acl ) { delete [] bf_acl; bf_acl=NULL; }
  if( st_buf ) { delete [] st_buf; st_buf=NULL; }
  if( ed_buf ) { delete [] ed_buf; ed_buf=NULL; }
}



// #################################################################
/* @brief 交点情報のグリフを生成する
 * @param [in] cut   カットの配列
 * @param [in] bid   境界IDの配列
 * @param [in] fp    file pointer
 * @param [in] m_st  範囲指定　デバッグ用
 * @param [in] m_ed  範囲指定　デバッグ用
 */
void FFV::generateGlyph(const long long* cut, const int* bid, FILE* fp, int* m_st, int* m_ed)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  // 交点の総数を求める
  unsigned global_cut=0;  /// 全カット数
  unsigned local_cut=0;   /// 担当プロセスのカット数
  unsigned g=0;

  int st[3], ed[3];

  if ( m_st == NULL )
  {
    st[0] = 1;
    st[1] = 1;
    st[2] = 1;
    ed[0] = ix;
    ed[1] = jx;
    ed[2] = kx;
  }
  else
  {
    st[0] = m_st[0];
    st[1] = m_st[1];
    st[2] = m_st[2];
    ed[0] = m_ed[0];
    ed[1] = m_ed[1];
    ed[2] = m_ed[2];
  }


  for (int k=st[2]; k<=ed[2]; k++) {
    for (int j=st[1]; j<=ed[1]; j++) {
      for (int i=st[0]; i<=ed[0]; i++) {

        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        int qq = bid[m];

        if ( IS_CUT(qq) ) // カットがあるか，IDによる判定
        {
          if (getBit5(qq, 0) != 0) g++;
          if (getBit5(qq, 1) != 0) g++;
          if (getBit5(qq, 2) != 0) g++;
          if (getBit5(qq, 3) != 0) g++;
          if (getBit5(qq, 4) != 0) g++;
          if (getBit5(qq, 5) != 0) g++;
        }

      }
    }
  }

  global_cut = local_cut = g;

  if ( numProc > 1 )
  {
    unsigned tmp = global_cut;
    if ( paraMngr->Allreduce(&tmp, &global_cut, 1, MPI_SUM, procGrp) != CPM_SUCCESS ) Exit(0);
  }

  Hostonly_
  {
    printf("\tGlyph : Number of Cut points = %u\n", global_cut);
    fprintf(fp, "\tGlyph : Number of Cut points = %u\n", global_cut);
  }


  // 格子幅（有次元）
  float m_pch[3] = {
    (float)pitch[0]*C.RefLength,
    (float)pitch[1]*C.RefLength,
    (float)pitch[2]*C.RefLength
  };

  // サブドメインの起点座標（有次元）
  float m_org[3] = {
    (float)origin[0]*C.RefLength,
    (float)origin[1]*C.RefLength,
    (float)origin[2]*C.RefLength
  };


  // ポリゴンをストアする配列を確保
  Glyph glyph(m_pch, m_org, local_cut, myRank);


  // グリフの生成モード
  bool inner_only = false;
  if (C.Hide.GlyphOutput == 2) inner_only=true;

  // カット点毎にグリフのポリゴン要素を生成し，配列にストア
  Vec3i idx;

  for (int k=st[2]; k<=ed[2]; k++) {
    for (int j=st[1]; j<=ed[1]; j++) {
      for (int i=st[0]; i<=ed[0]; i++) {

        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        int qq = bid[m];

        if ( IS_CUT(qq) ) // カットがあるか，IDによる判定
        {
          const long long pos = cut[m];

          idx.assign(i, j, k);

          if ( inner_only )
          {
            if (i != 1 )
            {
              if (getBit5(qq, 0) != 0) glyph.generateVertex(idx, pos, X_minus, qq);
            }
            if ( i != ix )
            {
              if (getBit5(qq, 1) != 0) glyph.generateVertex(idx, pos, X_plus, qq);
            }
            if ( j != 1 )
            {
              if (getBit5(qq, 2) != 0) glyph.generateVertex(idx, pos, Y_minus, qq);
            }
            if ( j != jx )
            {
              if (getBit5(qq, 3) != 0) glyph.generateVertex(idx, pos, Y_plus, qq);
            }
            if ( k != 1 )
            {
              if (getBit5(qq, 4) != 0) glyph.generateVertex(idx, pos, Z_minus, qq);
            }
            if ( k != kx )
            {
              if (getBit5(qq, 5) != 0) glyph.generateVertex(idx, pos, Z_plus, qq);
            }
          }
          else
          {
            if (getBit5(qq, 0) != 0) glyph.generateVertex(idx, pos, X_minus, qq);
            if (getBit5(qq, 1) != 0) glyph.generateVertex(idx, pos, X_plus , qq);
            if (getBit5(qq, 2) != 0) glyph.generateVertex(idx, pos, Y_minus, qq);
            if (getBit5(qq, 3) != 0) glyph.generateVertex(idx, pos, Y_plus , qq);
            if (getBit5(qq, 4) != 0) glyph.generateVertex(idx, pos, Z_minus, qq);
            if (getBit5(qq, 5) != 0) glyph.generateVertex(idx, pos, Z_plus , qq);
          }
        }

      }
    }
  }


  // ポリゴンの出力
  glyph.writeBinary("CutGlyph");

}




// #################################################################
// * @brief FielIO classの同定
// * @param [in] tpCntl  テキストパーサー
void FFV::identifyFIO(TextParser* tpCntl)
{
  string str;
  string label, leaf;
  int Format=0;


  // フォーマットパラメータの取得
  label = "/Output/Data/Format";

  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }

  if     ( !strcasecmp(str.c_str(), "sph") )    Format = sph_fmt;
  else if( !strcasecmp(str.c_str(), "bov") )    Format = bov_fmt;
  else if( !strcasecmp(str.c_str(), "plot3d") ) Format = plt3d_fun_fmt;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
    Exit(0);
  }


  // インスタンス
  if      ( Format == sph_fmt )
  {
    F = dynamic_cast<IO_BASE*>(new SPH);
    F->setFormat(sph_fmt);
  }
  else if ( Format == bov_fmt )
  {
    F->setFormat(bov_fmt);
  }
  else if ( Format == plt3d_fun_fmt )
  {
    F = dynamic_cast<IO_BASE*>(new PLT3D);
    F->setFormat(plt3d_fun_fmt);
  }

}


// #################################################################
/* @brief Intrisic classの同定
 * @param [in] fp  ファイル出力ポインタ
 */
void FFV::identifyExample(FILE* fp)
{

  // 例題クラスの実体をインスタンスし，Exにポイントする
  if      ( C.Mode.Example == id_PPLT2D)   Ex = dynamic_cast<Intrinsic*>(new IP_PPLT2D);
  else if ( C.Mode.Example == id_Duct )    Ex = dynamic_cast<Intrinsic*>(new IP_Duct);
  else if ( C.Mode.Example == id_PMT )     Ex = dynamic_cast<Intrinsic*>(new IP_PMT);
  else if ( C.Mode.Example == id_Rect )    Ex = dynamic_cast<Intrinsic*>(new IP_Rect);
  else if ( C.Mode.Example == id_Cylinder) Ex = dynamic_cast<Intrinsic*>(new IP_Cylinder);
  else if ( C.Mode.Example == id_Step )    Ex = dynamic_cast<Intrinsic*>(new IP_Step);
  else if ( C.Mode.Example == id_Polygon ) Ex = new Intrinsic;
  else if ( C.Mode.Example == id_Sphere )  Ex = dynamic_cast<Intrinsic*>(new IP_Sphere);
  else if ( C.Mode.Example == id_Jet )     Ex = dynamic_cast<Intrinsic*>(new IP_Jet);
  else
  {
    Hostonly_
    {
      stamped_printf(    "\tInvalid keyword is described for Exmple definition\n");
      stamped_fprintf(fp,"\tInvalid keyword is described for Exmple definition\n");
    }
    Exit(0);
  }


  // 組み込み例題クラス名を表示
  Hostonly_ Ex->printExample(fp, C.Mode.Example);

}


// #################################################################
/**
 * @brief 線形ソルバを特定
 * @param [in] tpCntl  テキストパーサー
 */
void FFV::identifyLinearSolver(TextParser* tpCntl)
{
  string str, label;
  double tmp=0.0;
  int inp=0;

  label = "/Iteration/DivMaxIteration";

  if ( !(tpCntl->getInspectedValue(label, inp)) )
  {
    Exit(0);
  }
  DivC.MaxIteration = inp;


  label = "/Iteration/DivCriterion";

  if ( !(tpCntl->getInspectedValue(label, tmp)) )
  {
    Exit(0);
  }
  DivC.divEPS = tmp;


  label = "/Iteration/DivNorm";

  if ( !(tpCntl->getInspectedValue(label, str)) )
  {
    Exit(0);
  }

  if ( !strcasecmp(str.c_str(), "L2") )
  {
    DivC.divType = nrm_div_l2;
  }
  else if ( !strcasecmp(str.c_str(), "max") )
  {
    DivC.divType = nrm_div_max;
  }
  else
  {
    Exit(0);
  }


  switch (C.AlgorithmF)
  {
    case Flow_FS_EE_EE:
    case Flow_FS_AB2:
      LS_setParameter(tpCntl, ic_prs1, "/Iteration/Pressure");
      break;

    case Flow_FS_AB_CN:
      LS_setParameter(tpCntl, ic_prs1, "/Iteration/Pressure");
      LS_setParameter(tpCntl, ic_vel1, "/Iteration/Velocity");
      break;

    case Flow_FS_RK_CN:
      LS_setParameter(tpCntl, ic_prs1, "/Iteration/Pressure");
      LS_setParameter(tpCntl, ic_prs2, "/Iteration/Pressure2nd");
      LS_setParameter(tpCntl, ic_vel1, "/Iteration/Velocity");
      break;

    default:
      Exit(0);
  }


  if ( !C.isHeatProblem() ) return;


  switch (C.AlgorithmH)
  {
    case Heat_EE_EE:
      break;

    case Heat_EE_EI:
      LS_setParameter(tpCntl, ic_tmp1, "/Iteration/Temperature");
      break;

    default:
      Exit(0);
  }

}



// #################################################################
/* @brief 距離情報の初期化
 */
void FFV::initCutInfo()
{
  size_t n_cell[3];

  for (int i=0; i<3; i++)
  {
    n_cell[i] = (size_t)(size[i] + 2*guide); // 分割数+ガイドセル
  }
  size_t size_n_cell = n_cell[0] * n_cell[1] * n_cell[2];


  // 初期値のセット
  for (size_t i=0; i<size_n_cell; i++)
  {
    for (int dir=0; dir<6; dir++)
    {
      initBit9(d_cut[i], dir);
    }
  }

}





// #################################################################
/* @brief インターバルの初期化
 */
void FFV::initInterval()
{
  double m_dt = DT.get_DT();
  unsigned m_Session_StartStep;   ///< セッションの開始ステップ

  if ( C.Interval[Control::tg_compute].getMode() == IntervalManager::By_step )
  {
    m_Session_StartStep = C.Interval[Control::tg_compute].getStartStep();
  }
  else
  {
    m_Session_StartStep = (unsigned)ceil( C.Interval[Control::tg_compute].getStartTime()  / (m_dt*C.Tscale) );
  }


  // セッションの最終ステップ
  if ( C.Interval[Control::tg_compute].getMode() == IntervalManager::By_step )
  {
    Session_LastStep = C.Interval[Control::tg_compute].getLastStep();
  }
  else
  {
    Session_LastStep = (unsigned)ceil( C.Interval[Control::tg_compute].getLastTime() / (m_dt*C.Tscale) );
  }



  // セッションの開始・終了時刻をセット >> @see Control::getTimeControl()
  for (int i=0; i<Control::tg_END; i++)
  {
    if ( (i != Control::tg_statistic) && (i != Control::tg_compute) )
    {
      C.Interval[i].setStart(m_Session_StartStep);
      C.Interval[i].setLast(Session_LastStep);
    }
  }



  // 入力モードが有次元の場合に，無次元に変換 >> 時制がBy_timeの場合のみ
  if ( C.Unit.Param == DIMENSIONAL )
  {
    for (int i=Control::tg_compute; i<=Control::tg_accelra; i++)
    {
      C.Interval[i].normalizeTime(C.Tscale);
    }

    if ( C.SamplingMode == ON )
    {
      C.Interval[Control::tg_sampled].normalizeTime(C.Tscale);
    }
  }

  // Reference frame
  RF.setAccel( C.Interval[Control::tg_accelra].getIntervalTime() );



  // インターバルの初期化
  double m_tm    = CurrentTime;  // 積算時間 Restart()で設定
  unsigned m_stp = CurrentStep;  // 積算ステップ数

  for (int i=Control::tg_compute; i<=Control::tg_accelra; i++)
  {
    if ( !C.Interval[i].initTrigger(m_stp, m_tm, m_dt) )
    {
      Hostonly_ printf("\t Error : initialize timing trigger [no=%d].\n", i);
      Exit(0);
    }
  }

  if ( C.SamplingMode == ON )
  {
    if ( !C.Interval[Control::tg_sampled].initTrigger(m_stp, m_tm, m_dt) )
    {
      Hostonly_ printf("\t Error : initialize timing trigger [tg_sampled].\n");
      Exit(0);
    }
  }

}


// #################################################################
/**
 * @brief 線形ソルバー関連の初期化
 * @param [in,out] m_total  本計算用のメモリリサイズ
 * @param [in]     tpCntl   テキストパーサーのツリーポインタ
 */
void FFV::LS_initialize(double& TotalMemory, TextParser* tpCntl)
{
  // communication buffer
  switch (LS[ic_prs1].getLS())
  {
    case SOR2SMA:
    case GMRES:
    case PCG:
    case BiCGSTAB:
      allocate_SOR2SMA_buffer(TotalMemory);
      break;
  }


  // extra arrays for Krylov subspace
  switch (LS[ic_prs1].getLS())
  {
    case GMRES:
      allocArray_Krylov(TotalMemory);
      break;

    case PCG:
      allocArray_PCG(TotalMemory);
      break;

    case BiCGSTAB:
      allocArray_BiCGstab(TotalMemory);
      if ( LS[ic_prs1].isPreconditioned() )
      {
        allocArray_BiCGSTABwithPreconditioning(TotalMemory);
      }
      break;
  }


  // Initialize

  for (int i=0; i<ic_END; i++)
  {
    if ( LS[i].getLS() != 0)
    {
      LS[i].Initialize(&C,
                       &BC,
                       ModeTiming,
                       face_comm_size,
                       &PM,
                       d_bcp,
                       d_bcd,
                       d_pcg_p,
                       d_pcg_p_,
                       d_pcg_r,
                       d_pcg_r0,
                       d_pcg_q,
                       d_pcg_s,
                       d_pcg_s_,
                       d_pcg_t,
                       d_pcg_t_,
                       ensPeriodic,
                       cf_sz,
                       cf_x,
                       cf_y,
                       cf_z);
    }
  }

}


// #################################################################
/**
 * @brief 線形ソルバを特定し，パラメータをセットする
 * @param [in] tpCntl テキストパーサーのポインタ
 * @param [in] odr    制御クラス配列の番号
 * @param [in] label  探索ラベル
 */
void FFV::LS_setParameter(TextParser* tpCntl, const int odr, const string label)
{
  string str;

  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
    Hostonly_ printf("\tParsing error : No '%s'\n", label.c_str());
    Exit(0);
  }

  //C.copyCriteria(dynamic_cast<IterationCtl*>(&LS[odr]), str);

  for (int i=0; i<C.NoBaseLS; i++)
  {
    if ( !strcasecmp( str.c_str(), C.Criteria[i].getAlias().c_str() ))
    {
      LS[odr].copy(&C.Criteria[i]);
    }
  }

}



// #################################################################
// 初期擾乱
void FFV::perturbation()
{
  REAL_TYPE width  = C.LES.ChannelWidth;
  REAL_TYPE Re_tau = C.LES.TurbulentReynoldsNum;
  REAL_TYPE Ubar   = C.LES.BulkVelocity;
  REAL_TYPE visc   = C.RefKviscosity;
  int mode;

  switch (C.LES.ChannelDir)
  {
    case X_minus:
    case X_plus:
      Hostonly_ printf("Currently, not supported.\n");
      Exit(0);
      break;

    case Y_minus:
    case Y_plus:
      mode = 1;
      perturb_u_y_ (d_v, size, &guide, pitch, origin, &width, &Re_tau, &Ubar, &visc, &mode);

      mode = 2;
      perturb_u_y_ (d_vf, size, &guide, pitch, origin, &width, &Re_tau, &Ubar, &visc, &mode);
      break;

    case Z_minus:
    case Z_plus:
      mode = 1;
      perturb_u_z_ (d_v, size, &guide, pitch, origin, &width, &Re_tau, &Ubar, &visc, &mode);

      mode = 2;
      perturb_u_z_ (d_vf, size, &guide, pitch, origin, &width, &Re_tau, &Ubar, &visc, &mode);
      break;

    default:
      Exit(0);
      break;
  }

}


// #################################################################
/* @brief 履歴の出力準備
 */
void FFV::prepHistoryOutput()
{
  // マスターノードでの履歴出力準備
  H = new History(&C);

  // gnu compilerでdynamic_castでエラーがでるため，一度配列に入れて渡す
  int container[2*ic_END];
  for (int i=0; i<ic_END; i++)
  {
    container[2*i+0] = LS[i].getResType();
    container[2*i+1] = LS[i].getErrType();
  }

  Hostonly_
  {
    H->printHistoryTitle(stdout, container, &C, &DivC, true);

    // コンポーネント情報
    if ( C.Mode.Log_Base == ON )
    {
      // 基本情報
      if ( !(fp_b=fopen("history_base.txt", "w")) )
      {
        stamped_printf("\tSorry, can't open 'history_base.txt' file. Write failed.\n");
        Exit(0);
      }
      H->printHistoryTitle(fp_b, container, &C, &DivC, true);

      // コンポーネント履歴情報
      if ( C.EnsCompo.monitor )
      {
        if ( !(fp_c=fopen("history_compo.txt", "w")) )
        {
          stamped_printf("\tSorry, can't open 'history_compo.txt' file. Write failed.\n");
          Exit(0);
        }
        H->printHistoryCompoTitle(fp_c, cmp, &C);
      }

      // 流量収支情報　
      if ( !(fp_d=fopen("history_domainflux.txt", "w")) )
      {
        stamped_printf("\tSorry, can't open 'history_domainflux.txt' file. Write failed.\n");
        Exit(0);
      }
      H->printHistoryDomfxTitle(fp_d, &C);


      // 力の履歴情報（コンポーネント毎）
      H->printHistoryForceTitle(cmp);
    }


    // 壁面情報　history_wall.log
    if ( C.Mode.Log_Wall == ON )
    {
      if ( !(fp_w=fopen("history_log_wall.txt", "w")) )
      {
				stamped_printf("\tSorry, can't open 'history_log_wall.txt' file.\n");
        Exit(0);
      }
      H->printHistoryWallTitle(fp_w);
    }

    // CCNVfile
    if ( C.Mode.CCNV == ON )
    {
      H->printCCNVtitle(container, &C);
    }
  }

}


// #################################################################
// グローバルな領域情報を表示する
void FFV::printGlobalDomain(FILE* fp)
{
  REAL_TYPE PB=0.0, TB=0.0, GB=0.0, MB=0.0, KB=0.0, total=0.0;
  KB = 1000.0;
  MB = 1000.0*KB;
  GB = 1000.0*MB;
  TB = 1000.0*GB;
  PB = 1000.0*TB;

  fprintf(fp,"\timax, jmax, kmax    = %13d %13d %13d     >> ",
          G_size[0],
          G_size[1],
          G_size[2]);

  total = (REAL_TYPE)G_size[0] * (REAL_TYPE)G_size[1] * (REAL_TYPE)G_size[2];

  if ( total > PB ) {
    fprintf (fp,"%6.2f (P cells)\n", total / PB);
  }
  else if ( total > TB ) {
    fprintf (fp,"%6.2f (T cells)\n", total / TB);
  }
  else if ( total > GB ) {
    fprintf (fp,"%6.2f (G cells)\n", total / GB);
  }
  else if ( total > MB ) {
    fprintf (fp,"%6.2f (M cells)\n", total / MB);
  }
  else if ( total > KB ) {
    fprintf (fp,"%6.2f (K cells)\n", total / KB);
  }
  else if ( total <= KB ){
    fprintf (fp,"%6.2f (cells)\n", total);
  }
  fprintf(fp,"\n");

  fprintf(fp,"\t(dx, dy, dz)  [m] / [-] = (%13.6e %13.6e %13.6e)  /  (%13.6e %13.6e %13.6e)\n",
          pitchD[0],
          pitchD[1],
          pitchD[2],
          pitch[0],
          pitch[1],
          pitch[2]);

  fprintf(fp,"\t(ox, oy, oz)  [m] / [-] = (%13.6e %13.6e %13.6e)  /  (%13.6e %13.6e %13.6e)\n",
          G_originD[0],
          G_originD[1],
          G_originD[2],
          G_origin[0],
          G_origin[1],
          G_origin[2]);

  fprintf(fp,"\t(Lx, Ly, Lz)  [m] / [-] = (%13.6e %13.6e %13.6e)  /  (%13.6e %13.6e %13.6e)\n",
          G_regionD[0],
          G_regionD[1],
          G_regionD[2],
          G_region[0],
          G_region[1],
          G_region[2]);
  fprintf(fp,"\n");

  fflush(fp);
}


// #################################################################
/* @brief 境界条件を読み込み，Controlクラスに保持する
 */
void FFV::setBCinfo()
{

  B.setControlVars(&C);


  // KOSと媒質の状態の整合性をチェックし，流体と固体の媒質数をカウント
  B.countMedium(&C, mat);


  // パラメータファイルの情報を元にCompoListの情報を設定する
  B.loadBCs(&C, mat, cmp);


  // 外部境界条件をBCクラスに保持する
  B.setOuterBC(BC.exportOBC(), cmp, ensPeriodic);



#if 0
  Hostonly_ B.checkList(mat, cmp);
#endif


  // 各コンポーネントが存在するかどうかを保持し, OBSTACLEの個数を返す
  num_obstacle = C.setExistComponent(cmp, BC.exportOBC(), global_obstacle);


  // KOSと境界条件種類の整合性をチェック
  B.chkBCconsistency(C.KindOfSolver, cmp);


  // RefMediumがMediumList中にあるかどうかをチェックし、RefMatを設定
  if ( (C.RefMat = U.findIDfromLabel(mat, C.NoMedium, C.RefMedium)) == 0 )
  {
    Hostonly_
    {
      printf("/Referece/Medium = \"%s\" is not listed in MediumTable.\n", C.RefMedium.c_str());
    }
    Exit(0);
  }



  // 各軸方向のフィル抑止モード（Periodic, Symmetric時の対策）
  for (int face=0; face<NOFACE; face++)
  {
    if( nID[face] >= 0 ) continue;

    BoundaryOuter* m_obc = BC.exportOBC(face);

    switch (m_obc->getClass())
    {
      case OBC_SYMMETRIC:
      case OBC_PERIODIC:
        if (face==X_minus || face==X_plus)
        {
          GM.FillSuppress[0] = OFF;
        }
        else if (face==Y_minus || face==Y_plus)
        {
          GM.FillSuppress[1] = OFF;
        }
        else
        {
          GM.FillSuppress[2] = OFF;
        }
        break;

      default:
        break;
    }
  }


  // 無次元パラメータの設定
  C.setRefParameters(mat, &RF);


  // パラメータの無次元化（正規化）に必要な参照物理量の設定
  B.setRefMediumProperty(C.RefDensity, C.RefSpecificHeat);


  // 演算用のパラメータ配列にコピーし，無次元化、Fortranからアクセスできるように一次元配列
  REAL_TYPE lmd0 = C.RefDensity * C.RefSpecificHeat * C.RefVelocity * C.RefLength;

  for (int n=1; n<=C.NoCompo; n++)
  {
    mat_tbl[n*3+0] = mat[n].P[p_density] / C.RefDensity;
    mat_tbl[n*3+1] = mat[n].P[p_specific_heat] / C.RefSpecificHeat;
    mat_tbl[n*3+2] = mat[n].P[p_thermal_conductivity] / lmd0;
  }


  // コンポーネントの無次元速度パラメータ
  const REAL_TYPE c_pai = (REAL_TYPE)(2.0*asin(1.0));

  for (int n=1; n<=C.NoCompo; n++)
  {
    if ( cmp[n].getType() == SPEC_VEL )
    {
      vec_tbl[n*7+0] = cmp[n].nv[0]; // normal vector
      vec_tbl[n*7+1] = cmp[n].nv[1];
      vec_tbl[n*7+2] = cmp[n].nv[2];
      vec_tbl[n*7+3] = cmp[n].ca[CompoList::amplitude] / C.RefVelocity;
      vec_tbl[n*7+4] = cmp[n].ca[CompoList::frequency] * C.RefLength / C.RefVelocity * 2.0*c_pai;
      vec_tbl[n*7+5] = cmp[n].ca[CompoList::initphase];
      vec_tbl[n*7+6] = cmp[n].ca[CompoList::bias] / C.RefVelocity;
    }
    else if ( cmp[n].getType() == SOLIDREV )
    {
      REAL_TYPE omg = cmp[n].ca[0] * C.RefLength / C.RefVelocity; // [rad]
      vec_tbl[n*7+0] = cmp[n].nv[0] * omg; // normal vector x angular velocity
      vec_tbl[n*7+1] = cmp[n].nv[1] * omg;
      vec_tbl[n*7+2] = cmp[n].nv[2] * omg;
      vec_tbl[n*7+3] = cmp[n].oc[0] / C.RefLength;
      vec_tbl[n*7+4] = cmp[n].oc[1] / C.RefLength;
      vec_tbl[n*7+5] = cmp[n].oc[2] / C.RefLength;
      vec_tbl[n*7+6] = 0.0;
    }
    else
    {
      vec_tbl[n*7+0] = 0.0;
      vec_tbl[n*7+1] = 0.0;
      vec_tbl[n*7+2] = 0.0;
      vec_tbl[n*7+3] = 0.0;
      vec_tbl[n*7+4] = 0.0;
      vec_tbl[n*7+5] = 0.0;
      vec_tbl[n*7+6] = 0.0;
    }

  }


  // 無次元媒質情報をコピー
  BC.copyNDmatTable(C.NoCompo, mat_tbl, vec_tbl);


  // 温度計算の場合の媒質毎の初期値指定オプション
  if ( C.isHeatProblem() )
  {
    B.getInitTempOfMedium(cmp, &C);
  }


  // 境界条件の媒質に初期温度にコピー
  for (int n=1; n<=C.NoCompo; n++)
  {
    // OBSTACLEより上のコンポーネントは境界条件
    if ( cmp[n].getType() > OBSTACLE )
    {
      for (int j=1; j<=C.NoMedium; j++)
      {
        if ( !strcasecmp(cmp[j].alias.c_str(), cmp[n].medium.c_str()) )
        {
          cmp[n].setInitTemp( cmp[j].getInitTemp() );
          //printf("init[%d] %f %s\n", n, cmp[n].getInitTemp(), cmp[n].medium.c_str());
        }
      }
    }
  }

}


// #################################################################
/* @brief 回転体をセット
 */
void FFV::setComponentSR()
{
  // 注意！！
  // ###   無次元パラメータを渡す  ###
  CompoFraction CF(size, guide, myRank, pitch, origin, 20);

  for (int n=1; n<=C.NoCompo; n++)
  {
    REAL_TYPE center[3];
    REAL_TYPE dirvec[3];
    REAL_TYPE depth, p1, p2;

    // 形状パラメータのセット
    int Ctype = cmp[n].getType();

    if ( Ctype==SOLIDREV )
    {
      center[0] = cmp[n].oc[0] / C.RefLength;
      center[1] = cmp[n].oc[1] / C.RefLength;
      center[2] = cmp[n].oc[2] / C.RefLength;
      dirvec[0] = cmp[n].dr[0] / C.RefLength;
      dirvec[1] = cmp[n].dr[1] / C.RefLength;
      dirvec[2] = cmp[n].dr[2] / C.RefLength;
      depth = cmp[n].depth  / C.RefLength;
      p1    = cmp[n].shp_p1 / C.RefLength;
      p2    = cmp[n].shp_p2 / C.RefLength;


      // getAttrb()にはオプションが入っている "palte" or "cylinder"
      CF.setShapeParam(cmp[n].nv, center, depth, p1, p2, cmp[n].getAttrb());


      // 回転角度の計算
      CF.getAngle();


      // bboxと投影面積、インデクスの計算
      int f_st[3], f_ed[3];
      cmp[n].area = CF.getBboxArea(f_st, f_ed);


      // インデクスのサイズ登録と存在フラグ
      cmp[n].setBbox(f_st, f_ed);
      cmp[n].setEnsLocal(ON);
      // ここで得られたst[],ed[]の値は、VoxInfo.CのencVIBCrev()で書き換えられる


      double flop=0.0;

      // 交点計算
      flop = 0.0;
      CF.intersectCylinder(f_st, f_ed, d_bid, d_cut, n);


      // ########## 確認のための出力
#if 0
      REAL_TYPE org[3], pit[3];

      //  ガイドセルがある場合(GuideOut != 0)にオリジナルポイントを調整
      for (int i=0; i<3; i++)
      {
        org[i] = C.org[i] - C.dx[i]*(REAL_TYPE)C.GuideOut;
        pit[i] = C.dx[i];
      }

      // 出力ファイルの指定が有次元の場合
      if ( C.Unit.File == DIMENSIONAL )
      {
        for (int i=0; i<3; i++)
        {
          org[i] *= C.RefLength;
          pit[i] *= C.RefLength;
        }
      }
      writeRawSPH(cvf, size, guide, 0, org, pit, sizeof(REAL_TYPE));
#endif
      // ##########
    }

  } // loop - NoCmpo

}



// #################################################################
/* @brief HEX,FANコンポーネントなどの体積率とbboxなどをセット
 */
void FFV::setComponentVF()
{
  // 有次元パラメータを渡す
  CompoFraction CF(size, guide, myRank, pitchD, originD, 20);

  for (int n=1; n<=C.NoCompo; n++)
  {
    REAL_TYPE center[3];
    REAL_TYPE dirvec[3];
    REAL_TYPE depth, p1, p2;

    // 形状パラメータのセット
    int Ctype = cmp[n].getType();

    if ( Ctype==HEX || Ctype==FAN )
    {
      switch ( Ctype )
      {
        case HEX:
          center[0] = cmp[n].oc[0] / C.RefLength;
          center[1] = cmp[n].oc[1] / C.RefLength;
          center[2] = cmp[n].oc[2] / C.RefLength;
          dirvec[0] = cmp[n].dr[0] / C.RefLength;
          dirvec[1] = cmp[n].dr[1] / C.RefLength;
          dirvec[2] = cmp[n].dr[2] / C.RefLength;
          depth = cmp[n].depth  / C.RefLength;
          p1    = cmp[n].shp_p1 / C.RefLength;
          p2    = cmp[n].shp_p2 / C.RefLength;
          CF.setShapeParam(cmp[n].nv, center, dirvec, depth, p1, p2);
          break;

        case FAN:
          center[0] = cmp[n].oc[0] / C.RefLength;
          center[1] = cmp[n].oc[1] / C.RefLength;
          center[2] = cmp[n].oc[2] / C.RefLength;
          depth = cmp[n].depth  / C.RefLength;
          p1    = cmp[n].shp_p1 / C.RefLength;
          p2    = cmp[n].shp_p2 / C.RefLength;
          CF.setShapeParam(cmp[n].nv, center, depth, p1, p2);
          break;

        case DARCY:
          Exit(0);
          break;
      }

      // 回転角度の計算
      CF.getAngle();


      // bboxと投影面積、インデクスの計算
      int f_st[3], f_ed[3];
      cmp[n].area = CF.getBboxArea(f_st, f_ed);


      // インデクスのサイズ登録と存在フラグ
      cmp[n].setBbox(f_st, f_ed);
      cmp[n].setEnsLocal(ON);


      double flop=0.0;

      TIMING_start("Compo_Vertex8");
      flop = 0.0;
      CF.vertex8(f_st, f_ed, d_cvf, flop);
      TIMING_stop("Compo_Vertex8", flop);

      TIMING_start("Compo_Subdivision");
      flop = 0.0;
      CF.subdivision(f_st, f_ed, d_cvf, flop);
      TIMING_stop("Compo_Subdivision", flop);

      // ########## 確認のための出力
#if 0
      REAL_TYPE org[3], pit[3];

      //  ガイドセルがある場合(GuideOut != 0)にオリジナルポイントを調整
      for (int i=0; i<3; i++)
      {
        org[i] = C.org[i] - C.dx[i]*(REAL_TYPE)C.GuideOut;
        pit[i] = C.dx[i];
      }

      // 出力ファイルの指定が有次元の場合
      if ( C.Unit.File == DIMENSIONAL )
      {
        for (int i=0; i<3; i++)
        {
          org[i] *= C.RefLength;
          pit[i] *= C.RefLength;
        }
      }
      writeRawSPH(cvf, size, guide, 0, org, pit, sizeof(REAL_TYPE));
#endif
      // ##########
    }

  } // loop - NoCmpo

}



// #################################################################
/* @brief パラメータのロードと計算領域を初期化し，並列モードを返す
 * @param [in] tpf ffvのパラメータを保持するTextParserインスタンス
 * @retval 並列モードの文字列
 */
string FFV::setDomain(TextParser* tpf)
{

  // ランク情報をセット >> 各クラスでランク情報メンバ変数を利用する前にセットすること
  C.setRankInfo    (paraMngr, procGrp);
  B.setRankInfo    (paraMngr, procGrp);
  V.setRankInfo    (paraMngr, procGrp);
  BC.setRankInfo   (paraMngr, procGrp);
  Ex->setRankInfo  (paraMngr, procGrp);
  MO.setRankInfo   (paraMngr, procGrp);
  F->setRankInfo   (paraMngr, procGrp);
  GM.setRankInfo   (paraMngr, procGrp);

  for (int i=0; i<ic_END; i++)
  {
    LS[i].setRankInfo(paraMngr, procGrp);
  }


  // 並列モードの取得
  string str = setParallelism();


  // 最初のパラメータの取得 >> C.guide
  C.get1stParameter(&DT);


  // 代表パラメータをコピー
  Ex->setRefParameter(&C);


  // 例題クラス固有のパラメータを取得
  if ( !Ex->getTP(&C, tpf) ) Exit(0);


  // 領域分割パラメータをロード : 分割指示 (divtype = 1-with / 2-without)
  // pitch[], G_region[], G_origin[]が確定する >> パラメータファイルで指定した値
  int div_type = SD_getParameter(tpf);


  // CPMlibの機能を用いて、計算領域分割情報を設定する
  SD_Initialize(div_type, tpf);


  // 各例題固有の領域パラメータを設定する
  Ex->setDomainParameter(&C, size, origin, region, pitch);


  // 各クラスで領域情報を保持
  setDomainInfo      (C.guide, C.RefLength);
  C.setDomainInfo    (C.guide, C.RefLength);
  B.setDomainInfo    (C.guide, C.RefLength);
  V.setDomainInfo    (C.guide, C.RefLength);
  BC.setDomainInfo   (C.guide, C.RefLength);
  Ex->setDomainInfo  (C.guide, C.RefLength);
  MO.setDomainInfo   (C.guide, C.RefLength);
  F->setDomainInfo   (C.guide, C.RefLength);
  GM.setDomainInfo   (C.guide, C.RefLength);


  for (int i=0; i<ic_END; i++)
  {
    LS[i].setDomainInfo(C.guide, C.RefLength);
  }

  // 境界条件のドライランの指定 >> F->getFIOparams()でパラメータを書き換えるので先にコール
  C.getDryRun();


  // ファイルIOパラメータ << get1stParameter()でgetTurbulenceModel()を呼んだあと
  F->getFIOparams();


  // 従属的なパラメータの取得
  C.get2ndParameter(&RF);

  // 係数行列の書き出しモード
  C.getAXB();

  // 粒子追跡
  C.getParticleTracking();

  // 全ノードについて，ローカルノード1面・一層あたりの通信量の和（要素数）を計算
  double c = 0.0;

  // 内部面のみをカウントする
  for (int n=0; n<6; n++)
  {
    if ( nID[n] >= 0 ) {

      switch (n)
      {
        case X_minus:
        case X_plus:
          c += (double)(size[1]*size[2]);
          break;

        case Y_minus:
        case Y_plus:
          c += (double)(size[0]*size[2]);
          break;

        case Z_minus:
        case Z_plus:
          c += (double)(size[0]*size[1]);
          break;
      }
    }
  }

  if ( numProc > 1 )
  {
    double tmp = c;
    if ( paraMngr->Allreduce(&tmp, &c, 1, MPI_SUM, procGrp) != CPM_SUCCESS ) Exit(0);
  }

  face_comm_size = c;


  return str;
}




// #################################################################
// @brief 初期条件の設定
void FFV::setInitialCondition()
{
  double flop_task;
  Gemini_R* m_buf = new Gemini_R [C.NoCompo+1];

  double tm = CurrentTime * C.Tscale;


  if ( C.Start == initial_start )
  {
		REAL_TYPE U0[3];

		// 速度の初期条件の設定
    if (C.Unit.Param == DIMENSIONAL)
    {
      U0[0] = C.iv.VecU/C.RefVelocity;
      U0[1] = C.iv.VecV/C.RefVelocity;
      U0[2] = C.iv.VecW/C.RefVelocity;
    }
    else
    {
      U0[0] = C.iv.VecU;
      U0[1] = C.iv.VecV;
      U0[2] = C.iv.VecW;
    }
		fb_set_vector_(d_v, size, &guide, U0, d_bcd);
    fb_set_fvector_(d_vf, size, &guide, U0, d_bcd);


    // LESの初期擾乱
    if (C.LES.InitialPerturbation == ON)
    {
      perturbation();
    }


    // セルフェイスの設定　発散値は関係なし
    BC.modDivergence(d_dv, d_cdf, CurrentTime, &C, v00, m_buf, flop_task);



		// 外部境界面の移流速度を計算し，外部境界条件を設定
		BC.OuterVBC(d_v, d_vf, d_cdf, tm, &C, v00, ensPeriodic);
    BC.InnerVBCperiodic(d_v, d_bcd);


		// 圧力
    REAL_TYPE ip;
    if (C.Unit.Param == DIMENSIONAL)
    {
      ip = FBUtility::convPrsD2ND(C.iv.Pressure, C.BasePrs, C.RefDensity, C.RefVelocity, C.Unit.Prs);
    }
    else
    {
      ip = C.iv.Pressure;
    }

    U.initS3D(d_p, size, guide, ip);
		BC.OuterPBC(d_p, ensPeriodic);

		// 温度　コンポーネントの初期値
		if ( C.isHeatProblem() )
    {
      for (int m=1; m<=C.NoCompo; m++) // Mediumでまわすと，オーダーがエンコードされていない場合もある
      {
        BC.setInitialTempCompo(m, d_bcd, d_ie);
      }

			BC.OuterTBCperiodic(d_ie, ensPeriodic);
		}

  }
  else // リスタート時
  {
    // 内部境界条件
    BC.InnerVBCperiodic(d_v, d_bcd);
    BC.InnerPBCperiodic(d_p, d_bcd);

    // 外部境界条件
    BC.OuterVBC(d_v, d_vf, d_cdf, tm, &C, v00, ensPeriodic);

    // 流出境界の流出速度の算出
    BC.modDivergence(d_ws, d_cdf, CurrentTime, &C, v00, m_buf, flop_task);

    //if ( C.isHeatProblem() ) BC.InnerTBC_Periodic()

  }



  // 外部境界面の流出流量と移流速度
  DomainMonitor( BC.exportOBC(), &C);



  // 初期解およびリスタート解の同期
  if ( numProc > 1 )
  {
    if ( paraMngr->BndCommV3D(d_v,  size[0], size[1], size[2], guide, guide, procGrp) != CPM_SUCCESS ) Exit(0);
    if ( paraMngr->BndCommV3D(d_vf, size[0], size[1], size[2], guide, guide, procGrp) != CPM_SUCCESS ) Exit(0);
    if ( paraMngr->BndCommS3D(d_p,  size[0], size[1], size[2], guide, guide, procGrp) != CPM_SUCCESS ) Exit(0);

    if ( C.isHeatProblem() )
    {
      if ( paraMngr->BndCommS3D(d_p, size[0], size[1], size[2], guide, guide, procGrp) != CPM_SUCCESS ) Exit(0);
    }
  }

  // VOF
  if ( C.BasicEqs == INCMP_2PHASE )
  {
    setVOF();

    if ( numProc > 1 )
    {
      if ( paraMngr->BndCommS3D(d_vof, size[0], size[1], size[2], guide, guide, procGrp) != CPM_SUCCESS ) Exit(0);
    }
  }

  // 後始末
  if ( m_buf ) { delete [] m_buf; m_buf=NULL; }
}



// #################################################################
/* @brief ParseMatクラスをセットアップし，媒質情報を入力ファイルから読み込み，媒質リストを作成する
 * @param [in] fp  ファイルポインタ
 */
void FFV::setMediumList(FILE* fp)
{
  Hostonly_
  {
    fprintf(fp,"\n----------\n\n");
    fprintf(fp,"\n\t>> Medium List\n\n");
  }


  // 媒質情報をmat[]にロード
  M.getMediumTable(C.NoMedium, mat);


  if ( !M.check(mat) )
  {
    Hostonly_ stamped_printf("Error : Duplicate label in Material Table\n");
    Exit(0);
  }

  // 媒質テーブルの表示
  Hostonly_
  {
    M.printMatList(fp, mat);
  }

}




// #################################################################
/* @brief 各種例題のモデルをセットアップ
 * @param [in] PrepMemory  前処理に必要なメモリ
 * @param [in] TotalMemory ソルバー実行に必要なメモリ
 * @param [in] fp          ファイルポインタ
 */
void FFV::setModel(double& PrepMemory, double& TotalMemory, FILE* fp)
{
  // 内部領域の交点と境界ID
  switch (C.Mode.Example)
  {
    case id_Polygon: // ユーザ例題
      SM_Polygon2Cut(PrepMemory, TotalMemory, fp);
      break;

    case id_Sphere:
    case id_Step:
    case id_Cylinder:
    case id_Duct:
      Ex->setup(d_bcd, &C, C.NoCompo, mat, C.NoCompo, cmp, d_cut, d_bid);
      break;

    default: // ほかのIntrinsic problems
      break;
  }



  /* 外部境界方向にカットがあるセルには、ガイドセルをCutIDの媒質でペイント
  unsigned long painted[6];

  V.paintCutIDonGC(d_bcd, d_bid, painted, cmp);

  Hostonly_ {
    fprintf(fp, "\n\tPainted guide cell by cut\n");
    fprintf(fp, "\t\t X minus = %10ld\n", painted[0]);
    fprintf(fp, "\t\t X plus  = %10ld\n", painted[1]);
    fprintf(fp, "\t\t Y minus = %10ld\n", painted[2]);
    fprintf(fp, "\t\t Y plus  = %10ld\n", painted[3]);
    fprintf(fp, "\t\t Z minus = %10ld\n", painted[4]);
    fprintf(fp, "\t\t Z plus  = %10ld\n\n", painted[5]);
  }
  */


  // 外部境界面の処理　ここで外部境界面の交点距離と交点ID、媒質IDをセット

  for (int face=0; face<NOFACE; face++)
  {
    if( nID[face] >= 0 ) continue;

    BoundaryOuter* m_obc = BC.exportOBC(face);
    int id = m_obc->getGuideMedium();
    int cls= m_obc->getClass();
    int ptr_cmp = m_obc->getPtr2cmp();

    // 周期境界以外
    switch (cls)
    {
      case OBC_SYMMETRIC: // >> 対称条件は流体
        if (mat[id].getState() != FLUID)
        {
          Hostonly_ printf("Specified medium in '%s' is not FLUID or not listed.\n", m_obc->alias.c_str());
          Exit(0);
        }
        V.setOBC(face, id, ptr_cmp, "fluid", d_bcd, d_cut, d_bid);
        break;

      case OBC_WALL:
        if (mat[id].getState() != SOLID)
        {
          Hostonly_ printf("Specified medium in '%s' is not SOLID or not listed.\n", m_obc->alias.c_str());
          Exit(0);
        }
        V.setOBC(face, id, ptr_cmp, "solid", d_bcd, d_cut, d_bid);
        break;

      case OBC_INTRINSIC:
        Ex->setOBC(face, d_bcd, &C, G_origin, C.NoCompo, mat, d_cut, d_bid);
        break;

      case OBC_PERIODIC: // nothing
        break;

      default:
        if (mat[id].getState() != FLUID) {
          Exit(0);
        }
        V.setOBC(face, id, ptr_cmp, "fluid", d_bcd, d_cut, d_bid);
        break;
    }
  }


  // 周期境界の場合の外部境界面処理
  V.setOBCperiodic(d_bcd, ensPeriodic);


  // ガイドセル同期
  if ( numProc > 1 )
  {
    if ( paraMngr->BndCommS3D(d_bcd, size[0], size[1], size[2], guide, 1, procGrp) != CPM_SUCCESS ) Exit(0);
    if ( paraMngr->BndCommS3D(d_bid, size[0], size[1], size[2], guide, 1, procGrp) != CPM_SUCCESS ) Exit(0);
    if ( paraMngr->BndCommS3D(d_cut, size[0], size[1], size[2], guide, 1, procGrp) != CPM_SUCCESS ) Exit(0);
  }

}




// #################################################################
/**
 * @brief MonitorListの処理
 */
void FFV::setMonitorList()
{
  MO.setControlVars(d_bid,
                    d_cut,
                    d_bcd,
                    C.RefVelocity,
                    C.BaseTemp,
                    C.DiffTemp,
                    C.RefDensity,
                    C.RefLength,
                    C.BasePrs,
                    C.RefLength,
                    C.Mode.Precision,
                    C.Unit.Prs,
                    C.num_process,
                    C.NoCompo,
                    mat_tbl);

  // パラメータを取得し，サンプリングの座標値をセットする >> サンプリング指定がなければ，return
  if ( !MO.getMonitor(&C, cmp) ) return;

  // モニター指定変数とKOSの整合性チェック
  if ( !MO.checkConsistency(C.KindOfSolver) )  Exit(0);

  // Polygon指定のセルモニターの場合に交点と境界IDを除去
  MO.clearCut();

  // モニタ結果出力ファイル群のオープン
  MO.openFile();


  // ########## 確認のための出力
  //Ex->writeSVX(d_bcd, &C);
#if 0
  REAL_TYPE org[3], pit[3];

  //  ガイドセルがある場合(GuideOut != 0)にオリジナルポイントを調整
  for (int i=0; i<3; i++)
  {
    org[i] = C.org[i] - C.dx[i]*(REAL_TYPE)C.GuideOut;
    pit[i] = C.dx[i];
  }

  // 出力ファイルの指定が有次元の場合
  if ( C.Unit.File == DIMENSIONAL )
  {
    for (int i=0; i<3; i++)
    {
      org[i] *= C.RefLength;
      pit[i] *= C.RefLength;
    }
  }
  writeRawSPH(d_bcd, size, guide, 0, org, pit, sizeof(REAL_TYPE));
#endif
  // ##########
}


// #################################################################
/* @brief 並列化と分割の方法を保持
 * @return 並列モード
 */
string FFV::setParallelism()
{
  string para_mode;

  C.num_thread  = omp_get_max_threads();

  // Serial or Parallel environment
  if( paraMngr->IsParallel() )
  {
    C.num_process = paraMngr->GetNumRank(procGrp);

    if ( C.num_thread > 1 )
    {
      C.Parallelism = Control::Hybrid;
      para_mode = "Hybrid";
    }
    else
    {
      C.Parallelism = Control::FlatMPI;
      para_mode = "FlatMPI";
    }
  }
  else
  {
    C.num_process = 1;

    if ( C.num_thread > 1 )
    {
      C.Parallelism = Control::OpenMP;
      para_mode = "OpenMP";
    }
    else
    {
      C.Parallelism = Control::Serial;
      para_mode = "Serial";
    }
  }
  return para_mode;
}


// #################################################################
/* @brief 時間積分幅や物理パラメータの設定
 */
void FFV::setParameters()
{
  // 無次元数などの計算パラメータを設定する．MediumListを決定した後，かつ，SetBC3Dクラスの初期化前に実施すること
  // 代表物性値をRefMatの示す媒質から取得
  // Δt=constとして，無次元の時間積分幅 deltaTを計算する

  double min_dx = std::min(pitch[0], std::min(pitch[1], pitch[2]));
  DT.set_Vars(C.KindOfSolver, C.Unit.Param, min_dx, (double)C.Reynolds, (double)C.Peclet);


  // 無次元速度1.0を与えてdeltaTをセットし，エラーチェック
  switch ( DT.set_DT(1.0) )
  {
    case 0: // 成功
      break;

    case 1:
      Hostonly_ stamped_printf("\tdt selection error(1) : 'Kind of Solver' is solid conduction. Consider to specify other dt scheme or confirm 'Kind of Solver'.\n");
      Exit(0);
      break;

    case 2:
      Hostonly_ stamped_printf("\tdt selection error(2) : 'Kind of Solver' is solid conduction. Consider to specify other dt scheme or confirm 'Kind of Solver'.\n");
      Exit(0);
      break;

    case 3:
      Hostonly_ stamped_printf("\tdt selection error(3) : 'Kind of Solver' includes flow effect. Consider to specify other dt scheme or confirm 'Kind of Solver'.\n");
      Exit(0);
      break;

    case 4:
      Hostonly_ stamped_printf("\tdt selection error(4) : 'Kind of Solver' is solid conduction. Consider to specify other dt scheme or confirm 'Kind of Solver'.\n");
      Exit(0);
      break;

    case 5:
      Hostonly_ stamped_printf("\tdt selection error(5) : 'Kind of Solver' is solid conduction. Consider to specify other dt scheme or confirm 'Kind of Solver'.\n");
      Exit(0);
      break;

    case 6:
      Hostonly_ stamped_printf("\tdt selection error(6) : CFL max V and Diffusion ; Not implemented yet\n");
      Exit(0);
      break;

    case 7:
      Hostonly_ stamped_printf("\tdt selection error(7) : CFL max V for cp ;  Not implemented yet\n");
      Exit(0);
      break;

    default:
      Exit(0);
      break;
  }


  // 無次元時間積分幅
  deltaT = DT.get_DT();


  // コンポーネントと外部境界のパラメータを有次元化
  C.setCmpParameters(mat, cmp, BC.exportOBC());
}



// #################################################################
// VOF値を気体(0.0)と液体(1.0)で初期化
void FFV::setVOF()
{
  size_t m;
  int s, odr;

  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;

  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        s = d_bcd[m];
        odr = DECODE_CMP(s);

        if ( cmp[odr].getState() == FLUID )
        {
          d_vof[m] = ( cmp[odr].getPhase() == GAS ) ? 0.0 : 1.0;
        }
      }
    }
  }
}




// #################################################################
/* @brief グローバルな領域情報を取得し、無次元の領域基本パラメータを返す
 * @param [in] tp_dom   TextParserクラス
 * @return 分割指示 (1-with / 2-without)
 * @note この関数では、pitch[], G_region[], G_origin[]を確定する
 */
int FFV::SD_getParameter(TextParser* tp_dom)
{
  // 領域分割モードのパターン
  //     分割指定(G_div指定)         domain.txt
  // 1)  G_div指定なし |  G_orign + G_region + (G_pitch || G_voxel)
  // 2)  G_div指定あり |  G_orign + G_region + (G_pitch || G_voxel)
  // 3)  G_div指定なし |          + ActiveDomainInfo
  // 4)  G_div指定あり |          + ActiveDomainInfo

  /* Policy
   * G_regionは必須
   * G_voxelとG_pitchでは，G_voxelが優先．両方が指定されている場合はエラー
   * G_pitchが指定されており，割り切れない場合は，領域を拡大
   *        G_voxel = (int)ceil(G_region/G_pitch)
   *        G_region = G_pitch * G_voxel
   *
   * CPMlibに渡す場合にG_voxelモードにする
   */

  string label, str;
  int div_type = 1; // 指定分割 => 1


  // 長さの単位
  label = "/DomainInfo/UnitOfLength";

  if ( !tp_dom->getInspectedValue(label, str) )
  {
    Hostonly_ stamped_printf("\tParsing error : Invalid string for '%s'\n", label.c_str());
    Exit(0);
  }

  if     ( !strcasecmp(str.c_str(), "NonDimensional") )  C.Unit.Length = LTH_ND;
  else if( !strcasecmp(str.c_str(), "M") )               C.Unit.Length = LTH_m;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described at '%s'\n", label.c_str());
    Exit(0);
  }


  // G_origin　必須
  label = "/DomainInfo/GlobalOrigin";

  if ( !tp_dom->getInspectedVector(label, G_origin, 3) )
  {
    cout << "ERROR : in parsing [" << label << "]" << endl;
    Exit(0);
  }


  // G_region 必須
  label = "/DomainInfo/GlobalRegion";

  if ( !tp_dom->getInspectedVector(label, G_region, 3) )
  {
    Hostonly_ cout << "ERROR : in parsing [" << label << "]" << endl;
    Exit(0);
  }

  if ( (G_region[0]>0.0) && (G_region[1]>0.0) && (G_region[2]>0.0) )
  {
    ; // skip
  }
  else
  {
    Hostonly_ cout << "ERROR : in parsing [" << label << "]" << endl;
    Exit(0);
  }

  // 排他チェックフラグ
  bool g_flag = true;
  bool p_flag = true;


  // G_voxel
  label = "/DomainInfo/GlobalVoxel";

  if ( !tp_dom->getInspectedVector(label, G_size, 3) ) g_flag = false;



  // pitch
  label = "/DomainInfo/GlobalPitch";

  if ( !tp_dom->getInspectedVector(label, pitch, 3) ) p_flag = false;


  // 排他チェック
  if ( !g_flag && !p_flag ) // 両方とも指定されていない
  {
    Hostonly_ cout << "\tNeither GlobalVoxel nor GlobalPitch is specified. You need to choose either parameter.\n" << endl;
    Exit(0);
  }
  else if ( g_flag ) // G_voxelの指定がある >> G_voxelが優先
  {
    if ( g_flag )
    {
      if ( (G_size[0]>0) && (G_size[1]>0) && (G_size[2]>0) )
      {
        pitch[0] = G_region[0] / G_size[0];
        pitch[1] = G_region[1] / G_size[1];
        pitch[2] = G_region[2] / G_size[2];

        if ( Ex->mode == Intrinsic::dim_2d )
        {
          pitch[2] = pitch[0];
          G_region[2] = (REAL_TYPE)G_size[2] * pitch[2];
          G_origin[2] = -0.5 * pitch[2];
        }

      }
      else
      {
        Hostonly_ printf("ERROR : in parsing [%s] >> (%d, %d, %d)\n", label.c_str(), G_size[0], G_size[1], G_size[2] );
        Exit(0);
      }
    }
  }
  else if ( !g_flag && p_flag ) // pitchのみ
  {
    if ( (pitch[0]>0.0) && (pitch[1]>0.0) && (pitch[2]>0.0) )
    {

      // pitchを基準にして、全計算領域の分割数を計算する。四捨五入
      G_size[0] = (int)floor(G_region[0]/pitch[0] + 0.5);
      G_size[1] = (int)floor(G_region[1]/pitch[1] + 0.5);
      G_size[2] = (int)floor(G_region[2]/pitch[2] + 0.5);


      // 再計算された全計算領域の分割数から、全計算領域の大きさを再計算
      double gr[3];
      gr[0] = (double)G_size[0] * (double)pitch[0];
      gr[1] = (double)G_size[1] * (double)pitch[1];
      gr[2] = (double)G_size[2] * (double)pitch[2];

      // 整合性チェック
      if ( Ex->mode == Intrinsic::dim_3d )
      {
        if ( (G_region[0]-gr[0]>ROUND_EPS) || (G_region[1]-gr[1]>ROUND_EPS) || (G_region[2]-gr[2]>ROUND_EPS) )
        {
          Hostonly_ {
            printf("\tGlobal Region is modified due to maintain the consistency of domain parameters.\n\n");
            printf("\t[%12.6e  %12.6e  %12.6e] >> [%12.6e  %12.6e  %12.6e]\n\n",
                   G_region[0], G_region[1], G_region[2], gr[0], gr[1], gr[2]);
          }
          G_region[0] = gr[0];
          G_region[1] = gr[1];
          G_region[2] = gr[2];
        }
      }
      else if ( Ex->mode == Intrinsic::dim_2d )
      {
        if ( (G_region[0]-gr[0]>ROUND_EPS) || (G_region[1]-gr[1]>ROUND_EPS) )
        {
          Hostonly_ {
            printf("\tGlobal Region is modified due to maintain the consistency of domain parameters.\n\n");
            printf("\t[%12.6e  %12.6e  %12.6e] >> [%12.6e  %12.6e  %12.6e]\n\n",
                   G_region[0], G_region[1], G_region[2], gr[0], gr[1], gr[2]);
          }
          G_region[0] = gr[0];
          G_region[1] = gr[1];
          G_region[2] = gr[2];
          G_origin[2] = -0.5 * pitch[2];
        }
      }
      else
      {
        Exit(0);
      }

    }
    else // パラメータが無効の場合
    {
      Hostonly_ printf("ERROR : in parsing [%s] >> (%e, %e, %e)\n", label.c_str(), pitch[0], pitch[1], pitch[2] );
      Exit(0);
    }
  }
  else
  {
    Hostonly_ cout << "\tError : unknown\n" << endl;
    Exit(0);
  }




  // 有次元の場合に無次元化する　paraMngr->VoxelInit()で計算する基本量
  if (C.Unit.Param == DIMENSIONAL )
  {
    // 無次元化
    for (int i=0; i<3; i++) {
      pitch[i]    /= C.RefLength;
      G_origin[i] /= C.RefLength;
      G_region[i] /= C.RefLength;
    }
  }



  // 2D check
  if ( (Ex->mode == Intrinsic::dim_2d) && (G_size[2] != 1) )
  {
    Hostonly_ {
      printf("\tError : In case of 2 dimensional problem, kmax must be 1 >> %d.\n", G_size[2]);
      Exit(0);
    }
  }


  // 偶数のチェック
  if ( Ex->even == ON )
  {
    if ( G_size[0]/2*2 != G_size[0] )
    {
      printf("\tDimension size must be even for x direction (%d %d %d)\n", G_size[0], G_size[1], G_size[2]);
      Exit(0);
    }
    if ( G_size[1]/2*2 != G_size[1] )
    {
      printf("\tDimension size must be even for y direction (%d %d %d)\n", G_size[0], G_size[1], G_size[2]);
      Exit(0);
    }
    if ( (Ex->mode == Intrinsic::dim_3d) && (G_size[2]/2*2 != G_size[2]) )
    {
      printf("\tDimension size must be even for z direction (%d %d %d)\n", G_size[0], G_size[1], G_size[2]);
      Exit(0);
    }
  }


  // G_division オプション
  label = "/DomainInfo/GlobalDivision";

  if ( !tp_dom->getInspectedVector(label, G_division, 3) )
  {
    Hostonly_ cout << "\tAutomatic domain division is selected." << endl;
    div_type = 2; // 自動分割
  }

  // プロセス分割数が指定されている場合のチェック
  if ( div_type == 1 )
  {
    if ( (G_division[0]>0) && (G_division[1]>0) && (G_division[2]>0) )
    {
      Hostonly_ printf("\tManual domain division is selected.\n");
    }
    else
    {
      Hostonly_ cout << "ERROR : in parsing [" << label << "]" << endl;
      Exit(0);
    }
  }



  // ActiveSubdomainファイル名の取得 >> ファイル名が指定されている場合はASモード
  label = "/DomainInfo/ActiveSubDomainFile";
  if ( tp_dom->chkLabel(label) )
  {
    if ( !tp_dom->getInspectedValue(label, str ) )
    {
      Hostonly_ cout << "\tNo option : in parsing [" << label << "]" << endl;
      Exit(0);
    }
    else
    {
      if ( str.empty() == true )
      {
        EXEC_MODE = ffvc_solver;
      }
      else
      {
        active_fname = str;
        EXEC_MODE = ffvc_solverAS;
        Hostonly_ printf("\n\tActive subdomain mode\n");
      }
    }
  }


  return div_type;
}


// #################################################################
/* @brief 計算領域情報を設定する
 * @param [in] div_type 分割数モード
 * @param [in] tp_dom   TextParserクラス
 */
void FFV::SD_Initialize(const int div_type, TextParser* tp_dom)
{
  // 袖通信の最大数
  size_t Nvc  = (size_t)C.guide;
  size_t Ncmp = 3; // 最大はベクトル3成分

  int m_sz[3]  = {G_size[0], G_size[1], G_size[2]};
  int m_div[3] = {G_division[0], G_division[1], G_division[2]};
  double m_org[3] = {(double)G_origin[0], (double)G_origin[1], (double)G_origin[2]};
  double m_reg[3] = {(double)G_region[0], (double)G_region[1], (double)G_region[2]};



  // 領域分割モードのパターン
  //     分割指定(G_div指定)         domain.txt
  // 1)  G_div指定なし |  G_orign + G_region + (pitch || G_voxel)
  // 2)  G_div指定あり |  G_orign + G_region + (pitch || G_voxel)
  // 3)  G_div指定なし |          + ActiveDomainInfo
  // 4)  G_div指定あり |          + ActiveDomainInfo

  /* Policy
   * G_regionは必須
   * G_voxelとpitchでは，G_voxelが優先．両方が指定されている場合はエラー
   * G_pitchが指定されており，割り切れない場合は，領域を拡大
   *        G_voxel = (int)ceil(G_region/pitch)
   *        G_region = pitch * G_voxel
   */


  // ノーマルモードとActive subdomainモードの切り替え
  if ( EXEC_MODE == ffvc_solver )
  {
    // 分割数を元に分割する >> CPMlibの仕様
    switch (div_type)
    {
      case 1: // 分割数が指示されている場合
        if ( paraMngr->VoxelInit(m_div, m_sz, m_org, m_reg, Nvc, Ncmp, DIV_COMM_SIZE, procGrp) != CPM_SUCCESS )
        {
          cout << "Domain decomposition error : " << endl;
          Exit(0);
        }
        break;


      case 2: // 分割数が指示されていない場合
        if ( paraMngr->VoxelInit(m_sz, m_org, m_reg, Nvc, Ncmp, DIV_COMM_SIZE, procGrp) != CPM_SUCCESS )
        {
          cout << "Domain decomposition error : " << endl;
          Exit(0);
        }
        break;


      default:
        Exit(0);
        break;
    }
  }
  else if ( EXEC_MODE == ffvc_solverAS )
  {
    if ( paraMngr->VoxelInit_Subdomain(m_div, m_sz, m_org, m_reg, active_fname, Nvc, Ncmp, procGrp) != CPM_SUCCESS )
    {
      cout << "Domain decomposition error : " << endl;
      Exit(0);
    }
  }
  else
  {
    Exit(0);
  }



  // アドレスサイズのチェック
  size_t n_cell[3];

  for ( int i=0; i<3; i++)
  {
    n_cell[i] = (size_t)(size[i] + 2*guide);  // 分割数+ガイドセル両側
  }

  size_t size_n_cell = n_cell[0] * n_cell[1] * n_cell[2];

  if (size_n_cell*6  > UINT_MAX)
  {
    Hostonly_
    {
      stamped_printf("\n\tError : Product of size[]*6 exceeds UINT_MAX\n\n");
    }
    Exit(0);
  }

}


// #################################################################
/* @brief 境界条件ポリゴンの正しい面積を保持する
 * @param [in]     fp       ファイルポインタ
 * @note マスターのみで実行するメソッド
 */
void FFV::SM_getVspecArea(FILE* fp)
{

  // ポリゴン情報へのアクセス
  vector<PolygonGroup*>* pg_roots = PL->get_root_groups();

  // Polygon Groupの数
  C.num_of_polygrp = pg_roots->size();

  if ( (C.num_of_polygrp < 1) || (C.num_of_polygrp > C.NoCompo) )
  {
    printf (   "\tError : Number of polygon group must be greater than 1 and less than NoCompo.\n");
    fprintf(fp,"\tError : Number of polygon group must be greater than 1 and less than NoCompo.\n");
    Exit(0);
  }

  fprintf(fp, "\t<< MASTER PROCESS ONLY >>\n");
  fprintf(fp, "\t   Polygon Group Label       Medium Alias              Local BC     Polygons          Area\n");
  fprintf(fp, "\t   ---------------------------------------------------------------------------------------\n");


  vector<PolygonGroup*>::iterator it;

  // ポリゴンのグループラベルを照合し、cmp[]の格納インデクスから境界条件を割り出す
  for (it = pg_roots->begin(); it != pg_roots->end(); it++)
  {
    std::string m_pg = (*it)->get_name();     // グループラベル
    std::string m_mat = (*it)->get_label();   // 媒質ラベル
    std::string m_bc = (*it)->get_type();     // 境界条件ラベル
    int ntria= (*it)->get_group_num_tria();   // ローカルのポリゴン数
    REAL_TYPE area = (*it)->get_group_area(); // ローカルのポリゴン面積

    // cmp[]の格納順を探す
    int c_id;

    for (int i=1; i<=C.NoCompo; i++)
    {
      if ( FBUtility::compare(m_pg, cmp[i].alias) )
      {
        c_id = i;
        break;
      }
    }

    int typ = cmp[c_id].getType();

    if ( typ==SPEC_VEL || typ==OUTFLOW )
    {
      cmp[c_id].area = area;
      fprintf(fp,"\t  %20s %18s  %20s %12d  %e\n",
              m_pg.c_str(),
              m_mat.c_str(),
              m_bc.c_str(),
              ntria,
              area);
    }
  }

  delete pg_roots;

  fprintf(fp, "\n\n");
}




// #################################################################
/* @brief 幾何形状情報を準備し，交点計算を行う
 * @param [in,out] m_prep   前処理用のメモリサイズ
 * @param [in,out] m_total  本計算用のメモリリサイズ
 * @param [in]     fp       ファイルポインタ
 */
void FFV::SM_Polygon2Cut(double& m_prep, double& m_total, FILE* fp)
{
  TIMING_start("Polylib_Section");


  Hostonly_
  {
    fprintf(fp,"\n----------\n\n");
    fprintf(fp,"\t>> Polylib configuration\n\n");
  }

  // Polylibファイルをテンポラリに出力
  if ( !F->writePolylibFile(cmp) )
  {
    Hostonly_
    {
      fprintf(fp,"\tError : writing polylib.tp\n");
      printf    ("\tError : writing polylib.tp\n");
    }
    Exit(0);
  }


  // Polylib: インスタンス取得
  PL = MPIPolylib::get_instance();


  // Polylib: 並列計算領域情報を設定 >> 有次元
  unsigned  poly_gc[3] = {guide, guide, guide};


  // 読み込むデータはオリジナルのまま >> Polylib管理のメソッドで検索などを行うため
  poly_stat = PL->init_parallel_info(paraMngr->GetMPI_Comm(procGrp),
                                     originD,           // 自ランクの基点座標
                                     (unsigned*)size,   // 自ランクの分割数
                                     poly_gc,           // ガイドセル数
                                     pitchD             // 格子幅
                                     );
  if ( poly_stat != PLSTAT_OK )
  {
    Hostonly_
    {
      fprintf(fp,"\tRank [%6d]: p_polylib->init_parallel_info() failed.\n", myRank);
      printf    ("\tRank [%6d]: p_polylib->init_parallel_info() failed.\n", myRank);
    }
    Exit(0);
  }



  // Polylib: STLデータ読み込み
  TIMING_start("Loading_Polygon_File");

  // ロード
  //poly_stat = PL->load_rank0( "polylib.tp");
  poly_stat = PL->load_only_in_rank0( "polylib.tp");

  if( poly_stat != PLSTAT_OK )
  {
    Hostonly_
    {
      printf    ("\tRank [%6d]: p_polylib->load_only_in_rank0() failed.", myRank);
      fprintf(fp,"\tRank [%6d]: p_polylib->load_only_in_rank0() failed.", myRank);
    }
    Exit(0);
  }

  TIMING_stop("Loading_Polygon_File");



  // 内部の速度境界条件を与えるポリゴンの面積をマスターランクのみで正しく求める
  TIMING_start("Get_Accurate_BC_Area");
  Hostonly_
  {
    SM_getVspecArea(fp);
  }
  TIMING_stop("Get_Accurate_BC_Area");


  // 各ランクに分配
  TIMING_start("Distributing_Polygon");

  poly_stat = PL->distribute_only_from_rank0();

  if( poly_stat != PLSTAT_OK )
  {
    Hostonly_
    {
      printf    ("\tRank [%6d]: p_polylib->distribute_only_from_rank0() failed.", myRank);
      fprintf(fp,"\tRank [%6d]: p_polylib->distribute_only_from_rank0() failed.", myRank);
    }
    Exit(0);
  }

  TIMING_stop("Distributing_Polygon");



  // 階層情報表示 debug brief hierarchy
  // ##########
#if 0
  PL->show_group_hierarchy();
  PL->show_group_hierarchy(fp);
#endif
  // ##########


  // ポリゴン情報へのアクセス
  vector<PolygonGroup*>* pg_roots = PL->get_root_groups();


  // Polygon Groupの数
  C.num_of_polygrp = pg_roots->size();

  if ( (C.num_of_polygrp < 1) || (C.num_of_polygrp > C.NoCompo) )  Exit(0);


  // PolygonPropertyの配列
  PG = new PolygonProperty[C.num_of_polygrp];

  Hostonly_ {
    fprintf(fp, "\t<< ALL PROCESSES >>\n");
    fprintf(fp, "\t   Polygon Group Label       Medium Alias              Local BC     Polygons          Area\n");
    fprintf(fp, "\t   ---------------------------------------------------------------------------------------\n");
  }

  vector<PolygonGroup*>::iterator it;

  // ポリゴンのグループラベルを照合し、cmp[]の格納インデクスをIDとして割り当てる
  int c = 0;
  for (it = pg_roots->begin(); it != pg_roots->end(); it++)
  {
    std::string m_pg = (*it)->get_name();     // グループラベル
    std::string m_mat = (*it)->get_label();   // 媒質ラベル
    std::string m_bc = (*it)->get_type();     // 境界条件ラベル
    int ntria= (*it)->get_group_num_tria();   // ローカルのポリゴン数
    REAL_TYPE area = (*it)->get_group_area(); // ローカルのポリゴン面積

    // 各ランクの保持するポリゴン数を設定
    PG[c].setLntria(ntria);


    // cmp[]の格納順を探す
    int c_id;

    for (int i=1; i<=C.NoCompo; i++)
    {
      if ( FBUtility::compare(m_pg, cmp[i].alias) )
      {
        c_id = i;
        break;
      }
    }

    // PolygonにIDを割り当てる
    poly_stat = (*it)->set_all_exid_of_trias(c_id);

    if ( poly_stat != PLSTAT_OK )
    {
      Hostonly_
      {
        printf(     "\tError : Polylib::set_all_exid_of_trias()\n");
        fprintf(fp, "\tError : Polylib::set_all_exid_of_trias()\n");
        Exit(0);
      }
    }


    if ( numProc > 1 )
    {
      int tmp = ntria;
      if ( paraMngr->Allreduce(&tmp, &ntria, 1, MPI_SUM, procGrp) != CPM_SUCCESS ) Exit(0);

      REAL_TYPE ta = area;
      if ( paraMngr->Allreduce(&ta, &area, 1, MPI_SUM, procGrp) != CPM_SUCCESS ) Exit(0);
    }


    // SPEC_VEL, OUTFLOW以外 >> SM_getVspecArea()
    int typ = cmp[c_id].getType();
    if ( typ!=SPEC_VEL && typ!=OUTFLOW )  cmp[c_id].area = area;


    PG[c].setID(c_id);
    PG[c].setGroup(m_pg);
    PG[c].setBClabel(m_bc);
    PG[c].setMaterial(m_mat);
    PG[c].setGntria(ntria);


    Hostonly_
    {
      fprintf(fp,"\t  %20s %18s  %20s %12d  %e\n",
              m_pg.c_str(),
              m_mat.c_str(),
              m_bc.c_str(),
              ntria,
              cmp[c_id].area);
    }

    c++;

    // ########## show corrdinates and area
#if 0
    PL->show_group_info(m_pg); //debug
#endif
    // ##########

  }


  delete pg_roots;

  Hostonly_
  {
    fprintf(fp, "\n");
  }



  // 使用メモリ量　基本クラスのみ
  double poly_mem, G_poly_mem;
  G_poly_mem = poly_mem = (double)PL->used_memory_size();
  m_prep += poly_mem;
  m_total+= poly_mem;

  displayMemoryInfo(fp, G_poly_mem, poly_mem, "Polygon");


  // ポリゴンのBbox
  GM.calcBboxFromPolygonGroup(PL, PG, C.num_of_polygrp);


  // Triangle display >> Debug
  // ##########
#if 0
  Vec3r m_min, m_max, t1(origin), t2(pitch), t3;
  t3.assign((REAL_TYPE)size[0]*t2.x, (REAL_TYPE)size[1]*t2.y, (REAL_TYPE)size[2]*t2.z);
  m_min = t1 - t2;      // 1層外側まで
  m_max = t1 + t3 + t2; //
  printf("min : %f %f %f\n", m_min.x, m_min.y, m_min.z);
  printf("max : %f %f %f\n", m_max.x, m_max.y, m_max.z);
  vector<Triangle*>* trias = PL->search_polygons("Ducky", m_min, m_max, false); // false; ポリゴンが一部でもかかる場合

  //Vec3r *p, nrl, n;
  Vec3r n;
  Vertex** p;
  c=0;
  vector<Triangle*>::iterator it3;
  for (it3 = trias->begin(); it3 != trias->end(); it3++) {
    p = (*it3)->get_vertex();
    n = (*it3)->get_normal();
    printf("%d : p0=(%6.3e %6.3e %6.3e)  p1=(%6.3e %6.3e %6.3e) p2=(%6.3e %6.3e %6.3e) n=(%6.3e %6.3e %6.3e)\n", c++,
           (*(p[0]))[0], (*(p[0]))[1], (*(p[0]))[2],
           (*(p[1]))[0], (*(p[1]))[1], (*(p[1]))[2],
           (*(p[2]))[0], (*(p[2]))[1], (*(p[2]))[2],
           n.x, n.y, n.z);
  }

  delete trias;  //後始末
#endif


  // ##########


  // Polylib: STLデータ書き出し
  // DomainInfo/UnitOfLength={"mm", "cm"}のとき，"m"でスケーリングされたポリゴンが出力される
  if ( C.Hide.GeomOutput == ON )
  {
    TIMING_start("Write_Polygon_File");

    unsigned poly_out_para = IO_DISTRIBUTE;
    if ( (C.Parallelism == Control::Serial) || (C.Parallelism == Control::OpenMP) ) poly_out_para = IO_GATHER;

    string fname;

    if ( poly_out_para == IO_GATHER )
    {
      poly_stat = PL->save_rank0( &fname, "stl_b" );

      if ( poly_stat != PLSTAT_OK )
      {
        Hostonly_
        {
          printf(    "Rank [%d]: p_polylib->save_rank0() failed to write into '%s'.", myRank, fname.c_str());
          fprintf(fp,"Rank [%d]: p_polylib->save_rank0() failed to write into '%s'.", myRank, fname.c_str());
        }
        Exit(0);
      }
    }
    else
    {
      poly_stat = PL->save_parallel( &fname, "stl_b" );

      if ( poly_stat != PLSTAT_OK )
      {
        Hostonly_
        {
          printf(    "Rank [%d]: p_polylib->save_parallel() failed to write into '%s'.", myRank, fname.c_str());
          fprintf(fp,"Rank [%d]: p_polylib->save_parallel() failed to write into '%s'.", myRank, fname.c_str());
        }
        Exit(0);
      }
    }

    TIMING_stop("Write_Polygon_File");
  }

  TIMING_stop("Polylib_Section");





  TIMING_start("Cut_Section");

  Hostonly_
  {
    fprintf(fp,"\n----------\n\n");
    fprintf(fp,"\t>> Calculate cut and quantize\n\n");
  }


  // 交点計算
  TIMING_start("Cut_Information");
  GM.quantizeCut(d_cut, d_bid, d_bcd, PL, PG);
  TIMING_stop("Cut_Information");


#if 0
  displayCutInfo(d_cut, d_bid);
#endif

  TIMING_stop("Cut_Section");

}
