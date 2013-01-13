// #################################################################
//
// FFV : Frontflow / violet Cartesian
//
// Copyright (c) 2012-2013  All right reserved.
//
// Institute of Industrial Science, University of Tokyo, Japan. 
//
// #################################################################

/** 
 * @file   ffv_Initialize.C
 * @brief  FFV Class
 * @author kero
 */

#include "ffv.h"
#include "limits.h"


int FFV::Initialize(int argc, char **argv)
{
  double TotalMemory   = 0.0;  ///< 計算に必要なメモリ量（ローカル）
  double PrepMemory    = 0.0;  ///< 初期化に必要なメモリ量（ローカル）
  double G_TotalMemory = 0.0;  ///< 計算に必要なメモリ量（グローバル）
  double G_PrepMemory  = 0.0;  ///< 初期化に必要なメモリ量（グローバル）
  double tmp_memory    = 0.0;  ///< 計算に必要なメモリ量（グローバル）？
  
  double flop_task     = 0.0;  ///< flops計算用
  
  // CPMバージョン表示
  Hostonly_
  {
    cpm_Base::VersionInfo();
  }
  
  
  // 固定パラメータ
  fixed_parameters();
  

  
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
    FBUtility::printVersion(fp,     "Welcome to FFV  ", FFV_VERS);
    FBUtility::printVersion(stdout, "Welcome to FFV  ", FFV_VERS);
    
    FBUtility::printVersion(fp,     "FlowBase        ", FB_VERS);
    FBUtility::printVersion(stdout, "FlowBase        ", FB_VERS);
  }
  
  // 入力ファイルの指定
  string input_file = argv[1];

  
  // ffvのパラメータローダのインスタンス生成
  TPControl tp_ffv;
  
  tp_ffv.getTPinstance();
  
  int ierror = tp_ffv.readTPfile(input_file);

  
  // TPControlクラスのポインタを各クラスに渡す
  C.importTP(&tp_ffv);
  B.importTP(&tp_ffv);
  M.importTP(&tp_ffv);
  MO.importTP(&tp_ffv);
  
  
  // 例題の種類を取得し，C.Mode.Exampleにフラグをセットする
  getExample(&C, &tp_ffv);
  
  // 組み込み例題クラスの実体をインスタンスし，*Exにポイントする
  connectExample(&C);
  
  
  // 組み込み例題クラス名を表示
  Hostonly_
  {
    Ex->printExample(fp, Ex->getExampleName());
  }
  
  
  // ランク情報をセット >> 各クラスでランク情報メンバ変数を利用する前にセットすること
  C.setRankInfo(paraMngr, procGrp);
  B.setRankInfo(paraMngr, procGrp);
  V.setRankInfo(paraMngr, procGrp);
  F.setRankInfo(paraMngr, procGrp);
  BC.setRankInfo(paraMngr, procGrp);
  Ex->setRankInfo(paraMngr, procGrp);
  MO.setRankInfo(paraMngr, procGrp);
  FP3DR.setRankInfo(paraMngr, procGrp);
  FP3DW.setRankInfo(paraMngr, procGrp);
  
  
  // 最初のパラメータの取得
  C.get_Steer_1(&DT, &FP3DR, &FP3DW);
  
  

  
  // 領域設定 計算領域全体のサイズ，並列計算時のローカルのサイズ，コンポーネントのサイズなどを設定
  DomainInitialize(&tp_ffv);

  
  // 各クラスで領域情報を保持
  C.setNeighborInfo(C.guide);
  B.setNeighborInfo(C.guide);
  V.setNeighborInfo(C.guide);
  F.setNeighborInfo(C.guide);
  BC.setNeighborInfo(C.guide);
  Ex->setNeighborInfo(C.guide);
  MO.setNeighborInfo(C.guide);
  FP3DR.setNeighborInfo(C.guide);
  FP3DW.setNeighborInfo(C.guide);
  
  
  // 各例題のパラメータ設定 -----------------------------------------------------
  Ex->setDomain(&C, size, origin, region, pitch);

  
  // パラメータを取得
  C.get_Steer_2(IC, &RF);
  
  
  // 組み込み例題の固有パラメータ
  if ( !Ex->getTP(&C, &tp_ffv) ) Exit(0);
  
  
  // 媒質情報をパラメータファイルから読み込み，媒質リストを作成する
  Hostonly_
  {
    printf(    "\n---------------------------------------------------------------------------\n\n");
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
    printf(    "\n\t>> Medium List\n\n");
    fprintf(fp,"\n\t>> Medium List\n\n");
  }
  
  // 媒質情報をロードし、 MediumTableタグ内の媒質数を保持
  C.NoMedium = M.get_MediumTable();
  
  // 媒質リストをインスタンス
  mat = new MediumList[C.NoMedium+1];
  
  
  // 媒質情報を設定
  setMediumList(fp);
  
  
  // パラメータファイルから C.NoBC, C.NoCompoを取得
  C.NoBC    = B.getNoLocalBC();    // LocalBoundaryタグ内の境界条件の個数
  C.NoCompo = C.NoBC + C.NoMedium; // コンポーネントの数の定義

  // ParseMatクラスの環境設定 
  M.setControlVars(C.NoCompo, C.NoBC, C.Unit.Temp, C.KindOfSolver);
  
  V.setNoCompo_BC(C.NoBC, C.NoCompo);
  
  B.setControlVars(&C);
  
  B.countMedium(&C, mat);
  
  // CompoListクラスをインスタンス．[0]はダミーとして利用しないので，配列の大きさはプラス１する
  cmp = new CompoList[C.NoCompo+1];

  

  
  // タイミング測定の初期化
  if ( C.Mode.Profiling != OFF )
  {
    ModeTiming = ON;
    TIMING__ PM.initialize( tm_END );
    TIMING__ PM.setRankInfo( paraMngr->GetMyRankID() );
    TIMING__ PM.setParallelMode(setParallelism(), C.num_thread, C.num_process);
    set_timing_label();
  }
  
  
  // タイミング測定開始
  TIMING_start(tm_init_sct); 
  
  
  
  // 前処理に用いるデータクラスのアロケート -----------------------------------------------------
  TIMING_start(tm_init_alloc); 
  allocArray_Prep(PrepMemory, TotalMemory);
  TIMING_stop(tm_init_alloc);
  
  
  // ファイルからIDを読み込む，または組み込み例題クラスでID情報を作成
  Hostonly_
  {
    printf(    "\n---------------------------------------------------------------------------\n\n");
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
    printf(    "\t>> Analysis Model Information\n\n");
    fprintf(fp,"\t>> Analysis Model Information\n\n");
  }
  
  TIMING_start(tm_voxel_prep_sct);
  
  
  
  // 各問題に応じてモデルを設定
  setModel(PrepMemory, TotalMemory, fp);
  
  
  
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
  
  display_memory_info(fp, G_PrepMemory, PrepMemory, "Preprocessor");
  
  
  
  // CompoList, MediumListのポインタをセット
  BC.importCMP_MAT(cmp, mat);
  
  
  // CompoListの設定，外部境界条件の読み込み保持
  setBCinfo();
  

  
  // Binaryの場合に，SOLIDセルを生成
  if ( !C.isCDS() && (C.Mode.Example == id_Polygon) )
  {
    generate_Solid(fp);
  }
  
  
  
  // ガイドセル上にパラメータファイルで指定する媒質IDを代入する．周期境界の場合の処理も含む．
  for (int face=0; face<NOFACE; face++)
  {
    V.adjMedium_on_GC(face, d_mid, BC.export_OBC(face)->get_Class(),
                      BC.export_OBC(face)->get_GuideMedium(), BC.export_OBC(face)->get_PrdcMode());
  }

  
  
  // Fill
  if ( (C.Mode.Example == id_Polygon) )
  {
    Hostonly_
    {
      printf(    "\n---------------------------------------------------------------------------\n\n");
      fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
      printf(    "\t>> Fill\n\n");
      fprintf(fp,"\t>> Fill\n\n");
    }
    
    fill(fp);
  }
  
  
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
  
  
  // ボクセルのスキャン
  VoxScan(fp);
  
  
  // スキャンしたセルIDの情報を表示する
  Hostonly_
  {
    printf("\n---------------------------------------------------------------------------\n\n");
    printf("\t>> Information of Scanned Voxel\n\n");
    V.printScannedCell(stdout);
		fflush(stdout);
    
    fprintf(fp, "\n---------------------------------------------------------------------------\n\n");
    fprintf(fp,"\t>> Information of Scanned Voxel\n\n");
		V.printScannedCell(fp);
		fflush(fp);
  }
  
  // ボクセルモデルの媒質インデクスがパラメータファイルに記述された媒質インデクスに含まれていること
  if ( !V.chkIDconsistency(C.NoMedium) )
  {
    Hostonly_
    {
			stamped_printf(    "\tIDs in between parameter file and scanned model are not consistent.\n");
      stamped_fprintf(fp,"\tIDs in between parameter file and scanned model are not consistent.\n");
		}
    return -1;
	}


  
  
  // Cell_Monitorの指定がある場合，モニタ位置をセット
  if ( (C.Sampling.log == ON) && (C.isMonitor() == ON) ) 
  {
    // ShapeMonitorのインスタンス
    ShapeMonitor SM(size, guide, (float*)pitch, (float*)origin);
    
    V.setShapeMonitor(d_mid, &SM, cmp, C.RefLength);
  }
  
  
  // CDSの場合，WALLとSYMMETRICのときに，カットを外部境界に接する内部セルに設定
  if ( C.isCDS() ) 
  {
    V.setOBC_Cut(&BC, d_cut);
  }
  
  
  // セルIDのノード間同期
  if ( paraMngr->BndCommS3D(d_mid, size[0], size[1], size[2], guide, 1) != CPM_SUCCESS ) Exit(0);
  

  // HEX/FANコンポーネントの形状情報からBboxと体積率を計算
  if ( C.isVfraction() )
  {
    TIMING_start(tm_init_alloc);
    allocArray_CompoVF(PrepMemory, TotalMemory);
    TIMING_stop(tm_init_alloc); 
    
    setComponentVF();
  }

  
  
  // 内部周期境界の場合のガイドセルのコピー処理
  V.adjMediumPrdc_Inner(d_mid, cmp);
  
  
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
  

  // BCIndexへのエンコード処理
  Hostonly_
  {
    printf(    "\n---------------------------------------------------------------------------\n\n");
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
  }
  
  
  // コンポーネント情報を保存
  //setLocalCmpIdx_Binary();
  
  // ポリゴンからBCのコンポーネント情報を設定
  Bbox_IBC();
  
  
  // BCIndexにビット情報をエンコードとコンポーネントインデクスの再構築
  VoxEncode();
  
  
  // 体積力を使う場合のコンポーネント配列の確保
  TIMING_start(tm_init_alloc);
  allocArray_Forcing(PrepMemory, TotalMemory, fp);
  TIMING_stop(tm_init_alloc); 
  

  // コンポーネントの体積率を8bitで量子化し，圧力損失コンポの場合にはFORCING_BITをON > bcdにエンコード
  V.setCmpFraction(cmp, d_bcd, d_cvf);

// ########## 
#if 1
  // CompoListとMediumListの関連を表示
  Hostonly_ 
  {
    M.printRelation(stdout, cmp, mat);
    M.printRelation(fp, cmp, mat);
  }
#endif
// ########## 

  
  // Ref_MediumがMediumList中にあるかどうかをチェックし、RefMatを設定
  if ( (C.RefMat = C.find_ID_from_Label(mat, C.NoMedium, C.Ref_Medium)) == 0 )
  {
    Hostonly_
    {
      printf(     "RefMat[%d] is not listed in MediumTable.\n", C.RefMat);
      fprintf(fp, "RefMat[%d] is not listed in MediumTable.\n", C.RefMat);
    }
    return -1;
  }


  // 周期境界条件が設定されている場合のBCIndexの周期条件の強制同期
  BC.setBCIperiodic(d_bcd);
  BC.setBCIperiodic(d_bcp);
  BC.setBCIperiodic(d_bcv);

  if ( C.isHeatProblem() )
  {
    BC.setBCIperiodic(d_bh1);
    BC.setBCIperiodic(d_bh2);
  }

  // bcd/bcp/bcv/bchの同期
  if ( numProc > 1 )
  {
    if ( paraMngr->BndCommS3D(d_bcd, size[0], size[1], size[2], guide, 1) != CPM_SUCCESS ) Exit(0);
    if ( paraMngr->BndCommS3D(d_bcp, size[0], size[1], size[2], guide, 1) != CPM_SUCCESS ) Exit(0);
    if ( paraMngr->BndCommS3D(d_bcv, size[0], size[1], size[2], guide, 1) != CPM_SUCCESS ) Exit(0);
    
    if ( C.isHeatProblem() )
    {
      if ( paraMngr->BndCommS3D(d_bh1, size[0], size[1], size[2], guide, 1) != CPM_SUCCESS ) Exit(0);
      if ( paraMngr->BndCommS3D(d_bh2, size[0], size[1], size[2], guide, 1) != CPM_SUCCESS ) Exit(0);
    }
  }

  
  // 法線計算
  get_Compo_Area();
  
  
  
  // 時間積分幅 deltaT や物理パラメータの設定
  setParameters();


  
  // 必要なパラメータをSetBC3Dクラスオブジェクトにコピーする >> setParameters()の後
  BC.setControlVars(&C, mat, cmp, &RF, Ex);
  
  
// ##########
#if 0
  // チェックのため，全計算セルのBCIndexの内容を表示する
  if ( !V.dbg_chkBCIndexP(bcd, bcp, "BCindex.txt") )
  {
    Hostonly_
    {
      stamped_printf("\tVoxInfo::DbgChkBCIndexP()\n");
    }
    return -1;
  }
#endif
// ##########
  
  
  
  // 温度計算の場合の初期値指定
  if ( C.isHeatProblem() )
  {
    B.get_Medium_InitTemp(cmp);
  }
  
  // set phase 
  if ( C.BasicEqs == INCMP_2PHASE )
  {
    B.get_Phase(cmp);
  }
  

  
  // CompoListの内容とセル数の情報を表示する
  Hostonly_
  {
    printf(    "\n---------------------------------------------------------------------------\n\n");
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
    printf(    "\t>> Component List\n\n");
    fprintf(fp,"\t>> Component List\n\n");
  }
  display_CompoList(fp);
  
  
  // 外部境界面の開口率を計算する
  V.countOpenAreaOfDomain(d_bcd, C.OpenDomain);
  
  
  Hostonly_
  {
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n\n");
    printf(    "\n---------------------------------------------------------------------------\n\n\n");
  }
  

  
  // Monitor Listの処理 --------------------------------------------
  MO.setControlVars(d_bcd,
                    C.RefVelocity, C.BaseTemp, C.DiffTemp, C.RefDensity, C.RefLength, C.BasePrs,
                    C.Unit.Temp, C.Mode.Precision, C.Unit.Prs, C.num_process);
  
  
  // モニタ機能がONの場合
  if ( C.Sampling.log == ON )
  {
    //パラメータを取得し，セットの配列を確保する
	  MO.get_Monitor(&C);
    
	  //プローブ位置をID=255としてボクセルファイルに書き込む
	  MO.write_ID(d_mid);
  }
  
  
  // 内部境界条件として指定されたモニタ設定を登録
  if ( (C.Sampling.log == ON) && (C.isMonitor() == ON) ) MO.setInnerBoundary(cmp, C.NoBC);
  
  
  // 性能測定モードがオフのときのみ，出力指定，あるいはMonitorListの場合に，svxファイルを出力する．
  if ( C.Hide.PM_Test == OFF )
  {
    if ( (C.Sampling.log == ON) || (C.FIO.IO_Voxel == Control::Sphere_SVX) )
    {
      Ex->writeSVX(d_mid, &C);
    }
  }
  
  
  // mid[]を解放する  ---------------------------
  if ( d_mid ) delete [] d_mid;
  
  
  
  // コンポーネントのグローバルインデクス情報を取得
  setGlobalCmpIdx();

  
  
  // コンポーネントの内容リストを表示し、コンポーネント数がゼロの場合と境界条件との整合性をチェック
  Hostonly_
  {
    printf(    "\n---------------------------------------------------------------------------\n\n");
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
    printf(    "\t>> Component Information\n");
    fprintf(fp,"\t>> Component Information\n");
  }
  display_Compo_Info(fp);

  
  // 各ノードの領域情報をファイル出力
  gather_DomainInfo();
  

  
  TIMING_stop(tm_voxel_prep_sct);
  // ここまでがボクセル準備の時間セクション
  
  
// ##########  
#if 0
  write_distance(cut);
#endif
// ##########  
  

  
  // 計算に用いる配列のアロケート ----------------------------------------------------------------------------------
  allocate_Main(TotalMemory);
  
  
  
  // 分散時のインデクスファイル生成
  setDFI();
  
  
  
  // スタート処理 瞬時値と平均値に分けて処理　------------------
  Hostonly_
  {
    printf(    "\n---------------------------------------------------------------------------\n\n");
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
  }
  
  
  
  // リスタート処理
  Restart(fp);

  
  
  // 制御インターバルの初期化
  init_Interval();
  
  
  // 平均値のロード
  if ( C.Start == restart )
  {
    TIMING_start(tm_restart);
    if ( C.Mode.Average == ON ) Restart_avrerage(fp, flop_task);
    TIMING_stop(tm_restart);
  }
  
  
  // リスタートの最大値と最小値の表示
  Restart_display_minmax(fp, flop_task);
  
  
  
  // 制御パラメータ，物理パラメータの表示
  Hostonly_ 
  {
    printf(    "\n---------------------------------------------------------------------------\n\n");
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
    printf(    "\t>> Outer Boundary Conditions\n\n");
    fprintf(fp,"\t>> Outer Boundary Conditions\n\n");
    
    display_Parameters(fp);
  }
  

  
  // ドライバ条件のチェック
  BC.checkDriver(fp);
  
  
  // 初期条件の条件設定
  setInitialCondition();
  


  if ( C.Sampling.log == ON ) 
  {
    // サンプリング指定がある場合，モニタ結果出力ファイル群のオープン
    MO.openFile(C.HistoryMonitorName.c_str());
    
    // サンプリング元となるデータ配列の登録
    if ( C.isHeatProblem() ) 
    {
      MO.setDataPtrs(d_v, d_p, d_t);
    }
    else 
    {
      MO.setDataPtrs(d_v, d_p);
    }
  }
  
  
  // PLOT3D形状データの書き出し
  PLT3D.Initialize(size, guide, deltaX, dfi_mng[var_Plot3D], &C, &FP3DW, &DFI, d_ws, d_p, d_wo, d_v, d_t, d_p0, d_wv, d_bcv, d_bcd);
  
  if (C.FIO.PLOT3D_OUT == ON)
  {
    PLT3D.setValuePlot3D();
    if (C.P3Op.IS_xyz == ON) PLT3D.OutputPlot3D_xyz(CurrentStep, origin, pitch);// ---> moving grid を考慮したときOutputPlot3D_postに組み込む
    if (C.P3Op.IS_DivideFunc == ON)
    {
      if (C.P3Op.IS_function_name == ON) Hostonly_ PLT3D.OutputPlot3D_function_name_divide();
    }
    else
    {
      if (C.P3Op.IS_function_name == ON) Hostonly_ PLT3D.OutputPlot3D_function_name();
    }
    if (C.P3Op.IS_fvbnd == ON) PLT3D.OutputPlot3D_fvbnd();
  }

  
  // セッションを開始したときに、初期値をファイル出力  リスタートと性能測定モードのときには出力しない
  if ( (C.Hide.PM_Test == OFF) && (0 == CurrentStep) )
  {
    flop_task = 0.0;
    FileOutput(flop_task);
    if (C.FIO.PLOT3D_OUT == ON) PLT3D.OutputPlot3D_post(CurrentStep, CurrentTime, v00, origin, pitch, flop_task);
  }

  
  // 粗い格子を用いたリスタート時には出力
  if ( C.Start == coarse_restart )
  {
    flop_task = 0.0;
    FileOutput(flop_task, true);
    if (C.FIO.PLOT3D_OUT == ON) PLT3D.OutputPlot3D_post(CurrentStep, CurrentTime, v00, origin, pitch, flop_task);
  }
  

  
  // SOR2SMA
  switch (IC[ItrCtl::ic_prs_pr].get_LS())
  {
    case SOR2SMA:
    case SOR2CMA:
    case GMRES:
    case RBGS:
		case PCG:
		case PBiCGSTAB:
      allocate_SOR2SMA_buffer(TotalMemory);
      break;
  }
  
  // Krylov subspace
  switch (IC[ItrCtl::ic_prs_pr].get_LS())
  {
    case GMRES:
      allocArray_Krylov(TotalMemory);
      break;
      
    case PCG:
			allocArray_PCG(TotalMemory);
			break;
      
		case PBiCGSTAB:
			allocArray_PBiCGSTAB(TotalMemory);
			break;
  }
  
  
  // メモリ使用量の表示
  Hostonly_
  {
    printf(    "\n---------------------------------------------------------------------------\n\n");
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
  }
  
  G_TotalMemory = TotalMemory;
  
  display_memory_info(fp, G_TotalMemory, TotalMemory, "Solver");

  
  
  // 履歴出力準備
  prep_HistoryOutput();
  
  
  
  Hostonly_ if ( fp ) fclose(fp);
  
  TIMING_stop(tm_init_sct);
  
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

  return 1;
}



// #################################################################
// ポリゴンのカット情報からVBCのboxをセット
void FFV::Bbox_IBC()
{
  int f_st[3], f_ed[3], len[3];
  
  for (int n=1; n<=C.NoBC; n++) {
    
    int tg = cmp[n].getMatOdr();
    
    // インデクスの計算 > インデクスの登録はVoxEncode()で、コンポーネント領域のリサイズ後に行う
    if ( V.find_IBC_bbox(tg, d_bid, d_cut, f_st, f_ed) )
    {
      len[0] = f_ed[0] - f_st[0] + 1;
      len[1] = f_ed[1] - f_st[1] + 1;
      len[2] = f_ed[2] - f_st[2] + 1;
      
      for (int d=0; d<3; d++)
      {
        int tmp_st=0;
        int tmp_ed=0;
        
        EnlargeIndex(tmp_st, tmp_ed, f_st[d], len[d], size[d], d, tg);
        
        f_st[d] = tmp_st;
        f_ed[d] = tmp_ed;
      }
      
      cmp[n].setBbox(f_st, f_ed);
      cmp[n].setEns(ON);
    }
    else
    {
      cmp[n].setEns(OFF);
    }
  }
  
}


// #################################################################
// 全Voxelモデルの媒質数とKOSの整合性をチェック
bool FFV::chkMediumConsistency()
{
  int nmSolid = C.NoMediumSolid;
  int nmFluid = C.NoMediumFluid;
  
  if ( numProc > 1 )
  {
    int nms = nmSolid;
    int nmf = nmFluid;
    paraMngr->Allreduce(&nms, &nmSolid, 1, MPI_SUM);
    paraMngr->Allreduce(&nmf, &nmFluid, 1, MPI_SUM);
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
      
    case CONJUGATE_HEAT_TRANSFER:
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
// 組み込み例題のインスタンス
void FFV::connectExample(Control* Cref)
{
  if      ( Cref->Mode.Example == id_PPLT2D)   Ex = dynamic_cast<Intrinsic*>(new IP_PPLT2D);
  else if ( Cref->Mode.Example == id_SHC1D)    Ex = dynamic_cast<Intrinsic*>(new IP_SHC1D);
  else if ( Cref->Mode.Example == id_Duct )    Ex = dynamic_cast<Intrinsic*>(new IP_Duct);
  else if ( Cref->Mode.Example == id_PMT )     Ex = dynamic_cast<Intrinsic*>(new IP_PMT);
  else if ( Cref->Mode.Example == id_Rect )    Ex = dynamic_cast<Intrinsic*>(new IP_Rect);
  else if ( Cref->Mode.Example == id_Cylinder) Ex = dynamic_cast<Intrinsic*>(new IP_Cylinder);
  else if ( Cref->Mode.Example == id_Step )    Ex = dynamic_cast<Intrinsic*>(new IP_Step);
  else if ( Cref->Mode.Example == id_Polygon ) Ex = dynamic_cast<Intrinsic*>(new IP_Polygon);
  else if ( Cref->Mode.Example == id_Sphere )  Ex = dynamic_cast<Intrinsic*>(new IP_Sphere);
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for Exmple definition\n");
    Exit(0);
  }
}


// #################################################################
// 時刻をRFクラスからv00[4]にコピーする
void FFV::copyV00fromRF(double m_time) 
{
  RF.setV00(m_time);
  
  double g[4];
  RF.copyV00(g);
  for (int i=0; i<4; i++) v00[i]=(REAL_TYPE)g[i];
}


// #################################################################
// コンポーネントの内容リストを表示する
void FFV::display_Compo_Info(FILE* fp)
{
  if ( C.NoBC >0 )
  {
    Hostonly_
    {
      B.printCompo( stdout, compo_global_bbox, mat, cmp, BC.export_OBC() );
      B.printCompo( fp, compo_global_bbox, mat, cmp, BC.export_OBC() );
    }
  }
  
  // コンポーネント数がゼロの場合のチェック
  for (int n=1; n<=C.NoBC; n++) {
    if ( cmp[n].getElement() == 0 )
    {
      Hostonly_ printf("\tError : No element was found in Component[%d]\n", n);
      fflush(stdout);
      Exit(0);
    }
  }
  
  // Check consistency of boundary condition
  for (int n=1; n<=C.NoBC; n++) {
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
// CompoListの内容とセル数の情報を表示する
void FFV::display_CompoList(FILE* fp)
{
  Hostonly_
  {
    M.chkList(stdout, cmp, C.BasicEqs);
    M.chkList(fp, cmp, C.BasicEqs);
  }
  
  // セル数の情報を表示する
  Hostonly_ {
    double cr = (double)G_Wcell/ ( (double)G_size[0] * (double)G_size[1] * (double)G_size[2]) *100.0;
    fprintf(stdout, "\tThis model includes %4d solid %s  [Solid cell ratio inside computational domain : %9.5f percent]\n\n", 
            C.NoMediumSolid, (C.NoMediumSolid>1) ? "IDs" : "ID", cr);
    fprintf(fp, "\tThis model includes %4d solid %s  [Solid cell ratio inside computational domain : %9.5f percent]\n\n", 
            C.NoMediumSolid, (C.NoMediumSolid>1) ? "IDs" : "ID", cr);
  }
  
}


// #################################################################
// 制御パラメータ，物理パラメータの表示
void FFV::display_Parameters(FILE* fp)
{
  C.displayParams(stdout, fp, IC, &DT, &RF, mat, &FP3DW);
  Ex->printPara(stdout, &C);
  Ex->printPara(fp, &C);
  
  // 外部境界面の開口率を表示
  C.printOuterArea(stdout, G_Fcell, G_Acell, G_size);
  C.printOuterArea(fp, G_Fcell, G_Acell, G_size);
  
  // 境界条件のリストと外部境界面のBC設定を表示

  B.printFaceOBC(stdout, G_region, BC.export_OBC(), mat);
  B.printFaceOBC(fp, G_region, BC.export_OBC(), mat);

  
  // モニタ情報の表示
  if ( C.Sampling.log == ON ) {
    
    MO.printMonitorInfo(stdout, C.HistoryMonitorName.c_str(), false); // ヘッダのみ
    
    FILE *fp_mon=NULL;
    Hostonly_
    {
      if ( !(fp_mon=fopen("sampling_info.txt", "w")) )
      {
        stamped_printf("\tSorry, can't open 'sampling_info.txt' file. Write failed.\n");
        //return -1;
        Exit(0);
      }
    }
    
    MO.printMonitorInfo(fp_mon, C.HistoryMonitorName.c_str(), true);  // 詳細モード
    Hostonly_ if ( fp_mon ) fclose(fp_mon);
  }
}


// #################################################################
// 計算領域情報を設定する
void FFV::DomainInitialize(TPControl* tp_dom)
{
  // メンバ変数にパラメータをロード : 分割指示 (1-with / 2-without)
  int div_type = get_DomainInfo(tp_dom);

  
// ##########  
# if 0
  printDomainInfo();
  fflush(stdout);
#endif
// ##########
  

  // 袖通信の最大数
  size_t Nvc  = (size_t)C.guide;
  size_t Ncmp = 6; // カットを使うので6成分（面）
  
  int m_sz[3]  = {G_size[0], G_size[1], G_size[2]};
  int m_div[3] = {G_division[0], G_division[1], G_division[2]};
  

  // 有次元の場合に無次元化する　paraMngr->VoxelInit()で計算する基本量
  if (C.Unit.Param == DIMENSIONAL )
  {
    for (int i=0; i<3; i++) {
      pitch[i]    /= C.RefLength;
      G_origin[i] /= C.RefLength;
      G_region[i] /= C.RefLength;
    }
  }

  
  double m_org[3] = {(double)G_origin[0], (double)G_origin[1], (double)G_origin[2]};
  double m_reg[3] = {(double)G_region[0], (double)G_region[1], (double)G_region[2]};
  
  
  // 領域分割モードのパターン
  //      分割指定(G_div指定)    |     domain.txt 
  // 1)  G_divなし >> 自動分割   |  G_orign + G_region + (G_pitch || G_voxel)
  // 2)  G_div指定あり          |  G_orign + G_region + (G_pitch || G_voxel)
  // 3)  G_divなし >> 自動分割   |   + ActiveDomainInfo
  // 4)  G_div指定あり          |   + ActiveDomainInfo
  
  // 分割数を元に分割する。pitchが指定値とならないこともある。
  switch (div_type) 
  {
    case 1: // 分割数が指示されている場合
      if ( paraMngr->VoxelInit(m_div, m_sz, m_org, m_reg, Nvc, Ncmp) != CPM_SUCCESS )
      {
        cout << "Domain decomposition error : " << endl;
        Exit(0);
      }
      break;
      
    case 2: // 分割数が指示されていない場
      if ( paraMngr->VoxelInit(m_sz, m_org, m_reg, Nvc, Ncmp) != CPM_SUCCESS )
      {
        cout << "Domain decomposition error : " << endl;
        Exit(0);
      }
      break;
      
    default:
      Exit(0);
      break;
  }


  // 分割後のパラメータをDomainInfoクラスメンバ変数に保持
  setNeighborInfo(C.guide);

  
  
  // チェック
  unsigned long tz = (unsigned long)size[0] * (unsigned long)size[1] * (unsigned long)size[2];
  if ( tz >= UINT_MAX)
  {
    Hostonly_ stamped_printf("\n\tError : Product of size[] exceeds UINT_MAX\n\n");
    Exit(0);
  }

}


// #################################################################
//初期インデクスの情報を元に，一層拡大したインデクス値を返す
void FFV::EnlargeIndex(int& m_st, int& m_ed, const int st_i, const int len, const int m_x, const int dir, const int m_id)
{
  int ed_i = st_i + len - 1;
  int n_st = st_i - 1;
  int n_ed = ed_i + 1;
  int max_c1 = m_x + guide;
  
  int label_st, label_ed;
  
  switch (dir) 
  {
    case 0:
      label_st = X_MINUS;
      label_ed = X_PLUS;
      break;
      
    case 1:
      label_st = Y_MINUS;
      label_ed = Y_PLUS;
      break;
      
    case 2:
      label_st = Z_MINUS;
      label_ed = Z_PLUS;
      break;
      
    default:
      Hostonly_ stamped_printf("\tError : DIRECTION\n");
      Exit(0);
  }
  
  // BVが-方向のガイドセル内のみにある場合
  if ( ed_i < guide )
  {
    if( nID[label_st] < 0 ) // 計算領域の外部面に接する場合は，対象外
    {
      m_st = 0;
      m_ed = 0;
    }
    else // 計算領域内部にある場合（並列時）
    {
      if ( n_ed == guide ) // ガイドセル1層外側の場合
      {
        m_st = 1; // F index
        m_ed = 1; // F index
      }
      else
      {
        m_st = 0;
        m_ed = 0;
      }
    }
  }
  
  // BVが+方向のガイドセル内のみにある場合
  else if ( st_i >= max_c1 )
  {
    if( nID[label_ed] < 0 ) // 計算領域の外部面に接する場合は，対象外
    {
      m_st = 0;
      m_ed = 0;
    }
    else
    {
      if ( n_st == (max_c1 - 1) ) // ガイドセル1層外側の場合
      {
        m_st = m_x; // F index
        m_ed = m_x; // F index
      }
      else
      {
        m_st = 0;
        m_ed = 0;
      }
    }
    //debug; Hostonly_ printf("(2) dir=%d : st=%d ed=%d\n", dir, m_st, m_ed);
  }
  
  // BVが内部領域のみにある場合（逐次・並列で同じ処理）
  else if ( (st_i >= guide) && (ed_i < max_c1) )
  {
    if ( st_i == guide ) // 最外層セル
    {
      m_st = 1; // F index
    }
    else
    {
      m_st = n_st + 1 - guide; // F index
    }
    
    if ( ed_i == (max_c1 - 1) ) // 最外層セル
    {
      m_ed = m_x; // F index
    }
    else // 内部
    {
      m_ed = n_ed + 1 - guide; // F index
    }
    //debug; Hostonly_ printf("(3) dir=%d : st=%d ed=%d\n", dir, m_st, m_ed);
  }
  
  // BVが両方向のガイドセルにまたがる場合（逐次・並列で同じ処理）
  else if ( (st_i < guide) && (ed_i >= max_c1) )
  {
    m_st = 1; // F index
    m_ed = m_x; // F index
    //debug; Hostonly_ printf("(4) dir=%d : st=%d ed=%d\n", dir, m_st, m_ed);
  }
  
  // BVが-方向のガイドセルから内部領域にある場合
  else if ( (st_i < guide) && (ed_i < max_c1) )
  {
    m_st = 1; // F index
    
    if ( ed_i == (max_c1 - 1) ) // 最外層セル
    {
      m_ed = m_x; // F index
    }
    else // 内部
    {
      m_ed = n_ed + 1 - guide; // F index
    }
    //debug; Hostonly_ printf("(5) dir=%d : st=%d ed=%d\n", dir, m_st, m_ed);
  }
  
  // BVが+方向のガイドセルから内部領域にある場合
  else if ( (st_i < max_c1) && (ed_i >= max_c1) )
  {
    m_ed = m_x; // F index
    
    if ( st_i == guide ) // 端点
    {
      m_st = 1; // F index
    }
    else // 内部
    {
      m_st = n_st + 1 - guide; // F index
    }
    //debug; Hostonly_ printf("(6) dir=%d : st=%d ed=%d\n", dir, m_st, m_ed);
  }
  
  else
  {
    string m_dir;
    if      ( dir == 0 ) m_dir = "X";
    else if ( dir == 1 ) m_dir = "Y";
    else                 m_dir = "Z";
    
    Hostonly_
    {
      stamped_printf("\tError : Unexpected case for ID=%d, (%d - %d): %s\n", m_id, st_i, ed_i, m_dir.c_str());
    }
    Exit(0);
  }
  
}



// #################################################################
// ポリゴンの場合のフィル操作
void FFV::fill(FILE* fp)
{
  
  // 指定媒質の属性をチェック
  bool flag = false;
  
  for (int i=C.NoBC+1; i<=C.NoCompo; i++) {
    if ( (cmp[i].getMatOdr() == C.Fill_Fluid) && (cmp[i].getState() == FLUID) )
    {
      flag = true;
    }
  }
  if ( !flag )
  {
    Hostonly_ printf("\tSpecified Medium of filling fluid is not FLUID\n");
    Exit(0);
  }
  
  
  flag = false;
  for (int i=C.NoBC+1; i<=C.NoCompo; i++) {
    if ( (cmp[i].getMatOdr() == C.Fill_Solid) && (cmp[i].getState() == SOLID) )
    {
      flag = true;
    }
  }
  if ( !flag )
  {
    Hostonly_ printf("\tSpecified Medium of filling solid is not SOLID\n");
    Exit(0);
  }

  
  unsigned long fc;
  
  // 最初にフィル対象のセル数を求める
  unsigned long fill_count = (unsigned long)size[0] * (unsigned long)size[1] * (unsigned long)size[2];
  unsigned long fs = 0;
  
  if ( numProc > 1 )
  {
    unsigned long tmp_fc = fill_count;
    if ( paraMngr->Allreduce(&tmp_fc, &fill_count, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
  }
  
  Hostonly_
  {
    printf(    "\t1st Fill -----\n\n\tFLUID\n");
    fprintf(fp,"\t1st Fill -----\n\n\tFLUID\n");

    printf    ("\t\tInitial target count : %15ld\n", fill_count);
    fprintf(fp,"\t\tInitial target count : %15ld\n", fill_count);
  }
  
  
  // 指定された媒質を使って、指定シード点を与える
  Hostonly_
  {
    printf    ("\t\tFilled by medium     : %s\n", mat[C.Fill_Fluid].getLabel().c_str());
    fprintf(fp,"\t\tFilled by medium     : %s\n", mat[C.Fill_Fluid].getLabel().c_str());
  }
  
  
  
  // 1st pass
  
  // ヒントが与えられている場合
  // 確実に流体のセルのみをペイントする
  
  if ( C.Fill_Hint >= 0 )
  {
    fs = V.fill_seed(d_mid, C.Fill_Hint, C.Fill_Fluid, d_cut);

    if ( numProc > 1 )
    {
      unsigned long tmp_fs = fs;
      if ( paraMngr->Allreduce(&tmp_fs, &fs, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
    }
    
    if ( fs == 0 )
    {
      Hostonly_
      {
        printf(    "Failed first painting (%s includes solid cell)\n", FBUtility::getDirection(C.Fill_Hint).c_str());
        fprintf(fp,"Failed first painting (%s includes solid cell)\n", FBUtility::getDirection(C.Fill_Hint).c_str());
      }
      Exit(0);
    }
    
    Hostonly_
    {
      printf(    "\t\tPainted %ld by Hint (%s includes solid cell)\n", fs, FBUtility::getDirection(C.Fill_Hint).c_str());
      fprintf(fp,"\t\tPainted %ld by Hint (%s includes solid cell)\n", fs, FBUtility::getDirection(C.Fill_Hint).c_str());
    }
  }

  // シード分のカウントデクリメント
  fill_count -= fs;
  
  Hostonly_
  {
    printf("\n");
    fprintf(fp,"\n");
  }



  // BIDによるフィル
  // 隣接する流体セルと接続しており，かつ固体セルに挟まれていないセルのみペイントする
  int c=0;
  while (fill_count > 0) {
    
    fc = (unsigned long)V.fill_by_bid(d_bid, d_mid, d_cut, C.Fill_Fluid, C.Fill_Solid);
    
    if ( numProc > 1 )
    {
      unsigned long t_fc = fc;
      if ( paraMngr->Allreduce(&t_fc, &fc, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
    }
    
    if ( fc == 0 ) break; // フィル対象がなくなったら終了
    
    fill_count -= fc;
    c++;
    
    // 同期
    if ( numProc > 1 )
    {
      if ( paraMngr->BndCommS4DEx(d_cut, 6, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
      if ( paraMngr->BndCommS3D(d_mid, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
      if ( paraMngr->BndCommS3D(d_bid, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
    }
  }
  
  Hostonly_
  {
    printf(    "\t\tIteration = %5d : FLUID filling by BID = %15ld\n", c+1, fill_count);
    fprintf(fp,"\t\tIteration = %5d : FLUID filling by BID = %15ld\n", c+1, fill_count);
  }
  
  
  
  // midによる穴埋め
  c = 0;
  unsigned long z1 = 0;
  unsigned long z2 = 0;
  
  while (fill_count > 0) {
    
    z1 = (unsigned long)V.fill_by_mid(d_bid, d_mid, d_cut, C.Fill_Fluid, C.Fill_Solid);
    
    if ( numProc > 1 )
    {
      unsigned long t_fc = z1;
      if ( paraMngr->Allreduce(&t_fc, &z1, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
    }
    
    z2 += z1;
    
    if ( z1 == 0 ) break; // フィル対象がなくなったら終了

    c++;
    
    // 同期
    if ( numProc > 1 )
    {
      if ( paraMngr->BndCommS4DEx(d_cut, 6, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
      if ( paraMngr->BndCommS3D(d_mid, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
      if ( paraMngr->BndCommS3D(d_bid, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
    }
  }
  
  Hostonly_
  {
    printf(    "\t\tIteration = %5d : Hole  filling by MID = %15ld\n", c+1, z2);
    fprintf(fp,"\t\tIteration = %5d : Hole  filling by MID = %15ld\n", c+1, z2);
  }
  
  
  
  // 2nd pass
  
  Hostonly_
  {
    printf(    "\n\t2nd Fill -----\n\n\tFLUID\n");
    fprintf(fp,"\n\t2nd Fill -----\n\n\tFLUID\n");
  }
  
  // 既にペイントした流体セルをクリア
  unsigned long fz = (unsigned long)V.fill_inside(d_mid, C.Fill_Fluid, 0);
  
  if ( numProc > 1 )
  {
    unsigned long t_fc = fz;
    if ( paraMngr->Allreduce(&t_fc, &fz, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
  }
  
  fill_count += fz;
  
  
  
  if ( C.Fill_Hint >= 0 ) // ヒントが与えられている場合
  {
    fs = V.fill_seed(d_mid, C.Fill_Hint, C.Fill_Fluid, d_cut);
    
    if ( numProc > 1 )
    {
      unsigned long tmp_fs = fs;
      if ( paraMngr->Allreduce(&tmp_fs, &fs, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
    }
    
    if ( fs == 0 )
    {
      Hostonly_
      {
        printf(    "Failed second painting (%s includes solid cell)\n", FBUtility::getDirection(C.Fill_Hint).c_str());
        fprintf(fp,"Failed second painting (%s includes solid cell)\n", FBUtility::getDirection(C.Fill_Hint).c_str());
      }
      Exit(0);
    }
  }
  
  
  
  // BIDによるフィル
  c = 0;
  while (fill_count > 0) {
    
    fc = (unsigned long)V.fill_by_bid(d_bid, d_mid, d_cut, C.Fill_Fluid, C.Fill_Solid);
    
    if ( numProc > 1 )
    {
      unsigned long t_fc = fc;
      if ( paraMngr->Allreduce(&t_fc, &fc, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
    }
    
    if ( fc == 0 ) break; // フィル対象がなくなったら終了
    
    fill_count -= fc;
    c++;
    
    // 同期
    if ( numProc > 1 )
    {
      if ( paraMngr->BndCommS4DEx(d_cut, 6, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
      if ( paraMngr->BndCommS3D(d_mid, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
      if ( paraMngr->BndCommS3D(d_bid, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
    }
  }
  
  Hostonly_
  {
    printf(    "\t\tIteration = %5d : FLUID filling by BID = %15ld\n\n", c+1, fill_count);
    fprintf(fp,"\t\tIteration = %5d : FLUID filling by BID = %15ld\n\n", c+1, fill_count);
  }
  
  
  
  // 固体に変更
  // Allreduce時の桁あふれ対策のため、unsigned long で集約
  c = 0;
  while ( fill_count > 0 ) {
    
    // 未ペイントのセルに対して、固体IDを与える
    unsigned long fc = (unsigned long)V.fill_inside(d_mid, 0, C.Fill_Solid);
    unsigned long t_fc = fc;
    
    if ( numProc > 1 )
    {
      if ( paraMngr->Allreduce(&t_fc, &fc, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
    }
    
    if ( fc == 0 ) break;
    
    fill_count -= fc;
    c++;
  }
  
  Hostonly_
  {
    printf(    "\tSOLID\n");
    fprintf(fp,"\tSOLID\n");

    printf(    "\t\tIteration = %5d : Filled cell =          %15ld\n\n", c, fill_count);
    fprintf(fp,"\t\tIteration = %5d : Filled cell =          %15ld\n\n", c, fill_count);
  }
  

  
  // 確認 paintedは未ペイントセルがある場合に1
  unsigned long painted = (unsigned long)V.fill_check(d_mid);
  
  if ( numProc > 1 )
  {
    unsigned long t_painted = painted;
    if ( paraMngr->Allreduce(&t_painted, &painted, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
  }
  
  if ( painted  != 0 )
  {
    Hostonly_
    {
      printf(    "\tFill operation is done, but still remains unpainted cells.\n");
      fprintf(fp,"\tFill operation is done, but still remains unpainted cells.\n");
    }
    Exit(0);
  }

}



// #################################################################
// 固定パラメータの設定
void FFV::fixed_parameters()
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
  
  // ログファイル名
  C.HistoryName        = "history_base.txt";
  C.HistoryCompoName   = "history_compo.txt";
  C.HistoryDomfxName   = "history_domainflux.txt";
  C.HistoryForceName   = "history_force.txt";
  C.HistoryWallName    = "history_log_wall.txt";
  C.HistoryItrName     = "history_iteration.txt";
  C.HistoryMonitorName = "sampling.txt";
  
  C.f_Pressure       = "prs";
  C.f_Velocity       = "vel";
  C.f_Fvelocity      = "fvel";
  C.f_Temperature    = "tmp";
  C.f_AvrPressure    = "prsa";
  C.f_AvrVelocity    = "vela";
  C.f_AvrTemperature = "tmpa";
  C.f_DivDebug       = "div";
  C.f_Helicity       = "hlt";
  C.f_TotalP         = "tp";
  C.f_I2VGT          = "i2vgt";
  C.f_Vorticity      = "vrt";
  
  C.FIO.IO_Input  = IO_DISTRIBUTE;
  C.FIO.IO_Output = IO_DISTRIBUTE;
  
}


// #################################################################
// IOモードに対応したディレクトリパスを返す
string FFV::directory_prefix(string path, const string fname, const int io_mode, const int para_mode)
{
  string tmp;
  
  switch (io_mode)
  {
    case Control::io_current:
      tmp = fname;
      break;
      
      
    case Control::io_specified:
      
      if ( !FBUtility::c_mkdir(path) )
      {
        Hostonly_ printf("Failed to create directory \"%s\"\n", path.c_str() );
        Exit(0);
      }
      tmp = path + "/" + fname;
      break;
      
      
    case Control::io_time_slice:
      
      // 1プロセスの場合にはランク番号がないので、タイムスライス毎のディレクトリは作らない
      if ( (para_mode == Control::Serial) || (para_mode == Control::OpenMP) )
      {
        return fname;
      }
      else
      {
        if ( !FBUtility::c_mkdir(path) )
        {
          Hostonly_ printf("Failed to create directory \"%s\"\n", path.c_str() );
          Exit(0);
        }
        tmp = path + "/" + fname;
      }
      break;
  }
  
  return tmp;
}


// #################################################################
// メモリ使用量の表示
void FFV::display_memory_info(FILE* fp, double G_mem, double L_mem, const char* str)
{
  if ( numProc > 1 )
  {
    double tmp_memory = G_mem;
    if ( paraMngr->Allreduce(&tmp_memory, &G_mem, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
  }
  
  Hostonly_
  {
    FBUtility::MemoryRequirement(str, G_mem, L_mem, stdout);
    FBUtility::MemoryRequirement(str, G_mem, L_mem, fp);
  }
  
  Hostonly_
  {
    printf("\n\n");
    fprintf(fp, "\n\n");
  }
}


// #################################################################
// 並列処理時の各ノードの分割数を集めてファイルに保存する
void FFV::gather_DomainInfo()
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
    if ( paraMngr->Gather(size, 3, m_size, 3, 0) != CPM_SUCCESS ) Exit(0);
    if ( paraMngr->Gather(origin, 3, m_org, 3, 0) != CPM_SUCCESS ) Exit(0);
    if ( paraMngr->Gather(region, 3, m_reg, 3, 0) != CPM_SUCCESS ) Exit(0);
    if ( paraMngr->Gather(&L_Fcell, 1, bf_fcl, 1, 0) != CPM_SUCCESS ) Exit(0);
    if ( paraMngr->Gather(&L_Wcell, 1, bf_wcl, 1, 0) != CPM_SUCCESS ) Exit(0);
    if ( paraMngr->Gather(&L_Acell, 1, bf_acl, 1, 0) != CPM_SUCCESS ) Exit(0);
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
  
  if ( nID[X_MINUS] < 0 ) m_srf -= (double)(jx*kx);  // remove face which does not join communication
  if ( nID[Y_MINUS] < 0 ) m_srf -= (double)(ix*kx);
  if ( nID[Z_MINUS] < 0 ) m_srf -= (double)(ix*jx);
  if ( nID[X_PLUS]  < 0 ) m_srf -= (double)(jx*kx);
  if ( nID[Y_PLUS]  < 0 ) m_srf -= (double)(ix*kx);
  if ( nID[Z_PLUS]  < 0 ) m_srf -= (double)(ix*jx);
  
  if ( numProc > 1 )
  {
    if ( paraMngr->Gather(&m_srf, 1, bf_srf, 1, 0) != CPM_SUCCESS ) Exit(0);
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
  
  for (int i=0; i<numProc; i++) {
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
  C.printGlobalDomain(fp, G_size, G_origin, G_region, pitch);
  
  // ローカルノードの情報を表示
  for (int i=0; i<numProc; i++) {
    Hostonly_ {
      fprintf(fp,"Domain %4d\n", i);
      fprintf(fp,"\t ix, jx,  kx        [-] =  %13ld %13ld %13ld\n",  m_size[i*3], m_size[i*3+1], m_size[i*3+2]);
      fprintf(fp,"\t(ox, oy, oz)  [m] / [-] = (%13.6e %13.6e %13.6e)  /  (%13.6e %13.6e %13.6e)\n", 
              m_org[i*3]*C.RefLength,  m_org[i*3+1]*C.RefLength,  m_org[i*3+2]*C.RefLength, m_org[i*3],  m_org[i*3+1],  m_org[i*3+2]);
      fprintf(fp,"\t(Lx, Ly, Lz)  [m] / [-] = (%13.6e %13.6e %13.6e)  /  (%13.6e %13.6e %13.6e)\n", 
              m_reg[i*3]*C.RefLength,  m_reg[i*3+1]*C.RefLength,  m_reg[i*3+2]*C.RefLength, m_reg[i*3],  m_reg[i*3+1],  m_reg[i*3+2]);
      
      if (C.NoBC != 0) fprintf(fp, "\t no            Label   Mat    i_st    i_ed    j_st    j_ed    k_st    k_ed\n");
    }
    
    if ( numProc > 1 )
    {
      for (int n=1; n<=C.NoBC; n++) {
        if( paraMngr->Gather(cmp[n].getBbox_st(), 3, st_buf, 3, 0) != CPM_SUCCESS ) Exit(0);
        if( paraMngr->Gather(cmp[n].getBbox_ed(), 3, ed_buf, 3, 0) != CPM_SUCCESS ) Exit(0);
        
        Hostonly_
        {
          fprintf(fp,"\t%3d %16s %5d %7d %7d %7d %7d %7d %7d\n",
                  n, cmp[n].getLabel().c_str(), cmp[n].getMatOdr(), st_buf[i*3], ed_buf[i*3], st_buf[i*3+1], ed_buf[i*3+1], st_buf[i*3+2], ed_buf[i*3+2]);
        }
      }
    }
    else // serial
    {
      int *st, *ed;
      for (int n=1; n<=C.NoBC; n++) {
        st = cmp[n].getBbox_st();
        ed = cmp[n].getBbox_ed();
        
        Hostonly_
        {
          fprintf(fp,"\t%3d %16s %5d %7d %7d %7d %7d %7d %7d\n",
                  n, cmp[n].getLabel().c_str(), cmp[n].getMatOdr(), st[0], ed[0], st[1], ed[1], st[2], ed[2]);
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
    fprintf(fp,"\t---------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
    
    double tmp_vol, tmp_acl, tmp_fcl, tmp_wcl;
    for (int i=0; i<numProc; i++) {
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
// Binary voxelをカット情報から生成
void FFV::generate_Solid(FILE* fp)
{
  unsigned long zc=0;
  

  for (int m=1; m<=C.NoCompo; m++) {
    
    if ( cmp[m].getState() == SOLID )
    {
      int tgt = cmp[m].getMatOdr();
      
      zc += V.Solid_from_Cut(d_mid, d_bid, d_cut, tgt);
    }
  }
  
  /* BC
  for (int m=1; m<=C.NoBC; m++) {
    
    int target = cmp[m].getMatOdr();
    int m_dir  = cmp[m].getBClocation();
    float vec[3] = { (float)cmp[m].nv[0], (float)cmp[m].nv[1], (float)cmp[m].nv[2] };
    
    zc += V.Solid_from_Cut_VBC(d_mid, target, C.Fill_Solid, vec, m_dir);
  }*/
      
  Hostonly_
  {
    printf(    "\n\tGenerated Solid cell from cut = %ld\n", zc);
    fprintf(fp,"\n\tGenerated Solid cell from cut = %ld\n", zc);
  }
  
  // midのガイドセル同期
  if ( paraMngr->BndCommS3D(d_mid, size[0], size[1], size[2], guide, 1) != CPM_SUCCESS ) Exit(0);
  
}



// #################################################################
// コンポーネントの面積を計算
void FFV::get_Compo_Area()
{

  if ( C.NoBC == 0 ) return;
  if ( C.Mode.Example != id_Polygon ) return;
  
  
  float area=0.0;
  
  // コンポーネントで指定されるID面の法線を計算
  for (int n=1; n<=C.NoBC; n++) {
    int type = cmp[n].getType();
    
    if ( (type==SPEC_VEL) || (type==SPEC_VEL_WH) || (type==OUTFLOW) )
    {
      string label = cmp[n].getLabel();
      
      for (int i=0; i<C.num_of_polygrp; i++) {
        
        if ( FBUtility::compare(poly_prop[i].label_grp, label) )
        {
          area = poly_prop[i].area;
          
          if ( numProc > 1 )
          {
            float ta = area;
            if ( paraMngr->Allreduce(&ta, &area, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
          }
          
          cmp[n].area = area;
        }
      }

    }
  }
}


// #################################################################
// グローバルな領域情報を取得
int FFV::get_DomainInfo(TPControl* tp_dom)
{
  // 領域分割モードのパターン
  //      分割指定(G_div指定)    |     domain.txt 
  // 1)  G_divなし >> 自動分割   |  G_orign + G_region + (G_pitch || G_voxel)
  // 2)  G_div指定あり          |  G_orign + G_region + (G_pitch || G_voxel)
  // 3)  G_divなし >> 自動分割   |   + ActiveDomainInfo
  // 4)  G_div指定あり          |   + ActiveDomainInfo
  
  
  string label, str;
  REAL_TYPE *rvec;
  int *ivec;
  int div_type = 1; // 指定分割 => 1
  
  
  // 長さの単位
  label = "/DomainInfo/UnitOfLength";
  
  if ( !tp_dom->GetValue(label, &str) )
  {
		Hostonly_ stamped_printf("\tParsing error : Invalid string for '%s'\n", label.c_str());
	  Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "NonDimensional") )  C.Unit.Length = LTH_ND;
  else if( !strcasecmp(str.c_str(), "M") )               C.Unit.Length = LTH_m;
  else if( !strcasecmp(str.c_str(), "cm") )              C.Unit.Length = LTH_cm;
  else if( !strcasecmp(str.c_str(), "mm") )              C.Unit.Length = LTH_mm;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described at '%s'\n", label.c_str());
    Exit(0);
  }
  
  
  // G_origin　必須
  rvec  = G_origin;
  label = "/DomainInfo/GlobalOrigin";
  
  if ( !tp_dom->GetVector(label, rvec, 3) )
  {
    cout << "ERROR : in parsing [" << label << "]" << endl;
    Exit(0);
  }
  
  // G_region 必須
  rvec  = G_region;
  label = "/DomainInfo/GlobalRegion";
  
  if ( !tp_dom->GetVector(label, rvec, 3) )
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
  
  
  // pitch オプション
  bool flag = true; // 排他チェック（voxel, pitch）
  rvec  = pitch;
  label = "/DomainInfo/GlobalPitch";
  
  if ( !tp_dom->GetVector(label, rvec, 3) )
  {
    Hostonly_ cout << "\tNo option : in parsing [" << label << "]" << endl;
    flag = false;
  }
  
  // pitchが入力されている場合のみチェック
  if ( flag )
  {
    if ( (pitch[0]>0.0) && (pitch[1]>0.0) && (pitch[2]>0.0) )
    {
      ; // skip
    }
    else
    {
      Hostonly_ printf("ERROR : in parsing [%s] >> (%e, %e, %e)\n", label.c_str(), pitch[0], pitch[1], pitch[2] );
      Exit(0);
    }
  }  
  
  
  // G_size オプション
  if ( flag ) // pitchが指定されている場合が優先
  {
    G_size[0] = (int)(G_region[0]/pitch[0]);
    G_size[1] = (int)(G_region[1]/pitch[1]);
    G_size[2] = (int)(G_region[2]/pitch[2]);
  }
  else
  {
    ivec  = G_size;
    label = "/DomainInfo/GlobalVoxel";
    
    if ( !tp_dom->GetVector(label, ivec, 3) )
    {
      Hostonly_ cout << "ERROR : Neither GlobalPitch nor GlobalVoxel is specified." << endl;
      Exit(0); // pitchもvoxelも有効でない
    }
    
    if ( (G_size[0]>0) && (G_size[1]>0) && (G_size[2]>0) )
    {
      pitch[0] = G_region[0] / (REAL_TYPE)G_size[0];
      pitch[1] = G_region[1] / (REAL_TYPE)G_size[1];
      pitch[2] = G_region[2] / (REAL_TYPE)G_size[2];
      
      pitch[1] = pitch[0];
      pitch[2] = pitch[0];
    }
    else
    {
      Hostonly_ cout << "ERROR : in parsing [" << label << "]" << endl;
      Exit(0);
    }
  }
  
  
  // G_division オプション
  ivec  = G_division;
  label = "/DomainInfo/GlobalDivision";
  
  if ( !tp_dom->GetVector(label, ivec, 3) )
  {
    Hostonly_ cout << "\tNo option : in parsing [" << label << "]" << endl;
    div_type = 2; // 自動分割
  }
  
  // プロセス分割数が指定されている場合のチェック
  if ( div_type == 1 )
  {
    if ( (G_division[0]>0) && (G_division[1]>0) && (G_division[2]>0) ) 
    {
      ; // skip
    }
    else
    {
      Hostonly_ cout << "ERROR : in parsing [" << label << "]" << endl;
      Exit(0);
    }
  }

  
  // 有次元の場合には、単位をメートルに変換
  REAL_TYPE factor;
  switch ( C.Unit.Length )
  {
    case LTH_m:
    case LTH_ND:
      factor = 1.0; // 変更なし
      break;
      
    case LTH_cm:
      factor = 0.01;
      break;
      
    case LTH_mm:
      factor = 0.001;
      break;
      
    default:
      Exit(0);
  }
  
  for (int i=0; i<3; i++) {
    pitch[i]    *= factor;
    G_origin[i] *= factor;
    G_region[i] *= factor;
  }
  
  
  // ActiveSubdomainファイル名の取得
  label = "/DomainInfo/ActiveSubDomainFile";
  
  if ( !tp_dom->GetValue(label, &str ) )
  {
    Hostonly_ cout << "\tNo option : in parsing [" << label << "]" << endl;
  }
  //@todo  string hoge = str;
  
  
  return div_type;
}


// #################################################################
// 組み込み例題の表示
void FFV::getExample(Control* Cref, TPControl* tpCntl)
{
  string keyword;
  string label;
  
  label = "/Steer/Example";
  
  if ( !(tpCntl->GetValue(label, &keyword )) )
  {
    Hostonly_ stamped_printf("\tError : '%s'\n", label.c_str());
    Exit(0);
  }
  
  if     ( FBUtility::compare(keyword, "ParallelPlate2D") )   Cref->Mode.Example = id_PPLT2D;
  else if( FBUtility::compare(keyword, "Duct") )              Cref->Mode.Example = id_Duct;
  else if( FBUtility::compare(keyword, "SHC1D") )             Cref->Mode.Example = id_SHC1D;
  else if( FBUtility::compare(keyword, "PerformanceTest") )   Cref->Mode.Example = id_PMT;
  else if( FBUtility::compare(keyword, "Rectangular") )       Cref->Mode.Example = id_Rect;
  else if( FBUtility::compare(keyword, "Cylinder") )          Cref->Mode.Example = id_Cylinder;
  else if( FBUtility::compare(keyword, "BackStep") )          Cref->Mode.Example = id_Step;
  else if( FBUtility::compare(keyword, "Users") )             Cref->Mode.Example = id_Polygon;
  else if( FBUtility::compare(keyword, "Sphere") )            Cref->Mode.Example = id_Sphere;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
    Exit(0);
  }
}


// #################################################################
// インターバルの初期化
void FFV::init_Interval()
{
  
  // セッションの初期時刻をセット
  for (int i=0; i<Interval_Manager::tg_END; i++) {
    C.Interval[i].setTime_init(Session_StartTime);
  }
  
  
  // 入力モードが有次元の場合に，無次元に変換
  if ( C.Unit.Param == DIMENSIONAL )
  {
    C.Interval[Interval_Manager::tg_compute].normalizeInterval(C.Tscale);
    C.Interval[Interval_Manager::tg_console].normalizeInterval(C.Tscale);
    C.Interval[Interval_Manager::tg_history].normalizeInterval(C.Tscale);
    C.Interval[Interval_Manager::tg_instant].normalizeInterval(C.Tscale);
    C.Interval[Interval_Manager::tg_average].normalizeInterval(C.Tscale);
    C.Interval[Interval_Manager::tg_accelra].normalizeInterval(C.Tscale);
    C.Interval[Interval_Manager::tg_avstart].normalizeInterval(C.Tscale);
    C.Interval[Interval_Manager::tg_sampled].normalizeInterval(C.Tscale);
    C.Interval[Interval_Manager::tg_plot3d].normalizeInterval(C.Tscale);
  }
  
  // Reference frame
  RF.setAccel( C.Interval[Interval_Manager::tg_accelra].getIntervalTime() );
  
  
  
  // インターバルの初期化
  double m_dt    = DT.get_DT();
  double m_tm    = CurrentTime;  // Restart()で設定
  unsigned m_stp = CurrentStep;
  
  if ( !C.Interval[Interval_Manager::tg_compute].initTrigger(m_stp, m_tm, m_dt, Interval_Manager::tg_compute, C.Tscale) )
  {
    Hostonly_ printf("\t Error : Computation Period is asigned to zero.\n");
    Exit(0);
  }
  
  if ( !C.Interval[Interval_Manager::tg_console].initTrigger(m_stp, m_tm, m_dt, Interval_Manager::tg_console, C.Tscale) )  // 基本履歴のコンソールへの出力
  {
    Hostonly_ printf("\t Error : Interval for Console output is asigned to zero.\n");
    Exit(0);
  }
  if ( !C.Interval[Interval_Manager::tg_history].initTrigger(m_stp, m_tm, m_dt, Interval_Manager::tg_history, C.Tscale) )  // 履歴のファイルへの出力
  {
    Hostonly_ printf("\t Error : Interval for History output is asigned to zero.\n");
    Exit(0);
  }
  if ( !C.Interval[Interval_Manager::tg_instant].initTrigger(m_stp, m_tm, m_dt, Interval_Manager::tg_instant, C.Tscale) )  // 瞬時値ファイル
  {
    Hostonly_ printf("\t Error : Interval for Instantaneous output is asigned to zero.\n");
    Exit(0);
  }
  if ( C.Mode.Average == ON )
  {
    // tg_averageの初期化はLoop中で行う．平均値開始時刻が未知のため．
    
    if ( !C.Interval[Interval_Manager::tg_avstart].initTrigger(m_stp, m_tm, m_dt, Interval_Manager::tg_avstart, C.Tscale) )  // 平均値開始
    {
      Hostonly_ printf("\t Error : Interval for Average start is asigned to zero.\n");
      Exit(0);
    }
  }
    
  if ( C.Sampling.log == ON )
  {
    if ( !C.Interval[Interval_Manager::tg_sampled].initTrigger(m_stp, m_tm, m_dt, Interval_Manager::tg_sampled, C.Tscale) )  // サンプリング履歴
    {
      Hostonly_ printf("\t Error : Interval for Sampling output is asigned to zero.\n");
      Exit(0);
    }    
  }
  
  if (C.FIO.PLOT3D_OUT == ON)
  {
    if ( !C.Interval[Interval_Manager::tg_plot3d].initTrigger(m_stp, m_tm, m_dt, Interval_Manager::tg_plot3d, C.Tscale) )  // 瞬時値ファイル
    {
      Hostonly_ printf("\t Error : Interval for plot3d output is asigned to zero.\n");
      Exit(0);
    }
  }
  
  Session_LastStep = C.Interval[Interval_Manager::tg_compute].getIntervalStep();
  
}


// #################################################################
// 距離の最小値を求める
void FFV::min_distance(float* cut, FILE* fp)
{
  float global_min;
  float local_min = 1.0;
  float eps = 1.0/255.0; // 0.000392
  unsigned g=0;
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
#pragma omp parallel firstprivate(ix, jx, kx, eps, gd)
  {
    float th_min = 1.0;
    
#pragma omp for schedule(static) reduction(+:g)
    for (int k=1; k<=kx; k++) {
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          
          for (int l=0; l<6; l++) {
            size_t m = _F_IDX_S4DEX(l, i, j, k, 6, ix, jx, kx, gd);
            float c = cut[m];
            
            th_min = min(th_min, c); //if ( local_min > c ) local_min = c;
            
            if ( (c > 0.0) && (c <= eps) )
            {
              cut[m] = eps;
              g++;
            }
          }
          
        }
      }
    }
    
#pragma omp critical
    {
      local_min = min(local_min, th_min);
    }
  }
    
  global_min = local_min;
  
  // Allreduce時の桁あふれ対策のため、unsigned long で集約
  unsigned long gl = (unsigned long)g;
  
  if ( numProc > 1 )
  {
    float tmp = global_min;
    if ( paraMngr->Allreduce(&tmp, &global_min, 1, MPI_MIN) != CPM_SUCCESS ) Exit(0);
    
    unsigned long tmp_g = gl;
    if ( paraMngr->Allreduce(&tmp_g, &gl, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
  }

  if ( gl > 0 )
  {
    Hostonly_
    {
      fprintf(fp, "\n\tMinimum non-dimnensional distance is %e and replaced to %e : num = %ld\n\n", global_min, eps, gl);
      printf     ("\n\tMinimum non-dimnensional distance is %e and replaced to %e : num = %ld\n\n", global_min, eps, gl);
    }
  }

}


// #################################################################
// 履歴の出力準備
void FFV::prep_HistoryOutput()
{
  // マスターノードでの履歴出力準備
  H = new History(&C);
  
  Hostonly_ {
    H->printHistoryTitle(stdout, IC, &C, true);
    
    // コンポーネント情報
    if ( C.Mode.Log_Base == ON ) 
    {
      // 基本情報　history.log, history_compo.log, history_domfx.log
      if ( !(fp_b=fopen(C.HistoryName.c_str(), "w")) ) 
      {
        stamped_printf("\tSorry, can't open '%s' file. Write failed.\n", C.HistoryName.c_str());
        Exit(0);
      }
      H->printHistoryTitle(fp_b, IC, &C, false);
      
      // コンポーネント履歴情報
      if ( !(fp_c=fopen(C.HistoryCompoName.c_str(), "w")) ) 
      {
        stamped_printf("\tSorry, can't open '%s' file. Write failed.\n", C.HistoryCompoName.c_str());
        Exit(0);
      }
      H->printHistoryCompoTitle(fp_c, cmp, &C);
      
      // 流量収支情報　
      if ( !(fp_d=fopen(C.HistoryDomfxName.c_str(), "w")) ) 
      {
        stamped_printf("\tSorry, can't open '%s' file. Write failed.\n", C.HistoryDomfxName.c_str());
        Exit(0);
      }
      H->printHistoryDomfxTitle(fp_d, &C);
      
      // 力の履歴情報　
      if ( !(fp_f=fopen(C.HistoryForceName.c_str(), "w")) ) 
      {
        stamped_printf("\tSorry, can't open '%s' file. Write failed.\n", C.HistoryForceName.c_str());
        Exit(0);
      }
      H->printHistoryForceTitle(fp_f);
    }
    
    // 反復履歴情報　history_itr.log
    if ( C.Mode.Log_Itr == ON ) 
    {
      if ( !(fp_i=fopen(C.HistoryItrName.c_str(), "w")) ) 
      {
				stamped_printf("\tSorry, can't open '%s' file.\n", C.HistoryItrName.c_str());
        Exit(0);
      }
    }
    
    // 壁面情報　history_wall.log
    if ( C.Mode.Log_Wall == ON ) 
    {
      if ( !(fp_w=fopen(C.HistoryWallName.c_str(), "w")) ) 
      {
				stamped_printf("\tSorry, can't open '%s' file.\n", C.HistoryWallName.c_str());
        Exit(0);
      }
      H->printHistoryWallTitle(fp_w);
    }
  }
  
  // サンプリング指定がある場合，モニタ結果出力ファイル群のオープン
  if ( C.Sampling.log == ON ) MO.openFile(C.HistoryMonitorName.c_str());
}


// #################################################################
// 読み込んだ領域情報のデバッグライト
void FFV::printDomainInfo()
{
  cout << "\n####### read parameters ########" << endl;
  cout << " Gorg      = " << G_origin[0] << "," << G_origin[1] << "," << G_origin[2] << endl;
  cout << " Gvoxel    = " << G_size[0]   << "," << G_size[1]   << "," << G_size[2]   << endl;
  cout << " Gpitch    = " << pitch[0]    << "," << pitch[1]    << "," << pitch[2]    << endl;
  cout << " Gregion   = " << G_region[0] << "," << G_region[1] << "," << G_region[2] << endl;
  cout << " Gdiv      = " << G_division[0]  << "," << G_division[1]  << "," << G_division[2]  << endl;
}


// #################################################################
// コンポーネントリストに登録されたセル要素BCのBV情報をリサイズする
void FFV::resizeBVface(const int* st, const int* ed, const int n, const int* bx)
{
  int s;
  size_t m;
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  // 初期値はローカルノードの大きさ
  int nst[3] = {ix, jx, kx};
  int ned[3] = {0, 0, 0};
  
  for (int k=st[2]; k<=ed[2]; k++) {
    for (int j=st[1]; j<=ed[1]; j++) {
      for (int i=st[0]; i<=ed[0]; i++) {
        
        //m = FBUtility::getFindexS3D(size, guide, i, j, k);
        m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        s = bx[m];
        
        if ( GET_FACE_BC(s, BC_FACE_W) == n ) 
        {
          if( i < nst[0] ) { nst[0] = i; }
          if( i > ned[0] ) { ned[0] = i; }
          if( j < nst[1] ) { nst[1] = j; }
          if( j > ned[1] ) { ned[1] = j; }
          if( k < nst[2] ) { nst[2] = k; }
          if( k > ned[2] ) { ned[2] = k; }
        }
        
        if ( GET_FACE_BC(s, BC_FACE_E) == n ) 
        {
          if( i < nst[0] ) { nst[0] = i; }
          if( i > ned[0] ) { ned[0] = i; }
          if( j < nst[1] ) { nst[1] = j; }
          if( j > ned[1] ) { ned[1] = j; }
          if( k < nst[2] ) { nst[2] = k; }
          if( k > ned[2] ) { ned[2] = k; }
        }
        
        if ( GET_FACE_BC(s, BC_FACE_S) == n ) 
        {
          if( i < nst[0] ) { nst[0] = i; }
          if( i > ned[0] ) { ned[0] = i; }
          if( j < nst[1] ) { nst[1] = j; }
          if( j > ned[1] ) { ned[1] = j; }
          if( k < nst[2] ) { nst[2] = k; }
          if( k > ned[2] ) { ned[2] = k; }
        }
        
        if ( GET_FACE_BC(s, BC_FACE_N) == n ) 
        {
          if( i < nst[0] ) { nst[0] = i; }
          if( i > ned[0] ) { ned[0] = i; }
          if( j < nst[1] ) { nst[1] = j; }
          if( j > ned[1] ) { ned[1] = j; }
          if( k < nst[2] ) { nst[2] = k; }
          if( k > ned[2] ) { ned[2] = k; }
        }
        
        if ( GET_FACE_BC(s, BC_FACE_B) == n ) 
        {
          if( i < nst[0] ) { nst[0] = i; }
          if( i > ned[0] ) { ned[0] = i; }
          if( j < nst[1] ) { nst[1] = j; }
          if( j > ned[1] ) { ned[1] = j; }
          if( k < nst[2] ) { nst[2] = k; }
          if( k > ned[2] ) { ned[2] = k; }
        }
        
        if ( GET_FACE_BC(s, BC_FACE_T) == n ) 
        {
          if( i < nst[0] ) { nst[0] = i; }
          if( i > ned[0] ) { ned[0] = i; }
          if( j < nst[1] ) { nst[1] = j; }
          if( j > ned[1] ) { ned[1] = j; }
          if( k < nst[2] ) { nst[2] = k; }
          if( k > ned[2] ) { ned[2] = k; }
        }
      }
    }
  }
  
  // replace
  cmp[n].setBbox(nst, ned);
}


// #################################################################
// コンポーネントリストに登録されたセル要素BCのBV情報をリサイズする
void FFV::resizeBVcell(const int* st, const int* ed, const int n, const int* bx)
{
  int s;
  size_t m;
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  // 初期値はローカルノードの大きさ
  int nst[3] = {ix, jx, kx};
  int ned[3] = {0, 0, 0};

  
  for (int k=st[2]; k<=ed[2]; k++) {
    for (int j=st[1]; j<=ed[1]; j++) {
      for (int i=st[0]; i<=ed[0]; i++) {
        
        //m = FBUtility::getFindexS3D(size, guide, i, j, k);
        m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        s = bx[m];
        
        if ( ( s & MASK_6) == n ) 
        {
          if( i < nst[0] ) { nst[0] = i; }
          if( i > ned[0] ) { ned[0] = i; }
          if( j < nst[1] ) { nst[1] = j; }
          if( j > ned[1] ) { ned[1] = j; }
          if( k < nst[2] ) { nst[2] = k; }
          if( k > ned[2] ) { ned[2] = k; }
        }
      }
    }
  }
  
  // replace
  cmp[n].setBbox(nst, ned);
}


// #################################################################
// コンポーネントリストに登録されたBV情報をリサイズする
void FFV::resizeCompoBV(const int kos, const bool isHeat)
{
  int st[3], ed[3];
  int typ;
  
  for (int n=1; n<=C.NoBC; n++) 
  {
    cmp[n].getBbox(st, ed);
    typ = cmp[n].getType();
    
    // デフォルトで流れ計算パートのBC
    switch ( typ ) 
    {
      case SPEC_VEL:
      case SPEC_VEL_WH:
      case OUTFLOW:
        resizeBVface(st, ed, n, d_bcv); // 速度のBCindex　セルフェイス位置でのflux型BC
        break;
        
      case PERIODIC:
      case IBM_DF:
      case HEX:
      case FAN:
      case DARCY:
        resizeBVcell(st, ed, n, d_bcd); // セルセンタ位置でのBC
        break;
    }
    
    // 熱計算のパート
    if ( isHeat ) 
    {
      switch ( typ ) 
      {
        case ADIABATIC:
        case HEATFLUX:
        case SPEC_VEL_WH:
        case OUTFLOW:
        case TRANSFER:
        case ISOTHERMAL:
          resizeBVface(st, ed, n, d_bh1);
          break;
          
        case RADIANT:
          break;
          
        case HEAT_SRC:
        case CNST_TEMP:
        case PERIODIC:
          resizeBVcell(st, ed, n, d_bh2);
          break;
      }            
    }
    
  }
}



// #################################################################
// 外部境界条件を読み込み，Controlクラスに保持する
void FFV::setBCinfo()
{
  // パラメータファイルをパースして，外部境界条件を保持する　>> VoxScan()以前に処理
  B.loadBC_Outer( BC.export_OBC(), M.export_MTI(), cmp );
  
  
  // パラメータファイルの情報を元にCompoListの情報を設定する
  B.loadBC_Local(&C, mat, cmp, poly_prop);
  
  
  // 各コンポーネントが存在するかどうかを保持しておく
  setEnsComponent();
  
  // KOSと境界条件種類の整合性をチェック
  B.chkBCconsistency(C.KindOfSolver, cmp);
  
}



// #################################################################
// HEX,FANコンポーネントなどの体積率とbboxなどをセット
// インデクスの登録と配列確保はVoxEncode()で、コンポーネント領域のリサイズ後に行う
void FFV::setComponentVF()
{
  const int subsampling = 20; // 体積率のサブサンプリングの基数
  int f_st[3], f_ed[3];
  double flop;
  
  CompoFraction CF(size, guide, (float*)pitch, (float*)origin, subsampling);
  
  for (int n=1; n<=C.NoBC; n++) 
  {
    
    if ( cmp[n].isFORCING() ) 
    {
      // 形状パラメータのセット
      switch ( cmp[n].getType() ) 
      {
        case HEX:
          CF.setShapeParam((float*)cmp[n].nv, (float*)cmp[n].oc, (float*)cmp[n].dr, (float)cmp[n].depth, (float)cmp[n].shp_p1, (float)cmp[n].shp_p2);
          break;
          
        case FAN:
          CF.setShapeParam((float*)cmp[n].nv, (float*)cmp[n].oc, (float)cmp[n].depth, (float)cmp[n].shp_p1, (float)cmp[n].shp_p2);
          break;
          
        case DARCY:
          Exit(0);
          break;
          
        default:
          Exit(0);
          break;
      }
      
      // 回転角度の計算
      CF.get_angle(); 
      
      // bboxと投影面積の計算
      cmp[n].area = CF.get_BboxArea();
      
      // インデクスの計算 > あとで，VoxEncode()でresize
      CF.bbox_index(f_st, f_ed);
      
      // インデクスのサイズ登録と存在フラグ
      cmp[n].setBbox(f_st, f_ed);
      cmp[n].setEns(ON);
      
      // 体積率
      TIMING_start(tm_cmp_vertex8);
      flop = 0.0;
      CF.vertex8(f_st, f_ed, d_cvf, flop);
      TIMING_stop(tm_cmp_vertex8, flop);
      
      TIMING_start(tm_cmp_subdivision);
      flop = 0.0;
      CF.subdivision(f_st, f_ed, d_cvf, flop);
      TIMING_stop(tm_cmp_subdivision, flop);
    }
  }
  
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
  F.writeRawSPH(cvf, size, guide, org, pit, sizeof(float));
#endif
// ##########
  
}



// #################################################################
// コンポーネントが存在するかを保持しておく
void FFV::setEnsComponent()
{
  int c;
  
  // Forcing > HEX, FAN, DARCY
  c = 0;
  for (int n=1; n<=C.NoBC; n++) 
  {
    if ( cmp[n].isFORCING() ) c++;
  }
  if ( c>0 ) C.EnsCompo.forcing = ON;
  
  // Heat source > HEAT_SRC, CNST_TEMP
  c = 0;
  for (int n=1; n<=C.NoBC; n++) 
  {
    if ( cmp[n].isHsrc() ) c++;
  }
  if ( c>0 ) C.EnsCompo.hsrc = ON;
  
  // 周期境界 > PERIODIC
  c = 0;
  for (int n=1; n<=C.NoBC; n++) 
  {
    if ( cmp[n].getType() == PERIODIC ) c++;
  }
  if ( c>0 ) C.EnsCompo.periodic = ON;
  
  // 流出境界 > OUTFLOW
  c = 0;
  for (int n=1; n<=C.NoBC; n++) 
  {
    if ( cmp[n].getType() == OUTFLOW ) c++;
  }
  if ( c>0 ) C.EnsCompo.outflow = ON;
  
  // 体積率コンポーネント
  c = 0;
  for (int n=1; n<=C.NoBC; n++) 
  {
    if ( cmp[n].isVFraction() ) c++;
  }
  if ( c>0 ) C.EnsCompo.fraction = ON;
  
  // モニタ
  c = 0;
  for (int n=1; n<=C.NoBC; n++) 
  {
    if ( cmp[n].isMONITOR() ) c++;
  }
  
  // MONITOR_LISTでCELL_MONITORが指定されている場合，C.EnsCompo.monitor==ON
  if ( (C.isMonitor() == ON) && (c < 1) ) 
  {
    Hostonly_ stamped_printf("\tError : CellMonitor in MonitorList is specified, however any Monitor can not be found.\n");
    Exit(0);
  }
  
  if ( (C.isMonitor() == OFF) && (c > 0) ) 
  {
    Hostonly_ stamped_printf("\tError : CellMonitor in MonitorList is NOT specified, however Monitor section is found in LocalBoundary.\n");
    Exit(0);
  }
  
}



// #################################################################
// コンポーネントのローカルなBbox情報からグローバルなBbox情報を求める
void FFV::setGlobalCmpIdx()
{
  int st_i, st_j, st_k, ed_i, ed_j, ed_k;
  int node_st_i, node_st_j, node_st_k;
  int st[3], ed[3];
  
  // グローバルインデクスの配列インスタンス
  compo_global_bbox = new int[6*(C.NoCompo+1)];
  int* cgb = compo_global_bbox;
  
  // ローカルインデクスからグローバルインデクスに変換
  for (int m=1; m<=C.NoCompo; m++) 
  {
    
    if ( !cmp[m].isEns() ) // コンポーネントが存在しないノードはゼロを代入
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
        
        cgb[6*m+0] = node_st_i + st_i;
        cgb[6*m+1] = node_st_j + st_j;
        cgb[6*m+2] = node_st_k + st_k;
        cgb[6*m+3] = node_st_i + ed_i;
        cgb[6*m+4] = node_st_j + ed_j;
        cgb[6*m+5] = node_st_k + ed_k;
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
  int* m_gArray = NULL;
  int* m_eArray = NULL;
  int array_size = 6*(C.NoCompo+1);
  int st_x, st_y, st_z, ed_x, ed_y, ed_z, es;
  
  if ( !(m_gArray = new int[numProc*6]) ) Exit(0);
  if ( !(m_eArray = new int[numProc]  ) ) Exit(0);
  
  for (int n=1; n<=C.NoBC; n++) 
  {
    if ( numProc > 1 ) 
    {
      es = ( cmp[n].isEns() ) ? 1 : 0;
      if ( paraMngr->Gather(&es, 1, m_eArray, 1, 0) != CPM_SUCCESS ) Exit(0);
      if ( paraMngr->Gather(&cgb[6*n], 6, m_gArray, 6, 0) != CPM_SUCCESS ) Exit(0);
      
      
      if (myRank == 0) { // マスターノードのみ
        
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



// #################################################################
// 初期条件の設定
void FFV::setInitialCondition()
{
  double flop_task;
  
  REAL_TYPE tm = CurrentTime * C.Tscale;
  
  if ( C.Start == initial_start ) {
		REAL_TYPE U0[3];
    
		// 速度の初期条件の設定
    if (C.Unit.Param == DIMENSIONAL) {
      U0[0] = C.iv.VecU/C.RefVelocity;
      U0[1] = C.iv.VecV/C.RefVelocity;
      U0[2] = C.iv.VecW/C.RefVelocity;
    }
    else {
      U0[0] = C.iv.VecU;
      U0[1] = C.iv.VecV;
      U0[2] = C.iv.VecW;
    }
		fb_set_vector_(d_v, size, &guide, U0, d_bcd);
    
		// 外部境界面の流出流量と移流速度
    DomainMonitor( BC.export_OBC(), &C);
    
		// 外部境界面の移流速度を計算し，外部境界条件を設定
    BC.OuterVBC_Periodic(d_v);
		BC.OuterVBC(d_v, d_v, d_bcv, tm, deltaT, &C, v00, flop_task);
    BC.InnerVBC(d_v, d_bcv, tm, v00, flop_task);
    BC.InnerVBC_Periodic(d_v, d_bcd);
    
		// 圧力
    REAL_TYPE ip;
    if (C.Unit.Param == DIMENSIONAL) {
      ip = FBUtility::convD2ND_P(C.iv.Pressure, C.BasePrs, C.RefDensity, C.RefVelocity, C.Unit.Prs);
    }
    else {
      ip = C.iv.Pressure;
    }
    
    U.xset(d_p, size, guide, ip, kind_scalar);
		BC.OuterPBC(d_p);
    
		// 温度
		if ( C.isHeatProblem() ) {
      REAL_TYPE it;
      if (C.Unit.Param == DIMENSIONAL) {
        it = FBUtility::convK2ND(C.iv.Temperature, C.BaseTemp, C.DiffTemp);
      }
      else {
        it = C.iv.Temperature;
      }
      
      U.xset(d_t, size, guide, it, kind_scalar);
      
      // コンポーネントの初期値
      for (int m=C.NoBC+1; m<=C.NoCompo; m++) {
        BC.setInitialTemp_Compo(m, d_bcd, d_t);
      }
      
			BC.OuterTBC(d_t);
		}
    
    // ユーザ例題のときに，速度の内部境界条件を設定する
    Ex->initCond(d_v, d_p);
    
  }
  else // リスタート時
  { 
    // 内部境界条件
    BC.InnerVBC(d_v, d_bcv, tm, v00, flop_task);
    BC.InnerVBC_Periodic(d_v, d_bcd);
    BC.InnerPBC_Periodic(d_p, d_bcd);
    
    // 外部境界条件
    BC.OuterVBC(d_v, d_v, d_bcv, tm, deltaT, &C, v00, flop_task);
    BC.OuterVBC_Periodic(d_v);
    
    // 流出境界の流出速度の算出
    // dummy
    Gemini_R* m_buf = new Gemini_R [C.NoBC];
    BC.mod_div(d_ws, d_bcv, tm, v00, m_buf, d_vf, d_v, flop_task);
    if ( m_buf ) { delete [] m_buf; m_buf=NULL; }
    
    DomainMonitor(BC.export_OBC(), &C);
    
    //if ( C.isHeatProblem() ) BC.InnerTBC_Periodic()
    
  }

  
  // 初期解およびリスタート解の同期
  if ( numProc > 1 )
  {
    if ( paraMngr->BndCommV3D(d_v, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
    if ( paraMngr->BndCommS3D(d_p, size[0], size[1], size[2], guide, 1    ) != CPM_SUCCESS ) Exit(0);
    
    if ( C.isHeatProblem() ) 
    {
      if ( paraMngr->BndCommS3D(d_p, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
    }
  }

  // VOF
  if ( C.BasicEqs == INCMP_2PHASE )
  {
    setVOF();
    if ( numProc > 1 )
    {
      if ( paraMngr->BndCommS3D(d_vof, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
    }
  }
  
}



// #################################################################
// midの情報から各BCコンポーネントのローカルなインデクスを取得する
// 計算内部領域の境界と外部境界とでは，ガイドセル部分にあるコンポーネントIDの取り扱いが異なる
// 外部境界に接する面では，幅はそのまま，始点はガイドセル部分を含む
// 内部境界に接する面では，始点と幅はローカルノード内の計算内部領域に含まれるように調整
void FFV::setLocalCmpIdx_Binary()
{
  int st_i[3], len[3];
  int id;
  int m_st[3], m_ed[3];
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  for (int m=1; m<=C.NoBC; m++) {
    id = cmp[m].getMatOdr();
    
    switch ( cmp[m].getType() ) 
    {
      case HEX:
      case FAN:
        break; // 体積率でエンコード済み
        
      default:
        
        for (int d=0; d<3; d++) 
        {
          st_i[d] = 0;
          len[d] = 0;
        }
        
        // コンポーネント範囲
        //GetBndIndexExtGc()は自ノード内でのidのバウンディングボックスを取得．インデクスはローカルインデクスで，ガイドセルを含む配列の基点をゼロとするCのインデクス
        if ( !paraMngr->GetBndIndexExtGc(id, d_mid, ix, jx, kx, gd, st_i[0], st_i[1], st_i[2], len[0], len[1], len[2]) ) 
        {
          Hostonly_ stamped_printf("\tError : can not get component local index for ID[%d]\n", id);
          Exit(0);
        }
        
        // ノード内にコンポーネントがあるかどうかをチェック
        if ( (len[0]==0) || (len[1]==0) || (len[2]==0) ) { // コンポーネントなし
          cmp[m].setEns(OFF);
          // BV情報はCompoListのコンストラクタでゼロに初期化されているので，すべてゼロ
        }
        else 
        {
          
          for (int d=0; d<3; d++) 
          {
            int tmp_st=0;
            int tmp_ed=0;
            
            EnlargeIndex(tmp_st, tmp_ed, st_i[d], len[d], size[d], d, id);
            
            m_st[d] = tmp_st;
            m_ed[d] = tmp_ed;
          }
          cmp[m].setBbox(m_st, m_ed);
          cmp[m].setEns(ON); // コンポーネントあり
        }
        
        break;
    }
    
  }
}


// #################################################################
// ParseMatクラスをセットアップし，媒質情報を入力ファイルから読み込み，媒質リストを作成する
void FFV::setMediumList(FILE* fp)
{
  if ( !mat ) Exit(0);
  
  if ( !M.makeMediumList(mat, C.NoMedium) ) {
    Hostonly_ stamped_printf("Error : Duplicate label in Material Table\n");
  }
  
  // 媒質テーブルの表示
  Hostonly_ {
    M.printMatList(stdout, mat, C.NoMedium);
    M.printMatList(fp, mat, C.NoMedium);
  }
}




// #################################################################
// 各種例題のモデルをセット
void FFV::setModel(double& PrepMemory, double& TotalMemory, FILE* fp)
{
  // d_midをゼロで初期化
  size_t mt = (size[0]+2*guide) * (size[1]+2*guide) *(size[2]+2*guide) * sizeof(int);
  memset(d_mid, 0, mt);
  
  
  switch (C.Mode.Example)
  {
    case id_Polygon: // ユーザ例題
      
      C.get_Geometry( M.export_MTI() );
      
      
      // PolylibとCutlibのセットアップ
      setup_Polygon2CutInfo(PrepMemory, TotalMemory, fp);
      break;
      
    case id_Sphere:
      if ( !C.isCDS() ) // binary
      {
        Ex->setup(d_mid, &C, G_origin, C.NoMedium, mat);
      }
      else // cut
      {
        // 初期値1.0をセット
        setup_CutInfo4IP(PrepMemory, TotalMemory, fp);
        Ex->setup_cut(d_mid, &C, G_origin, C.NoMedium, mat, d_cut);
      }
      break;
      
    case id_SHC1D:
      setup_CutInfo4IP(PrepMemory, TotalMemory, fp);
      Ex->setup(d_mid, &C, G_origin, C.NoMedium, mat);
      Ex->setup_bc(d_bid);
      break;
      
    default: // ほかのIntrinsic problems
      if ( C.isCDS() ) // カットの場合
      {
        setup_CutInfo4IP(PrepMemory, TotalMemory, fp);
      }
      Ex->setup(d_mid, &C, G_origin, C.NoMedium, mat);
      break;
  }
  
  // midのガイドセル同期
  if ( paraMngr->BndCommS3D(d_mid, size[0], size[1], size[2], guide, 1) != CPM_SUCCESS ) Exit(0);
  
}



// #################################################################
// 並列化と分割の方法を保持
string FFV::setParallelism()
{
  string para_mode;
  
  C.num_thread  = omp_get_max_threads();
  
  // Serial or Parallel environment
  if( paraMngr->IsParallel() )
  {
    C.num_process = paraMngr->GetNumRank();
    
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
// 時間積分幅や物理パラメータの設定
void FFV::setParameters()
{
  // 無次元数などの計算パラメータを設定する．MediumListを決定した後，かつ，SetBC3Dクラスの初期化前に実施すること
  // 代表物性値をRefMatの示す媒質から取得
  // Δt=constとして，無次元の時間積分幅 deltaTを計算する．ただし，一定幅の場合に限られる．不定幅の場合には別途考慮の必要
  DT.set_Vars(C.KindOfSolver, C.Unit.Param, (double)deltaX, (double)C.Reynolds, (double)C.Peclet);
  
  
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
  
  
  /* Interval Managerの計算の主管理タグ[tg_compute]に値を初期値を設定
  if ( !C.Interval[Interval_Manager::tg_compute].initTrigger(0, 0.0, DT.get_DT(), Interval_Manager::tg_compute, 
                                                             (double)(C.RefLength/C.RefVelocity)) ) 
  {
    Hostonly_ printf("\t Error : Computation Period is asigned to zero.\n");
    Exit(0);
  }
  
  Session_LastStep = C.Interval[Interval_Manager::tg_compute].getIntervalStep();
  */
  
  C.setParameters(mat, cmp, &RF, BC.export_OBC());
  
  
  // パラメータの無次元化（正規化）に必要な参照物理量の設定
  B.setRefMediumProperty(mat, cmp, C.RefMat);
}



// #################################################################
// IP用にカット領域をアロケートする
void FFV::setup_CutInfo4IP(double& m_prep, double& m_total, FILE* fp)
{
  Hostonly_ 
  {
    fprintf(fp, "\n---------------------------------------------------------------------------\n\n");
    fprintf(fp, "\t>> Cut Info\n\n");
    printf("\n---------------------------------------------------------------------------\n\n");
    printf("\t>> Cut Info\n\n");
  }
  
  size_t n_cell[3];
  
  for (int i=0; i<3; i++) 
  {
    n_cell[i] = (size_t)(size[i] + 2*guide); // 分割数+ガイドセル
  }
  size_t size_n_cell = n_cell[0] * n_cell[1] * n_cell[2];
  
  
  // カット情報保持領域のアロケート
  TIMING_start(tm_init_alloc);
  allocArray_Cut(m_total);
  TIMING_stop(tm_init_alloc);
  
  // 使用メモリ量　
  double cut_mem, G_cut_mem;
  G_cut_mem = cut_mem = (double)size_n_cell * (double)(6*sizeof(float) + sizeof(int));
  m_prep += cut_mem;
  m_total+= cut_mem;
  
  
  display_memory_info(fp, G_cut_mem, cut_mem, "Cut");
  

  
  // 初期値のセット
  for (size_t i=0; i<size_n_cell*6; i++) {
    d_cut[i] = 1.0f;
  }
  
  for (size_t i=0; i<size_n_cell; i++) {
    d_bid[i] = 0;
  }
  
}


// #################################################################
// 幾何形状情報を準備し，交点計算を行う
void FFV::setup_Polygon2CutInfo(double& m_prep, double& m_total, FILE* fp)
{
  unsigned poly_gc[3];
  float poly_org[3];
  float poly_dx[3];
  poly_gc[0]  = poly_gc[1] = poly_gc[2] = (unsigned)guide;
  
  // 有次元に変換 Polylib: 並列計算領域情報　ポリゴンは実スケールで，ガイドセル領域部分も含めて指定する
  poly_dx[0]  = (float)pitch[0] *C.RefLength;
  poly_dx[1]  = (float)pitch[1] *C.RefLength;
  poly_dx[2]  = (float)pitch[2] *C.RefLength;
  poly_org[0] = (float)origin[0]*C.RefLength;
  poly_org[1] = (float)origin[1]*C.RefLength;
  poly_org[2] = (float)origin[2]*C.RefLength;
  
  Hostonly_
  {
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
    fprintf(fp,"\t>> Polylib configuration\n\n");
    fprintf(fp,"\tfile name = %s\n\n", C.PolylibConfigName.c_str());
    
    printf("\n---------------------------------------------------------------------------\n\n");
    printf("\t>> Polylib configuration\n\n");
    printf("\tfile name = %s\n\n", C.PolylibConfigName.c_str());
  }
  
  // Polylib: インスタンス取得
  PL = MPIPolylib::get_instance();
  
  // Polylib: ポリゴングループのFactoryクラスを登録
  //PL->set_factory( new MyPolygonGroupFactory() );
  
  // Polylib: 並列計算領域情報を設定
  poly_stat = PL->init_parallel_info( MPI_COMM_WORLD,
                                     poly_org,        // 自ランクの基点座標
                                     (unsigned*)size, // 自ランクの分割数
                                     poly_gc,         // ガイドセル数
                                     poly_dx          // 格子幅
                                     );
  if ( poly_stat != PLSTAT_OK )
  {
    Hostonly_
    {
      fprintf(fp,"\tRank [%6d]: p_polylib->init_parallel_info() failed.", myRank);
      printf    ("\tRank [%6d]: p_polylib->init_parallel_info() failed.", myRank);
    }
    Exit(0);
  }
  
  
  
  // Polylib: STLデータ読み込みとスケーリング
  TIMING_start(tm_polygon_load);
  
  // ロード
  poly_stat = PL->load_rank0( C.PolylibConfigName );
  
  if( poly_stat != PLSTAT_OK )
  {
    Hostonly_
    {
      printf    ("\tRank [%6d]: p_polylib->load_rank0() failed.", myRank);
      fprintf(fp,"\tRank [%6d]: p_polylib->load_rank0() failed.", myRank);
    }
    Exit(0);
  }
  
  /* スケーリングする場合のみ表示
  if ( C.Scaling_Factor != 1.0 )
  {
    fprintf(fp,"\n\tScaling Factor           :   %f\n", C.Scaling_Factor);
  }
  
  // スケーリング
  poly_stat = PL->rescale_polygons(C.Scaling_Factor);
  
  if( poly_stat != PLSTAT_OK )
  {
    Hostonly_
    {
      printf    ("\tRank [%6d]: p_polylib->rescale_polygons() failed.", myRank);
      fprintf(fp,"\tRank [%6d]: p_polylib->rescale_polygons() failed.", myRank);
    }
    Exit(0);
  }*/
  
  
  TIMING_stop(tm_polygon_load);
  
  // 階層情報表示 debug
// ##########
#if 0
  PL->show_group_hierarchy();
  PL->show_group_hierarchy(fp);
#endif
// ##########
  
  
  // ポリゴン情報へのアクセス
  vector<PolygonGroup*>* pg_roots = PL->get_root_groups();
  vector<PolygonGroup*>::iterator it;
  
  // Polygon Groupの数
  C.num_of_polygrp = pg_roots->size();
  
  // ポリゴングループの属性保持のstruct array
  poly_prop = new Control::Polygon_property[C.num_of_polygrp];
  
  Hostonly_
  {
    printf(     "\tNumber of Polygon Group = %d\n\n", C.num_of_polygrp);
    fprintf(fp, "\tNumber of Polygon Group = %d\n\n", C.num_of_polygrp);
    
    printf(     "\t Medium ID         Material :          No. :  Polygon Group Label :         Area\n");
    fprintf(fp, "\t        ID         Material :          No. :  Polygon Group Label :         Area\n");
  }
  
  
  // ポリゴン情報の表示
  int c = 0;
  for (it = pg_roots->begin(); it != pg_roots->end(); it++) {
    std::string m_pg = (*it)->get_name();
    int m_id = (*it)->get_id();
    int ntria= (*it)->get_group_num_tria();
    std::string m_mat = (*it)->get_label();
    float area = (*it)->get_group_area();
    
    poly_prop[c].mat_id    = m_id;  // polylib.tpのID
    poly_prop[c].label_grp = m_pg;  // グループラベル
    poly_prop[c].label_mat = m_mat; // 媒質ラベル
    poly_prop[c].area      = area;  // ポリゴンの総面積
    c++;

    if ( numProc > 1 )
    {
      int tmp = ntria;
      if ( paraMngr->Allreduce(&tmp, &ntria, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
      
      float ta = area;
      if ( paraMngr->Allreduce(&ta, &area, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
    }
    
    Hostonly_
    {
      printf(    "\t %9d %16s : %12d : %20s : %e\n", m_id, m_mat.c_str(), ntria, m_pg.c_str(), area);
      fprintf(fp,"\t %9d %16s : %12d : %20s : %e\n", m_id, m_mat.c_str(), ntria, m_pg.c_str(), area);
    }
    
// ##########
#if 0
    PL->show_group_info(m_pg); //debug
#endif
// ##########
    
  }
  
  delete pg_roots;
  
  Hostonly_
  {
    printf("\n");
    fprintf(fp, "\n");
  }
  
  
  /* PolygonGroupの媒質数は，1以上、MediumTableの数以下であること
  if ( (C.num_of_polygrp < 1) || (C.num_of_polygrp > C.NoMedium) )
  {
    Hostonly_
    {
      printf("\tError : No of PolygonGroup must be less than one of Medium.\n");
      Exit(0);
    }
  }*/
  
  
  
  // 使用メモリ量　基本クラスのみ
  double poly_mem, G_poly_mem;
  G_poly_mem = poly_mem = (double)PL->used_memory_size();
  m_prep += poly_mem;
  m_total+= poly_mem;
  
  display_memory_info(fp, G_poly_mem, poly_mem, "Polygon");
  

  
  
  // Triangle display >> Debug
// ##########
#if 0
  PolylibNS::Vec3f m_min, m_max, t1(poly_org), t2(poly_dx), t3;
  t3.assign((float)size[0]*t2.t[0], (float)size[1]*t2.t[1], (float)size[2]*t2.t[2]);
  m_min = t1 - t2;      // 1層外側まで
  m_max = t1 + t3 + t2; //
  printf("min : %f %f %f\n", m_min.t[0], m_min.t[1], m_min.t[2]);
  printf("max : %f %f %f\n", m_max.t[0], m_max.t[1], m_max.t[2]);
  vector<Triangle*>* trias = PL->search_polygons("Tube", m_min, m_max, false); // false; ポリゴンが一部でもかかる場合
  
  PolylibNS::Vec3f *p, nrl, n;
  int c=0;
  vector<Triangle*>::iterator it2;
  for (it2 = trias->begin(); it2 != trias->end(); it2++) {
    p = (*it2)->get_vertex();
    n = (*it2)->get_normal();
    printf("%d : p0=(%6.3e %6.3e %6.3e)  p1=(%6.3e %6.3e %6.3e) p2=(%6.3e %6.3e %6.3e) n=(%6.3e %6.3e %6.3e)\n", c++,
           p[0].t[0], p[0].t[1], p[0].t[2],
           p[1].t[0], p[1].t[1], p[1].t[2],
           p[2].t[0], p[2].t[1], p[2].t[2],
           n.t[0], n.t[1], n.t[2]);
  }
  
  delete trias;  //後始末
#endif
// ##########
  
  // Polylib: STLデータ書き出しテスト >> Debug
// ##########
#if 0
  unsigned poly_out_para = IO_DISTRIBUTE; // 逐次の場合と並列の場合で明示的に切り分けている．あとで，考慮
  string fname;
  
  if ( poly_out_para == IO_GATHER )
  {
    poly_stat = PL->save_rank0( &fname, "stl_b" );
    if ( poly_stat != PLSTAT_OK )
    {
      Hostonly_
      {
        printf("Rank [%d]: p_polylib->save_rank0() failed to write into '%s'.", myRank, fname.c_str());
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
        printf("Rank [%d]: p_polylib->save_parallel() failed to write into '%s'.", myRank, fname.c_str());
        fprintf(fp,"Rank [%d]: p_polylib->save_parallel() failed to write into '%s'.", myRank, fname.c_str());
      }
      Exit(0);
    }
  }
#endif
// ##########  
  
  
  // Cutlib
  Hostonly_
  {
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
    fprintf(fp,"\t>> Cut Info\n\n");
    printf("\n---------------------------------------------------------------------------\n\n");
    printf("\t>> Cut Info\n\n");
  }
  
  size_t n_cell[3];

  
  for ( int i=0; i<3; i++) {
    n_cell[i] = (size_t)(size[i] + 2*guide); // 分割数+ガイドセル両側
    poly_org[i] -= poly_dx[i]*(float)guide;  // ガイドセル分だけシフト
  }
  
  size_t size_n_cell = n_cell[0] * n_cell[1] * n_cell[2];
  
  if (size_n_cell*6  > UINT_MAX)
  {
    Hostonly_
    {
      stamped_printf("\n\tError : Product of size[]*6 exceeds UINT_MAX\n\n");
      stamped_fprintf(fp,"\n\tError : Product of size[]*6 exceeds UINT_MAX\n\n");
    }
    Exit(0);
  }
  
  
  // Cutlibの配列は各方向(引数)のサイズ
  TIMING_start(tm_init_alloc);
  cutPos = new CutPos32Array(n_cell); // 6*(float)
  cutBid = new CutBid5Array(n_cell);  // (int32_t)
  TIMING_stop(tm_init_alloc);
  
  TIMING_start(tm_cutinfo);
  CutInfoCell(poly_org, poly_dx, PL, cutPos, cutBid); // ガイドセルを含む全領域で交点を計算
  TIMING_stop(tm_cutinfo);
  
  
  // 使用メモリ量
  double cut_mem, G_cut_mem;
  G_cut_mem = cut_mem = (double)size_n_cell * (double)(6*sizeof(float) + sizeof(int)); // float
  m_prep += cut_mem;
  m_total+= cut_mem;
  
  
  display_memory_info(fp, G_cut_mem, cut_mem, "Cut");
  
  
  // カットとID情報をポイント
  d_cut = (float*)cutPos->getDataPointer();
  d_bid = (int*)cutBid->getDataPointer();
  

  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  // カットの最小値を求める
  float f_min=1.0;
  
#pragma omp parallel firstprivate(ix, jx, kx, gd)
  {
    float th_min = 1.0;
    
#pragma omp parallel for schedule(static)
    for (int k=1; k<=kx; k++) {
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          
          size_t mp = _F_IDX_S4DEX(0, i, j, k, 6, ix, jx, kx, gd);
          size_t mb = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
          int bd = d_bid[mb];
          
          if ( TEST_BC(bd) ) // カットがあるか，IDによる判定
          //if ( (d_cut[mp+0]+d_cut[mp+1]+d_cut[mp+2]+d_cut[mp+3]+d_cut[mp+4]+d_cut[mp+5]) < 6.0 ) // 距離による判定
          {
            for (size_t n=0; n<6; n++) {
              
              th_min = min(th_min, d_cut[mp+n]);
              
            }
            
// ##########            
#if 0 // debug
            int b0 = (bd >> 0)  & MASK_5;
            int b1 = (bd >> 5)  & MASK_5;
            int b2 = (bd >> 10) & MASK_5;
            int b3 = (bd >> 15) & MASK_5;
            int b4 = (bd >> 20) & MASK_5;
            int b5 = (bd >> 25) & MASK_5;
            printf("%3d %3d %3d : %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %d %d %d %d %d %d\n", i,j,k,
                   d_cut[mp+0],
                   d_cut[mp+1],
                   d_cut[mp+2],
                   d_cut[mp+3],
                   d_cut[mp+4],
                   d_cut[mp+5],
                   b0, b1, b2, b3, b4, b5);
#endif
// ##########            
          }
        }
      }
    }

#pragma omp critical
    {
      f_min = min(f_min, th_min);
    }
  }

  
  if ( numProc > 1 )
  {
    float tmp = f_min;
    if ( paraMngr->Allreduce(&tmp, &f_min, 1, MPI_MIN) != CPM_SUCCESS ) Exit(0);
  }
  
  Hostonly_
  {
    printf(    "\n\tMinimum dist. = %5.3e  \n", f_min);
    fprintf(fp,"\n\tMinimum dist. = %5.3e  \n", f_min);
  }
  
  // カットの最小値
  min_distance(d_cut, fp);
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
// BCIndexにビット情報をエンコードする
void FFV::VoxEncode()
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  // 基本ビット情報（Active, State, コンポ，媒質情報）を全領域についてエンコードする
  V.setBCIndex_base1(d_bcd, d_mid, d_cvf, mat, cmp);

  // bcdの同期
  if ( paraMngr->BndCommS3D(d_bcd, ix, jx, kx, gd, 1) != CPM_SUCCESS ) Exit(0);
  


  V.setBCIndex_base2(d_bcd, d_mid, L_Acell, G_Acell, C.KindOfSolver, cmp);

  // STATEとACTIVEビットのコピー
  V.copyBCIbase(d_bcp, d_bcd);
  V.copyBCIbase(d_bcv, d_bcd);
  
  if ( C.isHeatProblem() ) 
  {
    V.copyBCIbase(d_bh1, d_bcd);
    V.copyBCIbase(d_bh2, d_bcd);
  }
  


  // BCIndexP に圧力計算のビット情報をエンコードする -----
  if ( C.isCDS() ) 
  {
    C.NoWallSurface = V.setBCIndexP(d_bcd, d_bcp, d_mid, &BC, cmp, d_cut, d_bid, true);
  }
  else // binary
  {
    C.NoWallSurface = V.setBCIndexP(d_bcd, d_bcp, d_mid, &BC, cmp, d_cut, d_bid, false);
  }

#if 0
  V.dbg_chkBCIndexP(d_bcd, d_bcp, "BCindexP.txt", cmp);
#endif
  
  
  
  // BCIndexV に速度計算のビット情報をエンコードする -----
  if ( C.isCDS() ) 
  {
    V.setBCIndexV(d_bcv, d_mid, d_bcp, &BC, cmp, true, d_cut, d_bid);
  }
  else // binary
  {
    V.setBCIndexV(d_bcv, d_mid, d_bcp, &BC, cmp);
  }
  


// ##########
#if 0
  V.dbg_chkBCIndexV(d_bcv, "BCindexV.txt");
#endif
// ##########
  
  // BCIndexT に温度計算のビット情報をエンコードする -----
  if ( C.isHeatProblem() )
  {
    if ( C.isCDS() )
    {
      V.setBCIndexH(d_bcd, d_bh1, d_bh2, d_mid, &BC, C.KindOfSolver, cmp, true, d_cut, d_bid);
    }
    else // binary
    {
      V.setBCIndexH(d_bcd, d_bh1, d_bh2, d_mid, &BC, C.KindOfSolver, cmp);
    }
    
    
// ##########
#if 0
    V.dbg_chkBCIndexH(d_bcv, "BCindexH.txt");
#endif
// ##########
    
  }

  // 内部領域のFluid, Solidのセル数を数える C.Wcell(Local), G_Wcell(global)
  V.countCellState(L_Wcell, G_Wcell, d_bcd, SOLID);
  V.countCellState(L_Fcell, G_Fcell, d_bcd, FLUID);
  

  // getLocalCmpIdx()などで作成したコンポーネントのインデクスの再構築
  resizeCompoBV(C.KindOfSolver, C.isHeatProblem());
}



// #################################################################
// ボクセルをスキャンし情報を表示する
void FFV::VoxScan(FILE* fp)
{
  // 外部境界面の媒質IDとその個数を取得
  int cell_id[NOFACE];
  
  for (int i=0; i<NOFACE; i++) {
    cell_id[i] = BC.export_OBC(i)->get_GuideMedium();
  }

// ##########
#if 0
  Hostonly_ 
  {
    fprintf(fp, "\tCell IDs on Guide cell region\n");
    
    for ( int i=0; i<NOFACE; i++) {
      fprintf(fp, "\t\t%s = %d\n", FBUtility::getDirection(i).c_str(), cell_id[i]);
    }
    fprintf(stdout, "\tCell IDs on Guide cell region\n");
    
    for ( int i=0; i<NOFACE; i++) {
      fprintf(stdout, "\t\t%s = %d\n", FBUtility::getDirection(i).c_str(), cell_id[i]);
    }
  }
#endif
// ##########
  
  // midに設定されたIDをスキャンし，IDの個数を返し，作業用のcolorList配列にIDを保持，midに含まれるIDの数をチェック
  int sc=0;
  if ( (sc = V.scanCell(d_mid, cell_id, C.Hide.Change_ID)) > C.NoMedium )
  {
    Hostonly_
    {
      stamped_printf(    "A number of Mediums included in model(%d) is greater than the one in 'Medium_List'(%d)\n", sc, C.NoMedium);
      stamped_fprintf(fp,"A number of Mediums included in model(%d) is greater than the one in 'Medium_List'(%d)\n", sc, C.NoMedium);
    }
    Exit(0);
  }
}
