// #################################################################
//
// FFV : Frontflow / violet
//
// Copyright (c) All right reserved. 2012
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
//
// #################################################################

/** 
 * @file   ffv_Initialize.C
 * @brief  FFV Class
 * @author kero
 */

#include "ffv.h"


int FFV::Initialize(int argc, char **argv)
{
  double TotalMemory   = 0.0;  ///< 計算に必要なメモリ量（ローカル）
  double PrepMemory    = 0.0;  ///< 初期化に必要なメモリ量（ローカル）
  double G_TotalMemory = 0.0;  ///< 計算に必要なメモリ量（グローバル）
  double G_PrepMemory  = 0.0;  ///< 初期化に必要なメモリ量（グローバル）
  double tmp_memory    = 0.0;  ///< 計算に必要なメモリ量（グローバル）？
  
  double flop_task     = 0.0;  ///< flops計算用
  
  // CPMバージョン表示
  if ( paraMngr->GetMyRankID() == 0 )
  {
    cpm_Base::VersionInfo();
  }
  
  
  // 固定パラメータ
  fixed_parameters();
  

  
  // ------------------------------------
  FILE* fp = NULL;
  
  // condition fileのオープン
  Hostonly_ {
    if ( !(fp=fopen("condition.txt", "w")) ) {
      stamped_printf("\tSorry, can't open 'condition.txt' file. Write failed.\n");
      return -1;
    }
  }
  
  // メッセージ表示
  Hostonly_ {
    FBUtility::printVersion(fp, "Welcome to FFV  ", FFV_VERS);
    FBUtility::printVersion(stdout, "Welcome to FFV  ", FFV_VERS);
    
    FBUtility::printVersion(fp, "FlowBase        ", FB_VERS);
    FBUtility::printVersion(stdout, "FlowBase        ", FB_VERS);
  }
  
  // 入力ファイルの指定
  string input_file = argv[1];
  
  int ierror;
  
  // TextParserのインスタンス生成
  ierror = tpCntl.getTPinstance();
  
  // TextParserのファイル読み込み
  ierror = tpCntl.readTPfile(input_file);

  
  // TPControlクラスのポインタをControlクラスに渡す
  C.importTP(&tpCntl);
  
  
  // 例題の種類を取得し，C.Mode.Exampleにフラグをセットする
  getExample(&C, &tpCntl);
  
  // 組み込み例題クラスの実体をインスタンスし，*Exにポイントする
  connectExample(&C);
  
  
  // 組み込み例題クラス名を表示
  Hostonly_ {
    Ex->printExample(fp, Ex->getExampleName());
  }
  
  // バージョン情報を取得し，ソルバークラスのバージョンと一致するかをチェックする
  C.get_Version();
  
  if ( C.version != FFV_VERS ) {
    Hostonly_ {
      fprintf(stdout, "\t##### Version of Input file (%d) is NOT compliant with FFV ver. %d #####\n", C.version, FFV_VERS);
      fprintf(fp, "\t##### Version of Input file (%d) is NOT compliant with FFV ver. %d #####\n", C.version, FFV_VERS);
    }
    return -1;
  }
  if ( C.FB_version != FB_VERS ) {
    Hostonly_ {
      fprintf(stdout, "\t##### Version of Input file (%d) is NOT compliant with FB ver. %d #####\n", C.FB_version, FB_VERS);
      fprintf(fp, "\t##### Version of Input file (%d) is NOT compliant with FB ver. %d #####\n", C.FB_version, FB_VERS);
    }
    return -1;
  }
  
  
  // 最初のパラメータの取得
  C.get_Steer_1(&DT);

  
  // 一度、テキストパーサーのDBを破棄
  if (tpCntl.remove() != TP_NO_ERROR ) 
  {
    Hostonly_ printf("Error : delete textparser\n");
    Exit(0);
  }
  
  
  // FiliIOのモードを修正
  if ( myRank == 1 ) {
    C.FIO.IO_Input  = IO_GATHER;
    C.FIO.IO_Output = IO_GATHER;
    Hostonly_ printf("\tMode of Parallel_Input/_Output was changed because of serial execution.\n");
  }
  
  
  // ###########################
  // 計算領域全体のサイズ，並列計算時のローカルのサイズ，コンポーネントのサイズなどを設定
  
  string dom_file = argv[2]; // 領域情報を記述したファイル名の取得 >> テキストパーサのDB切り替え
  
  // 領域設定
  DomainInitialize(dom_file);
  
  // 各クラスで領域情報を保持
  C.setDomainInfo(paraMngr, procGrp);   C.setNeighborInfo(C.guide);
  B.setDomainInfo(paraMngr, procGrp);   B.setNeighborInfo(C.guide);
  V.setDomainInfo(paraMngr, procGrp);   V.setNeighborInfo(C.guide);
  F.setDomainInfo(paraMngr, procGrp);   F.setNeighborInfo(C.guide);
  BC.setDomainInfo(paraMngr, procGrp);  BC.setNeighborInfo(C.guide);
  Ex->setDomainInfo(paraMngr, procGrp); Ex->setNeighborInfo(C.guide);
  
  // ###########################

  
  // 各例題のパラメータ設定 -----------------------------------------------------
  Ex->setDomain(&C, size, origin, region, pitch);
  
  
  // 　再度、入力ファイルをオープン
  ierror = tpCntl.readTPfile(input_file);

  // TPControlクラスのポインタを各クラスに渡す(2回目)
  C.importTP(&tpCntl);
  B.importTP(&tpCntl);
  M.importTP(&tpCntl);
  
  // パラメータを取得
  C.get_Steer_2(IC, &RF);
  
  
  // 組み込み例題の固有パラメータ
  if ( !Ex->getTP(&C, &tpCntl) ) Exit(0);
  
  
  // 媒質情報をパラメータファイルから読み込み，媒質リストを作成する
  Hostonly_  {
    fprintf(stdout,"\n---------------------------------------------------------------------------\n\n");
    fprintf(stdout,"\n\t>> Medium List\n\n");
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
    fprintf(fp,"\n\t>> Medium List\n\n");
  }
  
  // 媒質情報をロードし、 Medium_Tableタグ内の媒質数を保持
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
  
  B.setControlVars(&C, BC.export_OBC(), mat);
  
  B.setMediumTI( M.export_MTI() );
  
  B.countMedium(&C);
  
  // CompoListクラスをインスタンス．[0]はダミーとして利用しないので，配列の大きさはプラス１する
  cmp = new CompoList[C.NoCompo+1];

  

  
  // 並列処理モード
  string para_label = setParallelism();
  
  // タイミング測定の初期化
  if ( C.Mode.Profiling != OFF ) {
    ModeTiming = ON;
    TIMING__ PM.initialize( tm_END );
    TIMING__ PM.setRankInfo( paraMngr->GetMyRankID() );
    TIMING__ PM.setParallelMode(para_label, C.num_thread, C.num_process);
    set_timing_label();
  }
  
  // タイミング測定開始
  TIMING_start(tm_init_sct); 
  
  // 前処理に用いるデータクラスのアロケート -----------------------------------------------------
  TIMING_start(tm_init_alloc); 
  allocArray_Prep(PrepMemory, TotalMemory);
  TIMING_stop(tm_init_alloc);
  
  
  // ファイルからIDを読み込む，または組み込み例題クラスでID情報を作成
  Hostonly_ {
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
    fprintf(fp,"\t>> Voxel file information\n\n");
    fprintf(stdout,"\n---------------------------------------------------------------------------\n\n");
    fprintf(stdout,"\t>> Voxel file information\n\n");
  }
  
  TIMING_start(tm_voxel_prep_sct);
  
  
  // 各問題に応じてモデルを設定
  setModel(PrepMemory, TotalMemory, fp);
  
  
  
  // 領域情報の表示
  Hostonly_ {
    fprintf(stdout,"\n---------------------------------------------------------------------------\n");
    fprintf(stdout,"\n\t>> Global Domain Information\n\n");
    C.printGlobalDomain(stdout, G_size, G_origin, G_region, pitch);
    
    fprintf(fp,"\n---------------------------------------------------------------------------\n");
    fprintf(fp,"\n\t>> Global Domain Information\n\n");
    C.printGlobalDomain(fp, G_size, G_origin, G_region, pitch);
  }
  
  // メモリ消費量の情報を表示
  Hostonly_ {
    fprintf(stdout,"\n---------------------------------------------------------------------------\n\n");
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
  }
  G_PrepMemory = PrepMemory;
  
  if( numProc > 1 ) {
    tmp_memory = G_PrepMemory;
    if ( paraMngr->Allreduce(&tmp_memory, &G_PrepMemory, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
  }
  Hostonly_  {
    FBUtility::MemoryRequirement("prep", G_PrepMemory, PrepMemory, stdout);
    FBUtility::MemoryRequirement("prep", G_PrepMemory, PrepMemory, fp);
  }
  
  
  // Fill
  if ( (C.Mode.Example == id_Polygon) ) { //&& C.isCDS() ) {
    
    Hostonly_ {
      fprintf(stdout,"\n---------------------------------------------------------------------------\n\n");
      fprintf(stdout,"\t>> Fill\n\n");
      fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
      fprintf(fp,"\t>> Fill\n\n");
    }
    
    // Fill(&Vinfo);
  }
  
  
  
  
  // ∆tの決め方とKindOfSolverの組み合わせで無効なものをはねる
  if ( !DT.chkDtSelect() ) {
    Hostonly_ printf("\tCombination of specified 'Time_Increment' and 'Kind_of_Solver' is not permitted.\n");
    return -1;
  }
  
  
  // パラメータファイルから得られた内部BCコンポーネント数を表示
  Hostonly_ {
    fprintf(stdout,"\n---------------------------------------------------------------------------\n\n");
    fprintf(stdout,"\t>> Components\n\n");
    C.printNoCompo(stdout);
    fprintf(stdout,"\n"); fflush(stdout);
    
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
    fprintf(fp,"\t>> Components\n\n");
    C.printNoCompo(fp);
    fprintf(fp,"\n"); fflush(fp);
  }
  
  
  // パラメータファイルをパースして，外部境界条件を保持する　>> VoxScan()につづく
  B.loadBC_Outer();
  
  // ボクセルのスキャン
  VoxScan(fp);
  
  
  // スキャンしたセルIDの情報を表示する
  Hostonly_ {
    fprintf(stdout, "\n---------------------------------------------------------------------------\n\n");
    fprintf(stdout,"\t>> Information of Scaned Voxel\n\n");
    V.printScanedCell(stdout);
		fflush(stdout);
    
    fprintf(fp, "\n---------------------------------------------------------------------------\n\n");
    fprintf(fp,"\t>> Information of Scaned Voxel\n\n");
		V.printScanedCell(fp);
		fflush(fp);
  }
  
  // ボクセルモデルの媒質インデクスがパラメータファイルに記述された媒質インデクスに含まれていること
	Hostonly_ {
		if ( !V.chkIDconsistency(C.NoMedium) ) {
			stamped_printf("\tID in between XML and Voxel scaned is not consistent\n");
			return -1;
		}
	}

  
  // CompoList, MediumListのポインタをセット
  BC.importCMP_MAT(cmp, mat);
  B.importCompoPtr(cmp);
  
  // CompoListの設定，外部境界条件の読み込み保持、ガイドセル上にパラメータファイルで指定する媒質インデクスを代入
  setBCinfo();
  

  
  
  // Cell_Monitorの指定がある場合，モニタ位置をセット
  if ( (C.Sampling.log == ON) && (C.isMonitor() == ON) ) 
  {
    // ShapeMonitorのインスタンス
    ShapeMonitor SM(size, guide, pitch, origin);
    
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
  if ( C.isVfraction() ) {
    TIMING_start(tm_init_alloc);
    allocArray_CompoVF(PrepMemory, TotalMemory);
    TIMING_stop(tm_init_alloc); 
    
    setComponentVF();
  }
  
  // コンポーネントのローカルインデクスを保存
  setLocalCmpIdx_Binary();

  
  // 内部周期境界の場合のガイドセルのコピー処理
  V.adjMediumPrdc_Inner(d_mid, cmp);
  
  // 媒質数とKindOfSolverの整合性をチェックする
  if ( !chkMediumConsistency() ) {
    Hostonly_ stamped_printf("\tchkMediumConsistency()\n");
    return -1;
  }
  

  // BCIndexへのエンコード処理
  Hostonly_  {
    fprintf(stdout,"\n---------------------------------------------------------------------------\n\n");
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
  }
  
  // CDSのとき，ポリゴンからVBCのコンポーネント情報を設定
  if ( C.isCDS() ) {
    VIBC_Bbox_from_Cut();
  }
  
  // BCIndexにビット情報をエンコードとコンポーネントインデクスの再構築
  VoxEncode();
  

  
  // 体積力を使う場合のコンポーネント配列の確保
  TIMING_start(tm_init_alloc);
  allocArray_Forcing(PrepMemory, TotalMemory, fp);
  TIMING_stop(tm_init_alloc); 
  

  // コンポーネントの体積率を8bitで量子化し，圧力損失コンポの場合にはFORCING_BITをON > bcdにエンコード
  V.setCmpFraction(cmp, d_bcd, d_cvf);

// ########## 
#if 0
  // CompoListとMediumListの関連を表示
  Hostonly_ 
  {
    M.printRelation(stdout, cmp, mat);
    M.printRelation(fp, cmp, mat);
  }
#endif
// ########## 
  
  
  // Ref_MediumがCompoList中にあるかどうかをチェックし、RefMatを設定
  if ( (C.RefMat = C.find_ID_from_Label(mat, C.NoCompo, C.Ref_Medium)) == 0 )
  {
    Hostonly_ {
      fprintf(stdout, "RefMat[%d] is not listed in Medium_Table.\n");
      fprintf(fp, "RefMat[%d] is not listed in Medium_Table.\n");
    }
    return -1;
  }
  

  // 周期境界条件が設定されている場合のBCIndexの周期条件の強制同期
  BC.setBCIperiodic(d_bcd);
  BC.setBCIperiodic(d_bcp);
  BC.setBCIperiodic(d_bcv);
  if ( C.isHeatProblem() ) {
    BC.setBCIperiodic(d_bh1);
    BC.setBCIperiodic(d_bh2);
  }
  
  // bcd/bcp/bcv/bchの同期
  if ( paraMngr->BndCommS3D(d_bcd, size[0], size[1], size[2], guide, 1) != CPM_SUCCESS ) Exit(0);
  if ( paraMngr->BndCommS3D(d_bcp, size[0], size[1], size[2], guide, 1) != CPM_SUCCESS ) Exit(0);
  if ( paraMngr->BndCommS3D(d_bcv, size[0], size[1], size[2], guide, 1) != CPM_SUCCESS ) Exit(0);

  if ( C.isHeatProblem() ) {
    if ( paraMngr->BndCommS3D(d_bh1, size[0], size[1], size[2], guide, 1) != CPM_SUCCESS ) Exit(0);
    if ( paraMngr->BndCommS3D(d_bh2, size[0], size[1], size[2], guide, 1) != CPM_SUCCESS ) Exit(0);
  }
  
  
  // 法線計算
  if ( C.NoBC != 0 ) {
    
    // コンポーネントで指定されるID面の法線を計算，向きはblowing/suctionにより決まる．　bcdをセットしたあとに処理
    for (int n=1; n<=C.NoBC; n++) {
      if ( C.Mode.Example == id_Polygon ) {
        V.get_Compo_Area_Cut(n, cmp, PL);
      }
      else {
        //Vinfo.cal_Compo_Area_Normal(n, bcd, bcv, bh1, C.dh*C.RefLength, &compo_global_bbox[n*6]);
      }
    }
  }
  
  
  // 時間積分幅 deltaT や物理パラメータの設定
  setParameters();

  
  
  // 必要なパラメータをSetBC3Dクラスオブジェクトにコピーする >> setParameters()の後
  BC.setControlVars(&C, mat, cmp, &RF, Ex);
  
  
// ##########
#if 0
  // チェックのため，全計算セルのBCIndexの内容を表示する
  if ( !V.dbg_chkBCIndexP(bcd, bcp, "BCindex.txt") ) {
    Hostonly_ stamped_printf("\tVoxInfo::dbg_chkBCIndexP()\n");
    return -1;
  }
#endif
// ##########
  
  
  
  // 温度計算の場合の初期値指定
  if ( C.isHeatProblem() ) {
    cout <<  "now heat probrem *** initialize" << endl;
    B.get_Medium_InitTemp();
  }
  
  // set phase 
  if ( C.BasicEqs == INCMP_2PHASE ) {
    B.get_Phase();
  }
  
  
  // CompoListの内容とセル数の情報を表示する
  Hostonly_  {
    fprintf(stdout,"\n---------------------------------------------------------------------------\n\n");
    fprintf(stdout,"\t>> Component List\n\n");
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
    fprintf(fp,"\t>> Component List\n\n");
  }
  display_CompoList(fp);
  
  
  // 外部境界面の開口率を計算する
  V.countOpenAreaOfDomain(d_bcd, C.OpenDomain);
  
  
  Hostonly_ {
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n\n");
    fprintf(stdout,"\n---------------------------------------------------------------------------\n\n\n");
  }
  
  
  
  // Monitor Listの処理 --------------------------------------------
  //MO.setControlVars(bcd, G_org, G_reg, C.org, C.dx, C.Lbx, size, guide,
  //                  C.RefVelocity, C.BaseTemp, C.DiffTemp, C.RefDensity, C.RefLength, C.BasePrs,
  //                  C.Unit.Temp, C.Mode.Precision, C.Unit.Prs);
  
  
  // モニタ機能がONの場合に，パラメータを取得し，セットの配列を確保する
  //if ( C.Sampling.log == ON ) getXML_Monitor(m_solvCfg, &MO);
  //if ( C.Sampling.log == ON ) getTP_Monitor(m_solvCfg, &MO);//未完成
  
  
  // モニタリストが指定されている場合に，プローブ位置をID=255としてボクセルファイルに書き込む
  //if (C.Sampling.log == ON ) {
    //MO.write_ID(mid);
  //}
  
  // 内部境界条件として指定されたモニタ設定を登録
  //if ( (C.Sampling.log == ON) && (C.isMonitor() == ON) ) MO.setInnerBoundary(cmp, C.NoBC);
  
  
  // MonitorListの場合に，svxファイルを出力する．
  if ( C.Sampling.log == ON ) {
    Hostonly_ printf("\n\twrite ID which includes Monitor List ID\n\n");
    
    // 性能測定モードがオフのときのみ出力
    if ( C.Hide.PM_Test == OFF ) Ex->writeSVX(d_mid, &C); // writeSVX(); ユーザ問題の場合には，単にtrueを返す
  }
  
  
  
  // mid[]を解放する  ---------------------------
  if ( d_mid ) delete [] d_mid;
  
  
  
  // コンポーネントのグローバルインデクス情報を取得
  setGlobalCmpIdx();

  
  
  // コンポーネントの内容リストを表示し、コンポーネント数がゼロの場合と境界条件との整合性をチェック
  Hostonly_ {
    fprintf(stdout,"\n---------------------------------------------------------------------------\n\n");
    fprintf(stdout,"\t>> Component Information\n");
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
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
  Hostonly_ fprintf(stdout,"\n---------------------------------------------------------------------------\n\n");
  Hostonly_ fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
  
  
  // 制御インターバルの初期化
  init_Interval();
  
  
  
  // 平均値のロード
  if ( C.Start == restart ) {
    TIMING_start(tm_restart);
    if ( C.Mode.Average == ON ) Restart_avrerage(fp, flop_task);
    TIMING_stop(tm_restart);
  }
  
  
  // リスタートの最大値と最小値の表示
  Restart_display_minmax(fp, flop_task);
  
  
  
  // 制御パラメータ，物理パラメータの表示
  Hostonly_ 
  {
    fprintf(stdout,"\n---------------------------------------------------------------------------\n\n");
    fprintf(stdout,"\t>> Outer Boundary Conditions\n\n");
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
    fprintf(fp,"\t>> Outer Boundary Conditions\n\n");
    
    display_Parameters(fp);
  }
  
  
  
  // ドライバ条件のチェック
  BC.checkDriver(fp);
  
  
  // 初期条件の条件設定
  setInitialCondition();
  
  
  // 履歴出力準備
  prep_HistoryOutput();
  
  
  
  // サンプリング元となるデータ配列の登録
  if ( C.Sampling.log == ON ) {
    if ( C.isHeatProblem() ) {
      //MO.setDataPtrs(dc_v->GetData(), dc_p->GetData(), dc_t->GetData());
    }
    else {
      //MO.setDataPtrs(dc_v->GetData(), dc_p->GetData());
    }
  }
  
  
  // 初期状態のファイル出力  リスタート時と性能測定モードのときには出力しない
	if ( (C.Hide.PM_Test == OFF) && (0 == CurrentStep) ) FileOutput(flop_task);
  
  
  // 粗い格子を用いたリスタート時には出力
  if ( C.Start == coarse_restart ) FileOutput(flop_task, true);
  
  
  
  
  
  // SOR2SMA用のバッファ確保
  allocate_SOR2SMA_buffer(TotalMemory);
  
  
  // メモリ使用量の表示
  Hostonly_ {
    fprintf(stdout,"\n---------------------------------------------------------------------------\n\n");
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
  }
  G_TotalMemory = TotalMemory;
  if ( numProc > 1 ) {
    tmp_memory = G_TotalMemory;
    if ( paraMngr->Allreduce(&tmp_memory, &G_TotalMemory, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
  }
  Hostonly_ {
    double n = (double)( (size[0]+2*guide) * (size[1]+2*guide) * (size[2]+2*guide) * 3 );
    double mc = n * (double)sizeof(REAL_TYPE); // temporaty array for vector output, see prepOutput();
    FBUtility::MemoryRequirement("solver", G_TotalMemory, TotalMemory, stdout);
    FBUtility::MemoryRequirement("solver", G_TotalMemory, TotalMemory, fp);
  }
  
  
  
  // 初期化終了時に、入力パラメータのDBを破棄
  if (tpCntl.remove() != TP_NO_ERROR ) 
  {
    Hostonly_ printf("Error : delete textparser\n");
    Exit(0);
  }
  
  Hostonly_ if ( fp ) fclose(fp);
  
  
  // チェックモードの場合のコメント表示，前処理のみで中止---------------------------------------------------------
  if ( C.CheckParam == ON) {
		Hostonly_ fprintf(stdout, "\n\tCheck mode --- Only pre-process\n\n");
		Hostonly_ fprintf(fp, "\n\tCheck mode --- Only pre-process\n\n");
    return 0;
	}
  
  
  TIMING_stop(tm_init_sct);
  
  
  printf("CurrentTime         = %e\n",CurrentTime);
  printf("CurrentTime_Avr     = %e\n",CurrentTime_Avr);
  printf("Session_StartTime   = %e\n",Session_StartTime);
  printf("Session_CurrentTime = %e\n",Session_CurrentTime);

  printf("Session_LastStep    = %u\n",Session_LastStep);
  printf("Session_CurrentStep = %u\n",Session_CurrentStep);
  printf("Session_StartStep   = %u\n",Session_StartStep);
  printf("CurrentStep         = %u\n",CurrentStep);
  printf("CurrentStep_Avr     = %u\n",CurrentStep_Avr);
  
  Exit(0);
  return 1;
}



// 全Voxelモデルの媒質数とKOSの整合性をチェック
bool FFV::chkMediumConsistency()
{
  int nmSolid = C.NoMediumSolid;
  int nmFluid = C.NoMediumFluid;
  
  if ( numProc > 1 ) {
    int nms = nmSolid;
    int nmf = nmFluid;
    paraMngr->Allreduce(&nms, &nmSolid, 1, MPI_SUM);
    paraMngr->Allreduce(&nmf, &nmFluid, 1, MPI_SUM);
  }
  
  if ( (nmFluid == 0) && (nmSolid == 0) ) {
    Hostonly_ printf("\tError : No medium\n");
    return false;
  }
  
  switch (C.KindOfSolver) {
    case FLOW_ONLY:
    case THERMAL_FLOW:
    case THERMAL_FLOW_NATURAL:
      
      if ( nmFluid == 0 ) {
        Hostonly_ printf("\tError : No FLUID medium\n");
        return false;
      }
      break;
      
    case CONJUGATE_HEAT_TRANSFER:
      if ( ( nmFluid == 0 ) || ( nmSolid == 0 ) ) {
        Hostonly_ printf("\tError : Fluid/Solid should have at least one medium.\n");
        return false;
      }
      break;
      
    case SOLID_CONDUCTION:
      if ( nmSolid == 0 ) {
        Hostonly_ printf("\tError : No Solid medium\n");
        return false;
      }
      break;
  };
  
  return true;
}


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
  else {
    Hostonly_ stamped_printf("\tInvalid keyword is described for Exmple definition\n");
    Exit(0);
  }
}


/**
 * @brief 時刻をRFクラスからv00[4]にコピーする
 * @param [in] time 設定する時刻
 */
void FFV::copyV00fromRF(double m_time) 
{
  RF.setV00(m_time);
  
  double g[4];
  RF.copyV00(g);
  for (int i=0; i<4; i++) v00[i]=(REAL_TYPE)g[i];
}



// コンポーネントの内容リストを表示する
void FFV::display_Compo_Info(FILE* fp)
{
  if ( C.NoBC >0 ) {
    Hostonly_ {
      B.printCompo(stdout, compo_global_bbox, mat, cmp);
      B.printCompo(fp, compo_global_bbox, mat, cmp);
    }
  }
  
  // コンポーネント数がゼロの場合のチェック
  for (int n=1; n<=C.NoBC; n++) {
    if ( cmp[n].getElement() == 0 ) {
      Hostonly_ printf("\tError : No element was found in Component[%d]\n", n);
      fflush(stdout);
      Exit(0);
    }
  }
  
  // Check consistency of boundary condition
  for (int n=1; n<=C.NoBC; n++) {
    if ( cmp[n].getType() == HT_SN ) {
      if ( (C.KindOfSolver == FLOW_ONLY) || (C.KindOfSolver == THERMAL_FLOW) || (C.KindOfSolver == SOLID_CONDUCTION) ) {
        Hostonly_ printf("\tInconsistent parameters of combination between Kind of Solver and Heat Transfer type SN. Check QBCF\n");
        fflush(stdout);
        Exit(0);
      }
    }
    if ( cmp[n].getType() == HT_SF ) {
      if ( (C.KindOfSolver == FLOW_ONLY) || (C.KindOfSolver == THERMAL_FLOW_NATURAL) || (C.KindOfSolver == SOLID_CONDUCTION) ) {
        Hostonly_ printf("\tInconsistent parameters of combination between Kind of Solver and Heat Transfer type SF. Check QBCF\n");
        fflush(stdout);
        Exit(0);
      }
    }
  }
}


// CompoListの内容とセル数の情報を表示する
void FFV::display_CompoList(FILE* fp)
{
  Hostonly_  {
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


// 制御パラメータ，物理パラメータの表示
void FFV::display_Parameters(FILE* fp)
{
  C.displayParams(stdout, fp, IC, &DT, &RF, mat);
  Ex->printPara(stdout, &C);
  Ex->printPara(fp, &C);
  
  // 外部境界面の開口率を表示
  C.printOuterArea(stdout, G_Fcell, G_Acell, G_size);
  C.printOuterArea(fp, G_Fcell, G_Acell, G_size);
  
  // 境界条件のリストと外部境界面のBC設定を表示

  B.printFaceOBC(stdout, G_region);
  B.printFaceOBC(fp, G_region);

  
  /* モニタ情報の表示
   if ( C.Sampling.log == ON ) {
   
   MO.printMonitorInfo(mp, C.HistoryMonitorName, false); // ヘッダのみ
   
   FILE *fp_mon=NULL;
   Hostonly_ {
   if ( !(fp_mon=fopen("sampling_info.txt", "w")) ) {
   stamped_printf("\tSorry, can't open 'sampling_info.txt' file. Write failed.\n");
   return -1;
   }
   }
   
   MO.printMonitorInfo(fp_mon, C.HistoryMonitorName, true);  // 詳細モード
   Hostonly_ if ( fp_mon ) fclose(fp_mon);
   }*/
}

// 計算領域情報を設定する
void FFV::DomainInitialize(const string dom_file)
{
  // グローバルな領域情報のロード
  int ierror = tpCntl.readTPfile(dom_file);
  
  C.importTP(&tpCntl);

  // メンバ変数にパラメータをロード : 分割指示 (1-with / 2-without)
  int div_type = get_DomainInfo();

  
// ##########  
# if 0
  printDomainInfo();
  fflush(stdout);
#endif
// ##########
  
  
  // テキストパーサーのDBを破棄
  if (tpCntl.remove() != TP_NO_ERROR ) 
  {
    Hostonly_ printf("Error : delete textparser\n");
    Exit(0);
  }
  

  // 袖通信の成分数と袖の最大数
  // V(3), P(1), T(1), Cut(6)
  size_t Ncmp = ( C.isCDS() ) ? 6 : 3;
  size_t Nvc  = (size_t)C.guide;
  
  int m_sz[3]  = {G_size[0], G_size[1], G_size[2]};
  int m_div[3] = {G_division[0], G_division[1], G_division[2]};
  
  REAL_TYPE m_org[3] = {G_origin[0], G_origin[1], G_origin[2]};
  REAL_TYPE m_reg[3] = {G_region[0], G_region[1], G_region[2]};
  
  
  // 領域分割モードのパターン
  //      分割指定(G_div指定)    |     domain.txt 
  // 1)  G_divなし >> 自動分割   |  G_orign + G_region + (G_pitch || G_voxel)
  // 2)  G_div指定あり          |  G_orign + G_region + (G_pitch || G_voxel)
  // 3)  G_divなし >> 自動分割   |   + ActiveDomainInfo
  // 4)  G_div指定あり          |   + ActiveDomainInfo
  
  switch (div_type) 
  {
    case 1: // 分割数が指示されている場合
      if ( paraMngr->VoxelInit(m_sz, m_org, m_reg, Nvc, Ncmp) != CPM_SUCCESS )
      {
        cout << "Domain decomposition error : " << endl;
        Exit(0);
      }
      break;
      
    case 2: // 分割数が指示されていない場合
      if ( paraMngr->VoxelInit(m_div, m_sz, m_org, m_reg, Nvc, Ncmp) != CPM_SUCCESS )
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
  
  
  // 有次元の場合に無次元化 (Local)
  if (C.Unit.Param == DIMENSIONAL ) {
    for (int i=0; i<3; i++) {
      origin[i]   /= C.RefLength;
      pitch[i]    /= C.RefLength;
      region[i]   /= C.RefLength;
      G_origin[i] /= C.RefLength;
      G_region[i] /= C.RefLength;
    }
  }
  
  
  // チェック
  unsigned long tz = (unsigned long)size[0] * (unsigned long)size[1] * (unsigned long)size[2];
  if ( tz >= UINT_MAX) {
    Hostonly_ stamped_printf("\n\tError : Product of size[] exceeds UINT_MAX\n\n");
    Exit(0);
  }

}

/**
 @brief 初期インデクスの情報を元に，一層拡大したインデクス値を返す
 @param [in/out] m_st 拡大された開始点（Fortranインデクス）
 @param [in/out] m_ed 拡大された終了点（Fortranインデクス）
 @param [in]     st_i 開始点（Cインデクス）
 @param [in]     len  コンポーネントの存在長さ
 @param [in]     m_x  軸方向のサイズ
 @param [in]     dir  方向
 @param m_id キーID
 */
void FFV::EnlargeIndex(int& m_st, int& m_ed, const int st_i, const int len, const int m_x, const int dir, const int m_id)
{
  int ed_i = st_i + len - 1;
  int n_st = st_i - 1;
  int n_ed = ed_i + 1;
  int max_c1 = m_x + guide;
  
  int label_st, label_ed;
  
  switch (dir) {
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
  if ( ed_i < guide ) { 
    if( nID[label_st] < 0 ){ // 計算領域の外部面に接する場合は，対象外
      m_st = 0;
      m_ed = 0;
    }
    else { // 計算領域内部にある場合（並列時）
      if ( n_ed == guide ) { // ガイドセル1層外側の場合
        m_st = 1; // F index
        m_ed = 1; // F index
      }
      else {
        m_st = 0;
        m_ed = 0;
      }
    }
  }
  
  // BVが+方向のガイドセル内のみにある場合
  else if ( st_i >= max_c1 ) {
    if( nID[label_ed] < 0 ){ // 計算領域の外部面に接する場合は，対象外
      m_st = 0;
      m_ed = 0;
    }
    else {
      if ( n_st == (max_c1 - 1) ) { // ガイドセル1層外側の場合
        m_st = m_x; // F index
        m_ed = m_x; // F index
      }
      else {
        m_st = 0;
        m_ed = 0;
      }
    }
    //debug; Hostonly_ printf("(2) dir=%d : st=%d ed=%d\n", dir, m_st, m_ed);
  }
  
  // BVが内部領域のみにある場合（逐次・並列で同じ処理）
  else if ( (st_i >= guide) && (ed_i < max_c1) ) {
    if ( st_i == guide ) { // 最外層セル
      m_st = 1; // F index
    }
    else {
      m_st = n_st + 1 - guide; // F index
    }
    
    if ( ed_i == (max_c1 - 1) ) { // 最外層セル
      m_ed = m_x; // F index
    }
    else { // 内部
      m_ed = n_ed + 1 - guide; // F index
    }
    //debug; Hostonly_ printf("(3) dir=%d : st=%d ed=%d\n", dir, m_st, m_ed);
  }
  
  // BVが両方向のガイドセルにまたがる場合（逐次・並列で同じ処理）
  else if ( (st_i < guide) && (ed_i >= max_c1) ) {
    m_st = 1; // F index
    m_ed = m_x; // F index
    //debug; Hostonly_ printf("(4) dir=%d : st=%d ed=%d\n", dir, m_st, m_ed);
  }
  
  // BVが-方向のガイドセルから内部領域にある場合
  else if ( (st_i < guide) && (ed_i < max_c1) ) {
    m_st = 1; // F index
    
    if ( ed_i == (max_c1 - 1) ) { // 最外層セル
      m_ed = m_x; // F index
    }
    else { // 内部
      m_ed = n_ed + 1 - guide; // F index
    }
    //debug; Hostonly_ printf("(5) dir=%d : st=%d ed=%d\n", dir, m_st, m_ed);
  }
  
  // BVが+方向のガイドセルから内部領域にある場合
  else if ( (st_i < max_c1) && (ed_i >= max_c1) ) {
    m_ed = m_x; // F index
    
    if ( st_i == guide ) { // 端点
      m_st = 1; // F index
    }
    else { // 内部
      m_st = n_st + 1 - guide; // F index
    }
    //debug; Hostonly_ printf("(6) dir=%d : st=%d ed=%d\n", dir, m_st, m_ed);
  }
  
  else {
    string m_dir;
    if      ( dir == 0 ) m_dir = "X";
    else if ( dir == 1 ) m_dir = "Y";
    else                 m_dir = "Z";
    
    Hostonly_ stamped_printf("\tError : Unexpected case for ID=%d, (%d - %d): %s\n", m_id, st_i, ed_i, m_dir.c_str());
    Exit(0);
  }
  
}


/** 固定パラメータの設定 */
void FFV::fixed_parameters()
{
  // 精度
  if ( sizeof(REAL_TYPE) == sizeof(double) ) {
    C.Mode.Precision = FP_DOUBLE;
  }
  else {
    C.Mode.Precision = FP_SINGLE;
  }
  
  // ログファイル名
  C.HistoryName        = "history_base.txt";
  C.HistoryCompoName   = "history_compo.txt";
  C.HistoryDomfxName   = "history_domainflux.txt";
  C.HistoryForceName   = "history_force.txt";
  C.HistoryWallName    = "history_log_wall.txt";
  C.HistoryItrName     = "history_iteration.txt";
  C.HistoryMonitorName = "sample.log";
  
  C.f_Pressure       = "prs_";
  C.f_Velocity       = "vel_";
  C.f_Temperature    = "tmp_";
  C.f_AvrPressure    = "prsa_";
  C.f_AvrVelocity    = "vela_";
  C.f_AvrTemperature = "tmpa_";
  C.f_DivDebug       = "div_";
  C.f_Helicity       = "hlt_";
  C.f_TotalP         = "tp_";
  C.f_I2VGT          = "i2vgt_";
  C.f_Vorticity      = "vrt_";
  
}


/**
 * @brief 並列処理時の各ノードの分割数を集めてファイルに保存する
 */
void FFV::gather_DomainInfo()
{
  // 統計処理の母数
  double d = 1.0 /(double)numProc;
  double r;
  
  if ( numProc > 1 ) {
    r = 1.0 /(REAL_TYPE)(numProc-1);
  }
  else {
    r = 1.0;
  }
  
  int* m_size=NULL;           ///< 領域分割数
  unsigned long* bf_fcl=NULL; ///< Fluid cell
  unsigned long* bf_wcl=NULL; ///< Solid cell
  unsigned long* bf_acl=NULL; ///< Active cell
  
  REAL_TYPE* m_org =NULL;     ///< 基点
  REAL_TYPE* m_Lbx =NULL;     ///< 領域サイズ
  unsigned long* bf_srf=NULL; ///< 表面数
  
  int* st_buf=NULL; ///< 
  int* ed_buf=NULL; ///< 
  
  
  if( !(m_size = new int[numProc*3]) ) Exit(0);
  if( !(bf_fcl = new unsigned long[numProc]) )   Exit(0);
  if( !(bf_wcl = new unsigned long[numProc]) )   Exit(0);
  if( !(bf_acl = new unsigned long[numProc]) )   Exit(0);
  
  if( !(m_org  = new REAL_TYPE[numProc*3]) ) Exit(0);
  if( !(m_Lbx  = new REAL_TYPE[numProc*3]) ) Exit(0);
  if( !(bf_srf = new unsigned long[numProc]) )   Exit(0);
  
  if( !(st_buf = new int [numProc*3]) ) Exit(0);
  if( !(ed_buf = new int [numProc*3]) ) Exit(0);
  
  // 領域情報の収集
  if ( numProc > 1 ) {
    if ( paraMngr->Gather(size, 3, m_size, 3, 0) != CPM_SUCCESS ) Exit(0);
    if ( paraMngr->Gather(origin, 3, m_org, 3, 0) != CPM_SUCCESS ) Exit(0);
    if ( paraMngr->Gather(region, 3, m_Lbx, 3, 0) != CPM_SUCCESS ) Exit(0);
    if ( paraMngr->Gather(&L_Fcell, 1, bf_fcl, 1, 0) != CPM_SUCCESS ) Exit(0);
    if ( paraMngr->Gather(&L_Wcell, 1, bf_wcl, 1, 0) != CPM_SUCCESS ) Exit(0);
    if ( paraMngr->Gather(&L_Acell, 1, bf_acl, 1, 0) != CPM_SUCCESS ) Exit(0);
  }
  else { // serial
    memcpy(m_size, size, 3*sizeof(int));
    bf_fcl[0] = G_Fcell;
    bf_wcl[0] = G_Wcell;
    bf_acl[0] = G_Acell;
    memcpy(m_org, origin, 3*sizeof(REAL_TYPE));
    memcpy(m_Lbx, region, 3*sizeof(REAL_TYPE));
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
  unsigned long m_srf = (unsigned long)(2*(ix*jx + jx*kx + kx*ix));
  
  if ( nID[X_MINUS] < 0 ) m_srf -= (unsigned long)(jx*kx);  // remove face which does not join communication
  if ( nID[Y_MINUS] < 0 ) m_srf -= (unsigned long)(ix*kx);
  if ( nID[Z_MINUS] < 0 ) m_srf -= (unsigned long)(ix*jx);
  if ( nID[X_PLUS]  < 0 ) m_srf -= (unsigned long)(jx*kx);
  if ( nID[Y_PLUS]  < 0 ) m_srf -= (unsigned long)(ix*kx);
  if ( nID[Z_PLUS]  < 0 ) m_srf -= (unsigned long)(ix*jx);
  
  if ( numProc > 1 ) {
    if ( paraMngr->Gather(&m_srf, 1, bf_srf, 1, 0) != CPM_SUCCESS ) Exit(0);
  }
  else 
  {
    bf_srf[0] = m_srf;
  }
  
  // mean of domain
  unsigned long m_vol = 0;
  unsigned long m_efv = 0;
  m_srf = 0;
  
  for (int i=0; i<numProc; i++) {
    ix = m_size[3*i];
    jx = m_size[3*i+1];
    kx = m_size[3*i+2];
    m_vol += (unsigned long)(ix*jx*kx);
    m_srf += bf_srf[i];
    m_efv += bf_acl[i];
  }
  
  double d_vol = (double)m_vol * d;
  double d_srf = (double)m_srf * d;
  double d_efv = (double)m_efv * d;
  
  // std. deviation of domain
  double vol_dv = 0.0;
  double srf_dv = 0.0;
  double efv_dv = 0.0;
  
  unsigned long d1, d2, d3;
  
  for (int i=0; i<numProc; i++) {
    ix = m_size[3*i];
    jx = m_size[3*i+1];
    kx = m_size[3*i+2];
    d1 = (unsigned long)(ix*jx*kx) - m_vol;
    d2 = bf_srf[i] - m_srf;
    d3 = bf_acl[i] - m_efv;
    vol_dv += (double)d1*(double)d1;
    srf_dv += (double)d2*(double)d2;
    efv_dv += (double)d3*(double)d3;
  }
  vol_dv = sqrt(vol_dv*r);
  srf_dv = sqrt(srf_dv*r);
  efv_dv = sqrt(efv_dv*r);
  
  
  FILE *fp=NULL;
  if ( !(fp=fopen("DomainInfo.txt", "w")) ) {
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
              m_Lbx[i*3]*C.RefLength,  m_Lbx[i*3+1]*C.RefLength,  m_Lbx[i*3+2]*C.RefLength, m_Lbx[i*3],  m_Lbx[i*3+1],  m_Lbx[i*3+2]);
      
      if (C.NoBC != 0) fprintf(fp, "\t no            Label   Mat    i_st    i_ed    j_st    j_ed    k_st    k_ed\n");
    }
    
    if ( numProc > 1 ) {
      
      for (int n=1; n<=C.NoBC; n++) {
        if( paraMngr->Gather(cmp[n].getBbox_st(), 3, st_buf, 3, 0) != CPM_SUCCESS ) Exit(0);
        if( paraMngr->Gather(cmp[n].getBbox_ed(), 3, ed_buf, 3, 0) != CPM_SUCCESS ) Exit(0);
        
        
        Hostonly_ {
          fprintf(fp,"\t%3d %16s %5d %7d %7d %7d %7d %7d %7d\n",
                  n, cmp[n].getLabel().c_str(), cmp[n].getMatOdr(), st_buf[i*3], ed_buf[i*3], st_buf[i*3+1], ed_buf[i*3+1], st_buf[i*3+2], ed_buf[i*3+2]);
        }
      }
    }
  }
  
  Hostonly_ {
    fprintf(fp, "\n");
    fprintf(fp,"\n\t--------------------------------------------------\n");
    fprintf(fp,"\tReport of Whole Domain Statistics\n");
    fprintf(fp,"\tDomain size         = %7d %7d %7d\n", G_size[0], G_size[1], G_size[2]);
    fprintf(fp,"\tNumber of voxels    = %12.6e\n", vol);
    fprintf(fp,"\tNumber of surface   = %12.6e\n", srf);
    fprintf(fp,"\tEffective voxels    = %12.6e (%6.2f%%)\n", (REAL_TYPE)G_Acell, 100.0*(REAL_TYPE)G_Acell/vol);
    fprintf(fp,"\tFluid voxels        = %12.6e (%6.2f%%)\n", (REAL_TYPE)G_Fcell, 100.0*(REAL_TYPE)G_Fcell/vol);
    fprintf(fp,"\tWall  voxels        = %12.6e (%6.2f%%)\n", (REAL_TYPE)G_Wcell, 100.0*(REAL_TYPE)G_Wcell/vol);
    if ( numProc == 1 ) {
      fprintf(fp,"\tDivision :          = %d : %s\n", numProc, "Serial");
    }
    else {
      fprintf(fp,"\tDivision :          = %d : %s\n", numProc, "Equal segregation");
    }
    fprintf(fp,"\n\t--------------------------------------------------\n");
    fprintf(fp,"\tDomain Statistics per MPI process\n");
    fprintf(fp,"\tMean volume in each domain           = %12.6e\n", m_vol);
    fprintf(fp,"\tStd. deviation of domain             = %12.6e\n", vol_dv);
    fprintf(fp,"\tMean comm. in each domain            = %12.6e\n", m_srf);
    fprintf(fp,"\tStd. deviation of surface            = %12.6e\n", srf_dv);
    fprintf(fp,"\tMean effective volume in each domain = %12.6e\n", m_efv);
    fprintf(fp,"\tStd. deviation of effective volume   = %12.6e\n", efv_dv);
    fprintf(fp,"\n");
    
    fprintf(fp,"\tDomain :     ix     jx     kx       Volume Vol_dv[%%]      Surface Srf_dv[%%] Fluid[%%] Solid[%%]      Eff_Vol Eff_Vol_dv[%%]      Eff_Srf Eff_srf_dv[%%]  Itr_scheme\n");
    fprintf(fp,"\t---------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
    
    REAL_TYPE tmp_vol, tmp_acl, tmp_fcl, tmp_wcl;
    for (int i=0; i<numProc; i++) {
      ix = m_size[3*i];
      jx = m_size[3*i+1];
      kx = m_size[3*i+2];
      tmp_vol = (REAL_TYPE)(ix*jx*kx);
      tmp_acl = (REAL_TYPE)bf_acl[i];
      tmp_fcl = (REAL_TYPE)bf_fcl[i];
      tmp_wcl = (REAL_TYPE)bf_wcl[i];
      fprintf(fp,"\t%6d : %6d %6d %6d ", i, ix, jx, kx);
      fprintf(fp,"%12.4e  %8.3f ", tmp_vol, 100.0*(tmp_vol-m_vol)/m_vol);
      fprintf(fp,"%12.4e  %8.3f ", bf_srf[i], (m_srf == 0.0) ? 0.0 : 100.0*(bf_srf[i]-m_srf)/m_srf);
      fprintf(fp,"%8.3f %8.3f ", 100.0*tmp_fcl/tmp_vol, 100.0*tmp_wcl/tmp_vol);
      fprintf(fp,"%12.4e      %8.3f \n", tmp_acl, 100.0*(tmp_acl-m_efv)/m_efv);
    }
    fprintf(fp,"\t---------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
  }
  
  if (fp) fclose(fp);
  
  
  if( m_size ) { delete [] m_size; m_size=NULL; }
  if( m_org  ) { delete [] m_org;  m_org =NULL; }
  if( m_Lbx  ) { delete [] m_Lbx;  m_Lbx =NULL; }
  if( bf_srf ) { delete [] bf_srf; bf_srf=NULL; }
  if( bf_fcl ) { delete [] bf_fcl; bf_fcl=NULL; }
  if( bf_wcl ) { delete [] bf_wcl; bf_wcl=NULL; }
  if( bf_acl ) { delete [] bf_acl; bf_acl=NULL; }
  if( st_buf ) { delete [] st_buf; st_buf=NULL; }
  if( ed_buf ) { delete [] ed_buf; ed_buf=NULL; }
}



// 組み込み例題の設定
void FFV::getExample(Control* Cref, TPControl* tpCntl)
{
  string keyword;
  string label;
  
  label = "/Steer/Example";
  
  if ( !(tpCntl->GetValue(label, &keyword )) ) {
    Hostonly_ stamped_printf("\tExample error\n");
    Exit(0);
  }
  
  if     ( FBUtility::compare(keyword, "Parallel_Plate_2D") ) Cref->Mode.Example = id_PPLT2D;
  else if( FBUtility::compare(keyword, "Duct") )              Cref->Mode.Example = id_Duct;
  else if( FBUtility::compare(keyword, "SHC1D") )             Cref->Mode.Example = id_SHC1D;
  else if( FBUtility::compare(keyword, "Performance_Test") )  Cref->Mode.Example = id_PMT;
  else if( FBUtility::compare(keyword, "Rectangular") )       Cref->Mode.Example = id_Rect;
  else if( FBUtility::compare(keyword, "Cylinder") )          Cref->Mode.Example = id_Cylinder;
  else if( FBUtility::compare(keyword, "Back_Step") )         Cref->Mode.Example = id_Step;
  else if( FBUtility::compare(keyword, "Polygon") )           Cref->Mode.Example = id_Polygon;
  else if( FBUtility::compare(keyword, "Sphere") )            Cref->Mode.Example = id_Sphere;
  else {
    Hostonly_ stamped_printf("\tInvalid keyword is described for Example definition\n");
    Exit(0);
  }
}



// インターバルの初期化
void FFV::init_Interval()
{
  
  // セッションの初期時刻をセット
  for (int i=0; i<Interval_Manager::tg_END; i++) {
    C.Interval[i].setTime_init(Session_StartTime);
  }
  
  // インターバルの初期化
  double m_dt    = DT.get_DT();
  double m_tm    = CurrentTime;  // 設定した？
  unsigned m_stp = CurrentStep;
  
  if ( !C.Interval[Interval_Manager::tg_console].initTrigger(m_stp, m_tm, m_dt, Interval_Manager::tg_console) ) {  // 基本履歴のコンソールへの出力
    Hostonly_ printf("\t Error : Interval for Console output is asigned to zero.\n");
    Exit(0);
  }
  if ( !C.Interval[Interval_Manager::tg_history].initTrigger(m_stp, m_tm, m_dt, Interval_Manager::tg_history) ) {  // 履歴のファイルへの出力
    Hostonly_ printf("\t Error : Interval for History output is asigned to zero.\n");
    Exit(0);
  }
  if ( !C.Interval[Interval_Manager::tg_instant].initTrigger(m_stp, m_tm, m_dt, Interval_Manager::tg_instant) ) {  // 瞬時値ファイル
    Hostonly_ printf("\t Error : Interval for Instantaneous output is asigned to zero.\n");
    Exit(0);
  }
  if ( C.Mode.Average == ON ) {
    //if ( !C.Interval[Interval_Manager::tg_average].initTrigger(m_stp, m_tm, m_dt, Interval_Manager::tg_average) ) {  // 平均値ファイル
    //  Hostonly_ printf("\t Error : Interval for Average output is asigned to zero.\n");
    //  Exit(0);
    //}
    // tg_averageの初期化はLoop中で行う．平均値開始時刻が未知のため．
    if ( !C.Interval[Interval_Manager::tg_avstart].initTrigger(m_stp, m_tm, m_dt, Interval_Manager::tg_avstart) ) {  // 平均値開始
      Hostonly_ printf("\t Error : Interval for Average start is asigned to zero.\n");
      Exit(0);
    }
  }
  if ( C.Sampling.log == ON ) {
    if ( !C.Interval[Interval_Manager::tg_sampled].initTrigger(m_stp, m_tm, m_dt, Interval_Manager::tg_sampled) ) {  // サンプリング履歴
      Hostonly_ printf("\t Error : Interval for Sampling output is asigned to zero.\n");
      Exit(0);
    }    
  }
  
}


// 履歴の出力準備
void FFV::prep_HistoryOutput()
{
  // マスターノードでの履歴出力準備
  H = new History(&C);
  
  Hostonly_ {
    H->printHistoryTitle(stdout, IC, &C);
    
    // コンポーネント情報
    if ( C.Mode.Log_Base == ON ) {
      // 基本情報　history.log, history_compo.log, history_domfx.log
      if ( !(fp_b=fopen(C.HistoryName.c_str(), "w")) ) {
        stamped_printf("\tSorry, can't open '%s' file. Write failed.\n", C.HistoryName.c_str());
        Exit(0);
      }
      H->printHistoryTitle(fp_b, IC, &C);
      
      // コンポーネント履歴情報
      if ( !(fp_c=fopen(C.HistoryCompoName.c_str(), "w")) ) {
        stamped_printf("\tSorry, can't open '%s' file. Write failed.\n", C.HistoryCompoName.c_str());
        Exit(0);
      }
      H->printHistoryCompoTitle(fp_c, cmp, &C);
      
      // 流量収支情報　
      if ( !(fp_d=fopen(C.HistoryDomfxName.c_str(), "w")) ) {
        stamped_printf("\tSorry, can't open '%s' file. Write failed.\n", C.HistoryDomfxName.c_str());
        Exit(0);
      }
      H->printHistoryDomfxTitle(fp_d, &C);
      
      // 力の履歴情報　
      if ( !(fp_f=fopen(C.HistoryForceName.c_str(), "w")) ) {
        stamped_printf("\tSorry, can't open '%s' file. Write failed.\n", C.HistoryForceName.c_str());
        Exit(0);
      }
      H->printHistoryForceTitle(fp_f);
    }
    
    // 反復履歴情報　history_itr.log
    if ( C.Mode.Log_Itr == ON ) {
      if ( !(fp_i=fopen(C.HistoryItrName.c_str(), "w")) ) {
				stamped_printf("\tSorry, can't open '%s' file.\n", C.HistoryItrName.c_str());
        Exit(0);
      }
    }
    
    // 壁面情報　history_wall.log
    if ( C.Mode.Log_Wall == ON ) {
      if ( !(fp_w=fopen(C.HistoryWallName.c_str(), "w")) ) {
				stamped_printf("\tSorry, can't open '%s' file.\n", C.HistoryWallName.c_str());
        Exit(0);
      }
      H->printHistoryWallTitle(fp_w);
    }
  }
  
  // サンプリング指定がある場合，モニタ結果出力ファイル群のオープン
  //if ( C.Sampling.log == ON ) MO.openFile(C.HistoryMonitorName);
}



// 読み込んだ領域情報のデバッグライト
void FFV::printDomainInfo()
{
  cout << "\n####### read parameters ########" << endl;
  cout << " G_org      = " << G_origin[0] << "," << G_origin[1] << "," << G_origin[2] << endl;
  cout << " G_voxel    = " << G_size[0]   << "," << G_size[1]   << "," << G_size[2]   << endl;
  cout << " G_pitch    = " << pitch[0]    << "," << pitch[1]    << "," << pitch[2]    << endl;
  cout << " G_region   = " << G_region[0] << "," << G_region[1] << "," << G_region[2] << endl;
  cout << " G_div      = " << G_division[0]  << "," << G_division[1]  << "," << G_division[2]  << endl;
}



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
        
        if ( GET_FACE_BC(s, BC_FACE_W) == n ) {
          if( i < nst[0] ) { nst[0] = i; }
          if( i > ned[0] ) { ned[0] = i; }
          if( j < nst[1] ) { nst[1] = j; }
          if( j > ned[1] ) { ned[1] = j; }
          if( k < nst[2] ) { nst[2] = k; }
          if( k > ned[2] ) { ned[2] = k; }
        }
        
        if ( GET_FACE_BC(s, BC_FACE_E) == n ) {
          if( i < nst[0] ) { nst[0] = i; }
          if( i > ned[0] ) { ned[0] = i; }
          if( j < nst[1] ) { nst[1] = j; }
          if( j > ned[1] ) { ned[1] = j; }
          if( k < nst[2] ) { nst[2] = k; }
          if( k > ned[2] ) { ned[2] = k; }
        }
        
        if ( GET_FACE_BC(s, BC_FACE_S) == n ) {
          if( i < nst[0] ) { nst[0] = i; }
          if( i > ned[0] ) { ned[0] = i; }
          if( j < nst[1] ) { nst[1] = j; }
          if( j > ned[1] ) { ned[1] = j; }
          if( k < nst[2] ) { nst[2] = k; }
          if( k > ned[2] ) { ned[2] = k; }
        }
        
        if ( GET_FACE_BC(s, BC_FACE_N) == n ) {
          if( i < nst[0] ) { nst[0] = i; }
          if( i > ned[0] ) { ned[0] = i; }
          if( j < nst[1] ) { nst[1] = j; }
          if( j > ned[1] ) { ned[1] = j; }
          if( k < nst[2] ) { nst[2] = k; }
          if( k > ned[2] ) { ned[2] = k; }
        }
        
        if ( GET_FACE_BC(s, BC_FACE_B) == n ) {
          if( i < nst[0] ) { nst[0] = i; }
          if( i > ned[0] ) { ned[0] = i; }
          if( j < nst[1] ) { nst[1] = j; }
          if( j > ned[1] ) { ned[1] = j; }
          if( k < nst[2] ) { nst[2] = k; }
          if( k > ned[2] ) { ned[2] = k; }
        }
        
        if ( GET_FACE_BC(s, BC_FACE_T) == n ) {
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
        
        if ( ( s & MASK_6) == n ) {
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




// 外部境界条件を読み込み，Controlクラスに保持する
void FFV::setBCinfo()
{
  // パラメータファイルの情報を元にCompoListの情報を設定する
  B.loadBC_Local(&C);
  
  // 各コンポーネントが存在するかどうかを保持しておく
  setEnsComponent();
  
  // KOSと境界条件種類の整合性をチェック
  B.chkBCconsistency(C.KindOfSolver);
  
  // ガイドセル上にパラメータファイルで指定する媒質インデクスを代入する．周期境界の場合の処理も含む．
  for (int face=0; face<NOFACE; face++) {
    V.adjMedium_on_GC(face, d_mid, BC.export_OBC(face)->get_Class(), 
                      BC.export_OBC(face)->get_GuideMedium(), BC.export_OBC(face)->get_PrdcMode());
  }
}



// HEX,FANコンポーネントなどの体積率とbboxなどをセット
// インデクスの登録と配列確保はVoxEncode()で、コンポーネント領域のリサイズ後に行う
void FFV::setComponentVF()
{
  const int subsampling = 20; // 体積率のサブサンプリングの基数
  int f_st[3], f_ed[3];
  double flop;
  
  CompoFraction CF(size, guide, pitch, origin, subsampling);
  
  for (int n=1; n<=C.NoBC; n++) {
    
    if ( cmp[n].isFORCING() ) {
      // 形状パラメータのセット
      switch ( cmp[n].getType() ) {
        case HEX:
          CF.setShapeParam(cmp[n].nv, cmp[n].oc, cmp[n].dr, cmp[n].depth, cmp[n].shp_p1, cmp[n].shp_p2);
          break;
          
        case FAN:
          CF.setShapeParam(cmp[n].nv, cmp[n].oc, cmp[n].depth, cmp[n].shp_p1, cmp[n].shp_p2);
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
  for (int i=0; i<3; i++) {
    org[i] = C.org[i] - C.dx[i]*(REAL_TYPE)C.GuideOut;
    pit[i] = C.dx[i];
  }
  
  // 出力ファイルの指定が有次元の場合
  if ( C.Unit.File == DIMENSIONAL ) {
    for (int i=0; i<3; i++) {
      org[i] *= C.RefLength;
      pit[i] *= C.RefLength;
    }
  }
  F.writeRawSPH(cvf, size, guide, org, pit, FP_SINGLE);
#endif
// ##########
  
}



// コンポーネントが存在するかを保持しておく
void FFV::setEnsComponent()
{
  int c;
  
  // Forcing > HEX, FAN, DARCY
  c = 0;
  for (int n=1; n<=C.NoBC; n++) {
    if ( cmp[n].isFORCING() ) c++;
  }
  if ( c>0 ) C.EnsCompo.forcing = ON;
  
  // Heat source > HEAT_SRC, CNST_TEMP
  c = 0;
  for (int n=1; n<=C.NoBC; n++) {
    if ( cmp[n].isHsrc() ) c++;
  }
  if ( c>0 ) C.EnsCompo.hsrc = ON;
  
  // 周期境界 > PERIODIC
  c = 0;
  for (int n=1; n<=C.NoBC; n++) {
    if ( cmp[n].getType() == PERIODIC ) c++;
  }
  if ( c>0 ) C.EnsCompo.periodic = ON;
  
  // 流出境界 > OUTFLOW
  c = 0;
  for (int n=1; n<=C.NoBC; n++) {
    if ( cmp[n].getType() == OUTFLOW ) c++;
  }
  if ( c>0 ) C.EnsCompo.outflow = ON;
  
  // 体積率コンポーネント
  c = 0;
  for (int n=1; n<=C.NoBC; n++) {
    if ( cmp[n].isVFraction() ) c++;
  }
  if ( c>0 ) C.EnsCompo.fraction = ON;
  
  // モニタ
  c = 0;
  for (int n=1; n<=C.NoBC; n++) {
    if ( cmp[n].isMONITOR() ) c++;
  }
  // MONITOR_LISTでCELL_MONITORが指定されている場合，C.EnsCompo.monitor==ON
  if ( (C.isMonitor() == ON) && (c < 1) ) {
    Hostonly_ stamped_printf("\tError : Cell_Monitor in MONITOR_LIST is specified, however any MONITOR can not be found.\n");
    Exit(0);
  }
  if ( (C.isMonitor() == OFF) && (c > 0) ) {
    Hostonly_ stamped_printf("\tError : Cell_Monitor in MONITOR_LIST is NOT specified, however MONITOR section is found in LocalBoundary.\n");
    Exit(0);
  }
  
}


/**
 * @brief コンポーネントのローカルなBbox情報からグローバルなBbox情報を求める
 */
void FFV::setGlobalCmpIdx()
{
  int st_i, st_j, st_k, ed_i, ed_j, ed_k;
  int node_st_i, node_st_j, node_st_k;
  int st[3], ed[3];
  
  // グローバルインデクスの配列インスタンス
  compo_global_bbox = new int[6*(C.NoCompo+1)];
  int* cgb = compo_global_bbox;
  
  // ローカルインデクスからグローバルインデクスに変換
  for (int m=1; m<=C.NoCompo; m++) {
    
    if ( !cmp[m].isEns() ) { // コンポーネントが存在しないノードはゼロを代入
      cgb[6*m+0] = 0;
      cgb[6*m+1] = 0;
      cgb[6*m+2] = 0;
      cgb[6*m+3] = 0;
      cgb[6*m+4] = 0;
      cgb[6*m+5] = 0;
    }
    else { // コンポーネントが存在する場合
      cmp[m].getBbox(st, ed);
      st_i = st[0];
      st_j = st[1];
      st_k = st[2];
      ed_i = ed[0];
      ed_j = ed[1];
      ed_k = ed[2];
      
      if ( numProc > 1 ) {
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
      else {
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
  
  for (int n=1; n<=C.NoBC; n++) {
    if ( numProc > 1 ) {
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
        
        for (int m=0; m<numProc; m++) {
          if ( m_eArray[m]==1 ) { // コンポーネントの存在ランクのみを対象とする
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
    DomainMonitor( BC.export_OBC(), &C, flop_task);
    
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
    
    fb_set_real_s_(d_p, size, &guide, &ip);
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
      
      fb_set_real_s_(d_t, size, &guide, &it);
      
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
    
    //流出境界の流出速度の算出
    REAL_TYPE coef = deltaX/deltaT;
    REAL_TYPE m_av[2];
    BC.mod_div(d_ws, d_bcv, coef, tm, v00, m_av, flop_task);
    DomainMonitor(BC.export_OBC(), &C, flop_task);
    
    //if ( C.isHeatProblem() ) BC.InnerTBC_Periodic()
    
  }

  
  // 初期解およびリスタート解の同期
  if ( paraMngr->BndCommV3DEx(d_v, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
  if ( paraMngr->BndCommS3D  (d_p, size[0], size[1], size[2], guide, 1    ) != CPM_SUCCESS ) Exit(0);
  
  if ( C.isHeatProblem() ) {
    if ( paraMngr->BndCommS3D  (d_p, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
  }
  
  
  // VOF
  if ( C.BasicEqs == INCMP_2PHASE ) {
    setVOF();
    if ( paraMngr->BndCommS3D  (d_vof, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
  }
  
}



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
  
  // Interval Managerの計算の主管理タグ[tg_compute]に値を初期値を設定
  if ( !C.Interval[Interval_Manager::tg_compute].initTrigger(0, 0.0, DT.get_DT(), Interval_Manager::tg_compute, 
                                                             (double)(C.RefLength/C.RefVelocity)) ) 
  {
    Hostonly_ printf("\t Error : Computation Period is asigned to zero.\n");
    Exit(0);
  }
  C.LastStep = C.Interval[Interval_Manager::tg_compute].getIntervalStep();
  
  
  // C.Interval[Interval_Manager::tg_compute].initTrigger()で初期化後
  C.setParameters(mat, cmp, &RF, BC.export_OBC());
  
  
  // 媒質による代表パラメータのコピー
  B.setRefValue(mat, cmp, &C);
  
  // パラメータの無次元化（正規化）に必要な参照物理量の設定
  B.setRefMedium(mat, C.RefMat);
}


// 各種例題のモデルをセット
void FFV::setModel(double& PrepMemory, double& TotalMemory, FILE* fp)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  
  switch (C.Mode.Example) {
      
    case id_Polygon:
      
      C.get_Polygon();
      
      // PolylibとCutlibのセットアップ
      //setup_Polygon2CutInfo(PrepMemory, TotalMemory, fp);
      
      
      if ( !C.isCDS() ) {
        unsigned long zc = V.Solid_from_Cut(d_mid, d_cut, id_of_solid);
        Hostonly_ printf("\tGenerated Solid cell from cut = %ld\n", zc);
      }
      break;
      
    case id_Sphere:
      if ( !C.isCDS() ) {
        Ex->setup(d_mid, &C, G_origin, C.NoMedium, mat);
      }
      else {
        // cutをアロケートし，初期値1.0をセット
        setup_CutInfo4IP(PrepMemory, TotalMemory, fp);
        Ex->setup_cut(d_mid, &C, G_origin, C.NoMedium, mat, d_cut);
      }
      break;
      
    default: // ほかのIntrinsic problems
      if ( C.isCDS() ) {
        setup_CutInfo4IP(PrepMemory, TotalMemory, fp);
      }
      Ex->setup(d_mid, &C, G_origin, C.NoMedium, mat);
      break;
  }
  
  // midのガイドセル同期
  if ( paraMngr->BndCommS3D(d_mid, ix, jx, kx, guide, 1) != CPM_SUCCESS ) Exit(0);
}



// IP用にカット領域をアロケートする
void FFV::setup_CutInfo4IP(double& m_prep, double& m_total, FILE* fp)
{
  Hostonly_ {
    fprintf(fp, "\n---------------------------------------------------------------------------\n\n");
    fprintf(fp, "\t>> Cut Info\n\n");
    fprintf(stdout, "\n---------------------------------------------------------------------------\n\n");
    fprintf(stdout, "\t>> Cut Info\n\n");
  }
  
  size_t n_cell[3];
  
  for (int i=0; i<3; i++) {
    n_cell[i] = (size_t)(size[i] + 2*guide); // 分割数+ガイドセル
  }
  size_t size_n_cell = n_cell[0] * n_cell[1] * n_cell[2];
  
  
  TIMING_start(tm_init_alloc);
  allocArray_Cut(m_total);
  TIMING_stop(tm_init_alloc);
  
  // 使用メモリ量　
  double cut_mem, G_cut_mem;
  G_cut_mem = cut_mem = (double)size_n_cell * (double)(6*sizeof(float) + sizeof(int));
  m_prep += cut_mem;
  m_total+= cut_mem;
  
  if ( numProc > 1 ) {
    if ( paraMngr->Allreduce(&cut_mem, &G_cut_mem, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
  }
  
  Hostonly_  {
    FBUtility::MemoryRequirement("Cut", G_cut_mem, cut_mem, stdout);
    FBUtility::MemoryRequirement("Cut", G_cut_mem, cut_mem, fp);
  }
  
  // 初期値のセット
  for (size_t i=0; i<size_n_cell*6; i++) {
    d_cut[i] = 1.0f;
  }
  
}



// VOF値を気体(0.0)と液体(1.0)で初期化
void FFV::setVOF()
{
  size_t m;
  int s, odr;
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  for (int k=1; k<=size[2]; k++) {
    for (int j=1; j<=size[1]; j++) {
      for (int i=1; i<=size[0]; i++) {
        m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        s = d_bcd[m];
        odr = DECODE_CMP(s);
        if ( cmp[odr].getState() == FLUID ) {
          d_vof[m] = ( cmp[odr].getPhase() == GAS ) ? 0.0 : 1.0;
        }
      }
    }
  }
}


// ポリゴンのカット情報からVBCのboxをセット
void FFV::VIBC_Bbox_from_Cut()
{
  int f_st[3], f_ed[3];
  
  for (int n=1; n<=C.NoBC; n++) {
    
    if ( cmp[n].isVBC_IO() ) { // SPEC_VEL || SPEC_VEL_WH || OUTFLOW
      
      // インデクスの計算 > インデクスの登録はVoxEncode()で、コンポーネント領域のリサイズ後に行う
      V.findVIBCbbox(n, d_bcv, f_st, f_ed);
      
      // インデクスのサイズ登録と存在フラグ
      cmp[n].setBbox(f_st, f_ed);
      cmp[n].setEns(ON);
    }
  }
  
}


/**
 * @brief BCIndexにビット情報をエンコードする
 */
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
  if ( C.isHeatProblem() ) {
    V.copyBCIbase(d_bh1, d_bcd);
    V.copyBCIbase(d_bh2, d_bcd);
  }

  // BCIndexP に圧力計算のビット情報をエンコードする -----
  if ( C.isCDS() ) {
    C.NoWallSurface = V.setBCIndexP(d_bcd, d_bcp, d_mid, &BC, cmp, true, d_cut);
  }
  else { // binary
    C.NoWallSurface = V.setBCIndexP(d_bcd, d_bcp, d_mid, &BC, cmp);
  }

#if 0
  V.dbg_chkBCIndexP(d_bcd, d_bcp, "BCindexP.txt", cmp);
#endif
  
  // BCIndexV に速度計算のビット情報をエンコードする -----
  if ( C.isCDS() ) {
    V.setBCIndexV(d_bcv, d_mid, d_bcp, &BC, cmp, true, d_cut, d_bid);
  }
  else { // binary
    V.setBCIndexV(d_bcv, d_mid, d_bcp, &BC, cmp);
  }


// ##########
#if 0
  V.dbg_chkBCIndexV(d_bcv, "BCindexV.txt");
#endif
// ##########
  
  // BCIndexT に温度計算のビット情報をエンコードする -----
  if ( C.isHeatProblem() ) {
    V.setBCIndexH(d_bcd, d_bh1, d_bh2, d_mid, &BC, C.KindOfSolver, cmp);
    
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
  Hostonly_ {
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
  
  // midにロードされたIDをスキャンし，IDの個数を返し，作業用のcolorList配列にIDを保持，midに含まれるIDの数をチェック
  int sc=0;
  if ( (sc = V.scanCell(d_mid, cell_id, C.Hide.Change_ID)) != C.NoMedium ) 
  {
    Hostonly_ stamped_printf("A number of IDs included in voxel model(%d) is not agree with the one in 'Model_Setting'(%d)\n", 
                             sc, C.NoMedium);
    Exit(0);
  }
}

