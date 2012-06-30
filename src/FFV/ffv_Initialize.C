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
  double TotalMemory;    ///< 計算に必要なメモリ量（ローカル）
  double PrepMemory;     ///< 初期化に必要なメモリ量（ローカル）
  double G_TotalMemory;  ///< 計算に必要なメモリ量（グローバル）
  double G_PrepMemory;   ///< 初期化に必要なメモリ量（グローバル）
  double tmp_memory;     ///< 計算に必要なメモリ量（グローバル）？
  
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
  C.setDomainInfo(paraMngr, procGrp);   C.setNeighborInfo();
  B.setDomainInfo(paraMngr, procGrp);   B.setNeighborInfo();
  V.setDomainInfo(paraMngr, procGrp);   V.setNeighborInfo();
  F.setDomainInfo(paraMngr, procGrp);   F.setNeighborInfo();
  BC.setDomainInfo(paraMngr, procGrp);  BC.setNeighborInfo();
  Ex->setDomainInfo(paraMngr, procGrp); Ex->setNeighborInfo();
  
  // ###########################

  
  // 各例題のパラメータ設定 -----------------------------------------------------
  Ex->setDomain(&C, size, org, reg, pch);
  
  
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
  
  // CompoList, MediumListのポインタをセット
  BC.importCMP_MAT(cmp, mat);
  V.importCMP_MAT(cmp, mat);
  B.importCompoPtr(cmp);
  
  
  // ソルバークラスのノードローカルな変数の設定 -----------------------------------------------------
  dh0     = &C.dh;
  dh      = &C.dh;
  
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
    C.printGlobalDomain(stdout, G_size, G_org, G_reg);
    
    fprintf(fp,"\n---------------------------------------------------------------------------\n");
    fprintf(fp,"\n\t>> Global Domain Information\n\n");
    C.printGlobalDomain(fp, G_size, G_org, G_reg);
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
  VoxScan(&V, fp);
  
  
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

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  // 初期化終了時に、入力パラメータのDBを破棄
  if (tpCntl.remove() != TP_NO_ERROR ) 
  {
    Hostonly_ printf("Error : delete textparser\n");
    Exit(0);
  }
  
  Exit(0);
  return 1;
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
  int m_div[3] = {G_div [0], G_div [1], G_div [2]};
  
  REAL_TYPE m_org[3] = {G_org[0], G_org[1], G_org[2]};
  REAL_TYPE m_reg[3] = {G_reg[0], G_reg[1], G_reg[2]};
  
  
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
  setNeighborInfo();
  
  
  // 有次元の場合に無次元化 (Local)
  if (C.Unit.Param == DIMENSIONAL ) {
    for (int i=0; i<3; i++) {
      org[i]   /= C.RefLength;
      pch[i]   /= C.RefLength;
      reg[i]   /= C.RefLength;
      G_org[i] /= C.RefLength;
      G_reg[i] /= C.RefLength;
    }
  }
  
  
  // チェック
  unsigned long tz = (unsigned long)size[0] * (unsigned long)size[1] * (unsigned long)size[2];
  if ( tz >= UINT_MAX) {
    Hostonly_ stamped_printf("\n\tError : Product of size[] exceeds UINT_MAX\n\n");
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


// 読み込んだ領域情報のデバッグライト
void FFV::printDomainInfo()
{
  cout << "\n####### read parameters ########" << endl;
  cout << " G_org      = " << G_org[0]  << "," << G_org[1]  << "," << G_org[2]  << endl;
  cout << " G_voxel    = " << G_size[0] << "," << G_size[1] << "," << G_size[2] << endl;
  cout << " G_pitch    = " << pch[0]    << "," << pch[1]    << "," << pch[2]    << endl;
  cout << " G_region   = " << G_reg[0]  << "," << G_reg[1]  << "," << G_reg[2]  << endl;
  cout << " G_div      = " << G_div[0]  << "," << G_div[1]  << "," << G_div[2]  << endl;
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


// 各種例題のモデルをセット
void FFV::setModel(double& PrepMemory, double& TotalMemory, FILE* fp)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int id_of_solid = 2;
  
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
        Ex->setup(d_mid, &C, G_org, C.NoMedium, mat);
      }
      else {
        // cutをアロケートし，初期値1.0をセット
        setup_CutInfo4IP(PrepMemory, TotalMemory, fp);
        Ex->setup_cut(d_mid, &C, G_org, C.NoMedium, mat, d_cut);
      }
      break;
      
    default: // ほかのIntrinsic problems
      if ( C.isCDS() ) {
        setup_CutInfo4IP(PrepMemory, TotalMemory, fp);
      }
      Ex->setup(d_mid, &C, G_org, C.NoMedium, mat);
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

/**
 @brief ボクセルをスキャンし情報を表示する
 @param [in] V  前処理クラス
 @param [in] fp ファイルポインタ 
 */
void FFV::VoxScan(VoxInfo* V, FILE* fp)
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
  if ( (sc = V->scanCell(d_mid, cell_id, C.Hide.Change_ID)) != C.NoMedium ) 
  {
    Hostonly_ stamped_printf("A number of IDs included in voxel model(%d) is not agree with the one in 'Model_Setting'(%d)\n", 
                             sc, C.NoMedium);
    Exit(0);
  }
}
