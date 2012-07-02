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
  
  double flop_count;     ///< flops計算用
  
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
    ShapeMonitor SM(size, guide, C.dx, C.org);
    
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
    
    flop_count = 0.0;
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
  
  // BCIndexにビット情報をエンコードとコンポーネントインデクスの再構築
  VoxEncode();
  

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  // 初期化終了時に、入力パラメータのDBを破棄
  if (tpCntl.remove() != TP_NO_ERROR ) 
  {
    Hostonly_ printf("Error : delete textparser\n");
    Exit(0);
  }
  
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
  setNeighborInfo(C.guide);
  
  
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
  
  CompoFraction CF(size, guide, C.dx, C.org, subsampling);
  
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



// midの情報から各BCコンポーネントのローカルなインデクスを取得する
// 計算内部領域の境界と外部境界とでは，ガイドセル部分にあるコンポーネントIDの取り扱いが異なる
// 外部境界に接する面では，幅はそのまま，始点はガイドセル部分を含む
// 内部境界に接する面では，始点と幅はローカルノード内の計算内部領域に含まれるように調整
void FFV::setLocalCmpIdx_Binary()
{
  int st_i[3], len[3];
  int id;
  int m_st[3], m_ed[3];
  
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
        // GetBndIndexExtGc()は自ノード内でのidのバウンディングボックスを取得．インデクスはローカルインデクスで，ガイドセルを含む配列の基点をゼロとするCのインデクス
        //if ( !dc_mid->GetBndIndexExtGc(id, st_i[0], st_i[1], st_i[2], len[0], len[1], len[2], 0) ) 
        //{
        //  Hostonly_ stamped_printf("\tError : can not get component local index for ID[%d]\n", id);
        //  Exit(0);
        //}
        
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
