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
 * @file ffv_Initialize.C
 * @brief FFV Class
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
  
  TPControl tpCntl;      ///< テキストパーサのラッパークラス
  
  
  // CPMバージョン表示
  if ( paraMngr->GetMyRankID() == 0 )
  {
    cpm_Base::VersionInfo();
  }
  
  
  // 固定パラメータ
  fixed_parameters();
  
  
  // 前処理段階のみに使用するオブジェクトをインスタンス
  //VoxInfo Vinfo;
  ParseBC B;
  
  
  // CPMのポインタをセット
  //BC.importCPM(paraMngr);
  //Vinfo.importCPM(paraMngr);
  B.importCPM(paraMngr);
  F.importCPM(paraMngr);
  

  
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
    char buf[LABEL];
    memset(buf, 0, sizeof(char)*LABEL);
    strcpy(buf, "Welcome to FFV  ");
    FBUtility::printVersion(fp, buf, VERS_FFV);
    FBUtility::printVersion(mp, buf, VERS_FFV);
    
		memset(buf, 0, sizeof(char)*LABEL);
    strcpy(buf, "FlowBase        ");
    FBUtility::printVersion(fp, buf, FB_VERS);
    FBUtility::printVersion(mp, buf, FB_VERS);
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
  
  if ( C.version != VERS_FFV ) {
    Hostonly_ {
      fprintf(mp, "\t##### Version of Input file (%d) is NOT compliant with FFV ver. %d #####\n", C.version, VERS_FFV);
      fprintf(fp, "\t##### Version of Input file (%d) is NOT compliant with FFV ver. %d #####\n", C.version, VERS_FFV);
    }
    return -1;
  }
  if ( C.FB_version != FB_VERS ) {
    Hostonly_ {
      fprintf(mp, "\t##### Version of Input file (%d) is NOT compliant with FB ver. %d #####\n", C.FB_version, FB_VERS);
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
  if ( paraMngr->GetNumRank() == 1 ) {
    C.FIO.IO_Input  = IO_GATHER;
    C.FIO.IO_Output = IO_GATHER;
    Hostonly_ printf("\tMode of Parallel_Input/_Output was changed because of serial execution.\n");
  }
  
  // 計算領域全体のサイズ，並列計算時のローカルのサイズ，コンポーネントのサイズなどを設定する -----------------------------------------------------
  // 領域情報を記述したファイル名の取得 >> テキストパーサのDB切り替え
  string dom_file = argv[2];
  DomainInitialize(dom_file);

  // DomainInitialize()で使用したテキストパーサーのDBを破棄
  if (tpCntl.remove() != TP_NO_ERROR ) 
  {
    Hostonly_ printf("Error : delete textparser\n");
    Exit(0);
  }
  
  // 組み込み例題クラスへ変数の値をコピー，固有のパラメータ取得，領域設定 -----------------------------------------------------
  Ex->setControlVars (&C);
  Ex->importCPM(paraMngr);
  
  
  // 　再度、入力ファイルをオープン
  ierror = tpCntl.readTPfile(input_file);

  
  // TPControlクラスのポインタを各クラスに渡す
  C.importTP(&tpCntl);
  B.importTP(&tpCntl);
  M.importTP(&tpCntl);
  
  // パラメータを取得
  C.get_Steer_2(IC, &RF);
  
  
  // 組み込み例題の固有パラメータ
  if ( !Ex->getTP(&C, &tpCntl) ) Exit(0);
  
  
  // 媒質情報をパラメータファイルから読み込み，媒質リストを作成する
  Hostonly_  {
    fprintf(mp,"\n---------------------------------------------------------------------------\n\n");
    fprintf(mp,"\n\t>> Medium List\n\n");
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
  
  //Vinfo.setNoCompo_BC(C.NoBC, C.NoCompo);
  
  //B.setControlVars(&C, BC.export_OBC(), mat);
  
  B.setMediumPoint( M.export_MTI() );
  
  B.countMedium(&C);
  
  // CompoListクラスをインスタンス．[0]はダミーとして利用しないので，配列の大きさはプラス１する
  cmp = new CompoList[C.NoCompo+1];
  
  // CompoList, MediumListのポインタをセット
  //BC.setWorkList(cmp, mat);
  //Vinfo.setWorkList(cmp, mat);
  B.receiveCompoPtr(cmp);
  
  
  // ソルバークラスのノードローカルな変数の設定 -----------------------------------------------------
  dh0     = &C.dh;
  dh      = &C.dh;
  gc      = C.guide;
  sz[0]   = C.imax;
  sz[1]   = C.jmax;
  sz[2]   = C.kmax;
  
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
    fprintf(mp,"\n---------------------------------------------------------------------------\n\n");
    fprintf(mp,"\t>> Voxel file information\n\n");
  }
  
  TIMING_start(tm_voxel_prep_sct);
  
  
  
  
  
  
  
  
  // 初期化終了時に、入力パラメータのDBを破棄
  if (tpCntl.remove() != TP_NO_ERROR ) 
  {
    Hostonly_ printf("Error : delete textparser\n");
    Exit(0);
  }
  
  return 1;
}



// 組み込み例題のインスタンス
void FFV::connectExample(Control* Cref)
{
  if      ( Cref->Mode.Example == id_Users )   Ex = dynamic_cast<Intrinsic*>(new IP_Users);
  else if ( Cref->Mode.Example == id_PPLT2D)   Ex = dynamic_cast<Intrinsic*>(new IP_PPLT2D);
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
  const REAL_TYPE *m_pch; ///< 格子幅
  const int *m_sz;        ///< サブドメインの分割数
  const REAL_TYPE *m_org; ///< サブドメインの基点
  const REAL_TYPE *m_reg; ///< サブドメインのサイズ
  const int *g_sz;        ///< 全ドメインの分割数
  const REAL_TYPE *g_org; ///< 全ドメインの基点
  const REAL_TYPE *g_reg; ///< 全ドメインのサイズ
  
  // グローバルな領域情報のロード
  int ret = 0;
  cpm_GlobalDomainInfo *dInfo = cpm_TextParserDomain::Read(dom_file, ret);
  
  
# if 0
  printDomainInfo(dInfo);
  paraMngr->flush(cout);
#endif
  

  // 領域分割
  if( (ret = paraMngr->VoxelInit( dInfo, C.guide, 6 )) != CPM_SUCCESS )
  {
    cout << "Domain decomposition error : " << ret << endl;
    Exit(0);
  }

  // 分割後のパラメータ取得
  m_pch = paraMngr->GetPitch();
  m_sz  = paraMngr->GetLocalVoxelSize();
  m_org = paraMngr->GetLocalOrigin();
  m_reg = paraMngr->GetLocalRegion();
  
  g_sz  = paraMngr->GetGlobalVoxelSize();
  g_org = paraMngr->GetGlobalOrigin();
  g_reg = paraMngr->GetGlobalRegion();
  
  
  // 有次元の場合に無次元化
  REAL_TYPE org[3];
  REAL_TYPE pch[3];
  REAL_TYPE reg[3];
  for (int i=0; i<3; i++) 
  {
    org[i] = m_org[i];
    pch[i] = m_pch[i];
    reg[i] = m_reg[i];
  }
  
  if (C.Unit.Param == DIMENSIONAL ) {
    for (int i=0; i<3; i++) {
      org[i] /= C.RefLength;
      pch[i] /= C.RefLength;
      reg[i] /= C.RefLength;
    }
  }
  
  // 各例題のパラメータ設定
  Ex->setDomain(&C, m_sz, org, reg, pch);
  

  // グローバルな値の保持 （無次元値）
  for (int i=0; i<3; i++) {
    G_size[i] = g_sz[i];
    G_org[i]  = g_org[i];
    G_reg[i]  = g_reg[i];
  }
  
  // チェック
  unsigned long tz = (unsigned long)m_sz[0] * (unsigned long)m_sz[1] * (unsigned long)m_sz[2];
  if ( tz >= UINT_MAX) {
    Hostonly_ stamped_printf("\n\tError : Product of size[] exceeds UINT_MAX\n\n");
    Exit(0);
  }
  
  // コントロールクラスのメンバ変数で値を保持
  C.setDomainInfo(m_sz, m_org, m_pch, m_reg);
  
  if ( dInfo ) delete dInfo;

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
  
  if     ( FBUtility::compare(keyword, "Users") )             Cref->Mode.Example = id_Users;
  else if( FBUtility::compare(keyword, "Parallel_Plate_2D") ) Cref->Mode.Example = id_PPLT2D;
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
void FFV::printDomainInfo( cpm_GlobalDomainInfo* dInfo )
{
  const REAL_TYPE *org = dInfo->GetOrigin();
  const REAL_TYPE *pch = dInfo->GetPitch();
  const REAL_TYPE *rgn = dInfo->GetRegion();
  const int       *vox = dInfo->GetVoxNum();
  const int       *div = dInfo->GetDivNum();
  
  cout << "\n####### read parameters ########" << endl;
  cout << " G_org      = " << org[0] << "," << org[1] << "," << org[2] << endl;
  cout << " G_voxel    = " << vox[0] << "," << vox[1] << "," << vox[2] << endl;
  cout << " G_pitch    = " << pch[0] << "," << pch[1] << "," << pch[2] << endl;
  cout << " G_region   = " << rgn[0] << "," << rgn[1] << "," << rgn[2] << endl;
  cout << " G_div      = " << div[0] << "," << div[1] << "," << div[2] << endl;
  cout << " #subdomain = " << dInfo->GetSubDomainNum() << endl;
  for( size_t i=0;i<dInfo->GetSubDomainNum();i++ )
  {
    const cpm_ActiveSubDomainInfo *dom = dInfo->GetSubDomainInfo(i);
    const int *pos  = dom->GetPos();
    const int *bcid = dom->GetBCID();
    cout << "  domain" << i
    << "  pos="   << pos[0]  << "," << pos[1]  << "," << pos[2]
    << "  bcid="  << bcid[0] << "," << bcid[1] << "," << bcid[2]
    <<       ","  << bcid[3] << "," << bcid[4] << "," << bcid[5] << endl;
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
    M.printMatList(mp, mat, C.NoMedium);
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
