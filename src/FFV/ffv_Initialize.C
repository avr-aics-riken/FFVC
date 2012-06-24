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
  unsigned long TotalMemory;    ///< 計算に必要なメモリ量（ローカル）
  unsigned long PrepMemory;     ///< 初期化に必要なメモリ量（ローカル）
  unsigned long G_TotalMemory;  ///< 計算に必要なメモリ量（グローバル）
  unsigned long G_PrepMemory;   ///< 初期化に必要なメモリ量（グローバル）
  unsigned long tmp_memory;     ///< 計算に必要なメモリ量（グローバル）？
  
  int ret = 0;
  
  // CPMバージョン表示
  if ( paraMngr->GetMyRankID() == 0 )
  {
    cpm_Base::VersionInfo();
  }
  
  
  // 固定パラメータ
  fixed_parameters();
  
  // 前処理段階のみに使用するオブジェクトをインスタンス
  //VoxInfo Vinfo;
  ParseBC     B;
  ParseMat    M;
  
  
  // CPMのセット
  //BC.importCPM(paraMngr);
  //Vinfo.importCPM(paraMngr);
  B.importCPM(paraMngr);
  F.importCPM(paraMngr);
  
  
  // 並列処理モード
  setParallelism();

  
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
  

  int ierror;
  std::string inputfile("cavity.tp");
  
  // TextParserのインスタンス生成
  ierror = tpCntl.getTPinstance();
  
  // TextParserのファイル読み込み
  ierror = tpCntl.readTPfile(inputfile);
  
  // TPControlクラスのポインタを各クラスに渡す
  if ( !C.importTP(&tpCntl) ) {
    Hostonly_ stamped_printf("\tError during sending an object pointer of TPControl to Control class\n");
    return -1;
  }
  if ( !B.importTP(&tpCntl) ) {
    Hostonly_ stamped_printf("\tError during sending an object pointer of TPControl to ParseBC class\n");
    return -1;
  }
  if ( !M.importTP(&tpCntl) ) {
    Hostonly_ stamped_printf("\tError during sending an object pointer of TPControl to ParseMat class\n");
    return -1;
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  /*
  
  // 入力ファイルリスト
  vector<const char*> ifname;
  for (int i=1; i<argc; i++) 
  {
    ifname.push_back(argv[i]);
  }
  
  // 情報のプリント
  paraMngr->flush(cout);
  if ( paraMngr->GetMyRankID()==0 )
  {
    //入力となる領域分割ファイル
    cout << "Number of input files = " << ifname.size() << endl;
    for (int i=0; i<ifname.size(); i++)
    {
      cout << "  " << i << " : " << ifname[i] << endl;
    }
    cout << endl;
    
    //実数のサイズ
    if( cpm_Base::RealIsDouble() )
      cout << "REAL_TYPE=double" << endl;
    else
      cout << "REAL_TYPE=float" << endl;
  }
  paraMngr->flush(cout);
  
  
  // 領域分割情報の読み込み
  vector<cpm_GlobalDomainInfo*> domainInfo;
  for( int i=0;i<ifname.size();i++ )
  {
    // 領域分割情報の読み込み
    cpm_GlobalDomainInfo *dInfo = cpm_TextParserDomain::Read(ifname[i], ret);
    if( !dInfo && ret == TP_NO_ERROR )
    {
      delete dInfo;
      ret = CPM_ERROR_INVALID_PTR;
    }
    if( ret != TP_NO_ERROR )
    {
      for( int j=0;j<domainInfo.size();j++ )
      {
        delete domainInfo[j];
      }
      cerr << "TextParser error : " << ret << endl;
      return ret;
    }
    
    //リストに追加
    domainInfo.push_back(dInfo);
    
#ifdef _DEBUG
    printDomainInfo(dInfo);
    paraMngr->flush(cout);
#endif
  }
  
  // 領域分割
  for( int i=0;i<domainInfo.size();i++ )
  {
    if( i==0 )
    {
      // プロセスグループ0で領域分割
      if( (ret = paraMngr->VoxelInit( domainInfo[i], 3, 6 )) != CPM_SUCCESS )
      {
        cerr << "Domain decomposition error : " << ret << endl;
        return ret;
      }
    }
    else
    {
      // プロセスグループを生成
      int nproc = domainInfo[i]->GetSubDomainNum();
      int *proclist = new int[nproc];
      for( int j=0;j<nproc;j++ ) proclist[j] = j*2;
      int procGrpNo = paraMngr->CreateProcessGroup( nproc, proclist );
      
      // 領域分割
      if( procGrpNo >= 0 )
      {
        if( (ret = paraMngr->VoxelInit( domainInfo[i], 3, 3, procGrpNo )) != CPM_SUCCESS )
        {
          cerr << "Domain decomposition error : " << ret << endl;
          return ret;
        }
      }
    }
  }
#ifdef _DEBUG
  paraMngr->printVoxelInfo();
  paraMngr->printVoxelInfo(paraMngr->GetMyRankID());
  paraMngr->flush(cout);
#endif
  
  ////// MPI test //////
#include "include/_mpi_test.h"
  
  paraMngr->Barrier();
  
  return CPM_SUCCESS;
  */
  
  
  return 1;
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
  strcpy(C.HistoryName,        "history_base.txt");
  strcpy(C.HistoryCompoName,   "history_compo.txt");
  strcpy(C.HistoryDomfxName,   "history_domainflux.txt");
  strcpy(C.HistoryForceName,   "history_force.txt");
  strcpy(C.HistoryWallName,    "history_log_wall.txt");
  strcpy(C.HistoryItrName,     "history_iteration.txt");
  strcpy(C.HistoryMonitorName, "sample.log");
  
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



// 並列化と分割の方法を保持
void FFV::setParallelism()
{
  std::string para_mode;
  
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
  
}
