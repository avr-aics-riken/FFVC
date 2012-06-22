// #################################################################
//
// FFV : Frontflow / violet
//
// Copyright (c) All right reserved. 2012
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
//
// #################################################################

#include "ffv.h"

#ifdef _DEBUG
  #include "include/_debug.h"
#endif






int main( int argc, char **argv )
{

#ifdef TIME_MEASURE
  double init_str, init_end;
  double main_str, main_end;
  double post_str, post_end;
#endif 
  
  
  int ret = 0;


  // ##################################################################
  // 初期化
#ifdef TIME_MEASURE
  init_str = cpm_Base::GetWTime();
#endif
  
  
  // 並列管理クラスのインスタンスと初期化
  // ここでMPI_Initも行う
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance(argc, argv);
  if( !paraMngr ) return CPM_ERROR_PM_INSTANCE;
  
  
  if( paraMngr->GetMyRankID() == 0 )
  {
    cpm_Base::VersionInfo();
  }
  
  // FFV classのインスタンス
  FFV ffv;
  
  int init_ret = ffv.Initialize();
  
  switch( init_ret ){
    case -1:
      fprintf(mp, "\n\tSolver initialize error.\n\n");
      return 0;

    case 0:
      fprintf(mp, "\n\tForced termination during initialization.\n\n");
      return 1;
      
    case 1:
      // keep going on processing
      break;
      
    default:
      fprintf(mp, "\n\tSolver initialize error.\n\n");
      return 0;
  }
  
#ifdef TIME_MEASURE
  init_end = cpm_Base::GetWTime();
#endif
  
  
  // ##################################################################
  // タイムステップループ
#ifdef TIME_MEASURE
  main_str = cpm_Base::GetWTime();
#endif
  
  int loop_ret = ffv.MainLoop();
  
  switch (loop_ret) {
    case -1:
      fprintf(mp, "\n\tSolver error.\n\n");
      return 0;

    case 0:
      fprintf(mp, "\n\tSolver forced termination time-step loop.\n\n");
      break;
      
    case 1:
      if ( ffv.IsMaster(paraMngr) ) {
        fprintf(mp, "\n\tSolver finished.\n\n");
      }
      break;
  }
  
  
  
#ifdef TIME_MEASURE
  main_end = cpm_Base::GetWTime();
#endif
  
  
  // ##################################################################
  // ポスト処理
#ifdef TIME_MEASURE
  post_str = cpm_Base::GetWTime();
#endif
  
  if( !ffv.Post() ){
    SklErrMessage("solver post error.\n");
    return false;
  }
  
#ifdef TIME_MEASURE
  post_end = cpm_Base::GetWTime();
#endif
  
  
  
#ifdef TIME_MEASURE
  if( !ParaCmpo->IsParallel() || (ParaCmpo->GetMyID() == 0) ){
    SklMessage("TIME : Solver Init  %10.3f sec.\n", (init_end - init_str));
    SklMessage("TIME : Solver Main  %10.3f sec.\n", (main_end - main_str));
    SklMessage("TIME : Solver Post  %10.3f sec.\n", (post_end - post_str));
    SklMessage("TIME : Solver Total %10.3f sec.\n",
               (init_end-init_str)+(main_end-main_str)+(post_end-post_str));
  }
#endif // SKL_TIME_MEASURED
  
  
  
  
  
  
  
  // 入力ファイルリスト
  vector<const char*> ifname;
  for( int i=1;i<argc;i++ )
  {
    ifname.push_back(argv[i]);
  }

  // ちょっとした情報のプリント
  paraMngr->flush(cout);
  if( paraMngr->GetMyRankID()==0 )
  {
    //入力となる領域分割ファイル
    cout << "number of input files=" << ifname.size() << endl;
    for( int i=0;i<ifname.size();i++ )
    {
      cout << "  " << i << " : " << ifname[i] << endl;
    }

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
    cpm_GlobalDomainInfo *dInfo = cpm_TextParserDomain::Read(ifname[i],ret);
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
        cerr << "VoxelInit error : " << ret << endl;
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
          cerr << "VoxelInit error : " << ret << endl;
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
  double elapse = cpm_Base::GetWSpanTime(ts);
  if( paraMngr->GetMyRankID()==0 ) cout << endl << "### elapse = "<< elapse << " [sec]" << endl;

  return CPM_SUCCESS;
}

