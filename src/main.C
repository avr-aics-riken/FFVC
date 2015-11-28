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
// Copyright (c) 2012-2015 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

/** 
 * @file   main.C
 * @brief  ffvのmain関数
 * @author aics
 */

#include "FFV/ffv.h"
#include "ASD/ASDmodule.h"


// HPCPF status
FILE* fp_hpcpf = NULL;

#define hpcpf_status(x) \
(fprintf(fp_hpcpf, "status code = %d\nexit at %s:%u\n", x, __FILE__, __LINE__))


// return; 0 - normal
//         1 - others
int main( int argc, char **argv )
{
  // Version info
  
  if ( !strcasecmp(argv[1], "--version"))
  {
    printf("FFV-C : Frontflow / violet Cartesian  Version : %s : %s\n", FFVC_VERSION_NO, FFVC_REVISION);
    exit(0);
  }
  
  
  // タイミング用変数
  double init_str, init_end;
  double main_str, main_end;
  double post_str, post_end;
  
  // classのインスタンス
  FFV ffv;
  ASD asd;
  
  
  // FFVCのモード
  if ( !strcasecmp(argv[1], "--filter") )
  {
    ffv.EXEC_MODE = ffvc_filter;
  }
  else if ( !strcasecmp(argv[1], "--asd") )
  {
    ffv.EXEC_MODE = ffvc_asd;
  }
  else
  {
    ffv.EXEC_MODE = ffvc_solver; // ffv_Initialize() でffv_solver, ffv_solverASを選別
  }

  
  
  // 並列管理クラスのインスタンスと初期化
  // ここでMPI_Initも行う

  if ( (ffv.EXEC_MODE == ffvc_solver) || (ffv.EXEC_MODE == ffvc_filter) )
  {
    if ( !ffv.importCPM(cpm_ParaManager::get_instance(argc, argv)) ) return 1;
  }
  else //  ffvc_asd
  {
    if ( !asd.importCPM(cpm_ParaManager::get_instance(argc, argv)) ) return 1;
  }
  
  
  
  // 引数チェック
  bool usage_flag = true;
  if ( ffv.EXEC_MODE == ffvc_solver )
  {
    if ( argc != 2 ) usage_flag = false;
  }
  else if ( ffv.EXEC_MODE == ffvc_filter )
  {
    if ( argc != 3 ) usage_flag = false;
  }
  else // ffvc_asd
  {
    if ( (argc != 3) && (argc != 4) && (argc != 6) ) usage_flag = false;
  }
  
  if ( !usage_flag )
  {
    if ( ffv.IsMaster() )
    {
      printf("\n\tusage\n");
      printf("\n\t$ ffvc <parameter_file>\n");
      printf("\n\t$ ffvc --version\n");
      printf("\n\t$ ffvc --filter <parameter_file>\n");
      printf("\n\t$ ffvc --asd hoge.tp         : In case of global division is described in hoge.tp\n");
      printf("\n\t$ ffvc --asd hoge.tp 300     : command line mode with total number of division\n");
      printf("\n\t$ ffvc --asd hoge.tp 10 5 32 : command line mode with number of division for each direction\n");
    }
    
    if (cpm_ParaManager::get_instance()->GetMyRankID()==0) hpcpf_status(1);
    return 1;
  }

  // ##################################################################
  // 初期化
  init_str = cpm_Base::GetWTime();
  
  
  // Open HPCPF_STATUS file
  if (cpm_ParaManager::get_instance()->GetMyRankID()==0)
  {
    if ( !(fp_hpcpf=fopen("HPCPF_EXIT_STATUS", "w")) )
    {
      hpcpf_status(1);
      return 1;
    }
  }
  
  
  int init_ret;
  
  if ( ffv.EXEC_MODE == ffvc_solver )
  {
    init_ret = ffv.Initialize(argc, argv);
  }
  else if ( ffv.EXEC_MODE == ffvc_filter )
  {
    init_ret = ffv.FilterInitialize(argc, argv);
  }
  else // ffvc_asd
  {
    asd.evaluateASD(argc, argv);
    
    return 0;
  }
  
  switch( init_ret )
  {
    case -1:
      if ( ffv.IsMaster() ) printf("\n\tSolver initialize error.\n\n");
      if (cpm_ParaManager::get_instance()->GetMyRankID()==0) hpcpf_status(1);
      return 1;
      break;

    case 0:
      if ( ffv.IsMaster() ) printf("\n\tForced termination during initialization.\n\n");
      if (cpm_ParaManager::get_instance()->GetMyRankID()==0) hpcpf_status(1);
      return 1;
      break;
      
    case 1:
      // keep going on processing
      break;
      
    default:
      if ( ffv.IsMaster() ) printf("\n\tSolver initialize error.\n\n");
      if (cpm_ParaManager::get_instance()->GetMyRankID()==0) hpcpf_status(1);
      return 1;
      break;
  }
  
  init_end = cpm_Base::GetWTime();
  
  // シグナルハンドラの初期化
  FFV_TerminateCtrl::initialize();
  
  
  
  // ##################################################################
  // タイムステップループ
  main_str = cpm_Base::GetWTime();
  
  int loop_ret;
  
  if ( ffv.EXEC_MODE == ffvc_solver )
  {
    loop_ret = ffv.MainLoop();
    
    switch (loop_ret)
    {
      case -1:
        if ( ffv.IsMaster() ) printf("\n\tSolver error.\n\n");
        if (cpm_ParaManager::get_instance()->GetMyRankID()==0) hpcpf_status(1);
        break;
        
      case 0:
        if ( ffv.IsMaster() ) printf("\n\tSolver forced termination time-step loop.\n\n");
        if (cpm_ParaManager::get_instance()->GetMyRankID()==0) hpcpf_status(1);
        break;
        
      case 1:
        if ( ffv.IsMaster() ) printf("\n\tSolver finished.\n\n");
        if (cpm_ParaManager::get_instance()->GetMyRankID()==0) fprintf(fp_hpcpf, "status code = 0\n");
        break;
    }
  }
  else
  {
    loop_ret = ffv.FilterLoop();
  }

  
  main_end = cpm_Base::GetWTime();

  
  
  // ##################################################################
  // ポスト処理
  post_str = cpm_Base::GetWTime();
  
  if( !ffv.Post() )
  {
    printf("solver post error.\n");
    if (cpm_ParaManager::get_instance()->GetMyRankID()==0) hpcpf_status(1);
    return 1;
  }
  
  post_end = cpm_Base::GetWTime();

  if ( ffv.IsMaster() ) 
  {
    double init = init_end - init_str;
    double main = main_end - main_str;
    double post = post_end - post_str;
    
    printf("\n\n");
    printf("TIME : Solver Init  %10.3f sec.\n", init);
    printf("TIME : Solver Main  %10.3f sec.\n", main);
    printf("TIME : Solver Post  %10.3f sec.\n", post);
    printf("TIME : Solver Total %10.3f sec.\n", init+main+post);
  }
  
  if ( loop_ret != 1 )
  {
    if (cpm_ParaManager::get_instance()->GetMyRankID()==0) hpcpf_status(1);
    return 1;
  }

  // Normal return code
  if (cpm_ParaManager::get_instance()->GetMyRankID()==0) fprintf(fp_hpcpf, "status code = 0\n");
  
  return 0;
}
