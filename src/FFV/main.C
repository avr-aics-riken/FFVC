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
 * @file main.C
 * @brief ffvのmain関数
 * @author kero
 */

#include "ffv.h"

#ifdef _DEBUG
  #include "include/_debug.h"
#endif



int main( int argc, char **argv )
{

  double init_str, init_end;
  double main_str, main_end;
  double post_str, post_end;
  
  
  int ret = 0;

  // FFV classのインスタンス
  FFV ffv;

  // ##################################################################
  // 初期化
  init_str = cpm_Base::GetWTime();
  
  // 並列管理クラスのインスタンスと初期化
  // ここでMPI_Initも行う
  ffv.paraMngr = cpm_ParaManager::get_instance(argc, argv);
  if ( !ffv.paraMngr ) return CPM_ERROR_PM_INSTANCE;
  
  
  if ( ffv.paraMngr->GetMyRankID() == 0 )
  {
    cpm_Base::VersionInfo();
  }
  
  
  int init_ret = ffv.Initialize(argc, argv);
  
  switch( init_ret )
  {
    case -1:
      printf("\n\tSolver initialize error.\n\n");
      return 0;

    case 0:
      printf("\n\tForced termination during initialization.\n\n");
      return 1;
      
    case 1:
      // keep going on processing
      break;
      
    default:
      printf("\n\tSolver initialize error.\n\n");
      return 0;
  }
  
  init_end = cpm_Base::GetWTime();
  
  
  // ##################################################################
  // タイムステップループ
  main_str = cpm_Base::GetWTime();
  
  int loop_ret = ffv.MainLoop();
  
  switch (loop_ret) 
  {
    case -1:
      printf("\n\tSolver error.\n\n");
      return 0;

    case 0:
      printf("\n\tSolver forced termination time-step loop.\n\n");
      break;
      
    case 1:
      if ( ffv.IsMaster() ) 
      {
        printf("\n\tSolver finished.\n\n");
      }
      break;
  }
  
  main_end = cpm_Base::GetWTime();

  
  
  // ##################################################################
  // ポスト処理
  post_str = cpm_Base::GetWTime();
  
  if( !ffv.Post() )
  {
    printf("solver post error.\n");
    return false;
  }
  
  post_end = cpm_Base::GetWTime();

  if ( ffv.IsMaster() ) 
  {
    double init = init_end - init_str;
    double main = main_end - main_str;
    double post = post_end - post_str;
    
    printf("TIME : Solver Init  %10.3f sec.\n", init);
    printf("TIME : Solver Main  %10.3f sec.\n", main);
    printf("TIME : Solver Post  %10.3f sec.\n", post);
    printf("TIME : Solver Total %10.3f sec.\n", init+main+post);
  }
  
  return 0;
}
