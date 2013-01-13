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
 * @file   ffv_Post.C
 * @brief  FFV Class
 * @author kero
 */

#include "ffv.h"



// 終了時の処理
bool FFV::Post() 
{
  
  TIMING__ 
  { 
    FILE* fp = NULL;
    
    if ( IsMaster() ) 
    {
      if ( !(fp=fopen("profiling.txt", "w")) ) 
      {
        stamped_printf("\tSorry, can't open 'profiling.txt' file. Write failed.\n");
        Exit(0);
      }
    }
    
    // 測定結果の集計(gathreメソッドは全ノードで呼ぶこと)
    TIMING_start(tm_statistic);
    PM.gather();
    TIMING_stop(tm_statistic, 0.0);
    
    if ( IsMaster() ) 
    {
      // 結果出力(排他測定のみ)
      PM.print(stdout, HostName, C.OperatorName);
      PM.print(fp, HostName, C.OperatorName);
      
      // 結果出力(非排他測定も)
      if ( C.Mode.Profiling == DETAIL) 
      {
        PM.printDetail(stdout);
        PM.printDetail(fp);
      }
      
      if ( !fp ) fclose(fp);
    }
  }
  
  
  return true;
}

