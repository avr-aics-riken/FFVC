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
// Copyright (c) 2012-2016 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

/**
 * @file   ffv_Post.C
 * @brief  FFV Class
 * @author aics
 */

#include "ffv.h"



// 終了時の処理
bool FFV::Post()
{
  FILE* fp = NULL;


  // 統計情報
  if (C.Mode.Statistic == ON)
  {
    Hostonly_
    {
      if ( !H->printCompoStatistics(cmp, cmp_force_avr) ) Exit(0);
    }
  }


  TIMING__
  {
    fp = NULL;

    Hostonly_
    {
      if ( !(fp=fopen("profiling.txt", "w")) )
      {
        stamped_printf("\tSorry, can't open 'profiling.txt' file. Write failed.\n");
        Exit(0);
      }
    }

    // 測定結果の集計(gathreメソッドは全ノードで呼ぶこと)
    TIMING_start("Statistic");
    PM.gather();
    TIMING_stop("Statistic", 0.0);

    Hostonly_
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
