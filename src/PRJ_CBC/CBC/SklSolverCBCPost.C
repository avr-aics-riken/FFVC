/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2010
 *
 */

//@file SklSolverCBCPost.C
//@brief SklSolverCBC class
//@author keno, FSI Team, VCAD, RIKEN

#include "SklSolverCBC.h"

//@fn bool SklSolverCBC::SklSolverPost()
//@brief タイムステップループの後の処理
bool
SklSolverCBC::SklSolverPost() {

  TIMING__ { 
    FILE* fp = NULL;
    
    Hostonly_ {
      if ( !(fp=fopen("profiling.txt", "w")) ) {
        stamped_printf("\tSorry, can't open 'profiling.txt' file. Write failed.\n");
        assert(0);
      }
    }
    
    // 測定結果の集計(gathreメソッドは全ノードで呼ぶこと)
    PM.gather();
    
    Hostonly_ {
      // 結果出力(排他測定のみ)
      PM.print(stdout);
      PM.print(fp);
    
      // 結果出力(非排他測定も)
      if ( C.Mode.Profiling == DETAIL) {
        PM.printDetail(stdout);
        PM.printDetail(fp);
      }
    
     if ( !fp ) fclose(fp);
    }
  }
  
  Hostonly_ {
    if( cm_mode == 0 ){
      printf( "Communication Mode = CommBndCell\n" );
    } else if( cm_mode == 1 ){
      printf( "Communication Mode = CommBndCell2 or cbs3d_commface(no hide)\n" );
    } else {
      printf( "Communication Mode = CommBndCell2 or cbs3d_commface(hide)\n" );
    }
    fflush(stdout);
  }

  return true;
}
