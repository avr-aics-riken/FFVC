// #################################################################
//
// FFV : Frontflow / violet
//
// Copyright (c) All right reserved. 2012
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
//
// #################################################################


//@file ffv.C
//@brief ffv class
//@author kero

#include "ffv.h"

//@brief コンストラクタ
FFV::FFV()
{
  
  mp = stdout;
  
  
  
  
}

//@brief デストラクタ
FFV::~FFV()
{
  
  if ( mp ) fclose(mp);
  
}


//@brief タイムステップループ
int FFV::Loop(int m_step)
{
  
  
  
}


//@brief タイムステップループ
int FFV::MainLoop(void)
{
  int ret = 1;
  
  for (int i=1; i<=session_maxStep; i++){
    
    session_currentStep = i;
    
    int loop_ret = Loop(i);
    
    switch (loop_ret) {
      case -1: // error
        return -1;
        
      case 0: // forced terminated
        ret = 1;
        break;
        
      case 1: // normal
        ret = 1;
        break;
        
      default:
        return -1;
    }
    
    if( loop_ret == 0 ) break;
  }
  
  if ( !stepPost() ) return -1;
  
  return ret;
}


//@brief 終了時の処理
bool FFV::Post(cpm_ParaManager *paraMngr) 
{

  TIMING__ { 
    FILE* fp = NULL;
    
    if ( IsMaster(paraMngr) ) {
      if ( !(fp=fopen("profiling.txt", "w")) ) {
        stamped_printf("\tSorry, can't open 'profiling.txt' file. Write failed.\n");
        Exit(0);
      }
    }
    
    // 測定結果の集計(gathreメソッドは全ノードで呼ぶこと)
    TIMING_start(tm_statistic);
    PM.gather();
    TIMING_stop(tm_statistic, 0.0);
    
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
  
  if ( IsMaster(paraMngr) ) {
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


//@brief タイムステップループの後の処理
bool FFV::stepPost(void) 
{
  return true;
}

//@brief 利用例の表示
void FFV::Usage(void)
{
  char version[8];
  sprintf(version, "%s", VERS_FFV);
  
  cout << endl;
  cout << "\tFrontflow / violet" << endl;
  cout << "\t=======================================================" << endl;
  cout << "\tversion " << version << endl;
  cout << "\t==============" << endl;
  cout << endl;
  
  cout << " Usage : ";
  cout << "ffv"
  << " parameter_file" << endl;
  cout << endl;
  
  cout << " \tparameter_file includes all parameters for simulation." << endl;
  cout << endl;

}