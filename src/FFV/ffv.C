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
 * @file ffv.C
 * @brief FFV Class
 * @author kero
 */

#include "ffv.h"

// コンストラクタ
FFV::FFV()
{
  session_maxStep = 0;
  session_currentStep = 0;
  
  mp = stdout;
  
  
  paraMngr = NULL;
  
}

// デストラクタ
FFV::~FFV()
{
  
  if ( mp ) fclose(mp);
  
}


// タイムステップループ
int FFV::Loop(int m_step)
{
  return 1;
}


// タイムステップループ
int FFV::MainLoop()
{
  int ret = 1;
  
  for (int i=1; i<=session_maxStep; i++)
  {
    session_currentStep = i;
    
    int loop_ret = Loop(i);
    
    switch (loop_ret) 
    {
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


// 終了時の処理
bool FFV::Post() 
{
/*
  TIMING__ 
  { 
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
    
    if ( IsMaster(paraMngr) ) 
    {
      // 結果出力(排他測定のみ)
      PM.print(stdout);
      PM.print(fp);
      
      // 結果出力(非排他測定も)
      if ( C.Mode.Profiling == DETAIL) 
      {
        PM.printDetail(stdout);
        PM.printDetail(fp);
      }
      
      if ( !fp ) fclose(fp);
    }
  }
  
  if ( IsMaster(paraMngr) ) 
  {
    if( cm_mode == 0 )
    {
      printf( "Communication Mode = CommBndCell\n" );
    } 
    else if( cm_mode == 1 )
    {
      printf( "Communication Mode = CommBndCell2 or cbs3d_commface(no hide)\n" );
    } 
    else 
    {
      printf( "Communication Mode = CommBndCell2 or cbs3d_commface(hide)\n" );
    }
    paraMngr->flush(stdout);
  }
  */
  return true;
}


// タイムステップループの後の処理
bool FFV::stepPost() 
{
  return true;
}


// 利用例の表示
void FFV::Usage(void)
{
  FBUtility::printVersion(mp, "Frontflow/violet", VERS_FFV);
  
  cout << " Usage : ";
  cout << "ffv"
  << " parameter_file" << endl;
  cout << endl;
  
  cout << " \tparameter_file includes all parameters for simulation." << endl;
  cout << endl;

}