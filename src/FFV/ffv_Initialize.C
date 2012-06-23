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
  int ret = 0;
  
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
  
  
  
  
}