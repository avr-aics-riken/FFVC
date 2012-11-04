#ifndef _DFIINFO_H_
#define _DFIINFO_H_

// #################################################################
//
// Combine sph files and output 
//
// Copyright (c) All right reserved. 2012
//
// Institute of Industrial Science, University of Tokyo, Japan. 
//
// #################################################################

/** 
 * @file   dfiinfo.h
 * @brief  DfiInfo Class Header
 * @author kero
 */

//////#include <stdio.h>
//////#include <stdlib.h>
//////#include <string.h>
//////#include <string>
//////#include <iostream>
//////#include <fstream>

#include "cpm_Define.h"

#include "DomainInfo.h"
#include "FB_Define.h"
#include "mydebug.h"
#include "TPControl.h"

using namespace std;


class DifferentRestartInfo { // 異なる並列数からのリスタートのための情報保持クラス

public:

  int rank;
  int overlap_head[3];
  int overlap_tail[3];

  int read_file_voxel_head[3];
  int read_file_voxel_tail[3];
  int read_file_voxel_size[3];
  string f_prs;
  string f_vel;
  string f_temp;

public:

  /** コンストラクタ */
  DifferentRestartInfo(){
    rank=-1;
    for(int i=0;i<3;i++) overlap_head[i]=0;
    for(int i=0;i<3;i++) overlap_tail[i]=0;
    for(int i=0;i<3;i++) read_file_voxel_head[i]=0;
    for(int i=0;i<3;i++) read_file_voxel_tail[i]=0;
    for(int i=0;i<3;i++) read_file_voxel_size[i]=0;
  };

  /**　デストラクタ */
  ~DifferentRestartInfo(){};


};


class DfiInfo : public DomainInfo {

public:

  typedef struct 
  {
    int RankID;
    string HostName;
    int VoxelSize[3];
    int HeadIndex[3];
    int TailIndex[3];
    long IJK;    //(I,J,K)一次元アドレス変換
    long IJK_JK; //一次元アドレスY面
	  long IJK_K;  //一次元アドレスZ面
  } NodeInfo;

  string Prefix;
  int RankID_in_MPIworld;
  int GroupID_in_MPIworld;
  int Number_of_Rank_in_MPIworld;
  int Number_of_Group_in_MPIworld;
  int Global_Voxel[3];
  int Global_Division[3];
  string FileFormat;
  int GuideCell;
  NodeInfo *Node;
  vector<int> step;

  int NodeInfoSize;//=RankID_in_MPIworld
  vector<int> index_y;
  vector<int> index_z;

public:
  /** コンストラクタ */
  DfiInfo();
  
  /**　デストラクタ */
  ~DfiInfo();
  
  /**
   * @brief dfiファイルの読み込み
   * @param [in] fname  dfiファイル名
   * @param [in] D  DfiInfoクラスポインタ
   */
  void ReadDfiFile(string fname);

  /**
   * @brief 内部変数の計算
   */
  void SetValue();

};

#endif // _DFIINFO_H_
