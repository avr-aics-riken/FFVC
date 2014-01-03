#ifndef _FB_DFIINFO_H_
#define _FB_DFIINFO_H_

//##################################################################################
//
// Flow Base class
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2014 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

/** 
 * @file   dfiinfo.h
 * @brief  DfiInfo Class Header
 * @author kero
 */

#include "cpm_Define.h"

#include "DomainInfo.h"
#include "FB_Define.h"
#include "mydebug.h"
#include "TextParser.h"

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
  string f_fvel;
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

  typedef struct
  {
    int step;
    REAL_TYPE time;
    //MinMax...
  } Slice;
  
  //FileInfo
  string Prefix;
  string FileFormat;
  int GuideCell;
  string DataType;
  string Endian;
  string ArrayShape;
  int Component;
  
  //Unit
  string Length;
  REAL_TYPE L0;
  string Velocity;
  REAL_TYPE V0;
  string Pressure;
  REAL_TYPE P0;
  REAL_TYPE DiffPrs;
  
  //FilePath
  string dfi_proc;
  
  
  //Process
  int Global_Origin[3];
  int Global_Region[3];
  int Global_Voxel[3];
  int Global_Division[3];
  //int RankID_in_MPIworld;
  //int GroupID_in_MPIworld;
  int NumberOfRank;
  int NumberOfGroup;
  NodeInfo *Node;
  

  //TimeSlice
  vector<Slice*> Sc;
  
  //others...
  int dim;
  int NodeInfoSize;
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
   */
  void ReadDfiFile(string fname);
  
  /**
   * @brief Sliceのセット
   * @param [in] m_step
   * @param [in] m_time
   */
  void SetSlice(int m_step, REAL_TYPE m_time);
  
  /**
   * @brief dfi_procファイルの読み込み
   * @param [in] fname  dfiファイル名
   */
  void ReadDfiProc(string fname);

  /**
   * @brief 内部変数の計算
   */
  void SetValue();

};

#endif // _FB_DFIINFO_H_
