#ifndef _FB_DOMAIN_INFO_H_
#define _FB_DOMAIN_INFO_H_

//##################################################################################
//
// Flow Base class
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
 * @file   DomainInfo.h
 * @brief  FlowBase DomainInfo class Header
 * @author aics
 */

#include "cpm_ParaManager.h"
#include <string>
#include "FB_Define.h"


class DomainInfo {
  
public:
  cpm_ParaManager* paraMngr; ///< Cartesian Partition Manager
  

  int procGrp;         ///< プロセスグループ番号
  int myRank;          ///< 自ノードのランク番号
  int numProc;         ///< 全ランク数
  
  int nID[6];          ///< 隣接ブロックのランク番号
  int head[3];         ///< 開始インデクス（グローバルインデクス, Fortran）
  
  int guide;           ///< ガイドセル数
  int G_division[3];   ///< プロセス分割数
  
  REAL_TYPE pitch[3];     ///< 格子幅 (Non-dimensional)
  REAL_TYPE pitchD[3];    ///< 格子幅 (有次元)
  
  int size[3];            ///< 領域分割数 (Local, Non-dimensional)
  REAL_TYPE origin[3];    ///< 領域基点   (Local, Non-dimensional)
  REAL_TYPE region[3];    ///< 領域サイズ (Local, Non-dimensional)
  REAL_TYPE originD[3];   ///< 領域基点   (Local, 有次元)
  REAL_TYPE regionD[3];   ///< 領域サイズ (Local, 有次元)
  
  int G_size[3];          ///< 領域分割数 (Global, Non-dimensional)
  REAL_TYPE G_origin[3];  ///< 領域基点   (Global, Non-dimensional)
  REAL_TYPE G_region[3];  ///< 領域サイズ (Global, Non-dimensional)
  REAL_TYPE G_originD[3]; ///< 領域基点   (Global, 有次元)
  REAL_TYPE G_regionD[3]; ///< 領域サイズ (Global, 有次元)
  
  std::string HostName;   ///< ホスト名

  
  /** コンストラクタ */
  DomainInfo() {
    procGrp = 0;
    myRank  = -1;
    numProc = 0;
    for (int i=0; i<NOFACE; i++) nID[i] = -1;
    
    for (int i=0; i<3; i++)
    {
      head[i]       = 0;
      size[i]       = 0;
      G_size[i]     = 0;
      G_division[i] = 0;
      pitch[i]      = 0.0;
      origin[i]     = 0.0;
      region[i]     = 0.0;
      G_origin[i]   = 0.0;
      G_region[i]   = 0.0;
      pitchD[i]     = 0.0;
      originD[i]    = 0.0;
      regionD[i]    = 0.0;
      G_originD[i]  = 0.0;
      G_regionD[i]  = 0.0;
    }
    
    guide = 0;
    paraMngr = NULL;
  }
  
  /**　デストラクタ */
  virtual ~DomainInfo() {}
  
  
  /**
   * @brief CPMのポインタをコピーし、ランク情報を設定
   * @param [in] m_paraMngr  cpm_ParaManagerクラス
   * @return  エラーコード
   */
  bool importCPM(cpm_ParaManager* m_paraMngr)
  {
    if ( !m_paraMngr ) return false;
    paraMngr = m_paraMngr;
    
    return true;
  }
  
  /**
   * @brief ランク情報をセットする
   * @param [in] m_paraMngr  CPMlibポインタ
   * @param [in] m_proGrp    プロセスグループ番号
   */
  void setRankInfo(cpm_ParaManager* m_paraMngr, const int m_procGrp)
  {
    // ポインタコピー
    paraMngr = m_paraMngr;
    
    procGrp = m_procGrp;
    myRank  = paraMngr->GetMyRankID(procGrp);
    numProc = paraMngr->GetNumRank(procGrp);
    HostName= paraMngr->GetHostName();
  }
  
  
  /**
   * @brief 隣接情報と領域情報を設定する
   * @param [in] m_guide ガイドセル数
   * @param [in] RefL    代表長
   * @note cpm_ParaManager::VoxelInit()のあとにコール
   *       各パラメータは無次元値
   */
  void setDomainInfo(const int m_guide, const REAL_TYPE RefL)
  {
    // guide cell
    guide = m_guide;
    
    // 領域分割数 (Local)
    const int* l_sz = paraMngr->GetLocalVoxelSize(procGrp);
    size[0] = l_sz[0];
    size[1] = l_sz[1];
    size[2] = l_sz[2];
    
    // 領域基点 (Local)
    const double* l_org = paraMngr->GetLocalOrigin(procGrp);
    origin[0] = (REAL_TYPE)l_org[0];
    origin[1] = (REAL_TYPE)l_org[1];
    origin[2] = (REAL_TYPE)l_org[2];
    
    // 領域サイズ (Local)
    const double* l_reg = paraMngr->GetLocalRegion(procGrp);
    region[0] = (REAL_TYPE)l_reg[0];
    region[1] = (REAL_TYPE)l_reg[1];
    region[2] = (REAL_TYPE)l_reg[2];
    
    // 領域分割数 (Global)
    const int* g_sz  = paraMngr->GetGlobalVoxelSize(procGrp);
    G_size[0] = g_sz[0];
    G_size[1] = g_sz[1];
    G_size[2] = g_sz[2];
    
    // 領域基点   (Global)
    const double* g_org = paraMngr->GetGlobalOrigin(procGrp);
    G_origin[0] = (REAL_TYPE)g_org[0];
    G_origin[1] = (REAL_TYPE)g_org[1];
    G_origin[2] = (REAL_TYPE)g_org[2];
    
    // 領域サイズ (Global)
    const double* g_reg = paraMngr->GetGlobalRegion(procGrp);
    G_region[0] = (REAL_TYPE)g_reg[0];
    G_region[1] = (REAL_TYPE)g_reg[1];
    G_region[2] = (REAL_TYPE)g_reg[2];
    
    // 格子幅
    const double* m_pch = paraMngr->GetPitch(procGrp);
    pitch[0] = (REAL_TYPE)m_pch[0];
    pitch[1] = (REAL_TYPE)m_pch[1];
    pitch[2] = (REAL_TYPE)m_pch[2];
    
    // プロセス分割数
    const int* g_div = paraMngr->GetDivNum(procGrp);
    
    // 隣接情報
    const int* neighbour = paraMngr->GetNeighborRankID(procGrp);
    for (int i=0; i<6; i++)
    {
      nID[i] = neighbour[i];
    }
    
    // 開始インデクス > Fortran index
    const int* hdx = paraMngr->GetVoxelHeadIndex(procGrp);
    head[0] = hdx[0] + 1;
    head[1] = hdx[1] + 1;
    head[2] = hdx[2] + 1;
    
    // 無次元値をもとに有次元値を計算
    for (int i=0; i<3; i++) {
      pitchD[i]    = pitch[i]    * RefL;
      originD[i]   = origin[i]   * RefL;
      regionD[i]   = region[i]   * RefL;
      G_originD[i] = G_origin[i] * RefL;
      G_regionD[i] = G_region[i] * RefL;
    }
    
  }
  
  
  // for Graph Ploter
  const char * _ss(const char *label, double *v, int n=3)
  {
    char tmp_str[256]="";
    if( n == 3 ) sprintf(tmp_str,"%s=%lf %lf %lf", label, v[0], v[1], v[2]);
    if( n == 2 ) sprintf(tmp_str,"%s=%lf %lf", label, v[0], v[1]);
    if( n == 1 ) sprintf(tmp_str,"%s=%lf", label, v[0]);
    if( n == 0 ) sprintf(tmp_str,"%s=", label);
    std::string str = tmp_str;
    return str.c_str();
  }
  
  const char * _ss(const char *label, float *v, int n=3)
  {
    char tmp_str[256]="";
    if( n == 3 ) sprintf(tmp_str,"%s=%f %f %f", label, v[0], v[1], v[2]);
    if( n == 2 ) sprintf(tmp_str,"%s=%f %f", label, v[0], v[1]);
    if( n == 1 ) sprintf(tmp_str,"%s=%f", label, v[0]);
    if( n == 0 ) sprintf(tmp_str,"%s=", label);
    std::string str = tmp_str;
    return str.c_str();
  }
  
};

#endif // _FB_DOMAIN_INFO_H_
