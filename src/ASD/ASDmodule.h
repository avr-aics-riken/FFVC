#ifndef _ASD_MODULE_H_
#define _ASD_MODULE_H_

//##################################################################################
//
// FFV-C ASD module : Frontflow / violet Cartesian Active SubDomain
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
 * @file   SubDomain.h
 * @brief  ASD Class Header
 * @author aics
 * @note ASD moduleは逐次動作であるが，MPIPolylibを利用．対象パラメータは以下
 *
 * DomainInfo {
 *   GlobalOrigin   = (-8.0e+00, -7.0e-01, -3.0e-01)
 *   GlobalRegion   = (1.12e+01, 8.2e+00, 1.6e+00)
 *   GlobalVoxel    = (1120   , 820   , 160   )
 *   GlobalDivision = (10    , 6    , 4    )
 *   ActiveSubDomainFile = "curved.svx"
 *   Source = "polylib.tp"
 * }
 */

#include "../FB/DomainInfo.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <iostream>
#include <math.h>
#include "../FB/mydebug.h"
#include "SubDomain.h"
#include "../FB/VoxInfo.h"

#include "TextParser.h"
#include "Polylib.h"
#include "MPIPolylib.h"

#include "Cutlib.h"
#include "GridAccessor/Cell.h"

using namespace std;
using namespace PolylibNS;
using namespace ASubDomain;
using namespace cutlib;

class ASD : public DomainInfo {
  
  enum type_medium {
    md_fluid = 1,
    md_solid
  };
  
private:
  ///< int G_division[3];     プロセス分割数   GlobalDivision
  ///< int size[3];           全領域ボクセル数 GlobalVoxel
  ///< REAL_TYPE G_origin[3]; 領域基点        GlobalOrigin
  ///< REAL_TYPE G_region[3]; 全領域サイズ     GlobalRegion
  ///< REAL_TYPE pitch[3];    格子幅
  /// 以上はDomainInfoで定義済み
  
  float sd_rgn[3];   ///< subdomainサイズ
  unsigned G_Acell;  ///< グローバルなActive cell
  
  MPIPolylib* PL;     ///< Polylibクラス
  
  int FillSeedDir;
  int guide;
  
  // カット
  CutPos32Array *cutPos;
  CutBid5Array  *cutBid;
  
  float  *d_cut; ///< 距離情報
  int    *d_bid; ///< BC
  int    *d_bcd;
  
public:
  
  // default constructor
  ASD() {
    G_Acell = 0;
    FillSeedDir=-1;
    guide = 1;
    
    for (int i=0; i<3; i++)
    {
      sd_rgn[i] = 0.0;
    }
  };
  
  ~ASD() {
  };
  
  
public:
  // @brief Active SubDomainを作成し，統計情報を表示する
  // @param [in] argc  引数の数
  // @param [in] argv  引数文字列
  void evaluateASD(int argc, char **argv);
  
  
  /**
   * @brief CPMのポインタをコピーし、ランク情報を設定
   * @param [in] m_paraMngr  cpm_ParaManagerクラス
   * @return  エラーコード
   */
  bool importCPM(cpm_ParaManager* m_paraMngr)
  {
    if ( !m_paraMngr ) return false;
    paraMngr = m_paraMngr;
    
    setRankInfo(paraMngr, procGrp);
    
    return true;
  }
  
  
private:
  // @brief active subdomain flag
  int active(const float* px,
             const float* py,
             const float* pz,
             unsigned char* sd);
  
  
  // @brief カット計算
  void CalculateCut();
  
  
  // @brief position of min/max for each subdomain
  void createSubdomainTable(float* p_x, float* p_y, float* p_z);
  
  
  // @brief fill
  void fill(bool disp_flag);
  
  
  // @brief サブドメイン内に含まれるポリゴンリストを検索し，フラグを立てる
  void findPolygon(const float* px,
                   const float* py,
                   const float* pz,
                   unsigned char* sd,
                   const string label);
  
  
  // @brief グローバルな領域情報を取得
  void getDomainInfo(TextParser* tpCntl, bool flag);
  
  
  // @brief 幾何形状情報のロード
  void setupPolygonASD(const string fname, bool flag);
};

#endif // _ASD_MODULE_H_