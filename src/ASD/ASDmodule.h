#ifndef _ASD_MODULE_H_
#define _ASD_MODULE_H_

//##################################################################################
//
// FFV-C ASD module : Frontflow / violet Cartesian Active SubDomain
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
 *   Source = "polylib.tp"
 *   outputSVX = "hoge.svx"
 *   outputSubdomain = "no"
 * }
 */

#include "DomainInfo.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <iostream>
#include <math.h>
#include "mydebug.h"
#include "SubDomain.h"
#include "Geometry.h"
#include "PolyProperty.h"

#include "TextParser.h"
#include "Polylib.h"
#include "MPIPolylib.h"

#include "Vec3.h"

using namespace std;
using namespace PolylibNS;
using namespace ASubDomain;
using namespace Vec3class;

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
  
  REAL_TYPE sd_rgn[3];   ///< subdomainサイズ
  unsigned G_Acell;  ///< グローバルなActive cell
  
  string out_svx;   ///< SVX 形式出力 >> VXgen
  string out_sub;   ///< char形式出力 >> FXgen
  
  MPIPolylib* PL;     ///< Polylibクラス
  
  int FillSeedDir;
  int guide;
  int divPolicy;
  
  
  long long *d_cut; ///< 距離情報
  int    *d_bid; ///< BC
  int    *d_bcd;
  
  PolygonProperty* PG;
  Geometry GM;  ///< Geometry class
  
public:
  
  // default constructor
  ASD() {
    G_Acell = 0;
    FillSeedDir=-1;
    guide = 1;
    divPolicy = -1;
    PG = NULL;
    
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
  // active subdomain flag
  int active(const REAL_TYPE* px,
             const REAL_TYPE* py,
             const REAL_TYPE* pz,
             unsigned char* sd);
  
  
  // カット計算
  void CalculateCut();
  
  
  // position of min/max for each subdomain
  void createSubdomainTable(REAL_TYPE* p_x, REAL_TYPE* p_y, REAL_TYPE* p_z);
  
  
  // フィル　Geometryクラスを適用
  void fill(bool disp_flag, Geometry* GM);
  
  
  // サブドメイン内に含まれるポリゴンリストを検索し，フラグを立てる
  void findPolygon(const REAL_TYPE* px,
                   const REAL_TYPE* py,
                   const REAL_TYPE* pz,
                   unsigned char* sd,
                   const string label);
  
  
  // グローバルな領域情報を取得
  void getDomainInfo(TextParser* tpCntl, bool flag);
  
  
  // 幾何形状情報のロード
  void setupPolygonASD(const string fname, bool flag);
  
  
  // FXgenのソースより移動
  //
  // 通信面コストの計算 I,J,K分割を行った時の通信点数の総数を取得する
  inline
  unsigned long long
  CalcCommSize(const unsigned long long iDiv,
               const unsigned long long jDiv,
               const unsigned long long kDiv,
               const unsigned long long voxSize[3]);
  
  
  // 並列プロセス数からI,J,K方向の分割数を取得する
  // 通信面のトータルサイズが小さい分割パターンを採用する
  bool DecideDivPatternCommSize(const unsigned int divNum,
                                const unsigned int voxSize[3],
                                unsigned int divPttn[3]);
  
  
  // 並列プロセス数からI,J,K方向の分割数を取得する
  // １つのサブドメインが立方体に一番近い分割パターンを採用する
  bool DecideDivPatternCube(const unsigned int divNum,
                            const unsigned int voxSize[3],
                            unsigned int divPttn[3]);
  
  
  // I,J,K分割を行った時のI,J,Kボクセル数の最大/最小の差を取得する
  long long CheckCube(const unsigned long long iDiv,
                      const unsigned long long jDiv,
                      const unsigned long long kDiv,
                      const unsigned long long voxSize[3]);
};

#endif // _ASD_MODULE_H_
