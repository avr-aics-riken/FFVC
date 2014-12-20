#ifndef _FFV_GEOM_H_
#define _FFV_GEOM_H_
//##################################################################################
//
// FFV-C : Frontflow / violet Cartesian
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
 * @file   Geometry.h
 * @brief  Geometry Class header
 * @author aics
 */

#include "cpm_ParaManager.h"

#include "FB_Define.h"
#include "DomainInfo.h"
#include "Control.h"

#include "TextParser.h"

#include "Polylib.h"
#include "MPIPolylib.h"


using namespace std;
using namespace PolylibNS;
using namespace Vec3class;

class Geometry : public DomainInfo {
  
private:
  
  REAL_TYPE m_poly_org[3]; ///< ポリゴンの基点
  REAL_TYPE m_poly_dx[3];  ///< ポリゴンのピッチ

  int NumSuvDiv;       ///< 再分割数
  int FillSeedDir;     ///< フィルのヒント {x_minux | x_plus |...}
  string FillMedium;   ///< フィルに使う媒質 -> int FillID
  string SeedMedium;   ///< ヒントに使う媒質 -> int SeedID
  
  
public:
  
  int FillID;          ///< フィル媒質ID
  int SeedID;          ///< フィルシード媒質ID
  int FillSuppress[3]; ///< PeriodicとSymmetricの外部境界面フィル抑制
  
  
public:
  
  /** コンストラクタ */
  Geometry() {

    FillID = -1;
    SeedID = -1;
    NumSuvDiv = 0;
    FillSeedDir = -1;
    
    for (int i=0; i<3; i++) {
      m_poly_org[i] = 0.0;
      m_poly_dx[i] = 0.0;
      FillSuppress[i] = ON; // default is "fill"
    }
  }
  
  /**　デストラクタ */
  ~Geometry() {}
  

private:
  
  // 点pの属するセルインデクスを求める
  // @param [in]  pt 無次元座標
  // @param [out] w  インデクス
  void findIndex(const Vec3r pt, int* w) const
  {
    REAL_TYPE p[3], q[3];
    p[0] = (REAL_TYPE)pt.x;
    p[1] = (REAL_TYPE)pt.y;
    p[2] = (REAL_TYPE)pt.z;
    
    q[0] = (p[0]-origin[0])/pitch[0];
    q[1] = (p[1]-origin[1])/pitch[1];
    q[2] = (p[2]-origin[2])/pitch[2];
    
    w[0] = (int)ceil(q[0]);
    w[1] = (int)ceil(q[1]);
    w[2] = (int)ceil(q[2]);
  }
  
  
  
protected:
  
  // d_mid[]がtargetであるセルに対して、d_pvf[]に指定値valueを代入する
  unsigned long assignVF(const int target,
                         const REAL_TYPE value,
                         const int* d_mid,
                         REAL_TYPE* d_pvf);
  
  
  // mid[]内にあるm_idのセルを数える
  unsigned long countCellM(const int* mid,
                           const int m_id,
                           const bool painted=true,
                           const int* Dsize=NULL);
  
  
  // 未ペイントセルをtargetでフィル
  unsigned long fillByID(int* mid,
                         const int target,
                         const int* Dsize=NULL);
  
  
  // 流体媒質のフィルをbid情報を元に実行
  unsigned long fillByMid(int* mid,
                          const int tgt_id,
                          const int* Dsize=NULL);
  
  
  // 未ペイントセルを周囲のbidの固体最頻値でフィル
  unsigned long fillByModalSolid(int* bcd,
                                 const int fluid_id,
                                 const int* bid,
                                 const int m_NoCompo);
  
  
  // 未ペイントセルを周囲のmidの固体最頻値でフィル
  unsigned long fillByModalSolid(int* mid, const int fluid_id, const int m_NoCompo);
  
  
  // サブセルのSolid部分の値を代入
  unsigned long fillSubCellSolid(int* smd, REAL_TYPE* svf);
  
  
  // サブセルの未ペイントセルを周囲の媒質IDの固体最頻値でフィル
  unsigned long fillSubCellByModalSolid(int* smd,
                                        const int m_NoCompo,
                                        REAL_TYPE* svf,
                                        const MediumList* mat);
  
  
  // シード点をmid[]にペイントする
  unsigned long fillSeedMid(int* mid,
                            const int face,
                            const int target,
                            const int* Dsize=NULL);
  
  
  // list[]内の最頻値IDを求める
  int find_mode(const int m_sz,
                const int* list,
                const int m_NoCompo);
  
  
  // サブセル内の最頻値IDを求める
  int find_mode_smd(const int* smd, const int m_NoCompo);
  
  
  // セルに含まれるポリゴンを探索し、d_midに記録
  unsigned long findPolygonInCell(int* d_mid,
                                  MPIPolylib* PL,
                                  PolygonProperty* PG,
                                  const int m_NoCompo);
  
  
  // フィルパラメータを取得
  void getFillParam();
  
  
  // サブセルのペイント
  int SubCellFill(REAL_TYPE* svf,
                  int* smd,
                  const int dir,
                  const int refID,
                  const REAL_TYPE refVf);
  
  
  // サブセルのポリゴン含有テスト
  int SubCellIncTest(REAL_TYPE* svf,
                     int* smd,
                     const int ip,
                     const int jp,
                     const int kp,
                     const Vec3r pch,
                     const string m_pg,
                     MPIPolylib* PL,
                     const int m_NoCompo);
  
  
  // sub-division
  void SubDivision(REAL_TYPE* svf,
                   int* smd,
                   const int ip,
                   const int jp,
                   const int kp,
                   int* d_mid,
                   const MediumList* mat,
                   REAL_TYPE* d_pvf,
                   const int m_NoCompo);
  
  
public:
  
  // ポリゴングループの座標値からboxを計算する
  void calcBboxFromPolygonGroup(MPIPolylib* PL,
                                PolygonProperty* PG,
                                const int m_NoPolyGrp,
                                const REAL_TYPE* poly_org,
                                const REAL_TYPE* poly_dx);
  
  
  // ガイドセルのIDをbcdからmidに転写
  void copyIDonGuide(const int face, const int* bcd, int* mid);
  
  
  // bcd[]内にあるm_idのセルを数える
  unsigned long countCellB(const int* bcd,
                           const int m_id,
                           const bool painted=true,
                           const int* Dsize=NULL);
  
  
  // フィル操作
  void fill(FILE* fp,
            CompoList* cmp,
            MediumList* mat,
            int* d_bcd,
            float* d_cut,
            int* d_bid,
            const int m_NoCompo);

  
  // カットID情報に基づく流体媒質のフィルを実行
  unsigned long fillByBid(int* bid,
                          int* bcd,
                          float* cut,
                          const int tgt_id,
                          unsigned long& substituted,
                          const int m_NoCompo,
                          const int* Dsize=NULL);
  
  
  // 未ペイントセルをFLUIDでフィル
  unsigned long fillByFluid(int* bcd,
                            const int fluid_id,
                            const int* bid,
                            const int* Dsize=NULL);
  
  
  // 交点が定義点にある場合にそのポリゴンのエントリ番号でフィルする
  unsigned long fillCutOnCellCenter(int* bcd,
                                    const int* bid,
                                    const float* cut,
                                    const int* Dsize=NULL);
  
  
  // シード点をbcd[]にペイントする
  unsigned long fillSeedBcd(int* bcd,
                            const int face,
                            const int target,
                            const int* bid,
                            const int* Dsize=NULL);
  
  
  // フィルパラメータを取得
  void getFillParam(TextParser* tpCntl);
  
  
  // FIllIDとSeedIDをセット
  void setFillMedium(MediumList* mat, const int m_NoMedium);
  
  
  // サブサンプリング
  void SubSampling(FILE* fp,
                   MediumList* mat,
                   int* d_mid,
                   REAL_TYPE* d_pvf,
                   MPIPolylib* PL,
                   const int m_NoCompo);
  
  
  // ポリゴンの水密化
  void SeedFilling(FILE* fp,
                   CompoList* cmp,
                   MediumList* mat,
                   int* d_mid,
                   MPIPolylib* PL,
                   PolygonProperty* PG,
                   const int m_NoCompo);
  
  
  // @brief 再分割数を設定
  // @param [in] num 再分割数
  void setSubDivision(int num)
  {
    if ( num < 1 ) Exit(0);
    NumSuvDiv = num;
  }

};

#endif // _FFV_GEOM_H_