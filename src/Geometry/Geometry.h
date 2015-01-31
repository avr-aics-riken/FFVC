#ifndef _FFV_GEOM_H_
#define _FFV_GEOM_H_
//##################################################################################
//
// FFV-C : Frontflow / violet Cartesian
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
 * @file   Geometry.h
 * @brief  Geometry Class header
 * @author aics
 */

#include "DomainInfo.h"
#include "FB_Define.h"
#include "Medium.h"
#include "PolyProperty.h"
#include "Component.h"
#include "TextParser.h"

#include "Polylib.h"
#include "MPIPolylib.h"

using namespace std;
using namespace PolylibNS;
using namespace Vec3class;

class Geometry : public DomainInfo {
  
private:
  int NumSuvDiv;       ///< 再分割数
  int FillSeedDir;     ///< フィルのヒント {x_minux | x_plus |...}
  string FillMedium;   ///< フィルに使う媒質 -> int FillID
  string SeedMedium;   ///< ヒントに使う媒質 -> int SeedID
  unsigned temporary;
  
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
    temporary = 0;
    
    for (int i=0; i<3; i++) {
      FillSuppress[i] = ON; // default is "fill"
    }
  }
  
  /**　デストラクタ */
  ~Geometry() {}
  

private:
  
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
  
  
  // 点pの属するセルインデクスを求める
  // @param [in]  pt 無次元座標
  // @param [out] w  インデクス
  inline void findIndex(const Vec3r pt, int* w) const
  {
    REAL_TYPE p[3], q[3];
    p[0] = (REAL_TYPE)pt.x;
    p[1] = (REAL_TYPE)pt.y;
    p[2] = (REAL_TYPE)pt.z;
    
    q[0] = (p[0]-originD[0])/pitchD[0];
    q[1] = (p[1]-originD[1])/pitchD[1];
    q[2] = (p[2]-originD[2])/pitchD[2];
    
    w[0] = (int)ceil(q[0]);
    w[1] = (int)ceil(q[1]);
    w[2] = (int)ceil(q[2]);
  }
  
  
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
  
  
  // 交点が定義点にある場合の処理をした場合に、反対側のセルの状態を修正
  unsigned long modifyCutOnPoint(int* bid, long long* cut, const int* bcd, const int* Dsize=NULL);
  
  
  /**
   * @brief インデックスを(1,0,0)シフト
   * @param [in] index 元のインデクス
   * @param [in] h     シフト幅
   */
  inline Vec3r shift_E(const Vec3r index, const Vec3r h)
  {
    return Vec3r(index.x+h.x, index.y, index.z);
  }
  
  
  /**
   * @brief インデックスを(-1,0,0)シフト
   * @param [in] index 元のインデクス
   * @param [in] h     シフト幅
   */
  inline Vec3r shift_W(const Vec3r index, const Vec3r h)
  {
    return Vec3r(index.x-h.x, index.y, index.z  );
  }
  
  
  /**
   * @brief インデックスを(0,1,0)シフト
   * @param [in] index 元のインデクス
   * @param [in] h     シフト幅
   */
  inline Vec3r shift_N(const Vec3r index, const Vec3r h)
  {
    return Vec3r(index.x, index.y+h.y, index.z);
  }
  
  
  /**
   * @brief インデックスを(0,-1,0)シフト
   * @param [in] index 元のインデクス
   * @param [in] h     シフト幅
   */
  inline Vec3r shift_S(const Vec3r index, const Vec3r h)
  {
    return Vec3r(index.x, index.y-h.y, index.z);
  }
  
  
  /**
   * @brief インデックスを(0,0,1)シフト
   * @param [in] index 元のインデクス
   * @param [in] h     シフト幅
   */
  inline Vec3r shift_T(const Vec3r index, const Vec3r h)
  {
    return Vec3r(index.x, index.y, index.z+h.z);
  }
  
  
  /**
   * @brief インデックスを(0,0,-1)シフト
   * @param [in] index 元のインデクス
   * @param [in] h     シフト幅
   */
  inline Vec3r shift_B(const Vec3r index, const Vec3r h)
  {
    return Vec3r(index.x, index.y, index.z-h.z);
  }
  
  
  
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
  
  
  // ポリゴンと線分の交点計算
  bool TriangleIntersect(const Vec3r ray_o,
                         const Vec3r dir,
                         const Vec3r v0,
                         const Vec3r v1,
                         const Vec3r v2,
                         REAL_TYPE& pRetT,
                         REAL_TYPE& pRetU,
                         REAL_TYPE& pRetV);
  
  
  // 交点情報をアップデート
  unsigned updateCut(const Vec3r ray_o,
                     const int dir,
                     const Vec3r v0,
                     const Vec3r v1,
                     const Vec3r v2,
                     long long& cut,
                     int& bid,
                     const int pid);
  
  
  
  
public:
  
  // ポリゴングループの座標値からboxを計算する
  void calcBboxFromPolygonGroup(MPIPolylib* PL,
                                PolygonProperty* PG,
                                const int m_NoPolyGrp);
  
  
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
            long long* d_cut,
            int* d_bid,
            const int m_NoCompo);

  
  // カットID情報に基づく流体媒質のフィルを実行
  unsigned long fillByBid(int* bid,
                          int* bcd,
                          long long* cut,
                          const int tgt_id,
                          unsigned long& substituted,
                          const int m_NoCompo,
                          const int* Dsize=NULL);
  
  
  // 未ペイントセルをFLUIDでフィル
  unsigned long fillByFluid(int* bcd,
                            const int fluid_id,
                            const int* bid,
                            const int* Dsize=NULL);
  
  
  // シード点をbcd[]にペイントする
  unsigned long fillSeedBcd(int* bcd,
                            const int face,
                            const int target,
                            const int* bid,
                            const int* Dsize=NULL);
  
  
  // フィルパラメータを取得
  void getFillParam(TextParser* tpCntl);
  
  
  /**
   * @brief ベクトルの最小成分
   * @param [in,out] mn 比較して小さい成分
   * @param [in]     p  参照ベクトル
   */
  static inline void get_min(Vec3r& mn, const Vec3r p)
  {
    mn.x = (mn.x < p.x) ? mn.x : p.x;
    mn.y = (mn.y < p.y) ? mn.y : p.y;
    mn.z = (mn.z < p.z) ? mn.z : p.z;
  }
  
  
  /**
   * @brief ベクトルの最大値成分
   * @param [in,out] mx 比較して大きい成分
   * @param [in]     p  参照ベクトル
   */
  static inline void get_max(Vec3r& mx, const Vec3r p)
  {
    mx.x = (mx.x > p.x) ? mx.x : p.x;
    mx.y = (mx.y > p.y) ? mx.y : p.y;
    mx.z = (mx.z > p.z) ? mx.z : p.z;
  }
  
  
  // 交点が定義点にある場合にそのポリゴンのエントリ番号でフィルする
  unsigned long paintCutOnPoint(int* bcd,
                                int* bid,
                                long long* cut,
                                const int m_NoCompo,
                                const int* Dsize=NULL);
  
  
  // 交点計算を行い、量子化
  void quantizeCut(FILE* fp,
                   long long* cut,
                   int* bid,
                   int* bcd,
                   const int m_NoCompo,
                   MPIPolylib* PL,
                   PolygonProperty* PG);
  
  
  // ポリゴンの水密化
  void SeedFilling(FILE* fp,
                   CompoList* cmp,
                   MediumList* mat,
                   int* d_mid,
                   MPIPolylib* PL,
                   PolygonProperty* PG,
                   const int m_NoCompo);
  
  
  // FIllIDとSeedIDをセット
  void setFillMedium(MediumList* mat, const int m_NoMedium);
  
  
  // @brief 再分割数を設定
  // @param [in] num 再分割数
  void setSubDivision(int num)
  {
    if ( num < 1 ) Exit(0);
    NumSuvDiv = num;
  }
  
  
  // サブサンプリング
  void SubSampling(FILE* fp,
                   MediumList* mat,
                   int* d_mid,
                   REAL_TYPE* d_pvf,
                   MPIPolylib* PL,
                   const int m_NoCompo);

};

#endif // _FFV_GEOM_H_