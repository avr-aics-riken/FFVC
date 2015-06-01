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
#include "Alloc.h"

#include "Polylib.h"
#include "MPIPolylib.h"
#include <vector>
#include <algorithm>

using namespace std;
using namespace PolylibNS;
using namespace Vec3class;


class Geometry : public DomainInfo {
  
public:
  enum kind_fill_hint {
    kind_outerface=1,
    kind_point
  };
  
  
private:
  int NumSuvDiv;       ///< 再分割数
  int FillSeedDir;     ///< フィルのヒント {x_minux | x_plus |...}
  string SeedMedium;   ///< ヒントに使う媒質 -> int SeedID
  unsigned temporary;
  int NoHint;
  
  typedef struct {
    string identifier;
    kind_fill_hint kind;
    int dir;
    string medium;
    REAL_TYPE point[3];  ///< 有次元座標
  } KindFill;
  
  KindFill* fill_table;
  
  
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
    NoHint=0;
    
    for (int i=0; i<3; i++) {
      FillSuppress[i] = ON; // default is "fill"
    }
    
    fill_table = NULL;
  }
  
  /** デストラクタ */
  ~Geometry() {}
  

private:
  
  /* @brief 重複の無いようにリストにラベルを追加する
   * @param [in,out] tbl  ラベルリスト
   * @param [in]     dd   検査するラベル
   */
  void addLabel2List(vector<int>& tbl, const int dd)
  {
    vector<int>::iterator it;
    it = find(tbl.begin(), tbl.end(), dd);
    
    // labelTable内にddがないとき、末尾に追加
    if ( it == tbl.end() ) tbl.push_back(dd);
  }
  
  
  // 各ノードのラベルのユニーク性を担保する
  void assureUniqueLabel(vector<int> tbl, int* mid, const int* Dsize=NULL);
  
  
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
  
  
  // セル数をカウント（デバッグ用）
  unsigned long debug_countCellB(const int* bcd,
                                 const int m_id,
                                 const int* bid,
                                 const int* Dsize);
  
  
  // 未ペイントセルをtargetでフィル
  unsigned long fillByID(int* mid,
                         const int target,
                         const int* Dsize=NULL);
  
  
  // 流体媒質のフィルをbid情報を元に実行
  unsigned long fillByMid(int* mid,
                          const int tgt_id,
                          const int* Dsize=NULL);
  
  
  // 未ペイントセルを周囲の交点IDの最頻値でフィル
  bool fillByModalCutID(int* bcd,
                        const int* bid,
                        const int m_NoCompo,
                        const CompoList* cmp,
                        unsigned long& replaced);
  
  
  // 未ペイントセルを周囲のbidの固体最頻値でフィル
  unsigned long fillByModalSolid(int* bcd,
                                 const int* bid,
                                 const int m_NoCompo,
                                 const CompoList* cmp);
  
  
  // bid情報を元にフラッドフィル
  bool fill_connected(FILE* fp,
                      int* d_bcd,
                      const int* d_bid,
                      const MediumList* mat,
                      const int m_NoMedium,
                      const int fill_mode,
                      unsigned long& target_count);
  
  
  // サブセルのSolid部分の値を代入
  unsigned long fillSubCellSolid(int* smd, REAL_TYPE* svf);
  
  
  // 流体セルをマイナス値でワーク配列にコピー
  unsigned long fillFluidRegion(int* mid, const int* bcd);
  
  
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
  
  
  // ペイント対象のセルを含むbboxを探す
  unsigned long findBboxforSeeding(int* mid, int* bbox, const int* Dsize=NULL);
  
  
  // 点pの属するセルインデクスを求める
  // @param [in]  pt 有次元座標
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
  
  
  // 固体領域をペイントするシードセルを探す
  int findSeedCells(int* mid,
                    const int* bbox,
                    const int counter,
                    const int* Dsize=NULL);
  
  

  
  
  // 連結領域を同定する
  bool identifyConnectRegion(int* mid,
                             const int* bcd,
                             const int* bid,
                             const unsigned long paintable,
                             const unsigned long filled_fluid,
                             const int* Dsize=NULL);
  
  
  // 各ラベルの接続リストを作成
  int makeConnectList(vector< vector<int> >& cnct,
                      vector<int>& label,
                      const int* mid,
                      const int* bid,
                      const int* Dsize=NULL);
  
  
  // 距離の最小値を求める
  void minDistance(const long long* cut, const int* bid, FILE* fp);
  
  
  // 固体領域をスレッド毎にIDでペイント
  void paintSolidRegion(int* mid, const int* Dsize=NULL);
  
  
  // 接続している領域を指定ラベルでペイント
  void paintConnectedLabel(int* mid,
                           const vector<int>& label,
                           const int replace,
                           const int* Dsize=NULL);
  
  
  // 6方向にカットのあるセルを交点媒質でフィルする
  unsigned long replaceIsolatedCell(int* bcd,
                                    const int* bid,
                                    const int m_NoCompo,
                                    const CompoList* cmp);
  
  
  // 探査しながらペイント
  unsigned long searchPaint(int* mid,
                            const int* bbox,
                            const int* bid,
                            const int* Dsize=NULL);
  
  
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
  bool fill(FILE* fp,
            int* d_bcd,
            const int* d_bid,
            const int m_NoMedium,
            const MediumList* mat,
            const int m_NoCompo,
            const CompoList* cmp,
            int* d_mid);

  
  // bid情報を元にフラッドフィル
  unsigned long fillByBid(int* bcd,
                          const int* bid,
                          const MediumList* mat,
                          const int mode,
                          const int* Dsize=NULL);
  
  
  // 未ペイントセルをFLUIDでフィル
  unsigned long fillByFluid(int* bcd,
                            const int fluid_id,
                            const int* bid,
                            const int* Dsize=NULL);
  
  
  // bcd[]の内部セルにシードIDをペイントする
  unsigned long fillSeedBcdInner(int* bcd,
                                 const REAL_TYPE p[3],
                                 const int target,
                                 const int* Dsize=NULL);
  
  // bcd[]の外層にシードIDをペイントする
  unsigned long fillSeedBcdOuter(int* bcd,
                                 const int face,
                                 const int target,
                                 const int* bid,
                                 const int* Dsize=NULL);
  
  
  // フィルパラメータを取得
  void getFillParam(TextParser* tpCntl,
                    FILE* fp,
                    const int Unit,
                    const REAL_TYPE RefL,
                    const int m_NoMedium,
                    const MediumList* mat);
  
  
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
  void paintCutOnPoint(int* bcd,
                       int* bid,
                       long long* cut,
                       const int m_NoCompo,
                       const CompoList* cmp,
                       unsigned long& fillcut,
                       unsigned long& modopp,
                       const int* Dsize=NULL);
  
  
  // 交点計算を行い、量子化
  void quantizeCut(FILE* fp,
                   long long* cut,
                   int* bid,
                   int* bcd,
                   const int m_NoCompo,
                   const CompoList* cmp,
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