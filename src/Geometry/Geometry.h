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
  int NoMedium;        ///< 媒質数
  int NoCompo;         ///< コンポーネント数
  
  typedef struct {
    string identifier;
    kind_fill_hint kind;
    int dir;
    string medium;
    REAL_TYPE point[3];  ///< 有次元座標
  } KindFill;
  
  KindFill* fill_table;
  
  const MediumList* mat;
  const CompoList* cmp;
  FILE* fpc;            ///< condition.txtへのファイルポインタ
  
  
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
    NoHint = 0;
    NoMedium = 0;
    NoCompo  = 0;
    
    for (int i=0; i<3; i++) {
      FillSuppress[i] = ON; // default is "fill"
    }
    
    fill_table = NULL;
    mat = NULL;
    cmp = NULL;
    fpc = NULL;
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
  void assureUniqueLabel(int* labelTop,
                         int* labelsz,
                         vector<int>& tbl,
                         int* mid,
                         const int* Dsize=NULL);
  
  
  // d_mid[]がtargetであるセルに対して、d_pvf[]に指定値valueを代入する
  unsigned long assignVF(const int target,
                         const REAL_TYPE value,
                         const int* d_mid,
                         REAL_TYPE* d_pvf);
  
  
  // mid[]内にあるm_idのセルを数える
  unsigned long countCellM(const int* mid,
                           const int m_id,
                           const string mode,
                           const bool painted=true,
                           const int* Dsize=NULL);
  
  
  // セル数をカウント（デバッグ用）
  unsigned long debug_countCellB(const int* bcd,
                                 const int m_id,
                                 const int* bid,
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
  bool fillByModalSolid(int* bcd,
                        const int* bid,
                        unsigned long& replaced);
  
  
  // bid情報を元にフラッドフィル
  bool fillConnected(int* d_bcd,
                     const int* d_bid,
                     const int fill_mode,
                     unsigned long& target_count);
  
  
  // 未ペイントセルをpaintIDで連結フィルする
  unsigned long  fillConnected4ID(int* mid,
                                  const int* bid,
                                  const int paintID,
                                  const int* Dsize=NULL);
  
  
  // サブセルのSolid部分の値を代入
  unsigned long fillSubCellSolid(int* smd, REAL_TYPE* svf);
  
  
  // 流体セルをマイナス値でワーク配列にコピー
  unsigned long fillFluidRegion(int* mid, const int* bcd);
  
  
  // サブセルの未ペイントセルを周囲の媒質IDの固体最頻値でフィル
  unsigned long fillSubCellByModalSolid(int* smd, REAL_TYPE* svf);
  
  
  // シード点をmid[]にペイントする
  unsigned long fillSeedMid(int* mid,
                            const int face,
                            const int target,
                            const int* Dsize=NULL);
  
  
  // ペイント対象のセルを含むbboxを探す
  unsigned long findBboxforSeeding(const int* mid, int* bbox, const int* Dsize=NULL);
  
  
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
  
  
  // 未ペイントセルのシードセルをひとつ探す
  bool findKeyCell(int* mid,
                   const int counter,
                   const int* Dsize=NULL);
  
  
  // list[]内の最頻値IDを求める
  int find_mode(const int m_sz, const int* list);
  
  
  // サブセル内の最頻値IDを求める
  int find_mode_smd(const int* smd);
  
  
  // セルに含まれるポリゴンを探索し、d_midに記録
  unsigned long findPolygonInCell(int* d_mid,
                                  MPIPolylib* PL,
                                  PolygonProperty* PG);
  
  
  // 各ランクの最頻値情報を集める
  void gatherModes(vector<unsigned long>& mode);
  
                             
  // 各ランクの接続ルール数の情報を集める
  bool gatherRules(int* buffer,
                   const int maxRule,
                   const vector< vector<int> > cnct);
  
  
  // ラベルを集める
  bool gatherLabels(int* buffer,
                    const int maxSize,
                    const vector< vector<int> > cnct);
  
  
  // ラベル情報集約のため、バッファサイズを計算する
  int getNumBuffer(const vector< vector<int> > cnct);
  
  
  // 連結領域を同定する
  bool identifyConnectedRegion(int* mid,
                               int* bcd,
                               const int* bid,
                               const unsigned long paintable,
                               const unsigned long filled_fluid);

  
  // 各ラベルの接続リストを作成
  void makeConnectList(vector< vector<int> >& cnct,
                       const int* ltop,
                       const vector<int> localTbl,
                       const int* mid,
                       const int* bid,
                       const int* Dsize=NULL);
  
  
  // ラベルとmat[]インデクスの対応を決める
  bool makeHistgramByModalCut(vector<int>& matIndex,
                              const vector<int> labels,
                              const int* mid,
                              const int* bid,
                              const int* Dsize=NULL);
  
  
  // 全接続ルールを作成
  bool makeWholeRules(vector< vector<int> >& connectRules,
                      const int* labelSize,
                      const int* pckdLabel,
                      const int width_label,
                      const int* pckdRules,
                      const int width_rule);
  
  
  // 同種のセルでカットをもつ境界の修正
  void mergeSameSolidCell(int* bid,
                          long long* cut,
                          const int* bcd,
                          const int* Dsize=NULL);
  
  
  // 距離の最小値を求める
  void minDistance(const long long* cut, const int* bid);
  
  
  // 固体領域をスレッド毎にIDでペイント
  void paintSolidRegion(int* mid, const int* Dsize=NULL);
  
  
  // 対応づけしたmat[]インデクスによりラベル部分をペイントする
  void paintCellByUniqueLabel(int* bcd,
                              const int* mid,
                              const vector<int> labels,
                              const vector<int> matIndex,
                              const int* Dsize=NULL);
  
  
  // 接続している領域を指定ラベルでペイント
  void paintConnectedLabel(int* mid,
                           const vector<int> labels,
                           const vector< vector<int> > rules,
                           const int* Dsize=NULL);
  
  
  // ラベルを縮約する
  void reduceLabels(vector<int>& finalLabels,
                    vector< vector<int> >& finalRules,
                    vector< vector<int> > connectRules);
  
  
  // ラベルを登録しながら、ルールを縮約する
  void registerRule(vector< vector<int> >& finalRules,
                    const int reg_keys,
                    const int test_key,
                    vector< vector<int> >& connectRules,
                    vector<int>& valid);
  
  
  // 6方向にカットのあるセルを交点媒質でフィルする
  unsigned long replaceIsolatedCell(int* bcd, const int* bid);
  
  
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
                     MPIPolylib* PL);
  
  
  // sub-division
  void SubDivision(REAL_TYPE* svf,
                   int* smd,
                   const int ip,
                   const int jp,
                   const int kp,
                   int* d_mid,
                   REAL_TYPE* d_pvf);
  
  
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
  bool fill(int* d_bcd,
            int* d_bid,
            int* d_mid,
            long long* d_cut);

  
  // bid情報を元にフラッドフィル
  unsigned long fillByBid(int* bcd,
                          const int* bid,
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
                    const int Unit,
                    const REAL_TYPE RefL,
                    const int m_NoMedium,
                    const MediumList* m_mat,
                    FILE* m_fp);
  
  
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
                       unsigned long& fillcut,
                       unsigned long& modopp,
                       const int* Dsize=NULL);
  
  
  // 交点計算を行い、量子化
  void quantizeCut(long long* cut,
                   int* bid,
                   int* bcd,
                   MPIPolylib* PL,
                   PolygonProperty* PG,
                   const int* Dsize=NULL);
  
  
  // ポリゴンの水密化
  void SeedFilling(int* d_mid,
                   MPIPolylib* PL,
                   PolygonProperty* PG);
  
  
  // コンポーネント数をセット
  void setCompoPtr(const int m_NoCompo, const CompoList* m_cmp)
  {
    NoCompo = m_NoCompo;
    cmp = m_cmp;
  }
  
  
  // @brief 再分割数を設定
  // @param [in] num 再分割数
  void setSubDivision(int num)
  {
    if ( num < 1 ) Exit(0);
    NumSuvDiv = num;
  }
  
  
  // サブサンプリング
  void SubSampling(int* d_mid,
                   REAL_TYPE* d_pvf,
                   MPIPolylib* PL);

};

#endif // _FFV_GEOM_H_
