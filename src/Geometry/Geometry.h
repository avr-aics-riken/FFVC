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
// Copyright (c) 2012-2016 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
// Copyright (c) 2016-2020 Research Institute for Information Research, Kyushu University.
// All rights reserved.
//
//##################################################################################

/*
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
#ifndef DISABLE_MPI
  #include "MPIPolylib.h"
#endif

#include <vector>

#define SEARCH_WIDTH 6
#define DIST_AB 2.0

using namespace std;
using namespace PolylibNS;
using namespace Vec3class;


// 交点の法線方向探索結果を整列するためのテンポラリデータを格納

struct cut_probe {
  REAL_TYPE dist; // 距離
  Vec3r refp;     // 参照点
  Vec3r cp;       // 交点
  Vec3r nv;       // 交点からの法線
};
  

// 比較関数
static bool cmp_dist(const cut_probe& a, const cut_probe& b) {
  return a.dist < b.dist;
}
static bool cmp_i(const cut_probe& a, const cut_probe& b) {
  return a.refp.x < b.refp.x;
}
static bool cmp_j(const cut_probe& a, const cut_probe& b) {
  return a.refp.y < b.refp.y;
}
static bool cmp_k(const cut_probe& a, const cut_probe& b) {
  return a.refp.z < b.refp.z;
}




class Geometry : public DomainInfo {

private:
  int NumSuvDiv;       ///< 再分割数
  int FillSeedDir;     ///< フィルのヒント {x_minux | x_plus |...}
  string SeedMedium;   ///< ヒントに使う媒質 -> int SeedID
  unsigned temporary;
  int NoHint;
  int NoMedium;        ///< 媒質数
  int NoCompo;         ///< コンポーネント数
  
  Vec3r my_reg_min;    ///<  自領域の最小値　参照点チェック用
  Vec3r my_reg_max;    ///<  自領域の最大値

  typedef struct {
    string identifier;
    int dir;
    string medium;
  } KindFill;

  KindFill* fill_table;

  const MediumList* mat;
  const CompoList* cmp;
  FILE* fpc;            ///< condition.txtへのファイルポインタ
  
  // 各ランクごとのデータ保持
  vector<cut_probe> prb;    ///< 交点の参照座標と距離
  
  // debug
  vector<cut_probe> prb_dbg;


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

  
  // mid[]内にあるm_idのセルを数える
  unsigned long countCellM(const int* mid,
                           const int m_id,
                           const string mode,
                           const bool painted=true,
                           const int* Dsize=NULL);

  
  // bid情報を元にフラッドフィル
  unsigned long fillbyBid(int* bcd,
                          const int* bid,
                          const int mode,
                          const int* Dsize=NULL);
  
  
  // 未ペイントセルをFLUIDでフィル
  unsigned long fillByFluid(int* bcd,
                            const int fluid_id,
                            const int* bid,
                            const int* Dsize=NULL);
  

  // bid情報を元にフラッドフィル
  bool fillConnected(int* d_bcd,
                     const int* d_bid,
                     const int fill_mode,
                     unsigned long& target_count);


  // 流体セルをマイナス値でワーク配列にコピー
  unsigned long fillFluidRegion(int* mid, const int* bcd);


  // bcd[]の外層にシードIDをペイントする
  unsigned long fillSeedBcdOuter(int* bcd,
                                 const int face,
                                 const int target,
                                 const int* bid,
                                 const int* Dsize=NULL);

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
  unsigned mergeSameSolidCell(int* bid,
                              int* cutL,
                              int* cutU,
                              const int* bcd,
                              const int* Dsize=NULL);


  // 距離の最小値を求める
  void minDistance(const int* cutL, const int* cutU, const int* bid);


  


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


  // 6方向にカットのあるセルをSOLIDに変更
  void paintIsolatedCell(int* bid, int*cutL, int* cutU);
  
  
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
                     int& cutL,
                     int& cutU,
                     int& bid,
                     const int pid);

  // 未ペイントセルをpaintIDで連結フィルする
  unsigned long  fillConnected4ID(int* mid,
                                  const int* bid,
                                  const int paintID,
                                  const int* Dsize=NULL);
  
  

  // @brief 参照点の計算領域内のチェック
  // @param [in]   p   参照点の計算空間座標
  inline bool check_region(const Vec3r p)
  {
    if ( p.x < my_reg_min.x || p.x > my_reg_max.x
      || p.y < my_reg_min.y || p.y > my_reg_max.y
      || p.z < my_reg_min.z || p.z > my_reg_max.z ) return false;

    return true;
  }
  
  // @brief インデクスの範囲チェック
  inline bool check_index_region(const Vec3r t) {
    int i = (int)t.x;
    int j = (int)t.y;
    int k = (int)t.z;
    if ( i<0 || i>size[0]+1 || j<0 || j>size[1]+1 || k<0 || k>size[2]+1 ) return false;
    return true;
  }
  
  
  // SDF
  
  void estimateNV(int* nrm, const REAL_TYPE* fo);
  
  void smoothingNV(int* nrm, REAL_TYPE* fv);
  
  void tracingSDF(REAL_TYPE* fn, int* nrm, const int* wk, const int nLayer);
  
  void setLayerInit(REAL_TYPE* fn, REAL_TYPE* fo, const int* nrm, const int* wk);
  
  void smoothingScalar(int* nrm, REAL_TYPE* fo, const string str);
  
  int generateLayer(int* wk, const int* nrm);
  
  bool getDistance(const Vec3r o,
                   const Vec3r n,
                   const int* nrm,
                   const REAL_TYPE* f,
                   REAL_TYPE& ds,
                   REAL_TYPE& s,
                   Vec3r& p,
                   Vec3r& v);
  
  inline REAL_TYPE getInterp(const Vec3r g, const REAL_TYPE* f)
  {
    return (1.0-g.x)*(1.0-g.y)*(1.0-g.z)* f[0]
         +      g.x *(1.0-g.y)*(1.0-g.z)* f[1]
         + (1.0-g.x)*     g.y *(1.0-g.z)* f[2]
         +      g.x *     g.y *(1.0-g.z)* f[3]
         + (1.0-g.x)*(1.0-g.y)*     g.z * f[4]
         +      g.x *(1.0-g.y)*     g.z * f[5]
         + (1.0-g.x)*     g.y *     g.z * f[6]
         +      g.x *     g.y *     g.z * f[7];
  }
  
  inline Vec3r getInterp(const Vec3r g, const Vec3r f[8])
  {
    Vec3r r;
    r.x = (1.0-g.x)*(1.0-g.y)*(1.0-g.z)* f[0].x
        +      g.x *(1.0-g.y)*(1.0-g.z)* f[1].x
        + (1.0-g.x)*     g.y *(1.0-g.z)* f[2].x
        +      g.x *     g.y *(1.0-g.z)* f[3].x
        + (1.0-g.x)*(1.0-g.y)*     g.z * f[4].x
        +      g.x *(1.0-g.y)*     g.z * f[5].x
        + (1.0-g.x)*     g.y *     g.z * f[6].x
        +      g.x *     g.y *     g.z * f[7].x;
    r.y = (1.0-g.x)*(1.0-g.y)*(1.0-g.z)* f[0].y
        +      g.x *(1.0-g.y)*(1.0-g.z)* f[1].y
        + (1.0-g.x)*     g.y *(1.0-g.z)* f[2].y
        +      g.x *     g.y *(1.0-g.z)* f[3].y
        + (1.0-g.x)*(1.0-g.y)*     g.z * f[4].y
        +      g.x *(1.0-g.y)*     g.z * f[5].y
        + (1.0-g.x)*     g.y *     g.z * f[6].y
        +      g.x *     g.y *     g.z * f[7].y;
    r.z = (1.0-g.x)*(1.0-g.y)*(1.0-g.z)* f[0].z
        +      g.x *(1.0-g.y)*(1.0-g.z)* f[1].z
        + (1.0-g.x)*     g.y *(1.0-g.z)* f[2].z
        +      g.x *     g.y *(1.0-g.z)* f[3].z
        + (1.0-g.x)*(1.0-g.y)*     g.z * f[4].z
        +      g.x *(1.0-g.y)*     g.z * f[5].z
        + (1.0-g.x)*     g.y *     g.z * f[6].z
        +      g.x *     g.y *     g.z * f[7].z;
    return r;
  }
  
  
  // 参照点補間
  void getBboxCutDir(const int i,
                     const int j,
                     const int k,
                     const int dir,
                     const REAL_TYPE dh,
                     Vec3r& bx_min,
                     Vec3r& bx_max);
  
  bool getNV4CalculatedCut(vector<Triangle*>* trias,
                           const Vec3r ray_o,
                           const int dir,
                           const Vec3r d,
                           const REAL_TYPE dh,
                           const int rc,
                           Vec3r& nv);
  
  bool getMinRefPoint(const Vec3r o,
                      const int dir,
                      const Vec3r d,
                      const Vec3r nv,
                      const REAL_TYPE r,
                      const int* cutL,
                      const int* cutU,
                      cut_probe& cpr,
                      cut_probe& nonref);
  
  bool divideLineByPlane(const Vec3r o,
                         const Vec3r d,
                         const Vec3r a,
                         const Vec3r b,
                         REAL_TYPE& r);
  
  bool decideInterp(const Vec3r o,
                    const Vec3r aa,
                    const Vec3r bb,
                    const Vec3r nv,
                    const int face,
                    const int* cutL,
                    const int* cutU,
                    cut_probe& cpr);
  
  
  inline Vec3r floorVec(Vec3r p) {
    Vec3r q( (int)p.x, (int)p.y, (int)p.z );
    return q;
  }
  

  
  bool judgeInterp(const Vec3r p, const int dir, const int* cutL, const int* cutU);
  
  /* @brief nrmから指定方向の量子化値をとりだす
   * @param [in] n    nrm
   * @note 符号はマイナスの場合，符号ビットが立つ
   */
  inline Vec3r getNrm9 (const int n)
  {
    Vec3r v;
    
    v.x = (REAL_TYPE)getQuantized9(  n        & MASK_9);
    v.y = (REAL_TYPE)getQuantized9( (n >> 10) & MASK_9);
    v.z = (REAL_TYPE)getQuantized9( (n >> 20) & MASK_9);
    if (TEST_BIT(n, 9)) v.x = -v.x;
    if (TEST_BIT(n,19)) v.y = -v.y;
    if (TEST_BIT(n,29)) v.z = -v.z;
    
    return v;
  }


  /*
   * @brief 9bit幅の値の設定
   * @param [in,out] b   int 変数
   * @param [in]     v   法線の値
   */
  inline void setNrm9 (int& b, const Vec3r v)
  {
    b &= (~(MASK_30) ); // 対象30bitをゼロにする
    b |= ( quantize9( fabs(v.x) ) <<  0);          // 書き込む
    b |= ( quantize9( fabs(v.y) ) << 10);
    b |= ( quantize9( fabs(v.z) ) << 20);
    
    // 符号 マイナスのときに該当ビットを立てる
    if (v.x < 0.0)  b |= (0x1<<9);
    if (v.y < 0.0)  b |= (0x1<<19);
    if (v.z < 0.0)  b |= (0x1<<29);
  }
  
  
  /*
   * @brief 法線のコピー
   * @param [out] dst   コピー先
   * @param [in]  src   コピー元
   */
  inline void cpyNrm9 (int& dst, const int src)
  {
    dst &= (~(MASK_30) ); // 対象30bitをゼロにする
    int s = src & MASK_30;
    dst |= s;          // 書き込む
  }
  


  
  
  
// #################################################################################################
public:
  
  static void getDirVec(Vec3r& d, const int dir);
  
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
  
  
  
  
  // 参照点チェック用に自領域の大きさを設定
  void set_my_region()
  {
    // デフォルトサイズは計算領域内のランクを想定し、袖領域の1層外側まで
    my_reg_min.assign(0.0, 0.0, 0.0);
    my_reg_max.assign((REAL_TYPE)size[0] + 1.0,
                      (REAL_TYPE)size[1] + 1.0,
                      (REAL_TYPE)size[2] + 1.0);
    
    // 外部境界面は計算領域内
    if ( nID[X_minus] < 0 ) my_reg_min.x += 1.0;
    if ( nID[X_plus]  < 0 ) my_reg_max.x -= 1.0;
    if ( nID[Y_minus] < 0 ) my_reg_min.y += 1.0;
    if ( nID[Y_plus]  < 0 ) my_reg_max.y -= 1.0;
    if ( nID[Z_minus] < 0 ) my_reg_min.z += 1.0;
    if ( nID[Z_plus]  < 0 ) my_reg_max.z -= 1.0;
    
    //printf("%f %f %f : %f %f %f\n", my_reg_min.x, my_reg_min.y, my_reg_min.z, my_reg_max.x, my_reg_max.y, my_reg_max.z);
  }
  
    // ポリゴングループの座標値からboxを計算する
  #ifndef DISABLE_MPI
    void calcBboxFromPolygonGroup(MPIPolylib* PL,
                                  PolygonProperty* PG,
                                  const int m_NoPolyGrp);
  #else
    void calcBboxFromPolygonGroup(Polylib* PL,
                                  PolygonProperty* PG,
                                  const int m_NoPolyGrp);
  #endif
  

  // bcd[]内にあるm_idのセルを数える
  unsigned long countCellB(const int* bcd,
                           const int m_id,
                           const bool painted=true,
                           const int* Dsize=NULL);
  
  // フィル操作
  bool fillbyCut(int* d_bcd,
                 int* d_bid,
                 int* d_mid,
                 int* d_cutL,
                 int* d_cutU);


  // フィルパラメータを取得
  void getFillParam(TextParser* tpCntl,
                    const int Unit,
                    const REAL_TYPE m_RefL,
                    const int m_NoMedium,
                    const MediumList* m_mat,
                    FILE* m_fp);


  
    // 交点計算を行い、量子化
  #ifndef DISABLE_MPI
    void quantizeCut(int* cutL,
                     int* cutU,
                     int* bid,
                     MPIPolylib* PL,
                     PolygonProperty* PG,
                     const int* Dsize=NULL);
  #else
    void quantizeCut(int* cutL,
                     int* cutU,
                     int* bid,
                     Polylib* PL,
                     PolygonProperty* PG,
                     const int* Dsize=NULL);
  #endif
    



    // コンポーネント数をセット
  void setCompoPtr(const int m_NoCompo, const CompoList* m_cmp)
  {
    NoCompo = m_NoCompo;
    cmp = m_cmp;
  }



  
  // SDF
  void polygon2sdf(REAL_TYPE* sdf,
                   MPIPolylib* PL,
                   int* nrm,
                   const int* Dsize=NULL);
  
  void outerSDF(REAL_TYPE* sdf, int* nrm, const int* oflag);
  
  void trimNarrowBand(REAL_TYPE* sdf, int* nrm, const REAL_TYPE delta);
  
  bool marchingSDF(REAL_TYPE* fn, REAL_TYPE* fo, int* nrm, int* wk, REAL_TYPE* fv);
  
  void getNVfromIdx(const int* nrm, REAL_TYPE* f);
  
  
  
  // 参照点補間
  int getRefPointOnCut(const int* cutL,
                       const int* cutU,
                       MPIPolylib* PL,
                       int* bx,
                       bool inner=false);
  
  
  void writeProbe();

  void sortProbe(REAL_TYPE* cf, REAL_TYPE* ds);
  
};
#endif // _FFV_GEOM_H_
