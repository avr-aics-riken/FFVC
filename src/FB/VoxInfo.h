#ifndef _FB_BINVOX_H_
#define _FB_BINVOX_H_

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
 * @file   VoxInfo.h
 * @brief  FlowBase VoxInfo class Header
 * @author aics
 */

#include "cpm_ParaManager.h"
#include <math.h>
#include <cstdlib>
#include "DomainInfo.h"
#include "FBUtility.h"
#include "Component.h"
#include "Medium.h"
#include "SetBC.h"
#include "BndOuter.h"
#include "Vec3.h"
#include "Intrinsic.h"
#include "limits.h"
#include "omp.h"

using namespace Vec3class;

class VoxInfo : public DomainInfo {
  
private:
  int NoCompo;       ///< コンポーネントの総数
  int NoMedium;      ///< 媒質数
  Intrinsic *Ex;     ///< 例題クラスのポインタ

public:
  /** コンストラクタ */
  VoxInfo() {
    NoCompo = 0;
    NoMedium = 0;
    Ex = NULL;
  }
  
  /**　デストラクタ */
  ~VoxInfo() {}
  
private:
  
  // 外部境界に接するガイドセルのbcd[]にIDを内部周期境界からコピーする
  void copyIdPrdcInner (int* bcd, const int* m_st, const int* m_ed, const int m_id, const int m_dir);
  
  
  // Naive実装に使う係数を保持
  void copyCoefNaive(int* bx, REAL_TYPE* pn);
  
  
  // セルセンターのIDから対象セル数をカウントし，サブドメイン内にコンポーネントがあれば存在フラグを立てる
  unsigned long countCC (const int order, const int* bx, CompoList* cmp);
  
  
  // セルフェイスの交点IDから対象IDのセル数をカウントし，サブドメイン内にコンポーネントがあれば存在フラグを立てる
  unsigned long countCF (const int key, const int* bx, const int* bid, CompoList* cmp, const string attrb);
  
  
  // 外部境界面の有効セル数をカウントする
  unsigned long countValidCellOBC (const int face, const int* cdf, const int typ);

  
  // BCindexにそのセルが計算に有効(active)かどうかをエンコードする
  void encActive (unsigned long& Lcell, unsigned long& Gcell, int* bx, const int KOS);
  
  
  // 断熱マスクのエンコード
  void encAdiabatic (int* bd, const string target, const int* bid, int face=-1);
  
  
  // セルの各面を調べ，境界条件が設定されていれば，ビットをON
  void encHbit (const int* cdf, int* bd);
  
  
  // ディリクレ条件とノイマン条件の排他性をチェックし，反復行列の非対角要素/対角要素の係数をエンコードする
  void encPbit (int* bx);
  
  
  // 圧力のノイマン境界ビットをエンコードする（カット）
  unsigned long encPbitN (int* bx, const int* bid, const float* cut, const bool convergence);
  
  
  // 計算領域内部のコンポーネントの圧力境界条件フラグをbcp[]にエンコードする
  unsigned long encPbitIBC (const int order,
                            int* bcd,
                            int* bcp,
                            const int* bid,
                            const REAL_TYPE* vec,
                            const string condition,
                            const int bc_dir);
  
  
  // 外部境界に接するセルにおいて，bx[]に圧力境界条件keyに対応するビットフラグを設定する
  void encPbitOBC (int face, int* bx, string key, bool dir);
  
  
  // 熱境界条件のBCエントリをエンコードする
  unsigned long encQface (const int order,
                          const int* bid,
                          int* cdf,
                          int* bd,
                          const bool flag,
                          const int target,
                          CompoList* cmp);
  
  
  // cdf[]にVBCの境界条件に必要な情報をエンコードし，流入流出の場合にbp[]の隣接壁の方向フラグを除く
  // 境界条件指定キーセルのSTATEを流体に変更する
  unsigned long encVbitIBC (const int order,
                            int* cdf,
                            int* bp,
                            const int* cut_id,
                            const REAL_TYPE* vec,
                            const int bc_dir,
                            CompoList* cmp);
  
  
  // cdf[]に境界条件のビット情報をエンコードする
  void encVbitOBC (int face, int* cdf, string key, const bool enc_sw, string chk, int* bp, bool enc_uwd=false);
  
  
  /**
   * @brief 隣接6方向の最頻値IDを求める（fid以外）
   * @param [in] fid 流体のID
   * @param [in] qw  w方向のID
   * @param [in] qe  e方向のID
   * @param [in] qs  s方向のID
   * @param [in] qn  n方向のID
   * @param [in] qb  b方向のID
   * @param [in] qt  t方向のID
   * @note 候補なしの場合には、0が戻り値
   */
  inline int find_mode_id (const int fid, const int qw, const int qe, const int qs, const int qn, const int qb, const int qt)
  {
    unsigned key[CMP_BIT_W]; ///< ID毎の頻度 @note ffv_Initialize() >> fill()でif ( C.NoCompo+1 > CMP_BIT_W )をチェック
    int val[6];              ///< ID
    
    memset(key, 0, sizeof(unsigned)*CMP_BIT_W);
    
    // 0<= q <= Ncompo のはず
    val[0] = qw;
    val[1] = qe;
    val[2] = qs;
    val[3] = qn;
    val[4] = qb;
    val[5] = qt;
    
    
    for (int l=0; l<6; l++) key[ val[l] ]++;
    
    
    int mode = 0;
    int z = 0;
    
    for (int l=NoCompo; l>=1; l--)
    {
      if ( (l != fid) && (key[l] > mode) ) // 流体IDでなく、最大の頻度
      {
        mode = key[l];
        z = l;
      }
    }
    
    return z;
  }
  
  
  // MediumTable内の文字列格納番号を返す
  int getMediumOrder(const MediumList* mat, const std::string key);
  
  
  //@brief idxの第shiftビットをOFFにする
  inline int offBit (int idx, const int shift)
  {
    return ( idx & (~(0x1<<shift)) );
  }
  
  
  //@brief idxの第shiftビットをONにする
  inline int onBit (int idx, const int shift)
  {
    return ( idx | (0x1<<shift) );
  }
  
  
  /*
   * @brief 5bit幅の値の設定
   * @param [in,out] b   int 変数
   * @param [in]     q   5-bit幅のID (1-31)
   * @param [in]     sht シフト量
   */
  inline void setBit5raw (int& b, const int q, const int sht)
  {
    b &= (~(0x1f << sht) ); // 対象5bitをゼロにする
    b |= (q << sht);        // 書き込む
  }
  

  
public:
  
  /**
   * @brief 外部境界に接するガイドセルのbcd[]にIDをエンコードする（内部周期境界の場合）
   * @param [in,out] bcd  BCindex B
   * @param [in]     cmp  CompoList
   */
  void adjMediumPrdcInner (int* bcd, CompoList* cmp);
  
  
  /**
   * @brief BCindex Bのエンコードを確認
   * @param [in] bcd BCindex B
   */
  void chkOrder (const int* bcd);
  
  
  /**
   * @brief dst[]にsrc[]のstate, activeビットの情報をコピーする
   * @param [out] dst コピー先
   * @param [in]  src コピー元
   */
  void copyBCIbase (int* dst, const int* src);
  
  
  /**
   * @brief セル数をカウント
   * @param [in] bcd     BCindex B
   * @param [in] m_id    検査するID
   * @param [in] painted m_id以外のセル数をカウント(false), m_idのセル数をカウント(true)
   * @param [in] Dsize   サイズ
   */
  unsigned long countCellB(const int* bcd, const int m_id, const bool painted=true, const int* Dsize=NULL);
  
  
  /**
   * @brief セル数をカウント
   * @param [in] mid     work array
   * @param [in] m_id    検査するID
   * @param [in] painted m_id以外のセル数をカウント(false), m_idのセル数をカウント(true)
   * @param [in] Dsize   サイズ
   */
  unsigned long countCellM(const int* mid, const int m_id, const bool painted=true, const int* Dsize=NULL);
  
  
  /**
   * @brief セルの状態をカウントして，その個数をLcell, Gcellに保持する
   * @param [out] Lcell ノードローカルの値（戻り値）
   * @param [out] Gcell グローバルの値（戻り値）
   * @param [in]  bx    BCindex B
   * @param [in]  state カウントするセルの状態
   */
  void countCellState (unsigned long& Lcell, unsigned long& Gcell, const int* bx, const int state);
  
  
  /**
   @brief  計算領域の外部境界で外側1層と内側の両方が流体セル数の場合にカウントする
   @param [in]  bx       BCindex B
   @param [out] OpenArea 開口セル数
   */
  void countOpenAreaOfDomain (const int* bx, REAL_TYPE* OpenArea);
  
  
  /**
   * @brief BCindex を表示する（デバッグ用）
   * @param [in] bcd   BCindex B
   * @param [in] cdf   BCindex C
   * @param [in] bcp   BCindex P
   * @param [in] fname 出力用ファイル名
   */
  void dbg_chkBCIndex (const int* bcd, const int* cdf, const int* bcp, const char* fname);
  
  
  /**
   * @brief 流体媒質のフィルをbid情報を元に実行
   * @param [in,out] bid         境界ID（5ビット幅x6方向）
   * @param [in,out] bcd         BCindex B
   * @param [in,out] cut         カット情報
   * @param [in]     tgt_id      フィルする流体IDのエントリ
   * @param [in]     suppress    各軸方向のフィル抑止モード（Periodic, Symmetric時の対策）
   * @param [out]    substituted 固体IDに置換された数
   * @param [in]     Dsize       サイズ
   */
  unsigned long fillByBid (int* bid, int* bcd, float* cut, const int tgt_id, const int* suppress, unsigned long& substituted, const int* Dsize=NULL);
  
  
  /* @brief 未ペイントセルをフィル
   * @param [in,out] bcd      BCindex B
   * @param [in]     fluid_id フィルをする流体ID
   * @param [in]     bid      境界ID
   * @param [in]     Dsize    サイズ
   * @retval 置換されたセル数
   */
  unsigned long fillByFluid (int* bcd, const int fluid_id, const int* bid, const int* Dsize=NULL);
  
  
  /* @brief 未ペイントセルをフィル
   * @param [in,out] mid      work array
   * @param [in]     target   フィルをするID
   * @param [in]     Dsize    サイズ
   * @retval 置換されたセル数
   */
  unsigned long fillByID(int* mid, const int target, const int* Dsize=NULL);
  
  
  /**
   * @brief 流体媒質のフィルをbid情報を元に実行
   * @param [in,out] mid         work array
   * @param [in]     tgt_id      フィルする流体IDのエントリ
   * @param [in]     Dsize       サイズ
   */
  unsigned long fillByMid(int* mid, const int tgt_id, const int* Dsize=NULL);
  
  
  /* @brief 未ペイントセルをフィル
   * @param [in,out] bcd      BCindex B
   * @param [in]     fluid_id フィルをする流体ID
   * @param [in]     bid      境界ID
   * @retval 置換されたセル数
   */
  unsigned long fillByModalSolid (int* bcd, const int fluid_id, const int* bid);
  
  
  /**
   * @brief 交点が定義点にある場合にそのポリゴンのエントリ番号でフィルする
   * @param [in,out] bcd    BCindex B
   * @param [in]     bid    境界ID（5ビット幅x6方向）
   * @param [in]     cut    カット情報
   * @param [in]     Dsize  サイズ
   * @retval フィルされたセル数
   */
  unsigned long fillCutOnCellCenter (int* bcd, const int* bid, const float* cut, const int* Dsize=NULL);
  
  
  /**
   * @brief シード点をペイントする
   * @param [in,out] bcd    BCindex B
   * @param [in]     face   ヒント面
   * @param [in]     target ペイントするIDのエントリ
   * @param [in]     bid    境界ID
   * @param [in]     Dsize  サイズ
   */
  unsigned long fillSeedBcd(int* bcd, const int face, const int target, const int* bid, const int* Dsize=NULL);
  
  
  /**
   * @brief シード点をペイントする
   * @param [in,out] mid    work array
   * @param [in]     face   ヒント面
   * @param [in]     target ペイントするIDのエントリ
   * @param [in]     Dsize  サイズ
   */
  unsigned long fillSeedMid(int* mid, const int face, const int target, const int* Dsize=NULL);
  
  
  /* @brief 交点が定義点にある場合の修正
   * @param [in,out] bid      境界ID
   * @param [in]     cut      カット情報
   * @param [in]     fluid_id フィルをする流体ID
   * @retval 置換されたセル数
   */
  unsigned long modifyCutOnCellCenter (int* bid, const float* cut, const int fluid_id);
  
  
  /* @brief 領域境界のセルで境界方向にカットがある場合にガイドセルをそのIDで埋める
   * @param [in,out] bcd      BCindex B
   * @param [in]     bid      境界ID
   * @param [out]    painted  各外部領域面でペイントされた個数
   */
  void paintSolidGC (int* bcd, const int* bid, unsigned long* painted);
  
  
  /* @brief 孤立したゼロIDのセルを隣接セルIDで埋める
   * @param [in,out] bcd      BCindex B
   * @param [in]     fluid_id フィルする流体ID
   * @param [in]     bid      境界ID
   */
  unsigned long replaceIsolatedFcell (int* bcd, const int fluid_id, const int* bid);
  
  
  /**
   * @brief bx[]に各境界条件の共通のビット情報をエンコードする
   * @param [in,out] bx    BCindex B
   * @param [in]     bid   交点ID
   * @param [in]     mat   MediumList
   * @param [in,out] cmp   CompoList
   * @param [in,out] Lcell ノードローカルの有効セル数
   * @param [in,out] Gcell グローバルの有効セル数
   * @param [in]     KOS   解くべき方程式の種類 KIND_OF_SOLVER
   */
  void setBCIndexBase (int* bx,
                       const int* bid,
                       const MediumList* mat,
                       CompoList* cmp,
                       unsigned long& Lcell,
                       unsigned long& Gcell,
                       const int KOS);
  
  
  /**
   * @brief 温度境界条件のビット情報をエンコードする
   * @param [in,out] cdf    BCindex C
   * @param [in,out] bd     BCindex B
   * @param [in]     BC     SetBCクラスのポインタ
   * @param [in]     kos    KindOfSolver
   * @param [in,out] cmp    CompoList
   * @param [in]     cut    距離情報
   * @param [in]     cut_id カット点ID
   */
  void setBCIndexH (int* cdf,
                    int* bd,
                    SetBC* BC,
                    const int kos,
                    CompoList* cmp,
                    float* cut,
                    int* cut_id);
  
  
  /**
   * @brief 圧力境界条件のビット情報をエンコードする
   * @param [in,out] bcd      BCindex B
   * @param [in,out] bcp      BCindex P
   * @param [in]     BC       SetBCクラスのポインタ
   * @param [in,out] cmp      CompoList
   * @param [in]     icls     Intrinsic class
   * @param [in]     cut      距離情報
   * @param [in]     bid      カットID情報
   * @param [in]     naive    ナイーブ実装のON/OFF
   * @param [out]    pni      Naive実装の係数
   * @retval 表面セル数
   */
  unsigned long setBCIndexP (int* bcd,
                             int* bcp,
                             SetBC* BC,
                             CompoList* cmp,
                             int icls,
                             const float* cut,
                             const int* bid,
                             const int naive,
                             REAL_TYPE* pni);
  
  
  /**
   * @brief cdf[]に境界条件のビット情報をエンコードする
   * @param [in,out] cdf    BCindex C
   * @param [in,out] bp     BCindex P
   * @param [in]     BC     SetBCクラスのポインタ
   * @param [in]     cmp    CompoListクラスのポインタ
   * @param [in]     icls   Intrinsic class
   * @param [in]     cut    カット配列
   * @param [in]     cut_id BID配列
   */
  void setBCIndexV (int* cdf,
                    int* bp,
                    SetBC* BC,
                    CompoList* cmp,
                    int icls,
                    float* cut,
                    int* cut_id);
  
  
  /**
   * @brief bx[]のコンポーネントエントリを参照して体積率を計算し，圧力損失コンポーネントの場合にはビットを立てる
   * @param [in]     cmp コンポーネントリスト
   * @param [in,out] bx  BCindex B
   * @param [in]     vf  体積率
   */
  void setCmpFraction (CompoList* cmp, int* bx, const REAL_TYPE* vf);
  
  
  /**
   * @brief コンポーネントの操作に必要な定数の設定
   * @param [in] m_NoCompo  コンポーネントの総数
   * @param [in] m_NoMedium 媒質数
   * @param [in] ExRef      組み込み例題クラス
   */
  void setControlVars (const int m_NoCompo, const int m_NoMedium, Intrinsic* ExRef=NULL);
  
  
  /**
   * @brief 外部境界が周期境界の場合の距離情報・境界ID・媒質エントリをセット
   * @param [in]     face      外部境界面番号
   * @param [in]     prdc_mode 周期境界条件のモード
   * @param [in,out] bcd       BCindex B
   * @note 領域境界面は全て流体を想定
   */
  void setOBCperiodic (const int face, const int prdc_mode, int* bcd);
  
  
  /**
   * @brief 外部境界の距離情報・境界ID・媒質エントリをセット
   * @param [in]     face  外部境界面番号
   * @param [in]     c_id  媒質IDエントリ番号
   * @param [in]     str   "SOLID" or "FLUID"
   * @param [in,out] bcd   BCindex B
   * @param [in,out] cut   距離情報
   * @param [in,out] bid   カットID情報
   */
  void setOBC (const int face, const int c_id, const char* str, int* bcd, float* cut, int* bid);
  
};

#endif // _FB_BINVOX_H_
