#ifndef _FB_BINVOX_H_
#define _FB_BINVOX_H_

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
// Copyright (c) 2012-2016 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

/**
 * @file   VoxInfo.h
 * @brief  FlowBase VoxInfo class Header
 * @author aics
 */

#include "DomainInfo.h"
#include <math.h>
#include <cstdlib>
#include "FBUtility.h"
#include "Component.h"
#include "Medium.h"
#include "SetBC.h"
#include "BndOuter.h"
#include "common/Vec3.h" // defined in Polylib
#include "Intrinsic.h"
#include "limits.h"
#include "omp.h"

using namespace Vec3class;

class VoxInfo : public DomainInfo {

private:
  Intrinsic *Ex;       ///< 例題クラスのポインタ

public:
  /** コンストラクタ */
  VoxInfo() {
    Ex = NULL;
  }

  /** デストラクタ */
  ~VoxInfo() {}

private:

  // 外部境界に接するガイドセルのbcd[]にIDを内部周期境界からコピーする
  void copyIdPrdcInner (int* bcd, const int* m_st, const int* m_ed, const int m_id, const int m_dir);


  // セルセンターのIDから対象セル数をカウントし，サブドメイン内にコンポーネントがあれば存在フラグを立てる
  unsigned long countCC (const int order, const int* bcd, CompoList* cmp);


  // セルフェイスの交点IDから対象IDのセル数をカウントし，サブドメイン内にコンポーネントがあれば存在フラグを立てる
  unsigned long countCF (const int key, const int* bcd, const int* bid, CompoList* cmp, const string attrb);


  // 外部境界面の有効セル数をカウントする
  unsigned long countValidCellOBC (const int face, const int* cdf, const int typ);


  // BCindexにそのセルが計算に有効(active)かどうかをエンコードする
  void encActive (unsigned long& Lcell, unsigned long& Gcell, int* bx, const int KOS);


  // 断熱マスクのエンコード
  void encAdiabatic (int* bd, const string target, const int* bid, int face=-1);


  // セルの各面を調べ，境界条件が設定されていれば，ビットをON
  void encHbit (const int* cdf, int* bd);


  // ディリクレ条件とノイマン条件の排他性をチェックし，反復行列の非対角要素/対角要素の係数をエンコードする
  void encPbit (int* bx, const int* bid);


  // 圧力のノイマン境界ビットをエンコードする
  void encPbitN (int* bx, const int* bid, const int* cut_l, const int* cut_u, const bool convergence);


  // 計算領域内部のコンポーネントの圧力境界条件フラグをbcp[]にエンコードする
  unsigned long encPbitIBC (const int order,
                            int* bcd,
                            int* bcp,
                            const int* bid,
                            const REAL_TYPE* vec,
                            const string condition);


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
                            int* bid,
                            int* bcd,
                            CompoList* cmp,
                            const int solid_id);


  unsigned long encVbitIBCrev (const int order,
                               int* cdf,
                               int* bid,
                               int* bcd,
                               CompoList* cmp,
                               const int solid_id);


  // cdf[]に境界条件のビット情報をエンコードする
  void encVbitOBC (int face, int* cdf, string key, const bool enc_sw, string chk, int* bid, int* bcd);


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
   * @param [in,out] bcd        BCindex B
   * @param [in]     cmp        CompoList
   * @param [in]     m_NoCompo  コンポーネントの総数
   */
  void adjMediumPrdcInner(int* bcd, CompoList* cmp, const int m_NoCompo);


  /**
   * @brief BCindex Bのエンコードを確認
   * @param [in] bcd BCindex B
   */
  void chkOrder (const int* bcd);


  // フラグのチェック
  void chkFlag(const int* bx, const int* bid, const int* cut_l, const int* cut_u, const int* bcd);


  /**
   * @brief dst[]にsrc[]のstate, activeビットの情報をコピーする
   * @param [out] dst コピー先
   * @param [in]  src コピー元
   */
  void copyBCIbase (int* dst, const int* src);


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


  /* @brief 外部境界方向にカットがあるセルにはガイドセルをCutIDの媒質でペイント
   * @param [in,out] bcd      BCindex B
   * @param [in]     bid      境界ID
   * @param [out]    painted  各外部領域面でペイントされた個数
   * @param [in]     cmp      CompoList
   */
  void paintCutIDonGC (int* bcd, const int* bid, unsigned long* painted, const CompoList* cmp);


  /**
   * @brief bx[]に各境界条件の共通のビット情報をエンコードする
   * @param [in,out] bcd        BCindex B
   * @param [in]     bid        交点ID
   * @param [in]     mat        MediumList
   * @param [in,out] cmp        CompoList
   * @param [in,out] Lcell      ノードローカルの有効セル数
   * @param [in,out] Gcell      グローバルの有効セル数
   * @param [in]     KOS        解くべき方程式の種類 KIND_OF_SOLVER
   * @param [in]     m_NoMedium 媒質数
   * @param [in]     m_NoCompo  コンポーネント数
   */
  void setBCIndexBase(int* bcd,
                      const int* bid,
                      const MediumList* mat,
                      CompoList* cmp,
                      unsigned long& Lcell,
                      unsigned long& Gcell,
                      const int KOS,
                      const int m_NoMedium,
                      const int m_NoCompo);


  /**
   * @brief 温度境界条件のビット情報をエンコードする
   * @param [in,out] cdf        BCindex C
   * @param [in,out] bcd        BCindex B
   * @param [in]     BC         SetBCクラスのポインタ
   * @param [in]     kos        KindOfSolver
   * @param [in,out] cmp        CompoList
   * @param [in]     cut        距離情報
   * @param [in]     cut_id     カット点ID
   * @param [in]     m_NoCompo  コンポーネント数
   */
  void setBCIndexH(int* cdf,
                   int* bcd,
                   SetBC* BC,
                   const int kos,
                   CompoList* cmp,
                   int* cut_l,
                   int* cut_u,
                   int* cut_id,
                   const int m_NoCompo);


  /**
   * @brief 圧力境界条件のビット情報をエンコードする
   * @param [in,out] bcd        BCindex B
   * @param [in,out] bcp        BCindex P
   * @param [in]     BC         SetBCクラスのポインタ
   * @param [in,out] cmp        CompoList
   * @param [in]     icls       Intrinsic class
   * @param [in]     cut        距離情報
   * @param [in]     bid        カットID情報
   * @param [in]     m_NoCompo  コンポーネント数
   */
  void setBCIndexP(int* bcd,
                   int* bcp,
                   SetBC* BC,
                   CompoList* cmp,
                   int icls,
                   const int* cut_l,
                   const int* cut_u,
                   const int* bid,
                   const int m_NoCompo);


  /**
   * @brief cdf[]に境界条件のビット情報をエンコードする
   * @param [in,out] cdf        BCindex C
   * @param [in]     BC         SetBCクラスのポインタ
   * @param [in]     cmp        CompoListクラスのポインタ
   * @param [in]     icls       Intrinsic class
   * @param [in]     cut        カット配列
   * @param [in]     bid        BID配列
   * @param [in]     bcd        BCindex B
   * @param [in]     m_NoCompo  コンポーネント数
   * @param [in]     m_NoMedium 媒質数
   * @param [in]     mat        MediumList
   */
  void setBCIndexV(int* cdf,
                   SetBC* BC,
                   CompoList* cmp,
                   int icls,
                   int* cut_l,
                   int* cut_u,
                   int* bid,
                   int* bcd,
                   const int m_NoCompo,
                   const int m_NoMedium,
                   MediumList* mat);


  /**
   * @brief bx[]のコンポーネントエントリを参照して体積率を計算し，圧力損失コンポーネントの場合にはビットを立てる
   * @param [in]     cmp        コンポーネントリスト
   * @param [in,out] bx         BCindex B
   * @param [in]     vf         体積率
   * @param [in]     m_NoCompo  コンポーネント数
   */
  void setCmpFraction(CompoList* cmp,
                      int* bx,
                      const REAL_TYPE* vf,
                      const int m_NoCompo);


  /**
   * @brief コンポーネントの操作に必要な定数の設定
   * @param [in] ExRef      組み込み例題クラス
   */
  void setControlVars(Intrinsic* ExRef=NULL);


  /**
   * @brief 外部境界が周期境界の場合の距離情報・境界ID・媒質エントリをセット
   * @param [in,out] bcd       BCindex B
   * @param [in]     ens       周期境界方向フラグ
   * @note 領域境界面は全て流体を想定
   */
  void setOBCperiodic (int* bcd, const int* ens);


  /**
   * @brief 外部境界の距離情報・境界ID・媒質エントリをセット
   * @param [in]     face    外部境界面番号
   * @param [in]     c_id    媒質IDエントリ番号
   * @param [in]     ptr_cmp CompoListへのポインタ
   * @param [in]     str     "SOLID" or "FLUID"
   * @param [in,out] bcd     BCindex B
   * @param [in,out] cut     距離情報
   * @param [in,out] bid     カットID情報
   */
  void setOBC (const int face,
               const int c_id,
               const int ptr_cmp,
               const char* str,
               int* bcd,
               int* cut_l,
               int* cut_u,
               int* bid);

};

#endif // _FB_BINVOX_H_
