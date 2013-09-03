#ifndef _FB_BINVOX_H_
#define _FB_BINVOX_H_

//##################################################################################
//
// Flow Base class
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2013 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2013 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

/** 
 * @file   VoxInfo.h
 * @brief  FlowBase VoxInfo class Header
 * @author kero
 */

#include <math.h>
#include <cstdlib>

#include "cpm_ParaManager.h"

#include "DomainInfo.h"
#include "FBUtility.h"
#include "Component.h"
#include "Medium.h"
#include "SetBC.h"
#include "BndOuter.h"
#include "vec3.h"
#include "CompoFraction.h"
#include "Intrinsic.h"

#include "mpi.h"
#include "limits.h"
#include "omp.h"

class VoxInfo : public DomainInfo {
  
private:
  int NoCompo;       ///< コンポーネントの総数
  Intrinsic *Ex;     ///< 例題クラスのポインタ

public:
  /** コンストラクタ */
  VoxInfo() {
    NoCompo = 0;
    Ex = NULL;
  }
  
  /**　デストラクタ */
  ~VoxInfo() {}
  
private:
  
  
  // 外部境界に接するガイドセルのmid[]にIDを内部周期境界からコピーする
  void copyID_Prdc_Inner(int* mid, const int* m_st, const int* m_ed, const int m_id, const int m_dir);
  
  
  // 外部境界面の有効セル数をカウントする
  unsigned long countValidCellOBC (const int face, const int* bv, const int typ);

  
  // BCindexにそのセルが計算に有効(active)かどうかをエンコードする
  void encActive(unsigned long& Lcell, unsigned long& Gcell, int* bx, const int KOS);
  
  
  // 外部境界に接するセルに対称境界面の断熱マスクをセットする
  void encAmaskSymtrcBC(const int face, int* bh2);
  
  
  // セルの各面を調べ，境界条件が設定されていれば，ビットをON
  void encHbit(int* bh1, int* bh2);

  
  // CompoListのエントリをbx[]へエンコードする
  unsigned long encodeOrder(const int order,
                            const int* mid,
                            int* bx,
                            CompoList* cmp);
  
  
  // ディリクレ条件とノイマン条件の排他性をチェックし，反復行列の非対角要素/対角要素の係数をエンコードする
  void encPbit(int* bx);
  
  
  // 圧力のディリクレ境界ビットをエンコードする
  unsigned long encPbit_D_IBC(const int order, 
                              const int id, 
                              const int* mid, 
                              int* bcd, 
                              int* bcp, 
                              const int deface);
  
  
  // 圧力のノイマン境界ビットをエンコードする（バイナリボクセル）
  unsigned long encPbit_N_Binary(int* bx);
  
  
  // 圧力のノイマン境界ビットをエンコードする（カット）
  unsigned long encPbit_N_Cut(int* bx, const int* bid, const float* cut, const bool convergence);
  
  
  /**
   * @brief 計算領域内部のコンポーネントのNeumannフラグをbcp[]にエンコードする
   * @retval エンコードしたセル数
   * @param [in]     order  cmp[]のエントリ番号
   * @param [in]     mid    媒質ID配列
   * @param [in,out] bcd    BCindex ID
   * @param [in,out] bcp    BCindex P
   * @param [in]     vec    法線ベクトル
   * @param [in]     bc_dir 境界条件の方向
   */
  unsigned long encPbit_N_IBC(const int order,
                              const int* mid,
                              int* bcd,
                              int* bcp,
                              const float* vec,
                              const int bc_dir);
  
  
  // 外部境界に接するセルにおいて，bx[]に圧力境界条件keyに対応するビットフラグを設定する
  void encPbitOBC(int face, int* bx, string key, bool dir);
  
  
  // 熱境界条件のBCエントリをエンコードする
  unsigned long encQface(const int order,
                         const int* bid,
                         int* bh1,
                         int* bh2,
                         const bool flag,
                         const int target,
                         CompoList* cmp);
  
  
  void encQfaceSVO           (int order, int id, int* mid, int* bcd, int* bh1, int* bh2, int deface);
  
  
  // bv[]にVBCの境界条件に必要な情報をエンコードし，流入流出の場合にbp[]の隣接壁の方向フラグを除く
  // 境界条件指定キーセルのSTATEを流体に変更する
  unsigned long encVbitIBC(const int order,
                           int* bv,
                           int* bp,
                           const int* cut_id,
                           const float* vec,
                           const int bc_dir,
                           CompoList* cmp);
  
  
  // bv[]に境界条件のビット情報をエンコードする
  void encVbitOBC(int face, int* bv, string key, const bool enc_sw, string chk, int* bp, bool enc_uwd=false);
  
  
  /**
   * @brief 隣接６方向の最頻値IDを求める
   * @param [in] fid 流体のID
   * @param [in] qw  w方向のID
   * @param [in] qe  e方向のID
   * @param [in] qs  s方向のID
   * @param [in] qn  n方向のID
   * @param [in] qb  b方向のID
   * @param [in] qt  t方向のID
   */
  inline int find_mode_id(const int fid, const int qw, const int qe, const int qs, const int qn, const int qb, const int qt)
  {
    unsigned key[64]; ///< ID毎の頻度 @note ffv_Initialize() >> fill()でif ( C.NoCompo+1 > 64 )をチェック
    int val[6];       ///< ID
    
    memset(key, 0, sizeof(unsigned)*64);
    
    val[0] = qw;
    val[1] = qe;
    val[2] = qs;
    val[3] = qn;
    val[4] = qb;
    val[5] = qt;
    
    
    // 周囲6方向をテスト
    for (int l=0; l<6; l++)
    {
      if ( val[l] != fid ) // フィルする流体IDでない
      {
        key[ val[l] ]++;
      }
    }
    
    
    int mode = key[NoCompo]; // サーチの初期値，IDの大きい方から
    int z = NoCompo;         // 最頻値のID
    
    for (int l=NoCompo-1; l>=1; l--)
    {
      if ( key[l] > mode )
      {
        mode = key[l];
        z = l;
      }
    }
    
    return z;
  }
  
  
  /*
   * @brief 指定面方向のカットIDをとりだす
   * @param [in] dir  方向コード (w/X_MINUS=0, e/X_PLUS=1, s/2, n/3, b/4, t/5)
   * @param [in] bid  CutBid5のBoundrary ID
   */
  inline int getFaceBID(const int dir, const int bid) const
  {
    return ( (bid >> dir*5) & MASK_5 );
  }
  
  
  //@brief idxの第shiftビットをOFFにする
  inline int offBit(int idx, const int shift)
  {
    return ( idx & (~(0x1<<shift)) );
  }
  
  
  //@brief idxの第shiftビットをONにする
  inline int onBit(int idx, const int shift)
  {
    return ( idx | (0x1<<shift) );
  }
  
  
  // 不活性セルに対する断熱マスクをエンコード
  void setAmaskInactive(int* bh);
  
  
  //各方向の張るスクリーンにコンポーネント要素を投影し，面積を求める
  void projectCompo(const int* st, const int* ed, const int target, const int* bid, int* ss);
  
  
  // KOSがSOLID_CONDUCTIONの場合の断熱マスクの処理
  void setAmaskSolid(int* bh);
  
  
  // KOSがTHERMAL_FLOW, THERMAL_FLOW_NATURALの場合の断熱マスクの処理
  void setAmaskThermal(int* bh);
  
  
  /*
   * @brief CutBid5のBoundrary ID設定
   * @param [in,out] bid  CutBid5のBoundrary ID
   * @param [in]     dir  方向コード (w/X_MINUS=0, e/X_PLUS=1, s/2, n/3, b/4, t/5)
   * @param [in]     s_id ID (1-31)
   */
  inline void setFaceBID(int& bid, const int dir, const int s_id)
  {
    bid |= (s_id << (dir*5));
  }
  
  
  void setInactive_Compo(int id, int def, int* mid, int* bh1, int* bh2);


  
public:
  
  /**
   * @brief 外部境界に接するガイドセルのmid[]にIDをエンコードする（内部周期境界の場合）
   * @param [in,out] mid  ID配列のデータ
   * @param [in]     cmp  CompoList
   */
  void adjMediumPrdc_Inner(int* mid, CompoList* cmp);

  
  // コンポーネントのbid情報から近似的な面積と法線を計算する
  void calCompoArea(const float dhd, const int* bid, CompoList* cmp, FILE* fp);
  
  
  /**
   * @brief dst[]にsrc[]のstate, activeビットの情報をコピーする
   * @param [out] dst コピー先
   * @param [in]  src コピー元
   */
  void copyBCIbase (int* dst, const int* src);
  
  
  /**
   * @brief セル数をカウント
   * @param [in] mid     ID配列
   * @param [in] painted ID=0以外でペイント済みを求める(true), m_idのセルをカウント(false)
   * @param [in] m_id    検査するID
   */
  unsigned long countCell(const int* mid, bool painted=true, int m_id=0);
  
  
  /**
   * @brief セルの状態をカウントして，その個数をLcell, Gcellに保持する
   * @param [out] Lcell ノードローカルの値（戻り値）
   * @param [out] Gcell グローバルの値（戻り値）
   * @param [in]  bx    BCindex ID
   * @param [in]  state カウントするセルの状態
   */
  void countCellState(unsigned long& Lcell, unsigned long& Gcell, const int* bx, const int state);
  
  
  /**
   @brief  計算領域の外部境界で外側1層と内側の両方が流体セル数の場合にカウントする
   @param [in]  bx       BCindex ID
   @param [out] OpenArea 開口セル数
   */
  void countOpenAreaOfDomain (const int* bx, REAL_TYPE* OpenArea);
  
  
  /**
   * @brief BCindexIDを表示する（デバッグ用）
   * @param [in] bcd   BCindex ID
   * @param [in] fname 出力用ファイル名
   */
  void dbg_chkBCIndexD (const int* bcd, const char* fname);
  
  
  /**
   * @brief BCindexH1を表示する（デバッグ用）
   * @param [in] bh    BCindex H1
   * @param [in] fname 出力用ファイル名
   */
  void dbg_chkBCIndexH1 (const int* bh, const char* fname);
  
  
  /**
   * @brief BCindexH2を表示する（デバッグ用）
   * @param [in] bh    BCindex H2
   * @param [in] fname 出力用ファイル名
   */
  void dbg_chkBCIndexH2 (const int* bh, const char* fname);
  
  
  /**
   * @brief BCindexPを表示する（デバッグ用）
   * @param [in] bcd   BCindex ID
   * @param [in] bcp   BCindex P
   * @param [in] fname 出力用ファイル名
   */
  void dbg_chkBCIndexP (const int* bcd, const int* bcp, const char* fname);
  
  
  /**
   * @brief BCindexVを表示する（デバッグ用）
   * @param [in] bcv   BCindex V
   * @param [in] fname 出力用ファイル名
   */
  void dbg_chkBCIndexV (const int* bcv, const char* fname);
  
  
  /**
   * @brief 流体媒質のフィルをbid情報を元に実行
   * @param [in,out] bid         カット点のID配列（5ビット幅x6方向）
   * @param [in,out] mid         ID配列
   * @param [in,out] cut         カット情報
   * @param [in]     tgt_id      フィルする流体ID
   * @param [out]    substituted 固体IDに置換された数
   * @param [in]     m_list      CellMonitorの格納順がリストアップされた配列
   */
  unsigned long fillByBid (int* bid, int* mid, float* cut, const int tgt_id, unsigned long& substituted, const int* m_list);
  
  
  /**
   * @brief 流体媒質のフィルをmid情報を元に実行
   * @param [in,out] bid      カット点のID配列（5ビット幅x6方向）
   * @param [in,out] mid      ID配列
   * @param [in,out] cut      カット情報
   * @param [in]     tgt_id   フィルする流体ID
   * @param [in]     m_list   CellMonitorの格納順がリストアップされた配列
   */
  unsigned long fillByMid (int* bid, int* mid, float* cut, const int tgt_id, const int* m_list);
  
  
  /* @brief targetセルの周囲の固体最頻値でフィルを実行
   * @param [in,out] mid      ID配列
   * @param [in]     target   置換対象セルID
   * @param [in]     fluid_id フィルをする流体ID
   * @retval 置換されたセル数
   */
  unsigned long fillByModalSolid (int* mid, const int target, const int fluid_id);
  
  
  /**
   * @brief ID置換を実行
   * @param [in,out] mid     ID配列
   * @param [in]     target  置換対象ID
   * @param [in]     fill_id 置換ID
   * @retval 置換されたセル数
   */
  unsigned long fillReplace (int* mid, const int target, const int fill_id);
  
  /**
   * @brief シード点をペイントする
   * @param [in,out] mid    ID配列
   * @param [in]     face   ヒント面
   * @param [in]     target ペイントするID
   * @param [in]     cut    距離情報
   */
  unsigned long fillSeed (int* mid, const int face, const int target, const float* cut);
  
  
  /* @brief 孤立したゼロIDのセルを隣接セルIDで埋める
   * @param [in,out] mid      ID配列
   * @param [in]     fluid_id フィルする流体ID
   */
  unsigned long fillIsolatedEmptyCell (int* mid, const int fluid_id);
  
  
  /**
   * @brief 孤立した流体セルを探し，周囲の個体媒質で置換，BCindexを修正する
   * @param [in,out] bx       BCindex ID
   * @param [in]     fluid_id フィルする流体ID
   * @attention 事前にbx[]の同期が必要 >> 隣接セルがすべて固体の場合をチェックするため
   */
  unsigned long findIsolatedFcell (int* bx, const int fluid_id);
  
  
  /**
   * @brief BCindexにActive/Inactiveをエンコード
   * @retval 不活性化した数
   * @param [out]    L   ノードローカルの不活性セル数
   * @param [out]    G   グローバルの不活性セル数
   * @param [in]     id  セルID
   * @param [in]     mid ボクセル配列
   * @param [in,out] bx  BCindex ID
   * @param [in,out] cmp CompoList
   */
  unsigned long flipInactive(unsigned long& L,
                             unsigned long& G,
                             const int id,
                             const int* mid,
                             int* bx,
                             CompoList* cmp);
  
  
  /**
   * @brief IBCのbboxを取得する
   * @param [in]  tgt    コンポーネント配列のID
   * @param [in]  bid    カットID情報
   * @param [in]  cut    カット情報
   * @param [out] st     コンポーネントbboxの開始セル
   * @param [out] ed     コンポーネントbboxの終端セル
   * @param [in]  policy モニターのセル幅
   */
  bool findLBCbbox(const int tgt, const int* bid, const float* cut, int* st, int* ed, const int policy);
  
  
  /**
   * @brief cellで保持されるボクセルid配列をスキャンする
   * @param [in,out] mid        ボクセルIDを保持する配列
   * @param [in]     colorList  IDリスト
   * @param [in]     ID_replace ID=0を置換するID
   * @param [in]     fp         file pointer
   */
  void scanCell(int *mid, const int* colorList, const int ID_replace, FILE* fp);
  
  
  // bx[]に各境界条件の共通のビット情報をエンコードする
  void setBCIndexBase (int* bd,
                       int* mid,
                       const float* cvf,
                       const MediumList* mat,
                       CompoList* cmp,
                       unsigned long& Lcell,
                       unsigned long& Gcell,
                       const int KOS);
  
  
  /**
   * @brief 温度境界条件のビット情報をエンコードする
   * @param [in,out] bh1    BCindex H1
   * @param [in,out] bh1    BCindex H2
   * @param [in,out] mid    ID配列
   * @param [in]     BC     SetBCクラスのポインタ
   * @param [in]     kos    KindOfSolver
   * @param [in,out] cmp    CompoList
   * @param [in]     cut    距離情報
   * @param [in]     cut_id カット点ID
   */
  void setBCIndexH(int* bh1,
                   int* bh2,
                   int* mid,
                   SetBC* BC,
                   const int kos,
                   CompoList* cmp,
                   float* cut,
                   int* cut_id);
  
  
  /**
   * @brief 圧力境界条件のビット情報をエンコードする
   * @param [in,out] bcd      BCindex ID
   * @param [in,out] bcp      BCindex P
   * @param [in,out] mid      ID配列
   * @param [in]     BC       SetBCクラスのポインタ
   * @param [in,out] cmp      CompoList
   * @param [in]     icls     Intrinsic class
   * @param [in]     cut      距離情報
   * @param [in]     bid      カットID情報
   * @param [in]     isBinary バイナリの場合true
   * @retval 表面セル数
   */
  unsigned long setBCIndexP(int* bcd,
                            int* bcp,
                            int* mid,
                            SetBC* BC,
                            CompoList* cmp,
                            int icls,
                            const float* cut,
                            const int* bid,
                            const bool isBinary);
  
  
  /**
   * @brief bv[]に境界条件のビット情報をエンコードする
   * @param [in,out] bv     BCindex V
   * @param [in]     mid    ボクセルID配列
   * @param [in,out] bp     BCindex P
   * @param [in]     BC     SetBCクラスのポインタ
   * @param [in]     cmp    CompoListクラスのポインタ
   * @param [in]     icls   Intrinsic class
   * @param [in]     cut    カット配列
   * @param [in]     cut_id BID配列
   */
  void setBCIndexV(int* bv,
                   const int* mid,
                   int* bp,
                   SetBC* BC,
                   CompoList* cmp,
                   int icls,
                   float* cut,
                   int* cut_id);
  
  
  /**
   * @brief bx[]のコンポーネントエントリを参照して体積率を計算し，圧力損失コンポーネントの場合にはビットを立てる
   * @param [in]     cmp コンポーネントリスト
   * @param [in,out] bx  BCindex ID
   * @param [in]     vf  体積率
   */
  void setCmpFraction(CompoList* cmp, int* bx, const float* vf);
  
  
  /**
   * @brief コンポーネントの操作に必要な定数の設定
   * @param [in] m_NoCompo コンポーネントの総数
   * @param [in] ExRef  組み込み例題クラス
   */
  void setControlVars (const int m_NoCompo, Intrinsic* ExRef);
  
  
  /**
   * @brief 計算領域外部のガイドセルに媒質IDをエンコードする（周期境界以外の場合）
   * @param [in]     face      外部境界面番号
   * @param [in,out] mid       ID配列のデータクラス
   * @param [in]     c_id      媒質格納番号
   */
  void setMediumOnGC (const int face, int* mid, const int c_id);
  
  
  /**
   * @brief 計算領域外部のガイドセルに媒質IDをエンコードする（周期境界の場合）
   * @param [in]     face      外部境界面番号
   * @param [in,out] mid       ID配列のデータクラス
   * @param [in]     prdc_mode 周期境界条件のモード
   */
  void setMediumOnGCperiodic (const int face, int* mid, const int prdc_mode);
  
  /**
   * @brief セルモニターのIDをmid[]にセットする
   * @param [in,out] mid    ID配列
   * @param [in]     bid    カットID情報
   * @param [in]     cut    距離情報
   * @param [in]     target 指定ID
   * @param [in]     fluid  流体のID
   * @param [in]     policy セルモニター幅
   */
  unsigned long setMonitorCellID (int* mid, const int* bid, const float* cut, const int target, const int fluid, const int policy);
  
  
  /**
   * @brief Cell_Monitorで指定するIDでモニタ部分を指定するためのしかけ (SHAPE_BOX, SHAPE_CYLINDER)
   * @param [in] mid  ID配列
   * @param [in] n    BC格納番号
   * @param [in] SM   ShapeMonitorクラス
   * @param [in] cmp  CompoListクラス
   * @param [in] RefL 代表長さ
   */
  void setMonitorShape (int* mid, const int n, ShapeMonitor* SM, CompoList* cmp, const REAL_TYPE RefL);
  
  
  /**
   * @brief 外部境界のガイドセルが固体の場合に距離情報をセット
   * @param [in]     BC  SetBCクラスのポインタ
   * @param [in,out] cut 距離情報
   * @param [in]     bid カットID情報
   */
  void setOBCcut (SetBC* BC, float* cut, int* bid);
  
  
  /**
   * @brief ボクセルモデルにカット情報から得られた固体情報を転写する
   * @param [in,out] mid セルID
   * @param [in]     bid カットID情報
   * @param [in]     cut 距離情報
   * @param [in]     cmp CompoListクラス
   * @retval 固体セル数
   * @attention 境界条件ポリゴンのIDへペイントしないので，予めチェック用のリストをつくっておき，ペイント時に確認
   */
  unsigned long SolidFromCut (int* mid, const int* bid, const float* cut, CompoList* cmp);
  
};

#endif // _FB_BINVOX_H_
