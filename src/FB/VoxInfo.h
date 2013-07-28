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
  unsigned long count_ValidCell_OBC (const int face, const int* bv, const int typ);

  
  
  void encActive             (unsigned long& Lcell, unsigned long& Gcell, int* bx, const int KOS);
  void encAmask_SymtrcBC     (int face, int* bh2);
  void encHbit               (int* bh1, int* bh2);
  void encPbit               (int* bx);
  void encPbit_OBC           (int face, int* bx, string key, bool dir);
  
  
  // CompoListのエントリをbx[]へエンコードする
  unsigned long encodeOrder(const int order,
                            const int* mid,
                            int* bx,
                            CompoList* cmp);
  
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
  
  
  // 圧力のノイマン境界ビットをエンコードする（カット）
  unsigned long encPbit_N_Cut(int* bx, const float* cut, const bool convergence);
  
  
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
  
  
  
  /**
   * @brief 熱境界条件のBCエントリをエンコードする
   * @retval カウントしたセル数
   * @param [in]     order  CompoListのエントリ
   * @param [in]     target 対象セルID
   * @param [in]     mid    セルID配列
   * @param [in,out] bcd    BCindex ID
   * @param [in,out] bh1    BCindex H1
   * @param [in,out] bh2    BCindex H2
   * @param [in]     deface フラグエンコードの対象セルID
   * @param [in]     flag   断熱ビットのon(true)/off(false)
   * @note
   *   - targetとdefaceで挟まれる面に対して適用
   *   - 対象面の断熱ビットはflagで判断
   */
  unsigned long encQface(const int order,
                         const int target,
                         const int* mid, 
                         int* bcd, 
                         int* bh1, 
                         int* bh2, 
                         const int deface, 
                         const bool flag);
  
  /**
   * @brief 熱境界条件のBCエントリをエンコードする
   * @retval カウントしたセル数
   * @param [in]     order  CompoListのエントリ
   * @param [in]     target テスト対象セルID
   * @param [in]     mid    セルID配列
   * @param [in,out] bcd    BCindex ID
   * @param [in,out] bh1    BCindex H1
   * @param [in,out] bh2    BCindex H2
   * @param [in]     cutid  カットID
   * @param [in]     flag   断熱ビットのon(true)/off(false)
   * @note
   *   - 対象面の断熱ビットはflagで判断
   */
  unsigned long encQface_cut(const int order,
                         const int target,
                         const int* mid,
                         int* bcd,
                         int* bh1,
                         int* bh2,
                         const int cutid,
                         const bool flag);
  
  
  unsigned long encQfaceHT_B(const int order, 
                             const int id, 
                             const int* mid, 
                             int* bcd, 
                             int* bh1, 
                             int* bh2, 
                             const int deface);
  
  unsigned long encQfaceHT_S(const int order, 
                             const int id, 
                             const int* mid, 
                             int* bcd,
                             int* bh1, 
                             int* bh2, 
                             const int deface);
  
  unsigned long encQfaceISO_SF(const int order, 
                               const int id, 
                               const int* mid, 
                               int* bcd, 
                               int* bh1, 
                               int* bh2, 
                               const int deface);
  
  unsigned long encQfaceISO_SS(const int order, 
                               const int id, 
                               const int* mid, 
                               int* bcd, 
                               int* bh1, 
                               int* bh2, 
                               const int deface);
  
  void encQfaceSVO           (int order, int id, int* mid, int* bcd, int* bh1, int* bh2, int deface);
  
  
  // bv[]にVBCの境界条件に必要な情報をエンコードし，流入流出の場合にbp[]の隣接壁の方向フラグを除く
  // 境界条件指定キーセルのSTATEを流体に変更する
  unsigned long encVbitIBC(const int order,
                           const int* mid,
                           int* bv,
                           int* bp,
                           const float* vec,
                           const int bc_dir);
  
  
  // bv[]にVBCの境界条件に必要な情報をエンコード
  unsigned long encVbitIBCcut(const int order,
                              int* bv,
                              int* bp,
                              const float* cut,
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
    unsigned key[NoCompo+1]; // ID毎の頻度
    int val[6];              // ID
    
    
    for (int l=1; l<=NoCompo; l++) key[l]=0;
    
    
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
  inline int get_BID5(const int dir, const int bid) const
  {
    return ( (bid >> dir*5) & MASK_5 );
  }
  
  
  bool isInTable             (int MaxSize, int* cList, int target);
  
  
  
  //@brief idxの第shiftビットをOFFにする
  inline int offBit(int idx, const int shift) {
    return ( idx & (~(0x1<<shift)) );
  }
  
  
  //@brief idxの第shiftビットをONにする
  inline int onBit(int idx, const int shift) {
    return ( idx | (0x1<<shift) );
  }
  
  
  void setAmask_InActive     (int id, int* mid, int* bh);
  void setAmask_Solid        (int* bh);
  void setAmask_Thermal      (int* bh);
  
  
  /*
   * @brief CutBid5のBoundrary ID設定
   * @param [in,out] bid  CutBid5のBoundrary ID
   * @param [in]     dir  方向コード (w/X_MINUS=0, e/X_PLUS=1, s/2, n/3, b/4, t/5)
   * @param [in]     s_id 固体のID (1-31)
   */
  inline void set_BID5(int& bid, const int dir, const int s_id)
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

  
  /**
   * @brief dst[]にsrc[]のstate, activeビットの情報をコピーする
   * @param [out] dst コピー先
   * @param [in]  src コピー元
   */
  void copyBCIbase (int* dst, const int* src);
  
  
  /**
   * @brief ペイント済みの個数を返す
   * @retval 計算空間内における非empty(0)セル数
   * @param [in] mid ボクセルID配列
   */
  unsigned long countPainted(const int* mid);
  
  
  /**
   * @brief セルの状態をカウントして，その個数をLcell, Gcellに保持する
   * @param Lcell ノードローカルの値（戻り値）
   * @param Gcell グローバルの値（戻り値）
   * @param bx BCindex ID
   * @param state カウントするセルの状態
   */
  void countCellState(unsigned long& Lcell, unsigned long& Gcell, int* bx, const int state);
  
  
  /**
   @brief  計算領域の外部境界で外側1層と内側の両方が流体セル数の場合にカウントする
   @param [in]  bx       BCindex ID
   @param [out] OpenArea 開口セル数
   */
  void countOpenAreaOfDomain (int* bx, REAL_TYPE* OpenArea);
  
  
  /**
   * @brief BCindexIDを表示する（デバッグ用）
   * @param [in] bcd   BCindex ID
   * @param [in] fname 出力用ファイル名
   */
  void dbg_chkBCIndexD  (const int* bcd, const char* fname);
  
  
  /**
   * @brief BCindexPを表示する（デバッグ用）
   * @param [in] bcd   BCindex ID
   * @param [in] bcp   BCindex P
   * @param [in] fname 出力用ファイル名
   */
  void dbg_chkBCIndexP  (const int* bcd, const int* bcp, const char* fname);
  
  
  /**
   * @brief BCindexVを表示する（デバッグ用）
   * @param [in] bcv   BCindex V
   * @param [in] fname 出力用ファイル名
   */
  void dbg_chkBCIndexV  (const int* bcv, const char* fname);
  
  
  /**
   * @brief ペイント済みかどうかをチェックし、未ペイントセルがあれば1を返す
   * @param [in] mid ID配列
   */
  int fill_check(const int* mid);
  
  
  /**
   * @brief 流体媒質のフィルをbid情報を元に実行
   * @param [in,out] bid         カット点のID配列（5ビット幅x6方向）
   * @param [in,out] mid         ID配列
   * @param [in,out] cut         カット情報
   * @param [in]     tgt_id      フィルする流体ID
   * @param [out]    substituted 固体IDに置換された数
   * @param [in]     m_list      CellMonitorの格納順がリストアップされた配列
   */
  unsigned fill_by_bid(int* bid, int* mid, float* cut, const int tgt_id, unsigned& substituted, int* m_list);
  
  
  /**
   * @brief 流体媒質のフィルをmid情報を元に実行
   * @param [in,out] bid      カット点のID配列（5ビット幅x6方向）
   * @param [in,out] mid      ID配列
   * @param [in,out] cut      カット情報
   * @param [in]     tgt_id   フィルする流体ID
   * @param [in]     m_list   CellMonitorの格納順がリストアップされた配列
   */
  unsigned fill_by_mid(int* bid, int* mid, float* cut, const int tgt_id, int* m_list);
  
  
  /**
   * @brief ID置換を実行
   * @param [in,out] mid     ID配列
   * @param [in]     target  置換対象ID
   * @param [in]     fill_id 置換ID
   */
  unsigned fillReplace(int* mid, const int target, const int fill_id);
  
  /**
   * @brief シード点をペイントする
   * @param [in,out] mid    ID配列
   * @param [in]     face   ヒント面
   * @param [in]     target ペイントするID
   * @param [in]     cut    距離情報
   */
  unsigned fillSeed(int* mid, const int face, const int target, const float* cut);
  
  
  /**
   * @brief 孤立した流体セルを探し，周囲の個体媒質で置換，BCindexを修正する
   * @param [in]     bx    BCindex ID
   * @attention 事前にbx[]の同期が必要 >> 隣接セルがすべて固体の場合をチェックするため
   */
  void findIsolatedFcell(int* bx);
  
  
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
  unsigned long flipInActive(unsigned long& L,
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
   * @brief 内部フィルを実行する場合の固体IDを求める
   * @param [in] mid ボクセルIDを保持する配列
   */
  int getFillSolidID(const int* mid);
  
  
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
  
  
  // 境界条件のビット情報をエンコードする
  void setBCIndexH(int* bcd, int* bh1, int* bh2, int* mid, SetBC* BC, const int kos, CompoList* cmp,
                   bool isCDS=false, float* cut=NULL, int* cut_id=NULL);
  
  
  /**
   * @brief 圧力境界条件のビット情報をエンコードする
   * @param [in,out] bcd   BCindex ID
   * @param [in,out] bcp   BCindex P
   * @param [in,out] mid   ID配列
   * @param [in]     BC    SetBCクラスのポインタ
   * @param [in,out] cmp   CompoList
   * @param [in]     icls  Intrinsic class
   * @param [in]     cut   距離情報
   * @param [in]     bid   カットID情報
   * @param [in]     isCDS CDS->true
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
                            const bool isCDS);
  
  
  /**
   * @brief bv[]に境界条件のビット情報をエンコードする
   * @param [in,out] bv    BCindex V
   * @param [in]     mid   ボクセルID配列
   * @param [in,out] bp    BCindex P
   * @param [in]     BC    SetBCクラスのポインタ
   * @param [in]     cmp   CompoListクラスのポインタ
   * @param [in]     icls  Intrinsic class
   * @param [in]     isCDS カットかどうか
   * @param [in]     cut   カット配列
   * @param [in]     bid   BID配列
   */
  void setBCIndexV(int* bv,
                   const int* mid,
                   int* bp,
                   SetBC* BC,
                   CompoList* cmp,
                   int icls,
                   bool isCDS=false,
                   float* cut=NULL,
                   int* cut_id=NULL);
  
  
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
   * @brief 計算領域外部のガイドセルに媒質IDをエンコードする
   * @param [in]     face      外部境界面番号
   * @param [in,out] mid       ID配列のデータクラス
   * @param [in]     BCtype    外部境界面の境界条件の種類
   * @param [in]     c_id      媒質格納番号
   * @param [in]     prdc_mode 周期境界条件のモード
   */
  void setMediumOnGC(const int face, int* mid, const int BCtype, const int c_id, const int prdc_mode);
  
  
  /**
   * @brief セルモニターのIDをmid[]にセットする
   * @param [in,out] mid    ID配列
   * @param [in]     bid    カットID情報
   * @param [in]     cut    距離情報
   * @param [in]     target 指定ID
   * @param [in]     fluid  流体のID
   * @param [in]     policy セルモニター幅
   */
  unsigned long setMonitorCellID(int* mid, const int* bid, const float* cut, const int target, const int fluid, const int policy);
  
  
  /**
   * @brief Cell_Monitorで指定するIDでモニタ部分を指定するためのしかけ (SHAPE_BOX, SHAPE_CYLINDER)
   * @param [in] mid  ID配列
   * @param [in] n    BC格納番号
   * @param [in] SM   ShapeMonitorクラス
   * @param [in] cmp  CompoListクラス
   * @param [in] RefL 代表長さ
   */
  void setMonitorShape(int* mid, const int n, ShapeMonitor* SM, CompoList* cmp, const REAL_TYPE RefL);
  
  
  /**
   @brief 外部境界のガイドセルが固体の場合に距離情報をセット
   @param [in]     BC  SetBCクラスのポインタ
   @param [in,out] cut 距離情報
   */
  void setOBCcut(SetBC* BC, float* cut);
  
  
  /**
   * @brief ボクセルモデルにカット情報から得られた固体情報を転写する
   * @param [in,out] mid セルID
   * @param [in]     bid カットID情報
   * @param [in]     cut 距離情報
   * @param [in]     cmp CompoListクラス
   * @retval 固体セル数
   * @attention 境界条件ポリゴンのIDへペイントしないので，予めチェック用のリストをつくっておき，ペイント時に確認
   */
  unsigned long SolidFromCut(int* mid, const int* bid, const float* cut, CompoList* cmp);
  
};

#endif // _FB_BINVOX_H_
