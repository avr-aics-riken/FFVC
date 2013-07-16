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
  int NoBC;                      ///< 境界条件数
  int NoCompo;                   ///< コンポーネントの総数
  int NoVoxID;                   ///< 含まれるIDの数(Local/Global)
  int colorList[MODEL_ID_MAX+1]; ///< ボクセルモデルに含まれるIDのリスト(Global)
  
  Intrinsic *Ex;                 ///< 例題クラスのポインタ

public:
  /** コンストラクタ */
  VoxInfo() {
    NoVoxID = 0;
    NoBC = 0;
    NoCompo = 0;
    
    for (int i=0; i<MODEL_ID_MAX+1; i++) colorList[i]=0;
    
    Ex = NULL;
  }
  
  /**　デストラクタ */
  ~VoxInfo() {}
  
private:
  
  // モデル情報テーブルを作成，0で初期化する
  int* allocTable (const int size);
  
  
  // IDテーブルを表示（デバッグ用）
  void checkColorTable (FILE* fp, const int size, const int* table);
  
  
  // 指定されたIDが計算領域内部にあるかどうかを判定する
  bool chkIDinside(const int id, const int* mid);
  
  
  // 外部境界に接するガイドセルのmid[]にIDを内部周期境界からコピーする
  void copyID_Prdc_Inner(int* mid, const int* m_st, const int* m_ed, const int m_id, const int m_dir);
  
  
  // 媒質idの数を数え，値を返す
  unsigned long countState(const int id, const int* mid);
  
  
  // 外部境界面の有効セル数をカウントする
  unsigned long count_ValidCell_OBC (const int face, const int* bv, const int typ);
  
  
  // CompoListのエントリをbx[]へエンコードする
  unsigned long encodeOrder(const int order,
                            const int id, 
                            const int* mid, 
                            int* bx);
  
  
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
   * @param [in]     target BCのID
   * @param [in]     mid    媒質ID配列
   * @param [in,out] bcd    BCindex ID
   * @param [in,out] bcp    BCindex P
   * @param [in]     vec    法線ベクトル
   * @param [in]     bc_dir 境界条件の方向
   */
  unsigned long encPbit_N_IBC(const int order,
                              const int target,
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
  
  
  // bv[]にVBCの境界条件に必要な情報をエンコードし，流入流出の場合にbp[]の隣接壁の方向フラグを除く
  // 境界条件指定キーセルのSTATEを流体に変更する
  unsigned long encVbit_IBC(const int order,
                            const int target,
                            const int* mid,
                            int* bv,
                            int* bp,
                            const float* vec,
                            const int bc_dir);
  
  
  unsigned long encVbit_IBC_Cut(const int order, 
                                const int id, 
                                int* bv, 
                                int* bp, 
                                const float* cut, 
                                const int* cut_id, 
                                const float* vec, 
                                const int bc_dir);
  void encActive             (unsigned long& Lcell, unsigned long& Gcell, int* bx, const int KOS);
  void encAmask_SymtrcBC     (int face, int* bh2);
  void encHbit               (int* bh1, int* bh2);
  void encPbit               (int* bx);
  void encPbit_OBC           (int face, int* bx, string key, bool dir);
  void encQfaceSVO           (int order, int id, int* mid, int* bcd, int* bh1, int* bh2, int deface);
  
  
  // bv[]に境界条件のビット情報をエンコードする
  void encVbit_OBC(int face, int* bv, string key, const bool enc_sw, string chk, int* bp, bool enc_uwd=false);
  
  
  /**
   * @brief 孤立した流体セルを探し，周囲の個体媒質で置換，BCindexを修正する
   * @param [in]     order cmp[]に登録されたMediumListへのエントリ番号
   * @param [in,out] mid   ボクセルID配列
   * @param [in]     bx    BCindex ID
   * @attention 事前にbx[]の同期が必要 >> 隣接セルがすべて固体の場合をチェックするため
   */
  void find_isolated_Fcell(int order, int* mid, int* bx);
  
  
  /**
   * @brief cmp[]にエンコードされた媒質IDの中から対象となるIDのエントリを探す
   * @param [in] mat_id 対象とする媒質ID
   */
  int find_mat_odr(const int mat_id, CompoList* cmp);
  
  
  void getOffset             (int* st, int* ofst);
  
  
  bool isInTable             (int MaxSize, int* cList, int target);
  
  
  
  //@fn inline int offBit(int idx, const int shift)
  //@brief idxの第shiftビットをOFFにする
  inline int offBit(int idx, const int shift) {
    return ( idx & (~(0x1<<shift)) );
  }
  
  //@fn inline int onBit(int idx, const int shift)
  //@brief idxの第shiftビットをONにする
  inline int onBit(int idx, const int shift) {
    return ( idx | (0x1<<shift) );
  }
  
  
  
  
  void resizeBVcell          (const int* st, const int* ed, int n, int* bx, int* gcbv);
  void resizeBVface          (const int* st, const int* ed, int n, int* bx, int* gcbv);
  void setInactive_Compo     (int id, int def, int* mid, int* bh1, int* bh2);
  void setAmask_InActive     (int id, int* mid, int* bh);
  void setAmask_Solid        (int* bh);
  void setAmask_Thermal      (int* bh);
  void updateGlobalIndex     (const int* st, const int* ed, int n, int* gcbv);
  


  
  
public:
  
  // 計算領域外部のガイドセルに媒質IDをエンコードする
  void adjMedium_on_GC(const int face, int* mid, const int BCtype, const int c_id, const int prdc_mode);
  
  
  /**
   * @brief 外部境界に接するガイドセルのmid[]にIDをエンコードする（内部周期境界の場合）
   * @param [in,out] mid  ID配列のデータ
   * @param [in]     cmp  CompoList
   */
  void adjMediumPrdc_Inner(int* mid, CompoList* cmp);
  
  
  
  /**
   * @brief パラメータファイルとスキャンしたIDの同一性をチェック
   * @param [in] m_NoMedium  MediumTableに記述されたIDの個数
   */
  bool chkIDconsistency      (const int m_NoMedium);
  
  
  void copyBCIbase           (int* dst, int* src);
  
  
  /**
   * @brief セルの状態をカウントして，その個数をLcell, Gcellに保持する
   * @param Lcell ノードローカルの値（戻り値）
   * @param Gcell グローバルの値（戻り値）
   * @param bx BCindex ID
   * @param state カウントするセルの状態
   */
  void countCellState(unsigned long& Lcell, unsigned long& Gcell, int* bx, const int state);
  
  
  void countOpenAreaOfDomain (int* bx, REAL_TYPE* OpenArea);
  
  
  /**
   * @brief ペイント済みかどうかをチェックし、未ペイントセルがあれば1を返す
   * @param [in] mid ID配列
   */
  int fill_check(const int* mid);
  
  
  
  /**
   * @brief 流体媒質のフィルをbid情報を元に実行
   * @param [in,out] bid      カット点のID配列（5ビット幅x6方向）
   * @param [in,out] mid      ID配列
   * @param [in,out] cut      カット情報
   * @param [in]     tgt_id   フィルする流体ID
   * @param [in]     solid_id 固体ID
   * @param [in]     cmp      CompoList
   */
  unsigned fill_by_bid(int* bid, int* mid, float* cut, const int tgt_id, const int solid_id, CompoList* cmp);
  
  
  /**
   * @brief 流体媒質のフィルをmid情報を元に実行
   * @param [in,out] bid      カット点のID配列（5ビット幅x6方向）
   * @param [in,out] mid      ID配列
   * @param [in,out] cut      カット情報
   * @param [in]     tgt_id   フィルする流体ID
   * @param [in]     solid_id 固体ID
   * @param [in]     cmp      CompoList
   */
  unsigned fill_by_mid(int* bid, int* mid, float* cut, const int tgt_id, const int solid_id, CompoList* cmp);
  
  
  /**
   * @brief 内部フィルを実行
   * @param [in,out] mid     ID配列
   * @param [in]     target  置換対象ID
   * @param [in]     fill_id 置換ID
   */
  unsigned fill_inside(int* mid, const int target, const int fill_id);
  
  /**
   * @brief シード点をペイントする
   * @param [in,out] mid    ID配列
   * @param [in]     face   ヒント面
   * @param [in]     target ペイントするID
   * @param [in]     cut    距離情報
   * @param [in]     cmp      CompoList
   */
  unsigned long fill_seed(int* mid, const int face, const int target, const float* cut, CompoList* cmp);
  
  
  unsigned long flip_InActive(unsigned long& L,
                              unsigned long& G, 
                              const int id, 
                              const int* mid, 
                              int* bx);
  
  
  /**
   * @brief IBCのbboxを取得する
   * @param [in]  tgt コンポーネント配列のID
   * @param [in]  bid カットID情報
   * @param [in]  cut カット情報
   * @param [out] st  コンポーネントbboxの開始セル
   * @param [out] ed  コンポーネントbboxの終端セル
   */
  bool find_IBC_bbox(const int tgt, const int* bid, const float* cut, int* st, int* ed);
  
  
  //@fn const int* getColorList() const
  //@retval colorListのポインタ
	const int* getColorList() const { return colorList; }
	
  
  /*
   * @brief 指定面方向のカットIDをとりだす
   * @param [in] dir  方向コード (w/X_MINUS=0, e/X_PLUS=1, s/2, n/3, b/4, t/5)
   * @param [in] bid  CutBid5のBoundrary ID
   */
  inline int get_BID5(const int dir, const int bid) {
    return ( (bid >> dir*5) & MASK_5 );
  }
  
  

  // ボクセルをスキャンしたIDの数と境界条件の数，含まれるIDのリストを表示する
  void printScannedCell(FILE* fp);
  
  
  void resizeCompoBV         (int* bd, int* bv, int* bh1, int* bh2, int kos, bool isHeat, int* gcbv);
  
  

  // cellで保持されるボクセルid配列をスキャンし，coloList[]に登録する
  int scanCell(int *cell, const int* cid, const int ID_replace);
  
  
  /*
   * @brief CutBid5のBoundrary ID設定
   * @param [in,out] bid  CutBid5のBoundrary ID
   * @param [in]     dir  方向コード (w/X_MINUS=0, e/X_PLUS=1, s/2, n/3, b/4, t/5)
   * @param [in]     s_id 固体のID (1-31)
   */
  inline void set_BID5(int& bid, const int dir, const int s_id) {
    bid |= (s_id << (dir*5));
  }
  
  // bx[]に各境界条件の共通のビット情報をエンコードする（その1）
  void setBCIndex_base1(int* bd, int* mid, const float* cvf, const MediumList* mat, CompoList* cmp);
  
  
  // bx[]に各境界条件の共通のビット情報をエンコードする（その2）
  void setBCIndex_base2(int* bx, int* mid, unsigned long& Lcell, unsigned long& Gcell, const int KOS, CompoList* cmp);
  
  
  // 境界条件のビット情報をエンコードする
  void setBCIndexH(int* bcd, int* bh1, int* bh2, int* mid, SetBC* BC, const int kos, CompoList* cmp,
                   bool isCDS=false, float* cut=NULL, int* cut_id=NULL);
  
  
  // 圧力境界条件のビット情報をエンコードする
  unsigned long setBCIndexP(int* bcd,
                            int* bcp,
                            int* mid,
                            SetBC* BC,
                            CompoList* cmp,
                            int icls,
                            const float* cut,
                            const int* bid,
                            const bool isCDS);
  
  
  // bv[]に境界条件のビット情報をエンコードする
  void setBCIndexV(int* bv,
                   const int* mid,
                   int* bp,
                   SetBC* BC,
                   CompoList* cmp,
                   int icls,
                   bool isCDS=false,
                   float* cut=NULL,
                   int* cut_id=NULL);
  
  
  void setCmpFraction        (CompoList* compo, int* bx, float* vf);
  
  void setAdiabatic4SF       (int* bh);
  
  
  /*
   * @brief 組み込み例題クラスをアサイン
   * @param [in] ExRef  組み込み例題クラス
   */
  void setIntrinsic(Intrinsic* ExRef)
  {
    Ex = ExRef;
  }
  
  
  void setNoCompo_BC         (int m_NoBC, int m_NoCompo);
  void setOBC_Cut            (SetBC* BC, float* cut);
  
  
  // Cell_Monitorで指定するIDでモニタ部分を指定するための準備
  void setShapeMonitor(int* mid, ShapeMonitor* SM, CompoList* cmp, const REAL_TYPE RefL);
  
  
  // カット情報を用いて，指定IDからバイナリボクセルを作成する
  //unsigned long Solid_from_Cut(int* mid, const int* bid, const float* cut, const int target);
  
  
  // ボクセルモデルにカット情報から得られた固体情報を転写する
  unsigned long Solid_from_Cut(int* mid, const int* bid, const float* cut, CompoList* cmp);
  
  
  // // VBCコンポーネントのバックリング
  unsigned long Solid_from_Cut_VBC(int* mid,
                                   const int target,
                                   const int solid_id,
                                   const float* vec,
                                   const int bc_dir);
  
  // ----> debug function
  
  /**
   * @brief BCindexを表示する（デバッグ用）
   * @param [in] bcd   BCindex D
   * @param [in] fname 出力用ファイル名
   */
  void dbg_chkBCIndexD  (const int* bcd, const char* fname);
  
  /**
   * @brief BCindexを表示する（デバッグ用）
   * @param [in] bcd   BCindex ID
   * @param [in] bcp   BCindex P
   * @param [in] fname 出力用ファイル名
   * @param [in] cmp CompoListクラス
   */
  void dbg_chkBCIndexP  (const int* bcd, const int* bcp, const char* fname, CompoList* cmp);
  
  void dbg_chkBCIndexV  (int* bcv, const char* fname);
};

#endif // _FB_BINVOX_H_
