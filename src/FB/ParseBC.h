#ifndef _FB_PARA_BC_H_
#define _FB_PARA_BC_H_

// #################################################################
//
// CAERU Library
//
// Copyright (c) All right reserved. 2012
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
//
// #################################################################

//@file   ParseBC.h
//@brief  FlowBase ParseBC class Header
//@author kero

#include "string.h"

#include "cpm_Define.h"
#include "cpm_ParaManager.h"

#include "DomainInfo.h"
#include "FB_Define.h"
#include "FBUtility.h"
#include "BndOuter.h"
#include "Control.h"
#include "Component.h"
#include "Medium.h"
#include "vec3.h"
#include "Intrinsic.h"
#include "TPControl.h"


class ParseBC : public DomainInfo {
private:
  
  // NoCompo = NoBC + NoMedium
  int NoBC;        ///< LocalBCの数
  int NoBaseBC;    ///< 入力ファイルに記載された外部境界条件の指定種類数
  int NoMedium;    ///< 媒質数
  int NoCompo;     ///< コンポーネントリストの全登録数
  
  REAL_TYPE RefVelocity, BaseTemp, DiffTemp, RefDensity, RefSpecificHeat;
  REAL_TYPE RefLength, BasePrs;
  REAL_TYPE rho, nyu, cp, lambda, beta; // 無次元化の参照値

  int KindOfSolver;
  int Unit_Param;
  int monitor;
  int Unit_Temp;
  int Unit_Prs;
  int Mode_Gradp;
  bool HeatProblem;
  bool isCDS;
  
  TPControl* tpCntl;      ///< テキストパーサーのラッパークラス
  CompoList* compo;       ///< コンポーネントテーブル
  BoundaryOuter* BaseBc;  ///< テンポラリのテーブル


public:
  /** コンストラクタ */
  ParseBC(){
    KindOfSolver = 0;
    BaseTemp = DiffTemp = BasePrs = 0.0;
    RefVelocity = RefDensity = RefSpecificHeat = RefLength = 0.0;
    monitor = 0;
    Unit_Param = 0;
    Unit_Temp = 0;
    Unit_Prs = 0;
    Mode_Gradp = 0;
    HeatProblem = isCDS = false;
    NoBaseBC = 0;
    NoBC     = 0;
    NoCompo  = 0; 
    NoMedium = 0;
    
    BaseBc = NULL;
    compo = NULL;
  }
  
  /**　デストラクタ */
  ~ParseBC() {
    if (BaseBc) delete [] BaseBc;
  }
  
private:
  
  /**
   * @brief ラベルの重複を調べる
   * @param [in] n       テストするBaseBcの格納番号の最大値
   * @param [in] m_label テストラベル
   */
  bool chkDuplicate(const int n, const std::string m_label);
  
  
  /**
   * @brief ベクトルのコピー
   * @param [out] dst コピー先
   * @param [in]  src コピー元
   */
  void copyVec(REAL_TYPE* dst, const REAL_TYPE* src) 
  {
    dst[0] = src[0];
    dst[1] = src[1];
    dst[2] = src[2];
  }
  
  
  /**
   * @brief 境界条件の値(REAL_TYPE型)を取得し，返す
   * @param [in] label テストラベル
   */
  REAL_TYPE get_BCval_real(const std::string label);
  
  
  /**
   * @brief 外部境界条件のキーワードを照合し， BCの文字列を返す
   * @param [in] id 
   */
  std::string getOBCstr(const int id);
  

  int get_Vel_profile     (const std::string label_base);
  
  
  
  // コンポーネントのBbox情報st_xを返す
  int getCmpGbbox_st_x(const int odr, const int* gci) 
  {
    return ( gci[6*odr+0] );
  }
  
  
  // コンポーネントのBbox情報st_yを返す
  int getCmpGbbox_st_y(const int odr, const int* gci) 
  {
    return ( gci[6*odr+1] );
  }
  
  
  // コンポーネントのBbox情報st_zを返す
  int getCmpGbbox_st_z(const int odr, const int* gci) 
  {
    return ( gci[6*odr+2] );
  }
  
  // コンポーネントのBbox情報st_xを返す
  int getCmpGbbox_ed_x(const int odr, const int* gci) 
  {
    return ( gci[6*odr+3] );
  }
  
  
  // コンポーネントのBbox情報ed_yを返す
  int getCmpGbbox_ed_y(const int odr, const int* gci) 
  {
    return ( gci[6*odr+4] );
  }
  
  
  // コンポーネントのBbox情報ed_zを返す
  int getCmpGbbox_ed_z(const int odr, const int* gci) 
  {
    return ( gci[6*odr+5] );
  }
  
  
  /**
   * @brief 内部境界条件の座標値を取得し，登録する
   * @param [in]  label_base パラメータパス
   * @param [in]  n          オーダー
   * @param [out] v          ベクトルパラメータ
   */
  void get_Center(const std::string label_base, const int n, REAL_TYPE* v);
  
  
  
  void get_Dir            (const std::string label_base, const int n, REAL_TYPE* v);
  void get_NV             (const std::string label_base, const int n, REAL_TYPE* v);
  void get_IBC_Adiabatic  (const std::string label_base, const int n);
  void get_IBC_CnstTemp   (const std::string label_base, const int n);
  void get_IBC_Fan        (const std::string label_base, const int n);
  void get_IBC_IBM_DF     (const std::string label_base, const int n);
  void get_IBC_HeatFlux   (const std::string label_base, const int n);
  void get_IBC_HeatSrc    (const std::string label_base, const int n);
  void get_IBC_HT_B       (const std::string label_base, const int n);
  void get_IBC_HT_N       (const std::string label_base, const int n);
  void get_IBC_HT_S       (const std::string label_base, const int n);
  void get_IBC_HT_SF      (const std::string label_base, const int n);
  void get_IBC_HT_SN      (const std::string label_base, const int n);
  void get_IBC_IsoTherm   (const std::string label_base, const int n);
  void get_IBC_Monitor    (const std::string label_base, const int n, Control* C);
  void get_IBC_Outflow    (const std::string label_base, const int n);
  void get_IBC_Periodic   (const std::string label_base, const int n);
  void get_IBC_PrsLoss    (const std::string label_base, const int n);
  void get_IBC_Radiant    (const std::string label_base, const int n);
  void get_IBC_SpecVel    (const std::string label_base, const int n);
  void get_Darcy          (const std::string label_base, const int n);
  void get_OBC_FarField   (const std::string label_base, const int n);
  void get_OBC_HT         (const std::string label_base, const int n, const std::string kind);
  void get_OBC_Outflow    (const std::string label_base, const int n);
  void get_OBC_Periodic   (const std::string label_base, const int n);
  void get_OBC_SpecVH     (const std::string label_base, const int n);
  void get_OBC_Trcfree    (const std::string label_base, const int n);
  void get_OBC_Wall       (const std::string label_base, const int n);
  
  
  /**
   * @brief 単位ベクトルを計算して戻す
   * @param [in/out] v
   */
  void getUnitVec(REAL_TYPE* v);
  
  
  void get_Vel_Params     (const std::string label_base, const int prof, REAL_TYPE* ca, const char* str, const bool policy=false);
  
  
  /**
   * @brief コンポーネントが存在するかどうかを調べる
   * @param [in] label  テストするラベル
   * @retval bool値
   */
  bool isComponent(const int label);
  
  
  /**
   * @brief HTコンポーネントが存在するかどうかを調べる
   * @param [in] label  テストするラベル
   * @retval bool値
   */
  bool isCompoTransfer(const int label);
  
  
  /**
   * @brief 外部境界面の反対方向を返す
   * @param [in] dir 評価する方向
   * @return dirと反対方向
   */
  int oppositeDir(const int dir);
  
  
  /**
   * @brief 速度の外部境界条件処理の表示
   * @param [in] fp    ファイルポインタ
   * @param [in] ref   BoundaryOuter
   * @param [in] mat   MediumList
   * @param [in] G_reg グローバルの領域の大きさ
   * @param [in] face  面番号
   */
  void printOBC(FILE* fp, const BoundaryOuter* ref, const MediumList* mat, const REAL_TYPE* G_reg, const int face);
  
  
  /**
   * @brief 内部境界条件の照合を行う
   * @param [in] keyword テストキーワード
   * @param [in] m       BaseBcの格納番号
   */
  void setKeywordLBC(const std::string keyword, const int m);
  
  
  
  /**
   * @brief 内部境界条件の照合を行う
   * @param [in] keyword テストキーワード
   * @param [in] m       コンポーネントの格納番号
   */
  void setKeywordOBC(const std::string keyword, const int m);
  

  void set_Deface(const std::string label_base, const int n);
  
  
public:

  /**
   * @brief KOSと境界条件数の整合性をチェックする
   * @param [in] kos KindOfSolver
   */  
  void chkBCconsistency(const int kos);
  
  
  /**
   * @brief KOSと媒質の状態の整合性をチェックし，媒質数をカウント，C.NoMediumFluid, C.NoMediumSolidをセット
   * @param [in] Cref Controlクラス
   * @param [in] mat  MediumList
   */
  void countMedium(Control* Cref, const MediumList* mat);
  
  
  /**
   * @brief LocalBoundaryタグ直下のBCの個数（内部境界条件数）を返す
   */
  int getNoLocalBC();
  
  
  void get_Phase          ();
  void get_Medium_InitTemp();
  
  
  /**
   * @brief CMPのポインタを受け取る
   * @param [in] CMP  CompoListクラスのポインタ
   */
  void importCompoPtr(CompoList* CMP);
  
  
  /**
   * @brief TPのポインタを受け取る
   * @param [in] tp  TPControlクラスのポインタ
   */
  void importTP(TPControl* tp);
  
  
  /**
   * @brief 同じラベルが既にコンポーネントに登録されているかを調べる
   * @retval 重複していればfalseを返す
   * @param [in] candidate テストするラベル
   * @param [in] now       コンポーネントリストの現在までのエントリ番号
   */
  bool isLabelinCompo(const std::string candidate, const int now);
  
  
  
  /**
   * @brief CompoListに内部境界条件の情報を設定する
   * @param [in] C     Control
   * @param [in] mat   MediumList
   * @param [in] MTITP MediumTableInfo
   */
  void loadBC_Local(Control* C, const MediumList* mat, const MediumTableInfo *MTITP);
  
  
  /**
   * @brief パラメータファイルをパースして，外部境界条件を取得，保持する
   * @param [in/out] bc BoundaryOuter
   * @param [in] MTITP MediumTableInfo
   */
  void loadBC_Outer(BoundaryOuter* bc, const MediumTableInfo *MTITP);
  
  
  /**
   * @brief コンポーネントの情報を表示する
   * @param [in] fp  ファイルポインタ
   * @param [in] gci グローバルなコンポーネントのインデクス
   * @param [in] mat MediumList
   * @param [in] cmp CompoList
   * @param [in] bc  BoundaryOuter
   */
  void printCompo(FILE* fp, const int* gci, const MediumList* mat, CompoList* cmp, const BoundaryOuter* bc);
  
  
  /**
   * @brief 外部境界条件の各面の情報を表示する
   * @param [in] fp    ファイルポインタ
   * @param [in] G_reg グローバルの領域の大きさ
   * @param [in] bc    BoundaryOuter
   * @param [in] mat   MediumList
   */
  void printFaceOBC(FILE* fp, const REAL_TYPE* G_reg, const BoundaryOuter* bc, const MediumList* mat);
  
  
  void setControlVars(Control* Cref);
  
  
  /**
   * @brief 指定した媒質IDから参照物理量を設定する
   * @param [in] mat MediumList
   * @param [in] Ref 参照媒質番号
   */
  void setRefMedium(MediumList* mat, const int Ref);
  
  
  void setRefValue(MediumList* mat, CompoList* cmp, Control* C);
  
};

#endif // _FB_PARA_BC_H_
