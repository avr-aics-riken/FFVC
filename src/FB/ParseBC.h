#ifndef _FB_PARA_BC_H_
#define _FB_PARA_BC_H_

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

//@file   ParseBC.h
//@brief  FlowBase ParseBC class Header
//@author kero

#include "string.h"

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

using namespace std;


class ParseBC : public DomainInfo {
private:
  
  int NoBC;        ///< LocalBCの数
  int NoBaseBC;    ///< 入力ファイルに記載された外部境界条件の指定種類数
  int NoCompo;     ///< コンポーネントリストの全登録数
  int NoMedium;    ///< 媒質数
  int Unit_Temp;   ///< 温度単位
  int Unit_Prs;    ///< 圧力単位
  int globalBC[6]; ///< グローバルな外部境界条件を保持するBasicBc[]への格納番号ポインタ
  
  REAL_TYPE RefVelocity, BaseTemp, DiffTemp, RefDensity, RefSpecificHeat;
  REAL_TYPE RefLength, BasePrs;
  REAL_TYPE rho, nyu, cp, lambda, beta; // 無次元化の参照値

  int KindOfSolver; ///< 支配方程式の種類
  int Unit_Param;
  int monitor;
  int Mode_Gradp;
  bool HeatProblem;
  bool isCDS;
  
  TextParser* tpCntl;      ///< テキストパーサーのラッパークラス
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
    
    tpCntl = NULL;
    BaseBc = NULL;
  }
  
  /**　デストラクタ */
  ~ParseBC() {
    if (BaseBc) delete [] BaseBc;
  }
  
private:
  
  // ラベルの重複を調べる
  bool chkDuplicate(const int n, const string m_label);
  
  
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
   * @brief コンポーネントが存在するかどうかを調べる
   * @param [in] label  テストするラベル
   * @param [in] cmp    CompoList
   * @retval bool値
   */
  bool existComponent(const int label, const CompoList* cmp) const
  {
    for (int n=1; n<=NoCompo; n++)
    {
      if ( cmp[n].getType() == label ) return true;
    }
    return false;
  }
  
  
  /**
   * @brief HTコンポーネントが存在するかどうかを調べる
   * @param [in] label  テストするラベル
   * @param [in] cmp    CompoList
   * @retval bool値
   */
  bool existCompoTransfer(const int label, const CompoList* cmp) const
  {
    for (int n=1; n<=NoCompo; n++)
    {
      if ( cmp[n].getHtype() == label ) return true;
    }
    return false;
  }
  
  
  // 境界条件の値(REAL_TYPE型)を取得し，返す
  REAL_TYPE get_BCval_real(const string label);
  
  
  // 外部境界条件のキーワードを照合し， BCの文字列を返す
  string getOBCstr(const int id);
  

  // 外部境界の速度境界条件のタイプを取得し，返す
  int getVprofile(const string label_base);
  
  
  // コンポーネントのBbox始点情報を返す
  FB::Vec3i getCmpGbbox_st(const int odr, const int* gci) const
  {
    FB::Vec3i st;
    st.x = gci[6*odr+0];
    st.y = gci[6*odr+1];
    st.z = gci[6*odr+2];
    return st;
  }
  
  /*
   st_x ; gci[6*m+0];
   st_y ; gci[6*m+1];
   st_z ; gci[6*m+2];
   ed_x ; gci[6*m+3];
   ed_y ; gci[6*m+4];
   ed_z ; gci[6*m+5];
   */
  
  // コンポーネントのBbox終点情報を返す
  FB::Vec3i getCmpGbbox_ed(const int odr, const int* gci) const
  {
    FB::Vec3i ed;
    ed.x = gci[6*odr+3];
    ed.y = gci[6*odr+4];
    ed.z = gci[6*odr+5];
    return ed;
  }
  
  
  // 内部境界条件の座標値を取得し，登録する
  void get_Center(const string label_base, REAL_TYPE* v);
  
  
  // 内部境界条件の方向ベクトル値を取得し，登録する
  void get_Dir(const string label_base, REAL_TYPE* v);
  
  
  // 内部境界条件の法線ベクトル値を取得し，登録する
  void get_NV(const string label_base, REAL_TYPE* v);
  
  
  // Const_Temperatureのパラメータを取得する
  void get_IBC_CnstTemp(const string label_base, const int n, CompoList* cmp);
  
  
  // Fanのパラメータを取得する
  void get_IBC_Fan(const string label_base, const int n, CompoList* cmp);
  
  
  // Direct Forcingのパラメータを取得する
  void get_IBC_IBM_DF(const string label_base, const int n, CompoList* cmp);
  
  
  // Direct_Fluxのパラメータを取得する
  void get_IBC_HeatFlux(const string label_base, const int n, CompoList* cmp);
  
  
  // Heat_Generationのパラメータを取得する
  void get_IBC_HeatSrc(const string label_base, const int n, CompoList* cmp);
  
  
  // HeatTransfer_Bのパラメータを取得する
  void get_IBC_HT_B(const string label_base, const int n, CompoList* cmp);
  
  
  // HeatTransfer_Nのパラメータを取得する
  void get_IBC_HT_N(const string label_base, const int n, CompoList* cmp);
  
  
  // HeatTransfer_Sのパラメータを取得
  void get_IBC_HT_S(const string label_base, const int n, CompoList* cmp, const MediumList* mat);
  
  
  // HeatTransfer_SFのパラメータを取得する
  void get_IBC_HT_SF(const string label_base, const int n, CompoList* cmp);
  
  
  // HeatTransfer_SNのパラメータを取得する
  void get_IBC_HT_SN(const string label_base, const int n, CompoList* cmp, const MediumList* mat);
  
  
  // 境界条件IsoThermalのパラメータを取得し保持する
  void get_IBC_IsoTherm(const string label_base, const int n, CompoList* cmp);
  
  
  // Monitorの設定内容をパースし，パラメータを保持する
  void get_IBC_Monitor(const string label_base, const int n, CompoList* cmp);
  
  
  // 内部の流出境界のパラメータを取得する
  void get_IBC_Outflow(const string label_base, const int n, CompoList* cmp);
  
  
  // 内部の周期境界のパラメータを取得する
  void get_IBC_Periodic(const string label_base, const int n, CompoList* cmp);
  
  
  // HeatExchangerのパラメータを取得する
  void get_IBC_PrsLoss(const string label_base, const int n, CompoList* cmp);
  
  
  // 境界条件Radiantのパラメータを取得し保持する
  void get_IBC_Radiant(const string label_base, const int n, CompoList* cmp);
  
  
  // 内部の流入境界のパラメータを取得する
  void get_IBC_SpecVel(const string label_base, const int n, CompoList* cmp);
  
  
  // Darcyのパラメータを取得する
  void get_Darcy(const string label_base, const int n, CompoList* cmp);
  
  
  // 外部境界の遠方境界のパラメータを取得する
  void getObcFarField (const string label_base, const int n);
  
  
  // 外部の壁面熱伝達境界のパラメータを取得する
  void getObcHeatTransfer(const string label_base, const int n, const string kind);
  
  
  // 外部境界の流出条件のパラメータを取得する
  void getObcOutflow(const string label_base, const int n);
  
  
  // 外部境界の周期条件のパラメータを取得する
  void getObcPeriodic(const string label_base, const int n);
  
  
  // 外部境界の流入条件のパラメータを取得する
  void getObcSpecVH(const string label_base, const int n);
  
  
  // トラクションフリーの外部境界のパラメータを取得する
  void getObcTrcfree(const string label_base, const int n);
  
  
  // 外部境界の壁面条件のパラメータを取得する
  void getObcWall(const string label_base, const int n);
  
  
  // 隣接セル媒質のパラメータを取得する
  void get_Neighbor(const string label_base, const int n, CompoList* cmp, const MediumList* mat);
  
  
  // 単位ベクトルを計算して戻す
  void getUnitVec(REAL_TYPE* v);
  
  
  // 速度のパラメータを取得する
  void getVelocity(const string label_base, const int prof, REAL_TYPE* ca, const char* str, const bool policy=false);

  
  
  // 外部境界面の反対方向を返す
  int oppositeDir(const int dir);
  
  
  // 速度の外部境界条件処理の表示
  void printOBC(FILE* fp, const BoundaryOuter* ref, const MediumList* mat, const REAL_TYPE* G_reg, const int face);
  
  
  // 内部境界条件の照合を行う
  void setKeywordLBC(const string keyword, const int m, CompoList* cmp);
  
  
  // 外部境界条件のキーワードを照合し，コードを登録する
  void setKeywordOBC(const string keyword, const int m);
  
  
  
public:

  /**
   * @brief テーブルのチェック（デバッグ用）
   * @param [in] mat  MediumList
   * @param [out] cmp   CompoList
   */
  void checkList(const MediumList* mat, const CompoList* cmp);
  
  /**
   * @brief KOSと境界条件数の整合性をチェックする
   * @param [in] kos KindOfSolver
   * @param [in] cmp CompListクラスのポインタ
   */
  void chkBCconsistency(const int kos, const CompoList* cmp);
  
  
  /**
   * @brief KOSと媒質の状態の整合性をチェックし，媒質数をカウント，C.NoMediumFluid, C.NoMediumSolidをセット
   * @param [in] Cref Controlクラス
   * @param [in] mat  MediumList
   */
  void countMedium(Control* Cref, const MediumList* mat);
  
  
  /**
   * @brief 温度計算の場合の各媒質の初期値を取得する
   * @param [in,out] cmp    CompoList
   */
  void getInitTempOfMedium(CompoList* cmp);
  
  
  /**
   * @brief 2相流問題で気相か液相かを取得する
   * @param [out] cmp   CompoList
   */
  void get_Phase(CompoList* cmp);
  
  
  /**
   * @brief TPのポインタを受け取る
   * @param [in] tp  TextParserクラスのポインタ
   */
  void importTP(TextParser* tp);

  
  
  /**
   * @brief CompoListに内部境界条件の情報を設定する
   * @param [in]      C     Control
   * @param [in,out]  mat   MediumList
   * @param [out]     cmp   CompoList
   * @param [in]      PP    ポリゴン属性管理クラス [0]-[n-1]
   * @note 格納番号は1からスタート
   */
  void loadLocalBC(Control* C, MediumList* mat, CompoList* cmp, PolygonProperty* polyPP);
  
  
  /**
   * @brief パラメータファイルをパースして，外部境界条件を取得，保持する
   * @param [in,out] bc   BoundaryOuter
   * @param [in]     mat  MediumList
   * @param [out]    cmp  CompoList
   */
  void loadOuterBC(BoundaryOuter* bc, const MediumList* mat, CompoList* cmp);
  
  
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
   * @brief 取得したCompoList[]の内容を表示する
   * @param [in] fp      ファイルポインタ
   * @param [in] compo   CompoList
   * @param [in] basicEq 基礎方程式の種類
   * @note Hostonly
   */
  void printCompoSummary(FILE* fp, CompoList* cmp, const int basicEq);
  
  
  /**
   * @brief 外部境界条件の各面の情報を表示する
   * @param [in] fp    ファイルポインタ
   * @param [in] G_reg グローバルの領域の大きさ
   * @param [in] bc    BoundaryOuter
   * @param [in] mat   MediumList
   */
  void printFaceOBC(FILE* fp, const REAL_TYPE* G_reg, const BoundaryOuter* bc, const MediumList* mat);
  
  
  /**
   * @brief 変数の初期化
   * param [in] Cref Control class
   */
  void setControlVars(Control* Cref);
  
  
  /**
   * @brief 指定した媒質IDから参照物理量を設定する
   * @param [in] mat MediumList
   * @param [in] Ref 参照媒質の格納順
   */
  void setRefMediumProperty(const MediumList* mat, const int Ref);
  
};

#endif // _FB_PARA_BC_H_
