#ifndef _FB_MONITOR_COMPO_H_
#define _FB_MONITOR_COMPO_H_

// #################################################################
//
// FFV : Frontflow / violet
//
// Copyright (c) All right reserved. 2012
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
//
// #################################################################

/** 
 * @file   MonCompo.h
 * @brief  FlowBase MonitorCompo class Header
 * @author kero
 */

#include <string>
#include <vector>
#include <cassert>

#include "DomainInfo.h"
#include "FB_Define.h"
#include "vec3.h"
#include "basic_func.h"
#include "Component.h"
#include "Sampling.h"

using namespace std;


/**
 *  モニタグループクラス.
 */
class MonitorCompo : public DomainInfo {
public:
  
  /// PointSet要素用構造体
  struct MonitorPoint {
    Vec3r crd;      ///< モニタ点座標
    string label;   ///< モニタ点ラベル
    MonitorPoint(const REAL_TYPE v[3], const char* str) : crd(v), label(str) {}
    ~MonitorPoint() {}
  };
  
  /// 参照用パラメータ構造体
  struct ReferenceVariables {
    int modeUnit;       /// 出力単位指定フラグ (有次元，無次元)
    int unitTemp;       /// 温度単位指定フラグ (Kelvin / Celsius)
    int modePrecision;  /// 出力精度指定フラグ (単精度，倍精度)
    int unitPrs;        /// 圧力単位指定フラグ (絶対値，ゲージ圧)
    REAL_TYPE refVelocity;   /// 代表速度
    REAL_TYPE baseTemp;      /// 基準温度
    REAL_TYPE diffTemp;      /// 代表温度差
    REAL_TYPE refDensity;    /// 基準密度
    REAL_TYPE refLength;     /// 代表長さ
    REAL_TYPE basePrs;       /// 基準圧力
    Vec3r v00;               /// 参照（座標系移動）速度
  };
  
  /// モニタ点指定タイプ型
  enum Type {
    POINT_SET,
    LINE,
    INNER_BOUNDARY,
  };
  
  /// モニタ変数識別定数
  enum Var {
    VELOCITY,
    PRESSURE,
    TEMPERATURE,
    TOTAL_PRESSURE,
    VORTICITY,
    NUM_VAR,
  };
  
  
protected:
  int nPoint;             ///< モニタ点数
  bool variable[NUM_VAR]; ///< モニタ変数フラグ
  string label;           ///< グループのラベル
  Type type;              ///< モニタ点指定タイプ
  int method;             ///< サンプリング方法
  int mode;               ///< サンプリングモード
  
  Sampling** mon;    ///< 「モニタ点毎のSampligクラスへのポインタ」の配列
  
  Vec3r org;         ///< ローカル基点座標
  Vec3r pch;         ///< セル幅
  Vec3r box;         ///< ローカル領域サイズ
  
  Vec3r g_org;       ///< グローバル基点座標
  Vec3r g_box;       ///< グローバル領域サイズ
  
  ReferenceVariables refVar;  ///< 参照用パラメータ変数
  
  int* bcd;     ///< BCindex ID
  
  FILE* fp;          ///< 出力ファイルポインタ
  
  // サンプリング元データ
  REAL_TYPE* vSource;  ///< 速度サンプリング元データ
  REAL_TYPE* pSource;  ///< 圧力サンプリング元データ
  REAL_TYPE* tSource;  ///< 温度サンプリング元データ
  
  
  Vec3r* crd;      ///< モニタ点座標配列
  int* rank;       ///< モニタ点担当ランク番号配列
  string* comment; ///< モニタ点コメント配列
  int* pointStatus;///< 不正モニタ点フラグ配列
  
  Vec3r* vel;      ///< 速度サンプリング結果配列
  REAL_TYPE* prs;  ///< 圧力サンプリング結果配列
  REAL_TYPE* tmp;  ///< 温度サンプリング結果配列
  REAL_TYPE* tp;   ///< 全圧サンプリング結果配列
  Vec3r* vor;      ///< 渦度サンプリング結果配列
  
  CompoList* cmp;  ///< 内部境界条件指定の場合の対応するコンポーネントへのポインタ
  
  int num_process;
  
public:
  /// ディフォルトコンストラクタ.
  MonitorCompo() {
    nPoint = 0;
    for (int i = 0; i < NUM_VAR; i++) variable[i] = false;
    vel = vor  = NULL;
    prs = tmp = tp  = NULL;
    crd = NULL;
    rank = NULL;
    comment = NULL;
    pointStatus = NULL;
    vSource = pSource = tSource = NULL;
  }
  
  /// コンストラクタ.
  ///
  ///   param[in] pn Parallel_Info
  ///   param[in] org,pch,box ローカル領域基点座標，セル幅，領域サイズ
  ///   param[in] g_org,g_box グローバル領域基点座標，領域サイズ
  ///   param[in] size,guide  ローカルセルサイズ，ガイドセル数
  ///   param[in] refVar  参照パラメータ
  ///   param[in] bcd  BCindex ID
  ///
  MonitorCompo(Vec3r org, Vec3r pch, Vec3r box, Vec3r g_org, Vec3r g_box, 
               ReferenceVariables refVar,
               int* bcd, int num_process) {
    nPoint = 0;
    for (int i = 0; i < NUM_VAR; i++) variable[i] = false;
    vel = vor  = NULL;
    prs = tmp = tp  = NULL;
    crd = NULL;
    rank = NULL;
    comment = NULL;
    pointStatus = NULL;
    vSource = pSource = tSource = NULL;
    
    this->org = org;
    this->pch = pch;
    this->box = box;
    this->g_org = g_org;
    this->g_box = g_box;
    this->refVar = refVar;
    this->bcd = bcd;
    this->num_process = num_process;
  }
  
  /// デストラクタ.
  ~MonitorCompo() {
    if (crd) delete[] crd;
    if (vel) delete[] vel;
    if (prs) delete[] prs;
    if (tmp) delete[] tmp;
    if (tp) delete[] tp;
    if (rank) delete[] rank;
    if (comment) delete[] comment;
    if (pointStatus) delete[] pointStatus;
    
    if (mon) {
      for (int i = 0; i < nPoint; i++) if (mon[i]) delete mon[i];
      delete[] mon;
    }
  }
  
  /// モニタ点が指定領域内にあるかを判定.
  ///
  ///   @param[in] m モニタ点番号
  ///   @param[in] org 調査領域の基点
  ///   @param[in] box 調査領域のサイズ
  ///   @param[in] flag メッセージ出力フラグ(trueの時出力)
  ///   @return true=領域内/false=領域外
  ///
  bool check_region(int m, Vec3r org, Vec3r box, bool flag=false);
  
  
  /// PointSet登録.
  ///
  ///   @param[in] labelStr ラベル文字列
  ///   @param[in] variables モニタ変数vector
  ///   @param[in] methodStr method文字列
  ///   @param[in] modeStr   mode文字列
  ///   @param[in] pointSet  PointSet
  ///
  void setPointSet(const char* lableStr, vector<string>& variables, 
                   const char* methodStr, const char* modeStr,
                   vector<MonitorPoint>& pointSet);
  
  /// Line登録.
  ///
  ///   @param[in] labelStr ラベル文字列
  ///   @param[in] variables モニタ変数vector
  ///   @param[in] methodStr method文字列
  ///   @param[in] modeStr   mode文字列
  ///   @param[in] from Line始点
  ///   @param[in] to   Line終点
  ///   @param[in] nDivision 分割数(モニタ点数-1)
  ///
  void setLine(const char* labelStr, vector<string>& variables,
               const char* methodStr, const char* modeStr,
               REAL_TYPE from[3], REAL_TYPE to[3], int nDivision);
  
  /// 内部境界条件としてモニタ点を登録.
  ///
  ///   @param[in] n コンポーネントエントリ
  ///   @param[in] cmp コンポーネント
  ///
  void setInnerBoundary(int n, CompoList& cmp);
  
  /// サンプリング元データの登録.
  ///
  ///   @param[in] v 速度変数配列
  ///   @param[in] p 圧力変数配列
  ///   @param[in] t 温度変数配列
  ///
  void setDataPtrs(REAL_TYPE* v, REAL_TYPE* p, REAL_TYPE* t) {
    vSource = v;
    pSource = p;
    tSource = t;
    if (variable[TEMPERATURE] && !tSource) {
      Hostonly_ stamped_printf("\tError : Temperature monitoring source not assigned.\n");
      Exit(0);
    }
  }
  
  /// サンプリング(Line, PointSet).
  void sampling();
  
  /// 内部境界条件モニタ点でのサンプリング.
  ///
  ///   サンプリング結果を集計，コンポーネント領域での平均値を計算.
  ///   速度は法線ベクトルとの内積をとる．
  ///   結果はコンポーネントcmpに格納.
  ///
  void samplingInnerBoundary();
  
  /// モニタ結果出力(gather).
  ///
  ///   @param[in] step サンプリング時の計算ステップ
  ///   @param[in] tm   サンプリング時の計算時刻
  ///
  ///   @note 内部でノード0に集計するため全プロセスから呼ぶこと
  ///
  void print_gather(unsigned step, REAL_TYPE tm) { 
    gatherSampled();
    if (myRank == 0) print(step, tm, true);
  }
  
  /// モニタ結果出力(distribute).
  ///
  ///   @param[in] step サンプリング時の計算ステップ
  ///   @param[in] tm   サンプリング時の計算時刻
  ///
  ///   @note 全プロセスから呼ぶこと
  ///
  void print_distribute(unsigned step, REAL_TYPE tm) {
    print(step, tm, false);
  }
  
  /// グループラベルを返す.
  const char* getLabel() { return label.c_str(); }
  
  /// 登録タイプを返す.
  Type getType() { return type; }
  
  /// m番目のモニタ点を含むセルインデクスを返す.
  Vec3i getSamplingCellIndex(int m) {
    Vec3i index;
    Vec3r c = (crd[m] - org) / pch;
    index.x = int(c.x) + 1;
    index.y = int(c.y) + 1;
    index.z = int(c.z) + 1;
    return index;
  }
  
  /// モニタ点数を返す.
  int getSize() { return nPoint; }
  
  /// 出力ファイルオープン.
  ///
  ///    @param str ファイル名テンプレート
  ///    @param gathered true=gather出力/false=disutribute出力
  ///
  void openFile(const char* str, bool gathered);
  
  /// 出力ファイルクローズ.
  void closeFile();
  
  /// モニタ情報を出力.
  ///
  ///    @param[in] fp 出力ファイルポインタ
  ///    @param[in] no モニタグループ通し番号
  //
  void printInfo(FILE* fp, int no);
  
protected:
  
  /// サンプリングした変数をノード0に集約.
  void gatherSampled();
  
  /// サンプリングしたスカラー変数をノード0に集約.
  ///
  ///   @param[in,out] s スカラー変数配列
  ///   @param  sRecvBuf  通信用work領域
  ///
  void gatherSampledScalar(REAL_TYPE* s, REAL_TYPE* sRecvBuf);
  
  /// サンプリングしたベクトル変数をノード0に集約.
  ///
  ///   @param[in,out] v ベクトル変数配列
  ///   @param  vSendBuf,vRecvBuf  通信用work領域
  ///
  void gatherSampledVector(Vec3r* v, REAL_TYPE* vSendBuf, REAL_TYPE* vRecvBuf);
  
  /// モニタ対象物理量の設定.
  ///
  ///    @param[in] str モニタ対象物理量文字列
  ///
  void setMonitorVar(const char* str);
  
  /// サンプリング方法の設定.
  ///
  ///    @param[in] str サンプリング方法文字列
  ///
  void setSamplingMethod(const char* str);
  
  /// サンプリングモードの設定.
  ///
  ///    @param[in] str サンプリングモード文字列
  ///
  void setSamplingMode(const char* str);
  
  /// モニタリング管理用配列の確保.
  void allocArray();
  
  /// サンプリング値を格納する配列の確保.
  void allocSamplingArray();
  
  /// サンプリング方法文字列の取得.
  string getMethodStr();
  
  /// サンプリングモード文字列の取得.
  string getModeStr();
  
  /// モニタ点指定方法文字列の取得.
  string getTypeStr();
  
  /// モニタ変数を結合した文字列の取得.
  string getVarStr();
  
  /// モニタ点の状態を調べ，不正モニタ点フラグ配列pointStatusを設定.
  void checkMonitorPoints();
  
  /// 内部境界条件として指定されたモニタ領域のセル中心座標をcrd[]に設定.
  ///
  ///   @param[in] n コンポーネントエントリ
  ///   @param[in] cmp コンポーネント
  ///
  void setIBPoints(int n, CompoList& cmp);
  
  /// 内部境界条件として指定されたモニタ領域内でスカラー変数を平均.
  ///
  ///   @param[in] s スカラー変数配列
  ///   @return モニタ領域内平均値
  ///
  REAL_TYPE averageScalar(REAL_TYPE* s);
  
  /// 内部境界条件として指定されたモニタ領域内でベクトル変数を平均.
  ///
  ///   @param[in] v ベクトル変数配列
  ///   @return モニタ領域内平均値
  ///
  Vec3r averageVector(Vec3r* v);
  
  
  /// 各モニタ点を担当するランク番号を配列rank[]にセット.
  ///
  ///   @note 領域境界上のモニタ点は，ランク番号の大きい方の領域が担当
  ///
  void setRankArray();
  
  /// モニタ結果出力.
  ///
  ///   @param[in] step サンプリング時の計算ステップ
  ///   @param[in] tm   サンプリング時の計算時刻
  ///   @param[in] gathered 出力モードフラグ(true=gather出力/false=disutribute出力)
  ///
  void print(unsigned step, REAL_TYPE tm, bool gathered);
  
  /// モニタ結果出力ファイルにヘッダ部を出力.
  ///
  ///   @param[in] gathered 出力モードフラグ(true=gather出力/false=disutribute出力)
  ///
  void writeHeader(bool gathered);
  
  /// >>> 出力変換 --------------------------------
  /// 座標の単位変換
  REAL_TYPE convCrd(REAL_TYPE xyz) const {
    return ( (refVar.modeUnit==DIMENSIONAL) ? xyz*refVar.refLength : xyz );
  }
  
  /// 時間の単位変換
  REAL_TYPE convTime(REAL_TYPE tm) const {
    return ( (refVar.modeUnit==DIMENSIONAL) ? tm*refVar.refLength/refVar.refVelocity : tm );
  }
  
  /// 速度成分の単位変換
  REAL_TYPE convVel(REAL_TYPE vel) const {
    return ( (refVar.modeUnit==DIMENSIONAL) ? vel*refVar.refVelocity : vel ); 
  }
  
  /// 圧力の単位変換
  REAL_TYPE convPrs(REAL_TYPE prs) const {
    return ( (refVar.modeUnit==DIMENSIONAL) ? FBUtility::convND2D_P(prs, refVar.basePrs, refVar.refDensity, refVar.refVelocity, refVar.unitPrs) : prs );
  }
  
  /// 温度の単位変換
  REAL_TYPE convTmp(REAL_TYPE tmp) const {
    return ( (refVar.modeUnit==DIMENSIONAL) ? FBUtility::convK2Temp(FBUtility::convND2Kelvin(tmp, refVar.baseTemp, refVar.diffTemp), refVar.unitTemp) : tmp);
  }
  
  /// 全圧の単位変換
  REAL_TYPE convTP(REAL_TYPE tp) const {
    return ( (refVar.modeUnit==DIMENSIONAL) ? tp * (refVar.refVelocity * refVar.refVelocity * refVar.refDensity) : tp );
  }
  
  /// 渦度成分の単位変換
  REAL_TYPE convVor(REAL_TYPE vor) const {
    return ( (refVar.modeUnit==DIMENSIONAL) ? vor*refVar.refVelocity/refVar.refLength : vor ); 
  }
  
  /// >>> 集約 --------------------------------
  /// Allreduceによる総和(実数配列上書き，work配列指定).
  bool allReduceSum(REAL_TYPE* array, int n, REAL_TYPE* sendBuf);
  
  /// Allreduceによる総和(実数配列上書き).
  bool allReduceSum(REAL_TYPE* array, int n);
  
  /// Allreduceによる総和(整数配列上書き，work配列指定).
  bool allReduceSum(int* array, int n, int* sendBuf);
  
  /// Allreduceによる総和(整数配列上書き).
  bool allReduceSum(int* array, int n);
  
};

#endif // _FB_MONITOR_COMPO_H_
