#ifndef _FB_MONITOR_H_
#define _FB_MONITOR_H_

// #################################################################
//
// CAERU Library
//
// Copyright (c) All right reserved. 2012
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
//
// #################################################################

//@file   Monitor.h
//@brief  FlowBase MonitorList class Header
//@author kero

#include <string>
#include <vector>

#include "DomainInfo.h"
#include "FB_Define.h"
#include "vec3.h"
#include "MonCompo.h"
#include "basic_func.h"
#include "Component.h"
#include "TPControl.h"
#include "Control.h"

using namespace std;

/**
 *  モニタグループ管理クラス.
 */
class MonitorList : public DomainInfo {
  
public:
  /// 出力タイプ型
  enum OutputType {
    GATHER,      ///< 単一ファイル出力
    DISTRIBUTE,  ///< 分散出力
    NONE,
  };
  
protected:
  OutputType outputType;  ///< 出力タイプ
  
  int nGroup;    ///< モニタリンググループ数
  vector<MonitorCompo*> monGroup;  ///< モニタリンググループ配列
  
  Vec3r org;         ///< ローカル基点座標
  Vec3r pch;         ///< セル幅
  Vec3r box;         ///< ローカル領域サイズ
  
  Vec3r g_org;       ///< グローバル基点座標
  Vec3r g_box;       ///< グローバル領域サイズ
  
  MonitorCompo:: ReferenceVariables refVar;  ///< 参照用パラメータ変数
  
  int* bcd;     ///< BCindex ID
  int num_process;
  
  TPControl* tpCntl;   ///< テキストパーサへのポインタ
  
public:
  /// コンストラクタ
  MonitorList() {
    nGroup = 0;
    outputType = NONE;
  }
  
  /// デストラクタ
  ~MonitorList() { for (int i = 0; i < nGroup; i++) delete monGroup[i]; }
  
  
  /**
   @brief TPに記述されたモニタ座標情報を取得し，リストに保持する
   @param [in] C  Control クラスオブジェクトのポインタ
   */
  void get_Monitor(Control* C); //未完成
  
  /**
   @brief TPに記述されたモニタ座標情報(Line)を取得
   @param C  Control クラスオブジェクトのポインタ
   @param from Line始点座標
   @param to   Line終点座標
   @param nDivision Line分割数
   @note データは無次元化して保持
   */
  void get_Mon_Line(
                    Control* C,
                    const string label_base,
                    REAL_TYPE from[3],
                    REAL_TYPE to[3],
                    int& nDivision);
  
  /**
   @brief TPに記述されたモニタ座標情報を取得(PointSet)
   @param C  Control クラスオブジェクトのポインタ
   @param pointSet PointSet配列
   @note データは無次元化して保持
   */
  void get_Mon_Pointset(//未修正
                        Control* C,
                        const string label_base,
                        vector<MonitorCompo::MonitorPoint>& pointSet);
  
  
  /**
   * @brief TPのポインタを受け取る
   * @param [in] tp TPControl
   */
  void importTP(TPControl* tp);
  
  
  /// PointSet登録.
  ///
  ///   @param[in] str ラベル文字列
  ///   @param[in] variables モニタ変数vector
  ///   @param[in] method method文字列
  ///   @param[in] mode   mode文字列
  ///   @param[in] pointSet  PointSet
  ///
  void setPointSet(const char* labelStr, vector<string>& variables,
                   const char* methodStr, const char* modeStr,
                   vector<MonitorCompo::MonitorPoint>& pointSet);
  
  /// Line登録.
  ///
  ///   @param[in] str ラベル文字列
  ///   @param[in] variables モニタ変数vector
  ///   @param[in] method method文字列
  ///   @param[in] mode   mode文字列
  ///   @param[in] from Line始点
  ///   @param[in] to   Line終点
  ///   @param[in] nDivision 分割数(モニタ点数-1)
  ///
  void setLine(const char* labelStr, vector<string>& variables,
               const char* methodStr, const char* modeStr,
               REAL_TYPE from[3], REAL_TYPE to[3], int nDivision);
  
  /// 内部境界条件としてモニタ点を登録.
  ///
  ///   @param[in] cmp コンポーネント配列
  ///   @param[in] nBC コンポーネント数
  ///
  void setInnerBoundary(CompoList* cmp, int nBC);
  
  /// モニタ情報を出力.
  ///   @param[in] str 出力ファイルの基本名
  void printMonitorInfo(FILE* fp, const char* str, const bool verbose);
  
  /// 必要なパラメータのコピー.
  ///
  ///   @param[in] bcd BCindex ID
  ///   @param[in] refVelocity    代表速度
  ///   @param[in] baseTemp       基準温度
  ///   @param[in] diffTemp       代表温度差
  ///   @param[in] refDensity     基準密度
  ///   @param[in] refLength      代表長さ
  ///   @param[in] basePrs        基準圧力
  ///   @param[in] unitTemp       温度単位指定フラグ (Kelvin / Celsius)
  ///   @param[in] modePrecision  出力精度指定フラグ (単精度，倍精度)
  ///   @param[in] unitPrs        圧力単位指定フラグ (絶対値，ゲージ圧)
  ///
  void setControlVars(int* bcd,
                      REAL_TYPE refVelocity, REAL_TYPE baseTemp, REAL_TYPE diffTemp,
                      REAL_TYPE refDensity, REAL_TYPE refLength, REAL_TYPE basePrs,
                      int unitTemp, int modePrecision, int unitPrs,int num_process);
  
  /// 参照速度のコピー
  ///   @param[in] v00 座標系移動速度
  void set_V00(REAL_TYPE v00[4]);
  
  /// サンプリング元データの登録.
  ///
  ///   @param[in] v 速度変数配列
  ///   @param[in] p 圧力変数配列
  ///   @param[in] t 温度変数配列
  ///
  void setDataPtrs(REAL_TYPE* v, REAL_TYPE* p, REAL_TYPE* t=NULL) {
    for (int i = 0; i < nGroup; i++) monGroup[i]->setDataPtrs(v, p, t);
  }
  
  /// 出力タイプの設定.
  void setOutputType(OutputType type) { outputType = type; }
  
  /// サンプリングと出力の次元の設定
  /// @param[in] modeUnit       出力単位指定フラグ (有次元，無次元)
  void setSamplingUnit(int unit) { refVar.modeUnit = unit; }
  
  /// 出力ファイルオープン.
  ///
  ///    @param str ファイル名テンプレート
  ///
  ///    @note 出力タイプによらず全プロセスから呼んでも問題ない
  ///
  void openFile(const char* str);
  
  /// 出力ファイルクローズ.
  ///
  ///    @note 出力タイプによらず全プロセスから呼んでも問題ない
  ///
  void closeFile();
  
  /// サンプリング(Line, PointSet指定).
  void sampling() {
    for (int i = 0; i < nGroup; i++) {
      if (monGroup[i]->getType() != MonitorCompo::INNER_BOUNDARY) {
        monGroup[i]->sampling();
      }
    }
  }
  
  /// サンプリング(内部境界条件指定).
  ///
  ///   サンプリング結果を集計，コンポーネント領域での平均値を計算.
  ///   速度は法線ベクトルとの内積をとる．
  ///   結果はコンポーネントcmpに格納.
  ///
  void samplingInnerBoundary() {
    for (int i = 0; i < nGroup; i++) {
      if (monGroup[i]->getType() == MonitorCompo::INNER_BOUNDARY) {
        monGroup[i]->samplingInnerBoundary();
      }
    }
  }
  
  /// モニタ結果出力(Line, PointSet指定).
  ///
  ///   @param[in] step サンプリング時の計算ステップ
  ///   @param[in] tm   サンプリング時の計算時刻
  ///
  ///   @note gatherの場合も全プロセスから呼ぶこと
  ///
  void print(unsigned step, REAL_TYPE tm);
  
  /// 指定されるモニタ点位置にID=255を書き込む.
  ///
  ///   @param[in] id セルID配列
  ///
  void write_ID(int* id);
  
protected:
  
  /// Line指定の端点座標をグローバル計算領域内にクリッピング.
  ///
  ///   @param[in,out] from Line始点
  ///   @param[in,out] to   Line終点
  ///
  void clipLine(REAL_TYPE from[3], REAL_TYPE to[3]);
  
  /// 出力タイプ文字列の取得.
  string getOutputTypeStr();
  

  //@brief 座標値を無次元化する
  void normalizeCord(REAL_TYPE RefLength, REAL_TYPE x[3]) {
    x[0] /= RefLength;
    x[1] /= RefLength;
    x[2] /= RefLength;
  }
  
};

#endif // _FB_MONITOR_H_
