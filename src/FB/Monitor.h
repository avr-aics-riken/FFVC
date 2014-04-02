#ifndef _FB_MONITOR_H_
#define _FB_MONITOR_H_

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

//@file   Monitor.h
//@brief  FlowBase MonitorList class Header
//@author aics

#include "DomainInfo.h"
#include <string>
#include <vector>
#include "FB_Define.h"
#include "Vec3.h"
#include "MonCompo.h"
#include "Component.h"
#include "Control.h"
#include "CompoFraction.h"

using namespace std;
using namespace Vec3class;

/**
 *  モニタグループ管理クラス
 */
class MonitorList : public DomainInfo {
  
public:
  /// 出力タイプ型
  enum OutputType {
    GATHER,      ///< 単一ファイル出力
    DISTRIBUTE,  ///< 分散出力
    NONE
  };
  
protected:
  
  int num_process;       ///< プロセス数
  Vec3d org;             ///< ローカル基点座標
  Vec3d pch;             ///< セル幅
  Vec3d box;             ///< ローカル領域サイズ
  Vec3d g_org;           ///< グローバル基点座標
  Vec3d g_box;           ///< グローバル領域サイズ
  int* bid;              ///< 境界ID
  int* bcd;              ///< BCindex B
  float* cut;            ///< 交点情報
  TextParser* tpCntl;    ///< テキストパーサへのポインタ
  float area;            ///< 断面積 [m^2]
  int NoCompo;           ///< 物性テーブルの個数 >> 配列の大きさは[NoCompo+1]
  REAL_TYPE* mtbl;       ///< 物性テーブルへのポインタ
  
  OutputType outputType; ///< 出力タイプ
  vector<MonitorCompo*> monGroup;  ///< モニタリンググループ配列
  MonitorCompo::ReferenceVariables refVar;  ///< 参照用パラメータ変数

  string fname_sampling; ///< 基本ファイル名
  
  
public:
  /// コンストラクタ
  MonitorList()
  {
    outputType = NONE;
    num_process = 0;
    fname_sampling = "sampling.txt";
    area = 0.0;
    NoCompo = 0;
    mtbl = NULL;
  }
  
  /// デストラクタ
  ~MonitorList()
  {
    for (int i = 0; i < monGroup.size(); i++) delete monGroup[i];
  }
  
  
  /// @brief サンプリングの変数指定と計算モードの整合性をチェックする
  /// @param [in] KOS  KindOfSolver
  bool checkConsistency(const int KOS);
  
  
  /// @brief Polygonモニターの交点の状況を確認
  void checkStatus();
  
  
  /// @brief Polygonモニターの交点と境界IDの除去
  void clearCut();
  
  
  /// 出力ファイルクローズ
  ///
  ///    @note 出力タイプによらず全プロセスから呼んでも問題ない
  ///
  void closeFile();
  
  
  /**
   * @brief モニタ座標情報を取得し，リストに保持する
   * @param [in,out] C    Controlクラスオブジェクトのポインタ
   * @param [in,out] cmp  コンポーネント配列
   * @retval サンプリング指定のときtrue
   */
  bool getMonitor(Control* C, CompoList* cmp);
  
  
  /// VorticityとHelicityのサンプリングが指定されている場合にtrueを返す
  bool getStateVorticity();
  
  
  /**
   * @brief TPのポインタを受け取る
   * @param [in] tp TextParser
   */
  void importTP(TextParser* tp);
  
  
  
  /// 出力ファイルオープン
  ///
  ///    @note 出力タイプによらず全プロセスから呼んでも問題ない
  ///
  void openFile();
  
  
  /// モニタ結果出力(Line, PointSet指定)
  ///
  ///   @param [in] step サンプリング時の計算ステップ
  ///   @param [in] tm   サンプリング時の計算時刻
  ///
  ///   @note gatherの場合も全プロセスから呼ぶこと
  ///
  void print(const unsigned step, const double tm);
  
  
  /// モニタ情報を出力
  ///   @param [in] str 出力ファイルの基本名
  void printMonitorInfo(FILE* fp, const char* str, const bool verbose);
  
  
  /// サンプリング
  void sampling()
  {
    for (int i = 0; i < monGroup.size(); i++)
    {
      switch ( monGroup[i]->getType() )
      {
        case mon_LINE:
        case mon_POINT_SET:
          monGroup[i]->sampling();
          break;
          
        case mon_CYLINDER:
        case mon_BOX:
        case mon_POLYGON:
          monGroup[i]->samplingAverage();
          break;
      }
    }
  }
  
  
  /// 必要なパラメータのコピー
  ///
  ///   @param [in] bid            境界ID
  ///   @param [in] cut            交点情報
  ///   @param [in] bcd            BCindex B
  ///   @param [in] refVelocity    代表速度
  ///   @param [in] baseTemp       基準温度
  ///   @param [in] diffTemp       代表温度差
  ///   @param [in] refDensity     基準密度
  ///   @param [in] refLength      代表長さ
  ///   @param [in] basePrs        基準圧力
  ///   @param [in] modePrecision  出力精度指定フラグ (単精度，倍精度)
  ///   @param [in] unitPrs        圧力単位指定フラグ (絶対値，ゲージ圧)
  ///   @param [in] num_process    プロセス数
  ///   @param [in] m_NoCompo      物性テーブルの個数
  ///   @param [in] tbl            物性テーブル
  ///
  void setControlVars(int* bid,
                      float* cut,
                      int* bcd,
                      const double refVelocity,
                      const double baseTemp,
                      const double diffTemp,
                      const double refDensity,
                      const double refLength,
                      const double basePrs,
                      const int modePrecision,
                      const int unitPrs,
                      const int num_process,
                      const int m_NoCompo,
                      REAL_TYPE* tbl);
  
  
  /// サンプリング元データの登録
  ///
  ///   @param [in] v  速度変数配列
  ///   @param [in] p  圧力変数配列
  ///   @param [in] t  温度変数配列
  ///   @param [in] vr 渦度変数配列
  ///
  void setDataPtrs(REAL_TYPE* v, REAL_TYPE* p, REAL_TYPE* t, REAL_TYPE* vr)
  {
    for (int i = 0; i < monGroup.size(); i++) monGroup[i]->setDataPtrs(v, p, t, vr);
  }
  
  
  /**
   * @brief Polygonモニターの交点の数をcmp[]にセット
   * @param [in,out] cmp     CompoList
   * @param [in]     NoCompo リストの配列長
   */
  void setMonitorNpoint(CompoList* cmp, const int NoCompo);
  
  
  /// 参照速度のコピー
  ///   @param [in] v00 座標系移動速度
  void setV00(REAL_TYPE v00[4]);
  
  
  
  
  
protected:
  
  /// Line指定の端点座標をグローバル計算領域内にクリッピング
  void clipLine(double from[3], double to[3]);
  
  
  /// プリミティブ形状を作成
  void generatePrimitiveShape(const int odr,
                              const Monitor_Type montyp,
                              double nv[3],
                              double ctr[3],
                              double depth,
                              double tmp,
                              double height,
                              double dr[3]);
  
  
  /// モニタ座標情報(Line)を取得
  void getLine(const Control* C,
               const string label_base,
               double from[3],
               double to[3],
               int& nDivision);
  
  
  /// 出力タイプ文字列の取得
  string getOutputTypeStr();
  

  /// モニタ座標情報を取得(PointSet)
  void getPointset(const Control* C,
                   const string label_base,
                   vector<MonitorCompo::MonitorPoint>& pointSet);
  
  
  /// プリミティブのモニタ形状情報を取得
  void getPrimitive(Monitor_Type mon_type,
                    const string label_base,
                    double nv[3],
                    double center[3],
                    double& depth,
                    double& tmp,
                    double& height,
                    double dir[3]);
  
  
  /// 座標値を無次元化する
  void normalizeCord(double RefLength, double x[3])
  {
    x[0] /= RefLength;
    x[1] /= RefLength;
    x[2] /= RefLength;
  }
  
  
  /// サンプリング対象変数の登録
  void registVars(const string label, vector<string>& variables, const string key);
  
  
  /// Line登録
  void setLine(const char* str,
               vector<string>& variables,
               const char* method,
               const char* mode,
               double from[3],
               double to[3],
               int nDivision,
               Monitor_Type mon_type);
  
  
  /// 出力タイプの設定
  void setOutputType(OutputType type)
  {
    outputType = type;
  }
  
  
  /// PointSet登録
  void setPointSet(const char* str,
                   vector<string>& variables,
                   const char* method,
                   const char* mode,
                   vector<MonitorCompo::MonitorPoint>& pointSet,
                   Monitor_Type mon_type);
  
  
  /// Polygon登録
  void setPolygon(const char* str,
                  vector<string>& variables,
                  const char* method,
                  const char* mode,
                  const int order,
                  const double nv[3],
                  Monitor_Type mon_type);
  
  
  /// Primitiveの登録
  void setPrimitive(const char* str,
                    vector<string>& variables,
                    const char* method,
                    const char* mode,
                    Monitor_Type mon_type,
                    const double nv[3],
                    const double center[3],
                    const double& depth,
                    const double& tmp,
                    const double& height,
                    const double dir[3],
                    const int order);
  
  
  /// サンプリングと出力の次元の設定
  /// @param [in] unitOutput      出力単位指定フラグ (有次元，無次元)
  void setSamplingUnit(const int unitOutput)
  {
    refVar.modeUnitOutput = unitOutput;
  }
  
};

#endif // _FB_MONITOR_H_
