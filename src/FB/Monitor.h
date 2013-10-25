#ifndef _FB_MONITOR_H_
#define _FB_MONITOR_H_

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

//@file   Monitor.h
//@brief  FlowBase MonitorList class Header
//@author kero

#include "DomainInfo.h"
#include <string>
#include <vector>
#include "FB_Define.h"
#include "vec3.h"
#include "MonCompo.h"
#include "basic_func.h"
#include "Component.h"
#include "Control.h"
#include "CompoFraction.h"

using namespace std;

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
  
  int nGroup;            ///< モニタリンググループ数
  int num_process;       ///< プロセス数
  FB::Vec3r org;         ///< ローカル基点座標
  FB::Vec3r pch;         ///< セル幅
  FB::Vec3r box;         ///< ローカル領域サイズ
  FB::Vec3r g_org;       ///< グローバル基点座標
  FB::Vec3r g_box;       ///< グローバル領域サイズ
  int* bid;              ///< 境界ID
  int* bcd;              ///< BCindex B
  float* cut;            ///< 交点情報
  TextParser* tpCntl;    ///< テキストパーサへのポインタ
  float area;            ///< 断面積 [m^2]
  
  OutputType outputType; ///< 出力タイプ
  vector<MonitorCompo*> monGroup;  ///< モニタリンググループ配列
  MonitorCompo:: ReferenceVariables refVar;  ///< 参照用パラメータ変数

  string fname_integral; ///< mon_POLYGON, mon_BOX, mon_CYLINDER
  string fname_detail;   ///< mon_POINT_SET, mon_LINE
  
  
public:
  /// コンストラクタ
  MonitorList()
  {
    nGroup = 0;
    outputType = NONE;
    num_process = 0;
    fname_integral = "sampling_compo.txt";
    fname_detail   = "sampling.txt";
    area = 0.0;
  }
  
  /// デストラクタ
  ~MonitorList()
  {
    for (int i = 0; i < nGroup; i++) delete monGroup[i];
  }
  
  
  /// Polygonモニターの交点の状況を確認
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
   * @param [in,out] C   Controlクラスオブジェクトのポインタ
   * @param [in,out] cmp コンポーネント配列
   * @retval サンプリング指定のときtrue
   */
  bool getMonitor(Control* C, CompoList* cmp);
  
  
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
  void print(const unsigned step, const REAL_TYPE tm);
  
  
  /// モニタ情報を出力
  ///   @param [in] str 出力ファイルの基本名
  void printMonitorInfo(FILE* fp, const char* str, const bool verbose);
  
  
  /// サンプリング
  void sampling()
  {
    for (int i = 0; i < nGroup; i++)
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
  ///
  void setControlVars(int* bid,
                      float* cut,
                      int* bcd,
                      const REAL_TYPE refVelocity,
                      const REAL_TYPE baseTemp,
                      const REAL_TYPE diffTemp,
                      const REAL_TYPE refDensity,
                      const REAL_TYPE refLength,
                      const REAL_TYPE basePrs,
                      const int modePrecision,
                      const int unitPrs,
                      const int num_process);
  
  
  /// サンプリング元データの登録
  ///
  ///   @param [in] v 速度変数配列
  ///   @param [in] p 圧力変数配列
  ///   @param [in] t 温度変数配列
  ///
  void setDataPtrs(REAL_TYPE* v, REAL_TYPE* p, REAL_TYPE* t=NULL)
  {
    for (int i = 0; i < nGroup; i++) monGroup[i]->setDataPtrs(v, p, t);
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
  void clipLine(REAL_TYPE from[3], REAL_TYPE to[3]);
  
  
  /// プリミティブ形状を作成
  void generatePrimitiveShape(const int odr,
                              const Monitor_Type montyp,
                              REAL_TYPE nv[3],
                              REAL_TYPE ctr[3],
                              REAL_TYPE depth,
                              REAL_TYPE tmp,
                              REAL_TYPE height,
                              REAL_TYPE dr[3]);
  
  
  /// モニタ座標情報(Line)を取得
  void getLine(const Control* C,
               const string label_base,
               REAL_TYPE from[3],
               REAL_TYPE to[3],
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
                    REAL_TYPE nv[3],
                    REAL_TYPE center[3],
                    REAL_TYPE& depth,
                    REAL_TYPE& tmp,
                    REAL_TYPE& height,
                    REAL_TYPE dir[3]);
  
  
  /// 座標値を無次元化する
  void normalizeCord(REAL_TYPE RefLength, REAL_TYPE x[3])
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
               REAL_TYPE from[3],
               REAL_TYPE to[3],
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
                  const REAL_TYPE nv[3],
                  Monitor_Type mon_type);
  
  
  /// Primitiveの登録
  void setPrimitive(const char* str,
                    vector<string>& variables,
                    const char* method,
                    const char* mode,
                    Monitor_Type mon_type,
                    const REAL_TYPE nv[3],
                    const REAL_TYPE center[3],
                    const REAL_TYPE& depth,
                    const REAL_TYPE& tmp,
                    const REAL_TYPE& height,
                    const REAL_TYPE dir[3],
                    const int order);
  
  
  /// サンプリングと出力の次元の設定
  /// @param [in] modeUnit       出力単位指定フラグ (有次元，無次元)
  void setSamplingUnit(const int unit)
  {
    refVar.modeUnit = unit;
  }
  
};

#endif // _FB_MONITOR_H_
