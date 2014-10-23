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

// @file   Monitor.h
// @brief  FlowBase MonitorList class Header
// @author aics
// @note 浮動小数点の型はREAL_TYPE, mtblのみdouble

#include "DomainInfo.h"
#include <vector>
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
  int nGroup;            ///< モニタリンググループ数
  int num_process;       ///< プロセス数
  Vec3r org;             ///< ローカル基点座標
  Vec3r pch;             ///< セル幅
  Vec3r box;             ///< ローカル領域サイズ
  Vec3r g_org;           ///< グローバル基点座標
  Vec3r g_box;           ///< グローバル領域サイズ
  int* bid;              ///< 境界ID
  int* bcd;              ///< BCindex B
  float* cut;            ///< 交点情報
  TextParser* tpCntl;    ///< テキストパーサへのポインタ
  REAL_TYPE area;        ///< 断面積 [m^2]
  int NoCompo;           ///< 物性テーブルの個数 >> 配列の大きさは[NoCompo+1]
  double* mtbl;          ///< 物性テーブルへのポインタ
  
  OutputType outputType; ///< 出力タイプ
  vector<MonitorCompo*> monGroup;  ///< モニタリンググループ配列
  MonitorCompo::ReferenceVariables refVar;  ///< 参照用パラメータ変数

  string fname_sampling; ///< 基本ファイル名
  
  
public:
  /// コンストラクタ
  MonitorList()
  {
    nGroup = 0;
    outputType = NONE;
    num_process = 0;
    fname_sampling = "sampling.txt";
    area = 0.0;
    NoCompo = 0;
    mtbl = NULL;
    tpCntl = NULL;
  }
  
  /// デストラクタ
  ~MonitorList()
  {
    for (int i = 0; i < nGroup; i++) delete monGroup[i];
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
  
  
  /// TextParserの中身を出力する
  ///
  void printTextParser(std::string given_label="");
  
  
  /// サンプリング
  void sampling()
  {
    for (int i = 0; i < nGroup; i++)
    {
      switch ( monGroup[i]->getType() )
      {
        case mon_LINE:
          monGroup[i]->sampling();
          break;
          
        case mon_POINT_SET:
          monGroup[i]->sampling();
          break;
          
        case mon_CYLINDER:
          break;
          
        case mon_BOX:
          break;
          
        case mon_POLYGON:
          monGroup[i]->samplingAverage();
          break;
          
        default:
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
  ///   @param [in] m_mtbl         物性テーブル
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
                      const int num_process,
                      const int m_NoCompo,
                      double* m_mtbl);
  
  
  /// サンプリング元データの登録
  ///
  ///   @param [in] v  速度変数配列
  ///   @param [in] p  圧力変数配列
  ///   @param [in] t  温度変数配列
  ///   @param [in] vr 渦度変数配列
  ///
  void setDataPtrs(REAL_TYPE* v, REAL_TYPE* p, REAL_TYPE* t, REAL_TYPE* vr)
  {
    for (int i = 0; i < nGroup; i++) monGroup[i]->setDataPtrs(v, p, t, vr);
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
  
  //>> Graph ploter
  void getBox(Monitor_Type mon_type,
              const string label_base,
              REAL_TYPE center[3],
              REAL_TYPE main_vec[3],
              REAL_TYPE ref_vec[3],
              REAL_TYPE dim_vec[3]);
  
  void getCylinder(Monitor_Type mon_type,
                   const string label_base,
                   REAL_TYPE center[3],
                   REAL_TYPE main_vec[3],
                   REAL_TYPE ref_vec[3],
                   REAL_TYPE dim_vec[3]);
  
  void getPlane(Monitor_Type mon_type,
                const string label_base,
                REAL_TYPE center[3],
                REAL_TYPE main_vec[3],
                REAL_TYPE ref_vec[3],
                REAL_TYPE dim_vec[2],
                REAL_TYPE div_vec[2]);
  
  void getPolygon(Monitor_Type mon_type,
                  const string label_base,
                  string & stl_filename );
  
  //指定点　gp　は、指定された座標系内部のローカル座標を求める
  Vec3r localPt( Vec3r orig, Vec3r unit_z, Vec3r unit_x, Vec3r unit_y, Vec3r global_pt )
  {
    Vec3r vec = global_pt - orig;
    Vec3r local_pt( dot(vec, unit_x), dot(vec, unit_y), dot(vec, unit_z) );
    //lp->setIt( vec*unitX, vec*unitY, vec*unitZ );
    return local_pt;
  }
  
  //指定された座標系内部のローカル点　lp　のグローバル座標を求める
  Vec3r globalPt( Vec3r orig, Vec3r unit_z, Vec3r unit_x, Vec3r unit_y, Vec3r local_pt )
  {
    Vec3r global_pt = orig + unit_x*local_pt.x + unit_y*local_pt.y + unit_z*local_pt.z;
    return global_pt;
  }
  
  void getPointsInBox( Vec3r orig_o, Vec3r axis_z, Vec3r axis_x, REAL_TYPE dim[3], vector<MonitorCompo::MonitorPoint> *pointSet );
  void getPointsInCyl( Vec3r orig_o, Vec3r axis_z, Vec3r axis_x, REAL_TYPE dim[3], vector<MonitorCompo::MonitorPoint> *pointSet );
  void getPointsOnPlane( Vec3r orig_o, Vec3r axis_z, Vec3r axis_x, REAL_TYPE dim[2], REAL_TYPE div[2],
                        vector<MonitorCompo::MonitorPoint> *pointSet, vector<int> *grid_flags );
  
  //ローカル座標系内部で判定する。原点は箱の中心
  bool isInBox(Vec3r pt, REAL_TYPE dim[3])
  {
    REAL_TYPE hdx = dim[0]*0.5;
    REAL_TYPE hdy = dim[1]*0.5;
    REAL_TYPE hdz = dim[2]*0.5;
    if( pt.x < -hdx || pt.x > hdx ) return false;
    if( pt.y < -hdy || pt.y > hdy ) return false;
    if( pt.z < -hdz || pt.z > hdz ) return false;
    return true;
  }
  
  //ローカル座標系内部で判定する。原点は円柱の底面中心
  bool isInCyl(Vec3r pt, REAL_TYPE dim[3])
  {
    static REAL_TYPE _tol = 1.0e-6;
    
    REAL_TYPE r1 = dim[0];
    REAL_TYPE r2 = dim[1];
    REAL_TYPE h  = dim[2];
    
    if( pt.z > h || pt.z < 0.0 ) return false;
    
    REAL_TYPE rz = r1 + (r2-r1) * pt.z / h;
    REAL_TYPE rr = sqrt( pt.x*pt.x + pt.y*pt.y );
    
    return( rr<rz ? true : false );
  }
  //<< Graph ploter
  
  
  /// 座標値を無次元化する
  /// @param [in,out] x  coordinate
  void normalizeCord(REAL_TYPE x[3])
  {
    x[0] /= refVar.refLength;
    x[1] /= refVar.refLength;
    x[2] /= refVar.refLength;
  }
  
  
  /// サンプリング対象変数の登録
  void registVars(const string label, vector<string>& variables, const string key);
  
  
  /// Line登録
  MonitorCompo* setLine(const char* str,
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
  MonitorCompo* setPointSet(const char* str,
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
                    const int order,
                    const REAL_TYPE nv[3],
                    Monitor_Type mon_type);
  
  
  /// サンプリングと出力の次元の設定
  /// @param [in] unitOutput      出力単位指定フラグ (有次元，無次元)
  void setSamplingUnit(const int unitOutput)
  {
    refVar.modeUnitOutput = unitOutput;
  }
  
};

#endif // _FB_MONITOR_H_
