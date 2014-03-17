#ifndef _FB_MONITOR_COMPO_H_
#define _FB_MONITOR_COMPO_H_

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

/** 
 * @file   MonCompo.h
 * @brief  FlowBase MonitorCompo class Header
 * @author aics
 */

#include "DomainInfo.h"
#include <string>
#include <vector>
#include <cassert>
#include "FB_Define.h"
#include "vec3.h"
#include "basic_func.h"
#include "Component.h"
#include "Sampling.h"
#include "FBUtility.h"
#include "limits.h" // for UBUNTU

using namespace std;
using namespace Vec3class;

/**
 *  モニタグループクラス
 */
class MonitorCompo : public DomainInfo {
public:
  
  /// PointSet要素用構造体
  struct MonitorPoint {
    Vec3<REAL_TYPE> crd;      ///< モニタ点座標
    string label;       ///< モニタ点ラベル
    MonitorPoint(const REAL_TYPE v[3], const char* str) : crd(v), label(str) {}
    ~MonitorPoint() {}
  };
  
  /// 参照用パラメータ構造体
  struct ReferenceVariables {
    int modeUnitOutput;      /// 出力単位指定フラグ (有次元，無次元)
    int modePrecision;       /// 出力精度指定フラグ (単精度，倍精度)
    int unitPrs;             /// 圧力単位指定フラグ (絶対値，ゲージ圧)
    REAL_TYPE refVelocity;   /// 代表速度
    REAL_TYPE baseTemp;      /// 基準温度
    REAL_TYPE diffTemp;      /// 代表温度差
    REAL_TYPE refDensity;    /// 基準密度
    REAL_TYPE refLength;     /// 代表長さ
    REAL_TYPE basePrs;       /// 基準圧力
    Vec3<REAL_TYPE> v00;           /// 参照（座標系移動）速度
  };
  
  
protected:
  int nPoint;                ///< モニタ点数(Local)
  bool variable[var_END];    ///< モニタ変数フラグ
  string label;              ///< グループのラベル
  Monitor_Type monitor_type; ///< モニタ点指定タイプ
  int method;                ///< サンプリング方法
  int mode;                  ///< サンプリングモード
  int num_process;           ///< プロセス数
  int polyID;                ///< PolygonモニタのエントリID
  
  Sampling** mon;            ///< 「モニタ点毎のSampligクラスへのポインタ」の配列
  
  Vec3<REAL_TYPE> org;             ///< ローカル基点座標
  Vec3<REAL_TYPE> pch;             ///< セル幅
  Vec3<REAL_TYPE> box;             ///< ローカル領域サイズ
  Vec3<REAL_TYPE> g_org;           ///< グローバル基点座標
  Vec3<REAL_TYPE> g_box;           ///< グローバル領域サイズ
  REAL_TYPE nv[3];           ///< 法線ベクトル
  REAL_TYPE val[var_END];    ///< サンプリング値
  int NoCompo;               ///< 物性テーブルの個数 >> 配列の大きさは[NoCompo+1]
  REAL_TYPE* mtbl;           ///< 物性テーブル
  
  ReferenceVariables refVar; ///< 参照用パラメータ変数
  
  int* bid;            ///< 境界ID
  int* bcd;            ///< BCindex B
  float* cut;          ///< 交点情報
  
  FILE* fp;            ///< 出力ファイルポインタ
  
  // サンプリング元データ
  REAL_TYPE* vSource;   ///< 速度サンプリング元データ
  REAL_TYPE* pSource;   ///< 圧力サンプリング元データ
  REAL_TYPE* tSource;   ///< 温度サンプリング元データ
  REAL_TYPE* vrSource;  ///< 渦度サンプリング元データ
  
  
  Vec3<REAL_TYPE>* crd;      ///< モニタ点座標配列
  int* rank;           ///< モニタ点担当ランク番号配列
  string* comment;     ///< モニタ点コメント配列
  int* pointStatus;    ///< 不正モニタ点フラグ配列
  
  Vec3<REAL_TYPE>* vel;      ///< 速度サンプリング結果配列
  REAL_TYPE* prs;      ///< 圧力サンプリング結果配列
  REAL_TYPE* tmp;      ///< 温度サンプリング結果配列
  REAL_TYPE* tp;       ///< 全圧サンプリング結果配列
  REAL_TYPE* hlt;      ///< Helicityサンプリング結果配列
  Vec3<REAL_TYPE>* vor;      ///< 渦度サンプリング結果配列

  
public:
  /// デフォルトコンストラクタ
  MonitorCompo() {
    nPoint = 0;
    polyID = 0;
    nv[0] = nv[1] = nv[2] = 0.0;
    for (int i = 0; i < var_END; i++)
    {
      variable[i] = false;
      val[i] = 0.0;
    }
    NoCompo = 0;
    vel = vor = NULL;
    prs = tmp = tp = hlt = NULL;
    crd = NULL;
    rank = NULL;
    comment = NULL;
    pointStatus = NULL;
    vSource = pSource = tSource = vrSource = NULL;
    bid = NULL;
    bcd = NULL;
    cut = NULL;
    mtbl= NULL;
  }
  
  /// コンストラクタ
  ///
  ///   @param [in] org,pch,box ローカル領域基点座標，セル幅，領域サイズ
  ///   @param [in] g_org,g_box グローバル領域基点座標，領域サイズ
  ///   @param [in] size,guide  ローカルセルサイズ，ガイドセル数
  ///   @param [in] refVar      参照パラメータ
  ///   @param [in] bid         境界ID
  ///   @param [in] bcd         BCindex B
  ///   @param [in] cut         交点情報
  ///   @param [in] num_process プロセス数
  ///   @param [in] num_compo   物性テーブルの個数
  ///   @param [in] m_tbl       物性テーブル
  ///
  MonitorCompo(Vec3<REAL_TYPE> org,
               Vec3<REAL_TYPE> pch,
               Vec3<REAL_TYPE> box,
               Vec3<REAL_TYPE> g_org,
               Vec3<REAL_TYPE> g_box,
               ReferenceVariables refVar,
               int* bid,
               int* bcd,
               float* cut,
               const int num_process,
               const int num_compo,
               REAL_TYPE* m_tbl) {
    nPoint = 0;
    for (int i = 0; i < var_END; i++)
    {
      variable[i] = false;
      val[i] = 0.0;
    }
    vel = vor  = NULL;
    prs = tmp = tp  = NULL;
    crd = NULL;
    rank = NULL;
    comment = NULL;
    pointStatus = NULL;
    vSource = pSource = tSource = vrSource = NULL;
    
    this->org = org;
    this->pch = pch;
    this->box = box;
    this->g_org = g_org;
    this->g_box = g_box;
    this->refVar = refVar;
    this->bid = bid;
    this->bcd = bcd;
    this->cut = cut;
    this->num_process = num_process;
    this->NoCompo = num_compo;
    this->mtbl = m_tbl;
  }
  
  /// デストラクタ
  ~MonitorCompo() {
    if (crd) delete[] crd;
    if (vel) delete[] vel;
    if (prs) delete[] prs;
    if (tmp) delete[] tmp;
    if (tp) delete[] tp;
    if (rank) delete[] rank;
    if (comment) delete[] comment;
    if (pointStatus) delete[] pointStatus;
    
    
    if (mon)
    {
      for (int i = 0; i < nPoint; i++)
      {
        if (mon[i]) delete mon[i];
      }
      
      delete[] mon;
    }
  }
  
  
  /// モニタ点の状態を調べ，不正モニタ点フラグ配列pointStatusを設定
  void checkMonitorPoints();
  
  
  /// モニタ点が指定領域内にあるかを判定
  ///
  ///   @param [in] m    モニタ点番号
  ///   @param [in] org  調査領域の基点
  ///   @param [in] box  調査領域のサイズ
  ///   @param [in] flag メッセージ出力フラグ(trueの時出力)
  ///   @return true=領域内/false=領域外
  ///
  bool checkRegion(const int m, const Vec3<REAL_TYPE> org, const Vec3<REAL_TYPE> box, bool flag=false) const
  {
    if ((crd[m].x< org.x)         ||
        (crd[m].x>(org.x+box.x))  ||
        (crd[m].y< org.y)         ||
        (crd[m].y>(org.y+box.y))  ||
        (crd[m].z< org.z)         ||
        (crd[m].z>(org.z+box.z)) )
    {
      if (flag)
      {
        stamped_printf("\trank=%d : no.=%d [%14.6e %14.6e %14.6e] is out of region\n", myRank, m,
                       convCrd(crd[m].x), convCrd(crd[m].y), convCrd(crd[m].z));
      }
      
      return false;
    }
    
    return true;
  }
  
  
  ///セルモニタの場合の交点情報をクリアする
  unsigned long clearMonitorCut();
  
  
  /// 出力ファイルクローズ
  void closeFile();
  
  
  /// グループラベルを返す
  const char* getLabel() const
  {
    return label.c_str();
  }
  
  
  /// m番目のモニタ点を含むセルインデクスを返す
  Vec3i getSamplingCellIndex(int m) const
  {
    Vec3i index;
    Vec3<REAL_TYPE> c = (crd[m] - org) / pch;
    index.x = int(c.x) + 1;
    index.y = int(c.y) + 1;
    index.z = int(c.z) + 1;
    return index;
  }
  
  
  /// エントリを返す
  int getPolyID() const
  {
    return polyID;
  }
  
  
  /// モニタ点数を返す
  int getSize() const
  {
    return nPoint;
  }
  
  
  /// モニタ変数の状態を返す
  bool getStateVariable(int var) const
  {
    return variable[var];
  }
  
  
  /// 登録タイプを返す
  Monitor_Type getType() const
  {
    return monitor_type;
  }
  
  
  /// 出力ファイルオープン
  ///
  ///    @param str ファイル名テンプレート
  ///    @param gathered true=gather出力/false=disutribute出力
  ///
  void openFile(const char* str, const bool gathered);
  
  
  /// モニタ結果出力(distribute)
  ///
  ///   @param [in] step サンプリング時の計算ステップ
  ///   @param [in] tm   サンプリング時の計算時刻
  ///
  ///   @note 全プロセスから呼ぶこと
  ///
  void print_distribute(const unsigned step, const REAL_TYPE tm)
  {
    print(step, tm, false);
  }
  
  
  /// モニタ結果出力(gather)
  ///
  ///   @param [in] step サンプリング時の計算ステップ
  ///   @param [in] tm   サンプリング時の計算時刻
  ///
  ///   @note 内部でノード0に集計するため全プロセスから呼ぶこと
  ///
  void print_gather(const unsigned step, const REAL_TYPE tm)
  {
    gatherSampled();
    
    if (myRank == 0) print(step, tm, true);
  }
  
  
  /// モニタ情報を出力
  ///
  ///    @param [in] fp 出力ファイルポインタ
  ///    @param [in] no モニタグループ通し番号
  //
  void printInfo(FILE* fp, int no);
  
  
  /// 詳細サンプリング(Line, PointSet)
  void sampling();
  
  
  /// 平均値サンプリング
  ///
  ///   サンプリング結果を集計，領域での平均値を計算
  ///   速度は法線ベクトルとの内積をとる
  ///
  void samplingAverage();
  
  
  /// サンプリング元データの登録
  ///
  ///   @param [in] v  速度変数配列
  ///   @param [in] p  圧力変数配列
  ///   @param [in] t  温度変数配列
  ///   @param [in] vr 渦度変数配列
  ///
  void setDataPtrs(REAL_TYPE* v, REAL_TYPE* p, REAL_TYPE* t, REAL_TYPE* vr)
  {
    vSource  = v;
    pSource  = p;
    tSource  = t;
    vrSource = vr;
    
    
    if (variable[var_Velocity] && !vSource)
    {
      Hostonly_ stamped_printf("\tError : Null monitoring source for Velocity.\n");
      Exit(0);
    }
    
    if (variable[var_Pressure] && !pSource)
    {
      Hostonly_ stamped_printf("\tError : Null monitoring source for Pressure.\n");
      Exit(0);
    }
    
    if (variable[var_Temperature] && !tSource)
    {
      Hostonly_ stamped_printf("\tError : Null monitoring source for Temperature.\n");
      Exit(0);
    }

    
    if (variable[var_Vorticity] && !vrSource)
    {
      Hostonly_ stamped_printf("\tError : Null monitoring source for Vorticity.\n");
      Exit(0);
    }
    
    if (variable[var_TotalP] && !(vSource && pSource) )
    {
      Hostonly_ stamped_printf("\tError : Null monitoring source for Vorticity.\n");
      Exit(0);
    }
    
    if (variable[var_Helicity] && !vSource )
    {
      Hostonly_ stamped_printf("\tError : Null monitoring source for Vorticity.\n");
      Exit(0);
    }
  }
  
  
  /// Line登録
  ///
  ///   @param [in] labelStr  ラベル文字列
  ///   @param [in] variables モニタ変数vector
  ///   @param [in] methodStr method文字列
  ///   @param [in] modeStr   mode文字列
  ///   @param [in] from      Line始点
  ///   @param [in] to        Line終点
  ///   @param [in] nDivision 分割数(モニタ点数-1)
  ///   @param [in] m_type    モニタタイプ
  ///
  void setLine(const char* labelStr,
               vector<string>& variables,
               const char* methodStr,
               const char* modeStr,
               REAL_TYPE from[3],
               REAL_TYPE to[3],
               int nDivision,
               Monitor_Type m_type);
  
  
  /// PointSet登録
  ///
  ///   @param [in] labelStr  ラベル文字列
  ///   @param [in] variables モニタ変数vector
  ///   @param [in] methodStr method文字列
  ///   @param [in] modeStr   mode文字列
  ///   @param [in] pointSet  PointSet
  ///   @param [in] m_type    モニタタイプ
  ///
  void setPointSet(const char* lableStr,
                   vector<string>& variables,
                   const char* methodStr,
                   const char* modeStr,
                   vector<MonitorPoint>& pointSet,
                   Monitor_Type m_type);

  
  /// Polygon登録
  ///
  ///   @param [in] labelStr  ラベル文字列
  ///   @param [in] variables モニタ変数vector
  ///   @param [in] methodStr method文字列
  ///   @param [in] modeStr   mode文字列
  ///   @param [in] order     エントリ番号
  ///   @param [in] m_nv      法線ベクトル
  ///   @param [in] m_type    モニタタイプ
  ///
  void setPolygon(const char* labelStr,
                  vector<string>& variables,
                  const char* methodStr,
                  const char* modeStr,
                  const int order,
                  const REAL_TYPE m_nv[3],
                  Monitor_Type m_type);
  
  
  /// Primitiveを登録し，bcd[]をゼロクリア
  ///
  ///   @param [in] labelStr  ラベル文字列
  ///   @param [in] variables モニタ変数vector
  ///   @param [in] methodStr method文字列
  ///   @param [in] modeStr   mode文字列
  ///   @param [in] order     エントリ番号
  ///   @param [in] m_nv      法線ベクトル
  ///   @param [in] m_type    モニタタイプ
  ///
  void setPrimitive(const char* labelStr,
                    vector<string>& variables,
                    const char* methodStr,
                    const char* modeStr,
                    const int order,
                    const REAL_TYPE m_nv[3],
                    Monitor_Type m_type);
  
  
  
protected:
  
  /// モニタリング管理用配列の確保
  void allocArray();
  
  
  /// サンプリング値を格納する配列の確保
  void allocSamplingArray();
  
  
  /// Allreduceによる総和(実数配列上書き，work配列指定)
  bool allReduceSum(REAL_TYPE* array, int n, REAL_TYPE* sendBuf);
  
  
  /// Allreduceによる総和(実数配列上書き)
  bool allReduceSum(REAL_TYPE* array, int n);
  
  
  /// Allreduceによる総和(整数配列上書き，work配列指定)
  bool allReduceSum(int* array, int n, unsigned long* sendBuf);
  
  
  /// Allreduceによる総和(整数配列上書き)
  bool allReduceSum(int* array, int n);
  
  
  /// 指定されたモニタ領域内でスカラー変数を平均
  REAL_TYPE averageScalar(REAL_TYPE* s);
  
  
  /// 指定されたモニタ領域内でベクトル変数を平均
  Vec3<REAL_TYPE> averageVector(Vec3<REAL_TYPE>* v);
  
  
  /// 座標の単位変換
  REAL_TYPE convCrd(REAL_TYPE xyz) const
  {
    return ( (refVar.modeUnitOutput==DIMENSIONAL) ? xyz*refVar.refLength : xyz );
  }
  
  /// 時間の単位変換
  REAL_TYPE convTime(REAL_TYPE tm) const
  {
    return ( (refVar.modeUnitOutput==DIMENSIONAL) ? tm*refVar.refLength/refVar.refVelocity : tm );
  }
  
  /// 速度成分の単位変換
  REAL_TYPE convVel(REAL_TYPE vel) const
  {
    return ( (refVar.modeUnitOutput==DIMENSIONAL) ? vel*refVar.refVelocity : vel );
  }
  
  /// 圧力の単位変換
  REAL_TYPE convPrs(REAL_TYPE prs) const
  {
    return ( (refVar.modeUnitOutput==DIMENSIONAL) ? FBUtility::convPrsND2D(prs, refVar.basePrs, refVar.refDensity, refVar.refVelocity, refVar.unitPrs) : prs );
  }
  
  /// 温度の単位変換
  REAL_TYPE convTmp(REAL_TYPE tmp) const
  {
    return ( (refVar.modeUnitOutput==DIMENSIONAL) ? FBUtility::convTempND2D(tmp, refVar.baseTemp, refVar.diffTemp) : tmp);
  }
  
  /// 全圧の単位変換
  REAL_TYPE convTP(REAL_TYPE tp) const
  {
    return ( (refVar.modeUnitOutput==DIMENSIONAL) ? tp * (refVar.refVelocity * refVar.refVelocity * refVar.refDensity) : tp );
  }
  
  /// 渦度成分の単位変換
  REAL_TYPE convVor(REAL_TYPE vor) const
  {
    return ( (refVar.modeUnitOutput==DIMENSIONAL) ? vor*refVar.refVelocity/refVar.refLength : vor );
  }
  
  /// Helicityの単位変換
  REAL_TYPE convHlt(REAL_TYPE hlt) const
  {
    return ( (refVar.modeUnitOutput==DIMENSIONAL) ? hlt*refVar.refVelocity*refVar.refVelocity/refVar.refLength : hlt );
  }
  
  /// サンプリングした変数をノード0に集約
  void gatherSampled();
  
  
  /// サンプリングしたスカラー変数をノード0に集約
  ///
  ///   @param [in,out] s スカラー変数配列
  ///   @param  sRecvBuf  通信用work領域
  ///
  void gatherSampledScalar(REAL_TYPE* s, REAL_TYPE* sRecvBuf);
  
  
  /// サンプリングしたベクトル変数をノード0に集約
  ///
  ///   @param [in,out] v ベクトル変数配列
  ///   @param  vSendBuf,vRecvBuf  通信用work領域
  ///
  void gatherSampledVector(Vec3<REAL_TYPE>* v, REAL_TYPE* vSendBuf, REAL_TYPE* vRecvBuf);

  
  
  /// サンプリング方法文字列の取得
  string getMethodStr();
  
  
  /// サンプリングモード文字列の取得
  string getModeStr();
  
  
  /// モニタ点指定方法文字列の取得
  string getTypeStr();
  
  
  /// モニタ変数を結合した文字列の取得
  string getVarStr();
  
  
  /// モニタ結果出力
  ///
  ///   @param [in] step サンプリング時の計算ステップ
  ///   @param [in] tm   サンプリング時の計算時刻
  ///   @param [in] gathered 出力モードフラグ(true=gather出力/false=disutribute出力)
  ///
  void print(unsigned step, REAL_TYPE tm, bool gathered);
  
  
  /*
   * @brief 5bit幅の値の設定
   * @param [in,out] b   int 変数
   * @param [in]     q   5-bit幅のID (1-31)
   * @param [in]     dir 方向コード (w/X_MINUS=0, e/X_PLUS=1, s/2, n/3, b/4, t/5)
   
   */
  inline void setBit5 (int& b, const int q, const int dir)
  {
    b &= (~(0x1f << (dir*5)) ); // 対象5bitをゼロにする
    b |= (q << (dir*5));        // 書き込む
  }
  
  
  /// モニタ対象物理量の設定
  ///
  ///    @param [in] str モニタ対象物理量文字列
  ///
  void setMonitorVar(const char* str);
  
  
  /// 各モニタ点を担当するランク番号を配列rank[]にセット
  ///
  ///   @note 領域境界上のモニタ点は，ランク番号の大きい方の領域が担当
  ///
  void setRankArray();
  
  
  /// サンプリング方法の設定
  ///
  ///    @param [in] str サンプリング方法文字列
  ///
  void setSamplingMethod(const char* str);
  
  
  /// サンプリングモードの設定
  ///
  ///    @param [in] str サンプリングモード文字列
  ///
  void setSamplingMode(const char* str);
  
  
  /// モニタ結果出力ファイルにヘッダ部を出力
  ///
  ///   @param [in] gathered 出力モードフラグ(true=gather出力/false=disutribute出力)
  ///
  void writeHeader(bool gathered);
  
  
  /**
   * @brief mon_POLYGON, mon_BOX, mon_CYLINDERのヘッダー出力
   */
  void writeHeaderCompo();
};

#endif // _FB_MONITOR_COMPO_H_
