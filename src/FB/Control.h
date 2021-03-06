#ifndef _FB_CONTROL_H_
#define _FB_CONTROL_H_

//##################################################################################
//
// Flow Base class
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2015 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2016 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
// Copyright (c) 2016-2020 Research Institute for Information Research, Kyushu University.
// All rights reserved.
//
//##################################################################################

/*
 * @file   Control.h
 * @brief  FlowBase Control class Header
 * @author aics
 */

#include "cpm_Define.h"
#include <math.h>
#include "DomainInfo.h"
#include "FB_Define.h"
#include "Medium.h"
#include "Component.h"
#include "FBUtility.h"
#include "BndOuter.h"
#include "common/Vec3.h" // defined in Polylib
#include "IterationControl.h"
#include "TextParser.h"
#include "IntervalManager.h"
#include "mydebug.h"

using namespace std;
using namespace Vec3class;



// #################################################################
class DTcntl {
public:
  /** 時間積分幅の指定種別 */
  enum dt_Type
  {
    dt_direct=1,      ///< 入力値がΔt
    dt_cfl_ref_v,     ///< dt < c dx/U0
    dt_cfl_max_v,     ///< dt < c dx/Umax
    dt_dfn,           ///< 拡散数制限
    dt_cfl_dfn_ref_v, ///< dt = min( c dx/U0, diffusion number )
    dt_cfl_dfn_max_v, ///< dt = min( c dx/Umax, diffusion number )
    dt_cfl_max_v_cp   ///< 圧縮性　dt < (cfl+soundSpeed) dx/Umax
  };

private:
  int scheme;      ///< Δtのスキーム種別
  int KOS;         ///< Kind of Solver
  int mode;        ///< 入力パラメータの次元モード（無次元/有次元）
  double CFL;      ///< Δtの決定に使うCFLなど
  double deltaT;   ///< Δt（無次元）
  double min_dx;   ///< 最小格子幅（無次元）
  double Reynolds; ///< レイノルズ数
  double Peclet;   ///< ペクレ数

public:
  /** コンストラクタ */
  DTcntl() {
    scheme = 0;
    CFL = 0.0;
    KOS = 0;
    mode = 0;
    deltaT = 0.0;
    min_dx = Reynolds = Peclet = 0.0;
  }

  /**　デストラクタ */
  ~DTcntl() {}

public:
  /** スキームの種類を返す */
  int get_Scheme() const
  {
    return scheme;
  }

  /** 時間積分幅を返す */
  double get_DT() const
  {
    return deltaT;
  }

  /** CFL数を返す */
  double get_CFL() const
  {
    return CFL;
  }

  /**
   * @brief CFL数で指定されるdtを計算
   * @param [in] Uref 速度の参照値（無次元）
   */
  double dtCFL(const double Uref) const
  {
    double v = (Uref < 1.0) ? 1.0 : Uref; // 1.0 is non-dimensional reference velocity
    return (min_dx*CFL / v);
  }


  /** 拡散数のdt制限
   * @param [in] coef 係数 (Reynolds number or Peclet number)
   */
  double dtDFN(const double coef) const
  {
    return coef * min_dx*min_dx/6.0;
  }


  /** 時間積分幅とKindOfSolver種別の整合性をチェック */
  bool chkDtSelect();


  /**
   * @brief Δtのスキームを設定する
   * @retval エラーコード
   * @param [in] str  キーワード
   * @param [in] val  値
   */
  bool setScheme(const char* str, const double val);


  /**
   * @brief 各種モードに対応する時間積分幅を設定する
   * @retval return code
   * @param [in] vRef 速度の参照値（無次元）
   * @note deltaTは無次元
   */
  int set_DT(const double vRef);


  /**
   * @brief 基本変数をコピー
   * @param [in] m_kos     ソルバーの種類
   * @param [in] m_mode    次元モード
   * @param [in] m_min_dx  無次元最小格子幅
   * @param [in] re        レイノルズ数
   * @param [in] pe        ペクレ数
   */
  void set_Vars(const int m_kos, const int m_mode, const double m_min_dx, const double re, const double pe);
};




// #################################################################
class Control : public DomainInfo {
protected:
  TextParser* tpCntl;   ///< テキストパーサへのポインタ

public:

  /** 各種モード　パラメータ */
  typedef struct
  {
    int Statistic;
    int StatisticRestart;
    int Buoyancy;
    int Example;
    int Log_Base;
    int Log_Wall;
    int PDE;
    int Precision;
    int Profiling;
    int PrsNeuamnnType;
    int ShapeAprx;
    int Steady;
    int Base_Medium;
    int CCNV;
    int StatVelocity;
    int StatPressure;
    int StatTemperature;
    int ReynoldsStress;
    int ChannelOutputIter;
    int ChannelOutputMean;
    int ParticleTracking;
  } Mode_set;

  /** 隠しパラメータ */
  typedef struct
  {
    int Range_Limit;
    int PM_Test;
    int GeomOutput;
    int GlyphOutput;
    int DryRun;
  } Hidden_Parameter;

  // AXB係数書き出しオプション
  typedef struct
  {
    int func;
    int interval;
    int threshold;
  } AXB_param;

  /** LESパラメータ */
  typedef struct
  {
    int Calc;
    int Model;
    int InitialPerturbation;
    int ChannelDir;
    int VelocityProfile;
    REAL_TYPE Cs;
    REAL_TYPE damping_factor;
    REAL_TYPE BulkVelocity;
    REAL_TYPE ChannelWidth;
    REAL_TYPE TurbulentReynoldsNum;
  } LES_Parameter;

  /** 初期値 */
  typedef struct
  {
    REAL_TYPE Density;
    REAL_TYPE Energy;
    REAL_TYPE Pressure;
    REAL_TYPE VecU;
    REAL_TYPE VecV;
    REAL_TYPE VecW;
  } Initial_Value;

  /** 単位 */
  typedef struct
  {
    int Param;  /// 入力パラメータ単位 (Dimensional/NonDimensional)
    int Output; /// 出力単位 (Dimensional/NonDimensional)
    int File;   /// ファイルの記述単位 (Dimensional/NonDimensional)
    int Log;    /// 出力ログの単位 (Dimensional/NonDimensional)
    int Prs;    /// 圧力単位 (Absolute/Gauge)
    int Length; /// 入力パラメータの長さの単位 (non_dimensional/m/cm/mm)
  } Unit_Def;


  /** コンポーネント/BCのグローバルな存在 */
  typedef struct
  {
    int forcing;
    int hsrc;
    int outflow;
    int periodic;
    int fraction;
    int tfree;
    int vspec;
    int monitor;
    int obstacle;
  } Ens_of_Compo;

  /** ドライバ **/
  typedef struct
  {
    REAL_TYPE length;      ///< ドライバの長さ
    std::string Label;      ///< ドライバ部分のラベル
    std::string faceLabel; ///< ドライバ指定面のラベル
  } Driver_Def;


  // 安定化のパラメータ
  typedef struct
  {
    int control;
    REAL_TYPE begin;
    REAL_TYPE end;
    REAL_TYPE penalty_number;
  } Stability_Control;


  /** 偏微分方程式の型 */
  enum PDE_type
   {
    PDE_NS=1,
    PDE_EULER
  };


  /** 値域の制御　*/
  enum Range_Cntl
  {
    Range_Normal=1,
    Range_Cutoff
  };

  /** 出力タイミングの指定 */
  enum output_mode
  {
    IO_normal,
    IO_forced,
    IO_everytime
  };

  /** 対流項スキーム */
  enum convection_scheme
  {
    O1_upwind=1,
    O2_central,
    O3_muscl,
    O4_central
  };

  /** TVD limiter */
  enum Limiter_function
  {
    No_Limiter,
    MinMod
  };

  /** LES Model */
  enum LES_model
  {
    LES_no=0,
    LES_Smagorinsky,
    LES_CSM,
    LES_WALE
  };

  /** 壁面の扱い */
  enum Wall_Profile
  {
    No_Slip=0,   // 0
    Slip,        // 1
    Law_of_Wall, // 2
    Van_Driest   // 3
  };


  /** 並列化モード */
  enum Parallel_mode
  {
    Serial=1,
    OpenMP,
    FlatMPI,
    Hybrid
  };

  /** ファイル出力モード */
  enum File_IO_Mode
  {
    io_current=1,
    io_specified
  };

  /** 管理対象のリスト */
  enum flush_trigger
  {
    tg_compute=0,  ///< セッションの計算時間
    tg_console,    ///< コンソール出力
    tg_history,    ///< ファイル出力
    tg_basic,      ///< 基本変数の瞬時値出力
    tg_statistic,  ///< 統計値出力
    tg_derived,    ///< 派生変数の出力
    tg_accelra,    ///< 加速時間
    tg_sampled,    ///< サンプリング出力
    tg_END
  };


  int AlgorithmF;
  int AlgorithmH;
  int BasicEqs;
  int CheckParam;
  int CnvScheme;
  int GuideOut;
  int KindOfSolver;
  int Limiter;
  int MarchingScheme;
  int NoBaseLS;       ///< リストアップされた線形ソルバーの数
  int NoBC;           ///< 境界条件数
  int NoBCinner;      ///< 内部境界条件数
  int NoBCouter;      ///< 外部境界条件数
  int NoCompo;        ///< コンポーネント数
  int NoMedium;       ///< 媒質数
  int NoMediumFluid;  ///< 流体の媒質数
  int NoMediumSolid;  ///< 固体の媒質数
  int num_process;    ///< プロセス数
  int num_thread;     ///< スレッド数
  int num_of_polygrp; ///< ポリゴングループ数
  int Parallelism;    ///< 並列モード
  int RefMat;         ///< 参照媒質インデクス
  int Start;
  int SamplingMode;   ///< サンプリング指定
  int NvarsIns_plt3d; ///< Plot3dのバッファサイズ　瞬時値
  int NvarsAvr_plt3d; ///< Plot3dのバッファサイズ　統計値

  unsigned Restart_staging;    ///< リスタート時にリスタートファイルがSTAGINGされているか
  unsigned long NoWallSurface; ///< 固体表面セル数

  double Tscale;

  bool varState[var_END]; ///< 変数のActive/Inactive

  REAL_TYPE BasePrs;
  REAL_TYPE BaseTemp;
  REAL_TYPE HighTemp;
  REAL_TYPE LowTemp;
  REAL_TYPE DiffTemp;
  REAL_TYPE Domain_p1;
  REAL_TYPE Domain_p2;
  REAL_TYPE Gravity;
  REAL_TYPE Grashof;
	REAL_TYPE GridVel[3];
  REAL_TYPE H_Dface[NOFACE];
  REAL_TYPE Mach;
  REAL_TYPE OpenDomain[NOFACE];
  REAL_TYPE Peclet;
  REAL_TYPE PlotIntvl;
  REAL_TYPE Prandtl;
  REAL_TYPE Q_Dface[NOFACE];
  REAL_TYPE Rayleigh;
  REAL_TYPE RefDensity;
  REAL_TYPE RefLength;
  REAL_TYPE RefSoundSpeed;
  REAL_TYPE RefSpecificHeat;
  REAL_TYPE RefKviscosity;
  REAL_TYPE RefVelocity;
  REAL_TYPE Reynolds;
  REAL_TYPE SpecificHeatRatio;
  REAL_TYPE timeflag;
  REAL_TYPE V_Dface[NOFACE];

  // struct
  Mode_set          Mode;
  Initial_Value     iv;
  LES_Parameter     LES;
  Hidden_Parameter  Hide;
  Unit_Def          Unit;
  Ens_of_Compo      EnsCompo;
  Driver_Def        drv;
  Stability_Control Stab;
  AXB_param         axb;

  // class
  IntervalManager Interval[tg_END];  ///< タイミング制御
  IterationCtl* Criteria;            ///< 反復解法の収束判定パラメータ

  string RefMedium;      ///< 参照媒質名 -> int RefMat
  string OperatorName;

  string ver_TP;   ///< TextPerser version no.
  string ver_CPM;  ///< CPMlib
  string ver_CDM;  ///< CDMlib
  string ver_PM;   ///< PMlib
  string ver_Poly; ///< Polylib


  /** コンストラクタ */
  Control(){
    AlgorithmF = 0;
    AlgorithmH = 0;
    BasicEqs = 0;
    CheckParam = 0;
    CnvScheme = 0;
    GuideOut = 0;
    KindOfSolver = 0;
    Limiter = 0;
    MarchingScheme = 0;
    NoBaseLS = 0;
    NoBC = 0;
    NoBCinner = 0;
    NoBCouter = 0;
    NoCompo = 0;
    NoMedium = 0;
    NoMediumFluid = 0;
    NoMediumSolid = 0;
    NoWallSurface = 0;
    num_process = 0;
    num_thread = 0;
    Parallelism = 0;
    RefMat = -1;
    Restart_staging = 0;
    Start = 0;
    NvarsIns_plt3d = 0;
    NvarsAvr_plt3d = 0;

    PlotIntvl = 0.0;
    Domain_p1 = Domain_p2 = 0.0;
    RefVelocity = RefLength = RefDensity = RefSoundSpeed = RefSpecificHeat = RefKviscosity = 0.0;
    DiffTemp = BaseTemp = HighTemp  = LowTemp = 0.0;
    BasePrs = 0.0;
    Gravity = Mach = SpecificHeatRatio =  0.0;
    timeflag = 0.0;
    Tscale = 0.0;

    for (int i=0; i<NOFACE; i++) {
      OpenDomain[i]=0.0;
      Q_Dface[i]=0.0;
      V_Dface[i]=0.0;
      H_Dface[i]=0.0;
    }
    for (int i=0; i<3; i++) {
			GridVel[i]=0.0;
    }
    Reynolds = Peclet = 1.0;
    Prandtl = Rayleigh = Grashof = 0.0;

    for (int i=0; i<var_END; i++) varState[i]=false;

    iv.Density     = 0.0;
    iv.Energy      = 0.0;
    iv.Pressure    = 0.0;
    
    axb.func = 0;
    axb.interval = -1;
    axb.threshold = -1;

    Mode.Statistic = 0;
    Mode.StatisticRestart = 0;
    Mode.Base_Medium = 0;
    Mode.Buoyancy = 0;
    Mode.Example = 0;
    Mode.Log_Base = 0;
    Mode.Log_Wall = 0;
    Mode.PDE = 0;
    Mode.Precision = 0;
    Mode.Profiling = 0;
    Mode.PrsNeuamnnType = 0;
    Mode.ShapeAprx = 0;
    Mode.Steady = 0;
    Mode.CCNV = 0;
    Mode.StatVelocity = 0;
    Mode.StatPressure = 0;
    Mode.StatTemperature = 0;
    Mode.ReynoldsStress = 0;
    Mode.ChannelOutputIter = 0;
    Mode.ChannelOutputMean = 0;
    Mode.ParticleTracking = 0;

    LES.Calc=0;
    LES.Model=0;
    LES.InitialPerturbation = 0;
    LES.Cs = 0.0;
    LES.damping_factor = 0.0;
    LES.ChannelDir = 0;
    LES.BulkVelocity = 0.0;
    LES.ChannelWidth = 0.0;
    LES.TurbulentReynoldsNum = 0.0;
    LES.VelocityProfile = 0;


    Hide.Range_Limit = 0;
    Hide.PM_Test = 0;
    Hide.GeomOutput = OFF;
    Hide.GlyphOutput = OFF;
    Hide.DryRun = OFF;


    Stab.control = OFF;
    Stab.begin = 0.0;
    Stab.end   = 0.0;
    Stab.penalty_number = 0.0;

    Unit.Param  = 0;
    Unit.Output = 0;
    Unit.Prs    = 0;
    Unit.Log    = 0;
    Unit.File   = 0;

    SamplingMode = 0;

    EnsCompo.hsrc    = 0;
    EnsCompo.forcing = 0;
    EnsCompo.outflow = 0;
    EnsCompo.periodic= 0;
    EnsCompo.fraction= 0;
    EnsCompo.tfree   = 0;
    EnsCompo.vspec   = 0;
    EnsCompo.monitor = 0;
    EnsCompo.obstacle= 0;

    drv.length = 0.0;

    Criteria = NULL;
  }

  /**　デストラクタ */
  virtual ~Control()
  {
    if (Criteria) delete [] Criteria;
  }

protected:

  // 熱交換器パラメータの変換（Pa）
  void convertHexCoef(REAL_TYPE* cf);


  // 熱交換器パラメータの変換（水と水銀）
  void convertHexCoef(REAL_TYPE* cf, const REAL_TYPE DensityMode);


  // アプリケーションのパラメータを取得する
  void getApplicationControl();


  // 無次元パラメータを各種モードに応じて設定する
  void getDimensionlessParameter();


  // ドライバー情報を取得
  void getDriver();


  // 初期値擾乱
  void getInitialPerturbation();


  // ログ出力のパラメータを取得
  void getLog();


  // 参照パラメータを取得
  void getReference();


  // ソルバーの種類を特定するパラメータを取得
  void getShapeApproximation();


  // ソルバーの種類を特定するパラメータを取得し，ガイドセルの値を決定する
  void getSolverProperties ();


  // 時間制御に関するパラメータを取得する
  void getTimeControl(DTcntl* DT);


  // 乱流計算のオプションを取得する
  void getTurbulenceModel();


  // 入力ファイルに記述するパラメータとファイルの有次元・無次元の指定を取得する
  void getUnit();




public:

  /**
   * @brief 反復の収束判定パラメータをコピーする
   * @param [in,out] IC   反復制御用クラスの基底へのポインタ
   * @param [in]     name ラベル
   */
  void copyCriteria(IterationCtl* IC, const string name);


  // @brief 反復関連の情報を取得する
  void getIteration();


  // @brief 解法アルゴリズムを選択する
  void getSolvingMethod4Flow();


  // 制御，計算パラメータ群の取得
  void get1stParameter(DTcntl* DT);


  // 制御，計算パラメータ群の取得
  void get2ndParameter();


  /**
   * @brief 計算内部領域の全セル数を返す
   * @param [in] G_size 計算領域全体の分割数
   */
  REAL_TYPE getCellSize(const int* m_size);


  // DryRun parameter
  void getDryRun();


  // 係数行列の書き出し
  void getAXB();
  bool isAXB() {return (axb.func==ON)?true:false;}

  // 粒子追跡のオプション
  void getParticleTracking();

  // @brief 計算モデルの入力ソース情報を取得
  void getGeometryModel();


  /**
   * @brief コンポーネント数，媒質数，境界条件数を取得
   */
  void getNoOfComponent();


  /**
   * @brief ペクレ数の逆数を計算
   * @note Eulerの時にはゼロ
   */
  REAL_TYPE getRcpPeclet() const
  {
    return ( (Mode.PDE == PDE_NS) ? (1.0 / Peclet) : 0.0 );
  }


  /**
   * @brief レイノルズ数の逆数を計算
   * @note Eulerの時にはゼロ
   */
  REAL_TYPE getRcpReynolds() const
  {
    return ( (Mode.PDE == PDE_NS) ? (1.0 / Reynolds) : 0.0 );
  }


  /**
   * @brief 値(REAL_TYPE型)を取得し，返す
   * @param [in] label テストラベル
   * @param [in] tpc   TextParser pointer
   */
  static REAL_TYPE getValueReal(const std::string label, TextParser* tpc);


  /**
   * @brief ベクトル値を取得する
   * @param [in]  label     ラベルディレクトリ
   * @param [out] v         ベクトル値
   * @param [in]  tpc       TextParser pointer
   * @param [in]  normalize trueのとき無次元化
   */
  static bool getVec(const std::string label, REAL_TYPE* v, TextParser* tpc, bool normalize=false);


  /**
   * @brief ベクトル値(2D)を取得する
   * @param [in]  label     ラベルディレクトリ
   * @param [out] v         ベクトル値
   * @param [in]  tpc       TextParser pointer
   */
  static bool getVec2(const std::string label, REAL_TYPE* v, TextParser* tpc);


  /**
   * @brief TPのポインタを受け取る
   * @param [in] tp TextParser
   */
  void importTP(TextParser* tp);


  /**
   * @brief 形状近似がバイナリかどうか
   * @retval Binaryであればtrue
   */
  bool isBinary() const
  {
    return ( BINARY == Mode.ShapeAprx ) ? true : false;
  }


  //@brief 熱問題かどうかを返す
  bool isHeatProblem() const
  {
    return ( ( KindOfSolver != COLD_FLOW ) ? true : false );
  }


  /**
   * @brief 座標値を無次元化する
   * @param [in,out] x 座標値
   */
  void normalizeCord(REAL_TYPE x[3])
  {
    x[0] /= RefLength;
    x[1] /= RefLength;
    x[2] /= RefLength;
  }


  /**
   * @brief 外部境界の各方向の開口率（流体部分の比率）
   * @retval 開口率
   * @param [in] dir    方向
   * @param [in] area   計算外部領域の各面の開口率
   * @param [in] G_size 計算領域全体の分割数
   */
  REAL_TYPE OpenDomainRatio(const int dir, const REAL_TYPE area, const int* G_size);


  // 初期値の表示
  void printInitValues(FILE* fp, CompoList* cmp);


  // 計算パラメータの表示
  void printParaConditions(FILE* fp, const MediumList* mat);


  /*
   * @brief 制御パラメータSTEERの表示
   * @param [in] DT  DTcntl
   */
  void printSteerConditions(FILE* fp, const DTcntl* DT);


  /**
   * @brief 全計算領域の有効セル数と外部境界面の開口率を表示する
   * @param [in] fp 出力ファイルポインタ
   * @param [in] G_Fcell グローバルなFluid cell
   * @param [in] G_Acell グローバルなActive cell
   * @param [in] G_size  グローバルな分割数
   */
  void printOuterArea(FILE* fp, unsigned long G_Fcell, unsigned long G_Acell, int* G_size);



  /**
   * @brief 内部BCコンポーネントの数を表示する
   * @param [in] fp      ファイルポインタ
   */
  void printNoCompo(FILE* fp);


  /**
   * @brief コンポーネントが存在するかを保持しておく
   * @param [in,out] cmp        CompoList
   * @param [in,out] OBC        BoundaryOuter
   * @param [in,out] g_obstacle 各コンポーネント毎のOBSTACLEの有無
   * @retval OBSTACLEの個数
   */
  int setExistComponent(CompoList* cmp, BoundaryOuter* OBC, int* g_obstacle);


  /**
   * @brief コンポーネントと外部境界のパラメータを有次元に設定
   * @param [in,out] mat 媒質配列
   * @param [in,out] cmp コンポーネント配列
   * @param [in]     BO  外部境界条件パラメータ配列
   * @see setRefParameters()
   */
  void setCmpParameters(MediumList* mat, CompoList* cmp, BoundaryOuter* BO);


  /**
   * @brief 無次元パラメータを各種モードに応じて設定する
   * @param [in,out] mat 媒質配列
   * @param [in]     RF  参照フレーム
   * @note
   * - 代表長さと代表速度はパラメータで必ず与えること（読み込んだ値は変更しない）
   * - 純強制対流　有次元　（代表長さ，代表速度，動粘性係数，温度拡散係数）
   * -           無次元　（Pr, Re > RefV=RefL=1）
   * - 熱対流　　　有次元　（代表長さ，代表速度，温度差，体膨張係数，重力加速度，動粘性係数，温度拡散係数）
   * - 自然対流　　有次元　（代表長さ，代表速度，温度差，体膨張係数，重力加速度，動粘性係数，温度拡散係数）
   * - 固体熱伝導　有次元　（代表長さ，温度拡散係数)
   */
  void setRefParameters(MediumList* mat);


  /**
   * @brief 単位ベクトルを計算して戻す
   * @param [in,out] v ベクトル値
   */
  static void UnitVec(REAL_TYPE* v);

};

#endif // _FB_CONTROL_H_
