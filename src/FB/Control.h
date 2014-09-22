#ifndef _FB_CONTROL_H_
#define _FB_CONTROL_H_

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
#include "CompoFraction.h"
#include "IterationControl.h"
#include "TextParser.h"
#include "IntervalManager.h"

using namespace std;



// #################################################################
class PolygonProperty {
  
private:
  int l_ntria;       ///< ローカルなポリゴン数
  int g_ntria;       ///< グローバルなポリゴン数
  int m_id;          ///< ID
  REAL_TYPE l_area;  ///< ローカルな面積
  REAL_TYPE g_area;  ///< グローバルな面積
  string group;      ///< ポリゴングループ名
  string material;   ///< Mediumtable[@]のalias
  string bc;         ///< BCのラベル
  Vec3<REAL_TYPE> bx_min;      ///< Bboxの最小値
  Vec3<REAL_TYPE> bx_max;      ///< Bboxの最大値
  
public:
  PolygonProperty() {
    l_area = 0.0;
    l_ntria= 0;
    g_area = 0.0;
    g_ntria= 0;
  }
  
  ~PolygonProperty() {}
  
  string getGroup() const
  {
    return group;
  }
  
  int getID() const
  {
    return m_id;
  }
  
  string getMaterial() const
  {
    return material;
  }
  
  string getBClabel() const
  {
    return bc;
  }
  
  int getLntria() const
  {
    return l_ntria;
  }
  
  REAL_TYPE getLarea() const
  {
    return l_area;
  }
  
  int getGntria() const
  {
    return g_ntria;
  }
  
  REAL_TYPE getGarea() const
  {
    return g_area;
  }
  
  void setGroup(string key)
  {
    group = key;
  }
  
  void setID(int key)
  {
    m_id = key;
  }
  
  void setMaterial(string key)
  {
    material = key;
  }
  void setBClabel(string key)
  {
    bc = key;
  }
  
  void setLntria(int val)
  {
    l_ntria = val;
  }
  
  void setLarea(REAL_TYPE val)
  {
    l_area = val;
  }
  
  void setGntria(int val)
  {
    g_ntria = val;
  }
  
  void setGarea(REAL_TYPE val)
  {
    g_area = val;
  }
  
  Vec3<REAL_TYPE> getBboxMax() const
  {
    return bx_max;
  }
  
  Vec3<REAL_TYPE> getBboxMin() const
  {
    return bx_min;
  }
  
  void setBboxMax(const Vec3<REAL_TYPE> bmax)
  {
    bx_max = bmax;
  }
  
  void setBboxMin(const Vec3<REAL_TYPE> bmin)
  {
    bx_min = bmin;
  }
  
};



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
  double dh;       ///< 格子幅（無次元）
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
    dh = Reynolds = Peclet = 0.0;
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
    return (dh*CFL / v);
  }


  /** 拡散数のdt制限
   * @param [in] coef 係数 (Reynolds number or Peclet number)
   */
  double dtDFN(const double coef) const 
  { 
    return coef * dh*dh/6.0; 
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
   * @param [in] m_kos  ソルバーの種類
   * @param [in] m_mode 次元モード
   * @param [in] m_dh   無次元格子幅
   * @param [in] re     レイノルズ数
   * @param [in] pe     ペクレ数
   */
  void set_Vars(const int m_kos, const int m_mode, const double m_dh, const double re, const double pe);
};


// #################################################################
class ReferenceFrame {
  
protected:
  int Frame;         ///< 参照座標系
  double TimeAccel;  ///< 加速時間（無次元）
  double v00[4];     ///< 参照速度（無次元
  double GridVel[3]; ///< 座標系の移動速度（無次元）
  
public:
  /** 参照系の定義 */
  enum frame_type 
  {
    frm_static,
    frm_translation,
    frm_rotation
  };
  
  /** コンストラクタ */
  ReferenceFrame() 
  {
    Frame     = 0;
    TimeAccel = 0.0;
    v00[0] = v00[1] = v00[2] = v00[3] = 0.0;
    GridVel[0] = GridVel[1] = GridVel[2] = 0.0;
  }
  
  /**　デストラクタ */
  ~ReferenceFrame() {}
  
  
  /** 
   * @brief 加速時間をセット
   * @param [in] m_timeAccel 無次元の加速時間
   */
  void setAccel(const double m_timeAccel);
  
  
  /**
   * @brief 参照フレームの種類をセットする
   * @param [in] m_frame 参照フレームの種類
   */
  void setFrame(const int m_frame);
  
  
  /**
   * @brief 格子速度成分の単位方向ベクトルをセットする
   * @param [in] m_Gvel 格子速度成分の単位方向ベクトル
   */
  void setGridVel(const double* m_Gvel);
  
  
  /**
   * @brief 参照速度を設定する
   * @param [in] time 時刻（無次元）
   * @param [in] init フラグ
   */
  void setV00(const double time, const bool init=false);

  
  /**
   * @brief Frameを返す
   */
  int getFrame() const 
  {
    return Frame;
  }
  

  /**
   * @brief 加速時間を返す
   */
  double getAccel() const 
  {
    return TimeAccel;
  }
  

  /** 
   * @brief v00をコピーする
   * @param [out] m_v0 コピー元
   */
  void copyV00(double* m_v0)
  {
    for (int i=0; i<4; i++) m_v0[i] = v00[i];
  }
  

  /**
   * @brief GridVelocityをコピーする
   * @param [out] m_gv 格子速度
   */
  void copyGridVel(double* m_gv)
  {
    for (int i=0; i<3; i++) m_gv[i] = GridVel[i];
  }
};



// #################################################################
class Control : public DomainInfo {
protected:
  TextParser* tpCntl;   ///< テキストパーサへのポインタ
  
public:
  
  /** 各種モード　パラメータ */
  typedef struct 
  {
    int Average;
    int AverageRestart;
    int Buoyancy;
    int Example;
    int Log_Base;
    int Log_Itr;
    int Log_Wall;
    int PDE;
    int Precision;
    int Profiling;
    int PrsNeuamnnType;
    int ShapeAprx;
    int Steady;
    int Base_Medium;
    int CCNV;
  } Mode_set;
  
  /** 隠しパラメータ */
  typedef struct 
  {
    int Range_Limit;
    int PM_Test;
    int GeomOutput;
    int GlyphOutput;
    int Subdivision;
  } Hidden_Parameter;
  
  /** File IO control */
  typedef struct 
  {
    int IOmode;
    int Div_Debug;
    int IO_Voxel;
    int Format;
    int Slice;
    string OutDirPath;
    string InDirPath;
  } File_IO_Cntl;
  
  
  /** LESパラメータ */
  typedef struct 
  {
    int Calc;
    int Model;
    REAL_TYPE Cs;
    REAL_TYPE damping_factor;
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
    Smagorinsky=1,
    Low_Reynolds,
    Dynamic
  };
  
  /** 壁面の扱い */
  enum Wall_Condition 
  {
    No_Slip, // 0
    Slip,    // 1
    Log_Law  // 2
  };
  
  /** ボクセルファイルフォーマット */
  enum Voxel_Type 
  {
    Sphere_SVX=1,
    Sphere_SBX
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
    tg_average,    ///< 平均値出力
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
  int ExperimentNaive; ///< Naive実装のスイッチ
  int FillSeedDir;     ///< フィルのヒント {x_minux | x_plus |...}
  int GuideOut;
  int KindOfSolver;
  int Limiter;
  int MarchingScheme;
  int NoBaseLS;       ///< リストアップされた線形ソルバーの数
  int NoBC;           ///< 境界条件数
  int NoCompo;        ///< コンポーネント数
  int NoMedium;       ///< 媒質数
  int NoMediumFluid;  ///< 流体の媒質数
  int NoMediumSolid;  ///< 固体の媒質数
  int num_process;    ///< プロセス数
  int num_thread;     ///< スレッド数
  int num_of_polygrp; ///< ポリゴングループ数
  int Parallelism;    ///< 並列モード
  int FillID;         ///< フィル媒質ID
  int SeedID;         ///< フィルシード媒質ID
  int RefMat;         ///< 参照媒質インデクス
  int Start;
  int SamplingMode;   ///< サンプリング指定
  int FillSuppress[3];///< PeriodicとSymmetricの外部境界面フィル抑制
  
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
  REAL_TYPE RefVelocity;
  REAL_TYPE Reynolds;
  REAL_TYPE SpecificHeatRatio;
  REAL_TYPE timeflag;
  REAL_TYPE V_Dface[NOFACE];
  
  // struct
  Mode_set          Mode;
  Initial_Value     iv;
  LES_Parameter     LES;
  File_IO_Cntl      FIO;
  Hidden_Parameter  Hide;
  Unit_Def          Unit;
  Ens_of_Compo      EnsCompo;
  Driver_Def        drv;
  
  // class
  IntervalManager Interval[tg_END];  ///< タイミング制御
  IterationCtl* Criteria;            ///< 反復解法の収束判定パラメータ
  
  string file_fmt_ext;
  string PolylibConfigName;

  
  // 入力dfiファイルのプレフィックス
  string f_dfi_in_prs;
  string f_dfi_in_vel;
  string f_dfi_in_fvel;
  string f_dfi_in_temp;
  string f_dfi_in_prsa;
  string f_dfi_in_vela;
  string f_dfi_in_tempa;
  
  // 出力ファイルのプレフィックス
  string f_Velocity;
  string f_Pressure;
  string f_Temperature;
  string f_AvrPressure;
  string f_AvrVelocity;
  string f_AvrTemperature;
  string f_DivDebug;
  string f_Helicity;
  string f_TotalP;
  string f_I2VGT;
  string f_Vorticity;
  string f_Fvelocity;
  
  string RefMedium;      ///< 参照媒質名 -> int RefMat
  string FillMedium;     ///< フィルに使う媒質 -> int FillID
  string SeedMedium;     ///< ヒントに使う媒質 -> int SeedID
  string OperatorName;
  
  string ver_TP;   ///< TextPerser version no.
  string ver_CPM;  ///< CPMlib
  string ver_CIO;  ///< CIOlib
  string ver_PM;   ///< PMlib
  string ver_Poly; ///< Polylib
  string ver_CUT;  ///< Cutlib
  
  
  /** コンストラクタ */
  Control(){
    AlgorithmF = 0;
    AlgorithmH = 0;
    BasicEqs = 0;
    CheckParam = 0;
    CnvScheme = 0;
    ExperimentNaive = OFF;
    FillSeedDir = -1;
    GuideOut = 0;
    KindOfSolver = 0;
    Limiter = 0;
    MarchingScheme = 0;
    NoBaseLS = 0;
    NoBC = 0;
    NoCompo = 0;
    NoMedium = 0;
    NoMediumFluid = 0;
    NoMediumSolid = 0;
    NoWallSurface = 0;
    num_process = 0;
    num_thread = 0;
    Parallelism = 0;
    FillID = -1;
    SeedID = -1;
    RefMat = -1;
    Restart_staging = 0;
    Start = 0;
    
    PlotIntvl = 0.0;
    Domain_p1 = Domain_p2 = 0.0;
    RefVelocity = RefLength = RefDensity = RefSoundSpeed = RefSpecificHeat = 0.0;
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
    
    for (int i=0; i<3; i++)  FillSuppress[i] = -1;
    
    iv.Density     = 0.0;
    iv.Energy      = 0.0;
    iv.Pressure    = 0.0;
    
    Mode.Average = 0;
    Mode.AverageRestart = 0;
    Mode.Base_Medium = 0;
    Mode.Buoyancy = 0;
    Mode.Example = 0;
    Mode.Log_Base = 0;
    Mode.Log_Itr = 0;
    Mode.Log_Wall = 0;
    Mode.PDE = 0;
    Mode.Precision = 0;
    Mode.Profiling = 0;
    Mode.PrsNeuamnnType = 0;
    Mode.ShapeAprx = 0;
    Mode.Steady = 0;
    Mode.CCNV = 0;
    
    LES.Calc=0;
    LES.Model=0;
    LES.Cs = 0.0;
    LES.damping_factor=0.0;
    
    FIO.IOmode   = 0;
    FIO.Div_Debug  = 0;
    FIO.IO_Voxel   = 0;
    FIO.Format     = 0; // 0:sph, 1:BOV
    FIO.Slice      = 0;
    
    Hide.Range_Limit = 0;
    Hide.PM_Test = 0;
    Hide.GeomOutput = OFF;
    Hide.GlyphOutput = OFF;
    Hide.Subdivision = 0;
    
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
  
  
  // ファイル入出力に関するパラメータを取得し，sphフォーマットの出力の並列モードを指定する
  void getFieldData();
  
  
  // sphフォーマットのオプションを指定
  void getFormat_sph();
  

  // ログ出力のパラメータを取得
  void getLog();
  
  
  // 参照パラメータを取得
  void getReference();
  
  
  // 参照座標系を取得する
  void getReferenceFrame(ReferenceFrame* RF);
  
  
  // ソルバーの種類を特定するパラメータを取得
  void getShapeApproximation();
  
  
  // ソルバーの種類を特定するパラメータを取得し，ガイドセルの値を決定する
  void getSolverProperties ();
  
  
  // 初期値とリスタート条件
  void getStartCondition();
  
  
  // 時間制御に関するパラメータを取得する
  void getTimeControl(DTcntl* DT);
  
  
  // 乱流計算のオプションを取得する
  void getTurbulenceModel();
  
  
  // 入力ファイルに記述するパラメータとファイルの有次元・無次元の指定を取得する
  void getUnit();
  
  
  // 壁面上の扱いを指定する
  void getWallType();
  
  
  // 初期値の表示
  void printInitValues(FILE* fp, CompoList* cmp);
  
  
  // 線形ソルバー種別の表示
  void printLS(FILE* fp, const IterationCtl* IC);
  
  
  // 計算パラメータの表示
  void printParaConditions(FILE* fp, const MediumList* mat);
  
  
  /**
   * @brief 制御パラメータSTEERの表示
   * @param [in] IC  IterationCtl
   * @param [in] DT  DTcntl
   * @param [in] RF  ReferenceFrame
   * @param [in] em  ffvcの実行モード
   */
  void printSteerConditions(FILE* fp, IterationCtl* IC, const DTcntl* DT, const ReferenceFrame* RF, const int em);
  
  
  
  
public:
  
  /**
   * @brief 反復の収束判定パラメータをコピーする
   * @param [in,out] IC   反復制御用クラスの基底へのポインタ
   * @param [in]     name ラベル
   */
  void copyCriteria(IterationCtl& IC, const string name);
  
  
  /**
   * @brief 制御，計算パラメータ群の表示
   * @param [in] mp   ファイルポインタ（標準出力）
   * @param [in] fp   ファイルポインタ（ファイル出力）
   * @param [in] IC   IterationCtl
   * @param [in] DT   DTcntl
   * @param [in] RF   ReferenceFrame
   * @param [in] mat  MediumList
   * @param [in] cmp  CompoList
   * @param [in] em   ffvcの実行モード
   */
  void displayParams(FILE* mp, FILE* fp,
                     IterationCtl* IC,
                     DTcntl* DT,
                     ReferenceFrame* RF,
                     MediumList* mat,
                     CompoList* cmp,
                     const int em);
  
  
  /** MediumList中に登録されているkeyに対するIDを返す
   * @param [in] mat  MediumListクラス
   * @param [in] Namx リストの最大数
   * @param [in] key  探査するラベル
   * @return keyに対するIDを返す。発見できない場合はzero
   */
  int findIDfromLabel(const MediumList* mat, const int Nmax, const std::string key);
  
  
  // @brief 反復関連の情報を取得する
  void getIteration();
  
  
  // @brief 解法アルゴリズムを選択する
  void getSolvingMethod4Flow();
  
  
  // 制御，計算パラメータ群の取得
  void get1stParameter(DTcntl* DT);
  
  
  // 制御，計算パラメータ群の取得
  void get2ndParameter(ReferenceFrame* RF);
  
  
  /**
   * @brief 計算内部領域の全セル数を返す
   * @param [in] G_size 計算領域全体の分割数
   */
  REAL_TYPE getCellSize(const int* m_size);
  
  
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
    return ( ( KindOfSolver != FLOW_ONLY ) ? true : false );
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
	
  
  /**
   * @brief 全計算領域の有効セル数と外部境界面の開口率を表示する
   * @param [in] fp 出力ファイルポインタ
   * @param [in] G_Fcell グローバルなFluid cell
   * @param [in] G_Acell グローバルなActive cell
   * @param [in] G_size  グローバルな分割数
   */
  void printOuterArea(FILE* fp, unsigned long G_Fcell, unsigned long G_Acell, int* G_size);
  
  
  /**
   * @brief グローバルな領域情報を表示する
   * @param [in] fp      ファイルポインタ
   * @param [in] G_size  グローバルな分割数
   * @param [in] G_org   グローバルな領域基点
   * @param [in] G_reg   グローバルな領域サイズ
   * @param [in] pch     格子幅
   */
  void printGlobalDomain(FILE* fp, const int* G_size, const REAL_TYPE* G_org, const REAL_TYPE* G_reg, const REAL_TYPE* pch);
  
  
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
  void setRefParameters(MediumList* mat, ReferenceFrame* RF);

  
  /**
   * @brief 単位ベクトルを計算して戻す
   * @param [in,out] v ベクトル値
   */
  static void UnitVec(REAL_TYPE* v);
  
};

#endif // _FB_CONTROL_H_
