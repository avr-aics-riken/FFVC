#ifndef _FB_CONTROL_H_
#define _FB_CONTROL_H_

// #################################################################
//
// CAERU Library
//
// Copyright (c) 2012 All right reserved.
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
//
// #################################################################

/**
 @file   Control.h
 @brief  FlowBase Control class Header
 @author kero
 */

#include <math.h>

#include "cpm_Define.h"

#include "DomainInfo.h"
#include "FB_Define.h"
#include "Medium.h"
#include "Component.h"
#include "FBUtility.h"
#include "BndOuter.h"
#include "Interval_Mngr.h"
#include "TPControl.h"

#include "PLOT3D_read.h"
#include "PLOT3D_write.h"

using namespace std;

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
  bool set_Scheme(const char* str, const double val);
  
  
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



class ItrCtl {
private:
  int NormType;          ///< ノルムの種類
  int SubType;           ///< SKIP LOOP or MASK LOOP
  int ItrMax;            ///< 最大反復数
  int LinearSolver;      ///< 線形ソルバーの種類
  int SyncMode;          ///< 同期モード (comm_sync, comm_async)
  REAL_TYPE eps;         ///< 収束閾値
  REAL_TYPE omg;         ///< 加速/緩和係数
  REAL_TYPE NormValue;   ///< ノルムの値
  
public:
  /** 反復制御リスト */
  enum itr_cntl_key 
  {
    ic_prs_pr,
    ic_prs_cr,
    ic_vis_cn,
    ic_tdf_ei,
    ic_tdf_cn,
    ic_END
  };
  
  
  /** 反復法の収束基準種別 */
  enum norm_type 
  { 
    v_div_max,
    v_div_max_dbg,
    v_div_l2,
    v_res_l2_a,
    v_res_l2_r,
    p_res_l2_a,
    p_res_l2_r,
    t_res_l2_a,
    t_res_l2_r
  };
  
  int LoopCount;  ///< 反復回数
  
  /** コンストラクタ */
  ItrCtl() {
    NormType = 0;
    ItrMax = LoopCount = LinearSolver = SubType = 0;
    SyncMode = 0;
    eps = omg = 0.0;
  }
  
  /**　デストラクタ */
  ~ItrCtl() {}
  
  
public:
  /** @brief 線形ソルバの種類を返す */
  int get_LS() const 
  { 
    return LinearSolver; 
  }
  
  /** @brief 最大反復回数を返す */
  int get_ItrMax() const 
  { 
    return ItrMax; 
  }
  
  /** @brief ループ実装の種類を返す */
  int get_LoopType() const 
  { 
    return SubType; 
  }
  
  /** @brief 緩和/加速係数を返す */
  REAL_TYPE get_omg() const 
  { 
    return omg; 
  }
  
  /** @brief 収束閾値を返す */
  REAL_TYPE get_eps() const 
  { 
    return eps; 
  }
  
  /** @brief ノルムのタイプを返す */
  int get_normType() const 
  { 
    return NormType; 
  }
  
  /** @brief keyに対応するノルムの値を返す */
  REAL_TYPE get_normValue() const 
  { 
    return NormValue; 
  }
  
  /** @brief 同期モードを返す */
  int get_SyncMode() const 
  { 
    return SyncMode; 
  }
  
  
  /** @brief 線形ソルバの種類を設定する */
  void set_LS(const int key) 
  { 
    LinearSolver=key; 
  }
  
  /** @brief 最大反復回数を設定する */
  void set_ItrMax(const int key) 
  { 
    ItrMax=key; 
  }
  
  /** @brief ループ実装の種類を設定する */
  void set_LoopType(const int key) 
  { 
    SubType=key; 
  }

  /** @brief 緩和/加速係数を保持 */
  void set_omg(const REAL_TYPE r) 
  { 
    omg = r; 
  }
  
  /** @brief 収束閾値を保持 */
  void set_eps(const REAL_TYPE r) 
  { 
    eps = r; 
  }
  
  /** @brief ノルムのタイプを保持 */
  void set_normType(const int n) 
  { 
    NormType = n; 
  }
  
  /** @brief ノルム値を保持 */
  void set_normValue(const REAL_TYPE r) 
  { 
    NormValue = r; 
  }
  
  /** @brief 同期モードを保持 */
  void set_SyncMode(const Synch_Mode r) 
  { 
    SyncMode= r; 
  }
};



class Control : public DomainInfo {
protected:
  TPControl* tpCntl;   ///< テキストパーサへのポインタ
  
public:
  
  /** 各種モード　パラメータ */
  typedef struct 
  {
    int Average;
    int Buoyancy;
    int Example;
    int Helicity;
    int I2VGT;
    int Log_Base;
    int Log_Itr;
    int Log_Wall;
    int PDE;
    int Precision;
    int Profiling;
    int PrsNeuamnnType;
    int ShapeAprx;
    int Steady;
    int TP;
    int VRT;
    int Wall_profile;
    int Base_Medium;
    int Pshift;
  } Mode_set;
  
  /** 隠しパラメータ */
  typedef struct 
  {
    int Change_ID;
    int Range_Limit;
    int PM_Test;
  } Hidden_Parameter;
  
  /** File IO control */
  typedef struct 
  {
    int IO_Input;
    int IO_Output;
    int Div_Debug;
    int IO_Format; //=0 : *.sph , =1 : PLOT3D
    int IO_Voxel;
  } File_IO_Cntl;
  
  /** PLOT3D オプション */
  typedef struct 
  {
    string basename;
    int IS_xyz;
    int IS_q;
    int IS_funciton;
    int IS_function_name;
    int IS_fvbnd;
    //Initializeでセット
    int ngrid; //出力ブロック数
    int nvar;  //出力項目数
  } Plot3D_Option;
  
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
    REAL_TYPE Temperature;
    REAL_TYPE VecU;
    REAL_TYPE VecV;
    REAL_TYPE VecW;
  } Initial_Value;
  
  /** 単位 */
  typedef struct 
  {
    int Param;  /// 入力パラメータ単位 (Dimensional/NonDimensional)
    int File;   /// ファイルの記述単位 (Dimensional/NonDimensional)
    int Log;    /// 出力ログの単位 (Dimensional/NonDimensional)
    int Prs;    /// 圧力単位 (Absolute/Gauge)
    int Temp;   /// 温度単位 (Celsius/Kelvin)
    int Length; /// 入力パラメータの長さの単位 (non_dimensional/m/cm/mm)
  } Unit_Def;
  
  /** サンプリング機能 */
  typedef struct 
  {
    int log;       /// ログ出力 (ON / OFF)
    int out_mode;  /// 出力モード (Gather / Distribute)
    int unit;      /// 出力単位 (DImensional / NonDimensional)
  } Sampling_Def;
  
  /** コンポーネントの存在 */
  typedef struct 
  {
    int forcing;
    int hsrc;
    int outflow;
    int periodic;
    int fraction;
    int monitor;
  } Ens_of_Compo;
  
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
  
  /** 流れの計算アルゴリズム */
  enum Algorithm_Flow 
  {
    Flow_FS_EE_EE=0,
    Flow_FS_RK_CN,
    Flow_FS_AB2,
    Flow_FS_AB_CN
  };
  
  /** 温度計算アルゴリズム */
  enum Algorithm_Heat 
  {
    Heat_EE_EE=1, // Euler Explicit
    Heat_EE_EI    // Euler Explicit + Implicit
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
  
  
  int AlgorithmF;
  int AlgorithmH;
  int BasicEqs;
  int CheckParam;
  int CnvScheme;
  int FB_version;     ///< FlowBaseクラスのバージョン番号
  int Fill_Hint;      ///< フィルのヒント {no | x_minux | x_plus |...}
  int Fill_Fluid;     ///< 作動流体のフィル媒質
  int Fill_Solid;     ///< 固体のフィル媒質
  int GuideOut;
  int KindOfSolver;
  int Limiter;
  int MarchingScheme;
  int NoBC;
  int NoCompo;
  int NoMedium;       ///< 媒質数
  int NoMediumFluid;
  int NoMediumSolid;
  int num_process;
  int num_thread;
  int Parallelism;
  int RefMat;     ///< 参照媒質インデクス
  int Start;
  int version;    ///< FFVバージョン番号
  
  
  unsigned Restart_step;       ///< リスタートステップ
  unsigned long NoWallSurface; ///< 固体表面セル数
  
  float Scaling_Factor;  ///< ポリゴンのスケーリングファクター
	
  REAL_TYPE BasePrs;
  REAL_TYPE BaseTemp;
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
  REAL_TYPE RefKviscosity;
  REAL_TYPE RefLambda;
  REAL_TYPE RefLength;
  REAL_TYPE RefSoundSpeed;
  REAL_TYPE RefSpecificHeat;
  REAL_TYPE RefVelocity;
  REAL_TYPE RefViscosity;
  REAL_TYPE Reynolds;
  REAL_TYPE SpecificHeatRatio;
  REAL_TYPE timeflag;
  REAL_TYPE Tscale;
  REAL_TYPE V_Dface[NOFACE];
  
  // struct
  Mode_set          Mode;
  Initial_Value     iv;
  LES_Parameter     LES;
  File_IO_Cntl      FIO;
  Plot3D_Option     P3Op;
  Hidden_Parameter  Hide;
  Unit_Def          Unit;
  Sampling_Def      Sampling;
  Ens_of_Compo      EnsCompo;
  
  // class
  Interval_Manager  Interval[Interval_Manager::tg_END];
  
  string HistoryName;
  string HistoryCompoName;
  string HistoryDomfxName;
  string HistoryItrName;
  string HistoryMonitorName;
  string HistoryWallName;
  string PolylibConfigName;
  string HistoryForceName;
  
  string f_Coarse_pressure;
  string f_Coarse_velocity;
  string f_Coarse_temperature;
  string f_Coarse_dfi_prs;
  string f_Coarse_dfi_vel;
  string f_Coarse_dfi_temp;
  
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
  
  string Ref_Medium;
  string OperatorName;
  
  
  /** コンストラクタ */
  Control(){
    AlgorithmF = 0;
    AlgorithmH = 0;
    BasicEqs = 0;
    CheckParam = 0;
    FB_version = 0;
    CnvScheme = 0;
    Fill_Hint = -1;
    Fill_Fluid = 0;
    Fill_Solid = 0;
    GuideOut = 0;
    KindOfSolver = 0;
    Limiter = 0;
    MarchingScheme = 0;
    NoBC = 0;
    NoCompo = 0;
    NoMedium = 0;
    NoMediumFluid = 0;
    NoMediumSolid = 0;
    NoWallSurface = 0;
    num_process = 0;
    num_thread = 0;
    Parallelism = 0;
    RefMat = 0;
    Restart_step = 0;
    Start = 0;
    version = 0;
    
    Scaling_Factor = 1.0;
    
    PlotIntvl = 0.0;
    Domain_p1 = Domain_p2 = 0.0;
    RefVelocity = RefLength = RefDensity = RefSoundSpeed = RefSpecificHeat = RefKviscosity = RefLambda = RefViscosity = 0.0;
    DiffTemp = BaseTemp = BasePrs = 0.0;
    Gravity = Mach = SpecificHeatRatio =  0.0;
    timeflag = Tscale = 0.0;
    
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
    
    iv.Density     = 0.0;
    iv.Energy      = 0.0;
    iv.Pressure    = 0.0;
    iv.Temperature = 0.0;
    
    Mode.Average = 0;
    Mode.Base_Medium = 0;
    Mode.Buoyancy = 0;
    Mode.Example = 0;
    Mode.Helicity = 0;
    Mode.I2VGT = 0;
    Mode.Log_Base = 0;
    Mode.Log_Itr = 0;
    Mode.Log_Wall = 0;
    Mode.PDE = 0;
    Mode.Pshift = -1;
    Mode.Precision = 0;
    Mode.Profiling = 0;
    Mode.PrsNeuamnnType = 0;
    Mode.ShapeAprx = 0;
    Mode.Steady = 0;
    Mode.TP = 0;
    Mode.VRT = 0;
    Mode.Wall_profile = 0;
    
    LES.Calc=0;
    LES.Model=0;
    LES.Cs = 0.0;
    LES.damping_factor=0.0;
    
    FIO.IO_Input  = 0;
    FIO.IO_Output = 0;
    FIO.Div_Debug = 0;
    FIO.IO_Format = 0; //=0 : *.sph , =1 : PLOT3D
    FIO.IO_Voxel  = 0;
    
    P3Op.IS_xyz = ON;
    P3Op.IS_q = OFF;
    P3Op.IS_funciton = ON;
    P3Op.IS_function_name = ON;
    P3Op.IS_fvbnd = OFF;
    
    Hide.Change_ID = 0;
    Hide.Range_Limit = 0;
    Hide.PM_Test = 0;
    
    Unit.Param = 0;
    Unit.Prs   = 0;
    Unit.Log   = 0;
    Unit.File  = 0;
    Unit.Temp  = 0;
    
    Sampling.log = 0;
    Sampling.out_mode = 0;
    Sampling.unit = 0;
    
    EnsCompo.hsrc    = 0;
    EnsCompo.forcing = 0;
    EnsCompo.outflow = 0;
    EnsCompo.periodic= 0;
    EnsCompo.fraction= 0;
    EnsCompo.monitor = 0;
  }
  
  /**　デストラクタ */
  virtual ~Control() {}
  
protected:
  
  /**
   * @brief 熱交換器パラメータの変換（Pa）
   * @param [out] cf パラメータ値
   */
  void convertHexCoef(REAL_TYPE* cf);
  
  
  /**
   * @brief 熱交換器パラメータの変換（水と水銀）
   * @param [out] cf      パラメータ値
   * @param [in]  Density ヘッドの単位
   */
  void convertHexCoef(REAL_TYPE* cf, const REAL_TYPE DensityMode);
  
  
  /**
   * @brief 反復の収束判定パラメータを取得
   * @param [in]     label1  Nodeのラベル1
   * @param [in]     label2  Nodeのラベル2
   * @param [in]     order   ItrCtl配列の格納番号
   * @param [in,out] IC      反復制御用クラスの配列
   */
  void findCriteria(const string label1, const string label2, const int order, ItrCtl* IC);
  
  
  /** 
   * @brief 解法アルゴリズムを選択する 
   */
  void get_Algorithm();
  
  
  /**
   * @brief 平均値操作に関するパラメータを取得する
   * @note パラメータは，setParameters()で無次元して保持
   */
  void get_Average_option();
  
  
  /**
   * @brief Cell IDのゼロを指定IDに変更するオプションを取得する（隠しパラメータ）
   * @note 'Change_ID'の文字列チェックはしないので注意して使うこと
   */
  void get_ChangeID();
  
  
  /** 
   * @brief パラメータ入力チェックモードの取得 
   */
  void get_CheckParameter();
  
  
  /**
   * @brief 対流項スキームのパラメータを取得する
   */
  void get_Convection();
  
  
  /**
   * @brief 派生して計算する変数のオプションを取得する
   */
  void get_Derived();
  
  /**
   * @brief ファイル入出力に関するパラメータを取得し，sphフォーマットの出力の並列モードを指定する．
   * @note インターバルパラメータは，setParameters()で無次元して保持
   */
  void get_FileIO();
  
  
  /**
   * @brief 反復関連の情報を取得する
   */
  void get_Iteration(ItrCtl* IC);
  
  
  /**
   * @brief LES計算のオプションを取得する
   */
  void get_LES_option();
  
  
  
  /**
   * @brief PLOT3Dファイル入出力に関するパラメータ
   * @param [in]  FP3DR    PLOT3D READクラス ポインタ
   * @param [in]  FP3DW    PLOT3D WRITEクラス ポインタ
   */
  void get_PLOT3D(FileIO_PLOT3D_READ* FP3DR, FileIO_PLOT3D_WRITE* FP3DW);

  
  
  /**
   * @brief ソルバーの計算対象種別と浮力モードを取得
   */
  void get_KindOfSolver();
  

  void get_Log            ();
  
  /**
   * @brief 作業者情報を取得
   */
  void get_Operator();
  
  
  /**
   * @brief 無次元パラメータを各種モードに応じて設定する
   * @note
   *   - 純強制対流　有次元　（代表長さ，代表速度，動粘性係数，温度拡散係数）
   *   -           無次元　（Pr, Re > RefV=RefL=1）
   * @see bool setParameters(MediumList* mat, CompoList* cmp)
   */
  void get_Para_ND();
  
  
  /**
   * @brief 参照パラメータを取得
   * @note Ref_IDで指定される媒質を代表物性値とする
   */
  void get_Para_Ref();
  
  /**
   * @brief 温度の参照パラメータを取得
   */
  void get_Para_Temp();
  
  
  /**
   * @brief 性能試験モードを取得する（隠しパラメータ）
   */
  void get_PMtest();
  
  
  /**
   * @brief 参照座標系を取得する
   * @todo 回転は未
   */
  void get_ReferenceFrame (ReferenceFrame* RF);
  
  
  
  
  /** ソルバーの種類を特定するパラメータを取得し，ガイドセルの値を決定する */
  void get_Solver_Properties ();
  
  
  /**
   * @brief 初期値とリスタート条件
   */
  void get_start_condition();
  
  /**
   * @brief 時間制御に関するパラメータを取得する
   * @param [out] DT DTcntlクラス
   * @note パラメータは，setParameters()で無次元して保持
   */
  void get_Time_Control(DTcntl* DT);
  
  
  /** 入力ファイルに記述するパラメータとファイルの有次元・無次元の指定を取得する */
  void get_Unit();
  
  
  /**
   * @brief 変数の範囲制限モードを取得
   * @note 隠しパラメータ
   */
  void get_VarRange();
  
  
  /**
   * @brief 壁面上の扱いを指定する
   */
  void get_Wall_type();
  

  /*
   * @brief 初期値の表示
   * @note Init*には無次元値が保持されている
   * @see void Control::setInitialConditions(void)
   */
  void printInitValues(FILE* fp);
  
  
  /**
   * @brief 線形ソルバー種別の表示
   * @param [in] fp ファイルポインタ
   * @param [in] IC ItrCtl
   */
  void printLS(FILE* fp, const ItrCtl* IC);
  
  
  /** 
   * @brief 計算パラメータの表示
   * @param [in] fp  ファイルポインタ
   * @param [in] mat MediumListクラス
   */
  void printParaConditions(FILE* fp, const MediumList* mat);
  
  
  /** 
   * @brief 制御パラメータSTEERの表示
   * @param [in] IC ItrCtl
   * @param [in] DT DTcntl
   * @param [in] RF ReferenceFrame
   */
  void printSteerConditions(FILE* fp, const ItrCtl* IC, const DTcntl* DT, const ReferenceFrame* RF);

  
public:
  
  /** MediumList中に登録されているkeyに対するIDを返す
   * @param [in] mat  MediumListクラス
   * @param [in] Namx リストの最大数
   * @param [in] key  探査するラベル
   * @return keyに対するIDを返す。発見できない場合はzero
   */
  int find_ID_from_Label(MediumList* mat, const int Nmax, const std::string key);
  
  
  /**
   * @brief 計算内部領域の全セル数を返す
   * @param [in] G_size 計算領域全体の分割数
   */
  REAL_TYPE getCellSize(const int* m_size);
  
  
  /**
   * @brief モデルの形状情報
   * @param [in] MTITP 媒質情報テーブル
   */
  void get_Geometry(const MediumTableInfo *MTITP);
  
  
  /**
   * @brief 外部境界の各方向の開口率（流体部分の比率）
   * @retval 開口率
   * @param [in] dir    方向
   * @param [in] area   計算外部領域の各面の開口率
   * @param [in] G_size 計算領域全体の分割数
   */
  REAL_TYPE OpenDomainRatio(const int dir, const REAL_TYPE area, const int* G_size);
	
  
  /**
   * @brief 制御，計算パラメータ群の表示
   * @param [in] mp  ファイルポインタ（標準出力）
   * @param [in] fp  ファイルポインタ（ファイル出力）
   * @param [in] IC  ItrCtl
   * @param [in] DT  DTcntl
   * @param [in] RF  ReferenceFrame
   * @param [in] mat MediumList
   */
  void displayParams(FILE* mp, FILE* fp, ItrCtl* IC, DTcntl* DT, ReferenceFrame* RF, MediumList* mat);

  
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
   * @brief 無次元パラメータを各種モードに応じて設定する
   * @param [in,out] mat
   * @param [in,out] cmp
   * @param [in,out] RF
   * @param [in,out] BO
   * @note
   * - 代表長さと代表速度はパラメータで必ず与えること（読み込んだ値は変更しない）
   * - 純強制対流　有次元　（代表長さ，代表速度，動粘性係数，温度拡散係数）
   * -           無次元　（Pr, Re > RefV=RefL=1）
   * - 熱対流　　　有次元　（代表長さ，代表速度，温度差，体膨張係数，重力加速度，動粘性係数，温度拡散係数）
   * - 自然対流　　有次元　（代表長さ，代表速度，温度差，体膨張係数，重力加速度，動粘性係数，温度拡散係数）
   * - 固体熱伝導　有次元　（代表長さ，温度拡散係数 > Peclet=1）？
   * @see
   *  - bool Control::getXML_Para_ND(void)
   *  - void Control::getXML_Para_Init(void)
   */
  void setParameters(MediumList* mat, CompoList* cmp, ReferenceFrame* RF, BoundaryOuter* BO);
  
  
  /**
   * @brief ノルムのタイプを返す
   * @param [in] d ノルムの種類
   * @retval ノルムのラベル
   */
  string getNormString(const int d);
  
  
  /**
   * @brief labelのコンポーネント数を返す
   * @param [in] cmp   CompoListクラス
   * @param [in] label コンポーネントID
   */
  int countCompo(CompoList* cmp, const int label);
  

  //@brief Forcingコンポーネントがあれば1を返す
  int isForcing() const 
  {
    return EnsCompo.forcing;
  }
  

  //@brief Hsrcコンポーネントがあれば1を返す
  int isHsrc() const 
  {
    return EnsCompo.hsrc;
  }
  

  //@brief 部分周期境界コンポーネントがあれば1を返す
  int isPeriodic() const 
  {
    return EnsCompo.periodic;
  }


  //@brief モニタコンポーネントがあれば1を返す
  int isMonitor() const 
  {
    return EnsCompo.monitor;
  }
  

  //@brief 流出コンポーネントがあれば1を返す
  int isOutflow() const 
  {
    return EnsCompo.outflow;
  }
  

  //@brief 体積率コンポーネントがあれば1を返す
  int isVfraction() const 
  {
    return EnsCompo.fraction;
  }
  

  //@brief 熱問題かどうかを返す
  bool isHeatProblem() const 
  {
    return ( ( KindOfSolver != FLOW_ONLY ) ? true : false );
  }
  

  /**
   * @brief ソルバーがCDSタイプかどうかを返す
   * @retval CDSであればtrue
   */
  bool isCDS() const 
  {
    return ( CUT_INFO == Mode.ShapeAprx ) ? true : false;
  }
  

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
   * @brief TPのポインタを受け取る
   * @param [in] tp TPControl
   */
  void importTP(TPControl* tp);
  
  
  /**  モニタリングのON/OFFとセルモニタの有無のみを取得  */
  void get_Sampling();
  
  
  /**
   * @brief 制御，計算パラメータ群の取得
   * @param [in] DT     DTctlクラス ポインタ
   * @param [in] FP3DR  PLOT3D READクラス ポインタ
   * @param [in] FP3DW  PLOT3D WRITEクラス ポインタ
   * @note 他のパラメータ取得に先んじて処理しておくもの
   */
  void get_Steer_1(DTcntl* DT, FileIO_PLOT3D_READ* FP3DR, FileIO_PLOT3D_WRITE* FP3DW);
  

  /**
   * @brief 制御，計算パラメータ群の取得
   * @param [in] IC  反復制御クラス
   * @param [in] RF  ReferenceFrameクラス
   */
  void get_Steer_2(ItrCtl* IC, ReferenceFrame* RF);
  
  
  /** バージョン番号の取得 */
  void get_Version();
 

};

#endif // _FB_CONTROL_H_
