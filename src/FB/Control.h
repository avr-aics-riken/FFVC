#ifndef _FB_CONTROL_H_
#define _FB_CONTROL_H_

// #################################################################
//
// CAERU Library
//
// Copyright (c) All right reserved. 2012
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
//
// #################################################################

/**
 @file Control.h
 @brief FlowBase Control class Header
 @author kero
 */

#include <math.h>

#include "cpm_Define.h"

#include "FB_Define.h"
#include "Medium.h"
#include "Component.h"
#include "FBUtility.h"
//#include "Monitor.h"
#include "BndOuter.h"
#include "Interval_Mngr.h"
#include "TPControl.h"

using namespace std;

class DTcntl {
public:
  /** 時間積分幅の指定種別 */
  enum dt_Type {
    dt_direct=1,      ///< 入力値がΔt
    dt_cfl_ref_v,     ///< dt < c dx/U0
    dt_cfl_max_v,     ///< dt < c dx/Umax
    dt_dfn,           ///< 拡散数制限
    dt_cfl_dfn_ref_v, ///< dt = min( c dx/U0, diffusion number )
    dt_cfl_dfn_max_v, ///< dt = min( c dx/Umax, diffusion number )
    dt_cfl_max_v_cp   ///< 圧縮性　dt < (cfl+soundSpeed) dx/Umax
  };
  
private:
  int scheme;        ///< Δtのスキーム種別
  int KOS;           ///< Kind of Solver
  int mode;          ///< 入力パラメータの次元モード（無次元/有次元）
  double   CFL;      ///< Δtの決定に使うCFLなど
  double   deltaT;   ///< Δt（無次元）
  double   dh;       ///< 格子幅（無次元）
  double   Reynolds; ///< レイノルズ数
  double   Peclet;   ///< ペクレ数
  
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
   * @param[in] Uref 速度の参照値（無次元）
   */
  double dtCFL(const double Uref) const 
  {
    double v = (Uref < 1.0) ? 1.0 : Uref; // 1.0 is non-dimensional reference velocity
    return (dh*CFL / v);
  }


  /** 拡散数のdt制限
   * @param coef 係数 (Reynolds number or Peclet number)
   */
  double dtDFN(const double coef) const 
  { 
    return coef * dh*dh/6.0; 
  }
  
  
  /** 時間積分幅とKindOfSolver種別の整合性をチェック */
  bool chkDtSelect();
  
  
  /**
   * @brief Δtのスキームを設定する
   * @retval 設定の成否
   * @param str[in] str  キーワード
   * @param val[in] val  値
   */
  bool set_Scheme (const char* str, const double val);
  
  
  /**
   * @brief 各種モードに対応する時間積分幅を設定する
   * @retval return code
   * @param[in] vRef 速度の参照値（無次元）
   * @note deltaTは無次元
   */
  int set_DT(const double vRef);
  
  
  /**
   * @brief 基本変数をコピー
   * @param[in] m_kos  ソルバーの種類
   * @param[in] m_mode 次元モード
   * @param[in] m_dh   無次元格子幅
   * @param[in] re     レイノルズ数
   * @param[in] pe     ペクレ数
   */
  void set_Vars(const unsigned m_kos, const unsigned m_mode, const double m_dh, const double re, const double pe);
};



class ReferenceFrame {
  
protected:
  int Frame;         ///< 参照座標系
  double TimeAccel;  ///< 加速時間（無次元）
  double v00[4];     ///< 参照速度（無次元
  double GridVel[3]; ///< 座標系の移動速度（無次元）
  
public:
  /** 参照系の定義 */
  enum frame_type {
    frm_static,
    frm_translation,
    frm_rotation
  };
  
  /** コンストラクタ */
  ReferenceFrame(){
    Frame     = 0;
    TimeAccel = 0.0;
    v00[0] = v00[1] = v00[2] = v00[3] = 0.0;
    GridVel[0] = GridVel[1] = GridVel[2] = 0.0;
  }
  
  /**　デストラクタ */
  ~ReferenceFrame() {}
  
  
  /** 
   * @brief 加速時間をセット
   * @param[in] m_timeAccel 無次元の加速時間
   */
  void setAccel  (const double m_timeAccel);
  
  
  /**
   @brief 参照フレームの種類をセットする
   @param[in] m_frame 参照フレームの種類
   */
  void setFrame  (const int m_frame);
  
  
  /**
   @brief 格子速度成分の単位方向ベクトルをセットする
   @param[in] m_Gvel 格子速度成分の単位方向ベクトル
   */
  void setGridVel(const double* m_Gvel);
  
  
  /**
   @brief 参照速度を設定する
   @param[in] time 時刻（無次元）
   @param[in] init フラグ
   */
  void setV00(const double time, const bool init=false);

  
  //@brief Frameを返す
  int getFrame() const 
  {
    return Frame;
  }
  

  //@brief 加速時間を返す
  double getAccel() const 
  {
    return TimeAccel;
  }
  

  //@brief v00をコピーする
  void copyV00(double* m_v0)
  {
    for (int i=0; i<4; i++) m_v0[i] = v00[i];
  }
  

  //@brief GridVelocityをコピーする
  void copyGridVel(double* m_gv)
  {
    for (int i=0; i<3; i++) m_gv[i] = GridVel[i];
  }
};



class ItrCtl {
private:
  int NormType;          /// ノルムの種類
  int SubType;           /// SKIP LOOP or MASK LOOP
  int ItrMax;            /// 最大反復数
  int LinearSolver;      /// 線形ソルバーの種類
  REAL_TYPE eps;         /// 収束閾値
  REAL_TYPE omg;         /// 加速/緩和係数
  REAL_TYPE NormValue;   /// ノルムの値
  
public:
  /** 反復制御リスト */
  enum itr_cntl_key {
    ic_prs_pr,
    ic_prs_cr,
    ic_vis_cn,
    ic_tdf_ei,
    ic_tdf_cn,
    ic_END
  };
  
  
  /** 反復法の収束基準種別 */
  enum norm_type { 
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
};



class Control {
protected:
  TPControl* tpCntl;   ///< テキストパーサへのポインタ
  
public:
  
  // domain info
  DomainInfo dInfo;
  SubDomain dom;
  
  // 各種モード　パラメータ
  typedef struct {
    unsigned Average;
    unsigned Buoyancy;
    unsigned Example;
    unsigned Helicity;
    unsigned I2VGT;
    unsigned Log_Base;
    unsigned Log_Itr;
    unsigned Log_Wall;
    unsigned PDE;
    unsigned Precision;
    unsigned Profiling;
    unsigned PrsNeuamnnType;
    unsigned ShapeAprx;
    unsigned Steady;
    unsigned TP;
    unsigned VRT;
    unsigned Wall_profile;
    int Base_Medium;
    int Pshift;
  } Mode_set;
  
  // 隠しパラメータ
  typedef struct {
    unsigned Change_ID;
    unsigned Range_Limit;
    unsigned PM_Test;
    REAL_TYPE Scaling_Factor;
  } Hidden_Parameter;
  
  /// File IO control
  typedef struct {
    unsigned IO_Input;
    unsigned IO_Output;
    unsigned Div_Debug;
  } File_IO_Cntl;
  
  // LESパラメータ
  typedef struct {
    unsigned Calc;
    unsigned Model;
    REAL_TYPE Cs;
    REAL_TYPE damping_factor;
  } LES_Parameter;
  
  // 初期値
  typedef struct {
    REAL_TYPE Density;
    REAL_TYPE Energy;
    REAL_TYPE Pressure;
    REAL_TYPE Temperature;
    REAL_TYPE VecU;
    REAL_TYPE VecV;
    REAL_TYPE VecW;
  } Initial_Value;
  
  // 単位
  typedef struct {
    unsigned Param; /// 入力パラメータ単位 (Dimensional/NonDimensional)
    unsigned File;  /// ファイルの記述単位 (Dimensional/NonDimensional)
    unsigned Log;   /// 出力ログの単位 (Dimensional/NonDimensional)
    unsigned Prs;   /// 圧力単位 (Absolute/Gauge)
    unsigned Temp;  /// 温度単位 (Celsius/Kelvin)
  } Unit_Def;
  
  // サンプリング機能
  typedef struct {
    unsigned log;       /// ログ出力 (ON / OFF)
    unsigned out_mode;  /// 出力モード (Gather / Distribute)
    unsigned unit;      /// 出力単位 (DImensional / NonDimensional)
  } Sampling_Def;
  
  // コンポーネントの存在
  typedef struct {
    unsigned forcing;
    unsigned hsrc;
    unsigned outflow;
    unsigned periodic;
    unsigned fraction;
    unsigned monitor;
  } Ens_of_Compo;
  
  /// 偏微分方程式の型
  enum PDE_type {
    PDE_NS=1,
    PDE_EULER
  };
  
  /// 定常性
  enum Time_Variation {
    TV_Steady=1,
    TV_Unsteady
  };
  
  /// 値域の制御　
  enum Range_Cntl {
    Range_Normal=1,
    Range_Cutoff
  };
  
  /// 出力タイミングの指定
  enum output_mode {
    IO_normal,
    IO_forced,
    IO_everytime
  };
  
  /// 流れの計算アルゴリズム
  enum Algorithm_Flow {
    Flow_FS_EE_EE=0,
    Flow_FS_RK_CN,
    Flow_FS_AB2,
    Flow_FS_AB_CN
  };
  
  /// 温度計算アルゴリズム
  enum Algorithm_Heat {
    Heat_EE_EE=1, // Euler Explicit
    Heat_EE_EI    // Euler Explicit + Implicit
  };
  
  /// 対流項スキーム
  enum convection_scheme {
    O1_upwind=1,
    O2_central,
    O3_muscl,
    O4_central
  };
  
  /// TVD limiter
  enum Limiter_function {
    No_Limiter,
    MinMod
  };
  
  /// LES Model
  enum LES_model {
    Smagorinsky=1,
    Low_Reynolds,
    Dynamic
  };
  
  /// 壁面の扱い
  enum Wall_Condition {
    No_Slip, // 0
    Slip,    // 1
    Log_Law  // 2
  };
  
  /// ボクセルファイルフォーマット
  enum Voxel_Type {
    Sphere_SVX=1,
    Sphere_SBX
  };
  
  /// 並列化モード
  enum Parallel_mode {
    Serial=1,
    OpenMP,
    FlatMPI,
    Hybrid
  };
  
  unsigned long Acell,
                Fcell,
                Wcell;
  
  int NoMedium;   /// 媒質数
  int RefMat;     /// 参照媒質インデクス
  
  unsigned  AlgorithmF,
            AlgorithmH,
            BasicEqs,
            CheckParam,
            CnvScheme,
            FB_version,
            guide,
            GuideOut,
            imax, jmax, kmax,
            KindOfSolver,
            LastStep,
            Limiter,
            MarchingScheme,
            NoBC,
            NoCompo,
            NoMediumFluid,
            NoMediumSolid,
            NoWallSurface,
            num_process,
            num_thread,
            Parallelism,
            Restart_step,
            Start,
            version,
            vxFormat;
	
  REAL_TYPE BasePrs,
            BaseTemp,
            DiffTemp,
            dh,
            Domain_p1,
            Domain_p2,
            Eff_Cell_Ratio,
            dx[3],
            Gravity,
            Grashof,
	          GridVel[3],
            H_Dface[NOFACE],
            Lbx[3],
            Mach,
            OpenDomain[NOFACE],
            org[3],
            Peclet,
            pitch,
            PlotIntvl, 
            Prandtl,
            Q_Dface[NOFACE],
            Rayleigh,
            RefDensity,
            RefKviscosity,
            RefLambda,
            RefLength,
            RefSoundSpeed,
            RefSpecificHeat,
            RefVelocity,
            RefViscosity,
            Reynolds,
            roll,
            SpecificHeatRatio, 
            timeflag,
            Tscale,
            V_Dface[NOFACE],
            yaw;
  
  // struct
  Mode_set          Mode;
  Initial_Value     iv;
  LES_Parameter     LES;
  File_IO_Cntl      FIO;
  Hidden_Parameter  Hide;
  Unit_Def          Unit;
  Sampling_Def      Sampling;
  Ens_of_Compo      EnsCompo;
  
  // class
  Interval_Manager  Interval[Interval_Manager::tg_END];

  long TotalMemory;
  char HistoryName[LABEL], HistoryCompoName[LABEL], HistoryDomfxName[LABEL], HistoryItrName[LABEL], HistoryMonitorName[LABEL];
  char HistoryWallName[LABEL], PolylibConfigName[LABEL], HistoryForceName[LABEL];
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
  
  
  /** コンストラクタ */
  Control(){
    Acell = 0;
    AlgorithmF = 0;
    AlgorithmH = 0;
    BasicEqs = 0;
    CheckParam = 0;
    FB_version = 0;
    Fcell = 0;
    CnvScheme = 0;
    guide = 0;
    GuideOut = 0;
    imax = jmax = kmax = 0;
    KindOfSolver = 0;
    LastStep = 0;
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
    vxFormat = 0;
    Wcell = 0;
    
    dh = 0.0;
    PlotIntvl = 0.0;
    yaw = pitch = roll = 0.0;
    Domain_p1 = Domain_p2 = 0.0;
    RefVelocity = RefLength = RefDensity = RefSoundSpeed = RefSpecificHeat = RefKviscosity = RefLambda = RefViscosity = 0.0;
    Eff_Cell_Ratio = 0.0;
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
      org[i] = 0.0;
      dx[i]   = 0.0;
      Lbx[i] = 0.0;
			GridVel[i]=0.0;
    }
    Reynolds = Peclet = 1.0;
    Prandtl = Rayleigh = Grashof = 0.0;
    
    memset(HistoryName,  0, sizeof(char)*LABEL);
    memset(HistoryCompoName,  0, sizeof(char)*LABEL);
    memset(HistoryDomfxName,  0, sizeof(char)*LABEL);
    memset(HistoryForceName,  0, sizeof(char)*LABEL);
		memset(HistoryItrName,  0, sizeof(char)*LABEL);
    memset(HistoryMonitorName,  0, sizeof(char)*LABEL);
    memset(HistoryWallName,  0, sizeof(char)*LABEL);
    memset(PolylibConfigName,  0, sizeof(char)*LABEL);

    TotalMemory = 0;
    
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
    
    FIO.IO_Input = 0;
    FIO.IO_Output = 0;
    FIO.Div_Debug = 0;
    
    Hide.Change_ID = 0;
    Hide.Range_Limit = 0;
    Hide.PM_Test = 0;
    Hide.Scaling_Factor = 0.0;
    
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
  bool getXML_PrsAverage(void);
  
  void convertHexCoef        (REAL_TYPE* cf);
  void convertHexCoef        (REAL_TYPE* cf, REAL_TYPE DensityMode);
  
  
  /**
   @brief 反復の収束判定パラメータを取得
   @param[in]     label1  Nodeのラベル1
   @param[in]     label2  Nodeのラベル2
   @param[in]     order   ItrCtl配列の格納番号
   @param[in/out] IC      反復制御用クラスの配列
   */
  void findCriteria(const string label1, const string label2, const int order, ItrCtl* IC);
  
  
  /** 解法アルゴリズムを選択する */
  void get_Algorithm();
  
  
  void get_Average_option ();
  void get_ChangeID       ();
  
  
  /** パラメータ入力チェックモードの取得 */
  void get_CheckParameter ();
  
  
  void get_Convection     ();
  void get_Derived        ();
  
  /**
   @brief ファイル入出力に関するパラメータを取得し，sphフォーマットの出力の並列モードを指定する．
   @note インターバルパラメータは，setParameters()で無次元して保持
   */
  void get_FileIO();
  
  void get_Iteration      (ItrCtl* IC);
  void get_KindOfSolver   ();
  void get_LES_option     ();
  void get_Log            ();
  void get_Para_ND        ();
  void get_Para_Ref       ();
  void get_Para_Temp      ();
  void get_PMtest         ();
  void get_ReferenceFrame (ReferenceFrame* RF);
  ////void get_restart_rough  ();
  
  
  /**
   @brief スケーリングファクタを取得する（隠しパラメータ）
   @note 'Scaling_factor'の文字列チェックはしないので注意して使うこと
   */
  void get_Scaling();
  
  
  /** ソルバーの種類を特定するパラメータを取得し，ガイドセルの値を決定する */
  void get_Solver_Properties ();
  
  void get_start_condition();
  
  /**
   @brief 時間制御に関するパラメータを取得する
   @param[out] DT DTcntlクラス
   @note パラメータは，setParameters()で無次元して保持
   */
  void get_Time_Control(DTcntl* DT);
  
  
  /** 入力ファイルに記述するパラメータとファイルの有次元・無次元の指定を取得する */
  void get_Unit           ();
  
  
  void get_VarRange       ();
  void get_Wall_type      ();
  void printArea             (FILE* fp, unsigned G_Fcell, unsigned G_Acell, unsigned* G_size);
  void printVoxelSize        (unsigned* gs, FILE* fp);
  void printInitValues       (FILE* fp);
  void printLS               (FILE* fp, ItrCtl* IC);
  void printParaConditions   (FILE* fp);
  void printSteerConditions  (FILE* fp, ItrCtl* IC, DTcntl* DT, ReferenceFrame* RF);
  void setDimParameters      (void);

public:
  const char* get_VoxelFileName(void);
  
  REAL_TYPE getCellSize(unsigned* m_size);
  REAL_TYPE OpenDomainRatio(const unsigned dir, const REAL_TYPE area, const unsigned* G_size);
	
  void displayParams            (FILE* mp, FILE* fp, ItrCtl* IC, DTcntl* DT, ReferenceFrame* RF);

  void printAreaInfo            (FILE* fp, FILE* mp, unsigned G_Fcell, unsigned G_Acell, unsigned* G_size);
  void printDomain              (FILE* fp, unsigned* G_size, REAL_TYPE* G_org, REAL_TYPE* G_Lbx);
  void printDomainInfo          (FILE* mp, FILE* fp, unsigned* G_size, REAL_TYPE* G_org, REAL_TYPE* G_Lbx);
  void printNoCompo             (FILE* fp);
  void setDomainInfo            (unsigned* m_sz, REAL_TYPE* m_org, REAL_TYPE* m_pch, REAL_TYPE* m_wth);
  void setParameters            (MediumList* mat, CompoList* cmp, ReferenceFrame* RF, BoundaryOuter* BO);
  
  
  /**
   * @brief ノルムのタイプを返す
   * @param[in] d ノルムの種類
   * @retval ノルムのラベル
   */
  string getNormString(const int d);
  
  
  unsigned countCompo  (CompoList* cmp, unsigned label);
  
  //@fn unsigned isForcing(void) const
  //@brief Forcingコンポーネントがあれば1を返す
  unsigned isForcing(void) const {
    return EnsCompo.forcing;
  }
  
  //@fn unsigned isHsrc(void) const
  //@brief Hsrcコンポーネントがあれば1を返す
  unsigned isHsrc(void) const {
    return EnsCompo.hsrc;
  }
  
  //@fn unsigned isPeriodic(void) const
  //@brief 部分周期境界コンポーネントがあれば1を返す
  unsigned isPeriodic(void) const {
    return EnsCompo.periodic;
  }

  //@fn unsigned isMonitor(void) const
  //@brief モニタコンポーネントがあれば1を返す
  unsigned isMonitor(void) const {
    return EnsCompo.monitor;
  }
  
  //@fn unsigned isOutflow(void) const
  //@brief 流出コンポーネントがあれば1を返す
  unsigned isOutflow(void) const {
    return EnsCompo.outflow;
  }
  
  //@fn unsigned isVfraction(void) const
  //@brief 体積率コンポーネントがあれば1を返す
  unsigned isVfraction(void) const {
    return EnsCompo.fraction;
  }
  
  //@fn bool isHeatProblem(void) const
  //@brief 熱問題かどうかを返す
  bool isHeatProblem(void) const {
    return ( ( KindOfSolver != FLOW_ONLY ) ? true : false );
  }
  
  //@fn bool isCDS(void) const
  //@brief ソルバーがCDSタイプかどうかを返す
  //@retval CDSであればtrue
  bool isCDS(void) const {
    return ( CUT_INFO == Mode.ShapeAprx ) ? true : false;
  }
  
  //@fn REAL_TYPE getRcpPeclet(void) const
  //@brief ペクレ数の逆数を計算
  REAL_TYPE getRcpPeclet(void) const {
    return ( 1.0 / Peclet );
  }
  
  //@fn REAL_TYPE getRcpReynolds(void) const
  //@brief レイノルズ数の逆数を計算
  REAL_TYPE getRcpReynolds(void) const {
    return ( (Mode.PDE == PDE_NS) ? (1.0 / Reynolds) : 0.0 );
  }

  //@fn void normalizeCord(REAL_TYPE x[3])
  //@brief 座標値を無次元化する
  void normalizeCord(REAL_TYPE x[3]) {
    x[0] /= RefLength;
    x[1] /= RefLength;
    x[2] /= RefLength;
  }
  
  
public:
  bool set_DomainInfo   (unsigned* size,
                         REAL_TYPE* origin,
                         REAL_TYPE* pitch,
                         REAL_TYPE* width);
  bool get_DomainInfo   (unsigned* size);
  bool get_SubDomainInfo(unsigned* size);
  
  
  /**
   * @brief TPのポインタを受け取る
   */
  bool importTP(TPControl* tp);
  
  void get_Para_Init    (void);
  void get_Polygon      (void);
  
  
  /**  モニタリングのON/OFFとセルモニタの有無のみを取得  */
  void get_Sampling();
  
  void get_Steer_1      (DTcntl* DT);
  void get_Steer_2      (ItrCtl* IC, ReferenceFrame* RF);
  void get_Version      (void);
  
  //for text parser
protected:
  //bool get_PrsAverage(void);//未使用ルーチン
  ////未使用ルーチンget_Polygonからのみコールされている
  ////const CfgElem* get_Pointer(const char* key, string section);
 

};

#endif // _FB_CONTROL_H_
