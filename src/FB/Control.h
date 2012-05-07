#ifndef _SKL_FB_CONTROL_H_
#define _SKL_FB_CONTROL_H_

/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

/**
 @file Control.h
 @brief FlowBase Control class Header
 @author keno, FSI Team, VCAD, RIKEN
 */

#include <math.h>

#include "FBDefine.h"
#include "Skl.h"
#include "SklSolverBase.h"
#include "config/SklSolverConfig.h"
#include "SklUtil.h"
#include "Material.h"
#include "Component.h"
#include "FBUtility.h"
#include "Monitor.h"
#include "BndOuter.h"
#include "Interval_Mngr.h"

using namespace SklCfg;  // to use SklSolverConfig* cfg


class DTcntl {
public:
  enum dt_Type {
    dt_direct=1,      // 入力値がΔt
    dt_cfl_ref_v,     // dt < c dx/U0
    dt_cfl_max_v,     // dt < c dx/Umax
    dt_dfn,           // 拡散数制限
    dt_cfl_dfn_ref_v, // dt = min( c dx/U0, diffusion number )
    dt_cfl_dfn_max_v, // dt = min( c dx/Umax, diffusion number )
    dt_cfl_max_v_cp   // 圧縮性　dt < (cfl+soundSpeed) dx/Umax
  };
  
private:
  unsigned scheme;  // Δtのスキーム種別
  unsigned KOS;     // Kind of Solver
  unsigned mode;    // 入力パラメータの次元モード（無次元/有次元）
  double   CFL;     // Δtの決定に使うCFLなど
  double   deltaT;  // Δt（無次元）
  double   dh;      // 格子幅（無次元）
  double   Reynolds;// レイノルズ数
  double   Peclet;  // ペクレ数
  
public:
  DTcntl() {
    scheme = 0;
    CFL = 0.0;
    KOS = 0;
    mode = 0;
    deltaT = 0.0;
    dh = Reynolds = Peclet = 0.0;
  }
  ~DTcntl() {}
  
public:
  //@fn unsigned get_Scheme(void) const
  unsigned get_Scheme(void) const { return scheme; };
  
  //@fn double get_dt(void) const
  double get_DT(void) const { return deltaT; };
  
  //@fn Sdouble get_CFL(void) const
  double get_CFL(void) const { return CFL; };
  
  //@fn double DTcntl::dtCFL(const double Uref) const
  //@brief CFL数で指定されるdtを計算
  //@param Uref 速度の参照値（無次元）
  double dtCFL(const double Uref) const {
    double v = (Uref < 1.0) ? 1.0 : Uref; // 1.0 is non-dimensional reference velocity
    return (dh*CFL / v);
  }

  //@fn double DTcntl::dtDFN(const double coef)
  //@brief 拡散数のdt制限
  //@param coef 係数 (Reynolds number or Peclet number)
  double dtDFN(const double coef) const { return coef * dh*dh/6.0; }
  
  bool chkDtSelect(void);
  bool set_Scheme (const char* str, const double val);
  
  unsigned set_DT (const double vRef);
  void set_Vars   (const unsigned m_kos, const unsigned m_mode, const double m_dh, const double re, const double pe);
};



class ReferenceFrame {
  
protected:
  unsigned Frame;    /// 参照座標系
  double TimeAccel;  /// 加速時間（無次元）
  double v00[4];     /// 参照速度（無次元
  double GridVel[3]; /// 座標系の移動速度（無次元）
  
public:
  /// 参照系の定義
  enum frame_type {
    frm_static,
    frm_translation,
    frm_rotation
  };
  
  ReferenceFrame(){
    Frame     = 0;
    TimeAccel = 0.0;
    v00[0] = v00[1] = v00[2] = v00[3] = 0.0;
    GridVel[0] = GridVel[1] = GridVel[2] = 0.0;
  }
  ~ReferenceFrame() {}
  
  void setAccel  (const double m_timeAccel);
  void setFrame  (const unsigned m_frame);
  void setGridVel(const double* m_Gvel);
  void setV00    (const double time, const bool init=false);
  
  //@fn unsigned getFrame(void) const
  //@brief Frameを返す
  unsigned getFrame(void) const {
    return Frame;
  }
  
  //@fn double getAccel(void) const
  //@brief Frameを返す
  double getAccel(void) const {
    return TimeAccel;
  }
  
  //@fn void copyV00(double* m_v0) const
  //@brief v00をコピーする
  void copyV00(double* m_v0) const {
    for (int i=0; i<4; i++) m_v0[i] = v00[i];
  }
  
  //@fn void copyGridVel(double* m_gv) const
  //@brief GridVelをコピーする
  void copyGridVel(double* m_gv) const {
    for (int i=0; i<3; i++) m_gv[i] = GridVel[i];
  }
};



class ItrCtl {
private:
  unsigned NormType;     /// ノルムの種類
  unsigned SubType;      /// SKIP LOOP or MASK LOOP
  unsigned ItrMax;       /// 最大反復数
  unsigned LinearSolver; /// 線形ソルバーの種類
  REAL_TYPE eps;         /// 収束閾値
  REAL_TYPE omg;         /// 加速/緩和係数
  REAL_TYPE NormValue;   /// ノルムの値
  
public:
  /// 反復制御リスト
  enum itr_cntl_key {
    ic_prs_pr,
    ic_prs_cr,
    ic_vis_cn,
    ic_tdf_ei,
    ic_tdf_cn,
    ic_END
  };
  
  /// 反復法の収束基準種別
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
  
  unsigned LoopCount;     /// 反復回数
  
  ItrCtl() {
    NormType = 0;
    ItrMax = LoopCount = LinearSolver = SubType = 0;
    eps = omg = 0.0;
  }
  ~ItrCtl() {}
  
public:
  // @fn unsigned get_LS(void) const
  // @brief 線形ソルバの種類を返す
  unsigned get_LS(void) const { return LinearSolver; }
  
  // @fn unsigned get_ItrMax(void) const
  // @brief 最大反復回数を返す
  unsigned get_ItrMax(void) const { return ItrMax; }
  
  // @fn unsigned get_LoopType(void) const
  // @brief ループ実装の種類を返す
  unsigned get_LoopType(void) const { return SubType; }
  
  // @fn REAL_TYPE get_omg(void) const
  // @brief 緩和/加速係数を返す
  REAL_TYPE get_omg(void) const { return omg; }
  
  // @fn REAL_TYPE get_eps(void) const
  // @brief 収束閾値を返す
  REAL_TYPE get_eps(void) const { return eps; }
  
  // @fn unsigned get_normType(void) const
  // @brief ノルムのタイプを返す
  unsigned get_normType(void) const { return NormType; }
  
  // @fn REAL_TYPE get_normValue(void) const
  // @brief keyに対応するノルムの値を返す
  REAL_TYPE get_normValue(void) const { return NormValue; }
  
  // @fn void set_LS(unsigned key)
  // @brief 線形ソルバの種類を設定する
  void set_LS(unsigned key) { LinearSolver=key; }
  
  // @fn void set_ItrMax(unsigned key)
  // @brief 最大反復回数を設定する
  void set_ItrMax(unsigned key) { ItrMax=key; }
  
  // @fn void set_LoopType(unsigned key)
  // @brief ループ実装の種類を設定する
  void set_LoopType(unsigned key) { SubType=key; }
  
  // @fn void set_omg(REAL_TYPE r)
  // @brief 緩和/加速係数を保持
  void set_omg(REAL_TYPE r) { omg = r; }
  
  // @fn void set_eps(REAL_TYPE r)
  // @brief 収束閾値を保持
  void set_eps(REAL_TYPE r) { eps = r; }
  
  // @fn void set_normType(unsigned n) 
  // @brief ノルムのタイプを保持
  void set_normType(unsigned n) { NormType = n; }
  
  // @fn void set_normValue(REAL_TYPE r)
  // @brief ノルム値を保持
  void set_normValue(REAL_TYPE r) { NormValue = r; }
};



class Control {
protected:
  SklSolverConfig* CF;  // for XML parsing
  
public:
  
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
    int Rough_Initial;
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
    unsigned FileOut;
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
  
  /// リスタート指定
  enum start_type {
    initial_start,
    re_start
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
  
  /// 空間分割モード
  enum Partitioning {
    Equal=1,
    Mbx
  };
  
  unsigned  Acell,
            AlgorithmF,
            AlgorithmH,
            BasicEqs,
            CheckParam,
            CnvScheme,
            FB_version,
            Fcell,
            guide,
            GuideOut,
            imax, jmax, kmax,
            KindOfSolver,
            LastStep,
            Limiter,
            MarchingScheme,
            NoBC,
            NoCompo,
            NoDimension,
            NoID,
            NoMaterial,
            NoMediumFluid,
            NoMediumSolid,
            NoWallSurface,
            num_process,
            num_thread,
            Parallelism,
            Partition,
            RefID,
            Start,
            version,
            vxFormat,
            Wcell;
	
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
  char HistoryWallName[LABEL], PolylibConfigName[LABEL];
  char RoughInit_prs_file[LABEL];
  char RoughInit_vel_file[LABEL];
  char RoughInit_temp_file[LABEL];
  
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
    NoDimension = 0;
    NoID = 0;
    NoMaterial = 0;
    NoMediumFluid = 0;
    NoMediumSolid = 0;
    NoWallSurface = 0;
    num_process = 0;
    num_thread = 0;
    Parallelism = 0;
    Partition = 0;
    RefID = 0;
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
		memset(HistoryItrName,  0, sizeof(char)*LABEL);
    memset(HistoryMonitorName,  0, sizeof(char)*LABEL);
    memset(HistoryWallName,  0, sizeof(char)*LABEL);
    memset(PolylibConfigName,  0, sizeof(char)*LABEL);
    memset(RoughInit_prs_file,  0, sizeof(char)*LABEL);
    memset(RoughInit_vel_file,  0, sizeof(char)*LABEL);
    memset(RoughInit_temp_file,  0, sizeof(char)*LABEL);

    CF = NULL;
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
    Mode.Rough_Initial = 0;
    
    LES.Calc=0;
    LES.Model=0;
    LES.Cs = 0.0;
    LES.damping_factor=0.0;
    
    FIO.FileOut = 0;
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
  virtual ~Control() {}
  
protected:
  bool getXML_PrsAverage(void);
  
  const CfgElem* getXML_Pointer(const char* key, string section);
  
  void convertHexCoef        (REAL_TYPE* cf);
  void convertHexCoef        (REAL_TYPE* cf, REAL_TYPE DensityMode);
  void findXMLCriteria       (const CfgElem *elmL1, const char* key, unsigned order, ItrCtl* IC);
  void getXML_Algorithm      (void);
  void getXML_Average_option (void);
  void getXML_ChangeID       (void);
  void getXML_CheckParameter (void);
  void getXML_Convection     (void);
  void getXML_Derived        (void);
  void getXML_Dimension      (void);
  void getXML_FileIO         (void);
  void getXML_Iteration      (ItrCtl* IC);
  void getXML_KindOfSolver   (const CfgElem *elmL1);
  void getXML_LES_option     (void);
  void getXML_Log            (void);
  void getXML_Para_ND        (void);
  void getXML_Para_Ref       (void);
  void getXML_Para_Temp      (void);
  void getXML_Para_Wind      (void);
  void getXML_PMtest         (void);
  void getXML_ReferenceFrame (ReferenceFrame* RF);
  void getXML_restart_rough  (void);
  void getXML_Scaling        (void);
  void getXML_Solver_Properties (void);
  void getXML_Time_Control   (DTcntl* DT);
  void getXML_TimeMarching   (void);
  void getXML_Unit           (void);
  void getXML_VarRange       (void);
  void getXML_Wall_type      (void);
  void printArea             (FILE* fp, unsigned G_Fcell, unsigned G_Acell, unsigned* G_size);
  void printVoxelSize        (unsigned* gs, FILE* fp);
  void printInitValues       (FILE* fp);
  void printLS               (FILE* fp, ItrCtl* IC);
  void printParaConditions   (FILE* fp);
  void printSteerConditions  (FILE* fp, ItrCtl* IC, DTcntl* DT, ReferenceFrame* RF);
  void setDimParameters      (void);

public:
  bool receiveCfgPtr(SklSolverConfig* cfg);
  
  const char* getVoxelFileName(void);
  
  REAL_TYPE getCellSize(unsigned* m_size);
  REAL_TYPE OpenDomainRatio(unsigned dir, REAL_TYPE area, const unsigned Dims, unsigned* G_size);
	
  void displayParams            (FILE* mp, FILE* fp, ItrCtl* IC, DTcntl* DT, ReferenceFrame* RF);
  void getXML_Para_Init         (void);
  void getXML_Polygon           (void);
  void getXML_Sampling          (void);
  void getXML_Steer_1           (DTcntl* DT);
  void getXML_Steer_2           (ItrCtl* IC, ReferenceFrame* RF);
  void printAreaInfo            (FILE* fp, FILE* mp, unsigned G_Fcell, unsigned G_Acell, unsigned* G_size);
  void printDomain              (FILE* fp, unsigned* G_size, REAL_TYPE* G_org, REAL_TYPE* G_Lbx);
  void printDomainInfo          (FILE* mp, FILE* fp, unsigned* G_size, REAL_TYPE* G_org, REAL_TYPE* G_Lbx);
  void printNoCompo             (FILE* fp);
  void setDomainInfo            (unsigned* m_sz, REAL_TYPE* m_org, REAL_TYPE* m_pch, REAL_TYPE* m_wth);
  void setParameters            (MaterialList* mat, CompoList* cmp, unsigned NoBaseBC, BoundaryOuter* BO, ReferenceFrame* RF);
  void tell_Avr_Interval_2_Sphere(void);
  void tell_Interval_2_Sphere   (void);
  
  virtual void getXML_Version   (void);
  
  static string getDirection(unsigned dir);
	string getNormString(unsigned d);
  
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
	
  //@fn const SklSolverConfig* getSolverConfig(void)
  //@brief CFオブジェクトを返す
	const SklSolverConfig* getSolverConfig(void) const { return CF; };
  
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
  
};

#endif // _SKL_FB_CONTROL_H_
