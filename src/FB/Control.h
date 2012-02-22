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

#include "Skl.h"
#include "SklSolverBase.h"
#include "FBDefine.h"
#include "mydebug.h"
#include "config/SklSolverConfig.h"
#include "SklUtil.h"
#include "Material.h"
#include "Component.h"
#include "IterationCtl.h"
#include "Parallel_node.h"
#include "Ref_Frame.h"
#include "FBUtility.h"
#include "Monitor.h"
#include "BndOuter.h"
#include "DTcntl.h"
#include "Interval_Mngr.h"

using namespace SklCfg;  // to use SklSolverConfig* cfg

class Control : public Parallel_Node {
protected:
  SklSolverConfig* CF;  // for XML parsing
  
public:
  
  // 各種モード　パラメータ
  typedef struct {
    unsigned Average;
    unsigned Base_Medium;
    unsigned Buoyancy;
    unsigned Example;
    unsigned Helicity;
    unsigned I2VGT;
    unsigned Log_Base;
    unsigned Log_Itr;
    unsigned Log_Wall;
    unsigned Medium_Spec;
    unsigned PDE;
    unsigned Precision;
    unsigned Profiling;
    unsigned PrsNeuamnnType;
    unsigned ShapeAprx;
    unsigned Steady;
    unsigned TP;
    unsigned VarArrange;
    unsigned VRT;
    unsigned Wall_profile;
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
            InnerItr,
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
            RefID,
            Start,
            version,
            vxFormat,
            Wcell;
	
  REAL_TYPE  BasePrs,
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
    InnerItr = 0;
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

    CF = NULL;
    TotalMemory = 0;
    
    iv.Density     = 0.0;
    iv.Energy      = 0.0;
    iv.Pressure    = 0.0;
    iv.Temperature = 0.0;
    
    Mode.Average = 0;
    Mode.Buoyancy = 0;
    Mode.Example = 0;
    Mode.Helicity = 0;
    Mode.I2VGT = 0;
    Mode.Log_Base = 0;
    Mode.Log_Itr = 0;
    Mode.Log_Wall = 0;
    Mode.PDE = 0;
    Mode.Precision = 0;
    Mode.Profiling = 0;
    Mode.PrsNeuamnnType = 0;
    Mode.ShapeAprx = 0;
    Mode.Steady = 0;
    Mode.TP = 0;
    Mode.VarArrange = 0;
    Mode.VRT = 0;
    Mode.Wall_profile = 0;
    
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
  void getXML_Mon_Line       (MonitorList* M, const CfgElem *elmL2, REAL_TYPE from[3], REAL_TYPE to[3], int& nDivision);
  void getXML_Mon_Pointset   (MonitorList* M, const CfgElem *elmL2, vector<MonitorCompo::MonitorPoint>& pointSet);
  void getXML_Para_ND        (void);
  void getXML_Para_Ref       (void);
  void getXML_Para_Temp      (void);
  void getXML_Para_Wind      (void);
  void getXML_PMtest         (void);
  void getXML_ReferenceFrame (ReferenceFrame* RF);
  void getXML_Scaling        (void);
  void getXML_Solver_Properties (void);
  void getXML_Time_Control   (DTcntl* DT);
  void getXML_TimeMarching   (void);
  void getXML_Unit           (void);
  void getXML_VarArrangement (void);
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
  bool chkMediumConsistency(void);
  bool receiveCfgPtr(SklSolverConfig* cfg);
  
  const char* getVoxelFileName(void);
  
  REAL_TYPE getCellSize(unsigned* m_size);
  REAL_TYPE OpenDomainRatio(unsigned dir, REAL_TYPE area, const unsigned Dims, unsigned* G_size);
	
  void displayParams            (FILE* mp, FILE* fp, ItrCtl* IC, DTcntl* DT, ReferenceFrame* RF);
  void getXML_InnerItr          (void);
  void getXML_History           (void);
  void getXML_Monitor           (MonitorList* M);
  void getXML_Para_Init         (void);
  void getXML_Steer_1           (DTcntl* DT);
  void getXML_Steer_2           (ItrCtl* IC, ReferenceFrame* RF);
  void printAreaInfo            (FILE* fp, FILE* mp, unsigned G_Fcell, unsigned G_Acell, unsigned* G_size);
  void printDomain              (FILE* fp, unsigned* G_size, REAL_TYPE* G_org, REAL_TYPE* G_Lbx);
  void printDomain_debug        (void);
  void printDomainInfo          (FILE* mp, FILE* fp, unsigned* G_size, REAL_TYPE* G_org, REAL_TYPE* G_Lbx);
  void printNoCompo             (FILE* fp);
  void select_Itr_Impl          (ItrCtl* IC);
  void setDomainInfo            (unsigned* m_sz, REAL_TYPE* m_org, REAL_TYPE* m_pch, REAL_TYPE* m_wth);
  void setParameters            (MaterialList* mat, CompoList* cmp, unsigned NoBaseBC, BoundaryOuter* BO, ReferenceFrame* RF);
  void tell_Avr_Interval_2_Sphere(void);
  void tell_Interval_2_Sphere   (void);
  
  virtual void getXML_Version   (void);
  
  static string getDirection(unsigned dir);
	string getNormString(unsigned d);
  
  unsigned countCompo  (CompoList* cmp, unsigned label);
  unsigned getNoInFiles(const char* key);
  
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
  
  //@fn unsigned isOutflow(void) const
  //@brief 流出コンポーネントがあれば1を返す
  unsigned isOutflow(void) const {
    return EnsCompo.outflow;
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
