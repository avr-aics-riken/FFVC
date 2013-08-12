#ifndef _FB_CONTROL_H_
#define _FB_CONTROL_H_

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

/**
 * @file   Control.h
 * @brief  FlowBase Control class Header
 * @author kero
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
#include "CompoFraction.h"
#include "IterationControl.h"

/* 20130611 commentout
#include "PLOT3D_read.h"
#include "PLOT3D_write.h"
 */

using namespace std;



// #################################################################
class PolygonProperty {
  
private:
  int l_ntria;       ///< ローカルなポリゴン数
  int g_ntria;       ///< グローバルなポリゴン数
  REAL_TYPE l_area;  ///< ローカルな面積
  REAL_TYPE g_area;  ///< グローバルな面積
  string group;      ///< ポリゴングループ名
  string material;   ///< Mediumtable[@]のalias
  string bc;         ///< BCのラベル
  
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
  TPControl* tpCntl;   ///< テキストパーサへのポインタ
  
public:
    
  
  /** 各種モード　パラメータ */
  typedef struct 
  {
    int Average;
    int AverageRestart;
    int Buoyancy;
    int Example;
    int Helicity;
    int FaceV;
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
    int IOmode;
    int Div_Debug;
    int IO_Voxel;
    int Format;
    int Slice;
    string OutDirPath;
    string InDirPath;
  } File_IO_Cntl;
  
  /** PLOT3D オプション */
  typedef struct 
  {
    string basename_g;
    string basename_f;
    int IS_xyz;
    int IS_q;
    int IS_funciton;
    int IS_function_name;
    int IS_fvbnd;
    int IS_DivideFunc;
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
    int Output; /// 出力単位 (Dimensional/NonDimensional)
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
  
  /** コンポーネント/BCのグローバルな存在 */
  typedef struct 
  {
    int forcing;
    int hsrc;
    int outflow;
    int periodic;
    int fraction;
    int monitor;
    int tfree;
    int vspec;
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
  
  
  int AlgorithmF;
  int AlgorithmH;
  int BasicEqs;
  int CheckParam;
  int CnvScheme;
  int FB_version;     ///< FlowBaseクラスのバージョン番号
  int Fill_Hint;      ///< フィルのヒント {no | x_minux | x_plus |...}
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
  int RefMat;         ///< 参照媒質インデクス
  int Start;
  int version;        ///< FFVバージョン番号
  
  unsigned Restart_staging;    ///< リスタート時にリスタートファイルがSTAGINGされているか
  unsigned Restart_step;       ///< リスタートステップ
  unsigned Restart_stepAvr;    ///< 平均値のリスタートステップ
  unsigned long NoWallSurface; ///< 固体表面セル数
  
  float Scaling_Factor;  ///< ポリゴンのスケーリングファクター
	
  double Tscale;
  
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
  
  IterationCtl* Criteria; ///< 反復解法の収束判定パラメータ
  
  string file_fmt_ext;
  
  string HistoryName;
  string HistoryCompoName;
  string HistoryDomfxName;
  string HistoryItrName;
  string HistoryMonitorName;
  string HistoryWallName;
  string PolylibConfigName;
  string HistoryForceName;
  
  // 出力dfiファイルのプレフィックス
  string f_dfi_out_prs;
  string f_dfi_out_vel;
  string f_dfi_out_fvel;
  string f_dfi_out_temp;
  string f_dfi_out_prsa;
  string f_dfi_out_vela;
  string f_dfi_out_tempa;
  string f_dfi_out_hlt;
  string f_dfi_out_tp;
  string f_dfi_out_i2vgt;
  string f_dfi_out_vrt;
  string f_dfi_out_div;
  
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
  
  string RefMedium;    ///< 参照媒質名 -> int RefMat
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
    FB_version = 0;
    CnvScheme = 0;
    Fill_Hint = -1;
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
    RefMat = 0;
    Restart_staging = 0;    
    Restart_step = 0;
    Restart_stepAvr = 0;
    Start = 0;
    version = 0;
    
    Scaling_Factor = 1.0;
    
    PlotIntvl = 0.0;
    Domain_p1 = Domain_p2 = 0.0;
    RefVelocity = RefLength = RefDensity = RefSoundSpeed = RefSpecificHeat = RefKviscosity = RefLambda = RefViscosity = 0.0;
    DiffTemp = BaseTemp = BasePrs = 0.0;
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
    
    iv.Density     = 0.0;
    iv.Energy      = 0.0;
    iv.Pressure    = 0.0;
    iv.Temperature = 0.0;
    
    Mode.Average = 0;
    Mode.AverageRestart = 0;
    Mode.Base_Medium = 0;
    Mode.Buoyancy = 0;
    Mode.Example = 0;
    Mode.Helicity = 0;
    Mode.I2VGT = 0;
    Mode.FaceV = 0;
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
    Mode.VRT = 0;
    Mode.Wall_profile = 0;
    
    LES.Calc=0;
    LES.Model=0;
    LES.Cs = 0.0;
    LES.damping_factor=0.0;
    
    FIO.IOmode   = 0;
    FIO.Div_Debug  = 0;
    FIO.IO_Voxel   = 0;
    FIO.Format     = 0; // 0:sph, 1:PLOT3D, 2:BOV
    FIO.Slice      = 0;
    
    P3Op.IS_xyz = ON;
    P3Op.IS_q = OFF;
    P3Op.IS_funciton = ON;
    P3Op.IS_function_name = ON;
    P3Op.IS_fvbnd = OFF;
    
    Hide.Change_ID = 0;
    Hide.Range_Limit = 0;
    Hide.PM_Test = 0;
    
    Unit.Param  = 0;
    Unit.Output = 0;
    Unit.Prs    = 0;
    Unit.Log    = 0;
    Unit.File   = 0;
    Unit.Temp   = 0;
    
    Sampling.log = 0;
    Sampling.out_mode = 0;
    Sampling.unit = 0;
    
    EnsCompo.hsrc    = 0;
    EnsCompo.forcing = 0;
    EnsCompo.outflow = 0;
    EnsCompo.periodic= 0;
    EnsCompo.fraction= 0;
    EnsCompo.monitor = 0;
    EnsCompo.tfree   = 0;
    EnsCompo.vspec   = 0;
    
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
  
  
  // 平均値操作に関するパラメータを取得する
  void getAverageOption();

  
  // 対流項スキームのパラメータを取得する
  void getConvection();
  
  
  // 派生して計算する変数のオプションを取得する
  void getDerived();
  
  
  // 無次元パラメータを各種モードに応じて設定する
  void getDimensionlessParameter();
  
  
  // ファイル入出力に関するパラメータを取得し，sphフォーマットの出力の並列モードを指定する
  void getFieldData();
  
  
  // PLOt3Dフォーマットのオプションを指定する
  void getFormat_plot3d();
  
  
  // sphフォーマットのオプションを指定
  void getFormat_sph();
  

  // ログ出力のパラメータを取得
  void getLog();
  
  
  // Gmres反復固有のパラメータを指定する
  void getParaGmres(const string base, const int m);
  
  
  // Jacobi反復固有のパラメータを指定する
  void getParaJacobi(const string base, const int m);
  
  
  // PBiCGSTAB反復固有のパラメータを指定する
  void getParaPBiCGSTAB(const string base, const int m);
  
  
  // PCG反復固有のパラメータを指定する
  void getParaPCG(const string base, const int m);
  
  
  // RBGS反復固有のパラメータを指定する
  void getParaRBGS(const string base, const int m);
  
  
  // SOR反復固有のパラメータを指定する
  void getParaSOR(const string base, const int m);
  
  
  // RB-SOR反復固有のパラメータを指定する
  void getParaSOR2(const string base, const int m);
  
  
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
  void printInitValues(FILE* fp);
  
  
  // 線形ソルバー種別の表示
  void printLS(FILE* fp, const IterationCtl* IC);
  
  
  // 計算パラメータの表示
  void printParaConditions(FILE* fp, const MediumList* mat);
  
  
  /** 20130611
   * @brief 制御パラメータSTEERの表示
   * @param [in] IC IterationCtl
   * @param [in] DT DTcntl
   * @param [in] RF ReferenceFrame
   //* @param [in] FP3DW FileIO PLOT3D WRITE CLASS POINTER
   */
  void printSteerConditions(FILE* fp, IterationCtl* IC, const DTcntl* DT, const ReferenceFrame* RF);
  //void printSteerConditions(FILE* fp, const IterationCtl* IC, const DTcntl* DT, const ReferenceFrame* RF, FileIO_PLOT3D_WRITE* FP3DW);
  
  
  /**
   * @brief PLOT3Dファイル入出力に関するパラメータ
   * @param [in]  FP3DR    PLOT3D READクラス ポインタ
   * @param [in]  FP3DW    PLOT3D WRITEクラス ポインタ
   *
  void get_PLOT3D(FileIO_PLOT3D_READ* FP3DR, FileIO_PLOT3D_WRITE* FP3DW);
   */
  
  
  
  
public:
  
  /**
   * @brief 反復の収束判定パラメータをコピーする
   * @param [in,out] IC   反復制御用クラスの基底へのポインタ
   * @param [in]     name ラベル
   */
  void copyCriteria(IterationCtl& IC, const string name);
  
  
  /** 20130611
   * @brief 制御，計算パラメータ群の表示
   * @param [in] mp  ファイルポインタ（標準出力）
   * @param [in] fp  ファイルポインタ（ファイル出力）
   * @param [in] IC  IterationCtl
   * @param [in] DT  DTcntl
   * @param [in] RF  ReferenceFrame
   * @param [in] mat MediumList
   //* @param [in] FP3DW FileIO PLOT3D WRITE CLASS POINTER
   */
  //void displayParams(FILE* mp, FILE* fp, IterationCtl* IC, DTcntl* DT, ReferenceFrame* RF, MediumList* mat, FileIO_PLOT3D_WRITE* FP3DW);
  void displayParams(FILE* mp, FILE* fp, IterationCtl* IC, DTcntl* DT, ReferenceFrame* RF, MediumList* mat);
  
  
  //@brief Forcingコンポーネントが存在すれば1を返す
  int existForcing() const
  {
    return EnsCompo.forcing;
  }
  
  
  //@brief Hsrcコンポーネントが存在すれば1を返す
  int existHsrc() const
  {
    return EnsCompo.hsrc;
  }
  
  
  //@brief モニタコンポーネントが存在すれば1を返す
  int existMonitor() const
  {
    return EnsCompo.monitor;
  }
  
  
  //@brief 流出コンポーネントがあれば1を返す
  int existOutflow() const
  {
    return EnsCompo.outflow;
  }
  
  
  //@brief 部分周期境界コンポーネントが存在すれば1を返す
  int existPeriodic() const
  {
    return EnsCompo.periodic;
  }
  
  
  //@brief トラクションフリー境界が存在すれば1を返す
  int existTfree() const
  {
    return EnsCompo.tfree;
  }
  
  
  //@brief 体積率コンポーネントが存在すれば1を返す
  int existVfraction() const
  {
    return EnsCompo.fraction;
  }
  
  
  //@brief 体積率コンポーネントが存在すれば1を返す
  int existVspec() const
  {
    return EnsCompo.vspec;
  }
  
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
  void getSolvingMethod();
  
  
  // 制御，計算パラメータ群の取得 20130611
  //void get_Steer_1(DTcntl* DT, FileIO_PLOT3D_READ* FP3DR, FileIO_PLOT3D_WRITE* FP3DW);
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
  
  
  // @brief モニタリングのON/OFFとセルモニタの有無のみを取得
  void getMonitorList();
  
  
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
   * @brief TPのポインタを受け取る
   * @param [in] tp TPControl
   */
  void importTP(TPControl* tp);
  
  
  /**
   * @brief ソルバーがCDSタイプかどうかを返す
   * @retval CDSであればtrue
   */
  bool isCDS() const
  {
    return ( CUT_INFO == Mode.ShapeAprx ) ? true : false;
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
   * @param [in,out] cmp CompoList
   * @param [in,out] OBC BoundaryOuter
   */
  void setExistComponent(CompoList* cmp, BoundaryOuter* OBC);
  
  
  // @brief 無次元パラメータを各種モードに応じて設定する
  void setParameters(MediumList* mat, CompoList* cmp, ReferenceFrame* RF, BoundaryOuter* BO);

};

#endif // _FB_CONTROL_H_
