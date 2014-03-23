#ifndef _FB_DEFINE_H_
#define _FB_DEFINE_H_

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
 * @file   FB_Define.h
 * @brief  FlowBase Definition Header
 * @author aics
 * @attention このマクロを変更した場合には，対応するffv_f_param.hも変更すること
 */

#include "mydebug.h"

#define FB_VERS "1.4.8"


// precision
#ifdef _REAL_IS_DOUBLE_
#define REAL_TYPE double
#else
/** 実数型の指定
 * - デフォルトでは、REAL_TYPE=float
 * - コンパイル時オプション-D_REAL_IS_DOUBLE_を付与することで
 *   REAL_TYPE=doubleになる
 */
#define REAL_TYPE float
#endif

#define SINGLE_EPSILON 1.19e-7   // 5.96e-8  * 2
#define DOUBLE_EPSILON 4.44e-16  // 2.22e-16 * 2
#define ROUND_EPS 4.0e-7

#define KELVIN    273.15
#define BOLTZMAN  1.0

// general
#define NOFACE      6
#define ON          1
#define OFF         0
#define DETAIL      2
#define FB_FILE_PATH_LENGTH 64
#define FB_BUFF_LENGTH      256

// DIMENSITON
#define DIMENSION_2D 2
#define DIMENSION_3D 3

// 定常or非定常解析
#define FB_STEADY   0
#define FB_UNSTEADY 1

// IO mode
#define IO_GATHER     0 // bool => false
#define IO_DISTRIBUTE 1

// UNIT
#define NONDIMENSIONAL  0
#define DIMENSIONAL     1

// shape approximation
#define BINARY     1
#define CUT_INFO   2

// Time marching
#define TM_O_1ST    1
#define TM_O_2ND    2
#define TM_O_3RD    3

// Linear Solver, do not use zero
#define JACOBI        1
#define SOR           2
#define SOR2SMA       3
#define GMRES         4
#define RBGS          5
#define PCG           6
#define PBiCGSTAB     7
#define VP_ITERATION  8


// KindOfSolver
#define FLOW_ONLY               0
#define THERMAL_FLOW            1
#define THERMAL_FLOW_NATURAL    2
#define CONJUGATE_HT            3
#define CONJUGATE_HT_NATURAL    4
#define SOLID_CONDUCTION        5

// 基礎方程式
#define INCMP        0
#define LTDCMP       1
#define CMPRSS       2
#define INCMP_2PHASE 3

// 浮力計算のモード
#define NO_BUOYANCY 1
#define BOUSSINESQ  2
#define LOW_MACH    3


// 外部境界条件
#define OBC_MASK      31 // 外部境界と内部境界の識別子
#define OBC_WALL      1
#define OBC_SYMMETRIC 2
#define OBC_OUTFLOW   3
#define OBC_SPEC_VEL  4
#define OBC_PERIODIC  5
#define OBC_TRC_FREE  6
#define OBC_FAR_FIELD 7
#define OBC_INTRINSIC 8


// マスクのビット幅
#define MASK_8     0xff // 8 bit幅
#define MASK_6     0x3f // 6 bit幅
#define MASK_5     0x1f // 5 bit幅
#define CMP_BIT_W  32 // 5bit


// エンコードビット　共通
#define ACTIVE_BIT 31
#define STATE_BIT  30

// エンコードビット　B
#define FORCING_BIT   28 //  外力モデルの識別子
#define TOP_VF        20 //  Volume Fractionの先頭ビット

// エンコードビット　B : Heat part
#define ADIABATIC_T  19
#define ADIABATIC_B  18
#define ADIABATIC_N  17
#define ADIABATIC_S  16
#define ADIABATIC_E  15
#define ADIABATIC_W  14
#define GMA_T        13
#define GMA_B        12
#define GMA_N        11
#define GMA_S        10
#define GMA_E        9
#define GMA_W        8
#define H_DIAG       5
// 0-4 bitはコンポーネント番号

// エンコードビット　P
#define BC_D_T     29
#define BC_D_B     28
#define BC_D_N     27
#define BC_D_S     26
#define BC_D_E     25
#define BC_D_W     24

#define BC_N_T     23
#define BC_N_B     22
#define BC_N_N     21
#define BC_N_S     20
#define BC_N_E     19
#define BC_N_W     18

#define BC_NDAG_T  17
#define BC_NDAG_B  16
#define BC_NDAG_N  15
#define BC_NDAG_S  14
#define BC_NDAG_E  13
#define BC_NDAG_W  12

#define BC_DIAG    9

#define FACING_T   8
#define FACING_B   7
#define FACING_N   6
#define FACING_S   5
#define FACING_E   4
#define FACING_W   3
#define VLD_CNVG   2
#define VBC_UWD    1


// エンコードビット　C
#define BC_FACE_T  25
#define BC_FACE_B  20
#define BC_FACE_N  15
#define BC_FACE_S  10
#define BC_FACE_E  5
#define BC_FACE_W  0


// Component Type
#define OBSTACLE     1
#define ADIABATIC    2
#define HEATFLUX     3
#define TRANSFER     4
#define ISOTHERMAL   5
#define RADIANT      6
#define SPEC_VEL     7
#define OUTFLOW      8
#define IBM_DF       9
#define HEAT_SRC     10 // Hsrc
#define CNST_TEMP    11
#define HEX          12 // Forcing
#define FAN          13
#define DARCY        14
#define PERIODIC     15
#define MONITOR      16


// 熱伝達係数のモード コンポーネントタイプと並列
#define HT_S     21
#define HT_SN    22
#define HT_SF    23

// 圧力条件のタイプ
#define P_GRAD_ZERO  1
#define P_GRAD_NS    2
#define P_DIRICHLET  3

// 速度の流出境界における対流速度の評価法
#define V_MINMAX     1
#define V_AVERAGE    2

// state
#define SOLID      0
#define FLUID      1
#define GAS        2
#define LIQUID     3
#define ANY_STATE  4

// サンプリング方法
#define SAMPLING_NEAREST       1  /// モニタ点を含むセルでの値
#define SAMPLING_INTERPOLATION 2  /// 三重線形補間
#define SAMPLING_SMOOTHING     3  /// 局所平均

// サンプリングモード
#define SAMPLING_ALL        1 ///< 全タイプのセルを対象
#define SAMPLING_FLUID_ONLY 2 ///< 流体セルのみを対象
#define SAMPLING_SOLID_ONLY 3 ///< 固体セルのみを対象


// GRIDの種類
#define SINGLE_GRID 0
#define MULTI_GRID  1

// *.xyzファイルの出力を時刻暦にするかどうか
#define GRID_NOT_MOVE 0
#define GRID_MOVE     1

// IBLANKフラグのセットの有無
#define NOT_SET_IBLANK 0
#define SET_IBLANK     1


// 判定マクロ
// BCindex aの状態が流体であればtrueを返す (uint a)
#define IS_FLUID(a) ( ((a >> STATE_BIT) & 0x1) ? true : false )

// aをbだけ右シフトしてデコードする (uint a, b)
#define BIT_SHIFT(a,b) ( (a >> b) & 0x1 )

// コンポーネントエントリを返す (uint a)
#define DECODE_CMP(a) ( a & MASK_5 )

// Volume Fraction[0-255]を返す (uint a)
#define DECODE_VF(a) ( (a >> TOP_VF) & MASK_8 )

// BCindex aの第bビットがONかどうかを調べ，ONのとき，trueを返す
#define TEST_BIT(a,b) ( ( (a >> b) & 0x1 ) ? true : false )

// BCindex aの第bビットをREALにキャストして返す
#define GET_SHIFT_F(a,b) ( (REAL_TYPE)( (a>>b) & 0x1 ) )

// BCindexにエンコードされたFaceBCのインデクスを返す
#define GET_FACE_BC(a,b) ( (a>>b) & MASK_5 )

// 6面のいずれかにBCが設定されている場合，true
#define TEST_BC(s) ( (s & 0x3fffffff) != 0 )


/** 3次元インデクス(i,j,k) -> 1次元インデクス変換マクロ
 *  @note i,j,kインデクスはF表記
 *  @param [in] _I  i方向インデクス
 *  @param [in] _J  j方向インデクス
 *  @param [in] _K  k方向インデクス
 *  @param [in] _NI i方向インデクスサイズ
 *  @param [in] _NJ j方向インデクスサイズ
 *  @param [in] _NK k方向インデクスサイズ
 *  @param [in] _VC 仮想セル数
 *  @return 1次元インデクス
 */
#define _F_IDX_S3D(_I,_J,_K,_NI,_NJ,_NK,_VC) \
( (size_t)(_K+_VC-1) * (size_t)(_NI+2*_VC) * (size_t)(_NJ+2*_VC) \
+ (size_t)(_J+_VC-1) * (size_t)(_NI+2*_VC) \
+ (size_t)(_I+_VC-1) \
)


/** 4次元インデクス(n,i,j,k) -> 1次元インデクス変換マクロ
 *  @note i,j,kインデクスはF表記、NはC表記
 *  @param [in] _N  成分インデクス
 *  @param [in] _I  i方向インデクス
 *  @param [in] _J  j方向インデクス
 *  @param [in] _K  k方向インデクス
 *  @param [in] _NN 成分数
 *  @param [in] _NI i方向インデクスサイズ
 *  @param [in] _NJ j方向インデクスサイズ
 *  @param [in] _NK k方向インデクスサイズ
 *  @param [in] _VC 仮想セル数
 *  @return 1次元インデクス
 */
#define _F_IDX_S4DEX(_N,_I,_J,_K,_NN,_NI,_NJ,_NK,_VC) \
( (size_t)(_NN) * \
_F_IDX_S3D(_I,_J,_K,_NI,_NJ,_NK,_VC) \
+ (size_t)(_N) )


/** 3次元インデクス(3,i,j,k) -> 1次元インデクス変換マクロ
 *  @note i,j,kインデクスはF表記、NはC表記
 *  @param [in] _N  成分インデクス
 *  @param [in] _I  i方向インデクス
 *  @param [in] _J  j方向インデクス
 *  @param [in] _K  k方向インデクス
 *  @param [in] _NI i方向インデクスサイズ
 *  @param [in] _NJ j方向インデクスサイズ
 *  @param [in] _NK k方向インデクスサイズ
 *  @param [in] _VC 仮想セル数
 */
#define _F_IDX_V3DEX(_N,_I,_J,_K,_NI,_NJ,_NK,_VC) (_F_IDX_S4DEX(_N,_I,_J,_K,3,_NI,_NJ,_NK,_VC))


/** 4次元インデクス(i,j,k,n) -> 1次元インデクス変換マクロ
 *  @note i,j,kインデクスはF表記、NはC表記
 *  @param [in] _I  i方向インデクス
 *  @param [in] _J  j方向インデクス
 *  @param [in] _K  k方向インデクス
 *  @param [in] _N  成分インデクス
 *  @param [in] _NI i方向インデクスサイズ
 *  @param [in] _NJ j方向インデクスサイズ
 *  @param [in] _NK k方向インデクスサイズ
 *  @param [in] _VC 仮想セル数
 *  @return 1次元インデクス
 */
#define _F_IDX_S4D(_I,_J,_K,_N,_NI,_NJ,_NK,_VC) \
( (size_t)(_N) * (size_t)(_NI+2*_VC) * (size_t)(_NJ+2*_VC) * (size_t)(_NK+2*_VC) \
+ _F_IDX_S3D(_I,_J,_K,_NI,_NJ,_NK,_VC) \
)


/** 3次元インデクス(i,j,k,3) -> 1次元インデクス変換マクロ
 *  @note i,j,kインデクスはF表記、NはC表記
 *  @param[in] _I  i方向インデクス
 *  @param[in] _J  j方向インデクス
 *  @param[in] _K  k方向インデクス
 *  @param[in] _N  成分インデクス
 *  @param[in] _NI i方向インデクスサイズ
 *  @param[in] _NJ j方向インデクスサイズ
 *  @param[in] _NK k方向インデクスサイズ
 *  @param[in] _VC 仮想セル数
 *  @return 1次元インデクス
 */
#define _F_IDX_V3D(_I,_J,_K,_N,_NI,_NJ,_NK,_VC) (_F_IDX_S4D(_I,_J,_K,_N,_NI,_NJ,_NK,_VC))


/**
 * @brief 非同期通信のリクエストIDアクセス
 * @param [in] _DIR 送受信方向
 * @param [in] _KEY オペレーション(send/recv)
 */
#define _KEY_RQ(_DIR,_KEY) (_DIR*2+_KEY)


/**
 * @brief パック配列アクセス
 * @param [in] _SZ 配列サイズ（1D）
 * @param [in] _BF バッファ種類
 * @param [in] _L  アドレス（パックする変数カウント）
 * @note バッファ種類 >> enum Async_Buffer_Direction
 */
#define _PACK(_SZ,_BF,_L) (_BF*_SZ+_L)


// AoS
typedef struct {
  REAL_TYPE p0;
  REAL_TYPE p1;
} Gemini_R;

/// FFVC_EXEC_MODE
enum ffvc_execution {
  ffvc_solver=1,
  ffvc_solverAS,
  ffvc_filter,
  ffvc_asd
};

/// 変数の種類
enum Kind_of_vars {
  var_Velocity=0,
  var_Pressure,
  var_Temperature,
  var_Fvelocity,
  var_VelocityAvr,
  var_PressureAvr,
  var_TemperatureAvr,
  var_TotalP,
  var_Helicity,
  var_Vorticity,
  var_Qcr,
  var_END
};

/// スタート指定
enum start_type {
  initial_start=0,
  restart_sameDiv_sameRes,
  restart_sameDiv_refinement,
  restart_diffDiv_sameRes,
  restart_diffDiv_refinement
};

/// 圧力単位
enum Unit_Pressure {
  Unit_Gauge=1,
  Unit_Absolute
};

/// 長さ単位
enum Length_Unit {
  LTH_ND=0,
  LTH_m,
  LTH_cm,
  LTH_mm
};

/// 同期モード
enum Synch_Mode {
  comm_sync=1,
  comm_async
};

/// send/recv Key
enum CommKeys {
  key_send=0,
  key_recv
};

/// 非同期バッファの面と種類
enum Async_Buffer_Direction {
  face_m_send=0,
  face_p_send,
  face_p_recv,
  face_m_recv
};

/// スカラ/ベクトルの種別
enum sv_type {
  kind_scalar=1,
  kind_vector
};

/// File Format
enum File_format {
  sph_fmt=0,
  bov_fmt
};

/** 反復制御リスト */
enum itr_cntl_key
{
  ic_prs1,
  ic_prs2,
  ic_vel1,
  ic_tmp1,
  ic_div,
  ic_END
};


/** 反復法の収束基準種別 */
enum norm_type
{
  v_div_max=1,
  v_div_dbg,
  dx_b,
  r_b,
  r_r0
};

/** 組み込み例題のID */
enum Intrinsic_class
{
  id_Duct,
  id_PPLT2D,
  id_PMT,
  id_Rect,
  id_Cylinder,
  id_Step,
  id_Polygon,
  id_Sphere,
  id_Jet
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


/// モニタ点指定タイプ型
enum Monitor_Type {
  mon_POINT_SET,
  mon_LINE,
  mon_CYLINDER,
  mon_BOX,
  mon_POLYGON
};

enum DIRection {
  X_minus=0,
  X_plus,
  Y_minus,
  Y_plus,
  Z_minus,
  Z_plus
};

/*
 * @brief 指定面方向のカットIDをとりだす
 * @param [in] bid  int variable
 * @param [in] dir  方向コード (w/X_MINUS=0, e/X_PLUS=1, s/2, n/3, b/4, t/5)
 */
inline int getBit5 (const int bid, const int dir)
{
  return ( (bid >> dir*5) & MASK_5 );
}


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


/*
 * @brief 5bit幅のIDの設定
 * @param [in,out] b   int 変数
 * @param [in]     q   5-bit幅のID (1-31)
 */
inline void setBitID (int& b, const int q)
{
  b &= ( ~(0x1f) ); // 下位5bitをゼロにする
  b |= q;           // 書き込む
}


#endif // _FB_DEFINE_H_
