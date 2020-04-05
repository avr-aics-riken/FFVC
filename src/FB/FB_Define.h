#ifndef _FB_DEFINE_H_
#define _FB_DEFINE_H_

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
 * @file   FB_Define.h
 * @brief  FlowBase Definition Header
 * @author aics
 * @attention このマクロを変更した場合には，対応するffv_f_param.hも変更すること
 */

#include <float.h>
#include <math.h>
#include <stdio.h>

#define FB_VERS "1.5.4"

#define SINGLE_EPSILON 1.19e-7
#define DOUBLE_EPSILON 2.22e-16

#define DELTA_QT (1.0/511.0)

/** 実数型の指定
 * - デフォルトでは、REAL_TYPE=float
 * - コンパイル時オプション-D_REAL_IS_DOUBLE_を付与することで
 *   REAL_TYPE=doubleになる
 */
#ifdef _REAL_IS_DOUBLE_
#define REAL_TYPE double
#define REAL_TYPE_EPSILON DBL_MIN
#define ROUND_EPS SINGLE_EPSILON*5.0
#else
#define REAL_TYPE float
#define REAL_TYPE_EPSILON FLT_MIN
#define ROUND_EPS SINGLE_EPSILON*5.0
#endif


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

// Linear Solver, do not use zero. Zero indicates undefined.
#define JACOBI        1
#define SOR           2
#define SOR2SMA       3
#define GMRES         4
#define PCG           5
#define BiCGSTAB      6

#define FREQ_OF_RESTART 15 // リスタート周期

// KindOfSolver
#define COLD_FLOW               0
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


// マスクのビット幅
#define MASK_30    0x3fffffff // 30 bit幅
#define MASK_10    0x3ff // 10 bit幅
#define MASK_9     0x1ff // 9 bit幅
#define MASK_8     0xff  // 8 bit幅
#define MASK_6     0x3f  // 6 bit幅
#define MASK_5     0x1f  // 5 bit幅
#define CMP_BIT_W  32    // 5bit
#define QT_9       511   // 9bit幅の最大値

/*
 Bit   bcd[]            bid[]          bcp[]         cdf[]
 31    active           active         active        active
 30    state            stete          state         state
 29  + BC_D_T           +              BC_D_T        +
 28  | BC_D_B           |              BC_D_B        |
 27  | BC_D_N   VBC     |              BC_D_N        |
 26  | BC_D_S           |              BC_D_S        |
 25  | BC_D_E           + Z_plus       BC_D_E        + BC_FACE_T
 24  + BC_D_W           +              BC_D_W        +
 23                     |              BC_N_T        |
 22                     |              BC_N_B        |
 21                     |              BC_N_N        |
 20                     + Z_minus      BC_N_S        + BC_FACE_B
 19                     +              BC_N_E        +
 18                     |              BC_N_W        |
 17                     |              BC_DN_T       |
 16                     |              BC_DN_B       |
 15                     + Y_plus       BC_DN_N       + BC_FACE_N
 14                     +              BC_DN_S       +
 13                     |              BC_DN_E       |
 12                     |              BC_DN_W       |
 11                     |              BC_NDAG_T     |
 10                     + Y_minus      BC_NDAG_B     + BC_FACE_S
  9                     +              BC_NDAG_N     +
  8                     |              BC_NDAG_S     |
  7                     |              BC_NDAG_E     |
  6                     |              BC_NDAG_W     |
  5                     + X_plus       BC_DIAG       + BC_FACE_E
  4  +                  +                            +
  3  |                  |                            |
  2  | cmp[]/mat[]      |                            |
  1  |                  |                            |
  0  +                  + X_minus                    + BC_FACE_W


 
Bit    cut_l[]             cut_u[]             nrm[]
 31
 30     STATE_BIT            STATE_BIT          Initialized (STATE_BIT)
 29   + Ens Y_minus        + Ens Z_plus         + Z_nrm_Sign
 28   |                    |                    |
 27   |                    |                    |
 26   |                    |                    |
 25   |                    |                    |
 24   |                    |                    |
 23   |                    |                    |
 22   |                    |                    |
 21   |                    |                    |
 20   + Y_minus(2)         + Z_plus(5)          + Z_nrm
 19   + Ens X_plus         + Ens Z_minus        + Y_nrm_Sign
 18   |                    |                    |
 17   |                    |                    |
 16   |                    |                    |
 15   |                    |                    |
 14   |                    |                    |
 13   |                    |                    |
 12   |                    |                    |
 11   |                    |                    |
 10   + X_plus(1)          + Z_minus(4)         + Y_nrm
  9   + Ens X_minus        + Ens Y_plus         + X_nrm_Sign
  8   |                    |                    |
  7   |                    |                    |
  6   |                    |                    |
  5   |                    |                    |
  4   |                    |                    |
  3   |                    |                    |
  2   |                    |                    |
  1   |                    |                    |
  0   + X_minus(0)         + Y_plus(3)          + X_nrm

 */


// エンコードビット 共通
#define ACTIVE_BIT 31
#define STATE_BIT  30

// エンコードビット B
#define FORCING_BIT   28 //  外力モデルの識別子
#define TOP_VF        20 //  Volume Fractionの先頭ビット

// エンコードビット B : Heat part
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

// エンコードビット bcp[] for Poisson Eq
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

#define BC_DN_T    17
#define BC_DN_B    16
#define BC_DN_N    15
#define BC_DN_S    14
#define BC_DN_E    13
#define BC_DN_W    12

#define BC_NDAG_T  11
#define BC_NDAG_B  10
#define BC_NDAG_N  9
#define BC_NDAG_S  8
#define BC_NDAG_E  7
#define BC_NDAG_W  6
#define BC_DIAG    5


// エンコードビット cdf[] for veclocity BC
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
#define SOLIDREV     17
#define OUTER_BC     18


// 熱伝達係数のモード コンポーネントタイプと並列
#define HT_S     21
#define HT_SN    22
#define HT_SF    23


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


// ISNAN macro >> isnan() is defined as macro in math.h. If isnan() is undefined, use std::isnan()
#ifdef isnan
#define ISNAN(_X) isnan(_X)
#else
#define ISNAN(_X) std::isnan(_X)
#endif


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
#define IS_CUT(s) ( (s & 0x3fffffff) != 0 )


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
 * @brief 2D配列アクセス
 * @param [in] _I  最初の列
 * @param [in] _J  次の列
 * @param [in] _SZ 配列サイズ（1D方向）
 */
#define _IDX2D(_I,_J,_SZ) (_J*_SZ+_I)


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


// Divergence judgement
typedef struct {
  int MaxIteration;
  int Iteration;
  int divType;
  double divEPS;
  double divergence;
} DivConvergence;


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
  var_Div,
  var_RmsV,
  var_RmsMeanV,
  var_RmsP,
  var_RmsMeanP,
  var_RmsT,
  var_RmsMeanT,
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
  LTH_m
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
  bov_fmt,
  plt3d_fun_fmt
};

/** 反復制御リスト */
enum itr_cntl_key
{
  ic_prs1=0,
  ic_prs2,
  ic_vel1,
  ic_tmp1,
  ic_END
};


/** 反復法の収束基準種別 */
enum norm_type
{
  nrm_dx=0,
  nrm_dx_x,
  nrm_r_b,
  nrm_r_x,
  nrm_r_r0,
  nrm_div_l2,
  nrm_div_max
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
   mon_POINT_SET
  ,mon_LINE
  ,mon_CYLINDER
  ,mon_BOX
  ,mon_PLANE
  ,mon_POLYGON
  ,mon_DISC
  ,mon_UNKNOWN
};

enum DIRection {
  X_minus=0,
  X_plus,
  Y_minus,
  Y_plus,
  Z_minus,
  Z_plus
};


//@brief idxの第shiftビットをOFFにする
inline int offBit (int idx, const int shift)
{
  return ( idx & (~(0x1<<shift)) );
}


//@brief idxの第shiftビットをONにする
inline int onBit (int idx, const int shift)
{
  return ( idx | (0x1<<shift) );
}


/*
 * @brief 9bit幅の量子化
 * @param [in]  a  入力数値
 */
inline int quantize9(float a)
{
  int s = (int)((a + 0.5*DELTA_QT) / DELTA_QT);

  return s;
}

inline float getQuantized9(const int d)
{
  return (float)d*DELTA_QT;
}


/*
 * @brief 指定面方向のカットIDをとりだす
 * @param [in] bid  int variable
 * @param [in] dir  方向コード (w/X_MINUS=0, e/X_PLUS=1, s/2, n/3, b/4, t/5)
 */
inline int getBit5 (const int bid, const int dir)
{
  return ( bid >> (dir*5) & MASK_5 );
}


/*
 * @brief 5bit幅の値の設定
 * @param [in,out] b   int 変数
 * @param [in]     q   5-bit幅のID (1-31)
 * @param [in]     dir 方向コード (w/X_MINUS=0, e/X_PLUS=1, s/2, n/3, b/4, t/5)
 */
inline void setBit5 (int& b, const int q, const int dir)
{
  b &= (~(MASK_5 << (dir*5)) ); // 対象5bitをゼロにする
  b |= (q << (dir*5));          // 書き込む
}



/*
 * @brief 5bit幅の媒質IDの設定
 * @param [in,out] b   int 変数
 * @param [in]     q   5-bit幅のID (1-31)
 */
inline void setMediumID (int& b, const int q)
{
  b &= ( ~(0x1f) ); // 下位5bitをゼロにする
  b |= q;           // 書き込む
}


/*
 * @brief cut indexから指定方向交点の有無を返す
 * @param [in] cl   cut lower index
 * @param [in] dir  方向コード (w/X_MINUS=0, e/X_PLUS=1, s/2, n/3, b/4, t/5)
 * @retval 交点あり(1)、なし(0)
 */
inline int ensCutL(const int cl, const int dir)
{
  return cl >> (dir*10+9) & 0x1;
}

/*
 * @brief cut indexから指定方向交点の有無を返す
 * @param [in] cu   cut upper index
 * @param [in] dir  方向コード (w/X_MINUS=0, e/X_PLUS=1, s/2, n/3, b/4, t/5)
 * @retval 交点あり(1)、なし(0)
 */
inline int ensCutU(const int cu, const int dir)
{
  return cu >> ((dir-3)*10+9) & 0x1;
}


/*
 * @brief cut indexから指定方向の量子化値をとりだす
 * @param [in] cl   cut lower index
 * @param [in] dir  方向コード (w/X_MINUS=0, e/X_PLUS=1, s/2, n/3, b/4, t/5)
 */
inline int getBitL9(const int cl, const int dir)
{
  return cl >> (dir*10) & MASK_9;
}

/*
 * @brief cut indexから指定方向の量子化値をとりだす
 * @param [in] cu   cut upper index
 * @param [in] dir  方向コード (w/X_MINUS=0, e/X_PLUS=1, s/2, n/3, b/4, t/5)
 */
inline int getBitU9(const int cu, const int dir)
{
  return cu >> ((dir-3)*10) & MASK_9;
}



/*
 * @brief 指定方向交点の有無と距離0のチェック
 * @param [in] c    cut lower index
 * @param [in] dir  方向コード (w/X_MINUS=0, e/X_PLUS=1, s/2, n/3, b/4, t/5)
 * @retval 交点があり、かつ、距離がゼロのとき1
 */
inline int chkZeroCut(const int c, const int dir)
{
  // 量子化した距離 9ビット幅
  int d = (dir<3) ? getBitL9(c, dir) : getBitU9(c, dir);

  // 距離が記録されている
  if ( d > 0 ) return 0;

  // 距離ゼロ
  int e = (dir<3) ? ensCutL(c, dir) : ensCutU(c, dir);
  
  return (e == 1) ? 1 : 2; // 交点の記録あり(1) , なし(2)
}



/*
 * @brief cut indexから指定方向の距離を取り出す
 * @param [in] cl   cut lower index
 * @param [in] dir  方向コード (w/X_MINUS=0, e/X_PLUS=1, s/2, n/3, b/4, t/5)
 */
inline float getCutL9(const int cl, const int dir)
{
  return getQuantized9( cl >> (dir*10) & MASK_9 );
}

/*
 * @brief cut indexから指定方向の距離を取り出す
 * @param [in] cu   cut upper index
 * @param [in] dir  方向コード (w/X_MINUS=0, e/X_PLUS=1, s/2, n/3, b/4, t/5)
 */
inline float getCutU9(const int cu, const int dir)
{
  return getQuantized9( cu >> ((dir-3)*10) & MASK_9 );
}



/*
 * @brief cut indexの値の設定（9bit幅の値と交点フラグ）
 * @param [in,out] cl   cut lower index
 * @param [in]     q   9-bit幅の値
 * @param [in]     dir 方向コード (w/X_MINUS=0, e/X_PLUS=1, s/2, n/3, b/4, t/5)
 */
inline void setCutL9(int& cl, const int q, const int dir)
{
  int d = dir*10;
  cl &= ( ~(MASK_10 << d) );  // 対象10bitをゼロにする
  cl |= ( q<<d );             // 値を書き込む
  cl |= ( 0x1<<(d+9) );       // 交点フラグをON
}

/*
 * @brief cut indexの値の設定（9bit幅の値と交点フラグ）
 * @param [in,out] cl   cut lower index
 * @param [in,out] cu   cut upper index
 * @param [in]     q   9-bit幅の値
 * @param [in]     dir 方向コード (w/X_MINUS=0, e/X_PLUS=1, s/2, n/3, b/4, t/5)
 */
inline void setCutU9(int& cu, const int q, const int dir)
{
  int d = (dir-3)*10;
  cu &= ( ~(MASK_10 << d) );  // 対象10bitをゼロにする
  cu |= ( q<<d );             // 値を書き込む
  cu |= ( 0x1<<(d+9) );       // 交点フラグをON
}



/*
 * @brief cut indexの値を連結状態に変更（9bit幅の値と交点フラグ）
 * @param [in,out] cl   cut lower index
 * @param [in]     dir 方向コード (w/X_MINUS=0, e/X_PLUS=1, s/2, n/3, b/4, t/5)
 */
inline void setUncutL9(int& cl, const int dir)
{
  int d = dir*10;
  cl &= ( ~(MASK_10<<d) ); // 対象10bitをゼロにする
  cl |= ( QT_9<<d );       // 最大値を書き込む
  cl &= ( ~(0x1<<(d+9)) ); // 交点フラグをクリア
}

/*
 * @brief cut indexの値を連結状態に変更（9bit幅の値と交点フラグ）
 * @param [in,out] cu   cut upper index
 * @param [in]     dir 方向コード (w/X_MINUS=0, e/X_PLUS=1, s/2, n/3, b/4, t/5)
 */
inline void setUncutU9(int& cu, const int dir)
{
  int d = (dir-3)*10;
  cu &= ( ~(MASK_10<<d) ); // 対象10bitをゼロにする
  cu |= ( QT_9<<d );       // 最大値を書き込む
  cu &= ( ~(0x1<<(d+9)) ); // 交点フラグをクリア
}



/*
 * @brief cut indexを511で初期化
 * @param [in,out] cl   cut lower index
 * @param [in]     dir 方向コード (w/X_MINUS=0, e/X_PLUS=1, s/2, n/3, b/4, t/5)
 */
inline void initBitL9(int& cl, const int dir)
{
  int d = dir*10;
  cl &= ( ~( MASK_9<<d) );  // 対象9bitをゼロにする
  cl |= ( QT_9<<d );        // 値を書き込む
}


/*
 * @brief cut indexを511で初期化
 * @param [in,out] cu   cut upper index
 * @param [in]     dir 方向コード (w/X_MINUS=0, e/X_PLUS=1, s/2, n/3, b/4, t/5)
 */
inline void initBitU9(int& cu, const int dir)
{
  int d = (dir-3)*10;
  cu &= ( ~( MASK_9<<d) );  // 対象9bitをゼロにする
  cu |= ( QT_9<<d );        // 値を書き込む
}


#endif // _FB_DEFINE_H_
