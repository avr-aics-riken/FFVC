/*
 * CPM - Cartesian Partition Manager
 *
 * Copyright (C) 2012 University of Tokyo.
 *
 */

/**
 * @file   cpm_Define.h
 * CPMの定義マクロ記述ヘッダーファイル
 * @author University of Tokyo
 * @date   2012/05/31
 */
#ifndef _CPM_DEFINE_H_
#define _CPM_DEFINE_H_

#include "mpi.h"

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

#if defined(_BUFSIZE_LONG_DOUBLE_)
  #define REAL_BUF_TYPE long double
#elif defined(_BUFSIZE_DOUBLE_)
  #define REAL_BUF_TYPE double
#else
  /** 袖通信バッファの型指定
   *  - デフォルトでは、REAL_BUF_TYPE=REAL_TYPE
   *  - コンパイル時オプション-D_BUFSIZE_DOUBLE_を付与することで
   *    REAL_BUF_TYPE=doubleになる
   *  - コンパイル時オプション-D_BUFSIZE_LONG_DOUBLE_を付与することで
   *    REAL_BUF_TYPE=long doubleになる
   */
  #define REAL_BUF_TYPE REAL_TYPE
#endif


/** 3次元インデクス(i,j,k) -> 1次元インデクス変換マクロ
 *  @param[in] _I  i方向インデクス
 *  @param[in] _J  j方向インデクス
 *  @param[in] _K  k方向インデクス
 *  @param[in] _NI i方向インデクスサイズ
 *  @param[in] _NJ j方向インデクスサイズ
 *  @param[in] _NK k方向インデクスサイズ
 *  @param[in] _VC 仮想セル数
 *  @return 1次元インデクス
 */
#define _IDX_S3D(_I,_J,_K,_NI,_NJ,_NK,_VC) \
( size_t(_K+_VC) * size_t(_NI+2*_VC) * size_t(_NJ+2*_VC) \
+ size_t(_J+_VC) * size_t(_NI+2*_VC) \
+ size_t(_I+_VC) \
)

/** 4次元インデクス(i,j,k,n) -> 1次元インデクス変換マクロ
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
#define _IDX_S4D(_I,_J,_K,_N,_NI,_NJ,_NK,_VC) \
( size_t(_N) * size_t(_NI+2*_VC) * size_t(_NJ+2*_VC) * size_t(_NK+2*_VC) \
+ _IDX_S3D(_I,_J,_K,_NI,_NJ,_NK,_VC) \
)

/** 3次元インデクス(i,j,k,3) -> 1次元インデクス変換マクロ
 *  @param[in] _I  i方向インデクス
 *  @param[in] _J  j方向インデクス
 *  @param[in] _K  k方向インデクス
 *  @param[in] _N  成分インデクス
 *  @param[in] _NI i方向インデクスサイズ
 *  @param[in] _NJ j方向インデクスサイズ
 *  @param[in] _NK k方向インデクスサイズ
 *  @param[in] _VC 仮想セル数
 */
#define _IDX_V3D(_I,_J,_K,_N,_NI,_NJ,_NK,_VC) (_IDX_S4D(_I,_J,_K,_N,_NI,_NJ,_NK,_VC))

/** 4次元インデクス(n,i,j,k) -> 1次元インデクス変換マクロ
 *  @param[in] _N  成分インデクス
 *  @param[in] _I  i方向インデクス
 *  @param[in] _J  j方向インデクス
 *  @param[in] _K  k方向インデクス
 *  @param[in] _NN 成分数
 *  @param[in] _NI i方向インデクスサイズ
 *  @param[in] _NJ j方向インデクスサイズ
 *  @param[in] _NK k方向インデクスサイズ
 *  @param[in] _VC 仮想セル数
 *  @return 1次元インデクス
 */
#define _IDX_S4DEX(_N,_I,_J,_K,_NN,_NI,_NJ,_NK,_VC) \
( size_t(_NN) * _IDX_S3D(_I,_J,_K,_NI,_NJ,_NK,_VC) \
+ size_t(_N) )

/** 3次元インデクス(3,i,j,k) -> 1次元インデクス変換マクロ
 *  @param[in] _N  成分インデクス
 *  @param[in] _I  i方向インデクス
 *  @param[in] _J  j方向インデクス
 *  @param[in] _K  k方向インデクス
 *  @param[in] _NI i方向インデクスサイズ
 *  @param[in] _NJ j方向インデクスサイズ
 *  @param[in] _NK k方向インデクスサイズ
 *  @param[in] _VC 仮想セル数
 */
#define _IDX_V3DEX(_N,_I,_J,_K,_NI,_NJ,_NK,_VC) (_IDX_S4DEX(_N,_I,_J,_K,3,_NI,_NJ,_NK,_VC))

/** 面フラグ */
enum cpm_FaceFlag
{
  X_MINUS = 0  ///< -X face
, X_PLUS  = 1  ///< +X face
, Y_MINUS = 2  ///< -Y face
, Y_PLUS  = 3  ///< +Y face
, Z_MINUS = 4  ///< -Z face
, Z_PLUS  = 5  ///< +Z face
};

/** 軸方向フラグ */
enum cpm_DirFlag
{
  X_DIR = 0 ///< X direction
, Y_DIR = 1 ///< Y direction
, Z_DIR = 2 ///< Z direction
};

/** 方向フラグ */
enum cpm_PMFlag
{
  PLUS2MINUS = 0 ///< plus   -> minus direction
, MINUS2PLUS = 1 ///< minus  -> plus  direction
, BOTH       = 2 ///< plus  <-> minus direction
};

/** CPMのエラーコード */
enum cpm_ErrorCode
{
  CPM_SUCCESS                     = 0    ///< 正常終了

, CPM_ERROR                       = 1000 ///< その他のエラー
, CPM_ERROR_PM_INSTANCE           = 1001 ///< 並列管理クラスcpm_ParaManagerのインスタンス失敗
, CPM_ERROR_INVALID_PTR           = 1002 ///< ポインタのエラー
, CPM_ERROR_INVALID_DOMAIN_NO     = 1003 ///< 領域番号が不正
, CPM_ERROR_INVALID_OBJKEY        = 1004 ///< 指定登録番号のオブジェクトが存在しない
, CPM_ERROR_REGIST_OBJKEY         = 1005 ///< オブジェクト登録に失敗:

, CPM_ERROR_TEXTPARSER            = 2000 ///< テキストパーサーに関するエラー
, CPM_ERROR_NO_TEXTPARSER         = 2001 ///< テキストパーサーを組み込んでいない
, CPM_ERROR_TP_NOVECTOR           = 2002 ///< 領域分割情報ファイルのベクトルデータ読み込みエラー
, CPM_ERROR_TP_VECTOR_SIZE        = 2003 ///< 領域分割情報ファイルのベクトルデータのサイズが不正
, CPM_ERROR_TP_INVALID_G_ORG      = 2004 ///< 領域分割情報ファイルのドメイン原点情報が不正
, CPM_ERROR_TP_INVALID_G_VOXEL    = 2005 ///< 領域分割情報ファイルのドメインVOXEL数情報が不正
, CPM_ERROR_TP_INVALID_G_PITCH    = 2006 ///< 領域分割情報ファイルのドメインピッチ情報が不正
, CPM_ERROR_TP_INVALID_G_RGN      = 2007 ///< 領域分割情報ファイルのドメイン空間サイズ情報が不正
, CPM_ERROR_TP_INVALID_G_DIV      = 2008 ///< 領域分割情報ファイルのドメイン領域分割数情報が不正
, CPM_ERROR_TP_INVALID_POS        = 2009 ///< 領域分割情報ファイルのサブドメイン位置情報が不正

, CPM_ERROR_VOXELINIT             = 3000 ///< VoxelInitでエラー
, CPM_ERROR_NOT_IN_PROCGROUP      = 3001 ///< 自ランクがプロセスグループに含まれていない
, CPM_ERROR_ALREADY_VOXELINIIT    = 3002 ///< 指定されたプロセスグループが既に領域分割済み:
, CPM_ERROR_MISMATCH_NP_SUBDOMAIN = 3003 ///< 並列数とサブドメイン数が一致していない
, CPM_ERROR_CREATE_RANKMAP        = 3004 ///< ランクマップ生成に失敗
, CPM_ERROR_CREATE_NEIGHBOR       = 3005 ///< 隣接ランク情報生成に失敗
, CPM_ERROR_CREATE_LOCALDOMAIN    = 3006 ///< ローカル領域情報生成に失敗
, CPM_ERROR_INSERT_VOXELMAP       = 3007 ///< 領域情報のマップへの登録失敗
, CPM_ERROR_CREATE_PROCGROUP      = 3008 ///< プロセスグループ生成に失敗
, CPM_ERROR_INVALID_VOXELSIZE     = 3009 ///< VOXEL数が不正
, CPM_ERROR_INVALID_REGION        = 3010 ///< 全体空間サイズが不正
, CPM_ERROR_INVALID_DIVNUM        = 3011 ///< 領域分割数が不正

, CPM_ERROR_GET_INFO              = 4000 ///< 情報取得系関数でエラー
, CPM_ERROR_GET_DIVNUM            = 4001 ///< 領域分割数の取得エラー
, CPM_ERROR_GET_PITCH             = 4002 ///< ピッチの取得エラー 
, CPM_ERROR_GET_GLOBALVOXELSIZE   = 4003 ///< 全体ボクセル数の取得エラー
, CPM_ERROR_GET_GLOBALORIGIN      = 4004 ///< 全体空間の原点の取得エラー
, CPM_ERROR_GET_GLOBALREGION      = 4005 ///< 全体空間サイズの取得エラー
, CPM_ERROR_GET_LOCALVOXELSIZE    = 4006 ///< 自ランクのボクセル数の取得エラー
, CPM_ERROR_GET_LOCALORIGIN       = 4007 ///< 自ランクの空間原点の取得エラー
, CPM_ERROR_GET_LOCALREGION       = 4008 ///< 自ランクの空間サイズの取得エラー
, CPM_ERROR_GET_DIVPOS            = 4009 ///< 自ランクの領域分割位置の取得エラー
, CPM_ERROR_GET_HEADINDEX         = 4011 ///< 始点インデクスの取得エラー
, CPM_ERROR_GET_TAILINDEX         = 4012 ///< 終点インデクスの取得エラー
, CPM_ERROR_GET_NEIGHBOR_RANK     = 4013 ///< 隣接ランク番号の取得エラー
, CPM_ERROR_GET_PERIODIC_RANK     = 4014 ///< 周期境界位置の隣接ランク番号の取得エラー

, CPM_ERROR_GET_MYRANK            = 4015 ///< ランク番号の取得エラー
, CPM_ERROR_GET_NUMRANK           = 4016 ///< ランク数の取得エラー

, CPM_ERROR_MPI                   = 9000 ///< MPIのエラー
, CPM_ERROR_NO_MPI_INIT           = 9001 ///< MPI_Initがコールされていない
, CPM_ERROR_MPI_BARRIER           = 9003 ///< MPI_Barrierでエラー
, CPM_ERROR_MPI_BCAST             = 9004 ///< MPI_Bcastでエラー
, CPM_ERROR_MPI_SEND              = 9005 ///< MPI_Sendでエラー
, CPM_ERROR_MPI_RECV              = 9006 ///< MPI_Recvでエラー
, CPM_ERROR_MPI_ISEND             = 9007 ///< MPI_Isendでエラー
, CPM_ERROR_MPI_IRECV             = 9008 ///< MPI_Irecvでエラー
, CPM_ERROR_MPI_WAIT              = 9009 ///< MPI_Waitでエラー
, CPM_ERROR_MPI_WAITALL           = 9010 ///< MPI_Waitallでエラー
, CPM_ERROR_MPI_ALLREDUCE         = 9011 ///< MPI_Allreduceでエラー
, CPM_ERROR_MPI_GATHER            = 9012 ///< MPI_Gatherでエラー
, CPM_ERROR_MPI_ALLGATHER         = 9013 ///< MPI_Allgatherでエラー
, CPM_ERROR_MPI_GATHERV           = 9014 ///< MPI_Gathervでエラー
, CPM_ERROR_MPI_ALLGATHERV        = 9015 ///< MPI_Allgathervでエラー
, CPM_ERROR_MPI_DIMSCREATE        = 9016 ///< MPI_Dims_createでエラー

, CPM_ERROR_BNDCOMM               = 9500 ///< BndCommでエラー
, CPM_ERROR_BNDCOMM_VOXELSIZE     = 9501 ///< VoxelSize取得でエラー
, CPM_ERROR_BNDCOMM_BUFFER        = 9502 ///< 袖通信バッファ取得でエラー
, CPM_ERROR_BNDCOMM_BUFFERLENGTH  = 9503 ///< 袖通信バッファサイズが足りない

, CPM_ERROR_PERIODIC              = 9600 ///< PeriodicCommでエラー
, CPM_ERROR_PERIODIC_INVALID_DIR  = 9601 ///< 不正な軸方向フラグが指定された
, CPM_ERROR_PERIODIC_INVALID_PM   = 9602 ///< 不正な正負方向フラグが指定された

, CPM_ERROR_MPI_INVALID_COMM      = 9100 ///< MPIコミュニケータが不正
, CPM_ERROR_MPI_INVALID_DATATYPE  = 9101 ///< 対応しない型が指定された
, CPM_ERROR_MPI_INVALID_OPERATOR  = 9102 ///< 対応しないオペレータが指定された
, CPM_ERROR_MPI_INVALID_REQUEST   = 9103 ///< 不正なリクエストが指定された
};

/** fortran用のデータタイプ */
enum CPM_Datatype
{
  CPM_CHAR               =  1 ///< char
, CPM_UNSIGNED_CHAR      =  2 ///< unsigned char
, CPM_BYTE               =  3 ///< byte(not support)
, CPM_SHORT              =  4 ///< short
, CPM_UNSIGNED_SHORT     =  5 ///< unsigned short
, CPM_INT                =  6 ///< int
, CPM_UNSIGNED           =  7 ///< unsigned
, CPM_LONG               =  8 ///< long
, CPM_UNSIGNED_LONG      =  9 ///< unsigned long
, CPM_FLOAT              = 10 ///< float
, CPM_DOUBLE             = 11 ///< double
, CPM_LONG_DOUBLE        = 12 ///< long double
#ifdef MPI_LONG_LONG_INT
, CPM_LONG_LONG_INT      = 13 ///< long long int
#endif
#ifdef MPI_LONG_LONG
, CPM_LONG_LONG          = 13 ///< long long
#endif
#ifdef MPI_UNSIGNED_LONG_LONG
, CPM_UNSIGNED_LONG_LONG = 51 ///< unsigned long long
#endif
, CPM_REAL               = 52 ///< REAL_TYPE
};

/** fortran用のオペレータ */
enum CPM_Op
{
  CPM_MAX    = 100 ///< 最大値
, CPM_MIN    = 101 ///< 最小値
, CPM_SUM    = 102 ///< 和
, CPM_PROD   = 103 ///< 積
, CPM_LAND   = 104 ///< 論理積
, CPM_BAND   = 105 ///< ビット演算の積
, CPM_LOR    = 106 ///< 論理和
, CPM_BOR    = 107 ///< ビット演算の和
, CPM_LXOR   = 108 ///< 排他的論理和
, CPM_BXOR   = 109 ///< ビット演算の排他的論理和
, CPM_MINLOC = 110 ///< 最大値と位置(not support)
, CPM_MAXLOC = 111 ///< 最小値と位置(not support)
};

#endif /* _CPM_DEFINE_H_ */
