#ifndef _CIO_DEFINE_H_
#define _CIO_DEFINE_H_

/* #################################################################
 *
 * CIOlib - Cartesian Input / Output library
 *
 * Copyright (c) AICS, RIKEN. All right reserved. 2013
 *
 * #################################################################
 */

/**
 * @file   cio_Define.h
 * @brief  CIOの定義マクロ記述ヘッダーファイル
 * @author kero
 */


#include "mpi.h"


//data type
#ifdef _REAL_IS_DOUBLE_
#define REAL_TYPE double
#define INT_TYPE long long
#else
#define REAL_TYPE float
#define INT_TYPE int
#endif

/** CIOのエラーコード */
enum cio_ErrorCode
{
  CIO_SUCCESS  = 1      ///<正常終了
, CIO_ERROR    = -1     ///<エラー終了
};

/** Endian check */
enum cio_EMatchType
{
  CIO_UnKnown = 0       ///<未定 
, CIO_Match   = 1       ///<一致
, CIO_UnMatch = 2       ///<不一致
};

/** 粗密判定コード */
enum cio_EGlobalVoxel
{
  CIO_E_GV_SAME
, CIO_E_GVX2_SAME
, CIO_E_OTHER
};

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
#define _IDX_IJK(_I,_J,_K,_NI,_NJ,_NK,_VC) \
( (long long)(_K+_VC) * (long long)(_NI+2*_VC) * (long long)(_NJ+2*_VC) \
+ (long long)(_J+_VC) * (long long)(_NI+2*_VC) \
+ (long long)(_I+_VC) \
)

/** 2次元インデクス(i,j) -> 1次元インデクス変換マクロ
 *  @param[in] _I  i方向インデクス
 *  @param[in] _J  j方向インデクス
 *  @param[in] _NI i方向インデクスサイズ
 *  @param[in] _NJ j方向インデクスサイズ
 *  @param[in] _VC 仮想セル数
 *  @return 1次元インデクス
 */
#define _IDX_IJ(_I,_J,_NI,_NJ,_VC) \
( (long long)(_J+_VC) * (long long)(_NI+2*_VC) \
+ (long long)(_I+_VC) \
)

#define _IDX_NIJ(_N,_I,_J,_NI,_NJ,_NN,_VC) \
( (long long)(_NN)*_IDX_IJ(_I,_J,_NI,_NJ,_VC) \
 + (long long)(_N) \
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
#define _IDX_IJKN(_I,_J,_K,_N,_NI,_NJ,_NK,_VC) \
( (long long)(_N) * (long long)(_NI+2*_VC) * (long long)(_NJ+2*_VC) * (long long)(_NK+2*_VC) \
+ _IDX_IJK(_I,_J,_K,_NI,_NJ,_NK,_VC) \
)

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
#define _IDX_NIJK(_N,_I,_J,_K,_NN,_NI,_NJ,_NK,_VC) \
( (long long)(_NN) * _IDX_IJK(_I,_J,_K,_NI,_NJ,_NK,_VC) \
+ (long long)(_N) )



#endif /* _CIO_DEFINE_H_ */
