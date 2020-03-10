//##################################################################################
//
// FFV-C : Frontflow / violet Cartesian
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
//##################################################################################
//
///
/// @file  FileCommon.h
/// @brief ファイルIO用共通クラス群
///


#ifndef __FFV_FILE_HEADER_H__
#define __FFV_FILE_HEADER_H__

#include <limits.h>
#include <string>

#include "BitVoxel.h"

#include <stdint.h>

/// Octreeファイルのエンディアン識別子 (OC01)
#define OCTREE_FILE_IDENTIFIER    (('O' | ('C' << 8) | ('0' << 16) | ('1' << 24)))

/// LeafBlockファイルのエンディアン識別子 (LB01)
#define LEAFBLOCK_FILE_IDENTIFIER (('L' | ('B' << 8) | ('0' << 16) | ('1' << 24)))


namespace BVX_IO {

#ifdef __GNUC__
#pragma pack(push, 1)
#define ALIGNMENT __attribute__((packed))
#else
#pragma pack(1)
#define ALIGNMENT
#endif // __GNUC__
	
	/// Octreeファイルヘッダ構造体
	struct OctHeader
	{
		unsigned int identifier;   ///< エンディアン識別子
		double       org[3];       ///< 原点座標
		double       rgn[3];       ///< 領域サイズ
		unsigned int rootDims[3];  ///< ルート分割数
		unsigned int maxLevel;     ///< Octree最大分割レベル
		uint64_t     numLeaf;      ///< リーフノード数
		uint64_t     padding;      ///< 16バイトアライメント用パディング

		OctHeader() : padding(0) {}
	
	} ALIGNMENT;

	/// LeafBlockファイルヘッダ構造体
	struct LBHeader
	{
		unsigned int   identifier; ///< エンディアン識別子
		unsigned char  kind;       ///< ブロックファイル種類
		unsigned char  dataType;   ///< 1セルあたりのサイズ
		unsigned short bitWidth;   ///< 1セルあたりのビット幅
		unsigned int   vc;         ///< 仮想セルサイズ
		unsigned int   size[3];    ///< ブロックサイズ
		uint64_t       numBlock;   ///< ファイルに記載されている総ブロック数

	} ALIGNMENT;
	
	/// LeafBlockのCellIDヘッダ構造体
	struct LBCellIDHeader
	{
		uint64_t numBlock; ///< ブロック数
		uint64_t compSize; ///< 圧縮符号サイズ (バイト単位)

	} ALIGNMENT;

	/// RLE圧縮符号の走査用構造体
	struct GridRleCode
	{
		bitVoxelCell    c; ///< データ
		unsigned char len; ///< ラン長

	} ALIGNMENT;


#ifdef __GNUC__
#pragma pack(pop)
#else  // __GNUC__
#pragma pack()
#endif // __GNUC__

	/// ブロックデータタイプ
	enum LB_KIND
	{
		LB_CELLID = 0, ///< グリッド
		LB_SCALAR = 1, ///< スカラ
		LB_VECTOR = 3, ///< ベクター
		LB_TENSOR = 9, ///< テンソル
	};
	
	/// セルのデータ識別子
	enum LB_DATA_TYPE
	{
		LB_INT8    =  0, ///< 符号付き 8bit整数型
		LB_UINT8   =  1, ///< 符号なし 8bit整数型
		LB_INT16   =  2, ///< 符号付き16bit整数型
		LB_UINT16  =  3, ///< 符号なし16bit整数型
		LB_INT32   =  4, ///< 符号付き32bit整数型
		LB_UINT32  =  5, ///< 符号なし32bit整数型
		LB_INT64   =  6, ///< 符号付き64bit整数型
		LB_UINT64  =  7, ///< 符号なし64bit整数型
		LB_FLOAT32 =  8, ///< 32bit浮動小数点 (単精度浮動小数点)
		LB_FLOAT64 =  9  ///< 64bit浮動小数点 (倍精度浮動小数点)
	};

	
	/// 2byte用エンディアンスワップ
	static inline void BSwap16(void* a){
		unsigned short* x = (unsigned short*)a;
		*x = (unsigned short)( ((((*x) & 0xff00) >> 8 ) | (((*x) & 0x00ff) << 8)) );
	}

	/// 4byte用エンディアンスワップ
	static inline void BSwap32(void* a){
		unsigned int* x = (unsigned int*)a;
		*x = ( (((*x) & 0xff000000) >> 24) | (((*x) & 0x00ff0000) >> 8 ) |
		       (((*x) & 0x0000ff00) <<  8) | (((*x) & 0x000000ff) << 24) );
	}

	/// 8byte用エンディアンスワップ
	static inline void BSwap64(void* a){
		uint64_t* x = (uint64_t*)a;
		*x=( (((*x) & 0xff00000000000000ull) >> 56) | (((*x) & 0x00ff000000000000ull) >> 40) |
		     (((*x) & 0x0000ff0000000000ull) >> 24) | (((*x) & 0x000000ff00000000ull) >>  8) |
		     (((*x) & 0x00000000ff000000ull) <<  8) | (((*x) & 0x0000000000ff0000ull) << 24) |
		     (((*x) & 0x000000000000ff00ull) << 40) | (((*x) & 0x00000000000000ffull) << 56) );
	}


} // namespace BVX_IO

#endif // __FFV_FILE_HEADER_H__

