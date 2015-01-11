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
// Copyright (c) 2012-2015 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################
//
///
/// @file  RLE.h
/// @brief ランレングスによる圧縮/展開ライブラリ
///

#ifndef __FFV_RLE_H__
#define __FFV_RLE_H__

#include <cstdlib>
#include <cstring>

namespace BVX_IO {
	
	/// RLE圧縮
	/// 
	/// @param[in]  source     入力データの先頭ポインタ
	/// @param[in]  sourceSize 入力データのサイズ (Byte単位で指定)
	/// @param[out] destSize   出力データのサイズ (Byte単位で指定)
	/// @return     RLE圧縮符号の先頭ポインタ．エラーの場合NULLを返す
	///
	/// @note returnされたポインタは適宜解放(delete)してください．
	///
	template <typename rluint_t, typename runlen_t>
	unsigned char* rleEncode(const rluint_t* source, const size_t sourceSize, size_t* destSize)
	{
	#ifdef __GNUC__
	#pragma pack(push, 1)
	#define ALIGNMENT __attribute__((packed))
	#else
	#pragma pack(1)
	#define ALIGNMENT
	#endif // __GNUC__
		// RLE符号を簡単に走査するためアライメントを無効にした構造体を定義
		struct DR{
			rluint_t d;
			runlen_t len;
		} ALIGNMENT;
	#ifdef __GNUC__
	#pragma pack(pop)
	#else  // __GNUC__
	#pragma pack()
	#endif // __GNUC__
	
		const runlen_t maxCount = (runlen_t)~0;
		const size_t endData = sourceSize / sizeof(rluint_t);
		
		const rluint_t* pSrc = source;
	
		size_t maxSize = sourceSize * sizeof(rluint_t) + sourceSize * sizeof(runlen_t);
		unsigned char*  dest = new unsigned char[maxSize];
		memset(dest, 0, maxSize * sizeof(unsigned char));
	
		DR* pdr = reinterpret_cast<DR*>(dest);
	
		pdr[0].d   = pSrc[0];
		pdr[0].len = 1;
		size_t cnt = 0;
		for(size_t i = 1; i < endData; i++){
			const rluint_t d = pSrc[i];
			if( pdr[cnt].d != d || pdr[cnt].len == maxCount )
			{
				cnt++;
				pdr[cnt].d = d;
				pdr[cnt].len = 1;
			}
			else
			{
				pdr[cnt].len++;
			}
		}
	
		*destSize = (cnt + 1) * sizeof(DR);
	
		return dest;
	}
	
	
	/// RLE展開
	/// 
	/// @param[in]  source     入力データの先頭ポインタ (RLE圧縮符号)
	/// @param[in]  sourceSize 入力データのサイズ (Byte単位で指定)
	/// @param[out] destSize   出力データのサイズ (Byte単位で指定)
	/// @return     RLE圧縮符号を展開したデータのの先頭ポインタ
	///
	/// @note returnされたポインタは適宜解放(delete)してください．
	///
	template <typename rluint_t, typename runlen_t>
	rluint_t* rleDecode(const unsigned char* source, const size_t sourceSize, const size_t destSize)
	{
	
	#ifdef __GNUC__
	#pragma pack(push, 1)
	#define ALIGNMENT __attribute__((packed))
	#else
	#pragma pack(1)
	#define ALIGNMENT
	#endif // __GNUC__
		// RLE符号を簡単に走査するためアライメントを無効にした構造体を定義
		struct DR{
			rluint_t d;
			runlen_t len;
		} ALIGNMENT;
	#ifdef __GNUC__
	#pragma pack(pop)
	#else  // __GNUC__
	#pragma pack()
	#endif // __GNUC__
		
		size_t endData = destSize / sizeof(rluint_t);
	
		rluint_t* dest = new rluint_t[endData];
	
		const DR* pdr = reinterpret_cast<const DR*>(source);
		size_t num = sourceSize / sizeof(DR);
		
		size_t cnt = 0;
		for(size_t i = 0; i < num; i++){
			for(runlen_t l = 0; l < pdr->len; l++){
				if(cnt >= endData){
					delete [] dest;
					return NULL;
				}
				dest[cnt] = pdr->d;
				cnt++;
			}
			pdr++;
		}
	
		return dest;
	}

} // BVX_IO

#endif // __FFV_RLE_H__

