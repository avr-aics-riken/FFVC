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
/// @file  BlockSaver.h
/// @brief Blockファイルを出力する関数群
///

#ifndef __FFV_BLOCK_SAVER_H__
#define __FFV_BLOCK_SAVER_H__

#include "mpi.h"
#include <stdio.h>
#include <string>

#include "FileCommon.h"
#include "BitVoxel.h"
#include "RLE.h"
#include "FileSystemUtil.h"
#include "type.h"

namespace BVX_IO {
	
	/// Blockファイル(CellID)の出力
	/// 
	/// @param [in] size      格子サイズ
	/// @param [in] guide     ガイドセル
	/// @param [in] bitWidth  量子化するビット幅
	/// @param [in] rank      ランク番号
  /// @param [in] outDir    出力ディレクトリ
	/// @param [in] datas     CellIDが格納されたデータバッファ
	/// @param [in] rle       RLE圧縮フラグ (trueの場合RLE圧縮を行う)
	///
	/// @return 成功した場合true, 失敗した場合false
	/// 
  bool Save_Block_CellID(const int*            size,
                         const int             guide,
                         const unsigned        bitWidth,
                         const int             rank,
                         const std::string     outDir,
                         const unsigned char*  datas,
                         const bool            rle);
  
  
  /// Blockファイル(Scalar)の出力
  ///
  /// @param [in] size      格子サイズ
  /// @param [in] guide     ガイドセル
  /// @param [in] bitWidth  量子化するビット幅
  /// @param [in] rank      ランク番号
  /// @param [in] outDir    出力ディレクトリ
  /// @param [in] datas     BCflagが格納されたデータバッファ
  /// @param [in] rle       RLE圧縮フラグ (trueの場合RLE圧縮を行う)
  ///
  /// @return 成功した場合true, 失敗した場合false
  ///
  bool Save_Block_BCflag(const int*            size,
                         const int             guide,
                         const unsigned        bitWidth,
                         const int             rank,
                         const std::string     outDir,
                         const unsigned*       datas,
                         const bool            rle);
	
} // BVX_IO

#endif // __FFV_BLOCK_SAVER_H__

