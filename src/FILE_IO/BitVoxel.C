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
// Copyright (c) 2016-2020 Research Institute for Information Research, Kyushu University.
// All rights reserved.
//
//##################################################################################
//
///
/// @file  BitVoxel.C
/// @brief ビットボクセル圧縮/展開ライブラリ
///

#include "BitVoxel.h"
#include <cstring>

namespace BVX_IO {

  size_t GetBitVoxelSize(const size_t voxelSize, const unsigned char bitWidth)
  {
    const unsigned char vox_per_cell = (sizeof(bitVoxelCell) * 8) / bitWidth;
    return (voxelSize / vox_per_cell + (voxelSize % vox_per_cell == 0 ? 0 : 1 ));
  }
  
  
  bitVoxelCell* CompressBitVoxel( size_t* bitVoxelSize, const size_t voxelSize, const unsigned char* voxel, const unsigned char bitWidth)
  {
    const unsigned char vox_per_cell = (sizeof(bitVoxelCell) * 8) / bitWidth;
    size_t bsz  = voxelSize / vox_per_cell + (voxelSize % vox_per_cell == 0 ? 0 : 1);
    
    bitVoxelCell* bitVoxel = new bitVoxelCell[bsz];
    memset(bitVoxel, 0, sizeof(bitVoxelCell) * bsz);
    
    unsigned char mask = 0;
    for(int i = 0; i < bitWidth; i++) mask += (1 << i);
    
    for(size_t i = 0; i < voxelSize; i++){
      size_t       cellIdx =  i / vox_per_cell;
      unsigned int bitIdx  = (i % vox_per_cell) * bitWidth;
      
      unsigned char c = voxel[i];
      
      bitVoxel[cellIdx] += (c & mask) << bitIdx;
    }
    
    *bitVoxelSize = bsz;
    return bitVoxel;
    
  }
  
  
  unsigned char* DecompressBitVoxel( const size_t voxelSize, const bitVoxelCell* bitVoxel, const unsigned char bitWidth)
  {
    const unsigned char vox_per_cell = (sizeof(bitVoxelCell) * 8) / bitWidth;
    
    unsigned char* voxel = new unsigned char[voxelSize];
    memset(voxel, 0, sizeof(unsigned char) * voxelSize);
    
    unsigned char mask = 0;
    for(int i = 0; i < bitWidth; i++) mask += (1 << i);
    
    for(size_t i = 0; i < voxelSize; i++){
      size_t       cellIdx =  i / vox_per_cell;
      unsigned int bitIdx  = (i % vox_per_cell) * bitWidth;
      
      voxel[i] = (bitVoxel[cellIdx] >> bitIdx) & mask;
    }
    
    return voxel;
  }
}
