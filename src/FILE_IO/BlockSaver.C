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
/// @file  BlockSaver.C
/// @brief Blockファイルを出力する関数群
///

#include "BlockSaver.h"
#include "mydebug.h"
#include "FileSystemUtil.h"

namespace BVX_IO {
  
  //##################################################################################
  bool Save_Block_CellID(const int*            size,
                         const int             guide,
                         const unsigned        bitWidth,
                         const int             rank,
                         const std::string     outDir,
                         const unsigned char*  datas,
                         const bool            rle)
  {
    using namespace std;
    
    // ブロックヘッダを準備
    LBHeader header;
    header.identifier = LEAFBLOCK_FILE_IDENTIFIER;
    header.kind       = static_cast<unsigned char>(LB_CELLID);
    header.dataType   = static_cast<unsigned char>(LB_UINT8);
    header.bitWidth   = static_cast<unsigned short>(bitWidth);
    header.vc         = guide;
    header.size[0]    = size[0];
    header.size[1]    = size[1];
    header.size[2]    = size[2];
    header.numBlock   = 0; // temporary
    
    int vc = guide;
    const size_t numBlock = 1; // FFV-Cの場合には1ブロック
    
    // 自プロセスの担当ブロックの保存用一時バッファサイズを計算
    const size_t tsz = (size[0] + vc*2) * (size[1] + vc*2) * (size[2] + vc*2) * numBlock;
    
    // 自プロセスの担当ブロックをBitVoxel化 >> bitVoxelはCompressBitVoxel()でnew
    size_t bitVoxelSize = 0;
    bitVoxelCell* bitVoxel = CompressBitVoxel(&bitVoxelSize, tsz, datas, bitWidth);
    
    // ブロックのCellIDヘッダを準備
    LBCellIDHeader ch;
    ch.numBlock = numBlock;
    
    // BitVoxelバッファのサイズ(Byte単位)を計算
    size_t bvs = bitVoxelSize * sizeof(bitVoxelCell);
    unsigned char* rleBuf = NULL;
    unsigned char* dp     = NULL;
    
    if( rle ){ // RLE圧縮の場合
      size_t rleSize = 0;
      
      // BitVoxelをRLE圧縮
      rleBuf = rleEncode<bitVoxelCell, unsigned char>(bitVoxel, bvs, &rleSize);
      
      // ブロックのCellIDヘッダに圧縮サイズを記載
      ch.compSize = rleSize;
      dp = rleBuf;
      
    }
    else // RLEなしの場合
    {
      // ブロックのCellIDヘッダに圧縮サイズ(=0)を記載
      ch.compSize = 0;
      dp = reinterpret_cast<unsigned char*>(bitVoxel);
    }
    
    
    // オリジナルには集約モードがあるが、FFV-C用には分散のみ
    // write file
    char filename[128];
    sprintf(filename, "cid_%06d.bvx", rank);
    
    string file_dir = outDir + "/CID";
    CheckDir(file_dir);
    
    string filepath = outDir + "/CID/" + string(filename);
    FILE *fp = NULL;
    if( (fp = fopen(filepath.c_str(), "wb")) == NULL) {
      stamped_printf("fileopen error <%s>. [%s:%d]\n", filepath.c_str(), __FILE__, __LINE__);
      return false;
    }
    
    header.numBlock = numBlock;
    
    fwrite(&header, sizeof(LBHeader),     1, fp);
    fwrite(&ch,     sizeof(LBCellIDHeader), 1, fp);
    
    size_t dsz = ch.compSize == 0 ? bvs : ch.compSize;
    
    fwrite(dp, sizeof(unsigned char), dsz, fp);
    
    fclose(fp);
    
    delete [] bitVoxel;
    if( rleBuf != NULL ) delete [] rleBuf;
    
    return true;
  }
  
  
  
  //##################################################################################
  bool Save_Block_BCflag(const int*            size,
                         const int             guide,
                         const unsigned        bitWidth,
                         const int             rank,
                         const std::string     outDir,
                         const unsigned*       datas,
                         const bool            rle)
  {
    using namespace std;
    
    // ブロックヘッダを準備
    LBHeader header;
    header.identifier = LEAFBLOCK_FILE_IDENTIFIER;
    header.kind       = static_cast<unsigned char>(LB_CELLID);
    header.dataType   = static_cast<unsigned char>(LB_UINT8);
    header.bitWidth   = static_cast<unsigned short>(bitWidth);
    header.vc         = guide;
    header.size[0]    = size[0];
    header.size[1]    = size[1];
    header.size[2]    = size[2];
    header.numBlock   = 0; // temporary
    
    int vc = guide;
    const size_t numBlock = 1; // FFV-Cの場合には1ブロック
    
    // 自プロセスの担当ブロックの保存用一時バッファサイズを計算
    const size_t tsz = (size[0] + vc*2) * (size[1] + vc*2) * (size[2] + vc*2) * numBlock;
    
    // datasは既に5bit x 6個に詰められているのでint -> uintへキャストするのみ
    size_t bitVoxelSize = tsz;
    bitVoxelCell* bitVoxel = (bitVoxelCell*)datas;
    
    // ブロックのCellIDヘッダを準備
    LBCellIDHeader ch;
    ch.numBlock = numBlock;
    
    // BitVoxelバッファのサイズ(Byte単位)を計算
    size_t bvs = bitVoxelSize * sizeof(bitVoxelCell);
    unsigned char* rleBuf = NULL;
    unsigned char* dp     = NULL;
    
    if( rle ){ // RLE圧縮の場合
      size_t rleSize = 0;
      
      // BitVoxelをRLE圧縮
      rleBuf = rleEncode<bitVoxelCell, unsigned char>(bitVoxel, bvs, &rleSize);
      
      // ブロックのCellIDヘッダに圧縮サイズを記載
      ch.compSize = rleSize;
      dp = rleBuf;
      
    }
    else // RLEなしの場合
    {
      // ブロックのCellIDヘッダに圧縮サイズ(=0)を記載
      ch.compSize = 0;
      dp = reinterpret_cast<unsigned char*>(bitVoxel);
    }
    
    
    // オリジナルには集約モードがあるが、FFV-C用には分散のみ
    // write file
    char filename[128];
    sprintf(filename, "bcf_%06d.bvx", rank);
    
    string file_dir = outDir + "/BCflag";
    CheckDir(file_dir);
    
    string filepath = outDir + "/BCflag/" + string(filename);
    FILE *fp = NULL;
    if( (fp = fopen(filepath.c_str(), "wb")) == NULL) {
      stamped_printf("fileopen error <%s>. [%s:%d]\n", filepath.c_str(), __FILE__, __LINE__);
      return false;
    }
    
    header.numBlock = numBlock;
    
    fwrite(&header, sizeof(LBHeader),     1, fp);
    fwrite(&ch,     sizeof(LBCellIDHeader), 1, fp);
    
    size_t dsz = ch.compSize == 0 ? bvs : ch.compSize;
    
    fwrite(dp, sizeof(unsigned char), dsz, fp);
    
    fclose(fp);
    
    //delete [] bitVoxel; >> datasは借り物
    if( rleBuf != NULL ) delete [] rleBuf;
    
    return true;
  }

} // namespace BVX_IO
