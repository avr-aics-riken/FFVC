#ifndef _CIO_DFI_SPH_H_
#define _CIO_DFI_SPH_H_

/* #################################################################
 *
 * CIOlib - Cartesian Input / Output library
 *
 * Copyright (c) 2013 Advanced Institute for Computational Science, RIKEN.
 * All right reserved.
 *
 * #################################################################
 */

/** 
 * @file   cio_DFI_SPH.h
 * @brief  cio_DFI_SPH Class Header
 * @author kero    
 */

#include "cio_DFI.h"

using namespace std;


class cio_DFI_SPH : public cio_DFI {

public:
  /** data dims(scalar or vector) */
  typedef enum {_DATA_UNKNOWN=0, _SCALAR, _VECTOR} DataDims;

  /** data type(float or double) */
  typedef enum {_REAL_UNKNOWN=0, _FLOAT, _DOUBLE} RealType;

  /** コンストラクタ */
  cio_DFI_SPH();

  cio_DFI_SPH(cio_FileInfo F_Info, cio_FilePath F_Path, cio_Unit unit, cio_Domain domain, cio_MPI mpi,
              vector<cio_Slice> TSlice, vector<cio_Rank> RInfo)
  {
    DFI_Finfo  = F_Info; 
    DFI_Fpath  = F_Path;
    DFI_Unit   = unit;
    DFI_Domain = domain;
    DFI_MPI    = mpi;
    TimeSlice  = TSlice;
    RankInfo   = RInfo;
  };
  
  /**　デストラクタ */
  ~cio_DFI_SPH();

  /**
   * @brief read sph data
   * @param[in]  step         読込むstep番号
   * @param[in]  gc           仮想セル数
   * @param[in]  Gvoxel[3]    グローバルボクセルサイズ　
   * @param[in]  Gdivision[3] 領域分割数
   * @param[in]  head[3]      計算領域の開始位置
   * @param[in]  tail[3]      計算領域の終了位置
   * @param[out] val          フィールドデータポインタ
   */ 
  void ReadData(int step, int gc, 
                int Gvoxel[3], int Gdivision[3], int head[3], int tail[3],
                REAL_TYPE *val, REAL_TYPE &time);

  /**
   * @brief write sph data
   * @param[in]  step 出力step番号
   * @param[in]  time 出力step番号
   * @param[in]  val  フィールドデータポインタ
   */ 
  void WriteData(int step, int gc, REAL_TYPE time,
                REAL_TYPE *val, REAL_TYPE *minmax);

  /**
   * @brief sphファイルのヘッダーレコード読込み
   * @param[in]  fp         ファイルポインタ
   * @param[in]  Etype      Endian Type
   * @param[in]  step       ステップ番号
   * @param[out] voxsize[3] voxsize
   * @return true:出力成功 false:出力失敗
   */
  bool read_Head(FILE* fp, int Etype, int step, INT_TYPE voxsize[3] ,
                 REAL_TYPE &time);

  /**
   * @brief sph S3D ファイルの読込み 
   * @param[in]  fname sphファイル名
   * @param[in]  step  読込むstep番号
   * @param[in]  gc    仮想セル数
   * @param[out] val   フィールドデータポインタ
   * @return true:出力成功 false:出力失敗
   */
  bool read_S3D(const char* fname,
                int         step,
                int         gc,
                REAL_TYPE*  val,
                REAL_TYPE&  time
                );

  /**
   * @brief sph V3D ファイルの読込み 
   * @param[in]  fname sphファイル名
   * @param[in]  step  読込むstep番号
   * @param[in]  sz    サイズ
   * @param[in]  gc    仮想セル数
   * @param[out] val   フィールドデータポインタ
   * @return true:出力成功 false:出力失敗
   */
  bool read_V3D(const char* fname,
                int         step,
                int*        sz,
                int         gc,
                REAL_TYPE*  val
                );

  /**
   * @brief 密データの読込みコピー（同一粒子）
   * @param[in]  fname  sphファイル名
   * @param[in]  step   読込むstep番号
   * @param[in]  gc     仮想セル数
   * @param[in]  head   自領域の起点
   * @param[in]  tail   自領域の終点
   * @param[in]  R_head 読み込む領域の起点
   * @param[in]  R_tail 読み込む領域の終点
   * @param[in]  sta[3] start   
   * @param[in]  end[3] end     
   * @param[out] val    フィールドデータポインタ
   * @return true:出力成功 false:出力失敗
   */ 
  bool read_MxN(const char* fname,
                int         step,
                int         gc,
                int         head[3],
                int         tail[3],
                int         R_head[3],
                int         R_tail[3],
                int         sta[3],
                int         end[3],
                REAL_TYPE*  val,
                REAL_TYPE&  time,
                int         n);

  /**
   * @brief 粗い粒子、１対１の読込みコピー&補間
   * @param[in]  fname     sphファイル名
   * @param[in]  step      読込むstep番号
   * @param[in]  gc        仮想セル数（自）
   * @param[in]  head      開始インデックス（自）
   * @param[in]  tail      終了インデックス（自）
   * @param[in]  R_head    開始インデックス（dfi）
   * @param[in]  R_tail    終了インデックス（dif）
   * @param[in]  Voxelsize ボクセルサイズ（dif）
   * @param[in]  dfi_gc    仮想セル数（dfi）
   * @param[in]  ncomp     コンポーネント数
   * @param[in]  sta[3]    start
   * @param[in]  end[3]    end
   * @param[out] val       フィールドデータポインタ
   * @return true:出力成功 false:出力失敗
   */ 
  bool read_Coarse(const char* fname,
                int         step,
                int         gc,
                int         head[3],
                int         tail[3],
                int         G_Voxelsize[3],
                int         R_head[3],
                int         R_tail[3],
                int         Voxelsize[3],
                int         dfi_gc,
                int         ncomp,
                int         sta[3],
                int         end[3],
                REAL_TYPE*  val,
                REAL_TYPE&  time);

  /**
   * @brief 粗い粒子、１対多の読込みコピー&補間
   * @param[in]  fname     sphファイル名
   * @param[in]  step      読込むstep番号
   * @param[in]  gc        仮想セル数（自）
   * @param[in]  head      開始インデックス（自）
   * @param[in]  tail      終了インデックス（自）
   * @param[in]  R_head    開始インデックス（dfi）
   * @param[in]  R_tail    終了インデックス（dif）
   * @param[in]  Voxelsize ボクセルサイズ（dif）
   * @param[in]  dfi_gc    仮想セル数（dfi）
   * @param[in]  ncomp     コンポーネント数
   * @param[in]  sta[3]    start
   * @param[in]  end[3]    end
   * @param[out] val       フィールドデータポインタ
   * @return true:出力成功 false:出力失敗
   */ 
  bool read_Coarse_MxN(const char* fname,
                int         step,
                int         gc,
                int         head[3],
                int         tail[3],
                int         G_Voxelsize[3],
                int         R_head[3],
                int         R_tail[3],
                int         Voxelsize[3],
                int         dfi_gc,
                int         ncomp,
                int         sta[3],
                int         end[3],
                REAL_TYPE*  val,
                REAL_TYPE&  time,
                int         n);
  /**
   * @brief SPHヘッダファイルの出力
   * @param[in] fp     ファイルポインタ
   * @param[in] step   ステップ番号
   * @param[in] time   時刻
   * @param[in] RankID ランク番号
   * @return true:出力成功 false:出力失敗
   */
  //bool write_head(FILE* fp, int step, REAL_TYPE time, int gc, int RankID); 
  bool write_head(FILE* fp, int step, REAL_TYPE time, int RankID); 

  /**
   * @brief SPHデータ出力
   * @param[in]  fp ファイルポインタ
   * @param[in]  val データポインタ
   * @param[in]  RankID ランク番号
   * @return true:出力成功 false:出力失敗
   */
  bool write_data(FILE* fp, REAL_TYPE* val, int gc, int RankID); 

};

#endif // _cio_DFI_SPH_H_
