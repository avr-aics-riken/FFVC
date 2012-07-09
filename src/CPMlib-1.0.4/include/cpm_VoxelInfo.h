/*
 * CPM - Cartesian Partition Manager
 *
 * Copyright (C) 2012 University of Tokyo.
 *
 */

/**
 * @file   cpm_VoxelInfo.h
 * VOXEL空間情報クラスのヘッダーファイル
 * @author University of Tokyo
 * @date   2012/05/31
 */

#ifndef _CPM_VOXELINFO_H_
#define _CPM_VOXELINFO_H_

#include "cpm_Base.h"
#include "cpm_DomainInfo.h"

/** CPMのVOXEL空間情報管理クラス
 */
class cpm_VoxelInfo : public cpm_Base
{
friend class cpm_ParaManager;
////////////////////////////////////////////////////////////////////////////////
// メンバー関数 
////////////////////////////////////////////////////////////////////////////////
public:


private:
  /** コンストラクタ */
  cpm_VoxelInfo();

  /** デストラクタ */
  virtual ~cpm_VoxelInfo();

  /** CPM領域分割情報の生成
   *  - MPI_COMM_WORLDを使用した領域を生成する。
   *  @param[in]  comm  MPIコミュニケータ
   *  @param[in]  dInfo 領域分割情報
   *  @param[in]  maxVC 最大の袖数(袖通信用)
   *  @param[in]  maxN  最大の成分数(袖通信用)
   *  @return 終了コード(CPM_SUCCESS=正常終了)
   */
  cpm_ErrorCode Init( MPI_Comm comm, cpm_GlobalDomainInfo* dInfo );

  /** ランクマップを生成
   *  @retval true  正常終了
   *  @retval false エラー
   */
  bool CreateRankMap();

  /** 隣接ランク情報を生成
   *  @retval true  正常終了
   *  @retval false エラー
   */
  bool CreateNeighborRankInfo();

  /** ローカル領域情報を生成
   *  @retval true  正常終了
   *  @retval false エラー
   */
  bool CreateLocalDomainInfo();

  /** 領域分割数を取得
   *  @return 領域分割数整数配列のポインタ
   */
  const int* GetDivNum() const;

  /** ピッチを取得
   *  @return ピッチ実数配列のポインタ
   */
  const REAL_TYPE* GetPitch() const;

  /** 全体ボクセル数を取得
   *  @return 全体ボクセル数整数配列のポインタ
   */
  const int* GetGlobalVoxelSize() const;

  /** 全体空間の原点を取得
   *  @return 全体空間の原点実数配列のポインタ
   */
  const REAL_TYPE* GetGlobalOrigin() const;

  /** 全体空間サイズを取得
   *  @return 全体空間サイズ実数配列のポインタ
   */
  const REAL_TYPE* GetGlobalRegion() const;

  /** 自ランクのボクセル数を取得
   *  @return 自ランクのボクセル数整数配列のポインタ
   */
  const int* GetLocalVoxelSize() const;

  /** 自ランクの空間原点を取得
   *  @return 自ランクの空間原点実数配列のポインタ
   */
  const REAL_TYPE* GetLocalOrigin() const;

  /** 自ランクの空間サイズを取得
   *  @return 自ランクの空間サイズ実数配列のポインタ
   */
  const REAL_TYPE* GetLocalRegion() const;

  /** 自ランクの領域分割位置を取得
   *  @return 自ランクの領域分割位置整数配列のポインタ
   */
  const int* GetDivPos() const;

  /** 自ランクの始点VOXELの全体空間でのインデクスを取得
   *  @return 自ランクの始点インデクス整数配列のポインタ
   */
  const int* GetVoxelHeadIndex() const;

  /** 自ランクの終点VOXELの全体空間でのインデクスを取得
   *  @return 自ランクの終点インデクス整数配列のポインタ
   */
  const int* GetVoxelTailIndex() const;

  /** 自ランクの隣接ランク番号を取得
   *  @return 自ランクの隣接ランク番号整数配列のポインタ
   */
  const int* GetNeighborRankID() const;

  /** 自ランクの周期境界の隣接ランク番号を取得
   *  @return 自ランクの周期境界の隣接ランク番号整数配列のポインタ
   */
  const int* GetPeriodicRankID() const;





////////////////////////////////////////////////////////////////////////////////
// メンバー変数
////////////////////////////////////////////////////////////////////////////////
public:


private:
  /**** 全体空間の情報 ****/
  cpm_GlobalDomainInfo m_globalDomainInfo; ///< 空間全体の領域情報
 
  /**** 自ランクの情報 ****/
  cpm_LocalDomainInfo m_localDomainInfo;   ///< 自ランクの領域情報
  int m_voxelHeadIndex[3];                 ///< 自ランクの始点ボクセルインデックス
  int m_voxelTailIndex[3];                 ///< 自ランクの終点ボクセルインデックス

  /**** 並列情報 ****/
  MPI_Comm  m_comm;        ///< MPIコミュニケータ
  int m_nRank;             ///< コミュニケータ内のランク数(=プロセス並列数)
  int m_rankNo;            ///< コミュニケータ内でのランク番号
  int m_neighborRankID[6]; ///< 隣接ランク番号(外部境界は負の値)
  int m_periodicRankID[6]; ///< 周期境界の隣接ランク番号

  int *m_rankMap; ///< ランクマップ

};

#endif /* _CPM_VOXELINFO_H_ */

