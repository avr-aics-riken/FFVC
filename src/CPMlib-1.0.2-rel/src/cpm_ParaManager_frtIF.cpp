/*
 * CPM - Cartesian Partition Manager
 *
 * Copyright (C) 2012 University of Tokyo.
 *
 */

/**
 * @file   cpm_ParaManager_frtIF.cpp
 * パラレルマネージャクラスのFortranインターフェイスのソースファイル
 * @author University of Tokyo
 * @date   2012/05/31
 */
#include "cpm_ParaManager.h"

/** extern宣言 */
#define CPM_EXTERN extern "C"

#if 0
/** S3D,V3D通信でS4D版を使う */
#define _USE_S4D_
#endif

#ifndef CPM_WINDOWS
  #define cpm_Initialize_            cpm_initialize_
  #define cpm_VoxelInit_             cpm_voxelinit_
  #define cpm_VoxelInit_nodiv_       cpm_voxelinit_nodiv_
  #define cpm_IsParallel_            cpm_isparallel_
  #define cpm_GetDivNum_             cpm_getdivnum_
  #define cpm_GetPitch_              cpm_getpitch_
  #define cpm_GetGlobalVoxelSize_    cpm_getglobalvoxelsize_
  #define cpm_GetGlobalOrigin_       cpm_getglobalorigin_
  #define cpm_GetGlobalRegion_       cpm_getglobalregion_
  #define cpm_GetLocalVoxelSize_     cpm_getlocalvoxelsize_
  #define cpm_GetLocalOrigin_        cpm_getlocalorigin_
  #define cpm_GetLocalRegion_        cpm_getlocalregion_
  #define cpm_GetDivPos_             cpm_getdivpos_
  #define cpm_GetVoxelHeadIndex_     cpm_getvoxelheadindex_
  #define cpm_GetVoxelTailIndex_     cpm_getvoxeltailindex_
  #define cpm_GetNeighborRankID_     cpm_getneighborrankid_
  #define cpm_GetPeriodicRankID_     cpm_getperiodicrankid_
  #define cpm_GetMyRankID_           cpm_getmyrankid_
  #define cpm_GetNumRank_            cpm_getnumrank_
  #define cpm_Abort_                 cpm_abort_
  #define cpm_Barrier_               cpm_barrier_
  #define cpm_Wait_                  cpm_wait_
  #define cpm_Waitall_               cpm_waitall_
  #define cpm_Bcast_                 cpm_bcast_
  #define cpm_Send_                  cpm_send_
  #define cpm_Recv_                  cpm_recv_
  #define cpm_Isend_                 cpm_isend_
  #define cpm_Irecv_                 cpm_irecv_
  #define cpm_Allreduce_             cpm_allreduce_
  #define cpm_Gather_                cpm_gather_
  #define cpm_Allgather_             cpm_allgather_
  #define cpm_Gatherv_               cpm_gatherv_
  #define cpm_Allgatherv_            cpm_allgatherv_
  #define cpm_SetBndCommBuffer_      cpm_setbndcommbuffer_
  #define cpm_BndCommS3D_            cpm_bndcomms3d_
  #define cpm_BndCommV3D_            cpm_bndcommv3d_
  #define cpm_BndCommS4D_            cpm_bndcomms4d_
  #define cpm_BndCommS3D_nowait_     cpm_bndcomms3d_nowait_
  #define cpm_BndCommV3D_nowait_     cpm_bndcommv3d_nowait_
  #define cpm_BndCommS4D_nowait_     cpm_bndcomms4d_nowait_
  #define cpm_wait_BndCommS3D_       cpm_wait_bndcomms3d_
  #define cpm_wait_BndCommV3D_       cpm_wait_bndcommv3d_
  #define cpm_wait_BndCommS4D_       cpm_wait_bndcomms4d_
  #define cpm_BndCommV3DEx_          cpm_bndcommv3dex_
  #define cpm_BndCommS4DEx_          cpm_bndcomms4dex_
  #define cpm_BndCommV3DEx_nowait_   cpm_bndcommv3dex_nowait_
  #define cpm_BndCommS4DEx_nowait_   cpm_bndcomms4dex_nowait_
  #define cpm_wait_BndCommV3DEx_     cpm_wait_bndcommv3dex_
  #define cpm_wait_BndCommS4DEx_     cpm_wait_bndcomms4dex_
  #define cpm_PeriodicCommS3D        cpm_periodiccomms3d_
  #define cpm_PeriodicCommV3D        cpm_periodiccommv3d_
  #define cpm_PeriodicCommS4D        cpm_periodiccomms4d_
  #define cpm_PeriodicCommV3DEx      cpm_periodiccommv3dex_
  #define cpm_PeriodicCommS4DEx      cpm_periodiccomms4dex_
#else
  #define cpm_Initialize_            CPM_INITIALIZE
  #define cpm_VoxelInit_             CPM_VOXELINIT
  #define cpm_VoxelInit_nodiv_       CPM_VOXELINIT_NODIV
  #define cpm_IsParallel_            CPM_ISPARALLEL
  #define cpm_GetDivNum_             CPM_GETDIVNUM
  #define cpm_GetPitch_              CPM_GETPITCH
  #define cpm_GetGlobalVoxelSize_    CPM_GETGLOBALVOXELSIZE
  #define cpm_GetGlobalOrigin_       CPM_GETGLOBALORIGIN
  #define cpm_GetGlobalRegion_       CPM_GETGLOBALREGION
  #define cpm_GetLocalVoxelSize_     CPM_GETLOCALVOXELSIZE
  #define cpm_GetLocalOrigin_        CPM_GETLOCALORIGIN
  #define cpm_GetLocalRegion_        CPM_GETLOCALREGION
  #define cpm_GetDivPos_             CPM_GETDIVPOS
  #define cpm_GetVoxelHeadIndex_     CPM_GETVOXELHEADINDEX
  #define cpm_GetVoxelTailIndex_     CPM_GETVOXELTAILINDEX
  #define cpm_GetNeighborRankID_     CPM_GETNEIGHBORRANKID
  #define cpm_GetPeriodicRankID_     CPM_GETPERIODICRANKID
  #define cpm_GetMyRankID_           CPM_GETMYRANKID
  #define cpm_GetNumRank_            CPM_GETNUMRANK
  #define cpm_Abort_                 CPM_ABORT
  #define cpm_Barrier_               CPM_BARRIER
  #define cpm_Wait_                  CPM_WAIT
  #define cpm_Waitall_               CPM_WAITALL
  #define cpm_Bcast_                 CPM_BCAST
  #define cpm_Send_                  CPM_SEND
  #define cpm_Recv_                  CPM_RECV
  #define cpm_Isend_                 CPM_ISEND
  #define cpm_Irecv_                 CPM_IRECV
  #define cpm_Allreduce_             CPM_ALLREDUCE
  #define cpm_Gather_                CPM_GATHER
  #define cpm_Allgather_             CPM_ALLGATHER
  #define cpm_Gatherv_               CPM_GATHERV
  #define cpm_Allgatherv_            CPM_ALLGATHERV
  #define cpm_SetBndCommBuffer_      CPM_SETBNDCOMMBUFFER
  #define cpm_BndCommS3D_            CPM_BNDCOMMS3D
  #define cpm_BndCommV3D_            CPM_BNDCOMMV3D
  #define cpm_BndCommS4D_            CPM_BNDCOMMS4D
  #define cpm_BndCommS3D_nowait_     CPM_BNDCOMMS3D_NOWAIT
  #define cpm_BndCommV3D_nowait_     CPM_BNDCOMMV3D_NOWAIT
  #define cpm_BndCommS4D_nowait_     CPM_BNDCOMMS4D_NOWAIT
  #define cpm_wait_BndCommS3D_       CPM_WAIT_BNDCOMMS3D
  #define cpm_wait_BndCommV3D_       CPM_WAIT_BNDCOMMV3D
  #define cpm_wait_BndCommS4D_       CPM_WAIT_BNDCOMMS4D
  #define cpm_BndCommV3DEx_          CPM_BNDCOMMV3DEX
  #define cpm_BndCommS4DEx_          CPM_BNDCOMMS4DEX
  #define cpm_BndCommV3DEx_nowait_   CPM_BNDCOMMV3DEX_NOWAIT
  #define cpm_BndCommS4DEx_nowait_   CPM_BNDCOMMS4DEX_NOWAIT
  #define cpm_wait_BndCommV3DEx_     CPM_WAIT_BNDCOMMV3DEX
  #define cpm_wait_BndCommS4DEx_     CPM_WAIT_BNDCOMMS4DEX
  #define cpm_PeriodicCommS3D        CPM_PERIODICCOMMS3D
  #define cpm_PeriodicCommV3D        CPM_PERIODICCOMMV3D
  #define cpm_PeriodicCommS4D        CPM_PERIODICCOMMS4D
  #define cpm_PeriodicCommV3DEx      CPM_PERIODICCOMMV3DEX
  #define cpm_PeriodicCommS4DEx      CPM_PERIODICCOMMS4DEX
#endif

////////////////////////////////////////////////////////////////////////////////
/** 初期化処理(MPI_Initは実行済みの場合)
 *  - InitializeのFortranインターフェイス関数
 *  - FortranでMPI_Initがコールされている必要がある
 *  @param[out] ierr 終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_Initialize_( int *ierr )
{
  if( ierr )
  {
    return;
  }

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // Initialize
  *ierr = paraMngr->Initialize();
}

////////////////////////////////////////////////////////////////////////////////
/** 領域分割
 *  - VoxelInitのFortranインターフェイス関数
 *  - 領域分割の各種情報を引数で渡して領域分割を行う
 *  - プロセスグループの全てのランクが活性ドメインになる
 *  - 領域分割数を指定する
 *  @param[in]  div       領域分割数(サイズ3)
 *  @param[in]  vox       空間全体のボクセル数(サイズ3)
 *  @param[in]  origin    空間全体の原点(サイズ3)
 *  @param[in]  pitch     ボクセルピッチ(サイズ3)
 *  @param[in]  maxVC     最大の袖数(袖通信用)
 *  @param[in]  maxN      最大の成分数(袖通信用)
 *  @param[in]  procGrpNo 領域分割を行うプロセスグループ番号
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_VoxelInit_( int *div, int *vox, REAL_TYPE *origin, REAL_TYPE *pitch
              , int *maxVC, int *maxN, int *procGrpNo, int *ierr )
{
  if( !div || !vox || !origin || !pitch || !maxVC || !maxN ||
      !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // VoxelInit
  *ierr = paraMngr->VoxelInit( div, vox, origin, pitch
                             , size_t(*maxVC), size_t(*maxN), *procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
/** 領域分割
 *  - VoxelInitのFortranインターフェイス関数
 *  - 領域分割の各種情報を引数で渡して領域分割を行う
 *  - プロセスグループの全てのランクが活性ドメインになる
 *  - プロセスグループのランク数で自動領域分割
 *  @param[in]  vox       空間全体のボクセル数(サイズ3)
 *  @param[in]  origin    空間全体の原点(サイズ3)
 *  @param[in]  pitch     ボクセルピッチ(サイズ3)
 *  @param[in]  maxVC     最大の袖数(袖通信用)
 *  @param[in]  maxN      最大の成分数(袖通信用)
 *  @param[in]  procGrpNo 領域分割を行うプロセスグループ番号
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_VoxelInit_nodiv_( int *vox, REAL_TYPE *origin, REAL_TYPE *pitch
                    , int *maxVC, int *maxN, int *procGrpNo, int *ierr )
{
  if( !vox || !origin || !pitch || !maxVC || !maxN ||
      !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // VoxelInit
  *ierr = paraMngr->VoxelInit( vox, origin, pitch
                             , size_t(*maxVC), size_t(*maxN), *procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
/** 並列実行であるかチェックする
 *  - IsParallelのFortranインターフェイス関数
 *  @param[out] ipara 並列実行フラグ(1=並列実行、1以外=逐次実行)
 *  @param[out] ierr  終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_IsParallel_( int *ipara, int *ierr )
{
  if( !ipara || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  *ipara = 0;

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // IsParallel
  if( paraMngr->IsParallel() )
    *ipara = 1;
  else
    *ipara = 0;
}

////////////////////////////////////////////////////////////////////////////////
/** 領域分割数を取得
 *  - GetDivNumのFortranインターフェイス関数
 *  @param[out] div       領域分割数(3wordの整数配列)
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_GetDivNum_( int *div, int *procGrpNo, int *ierr )
{
  if( !procGrpNo || !div || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  div[0] = div[1] = div[2] = 0;

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  const int *dd = paraMngr->GetDivNum( *procGrpNo );
  if( !dd )
  {
    *ierr = CPM_ERROR_GET_DIVNUM;
    return;
  }

  div[0] = dd[0];
  div[1] = dd[1];
  div[2] = dd[2];

  *ierr = CPM_SUCCESS;
  return;
}

////////////////////////////////////////////////////////////////////////////////
/** ピッチを取得
 *  - GetPitchのFortranインターフェイス関数
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] pch       ピッチ(3wordの実数配列)
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_GetPitch_( REAL_TYPE *pch, int *procGrpNo, int *ierr )
{
  if( !procGrpNo || !pch || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  pch[0] = pch[1] = pch[2] = REAL_TYPE(0);

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  const REAL_TYPE *pp = paraMngr->GetPitch( *procGrpNo );
  if( !pp )
  {
    *ierr = CPM_ERROR_GET_PITCH;
    return;
  }

  pch[0] = pp[0];
  pch[1] = pp[1];
  pch[2] = pp[2];

  *ierr = CPM_SUCCESS;
  return;
}

////////////////////////////////////////////////////////////////////////////////
/** 全体ボクセル数を取得
 *  - GetGlobalVoxelSizeのFortranインターフェイス関数
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] wsz       全体ボクセル数(3wordの整数配列)
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_GetGlobalVoxelSize_( int *wsz, int *procGrpNo, int *ierr )
{
  if( !procGrpNo || !wsz || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  wsz[0] = wsz[1] = wsz[2] = 0;

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  const int *sz = paraMngr->GetGlobalVoxelSize( *procGrpNo );
  if( !sz )
  {
    *ierr = CPM_ERROR_GET_GLOBALVOXELSIZE;
    return;
  }

  wsz[0] = sz[0];
  wsz[1] = sz[1];
  wsz[2] = sz[2];

  *ierr = CPM_SUCCESS;
  return;
}

////////////////////////////////////////////////////////////////////////////////
/** 全体空間の原点を取得
 *  - GetGlobalOriginのFortranインターフェイス関数
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] worg      全体空間の原点(3wordの実数配列)
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_GetGlobalOrigin_( REAL_TYPE *worg, int *procGrpNo, int *ierr )
{
  if( !procGrpNo || !worg || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  worg[0] = worg[1] = worg[2] = REAL_TYPE(0);

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  const REAL_TYPE *org = paraMngr->GetGlobalOrigin( *procGrpNo );
  if( !org )
  {
    *ierr = CPM_ERROR_GET_GLOBALORIGIN;
    return;
  }

  worg[0] = org[0];
  worg[1] = org[1];
  worg[2] = org[2];

  *ierr = CPM_SUCCESS;
  return;
}

////////////////////////////////////////////////////////////////////////////////
/** 全体空間サイズを取得
 *  - GetGlobalRegionのFortranインターフェイス関数
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] wrgn      全体空間サイズ(3wordの実数配列)
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_GetGlobalRegion_( REAL_TYPE *wrgn, int *procGrpNo, int *ierr )
{
  if( !procGrpNo || !wrgn || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  wrgn[0] = wrgn[1] = wrgn[2] = REAL_TYPE(0);

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  const REAL_TYPE *rgn = paraMngr->GetGlobalRegion( *procGrpNo );
  if( !rgn )
  {
    *ierr = CPM_ERROR_GET_GLOBALREGION;
    return;
  }

  wrgn[0] = rgn[0];
  wrgn[1] = rgn[1];
  wrgn[2] = rgn[2];

  *ierr = CPM_SUCCESS;
  return;
}

////////////////////////////////////////////////////////////////////////////////
/** 自ランクのボクセル数を取得
 *  - GetLocalVoxelSizeのFortranインターフェイス関数
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] lsz       自ランクのボクセル数(3wordの整数配列)
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_GetLocalVoxelSize_( int *lsz, int *procGrpNo, int *ierr )
{
  if( !procGrpNo || !lsz || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  lsz[0] = lsz[1] = lsz[2] = 0;

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  const int *sz = paraMngr->GetLocalVoxelSize( *procGrpNo );
  if( !sz )
  {
    *ierr = CPM_ERROR_GET_LOCALVOXELSIZE;
    return;
  }

  lsz[0] = sz[0];
  lsz[1] = sz[1];
  lsz[2] = sz[2];

  *ierr = CPM_SUCCESS;
  return;
}

////////////////////////////////////////////////////////////////////////////////
/** 自ランクの空間原点を取得
 *  - GetLocalOriginのFortranインターフェイス関数
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] lorg      自ランクの空間原点(3wordの実数配列)
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_GetLocalOrigin_( REAL_TYPE *lorg, int *procGrpNo, int *ierr )
{
  if( !procGrpNo || !lorg || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  lorg[0] = lorg[1] = lorg[2] = REAL_TYPE(0);

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  const REAL_TYPE *org = paraMngr->GetLocalOrigin( *procGrpNo );
  if( !org )
  {
    *ierr = CPM_ERROR_GET_LOCALORIGIN;
    return;
  }

  lorg[0] = org[0];
  lorg[1] = org[1];
  lorg[2] = org[2];

  *ierr = CPM_SUCCESS;
  return;
}

////////////////////////////////////////////////////////////////////////////////
/** 自ランクの空間サイズを取得
 *  - GetLocalRegionのFortranインターフェイス関数
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] lrgn      自ランクの空間サイズ(3wordの実数配列)
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_GetLocalRegion_( REAL_TYPE *lrgn, int *procGrpNo, int *ierr )
{
  if( !procGrpNo || !lrgn || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  lrgn[0] = lrgn[1] = lrgn[2] = REAL_TYPE(0);

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  const REAL_TYPE *rgn = paraMngr->GetLocalOrigin( *procGrpNo );
  if( !rgn )
  {
    *ierr = CPM_ERROR_GET_LOCALREGION;
    return;
  }

  lrgn[0] = rgn[0];
  lrgn[1] = rgn[1];
  lrgn[2] = rgn[2];

  *ierr = CPM_SUCCESS;
  return;
}

////////////////////////////////////////////////////////////////////////////////
/** 自ランクの領域分割位置を取得
 *  - GetDivPosのFortranインターフェイス関数
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] pos       自ランクの領域分割位置(3wordの整数配列)
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_GetDivPos_( int *pos, int *procGrpNo, int *ierr )
{
  if( !procGrpNo || !pos || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  pos[0] = pos[1] = pos[2] = 0;

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  const int *pp = paraMngr->GetDivPos( *procGrpNo );
  if( !pp )
  {
    *ierr = CPM_ERROR_GET_DIVPOS;
    return;
  }

  pos[0] = pp[0];
  pos[1] = pp[1];
  pos[2] = pp[2];

  *ierr = CPM_SUCCESS;
  return;
}

////////////////////////////////////////////////////////////////////////////////
/** 自ランクの始点VOXELの全体空間でのインデクスを取得
 *  - GetVoxelHeadIndexのFortranインターフェイス関数
 *  - 全体空間の先頭インデクスを0としたC型のインデクス
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] idx       自ランクの始点VOXELインデクス(3wordの整数配列)
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_GetVoxelHeadIndex_( int *idx, int *procGrpNo, int *ierr )
{
  if( !procGrpNo || !idx || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

 idx[0] = idx[1] = idx[2] = 0;

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  const int *id = paraMngr->GetVoxelHeadIndex( *procGrpNo );
  if( !id )
  {
    *ierr = CPM_ERROR_GET_HEADINDEX;
    return;
  }

  idx[0] = id[0];
  idx[1] = id[1];
  idx[2] = id[2];

  *ierr = CPM_SUCCESS;
  return;
}

////////////////////////////////////////////////////////////////////////////////
/** 自ランクの終点VOXELの全体空間でのインデクスを取得
 *  - GetVoxelTailIndexのFortranインターフェイス関数
 *  - 全体空間の先頭インデクスを0としたC型のインデクス
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] idx       自ランクの終点VOXELインデクス(3wordの整数配列)
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_GetVoxelTailIndex_( int *idx, int *procGrpNo, int *ierr )
{
  if( !procGrpNo || !idx || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

 idx[0] = idx[1] = idx[2] = 0;

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  const int *id = paraMngr->GetVoxelTailIndex( *procGrpNo );
  if( !id )
  {
    *ierr = CPM_ERROR_GET_TAILINDEX;
    return;
  }

  idx[0] = id[0];
  idx[1] = id[1];
  idx[2] = id[2];

  *ierr = CPM_SUCCESS;
  return;
}

////////////////////////////////////////////////////////////////////////////////
/** 自ランクの隣接ランク番号を取得
 *  - GetNeighborRankIDのFortranインターフェイス関数
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] nID       自ランクの隣接ランク番号(6wordの整数配列)
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_GetNeighborRankID_( int *nID, int *procGrpNo, int *ierr )
{
  if( !procGrpNo || !nID || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

 nID[0] = nID[1] = nID[2] = nID[3] = nID[4] = nID[5] = 0;

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  const int *id = paraMngr->GetNeighborRankID( *procGrpNo );
  if( !id )
  {
    *ierr = CPM_ERROR_GET_NEIGHBOR_RANK;
    return;
  }

  nID[0] = id[0];
  nID[1] = id[1];
  nID[2] = id[2];
  nID[3] = id[3];
  nID[4] = id[4];
  nID[5] = id[5];

  *ierr = CPM_SUCCESS;
  return;
}

////////////////////////////////////////////////////////////////////////////////
/** 自ランクの周期境界の隣接ランク番号を取得
 *  - GetPeriodicRankIDのFortranインターフェイス関数
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] nID       自ランクの周期境界の隣接ランク番号6wordの整数配列)
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_GetPeriodicRankID_( int *nID, int *procGrpNo, int *ierr )
{
  if( !procGrpNo || !nID || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

 nID[0] = nID[1] = nID[2] = nID[3] = nID[4] = nID[5] = 0;

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  const int *id = paraMngr->GetPeriodicRankID( *procGrpNo );
  if( !id )
  {
    *ierr = CPM_ERROR_GET_PERIODIC_RANK;
    return;
  }

  nID[0] = id[0];
  nID[1] = id[1];
  nID[2] = id[2];
  nID[3] = id[3];
  nID[4] = id[4];
  nID[5] = id[5];

  *ierr = CPM_SUCCESS;
  return;
}

////////////////////////////////////////////////////////////////////////////////
/** ランク番号の取得
 *  - GetMyRankIDのFortranインターフェイス関数
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] id        ランク番号
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_GetMyRankID_( int *id, int *procGrpNo, int *ierr )
{
  if( !procGrpNo || !id || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  *id = 0;

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  *id = paraMngr->GetMyRankID( *procGrpNo );
  if( !id )
  {
    *ierr = CPM_ERROR_GET_MYRANK;
    return;
  }

  *ierr = CPM_SUCCESS;
  return;
}

////////////////////////////////////////////////////////////////////////////////
/** ランク数の取得
 *  - GetNumRankのFortranインターフェイス関数
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] nrank     ランク数
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_GetNumRank_( int *nrank, int *procGrpNo, int *ierr )
{
  if( !procGrpNo || !nrank || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  *nrank = 1;

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  *nrank = paraMngr->GetNumRank( *procGrpNo );
  if( !nrank )
  {
    *ierr = CPM_ERROR_GET_NUMRANK;
    return;
  }

  *ierr = CPM_SUCCESS;
  return;
}

////////////////////////////////////////////////////////////////////////////////
/** Abort
 *  - AbortのFortranインターフェイス関数
 *  @param[in]  errorcode MPI_Abortに渡すエラーコード
 */
CPM_EXTERN
void
cpm_Abort_( int *errorcode )
{
  int err = 0;
  if( errorcode )
  {
    err = *errorcode;
  }

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    MPI_Abort( MPI_COMM_WORLD, err );
    exit(err);
    return;
  }

  paraMngr->Abort( err );
}

////////////////////////////////////////////////////////////////////////////////
/** Barrier
 *  - BarrierのFortranインターフェイス関数
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_Barrier_( int *procGrpNo, int *ierr )
{
  if( !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // Barrier
  *ierr = paraMngr->Barrier( *procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
/** Wait
 *  - WaitのFortranインターフェイス関数
 *  @param[in]  reqNo     リクエスト番号(0以上の整数)
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_Wait_( int *reqNo, int *ierr )
{
  if( !reqNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // cpm_Wait
  *ierr = paraMngr->cpm_Wait( *reqNo );
}

////////////////////////////////////////////////////////////////////////////////
/** Waitall
 *  - WaitallのFortranインターフェイス関数
 *  @param[in]  count     リクエストの数
 *  @param[in]  reqlist   リクエスト番号のリスト(0以上の整数)
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_Waitall_( int *count, int *reqlist, int *ierr )
{
  if( !count || !reqlist || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // cpm_Waitall
  *ierr = paraMngr->cpm_Waitall( *count, reqlist );
}

////////////////////////////////////////////////////////////////////////////////
/** Bcast
 *  - BcastのFortranインターフェイス関数
 *  @param[inout] buf       送受信バッファ
 *  @param[in]    count     送信バッファのサイズ(ワード数)
 *  @param[in]    datatype  データタイプ(fparam.fiを参照)
 *  @param[in]    root      送信元のランク番号(procGrpNo内でのランク番号)
 *  @param[in]    procGrpNo プロセスグループ番号
 *  @param[out]   ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_Bcast_( void *buf, int *count, int *datatype, int *root, int *procGrpNo, int *ierr )
{
  if( !buf || !count || !datatype || !root || !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // MPI_Datatype
  MPI_Datatype dtype = cpm_ParaManager::GetMPI_Datatype( *datatype );
  if( dtype == MPI_DATATYPE_NULL )
  {
    *ierr = CPM_ERROR_MPI_INVALID_DATATYPE;
    return;
  }

  // Bcast
  *ierr = paraMngr->Bcast( dtype, buf, *count, *root, *procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
/** Send
 *  - SendのFortranインターフェイス関数
 *  @param[inout] buf       送信バッファ
 *  @param[in]    count     送信バッファのサイズ(ワード数)
 *  @param[in]    datatype  データタイプ(fparam.fiを参照)
 *  @param[in]    dest      送信先のランク番号(procGrpNo内でのランク番号)
 *  @param[in]    procGrpNo プロセスグループ番号
 *  @param[out]   ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_Send_( void *buf, int *count, int *datatype, int *dest, int *procGrpNo, int *ierr )
{
  if( !buf || !count || !datatype || !dest || !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // MPI_Datatype
  MPI_Datatype dtype = cpm_ParaManager::GetMPI_Datatype( *datatype );
  if( dtype == MPI_DATATYPE_NULL )
  {
    *ierr = CPM_ERROR_MPI_INVALID_DATATYPE;
    return;
  }

  // Send
  *ierr = paraMngr->Send( dtype, buf, *count, *dest, *procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
/** Recv
 *  - RecvのFortranインターフェイス関数
 *  @param[inout] buf       受信バッファ
 *  @param[in]    count     受信バッファのサイズ(ワード数)
 *  @param[in]    datatype  データタイプ(fparam.fiを参照)
 *  @param[in]    source    送信元のランク番号(procGrpNo内でのランク番号)
 *  @param[in]    procGrpNo プロセスグループ番号
 *  @param[out]   ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_Recv_( void *buf, int *count, int *datatype, int *source, int *procGrpNo, int *ierr )
{
  if( !buf || !count || !datatype || !source || !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // MPI_Datatype
  MPI_Datatype dtype = cpm_ParaManager::GetMPI_Datatype( *datatype );
  if( dtype == MPI_DATATYPE_NULL )
  {
    *ierr = CPM_ERROR_MPI_INVALID_DATATYPE;
    return;
  }

  // Recv
  *ierr = paraMngr->Recv( dtype, buf, *count, *source, *procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
/** Isend
 *  - IsendのFortranインターフェイス関数
 *  @param[inout] buf       送信バッファ
 *  @param[in]    count     送信バッファのサイズ(ワード数)
 *  @param[in]    datatype  データタイプ(fparam.fiを参照)
 *  @param[in]    dest      送信先のランク番号(procGrpNo内でのランク番号)
 *  @param[in]    procGrpNo プロセスグループ番号
 *  @param[out]   reqNo     リクエスト番号
 *  @param[out]   ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_Isend_( void *buf, int *count, int *datatype, int *dest, int *procGrpNo, int *reqNo, int *ierr )
{
  if( !buf || !count || !datatype || !dest || !procGrpNo || !reqNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // cpm_Isend
  if( (*ierr = paraMngr->cpm_Isend( buf, *count, *datatype, *dest, reqNo, *procGrpNo ) ) != CPM_SUCCESS )
  {
    return;
  }

  *ierr = CPM_SUCCESS;
  return;
}

////////////////////////////////////////////////////////////////////////////////
/** Irecv
 *  - IrecvのFortranインターフェイス関数
 *  @param[inout] buf       受信バッファ
 *  @param[in]    count     受信バッファのサイズ(ワード数)
 *  @param[in]    datatype  データタイプ(fparam.fiを参照)
 *  @param[in]    source    送信元先のランク番号(procGrpNo内でのランク番号)
 *  @param[in]    procGrpNo プロセスグループ番号
 *  @param[in]    reqNo     リクエスト番号
 *  @param[out]   ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_Irecv_( void *buf, int *count, int *datatype, int *source, int *procGrpNo, int *reqNo, int *ierr )
{
  if( !buf || !count || !datatype || !source || !procGrpNo || !reqNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // cpm_Irecv
  if( (*ierr = paraMngr->cpm_Irecv( buf, *count, *datatype, *source, reqNo, *procGrpNo ) ) != CPM_SUCCESS )
  {
    return;
  }

  *ierr = CPM_SUCCESS;
  return;
}

////////////////////////////////////////////////////////////////////////////////
/** MPI_AllreduceのFortranインターフェイス
 *  - MPI_AllreduceのFortranインターフェイス関数
 *  @param[in]  sendbuf   送信データ
 *  @param[out] recvbuf   受信データ
 *  @param[in]  count     送受信データのサイズ
 *  @param[in]  datatype  送受信データのデータタイプ(cpm_fparam.fiを参照)
 *  @param[in]  op        オペレータ
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_Allreduce_( void *sendbuf, void *recvbuf, int *count, int *datatype, int *op, int *procGrpNo, int *ierr )
{
  if( !sendbuf || !recvbuf || !count || !datatype || !op || !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // MPI_Datatype
  MPI_Datatype dtype = cpm_ParaManager::GetMPI_Datatype( *datatype );
  if( dtype == MPI_DATATYPE_NULL )
  {
    *ierr = CPM_ERROR_MPI_INVALID_DATATYPE;
    return;
  }

  // MPI_Op
  MPI_Op ope = cpm_ParaManager::GetMPI_Op( *op );
  if( ope == MPI_OP_NULL )
  {
    *ierr = CPM_ERROR_MPI_INVALID_OPERATOR;
    return;
  }

  // Allreduce
  *ierr = paraMngr->Allreduce( dtype, sendbuf, recvbuf, *count, ope, *procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
/** MPI_GatherのFortranインターフェイス
 *  - MPI_GatherのFortranインターフェイス関数
 *  @param[in]  sendbuf   送信データ
 *  @param[in]  sendcnt   送信データのサイズ
 *  @param[in]  sendtype  送信データのデータタイプ(cpm_fparam.fiを参照)
 *  @param[out] recvbuf   受信データ
 *  @param[in]  recvcnt   受信データのサイズ
 *  @param[in]  recvtype  受信データのデータタイプ(cpm_fparam.fiを参照)
 *  @param[in]  root      受信するランク番号(procGrpNo内でのランク番号)
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_Gather_( void *sendbuf, int *sendcnt, int *sendtype, void *recvbuf, int *recvcnt, int *recvtype
           , int *root, int *procGrpNo, int *ierr )
{
  if( !sendbuf || !sendcnt || !sendtype || !recvbuf || !recvcnt || !recvtype ||
      !root || !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // MPI_Datatype
  MPI_Datatype stype = cpm_ParaManager::GetMPI_Datatype( *sendtype );
  if( stype == MPI_DATATYPE_NULL )
  {
    *ierr = CPM_ERROR_MPI_INVALID_DATATYPE;
    return;
  }
  MPI_Datatype rtype = cpm_ParaManager::GetMPI_Datatype( *recvtype );
  if( rtype == MPI_DATATYPE_NULL )
  {
    *ierr = CPM_ERROR_MPI_INVALID_DATATYPE;
    return;
  }

  // Gather
  *ierr = paraMngr->Gather( stype, sendbuf, *sendcnt, rtype, recvbuf, *recvcnt, *root, *procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
/** MPI_AllgatherのFortranインターフェイス
 *  - MPI_AllgatherのFortranインターフェイス関数
 *  @param[in]  sendbuf   送信データ
 *  @param[in]  sendcnt   送信データのサイズ
 *  @param[in]  sendtype  送信データのデータタイプ(cpm_fparam.fiを参照)
 *  @param[out] recvbuf   受信データ
 *  @param[in]  recvcnt   受信データのサイズ
 *  @param[in]  recvtype  受信データのデータタイプ(cpm_fparam.fiを参照)
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_Allgather_( void *sendbuf, int *sendcnt, int *sendtype, void *recvbuf, int *recvcnt, int *recvtype
              , int *procGrpNo, int *ierr )
{
  if( !sendbuf || !sendcnt || !sendtype || !recvbuf || !recvcnt || !recvtype || !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // MPI_Datatype
  MPI_Datatype stype = cpm_ParaManager::GetMPI_Datatype( *sendtype );
  if( stype == MPI_DATATYPE_NULL )
  {
    *ierr = CPM_ERROR_MPI_INVALID_DATATYPE;
    return;
  }
  MPI_Datatype rtype = cpm_ParaManager::GetMPI_Datatype( *recvtype );
  if( rtype == MPI_DATATYPE_NULL )
  {
    *ierr = CPM_ERROR_MPI_INVALID_DATATYPE;
    return;
  }

  // Allgather
  *ierr = paraMngr->Allgather( stype, sendbuf, *sendcnt, rtype, recvbuf, *recvcnt, *procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
/** MPI_GathervのFortranインターフェイス
 *  - MPI_GathervのFortranインターフェイス関数
 *  @param[in]  sendbuf   送信データ
 *  @param[in]  sendcnt   送信データのサイズ
 *  @param[in]  sendtype  送信データのデータタイプ(cpm_fparam.fiを参照)
 *  @param[out] recvbuf   受信データ
 *  @param[in]  recvcnts  各ランクからの受信データサイズ
 *  @param[in]  displs    各ランクからの受信データ配置位置
 *  @param[in]  recvtype  受信データのデータタイプ(cpm_fparam.fiを参照)
 *  @param[in]  root      受信するランク番号(procGrpNo内でのランク番号)
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_Gatherv_( void *sendbuf, int *sendcnt, int *sendtype, void *recvbuf, int *recvcnts, int *displs, int *recvtype
            , int *root, int *procGrpNo, int *ierr )
{
  if( !sendbuf || !sendcnt || !sendtype || !recvbuf || !recvcnts || !displs || !recvtype ||
      !root || !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // MPI_Datatype
  MPI_Datatype stype = cpm_ParaManager::GetMPI_Datatype( *sendtype );
  if( stype == MPI_DATATYPE_NULL )
  {
    *ierr = CPM_ERROR_MPI_INVALID_DATATYPE;
    return;
  }
  MPI_Datatype rtype = cpm_ParaManager::GetMPI_Datatype( *recvtype );
  if( rtype == MPI_DATATYPE_NULL )
  {
    *ierr = CPM_ERROR_MPI_INVALID_DATATYPE;
    return;
  }

  // Gatherv
  *ierr = paraMngr->Gatherv( stype, sendbuf, *sendcnt, rtype, recvbuf, recvcnts, displs, *root, *procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
/** MPI_AllgathervのFortranインターフェイス
 *  - MPI_AllgathervのFortranインターフェイス関数
 *  @param[in]  sendbuf   送信データ
 *  @param[in]  sendcnt   送信データのサイズ
 *  @param[in]  sendtype  送信データのデータタイプ(cpm_fparam.fiを参照)
 *  @param[out] recvbuf   受信データ
 *  @param[in]  recvcnts  各ランクからの受信データサイズ
 *  @param[in]  displs    各ランクからの受信データ配置位置
 *  @param[in]  recvtype  受信データのデータタイプ(cpm_fparam.fiを参照)
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] ierr      終了コード(0=正常終了、0以外=cpm_ErrorCodeの値)
 */
CPM_EXTERN
void
cpm_Allgatherv_( void *sendbuf, int *sendcnt, int *sendtype, void *recvbuf, int *recvcnts, int *displs, int *recvtype
               , int *procGrpNo, int *ierr )
{
  if( !sendbuf || !sendcnt || !sendtype || !recvbuf || !recvcnts || !displs || !recvtype ||
      !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // MPI_Datatype
  MPI_Datatype stype = cpm_ParaManager::GetMPI_Datatype( *sendtype );
  if( stype == MPI_DATATYPE_NULL )
  {
    *ierr = CPM_ERROR_MPI_INVALID_DATATYPE;
    return;
  }
  MPI_Datatype rtype = cpm_ParaManager::GetMPI_Datatype( *recvtype );
  if( rtype == MPI_DATATYPE_NULL )
  {
    *ierr = CPM_ERROR_MPI_INVALID_DATATYPE;
    return;
  }

  // Gatherv
  *ierr = paraMngr->Allgatherv( stype, sendbuf, *sendcnt, rtype, recvbuf, recvcnts, displs, *procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
/** 袖通信バッファのセット(Fortranインターフェイス)
 *  - 袖通信バッファ確保処理のFortranインターフェイス関数
 *  @param[in]  maxVC     送受信バッファの最大袖数
 *  @param[in]  maxN      送受信バッファの最大成分数
 *  @param[in]  procGrpNo プロセスグループ番号
 *  @param[out] ierr      終了コード(CPM_SUCCESS=正常終了)
 */
CPM_EXTERN
void
cpm_SetBndCommBuffer_( int *maxVC, int *maxN, int *procGrpNo, int *ierr )
{
  if( !maxVC || !maxN || !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // SetBndCommBuffer
  *ierr = paraMngr->SetBndCommBuffer( size_t(*maxVC), size_t(*maxN), *procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
/** 袖通信(Scalar4D版)のFortranインターフェイス
 *  - (imax,jmax,kmax,nmax)の形式の配列の袖通信を行う
 *  - BndCommS4DのFortranインターフェイス関数
 *  @param[inout] array     袖通信をする配列の先頭ポインタ
 *  @param[in]    imax      配列サイズ(I方向)
 *  @param[in]    jmax      配列サイズ(J方向)
 *  @param[in]    kmax      配列サイズ(K方向)
 *  @param[in]    nmax      配列サイズ(成分数)
 *  @param[in]    vc        仮想セル数
 *  @param[in]    vc_comm   通信する仮想セル数
 *  @param[in]    datatype  袖通信データのデータタイプ(cpm_fparam.fiを参照)
 *  @param[in]    procGrpNo プロセスグループ番号
 *  @param[out]   ierr      終了コード(CPM_SUCCESS=正常終了)
 */
CPM_EXTERN
void
cpm_BndCommS4D_( void *array, int *imax, int *jmax, int *kmax, int *nmax, int *vc, int *vc_comm
               , int *datatype, int *procGrpNo, int *ierr )
{
  if( !array || !imax || !jmax || !kmax || !nmax || !vc || !vc_comm || !datatype || !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // MPI_Datatype
  MPI_Datatype dtype = cpm_ParaManager::GetMPI_Datatype( *datatype );
  if( dtype == MPI_DATATYPE_NULL )
  {
    *ierr = CPM_ERROR_MPI_INVALID_DATATYPE;
    return;
  }

  // BndCommS4D
   *ierr = paraMngr->BndCommS4D( dtype, array, *imax, *jmax, *kmax, *nmax, *vc, *vc_comm, *procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
/** 袖通信(Scalar3D版)のFortranインターフェイス
 *  - (imax,jmax,kmax)の形式の配列の袖通信を行う
 *  - BndCommS3DのFortranインターフェイス関数
 *  @param[inout] array     袖通信をする配列の先頭ポインタ
 *  @param[in]    imax      配列サイズ(I方向)
 *  @param[in]    jmax      配列サイズ(J方向)
 *  @param[in]    kmax      配列サイズ(K方向)
 *  @param[in]    vc        仮想セル数
 *  @param[in]    vc_comm   通信する仮想セル数
 *  @param[in]    datatype  袖通信データのデータタイプ(cpm_fparam.fiを参照)
 *  @param[in]    procGrpNo プロセスグループ番号
 *  @param[out]   ierr      終了コード(CPM_SUCCESS=正常終了)
 */
CPM_EXTERN
void
cpm_BndCommS3D_( void *array, int *imax, int *jmax, int *kmax, int *vc, int *vc_comm, int *datatype
               , int *procGrpNo, int *ierr )
{
#ifdef _USE_S4D_
  int nmax = 1;
  cpm_BndCommS4D_( array, imax, jmax, kmax, &nmax, vc, vc_comm, datatype, procGrpNo, ierr );
#else
  if( !array || !imax || !jmax || !kmax || !vc || !vc_comm || !datatype || !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // MPI_Datatype
  MPI_Datatype dtype = cpm_ParaManager::GetMPI_Datatype( *datatype );
  if( dtype == MPI_DATATYPE_NULL )
  {
    *ierr = CPM_ERROR_MPI_INVALID_DATATYPE;
    return;
  }

  // BndCommS3D
  *ierr = paraMngr->BndCommS3D( dtype, array, *imax, *jmax, *kmax, *vc, *vc_comm, *procGrpNo );
#endif
}

////////////////////////////////////////////////////////////////////////////////
/** 袖通信(Vector3D版)のFortranインターフェイス
 *  - (imax,jmax,kmax,3)の形式の配列の袖通信を行う
 *  - BndCommV3DのFortranインターフェイス関数
 *  @param[inout] array     袖通信をする配列の先頭ポインタ
 *  @param[in]    imax      配列サイズ(I方向)
 *  @param[in]    jmax      配列サイズ(J方向)
 *  @param[in]    kmax      配列サイズ(K方向)
 *  @param[in]    vc        仮想セル数
 *  @param[in]    vc_comm   通信する仮想セル数
 *  @param[in]    datatype  袖通信データのデータタイプ(cpm_fparam.fiを参照)
 *  @param[in]    procGrpNo プロセスグループ番号
 *  @param[out]   ierr      終了コード(CPM_SUCCESS=正常終了)
 */
CPM_EXTERN
void
cpm_BndCommV3D_( void *array, int *imax, int *jmax, int *kmax, int *vc, int *vc_comm, int *datatype
               , int *procGrpNo, int *ierr )
{
#ifdef _USE_S4D_
  int nmax = 3;
  cpm_BndCommS4D_( array, imax, jmax, kmax, &nmax, vc, vc_comm, datatype, procGrpNo, ierr );
#else
  if( !array || !imax || !jmax || !kmax || !vc || !vc_comm || !datatype || !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // MPI_Datatype
  MPI_Datatype dtype = cpm_ParaManager::GetMPI_Datatype( *datatype );
  if( dtype == MPI_DATATYPE_NULL )
  {
    *ierr = CPM_ERROR_MPI_INVALID_DATATYPE;
    return;
  }

  // BndCommV3D
  *ierr = paraMngr->BndCommV3D( dtype, array, *imax, *jmax, *kmax, *vc, *vc_comm, *procGrpNo );
#endif
}

////////////////////////////////////////////////////////////////////////////////
/** 非同期版袖通信(Scalar4D版)のFortranインターフェイス
 *  - (imax,jmax,kmax,nmax)の形式の配列の袖通信を行う
 *  - BndCommS4D_nowaitのFortranインターフェイス関数
 *  @param[inout] array     袖通信をする配列の先頭ポインタ
 *  @param[in]    imax      配列サイズ(I方向)
 *  @param[in]    jmax      配列サイズ(J方向)
 *  @param[in]    kmax      配列サイズ(K方向)
 *  @param[in]    nmax      配列サイズ(成分数)
 *  @param[in]    vc        仮想セル数
 *  @param[in]    vc_comm   通信する仮想セル数
 *  @param[in]    datatype  袖通信データのデータタイプ(cpm_fparam.fiを参照)
 *  @param[out]   reqlist   リクエスト番号のリスト(サイズ12)
 *  @param[in]    procGrpNo プロセスグループ番号
 *  @param[out]   ierr      終了コード(CPM_SUCCESS=正常終了)
 */
CPM_EXTERN
void
cpm_BndCommS4D_nowait_( void *array, int *imax, int *jmax, int *kmax, int *nmax
                      , int *vc, int *vc_comm, int *datatype, int *reqlist, int *procGrpNo, int *ierr )
{
  if( !array || !imax || !jmax || !kmax || !nmax || !vc || !vc_comm || !datatype ||
      !reqlist || !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // cpm_BndCommS4D_nowait
   *ierr = paraMngr->cpm_BndCommS4D_nowait( array, *imax, *jmax, *kmax, *nmax, *vc, *vc_comm
                                          , *datatype, reqlist, *procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
/** 非同期版袖通信(Scalar3D版)のFortranインターフェイス
 *  - (imax,jmax,kmax)の形式の配列の袖通信を行う
 *  - BndCommS3D_nowaitのFortranインターフェイス関数
 *  @param[inout] array     袖通信をする配列の先頭ポインタ
 *  @param[in]    imax      配列サイズ(I方向)
 *  @param[in]    jmax      配列サイズ(J方向)
 *  @param[in]    kmax      配列サイズ(K方向)
 *  @param[in]    vc        仮想セル数
 *  @param[in]    vc_comm   通信する仮想セル数
 *  @param[in]    datatype  袖通信データのデータタイプ(cpm_fparam.fiを参照)
 *  @param[out]   reqlist   リクエスト番号のリスト(サイズ12)
 *  @param[in]    procGrpNo プロセスグループ番号
 *  @param[out]   ierr      終了コード(CPM_SUCCESS=正常終了)
 */
CPM_EXTERN
void
cpm_BndCommS3D_nowait_( void *array, int *imax, int *jmax, int *kmax
                      , int *vc, int *vc_comm, int *datatype, int *reqlist, int *procGrpNo, int *ierr )
{
#ifdef _USE_S4D_
  int nmax = 1;
  cpm_BndCommS4D_nowait_( array, imax, jmax, kmax, &nmax, vc, vc_comm
                        , datatype, reqlist, procGrpNo, ierr );
#else
  if( !array || !imax || !jmax || !kmax || !vc || !vc_comm || !datatype ||
      !reqlist || !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // cpm_BndCommS3D_nowait
   *ierr = paraMngr->cpm_BndCommS3D_nowait( array, *imax, *jmax, *kmax, *vc, *vc_comm
                                          , *datatype, reqlist, *procGrpNo );
#endif
}

////////////////////////////////////////////////////////////////////////////////
/** 非同期版袖通信(Vector3D版)のFortranインターフェイス
 *  - (imax,jmax,kmax,3)の形式の配列の袖通信を行う
 *  - BndCommV3D_nowaitのFortranインターフェイス関数
 *  @param[inout] array     袖通信をする配列の先頭ポインタ
 *  @param[in]    imax      配列サイズ(I方向)
 *  @param[in]    jmax      配列サイズ(J方向)
 *  @param[in]    kmax      配列サイズ(K方向)
 *  @param[in]    vc        仮想セル数
 *  @param[in]    vc_comm   通信する仮想セル数
 *  @param[in]    datatype  袖通信データのデータタイプ(cpm_fparam.fiを参照)
 *  @param[out]   reqlist   リクエスト番号のリスト(サイズ12)
 *  @param[in]    procGrpNo プロセスグループ番号
 *  @param[out]   ierr      終了コード(CPM_SUCCESS=正常終了)
 */
CPM_EXTERN
void
cpm_BndCommV3D_nowait_( void *array, int *imax, int *jmax, int *kmax
                      , int *vc, int *vc_comm, int *datatype, int *reqlist, int *procGrpNo, int *ierr )
{
#ifdef _USE_S4D_
  int nmax = 3;
  cpm_BndCommS4D_nowait_( array, imax, jmax, kmax, &nmax, vc, vc_comm
                        , datatype, reqlist, procGrpNo, ierr );
#else
  if( !array || !imax || !jmax || !kmax || !vc || !vc_comm || !datatype ||
      !reqlist || !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // cpm_BndCommV3D_nowait
   *ierr = paraMngr->cpm_BndCommV3D_nowait( array, *imax, *jmax, *kmax, *vc, *vc_comm
                                          , *datatype, reqlist, *procGrpNo );
#endif
}

////////////////////////////////////////////////////////////////////////////////
/** 非同期版袖通信のwait、展開(Scalar4D版)のFortranインターフェイス
 *  - (imax,jmax,kmax,nmax)の形式の配列の非同期版袖通信のwaitと展開を行う
 *  - wait_BndCommS4DのFortranインターフェイス関数
 *  @param[inout] array     袖通信をする配列の先頭ポインタ
 *  @param[in]    imax      配列サイズ(I方向)
 *  @param[in]    jmax      配列サイズ(J方向)
 *  @param[in]    kmax      配列サイズ(K方向)
 *  @param[in]    nmax      配列サイズ(成分数)
 *  @param[in]    vc        仮想セル数
 *  @param[in]    vc_comm   通信する仮想セル数
 *  @param[in]    datatype  袖通信データのデータタイプ(cpm_fparam.fiを参照)
 *  @param[out]   reqlist   リクエスト番号のリスト(サイズ12)
 *  @param[in]    procGrpNo プロセスグループ番号
 *  @param[out]   ierr      終了コード(CPM_SUCCESS=正常終了)
 */
CPM_EXTERN
void
cpm_wait_BndCommS4D_( void *array, int *imax, int *jmax, int *kmax, int *nmax
                    , int *vc, int *vc_comm, int *datatype, int *reqlist, int *procGrpNo, int *ierr )
{
  if( !array || !imax || !jmax || !kmax || !nmax || !vc || !vc_comm || !datatype ||
      !reqlist || !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // cpm_wait_BndCommS4D
   *ierr = paraMngr->cpm_wait_BndCommS4D( array, *imax, *jmax, *kmax, *nmax, *vc, *vc_comm
                                        , *datatype, reqlist, *procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
/** 非同期版袖通信のwait、展開(Scalar3D版)のFortranインターフェイス
 *  - (imax,jmax,kmax)の形式の配列の非同期版袖通信のwaitと展開を行う
 *  - wait_BndCommS3DのFortranインターフェイス関数
 *  @param[inout] array     袖通信をする配列の先頭ポインタ
 *  @param[in]    imax      配列サイズ(I方向)
 *  @param[in]    jmax      配列サイズ(J方向)
 *  @param[in]    kmax      配列サイズ(K方向)
 *  @param[in]    vc        仮想セル数
 *  @param[in]    vc_comm   通信する仮想セル数
 *  @param[in]    datatype  袖通信データのデータタイプ(cpm_fparam.fiを参照)
 *  @param[out]   reqlist   リクエスト番号のリスト(サイズ12)
 *  @param[in]    procGrpNo プロセスグループ番号
 *  @param[out]   ierr      終了コード(CPM_SUCCESS=正常終了)
 */
CPM_EXTERN
void
cpm_wait_BndCommS3D_( void *array, int *imax, int *jmax, int *kmax
                    , int *vc, int *vc_comm, int *datatype, int *reqlist, int *procGrpNo, int *ierr )
{
#ifdef _USE_S4D_
  int nmax = 1;
  cpm_wait_BndCommS4D_( array, imax, jmax, kmax, &nmax, vc, vc_comm, datatype, reqlist, procGrpNo, ierr );
#else
  if( !array || !imax || !jmax || !kmax || !vc || !vc_comm || !datatype ||
      !reqlist || !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // cpm_wait_BndCommS3D
   *ierr = paraMngr->cpm_wait_BndCommS3D( array, *imax, *jmax, *kmax, *vc, *vc_comm
                                        , *datatype, reqlist, *procGrpNo );
#endif
}

////////////////////////////////////////////////////////////////////////////////
/** 非同期版袖通信のwait、展開(Vector3D版)のFortranインターフェイス
 *  - (imax,jmax,kmax,3)の形式の配列の非同期版袖通信のwaitと展開を行う
 *  - wait_BndCommV3DのFortranインターフェイス関数
 *  @param[inout] array     袖通信をする配列の先頭ポインタ
 *  @param[in]    imax      配列サイズ(I方向)
 *  @param[in]    jmax      配列サイズ(J方向)
 *  @param[in]    kmax      配列サイズ(K方向)
 *  @param[in]    vc        仮想セル数
 *  @param[in]    vc_comm   通信する仮想セル数
 *  @param[in]    datatype  袖通信データのデータタイプ(cpm_fparam.fiを参照)
 *  @param[out]   reqlist   リクエスト番号のリスト(サイズ12)
 *  @param[in]    procGrpNo プロセスグループ番号
 *  @param[out]   ierr      終了コード(CPM_SUCCESS=正常終了)
 */
CPM_EXTERN
void
cpm_wait_BndCommV3D_( void *array, int *imax, int *jmax, int *kmax
                    , int *vc, int *vc_comm, int *datatype, int *reqlist, int *procGrpNo, int *ierr )
{
#ifdef _USE_S4D_
  int nmax = 3;
  cpm_wait_BndCommS4D_( array, imax, jmax, kmax, &nmax, vc, vc_comm, datatype, reqlist, procGrpNo, ierr );
#else
  if( !array || !imax || !jmax || !kmax || !vc || !vc_comm || !datatype ||
      !reqlist || !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // cpm_wait_BndCommV3D
   *ierr = paraMngr->cpm_wait_BndCommV3D( array, *imax, *jmax, *kmax, *vc, *vc_comm
                                        , *datatype, reqlist, *procGrpNo );
#endif
}

////////////////////////////////////////////////////////////////////////////////
/** 袖通信(Scalar4DEx版)のFortranインターフェイス
 *  - (nmax,imax,jmax,kmax)の形式の配列の袖通信を行う
 *  - BndCommS4DExのFortranインターフェイス関数
 *  @param[inout] array     袖通信をする配列の先頭ポインタ
 *  @param[in]    nmax      配列サイズ(成分数)
 *  @param[in]    imax      配列サイズ(I方向)
 *  @param[in]    jmax      配列サイズ(J方向)
 *  @param[in]    kmax      配列サイズ(K方向)
 *  @param[in]    vc        仮想セル数
 *  @param[in]    vc_comm   通信する仮想セル数
 *  @param[in]    datatype  袖通信データのデータタイプ(cpm_fparam.fiを参照)
 *  @param[in]    procGrpNo プロセスグループ番号
 *  @param[out]   ierr      終了コード(CPM_SUCCESS=正常終了)
 */
CPM_EXTERN
void
cpm_BndCommS4DEx_( void *array, int *nmax, int *imax, int *jmax, int *kmax, int *vc, int *vc_comm
                 , int *datatype, int *procGrpNo, int *ierr )
{
  if( !array || !imax || !jmax || !kmax || !nmax || !vc || !vc_comm || !datatype || !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // MPI_Datatype
  MPI_Datatype dtype = cpm_ParaManager::GetMPI_Datatype( *datatype );
  if( dtype == MPI_DATATYPE_NULL )
  {
    *ierr = CPM_ERROR_MPI_INVALID_DATATYPE;
    return;
  }

  // BndCommS4DEx
   *ierr = paraMngr->BndCommS4DEx( dtype, array, *nmax, *imax, *jmax, *kmax, *vc, *vc_comm, *procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
/** 袖通信(Vector3DEx版)のFortranインターフェイス
 *  - (3,imax,jmax,kmax)の形式の配列の袖通信を行う
 *  - BndCommV3DのFortranインターフェイス関数
 *  @param[inout] array     袖通信をする配列の先頭ポインタ
 *  @param[in]    imax      配列サイズ(I方向)
 *  @param[in]    jmax      配列サイズ(J方向)
 *  @param[in]    kmax      配列サイズ(K方向)
 *  @param[in]    vc        仮想セル数
 *  @param[in]    vc_comm   通信する仮想セル数
 *  @param[in]    datatype  袖通信データのデータタイプ(cpm_fparam.fiを参照)
 *  @param[in]    procGrpNo プロセスグループ番号
 *  @param[out]   ierr      終了コード(CPM_SUCCESS=正常終了)
 */
CPM_EXTERN
void
cpm_BndCommV3DEx_( void *array, int *imax, int *jmax, int *kmax, int *vc, int *vc_comm, int *datatype
                 , int *procGrpNo, int *ierr )
{
#ifdef _USE_S4D_
  int nmax = 3;
  cpm_BndCommS4DEx_( array, &nmax, imax, jmax, kmax, vc, vc_comm, datatype, procGrpNo, ierr );
#else
  if( !array || !imax || !jmax || !kmax || !vc || !vc_comm || !datatype || !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // MPI_Datatype
  MPI_Datatype dtype = cpm_ParaManager::GetMPI_Datatype( *datatype );
  if( dtype == MPI_DATATYPE_NULL )
  {
    *ierr = CPM_ERROR_MPI_INVALID_DATATYPE;
    return;
  }

  // BndCommV3DEx
  *ierr = paraMngr->BndCommV3DEx( dtype, array, *imax, *jmax, *kmax, *vc, *vc_comm, *procGrpNo );
#endif
}

////////////////////////////////////////////////////////////////////////////////
/** 非同期版袖通信(Scalar4DEx版)のFortranインターフェイス
 *  - (nmax,imax,jmax,kmax)の形式の配列の袖通信を行う
 *  - BndCommS4DEx_nowaitのFortranインターフェイス関数
 *  @param[inout] array     袖通信をする配列の先頭ポインタ
 *  @param[in]    nmax      配列サイズ(成分数)
 *  @param[in]    imax      配列サイズ(I方向)
 *  @param[in]    jmax      配列サイズ(J方向)
 *  @param[in]    kmax      配列サイズ(K方向)
 *  @param[in]    vc        仮想セル数
 *  @param[in]    vc_comm   通信する仮想セル数
 *  @param[in]    datatype  袖通信データのデータタイプ(cpm_fparam.fiを参照)
 *  @param[out]   reqlist   リクエスト番号のリスト(サイズ12)
 *  @param[in]    procGrpNo プロセスグループ番号
 *  @param[out]   ierr      終了コード(CPM_SUCCESS=正常終了)
 */
CPM_EXTERN
void
cpm_BndCommS4DEx_nowait_( void *array, int *nmax, int *imax, int *jmax, int *kmax
                        , int *vc, int *vc_comm, int *datatype, int *reqlist, int *procGrpNo, int *ierr )
{
  if( !array || !imax || !jmax || !kmax || !nmax || !vc || !vc_comm || !datatype ||
      !reqlist || !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // cpm_BndCommS4DEx_nowait
   *ierr = paraMngr->cpm_BndCommS4DEx_nowait( array, *nmax, *imax, *jmax, *kmax, *vc, *vc_comm
                                            , *datatype, reqlist, *procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
/** 非同期版袖通信(Vector3DEx版)のFortranインターフェイス
 *  - (3,imax,jmax,kmax)の形式の配列の袖通信を行う
 *  - BndCommV3D_nowaitのFortranインターフェイス関数
 *  @param[inout] array     袖通信をする配列の先頭ポインタ
 *  @param[in]    imax      配列サイズ(I方向)
 *  @param[in]    jmax      配列サイズ(J方向)
 *  @param[in]    kmax      配列サイズ(K方向)
 *  @param[in]    vc        仮想セル数
 *  @param[in]    vc_comm   通信する仮想セル数
 *  @param[in]    datatype  袖通信データのデータタイプ(cpm_fparam.fiを参照)
 *  @param[out]   reqlist   リクエスト番号のリスト(サイズ12)
 *  @param[in]    procGrpNo プロセスグループ番号
 *  @param[out]   ierr      終了コード(CPM_SUCCESS=正常終了)
 */
CPM_EXTERN
void
cpm_BndCommV3DEx_nowait_( void *array, int *imax, int *jmax, int *kmax
                        , int *vc, int *vc_comm, int *datatype, int *reqlist, int *procGrpNo, int *ierr )
{
#ifdef _USE_S4D_
  int nmax = 3;
  cpm_BndCommS4DEx_nowait_( array, &nmax, imax, jmax, kmax, vc, vc_comm
                          , datatype, reqlist, procGrpNo, ierr );
#else
  if( !array || !imax || !jmax || !kmax || !vc || !vc_comm || !datatype ||
      !reqlist || !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // cpm_BndCommV3DEx_nowait
   *ierr = paraMngr->cpm_BndCommV3DEx_nowait( array, *imax, *jmax, *kmax, *vc, *vc_comm
                                            , *datatype, reqlist, *procGrpNo );
#endif
}

////////////////////////////////////////////////////////////////////////////////
/** 非同期版袖通信のwait、展開(Scalar4DEx版)のFortranインターフェイス
 *  - (nmax,imax,jmax,kmax)の形式の配列の非同期版袖通信のwaitと展開を行う
 *  - wait_BndCommS4DExのFortranインターフェイス関数
 *  @param[inout] array     袖通信をする配列の先頭ポインタ
 *  @param[in]    nmax      配列サイズ(成分数)
 *  @param[in]    imax      配列サイズ(I方向)
 *  @param[in]    jmax      配列サイズ(J方向)
 *  @param[in]    kmax      配列サイズ(K方向)
 *  @param[in]    vc        仮想セル数
 *  @param[in]    vc_comm   通信する仮想セル数
 *  @param[in]    datatype  袖通信データのデータタイプ(cpm_fparam.fiを参照)
 *  @param[out]   reqlist   リクエスト番号のリスト(サイズ12)
 *  @param[in]    procGrpNo プロセスグループ番号
 *  @param[out]   ierr      終了コード(CPM_SUCCESS=正常終了)
 */
CPM_EXTERN
void
cpm_wait_BndCommS4DEx_( void *array, int *nmax, int *imax, int *jmax, int *kmax
                      , int *vc, int *vc_comm, int *datatype, int *reqlist, int *procGrpNo, int *ierr )
{
  if( !array || !imax || !jmax || !kmax || !nmax || !vc || !vc_comm || !datatype ||
      !reqlist || !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // cpm_wait_BndCommS4DEx
   *ierr = paraMngr->cpm_wait_BndCommS4DEx( array, *nmax, *imax, *jmax, *kmax, *vc, *vc_comm
                                          , *datatype, reqlist, *procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
/** 非同期版袖通信のwait、展開(Vector3DEx版)のFortranインターフェイス
 *  - (3,imax,jmax,kmax)の形式の配列の非同期版袖通信のwaitと展開を行う
 *  - wait_BndCommV3DExのFortranインターフェイス関数
 *  @param[inout] array     袖通信をする配列の先頭ポインタ
 *  @param[in]    imax      配列サイズ(I方向)
 *  @param[in]    jmax      配列サイズ(J方向)
 *  @param[in]    kmax      配列サイズ(K方向)
 *  @param[in]    vc        仮想セル数
 *  @param[in]    vc_comm   通信する仮想セル数
 *  @param[in]    datatype  袖通信データのデータタイプ(cpm_fparam.fiを参照)
 *  @param[out]   reqlist   リクエスト番号のリスト(サイズ12)
 *  @param[in]    procGrpNo プロセスグループ番号
 *  @param[out]   ierr      終了コード(CPM_SUCCESS=正常終了)
 */
CPM_EXTERN
void
cpm_wait_BndCommV3DEx_( void *array, int *imax, int *jmax, int *kmax
                      , int *vc, int *vc_comm, int *datatype, int *reqlist, int *procGrpNo, int *ierr )
{
#ifdef _USE_S4D_
  int nmax = 3;
  cpm_wait_BndCommS4DEx_( array, &nmax, imax, jmax, kmax, vc, vc_comm, datatype, reqlist, procGrpNo, ierr );
#else
  if( !array || !imax || !jmax || !kmax || !vc || !vc_comm || !datatype ||
      !reqlist || !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // cpm_wait_BndCommV3DEx
   *ierr = paraMngr->cpm_wait_BndCommV3DEx( array, *imax, *jmax, *kmax, *vc, *vc_comm
                                          , *datatype, reqlist, *procGrpNo );
#endif
}

////////////////////////////////////////////////////////////////////////////////
/** 周期境界袖通信(Scalar4D版)のFortranインターフェイス
 *  - (imax,jmax,kmax,nmax)の形式の配列の周期境界方向の袖通信を行う
 *  - PeriodicCommS4DのFortranインターフェイス関数
 *  @param[inout] array     袖通信をする配列の先頭ポインタ
 *  @param[in]    imax      配列サイズ(I方向)
 *  @param[in]    jmax      配列サイズ(J方向)
 *  @param[in]    kmax      配列サイズ(K方向)
 *  @param[in]    nmax      配列サイズ(成分数)
 *  @param[in]    vc        仮想セル数
 *  @param[in]    vc_comm   通信する仮想セル数
 *  @param[in]    dir       通信する軸方向(X_DIR or Y_DIR or Z_DIR)
 *  @param[in]    pm        通信する正負方向(PLUS2MINUS or MINUS2PLUS or BOTH)
 *  @param[in]    datatype  袖通信データのデータタイプ(cpm_fparam.fiを参照)
 *  @param[in]    procGrpNo プロセスグループ番号
 *  @param[out]   ierr      終了コード(CPM_SUCCESS=正常終了)
 */
CPM_EXTERN
void
cpm_PeriodicCommS4D_( void *array, int *imax, int *jmax, int *kmax, int *nmax, int *vc, int *vc_comm
                    , int *dir, int *pm, int *datatype, int *procGrpNo, int *ierr )
{
  if( !array || !imax || !jmax || !kmax || !nmax || !vc || !vc_comm || !dir || !pm ||
      !datatype || !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // MPI_Datatype
  MPI_Datatype dtype = cpm_ParaManager::GetMPI_Datatype( *datatype );
  if( dtype == MPI_DATATYPE_NULL )
  {
    *ierr = CPM_ERROR_MPI_INVALID_DATATYPE;
    return;
  }

  // cpm_DirFlag
  if( *dir < 0 || *dir > 2 )
  {
    *ierr = CPM_ERROR_PERIODIC_INVALID_DIR;
    return;
  }

  // cpm_PMFlag
  if( *pm < 0 || *pm > 2 )
  {
    *ierr = CPM_ERROR_PERIODIC_INVALID_PM;
    return;
  }

  // PeriodicCommS4D
  *ierr = paraMngr->PeriodicCommS4D( dtype, array, *imax, *jmax, *kmax, *nmax, *vc, *vc_comm
                                   , (cpm_DirFlag)*dir, (cpm_PMFlag)*pm, *procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
/** 周期境界袖通信(Scalar3D版)のFortranインターフェイス
 *  - (imax,jmax,kmax)の形式の配列の周期境界方向の袖通信を行う
 *  - PeriodicCommS3DのFortranインターフェイス関数
 *  @param[inout] array     袖通信をする配列の先頭ポインタ
 *  @param[in]    imax      配列サイズ(I方向)
 *  @param[in]    jmax      配列サイズ(J方向)
 *  @param[in]    kmax      配列サイズ(K方向)
 *  @param[in]    vc        仮想セル数
 *  @param[in]    vc_comm   通信する仮想セル数
 *  @param[in]    dir       通信する軸方向(X_DIR or Y_DIR or Z_DIR)
 *  @param[in]    pm        通信する正負方向(PLUS2MINUS or MINUS2PLUS or BOTH)
 *  @param[in]    datatype  袖通信データのデータタイプ(cpm_fparam.fiを参照)
 *  @param[in]    procGrpNo プロセスグループ番号
 *  @param[out]   ierr      終了コード(CPM_SUCCESS=正常終了)
 */
CPM_EXTERN
void
cpm_PeriodicCommS3D_( void *array, int *imax, int *jmax, int *kmax, int *vc, int *vc_comm
                    , int *dir, int *pm, int *datatype, int *procGrpNo, int *ierr )
{
#ifdef _USE_S4D_
  int nmax = 1;
  cpm_PeriodicCommS4D_( array, imax, jmax, kmax, &nmax, vc, vc_comm, dir, pm
                      , datatype, procGrpNo, ierr );
#else
  if( !array || !imax || !jmax || !kmax || !vc || !vc_comm || !dir || !pm ||
      !datatype || !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // MPI_Datatype
  MPI_Datatype dtype = cpm_ParaManager::GetMPI_Datatype( *datatype );
  if( dtype == MPI_DATATYPE_NULL )
  {
    *ierr = CPM_ERROR_MPI_INVALID_DATATYPE;
    return;
  }

  // cpm_DirFlag
  if( *dir < 0 || *dir > 2 )
  {
    *ierr = CPM_ERROR_PERIODIC_INVALID_DIR;
    return;
  }

  // cpm_PMFlag
  if( *pm < 0 || *pm > 2 )
  {
    *ierr = CPM_ERROR_PERIODIC_INVALID_PM;
    return;
  }

  // PeriodicCommS3D
  *ierr = paraMngr->PeriodicCommS3D( dtype, array, *imax, *jmax, *kmax, *vc, *vc_comm
                                   , (cpm_DirFlag)*dir, (cpm_PMFlag)*pm, *procGrpNo );
#endif
}

////////////////////////////////////////////////////////////////////////////////
/** 周期境界袖通信(Vector3D版)のFortranインターフェイス
 *  - (imax,jmax,kmax,3)の形式の配列の周期境界方向の袖通信を行う
 *  - PeriodicCommV3DのFortranインターフェイス関数
 *  @param[inout] array     袖通信をする配列の先頭ポインタ
 *  @param[in]    imax      配列サイズ(I方向)
 *  @param[in]    jmax      配列サイズ(J方向)
 *  @param[in]    kmax      配列サイズ(K方向)
 *  @param[in]    vc        仮想セル数
 *  @param[in]    vc_comm   通信する仮想セル数
 *  @param[in]    dir       通信する軸方向(X_DIR or Y_DIR or Z_DIR)
 *  @param[in]    pm        通信する正負方向(PLUS2MINUS or MINUS2PLUS or BOTH)
 *  @param[in]    datatype  袖通信データのデータタイプ(cpm_fparam.fiを参照)
 *  @param[in]    procGrpNo プロセスグループ番号
 *  @param[out]   ierr      終了コード(CPM_SUCCESS=正常終了)
 */
CPM_EXTERN
void
cpm_PeriodicCommV3D_( void *array, int *imax, int *jmax, int *kmax, int *vc, int *vc_comm
                    , int *dir, int *pm, int *datatype, int *procGrpNo, int *ierr )
{
#ifdef _USE_S4D_
  int nmax = 3;
  cpm_PeriodicCommS4D_( array, imax, jmax, kmax, &nmax, vc, vc_comm, dir, pm
                      , datatype, procGrpNo, ierr );
#else
  if( !array || !imax || !jmax || !kmax || !vc || !vc_comm || !dir || !pm ||
      !datatype || !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // MPI_Datatype
  MPI_Datatype dtype = cpm_ParaManager::GetMPI_Datatype( *datatype );
  if( dtype == MPI_DATATYPE_NULL )
  {
    *ierr = CPM_ERROR_MPI_INVALID_DATATYPE;
    return;
  }

  // cpm_DirFlag
  if( *dir < 0 || *dir > 2 )
  {
    *ierr = CPM_ERROR_PERIODIC_INVALID_DIR;
    return;
  }

  // cpm_PMFlag
  if( *pm < 0 || *pm > 2 )
  {
    *ierr = CPM_ERROR_PERIODIC_INVALID_PM;
    return;
  }

  // PeriodicCommV3D
  *ierr = paraMngr->PeriodicCommV3D( dtype, array, *imax, *jmax, *kmax, *vc, *vc_comm
                                   , (cpm_DirFlag)*dir, (cpm_PMFlag)*pm, *procGrpNo );
#endif
}

////////////////////////////////////////////////////////////////////////////////
/** 周期境界袖通信(Scalar4DEx版)のFortranインターフェイス
 *  - (nmax,imax,jmax,kmax)の形式の配列の周期境界方向の袖通信を行う
 *  - PeriodicCommS4DExのFortranインターフェイス関数
 *  @param[inout] array     袖通信をする配列の先頭ポインタ
 *  @param[in]    nmax      配列サイズ(成分数)
 *  @param[in]    imax      配列サイズ(I方向)
 *  @param[in]    jmax      配列サイズ(J方向)
 *  @param[in]    kmax      配列サイズ(K方向)
 *  @param[in]    vc        仮想セル数
 *  @param[in]    vc_comm   通信する仮想セル数
 *  @param[in]    dir       通信する軸方向(X_DIR or Y_DIR or Z_DIR)
 *  @param[in]    pm        通信する正負方向(PLUS2MINUS or MINUS2PLUS or BOTH)
 *  @param[in]    datatype  袖通信データのデータタイプ(cpm_fparam.fiを参照)
 *  @param[in]    procGrpNo プロセスグループ番号
 *  @param[out]   ierr      終了コード(CPM_SUCCESS=正常終了)
 */
CPM_EXTERN
void
cpm_PeriodicCommS4DEx_( void *array, int *nmax, int *imax, int *jmax, int *kmax, int *vc, int *vc_comm
                      , int *dir, int *pm, int *datatype, int *procGrpNo, int *ierr )
{
  if( !array || !imax || !jmax || !kmax || !nmax || !vc || !vc_comm || !dir || !pm ||
      !datatype || !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // MPI_Datatype
  MPI_Datatype dtype = cpm_ParaManager::GetMPI_Datatype( *datatype );
  if( dtype == MPI_DATATYPE_NULL )
  {
    *ierr = CPM_ERROR_MPI_INVALID_DATATYPE;
    return;
  }

  // cpm_DirFlag
  if( *dir < 0 || *dir > 2 )
  {
    *ierr = CPM_ERROR_PERIODIC_INVALID_DIR;
    return;
  }

  // cpm_PMFlag
  if( *pm < 0 || *pm > 2 )
  {
    *ierr = CPM_ERROR_PERIODIC_INVALID_PM;
    return;
  }

  // PeriodicCommS4DEx
  *ierr = paraMngr->PeriodicCommS4DEx( dtype, array, *nmax, *imax, *jmax, *kmax, *vc, *vc_comm
                                     , (cpm_DirFlag)*dir, (cpm_PMFlag)*pm, *procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
/** 周期境界袖通信(Vector3DEx版)のFortranインターフェイス
 *  - (3,imax,jmax,kmax)の形式の配列の周期境界方向の袖通信を行う
 *  - PeriodicCommV3DExのFortranインターフェイス関数
 *  @param[inout] array     袖通信をする配列の先頭ポインタ
 *  @param[in]    imax      配列サイズ(I方向)
 *  @param[in]    jmax      配列サイズ(J方向)
 *  @param[in]    kmax      配列サイズ(K方向)
 *  @param[in]    vc        仮想セル数
 *  @param[in]    vc_comm   通信する仮想セル数
 *  @param[in]    dir       通信する軸方向(X_DIR or Y_DIR or Z_DIR)
 *  @param[in]    pm        通信する正負方向(PLUS2MINUS or MINUS2PLUS or BOTH)
 *  @param[in]    datatype  袖通信データのデータタイプ(cpm_fparam.fiを参照)
 *  @param[in]    procGrpNo プロセスグループ番号
 *  @param[out]   ierr      終了コード(CPM_SUCCESS=正常終了)
 */
CPM_EXTERN
void
cpm_PeriodicCommV3DEx_( void *array, int *imax, int *jmax, int *kmax, int *vc, int *vc_comm
                      , int *dir, int *pm, int *datatype, int *procGrpNo, int *ierr )
{
#ifdef _USE_S4D_
  int nmax = 3;
  cpm_PeriodicCommS4DEx_( array, &nmax, imax, jmax, kmax, vc, vc_comm, dir, pm
                        , datatype, procGrpNo, ierr );
#else
  if( !array || !imax || !jmax || !kmax || !vc || !vc_comm || !dir || !pm ||
      !datatype || !procGrpNo || !ierr )
  {
    if( ierr ) *ierr = CPM_ERROR_INVALID_PTR;
    return;
  }

  // インスタンス取得
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance();
  if( !paraMngr )
  {
    *ierr = CPM_ERROR_PM_INSTANCE;
    return;
  }

  // MPI_Datatype
  MPI_Datatype dtype = cpm_ParaManager::GetMPI_Datatype( *datatype );
  if( dtype == MPI_DATATYPE_NULL )
  {
    *ierr = CPM_ERROR_MPI_INVALID_DATATYPE;
    return;
  }

  // cpm_DirFlag
  if( *dir < 0 || *dir > 2 )
  {
    *ierr = CPM_ERROR_PERIODIC_INVALID_DIR;
    return;
  }

  // cpm_PMFlag
  if( *pm < 0 || *pm > 2 )
  {
    *ierr = CPM_ERROR_PERIODIC_INVALID_PM;
    return;
  }

  // PeriodicCommV3DEx
  *ierr = paraMngr->PeriodicCommV3DEx( dtype, array, *imax, *jmax, *kmax, *vc, *vc_comm
                                     , (cpm_DirFlag)*dir, (cpm_PMFlag)*pm, *procGrpNo );
#endif
}

////////////////////////////////////////////////////////////////////////////////
// cpm_Wait
cpm_ErrorCode
cpm_ParaManager::cpm_Wait( int reqNo )
{
  // MPI_Request
  MPI_Request *req = m_reqList.Get(reqNo);
  if( !req )
  {
    return CPM_ERROR_INVALID_OBJKEY;
  }

  // MPI_Wait
  MPI_Status status;
  if( MPI_Wait( req, &status ) != MPI_SUCCESS )
  {
    return CPM_ERROR_MPI_WAIT;
  }

  // 削除
  return m_reqList.Delete(reqNo);
}

////////////////////////////////////////////////////////////////////////////////
// cpm_Waitall
cpm_ErrorCode
cpm_ParaManager::cpm_Waitall( int count, int reqNoList[] )
{
  // MPI_Request
  MPI_Status* stat = new MPI_Status [count];
  MPI_Request* req = new MPI_Request[count];
  int cnt = 0;
  for( int i=0;i<count;i++ )
  {
    MPI_Request *r = m_reqList.Get( reqNoList[i] );
    if( r )
    {
      r[cnt++] = *req;
    }
  }
  if( cnt == 0 )
  {
    delete [] stat;
    delete [] req;
    return CPM_ERROR_INVALID_OBJKEY;
  }

  // MPI_Waitall
  if( MPI_Waitall( cnt, req, stat ) != MPI_SUCCESS )
  {
    delete [] stat;
    delete [] req;
    return CPM_ERROR_MPI_WAITALL;
  }

  // 削除
  for( int i=0;i<count;i++ )
  {
    m_reqList.Delete( reqNoList[i] );
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// cpm_Isend
cpm_ErrorCode
cpm_ParaManager::cpm_Isend( void *buf, int count, int datatype, int dest, int *reqNo, int procGrpNo )
{
  // MPI_Datatype
  MPI_Datatype dtype = cpm_ParaManager::GetMPI_Datatype( datatype );
  if( dtype == MPI_DATATYPE_NULL )
  {
    return CPM_ERROR_MPI_INVALID_DATATYPE;
  }

  // MPI_Request
  MPI_Request *req = m_reqList.Create();
  if( !req )
  {
    return CPM_ERROR_INVALID_PTR;
  }

  // Isend
  cpm_ErrorCode ret = Isend( dtype, buf, count, dest, req, procGrpNo );
  if( ret != MPI_SUCCESS )
  {
    delete req;
    return ret;
  }

  // MPI_Requestを登録
  if( (*reqNo = m_reqList.Add(req) ) < 0 )
  {
    delete req;
    return CPM_ERROR_REGIST_OBJKEY;
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// cpm_Irecv
cpm_ErrorCode
cpm_ParaManager::cpm_Irecv( void *buf, int count, int datatype, int source, int *reqNo, int procGrpNo )
{
  // MPI_Datatype
  MPI_Datatype dtype = cpm_ParaManager::GetMPI_Datatype( datatype );
  if( dtype == MPI_DATATYPE_NULL )
  {
    return CPM_ERROR_MPI_INVALID_DATATYPE;
  }

  // Irecv
  MPI_Request req;
  cpm_ErrorCode ret = Irecv( dtype, buf, count, source, &req, procGrpNo );
  if( ret != MPI_SUCCESS )
  {
    return ret;
  }

  // MPI_Requestを登録
  MPI_Request *r = m_reqList.Create();
  *r = req;
  if( (*reqNo = m_reqList.Add(r) ) < 0 )
  {
    delete r;
    return CPM_ERROR_REGIST_OBJKEY;
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// cpm_BndCommS3D_nowait
cpm_ErrorCode
cpm_ParaManager::cpm_BndCommS3D_nowait( void *array, int imax, int jmax, int kmax
                                      , int vc, int vc_comm, int datatype, int *reqNo, int procGrpNo )
{
#ifdef _USE_S4D_
  return cpm_BndCommS4D_nowait( array, imax, jmax, kmax, 1, vc, vc_comm, datatype, reqNo, procGrpNo );
#else
  // MPI_Datatype
  MPI_Datatype dtype = cpm_ParaManager::GetMPI_Datatype( datatype );
  if( dtype == MPI_DATATYPE_NULL )
  {
    return CPM_ERROR_MPI_INVALID_DATATYPE;
  }

  // MPI_Request
  MPI_Request req[12];

  // BndCommS3D_nowait
  cpm_ErrorCode ret = BndCommS3D_nowait( dtype, array, imax, jmax, kmax
                                       , vc, vc_comm, req, procGrpNo );
  if( ret != MPI_SUCCESS )
  {
    return ret;
  }

  // MPI_Requestを登録
  for( int i=0;i<12;i++ )
  {
    MPI_Request *r = m_reqList.Create();
    *r = req[i];
    if( (reqNo[i] = m_reqList.Add(r) ) < 0 )
    {
      delete r;
      return CPM_ERROR_REGIST_OBJKEY;
    }
  }
  return CPM_SUCCESS;
#endif
}

////////////////////////////////////////////////////////////////////////////////
// cpm_BndCommV3D_nowait
cpm_ErrorCode
cpm_ParaManager::cpm_BndCommV3D_nowait( void *array, int imax, int jmax, int kmax
                                      , int vc, int vc_comm, int datatype, int *reqNo, int procGrpNo )
{
#ifdef _USE_S4D_
  return cpm_BndCommS4D_nowait( array, imax, jmax, kmax, 3, vc, vc_comm, datatype, reqNo, procGrpNo );
#else
  // MPI_Datatype
  MPI_Datatype dtype = cpm_ParaManager::GetMPI_Datatype( datatype );
  if( dtype == MPI_DATATYPE_NULL )
  {
    return CPM_ERROR_MPI_INVALID_DATATYPE;
  }

  // MPI_Request
  MPI_Request req[12];

  // BndCommV3D_nowait
  cpm_ErrorCode ret = BndCommV3D_nowait( dtype, array, imax, jmax, kmax
                                       , vc, vc_comm, req, procGrpNo );
  if( ret != MPI_SUCCESS )
  {
    return ret;
  }

  // MPI_Requestを登録
  for( int i=0;i<12;i++ )
  {
    MPI_Request *r = m_reqList.Create();
    *r = req[i];
    if( (reqNo[i] = m_reqList.Add(r) ) < 0 )
    {
      delete r;
      return CPM_ERROR_REGIST_OBJKEY;
    }
  }
  return CPM_SUCCESS;
#endif
}

////////////////////////////////////////////////////////////////////////////////
// cpm_BndCommS4D_nowait
cpm_ErrorCode
cpm_ParaManager::cpm_BndCommS4D_nowait( void *array, int imax, int jmax, int kmax, int nmax
                                      , int vc, int vc_comm, int datatype, int *reqNo, int procGrpNo )
{
  // MPI_Datatype
  MPI_Datatype dtype = cpm_ParaManager::GetMPI_Datatype( datatype );
  if( dtype == MPI_DATATYPE_NULL )
  {
    return CPM_ERROR_MPI_INVALID_DATATYPE;
  }

  // MPI_Request
  MPI_Request req[12];

  // BndCommS4D_nowait
  cpm_ErrorCode ret = BndCommS4D_nowait( dtype, array, imax, jmax, kmax, nmax
                                       , vc, vc_comm, req, procGrpNo );
  if( ret != CPM_SUCCESS )
  {
    return ret;
  }

  // MPI_Requestを登録
  for( int i=0;i<12;i++ )
  {
    MPI_Request *r = m_reqList.Create();
    *r = req[i];
    if( (reqNo[i] = m_reqList.Add(r) ) < 0 )
    {
      delete r;
      return CPM_ERROR_REGIST_OBJKEY;
    }
  }
  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// cpm_wait_BndCommS3D
cpm_ErrorCode
cpm_ParaManager::cpm_wait_BndCommS3D( void *array, int imax, int jmax, int kmax
                                    , int vc, int vc_comm, int datatype, int *reqNo, int procGrpNo )
{
#ifdef _USE_S4D_
  return cpm_wait_BndCommS4D( array, imax, jmax, kmax, 1, vc, vc_comm, datatype, reqNo, procGrpNo );
#else
  // MPI_Datatype
  MPI_Datatype dtype = cpm_ParaManager::GetMPI_Datatype( datatype );
  if( dtype == MPI_DATATYPE_NULL )
  {
    return CPM_ERROR_MPI_INVALID_DATATYPE;
  }

  // MPI_Request
  MPI_Request req[12];
  for( int i=0;i<12;i++ )
  {
    MPI_Request *r = m_reqList.Get(reqNo[i]);
    if( !r )
    {
      return CPM_ERROR_INVALID_OBJKEY;
    }
    req[i] = *r;
  }

  // wait_BndCommS3D
  cpm_ErrorCode ret = wait_BndCommS3D( dtype, array, imax, jmax, kmax
                                     , vc, vc_comm, req, procGrpNo );
  if( ret != CPM_SUCCESS )
  {
    return ret;
  }

  // MPI_Requestを削除
  for( int i=0;i<12;i++ )
  {
    m_reqList.Delete(reqNo[i]);
  }

  return CPM_SUCCESS;
#endif
}

////////////////////////////////////////////////////////////////////////////////
// cpm_wait_BndCommV3D
cpm_ErrorCode
cpm_ParaManager::cpm_wait_BndCommV3D( void *array, int imax, int jmax, int kmax
                                    , int vc, int vc_comm, int datatype, int *reqNo, int procGrpNo )
{
#ifdef _USE_S4D_
  return cpm_wait_BndCommS4D( array, imax, jmax, kmax, 3, vc, vc_comm, datatype, reqNo, procGrpNo );
#else
  // MPI_Datatype
  MPI_Datatype dtype = cpm_ParaManager::GetMPI_Datatype( datatype );
  if( dtype == MPI_DATATYPE_NULL )
  {
    return CPM_ERROR_MPI_INVALID_DATATYPE;
  }

  // MPI_Request
  MPI_Request req[12];
  for( int i=0;i<12;i++ )
  {
    MPI_Request *r = m_reqList.Get(reqNo[i]);
    if( !r )
    {
      return CPM_ERROR_INVALID_OBJKEY;
    }
    req[i] = *r;
  }

  // wait_BndCommV3D
  cpm_ErrorCode ret = wait_BndCommV3D( dtype, array, imax, jmax, kmax
                                     , vc, vc_comm, req, procGrpNo );
  if( ret != CPM_SUCCESS )
  {
    return ret;
  }

  // MPI_Requestを削除
  for( int i=0;i<12;i++ )
  {
    m_reqList.Delete(reqNo[i]);
  }

  return CPM_SUCCESS;
#endif
}

////////////////////////////////////////////////////////////////////////////////
// cpm_wait_BndCommS4D
cpm_ErrorCode
cpm_ParaManager::cpm_wait_BndCommS4D( void *array, int imax, int jmax, int kmax, int nmax
                                    , int vc, int vc_comm, int datatype, int *reqNo, int procGrpNo )
{
  // MPI_Datatype
  MPI_Datatype dtype = cpm_ParaManager::GetMPI_Datatype( datatype );
  if( dtype == MPI_DATATYPE_NULL )
  {
    return CPM_ERROR_MPI_INVALID_DATATYPE;
  }

  // MPI_Request
  MPI_Request req[12];
  for( int i=0;i<12;i++ )
  {
    MPI_Request *r = m_reqList.Get(reqNo[i]);
    if( !r )
    {
      return CPM_ERROR_INVALID_OBJKEY;
    }
    req[i] = *r;
  }

  // wait_BndCommS4D
  cpm_ErrorCode ret = wait_BndCommS4D( dtype, array, imax, jmax, kmax, nmax
                                     , vc, vc_comm, req, procGrpNo );
  if( ret != CPM_SUCCESS )
  {
    return ret;
  }

  // MPI_Requestを削除
  for( int i=0;i<12;i++ )
  {
    m_reqList.Delete(reqNo[i]);
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// cpm_BndCommV3DEx_nowait
cpm_ErrorCode
cpm_ParaManager::cpm_BndCommV3DEx_nowait( void *array, int imax, int jmax, int kmax
                                        , int vc, int vc_comm, int datatype, int *reqNo, int procGrpNo )
{
#ifdef _USE_S4D_
  return cpm_BndCommS4DEx_nowait( array, 3, imax, jmax, kmax, vc, vc_comm, datatype, reqNo, procGrpNo );
#else
  // MPI_Datatype
  MPI_Datatype dtype = cpm_ParaManager::GetMPI_Datatype( datatype );
  if( dtype == MPI_DATATYPE_NULL )
  {
    return CPM_ERROR_MPI_INVALID_DATATYPE;
  }

  // MPI_Request
  MPI_Request req[12];

  // BndCommV3DEx_nowait
  cpm_ErrorCode ret = BndCommV3DEx_nowait( dtype, array, imax, jmax, kmax
                                         , vc, vc_comm, req, procGrpNo );
  if( ret != MPI_SUCCESS )
  {
    return ret;
  }

  // MPI_Requestを登録
  for( int i=0;i<12;i++ )
  {
    MPI_Request *r = m_reqList.Create();
    *r = req[i];
    if( (reqNo[i] = m_reqList.Add(r) ) < 0 )
    {
      delete r;
      return CPM_ERROR_REGIST_OBJKEY;
    }
  }
  return CPM_SUCCESS;
#endif
}

////////////////////////////////////////////////////////////////////////////////
// cpm_BndCommS4DEx_nowait
cpm_ErrorCode
cpm_ParaManager::cpm_BndCommS4DEx_nowait( void *array, int nmax, int imax, int jmax, int kmax
                                        , int vc, int vc_comm, int datatype, int *reqNo, int procGrpNo )
{
  // MPI_Datatype
  MPI_Datatype dtype = cpm_ParaManager::GetMPI_Datatype( datatype );
  if( dtype == MPI_DATATYPE_NULL )
  {
    return CPM_ERROR_MPI_INVALID_DATATYPE;
  }

  // MPI_Request
  MPI_Request req[12];

  // BndCommS4DEx_nowait
  cpm_ErrorCode ret = BndCommS4DEx_nowait( dtype, array, nmax, imax, jmax, kmax
                                         , vc, vc_comm, req, procGrpNo );
  if( ret != CPM_SUCCESS )
  {
    return ret;
  }

  // MPI_Requestを登録
  for( int i=0;i<12;i++ )
  {
    MPI_Request *r = m_reqList.Create();
    *r = req[i];
    if( (reqNo[i] = m_reqList.Add(r) ) < 0 )
    {
      delete r;
      return CPM_ERROR_REGIST_OBJKEY;
    }
  }
  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// cpm_wait_BndCommV3DEx
cpm_ErrorCode
cpm_ParaManager::cpm_wait_BndCommV3DEx( void *array, int imax, int jmax, int kmax
                                      , int vc, int vc_comm, int datatype, int *reqNo, int procGrpNo )
{
#ifdef _USE_S4D_
  return cpm_wait_BndCommS4DEx( array, 3, imax, jmax, kmax, vc, vc_comm, datatype, reqNo, procGrpNo );
#else
  // MPI_Datatype
  MPI_Datatype dtype = cpm_ParaManager::GetMPI_Datatype( datatype );
  if( dtype == MPI_DATATYPE_NULL )
  {
    return CPM_ERROR_MPI_INVALID_DATATYPE;
  }

  // MPI_Request
  MPI_Request req[12];
  for( int i=0;i<12;i++ )
  {
    MPI_Request *r = m_reqList.Get(reqNo[i]);
    if( !r )
    {
      return CPM_ERROR_INVALID_OBJKEY;
    }
    req[i] = *r;
  }

  // wait_BndCommV3DEx
  cpm_ErrorCode ret = wait_BndCommV3DEx( dtype, array, imax, jmax, kmax
                                       , vc, vc_comm, req, procGrpNo );
  if( ret != CPM_SUCCESS )
  {
    return ret;
  }

  // MPI_Requestを削除
  for( int i=0;i<12;i++ )
  {
    m_reqList.Delete(reqNo[i]);
  }

  return CPM_SUCCESS;
#endif
}

////////////////////////////////////////////////////////////////////////////////
// cpm_wait_BndCommS4DEx
cpm_ErrorCode
cpm_ParaManager::cpm_wait_BndCommS4DEx( void *array, int nmax, int imax, int jmax, int kmax
                                      , int vc, int vc_comm, int datatype, int *reqNo, int procGrpNo )
{
  // MPI_Datatype
  MPI_Datatype dtype = cpm_ParaManager::GetMPI_Datatype( datatype );
  if( dtype == MPI_DATATYPE_NULL )
  {
    return CPM_ERROR_MPI_INVALID_DATATYPE;
  }

  // MPI_Request
  MPI_Request req[12];
  for( int i=0;i<12;i++ )
  {
    MPI_Request *r = m_reqList.Get(reqNo[i]);
    if( !r )
    {
      return CPM_ERROR_INVALID_OBJKEY;
    }
    req[i] = *r;
  }

  // wait_BndCommS4DEx
  cpm_ErrorCode ret = wait_BndCommS4DEx( dtype, array, nmax, imax, jmax, kmax
                                       , vc, vc_comm, req, procGrpNo );
  if( ret != CPM_SUCCESS )
  {
    return ret;
  }

  // MPI_Requestを削除
  for( int i=0;i<12;i++ )
  {
    m_reqList.Delete(reqNo[i]);
  }

  return CPM_SUCCESS;
}





