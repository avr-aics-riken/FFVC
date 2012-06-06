/*
 * CPM - Cartesian Partition Manager
 *
 * Copyright (C) 2012 University of Tokyo.
 *
 */

/**
 * @file   cpm_ParaManager_BndCommEx.h
 * パラレルマネージャクラスのインラインヘッダーファイル
 * @author University of Tokyo
 * @date   2012/05/31
 */

#ifndef _CPM_PARAMANAGER_BNDCOMMEX_H_
#define _CPM_PARAMANAGER_BNDCOMMEX_H_

#define _IDXFX(_N,_I,_J,_K,_NN,_IS,_NJ,_NK,_VC) \
( size_t(_NN) \
* ( size_t(_K+_VC) * size_t(_VC) * size_t(_NJ+2*_VC) \
  + size_t(_J+_VC) * size_t(_VC) \
  + size_t(_I-(_IS)) \
  ) \
+ size_t(_N) \
)

#define _IDXFY(_N,_I,_J,_K,_NN,_NI,_JS,_NK,_VC) \
( size_t(_NN) \
* ( size_t(_K+_VC)   * size_t(_NI+2*_VC) * size_t(_VC) \
  + size_t(_J-(_JS)) * size_t(_NI+2*_VC) \
  + size_t(_I+_VC) \
  ) \
+ size_t(_N) \
)

#define _IDXFZ(_N,_I,_J,_K,_NN,_NI,_NJ,_KS,_VC) \
( size_t(_NN) \
* ( size_t(_K-(_KS)) * size_t(_NI+2*_VC) * size_t(_NJ+2*_VC) \
  + size_t(_J+_VC)   * size_t(_NI+2*_VC) \
  + size_t(_I+_VC) \
  ) \
+ size_t(_N) \
)

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Vector3DEx版)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::BndCommV3DEx( T *array, int imax, int jmax, int kmax, int vc, int vc_comm, int procGrpNo )
{
  return BndCommS4DEx( array, 3, imax, jmax, kmax, vc, vc_comm, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar4DEx版)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::BndCommS4DEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm, int procGrpNo )
{
  cpm_ErrorCode ret;

  if( !array )
  {
    return CPM_ERROR_INVALID_PTR;
  }

  // 通信バッファを取得
  S_BNDCOMM_BUFFER *bufInfo = GetBndCommBuffer(procGrpNo);
  if( !bufInfo )
  {
    return CPM_ERROR_BNDCOMM_BUFFER;
  }

  // 隣接ランクを取得
  const int *nID = GetNeighborRankID(procGrpNo);
  if( !nID )
  {
    return CPM_ERROR_GET_NEIGHBOR_RANK;
  }

  // 通信バッファサイズを計算
  size_t nwX = size_t(jmax+2*vc_comm) * size_t(kmax+2*vc_comm) * size_t(vc_comm) * size_t(nmax);
  size_t nwY = size_t(kmax+2*vc_comm) * size_t(imax+2*vc_comm) * size_t(vc_comm) * size_t(nmax);
  size_t nwZ = size_t(imax+2*vc_comm) * size_t(jmax+2*vc_comm) * size_t(vc_comm) * size_t(nmax);
  if( nwX > bufInfo->m_nwX || nwY > bufInfo->m_nwY || nwZ > bufInfo->m_nwZ )
  {
    return CPM_ERROR_BNDCOMM_BUFFERLENGTH;
  }

  // ポインタをセット
  T *sendmx = (T*)(bufInfo->m_bufX[0]);
  T *recvmx = (T*)(bufInfo->m_bufX[1]);
  T *sendpx = (T*)(bufInfo->m_bufX[2]);
  T *recvpx = (T*)(bufInfo->m_bufX[3]);
  T *sendmy = (T*)(bufInfo->m_bufY[0]);
  T *recvmy = (T*)(bufInfo->m_bufY[1]);
  T *sendpy = (T*)(bufInfo->m_bufY[2]);
  T *recvpy = (T*)(bufInfo->m_bufY[3]);
  T *sendmz = (T*)(bufInfo->m_bufZ[0]);
  T *recvmz = (T*)(bufInfo->m_bufZ[1]);
  T *sendpz = (T*)(bufInfo->m_bufZ[2]);
  T *recvpz = (T*)(bufInfo->m_bufZ[3]);

  MPI_Request req[12];

  //// X face ////
  int nIDmx = nID[X_MINUS];
  int nIDpx = nID[X_PLUS];

  // pack
  if( (ret = packXEx( array, nmax, imax, jmax, kmax, vc, vc_comm, sendmx, sendpx, nIDmx, nIDpx )) != CPM_SUCCESS ) return ret;

  // Isend/Irecv
  if( (ret = sendrecv( sendmx, recvmx, sendpx, recvpx, nwX, &req[0], nIDmx, nIDmx, nIDpx, nIDpx, procGrpNo )) != CPM_SUCCESS ) return ret;

  // wait
  if( (ret = Waitall( 4, &req[0] )) != CPM_SUCCESS ) return ret;

  // unpack
  if( (ret = unpackXEx( array, nmax, imax, jmax, kmax, vc, vc_comm, recvmx, recvpx, nIDmx, nIDpx )) != CPM_SUCCESS ) return ret;

  //// Y face ////
  int nIDmy = nID[Y_MINUS];
  int nIDpy = nID[Y_PLUS];

  // pack
  if( (ret = packYEx( array, nmax, imax, jmax, kmax, vc, vc_comm, sendmy, sendpy, nIDmy, nIDpy )) != CPM_SUCCESS ) return ret;

  // Isend/Irecv
  if( (ret = sendrecv( sendmy, recvmy, sendpy, recvpy, nwY, &req[4], nIDmy, nIDmy, nIDpy, nIDpy, procGrpNo )) != CPM_SUCCESS ) return ret;

  // wait
  if( (ret = Waitall( 4, &req[4] )) != CPM_SUCCESS ) return ret;

  // unpack
  if( (ret = unpackYEx( array, nmax, imax, jmax, kmax, vc, vc_comm, recvmy, recvpy, nIDmy, nIDpy )) != CPM_SUCCESS ) return ret;

  //// Z face ////
  int nIDmz = nID[Z_MINUS];
  int nIDpz = nID[Z_PLUS];

  // pack
  if( (ret = packZEx( array, nmax, imax, jmax, kmax, vc, vc_comm, sendmz, sendpz, nIDmz, nIDpz )) != CPM_SUCCESS ) return ret;

  // Isend/Irecv
  if( (ret = sendrecv( sendmz, recvmz, sendpz, recvpz, nwZ, &req[8], nIDmz, nIDmz, nIDpz, nIDpz, procGrpNo )) != CPM_SUCCESS ) return ret;

  // wait
  if( (ret = Waitall( 4, &req[8] )) != CPM_SUCCESS ) return ret;

  // unpack
  if( (ret = unpackZEx( array, nmax, imax, jmax, kmax, vc, vc_comm, recvmz, recvpz, nIDmz, nIDpz )) != CPM_SUCCESS ) return ret;

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 非同期版袖通信(Vector3DEx版)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::BndCommV3DEx_nowait( T *array, int imax, int jmax, int kmax, int vc, int vc_comm
                                    , MPI_Request req[12], int procGrpNo )
{
  return BndCommS4DEx_nowait( array, 3, imax, jmax, kmax, vc, vc_comm, req, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 非同期版袖通信(Scalar4DEx版)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::BndCommS4DEx_nowait( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                                    , MPI_Request req[12], int procGrpNo )
{
  cpm_ErrorCode ret;

  if( !array )
  {
    return CPM_ERROR_INVALID_PTR;
  }

  for( int i=0;i<12;i++ )
  {
    req[i] = MPI_REQUEST_NULL;
  }

  // 通信バッファを取得
  S_BNDCOMM_BUFFER *bufInfo = GetBndCommBuffer(procGrpNo);
  if( !bufInfo )
  {
    return CPM_ERROR_BNDCOMM_BUFFER;
  }

  // 隣接ランクを取得
  const int *nID = GetNeighborRankID(procGrpNo);
  if( !nID )
  {
    return CPM_ERROR_GET_NEIGHBOR_RANK;
  }

  // 通信バッファサイズを計算
  size_t nwX = size_t(jmax+2*vc_comm) * size_t(kmax+2*vc_comm) * size_t(vc_comm) * size_t(nmax);
  size_t nwY = size_t(kmax+2*vc_comm) * size_t(imax+2*vc_comm) * size_t(vc_comm) * size_t(nmax);
  size_t nwZ = size_t(imax+2*vc_comm) * size_t(jmax+2*vc_comm) * size_t(vc_comm) * size_t(nmax);
  if( nwX > bufInfo->m_nwX || nwY > bufInfo->m_nwY || nwZ > bufInfo->m_nwZ )
  {
    return CPM_ERROR_BNDCOMM_BUFFERLENGTH;
  }

  // ポインタをセット
  T *sendmx = (T*)(bufInfo->m_bufX[0]);
  T *recvmx = (T*)(bufInfo->m_bufX[1]);
  T *sendpx = (T*)(bufInfo->m_bufX[2]);
  T *recvpx = (T*)(bufInfo->m_bufX[3]);
  T *sendmy = (T*)(bufInfo->m_bufY[0]);
  T *recvmy = (T*)(bufInfo->m_bufY[1]);
  T *sendpy = (T*)(bufInfo->m_bufY[2]);
  T *recvpy = (T*)(bufInfo->m_bufY[3]);
  T *sendmz = (T*)(bufInfo->m_bufZ[0]);
  T *recvmz = (T*)(bufInfo->m_bufZ[1]);
  T *sendpz = (T*)(bufInfo->m_bufZ[2]);
  T *recvpz = (T*)(bufInfo->m_bufZ[3]);

  //// X face ////
  int nIDmx = nID[X_MINUS];
  int nIDpx = nID[X_PLUS];

  // pack
  if( (ret = packXEx( array, nmax, imax, jmax, kmax, vc, vc_comm, sendmx, sendpx, nIDmx, nIDpx )) != CPM_SUCCESS ) return ret;

  // Isend/Irecv
  if( (ret = sendrecv( sendmx, recvmx, sendpx, recvpx, nwX, &req[0], nIDmx, nIDmx, nIDpx, nIDpx, procGrpNo )) != CPM_SUCCESS ) return ret;

  //// Y face ////
  int nIDmy = nID[Y_MINUS];
  int nIDpy = nID[Y_PLUS];

  // pack
  if( (ret = packYEx( array, nmax, imax, jmax, kmax, vc, vc_comm, sendmy, sendpy, nIDmy, nIDpy )) != CPM_SUCCESS ) return ret;

  // Isend/Irecv
  if( (ret = sendrecv( sendmy, recvmy, sendpy, recvpy, nwY, &req[4], nIDmy, nIDmy, nIDpy, nIDpy, procGrpNo )) != CPM_SUCCESS ) return ret;

  //// Z face ////
  int nIDmz = nID[Z_MINUS];
  int nIDpz = nID[Z_PLUS];

  // pack
  if( (ret = packZEx( array, nmax, imax, jmax, kmax, vc, vc_comm, sendmz, sendpz, nIDmz, nIDpz )) != CPM_SUCCESS ) return ret;

  // Isend/Irecv
  if( (ret = sendrecv( sendmz, recvmz, sendpz, recvpz, nwZ, &req[8], nIDmz, nIDmz, nIDpz, nIDpz, procGrpNo )) != CPM_SUCCESS ) return ret;

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 非同期版袖通信のwait、展開(Vector3DEx版)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::wait_BndCommV3DEx( T *array, int imax, int jmax, int kmax, int vc, int vc_comm
                                  , MPI_Request req[12], int procGrpNo )
{
  return wait_BndCommS4DEx( array, 3, imax, jmax, kmax, vc, vc_comm, req, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 非同期版袖通信のwait、展開(Scalar4DEx版)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::wait_BndCommS4DEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                                  , MPI_Request req[12], int procGrpNo )
{
  cpm_ErrorCode ret;

  if( !array )
  {
    return CPM_ERROR_INVALID_PTR;
  }

  // 通信バッファを取得
  S_BNDCOMM_BUFFER *bufInfo = GetBndCommBuffer(procGrpNo);
  if( !bufInfo )
  {
    return CPM_ERROR_BNDCOMM_BUFFER;
  }

  // 隣接ランクを取得
  const int *nID = GetNeighborRankID(procGrpNo);
  if( !nID )
  {
    return CPM_ERROR_GET_NEIGHBOR_RANK;
  }

  // 通信バッファサイズを計算
  size_t nwX = size_t(jmax+2*vc_comm) * size_t(kmax+2*vc_comm) * size_t(vc_comm) * size_t(nmax);
  size_t nwY = size_t(kmax+2*vc_comm) * size_t(imax+2*vc_comm) * size_t(vc_comm) * size_t(nmax);
  size_t nwZ = size_t(imax+2*vc_comm) * size_t(jmax+2*vc_comm) * size_t(vc_comm) * size_t(nmax);
  if( nwX > bufInfo->m_nwX || nwY > bufInfo->m_nwY || nwZ > bufInfo->m_nwZ )
  {
    return CPM_ERROR_BNDCOMM_BUFFERLENGTH;
  }

  // ポインタをセット
  T *recvmx = (T*)(bufInfo->m_bufX[1]);
  T *recvpx = (T*)(bufInfo->m_bufX[3]);
  T *recvmy = (T*)(bufInfo->m_bufY[1]);
  T *recvpy = (T*)(bufInfo->m_bufY[3]);
  T *recvmz = (T*)(bufInfo->m_bufZ[1]);
  T *recvpz = (T*)(bufInfo->m_bufZ[3]);

  //// X face ////
  int nIDmx = nID[X_MINUS];
  int nIDpx = nID[X_PLUS];

  // wait
  if( (ret = Waitall( 4, &req[0] )) != CPM_SUCCESS ) return ret;

  // unpack
  if( (ret = unpackXEx( array, nmax, imax, jmax, kmax, vc, vc_comm, recvmx, recvpx, nIDmx, nIDpx )) != CPM_SUCCESS ) return ret;

  //// Y face ////
  int nIDmy = nID[Y_MINUS];
  int nIDpy = nID[Y_PLUS];

  // wait
  if( (ret = Waitall( 4, &req[4] )) != CPM_SUCCESS ) return ret;

  // unpack
  if( (ret = unpackYEx( array, nmax, imax, jmax, kmax, vc, vc_comm, recvmy, recvpy, nIDmy, nIDpy )) != CPM_SUCCESS ) return ret;

  //// Z face ////
  int nIDmz = nID[Z_MINUS];
  int nIDpz = nID[Z_PLUS];

  // wait
  if( (ret = Waitall( 4, &req[8] )) != CPM_SUCCESS ) return ret;

  // unpack
  if( (ret = unpackZEx( array, nmax, imax, jmax, kmax, vc, vc_comm, recvmz, recvpz, nIDmz, nIDpz )) != CPM_SUCCESS ) return ret;

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 周期境界袖通信(Vector3DEx版)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::PeriodicCommV3DEx( T *array, int imax, int jmax, int kmax, int vc, int vc_comm
                                  , cpm_DirFlag dir, cpm_PMFlag pm, int procGrpNo )
{
  return PeriodicCommS4DEx( array, 3, imax, jmax, kmax, vc, vc_comm, dir, pm, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 周期境界袖通信(Scalar4DEx版)
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::PeriodicCommS4DEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                                  , cpm_DirFlag dir, cpm_PMFlag pm, int procGrpNo )
{
  cpm_ErrorCode ret;

  if( !array )
  {
    return CPM_ERROR_INVALID_PTR;
  }

  // 通信バッファを取得
  S_BNDCOMM_BUFFER *bufInfo = GetBndCommBuffer(procGrpNo);
  if( !bufInfo )
  {
    return CPM_ERROR_BNDCOMM_BUFFER;
  }

  // 周期境界隣接ランクを取得
  const int *nID = GetPeriodicRankID(procGrpNo);
  if( !nID )
  {
    return CPM_ERROR_GET_PERIODIC_RANK;
  }

  // 通信バッファサイズを計算
  size_t nwX = size_t(jmax+2*vc_comm) * size_t(kmax+2*vc_comm) * size_t(vc_comm) * size_t(nmax);
  size_t nwY = size_t(kmax+2*vc_comm) * size_t(imax+2*vc_comm) * size_t(vc_comm) * size_t(nmax);
  size_t nwZ = size_t(imax+2*vc_comm) * size_t(jmax+2*vc_comm) * size_t(vc_comm) * size_t(nmax);
  if( nwX > bufInfo->m_nwX || nwY > bufInfo->m_nwY || nwZ > bufInfo->m_nwZ )
  {
    return CPM_ERROR_BNDCOMM_BUFFERLENGTH;
  }

  // ポインタをセット
  T *sendmx = (T*)(bufInfo->m_bufX[0]);
  T *recvmx = (T*)(bufInfo->m_bufX[1]);
  T *sendpx = (T*)(bufInfo->m_bufX[2]);
  T *recvpx = (T*)(bufInfo->m_bufX[3]);
  T *sendmy = (T*)(bufInfo->m_bufY[0]);
  T *recvmy = (T*)(bufInfo->m_bufY[1]);
  T *sendpy = (T*)(bufInfo->m_bufY[2]);
  T *recvpy = (T*)(bufInfo->m_bufY[3]);
  T *sendmz = (T*)(bufInfo->m_bufZ[0]);
  T *recvmz = (T*)(bufInfo->m_bufZ[1]);
  T *sendpz = (T*)(bufInfo->m_bufZ[2]);
  T *recvpz = (T*)(bufInfo->m_bufZ[3]);

  MPI_Request req[12];
  for( int i=0;i<12;i++ ) req[i] = MPI_REQUEST_NULL;

  //// X face ////
  int nIDsmx = cpm_ParaManager::getRankNull();
  int nIDspx = cpm_ParaManager::getRankNull();
  int nIDrmx = cpm_ParaManager::getRankNull();
  int nIDrpx = cpm_ParaManager::getRankNull();
  if( dir==X_DIR )
  {
    if( pm==MINUS2PLUS || pm==BOTH )
    {
      nIDsmx = nID[X_MINUS];
      nIDrpx = nID[X_PLUS];
    }
    if( pm==PLUS2MINUS || pm==BOTH )
    {
      nIDspx = nID[X_PLUS];
      nIDrmx = nID[X_MINUS];
    }
  }

  // pack
  if( (ret = packXEx( array, nmax, imax, jmax, kmax, vc, vc_comm, sendmx, sendpx, nIDsmx, nIDspx )) != CPM_SUCCESS ) return ret;

  // Isend/Irecv
  if( (ret = sendrecv( sendmx, recvmx, sendpx, recvpx, nwX, &req[0], nIDsmx, nIDrmx, nIDspx, nIDrpx, procGrpNo )) != CPM_SUCCESS ) return ret;

  // wait
  if( (ret = Waitall( 4, &req[0] )) != CPM_SUCCESS ) return ret;

  // unpack
  if( (ret = unpackXEx( array, nmax, imax, jmax, kmax, vc, vc_comm, recvmx, recvpx, nIDrmx, nIDrpx )) != CPM_SUCCESS ) return ret;

  //// Y face ////
  int nIDsmy = cpm_ParaManager::getRankNull();
  int nIDspy = cpm_ParaManager::getRankNull();
  int nIDrmy = cpm_ParaManager::getRankNull();
  int nIDrpy = cpm_ParaManager::getRankNull();
  if( dir==Y_DIR )
  {
    if( pm==MINUS2PLUS || pm==BOTH )
    {
      nIDsmy = nID[Y_MINUS];
      nIDrpy = nID[Y_PLUS];
    }
    if( pm==PLUS2MINUS || pm==BOTH )
    {
      nIDspy = nID[Y_PLUS];
      nIDrmy = nID[Y_MINUS];
    }
  }

  // pack
  if( (ret = packYEx( array, nmax, imax, jmax, kmax, vc, vc_comm, sendmy, sendpy, nIDsmy, nIDspy )) != CPM_SUCCESS ) return ret;

  // Isend/Irecv
  if( (ret = sendrecv( sendmy, recvmy, sendpy, recvpy, nwY, &req[4], nIDsmy, nIDrmy, nIDspy, nIDrpy, procGrpNo )) != CPM_SUCCESS ) return ret;

  // wait
  if( (ret = Waitall( 4, &req[4] )) != CPM_SUCCESS ) return ret;

  // unpack
  if( (ret = unpackYEx( array, nmax, imax, jmax, kmax, vc, vc_comm, recvmy, recvpy, nIDrmy, nIDrpy )) != CPM_SUCCESS ) return ret;

  //// Z face ////
  int nIDsmz = cpm_ParaManager::getRankNull();
  int nIDspz = cpm_ParaManager::getRankNull();
  int nIDrmz = cpm_ParaManager::getRankNull();
  int nIDrpz = cpm_ParaManager::getRankNull();
  if( dir==Z_DIR )
  {
    if( pm==MINUS2PLUS || pm==BOTH )
    {
      nIDsmz = nID[Z_MINUS];
      nIDrpz = nID[Z_PLUS];
    }
    if( pm==PLUS2MINUS || pm==BOTH )
    {
      nIDspz = nID[Z_PLUS];
      nIDrmz = nID[Z_MINUS];
    }
  }

  // pack
  if( (ret = packZEx( array, nmax, imax, jmax, kmax, vc, vc_comm, sendmz, sendpz, nIDsmz, nIDspz )) != CPM_SUCCESS ) return ret;

  // Isend/Irecv
  if( (ret = sendrecv( sendmz, recvmz, sendpz, recvpz, nwZ, &req[8], nIDsmz, nIDrmz, nIDspz, nIDrpz, procGrpNo )) != CPM_SUCCESS ) return ret;

  // wait
  if( (ret = Waitall( 4, &req[8] )) != CPM_SUCCESS ) return ret;

  // unpack
  if( (ret = unpackZEx( array, nmax, imax, jmax, kmax, vc, vc_comm, recvmz, recvpz, nIDrmz, nIDrpz )) != CPM_SUCCESS ) return ret;

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar4DEx,Vector3DEx版)のX方向送信バッファのセット
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::packXEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                        , T *sendm, T *sendp, int nIDm, int nIDp )
{
  if( !IsRankNull(nIDm) )
  {
    for( int k=0-vc_comm;k<kmax+vc_comm;k++ ){
    for( int j=0-vc_comm;j<jmax+vc_comm;j++ ){
    for( int i=0;i<vc_comm;i++ ){
    for( int n=0;n<nmax;n++){
      sendm[_IDXFX(n,i,j,k,nmax,0,jmax,kmax,vc_comm)] = array[_IDX_S4DEX(n,i,j,k,nmax,imax,jmax,kmax,vc)];
    }}}}
  }

  if( !IsRankNull(nIDp) )
  {
    for( int k=0-vc_comm;k<kmax+vc_comm;k++ ){
    for( int j=0-vc_comm;j<jmax+vc_comm;j++ ){
    for( int i=imax-vc_comm;i<imax;i++ ){
    for( int n=0;n<nmax;n++){
      sendp[_IDXFX(n,i,j,k,nmax,imax-vc_comm,jmax,kmax,vc_comm)] = array[_IDX_S4DEX(n,i,j,k,nmax,imax,jmax,kmax,vc)];
    }}}}
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar4DEx,Vector3DEx版)のX方向受信バッファを元に戻す
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::unpackXEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                          , T *recvm, T *recvp, int nIDm, int nIDp )
{
  if( !IsRankNull(nIDm) )
  {
    for( int k=0-vc_comm;k<kmax+vc_comm;k++ ){
    for( int j=0-vc_comm;j<jmax+vc_comm;j++ ){
    for( int i=0-vc_comm;i<0;i++ ){
    for( int n=0;n<nmax;n++){
      array[_IDX_S4DEX(n,i,j,k,nmax,imax,jmax,kmax,vc)] = recvm[_IDXFX(n,i,j,k,nmax,0-vc_comm,jmax,kmax,vc_comm)];
    }}}}
  }

  if( !IsRankNull(nIDp) )
  {
    for( int k=0-vc_comm;k<kmax+vc_comm;k++ ){
    for( int j=0-vc_comm;j<jmax+vc_comm;j++ ){
    for( int i=imax;i<imax+vc_comm;i++ ){
    for( int n=0;n<nmax;n++){
      array[_IDX_S4DEX(n,i,j,k,nmax,imax,jmax,kmax,vc)] = recvp[_IDXFX(n,i,j,k,nmax,imax,jmax,kmax,vc_comm)];
    }}}}
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar4DEx,Vector3DEx版)のY方向送信バッファのセット
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::packYEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                        , T *sendm, T *sendp, int nIDm, int nIDp )
{
  if( !IsRankNull(nIDm) )
  {
    for( int k=0-vc_comm;k<kmax+vc_comm;k++ ){
    for( int j=0;j<vc_comm;j++ ){
    for( int i=0-vc_comm;i<imax+vc_comm;i++ ){
    for( int n=0;n<nmax;n++){
      sendm[_IDXFY(n,i,j,k,nmax,imax,0,kmax,vc_comm)] = array[_IDX_S4DEX(n,i,j,k,nmax,imax,jmax,kmax,vc)];
    }}}}
  }

  if( !IsRankNull(nIDp) )
  {
    for( int k=0-vc_comm;k<kmax+vc_comm;k++ ){
    for( int j=jmax-vc_comm;j<jmax;j++ ){
    for( int i=0-vc_comm;i<imax+vc_comm;i++ ){
    for( int n=0;n<nmax;n++){
      sendp[_IDXFY(n,i,j,k,nmax,imax,jmax-vc_comm,kmax,vc_comm)] = array[_IDX_S4DEX(n,i,j,k,nmax,imax,jmax,kmax,vc)];
    }}}}
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar4DEx,Vector3DEx版)のY方向受信バッファを元に戻す
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::unpackYEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                          , T *recvm, T *recvp, int nIDm, int nIDp )
{
  if( !IsRankNull(nIDm) )
  {
    for( int k=0-vc_comm;k<kmax+vc_comm;k++ ){
    for( int j=0-vc_comm;j<0;j++ ){
    for( int i=0-vc_comm;i<imax+vc_comm;i++ ){
    for( int n=0;n<nmax;n++){
      array[_IDX_S4DEX(n,i,j,k,nmax,imax,jmax,kmax,vc)] = recvm[_IDXFY(n,i,j,k,nmax,imax,0-vc_comm,kmax,vc_comm)];
    }}}}
  }

  if( !IsRankNull(nIDp) )
  {
    for( int k=0-vc_comm;k<kmax+vc_comm;k++ ){
    for( int j=jmax;j<jmax+vc_comm;j++ ){
    for( int i=0-vc_comm;i<imax+vc_comm;i++ ){
    for( int n=0;n<nmax;n++){
      array[_IDX_S4DEX(n,i,j,k,nmax,imax,jmax,kmax,vc)] = recvp[_IDXFY(n,i,j,k,nmax,imax,jmax,kmax,vc_comm)];
    }}}}
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar4DEx,Vector3DEx版)のZ方向送信バッファのセット
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::packZEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                        , T *sendm, T *sendp, int nIDm, int nIDp )
{
  if( !IsRankNull(nIDm) )
  {
    for( int k=0;k<vc_comm;k++ ){
    for( int j=0-vc_comm;j<jmax+vc_comm;j++ ){
    for( int i=0-vc_comm;i<imax+vc_comm;i++ ){
    for( int n=0;n<nmax;n++){
      sendm[_IDXFZ(n,i,j,k,nmax,imax,jmax,0,vc_comm)] = array[_IDX_S4DEX(n,i,j,k,nmax,imax,jmax,kmax,vc)];
    }}}}
  }

  if( !IsRankNull(nIDp) )
  {
    for( int k=kmax-vc_comm;k<kmax;k++ ){
    for( int j=0-vc_comm;j<jmax+vc_comm;j++ ){
    for( int i=0-vc_comm;i<imax+vc_comm;i++ ){
    for( int n=0;n<nmax;n++){
      sendp[_IDXFZ(n,i,j,k,nmax,imax,jmax,kmax-vc_comm,vc_comm)] = array[_IDX_S4DEX(n,i,j,k,nmax,imax,jmax,kmax,vc)];
    }}}}
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar4DEx,Vector3DEx版)のZ方向受信バッファを元に戻す
template<class T> CPM_INLINE
cpm_ErrorCode
cpm_ParaManager::unpackZEx( T *array, int nmax, int imax, int jmax, int kmax, int vc, int vc_comm
                          , T *recvm, T *recvp, int nIDm, int nIDp )
{
  if( !IsRankNull(nIDm) )
  {
    for( int k=0-vc_comm;k<0;k++ ){
    for( int j=0-vc_comm;j<jmax+vc_comm;j++ ){
    for( int i=0-vc_comm;i<imax+vc_comm;i++ ){
    for( int n=0;n<nmax;n++){
      array[_IDX_S4DEX(n,i,j,k,nmax,imax,jmax,kmax,vc)] = recvm[_IDXFZ(n,i,j,k,nmax,imax,jmax,0-vc_comm,vc_comm)];
    }}}}
  }

  if( !IsRankNull(nIDp) )
  {
    for( int k=kmax;k<kmax+vc_comm;k++ ){
    for( int j=0-vc_comm;j<jmax+vc_comm;j++ ){
    for( int i=0-vc_comm;i<imax+vc_comm;i++ ){
    for( int n=0;n<nmax;n++){
      array[_IDX_S4DEX(n,i,j,k,nmax,imax,jmax,kmax,vc)] = recvp[_IDXFZ(n,i,j,k,nmax,imax,jmax,kmax,vc_comm)];
    }}}}
  }

  return CPM_SUCCESS;
}

#undef _IDXFX
#undef _IDXFY
#undef _IDXFZ

#endif /* _CPM_PARAMANAGER_BNDCOMMEX_H_ */

