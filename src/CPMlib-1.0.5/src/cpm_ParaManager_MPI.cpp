/*
 * CPM - Cartesian Partition Manager
 *
 * Copyright (C) 2012 University of Tokyo.
 *
 */

/**
 * @file   cpm_ParaManager_MPI.cpp
 * @brief  パラレルマネージャクラスのMPIインターフェイス関数ソースファイル
 * @author University of Tokyo
 * @date   2012/05/31
 */
#include "stdlib.h"
#include "cpm_ParaManager.h"
#include <unistd.h> // for gethostname()

////////////////////////////////////////////////////////////////////////////////
// MPI_Datatypeを取得
MPI_Datatype
cpm_ParaManager::GetMPI_Datatype(int datatype)
{
  if( datatype == CPM_REAL )
  {
    if( RealIsDouble() ) return MPI_DOUBLE;
    return MPI_FLOAT;
  }
  else if( datatype == CPM_CHAR )               return MPI_CHAR;
  else if( datatype == CPM_SHORT )              return MPI_SHORT;
  else if( datatype == CPM_INT )                return MPI_INT;
  else if( datatype == CPM_LONG )               return MPI_LONG;
  else if( datatype == CPM_FLOAT )              return MPI_FLOAT;
  else if( datatype == CPM_DOUBLE )             return MPI_DOUBLE;
  else if( datatype == CPM_LONG_DOUBLE )        return MPI_LONG_DOUBLE;
  else if( datatype == CPM_UNSIGNED_CHAR )      return MPI_UNSIGNED_CHAR;
  else if( datatype == CPM_UNSIGNED_SHORT )     return MPI_UNSIGNED_SHORT;
  else if( datatype == CPM_UNSIGNED )           return MPI_UNSIGNED;
  else if( datatype == CPM_UNSIGNED_LONG )      return MPI_UNSIGNED_LONG;
#ifdef MPI_LONG_LONG_INT
  else if( datatype == CPM_LONG_LONG_INT )      return MPI_LONG_LONG_INT;
#endif
#ifdef MPI_LONG_LONG
  else if( datatype == CPM_LONG_LONG )          return MPI_LONG_LONG;
#endif
#ifdef MPI_UNSIGNED_LONG_LONG
  else if( datatype == CPM_UNSIGNED_LONG_LONG ) return MPI_UNSIGNED_LONG_LONG;
#endif

  return MPI_DATATYPE_NULL;
}

////////////////////////////////////////////////////////////////////////////////
// MPI_Opを取得
MPI_Op
cpm_ParaManager::GetMPI_Op(int op)
{
  if     ( op == CPM_MAX    ) return MPI_MAX;
  else if( op == CPM_MIN    ) return MPI_MIN;
  else if( op == CPM_SUM    ) return MPI_SUM;
  else if( op == CPM_PROD   ) return MPI_PROD;
  else if( op == CPM_LAND   ) return MPI_LAND;
  else if( op == CPM_BAND   ) return MPI_BAND;
  else if( op == CPM_LOR    ) return MPI_LOR;
  else if( op == CPM_BOR    ) return MPI_BOR;
  else if( op == CPM_LXOR   ) return MPI_LXOR;
  else if( op == CPM_BXOR   ) return MPI_BXOR;
//  else if( op == CPM_MINLOC ) return CPM_MINLOC; // not support
//  else if( op == CPM_MAXLOC ) return CPM_MAXLOC; // not support

  return MPI_OP_NULL;
}

////////////////////////////////////////////////////////////////////////////////
// ランク番号の取得
int
cpm_ParaManager::GetMyRankID( int procGrpNo )
{
  // 不正なプロセスグループ番号
  if( procGrpNo < 0 || procGrpNo >= int(m_procGrpList.size()) )
  {
    // プロセスグループが存在しない
    return getRankNull();
  }

  // コミュニケータをチェック
  MPI_Comm comm = m_procGrpList[procGrpNo];
  if( IsCommNull(comm) )
  {
    // プロセスグループに自ランクが含まれない
   return getRankNull();
  }

  // ランク番号を取得
  int rankNo;
  MPI_Comm_rank( comm, &rankNo );

  // ランク番号
  return rankNo;
}

////////////////////////////////////////////////////////////////////////////////
// ランク数の取得
int
cpm_ParaManager::GetNumRank( int procGrpNo )
{
  // 不正なプロセスグループ番号
  if( procGrpNo < 0 || procGrpNo >= int(m_procGrpList.size()) )
  {
    // プロセスグループが存在しない
    return -1;
  }

  // コミュニケータをチェック
  MPI_Comm comm = m_procGrpList[procGrpNo];
  if( IsCommNull(comm) )
  {
    // プロセスグループに自ランクが含まれない
    return -1;
  }

  // ランク数を取得
  int nrank;
  MPI_Comm_size( comm, &nrank );

  // ランク数
  return nrank;
}

////////////////////////////////////////////////////////////////////////////////
// ホスト名の取得
std::string
cpm_ParaManager::GetHostName()
{
  char name[512];
  memset(name, 0x00, sizeof(char)*512);
  if( gethostname(name, 512) != 0 ) return std::string("");
  return std::string(name);
}

////////////////////////////////////////////////////////////////////////////////
// コミュニケータの取得
MPI_Comm
cpm_ParaManager::GetMPI_Comm( int procGrpNo )
{
  // 不正なプロセスグループ番号
  if( procGrpNo < 0 || procGrpNo >= int(m_procGrpList.size()) )
  {
    // プロセスグループが存在しない
    return getCommNull();
  }

  return m_procGrpList[procGrpNo];
}

////////////////////////////////////////////////////////////////////////////////
// Abort
void
cpm_ParaManager::Abort( int errorcode )
{
  // MPI_Abort
  MPI_Abort(MPI_COMM_WORLD, errorcode);
  exit(errorcode);
}

////////////////////////////////////////////////////////////////////////////////
// Barrier
cpm_ErrorCode
cpm_ParaManager::Barrier( int procGrpNo )
{
  // コミュニケータを取得
  MPI_Comm comm = GetMPI_Comm( procGrpNo );
  if( IsCommNull( comm ) )
  {
    // プロセスグループが存在しない
    return CPM_ERROR_NOT_IN_PROCGROUP;
  }

  // MPI_Barrier
  if( MPI_Barrier(comm) != MPI_SUCCESS )
  {
    return CPM_ERROR_MPI_BARRIER;
  }
  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// Wait
cpm_ErrorCode
cpm_ParaManager::Wait( MPI_Request *request )
{
  if( !request )
  {
    return CPM_ERROR_INVALID_PTR;
  }
  if( *request == MPI_REQUEST_NULL )
  {
    return CPM_ERROR_MPI_INVALID_REQUEST;
  }

  // MPI_Wait
  MPI_Status status;
  if( MPI_Wait( request, &status ) != MPI_SUCCESS )
  {
    return CPM_ERROR_MPI_WAIT;
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// Waitall
cpm_ErrorCode
cpm_ParaManager::Waitall( int count, MPI_Request requests[] )
{
  // status
  int cnt = 0;
  MPI_Status  *stat = new MPI_Status[count];
  MPI_Request *req  = new MPI_Request[count];
  for( int i=0;i<count;i++ )
  {
    if( requests[i] != MPI_REQUEST_NULL ) req[cnt++] = requests[i];
  }
  if( cnt == 0 )
  {
    delete [] stat;
    delete [] req;
    return CPM_SUCCESS;
  }

  // MPI_Waitall
  if( MPI_Waitall( cnt, req, stat ) != MPI_SUCCESS )
  {
    delete [] stat;
    delete [] req;
    return CPM_ERROR_MPI_WAITALL;
  }

  delete [] stat;
  delete [] req;
  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// Bcast(MPI_Datatype指定)
cpm_ErrorCode
cpm_ParaManager::Bcast( MPI_Datatype dtype, void *buf, int count, int root, int procGrpNo )
{
  if( !buf )
  {
    return CPM_ERROR_INVALID_PTR;
  }

  // コミュニケータを取得
  MPI_Comm comm = GetMPI_Comm(procGrpNo);
  if( IsCommNull(comm) )
  {
    // プロセスグループが存在しない
    return CPM_ERROR_NOT_IN_PROCGROUP;
  }

  // MPI_Bcast
  if( MPI_Bcast( buf, count, dtype, root, comm ) != MPI_SUCCESS )
  {
    return CPM_ERROR_MPI_BCAST;
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// Send(MPI_Datatype指定)
cpm_ErrorCode
cpm_ParaManager::Send( MPI_Datatype dtype, void *buf, int count, int dest, int procGrpNo )
{
  if( !buf )
  {
    return CPM_ERROR_INVALID_PTR;
  }

  // コミュニケータを取得
  MPI_Comm comm = GetMPI_Comm(procGrpNo);
  if( IsCommNull(comm) )
  {
    // プロセスグループが存在しない
    return CPM_ERROR_NOT_IN_PROCGROUP;
  }

  // MPI_Send
  int tag = 1;
  if( MPI_Send( buf, count, dtype, dest, tag, comm ) != MPI_SUCCESS )
  {
    return CPM_ERROR_MPI_SEND;
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// Recv(MPI_Datatype指定)
cpm_ErrorCode
cpm_ParaManager::Recv( MPI_Datatype dtype, void *buf, int count, int source, int procGrpNo ) 
{
  if( !buf )
  {
    return CPM_ERROR_INVALID_PTR;
  }

  // コミュニケータを取得
  MPI_Comm comm = GetMPI_Comm(procGrpNo);
  if( IsCommNull(comm) )
  {
    // プロセスグループが存在しない
    return CPM_ERROR_NOT_IN_PROCGROUP;
  }

  // MPI_Send
  int tag = 1;
  MPI_Status status;
  if( MPI_Recv( buf, count, dtype, source, tag, comm, &status ) != MPI_SUCCESS )
  {
    return CPM_ERROR_MPI_SEND;
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// Isend(MPI_Datatype指定)
cpm_ErrorCode
cpm_ParaManager::Isend( MPI_Datatype dtype, void *buf, int count, int dest, MPI_Request *request, int procGrpNo )
{
  if( !buf || !request )
  {
    return CPM_ERROR_INVALID_PTR;
  }
  *request = MPI_REQUEST_NULL;

  // コミュニケータを取得
  MPI_Comm comm = GetMPI_Comm(procGrpNo);
  if( IsCommNull(comm) )
  {
    // プロセスグループが存在しない
    return CPM_ERROR_NOT_IN_PROCGROUP;
  }

  // MPI_Isend
  int tag = 1;
  if( MPI_Isend( buf, count, dtype, dest, tag, comm, request ) != MPI_SUCCESS )
  {
    return CPM_ERROR_MPI_ISEND;
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// Irecv(MPI_Datatype指定)
cpm_ErrorCode
cpm_ParaManager::Irecv( MPI_Datatype dtype, void *buf, int count, int source, MPI_Request *request, int procGrpNo )
{
  if( !buf || !request )
  {
    return CPM_ERROR_INVALID_PTR;
  }
  *request = MPI_REQUEST_NULL;

  // コミュニケータを取得
  MPI_Comm comm = GetMPI_Comm(procGrpNo);
  if( IsCommNull(comm) )
  {
    // プロセスグループが存在しない
    return CPM_ERROR_NOT_IN_PROCGROUP;
  }

  // MPI_Irecv
  int tag = 1;
  if( MPI_Irecv( buf, count, dtype, source, tag, comm, request ) != MPI_SUCCESS )
  {
    return CPM_ERROR_MPI_IRECV;
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// Allreduce(MPI_Datatype指定)
cpm_ErrorCode
cpm_ParaManager::Allreduce( MPI_Datatype dtype, void *sendbuf, void *recvbuf, int count, MPI_Op op, int procGrpNo )
{
  if( !sendbuf || !recvbuf )
  {
    return CPM_ERROR_INVALID_PTR;
  }

  // コミュニケータを取得
  MPI_Comm comm = GetMPI_Comm(procGrpNo);
  if( IsCommNull(comm) )
  {
    // プロセスグループが存在しない
    return CPM_ERROR_NOT_IN_PROCGROUP;
  }

  // MPI_Allreduce
  if( MPI_Allreduce( sendbuf, recvbuf, count, dtype, op, comm ) != MPI_SUCCESS )
  {
    return CPM_ERROR_MPI_ALLREDUCE;
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// Gather(MPI_Datatype指定)
cpm_ErrorCode
cpm_ParaManager::Gather( MPI_Datatype stype, void *sendbuf, int sendcnt
                       , MPI_Datatype rtype, void *recvbuf, int recvcnt
                       , int root, int procGrpNo )
{
  if( !sendbuf || !recvbuf )
  {
    return CPM_ERROR_INVALID_PTR;
  }

  //コミュニケータを取得
  MPI_Comm comm = GetMPI_Comm(procGrpNo);
  if( IsCommNull(comm) )
  {
    // プロセスグループが存在しない
    return CPM_ERROR_NOT_IN_PROCGROUP;
  }

  // MPI_Gather
  if( MPI_Gather( sendbuf, sendcnt, stype, recvbuf, recvcnt, rtype, root, comm ) != MPI_SUCCESS )
  {
    return CPM_ERROR_MPI_GATHER;
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// Allgather(MPI_Datatype指定)
cpm_ErrorCode
cpm_ParaManager::Allgather( MPI_Datatype stype, void *sendbuf, int sendcnt
                          , MPI_Datatype rtype, void *recvbuf, int recvcnt
                          , int procGrpNo )
{
  if( !sendbuf || !recvbuf )
  {
    return CPM_ERROR_INVALID_PTR;
  }

  //コミュニケータを取得
  MPI_Comm comm = GetMPI_Comm(procGrpNo);
  if( IsCommNull(comm) )
  {
    // プロセスグループが存在しない
    return CPM_ERROR_NOT_IN_PROCGROUP;
  }

  // MPI_Allaather
  if( MPI_Allgather( sendbuf, sendcnt, stype, recvbuf, recvcnt, rtype, comm ) != MPI_SUCCESS )
  {
    return CPM_ERROR_MPI_ALLGATHER;
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// Gatherv(MPI_Datatype指定)
cpm_ErrorCode
cpm_ParaManager::Gatherv( MPI_Datatype stype, void *sendbuf, int sendcnt
                        , MPI_Datatype rtype, void *recvbuf, int *recvcnts
                        , int *displs, int root, int procGrpNo )
{
  if( !sendbuf || !recvbuf || !recvcnts || !displs )
  {
    return CPM_ERROR_INVALID_PTR;
  }

  //コミュニケータを取得
  MPI_Comm comm = GetMPI_Comm(procGrpNo);
  if( IsCommNull(comm) )
  {
    // プロセスグループが存在しない
    return CPM_ERROR_NOT_IN_PROCGROUP;
  }

  // MPI_Gather
  if( MPI_Gatherv( sendbuf, sendcnt, stype, recvbuf, recvcnts, displs
                 , rtype, root, comm ) != MPI_SUCCESS )
  {
    return CPM_ERROR_MPI_GATHERV;
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// Allgatherv(MPI_Datatype指定)
cpm_ErrorCode
cpm_ParaManager::Allgatherv( MPI_Datatype stype, void *sendbuf, int sendcnt
                           , MPI_Datatype rtype, void *recvbuf, int *recvcnts
                           , int *displs, int procGrpNo )
{
  if( !sendbuf || !recvbuf || !recvcnts || !displs )
  {
    return CPM_ERROR_INVALID_PTR;
  }

  // コミュニケータを取得
  MPI_Comm comm = GetMPI_Comm(procGrpNo);
  if( IsCommNull(comm) )
  {
    // プロセスグループが存在しない
    return CPM_ERROR_NOT_IN_PROCGROUP;
  }

  // MPI_Allaather
  if( MPI_Allgatherv( sendbuf, sendcnt, stype, recvbuf, recvcnts, displs, rtype, comm ) != MPI_SUCCESS )
  {
    return CPM_ERROR_MPI_ALLGATHERV;
  }

  return CPM_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar3D版, MPI_Datatype指定)
cpm_ErrorCode
cpm_ParaManager::BndCommS3D( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                           , int vc, int vc_comm, int procGrpNo ) 
{
  return BndCommS4D( dtype, array, imax, jmax, kmax, 1, vc, vc_comm, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Vector3D版, MPI_Datatype指定)
cpm_ErrorCode
cpm_ParaManager::BndCommV3D( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                           , int vc, int vc_comm, int procGrpNo ) 
{
  return BndCommS4D( dtype, array, imax, jmax, kmax, 3, vc, vc_comm, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar4D版, MPI_Datatype指定)
cpm_ErrorCode
cpm_ParaManager::BndCommS4D( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax, int nmax
                           , int vc, int vc_comm, int procGrpNo ) 
{
  if( dtype == MPI_CHAR )
    return BndCommS4D( (char*)array, imax, jmax, kmax, nmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_SHORT )
    return BndCommS4D( (short*)array, imax, jmax, kmax, nmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_INT )
    return BndCommS4D( (int*)array, imax, jmax, kmax, nmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_LONG )
    return BndCommS4D( (long*)array, imax, jmax, kmax, nmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_FLOAT )
    return BndCommS4D( (float*)array, imax, jmax, kmax, nmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_DOUBLE )
    return BndCommS4D( (double*)array, imax, jmax, kmax, nmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_LONG_DOUBLE )
    return BndCommS4D( (long double*)array, imax, jmax, kmax, nmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_UNSIGNED_CHAR )
    return BndCommS4D( (unsigned char*)array, imax, jmax, kmax, nmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_UNSIGNED_SHORT )
    return BndCommS4D( (unsigned short*)array, imax, jmax, kmax, nmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_UNSIGNED )
    return BndCommS4D( (unsigned*)array, imax, jmax, kmax, nmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_UNSIGNED_LONG )
    return BndCommS4D( (unsigned long*)array, imax, jmax, kmax, nmax, vc, vc_comm, procGrpNo );
#ifdef MPI_LONG_LONG_INT
  else if( dtype == MPI_LONG_LONG_INT )
    return BndCommS4D( (long long int*)array, imax, jmax, kmax, nmax, vc, vc_comm, procGrpNo );
#endif
#ifdef MPI_LONG_LONG
  else if( dtype == MPI_LONG_LONG )
    return BndCommS4D( (long long*)array, imax, jmax, kmax, nmax, vc, vc_comm, procGrpNo );
#endif
#ifdef MPI_UNSIGNED_LONG_LONG
  else if( dtype == MPI_UNSIGNED_LONG_LONG )
    return BndCommS4D( (unsigned long long*)array, imax, jmax, kmax, nmax, vc, vc_comm, procGrpNo );
#endif

  return CPM_ERROR_MPI_INVALID_DATATYPE;
}

////////////////////////////////////////////////////////////////////////////////
// 非同期版袖通信(Scalar3D版, MPI_Datatype指定)
cpm_ErrorCode
cpm_ParaManager::BndCommS3D_nowait( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                                   , int vc, int vc_comm, MPI_Request req[12], int procGrpNo ) 
{
  return BndCommS4D_nowait( dtype, array, imax, jmax, kmax, 1, vc, vc_comm, req, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 非同期版袖通信(Vector3D版, MPI_Datatype指定)
cpm_ErrorCode
cpm_ParaManager::BndCommV3D_nowait( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                                   , int vc, int vc_comm, MPI_Request req[12], int procGrpNo ) 
{
  return BndCommS4D_nowait( dtype, array, imax, jmax, kmax, 3, vc, vc_comm, req, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 非同期版袖通信(Scalar4D版, MPI_Datatype指定)
cpm_ErrorCode
cpm_ParaManager::BndCommS4D_nowait( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax, int nmax
                                  , int vc, int vc_comm, MPI_Request req[12], int procGrpNo ) 
{
  if( dtype == MPI_CHAR )
    return BndCommS4D_nowait( (char*)array, imax, jmax, kmax, nmax, vc, vc_comm, req, procGrpNo );
  else if( dtype == MPI_SHORT )
    return BndCommS4D_nowait( (short*)array, imax, jmax, kmax, nmax, vc, vc_comm, req, procGrpNo );
  else if( dtype == MPI_INT )
    return BndCommS4D_nowait( (int*)array, imax, jmax, kmax, nmax, vc, vc_comm, req, procGrpNo );
  else if( dtype == MPI_LONG )
    return BndCommS4D_nowait( (long*)array, imax, jmax, kmax, nmax, vc, vc_comm, req, procGrpNo );
  else if( dtype == MPI_FLOAT )
    return BndCommS4D_nowait( (float*)array, imax, jmax, kmax, nmax, vc, vc_comm, req, procGrpNo );
  else if( dtype == MPI_DOUBLE )
    return BndCommS4D_nowait( (double*)array, imax, jmax, kmax, nmax, vc, vc_comm, req, procGrpNo );
  else if( dtype == MPI_LONG_DOUBLE )
    return BndCommS4D_nowait( (long double*)array, imax, jmax, kmax, nmax, vc, vc_comm, req, procGrpNo );
  else if( dtype == MPI_UNSIGNED_CHAR )
    return BndCommS4D_nowait( (unsigned char*)array, imax, jmax, kmax, nmax, vc, vc_comm, req, procGrpNo );
  else if( dtype == MPI_UNSIGNED_SHORT )
    return BndCommS4D_nowait( (unsigned short*)array, imax, jmax, kmax, nmax, vc, vc_comm, req, procGrpNo );
  else if( dtype == MPI_UNSIGNED )
    return BndCommS4D_nowait( (unsigned*)array, imax, jmax, kmax, nmax, vc, vc_comm, req, procGrpNo );
  else if( dtype == MPI_UNSIGNED_LONG )
    return BndCommS4D_nowait( (unsigned long*)array, imax, jmax, kmax, nmax, vc, vc_comm, req, procGrpNo );
#ifdef MPI_LONG_LONG_INT
  else if( dtype == MPI_LONG_LONG_INT )
    return BndCommS4D_nowait( (long long int*)array, imax, jmax, kmax, nmax, vc, vc_comm, req, procGrpNo );
#endif
#ifdef MPI_LONG_LONG
  else if( dtype == MPI_LONG_LONG )
    return BndCommS4D_nowait( (long long*)array, imax, jmax, kmax, nmax, vc, vc_comm, req, procGrpNo );
#endif
#ifdef MPI_UNSIGNED_LONG_LONG
  else if( dtype == MPI_UNSIGNED_LONG_LONG )
    return BndCommS4D_nowait( (unsigned long long*)array, imax, jmax, kmax, nmax, vc, vc_comm, req, procGrpNo );
#endif

  return CPM_ERROR_MPI_INVALID_DATATYPE;
}

////////////////////////////////////////////////////////////////////////////////
// 非同期版袖通信のwait、展開(Scalar3D版, MPI_Datatype指定)
cpm_ErrorCode
cpm_ParaManager::wait_BndCommS3D( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                                , int vc, int vc_comm, MPI_Request req[12], int procGrpNo )
{
  return wait_BndCommS4D( dtype, array, imax, jmax, kmax, 1, vc, vc_comm, req, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 非同期版袖通信のwait、展開(Vector3D版, MPI_Datatype指定)
cpm_ErrorCode
cpm_ParaManager::wait_BndCommV3D( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                                , int vc, int vc_comm, MPI_Request req[12], int procGrpNo )
{
  return wait_BndCommS4D( dtype, array, imax, jmax, kmax, 3, vc, vc_comm, req, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 非同期版袖通信のwait、展開(Scalar4D版, MPI_Datatype指定)
cpm_ErrorCode
cpm_ParaManager::wait_BndCommS4D( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax, int nmax
                                , int vc, int vc_comm, MPI_Request req[12], int procGrpNo ) 
{
  if( dtype == MPI_CHAR )
    return wait_BndCommS4D( (char*)array, imax, jmax, kmax, nmax, vc, vc_comm, req, procGrpNo );
  else if( dtype == MPI_SHORT )
    return wait_BndCommS4D( (short*)array, imax, jmax, kmax, nmax, vc, vc_comm, req, procGrpNo );
  else if( dtype == MPI_INT )
    return wait_BndCommS4D( (int*)array, imax, jmax, kmax, nmax, vc, vc_comm, req, procGrpNo );
  else if( dtype == MPI_LONG )
    return wait_BndCommS4D( (long*)array, imax, jmax, kmax, nmax, vc, vc_comm, req, procGrpNo );
  else if( dtype == MPI_FLOAT )
    return wait_BndCommS4D( (float*)array, imax, jmax, kmax, nmax, vc, vc_comm, req, procGrpNo );
  else if( dtype == MPI_DOUBLE )
    return wait_BndCommS4D( (double*)array, imax, jmax, kmax, nmax, vc, vc_comm, req, procGrpNo );
  else if( dtype == MPI_LONG_DOUBLE )
    return wait_BndCommS4D( (long double*)array, imax, jmax, kmax, nmax, vc, vc_comm, req, procGrpNo );
  else if( dtype == MPI_UNSIGNED_CHAR )
    return wait_BndCommS4D( (unsigned char*)array, imax, jmax, kmax, nmax, vc, vc_comm, req, procGrpNo );
  else if( dtype == MPI_UNSIGNED_SHORT )
    return wait_BndCommS4D( (unsigned short*)array, imax, jmax, kmax, nmax, vc, vc_comm, req, procGrpNo );
  else if( dtype == MPI_UNSIGNED )
    return wait_BndCommS4D( (unsigned*)array, imax, jmax, kmax, nmax, vc, vc_comm, req, procGrpNo );
  else if( dtype == MPI_UNSIGNED_LONG )
    return wait_BndCommS4D( (unsigned long*)array, imax, jmax, kmax, nmax, vc, vc_comm, req, procGrpNo );
#ifdef MPI_LONG_LONG_INT
  else if( dtype == MPI_LONG_LONG_INT )
    return wait_BndCommS4D( (long long int*)array, imax, jmax, kmax, nmax, vc, vc_comm, req, procGrpNo );
#endif
#ifdef MPI_LONG_LONG
  else if( dtype == MPI_LONG_LONG )
    return wait_BndCommS4D( (long long*)array, imax, jmax, kmax, nmax, vc, vc_comm, req, procGrpNo );
#endif
#ifdef MPI_UNSIGNED_LONG_LONG
  else if( dtype == MPI_UNSIGNED_LONG_LONG )
    return wait_BndCommS4D( (unsigned long long*)array, imax, jmax, kmax, nmax, vc, vc_comm, req, procGrpNo );
#endif

  return CPM_ERROR_MPI_INVALID_DATATYPE;
}

////////////////////////////////////////////////////////////////////////////////
// 周期境界袖通信(Scalar3D版, MPI_Datatype指定)
cpm_ErrorCode
cpm_ParaManager::PeriodicCommS3D( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                                , int vc, int vc_comm, cpm_DirFlag dir, cpm_PMFlag pm, int procGrpNo )
{
  return PeriodicCommS4D( dtype, array, imax, jmax, kmax, 1, vc, vc_comm, dir, pm, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 周期境界袖通信(Vector3D版, MPI_Datatype指定)
cpm_ErrorCode
cpm_ParaManager::PeriodicCommV3D( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                                , int vc, int vc_comm, cpm_DirFlag dir, cpm_PMFlag pm, int procGrpNo )
{
  return PeriodicCommS4D( dtype, array, imax, jmax, kmax, 3, vc, vc_comm, dir, pm, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 周期境界袖通信(Scalar4D版, MPI_Datatype指定)
cpm_ErrorCode
cpm_ParaManager::PeriodicCommS4D( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax, int nmax
                                , int vc, int vc_comm, cpm_DirFlag dir, cpm_PMFlag pm, int procGrpNo )
{
  if( dtype == MPI_CHAR )
    return PeriodicCommS4D( (char*)array, imax, jmax, kmax, nmax, vc, vc_comm, dir, pm, procGrpNo );
  else if( dtype == MPI_SHORT )
    return PeriodicCommS4D( (short*)array, imax, jmax, kmax, nmax, vc, vc_comm, dir, pm, procGrpNo );
  else if( dtype == MPI_INT )
    return PeriodicCommS4D( (int*)array, imax, jmax, kmax, nmax, vc, vc_comm, dir, pm, procGrpNo );
  else if( dtype == MPI_LONG )
    return PeriodicCommS4D( (long*)array, imax, jmax, kmax, nmax, vc, vc_comm, dir, pm, procGrpNo );
  else if( dtype == MPI_FLOAT )
    return PeriodicCommS4D( (float*)array, imax, jmax, kmax, nmax, vc, vc_comm, dir, pm, procGrpNo );
  else if( dtype == MPI_DOUBLE )
    return PeriodicCommS4D( (double*)array, imax, jmax, kmax, nmax, vc, vc_comm, dir, pm, procGrpNo );
  else if( dtype == MPI_LONG_DOUBLE )
    return PeriodicCommS4D( (long double*)array, imax, jmax, kmax, nmax, vc, vc_comm, dir, pm, procGrpNo );
  else if( dtype == MPI_UNSIGNED_CHAR )
    return PeriodicCommS4D( (unsigned char*)array, imax, jmax, kmax, nmax, vc, vc_comm, dir, pm, procGrpNo );
  else if( dtype == MPI_UNSIGNED_SHORT )
    return PeriodicCommS4D( (unsigned short*)array, imax, jmax, kmax, nmax, vc, vc_comm, dir, pm, procGrpNo );
  else if( dtype == MPI_UNSIGNED )
    return PeriodicCommS4D( (unsigned*)array, imax, jmax, kmax, nmax, vc, vc_comm, dir, pm, procGrpNo );
  else if( dtype == MPI_UNSIGNED_LONG )
    return PeriodicCommS4D( (unsigned long*)array, imax, jmax, kmax, nmax, vc, vc_comm, dir, pm, procGrpNo );
#ifdef MPI_LONG_LONG_INT
  else if( dtype == MPI_LONG_LONG_INT )
    return PeriodicCommS4D( (long long int*)array, imax, jmax, kmax, nmax, vc, vc_comm, dir, pm, procGrpNo );
#endif
#ifdef MPI_LONG_LONG
  else if( dtype == MPI_LONG_LONG )
    return PeriodicCommS4D( (long long*)array, imax, jmax, kmax, nmax, vc, vc_comm, dir, pm, procGrpNo );
#endif
#ifdef MPI_UNSIGNED_LONG_LONG
  else if( dtype == MPI_UNSIGNED_LONG_LONG )
    return PeriodicCommS4D( (unsigned long long*)array, imax, jmax, kmax, nmax, vc, vc_comm, dir, pm, procGrpNo );
#endif

  return CPM_ERROR_MPI_INVALID_DATATYPE;
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Vector3DEx版, MPI_Datatype指定)
cpm_ErrorCode
cpm_ParaManager::BndCommV3DEx( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                             , int vc, int vc_comm, int procGrpNo )
{
  return BndCommS4DEx( dtype, array, 3, imax, jmax, kmax, vc, vc_comm, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 袖通信(Scalar4DEx版, MPI_Datatype指定)
cpm_ErrorCode
cpm_ParaManager::BndCommS4DEx( MPI_Datatype dtype, void *array, int nmax, int imax, int jmax, int kmax
                             , int vc, int vc_comm, int procGrpNo )
{
  if( dtype == MPI_CHAR )
    return BndCommS4DEx( (char*)array, nmax, imax, jmax, kmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_SHORT )
    return BndCommS4DEx( (short*)array, nmax, imax, jmax, kmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_INT )
    return BndCommS4DEx( (int*)array, nmax, imax, jmax, kmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_LONG )
    return BndCommS4DEx( (long*)array, nmax, imax, jmax, kmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_FLOAT )
    return BndCommS4DEx( (float*)array, nmax, imax, jmax, kmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_DOUBLE )
    return BndCommS4DEx( (double*)array, nmax, imax, jmax, kmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_LONG_DOUBLE )
    return BndCommS4DEx( (long double*)array, nmax, imax, jmax, kmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_UNSIGNED_CHAR )
    return BndCommS4DEx( (unsigned char*)array, nmax, imax, jmax, kmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_UNSIGNED_SHORT )
    return BndCommS4DEx( (unsigned short*)array, nmax, imax, jmax, kmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_UNSIGNED )
    return BndCommS4DEx( (unsigned*)array, nmax, imax, jmax, kmax, vc, vc_comm, procGrpNo );
  else if( dtype == MPI_UNSIGNED_LONG )
    return BndCommS4DEx( (unsigned long*)array, nmax, imax, jmax, kmax, vc, vc_comm, procGrpNo );
#ifdef MPI_LONG_LONG_INT
  else if( dtype == MPI_LONG_LONG_INT )
    return BndCommS4DEx( (long long int*)array, nmax, imax, jmax, kmax, vc, vc_comm, procGrpNo );
#endif
#ifdef MPI_LONG_LONG
  else if( dtype == MPI_LONG_LONG )
    return BndCommS4DEx( (long long*)array, nmax, imax, jmax, kmax, vc, vc_comm, procGrpNo );
#endif
#ifdef MPI_UNSIGNED_LONG_LONG
  else if( dtype == MPI_UNSIGNED_LONG_LONG )
    return BndCommS4DEx( (unsigned long long*)array, nmax, imax, jmax, kmax, vc, vc_comm, procGrpNo );
#endif

  return CPM_ERROR_MPI_INVALID_DATATYPE;
}

////////////////////////////////////////////////////////////////////////////////
// 非同期版袖通信(Vector3DEx版, MPI_Datatype指定)
cpm_ErrorCode
cpm_ParaManager::BndCommV3DEx_nowait( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                                    , int vc, int vc_comm, MPI_Request req[12], int procGrpNo ) 
{
  return BndCommS4DEx_nowait( dtype, array, 3, imax, jmax, kmax, vc, vc_comm, req, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 非同期版袖通信(Scalar4DEx版, MPI_Datatype指定)
cpm_ErrorCode
cpm_ParaManager::BndCommS4DEx_nowait( MPI_Datatype dtype, void *array, int nmax, int imax, int jmax, int kmax
                                    , int vc, int vc_comm, MPI_Request req[12], int procGrpNo ) 
{
  if( dtype == MPI_CHAR )
    return BndCommS4DEx_nowait( (char*)array, nmax, imax, jmax, kmax, vc, vc_comm, req, procGrpNo );
  else if( dtype == MPI_SHORT )
    return BndCommS4DEx_nowait( (short*)array, nmax, imax, jmax, kmax, vc, vc_comm, req, procGrpNo );
  else if( dtype == MPI_INT )
    return BndCommS4DEx_nowait( (int*)array, nmax, imax, jmax, kmax, vc, vc_comm, req, procGrpNo );
  else if( dtype == MPI_LONG )
    return BndCommS4DEx_nowait( (long*)array, nmax, imax, jmax, kmax, vc, vc_comm, req, procGrpNo );
  else if( dtype == MPI_FLOAT )
    return BndCommS4DEx_nowait( (float*)array, nmax, imax, jmax, kmax, vc, vc_comm, req, procGrpNo );
  else if( dtype == MPI_DOUBLE )
    return BndCommS4DEx_nowait( (double*)array, nmax, imax, jmax, kmax, vc, vc_comm, req, procGrpNo );
  else if( dtype == MPI_LONG_DOUBLE )
    return BndCommS4DEx_nowait( (long double*)array, nmax, imax, jmax, kmax, vc, vc_comm, req, procGrpNo );
  else if( dtype == MPI_UNSIGNED_CHAR )
    return BndCommS4DEx_nowait( (unsigned char*)array, nmax, imax, jmax, kmax, vc, vc_comm, req, procGrpNo );
  else if( dtype == MPI_UNSIGNED_SHORT )
    return BndCommS4DEx_nowait( (unsigned short*)array, nmax, imax, jmax, kmax, vc, vc_comm, req, procGrpNo );
  else if( dtype == MPI_UNSIGNED )
    return BndCommS4DEx_nowait( (unsigned*)array, nmax, imax, jmax, kmax, vc, vc_comm, req, procGrpNo );
  else if( dtype == MPI_UNSIGNED_LONG )
    return BndCommS4DEx_nowait( (unsigned long*)array, nmax, imax, jmax, kmax, vc, vc_comm, req, procGrpNo );
#ifdef MPI_LONG_LONG_INT
  else if( dtype == MPI_LONG_LONG_INT )
    return BndCommS4DEx_nowait( (long long int*)array, nmax, imax, jmax, kmax, vc, vc_comm, req, procGrpNo );
#endif
#ifdef MPI_LONG_LONG
  else if( dtype == MPI_LONG_LONG )
    return BndCommS4DEx_nowait( (long long*)array, nmax, imax, jmax, kmax, vc, vc_comm, req, procGrpNo );
#endif
#ifdef MPI_UNSIGNED_LONG_LONG
  else if( dtype == MPI_UNSIGNED_LONG_LONG )
    return BndCommS4DEx_nowait( (unsigned long long*)array, nmax, imax, jmax, kmax, vc, vc_comm, req, procGrpNo );
#endif

  return CPM_ERROR_MPI_INVALID_DATATYPE;
}

////////////////////////////////////////////////////////////////////////////////
// 非同期版袖通信のwait、展開(Vector3DEx版, MPI_Datatype指定)
cpm_ErrorCode
cpm_ParaManager::wait_BndCommV3DEx( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                                  , int vc, int vc_comm, MPI_Request req[12], int procGrpNo ) 
{
  return wait_BndCommS4DEx( dtype, array, 3, imax, jmax, kmax, vc, vc_comm, req, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 非同期版袖通信のwait、展開(Scalar4DEx版, MPI_Datatype指定)
cpm_ErrorCode
cpm_ParaManager::wait_BndCommS4DEx( MPI_Datatype dtype, void *array, int nmax, int imax, int jmax, int kmax
                                  , int vc, int vc_comm, MPI_Request req[12], int procGrpNo ) 
{
  if( dtype == MPI_CHAR )
    return wait_BndCommS4DEx( (char*)array, nmax, imax, jmax, kmax, vc, vc_comm, req, procGrpNo );
  else if( dtype == MPI_SHORT )
    return wait_BndCommS4DEx( (short*)array, nmax, imax, jmax, kmax, vc, vc_comm, req, procGrpNo );
  else if( dtype == MPI_INT )
    return wait_BndCommS4DEx( (int*)array, nmax, imax, jmax, kmax, vc, vc_comm, req, procGrpNo );
  else if( dtype == MPI_LONG )
    return wait_BndCommS4DEx( (long*)array, nmax, imax, jmax, kmax, vc, vc_comm, req, procGrpNo );
  else if( dtype == MPI_FLOAT )
    return wait_BndCommS4DEx( (float*)array, nmax, imax, jmax, kmax, vc, vc_comm, req, procGrpNo );
  else if( dtype == MPI_DOUBLE )
    return wait_BndCommS4DEx( (double*)array, nmax, imax, jmax, kmax, vc, vc_comm, req, procGrpNo );
  else if( dtype == MPI_LONG_DOUBLE )
    return wait_BndCommS4DEx( (long double*)array, nmax, imax, jmax, kmax, vc, vc_comm, req, procGrpNo );
  else if( dtype == MPI_UNSIGNED_CHAR )
    return wait_BndCommS4DEx( (unsigned char*)array, nmax, imax, jmax, kmax, vc, vc_comm, req, procGrpNo );
  else if( dtype == MPI_UNSIGNED_SHORT )
    return wait_BndCommS4DEx( (unsigned short*)array, nmax, imax, jmax, kmax, vc, vc_comm, req, procGrpNo );
  else if( dtype == MPI_UNSIGNED )
    return wait_BndCommS4DEx( (unsigned*)array, nmax, imax, jmax, kmax, vc, vc_comm, req, procGrpNo );
  else if( dtype == MPI_UNSIGNED_LONG )
    return wait_BndCommS4DEx( (unsigned long*)array, nmax, imax, jmax, kmax, vc, vc_comm, req, procGrpNo );
#ifdef MPI_LONG_LONG_INT
  else if( dtype == MPI_LONG_LONG_INT )
    return wait_BndCommS4DEx( (long long int*)array, nmax, imax, jmax, kmax, vc, vc_comm, req, procGrpNo );
#endif
#ifdef MPI_LONG_LONG
  else if( dtype == MPI_LONG_LONG )
    return wait_BndCommS4DEx( (long long*)array, nmax, imax, jmax, kmax, vc, vc_comm, req, procGrpNo );
#endif
#ifdef MPI_UNSIGNED_LONG_LONG
  else if( dtype == MPI_UNSIGNED_LONG_LONG )
    return wait_BndCommS4DEx( (unsigned long long*)array, nmax, imax, jmax, kmax, vc, vc_comm, req, procGrpNo );
#endif

  return CPM_ERROR_MPI_INVALID_DATATYPE;
}

////////////////////////////////////////////////////////////////////////////////
// 周期境界袖通信(Vector3DEx版, MPI_Datatype指定)
cpm_ErrorCode
cpm_ParaManager::PeriodicCommV3DEx( MPI_Datatype dtype, void *array, int imax, int jmax, int kmax
                                  , int vc, int vc_comm, cpm_DirFlag dir, cpm_PMFlag pm, int procGrpNo )
{
  return PeriodicCommS4DEx( dtype, array, 3, imax, jmax, kmax, vc, vc_comm, dir, pm, procGrpNo );
}

////////////////////////////////////////////////////////////////////////////////
// 周期境界袖通信(Scalar4DEx版, MPI_Datatype指定)
cpm_ErrorCode
cpm_ParaManager::PeriodicCommS4DEx( MPI_Datatype dtype, void *array, int nmax, int imax, int jmax, int kmax
                                  , int vc, int vc_comm, cpm_DirFlag dir, cpm_PMFlag pm, int procGrpNo )
{
  if( dtype == MPI_CHAR )
    return PeriodicCommS4DEx( (char*)array, nmax, imax, jmax, kmax, vc, vc_comm, dir, pm, procGrpNo );
  else if( dtype == MPI_SHORT )
    return PeriodicCommS4DEx( (short*)array, nmax, imax, jmax, kmax, vc, vc_comm, dir, pm, procGrpNo );
  else if( dtype == MPI_INT )
    return PeriodicCommS4DEx( (int*)array, nmax, imax, jmax, kmax, vc, vc_comm, dir, pm, procGrpNo );
  else if( dtype == MPI_LONG )
    return PeriodicCommS4DEx( (long*)array, nmax, imax, jmax, kmax, vc, vc_comm, dir, pm, procGrpNo );
  else if( dtype == MPI_FLOAT )
    return PeriodicCommS4DEx( (float*)array, nmax, imax, jmax, kmax, vc, vc_comm, dir, pm, procGrpNo );
  else if( dtype == MPI_DOUBLE )
    return PeriodicCommS4DEx( (double*)array, nmax, imax, jmax, kmax, vc, vc_comm, dir, pm, procGrpNo );
  else if( dtype == MPI_LONG_DOUBLE )
    return PeriodicCommS4DEx( (long double*)array, nmax, imax, jmax, kmax, vc, vc_comm, dir, pm, procGrpNo );
  else if( dtype == MPI_UNSIGNED_CHAR )
    return PeriodicCommS4DEx( (unsigned char*)array, nmax, imax, jmax, kmax, vc, vc_comm, dir, pm, procGrpNo );
  else if( dtype == MPI_UNSIGNED_SHORT )
    return PeriodicCommS4DEx( (unsigned short*)array, nmax, imax, jmax, kmax, vc, vc_comm, dir, pm, procGrpNo );
  else if( dtype == MPI_UNSIGNED )
    return PeriodicCommS4DEx( (unsigned*)array, nmax, imax, jmax, kmax, vc, vc_comm, dir, pm, procGrpNo );
   else if( dtype == MPI_UNSIGNED_LONG )
    return PeriodicCommS4DEx( (unsigned long*)array, nmax, imax, jmax, kmax, vc, vc_comm, dir, pm, procGrpNo );
#ifdef MPI_LONG_LONG_INT
  else if( dtype == MPI_LONG_LONG_INT )
    return PeriodicCommS4DEx( (long long int*)array, nmax, imax, jmax, kmax, vc, vc_comm, dir, pm, procGrpNo );
#endif
#ifdef MPI_LONG_LONG
  else if( dtype == MPI_LONG_LONG )
    return PeriodicCommS4DEx( (long long*)array, nmax, imax, jmax, kmax, vc, vc_comm, dir, pm, procGrpNo );
#endif
#ifdef MPI_UNSIGNED_LONG_LONG
  else if( dtype == MPI_UNSIGNED_LONG_LONG )
    return PeriodicCommS4DEx( (unsigned long long*)array, nmax, imax, jmax, kmax, vc, vc_comm, dir, pm, procGrpNo );
#endif

  return CPM_ERROR_MPI_INVALID_DATATYPE;
}

