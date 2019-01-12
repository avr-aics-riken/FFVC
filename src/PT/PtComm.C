/*
###################################################################################
#
# CBrick
#
# Copyright (c) 2017-2018 Research Institute for Information Technology(RIIT),
#                    Kyushu University.  All rights reserved.
#
####################################################################################
*/

/*
 * @file   PtComm.C
 * @brief  PtComm class
 */

#include "PtComm.h"


// @brief 経路確定のための通信
bool PtComm::establishCommPath()
{
  // Communication identifier
  for (int i=0; i<2*NDST; i++)
       mpi_req[i] = MPI_REQUEST_NULL;

  memset(is_buf, 0, sizeof(int)*NDST*3);
  memset(ir_buf, 0, sizeof(int)*NDST*3);


  for (int i=0; i<NDST; i++)
  {
    is_buf[3*i+0] = pInfo[4*i+0];  // 送信粒子数
    is_buf[3*i+1] = pInfo[4*i+1];  // グループID
    is_buf[3*i+2] = pInfo[4*i+2];  // 粒子ID
  }


  for (int i=0; i<NDST; i++)
  {
    if ( !commInfo(&is_buf[3*i],
                   &ir_buf[3*i],
                   3,
                   Rmap[i],
                   &mpi_req[2*i]) ) return false;
  }


  MPI_Status mpi_stat[NDST*2];

  for (int i=0; i<NDST; i++)
  {
    if ( !waitComm( &mpi_req[2*i],
                    &mpi_stat[2*i]) ) return false;
  }


  // グループIDと粒子IDは上書き
  for (int i=0; i<NDST; i++)
  {
    pInfo[4*i+1] = ir_buf[3*i+1];  // グループID
    pInfo[4*i+2] = ir_buf[3*i+2];  // 粒子ID
    pInfo[4*i+3] = ir_buf[3*i+0];  // 受信粒子数
  }

  return true;
}


// @brief 粒子データの通信
bool PtComm::commParticle()
{
  // Communication identifier
  for (int i=0; i<2*NDST; i++)
       mpi_req[i] = MPI_REQUEST_NULL;


  for (int i=0; i<NDST; i++)
  {
    if ( !SendRecvParticle(pInfo[4*i+0],
                           pInfo[4*i+3],
                           Rmap[i],
                           &is_buf[3*i],
                           &ir_buf[3*i],
                           &mpi_req[2*i]) ) return false;
  }


  MPI_Status mpi_stat[NDST*2];

  for (int i=0; i<NDST; i++)
  {
    if ( !waitComm( &mpi_req[2*i],
                   &mpi_stat[2*i]) ) return false;
  }


  return true;
}


/*
 * @brief 粒子データの送受信
 * @param [in]  s_msg 送信カウント
 * @param [in]  r_msg 受信カウント
 * @param [in]  nID   近傍のランク番号
 * @param [in]  sbuf  Send buffer to neighbor
 * @param [out] rbuf  Recieve buffer from neighbor
 * @param [out] req   Array of MPI request
 * @retval true-success, false-fail
 */
bool PtComm::SendRecvParticle(const int s_msg,
                              const int r_msg,
                              const int nID,
                              int* sbuf,
                              int* rbuf,
                              MPI_Request req[2])
{
  // Identifier
  MPI_Request r0 = MPI_REQUEST_NULL;
  MPI_Request r1 = MPI_REQUEST_NULL;

  int tag_r=0;
  int tag_s=0;

  // 近傍からの受信
  if ( r_msg > 0 ) {
    if ( MPI_SUCCESS != MPI_Irecv(rbuf,
                                  r_msg ,
                                  data_type,
                                  nID,
                                  tag_r,
                                  MPI_COMM_WORLD,
                                  &r1) ) return false;
  }

  // 近傍への送信
  if ( s_msg > 0 ) {
    if ( MPI_SUCCESS != MPI_Isend(sbuf,
                                  s_msg,
                                  data_type,
                                  nID,
                                  tag_s,
                                  MPI_COMM_WORLD,
                                  &r0) ) return false;
  }

  req[0] = r0; // Send
  req[1] = r1; // Recieve

  return true;
}


/*
 * @brief 通信経路確定のための事前情報の通信
 * @param [in]  sbuf Send buffer to neighbor
 * @param [out] rbuf Recieve buffer from neighbor
 * @param [in]  msz  send message size
 * @param [in]  nID  Neighbor ID to send/recv
 * @param [out] req  Array of MPI request
 * @retval true-success, false-fail
 */
bool PtComm::commInfo(int* sbuf,
                      int* rbuf,
                      const int msz,
                      const int nID,
                      MPI_Request req[2])
{
  // Identifier
  MPI_Request r0 = MPI_REQUEST_NULL;
  MPI_Request r1 = MPI_REQUEST_NULL;

  int tag_r=0;
  int tag_s=0;

  if ( nID >= 0 ) {

    // 近傍からの受信
    if ( MPI_SUCCESS != MPI_Irecv(rbuf,
                                  msz,
                                  data_type,
                                  nID,
                                  tag_r,
                                  MPI_COMM_WORLD,
                                  &r1) ) return false;

    // 近傍への送信
    if ( MPI_SUCCESS != MPI_Isend(sbuf,
                                  msz,
                                  data_type,
                                  nID,
                                  tag_s,
                                  MPI_COMM_WORLD,
                                  &r0) ) return false;
  }

  req[0] = r0; // Send
  req[1] = r1; // Recieve

  return true;
}


// @brief 通信の確定
bool PtComm::waitComm(MPI_Request* req,
                      MPI_Status* stat)
{

  if ( MPI_SUCCESS != MPI_Waitall(2, req, stat) ) return false;

  return true;
}
