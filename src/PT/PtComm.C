
//##################################################################################
//
// Copyright (c) 2019 Research Institute for Information Technology, Kyushu university
// All rights researved.
//
//##################################################################################

/*
 * @file   PtComm.C
 * @brief  PtComm class
 */

#include "PtComm.h"


//#############################################################################
// @brief 経路確定のための通信
// @note 送受信情報を周囲の26方向と通信する
bool PtComm::establishCommPath()
{
  // Communication identifier
  for (int i=0; i<2*NDST; i++)
       reqI[i] = MPI_REQUEST_NULL;

  memset(is_buf, 0, sizeof(int)*NDST);
  memset(ir_buf, 0, sizeof(int)*NDST);


  for (int i=0; i<NDST; i++)
  {
    is_buf[i] = pInfo[2*i];  // 送信粒子数
  }


  for (int i=0; i<NDST; i++)
  {
    int nID = Rmap[i];
    
    if ( nID >= 0 )
    {
      if ( i != myRank ) {
        if ( !commInfo(&is_buf[i],
                       &ir_buf[i],
                       1,
                       nID,
                       &reqI[2*i]) ) return false;
      }
    }
  }


  MPI_Status statI[NDST*2];

  for (int i=0; i<NDST; i++)
  {
    int nID = Rmap[i];
    
    if ( nID >= 0 )
    {
      if ( i != myRank ) {
        if ( !waitCommInfo( &reqI[2*i],
                           &statI[2*i]) ) return false;
      }
    }
  }


  // 受信粒子数
  for (int i=0; i<NDST; i++)
  {
    pInfo[2*i+1] = ir_buf[i];
  }
  
  
  // 経路の確定 送受信があるランク方向のみ
  memset(Cmap, 0, sizeof(int)*NDST*2);
  
  for (int i=0; i<NDST; i++)
  {
    int nID = Rmap[i];
    
    // 計算領域内のみ
    if ( nID >= 0 )
    {
      // 送信フラグ
      if ( pInfo[2*i+0] > 0 ) Cmap[2*i+0] = 1;
      
      // 受信フラグ
      if ( pInfo[2*i+1] > 0 ) Cmap[2*i+1] = 1;
    }
  }
  
  return true;
}



//#############################################################################
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
  int tag_r=0;
  int tag_s=0;
  
  // 近傍からの受信
  if ( MPI_SUCCESS != MPI_Irecv(rbuf,
                                msz,
                                MPI_INT,
                                nID,
                                tag_r,
                                MPI_COMM_WORLD,
                                &req[1]) ) return false;
  
  // 近傍への送信
  if ( MPI_SUCCESS != MPI_Isend(sbuf,
                                msz,
                                MPI_INT,
                                nID,
                                tag_s,
                                MPI_COMM_WORLD,
                                &req[0]) ) return false;
  
  return true;
}


/*#############################################################################
// @brief 送受信データ型を定義する
// @note
bool PtComm::defineDatatypes(const int i, const unsigned buf_length)
{
  int blockcount[2];
  MPI_Datatype types[2];
  MPI_Aint     displs[2];
  MPI_Datatype packedSend;
  
  blockcount[0] = buf_length * 6; // pos + vel
  blockcount[1] = buf_length * 2; // bf + foo
  
  types[0] = data_type;
  types[1] = MPI_INT;
  
  MPI_Get_address(&ps_buf[buf_length*i], &displs[0]);
  MPI_Get_address(&bs_buf[buf_length*i], &displs[1]);
}
*/

//#############################################################################
// @brief 粒子データの通信
// @param [in] buf_max_particle バッファ計算の元になる最大粒子数
bool PtComm::commParticle(const unsigned buf_max_particle)
{
  for (int i=0; i<2*NDST; i++) {
    reqR[i] = MPI_REQUEST_NULL;
    reqI[i] = MPI_REQUEST_NULL;
  }
  
  int ofst6 = buf_max_particle * 6;  // pos, velの6要素
  int ofst2 = buf_max_particle * 2;  // bf, fooの2要素
  
  for (int i=0; i<NDST; i++)
  {
    int nID = Rmap[i];
    int ms = pInfo[2*i+0];
    int mr = pInfo[2*i+1];
    
    if ( i != myRank )
    {
      if ( !SendRecvParticle(ms * 6,  // 送信数
                             mr * 6,  // 受信数
                             nID,
                             &ps_buf[ofst6 * i],
                             &pr_buf[ofst6 * i],
                             &reqR[2*i],
                             &Cmap[2*i]) ) return false;
      
      if ( !SendRecvParticle(ms * 2,
                             mr * 2,
                             nID,
                             &bs_buf[ofst2 * i],
                             &br_buf[ofst2 * i],
                             &reqI[2*i],
                             &Cmap[2*i]) ) return false;
    }
  }


  MPI_Status statR[NDST*2];
  MPI_Status statI[NDST*2];

  for (int i=0; i<NDST; i++)
  {
    int mr = pInfo[2*i+1];
    
    if ( i != myRank )
    {
      if ( !waitCommPart( &reqR[2*i],
                         &statR[2*i],
                          &Cmap[2*i],
                         mr * 6,
                         &pr_buf[ofst6 * i]) ) return false;
      
      if ( !waitCommPart( &reqI[2*i],
                         &statI[2*i],
                          &Cmap[2*i],
                         mr * 2,
                         &br_buf[ofst2 * i]) ) return false;
    }
  }


  return true;
}



//#############################################################################
/*
 * @brief 粒子データの送受信(REAL_TYPE)
 * @param [in]  s_msg 送信カウント
 * @param [in]  r_msg 受信カウント
 * @param [in]  nID   近傍のランク番号
 * @param [in]  sbuf  Send buffer to neighbor
 * @param [out] rbuf  Recieve buffer from neighbor
 * @param [out] req   Array of MPI request
 * @param [in]  flag  送受信実行フラグ
 * @retval true-success, false-fail
 */
bool PtComm::SendRecvParticle(const int s_msg,
                              const int r_msg,
                              const int nID,
                              REAL_TYPE* sbuf,
                              REAL_TYPE* rbuf,
                              MPI_Request req[2],
                              const int flag[2])
{
  int tag_r=0;
  int tag_s=0;
  

  // 近傍からの受信
  if ( flag[1] == 1 ) {
    if ( MPI_SUCCESS != MPI_Irecv(rbuf,
                                  r_msg,
                                  data_type,
                                  nID,
                                  tag_r,
                                  MPI_COMM_WORLD,
                                  &req[1]) ) return false;
  }

  // 近傍への送信
  if ( flag[0] == 1 ) {
    if ( MPI_SUCCESS != MPI_Isend(sbuf,
                                  s_msg,
                                  data_type,
                                  nID,
                                  tag_s,
                                  MPI_COMM_WORLD,
                                  &req[0]) ) return false;
#ifdef PT_DEBUG
    // 送信内容を表示、受信内容はwaitAll後
    //printf("snd[%d] : ", myRank);
    //for (int i=0; i<s_msg; i++) printf("%f ", sbuf[i]);
    //printf("\n");
#endif
  }

  return true;
}


//#############################################################################
/*
 * @brief 粒子データの送受信(int)
 * @param [in]  s_msg 送信カウント
 * @param [in]  r_msg 受信カウント
 * @param [in]  nID   近傍のランク番号
 * @param [in]  sbuf  Send buffer to neighbor
 * @param [out] rbuf  Recieve buffer from neighbor
 * @param [out] req   Array of MPI request
 * @param [in]  flag  送受信実行フラグ
 * @retval true-success, false-fail
 */
bool PtComm::SendRecvParticle(const int s_msg,
                              const int r_msg,
                              const int nID,
                              int* sbuf,
                              int* rbuf,
                              MPI_Request req[2],
                              const int flag[2])
{
  int tag_r=0;
  int tag_s=0;
  
  // 近傍からの受信
  if ( flag[1] == 1 ) {
    if ( MPI_SUCCESS != MPI_Irecv(rbuf,
                                  r_msg ,
                                  MPI_INT,
                                  nID,
                                  tag_r,
                                  MPI_COMM_WORLD,
                                  &req[1]) ) return false;
  }
  
  // 近傍への送信
  if ( flag[0] == 1 ) {
    if ( MPI_SUCCESS != MPI_Isend(sbuf,
                                  s_msg,
                                  MPI_INT,
                                  nID,
                                  tag_s,
                                  MPI_COMM_WORLD,
                                  &req[0]) ) return false;
#ifdef PT_DEBUG
    // 送信内容を表示、受信内容はwaitAll後
    //printf("snd[%d] : ", myRank);
    //for (int i=0; i<s_msg; i++) printf("%d ", sbuf[i]);
    //printf("\n");
#endif
  }
  
  return true;
}


//#############################################################################
// @brief 粒子情報通信の確定
bool PtComm::waitCommInfo(MPI_Request* req,
                          MPI_Status* stat)
{
  if ( MPI_SUCCESS != MPI_Waitall(2, req, stat) ) return false;

  return true;
}


//#############################################################################
// @brief 粒子データ通信の確定 (REAL_TYPE)
bool PtComm::waitCommPart(MPI_Request* req,
                          MPI_Status* stat,
                          const int flag[2],
                          const int r_msg,
                          REAL_TYPE* rbuf)
{
  // 近傍への送信
  if ( flag[0] == 1 )
  {
    if ( MPI_SUCCESS != MPI_Waitall(1, &req[0], &stat[0]) ) return false;
  }
  // 近傍からの受信
  if ( flag[1] == 1 )
  {
    if ( MPI_SUCCESS != MPI_Waitall(1, &req[1], &stat[1]) ) return false;
#ifdef PT_DEBUG
    //printf("rcv[%d] : ", myRank);
    //for (int i=0; i<r_msg; i++) printf("%f ", rbuf[i]);
    //printf("\n");
#endif
  }
  
  return true;
}


//#############################################################################
// @brief 粒子データ通信の確定 (int)
bool PtComm::waitCommPart(MPI_Request* req,
                          MPI_Status* stat,
                          const int flag[2],
                          const int r_msg,
                          int* rbuf)
{
  // 近傍への送信
  if ( flag[0] == 1 )
  {
    if ( MPI_SUCCESS != MPI_Waitall(1, &req[0], &stat[0]) ) return false;
  }
  // 近傍からの受信
  if ( flag[1] == 1 )
  {
    if ( MPI_SUCCESS != MPI_Waitall(1, &req[1], &stat[1]) ) return false;
#ifdef PT_DEBUG
    //printf("rcv[%d] : ", myRank);
    //for (int i=0; i<r_msg; i++) printf("%d ", rbuf[i]);
    //printf("\n");
#endif
  }
  
  return true;
}


//#############################################################################
// @brief 統計
// @param [out] nCommP   送受信する総粒子数（ランク0のみ）
// @param [in]  l_part   自ランクのもつ全粒子数
// @param [out] g_part   全粒子数（ランク0のみ）
bool PtComm::Statistics(int& nCommP,
                        const unsigned l_part,
                        unsigned& g_part)
{
  int sum = 0;
  for (int i=0; i<NDST; i++)
  {
    sum += pInfo[2*i+0] + pInfo[2*i+1]; // 送受信数
  }
  
  // 総受信数の総和
  int tmp = sum;
  if ( MPI_SUCCESS != MPI_Reduce(&tmp,
                                 &sum,
                                 1,
                                 MPI_INT,
                                 MPI_SUM,
                                 0,
                                 MPI_COMM_WORLD) ) return false;
  nCommP = sum / 2;
  
  unsigned tmp2 = l_part;
  if ( MPI_SUCCESS != MPI_Gather(&tmp2,
                                 1,
                                 MPI_UNSIGNED,
                                 nPart,
                                 1,
                                 MPI_UNSIGNED,
                                 0,
                                 MPI_COMM_WORLD) ) return false;
  
  unsigned sum2 = 0;
  for (int i=0; i<numProc; i++) {
    sum2 += nPart[i];
  }
  g_part = sum2;
  
  return true;
}


//#############################################################################
// @brief 最大数を求める
// @param [in,out]  npart  バッファ計算用の最大粒子数
bool PtComm::getMax(unsigned& npart)
{
  // 最大バッファ要素数
  unsigned tmp = npart;
  if ( MPI_SUCCESS != MPI_Allreduce(&tmp,
                                    &npart,
                                    1,
                                    MPI_UNSIGNED,
                                    MPI_MAX,
                                    MPI_COMM_WORLD) ) return false;
  
  return true;
}


//#############################################################################
// @brief 引数の和
// @param [in,out]  var  和をトル変数
bool PtComm::getSum(unsigned& var)
{
	// 最大バッファ要素数
	unsigned tmp = var;
	if ( MPI_SUCCESS != MPI_Allreduce(&tmp,
																		&var,
																		1,
																		MPI_UNSIGNED,
																		MPI_SUM,
																		MPI_COMM_WORLD) ) return false;
	
	return true;
}


//#############################################################################
/*
 * @brief 粒子座標のbroadcast
 * @param [in]  msg  送信カウント
 * @param [in]  buf  Send/Recv buffer
 * @retval true-success, false-fail
 */
bool PtComm::BcastParticles(const int msg, REAL_TYPE* buf)
{
	if (MPI_Bcast(buf, msg, data_type, 0, MPI_COMM_WORLD) != MPI_SUCCESS) {
		return false;
	}
	
	return true;
}

