#ifndef _PT_COMM_H_
#define _PT_COMM_H_

/*
//##################################################################################
//
// Copyright (c) 2019 Research Institute for Information Technology, Kyushu university
// All rights researved.
//
//##################################################################################
*/

/**
 * @file   PtComm.h
 * @brief  PtComm class Header
 */

#include "mpi.h"
#include <string.h>
#include <stdlib.h>
#include "FB_Define.h"
#include "mydebug.h"
#include "PtDefine.h"
#include <vector>

using std::vector;
using std::string;
using std::fill;


class PtComm {
private:
  int neighbor[NDST];        ///< 隣接ランク番号
  int pInfo[NDST*4];         ///< 粒子+情報保持
                             ///< pInfo[27][m]
                             ///<   m=0 送信粒子個数
                             ///<   m=1 粒子グループ
                             ///<   m-2 粒子ID
                             ///<   m=3 受信粒子個数
  vector<REAL_TYPE> ps_buf;  ///< 粒子用送信バッファ
  vector<REAL_TYPE> pr_buf;  ///< 粒子用受信バッファ
  vector<int> bs_buf;        ///< bf用送信バッファ
  vector<int> br_buf;        ///< bf用受信バッファ
  int is_buf[NDST*3];        ///< 情報送信用バッファ
  int ir_buf[NDST*3];        ///< 情報受信用バッファ

  int G_div[3];              ///< 領域分割数
  int myRank;                ///< 自ランク番号
  int Rmap[NDST];            ///< 3x3x3のランクマップ

#ifndef DISABLE_MPI
  MPI_Datatype data_type;    ///< float / doubleの切り替え
  MPI_Request mpi_req[2*NDST];
#endif


public:
  // デフォルト コンストラクタ
  PtComm() {
    myRank = -1;

    for (int i=0; i<3; i++) {
      G_div[i] = 0;
    }

    for (int i=0; i<NDST; i++) Rmap[i] = -1;

    ps_buf = vector<REAL_TYPE>(3*BUF_UNIT*NDST, 0);
    pr_buf = vector<REAL_TYPE>(3*BUF_UNIT*NDST, 0);
    bs_buf = vector<int>(BUF_UNIT*NDST, 0);
    br_buf = vector<int>(BUF_UNIT*NDST, 0);

    memset(neighbor, -1, sizeof(int)*NDST);
    memset(pInfo, 0, sizeof(unsigned)*NDST*2);

#ifndef DISABLE_MPI
    data_type = MPI_FLOAT;
    if( sizeof(REAL_TYPE) == 8 )
    {
      data_type = MPI_DOUBLE;
    }

    for (int i=0; i<2*NDST; i++)
       mpi_req[i] = MPI_REQUEST_NULL;
#endif
  }

  /** デストラクタ */
  ~PtComm() {}


public:

  REAL_TYPE* ps_ptr() {
    return &ps_buf[0];
  }

  REAL_TYPE* pr_ptr() {
    return &pr_buf[0];
  }

  int* bs_ptr() {
    return &bs_buf[0];
  }

  int* br_ptr() {
    return &br_buf[0];
  }

  int* pInfo_ptr() {
    return pInfo;
  }

  void setPtComm(const int m_div[],
                 const int m_rank)
  {
    this->G_div[0]    = m_div[0];
    this->G_div[1]    = m_div[1];
    this->G_div[2]    = m_div[2];
    this->myRank      = m_rank;

    int ib = G_div[0];
    int jb = G_div[1];
    int kb = G_div[2];

    // 自ランクのブロック番号を得る
    int w = myRank / (ib*jb);
    int v = (myRank - w*ib*jb) / ib;
    int u = myRank - w*ib*jb - v*ib;

    for (int k=w-1; k<w+2; k++) {
    for (int j=v-1; j<v+2; j++) {
    for (int i=u-1; i<u+2; i++) {
      if ( !(k<0 || k>kb-1
          || j<0 || j>jb-1
          || i<0 || i>ib-1) ) {
        int m = (k-w+1)*9 + (j-v+1)*3 + i-u+1;
        Rmap[m] = k*ib*jb + j*ib + i;
        //printf("%d %d %d : m=%d : rank=%d\n",i,j,k,m,Rmap[m]);
      }
    }}}
  }


  // @brief 送信バッファのリサイズと再初期化
  void resizeSendBuffer(const unsigned len) {
    ps_buf.clear();
    ps_buf.resize(len*3*NDST); // 初期値がフィルされる
    bs_buf.clear();
    bs_buf.resize(len*NDST);
  }


  // @brief 送信バッファの再初期化
  void initSendBuffer() {
    fill(ps_buf.begin(), ps_buf.end(), 0);
    fill(bs_buf.begin(), bs_buf.end(), 0);
  }


  // @brief 経路確定のための通信
  bool establishCommPath();


  // @brief 粒子データの通信
  bool commParticle();


private:

  int min1(int a, int b) const {
    return a > b ? b : a;
  }


  int max1(int a, int b) const {
    return a > b ? a : b;
  }


  // @brief 通信経路確定のための事前情報の通信
  bool commInfo(int* sbuf,
                int* rbuf,
                const int msz,
                const int nID,
                MPI_Request req[2]);


  // @brief 粒子データの送受信
  bool SendRecvParticle(const int s_msg,
                        const int r_msg,
                        const int nID,
                        int* sbuf,
                        int* rbuf,
                        MPI_Request req[2]);


  // @brief 通信の確定
  bool waitComm(MPI_Request* req,
                MPI_Status* stat);

};

#endif // _PT_COMM_H_
