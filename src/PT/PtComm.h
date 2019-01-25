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
public:
	int pInfo[NDST*2];         ///< 粒子+情報保持
														 ///< pInfo[27][m]
														 ///<   m=0 送信粒子個数
														 ///<   m=1 受信粒子個数
	vector<REAL_TYPE> ps_buf;  ///< 粒子用送信バッファ
	vector<REAL_TYPE> pr_buf;  ///< 粒子用受信バッファ
	vector<int> bs_buf;        ///< bf用送信バッファ
	vector<int> br_buf;        ///< bf用受信バッファ
	int is_buf[NDST];          ///< 情報送信用バッファ
	int ir_buf[NDST];          ///< 情報受信用バッファ
	unsigned* nPart;           ///< 全ランクの粒子数
	
	
protected:
  int neighbor[NDST];        ///< 隣接ランク番号
  int G_div[3];              ///< 領域分割数
  int myRank;                ///< 自ランク番号
  int numProc;               ///< プロセス数
  int Rmap[NDST];            ///< 3x3x3のランクマップ
  int Cmap[NDST*2];          ///< 送受信対象 > 1
	

#ifndef DISABLE_MPI
  MPI_Datatype data_type;    ///< float / doubleの切り替え
  MPI_Request reqR[2*NDST];
  MPI_Request reqI[2*NDST];
#endif


public:
  // デフォルト コンストラクタ
  PtComm() {
    myRank = -1;

    for (int i=0; i<3; i++) {
      G_div[i] = 0;
    }

    ps_buf = vector<REAL_TYPE>(6*BUF_UNIT*NDST, 0);
    pr_buf = vector<REAL_TYPE>(6*BUF_UNIT*NDST, 0);
    bs_buf = vector<int>(2*BUF_UNIT*NDST, 0);
    br_buf = vector<int>(2*BUF_UNIT*NDST, 0);

    memset(neighbor, -1, sizeof(int)*NDST);
    memset(pInfo,  0, sizeof(int)*NDST*2);
    memset(is_buf, 0, sizeof(int)*NDST);
    memset(ir_buf, 0, sizeof(int)*NDST);
    memset(Rmap, -1, sizeof(int)*NDST);
    memset(Cmap, 0, sizeof(int)*NDST*2);


#ifndef DISABLE_MPI
    data_type = MPI_FLOAT;
    if( sizeof(REAL_TYPE) == 8 )
    {
      data_type = MPI_DOUBLE;
    }

    for (int i=0; i<2*NDST; i++) {
      reqR[i] = MPI_REQUEST_NULL;
      reqI[i] = MPI_REQUEST_NULL;
    }
#endif
    
    nPart = NULL;
  }

  /** デストラクタ */
  ~PtComm() {
    if (nPart) delete [] nPart;
  }

  

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
  
  unsigned* nPart_ptr() {
    return nPart;
  }

  
  // @note Rmapのポインタを返す
  int* setPtComm(const int m_div[],
                 const int m_rank,
                 const int m_np)
  {
    this->G_div[0]    = m_div[0];
    this->G_div[1]    = m_div[1];
    this->G_div[2]    = m_div[2];
    this->myRank      = m_rank;
    this->numProc     = m_np;

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
    
    nPart = new unsigned[numProc];
    memset(nPart, 0, sizeof(unsigned)*numProc);
    
    return Rmap;
  }


  // @brief 送信バッファのリサイズと再初期化
  void resizeSendBuffer(const unsigned len) {
    ps_buf.clear();
    ps_buf.resize(len*6*NDST); // pos+vel=6要素　初期値がフィルされる
    bs_buf.clear();
    bs_buf.resize(len*2*NDST); // b + foo = 2要素
  }


  // @brief 送信バッファの再初期化
  void initSendBuffer() {
    fill(ps_buf.begin(), ps_buf.end(), 0);
    fill(bs_buf.begin(), bs_buf.end(), 0);
  }


  // @brief 経路確定のための通信
  bool establishCommPath();


  // @brief 粒子データの通信
  bool commParticle(const unsigned buf_length);


  // @brief 統計
  bool Statistics(int& nCommP,
                  const unsigned l_part,
                  unsigned& g_part);
  
  
  // @brief 引数のご全ランクの最大数を求める
  bool getMax(unsigned& var);
  
	
	// @brief 引数の和
	bool getSum(unsigned& var);
  
	
	
	
private:

  // @brief 通信経路確定のための事前情報の通信
  bool commInfo(int* sbuf,
                int* rbuf,
                const int msz,
                const int nID,
                MPI_Request req[2]);


  // @brief 粒子データの送受信(REAL_TYPE)
  bool SendRecvParticle(const int s_msg,
                        const int r_msg,
                        const int nID,
                        REAL_TYPE* sbuf,
                        REAL_TYPE* rbuf,
                        MPI_Request req[2],
                        const int flag[2]);
  
  
  // @brief 粒子データの送受信(int)
  bool SendRecvParticle(const int s_msg,
                        const int r_msg,
                        const int nID,
                        int* sbuf,
                        int* rbuf,
                        MPI_Request req[2],
                        const int flag[2]);

  
  // @brief 粒子情報通信の確定
  bool waitCommInfo(MPI_Request* req,
                    MPI_Status* stat);
  
  
  // @brief 粒子データ通信の確定 (REAL_TYPE)
  bool waitCommPart(MPI_Request* req,
                    MPI_Status* stat,
                    const int flag[2],
                    const int r_msg,
                    REAL_TYPE* rbuf);
  
  
  // @brief 粒子データ通信の確定 (int)
  bool waitCommPart(MPI_Request* req,
                    MPI_Status* stat,
                    const int flag[2],
                    const int r_msg,
                    int* rbuf);

};

#endif // _PT_COMM_H_
