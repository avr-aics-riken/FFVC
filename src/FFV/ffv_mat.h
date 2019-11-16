#ifndef _FFV_MAT_H_
#define _FFV_MAT_H_

/*
##################################################################################
#
# RIAM-COMPACT HPC version : RAinWATER
#
# Copyright (C) 2015-2018 Research Institute for Applied Mechanics(RIAM)
#                       / Research Institute for Information Technology(RIIT), Kyushu University.
# All rights reserved.
#
# Copyright (C) 2015-2018 Advanced Institute for Computational Science, RIKEN.
# All rights reserved.
#
##################################################################################
*/

/**
 * @file   ffv_mat.h
 * @brief  matA class Header
 * @author RIIT
 */

#include "ffv_Define.h"
#include "DomainInfo.h"
#include "ffv_Ffunc.h"
#include "mydebug.h"

using namespace std;


// #################################################################
class matA : public DomainInfo {

private:
  bool exist;           ///< 初回のみを保証するフラグ
  int szRcvBf;          ///< 受信バッファのサイズ
  int szSndBf;          ///< 送信バッファのサイズ
  REAL_TYPE* SndBfR;    ///< 送信バッファ
  REAL_TYPE* RcvBfR;    ///< 受信バッファ
  int* SndBfI;          ///< 送信バッファ
  int* RcvBfI;          ///< 受信バッファ
  int Maxlsz;           ///< ローカル配列の最大サイズ

  // ローカル
  int* l_bp;             ///< ビット行列
  REAL_TYPE* l_rhs;      ///< 右辺項

  // マスターのみ
  int* sd_sz;            ///< サブドメイン要素数
  int* sd_hd;            ///< サブドメイン先頭インデクス
  int* g_bp;             ///< ビット行列
  REAL_TYPE* g_rhs;      ///< 右辺項


public:
  /** コンストラクタ */
  matA() {
    exist = false;
    Maxlsz = 0;
    szSndBf = 0;
    szRcvBf = 0;

    SndBfI = NULL;
    RcvBfI = NULL;
    SndBfR = NULL;
    RcvBfR = NULL;

    sd_sz = NULL;
    sd_hd  = NULL;
    g_bp = NULL;
    g_rhs= NULL;

    l_bp = NULL;
    l_rhs= NULL;
  }


  /** デストラクタ */
  ~matA() { }



public:

  bool prep(double &memory);

  bool outArrayInt(int* src, const REAL_TYPE cs, char* fname);

  bool outArrayReal(REAL_TYPE* src, char* fname);
  
  void existOn() {
    exist = true;
  }
  
  bool isExist() {
    return exist;
  }


private:
  // #################################################################
  /**
   * @brief 配列の集約
   */
  template <typename T>
  bool gatherArray(T* sbuf,
                         T* src,
                         T* rbuf,
                         T* dst,
                         MPI_Datatype d_type,
                         int sw=0)
  {
    // 送信バッファにデータをコピー
    pack2SndBuf(sbuf, src);

    if ( MPI_SUCCESS != MPI_Gather(sbuf,
                                   szSndBf,
                                   d_type,
                                   rbuf,
                                   szSndBf,
                                   d_type,
                                   0,
                                   MPI_COMM_WORLD) ) return false;


    // バッファからグローバルにデータをコピー
    Hostonly_ copyRcvBuf2gl(rbuf, dst, sw);

    return true;
  }



  // #################################################################
  // @brief スカラーデータを送信バッファにコピー
  template <typename T>
  void pack2SndBuf(T* sbuf, T* src)
  {
    #pragma omp parallel for collapse(2)
    for (int k=0; k<size[2]; k++) {
      for (int j=0; j<size[1]; j++) {
        for (int i=0; i<size[0]; i++) {
          sbuf[_IDX_S3D(i,j,k,size[0],size[1],size[2],0)]
         = src[_IDX_S3D(i,j,k,size[0],size[1],size[2],2)];
        }
      }
    }
  }


  // #################################################################
  // @brief 受信バッファの内容をグローバル配列にコピー
  // @note ランク０のみ
  template <typename T>
  void copyRcvBuf2gl(T* rbuf, T* dst, int sw)
  {

    for (int n=0; n<numProc; n++)
    {
      int* msz = &sd_sz[n*3];
      int* mhd = &sd_hd[n*3];
      int ofst = Maxlsz*n;
      int u,v,w;

      #pragma omp parallel for collapse(2) private(u,v,w)
      for (int k=0; k<msz[2]; k++) {
        for (int j=0; j<msz[1]; j++) {
          for (int i=0; i<msz[0]; i++) {
            u = mhd[0] - 1 + i;
            v = mhd[1] - 1 + j;
            w = mhd[2] - 1 + k;
            dst[w*G_size[0]*G_size[1] + v*G_size[0] + u]
               = rbuf[ofst + k*msz[0]*msz[1] + j*msz[0] + i];
          }
        }
      }
    } // numProc

    // Debug
    if (sw==1)
    {
      for (int k=0; k<G_size[2]; k++) {
        for (int j=0; j<G_size[1]; j++) {
          for (int i=0; i<G_size[0]; i++) {
            T s = dst[k*G_size[0]*G_size[1] + j*G_size[0] + i];
            if ( s != 0.0 ) printf("%d %d %d %e\n", i,j,k,(double)s);
          }
        }
      }
    }

  }

};

#endif // _FFV_MAT_H_
