/*
###################################################################################
#
# RIAM-COMPACT HPC version : RAinWATER
#
# Copyright (C) 2015-2018 Research Institute for Applied Mechanics(RIAM)
#                       / Research Institute for Information Technology(RIIT), Kyushu University.
# All rights reserved.
#
# Copyright (C) 2015-2018 Advanced Institute for Computational Science, RIKEN.
# All rights reserved.
//
// Copyright (c) 2016-2020 Research Institute for Information Research, Kyushu University.
// All rights reserved.
#
###################################################################################
*/

#include "ffv_mat.h"
#include "ffv_LSfunc.h"

// #################################################################
bool matA::outArrayInt(int* src, const REAL_TYPE cs, char* fname)
{
  if (!gatherArray(SndBfI, src, RcvBfI, g_bp, MPI_INT)) return false;

  Hostonly_
  {
    REAL_TYPE cm = cs;
    write_mat_(G_size, pitch, g_bp, &cm, fname);
  }

  return true;
}


// #################################################################
bool matA::outArrayReal(REAL_TYPE* src, char* fname)
{
  MPI_Datatype data_type;
  data_type = MPI_FLOAT;
  if( sizeof(REAL_TYPE) == 8 )
  {
    data_type = MPI_DOUBLE;
  }

  if (!gatherArray(SndBfR, src, RcvBfR, g_rhs, data_type, 0)) return false;

  Hostonly_
  {
    write_rhs_(G_size, g_rhs, g_bp, fname);
  }

  return true;
}

// #################################################################
/**
 * @brief Gather
 */
bool matA::prep(double &memory)
{
  // グローバル配列 マスターのみ
  if (myRank == 0)
  {
    if ( !(sd_sz = new int[numProc*3]) ) return false;
    memory += (double)numProc * (double)sizeof(int) * 3.0;

    if ( !(sd_hd = new int[numProc*3]) ) return false;
    memory += (double)numProc * (double)sizeof(int) * 3.0;

    int nx = G_size[0] * G_size[1] * G_size[2];
    if ( !(g_bp = new int[nx]) ) return false;
    memory += (double)nx * (double)sizeof(int);

    if ( !(g_rhs = new REAL_TYPE[nx]) ) return false;
    memory += (double)nx * (double)sizeof(REAL_TYPE);
  }


  // 領域情報の収集
  if ( numProc > 1 )
  {
    if ( MPI_SUCCESS != MPI_Gather(size,
                                   3,
                                   MPI_INT,
                                   sd_sz,
                                   3,
                                   MPI_INT,
                                   0,
                                   MPI_COMM_WORLD) ) return false;

    if ( MPI_SUCCESS != MPI_Gather(head,
                                   3,
                                   MPI_INT,
                                   sd_hd,
                                   3,
                                   MPI_INT,
                                   0,
                                   MPI_COMM_WORLD) ) return false;

  }
  else
  {
    sd_sz = size;
    sd_hd = head;
  }


  // バッファの最大値を求める
  Maxlsz = 0;

  Hostonly_
  {
    for (int m=0; m<numProc; m++)
    {
      int* msz = &sd_sz[m*3];
      int tmp = msz[0] * msz[1] * msz[2];
      if (Maxlsz < tmp) Maxlsz = tmp;
    }
  }



  if ( MPI_SUCCESS != MPI_Bcast(&Maxlsz,
                                1,
                                MPI_INT,
                                0,
                                MPI_COMM_WORLD) ) return false;


  // ローカル送信バッファの確保 送信元のみ
  szSndBf = Maxlsz;
  SndBfR = new REAL_TYPE[szSndBf];
  memory += (double)szSndBf * (double)sizeof(REAL_TYPE);

  SndBfI = new int[szSndBf];
  memory += (double)szSndBf * (double)sizeof(int);


  // 受信バッファはマスターランクのみ
  if (myRank == 0)
  {
    szRcvBf = Maxlsz * numProc;
    RcvBfR = new REAL_TYPE[szRcvBf];
    memory += (double)szRcvBf * (double)sizeof(REAL_TYPE);
    RcvBfI = new int[szRcvBf];
    memory += (double)szRcvBf * (double)sizeof(int);
  }

  return true;
}
