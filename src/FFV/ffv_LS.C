//##################################################################################
//
// FFV-C : Frontflow / violet Cartesian
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2014 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

/** 
 * @file   ffv_LS.C
 * @brief  FFV Class
 * @author aics
 */

#include "ffv.h"

/**
 * @brief 2D配列アクセス
 * @param [in] _I  最初の列
 * @param [in] _J  次の列
 * @param [in] _SZ 配列サイズ（1D方向）
 */
#define _IDX2D(_I,_J,_SZ) (_J*_SZ+_I)


// #################################################################
// SOR2SMAの非同期通信処理
// 並列時のみコールされる
void FFV::comm_SOR2SMA(REAL_TYPE* d_x, const int col, const int ip, MPI_Request* key)
{
  // cf_sz バッファサイズ
  // cf_x x方向のバッファ (cf_sz[0]*4) 
  // cf_y y方向のバッファ (cf_sz[1]*4) 
  // cf_z z方向のバッファ (cf_sz[2]*4) 
  
  // key [12] >> [dir(0-5), send(0)/recv(1)]
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  int cx = cf_sz[0];
  int cy = cf_sz[1];
  int cz = cf_sz[2];
  
  int ic = (col + ip) % 2; // スタートインデクス
  

  // X_MINUS
  // send
  if ( nID[X_minus]>=0 )
  {
    int c = 0;
    int i = 1;
    
    for (int k=1; k<=kx; k++) {
      int js= (k + ic + i) % 2;
      
      for (int j=1+js; j<=jx; j+=2) {
        cf_x[_PACK(cx, face_m_send, c)] = d_p[ _F_IDX_S3D(i, j, k, ix, jx, kx, gd) ];
        c++;
      }
    }
    
    if ( paraMngr->Isend(&cf_x[_PACK(cx, face_m_send, 0)], cx, nID[X_minus], &key[_KEY_RQ(X_minus, key_send)]) != CPM_SUCCESS ) Exit(0);
  }
    
  // recieve
  if ( nID[X_plus]>=0 )
  {
    if ( paraMngr->Irecv(&cf_x[_PACK(cx, face_p_recv, 0)], cx, nID[X_plus], &key[_KEY_RQ(X_plus, key_recv)]) != CPM_SUCCESS ) Exit(0);
  }
    
    
  // X_PLUS
  // send
  if ( nID[X_plus]>=0 )
  {
    int c = 0;
    int i = ix;
    
    for (int k=1; k<=kx; k++) {
      int js= (k + ic + i) % 2;
      
      for (int j=1+js; j<=jx; j+=2) {
        cf_x[_PACK(cx, face_p_send, c)] = d_p[ _F_IDX_S3D(i, j, k, ix, jx, kx, gd) ];
        c++;
      }
    }
    
    if ( paraMngr->Isend(&cf_x[_PACK(cx, face_p_send, 0)], cx, nID[X_plus], &key[_KEY_RQ(X_plus, key_send)]) != CPM_SUCCESS ) Exit(0);
  }
  
  // recieve
  if ( nID[X_minus]>=0 )
  {
    if ( paraMngr->Irecv(&cf_x[_PACK(cx, face_m_recv, 0)], cx, nID[X_minus], &key[_KEY_RQ(X_minus, key_recv)]) != CPM_SUCCESS ) Exit(0);
  }
  
  
  // Y_MINUS
  // send
  if ( nID[Y_minus]>=0 )
  {
    int c = 0;
    int j = 1;
    
    for (int k=1; k<=kx; k++) {
      int is= (k + ic + j) % 2;
      
      for (int i=1+is; i<=ix; i+=2) {
        cf_y[_PACK(cy, face_m_send, c)] = d_p[ _F_IDX_S3D(i, j, k, ix, jx, kx, gd) ];
        c++;
      }
    }
    
    if ( paraMngr->Isend(&cf_y[_PACK(cy, face_m_send, 0)], cy, nID[Y_minus], &key[_KEY_RQ(Y_minus, key_send)]) != CPM_SUCCESS ) Exit(0);
  }
  
  // recv
  if ( nID[Y_plus]>=0 )
  {
    if ( paraMngr->Irecv(&cf_y[_PACK(cy, face_p_recv, 0)], cy, nID[Y_plus], &key[_KEY_RQ(Y_plus, key_recv)]) != CPM_SUCCESS ) Exit(0);
  }
  
  
  // Y_PLUS
  // send
  if ( nID[Y_plus]>=0 )
  {
    int c = 0;
    int j = jx;
    
    for (int k=1; k<=kx; k++) {
      int is= (k + ic + j) % 2;
      
      for (int i=1+is; i<=ix; i+=2) {
        cf_y[_PACK(cy, face_p_send, c)] = d_p[ _F_IDX_S3D(i, j, k, ix, jx, kx, gd) ];
        c++;
      }
    }
    
    if ( paraMngr->Isend(&cf_y[_PACK(cy, face_p_send, 0)], cy, nID[Y_plus], &key[_KEY_RQ(Y_plus, key_send)]) != CPM_SUCCESS ) Exit(0);
  }
  
  // recv
  if ( nID[Y_minus]>=0 )
  {
    if ( paraMngr->Irecv(&cf_y[_PACK(cy, face_m_recv, 0)], cy, nID[Y_minus], &key[_KEY_RQ(Y_minus, key_recv)]) != CPM_SUCCESS ) Exit(0);
  }
  
  
  // Z_MINUS
  // send
  if ( nID[Z_minus]>=0 )
  {
    int c = 0;
    int k = 1;
    
    for (int j=1; j<=jx; j++) {
      int is= (j + ic + k) % 2;
      
      for (int i=1+is; i<=ix; i+=2) {
        cf_z[_PACK(cz, face_m_send, c)] = d_p[ _F_IDX_S3D(i, j, k, ix, jx, kx, gd) ];
        c++;
      }
    }
    
    if ( paraMngr->Isend(&cf_z[_PACK(cz, face_m_send, 0)], cz, nID[Z_minus], &key[_KEY_RQ(Z_minus, key_send)]) != CPM_SUCCESS ) Exit(0);
  }
  
  // recv
  if ( nID[Z_plus]>=0 )
  {
    if ( paraMngr->Irecv(&cf_z[_PACK(cz, face_p_recv, 0)], cz, nID[Z_plus], &key[_KEY_RQ(Z_plus, key_recv)]) != CPM_SUCCESS ) Exit(0);
  }
  
  
  // Z_PLUS
  // send
  if ( nID[Z_plus]>=0 )
  {
    int c = 0;
    int k = kx;
    
    for (int j=1; j<=jx; j++) {
      int is= (j + ic + k) % 2;
      
      for (int i=1+is; i<=ix; i+=2) {
        cf_z[_PACK(cz, face_p_send, c)] = d_p[ _F_IDX_S3D(i, j, k, ix, jx, kx, gd) ];
        c++;
      }
    }
    
    if ( paraMngr->Isend(&cf_z[_PACK(cz, face_p_send, 0)], cz, nID[Z_plus], &key[_KEY_RQ(Z_plus, key_send)]) != CPM_SUCCESS ) Exit(0);
  }
  
  // recv
  if ( nID[Z_minus]>=0 )
  {
    if ( paraMngr->Irecv(&cf_z[_PACK(cz, face_m_recv, 0)], cz, nID[Z_minus], &key[_KEY_RQ(Z_minus, key_recv)]) != CPM_SUCCESS ) Exit(0);
  }
  
}


// #################################################################
// 種類Lの線形ソルバを利用する場合，trueを返す
bool FFV::hasLinearSolver(const int L)
{
  for (int i=0; i<ic_END; i++)
  {
    if ( IC[i].getLS() == L ) return true;
  }
  
  return false;
}


// #################################################################
// Point SOR
int FFV::Point_SOR(IterationCtl* IC, REAL_TYPE* x, REAL_TYPE* b, const double rhs_nrm, const double r0)
{
  double flop_count=0.0;         /// 浮動小数点演算数
  REAL_TYPE omg = IC->getOmega(); /// 加速係数
	double res = 0.0;              /// 残差
  int lc=0;                      /// ループカウント
  
  // x     圧力 p^{n+1}
	// b     RHS vector
	// d_bcp ビットフラグ
  
  TIMING_start(tm_poi_itr_sct_2); // >>> Poisson Iteration section 2
  
  for (lc=0; lc<IC->getMaxIteration(); lc++)
  {
    // 反復処理
    TIMING_start(tm_poi_PSOR);
    flop_count = 0.0;
    psor_(x, size, &guide, &omg, &res, b, d_bcp, &flop_count);
    TIMING_stop(tm_poi_PSOR, flop_count);
    
    // 境界条件
    TIMING_start(tm_poi_BC);
    BC.OuterPBC(x);
    if ( C.EnsCompo.periodic == ON ) BC.InnerPBCperiodic(x, d_bcd);
    TIMING_stop(tm_poi_BC, 0.0);
    
    // 同期処理
    Sync_Scalar(IC, x, 1);
    
    // 残差の集約
    if ( numProc > 1 )
    {
      TIMING_start(tm_poi_res_comm);
      double tmp = res;
      if ( paraMngr->Allreduce(&tmp, &res, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
      TIMING_stop(tm_poi_res_comm, 2.0*numProc*sizeof(double) ); // 双方向 x ノード数
    }
    
    res = sqrt(res); // L2 of (x^{m+1} - x^m)
    
    
    // Residual
    switch ( IC->getNormType() )
    {
      case r_b:
      case r_r0:
        
        TIMING_start(tm_poi_src_nrm);
        res = 0.0;
        flop_count = 0.0;
        poi_residual_(&res, size, &guide, x, b, d_bcp, &flop_count);
        TIMING_stop(tm_poi_src_nrm, flop_count);
        
        if ( numProc > 1 )
        {
          TIMING_start(tm_poi_src_comm);
          double m_tmp = res;
          if ( paraMngr->Allreduce(&m_tmp, &res, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
          TIMING_stop(tm_poi_src_comm, 2.0*numProc*sizeof(double) ); // 双方向 x ノード数
        }
        
        res = sqrt(res);
        break;
        
      case dx_b:
        // nothing to do, dx is already obtained in psor_(&res,...)
        break;
    }
    
    // 残差の保存
    switch ( IC->getNormType() )
    {
      case dx_b:
        IC->setNormValue( (rhs_nrm==0.0) ? res : res/rhs_nrm );
        break;
        
      case r_b:
        IC->setNormValue( (rhs_nrm==0.0) ? res : res/rhs_nrm );
        break;
        
      case r_r0:
        IC->setNormValue( (r0==0.0) ? res : res/r0 );
        break;
        
      default:
        printf("\tInvalid Linear Solver for Pressure\n");
        Exit(0);
        break;
    }

    // 収束判定　性能測定モードのときは収束判定を行わない
    if ( (C.Hide.PM_Test == OFF) && (IC->getNormValue() < IC->getCriterion()) ) break;
    
  }
  
  TIMING_stop(tm_poi_itr_sct_2, 0.0); // <<< Poisson Iteration subsection 2
  
  return lc;
}


// #################################################################
// 反復変数の同期処理
void FFV::Sync_Scalar(IterationCtl* IC, REAL_TYPE* d_class, const int num_layer)
{
  if ( numProc > 1 )
  {
    TIMING_start(tm_poi_comm);
    
    if (IC->getSyncMode() == comm_sync )
    {
      if ( paraMngr->BndCommS3D(d_class, size[0], size[1], size[2], guide, num_layer) != CPM_SUCCESS ) Exit(0);
    }
    else
    {
      MPI_Request req[12];
      for (int i=0; i<12; i++) req[i] = MPI_REQUEST_NULL;
      
      if ( paraMngr->BndCommS3D_nowait(d_class, size[0], size[1], size[2], guide, num_layer, req ) != CPM_SUCCESS ) Exit(0);
      if ( paraMngr->wait_BndCommS3D  (d_class, size[0], size[1], size[2], guide, num_layer, req ) != CPM_SUCCESS ) Exit(0);
    }
    TIMING_stop(tm_poi_comm, face_comm_size*(double)num_layer);
  }
}



// #################################################################
// SOR2SMA
int FFV::SOR_2_SMA(IterationCtl* IC, REAL_TYPE* x, REAL_TYPE* b, const double rhs_nrm, const double r0)
{
  int ip;                         /// ローカルノードの基点(1,1,1)のカラーを示すインデクス
                                  /// ip=0 > R, ip=1 > B
  double flop_count=0.0;          /// 浮動小数点演算数
  REAL_TYPE omg = IC->getOmega(); /// 加速係数
	double res = 0.0;               /// 残差
  int lc=0;                       /// ループカウント
  
  // x     圧力 p^{n+1}
  // b     RHS vector
	// d_bcp ビットフラグ

  
  TIMING_start(tm_poi_itr_sct_2); // >>> Poisson Iteration section 2
  
  for (lc=0; lc<IC->getMaxIteration(); lc++)
  {
    // 2色のマルチカラー(Red&Black)のセットアップ
    TIMING_start(tm_poi_setup);
    
    MPI_Request req[12]; /// 送信ID
    
    for (int i=0; i<12; i++) req[i] = MPI_REQUEST_NULL;
    
    // ip = 0 基点(1,1,1)がRからスタート
    //    = 1 基点(1,1,1)がBからスタート
    if ( numProc > 1 )
    {
      ip = (head[0]+head[1]+head[2]+1) % 2;
    }
    else
    {
      ip = 0;
    }
    TIMING_stop(tm_poi_setup, 0.0);
    
    
    // 各カラー毎の間に同期, 残差は色間で積算する
    // R - color=0 / B - color=1
    for (int color=0; color<2; color++) {
      
      TIMING_start(tm_poi_SOR2SMA);
      flop_count = 0.0; // 色間で積算しない
      if ( IC->getNaive()==OFF)
      {
        if (IC->getBit3() == OFF)
        {
          psor2sma_core_(x, size, &guide, &ip, &color, &omg, &res, b, d_bcp, &flop_count);
        }
        else
        {
          psor2sma_core_bit3_(x, size, &guide, &ip, &color, &omg, &res, b, d_bcp, &flop_count);
        }
      }
      else
      {
        psor2sma_naive_(x, size, &guide, &ip, &color, &omg, &res, b, d_bcp, d_pni, &flop_count);
      }
      
      TIMING_stop(tm_poi_SOR2SMA, flop_count);
      
      // 境界条件
      TIMING_start(tm_poi_BC);
      BC.OuterPBC(x);
      if ( C.EnsCompo.periodic == ON ) BC.InnerPBCperiodic(x, d_bcd);
      TIMING_stop(tm_poi_BC, 0.0);
      
      
      // 同期処理
      if ( numProc > 1 )
      {
        TIMING_start(tm_poi_comm);
        
        
        if (IC->getSyncMode() == comm_sync )
        {
          if ( paraMngr->BndCommS3D(x, size[0], size[1], size[2], guide, 1) != CPM_SUCCESS ) Exit(0); // 1 layer communication
        }
        else
        {
          //if ( paraMngr->BndCommS3D_nowait(x, size[0], size[1], size[2], guide, 1, req) != CPM_SUCCESS ) Exit(0);
          //if ( paraMngr->wait_BndCommS3D  (x, size[0], size[1], size[2], guide, 1, req) != CPM_SUCCESS ) Exit(0);
          //comm_SOR2SMA(x, color, ip, req);
          //wait_SOR2SMA(x, color, ip, req);
          int ireq[12];
          sma_comm_     (x, size, &guide, &color, &ip, cf_sz, cf_x, cf_y, cf_z, ireq, nID);
          sma_comm_wait_(x, size, &guide, &color, &ip, cf_sz, cf_x, cf_y, cf_z, ireq);
        }
        TIMING_stop(tm_poi_comm, face_comm_size*0.5);
      }
    }
    
    
    // Residual
    switch ( IC->getNormType() )
    {
      case r_b:
      case r_r0:
        
        TIMING_start(tm_poi_src_nrm);
        res = 0.0;
        poi_residual_(&res, size, &guide, x, b, d_bcp, &flop_count);
        TIMING_stop(tm_poi_src_nrm, flop_count);
        break;
        
      case dx_b:
        // nothing to do, dx is already obtained in psor_(&res,...)
        break;
    }
    
    if ( numProc > 1 )
    {
      TIMING_start(tm_poi_src_comm);
      double m_tmp = res;
      if ( paraMngr->Allreduce(&m_tmp, &res, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
      TIMING_stop(tm_poi_src_comm, 2.0*numProc*sizeof(double) ); // 双方向 x ノード数
    }
    
    res = sqrt(res);
    
    
    // 残差の保存
    switch ( IC->getNormType() )
    {
      case dx_b:
        IC->setNormValue( (rhs_nrm==0.0) ? res : res/rhs_nrm );
        break;
        
      case r_b:
        IC->setNormValue( (rhs_nrm==0.0) ? res : res/rhs_nrm );
        break;
        
      case r_r0:
        IC->setNormValue( (r0==0.0) ? res : res/r0 );
        break;
        
      default:
        printf("\tInvalid Linear Solver for Pressure\n");
        Exit(0);
        break;
    }
    
    // 収束判定　性能測定モードのときは収束判定を行わない
    if ( (C.Hide.PM_Test == OFF) && (IC->getNormValue() < IC->getCriterion()) ) break;
    
  }
  
  TIMING_stop(tm_poi_itr_sct_2, 0.0); // <<< Poisson Iteration subsection 2
  
  return lc;
}


// #################################################################
// SOR2SMAの非同期通信処理
// 並列時のみコールされる
void FFV::wait_SOR2SMA(REAL_TYPE* d_x, const int col, const int ip, MPI_Request* key)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  int cx = cf_sz[0];
  int cy = cf_sz[1];
  int cz = cf_sz[2];
  
  int ic = (col + ip) % 2; // スタートインデクス
  
  int rq;
  
  
  // X_PLUS面の受信バッファを展開
  rq = _KEY_RQ(X_plus, key_recv);
  
  if (key[rq]>=0 )
  {
    if ( paraMngr->Wait(&key[rq]) != CPM_SUCCESS ) Exit(0);
    
    int c = 0;
    int i = ix+1;
    
    for (int k=1; k<=kx; k++) {
      int js= (k + ic + i) % 2;
      
      for (int j=1+js; j<=jx; j+=2) {
        d_x[ _F_IDX_S3D(i, j, k, ix, jx, kx, gd) ] = cf_x[_PACK(cx, face_p_recv, c)];
        c++;
      }
    }
  }
  
  // X_MINUS面の受信バッファを展開
  rq = _KEY_RQ(X_minus, key_recv);
  
  if (key[rq]>=0 )
  {
    if ( paraMngr->Wait(&key[rq]) != CPM_SUCCESS ) Exit(0);
    
    int c = 0;
    int i = 0;
    
    for (int k=1; k<=kx; k++) {
      int js= (k + ic + i) % 2;
      
      for (int j=1+js; j<=jx; j+=2) {
        d_x[ _F_IDX_S3D(i, j, k, ix, jx, kx, gd) ] = cf_x[_PACK(cx, face_m_recv, c)];
        c++;
      }
    }
  }
          
  // Y_PLUS面の受信バッファを展開
  rq = _KEY_RQ(Y_plus, key_recv);
  
  if (key[rq]>=0 )
  {
    if ( paraMngr->Wait(&key[rq]) != CPM_SUCCESS ) Exit(0);
    
    int c = 0;
    int j = jx+1;
    
    for (int k=1; k<=kx; k++) {
      int is= (k + ic + j) % 2;
      
      for (int i=1+is; i<=ix; i+=2) {
        d_x[ _F_IDX_S3D(i, j, k, ix, jx, kx, gd) ] = cf_y[_PACK(cy, face_p_recv, c)];
        c++;
      }
    }
  }
  
  // Y_MINUS面の受信バッファを展開
  rq = _KEY_RQ(Y_minus, key_recv);
  
  if (key[rq]>=0 )
  {
    if ( paraMngr->Wait(&key[rq]) != CPM_SUCCESS ) Exit(0);
    
    int c = 0;
    int j = 0;
    
    for (int k=1; k<=kx; k++) {
      int is= (k + ic + j) % 2;
      
      for (int i=1+is; i<=ix; i+=2) {
        d_x[ _F_IDX_S3D(i, j, k, ix, jx, kx, gd) ] = cf_y[_PACK(cy, face_m_recv, c)];
        c++;
      }
    }
  }
  
  // Z_PLUS面の受信バッファを展開
  rq = _KEY_RQ(Z_plus, key_recv);
  
  if (key[rq]>=0 )
  {
    if ( paraMngr->Wait(&key[rq]) != CPM_SUCCESS ) Exit(0);
    
    int c = 0;
    int k = kx+1;
    
    for (int j=1; j<=jx; j++) {
      int is= (j + ic + k) % 2;
      
      for (int i=1+is; i<=ix; i+=2) {
        d_x[ _F_IDX_S3D(i, j, k, ix, jx, kx, gd) ] = cf_z[_PACK(cz, face_p_recv, c)];
        c++;
      }
    }
  }
  
  // Z_MINUS面の受信バッファを展開
  rq = _KEY_RQ(Z_minus, key_recv);
  
  if (key[rq]>=0 )
  {
    if ( paraMngr->Wait(&key[rq]) != CPM_SUCCESS ) Exit(0);
    
    int c = 0;
    int k = 0;
    
    for (int j=1; j<=jx; j++) {
      int is= (j + ic + k) % 2;
      
      for (int i=1+is; i<=ix; i+=2) {
        d_x[ _F_IDX_S3D(i, j, k, ix, jx, kx, gd) ] = cf_z[_PACK(cz, face_m_recv, c)];
        c++;
      }
    }
  }
  
  
  // Wait for send ------------------------
  
  // X_MINUS
  rq = _KEY_RQ(X_minus, key_send);
  
  if (key[rq]>=0 )
  {
    if ( paraMngr->Wait(&key[rq]) != CPM_SUCCESS ) Exit(0);
  }
  
  // X_PLUS
  rq = _KEY_RQ(X_plus, key_send);
  
  if (key[rq]>=0 )
  {
    if ( paraMngr->Wait(&key[rq]) != CPM_SUCCESS ) Exit(0);
  }
  
  // Y_MINUS
  rq = _KEY_RQ(Y_minus, key_send);
  
  if (key[rq]>=0 )
  {
    if ( paraMngr->Wait(&key[rq]) != CPM_SUCCESS ) Exit(0);
  }
  
  // Y_PLUS
  rq = _KEY_RQ(Y_plus, key_send);
  
  if (key[rq]>=0 )
  {
    if ( paraMngr->Wait(&key[rq]) != CPM_SUCCESS ) Exit(0);
  }
  
  // Z_MINUS
  rq = _KEY_RQ(Z_minus, key_send);
  
  if (key[rq]>=0 )
  {
    if ( paraMngr->Wait(&key[rq]) != CPM_SUCCESS ) Exit(0);
  }
  
  // Z_PLUS
  rq = _KEY_RQ(Z_plus, key_send);
  
  if (key[rq]>=0 )
  {
    if ( paraMngr->Wait(&key[rq]) != CPM_SUCCESS ) Exit(0);
  }
}



// #################################################################
// Flexible Gmres(m)
void FFV::Fgmres(IterationCtl* IC, const double rhs_nrm, const double r0)
{
  const double eps_1 = 1.0e-30;
  const double eps_2 = IC->getCriterion();
  int mode_precond = 1; // pre-conditioning
  
  // 残差収束チェック
  if ( rhs_nrm < eps_1 ) return;
  
  // d_wg : work
  // d_ws : h^2 \Psi
  // d_p  : x
  // 
  
  // >>> Gmres section
  TIMING_start(tm_gmres_sor_sct);
  
  double beta, beta1, res, eps_abs, res_abs, al;
  double r4;
  double flop=0.0;
  
  const int Iteration_Max = IC->getMaxIteration();
  int m = FREQ_OF_RESTART;
  int s_length = (size[0]+2*guide) * (size[1]+2*guide) * (size[2]+2*guide);
  
  int isfin  = 0;
  
  const int oki  = 2;
  const int step = 3;
  int nrm  = 0;
  
  
  int Nmax = m+1; // 配列確保用，ゼロはダミー
      
  double *hg = new double[(Nmax+1) * (Nmax+1)]; // Hessnberg matrix
  double *cs = new double[Nmax+1];              // Cosine for Givens rotations
  double *sn = new double[Nmax+1];              // Sine   for Givens rotations
  double *rm = new double[Nmax+1];              // residual vector for minimization problem
  
  memset(hg, 0, (Nmax+1)*(Nmax+1)*sizeof(double) );
  memset(cs, 0, (Nmax+1)*sizeof(double) );
  memset(sn, 0, (Nmax+1)*sizeof(double) );
  memset(rm, 0, (Nmax+1)*sizeof(double) );
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  int first = 1;

  
  // Outer iteration
  for (int j=1; j<=Iteration_Max; j++) {
    
    TIMING_start(tm_gmres_mvprod);
    flop = 0.0;
    matvec_p_(d_wg, size, &guide, d_p, d_bcp, &flop); // d_wg <- Ax
    TIMING_stop(tm_gmres_mvprod, flop);
    
    
    TIMING_start(tm_gmres_res_sample);
    flop = 0.0;
    residual_(d_wg, size, &guide, &res_abs, d_ws, d_sq, d_bcp, &flop); // d_wg <- b - Ax (= r); res_abs <- \sum{r^2}
    TIMING_stop(tm_gmres_res_sample, flop);
    
    
    TIMING_start(tm_gmres_comm);
    if ( numProc > 1 )
    {
      double tmp = res_abs;
      if ( paraMngr->Allreduce(&tmp, &res_abs, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
    }
    TIMING_stop(tm_gmres_comm, 2.0*numProc*sizeof(double)); // 双方向 x ノード数
    
    
    beta = sqrt(res_abs); // beta <- \|r\|_2
    
    // convergence check
    if (beta < eps_1) goto jump_2;
    
    
    
    // set initial value of residual vector
    rm[1] = beta; // rm <- {\beta, 0, 0, ..., 0}^T

    
    beta1 = 1.0 / beta;
    
    orth_basis_(d_vm, size, &guide, &Nmax, &first, &beta1, d_wg, &flop); // d_vm(1) <- v^1 
    
    
    
    // Orthogoanl basis
    for (int i=1; i<=m; i++) {
      
      // address
      int g = guide;
      size_t adrs = _F_IDX_S4D(1-g, 1-g, 1-g, i-1, size[0], size[1], size[2], guide);
      
      if ( mode_precond == 1)
      {
        // Decide the number of iteration for inner-iteration
        int n_inner = 3; //((im / oki) + 1) * step;
        
        // Inner-iteration
        for (int i_inner=1; i_inner<=n_inner; i_inner++) {
          SOR_2_SMA(IC, &d_zm[adrs], &d_vm[adrs], rhs_nrm, r0); // z^i <- K^{-1} v^i
        }
      }
      else
      {
        cp_orth_basis_(&d_zm[adrs], size, &guide, &d_vm[adrs]); // z^i <- v^i
      }
      
      TIMING_start(tm_gmres_mvprod);
      flop = 0.0;
      matvec_p_(d_wg, size, &guide, &d_zm[adrs], d_bcp, &flop); // d_wg <- Az^i
      TIMING_stop(tm_gmres_mvprod, flop);
      
      
      for (int km=1; km<=i; km++) {
        al = 0.0;
        ml_add_1_(&al, size, &guide, &Nmax, d_vm, d_wg, &km, &flop);
        
        
        TIMING_start(tm_gmres_comm);
        if ( numProc > 1 )
        {
          double tmp = al;
          if ( paraMngr->Allreduce(&tmp, &al, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
        }
        TIMING_stop(tm_gmres_comm, 2.0*numProc*sizeof(double)); // 双方向 x ノード数
        
        
        hg[_IDX2D(km, i, Nmax)] = al;
        
        double r_al = -al;
        ml_add_3_(d_wg, size, &guide, &Nmax, &r_al, d_vm, &km, &flop);
      }
      
      
      al =0.0;
      ml_add_2_(&al, size, &guide, d_wg, &flop);
      
      
      TIMING_start(tm_gmres_comm);
      if ( numProc > 1 )
      {
        double tmp = al;
        if ( paraMngr->Allreduce(&tmp, &al, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
      }
      TIMING_stop(tm_gmres_comm, 2.0*numProc*sizeof(double)); // 双方向 x ノード数
      
      
      TIMING_start(tm_gmres_others);
      flop = 0.0;
      
      hg[_IDX2D(i+1, i, Nmax)] = sqrt(al);
      flop += 20.0;
      
      if ( hg[_IDX2D(i+1, i, Nmax)] < eps_1 )
      {
        nrm = i-1;
        TIMING_stop(tm_gmres_others, flop);
        goto jump_1;
      }
      
      if (i < Nmax)
      {
        int idx = i + 1;
        al = 1.0 / hg[_IDX2D(idx, i, Nmax)];
        flop += 13.0;
        orth_basis_(d_vm, size, &guide, &Nmax, &idx, &al, d_wg, &flop);
      }
      
      rm[_IDX2D(1, i, Nmax)] = hg[_IDX2D(1, i, Nmax)];
      
      for (int km=2; km<=i; km++) {
        rm[_IDX2D(km  , i, Nmax)] = cs[km-1]*hg[_IDX2D(km, i, Nmax)] - sn[km-1]*rm[_IDX2D(km-1, i, Nmax)];
        rm[_IDX2D(km-1, i, Nmax)] = sn[km-1]*hg[_IDX2D(km, i, Nmax)] + cs[km-1]*rm[_IDX2D(km-1, i, Nmax)];
      }
      flop += (double)(i-2+1)*6.0;
      
      al = sqrt(rm[_IDX2D(i, i, Nmax)]*rm[_IDX2D(i, i, Nmax)] + hg[_IDX2D(i+1, i, Nmax)]*hg[_IDX2D(i+1, i, Nmax)]);
      flop += 23;
      
      if ( al < eps_1 )
      {
        nrm = i - 1;
        TIMING_stop(tm_gmres_others, flop);
        goto jump_1;
      }
      
      cs[i]    = rm[_IDX2D(i,  i, Nmax)] / al;
      sn[i]    = hg[_IDX2D(i+1,i, Nmax)] / al;
      
      rm[_IDX2D(i, i, Nmax)] = cs[i] * rm[_IDX2D(i, i, Nmax)] + sn[i] * hg[_IDX2D(i+1, i, Nmax)];
      
      rm[i+1]  = - sn[i]*rm[i];
      rm[i]    =   cs[i]*rm[i];
      
      al = rm[i+1] * rm[i+1];
      
      flop += 2.0*13.0 + 4.0;
      
      if (al < eps_abs)
      {
        nrm = i;
        TIMING_stop(tm_gmres_others, flop);
        goto jump_1;
      }
      
      TIMING_stop(tm_gmres_others, flop);
      
    } // loop; im
    
    nrm = Nmax;
    
  jump_1:
    
    rm[nrm] = rm[nrm] / rm[_IDX2D(nrm, nrm, Nmax)];
    
    for (int im=nrm-1; im>=1; im--) {
      al = rm[im];
      
      for (int jm=im+1; jm<=nrm; jm++) {
        al -= rm[_IDX2D(im, jm, Nmax)] * rm[jm];
      }
      
      rm[im] = al / rm[_IDX2D(im, im, Nmax)];
    }
    
    for (int im=1; im<=nrm; im++) {
      al = rm[im];
      
      ml_add_4_(d_p, size, &guide, &Nmax, &al, d_zm, &im, &flop);
    }
    
  } // loop; j
  
jump_2:
  
  SOR_2_SMA(IC, &d_zm[0], &d_vm[0], rhs_nrm, r0);
  
  
jump_4:
  
  if ( (isfin == 0) && (res < (rhs_nrm * eps_2 * eps_2)) )
  {
    res = 100.0 * rhs_nrm * eps_2 * eps_2;
  }
  
  if ( rm ) delete [] rm;
  if ( hg ) delete [] hg;
  if ( cs ) delete [] cs;
  if ( sn ) delete [] sn;
  
  TIMING_stop(tm_gmres_sor_sct, 0.0);
  // <<< Poisson Source section

}

// #################################################################
// RBGS
int FFV::Frbgs(IterationCtl* IC, REAL_TYPE* x, REAL_TYPE* b, const double rhs_nrm, const double r0)
{
  REAL_TYPE omg = IC->getOmega();
  int lc=0;                      /// ループカウント
  for (lc=0; lc<IC->getMaxIteration(); lc++)
  {
		Fsmoother(x, b, omg);
    
		REAL_TYPE rr = 0.0;
		blas_calcr2_(&rr, x, b, d_bcp, size, &guide);
    
    if ( numProc > 1 )
    {
      REAL_TYPE m_tmp = rr;
      if ( paraMngr->Allreduce(&m_tmp, &rr, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
    }
    
    REAL_TYPE res = sqrt(rr);
    
		if ( Fcheck(IC, res, rhs_nrm, r0) == true ) break;
    
	}
	Sync_Scalar(IC, x, 1);
  
	return lc;
}

// #################################################################
// PCG
int FFV::Fpcg(IterationCtl* IC, REAL_TYPE* x, REAL_TYPE* b, const double rhs_nrm, const double r0)
{
	REAL_TYPE res = 0.0;
  REAL_TYPE omg = IC->getOmega();
  
	blas_clear_(d_pcg_r , size, &guide);
	blas_clear_(d_pcg_p , size, &guide);
	blas_clear_(d_pcg_q , size, &guide);
	blas_clear_(d_pcg_z , size, &guide);
  
	blas_calcr_(d_pcg_r, x, b, d_bcp, size, &guide);
	Sync_Scalar(IC, d_pcg_r, 1);
  
	REAL_TYPE rr0 = 1.0;
	REAL_TYPE rr1 = 0.0;
  int lc=0;                      /// ループカウント
  
  for (lc=0; lc<IC->getMaxIteration(); lc++)
  {
		blas_clear_(d_pcg_z, size, &guide);
		Fpreconditioner(IC, d_pcg_z, d_pcg_r);
    
		rr1 = 0.0;
		Fdot(&rr1, d_pcg_r, d_pcg_z);
    
		if( fabs(rr1) < FLT_MIN ) {
			res = rr1;
			lc = 0;
			break;
		}
    
		REAL_TYPE beta = rr1/rr0;
    
		if ( lc == 1 )
    {
			blas_copy_(d_pcg_p, d_pcg_z, size, &guide);
		}
    else
    {
			blas_xpay_(d_pcg_p, d_pcg_z, &beta, size, &guide);
		}
		Sync_Scalar(IC, d_pcg_p, 1);
    
		blas_calcax_(d_pcg_q, d_pcg_p, d_bcp, size, &guide);
    
		REAL_TYPE qp = 0.0;
		Fdot(&qp, d_pcg_q, d_pcg_p);
    
		REAL_TYPE alpha  = rr1/qp;
		REAL_TYPE alphan = -alpha;
    
		blas_axpy_(x      , d_pcg_p, &alpha , size, &guide);
		blas_axpy_(d_pcg_r, d_pcg_q, &alphan, size, &guide);
		Sync_Scalar(IC, d_pcg_r, 1);
    
		REAL_TYPE rr = 0.0;
		Fdot(&rr, d_pcg_r, d_pcg_r);
    
    res = sqrt(rr);
    
		if( Fcheck(IC, res, rhs_nrm, r0) == true ) break;
    
		rr0 = rr1;
	}
  
	BC.OuterPBC(x);
  
	if ( C.EnsCompo.periodic == ON )
  {
		BC.InnerPBCperiodic(x, d_bcd);
	}
  
	Sync_Scalar(IC, x, 1);
  
	return lc;
}

// #################################################################
// PBiCBSTAB
int FFV::Fpbicgstab(IterationCtl* IC, REAL_TYPE* x, REAL_TYPE* b, const double rhs_nrm, const double r0)
{
	REAL_TYPE res = 0.0;
  REAL_TYPE omg = IC->getOmega();
  
	blas_clear_(d_pcg_r , size, &guide);
	blas_clear_(d_pcg_p , size, &guide);
  
	blas_clear_(d_pcg_r0, size, &guide);
	blas_clear_(d_pcg_p_, size, &guide);
	blas_clear_(d_pcg_q_, size, &guide);
	blas_clear_(d_pcg_s , size, &guide);
	blas_clear_(d_pcg_s_, size, &guide);
	blas_clear_(d_pcg_t_, size, &guide);
  
	blas_calcr_(d_pcg_r, x, b, d_bcp, size, &guide);
	Sync_Scalar(IC, d_pcg_r, 1);
  
	blas_copy_(d_pcg_r0, d_pcg_r, size, &guide);
  
	REAL_TYPE rr0 = 1.0;
	REAL_TYPE alpha = 0.0;
	REAL_TYPE gamma  = 1.0;
	REAL_TYPE gamman = -gamma;
  int lc=0;                      /// ループカウント
  
  for (lc=0; lc<IC->getMaxIteration(); lc++)
  {
		REAL_TYPE rr1 = 0.0;
		Fdot(&rr1, d_pcg_r, d_pcg_r0);
    
		if( fabs(rr1) < FLT_MIN )
    {
			res = rr1;
			lc = 0;
			break;
		}
    
		if( lc == 1 )
    {
			blas_copy_(d_pcg_p, d_pcg_r, size, &guide);
		}
    else
    {
			REAL_TYPE beta = rr1/rr0*alpha/gamma;
			blas_axpy_(d_pcg_p, d_pcg_q_, &gamman, size, &guide);
			blas_xpay_(d_pcg_p, d_pcg_r , &beta  , size, &guide);
		}
		Sync_Scalar(IC, d_pcg_p, 1);
    
		blas_clear_(d_pcg_p_, size, &guide);
		Fpreconditioner(IC, d_pcg_p_, d_pcg_p);
    
		blas_calcax_(d_pcg_q_, d_pcg_p_, d_bcp, size, &guide);
    
		REAL_TYPE q_r0 = 0.0;
		Fdot(&q_r0, d_pcg_q_, d_pcg_r0);
    
		REAL_TYPE alpha  = rr1/q_r0;
		REAL_TYPE alphan = -alpha;
		blas_axpyz_(d_pcg_s, d_pcg_q_, d_pcg_r, &alphan, size, &guide);
		Sync_Scalar(IC, d_pcg_s, 1);
    
		blas_clear_(d_pcg_s_, size, &guide);
		Fpreconditioner(IC, d_pcg_s_, d_pcg_s);
    
		blas_calcax_(d_pcg_t_, d_pcg_s_, d_bcp, size, &guide);
    
		REAL_TYPE t_s = 0.0;
		Fdot(&t_s, d_pcg_t_, d_pcg_s);
    
		REAL_TYPE t_t_ = 0.0;
		Fdot(&t_t_, d_pcg_t_, d_pcg_t_);
    
		gamma  = t_s/t_t_;
		gamman = -gamma;
    
		blas_axpbypz_(x      , d_pcg_p_, d_pcg_s_, &alpha , &gamma, size, &guide);
		blas_axpyz_  (d_pcg_r, d_pcg_t_, d_pcg_s , &gamman, size, &guide);
    
		REAL_TYPE rr = 0.0;
		Fdot(&rr, d_pcg_r, d_pcg_r);
    
    res = sqrt(rr);
    
		if ( Fcheck(IC, res, rhs_nrm, r0) == true ) break;
    
		rr0 = rr1;
	}
  
	BC.OuterPBC(x);
	if ( C.EnsCompo.periodic == ON )
  {
		BC.InnerPBCperiodic(x, d_bcd);
	}
  
	Sync_Scalar(IC, x, 1);
  
	return lc;
}

// #################################################################
// Check
bool FFV::Fcheck(IterationCtl* IC, REAL_TYPE res, const double rhs_nrm, const double r0)
{
	switch ( IC->getNormType() )
  {
		case dx_b:
			IC->setNormValue( (rhs_nrm==0.0) ? res : res/rhs_nrm );
			break;
		case r_b:
			IC->setNormValue( (rhs_nrm==0.0) ? res : res/rhs_nrm );
			break;
		case r_r0:
			IC->setNormValue( (r0==0.0) ? res : res/r0 );
			break;
		default:
			printf("\tInvalid Linear Solver for Pressure\n");
			Exit(0);
			break;
	}
  
	if ( (C.Hide.PM_Test == OFF) && (IC->getNormValue() < IC->getCriterion()) )
  {
		return true;
	}
  
	return false;
}

// #################################################################
// Preconditioner
int FFV::Fpreconditioner(IterationCtl* IC, REAL_TYPE* x, REAL_TYPE* b)
{
  REAL_TYPE omg = IC->getOmega();
  
  int lc=0;                      /// ループカウント
	int lc_max = 4;
  
  // 前処理なし(コピー)
	if( lc_max == 0 )
  {
		blas_copy_(x, b, size, &guide);
		return lc;
	}
  
  for (lc=0; lc<lc_max; lc++)
  {
		Fsmoother(x, b, omg);
	}
  
	return lc;
}

// #################################################################
// Smoother
void FFV::Fsmoother(REAL_TYPE* x, REAL_TYPE* b, REAL_TYPE omg)
{
	int ip = 0;
  
	if ( numProc > 1 )
  {
		ip = (head[0]+head[1]+head[2]+1) % 2;
	}
  else
  {
		ip = 0;
	}
  
	for (int color=0; color<2; color++)
  {
		blas_smoother_core_(x, b, d_bcp, &ip, &color, &omg, size, &guide);
    
		BC.OuterPBC(x);
    
		if ( C.EnsCompo.periodic == ON )
    {
			BC.InnerPBCperiodic(x, d_bcd);
		}
    
		if ( numProc > 1 )
    {
			if (IC->getSyncMode() == comm_sync )
      {
				if ( paraMngr->BndCommS3D(x, size[0], size[1], size[2], guide, 1) != CPM_SUCCESS ) {
					Exit(0);
				}
			}
      else
      {
				int ireq[12];
				sma_comm_     (x, size, &guide, &color, &ip, cf_sz, cf_x, cf_y, cf_z, ireq, nID);
				sma_comm_wait_(x, size, &guide, &color, &ip, cf_sz, cf_x, cf_y, cf_z, ireq);
			}
		}
	}
}

// #################################################################
// Dot
void FFV::Fdot(REAL_TYPE* xy, REAL_TYPE* x, REAL_TYPE* y)
{
	blas_dot_(xy, x, y, size, &guide);
  
	if ( numProc > 1 )
  {
		REAL_TYPE xy_tmp = *xy;
		if  ( paraMngr->Allreduce(&xy_tmp, xy, 1, MPI_SUM) != CPM_SUCCESS )Exit(0);
	}
}
