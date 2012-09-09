// #################################################################
//
// FFV : Frontflow / violet Cartesian
//
// Copyright (c) 2012 All right reserved.
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
//
// #################################################################

/** 
 * @file   ffv_LS.C
 * @brief  FFV Class
 * @author kero
 */

#include "ffv.h"


// SOR2SMAの非同期通信処理
// 並列時のみコールされる
void FFV::comm_SOR2SMA(const int col, const int ip, MPI_Request* key)
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
  if ( nID[X_MINUS]>=0 )
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
    
    if ( paraMngr->Isend(&cf_x[_PACK(cx, face_m_send, 0)], cx, nID[X_MINUS], &key[_KEY_RQ(X_MINUS, key_send)]) != CPM_SUCCESS ) Exit(0);
  }
    
  // recieve
  if ( nID[X_PLUS]>=0 )
  {
    if ( paraMngr->Irecv(&cf_x[_PACK(cx, face_p_recv, 0)], cx, nID[X_PLUS], &key[_KEY_RQ(X_PLUS, key_recv)]) != CPM_SUCCESS ) Exit(0);
  }
    
    
  // X_PLUS
  // send
  if ( nID[X_PLUS]>=0 )
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
    
    if ( paraMngr->Isend(&cf_x[_PACK(cx, face_p_send, 0)], cx, nID[X_PLUS], &key[_KEY_RQ(X_PLUS, key_send)]) != CPM_SUCCESS ) Exit(0);
  }
  
  // recieve
  if ( nID[X_MINUS]>=0 )
  {
    if ( paraMngr->Irecv(&cf_x[_PACK(cx, face_m_recv, 0)], cx, nID[X_MINUS], &key[_KEY_RQ(X_MINUS, key_recv)]) != CPM_SUCCESS ) Exit(0);
  }
  
  
  // Y_MINUS
  // send
  if ( nID[Y_MINUS]>=0 )
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
    
    if ( paraMngr->Isend(&cf_y[_PACK(cy, face_m_send, 0)], cy, nID[Y_MINUS], &key[_KEY_RQ(Y_MINUS, key_send)]) != CPM_SUCCESS ) Exit(0);
  }
  
  // recv
  if ( nID[Y_PLUS]>=0 )
  {
    if ( paraMngr->Irecv(&cf_y[_PACK(cy, face_p_recv, 0)], cy, nID[Y_PLUS], &key[_KEY_RQ(Y_PLUS, key_recv)]) != CPM_SUCCESS ) Exit(0);
  }
  
  
  // Y_PLUS
  // send
  if ( nID[Y_PLUS]>=0 )
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
    
    if ( paraMngr->Isend(&cf_y[_PACK(cy, face_p_send, 0)], cy, nID[Y_PLUS], &key[_KEY_RQ(Y_PLUS, key_send)]) != CPM_SUCCESS ) Exit(0);
  }
  
  // recv
  if ( nID[Y_MINUS]>=0 )
  {
    if ( paraMngr->Irecv(&cf_y[_PACK(cy, face_m_recv, 0)], cy, nID[Y_MINUS], &key[_KEY_RQ(Y_MINUS, key_recv)]) != CPM_SUCCESS ) Exit(0);
  }
  
  
  // Z_MINUS
  // send
  if ( nID[Z_MINUS]>=0 )
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
    
    if ( paraMngr->Isend(&cf_z[_PACK(cz, face_m_send, 0)], cz, nID[Z_MINUS], &key[_KEY_RQ(Z_MINUS, key_send)]) != CPM_SUCCESS ) Exit(0);
  }
  
  // recv
  if ( nID[Z_PLUS]>=0 )
  {
    if ( paraMngr->Irecv(&cf_z[_PACK(cz, face_p_recv, 0)], cz, nID[Z_PLUS], &key[_KEY_RQ(Z_PLUS, key_recv)]) != CPM_SUCCESS ) Exit(0);
  }
  
  
  // Z_PLUS
  // send
  if ( nID[Z_PLUS]>=0 )
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
    
    if ( paraMngr->Isend(&cf_z[_PACK(cz, face_p_send, 0)], cz, nID[Z_PLUS], &key[_KEY_RQ(Z_PLUS, key_send)]) != CPM_SUCCESS ) Exit(0);
  }
  
  // recv
  if ( nID[Z_MINUS]>=0 )
  {
    if ( paraMngr->Irecv(&cf_z[_PACK(cz, face_m_recv, 0)], cz, nID[Z_MINUS], &key[_KEY_RQ(Z_MINUS, key_recv)]) != CPM_SUCCESS ) Exit(0);
  }
  
}



// 種類Lの線形ソルバを利用する場合，trueを返す
bool FFV::hasLinearSolver(const int L)
{
  for (int i=0; i<ItrCtl::ic_END; i++)
    if ( IC[i].get_LS() == L ) return true;
  
  return false;
}


// 線形ソルバーの選択実行
void FFV::LS_Binary(ItrCtl* IC, REAL_TYPE b2)
{	
  double flop_count=0.0;         /// 浮動小数点演算数
  REAL_TYPE omg = IC->get_omg(); /// 加速係数
	REAL_TYPE res = 0.0;           /// 残差
  
	// d_p   圧力 p^{n+1}
  // d_p0  圧力 p^n
	// d_ws  非反復のソース項
  // d_sq  反復毎に変化するソース項
	// d_bcp ビットフラグ
  
  
  // 反復処理
  switch (IC->get_LS()) 
  {
    case SOR:
      TIMING_start(tm_poi_itr_sct_2); // >>> Poisson Iteration section 2
      
      // 反復処理
      TIMING_start(tm_poi_PSOR);
      flop_count = 0.0;
      psor_(d_p, size, &guide, &omg, &res, d_ws, d_sq, d_bcp, &flop_count);
      //r = PSOR(p, d_ws, d_sq, bcp, IC, flop_count); //実装速度比較
      TIMING_stop(tm_poi_PSOR, flop_count);
      
      // 境界条件
      TIMING_start(tm_poi_BC);
      BC.OuterPBC(d_p);
      if ( C.isPeriodic() == ON ) BC.InnerPBC_Periodic(d_p, d_bcd);
      TIMING_stop(tm_poi_BC, 0.0);
      
      TIMING_stop(tm_poi_itr_sct_2, 0.0); // <<< Poisson Iteration subsection 2
      
      
      TIMING_start(tm_poi_itr_sct_3); // >>> Poisson Iteration section 3
      
      // 同期処理
      if ( numProc > 1 ) 
      {
        TIMING_start(tm_poi_comm);
        if (IC->get_SyncMode() == comm_sync ) 
        {
          if ( paraMngr->BndCommS3D(d_p, size[0], size[1], size[2], guide, 1) != CPM_SUCCESS ) Exit(0); // 1 layer communication
        }
        else 
        {
          MPI_Request req[12];
          for (int i=0; i<12; i++) req[i] = MPI_REQUEST_NULL;
          
          if ( paraMngr->BndCommS3D_nowait(d_p, size[0], size[1], size[2], guide, 1, req ) != CPM_SUCCESS ) Exit(0); // 1 layer communication
          if ( paraMngr->wait_BndCommS3D  (d_p, size[0], size[1], size[2], guide, 1, req ) != CPM_SUCCESS ) Exit(0); // 1 layer communication
        }
        TIMING_stop(tm_poi_comm, comm_size);
      }
      
      // 残差の集約
      if ( numProc > 1 ) 
      {
        TIMING_start(tm_poi_res_comm);
        REAL_TYPE tmp = res;
        if ( paraMngr->Allreduce(&tmp, &res, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
        TIMING_stop(tm_poi_res_comm, 2.0*numProc*sizeof(REAL_TYPE) ); // 双方向 x ノード数
      }
      
      TIMING_stop(tm_poi_itr_sct_3, 0.0); // <<< Poisson Iteration subsection 3
      break;
      
  
    case SOR2SMA:
      res = SOR_2_SMA(IC);
      
      break;
      
    default:
      printf("\tInvalid Linear Solver for Pressure\n");
      Exit(0);
      break;
  }
  
  
  // 残差の保存 
  /// @note Control::p_res_l2_r/aのときのみ，それ以外はNorm_Poissonで計算
  if ( (IC->get_normType() == ItrCtl::p_res_l2_r) || (IC->get_normType() == ItrCtl::p_res_l2_a) ) 
  {
    REAL_TYPE res_a = sqrt( res /(REAL_TYPE)G_Acell ); // 残差のRMS
    REAL_TYPE bb  = sqrt( b2/(REAL_TYPE)G_Acell ); // ソースベクトルのRMS
    REAL_TYPE res_r = ( bb == 0.0 ) ? res_a : res_a/bb;
    REAL_TYPE nrm = ( IC->get_normType() == ItrCtl::p_res_l2_r ) ? res_r : res_a;
    IC->set_normValue( nrm );
  }
}



// SOR2SMA
REAL_TYPE FFV::SOR_2_SMA(ItrCtl* IC)
{
  int ip;                        /// ローカルノードの基点(1,1,1)のカラーを示すインデクス
                                 /// ip=0 > R, ip=1 > B
  double flop_count=0.0;         /// 浮動小数点演算数
  double comm_size;              /// 通信面1面あたりの通信量
  REAL_TYPE omg = IC->get_omg(); /// 加速係数
	REAL_TYPE res = 0.0;           /// 残差
  
  comm_size = count_comm_size(size, guide);
  
  
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
    TIMING_start(tm_poi_itr_sct_2); // >>> Poisson Iteration section 2
    
    TIMING_start(tm_poi_SOR2SMA);
    flop_count = 0.0; // 色間で積算しない
    psor2sma_core_(d_p, size, &guide, &ip, &color, &omg, &res, d_ws, d_sq, d_bcp, &flop_count);
    //PSOR2sma_core(p, ip, color, d_ws, d_sq, bcp, IC, flop_count);
    TIMING_stop(tm_poi_SOR2SMA, flop_count);
    
    // 境界条件
    TIMING_start(tm_poi_BC);
    BC.OuterPBC(d_p);
    if ( C.isPeriodic() == ON ) BC.InnerPBC_Periodic(d_p, d_bcd);
    TIMING_stop(tm_poi_BC, 0.0);
    
    TIMING_stop(tm_poi_itr_sct_2, 0.0); // <<< Poisson Iteration subsection 2
    
    
    TIMING_start(tm_poi_itr_sct_3); // >>> Poisson Iteration section 3
    
    // 同期処理
    if ( numProc > 1 )
    {
      TIMING_start(tm_poi_comm);
      if (IC->get_SyncMode() == comm_sync )
      {
        if ( paraMngr->BndCommS3D(d_p, size[0], size[1], size[2], guide, 1) != CPM_SUCCESS ) Exit(0); // 1 layer communication
      }
      else
      {
        //if ( paraMngr->BndCommS3D_nowait(d_p, size[0], size[1], size[2], guide, 1, req) != CPM_SUCCESS ) Exit(0);
        //if ( paraMngr->wait_BndCommS3D  (d_p, size[0], size[1], size[2], guide, 1, req) != CPM_SUCCESS ) Exit(0);
        //comm_SOR2SMA(color, ip, req);
        //wait_SOR2SMA(color, ip, req);
        int ireq[12];
        sma_comm_     (d_p, size, &guide, &color, &ip, cf_sz, cf_x, cf_y, cf_z, ireq, nID);
        sma_comm_wait_(d_p, size, &guide, &color, &ip, cf_sz, cf_x, cf_y, cf_z, ireq);
      }
      TIMING_stop(tm_poi_comm, comm_size*0.5);
    }
    
    TIMING_stop(tm_poi_itr_sct_3, 0.0); // <<< Poisson Iteration subsection 3
  }
  
  // 残差の集約
  if ( numProc > 1 )
  {
    TIMING_start(tm_poi_res_comm);
    REAL_TYPE tmp = res;
    if ( paraMngr->Allreduce(&tmp, &res, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
    TIMING_stop(tm_poi_res_comm, 2.0*numProc*sizeof(REAL_TYPE)*0.5 ); // 双方向 x ノード数 check
  }
  
  return res;
}



// SOR2SMAの非同期通信処理
// 並列時のみコールされる
void FFV::wait_SOR2SMA(const int col, const int ip, MPI_Request* key)
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
  rq = _KEY_RQ(X_PLUS, key_recv);
  
  if (key[rq]>=0 )
  {
    if ( paraMngr->Wait(&key[rq]) != CPM_SUCCESS ) Exit(0);
    
    int c = 0;
    int i = ix+1;
    
    for (int k=1; k<=kx; k++) {
      int js= (k + ic + i) % 2;
      
      for (int j=1+js; j<=jx; j+=2) {
        d_p[ _F_IDX_S3D(i, j, k, ix, jx, kx, gd) ] = cf_x[_PACK(cx, face_p_recv, c)];
        c++;
      }
    }
  }
  
  // X_MINUS面の受信バッファを展開
  rq = _KEY_RQ(X_MINUS, key_recv);
  
  if (key[rq]>=0 )
  {
    if ( paraMngr->Wait(&key[rq]) != CPM_SUCCESS ) Exit(0);
    
    int c = 0;
    int i = 0;
    
    for (int k=1; k<=kx; k++) {
      int js= (k + ic + i) % 2;
      
      for (int j=1+js; j<=jx; j+=2) {
        d_p[ _F_IDX_S3D(i, j, k, ix, jx, kx, gd) ] = cf_x[_PACK(cx, face_m_recv, c)];
        c++;
      }
    }
  }
          
  // Y_PLUS面の受信バッファを展開
  rq = _KEY_RQ(Y_PLUS, key_recv);
  
  if (key[rq]>=0 )
  {
    if ( paraMngr->Wait(&key[rq]) != CPM_SUCCESS ) Exit(0);
    
    int c = 0;
    int j = jx+1;
    
    for (int k=1; k<=kx; k++) {
      int is= (k + ic + j) % 2;
      
      for (int i=1+is; i<=ix; i+=2) {
        d_p[ _F_IDX_S3D(i, j, k, ix, jx, kx, gd) ] = cf_y[_PACK(cy, face_p_recv, c)];
        c++;
      }
    }
  }
  
  // Y_MINUS面の受信バッファを展開
  rq = _KEY_RQ(Y_MINUS, key_recv);
  
  if (key[rq]>=0 )
  {
    if ( paraMngr->Wait(&key[rq]) != CPM_SUCCESS ) Exit(0);
    
    int c = 0;
    int j = 0;
    
    for (int k=1; k<=kx; k++) {
      int is= (k + ic + j) % 2;
      
      for (int i=1+is; i<=ix; i+=2) {
        d_p[ _F_IDX_S3D(i, j, k, ix, jx, kx, gd) ] = cf_y[_PACK(cy, face_m_recv, c)];
        c++;
      }
    }
  }
  
  // Z_PLUS面の受信バッファを展開
  rq = _KEY_RQ(Z_PLUS, key_recv);
  
  if (key[rq]>=0 )
  {
    if ( paraMngr->Wait(&key[rq]) != CPM_SUCCESS ) Exit(0);
    
    int c = 0;
    int k = kx+1;
    
    for (int j=1; j<=jx; j++) {
      int is= (j + ic + k) % 2;
      
      for (int i=1+is; i<=ix; i+=2) {
        d_p[ _F_IDX_S3D(i, j, k, ix, jx, kx, gd) ] = cf_z[_PACK(cz, face_p_recv, c)];
        c++;
      }
    }
  }
  
  // Z_MINUS面の受信バッファを展開
  rq = _KEY_RQ(Z_MINUS, key_recv);
  
  if (key[rq]>=0 )
  {
    if ( paraMngr->Wait(&key[rq]) != CPM_SUCCESS ) Exit(0);
    
    int c = 0;
    int k = 0;
    
    for (int j=1; j<=jx; j++) {
      int is= (j + ic + k) % 2;
      
      for (int i=1+is; i<=ix; i+=2) {
        d_p[ _F_IDX_S3D(i, j, k, ix, jx, kx, gd) ] = cf_z[_PACK(cz, face_m_recv, c)];
        c++;
      }
    }
  }
  
  
  // Wait for send ------------------------
  
  // X_MINUS
  rq = _KEY_RQ(X_MINUS, key_send);
  
  if (key[rq]>=0 )
  {
    if ( paraMngr->Wait(&key[rq]) != CPM_SUCCESS ) Exit(0);
  }
  
  // X_PLUS
  rq = _KEY_RQ(X_PLUS, key_send);
  
  if (key[rq]>=0 )
  {
    if ( paraMngr->Wait(&key[rq]) != CPM_SUCCESS ) Exit(0);
  }
  
  // Y_MINUS
  rq = _KEY_RQ(Y_MINUS, key_send);
  
  if (key[rq]>=0 )
  {
    if ( paraMngr->Wait(&key[rq]) != CPM_SUCCESS ) Exit(0);
  }
  
  // Y_PLUS
  rq = _KEY_RQ(Y_PLUS, key_send);
  
  if (key[rq]>=0 )
  {
    if ( paraMngr->Wait(&key[rq]) != CPM_SUCCESS ) Exit(0);
  }
  
  // Z_MINUS
  rq = _KEY_RQ(Z_MINUS, key_send);
  
  if (key[rq]>=0 )
  {
    if ( paraMngr->Wait(&key[rq]) != CPM_SUCCESS ) Exit(0);
  }
  
  // Z_PLUS
  rq = _KEY_RQ(Z_PLUS, key_send);
  
  if (key[rq]>=0 )
  {
    if ( paraMngr->Wait(&key[rq]) != CPM_SUCCESS ) Exit(0);
  }
}
