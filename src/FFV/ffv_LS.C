// #################################################################
//
// FFV : Frontflow / violet Cartesian
//
// Copyright (c) All right reserved. 2012
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
void FFV::comm_SOR2SMA(const int col, const int ip, int* key)
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
  
  int a = col + ip;
  int ic = a - int(a/2)*2; // スタートインデクス
  
  // 送信IDを初期化
  for (int i=0; i<NOFACE*2; i++) key[i] = -1;
  

  // X_MINUS
  // send
  if ( nID[X_MINUS]>=0 )
  {
    int c = 0;
    int i = 1;
    
    for (int k=1; k<=kx; k++) {
      int b = k + ic + i;
      int js= b - int(b/2)*2;
      
      for (int j=1+js; j<=jx; j+=2) {
        cf_x[_PACK(cx, send_to_minus, c)] = d_p[ _F_IDX_S3D(i, j, k, ix, jx, kx, gd) ];
        c++;
      }
    }
    
    if ( paraMngr->Isend(cf_x[_PACK(cx, send_to_minus, 0)], cx, nID[X_MINUS], key[_ASYNC_RQ(X_MINUS, async_send)]) != CPM_SUCCESS ) Exit(0);
  }
    
  // recieve
  if ( nID[X_PLUS]>=0 )
  {
    if ( paraMngr->Irecv(cf_x[_PACK(cx, recv_from_plus, 0)], cx, nID[X_PLUS], key[_ASYNC_RQ(X_MINUS, async_recv)]) != CPM_SUCCESS ) Exit(0);
  }
    
    
  // X_PLUS
  // send
  if ( nID[X_PLUS]>=0 )
  {
    int c = 0;
    int i = ix;
    
    for (int k=1; k<=kx; k++) {
      int b = k + ic + i;
      int js= b - int(b/2)*2;
      
      for (int j=1+js; j<=jx; j+=2) {
        cf_x[_PACK(cx, send_to_plus, c)] = d_p[ _F_IDX_S3D(i, j, k, ix, jx, kx, gd) ];
        c++;
      }
    }
    
    if ( paraMngr->Isend(cf_x[_PACK(cx, send_to_plus, 0)], cx, nID[X_PLUS], key[_ASYNC_RQ(X_PLUS, async_send)]) != CPM_SUCCESS ) Exit(0);
  }
  
  // recieve
  if ( nID[X_MINUS]>=0 )
  {
    if ( paraMngr->Irecv(cf_x[_PACK(cx, recv_from_minus, 0)], cx, nID[X_MINUS], key[_ASYNC_RQ(X_PLUS, async_recv)]) != CPM_SUCCESS ) Exit(0);
  }
  
  
  // Y_MINUS
  // send
  if ( nID[Y_MINUS]>=0 )
  {
    int c = 0;
    int j = 1;
    
    for (int k=1; k<=kx; k++) {
      int b = k + ic + j;
      int is= b - int(b/2)*2;
      
      for (int i=1+is; i<=ix; i+=2) {
        cf_y[_PACK(cy, send_to_minus, c)] = d_p[ _F_IDX_S3D(i, j, k, ix, jx, kx, gd) ];
        c++;
      }
    }
    
    if ( paraMngr->Isend(cf_y[_PACK(cy, send_to_minus, 0)], cy, nID[Y_MINUS], key[_ASYNC_RQ(Y_MINUS, async_send)]) != CPM_SUCCESS ) Exit(0);
  }
  
  // recv
  if ( nID[Y_PLUS]>=0 )
  {
    if ( paraMngr->Irecv(cf_y[_PACK(cy, recv_from_plus, 0)], cy, nID[Y_PLUS], key[_ASYNC_RQ(Y_MINUS, async_recv)]) != CPM_SUCCESS ) Exit(0);
  }
  
  
  // Y_PLUS
  // send
  if ( nID[Y_PLUS]>=0 )
  {
    int c = 0;
    int j = jx;
    
    for (int k=1; k<=kx; k++) {
      int b = k + ic + j;
      int is= b - int(b/2)*2;
      
      for (int i=1+is; i<=ix; i+=2) {
        cf_y[_PACK(cy, send_to_plus, c)] = d_p[ _F_IDX_S3D(i, j, k, ix, jx, kx, gd) ];
        c++;
      }
    }
    
    if ( paraMngr->Isend(cf_y[_PACK(cy, send_to_plus, 0)], cy, nID[Y_PLUS], key[_ASYNC_RQ(Y_PLUS, async_send)]) != CPM_SUCCESS ) Exit(0);
  }
  
  // recv
  if ( nID[Y_MINUS]>=0 )
  {
    if ( paraMngr->Irecv(cf_y[_PACK(cy, recv_from_minus, 0)], cy, nID[Y_MINUS], key[_ASYNC_RQ(Y_PLUS, async_recv)]) != CPM_SUCCESS ) Exit(0);
  }
  
  
  // Z_MINUS
  // send
  if ( nID[Z_MINUS]>=0 )
  {
    int c = 0;
    int k = 1;
    
    for (int j=1; j<=jx; j++) {
      int b = j + ic + k;
      int is= b - int(b/2)*2;
      
      for (int i=1+is; i<=ix; i+=2) {
        cf_z[_PACK(cz, send_to_minus, c)] = d_p[ _F_IDX_S3D(i, j, k, ix, jx, kx, gd) ];
        c++;
      }
    }
    
    if ( paraMngr->Isend(cf_z[_PACK(cz, send_to_minus, 0)], cz, nID[Z_MINUS], key[_ASYNC_RQ(Z_MINUS, async_send)]) != CPM_SUCCESS ) Exit(0);
  }
  
  // recv
  if ( nID[Z_PLUS]>=0 )
  {
    if ( paraMngr->Irecv(cf_z[_PACK(cz, recv_from_plus, 0)], cz, nID[Z_PLUS], key[_ASYNC_RQ(Z_MINUS, async_recv)]) != CPM_SUCCESS ) Exit(0);
  }
  
  
  // Z_PLUS
  // send
  if ( nID[Z_PLUS]>=0 )
  {
    int c = 0;
    int k = kx;
    
    for (int j=1; j<=jx; j++) {
      int b = j + ic + k;
      int is= b - int(b/2)*2;
      
      for (int i=1+is; i<=ix; i+=2) {
        cf_z[_PACK(cz, send_to_plus, c)] = d_p[ _F_IDX_S3D(i, j, k, ix, jx, kx, gd) ];
        c++;
      }
    }
    
    if ( paraMngr->Isend(cf_z[_PACK(cz, send_to_plus, 0)], cz, nID[Z_PLUS], key[_ASYNC_RQ(Z_PLUS, async_send)]) != CPM_SUCCESS ) Exit(0);
  }
  
  // recv
  if ( nID[Z_MINUS]>=0 )
  {
    if ( paraMngr->Irecv(cf_z[_PACK(cz, recv_from_minus, 0)], cz, nID[Z_MINUS], key[_ASYNC_RQ(Z_PLUS, async_recv)]) != CPM_SUCCESS ) Exit(0);
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
  double comm_size;              /// 通信面1面あたりの通信量
  REAL_TYPE clear_value=0.0;	   /// 初期化の値
  REAL_TYPE omg = IC->get_omg(); /// 加速係数
	REAL_TYPE r = 0.0;             /// 残差
  
  comm_size = count_comm_size(size, guide);
  
	// c_p   圧力 p^{n+1}
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
      psor_(d_p, size, &guide, &omg, &r, d_ws, d_sq, d_bcp, &flop_count);
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
        if ( paraMngr->BndCommS3D(d_p, size[0], size[1], size[2], guide, 1) != CPM_SUCCESS ) Exit(0); // 1 layer communication
        TIMING_stop(tm_poi_comm, comm_size);
      }
      
      // 残差の集約
      if ( numProc > 1 ) 
      {
        TIMING_start(tm_poi_res_comm);
        REAL_TYPE tmp = r;
        if ( paraMngr->Allreduce(&tmp, &r, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
        TIMING_stop(tm_poi_res_comm, 2.0*numProc*sizeof(REAL_TYPE) ); // 双方向 x ノード数
      }
      
      TIMING_stop(tm_poi_itr_sct_3, 0.0); // <<< Poisson Iteration subsection 3
      break;
      
      /*
       case SOR2SMA:
       // 2色のマルチカラーのセットアップ
       TIMING_start(tm_poi_setup);
       
       int nID[6], sidx[3], ip;
       
       if ( numProc > 1 ) {
       sidx[0] = head[0];
       sidx[1] = head[1];
       sidx[2] = head[2];
       ip = (sidx[0]+sidx[1]+sidx[2]) % 2;
       } 
       else {
       sidx[0] = sidx[1] = sidx[2] = 1;
       ip = 0;
       }
       TIMING_stop(tm_poi_setup, 0.0);
       
       
       // 各カラー毎の間に同期
       r = 0.0;          // 色間で積算する
       for (int color=0; color<2; color++) 
       {
       TIMING_start(tm_poi_itr_sct_2); // >>> Poisson Iteration section 2
       
       TIMING_start(tm_poi_SOR2SMA);
       flop_count = 0.0; // 色間で積算しない
       psor2sma_core_(d_p, size, &guide, &ip, &color, &omg, &r, d_ws, d_sq, d_bcp, &flop_count);
       //PSOR2sma_core(p, ip, color, d_ws, d_sq, bcp, IC, flop_count);
       TIMING_stop(tm_poi_SOR2SMA, flop_count);
       
       // 境界条件
       TIMING_start(tm_poi_BC);
       BC.OuterPBC(d_p);
       if ( C.isPeriodic() == ON ) BC.InnerPBC_Periodic(d_p, d_bcd);
       TIMING_stop(tm_poi_BC, 0.0);
       
       TIMING_stop(tm_poi_itr_sct_2, 0.0); // <<< Poisson Iteration subsection 2
       
       
       TIMING_start(tm_poi_itr_sct_3); // >>> Poisson Iteration section 3
       
       // 残差の集約と同期処理
       if ( numProc > 1 ) {
       TIMING_start(tm_poi_res_comm);
       REAL_TYPE tmp = r;
       if ( paraMngr->Allreduce(&tmp, &r, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
       TIMING_stop(tm_poi_res_comm, 2.0*numProc*sizeof(REAL_TYPE)*0.5 ); // 双方向 x ノード数 check
       
       TIMING_start(tm_poi_comm);
       if (cm_mode == 0 ) 
       {
       if ( paraMngr->BndCommS3D(d_p, size[0], size[1], size[2], guide, 1) != CPM_SUCCESS ) Exit(0); // 1 layer communication
       }
       else {
       sma_comm_     (d_p, size, &guide, &color, &ip, cf_sz, cf_x, cf_y, cf_z, req, &para_key);
       sma_comm_wait_(d_p, size, &guide, &color, &ip, cf_sz, cf_x, cf_y, cf_z, req);
       }
       TIMING_stop(tm_poi_comm, comm_size*0.5);
       }
       
       TIMING_stop(tm_poi_itr_sct_3, 0.0); // <<< Poisson Iteration subsection 3
       }
       break; */
      
    default:
      printf("\tInvalid Linear Solver for Pressure\n");
      Exit(0);
      break;
  }
  
  // 残差の保存 
  /// @note Control::p_res_l2_r/aのときのみ，それ以外はNorm_Poissonで計算
  if ( (IC->get_normType() == ItrCtl::p_res_l2_r) || (IC->get_normType() == ItrCtl::p_res_l2_a) ) 
  {
    REAL_TYPE res_a = sqrt( r /(REAL_TYPE)G_Acell ); // 残差のRMS
    REAL_TYPE bb  = sqrt( b2/(REAL_TYPE)G_Acell ); // ソースベクトルのRMS
    REAL_TYPE res_r = ( bb == 0.0 ) ? res_a : res_a/bb;
    REAL_TYPE nrm = ( IC->get_normType() == ItrCtl::p_res_l2_r ) ? res_r : res_a;
    IC->set_normValue( nrm );
  }
}



// SOR2SMAの非同期通信処理
// 並列時のみコールされる
void FFV::wait_SOR2SMA(const int col, const int ip, int* key)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  int cx = cf_sz[0];
  int cy = cf_sz[1];
  int cz = cf_sz[2];
  
  int a = col + ip;
  int ic = a - int(a/2)*2; // スタートインデクス
  
  int rq;
  
  // Wait for recv ------------------------
  
  // from X_MINUS
  rq = _ASYNC_RQ(X_MINUS, async_recv);
  
  if (key[rq]>=0 )
  {
    if ( paraMngr->Wait(key[rq]) != CPM_SUCCESS ) Exit(0);
    
    int c = 0;
    int i = ix+1;
    
    for (int k=1; k<=kx; k++) {
      int b = k + ic + i;
      int js= b - int(b/2)*2;
      
      for (int j=1+js; j<=jx; j+=2) {
        d_p[ _F_IDX_S3D(i, j, k, ix, jx, kx, gd) ] = cf_x[_PACK(cx, recv_from_plus, c)];
        c++;
      }
    }
  }
  
  // from X_PLUS
  rq = _ASYNC_RQ(X_PLUS, async_recv);
  
  if (key[rq]>=0 )
  {
    if ( paraMngr->Wait(key[rq]) != CPM_SUCCESS ) Exit(0);
    
    int c = 0;
    int i = 0;
    
    for (int k=1; k<=kx; k++) {
      int b = k + ic + i;
      int js= b - int(b/2)*2;
      
      for (int j=1+js; j<=jx; j+=2) {
        d_p[ _F_IDX_S3D(i, j, k, ix, jx, kx, gd) ] = cf_x[_PACK(cx, recv_from_minus, c)];
        c++;
      }
    }
  }
          
  // from Y_MINUS
  rq = _ASYNC_RQ(Y_MINUS, async_recv);
  
  if (key[rq]>=0 )
  {
    if ( paraMngr->Wait(key[rq]) != CPM_SUCCESS ) Exit(0);
    
    int c = 0;
    int j = jx+1;
    
    for (int k=1; k<=kx; k++) {
      int b = k + ic + j;
      int is= b - int(b/2)*2;
      
      for (int i=1+is; i<=ix; i+=2) {
        d_p[ _F_IDX_S3D(i, j, k, ix, jx, kx, gd) ] = cf_y[_PACK(cy, recv_from_plus, c)];
        c++;
      }
    }
  }
  
  // from Y_PLUS
  rq = _ASYNC_RQ(Y_PLUS, async_recv);
  
  if (key[rq]>=0 )
  {
    if ( paraMngr->Wait(key[rq]) != CPM_SUCCESS ) Exit(0);
    
    int c = 0;
    int j = 0;
    
    for (int k=1; k<=kx; k++) {
      int b = k + ic + j;
      int is= b - int(b/2)*2;
      
      for (int i=1+is; i<=ix; i+=2) {
        d_p[ _F_IDX_S3D(i, j, k, ix, jx, kx, gd) ] = cf_y[_PACK(cy, recv_from_minus, c)];
        c++;
      }
    }
  }
  
  // from Z_MINUS
  rq = _ASYNC_RQ(Z_MINUS, async_recv);
  
  if (key[rq]>=0 )
  {
    if ( paraMngr->Wait(key[rq]) != CPM_SUCCESS ) Exit(0);
    
    int c = 0;
    int k = kx+1;
    
    for (int j=1; j<=jx; j++) {
      int b = j + ic + k;
      int is= b - int(b/2)*2;
      
      for (int i=1+is; i<=ix; i+=2) {
        d_p[ _F_IDX_S3D(i, j, k, ix, jx, kx, gd) ] = cf_z[_PACK(cz, recv_from_plus, c)];
        c++;
      }
    }
  }
  
  // from Z_PLUS
  rq = _ASYNC_RQ(Z_PLUS, async_recv);
  
  if (key[rq]>=0 )
  {
    if ( paraMngr->Wait(key[rq]) != CPM_SUCCESS ) Exit(0);
    
    int c = 0;
    int k = 0;
    
    for (int j=1; j<=jx; j++) {
      int b = j + ic + k;
      int is= b - int(b/2)*2;
      
      for (int i=1+is; i<=ix; i+=2) {
        d_p[ _F_IDX_S3D(i, j, k, ix, jx, kx, gd) ] = cf_z[_PACK(cz, recv_from_minus, c)];
        c++;
      }
    }
  }
  
  
  // Wait for send ------------------------
  
  // X_MINUS
  rq = _ASYNC_RQ(X_MINUS, async_send);
  
  if (key[rq]>=0 )
  {
    if ( paraMngr->Wait(key[rq]) != CPM_SUCCESS ) Exit(0);
  }
  
  // X_PLUS
  rq = _ASYNC_RQ(X_PLUS, async_send);
  
  if (key[rq]>=0 )
  {
    if ( paraMngr->Wait(key[rq]) != CPM_SUCCESS ) Exit(0);
  }
  
  // Y_MINUS
  rq = _ASYNC_RQ(Y_MINUS, async_send);
  
  if (key[rq]>=0 )
  {
    if ( paraMngr->Wait(key[rq]) != CPM_SUCCESS ) Exit(0);
  }
  
  // Y_PLUS
  rq = _ASYNC_RQ(Y_PLUS, async_send);
  
  if (key[rq]>=0 )
  {
    if ( paraMngr->Wait(key[rq]) != CPM_SUCCESS ) Exit(0);
  }
  
  // Z_MINUS
  rq = _ASYNC_RQ(Z_MINUS, async_send);
  
  if (key[rq]>=0 )
  {
    if ( paraMngr->Wait(key[rq]) != CPM_SUCCESS ) Exit(0);
  }
  
  // Z_PLUS
  rq = _ASYNC_RQ(Z_PLUS, async_send);
  
  if (key[rq]>=0 )
  {
    if ( paraMngr->Wait(key[rq]) != CPM_SUCCESS ) Exit(0);
  }
}
