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
