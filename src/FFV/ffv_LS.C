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


// #################################################################
// 種類Lの線形ソルバを利用する場合，trueを返す
bool FFV::hasLinearSolver(const int L)
{
  for (int i=0; i<ItrCtl::ic_END; i++)
    if ( IC[i].get_LS() == L ) return true;
  
  return false;
}


// #################################################################
// 線形ソルバーの選択実行
void FFV::LS_Binary(ItrCtl* IC, const double rhs_nrm)
{	
	double res = 0.0;    /// 残差

  
  // 反復処理
  switch (IC->get_LS()) 
  {
    case SOR:
      res = Point_SOR(IC);
      break;
      
    case SOR2SMA:
      res = SOR_2_SMA(IC);
      break;
      
    case GMRES_SOR:
      res = Gmres_SOR(IC, rhs_nrm);
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
    double res_a = sqrt( res / (double)G_Acell ); // 残差のRMS
    double bb  = sqrt( rhs_nrm / (double)G_Acell ); // ソースベクトルのRMS
    double res_r = ( bb == 0.0 ) ? res_a : res_a/bb;
    double nrm = ( IC->get_normType() == ItrCtl::p_res_l2_r ) ? res_r : res_a;
    IC->set_normValue( nrm );
  }
}


// #################################################################
// Point SOR
double FFV::Point_SOR(ItrCtl* IC)
{
  double flop_count=0.0;         /// 浮動小数点演算数
  REAL_TYPE omg = IC->get_omg(); /// 加速係数
	double res = 0.0;              /// 残差
  
  // d_p   圧力 p^{n+1}
  // d_p0  圧力 p^n
	// d_ws  非反復のソース項
  // d_sq  反復毎に変化するソース項
	// d_bcp ビットフラグ
  
  
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
  Sync_Scalar(IC, d_p, 1);
  
  TIMING_stop(tm_poi_itr_sct_3, 0.0); // <<< Poisson Iteration subsection 3
  
  
  // 残差の集約
  if ( numProc > 1 )
  {
    TIMING_start(tm_poi_res_comm);
    double tmp = res;
    if ( paraMngr->Allreduce(&tmp, &res, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
    TIMING_stop(tm_poi_res_comm, 2.0*numProc*sizeof(double) ); // 双方向 x ノード数
  }

  return res;
}


// #################################################################
// 反復変数の同期処理
void FFV::Sync_Scalar(ItrCtl* IC, REAL_TYPE* d_class, const int num_layer)
{
  if ( numProc > 1 )
  {
    TIMING_start(tm_poi_comm);
    
    /// 通信面1面あたりの通信量
    double comm_size = count_comm_size(size, guide);
    
    if (IC->get_SyncMode() == comm_sync )
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
    TIMING_stop(tm_poi_comm, comm_size*(double)num_layer);
  }
}



// #################################################################
// SOR2SMA
double FFV::SOR_2_SMA(ItrCtl* IC)
{
  int ip;                        /// ローカルノードの基点(1,1,1)のカラーを示すインデクス
                                 /// ip=0 > R, ip=1 > B
  double flop_count=0.0;         /// 浮動小数点演算数
  REAL_TYPE omg = IC->get_omg(); /// 加速係数
	double res = 0.0;              /// 残差
  
  
  // d_p   圧力 p^{n+1}
  // d_p0  圧力 p^n
	// d_ws  非反復のソース項
  // d_sq  反復毎に変化するソース項
	// d_bcp ビットフラグ
  
  
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
      
      /// 通信面1面あたりの通信量
      double comm_size = count_comm_size(size, guide);
      
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
    double tmp = res;
    if ( paraMngr->Allreduce(&tmp, &res, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
    TIMING_stop(tm_poi_res_comm, 2.0*numProc*sizeof(REAL_TYPE)*0.5 ); // 双方向 x ノード数 check
  }
  
  return res;
}


// #################################################################
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

// #################################################################
// Flexible gmres(m)
double FFV::Fgmres(ItrCtl* IC, double res_rhs)
{
  const double fct2  = 6.0;
  const double eps_1 = 1.0e-30;
  const double eps_2 = IC->get_eps();
  
  // 残差収束チェック
  if ( res_rhs < eps_1 ) return res_rhs;
  
  
  // >>> Gmres section
  TIMING_start(tm_gmres_sor_sct);
  
  double t_eps, beta, beta_1, res, eps_abs, res_abs, al;
  double r4;
  double flop=0.0;
  
  const int Iteration_Max = IC->get_ItrMax();
  int m_max = FREQ_OF_RESTART;
  int s_length = (size[0]+2*guide) * (size[1]+2*guide) * (size[2]+2*guide);
  
  int isfin  = 0;
  
  const int oki  = 2;
  const int step = 3;
  int nrm  = 0;
  
  if ( C.Mode.Precision == sizeof(float) )
  {
    t_eps = 0.995 * eps_2 * eps_2;
  }
  else
  {
    t_eps = 0.999 * eps_2 * eps_2;
  }
  
  const int Nmax = m_max+1; // 配列確保用，ゼロはダミー
  
  double *rgm = new double[Nmax * Nmax];     // (Nmax, Nmax)
  double *hgm = new double[(Nmax+1) * Nmax]; // (Nmax+1, Nmax)
  double *cgm = new double[Nmax];
  double *sgm = new double[Nmax];
  double *bgm = new double[Nmax];
  double *ygm = new double[Nmax];
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  int column = 1;
  
  t_eps = res_rhs * t_eps;
  
  res = SOR_2_SMA(IC);
  
  r4 = sqrt(res/res_rhs);
  
  
  if (res < t_eps)
  {
    //Hostonly_ printf("Final : %e\n", r4);
    isfin = 1;
    goto jump_4;
  }
  
  
  if ( res > (fct2*t_eps) ) t_eps = res/fct2;
  
  eps_abs = 0.15 * t_eps;
  
  
  for (int i_iter=1; i_iter<=Iteration_Max; i_iter++) {
    
    // テスト的にはずす >> 収束性良い？
    //res = SOR_2_SMA(IC);
    //if (res < t_eps) goto jump_3;
    
    TIMING_start(tm_gmres_mvprod);
    flop = 0.0;
    mv_prod_(d_yt, size, &guide, d_p, d_bcp, &flop);
    TIMING_stop(tm_gmres_mvprod, flop);
    
    
    TIMING_start(tm_gmres_res_sample);
    flop = 0.0;
    res_smpl_(d_rest, size, &guide, &res_abs, d_ws, d_yt, d_bcp, &flop);
    TIMING_stop(tm_gmres_res_sample, flop);
    
    
    TIMING_start(tm_gmres_comm);
    if ( numProc > 1 )
    {
      double tmp = res_abs;
      if ( paraMngr->Allreduce(&tmp, &res_abs, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
    }
    TIMING_stop(tm_gmres_comm, 2.0*numProc*sizeof(double)); // 双方向 x ノード数
    
    
    beta = sqrt(res_abs);
    
    if (beta < eps_1) goto jump_2;
    
    
    beta_1 = 1.0/beta;
    
    multiply_(d_vm, size, &guide, &m_max, &column, &beta_1, d_rest, &flop);
    
    
    bgm[1] = beta;
    
    for (int i=2; i<=m_max; i++) {
      bgm[i] = 0.0;
    }
    
    
    
    for (int im=1; im<=m_max; im++) {
      
      copy_1_(d_rest, size, &guide, &m_max, d_vm, &im);
      
      
#pragma omp parallel for firstprivate(s_length)
      for (int i=0; i<s_length; i++) {
        d_xt[i] = 0.0;
      }
      
      
      // Decide the number of iteration for inner-iteration
      int n_inner = 3; //((im / oki) + 1) * step;
      
      
      // Inner-iteration
      for (int i_inner=1; i_inner<=n_inner; i_inner++) {
        
        res = SOR_2_SMA(IC);
      }
      
      copy_2_(d_zm, size, &guide, &m_max, d_xt, &im);
      
      TIMING_start(tm_gmres_mvprod);
      flop = 0.0;
      mv_prod_(d_xt, size, &guide, d_p, d_bcp, &flop);
      TIMING_stop(tm_gmres_mvprod, flop);
      
      
      for (int km=1; km<=im; km++) {
        al = 0.0;
        ml_add_1_(&al, size, &guide, &m_max, d_vm, d_yt, &km, &flop);
        
        
        TIMING_start(tm_gmres_comm);
        if ( numProc > 1 )
        {
          double tmp = al;
          if ( paraMngr->Allreduce(&tmp, &al, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
        }
        TIMING_stop(tm_gmres_comm, 2.0*numProc*sizeof(double)); // 双方向 x ノード数
        
        
        hgm[_IDX2D(km, im, Nmax)] = al;
        
        double r_al = -al;
        ml_add_3_(d_yt, size, &guide, &m_max, &r_al, d_vm, &km, &flop);
      }
      
      
      al =0.0;
      ml_add_2_(&al, size, &guide, d_yt, &flop);
      
      
      TIMING_start(tm_gmres_comm);
      if ( numProc > 1 )
      {
        double tmp = al;
        if ( paraMngr->Allreduce(&tmp, &al, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
      }
      TIMING_stop(tm_gmres_comm, 2.0*numProc*sizeof(double)); // 双方向 x ノード数
      
      
      TIMING_start(tm_gmres_others);
      flop = 0.0;
      
      hgm[_IDX2D(im+1, im, Nmax)] = sqrt(al);
      flop += 20.0;
      
      if ( hgm[_IDX2D(im+1, im, Nmax)] < eps_1 )
      {
        nrm = im-1;
        TIMING_stop(tm_gmres_others, flop);
        goto jump_1;
      }
      
      if (im < m_max)
      {
        int idx = im + 1;
        al = 1.0 / hgm[_IDX2D(idx, im, Nmax)];
        flop += 13.0;
        multiply_(d_vm, size, &guide, &m_max, &idx, &al, d_yt, &flop);
      }
      
      rgm[_IDX2D(1, im, Nmax)] = hgm[_IDX2D(1, im, Nmax)];
      
      for (int km=2; km<=im; km++) {
        rgm[_IDX2D(km  , im, Nmax)] = cgm[km-1]*hgm[_IDX2D(km, im, Nmax)] - sgm[km-1]*rgm[_IDX2D(km-1, im, Nmax)];
        rgm[_IDX2D(km-1, im, Nmax)] = sgm[km-1]*hgm[_IDX2D(km, im, Nmax)] + cgm[km-1]*rgm[_IDX2D(km-1, im, Nmax)];
      }
      flop += (double)(im-2+1)*6.0;
      
      al = sqrt(rgm[_IDX2D(im, im, Nmax)]*rgm[_IDX2D(im, im, Nmax)] + hgm[_IDX2D(im+1, im, Nmax)]*hgm[_IDX2D(im+1, im, Nmax)]);
      flop += 23;
      
      if ( al < eps_1 )
      {
        nrm = im - 1;
        TIMING_stop(tm_gmres_others, flop);
        goto jump_1;
      }
      
      cgm[im]    = rgm[_IDX2D(im,  im, Nmax)] / al;
      sgm[im]    = hgm[_IDX2D(im+1,im, Nmax)] / al;
      
      rgm[_IDX2D(im, im, Nmax)] = cgm[im] * rgm[_IDX2D(im, im, Nmax)] + sgm[im] * hgm[_IDX2D(im+1, im, Nmax)];
      
      bgm[im+1]  = - sgm[im]*bgm[im];
      bgm[im]    =   cgm[im]*bgm[im];
      
      al = bgm[im+1] * bgm[im+1];
      
      flop += 2.0*13.0 + 4.0;
      
      if (al < eps_abs)
      {
        nrm = im;
        TIMING_stop(tm_gmres_others, flop);
        goto jump_1;
      }
      
      TIMING_stop(tm_gmres_others, flop);
      
    } // loop; im
    
    nrm = m_max;
    
  jump_1:
    
    ygm[nrm] = bgm[nrm] / rgm[_IDX2D(nrm, nrm, Nmax)];
    
    for (int im=nrm-1; im>=1; im--) {
      al = bgm[im];
      
      for (int jm=im+1; jm<=nrm; jm++) {
        al -= rgm[_IDX2D(im, jm, Nmax)] * ygm[jm];
      }
      
      ygm[im] = al / rgm[_IDX2D(im, im, Nmax)];
    }
    
    for (int im=1; im<=nrm; im++) {
      al = ygm[im];
      
      ml_add_4_(d_p, size, &guide, &m_max, &al, d_zm, &im, &flop);
    }
    
  } // loop; i_iter
  
jump_2:
  
  res = SOR_2_SMA(IC);
  
  
jump_3:
  
  
jump_4:
  
  if ( (isfin == 0) && (res < (res_rhs * eps_2 * eps_2)) )
  {
    res = 100.0 * res_rhs * eps_2 * eps_2;
  }
  
  if ( rgm ) delete [] rgm;
  if ( hgm ) delete [] hgm;
  if ( cgm ) delete [] cgm;
  if ( sgm ) delete [] sgm;
  if ( bgm ) delete [] bgm;
  if ( ygm ) delete [] ygm;
  
  TIMING_stop(tm_gmres_sor_sct, 0.0);
  // <<< Poisson Source section
  
  return res;
}


// #################################################################
// GMres SOR
double FFV::Gmres_SOR(ItrCtl* IC, double res_rhs)
{
  const double fct2  = 6.0;
  const double eps_1 = 1.0e-30;
  const double eps_2 = IC->get_eps();
  
  // 残差収束チェック
  if ( res_rhs < eps_1 ) return res_rhs;
  
  
  // >>> Gmres section
  TIMING_start(tm_gmres_sor_sct);
  
  double t_eps, beta, beta_1, res, eps_abs, res_abs, al;
  double r4;
  double flop=0.0;
  
  const int Iteration_Max = IC->get_ItrMax();
  int m_max = FREQ_OF_RESTART;
  int s_length = (size[0]+2*guide) * (size[1]+2*guide) * (size[2]+2*guide);
  
  int isfin  = 0;
  
  const int oki  = 2;
  const int step = 3;
  int nrm  = 0;

  if ( C.Mode.Precision == sizeof(float) )
  {
    t_eps = 0.995 * eps_2 * eps_2;
  }
  else
  {
    t_eps = 0.999 * eps_2 * eps_2;
  }

  const int Nmax = m_max+1; // 配列確保用，ゼロはダミー
  
  double *rgm = new double[Nmax * Nmax];     // (Nmax, Nmax)
  double *hgm = new double[(Nmax+1) * Nmax]; // (Nmax+1, Nmax)
  double *cgm = new double[Nmax];
  double *sgm = new double[Nmax];
  double *bgm = new double[Nmax];
  double *ygm = new double[Nmax];
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  int column = 1;
  
  t_eps = res_rhs * t_eps;

  res = SOR_2_SMA(IC);

  r4 = sqrt(res/res_rhs);

  
  if (res < t_eps)
  {
    //Hostonly_ printf("Final : %e\n", r4);
    isfin = 1;
    goto jump_4;
  }
  
  
  if ( res > (fct2*t_eps) ) t_eps = res/fct2;
  
  eps_abs = 0.15 * t_eps;
  
  
  for (int i_iter=1; i_iter<=Iteration_Max; i_iter++) {
    
    // テスト的にはずす >> 収束性良い？
    //res = SOR_2_SMA(IC);
    //if (res < t_eps) goto jump_3;

    TIMING_start(tm_gmres_mvprod);
    flop = 0.0;
    mv_prod_(d_yt, size, &guide, d_p, d_bcp, &flop);
    TIMING_stop(tm_gmres_mvprod, flop);
    
    
    TIMING_start(tm_gmres_res_sample);
    flop = 0.0;
    res_smpl_(d_rest, size, &guide, &res_abs, d_ws, d_yt, d_bcp, &flop);
    TIMING_stop(tm_gmres_res_sample, flop);
    
    
    TIMING_start(tm_gmres_comm);
    if ( numProc > 1 )
    {
      double tmp = res_abs;
      if ( paraMngr->Allreduce(&tmp, &res_abs, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
    }
    TIMING_stop(tm_gmres_comm, 2.0*numProc*sizeof(double)); // 双方向 x ノード数
    
    
    beta = sqrt(res_abs);
    
    if (beta < eps_1) goto jump_2;
    
    
    beta_1 = 1.0/beta;
    
    multiply_(d_vm, size, &guide, &m_max, &column, &beta_1, d_rest, &flop);
    
    
    bgm[1] = beta;
    
    for (int i=2; i<=m_max; i++) {
      bgm[i] = 0.0;
    }
    
    
    
    for (int im=1; im<=m_max; im++) {
      
      copy_1_(d_rest, size, &guide, &m_max, d_vm, &im);
      
      
#pragma omp parallel for firstprivate(s_length)
      for (int i=0; i<s_length; i++) {
        d_xt[i] = 0.0;
      }
      
      
      // Decide the number of iteration for inner-iteration
      int n_inner = 3; //((im / oki) + 1) * step;
      
      
      // Inner-iteration
      for (int i_inner=1; i_inner<=n_inner; i_inner++) {
        
        res = SOR_2_SMA(IC);
      }
      
      copy_2_(d_zm, size, &guide, &m_max, d_xt, &im);
      
      TIMING_start(tm_gmres_mvprod);
      flop = 0.0;
      mv_prod_(d_xt, size, &guide, d_p, d_bcp, &flop);
      TIMING_stop(tm_gmres_mvprod, flop);
      

      for (int km=1; km<=im; km++) {
        al = 0.0;
        ml_add_1_(&al, size, &guide, &m_max, d_vm, d_yt, &km, &flop);
        
        
        TIMING_start(tm_gmres_comm);
        if ( numProc > 1 )
        {
          double tmp = al;
          if ( paraMngr->Allreduce(&tmp, &al, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
        }
        TIMING_stop(tm_gmres_comm, 2.0*numProc*sizeof(double)); // 双方向 x ノード数
        
        
        hgm[_IDX2D(km, im, Nmax)] = al;
        
        double r_al = -al;
        ml_add_3_(d_yt, size, &guide, &m_max, &r_al, d_vm, &km, &flop);
      }

      
      al =0.0;
      ml_add_2_(&al, size, &guide, d_yt, &flop);
      
      
      TIMING_start(tm_gmres_comm);
      if ( numProc > 1 )
      {
        double tmp = al;
        if ( paraMngr->Allreduce(&tmp, &al, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
      }
      TIMING_stop(tm_gmres_comm, 2.0*numProc*sizeof(double)); // 双方向 x ノード数
      
      
      TIMING_start(tm_gmres_others);
      flop = 0.0;
      
      hgm[_IDX2D(im+1, im, Nmax)] = sqrt(al);
      flop += 20.0;
      
      if ( hgm[_IDX2D(im+1, im, Nmax)] < eps_1 )
      {
        nrm = im-1;
        TIMING_stop(tm_gmres_others, flop);
        goto jump_1;
      }

      if (im < m_max)
      {
        int idx = im + 1;
        al = 1.0 / hgm[_IDX2D(idx, im, Nmax)];
        flop += 13.0;
        multiply_(d_vm, size, &guide, &m_max, &idx, &al, d_yt, &flop);
      }

      rgm[_IDX2D(1, im, Nmax)] = hgm[_IDX2D(1, im, Nmax)];
      
      for (int km=2; km<=im; km++) {
        rgm[_IDX2D(km  , im, Nmax)] = cgm[km-1]*hgm[_IDX2D(km, im, Nmax)] - sgm[km-1]*rgm[_IDX2D(km-1, im, Nmax)];
        rgm[_IDX2D(km-1, im, Nmax)] = sgm[km-1]*hgm[_IDX2D(km, im, Nmax)] + cgm[km-1]*rgm[_IDX2D(km-1, im, Nmax)];
      }
      flop += (double)(im-2+1)*6.0;
      
      al = sqrt(rgm[_IDX2D(im, im, Nmax)]*rgm[_IDX2D(im, im, Nmax)] + hgm[_IDX2D(im+1, im, Nmax)]*hgm[_IDX2D(im+1, im, Nmax)]);
      flop += 23;
      
      if ( al < eps_1 )
      {
        nrm = im - 1;
        TIMING_stop(tm_gmres_others, flop);
        goto jump_1;
      }

      cgm[im]    = rgm[_IDX2D(im,  im, Nmax)] / al;
      sgm[im]    = hgm[_IDX2D(im+1,im, Nmax)] / al;
      
      rgm[_IDX2D(im, im, Nmax)] = cgm[im] * rgm[_IDX2D(im, im, Nmax)] + sgm[im] * hgm[_IDX2D(im+1, im, Nmax)];
      
      bgm[im+1]  = - sgm[im]*bgm[im];
      bgm[im]    =   cgm[im]*bgm[im];
      
      al = bgm[im+1] * bgm[im+1];

      flop += 2.0*13.0 + 4.0;

      if (al < eps_abs)
      {
        nrm = im;
        TIMING_stop(tm_gmres_others, flop);
        goto jump_1;
      }
      
      TIMING_stop(tm_gmres_others, flop);

    } // loop; im
    
    nrm = m_max;
    
  jump_1:
    
    ygm[nrm] = bgm[nrm] / rgm[_IDX2D(nrm, nrm, Nmax)];
    
    for (int im=nrm-1; im>=1; im--) {
      al = bgm[im];
      
      for (int jm=im+1; jm<=nrm; jm++) {
        al -= rgm[_IDX2D(im, jm, Nmax)] * ygm[jm];
      }
      
      ygm[im] = al / rgm[_IDX2D(im, im, Nmax)];
    }
    
    for (int im=1; im<=nrm; im++) {
      al = ygm[im];
      
      ml_add_4_(d_p, size, &guide, &m_max, &al, d_zm, &im, &flop);
    }
    
  } // loop; i_iter
  
jump_2:
  
  res = SOR_2_SMA(IC);
  
  
jump_3:
  
  
jump_4:
  
  if ( (isfin == 0) && (res < (res_rhs * eps_2 * eps_2)) )
  {
    res = 100.0 * res_rhs * eps_2 * eps_2;
  }
  
  if ( rgm ) delete [] rgm;
  if ( hgm ) delete [] hgm;
  if ( cgm ) delete [] cgm;
  if ( sgm ) delete [] sgm;
  if ( bgm ) delete [] bgm;
  if ( ygm ) delete [] ygm;
  
  TIMING_stop(tm_gmres_sor_sct, 0.0);
  // <<< Poisson Source section
  
  return res;
}
