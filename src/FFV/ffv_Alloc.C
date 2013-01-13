// #################################################################
//
// FFV : Frontflow / violet Cartesian
//
// Copyright (c) 2012-2013  All right reserved.
//
// Institute of Industrial Science, University of Tokyo, Japan. 
//
// #################################################################


/** 
 * @file   ffv_Alloc.C
 * @brief  FFV Class
 * @author kero
 */

#include "ffv.h"

// Adams-Bashforth法に用いる配列のアロケーション
void FFV::allocArray_AB2 (double &total)
{
  double mc = (double)(size[0] * size[1] * size[2]);
  
  // d_abf
  if ( !(d_abf = Alloc::Real_V3D(size, guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE) * 3.0;
}



// 平均値処理に用いる配列のアロケーション
void FFV::allocArray_Average (double &total)
{
  double mc = (double)(size[0] * size[1] * size[2]);
  
  // d_ap
  if ( !(d_ap = Alloc::Real_S3D(size, guide)) ) Exit(0);
  total += mc * (double)sizeof(REAL_TYPE);
  
  // d_av
  if ( !(d_av = Alloc::Real_V3D(size, guide)) ) Exit(0);
  total += mc * (double)sizeof(REAL_TYPE) * 3.0;
  
  if ( C.isHeatProblem() ) 
  {
    // d_at
    if ( !(d_at = Alloc::Real_S3D(size, guide)) ) Exit(0);
    total += mc * (double)sizeof(REAL_TYPE);
  }
}



// 粗格子読み込みに用いる配列のアロケーション
void FFV::allocArray_CoarseMesh(const int* r_size, double &prep)
{
  double mc = (double)(r_size[0] * r_size[1] * r_size[2]);
  
  // d_r_p
  if ( !(d_r_p = Alloc::Real_S3D(r_size, guide)) ) Exit(0);
  prep += mc * (double)sizeof(REAL_TYPE);
  
  // d_r_v
  if ( !(d_r_v = Alloc::Real_V3D(r_size, guide)) ) Exit(0);
  prep += mc * (double)sizeof(REAL_TYPE) * 3.0;
  
  
  if ( C.isHeatProblem() ) 
  {
    // d_r_t
    if ( !(d_r_t = Alloc::Real_S3D(r_size, guide)) ) Exit(0);
    prep += mc * (double)sizeof(REAL_TYPE);
    
  }
}



// コンポーネント体積率の配列のアロケーション
void FFV::allocArray_CompoVF(double &prep, double &total)
{
  double mc = (double)(size[0] * size[1] * size[2]);
  
  // d_cvf
  if ( !(d_cvf = paraMngr->AllocFloatS3D(guide)) ) Exit(0);
  prep += mc * (double)sizeof(float);
  total+= mc * (double)sizeof(float);
}



// カット情報の配列
void FFV::allocArray_Cut(double &total)
{
  double mc = (double)(size[0] * size[1] * size[2]);
  
  // d_cut
  if ( !(d_cut = paraMngr->AllocFloatS4DEx(6, guide)) ) Exit(0);
  total+= mc * (double)sizeof(float) * 6.0;
  
  // d_bid
  if ( !(d_bid = paraMngr->AllocIntS3D(guide)) ) Exit(0);
  total += mc * (double)sizeof(int);
}


// コンポーネントのワーク用配列のアロケート
void FFV::allocArray_Forcing(double& m_prep, double& m_total, FILE* fp)
{
  
  // 管理用のポインタ配列の確保
  component_array = new REAL_TYPE* [C.NoBC];
  
  for (int i=0; i<C.NoBC; i++) 
  {
    component_array[i] = NULL;
  }
  
  
  // リサイズ後のインデクスサイズの登録と配列領域の確保
  int c_sz[3];
  int gd=2;    // ガイドセルは，両側それぞれ2セル
  size_t m_cmp_size=0;
  
  for (int n=1; n<=C.NoBC; n++) 
  {
    
    if ( cmp[n].isFORCING() ) 
    {
      // インデクスサイズをリサイズしたst[], ed[]から計算
      cmp[n].set_cmp_sz();
      cmp[n].get_cmp_sz(c_sz);
      
      // ワーク用の配列を確保
      size_t array_size = (c_sz[0]+2*gd) * (c_sz[1]+2*gd) * (c_sz[2]+2*gd) * 3;
      component_array[n-1] = new REAL_TYPE[array_size];
      m_cmp_size += array_size;
    }
  }
  
  // 使用メモリ量　
  double cmp_mem, G_cmp_mem;
  G_cmp_mem = cmp_mem = (double)m_cmp_size * sizeof(double);
  m_prep += cmp_mem;
  m_total+= cmp_mem;
  
  if ( numProc > 1 ) 
  {
    paraMngr->Allreduce(&cmp_mem, &G_cmp_mem, 1, MPI_SUM);
  }
  
  if ( C.isForcing() )
  {
    Hostonly_  
    {
      FBUtility::MemoryRequirement("Component", G_cmp_mem, cmp_mem, stdout);
      FBUtility::MemoryRequirement("Component", G_cmp_mem, cmp_mem, fp);
    }   
  }
}


// Krylov-subspace法に用いる配列のアロケーション
void FFV::allocArray_Krylov(double &total)
{
  double mc = (double)(size[0] * size[1] * size[2]);
  
  if ( !(d_wg = Alloc::Real_S3D(size, guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE);
  
  if ( !(d_res = Alloc::Real_S3D(size, guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE);
  
  if ( !(d_vm = Alloc::Real_S4D(size, guide, FREQ_OF_RESTART+1)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE);
  
  if ( !(d_zm = Alloc::Real_S4D(size, guide, FREQ_OF_RESTART+1)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE);
  
}


// PCG法に用いる配列のアロケーション
void FFV::allocArray_PCG(double &total)
{
  double mc = (double)(size[0] * size[1] * size[2]);
  
  if ( !(d_pcg_r = Alloc::Real_S3D(size, guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE);
  
  if ( !(d_pcg_p = Alloc::Real_S3D(size, guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE);
  
  if ( !(d_pcg_q = Alloc::Real_S3D(size, guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE);
  
  if ( !(d_pcg_z = Alloc::Real_S3D(size, guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE);
}


// PBiCGSTAB法に用いる配列のアロケーション
void FFV::allocArray_PBiCGSTAB(double &total)
{
  double mc = (double)(size[0] * size[1] * size[2]);
  
  if ( !(d_pcg_r = Alloc::Real_S3D(size, guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE);
  
  if ( !(d_pcg_p = Alloc::Real_S3D(size, guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE);
  
  if ( !(d_pcg_r0 = Alloc::Real_S3D(size, guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE);
  
  if ( !(d_pcg_p_ = Alloc::Real_S3D(size, guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE);
  
  if ( !(d_pcg_q_ = Alloc::Real_S3D(size, guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE);
  
	if ( !(d_pcg_s = Alloc::Real_S3D(size, guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE);
  
	if ( !(d_pcg_s_ = Alloc::Real_S3D(size, guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE);
  
	if ( !(d_pcg_t_ = Alloc::Real_S3D(size, guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE);
}


// 熱の主計算部分に用いる配列のアロケーション
void FFV::allocArray_Heat(double &total)
{
  double mc = (double)(size[0] * size[1] * size[2]);
  
  // d_t
  if ( !(d_t = Alloc::Real_S3D(size, guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE);
  
  // d_t0
  if ( !(d_t0 = Alloc::Real_S3D(size, guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE);
  
  // d_qbc
  if ( !(d_qbc = Alloc::Real_V3D(size, guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE) * 3.0;
}



// 体積率の配列のアロケーション
void FFV::allocArray_Interface(double &total)
{
  double mc = (double)(size[0] * size[1] * size[2]);
  
  // d_vof
  if ( !(d_vof = Alloc::Real_S3D(size, guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE);
  
}



// LES計算に用いる配列のアロケーション
void FFV::allocArray_LES(double &total)
{
  double mc = (double)(size[0] * size[1] * size[2]);
  
  // d_vt
  if ( !(d_vt = Alloc::Real_S3D(size, guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE);
}



// 主計算部分に用いる配列のアロケーション
void FFV::allocArray_Main(double &total)
{
  double mc = (double)(size[0] * size[1] * size[2]);
  
  // d_v
  if ( !(d_v = Alloc::Real_V3D(size, guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE) * 3.0;
  
  // d_vf
  if ( !(d_vf = Alloc::Real_V3D(size, guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE) * 3.0;
  
  // d_vc
  if ( !(d_vc = Alloc::Real_V3D(size, guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE) * 3.0;
  
  // d_v0
  if ( !(d_v0 = Alloc::Real_V3D(size, guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE) * 3.0;
  
  // d_wv
  if ( !(d_wv = Alloc::Real_V3D(size, guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE) * 3.0;
  
  // d_wo
  if ( !(d_wo = Alloc::Real_V3D(size, guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE) * 3.0;
  
  // d_p
  if ( !(d_p = Alloc::Real_S3D(size, guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE);
  
  // d_p0
  if ( !(d_p0 = Alloc::Real_S3D(size, guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE);
  
  // d_sq
  if ( !(d_sq = Alloc::Real_S3D(size, guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE);
  
  // d_dv
  if ( !(d_dv = Alloc::Real_S3D(size, guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE);
  
  // d_b
  if ( !(d_b = Alloc::Real_S3D(size, guide)) ) Exit(0);
  total+= mc * (double)sizeof(REAL_TYPE);
}



// 前処理に用いる配列のアロケーション
void FFV::allocArray_Prep(double &prep, double &total)
{
  double mc = (double)(size[0] * size[1] * size[2]);
  
  // d_ws
  if ( !(d_ws = Alloc::Real_S3D(size, guide)) ) Exit(0);
  prep += mc * (double)sizeof(REAL_TYPE);
  total+= mc * (double)sizeof(REAL_TYPE);
  
  // d_mid
  if ( !(d_mid = paraMngr->AllocIntS3D(guide)) ) Exit(0);
  prep += mc * (double)sizeof(int);
  
  // d_bcd
  if ( !(d_bcd = paraMngr->AllocIntS3D(guide)) ) Exit(0);
  prep += mc * (double)sizeof(int);
  total+= mc * (double)sizeof(int);
  
  // d_bcp
  if ( !(d_bcp = paraMngr->AllocIntS3D(guide)) ) Exit(0);
  prep += mc * (double)sizeof(int);
  total+= mc * (double)sizeof(int);
  
  // d_bcv
  if ( !(d_bcv = paraMngr->AllocIntS3D(guide)) ) Exit(0);
  prep += mc * (double)sizeof(int);
  total+= mc * (double)sizeof(int);
  
  if ( C.isHeatProblem() ) 
  {
    // d_bh1
    if ( !(d_bh1 = paraMngr->AllocIntS3D(guide)) ) Exit(0);
    prep += mc * (double)sizeof(int);
    total+= mc * (double)sizeof(int);
    
    // d_bh2
    if ( !(d_bh2 = paraMngr->AllocIntS3D(guide)) ) Exit(0);
    prep += mc * (double)sizeof(int);
    total+= mc * (double)sizeof(int);
  }
}


// SOR2SMAのバッファ確保
void FFV::allocate_SOR2SMA_buffer(double &total)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  
  cf_sz[0] = (jx+1) * (kx+1) / 2; // バッファサイズ
  cf_sz[1] = (kx+1) * (ix+1) / 2; // +1はマージン
  cf_sz[2] = (ix+1) * (jx+1) / 2;
  
  size_t n1 = cf_sz[0]*4;
  size_t n2 = cf_sz[1]*4;
  size_t n3 = cf_sz[2]*4;
  
  if( (cf_x = new REAL_TYPE[n1]) == NULL ) Exit(0);
  if( (cf_y = new REAL_TYPE[n2]) == NULL ) Exit(0);
  if( (cf_z = new REAL_TYPE[n3]) == NULL ) Exit(0);
  
  memset(cf_x, 0, sizeof(REAL_TYPE)*n1);
  memset(cf_y, 0, sizeof(REAL_TYPE)*n2);
  memset(cf_z, 0, sizeof(REAL_TYPE)*n3);
  
  total += (double)( (n1+n2+n3)*sizeof(REAL_TYPE) );
}


// 主計算に用いる配列の確保
void FFV::allocate_Main(double &total)
{
  TIMING_start(tm_init_alloc);
  allocArray_Main(total);
  
  //allocArray_Collocate (total);
  
  if ( C.LES.Calc == ON ) 
  {
    allocArray_LES (total);
  }
  
  if ( C.isHeatProblem() ) 
  {
    allocArray_Heat(total);
  }
  
  if ( (C.AlgorithmF == Control::Flow_FS_AB2) || (C.AlgorithmF == Control::Flow_FS_AB_CN) ) 
  {
    allocArray_AB2(total);
  }
  
  if ( C.BasicEqs == INCMP_2PHASE ) 
  {
    allocArray_Interface(total);
  }
  
  // 時間平均用の配列をアロケート
  if ( C.Mode.Average == ON ) 
  {
    allocArray_Average(total);
  }
  
  TIMING_stop(tm_init_alloc);
}

