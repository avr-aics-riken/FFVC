//##################################################################################
//
// FFV-C : Frontflow / violet Cartesian
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2015 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2015 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################


/** 
 * @file   ffv_Alloc.C
 * @brief  FFV Class
 * @author aics
 */

#include "ffv_Alloc.h"
#include "mydebug.h"
#include "FBUtility.h"
#include "Alloc.h"
#include <string.h>

//#include <stdio.h>
//#include <stdlib.h>

//#include <string>
//#include <iostream>
//#include <fstream>

/*
#ifndef _WIN32
#include <unistd.h>
#include <strings.h>
#else
#include "sph_win32_util.h"
#endif
#include <sys/types.h>

#if defined(IA32_LINUX) || defined(IA64_LINUX) || defined(SGI_ALTIX)
#include <sys/stat.h>
#endif

#ifdef MacOSX
#include <sys/uio.h>
#endif
 */

// #################################################################
/* @brief Adams-Bashforth法に用いる配列のアロケーション
 * @param [in,out] total ソルバーに使用するメモリ量
 */
void FALLOC::allocArray_AB2 (double &total)
{
  // d_abf
  if ( !(d_abf = Alloc::Real_V3D(size, guide)) ) Exit(0);
  total+= array_size * (double)sizeof(REAL_TYPE) * 3.0;
}


// #################################################################
/* @brief 平均処理に用いる配列のアロケーション
 * @param [in,out] total  ソルバーに使用するメモリ量
 * @param [in]     C      Control class
 */
void FALLOC::allocArray_Average(double &total, Control* C)
{
  // d_ap
  if ( !(d_ap = Alloc::Real_S3D(size, guide)) ) Exit(0);
  total += array_size * (double)sizeof(REAL_TYPE);
  
  // d_av
  if ( !(d_av = Alloc::Real_V3D(size, guide)) ) Exit(0);
  total += array_size * (double)sizeof(REAL_TYPE) * 3.0;
  
  if ( C->isHeatProblem() )
  {
    // d_ae
    if ( !(d_ae = Alloc::Real_S3D(size, guide)) ) Exit(0);
    total += array_size * (double)sizeof(REAL_TYPE);
  }
}


// #################################################################
/* @brief 粗格子読み込みに用いる配列のアロケーション
 * @param [in]     r_size  粗格子の領域サイズ
 * @param [in,out] prep    前処理に使用するメモリ量
 * @param [in]     isHeat 熱問題のときtrue
 */
void FALLOC::allocArray_CoarseMesh(const int* r_size, double &prep, const bool isHeat)
{
  // d_r_p
  if ( !(d_r_p = Alloc::Real_S3D(r_size, guide)) ) Exit(0);
  prep += array_size * (double)sizeof(REAL_TYPE);
  
  // d_r_v
  if ( !(d_r_v = Alloc::Real_V3D(r_size, guide)) ) Exit(0);
  prep += array_size * (double)sizeof(REAL_TYPE) * 3.0;
  
  
  if ( isHeat )
  {
    // d_r_t
    if ( !(d_r_t = Alloc::Real_S3D(r_size, guide)) ) Exit(0);
    prep += array_size * (double)sizeof(REAL_TYPE);
    
  }
}


// #################################################################
/**
 * @brief コンポーネント体積率の配列のアロケーション
 * @param [in,out] prep  前処理に使用するメモリ量
 * @param [in,out] total ソルバーに使用するメモリ量
 */
void FALLOC::allocArray_CompoVF(double &prep, double &total)
{
  if ( !(d_cvf = Alloc::Real_S3D(size, guide)) ) Exit(0);
  prep += array_size * (double)sizeof(REAL_TYPE);
  total+= array_size * (double)sizeof(REAL_TYPE);
}


// #################################################################
/**
 * @brief カット情報の配列
 * @param [in,out] total ソルバーに使用するメモリ量
 */
void FALLOC::allocArray_Cut(double &prep, double &total)
{
  if ( !(d_cut = Alloc::LLong_S3D(size, guide)) ) Exit(0);
  prep  += array_size * (double)sizeof(long long);
  total += array_size * (double)sizeof(long long);
  
  
  if ( !(d_bid = Alloc::Int_S3D(size, guide)) ) Exit(0);
  prep  += array_size * (double)sizeof(int);
  total += array_size * (double)sizeof(int);
}


// #################################################################
/**
 * @brief コンポーネントのワーク用配列のアロケート
 * @param [in,out] m_prep  前処理用のメモリサイズ
 * @param [in,out] m_total 本計算用のメモリリサイズ
 * @param [in]     fp      ファイルポインタ
 * @param [in]     C       Control class
 * @param [in]     cmp     CompoList class
 */
void FALLOC::allocArray_Forcing(double& m_prep, double& m_total, FILE* fp, Control* C, CompoList* cmp)
{
  
  // 管理用のポインタ配列の確保
  component_array = new REAL_TYPE* [C->NoCompo];
  
  for (int i=0; i<C->NoCompo; i++)
  {
    component_array[i] = NULL;
  }
  
  
  // リサイズ後のインデクスサイズの登録と配列領域の確保
  int c_sz[3];
  int gd=2;    // ガイドセルは，両側それぞれ2セル
  size_t m_cmp_size=0;
  
  for (int n=1; n<=C->NoCompo; n++)
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
  
  if ( C->EnsCompo.forcing )
  {
    Hostonly_  
    {
      FBUtility::MemoryRequirement("Component", G_cmp_mem, cmp_mem, stdout);
      FBUtility::MemoryRequirement("Component", G_cmp_mem, cmp_mem, fp);
    }   
  }
}


// #################################################################
/**
 * @brief Krylov-subspace Iteration
 * @param [in,out] total ソルバーに使用するメモリ量
 */
void FALLOC::allocArray_Krylov(double &total)
{
  if ( !(d_wg = Alloc::Real_S3D(size, guide)) ) Exit(0);
  total+= array_size * (double)sizeof(REAL_TYPE);
  
  
  if ( !(d_res = Alloc::Real_S3D(size, guide)) ) Exit(0);
  total+= array_size * (double)sizeof(REAL_TYPE);
  
  
  if ( !(d_vm = Alloc::Real_S4D(size, guide, FREQ_OF_RESTART+1)) ) Exit(0);
  total+= array_size * (double)sizeof(REAL_TYPE);
  
  
  if ( !(d_zm = Alloc::Real_S4D(size, guide, FREQ_OF_RESTART+1)) ) Exit(0);
  total+= array_size * (double)sizeof(REAL_TYPE);
  
}



// #################################################################
/**
 * @brief 体積率の配列のアロケーション
 * @param [in,out] total ソルバーに使用するメモリ量
 */
void FALLOC::allocArray_Interface(double &total)
{
  if ( !(d_vof = Alloc::Real_S3D(size, guide)) ) Exit(0);
  total+= array_size * (double)sizeof(REAL_TYPE);
  
}


// #################################################################
/**
 * @brief LES計算に用いる配列のアロケーション
 * @param [in,out] total ソルバーに使用するメモリ量
 */
void FALLOC::allocArray_LES(double &total)
{
  if ( !(d_vt = Alloc::Real_S3D(size, guide)) ) Exit(0);
  total+= array_size * (double)sizeof(REAL_TYPE);
}


// #################################################################
/**
 * @brief 統計に用いる配列のアロケーション
 * @param [in,out] total ソルバーに使用するメモリ量
 */
void FALLOC::allocArray_Statistic(double &total, Control* C)
{
  // Velocity
  if ( C->Mode.StatVelocity == ON )
  {
    if ( !(d_rms_v = Alloc::Real_V3D(size, guide)) ) Exit(0);
    total+= array_size * (double)sizeof(REAL_TYPE) * 3.0;
    
    if ( !(d_rms_mean_v = Alloc::Real_V3D(size, guide)) ) Exit(0);
    total+= array_size * (double)sizeof(REAL_TYPE) * 3.0;
  }
  
  
  // Pressure
  if ( C->Mode.StatPressure == ON )
  {
    if ( !(d_rms_p = Alloc::Real_S3D(size, guide)) ) Exit(0);
    total+= array_size * (double)sizeof(REAL_TYPE);
    
    if ( !(d_rms_mean_p = Alloc::Real_S3D(size, guide)) ) Exit(0);
    total+= array_size * (double)sizeof(REAL_TYPE);
  }
  
  
  // Temperature
  if ( C->isHeatProblem() && C->Mode.StatTemperature==ON )
  {
    if ( !(d_rms_t = Alloc::Real_S3D(size, guide)) ) Exit(0);
    total+= array_size * (double)sizeof(REAL_TYPE);
    
    if ( !(d_rms_mean_t = Alloc::Real_S3D(size, guide)) ) Exit(0);
    total+= array_size * (double)sizeof(REAL_TYPE);
  }

  
  // Reynolds Stress
  if ( C->Mode.ReynoldsStress == ON )
  {
    // レイノルズ応力テンソル
    if ( !(d_R = Alloc::Real_T3D(size, guide)) ) Exit(0);
    total += array_size * (double)sizeof(REAL_TYPE) * 6.0;
    
    // レイノルズ応力テンソル (時間平均値)
    if ( !(d_aR = Alloc::Real_T3D(size, guide)) ) Exit(0);
    total += array_size * (double)sizeof(REAL_TYPE) * 6.0;
    
    // 生成項 (時間平均値)
    if ( !(d_aP = Alloc::Real_T3D(size, guide)) ) Exit(0);
    total += array_size * (double)sizeof(REAL_TYPE) * 6.0;
    
    // 散逸項 (時間平均値)
    if ( !(d_aE = Alloc::Real_T3D(size, guide)) ) Exit(0);
    total += array_size * (double)sizeof(REAL_TYPE) * 6.0;
    
    // 乱流拡散項 (時間平均値)
    if ( !(d_aT = Alloc::Real_T3D(size, guide)) ) Exit(0);
    total += array_size * (double)sizeof(REAL_TYPE) * 6.0;
    
    // 速度圧力勾配相関項 (時間平均値)
    if ( !(d_aPI = Alloc::Real_T3D(size, guide)) ) Exit(0);
    total += array_size * (double)sizeof(REAL_TYPE) * 6.0;
  }
  
  // Channel Mean
  if ( C->Mode.ChannelOutputMean == ON )
  {
    if ( !(d_av_mean = Alloc::Real_S3D(size, guide)) ) Exit(0);
    total += array_size * (double)sizeof(REAL_TYPE) * 3.0;
    
    if ( !(d_arms_mean = Alloc::Real_S3D(size, guide)) ) Exit(0);
    total += array_size * (double)sizeof(REAL_TYPE) * 3.0;
    
    if ( !(d_aR_mean = Alloc::Real_T3D(size, guide)) ) Exit(0);
    total += array_size * (double)sizeof(REAL_TYPE) * 6.0;
    
    if ( !(d_aP_mean = Alloc::Real_T3D(size, guide)) ) Exit(0);
    total += array_size * (double)sizeof(REAL_TYPE) * 6.0;
    
    if ( !(d_aE_mean = Alloc::Real_T3D(size, guide)) ) Exit(0);
    total += array_size * (double)sizeof(REAL_TYPE) * 6.0;
    
    if ( !(d_aT_mean = Alloc::Real_T3D(size, guide)) ) Exit(0);
    total += array_size * (double)sizeof(REAL_TYPE) * 6.0;
    
    if ( !(d_aPI_mean = Alloc::Real_T3D(size, guide)) ) Exit(0);
    total += array_size * (double)sizeof(REAL_TYPE) * 6.0;
  }
  
}


// #################################################################
/**
 * @brief 主計算部分に用いる配列のアロケーション
 * @param [in,out] total ソルバーに使用するメモリ量
 */
void FALLOC::allocArray_Main(double &total, Control* C)
{
  // [*] ステップ間でデータを保持
  if ( !(d_v = Alloc::Real_V3D(size, guide)) ) Exit(0);
  total+= array_size * (double)sizeof(REAL_TYPE) * 3.0;
  

  if ( !(d_vf = Alloc::Real_V3D(size, guide)) ) Exit(0);
  total+= array_size * (double)sizeof(REAL_TYPE) * 3.0;
  

  if ( !(d_p = Alloc::Real_S3D(size, guide)) ) Exit(0);
  total+= array_size * (double)sizeof(REAL_TYPE);
  
  
  if ( !(d_dv = Alloc::Real_S3D(size, guide)) ) Exit(0);
  total+= array_size * (double)sizeof(REAL_TYPE);
  
  
  if ( !(d_wv = Alloc::Real_V3D(size, guide)) ) Exit(0);
  total+= array_size * (double)sizeof(REAL_TYPE) * 3.0;
  
  
  if ( C->isHeatProblem() )
  {
    if ( !(d_ie = Alloc::Real_S3D(size, guide)) ) Exit(0);
    total+= array_size * (double)sizeof(REAL_TYPE);
  }
  
  
  // 大きなバッファを用意
  int dnum;
  
  if ( C->isHeatProblem() )
  {
    dnum = IO_BLOCK_SIZE_HEAT;
  }
  else
  {
    dnum = IO_BLOCK_SIZE_FLOW;
  }
  
  if ( !(d_io_buffer = Alloc::Real_S4D(size, guide, dnum)) ) Exit(0);
  total+= array_size * (double)sizeof(REAL_TYPE) * (double)dnum;
  
  
  // 割り当て
  size_t dims[3];
  dims[0] = (size_t)(size[0] + 2*guide);
  dims[1] = (size_t)(size[1] + 2*guide);
  dims[2] = (size_t)(size[2] + 2*guide);
  
  size_t nx = dims[0] * dims[1] * dims[2];
  
  
  d_vc = &d_io_buffer[0];
  d_v0 = &d_io_buffer[nx*3];
  d_p0 = &d_io_buffer[nx*6];
  d_sq = &d_io_buffer[nx*7];
  d_b  = &d_io_buffer[nx*8];
  
  if ( C->isHeatProblem() )
  {
    d_ie0= &d_io_buffer[nx*9];
    d_qbc= &d_io_buffer[nx*10];
  }
  
  
  
  // 渦度の出力指定がある，あるいは渦度関連のサンプリングがある場合にアロケート
  if ( C->varState[var_Vorticity] )
  {
    if ( !(d_vrt = Alloc::Real_V3D(size, guide)) ) Exit(0);
    total+= array_size * (double)sizeof(REAL_TYPE) * 3.0;
  }
}


// #################################################################
/**
 * @brief BiCGstab Iteration
 * @param [in,out] total ソルバーに使用するメモリ量
 */
void FALLOC::allocArray_BiCGstab(double &total)
{
  if ( !(d_pcg_r = Alloc::Real_S3D(size, guide)) ) Exit(0);
  total+= array_size * (double)sizeof(REAL_TYPE);
  
  
  if ( !(d_pcg_p = Alloc::Real_S3D(size, guide)) ) Exit(0);
  total+= array_size * (double)sizeof(REAL_TYPE);
  
  
  if ( !(d_pcg_r0 = Alloc::Real_S3D(size, guide)) ) Exit(0);
  total+= array_size * (double)sizeof(REAL_TYPE);
  
  
  if ( !(d_pcg_q = Alloc::Real_S3D(size, guide)) ) Exit(0);
  total+= array_size * (double)sizeof(REAL_TYPE);
  
  
  if ( !(d_pcg_s = Alloc::Real_S3D(size, guide)) ) Exit(0);
  total+= array_size * (double)sizeof(REAL_TYPE);
  
  
  if ( !(d_pcg_t = Alloc::Real_S3D(size, guide)) ) Exit(0);
  total+= array_size * (double)sizeof(REAL_TYPE);
}


// #################################################################
/**
 * @brief BiCGSTAB Iteration
 * @param [in,out] total ソルバーに使用するメモリ量
 */
void FALLOC::allocArray_BiCGSTABwithPreconditioning(double &total)
{
  if ( !(d_pcg_p_ = Alloc::Real_S3D(size, guide)) ) Exit(0);
  total+= array_size * (double)sizeof(REAL_TYPE);
  
	if ( !(d_pcg_s_ = Alloc::Real_S3D(size, guide)) ) Exit(0);
  total+= array_size * (double)sizeof(REAL_TYPE);
  
  if ( !(d_pcg_t_ = Alloc::Real_S3D(size, guide)) ) Exit(0);
  total+= array_size * (double)sizeof(REAL_TYPE);
}


// #################################################################
/**
 * @brief PCG Iteration
 * @param [in,out] total ソルバーに使用するメモリ量
 */
void FALLOC::allocArray_PCG(double &total)
{
  if ( !(d_pcg_r = Alloc::Real_S3D(size, guide)) ) Exit(0);
  total+= array_size * (double)sizeof(REAL_TYPE);
  
  
  if ( !(d_pcg_p = Alloc::Real_S3D(size, guide)) ) Exit(0);
  total+= array_size * (double)sizeof(REAL_TYPE);
  
  
  if ( !(d_pcg_q = Alloc::Real_S3D(size, guide)) ) Exit(0);
  total+= array_size * (double)sizeof(REAL_TYPE);
  
  
  if ( !(d_pcg_z = Alloc::Real_S3D(size, guide)) ) Exit(0);
  total+= array_size * (double)sizeof(REAL_TYPE);
}


// #################################################################
/**
 * @brief 前処理に用いる配列のアロケーション
 * @param [in,out] prep  前処理に使用するメモリ量
 * @param [in,out] total ソルバーに使用するメモリ量
 */
void FALLOC::allocArray_Prep(double &prep, double &total)
{
  if ( !(d_ws = Alloc::Real_S3D(size, guide)) ) Exit(0);
  prep += array_size * (double)sizeof(REAL_TYPE);
  total+= array_size * (double)sizeof(REAL_TYPE);
  

  if ( !(d_mid = Alloc::Int_S3D(size, guide)) ) Exit(0);
  prep += array_size * (double)sizeof(int);
  

  if ( !(d_bcd = Alloc::Int_S3D(size, guide)) ) Exit(0);
  prep += array_size * (double)sizeof(int);
  total+= array_size * (double)sizeof(int);
  

  if ( !(d_bcp = Alloc::Int_S3D(size, guide)) ) Exit(0);
  prep += array_size * (double)sizeof(int);
  total+= array_size * (double)sizeof(int);
  

  if ( !(d_cdf = Alloc::Int_S3D(size, guide)) ) Exit(0);
  prep += array_size * (double)sizeof(int);
  total+= array_size * (double)sizeof(int);
  
  if ( !(d_pvf = Alloc::Real_S3D(size, guide)) ) Exit(0);
  prep += array_size * (double)sizeof(REAL_TYPE);
  total+= array_size * (double)sizeof(REAL_TYPE);
}


// #################################################################
/**
 * @brief SOR2SMAのバッファ確保
 * @param [in,out] total ソルバーに使用するメモリ量
 */
void FALLOC::allocate_SOR2SMA_buffer(double &total)
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
