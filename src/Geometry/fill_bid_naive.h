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
 * @file   fill_bid_naive.h
<<<<<<< HEAD
 * @brief  bidによるフィルアルゴリズム、simpleなフラッドフィル
=======
 * @brief  bidによるフィルアルゴリズム、simpleなシードフィル
>>>>>>> origin/master
 * @author aics
 */


size_t m_p = _F_IDX_S3D(i  , j  , k  , ix, jx, kx, gd);
size_t m_e = _F_IDX_S3D(i+1, j,   k,   ix, jx, kx, gd);
size_t m_w = _F_IDX_S3D(i-1, j,   k,   ix, jx, kx, gd);
size_t m_n = _F_IDX_S3D(i,   j+1, k,   ix, jx, kx, gd);
size_t m_s = _F_IDX_S3D(i,   j-1, k,   ix, jx, kx, gd);
size_t m_t = _F_IDX_S3D(i,   j,   k+1, ix, jx, kx, gd);
size_t m_b = _F_IDX_S3D(i,   j,   k-1, ix, jx, kx, gd);

// bcd[]の0-4ビットにエンコードされたエントリの取り出し
int zp = DECODE_CMP( bcd[m_p] );
int zw = DECODE_CMP( bcd[m_w] );
int ze = DECODE_CMP( bcd[m_e] );
int zs = DECODE_CMP( bcd[m_s] );
int zn = DECODE_CMP( bcd[m_n] );
int zb = DECODE_CMP( bcd[m_b] );
int zt = DECODE_CMP( bcd[m_t] );


// 未ペイントの場合にテスト
if ( zp == 0 )
{
  int qq = bid[m_p];
  
  // 隣接セル方向に対する境界ID
  int qw = getBit5(qq, 0);
  int qe = getBit5(qq, 1);
  int qs = getBit5(qq, 2);
  int qn = getBit5(qq, 3);
  int qb = getBit5(qq, 4);
  int qt = getBit5(qq, 5);

  
  // 各方向のテスト
  // 対象セルは，未評価セル(zp == 0) 
  // 例えば、X-側からテストする場合は，X-側のセルがfill_modeのIDで既にペイントされている
  
  
  // フラッドフィル
  int tag = 0;
  int tg;
  
  // mode_x==0の時には，X方向の領域境界でフィルしない
  if ( (sdw < 0) && (i == 1) && !mode_x )
  {
    ; // skip
  }
  else
  {
    if ( mat[zw].getState()==fill_mode && qw==0 )
    {
      tg = zw;
      tag++;
    }
  }
  
  if ( (sde < 0) && (i == ix) && !mode_x )
  {
    ; // skip
  }
  else
  {
    if ( mat[ze].getState()==fill_mode && qe==0 )
    {
      tg = ze;
      tag++;
    }
  }
  
  // mode_y==0の時には，Y方向の領域境界でフィルしない
  if ( (sds < 0) && (j == 1) && !mode_y )
  {
    ; // skip
  }
  else
  {
    if ( mat[zs].getState()==fill_mode && qs==0 )
    {
      tg = zs;
      tag++;
    }
  }
  
  if ( (sdn < 0) && (j == jx) && !mode_y )
  {
    ;
  }
  else
  {
    if ( mat[zn].getState()==fill_mode && qn==0 )
    {
      tg = zn;
      tag++;
    }
  }
  
  // mode_z==0の時には，Z方向の領域境界でフィルしない
  if ( (sdb < 0) && (k == 1) && !mode_z )
  {
    ; // skip
  }
  else
  {
    if ( mat[zb].getState()==fill_mode && qb==0 )
    {
      tg = zb;
      tag++;
    }
  }
  
  if ( (sdt < 0) && (k == kx) && !mode_z )
  {
    ; // skip
  }
  else
  {
    if ( mat[zt].getState()==fill_mode && qt==0 )
    {
      tg = zt;
      tag++;
    }
  }
  

  if ( tag>0 )
  {
    if ( (DECODE_CMP(bcd[m_p]) == 0) )
    {
      setMediumID(bcd[m_p], tg);
      filled++;
    }
  }

  
} // target
