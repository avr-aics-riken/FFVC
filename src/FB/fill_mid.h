//##################################################################################
//
// Flow Base class
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2013 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2013 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

/**
 * @file   fill_mid.h
 * @brief  midによるフィルアルゴのコア
 * @author kero
 */

float cpos = 0.5;

size_t m_p = _F_IDX_S3D(i,   j,   k,   ix, jx, kx, gd);
size_t m_e = _F_IDX_S3D(i+1, j,   k,   ix, jx, kx, gd);
size_t m_w = _F_IDX_S3D(i-1, j,   k,   ix, jx, kx, gd);
size_t m_n = _F_IDX_S3D(i,   j+1, k,   ix, jx, kx, gd);
size_t m_s = _F_IDX_S3D(i,   j-1, k,   ix, jx, kx, gd);
size_t m_t = _F_IDX_S3D(i,   j,   k+1, ix, jx, kx, gd);
size_t m_b = _F_IDX_S3D(i,   j,   k-1, ix, jx, kx, gd);

int zp = mid[m_p];
int zw = mid[m_w];
int ze = mid[m_e];
int zs = mid[m_s];
int zn = mid[m_n];
int zb = mid[m_b];
int zt = mid[m_t];

// 未ペイントの場合にテスト
if ( zp == tg )
{
  int qq = bid[m_p];
  int sd;
  
  // 隣接セルの方向に対するカットIDの有無>> 0ならばカット無し、チェック半径は1
  int qw = getFaceBID(0, qq);
  int qe = getFaceBID(1, qq);
  int qs = getFaceBID(2, qq);
  int qn = getFaceBID(3, qq);
  int qb = getFaceBID(4, qq);
  int qt = getFaceBID(5, qq);
  

  int ff = 0;
  for (int n=1; n<=NoCompo; n++) {
    if ( list[n] == qw ) ff++; // CELL_MONITORであれば ff>0
  }
  if ( ff > 0 ) qw=0; // カットなし
    
  ff = 0;
  for (int n=1; n<=NoCompo; n++) {
    if ( list[n] == qe ) ff++;
  }
  if ( ff > 0 ) qe=0;
    
  ff = 0;
  for (int n=1; n<=NoCompo; n++) {
    if ( list[n] == qs ) ff++;
  }
  if ( ff > 0 ) qs=0;
    
  ff = 0;
  for (int n=1; n<=NoCompo; n++) {
    if ( list[n] == qn ) ff++;
  }
  if ( ff > 0 ) qn=0;
    
  ff = 0;
  for (int n=1; n<=NoCompo; n++) {
    if ( list[n] == qb ) ff++;
  }
  if ( ff > 0 ) qb=0;
    
  ff = 0;
  for (int n=1; n<=NoCompo; n++) {
    if ( list[n] == qt ) ff++;
  }
  if ( ff > 0 ) qt=0;
    
    
  
  if ( ((zs != 0) && (zs != tg)  &&  // X平面：ゼロでなく、かつ、targetでもない ==> 固体
        (zn != 0) && (zn != tg)) ||  // n-sセルの両方が固体の場合
       ((zb != 0) && (zb != tg)  &&
        (zt != 0) && (zt != tg)) )
  {
    
    if ( (zw == tg) && (qw == 0) ) // from X-
    {
      sd = find_mode_id(tgt_id, qw, qe, qs, qn, qb, qt);
      mid[m_p] = sd; // セルIDを固体に変更
      setFaceBID(bid[m_w], X_PLUS, sd); // カットIDを設定
      cut[_F_IDX_S4DEX(X_PLUS, i-1, j, k, 6, ix, jx, kx, gd)] = cpos; // カット位置をセット
      replaced++;
    }
    else if ( (ze == tg) && (qe == 0) ) // from X+
    {
      sd = find_mode_id(tgt_id, qw, qe, qs, qn, qb, qt);
      mid[m_p] = sd;
      setFaceBID(bid[m_e], X_MINUS, sd);
      cut[_F_IDX_S4DEX(X_MINUS, i+1, j, k, 6, ix, jx, kx, gd)] = cpos;
      replaced++;
    }
  }
  
  if ( ((zw != 0) && (zw != tg)  && // Y平面
        (ze != 0) && (ze != tg)) ||
       ((zb != 0) && (zb != tg)  &&
        (zt != 0) && (zt != tg)) )
  {
    
    if ( (zs == tg) && (qs == 0) ) // from Y-
    {
      sd = find_mode_id(tgt_id, qw, qe, qs, qn, qb, qt);
      mid[m_p] = sd;
      setFaceBID(bid[m_s], Y_PLUS, sd);
      cut[_F_IDX_S4DEX(Y_PLUS, i, j-1, k, 6, ix, jx, kx, gd)] = cpos;
      replaced++;
    }
    else if ( (zn == tg) && (qn == 0) ) // from Y+
    {
      sd = find_mode_id(tgt_id, qw, qe, qs, qn, qb, qt);
      mid[m_p] = sd;
      setFaceBID(bid[m_n], Y_MINUS, sd);
      cut[_F_IDX_S4DEX(Y_MINUS, i, j+1, k, 6, ix, jx, kx, gd)] = cpos;
      replaced++;
    }
  }
  
  if ( ((zw != 0) && (zw != tg)  && // Z平面
        (ze != 0) && (ze != tg)) ||
       ((zs != 0) && (zs != tg)  &&
        (zn != 0) && (zn != tg)) )
  {
    if ( (zb == tg) && (qb == 0) ) // from Z-
    {
      sd = find_mode_id(tgt_id, qw, qe, qs, qn, qb, qt);
      mid[m_p] = sd;
      setFaceBID(bid[m_b], Z_PLUS, sd);
      cut[_F_IDX_S4DEX(Z_PLUS, i, j, k-1, 6, ix, jx, kx, gd)] = cpos;
      replaced++;
    }
    else if ( (zt == tg) && (qt == 0) ) // from Z+
    {
      sd = find_mode_id(tgt_id, qw, qe, qs, qn, qb, qt);
      mid[m_p] = sd;
      setFaceBID(bid[m_t], Z_MINUS, sd);
      cut[_F_IDX_S4DEX(Z_MINUS, i, j, k+1, 6, ix, jx, kx, gd)] = cpos;
      replaced++;
    }
  }
  
} // target