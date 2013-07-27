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
 * @file   fill_bid.h
 * @brief  bidによるフィルアルゴのコア
 * @author kero
 */

float cpos = 0.5;

size_t m_p = _F_IDX_S3D(i  , j  , k  , ix, jx, kx, gd);
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
if ( zp == 0 )
{
  int qq = bid[m_p];
  int sd;
  
  // 隣接セルの方向に対するカットIDの有無>> 0ならばカット無し、チェック半径は1
  int qw = get_BID5(X_MINUS, qq);
  int qe = get_BID5(X_PLUS,  qq);
  int qs = get_BID5(Y_MINUS, qq);
  int qn = get_BID5(Y_PLUS,  qq);
  int qb = get_BID5(Z_MINUS, qq);
  int qt = get_BID5(Z_PLUS,  qq);
  
  // 流体属性をもつCELL_MONITORについては，カットがないものとしてフィル対象とする
  int ff = 0;
  for (int n=1; n<=NoCompo; n++) {
    if ( list[n] == qw ) ff++; // CELL_MONITORであれば ff>0
  }
  if ( ff > 0 ) qw=0; // カットなしとする
  
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
    
  
  // 各方向のテスト
  // 例えば、W側をみてテストするとき、次の2つを満たす場合のみがフィルの候補となる
  // 1) W側のセルがフィルを行うターゲットID（Fluid）で既にペイントされている(zw == tg)
  // 2) 対象セルのカットIDのW方向がゼロ(qw==0)、つまりカットがない
  
  if ( (zw == tg) && (qw == 0) )
  {
    // W方向からみる場合、面直平面のY, Z各方向の対辺の両方にカットがある場合は塞ぐ
    if ( (qs*qn != 0) || (qb*qt != 0) )
    {
      sd = find_mode_id(tgt_id, qw, qe, qs, qn, qb, qt);
      mid[m_p] =  sd;
      set_BID5(bid[m_w], X_PLUS, sd); // テストする方向からみて、カットIDを設定
      cut[_F_IDX_S4DEX(X_PLUS, i-1, j, k, 6, ix, jx, kx, gd)] = cpos; // カット位置をセット
      replaced++;
    }
    else // 流体IDでフィルする
    {
      mid[m_p] = tg;
      filled++;
    }
  }
  else if ( (zs == tg) && (qs == 0) )
  {
    if ( (qw*qe != 0) || (qb*qt != 0) )
    {
      sd = find_mode_id(tgt_id, qw, qe, qs, qn, qb, qt);
      mid[m_p] =  sd;
      set_BID5(bid[m_s], Y_PLUS, sd);
      cut[_F_IDX_S4DEX(Y_PLUS, i, j-1, k, 6, ix, jx, kx, gd)] = cpos;
      replaced++;
    }
    else
    {
      mid[m_p] = tg;
      filled++;
    }
  }
  else if ( (zb == tg) && (qb == 0) )
  {
    if ( (qw*qe != 0) || (qs*qn != 0) )
    {
      sd = find_mode_id(tgt_id, qw, qe, qs, qn, qb, qt);
      mid[m_p] = sd;
      set_BID5(bid[m_b], Z_PLUS, sd);
      cut[_F_IDX_S4DEX(Z_PLUS, i, j, k-1, 6, ix, jx, kx, gd)] = cpos;
      replaced++;
    }
    else
    {
      mid[m_p] = tg;
      filled++;
    }
  }
  else if ( (zt == tg) && (qt == 0) )
  {
    if ( (qw*qe != 0) || (qs*qn != 0) )
    {
      sd = find_mode_id(tgt_id, qw, qe, qs, qn, qb, qt);
      mid[m_p] = sd;
      set_BID5(bid[m_t], Z_MINUS, sd);
      cut[_F_IDX_S4DEX(Z_MINUS, i, j, k+1, 6, ix, jx, kx, gd)] = cpos;
      replaced++;
    }
    else
    {
      mid[m_p] = tg;
      filled++;
    }
  }
  else if ( (zn == tg) && (qn == 0) )
  {
    if ( (qw*qe != 0) || (qb*qt != 0) )
    {
      sd = find_mode_id(tgt_id, qw, qe, qs, qn, qb, qt);
      mid[m_p] = sd;
      set_BID5(bid[m_n], Y_MINUS, sd);
      cut[_F_IDX_S4DEX(Y_MINUS, i, j+1, k, 6, ix, jx, kx, gd)] = cpos;
      replaced++;
    }
    else
    {
      mid[m_p] = tg;
      filled++;
    }
  }
  else if ( (ze == tg) && (qe == 0) )
  {
    if ( (qs*qn != 0) || (qb*qt != 0) )
    {
      sd = find_mode_id(tgt_id, qw, qe, qs, qn, qb, qt);
      mid[m_p] = sd;
      set_BID5(bid[m_e], X_MINUS, sd);
      cut[_F_IDX_S4DEX(X_MINUS, i+1, j, k, 6, ix, jx, kx, gd)] = cpos;
      replaced++;
    }
    else
    {
      mid[m_p] = tg;
      filled++;
    }
  }
  
} // target