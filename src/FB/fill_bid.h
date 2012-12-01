// #################################################################
//
// CAERU Library
//
// Copyright (c) 2012 All right reserved.
//
// Institute of Industrial Science, University of Tokyo, Japan. 
//
// #################################################################

/**
 * @file   fill_bid.h
 * @brief  bidによるフィルアルゴのコア
 * @author kero
 */
 
m_p = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);

m_e = _F_IDX_S3D(i+1, j,   k,   ix, jx, kx, gd);
m_w = _F_IDX_S3D(i-1, j,   k,   ix, jx, kx, gd);
m_n = _F_IDX_S3D(i,   j+1, k,   ix, jx, kx, gd);
m_s = _F_IDX_S3D(i,   j-1, k,   ix, jx, kx, gd);
m_t = _F_IDX_S3D(i,   j,   k+1, ix, jx, kx, gd);
m_b = _F_IDX_S3D(i,   j,   k-1, ix, jx, kx, gd);

zp = mid[m_p];
zw = mid[m_w];
ze = mid[m_e];
zs = mid[m_s];
zn = mid[m_n];
zb = mid[m_b];
zt = mid[m_t];

// 未ペイントの場合にテスト
if ( zp == 0 )
{
  qq = bid[m_p];
  
  // 隣接セルの方向に対するカットIDの有無>> 0ならばカット無し、チェック半径は1
  qw = get_BID5(X_MINUS, qq);
  qe = get_BID5(X_PLUS,  qq);
  qs = get_BID5(Y_MINUS, qq);
  qn = get_BID5(Y_PLUS,  qq);
  qb = get_BID5(Z_MINUS, qq);
  qt = get_BID5(Z_PLUS,  qq);
  
  // 各方向のテスト
  // 例えば、W側をみてテストするとき、次の2つを満たす場合のみがフィルの候補となる
  // 1) W側のセルがフィルを行うターゲットID（Fluid）で既にペイントされている
  // 2) 対象セルのカットIDのW方向がゼロ、つまりカットがない
  
  if ( (zw == tg) && (qw == 0) )
  {
    // W方向からみる場合、面直平面のY, Z各方向の対辺の両方にカットがある場合は塞ぐ
    if ( (qs*qn != 0) || (qb*qt != 0) )
    { 
      mid[m_p] = sd; // セルIDを固体に変更
      set_BID5(bid[m_w], X_PLUS, sd); // テストする方向からみて、カットIDを設定
      cut[_F_IDX_S4DEX(X_PLUS, i-1, j, k, 6, ix, jx, kx, gd)] = cpos; // カット位置をセット
    }
    else // フィルする
    {
      mid[m_p] = tg;
      c++;
    }
  }
  else if ( (zs == tg) && (qs == 0) )
  {
    if ( (qw*qe != 0) || (qb*qt != 0) )
    {
      mid[m_p] = sd;
      set_BID5(bid[m_s], Y_PLUS, sd);
      cut[_F_IDX_S4DEX(Y_PLUS, i, j-1, k, 6, ix, jx, kx, gd)] = cpos;
    }
    else
    {
      mid[m_p] = tg;
      c++;
    }
  }
  else if ( (zb == tg) && (qb == 0) )
  {
    if ( (qw*qe != 0) || (qs*qn != 0) )
    {
      mid[m_p] = sd;
      set_BID5(bid[m_b], Z_PLUS, sd);
      cut[_F_IDX_S4DEX(Z_PLUS, i, j, k-1, 6, ix, jx, kx, gd)] = cpos;
    }
    else
    {
      mid[m_p] = tg;
      c++;
    }
  }
  else if ( (zt == tg) && (qt == 0) )
  {
    if ( (qw*qe != 0) || (qs*qn != 0) )
    {
      mid[m_p] = sd;
      set_BID5(bid[m_t], Z_MINUS, sd);
      cut[_F_IDX_S4DEX(Z_MINUS, i, j, k+1, 6, ix, jx, kx, gd)] = cpos;
    }
    else
    {
      mid[m_p] = tg;
      c++;
    }
  }
  else if ( (zn == tg) && (qn == 0) )
  {
    if ( (qw*qe != 0) || (qb*qt != 0) )
    {
      mid[m_p] = sd;
      set_BID5(bid[m_n], Y_MINUS, sd);
      cut[_F_IDX_S4DEX(Y_MINUS, i, j+1, k, 6, ix, jx, kx, gd)] = cpos;
    }
    else
    {
      mid[m_p] = tg;
      c++;
    }
  }
  else if ( (ze == tg) && (qe == 0) )
  {
    if ( (qs*qn != 0) || (qb*qt != 0) )
    {
      mid[m_p] = sd;
      set_BID5(bid[m_e], X_MINUS, sd);
      cut[_F_IDX_S4DEX(X_MINUS, i+1, j, k, 6, ix, jx, kx, gd)] = cpos;
    }
    else
    {
      mid[m_p] = tg;
      c++;
    }
  }
  
} // target