// #################################################################
//
// CAERU Library
//
// Copyright (c) 2012 All right reserved.
//
// Institute of Industrial Science, The University of Tokyo, Japan.
//
// #################################################################

/**
 * @file   fill_mid.h
 * @brief  midによるフィルアルゴのコア
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
if ( zp == tg )
{
  qq = bid[m_p];
  
  // 隣接セルの方向に対するカットIDの有無>> 0ならばカット無し、チェック半径は1
  qw = get_BID5(X_MINUS, qq);
  qe = get_BID5(X_PLUS,  qq);
  qs = get_BID5(Y_MINUS, qq);
  qn = get_BID5(Y_PLUS,  qq);
  qb = get_BID5(Z_MINUS, qq);
  qt = get_BID5(Z_PLUS,  qq);
  
  
  if ( ((zs != 0) && (zs != tg)  &&  // X平面：ゼロでなく、かつ、targetでもない ==> 固体
        (zn != 0) && (zn != tg)) ||
       ((zb != 0) && (zb != tg)  &&
        (zt != 0) && (zt != tg)) )
  {
    
    if ( (zw == tg) && (qw == 0) ) // from X-
    {
      mid[m_p] = sd; // セルIDを固体に変更
      set_BID5(bid[m_w], X_PLUS, sd); // カットIDを設定
      cut[_F_IDX_S4DEX(X_PLUS, i-1, j, k, 6, ix, jx, kx, gd)] = cpos; // カット位置をセット
      c++;
    }
    else if ( (ze == tg) && (qe == 0) ) // from X+
    {
      mid[m_p] = sd;
      set_BID5(bid[m_e], X_MINUS, sd);
      cut[_F_IDX_S4DEX(X_MINUS, i+1, j, k, 6, ix, jx, kx, gd)] = cpos;
      c++;
    }
  }
  
  if ( ((zw != 0) && (zw != tg)  && // Y平面
        (ze != 0) && (ze != tg)) ||
       ((zb != 0) && (zb != tg)  &&
        (zt != 0) && (zt != tg)) )
  {
    
    if ( (zs == tg) && (qs == 0) ) // from Y-
    {
      mid[m_p] = sd;
      set_BID5(bid[m_s], Y_PLUS, sd);
      cut[_F_IDX_S4DEX(Y_PLUS, i, j-1, k, 6, ix, jx, kx, gd)] = cpos;
      c++;
    }
    else if ( (zn == tg) && (qn == 0) ) // from Y+
    {
      mid[m_p] = sd;
      set_BID5(bid[m_n], Y_MINUS, sd);
      cut[_F_IDX_S4DEX(Y_MINUS, i, j+1, k, 6, ix, jx, kx, gd)] = cpos;
      c++;
    }
  }
  
  if ( ((zw != 0) && (zw != tg)  && // Z平面
        (ze != 0) && (ze != tg)) ||
       ((zs != 0) && (zs != tg)  &&
        (zn != 0) && (zn != tg)) )
  {
    if ( (zb == tg) && (qb == 0) ) // from Z-
    {
      mid[m_p] = sd;
      set_BID5(bid[m_b], Z_PLUS, sd);
      cut[_F_IDX_S4DEX(Z_PLUS, i, j, k-1, 6, ix, jx, kx, gd)] = cpos;
      c++;
    }
    else if ( (zt == tg) && (qt == 0) ) // from Z+
    {
      mid[m_p] = sd;
      set_BID5(bid[m_t], Z_MINUS, sd);
      cut[_F_IDX_S4DEX(Z_MINUS, i, j, k+1, 6, ix, jx, kx, gd)] = cpos;
      c++;
    }
  }
  
} // target