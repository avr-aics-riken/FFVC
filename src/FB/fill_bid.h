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
  
  // 隣接セルの方向に対する境界ID
  int qw = getFaceBID(0, qq);
  int qe = getFaceBID(1, qq);
  int qs = getFaceBID(2, qq);
  int qn = getFaceBID(3, qq);
  int qb = getFaceBID(4, qq);
  int qt = getFaceBID(5, qq);
  
  float* pos = &cut[ _F_IDX_S4DEX(0, i, j, k, 6, ix, jx, kx, gd) ];
  
  // 最頻値ID
  int sd = find_mode_id(tg, qw, qe, qs, qn, qb, qt);
  
  
  // 各方向のテスト
  // 対象セルは，未評価セル(zp == 0) 
  // 例えば、X-側からテストする場合は，X-側のセルがフィルを行うターゲットID（Fluid）で既にペイントされている(zw == tg)
  // 次のパターンではない場合にX-方向からみてi点上にsolidがあるとみなす
  // 1) X方向面直平面のYもしくはZ各方向の対辺の両方のみにカットがある
  // 2) YとZの4方向にカットがある

  // 固体IDでフィルする場合にインクリメント
  int sd_flag = 0;
  
  // 検査方向に貫通しない場合にインクリメント
  int fill_flag;
  
  // X-
  fill_flag = 0;
  if ( zw==tg && qw==0 )
  {
    if ( (qs+qn)==0 && (qb*qt)>0 ) fill_flag++;
    if ( (qb+qt)==0 && (qs*qn)>0 ) fill_flag++;
    if ( qs*qn*qb*qt > 0 )         fill_flag++;
  }
  
  if ( fill_flag > 0 )
  {
    cut[ _F_IDX_S4DEX(X_PLUS, i-1, j, k, 6, ix, jx, kx, gd) ] = 1.0;
    setFaceBID(bid[m_w], X_PLUS, sd);
    sd_flag++;
  }
  
  // X+
  fill_flag = 0;
  if ( ze==tg && qe==0 )
  {
    if ( (qs+qn)==0 && (qb*qt)>0 ) fill_flag++;
    if ( (qb+qt)==0 && (qs*qn)>0 ) fill_flag++;
    if ( qs*qn*qb*qt > 0 )         fill_flag++;
  }
  
  if ( fill_flag > 0 )
  {
    cut[ _F_IDX_S4DEX(X_MINUS, i+1, j, k, 6, ix, jx, kx, gd) ] = 1.0;
    setFaceBID(bid[m_e], X_MINUS, sd);
    sd_flag++;
  }
  
  // Y-
  fill_flag = 0;
  if ( zs==tg && qs==0 )
  {
    if ( (qw+qe)==0 && (qb*qt)>0 ) fill_flag++;
    if ( (qb+qt)==0 && (qw*qe)>0 ) fill_flag++;
    if ( qw*qe*qb*qt > 0 )         fill_flag++;
  }
  
  if ( fill_flag > 0 )
  {
    cut[ _F_IDX_S4DEX(Y_PLUS, i, j-1, k, 6, ix, jx, kx, gd) ] = 1.0;
    setFaceBID(bid[m_s], Y_PLUS, sd);
    sd_flag++;
  }

  // Y+
  fill_flag = 0;
  if ( zn==tg && qn==0 )
  {
    if ( (qw+qe)==0 && (qb*qt)>0 ) fill_flag++;
    if ( (qb+qt)==0 && (qw*qe)>0 ) fill_flag++;
    if ( qw*qe*qb*qt > 0 )         fill_flag++;
  }
  
  if ( fill_flag > 0 )
  {
    cut[ _F_IDX_S4DEX(Y_MINUS, i, j+1, k, 6, ix, jx, kx, gd) ] = 1.0;
    setFaceBID(bid[m_n], Y_MINUS, sd);
    sd_flag++;
  }
  
  // Z-
  fill_flag = 0;
  if ( zb==tg && qb==0 )
  {
    if ( (qw+qe)==0 && (qs*qn)>0 ) fill_flag++;
    if ( (qs+qn)==0 && (qw*qe)>0 ) fill_flag++;
    if ( qw*qe*qs*qn > 0 )         fill_flag++;
  }
  
  if ( fill_flag > 0 )
  {
    cut[ _F_IDX_S4DEX(Z_PLUS, i, j, k-1, 6, ix, jx, kx, gd) ] = 1.0;
    setFaceBID(bid[m_b], Z_PLUS, sd);
    sd_flag++;
  }
  
  // Z+
  fill_flag = 0;
  if ( zt==tg && qt==0 )
  {
    if ( (qw+qe)==0 && (qs*qn)>0 ) fill_flag++;
    if ( (qs+qn)==0 && (qw*qe)>0 ) fill_flag++;
    if ( qw*qe*qs*qn > 0 )         fill_flag++;
  }
  
  if ( fill_flag > 0 )
  {
    cut[ _F_IDX_S4DEX(Z_MINUS, i, j, k+1, 6, ix, jx, kx, gd) ] = 1.0;
    setFaceBID(bid[m_t], Z_MINUS, sd);
    sd_flag++;
  }
  
  if ( sd_flag > 0 )
  {
    bcd[m_p] |= sd;
    
    for (int l=0; l<6; l++)
    {
      pos[l] = 0.0;
      setFaceBID(bid[m_p], l, sd);
    }
    replaced++;
  }
  else
  {
    bcd[m_p] |= tg;
    filled++;
  }
  
  
  //if ( cut[_F_IDX_S4DEX(X_PLUS,  i, j, k, 6, ix, jx, kx, gd)] <= ROUND_EPS ) printf("X+ : %d %d %d\n",i,j,k);
  //if ( cut[_F_IDX_S4DEX(X_MINUS, i, j, k, 6, ix, jx, kx, gd)] <= ROUND_EPS ) printf("X- : %d %d %d\n",i,j,k);
  //if ( cut[_F_IDX_S4DEX(Y_PLUS,  i, j, k, 6, ix, jx, kx, gd)] <= ROUND_EPS ) printf("Y+ : %d %d %d\n",i,j,k);
  //if ( cut[_F_IDX_S4DEX(Y_MINUS, i, j, k, 6, ix, jx, kx, gd)] <= ROUND_EPS ) printf("Y- : %d %d %d\n",i,j,k);
  //if ( cut[_F_IDX_S4DEX(Z_PLUS,  i, j, k, 6, ix, jx, kx, gd)] <= ROUND_EPS ) printf("Z+ : %d %d %d\n",i,j,k);
  //if ( cut[_F_IDX_S4DEX(Z_MINUS, i, j, k, 6, ix, jx, kx, gd)] <= ROUND_EPS ) printf("Z- : %d %d %d\n",i,j,k);
  //if ( (i==26) && (j==14) && (k==239) )
    //printf("%d %d %d : %f %f %f %f %f %f\n",i,j,k,
      //     cut[_F_IDX_S4DEX(X_MINUS, i, j, k, 6, ix, jx, kx, gd)],
      //     cut[_F_IDX_S4DEX(X_PLUS, i, j, k, 6, ix, jx, kx, gd)],
      //     cut[_F_IDX_S4DEX(Y_MINUS, i, j, k, 6, ix, jx, kx, gd)],
      //     cut[_F_IDX_S4DEX(Y_PLUS, i, j, k, 6, ix, jx, kx, gd)],
      //     cut[_F_IDX_S4DEX(Z_MINUS, i, j, k, 6, ix, jx, kx, gd)],
      //     cut[_F_IDX_S4DEX(Z_PLUS, i, j, k, 6, ix, jx, kx, gd)]);
    //printf("%d %d %d : %d %d %d %d %d %d\n",i,j,k,
    //       getFaceBID(0, bid[m_p]),
    //       getFaceBID(1, bid[m_p]),
    //       getFaceBID(2, bid[m_p]),
    //       getFaceBID(3, bid[m_p]),
    //       getFaceBID(4, bid[m_p]),
    //       getFaceBID(5, bid[m_p]) );
  
} // target
