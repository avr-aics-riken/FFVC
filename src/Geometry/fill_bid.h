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
 * @file   fill_bid.h
 * @brief  bidによるフィルアルゴのコア
 * @author aics
 */


/// 最初に1セルの隙間の穴埋め処理を行い，その後，フラッドフィル

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
  
  
  // 最頻値ID >> tg以外
  int sd = FBUtility::find_mode_id(tg, qw, qe, qs, qn, qb, qt, m_NoCompo);

  
  // 各方向のテスト
  // 対象セルは，未評価セル(zp == 0) 
  // 例えば、X-側からテストする場合は，X-側のセルがフィルを行うターゲットIDで既にペイントされている(zw == tg)
  // 次のパターンの場合にX-方向からみてi点上にsolidがあるとみなす
  // 1) X方向面直平面のYもしくはZ各方向の対辺の両方のみにカットがある
  // 2) YとZの4方向にカットがある

  // 固体IDでフィルする場合にインクリメント
  int sd_flag = 0;
  
  
  
  // X-
  if ( zw==tg && qw==0 )
  {
    // 検査方向に壁をおく場合にインクリメント
    int fill_flag = 0;
    
    // qs+qn==0は両方とも流体、qb*qt>0は両方とも物体を意味する
    if ( ((qs+qn)==0) && ((qb*qt)>0) && (qe==0) ) fill_flag++; // Y方向の両側にはカットが無く，Z方向の両側にカットがある，X+方向にカットなし
    if ( ((qb+qt)==0) && ((qs*qn)>0) && (qe==0) ) fill_flag++; // Z方向の両側にはカットが無く，Y方向の両側にカットがある，X+方向にカットなし
    if ( ((qs*qn*qb*qt)>0)           && (qe==0) ) fill_flag++; // Y方向とZ方向の両側にカットがある，X+方向にカットなし
  
    if ( fill_flag > 0 )
    {
      if (sd == 0) Exit(0);
        
      cut[ _F_IDX_S4DEX(X_plus, i-1, j, k, 6, ix, jx, kx, gd) ] = 1.0;
      setBit5(bid[m_w], sd, X_plus);
      
      cut[ _F_IDX_S4DEX(X_minus, i, j, k, 6, ix, jx, kx, gd) ] = 0.0;
      setBit5(bid[m_p], sd, X_minus);
      sd_flag++;
    }
  }
  
  
  // X+
  if ( ze==tg && qe==0 )
  {
    int fill_flag = 0;
    
    if ( ((qs+qn)==0) && ((qb*qt)>0) && (qw==0) ) fill_flag++;
    if ( ((qb+qt)==0) && ((qs*qn)>0) && (qw==0) ) fill_flag++;
    if ( ((qs*qn*qb*qt)>0)           && (qw==0) ) fill_flag++;
    
    if ( fill_flag > 0 )
    {
      if (sd == 0) Exit(0);
        
      cut[ _F_IDX_S4DEX(X_minus, i+1, j, k, 6, ix, jx, kx, gd) ] = 1.0;
      setBit5(bid[m_e], sd, X_minus);
      
      cut[ _F_IDX_S4DEX(X_plus, i, j, k, 6, ix, jx, kx, gd) ] = 0.0;
      setBit5(bid[m_p], sd, X_plus);
      sd_flag++;
    }
  }
  
  
  // Y-
  if ( zs==tg && qs==0 )
  {
    int fill_flag = 0;
    
    if ( ((qw+qe)==0) && ((qb*qt)>0) && (qn==0) ) fill_flag++;
    if ( ((qb+qt)==0) && ((qw*qe)>0) && (qn==0) ) fill_flag++;
    if ( ((qw*qe*qb*qt)>0)           && (qn==0) ) fill_flag++;
    
    if ( fill_flag > 0 )
    {
      if (sd == 0) Exit(0);
        
        cut[ _F_IDX_S4DEX(Y_plus, i, j-1, k, 6, ix, jx, kx, gd) ] = 1.0;
        setBit5(bid[m_s], sd, Y_plus);
        
        cut[ _F_IDX_S4DEX(Y_minus, i, j, k, 6, ix, jx, kx, gd) ] = 0.0;
        setBit5(bid[m_p], sd, Y_minus);
        sd_flag++;
    }
  }

  
  // Y+
  if ( zn==tg && qn==0 )
  {
    int fill_flag = 0;
    
    if ( ((qw+qe)==0) && ((qb*qt)>0) && (qs==0) ) fill_flag++;
    if ( ((qb+qt)==0) && ((qw*qe)>0) && (qs==0) ) fill_flag++;
    if ( ((qw*qe*qb*qt)>0)           && (qs==0) ) fill_flag++;
    
    if ( fill_flag > 0 )
    {
      if (sd == 0) Exit(0);
        
        cut[ _F_IDX_S4DEX(Y_minus, i, j+1, k, 6, ix, jx, kx, gd) ] = 1.0;
        setBit5(bid[m_n], sd, Y_minus);
        
        cut[ _F_IDX_S4DEX(Y_plus, i, j, k, 6, ix, jx, kx, gd) ] = 0.0;
        setBit5(bid[m_p], sd, Y_plus);
        sd_flag++;
    }
  }
  
  
  // Z-
  if ( zb==tg && qb==0 )
  {
    int fill_flag = 0;
    
    if ( ((qw+qe)==0) && ((qs*qn)>0) && (qt==0) ) fill_flag++;
    if ( ((qs+qn)==0) && ((qw*qe)>0) && (qt==0) ) fill_flag++;
    if ( ((qw*qe*qs*qn)>0)           && (qt==0) ) fill_flag++;
    
    if ( fill_flag > 0 )
    {
      if (sd == 0) Exit(0);
        
      cut[ _F_IDX_S4DEX(Z_plus, i, j, k-1, 6, ix, jx, kx, gd) ] = 1.0;
      setBit5(bid[m_b], sd, Z_plus);
      
      cut[ _F_IDX_S4DEX(Z_minus, i, j, k, 6, ix, jx, kx, gd) ] = 0.0;
      setBit5(bid[m_p], sd, Z_minus);
      sd_flag++;
    }
  }
  
  
  // Z+
  if ( zt==tg && qt==0 )
  {
    int fill_flag = 0;
    
    if ( ((qw+qe)==0) && ((qs*qn)>0) && (qb==0) ) fill_flag++;
    if ( ((qs+qn)==0) && ((qw*qe)>0) && (qb==0) ) fill_flag++;
    if ( ((qw*qe*qs*qn)>0)           && (qb==0) ) fill_flag++;
    
    if ( fill_flag > 0 )
    {
      if (sd == 0) Exit(0);
        
      cut[ _F_IDX_S4DEX(Z_minus, i, j, k+1, 6, ix, jx, kx, gd) ] = 1.0;
      setBit5(bid[m_t], sd, Z_minus);
      
      cut[ _F_IDX_S4DEX(Z_plus, i, j, k, 6, ix, jx, kx, gd) ] = 0.0;
      setBit5(bid[m_p], sd, Z_plus);
      sd_flag++;
    }
  }

  
  /* ?
   long long pos = cut[m_p];
   
  if ( sd_flag > 0 )
  {
    setMediumID(bcd[m_p], sd);
    
    for (int l=0; l<6; l++)
    {
      pos[l] = 0.0;
      setBit5(bid[m_p], sd, l);
    }
    replaced++;
  }
   */
  
  
  // フラッドフィル
  int tag = 0;
  
  // mode_x==0の時には，X方向の領域境界でフィルしない
  if ( (sdw < 0) && (i == 1) && !mode_x )
  {
    ; // skip
  }
  else
  {
    if ( zw==tg && qw==0 ) tag++;
  }
  
  if ( (sde < 0) && (i == ix) && !mode_x )
  {
    ; // skip
  }
  else
  {
    if ( ze==tg && qe==0 ) tag++;
  }
  
  // mode_y==0の時には，Y方向の領域境界でフィルしない
  if ( (sds < 0) && (j == 1) && !mode_y )
  {
    ; // skip
  }
  else
  {
    if ( zs==tg && qs==0 ) tag++;
  }
  
  if ( (sdn < 0) && (j == jx) && !mode_y )
  {
    ;
  }
  else
  {
    if ( zn==tg && qn==0 ) tag++;
  }
  
  // mode_z==0の時には，Z方向の領域境界でフィルしない
  if ( (sdb < 0) && (k == 1) && !mode_z )
  {
    ; // skip
  }
  else
  {
    if ( zb==tg && qb==0 ) tag++;
  }
  
  if ( (sdt < 0) && (k == kx) && !mode_z )
  {
    ; // skip
  }
  else
  {
    if ( zt==tg && qt==0 ) tag++;
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
