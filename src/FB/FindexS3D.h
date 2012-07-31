// #################################################################
//
// CAERU Library
//
// Copyright (c) All right reserved. 2012
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
//
// #################################################################

/**
 * @file   FindexS3D.h
 * @brief  アドレス計算の共通ヘッダ
 * @author kero
 */

  // int ix, jx, kx, gd;
  // size_t m_p, m_e, m_w, m_n, m_s, m_t, m_b;

  m_e = _F_IDX_S3D(i+1, j,   k,   ix, jx, kx, gd);
  m_w = _F_IDX_S3D(i-1, j,   k,   ix, jx, kx, gd);
  m_n = _F_IDX_S3D(i,   j+1, k,   ix, jx, kx, gd);
  m_s = _F_IDX_S3D(i,   j-1, k,   ix, jx, kx, gd);
  m_t = _F_IDX_S3D(i,   j,   k+1, ix, jx, kx, gd);
  m_b = _F_IDX_S3D(i,   j,   k-1, ix, jx, kx, gd);