//##################################################################################
//
// Flow Base class
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
 * @file   FindexS3D.h
 * @brief  アドレス計算の共通ヘッダ
 * @author aics
 */

  // int ix, jx, kx, gd;

  size_t m_w = _F_IDX_S3D(i-1, j,   k,   ix, jx, kx, gd);
  size_t m_e = _F_IDX_S3D(i+1, j,   k,   ix, jx, kx, gd);
  size_t m_s = _F_IDX_S3D(i,   j-1, k,   ix, jx, kx, gd);
  size_t m_n = _F_IDX_S3D(i,   j+1, k,   ix, jx, kx, gd);
  size_t m_b = _F_IDX_S3D(i,   j,   k-1, ix, jx, kx, gd);
  size_t m_t = _F_IDX_S3D(i,   j,   k+1, ix, jx, kx, gd);
