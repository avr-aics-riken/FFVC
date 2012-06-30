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

  m_e = _F_IDX_S3D(i+1, j,   k,   ix, jx, kx, gd); //FBUtility::getFindexS3D(sz, gd, i+1, j  , k  );
  m_w = _F_IDX_S3D(i-1, j,   k,   ix, jx, kx, gd); //FBUtility::getFindexS3D(sz, gd, i-1, j  , k  );
  m_n = _F_IDX_S3D(i,   j+1, k,   ix, jx, kx, gd); //FBUtility::getFindexS3D(sz, gd, i  , j+1, k  );
  m_s = _F_IDX_S3D(i,   j-1, k,   ix, jx, kx, gd); //FBUtility::getFindexS3D(sz, gd, i  , j-1, k  );
  m_t = _F_IDX_S3D(i,   j,   k+1, ix, jx, kx, gd); //FBUtility::getFindexS3D(sz, gd, i  , j  , k+1);
  m_b = _F_IDX_S3D(i,   j,   k-1, ix, jx, kx, gd); //FBUtility::getFindexS3D(sz, gd, i  , j  , k-1);