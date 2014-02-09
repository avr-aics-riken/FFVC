//##################################################################################
//
// Flow Base class
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2014 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

/** 
 * @file   FileIO.C
 * @brief  FlowBase FileIO class Header
 * @author aics
 */

#include "FileIO.h"



// #################################################################
// sphファイルの書き出し（内部領域のみ）
void FileIO::writeRawSPH(const REAL_TYPE *vf, const int* sz, const int gc, const REAL_TYPE* org, const REAL_TYPE* ddx, const int m_ModePrecision)
{
  int pad, dType, stp, svType;
  int i, j, k;
  REAL_TYPE ox, oy, oz, dx, dy, dz, tm;
  long long szl[3], stpl;
  
  
  char sph_fname[512];
  
  if ( paraMngr->IsParallel() ) 
  {
    sprintf( sph_fname, "field%010d.sph", paraMngr->GetMyRankID() );
  } 
  else 
  {
    sprintf( sph_fname, "field.sph" );
  }
  
  ofstream ofs(sph_fname, ios::out | ios::binary);
  if (!ofs)
  {
    cout << "\tCan't open " << sph_fname << " file" << endl;
    Exit(0);
  }
  
  int ix = sz[0]; //+2*gc;
  int jx = sz[1]; //+2*gc;
  int kx = sz[2]; //+2*gc;
  int gd = gc;
  
  size_t nx = ix * jx * kx;
  
  ox = org[0]; //-ddx[0]*(REAL_TYPE)gc;
  oy = org[1]; //-ddx[1]*(REAL_TYPE)gc;
  oz = org[2]; //-ddx[2]*(REAL_TYPE)gc;
  dx = ddx[0];
  dy = ddx[1];
  dz = ddx[2];
  //printf("org: %f %f %f\n", ox, oy, oz);
  //printf("dx : %f %f %f\n", dx, dy, dz);
  
  svType = kind_scalar;
  if ( sizeof(REAL_TYPE) == sizeof(double) )
  {
    for (i=0; i<3; i++)   szl[i] = (long long)sz[i];
  }
  
  REAL_TYPE *f = new REAL_TYPE[nx];
  
  size_t m, l;
  
  for (k=1; k<=kx; k++) {
    for (j=1; j<=jx; j++) {
      for (i=1; i<=ix; i++) {
        l = _F_IDX_S3D(i, j, k, ix, jx, kx, 0);
        m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        f[l] = (REAL_TYPE)vf[m];
      }
    }
  }
  
  // data property
  ( m_ModePrecision == sizeof(float) ) ? dType=1 : dType=2;
  pad = sizeof(int)*2;
  ofs.write( (char*)&pad, sizeof(int) );
  ofs.write( (char*)&svType, sizeof(int) );
  ofs.write( (char*)&dType, sizeof(int) );
  ofs.write( (char*)&pad, sizeof(int) );
  
  // voxel size
  if (dType == 1) {
    pad = sizeof(int)*3;
    ofs.write( (char*)&pad, sizeof(int) );
    ofs.write( (char*)&ix, sizeof(int) );
    ofs.write( (char*)&jx, sizeof(int) );
    ofs.write( (char*)&kx, sizeof(int) );
    ofs.write( (char*)&pad, sizeof(int) );
  }
  else {
    pad = sizeof(long long)*3;
    ofs.write( (char*)&pad, sizeof(int) );
    ofs.write( (char*)&szl[0], sizeof(long long) );
    ofs.write( (char*)&szl[1], sizeof(long long) );
    ofs.write( (char*)&szl[2], sizeof(long long) );
    ofs.write( (char*)&pad, sizeof(int) );
  }
  
  // original point of domain
  if (dType == 1) {
    pad = sizeof(float)*3;
    ofs.write( (char*)&pad, sizeof(int) );
    ofs.write( (char*)&ox, sizeof(float) );
    ofs.write( (char*)&oy, sizeof(float) );
    ofs.write( (char*)&oz, sizeof(float) );
    ofs.write( (char*)&pad, sizeof(int) );
  }
  else {
    pad = sizeof(double)*3;
    ofs.write( (char*)&pad, sizeof(int) );
    ofs.write( (char*)&ox, sizeof(double) );
    ofs.write( (char*)&oy, sizeof(double) );
    ofs.write( (char*)&oz, sizeof(double) );
    ofs.write( (char*)&pad, sizeof(int) );
  }
  
  // pitch of voxel
  if (dType == 1) {
    pad = sizeof(float)*3;
    ofs.write( (char*)&pad, sizeof(int) );
    ofs.write( (char*)&dx, sizeof(float) );
    ofs.write( (char*)&dy, sizeof(float) );
    ofs.write( (char*)&dz, sizeof(float) );
    ofs.write( (char*)&pad, sizeof(int) );
  }
  else {
    pad = sizeof(double)*3;
    ofs.write( (char*)&pad, sizeof(int) );
    ofs.write( (char*)&dx, sizeof(double) );
    ofs.write( (char*)&dy, sizeof(double) );
    ofs.write( (char*)&dz, sizeof(double) );
    ofs.write( (char*)&pad, sizeof(int) );
  }
  
  // time stamp
  if (dType == 1) {
    stp = 0;
    tm = 0.0;
    pad = sizeof(int)+sizeof(float);
    ofs.write( (char*)&pad, sizeof(int) );
    ofs.write( (char*)&stp, sizeof(int) );
    ofs.write( (char*)&tm, sizeof(float) );
    ofs.write( (char*)&pad, sizeof(int) );
  }
  else {
    stpl =0;
    tm = 0.0;
    pad = sizeof(long long)+sizeof(double);
    ofs.write( (char*)&pad, sizeof(int) );
    ofs.write( (char*)&stpl, sizeof(long long) );
    ofs.write( (char*)&tm, sizeof(double) );
    ofs.write( (char*)&pad, sizeof(int) );
  }
  
  if (svType == kind_scalar) {
    pad = (m_ModePrecision == sizeof(float)) ? nx * sizeof(float) : nx * sizeof(double);
    ofs.write( (char*)&pad, sizeof(int) );
    ofs.write( (char*)f,   pad );
    ofs.write( (char*)&pad, sizeof(int) );
  }
  
  ofs.close();
  
  if (f) { delete [] f; f=NULL; }
}


// #################################################################
// スカラーファイルを出力する
void FileIO::writeScalar(const string fname,
                         int* sz, 
                         int gc,
                         REAL_TYPE* s, 
                         const unsigned step, 
                         const double time, 
                         const REAL_TYPE* org, 
                         const REAL_TYPE* pit, 
                         const int guide_out,
                         const bool mode,
                         const unsigned step_avr,
                         const double time_avr)
{
  if ( fname.empty() ) Exit(0);
  
  if ( fname.size() > FB_FILE_PATH_LENGTH ) {
    Hostonly_ printf ("Error : Length of file path is greater than %d\n", FB_FILE_PATH_LENGTH);
    Exit(0);
  }
  
  
  char tmp[FB_FILE_PATH_LENGTH];
  memset(tmp, 0, sizeof(char)*FB_FILE_PATH_LENGTH);
  strcpy(tmp, fname.c_str());

  int stp = (int)step;
  REAL_TYPE tm = (REAL_TYPE)time;
  int g = guide_out;

  REAL_TYPE o[3] = {org[0], org[1], org[2]};
  REAL_TYPE p[3] = {pit[0], pit[1], pit[2]};
  
  int avs = (mode == true) ? 1 : 0;
  int stp_a = (int)step_avr;
  REAL_TYPE tm_a = (REAL_TYPE)time_avr;
  int d_type = (sizeof(REAL_TYPE) == 4) ? 1 : 2;  // 1-float / 2-double
  
  fb_write_sph_s_ (s, sz, &gc, tmp, &stp, &tm, o, p, &d_type, &g, &avs, &stp_a, &tm_a);
  
}


// #################################################################
// ベクトルファイルを出力する
void FileIO::writeVector(const string fname,
                         int* sz, 
                         int gc, 
                         REAL_TYPE* v, 
                         const unsigned step, 
                         const double time, 
                         const REAL_TYPE* org, 
                         const REAL_TYPE* pit, 
                         const int guide_out,
                         const bool mode,
                         const unsigned step_avr,
                         const double time_avr)
{
  if ( fname.empty() ) Exit(0);
  
  if ( fname.size() > FB_FILE_PATH_LENGTH ) {
    Hostonly_ printf ("Error : Length of file path is greater than %d\n", FB_FILE_PATH_LENGTH);
    Exit(0);
  }
  
  char tmp[FB_FILE_PATH_LENGTH];
  memset(tmp, 0, sizeof(char)*FB_FILE_PATH_LENGTH);
  strcpy(tmp, fname.c_str());
  
  int stp = (int)step;
  REAL_TYPE tm = (REAL_TYPE)time;
  int g = guide_out;
  
  REAL_TYPE o[3] = {org[0], org[1], org[2]};
  REAL_TYPE p[3] = {pit[0], pit[1], pit[2]};
  
  int avs = (mode == true) ? 1 : 0;
  int stp_a = (int)step_avr;
  REAL_TYPE tm_a = (REAL_TYPE)time_avr;
  int d_type = (sizeof(REAL_TYPE) == 4) ? 1 : 2;  // 1-float / 2-double
  
  fb_write_sph_v_ (v, sz, &gc, tmp, &stp, &tm, o, p, &d_type, &g, &avs, &stp_a, &tm_a);

}
