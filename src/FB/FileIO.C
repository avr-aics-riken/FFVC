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
void FileIO::writeRawSPH(const REAL_TYPE *vf, const int* sz, const int gc, const int gc_out, const REAL_TYPE* org, const REAL_TYPE* ddx, const int m_ModePrecision)
{
  int pad, dType, stp, svType;
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
    printf("\tCan't open %s file\n", sph_fname);
    Exit(0);
  }
  
  int m_sz[3];
  m_sz[0] = sz[0]+2*gc_out;
  m_sz[1] = sz[1]+2*gc_out;
  m_sz[2] = sz[2]+2*gc_out;
  int gd = gc;
  
  size_t nx = m_sz[0] * m_sz[1] * m_sz[2];
  
  ox = org[0]-ddx[0]*(REAL_TYPE)gc_out;
  oy = org[1]-ddx[1]*(REAL_TYPE)gc_out;
  oz = org[2]-ddx[2]*(REAL_TYPE)gc_out;
  dx = ddx[0];
  dy = ddx[1];
  dz = ddx[2];
  //printf("org: %f %f %f\n", ox, oy, oz);
  //printf("dx : %f %f %f\n", dx, dy, dz);
  
  svType = kind_scalar;
  if ( sizeof(REAL_TYPE) == sizeof(double) )
  {
    for (int i=0; i<3; i++)   szl[i] = (long long)m_sz[i];
  }
  
  REAL_TYPE *f = new REAL_TYPE[nx];
  
  size_t m, l;
  
  for (int k=1; k<=m_sz[2]; k++) {
    for (int j=1; j<=m_sz[1]; j++) {
      for (int i=1; i<=m_sz[0]; i++) {
        l = _F_IDX_S3D(i, j, k, m_sz[0], m_sz[1], m_sz[2], gc_out);
        m = _F_IDX_S3D(i, j, k, m_sz[0], m_sz[1], m_sz[2], gd);
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
    ofs.write( (char*)&m_sz[0], sizeof(int) );
    ofs.write( (char*)&m_sz[1], sizeof(int) );
    ofs.write( (char*)&m_sz[2], sizeof(int) );
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
    int pad = (m_ModePrecision == sizeof(float)) ? nx * sizeof(float) : nx * sizeof(double);
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

// #################################################################
// svxフォーマットで出力する(ID)
void FileIO::writeSVX(int* mid,
                      const int ip,
                      const int jp,
                      const int kp,
                      const int m_sz[3],
                      const int m_gd,
                      const float m_pit[3],
                      const float m_org[3]
                      )
{
  char svx_fname[512];
  sprintf( svx_fname, "subc_%d_%d_%d.svx", ip,jp,kp );

  ofstream ofs(svx_fname, ios::out | ios::binary);
  if (!ofs) {
    cout << "\tCan't open " << svx_fname << " file" << endl;
    Exit(0);
  }
  
  int imax = m_sz[0];
  int jmax = m_sz[1];
  int kmax = m_sz[2];
  int gd = m_gd;
  
  int ix = imax+2*gd;  // guide cell
  int jx = jmax+2*gd;
  int kx = kmax+2*gd;
  
  size_t nx = (size_t)(ix*jx*kx);
  
  float dx = m_pit[0];
  float dy = m_pit[1];
  float dz = m_pit[2];
  float ox = m_org[0] - dx*gd;
  float oy = m_org[1] - dy*gd;
  float oz = m_org[2] - dz*gd;
  
  int* q = mid;
  
  /*
  int* q = new int[nx];
  
  for (int k=0; k<=(kmax+1); k++) {
    for (int j=0; j<=(jmax+1); j++) {
      for (int i=0; i<=(imax+1); i++) {
        l = _F_IDX_S3D(i, j, k, imax, jmax, kmax, gd);
        m = _F_IDX_S3D(i, j, k, imax, jmax, kmax, gd);
        q[l] = mid[m];
      }
    }
  }
   */
  
  // voxel size
  int sz = sizeof(int)*3;
  ofs.write( (char*)&sz, sizeof(int) );
  ofs.write( (char*)&ix, sizeof(int) );
  ofs.write( (char*)&jx, sizeof(int) );
  ofs.write( (char*)&kx, sizeof(int) );
  ofs.write( (char*)&sz, sizeof(int) );
  
  // original point of domain
  sz = sizeof(float)*3;
  ofs.write( (char*)&sz, sizeof(int) );
  ofs.write( (char*)&ox, sizeof(float) );
  ofs.write( (char*)&oy, sizeof(float) );
  ofs.write( (char*)&oz, sizeof(float) );
  ofs.write( (char*)&sz, sizeof(int) );
  
  // pitch of voxel
  ofs.write( (char*)&sz, sizeof(int) );
  ofs.write( (char*)&dx, sizeof(float) );
  ofs.write( (char*)&dy, sizeof(float) );
  ofs.write( (char*)&dz, sizeof(float) );
  ofs.write( (char*)&sz, sizeof(int) );
  
  // type of stored data
  sz = sizeof(int)*1;
  int dtype = 0;
  dtype |= ( 0x1<<2 );  // medium ID
  ofs.write( (char*)&sz,  sizeof(int) );
  ofs.write( (char*)&dtype, sizeof(int) );
  ofs.write( (char*)&sz,  sizeof(int) );
  
  // medium ID
  sz = nx * sizeof(int);
  ofs.write( (char*)&sz, sizeof(int) );
  ofs.write( (char*)q,   sz );
  ofs.write( (char*)&sz, sizeof(int) );
  
  ofs.close();
  
  //if (q) { delete [] q; q=NULL; }
}

