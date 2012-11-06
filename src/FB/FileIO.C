// #################################################################
//
// FFV : Frontflow / violet
//
// Copyright (c) 2012 All right reserved.
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
//
// #################################################################

/** 
 * @file   FileIO.C
 * @brief  FlowBase FileIO class Header
 * @author kero
 */

#include "FileIO.h"



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




// 圧力のファイルをロードする
void FileIO::readPressure(FILE* fp,
                          const string fname,
                          int* sz,
                          int gc,
                          REAL_TYPE* p,
                          unsigned& step,
                          double& time,
                          const int Dmode,
                          const REAL_TYPE BasePrs,
                          const REAL_TYPE RefDensity,
                          const REAL_TYPE RefVelocity,
                          double& flop,
                          const int guide_out,
                          const bool mode,
                          unsigned& step_avr,
                          double& time_avr
                          )
{
  
  if ( fname.empty() ) Exit(0);
  
  if ( fname.size() > FB_FILE_PATH_LENGTH )
  {
    Hostonly_ printf ("Error : Length of file path is greater than %d\n", FB_FILE_PATH_LENGTH);
    Exit(0);
  }
  

  char tmp[FB_FILE_PATH_LENGTH];
  memset(tmp, 0, sizeof(char)*FB_FILE_PATH_LENGTH);
  strcpy(tmp, fname.c_str());

  int g = guide_out;
  int avs = (mode == true) ? 1 : 0;
  int f_step, a_step;
  REAL_TYPE f_time, a_time;
  
  fb_read_sph_s_ (p, sz, &gc, tmp, &f_step, &f_time, &g, &avs, &a_step, &a_time);
  
  if ( !mode )
  {
    if ( (step_avr == 0) || (time_avr <= 0.0) )
    {
      Hostonly_ printf ("Error : restarted step[%d] or time[%e] is invalid\n", a_step, a_time);
      Exit(0);
    }
  }

  
  // 有次元ファイルの場合，無次元に変換する
  if ( Dmode == DIMENSIONAL )
  {
    REAL_TYPE scale = (mode == true) ? 1.0 : (REAL_TYPE)step_avr; // 瞬時値の時スケールは1.0、平均値の時は平均数
  
    U.prs_array_D2ND(p, sz, gc, BasePrs, RefDensity, RefVelocity, scale, flop);
  }
  
  if ( mode )
  {
    Hostonly_ printf     ("\t[%s] has read :\tstep=%d  time=%e [%s]\n", tmp, f_step, f_time, (Dmode==DIMENSIONAL)?"sec.":"-");
    Hostonly_ fprintf(fp, "\t[%s] has read :\tstep=%d  time=%e [%s]\n", tmp, f_step, f_time, (Dmode==DIMENSIONAL)?"sec.":"-");
  }
  else
  {
    Hostonly_ printf     ("\t[%s] has read :\tstep=%d  time=%e \t:Averaged step=%d  time=%e [%s]\n",
                          tmp, f_step, f_time, a_step, a_time, (Dmode==DIMENSIONAL)?"sec.":"-");
    Hostonly_ fprintf(fp, "\t[%s] has read :\tstep=%d  time=%e \t:Averaged step=%d  time=%e [%s]\n", 
                          tmp, f_step, f_time, a_step, a_time, (Dmode==DIMENSIONAL)?"sec.":"-");
  }
  
  step = (unsigned)f_step;
  time = (double)f_time;
  step_avr = (unsigned)a_step;
  time_avr = (double)a_time;

}



// 速度のファイルをロードする
void FileIO::readVelocity(FILE* fp, 
                          const string fname,
                          int* sz, 
                          int gc, 
                          REAL_TYPE* v, 
                          unsigned& step, 
                          double& time, 
                          const REAL_TYPE *v00, 
                          const int Dmode, 
                          const REAL_TYPE RefVelocity, 
                          double& flop, 
                          const int guide_out,
                          const bool mode,
                          unsigned& step_avr,
                          double& time_avr)
{
  
  if ( fname.empty() ) Exit(0);
  
  if ( fname.size() > FB_FILE_PATH_LENGTH )
  {
    Hostonly_ printf ("Error : Length of file path is greater than %d\n", FB_FILE_PATH_LENGTH);
    Exit(0);
  }
  
  char tmp[FB_FILE_PATH_LENGTH];
  memset(tmp, 0, sizeof(char)*FB_FILE_PATH_LENGTH);
  strcpy(tmp, fname.c_str());

  int g = guide_out;
  int avs = (mode == true) ? 1 : 0;
  int f_step, a_step;
  REAL_TYPE f_time, a_time;
  
  fb_read_sph_v_ (v, sz, &gc, tmp, &f_step, &f_time, &g, &avs, &a_step, &a_time);
  
  if ( !mode )
  {
    if ( (a_step == 0) || (a_time <= 0.0) ) {
      Hostonly_ printf ("Error : restarted step[%d] or time[%e] is invalid\n", a_step, a_time);
      Exit(0);
    }
  }

  REAL_TYPE refv = (Dmode == DIMENSIONAL) ? RefVelocity : 1.0;
  REAL_TYPE scale = (mode == true) ? 1.0 : (REAL_TYPE)step_avr; // 瞬時値の時スケールは1.0、平均値の時は平均数
  REAL_TYPE u0[4];
  u0[0] = v00[0];
  u0[1] = v00[1];
  u0[2] = v00[2];
  u0[3] = v00[3];

  fb_shift_refv_in_(v, sz, &gc, u0, &scale, &refv, &flop);

  if ( mode )
  {
    Hostonly_ printf     ("\t[%s] has read :\tstep=%d  time=%e [%s]\n", tmp, f_step, f_time, (Dmode==DIMENSIONAL)?"sec.":"-");
    Hostonly_ fprintf(fp, "\t[%s] has read :\tstep=%d  time=%e [%s]\n", tmp, f_step, f_time, (Dmode==DIMENSIONAL)?"sec.":"-");
  }
  else
  {
    Hostonly_ printf     ("\t[%s] has read :\tstep=%d  time=%e \t:Averaged step=%d  time=%e [%s]\n",
                          tmp, f_step, f_time, a_step, a_time, (Dmode==DIMENSIONAL)?"sec.":"-");
    Hostonly_ fprintf(fp, "\t[%s] has read :\tstep=%d  time=%e \t:Averaged step=%d  time=%e [%s]\n", 
                      tmp, f_step, f_time, a_step, a_time, (Dmode==DIMENSIONAL)?"sec.":"-");
  }

  step = (unsigned)f_step;
  time = (double)f_time;
  step_avr = (unsigned)a_step;
  time_avr = (double)a_time;
}


// 温度のファイルをロードする
void FileIO::readTemperature(FILE* fp, 
                             const string fname,
                             int* sz, 
                             int gc, 
                             REAL_TYPE* t, 
                             unsigned& step, 
                             double& time, 
                             const int Dmode, 
                             const REAL_TYPE Base_tmp, 
                             const REAL_TYPE Diff_tmp,
                             const REAL_TYPE Kelvin,
                             double& flop, 
                             const int guide_out,
                             const bool mode,
                             unsigned& step_avr,
                             double& time_avr)
{
  
  if ( fname.empty() ) Exit(0);
  
  if ( fname.size() > FB_FILE_PATH_LENGTH )
  {
    Hostonly_ printf ("Error : Length of file path is greater than %d\n", FB_FILE_PATH_LENGTH);
    Exit(0);
  }
  
  char tmp[FB_FILE_PATH_LENGTH];
  memset(tmp, 0, sizeof(char)*FB_FILE_PATH_LENGTH);
  strcpy(tmp, fname.c_str());
  
  int g = guide_out;
  int avs = (mode == true) ? 1 : 0;
  int f_step, a_step;
  REAL_TYPE f_time, a_time;
  
  fb_read_sph_s_ (t, sz, &gc, tmp, &f_step, &f_time, &g, &avs, &a_step, &a_time);
  if ( !mode )
  {
    if ( (a_step == 0) || (a_time <= 0.0) )
    {
      Hostonly_ printf ("Error : restarted step[%d] or time[%e] is invalid\n", a_step, a_time);
      Exit(0);
    }
  }
  
  // 有次元ファイルの場合，無次元に変換する
  REAL_TYPE scale = (mode == true) ? 1.0 : (REAL_TYPE)step_avr; // 瞬時値の時スケールは1.0、平均値の時は平均数
  
  if ( Dmode == DIMENSIONAL )
  {
    U.tmp_array_D2ND(t, sz, gc, Base_tmp, Diff_tmp, Kelvin, scale, flop);
  }
  
  if ( mode ) {
    Hostonly_ printf     ("\t[%s] has read :\tstep=%d  time=%e [%s]\n", tmp, f_step, f_time, (Dmode==DIMENSIONAL)?"sec.":"-");
    Hostonly_ fprintf(fp, "\t[%s] has read :\tstep=%d  time=%e [%s]\n", tmp, f_step, f_time, (Dmode==DIMENSIONAL)?"sec.":"-");
  }
  else {
    Hostonly_ printf     ("\t[%s] has read :\tstep=%d  time=%e \t:Averaged step=%d  time=%e [%s]\n",
                          tmp, f_step, f_time, a_step, a_time, (Dmode==DIMENSIONAL)?"sec.":"-");
    Hostonly_ fprintf(fp, "\t[%s] has read :\tstep=%d  time=%e \t:Averaged step=%d  time=%e [%s]\n", 
                      tmp, f_step, f_time, a_step, a_time, (Dmode==DIMENSIONAL)?"sec.":"-");
  }

  step = (unsigned)f_step;
  time = (double)f_time;
  step_avr = (unsigned)a_step;
  time_avr = (double)a_time;
}


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
  REAL_TYPE o[3], p[3];
  o[0] = org[0];
  o[1] = org[1];
  o[2] = org[2];
  p[0] = pit[0];
  p[1] = pit[1];
  p[2] = pit[2];
  
  int avs = (mode == true) ? 1 : 0;
  int stp_a = (int)step_avr;
  REAL_TYPE tm_a = (REAL_TYPE)time_avr;
  int d_type = (sizeof(REAL_TYPE) == 4) ? 1 : 2;  // 1-float / 2-double
  
  fb_write_sph_s_ (s, sz, &gc, tmp, &stp, &tm, o, p, &d_type, &g, &avs, &stp_a, &tm_a);
  
}


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
  REAL_TYPE o[3], p[3];
  o[0] = org[0];
  o[1] = org[1];
  o[2] = org[2];
  p[0] = pit[0];
  p[1] = pit[1];
  p[2] = pit[2];
  
  int avs = (mode == true) ? 1 : 0;
  int stp_a = (int)step_avr;
  REAL_TYPE tm_a = (REAL_TYPE)time_avr;
  int d_type = (sizeof(REAL_TYPE) == 4) ? 1 : 2;  // 1-float / 2-double
  
  fb_write_sph_v_ (v, sz, &gc, tmp, &stp, &tm, o, p, &d_type, &g, &avs, &stp_a, &tm_a);

}

