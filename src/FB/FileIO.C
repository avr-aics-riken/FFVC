// #################################################################
//
// FFV : Frontflow / violet
//
// Copyright (c) All right reserved. 2012
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
//
// #################################################################

/** 
 * @file FileIO.C
 * @brief FlowBase FileIO class Header
 * @author kero
 */

#include "FileIO.h"


// CPMクラスポインタのコピー
void FileIO::importCPM(cpm_ParaManager* m_paraMngr)
{
  if ( !m_paraMngr ) Exit(0);
  paraMngr = m_paraMngr;
}


// ファイル出力時，発散値を計算する
void FileIO::cnv_Div(REAL_TYPE* dst, REAL_TYPE* src, const int* size, const int guide, const REAL_TYPE coef, REAL_TYPE& flop)
{
  if( !dst || !src || !size ) Exit(0);
  
  int sz[3];
  sz[0] = size[0];
  sz[1] = size[1];
  sz[2] = size[2];
  int gc = guide;
  REAL_TYPE cf = coef;
  
  fb_mulcpy_ (dst, src, sz, &gc, &cf, &flop);
}



// 全圧データについて，無次元から有次元単位に変換する
void FileIO::cnv_TP_ND2D(REAL_TYPE* dst, REAL_TYPE* src, const int* size, const int guide, 
                         const REAL_TYPE Ref_rho, const REAL_TYPE Ref_v, REAL_TYPE& flop)
{
  if( !dst || !src || !size ) Exit(0);
  
  int sz[3];
  sz[0] = size[0];
  sz[1] = size[1];
  sz[2] = size[2];
  int gc = guide;
  
  REAL_TYPE cf = Ref_rho*Ref_v*Ref_v;
  
  fb_mulcpy_ (dst, src, sz, &gc, &cf, &flop);
}


// sphファイルの書き出し（内部領域のみ）
void FileIO::writeRawSPH(const REAL_TYPE *vf, const int* size, const int gc, const REAL_TYPE* org, const REAL_TYPE* ddx, const int m_ModePrecision)
{
  int sz, dType, stp, svType;
  int ix, jx, kx, i, j, k;
  REAL_TYPE ox, oy, oz, dx, dy, dz, tm;
  long long szl[3], stpl;
  
  
  char sph_fname[512];
  
  if ( paraMngr->IsParallel() ) {
    sprintf( sph_fname, "field%010d.sph", paraMngr->GetMyRankID() );
  } else {
    sprintf( sph_fname, "field.sph" );
  }
  
  ofstream ofs(sph_fname, ios::out | ios::binary);
  if (!ofs) {
    cout << "\tCan't open " << sph_fname << " file" << endl;
    Exit(0);
  }
  
  ix = size[0]; //+2*gc;
  jx = size[1]; //+2*gc;
  kx = size[2]; //+2*gc;
  
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
  if ( sizeof(REAL_TYPE) == sizeof(double) ) {
    for (i=0; i<3; i++)   szl[i] = (long long)size[i];
  }
  
  REAL_TYPE *f = new REAL_TYPE[nx];
  
  size_t m, l;
  
  for (k=1; k<=kx; k++) {
    for (j=1; j<=jx; j++) {
      for (i=1; i<=ix; i++) {
        l = ix*jx*(k-1) + ix*(j-1) + i-1;
        m = FBUtility::getFindexS3D(size, gc, i, j, k);
        f[l] = (REAL_TYPE)vf[m];
      }
    }
  }
  
  // data property
  ( m_ModePrecision == FP_SINGLE ) ? dType=1 : dType=2;
  sz = sizeof(int)*2;
  ofs.write( (char*)&sz, sizeof(int) );
  ofs.write( (char*)&svType, sizeof(int) );
  ofs.write( (char*)&dType, sizeof(int) );
  ofs.write( (char*)&sz, sizeof(int) );
  
  // voxel size
  if (dType == 1) {
    sz = sizeof(int)*3;
    ofs.write( (char*)&sz, sizeof(int) );
    ofs.write( (char*)&ix, sizeof(int) );
    ofs.write( (char*)&jx, sizeof(int) );
    ofs.write( (char*)&kx, sizeof(int) );
    ofs.write( (char*)&sz, sizeof(int) );
  }
  else {
    sz = sizeof(long long)*3;
    ofs.write( (char*)&sz, sizeof(int) );
    ofs.write( (char*)&szl[0], sizeof(long long) );
    ofs.write( (char*)&szl[1], sizeof(long long) );
    ofs.write( (char*)&szl[2], sizeof(long long) );
    ofs.write( (char*)&sz, sizeof(int) );
  }
  
  // original point of domain
  if (dType == 1) {
    sz = sizeof(float)*3;
    ofs.write( (char*)&sz, sizeof(int) );
    ofs.write( (char*)&ox, sizeof(float) );
    ofs.write( (char*)&oy, sizeof(float) );
    ofs.write( (char*)&oz, sizeof(float) );
    ofs.write( (char*)&sz, sizeof(int) );
  }
  else {
    sz = sizeof(double)*3;
    ofs.write( (char*)&sz, sizeof(int) );
    ofs.write( (char*)&ox, sizeof(double) );
    ofs.write( (char*)&oy, sizeof(double) );
    ofs.write( (char*)&oz, sizeof(double) );
    ofs.write( (char*)&sz, sizeof(int) );
  }
  
  // pitch of voxel
  if (dType == 1) {
    sz = sizeof(float)*3;
    ofs.write( (char*)&sz, sizeof(int) );
    ofs.write( (char*)&dx, sizeof(float) );
    ofs.write( (char*)&dy, sizeof(float) );
    ofs.write( (char*)&dz, sizeof(float) );
    ofs.write( (char*)&sz, sizeof(int) );
  }
  else {
    sz = sizeof(double)*3;
    ofs.write( (char*)&sz, sizeof(int) );
    ofs.write( (char*)&dx, sizeof(double) );
    ofs.write( (char*)&dy, sizeof(double) );
    ofs.write( (char*)&dz, sizeof(double) );
    ofs.write( (char*)&sz, sizeof(int) );
  }
  
  // time stamp
  if (dType == 1) {
    stp = 0;
    tm = 0.0;
    sz = sizeof(int)+sizeof(float);
    ofs.write( (char*)&sz, sizeof(int) );
    ofs.write( (char*)&stp, sizeof(int) );
    ofs.write( (char*)&tm, sizeof(float) );
    ofs.write( (char*)&sz, sizeof(int) );
  }
  else {
    stpl =0;
    tm = 0.0;
    sz = sizeof(long long)+sizeof(double);
    ofs.write( (char*)&sz, sizeof(int) );
    ofs.write( (char*)&stpl, sizeof(long long) );
    ofs.write( (char*)&tm, sizeof(double) );
    ofs.write( (char*)&sz, sizeof(int) );
  }
  
  if (svType == kind_scalar) {
    sz = (m_ModePrecision == FP_SINGLE) ? nx * sizeof(float) : nx * sizeof(double);
    ofs.write( (char*)&sz, sizeof(int) );
    ofs.write( (char*)f,   sz );
    ofs.write( (char*)&sz, sizeof(int) );
  }
  
  ofs.close();
  
  if (f) { delete [] f; f=NULL; }
}




// 圧力のファイルをロードする
void FileIO::readPressure(FILE* fp,
                          const std::string fname,
                          int* size,
                          int gc,
                          REAL_TYPE* p,
                          int& step,
                          REAL_TYPE& time,
                          const int Dmode,
                          const REAL_TYPE BasePrs,
                          const REAL_TYPE RefDensity,
                          const REAL_TYPE RefVelocity,
                          REAL_TYPE& flop,
                          const int guide_out,
                          const bool mode,
                          int& step_avr,
                          REAL_TYPE& time_avr
                          )
{
  if ( fname.empty() ) Exit(0);
  
  if ( fname.size() > FB_FILE_PATH_LENGTH ) {
    Hostonly_ printf ("Error : Length of file path is greater than %d\n", FB_FILE_PATH_LENGTH);
    Exit(0);
  }
  

  char tmp[FB_FILE_PATH_LENGTH];
  memset(tmp, 0, sizeof(char)*FB_FILE_PATH_LENGTH);
  strcpy(tmp, fname.c_str());

  int g = guide_out;
  int avs = (mode == true) ? 1 : 0;
  
  fb_read_sph_s_ (p, size, &gc, tmp, &step, &time, &g, &avs, &step_avr, &time_avr);
  
  if ( !mode ) {
    if ( (step_avr == 0) || (time_avr <= 0.0) ) {
      Hostonly_ printf ("Error : restarted step[%d] or time[%e] is invalid\n", step_avr, time_avr);
      Exit(0);
    }
  }

  
  // 有次元ファイルの場合，無次元に変換する
  if ( Dmode == DIMENSIONAL ) {
    int d_length = (size[0]+2*gc) * (size[1]+2*gc) * (size[2]+2*gc);
    REAL_TYPE scale = (mode == true) ? 1.0 : (REAL_TYPE)step_avr; // 瞬時値の時スケールは1.0、平均値の時は平均数
    REAL_TYPE basep = BasePrs;
    REAL_TYPE ref_d = RefDensity;
    REAL_TYPE ref_v = RefVelocity;
  
    fb_prs_d2nd_(p, &d_length, &basep, &ref_d, &ref_v, &scale, &flop);
  }
  
  if ( mode ) {
    Hostonly_ printf     ("\t[%s] has read :\tstep=%d  time=%e [%s]\n", tmp, step, time, (Dmode==DIMENSIONAL)?"sec.":"-");
    Hostonly_ fprintf(fp, "\t[%s] has read :\tstep=%d  time=%e [%s]\n", tmp, step, time, (Dmode==DIMENSIONAL)?"sec.":"-");
  }
  else {
    Hostonly_ printf     ("\t[%s] has read :\tstep=%d  time=%e \t:Averaged step=%d  time=%e [%s]\n",
                          tmp, step, time, step_avr, time_avr, (Dmode==DIMENSIONAL)?"sec.":"-");
    Hostonly_ fprintf(fp, "\t[%s] has read :\tstep=%d  time=%e \t:Averaged step=%d  time=%e [%s]\n", 
                          tmp, step, time, step_avr, time_avr, (Dmode==DIMENSIONAL)?"sec.":"-");
  }

}



// 速度のファイルをロードする
void FileIO::readVelocity(FILE* fp, 
                          const std::string fname,
                          int* size, 
                          int gc, 
                          REAL_TYPE* v, 
                          int& step, 
                          REAL_TYPE& time, 
                          const REAL_TYPE *v00, 
                          const int Dmode, 
                          const REAL_TYPE RefVelocity, 
                          REAL_TYPE& flop, 
                          const int guide_out,
                          const bool mode,
                          int& step_avr,
                          REAL_TYPE& time_avr)
{
  if ( fname.empty() ) Exit(0);
  
  if ( fname.size() > FB_FILE_PATH_LENGTH ) {
    Hostonly_ printf ("Error : Length of file path is greater than %d\n", FB_FILE_PATH_LENGTH);
    Exit(0);
  }
  
  char tmp[FB_FILE_PATH_LENGTH];
  memset(tmp, 0, sizeof(char)*FB_FILE_PATH_LENGTH);
  strcpy(tmp, fname.c_str());

  int g = guide_out;
  int avs = (mode == true) ? 1 : 0;
  
  fb_read_sph_v_ (v, size, &gc, tmp, &step, &time, &g, &avs, &step_avr, &time_avr);
  
  if ( !mode ) {
    if ( (step_avr == 0) || (time_avr <= 0.0) ) {
      Hostonly_ printf ("Error : restarted step[%d] or time[%e] is invalid\n", step_avr, time_avr);
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

  fb_shift_refv_in_(v, size, &gc, u0, &scale, &refv, &flop);

  if ( mode ) {
    Hostonly_ printf     ("\t[%s] has read :\tstep=%d  time=%e [%s]\n", tmp, step, time, (Dmode==DIMENSIONAL)?"sec.":"-");
    Hostonly_ fprintf(fp, "\t[%s] has read :\tstep=%d  time=%e [%s]\n", tmp, step, time, (Dmode==DIMENSIONAL)?"sec.":"-");
  }
  else {
    Hostonly_ printf     ("\t[%s] has read :\tstep=%d  time=%e \t:Averaged step=%d  time=%e [%s]\n",
                          tmp, step, time, step_avr, time_avr, (Dmode==DIMENSIONAL)?"sec.":"-");
    Hostonly_ fprintf(fp, "\t[%s] has read :\tstep=%d  time=%e \t:Averaged step=%d  time=%e [%s]\n", 
                      tmp, step, time, step_avr, time_avr, (Dmode==DIMENSIONAL)?"sec.":"-");
  }

}


// 温度のファイルをロードする
void FileIO::readTemperature(FILE* fp, 
                             const std::string fname,
                             int* size, 
                             int gc, 
                             REAL_TYPE* t, 
                             int& step, 
                             REAL_TYPE& time, 
                             const int Dmode, 
                             const REAL_TYPE Base_tmp, 
                             const REAL_TYPE Diff_tmp, 
                             const REAL_TYPE Kelvin, 
                             REAL_TYPE& flop, 
                             const int guide_out,
                             const bool mode,
                             int& step_avr,
                             REAL_TYPE& time_avr)
{
  if ( fname.empty() ) Exit(0);
  
  if ( fname.size() > FB_FILE_PATH_LENGTH ) {
    Hostonly_ printf ("Error : Length of file path is greater than %d\n", FB_FILE_PATH_LENGTH);
    Exit(0);
  }
  
  char tmp[FB_FILE_PATH_LENGTH];
  memset(tmp, 0, sizeof(char)*FB_FILE_PATH_LENGTH);
  strcpy(tmp, fname.c_str());
  
  int g = guide_out;
  int avs = (mode == true) ? 1 : 0;
  
  fb_read_sph_s_ (t, size, &gc, tmp, &step, &time, &g, &avs, &step_avr, &time_avr);
  if ( !mode ) {
    if ( (step_avr == 0) || (time_avr <= 0.0) ) {
      Hostonly_ printf ("Error : restarted step[%d] or time[%e] is invalid\n", step_avr, time_avr);
      Exit(0);
    }
  }
  
  // 有次元ファイルの場合，無次元に変換する
  int d_length = (size[0]+2*gc) * (size[1]+2*gc) * (size[2]+2*gc);
  REAL_TYPE scale = (mode == true) ? 1.0 : (REAL_TYPE)step_avr; // 瞬時値の時スケールは1.0、平均値の時は平均数
  REAL_TYPE base_t = Base_tmp;
  REAL_TYPE diff_t = Diff_tmp;
  REAL_TYPE klv    = Kelvin;
  
  if ( Dmode == DIMENSIONAL ) {
    fb_tmp_d2nd_(t, &d_length, &base_t, &diff_t, &klv, &scale, &flop);
  }
  
  if ( mode ) {
    Hostonly_ printf     ("\t[%s] has read :\tstep=%d  time=%e [%s]\n", tmp, step, time, (Dmode==DIMENSIONAL)?"sec.":"-");
    Hostonly_ fprintf(fp, "\t[%s] has read :\tstep=%d  time=%e [%s]\n", tmp, step, time, (Dmode==DIMENSIONAL)?"sec.":"-");
  }
  else {
    Hostonly_ printf     ("\t[%s] has read :\tstep=%d  time=%e \t:Averaged step=%d  time=%e [%s]\n",
                          tmp, step, time, step_avr, time_avr, (Dmode==DIMENSIONAL)?"sec.":"-");
    Hostonly_ fprintf(fp, "\t[%s] has read :\tstep=%d  time=%e \t:Averaged step=%d  time=%e [%s]\n", 
                      tmp, step, time, step_avr, time_avr, (Dmode==DIMENSIONAL)?"sec.":"-");
  }

}


// スカラーファイルを出力する
void FileIO::writeScalar(const std::string fname, 
                         int* size, 
                         int gc,
                         REAL_TYPE* s, 
                         const int step, 
                         const REAL_TYPE time, 
                         const REAL_TYPE* org, 
                         const REAL_TYPE* pit, 
                         const int guide_out,
                         const bool mode,
                         const int step_avr,
                         const REAL_TYPE time_avr)
{
  if ( fname.empty() ) Exit(0);
  
  if ( fname.size() > FB_FILE_PATH_LENGTH ) {
    Hostonly_ printf ("Error : Length of file path is greater than %d\n", FB_FILE_PATH_LENGTH);
    Exit(0);
  }
  
  
  char tmp[FB_FILE_PATH_LENGTH];
  memset(tmp, 0, sizeof(char)*FB_FILE_PATH_LENGTH);
  strcpy(tmp, fname.c_str());

  int stp = step;
  REAL_TYPE tm = time;
  int g = guide_out;
  REAL_TYPE o[3], p[3];
  o[0] = org[0];
  o[1] = org[1];
  o[2] = org[2];
  p[0] = pit[0];
  p[1] = pit[1];
  p[2] = pit[2];
  
  int avs = (mode == true) ? 1 : 0;
  int stp_a = step_avr;
  REAL_TYPE tm_a = time_avr;
  int d_type = (sizeof(REAL_TYPE) == 4) ? 1 : 2;  // 1-float / 2-double
  
  fb_write_sph_s_ (s, size, &gc, tmp, &stp, &tm, o, p, &d_type, &g, &avs, &stp_a, &tm_a);
  
}


// ベクトルファイルを出力する
void FileIO::writeVector(const std::string fname, 
                         int* size, 
                         int gc, 
                         REAL_TYPE* v, 
                         const int step, 
                         const REAL_TYPE time, 
                         const REAL_TYPE* org, 
                         const REAL_TYPE* pit, 
                         const int guide_out,
                         const bool mode,
                         const int step_avr,
                         const REAL_TYPE time_avr)
{
  if ( fname.empty() ) Exit(0);
  
  if ( fname.size() > FB_FILE_PATH_LENGTH ) {
    Hostonly_ printf ("Error : Length of file path is greater than %d\n", FB_FILE_PATH_LENGTH);
    Exit(0);
  }
  
  char tmp[FB_FILE_PATH_LENGTH];
  memset(tmp, 0, sizeof(char)*FB_FILE_PATH_LENGTH);
  strcpy(tmp, fname.c_str());
  
  int stp = step;
  REAL_TYPE tm = time;
  int g = guide_out;
  REAL_TYPE o[3], p[3];
  o[0] = org[0];
  o[1] = org[1];
  o[2] = org[2];
  p[0] = pit[0];
  p[1] = pit[1];
  p[2] = pit[2];
  
  int avs = (mode == true) ? 1 : 0;
  int stp_a = step_avr;
  REAL_TYPE tm_a = time_avr;
  int d_type = (sizeof(REAL_TYPE) == 4) ? 1 : 2;  // 1-float / 2-double
  
  fb_write_sph_v_ (v, size, &gc, tmp, &stp, &tm, o, p, &d_type, &g, &avs, &stp_a, &tm_a);

}

