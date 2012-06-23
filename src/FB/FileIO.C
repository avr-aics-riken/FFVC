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


// 作業用ポインタのコピー
void FileIO::setPartitionManager(cpm_ParaManager* m_paraMngr)
{
  if ( !m_paraMngr ) Exit(0);
  paraMngr = m_paraMngr;
}


/**
 @fn void FileIO::cnv_Div(SklScalar3D<REAL_TYPE>* dst, const SklScalar3D<REAL_TYPE>* src, const REAL_TYPE coef, REAL_TYPE& flop)
 @brief ファイル出力時，発散値を計算する
 @param dst 単位変換のデータクラス
 @param src 単位変換前のデータクラス
 @param coef 係数
 @param flop 浮動小数演算数
 @see SklUtil::cpyS3D()
 */
void FileIO::cnv_Div(SklScalar3D<REAL_TYPE>* dst, const SklScalar3D<REAL_TYPE>* src, const REAL_TYPE coef, REAL_TYPE& flop)
{
  if( !dst || !src ) Exit(0);
  
  const unsigned* dst_sz = dst->GetSize();
  const unsigned* src_sz = src->GetSize();
  if(   (src_sz[0] != dst_sz[0])
     || (src_sz[1] != dst_sz[1])
     || (src_sz[2] != dst_sz[2]) ) Exit(0);
  
  REAL_TYPE* dst_data = dst->GetData();
  const REAL_TYPE* src_data = src->GetData();
  if( !dst_data || !src_data ) Exit(0);
  
  dst_sz = dst->_GetSize();
  src_sz = src->_GetSize();
  unsigned dst_gc = dst->GetVCellSize();
  unsigned src_gc = src->GetVCellSize();
  unsigned sta, ix, jx, kx;
  int diff;
  register unsigned i, j, k;
  
  unsigned long dst_lsz[3];
  dst_lsz[0] = dst_sz[0];
  dst_lsz[1] = dst_sz[1];
  dst_lsz[2] = dst_sz[2];
  
  // compare whole size
  if(   (src_sz[0] == dst_sz[0])
     && (src_sz[1] == dst_sz[1])
     && (src_sz[2] == dst_sz[2]) ) {
    unsigned long idx;
    
    for(k=0; k<dst_sz[2]; k++){
      for(j=0; j<dst_sz[1]; j++){
        for(i=0; i<dst_sz[0]; i++){
          idx = dst_lsz[0]*dst_lsz[1]*k  + dst_lsz[0]*j  + i;
          dst_data[idx] = src_data[idx]*coef;
        }
      }
    }
    flop += (REAL_TYPE)(dst_sz[0]*dst_sz[1]*dst_sz[2]);
  }
  else {
    Exit(0);
  }
}

/**
 @fn void FileIO::cnv_TP_ND2D(SklScalar3D<REAL_TYPE>* dst, const SklScalar3D<REAL_TYPE>* src, const REAL_TYPE Ref_rho, const REAL_TYPE Ref_v, REAL_TYPE& flop)
 @brief 全圧データについて，無次元から有次元単位に変換する
 @param dst 単位変換のデータクラス
 @param src 単位変換前のデータクラス
 @param Ref_rho 代表密度(kg/m^3)
 @param Ref_v 代表速度(m/s)
 @param flop 浮動小数演算数
 @see SklUtil::cpyS3D()
 */
void FileIO::cnv_TP_ND2D(SklScalar3D<REAL_TYPE>* dst, const SklScalar3D<REAL_TYPE>* src, const REAL_TYPE Ref_rho, const REAL_TYPE Ref_v, REAL_TYPE& flop)
{
  if( !dst || !src ) Exit(0);
  
  const unsigned* dst_sz = dst->GetSize();
  const unsigned* src_sz = src->GetSize();
  if(   (src_sz[0] != dst_sz[0])
     || (src_sz[1] != dst_sz[1])
     || (src_sz[2] != dst_sz[2]) ) Exit(0);
  
  REAL_TYPE* dst_data = dst->GetData();
  const REAL_TYPE* src_data = src->GetData();
  if( !dst_data || !src_data ) Exit(0);
  
  dst_sz = dst->_GetSize();
  src_sz = src->_GetSize();
  unsigned dst_gc = dst->GetVCellSize();
  unsigned src_gc = src->GetVCellSize();
  unsigned sta, ix, jx, kx;
  int diff;
  register unsigned i, j, k;
  
  REAL_TYPE c = Ref_rho*Ref_v*Ref_v;
  flop += (REAL_TYPE)(dst_sz[0]*dst_sz[1]*dst_sz[2]*1);
  
  unsigned long src_idx, dst_idx, idx, src_lsz[3], dst_lsz[3];
  src_lsz[0] = src_sz[0];
  src_lsz[1] = src_sz[1];
  src_lsz[2] = src_sz[2];
  dst_lsz[0] = dst_sz[0];
  dst_lsz[1] = dst_sz[1];
  dst_lsz[2] = dst_sz[2];
  
  // compare whole size
  if(   (src_sz[0] == dst_sz[0])
     && (src_sz[1] == dst_sz[1])
     && (src_sz[2] == dst_sz[2]) ) {
    unsigned idx;
    
    for(k=0; k<dst_sz[2]; k++){
      for(j=0; j<dst_sz[1]; j++){
        for(i=0; i<dst_sz[0]; i++){
          idx = dst_lsz[0]*dst_lsz[1]*k  + dst_lsz[0]*j  + i;
          dst_data[idx] = src_data[idx]*c;
        }
      }
    }
  }
  else {
    CalcIndex(dst_sz[0], dst_sz[1], dst_sz[2], dst_gc, src_gc, ix, jx, kx, diff, sta);
    
    for(k=sta; k<kx; k++){
      unsigned kk = k+diff;
      for(j=sta; j<jx; j++){
        unsigned jj = j+diff;
        for(i=sta; i<ix; i++){
          src_idx = src_lsz[0]*src_lsz[1]*kk + src_lsz[0]*jj + i+diff;
          dst_idx = dst_lsz[0]*dst_lsz[1]*k  + dst_lsz[0]*j  + i;
          dst_data[dst_idx] = src_data[src_idx]*c;
        }
      }
    }
  }
}


/**
 @fn void FileIO::writeRawSPH(const REAL_TYPE *vf, const unsigned* size, const unsigned gc, const REAL_TYPE* org, const REAL_TYPE* ddx, const unsigned m_ModePrecision)
 @brief sphファイルの書き出し（内部領域のみ）
 @param vf スカラデータ
 @param size 配列サイズ
 @param gc ガイドセル
 @param org 基点
 @param ddx ピッチ
 @param m_ModePrecision 浮動小数点の精度
 @note 標記上，long 対応になっているが，ファイルフォーマットとの対応を確認のこと
 */
void FileIO::writeRawSPH(const REAL_TYPE *vf, const unsigned* size, const unsigned gc, const REAL_TYPE* org, const REAL_TYPE* ddx, const unsigned m_ModePrecision)
{
  //SklParaManager* para_mng = ParaCmpo->GetParaManager();
  
  int sz, dType, stp, svType;
  int ix, jx, kx, i, j, k;
  unsigned long l, nx;
  REAL_TYPE ox, oy, oz, dx, dy, dz, tm;
  long long szl[3], stpl;
  unsigned m;
  
  char sph_fname[512];
  
  if ( pn.numProc > 1 ) {
    sprintf( sph_fname, "field%010d.sph", pn.myrank );
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
  nx = ix * jx * kx;
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
  ( m_ModePrecision == SPH_SINGLE ) ? dType=1 : dType=2;
  sz = sizeof(unsigned)*2;
  ofs.write( (char*)&sz, sizeof(int) );
  ofs.write( (char*)&svType, sizeof(int) );
  ofs.write( (char*)&dType, sizeof(int) );
  ofs.write( (char*)&sz, sizeof(int) );
  
  // voxel size
  if (dType == 1) {
    sz = sizeof(unsigned)*3;
    ofs.write( (char*)&sz, sizeof(int) );
    ofs.write( (char*)&ix, sizeof(unsigned) );
    ofs.write( (char*)&jx, sizeof(unsigned) );
    ofs.write( (char*)&kx, sizeof(unsigned) );
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
    sz = (m_ModePrecision == SPH_SINGLE) ? nx * sizeof(float) : nx * sizeof(double);
    ofs.write( (char*)&sz, sizeof(int) );
    ofs.write( (char*)f,   sz );
    ofs.write( (char*)&sz, sizeof(int) );
  }
  
  ofs.close();
  
  if (f) { delete [] f; f=NULL; }
}



//@fn void FileIO::readPressure()
//@brief 圧力ファイルをロードする
void FileIO::readPressure(FILE* fp,                    /// @param fp ファイルポインタ（ファイル出力）
                          const std::string fname,     /// @param fname ファイル名
                          const unsigned* size,        /// @param size サイズ
                          const unsigned gc,           /// @param guide ガイドセルサイズ
                          REAL_TYPE* p,                /// @param p 圧力データ
                          int& step,                   /// @param step[out] ステップ
                          REAL_TYPE& time,             /// @param time[out] 時刻
                          const unsigned Dmode,        /// @param Dmode 次元（無次元-0 / 有次元-1）
                          const REAL_TYPE BasePrs,     /// @param BasePrs 基準圧力
                          const REAL_TYPE RefDensity,  /// @param RefDensity　代表密度
                          const REAL_TYPE RefVelocity, /// @param RefVelocity 代表速度
                          REAL_TYPE& flop,             /// @param flop
                          const int guide_out,         /// @param guide_out 出力ガイドセル数
                          const bool mode,             /// @param mode 平均値出力指示（瞬時値のときtrue，平均値のときfalse）
                          int& step_avr,               /// @param step_avr 平均操作したステップ数
                          REAL_TYPE& time_avr          /// @param time_avr 平均操作した時間
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
  
  fb_read_sph_s_ (p, (int*)size, (int*)&gc, tmp, &step, &time, &g, &avs, &step_avr, &time_avr);
  
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

/**
 @fn void FileIO::readVelocity(FILE* fp, 
 const std::string fname,
 const unsigned* size, 
 const unsigned gc, 
 REAL_TYPE* v, 
 int& step, 
 REAL_TYPE& time, 
 const REAL_TYPE *v00, 
 const unsigned Dmode, 
 const REAL_TYPE RefVelocity, 
 REAL_TYPE& flop, 
 const bool mode)
 @brief 速度をロードする
 @param fp ファイルポインタ（ファイル出力）
 @param fname ファイル名
 @param size サイズ
 @param guide ガイドセルサイズ
 @param block ブロック数（粗い格子のロードのときに指定、通常は1）
 @param v  結果を保持するデータ
 @param step ステップ
 @param time 時刻
 @param v00[4]
 @param Dmode 次元（無次元-0 / 有次元-1）
 @param RefVelocity 代表速度
 @param flop
 @param guide_out 出力ガイドセル数
 @param mode 平均値出力指示（瞬時値のときtrue，平均値のときfalse）
 @param step_avr 平均操作したステップ数
 @param time_avr 平均操作した時間
 */
void FileIO::readVelocity(FILE* fp, 
                          const std::string fname,
                          const unsigned* size, 
                          const unsigned gc, 
                          REAL_TYPE* v, 
                          int& step, 
                          REAL_TYPE& time, 
                          const REAL_TYPE *v00, 
                          const unsigned Dmode, 
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
  
  fb_read_sph_v_ (v, (int*)size, (int*)&gc, tmp, &step, &time, &g, &avs, &step_avr, &time_avr);
  
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

  fb_shift_refv_in_(v, (int*)size, (int*)&gc, u0, &scale, &refv, &flop);

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

/**
 @fn void FileIO::readTemperature(FILE* fp, 
 const std::string fname,
 const unsigned* size, 
 const unsigned gc, 
 REAL_TYPE* t, 
 int& step, 
 REAL_TYPE& time, 
 const unsigned Dmode, 
 const REAL_TYPE Base_tmp, 
 const REAL_TYPE Diff_tmp, 
 const REAL_TYPE Kelvin, 
 REAL_TYPE& flop, 
 const bool mode)
 @brief 温度をロードする
 @param fp ファイルポインタ（ファイル出力）
 @param fname ファイル名
 @param size グローバルなサイズ
 @param gc ガイドセルサイズ
 @param block ブロック数（粗い格子のロードのときに指定、通常は1）
 @param t  結果を保持するデータクラス
 @param step[out] ステップ
 @param time[out] 時刻
 @param Dmode 次元（無次元-0 / 有次元-1）
 @param Base_tmp 基準温度
 @param Diff_tmp　代表温度差
 @param Kelvin 定数
 @param flop
 @param guide_out 出力ガイドセル数
 @param mode 平均値出力指示（瞬時値のときtrue，平均値のときfalse）
 @param step_avr 平均操作したステップ数
 @param time_avr 平均操作した時間
 */
void FileIO::readTemperature(FILE* fp, 
                             const std::string fname,
                             const unsigned* size, 
                             const unsigned gc, 
                             REAL_TYPE* t, 
                             int& step, 
                             REAL_TYPE& time, 
                             const unsigned Dmode, 
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
  
  fb_read_sph_s_ (t, (int*)size, (int*)&gc, tmp, &step, &time, &g, &avs, &step_avr, &time_avr);
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

/**
 @fn void FileIO::writeScalar(const std::string fname, 
 const unsigned* size, 
 const unsigned gc,
 REAL_TYPE* s, 
 const int step, 
 const REAL_TYPE time, 
 const REAL_TYPE* org, 
 const REAL_TYPE* pit, 
 const int guide_out)
 @brief スカラー場を出力する
 @param fname ファイル名
 @param size
 @param gc
 @param s スカラー場
 @param step ステップ
 @param time 時刻
 @param org
 @param pit
 @param guide_out ガイドセル数
 @param mode 平均値出力指示（瞬時値のときtrue，平均値のときfalse）
 @param step_avr 平均操作したステップ数
 @param time_avr 平均操作した時間
 */
void FileIO::writeScalar(const std::string fname, 
                         const unsigned* size, 
                         const unsigned gc,
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
  
  fb_write_sph_s_ (s, (int*)size, (int*)&gc, tmp, &stp, &tm, o, p, &d_type, &g, &avs, &stp_a, &tm_a);
  
}

/**
 @fn void FileIO::writeVector(const std::string fname, 
 const unsigned* size, 
 const unsigned gc, 
 REAL_TYPE* v, 
 const int step, 
 const REAL_TYPE time, 
 const REAL_TYPE* org, 
 const REAL_TYPE* pit, 
 const int guide_out)
 @brief ベクトル場を出力する
 @param fname ファイル名
 @param size
 @param gc
 @param v ベクトル場
 @param step ステップ
 @param time 時刻
 @param org
 @param pit
 @param guide_out ガイドセル数
 @param mode 平均値出力指示（瞬時値のときtrue，平均値のときfalse）
 @param step_avr 平均操作したステップ数
 @param time_avr 平均操作した時間
 */
void FileIO::writeVector(const std::string fname, 
                         const unsigned* size, 
                         const unsigned gc, 
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
  
  fb_write_sph_v_ (v, (int*)size, (int*)&gc, tmp, &stp, &tm, o, p, &d_type, &g, &avs, &stp_a, &tm_a);

}

