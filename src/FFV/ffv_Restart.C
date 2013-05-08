// #################################################################
//
// FFV : Frontflow / violet Cartesian
//
// Copyright (c) 2012-2013  All right reserved.
//
// Institute of Industrial Science, University of Tokyo, Japan. 
//
// #################################################################

/** 
 * @file   ffv_Restart.C
 * @brief  FFV Class
 * @author kero
 */

#include "ffv.h"



// #################################################################
// ファイルのオープンチェック
bool FFV::checkFile(std::string fname)
{
  ifstream ifs(fname.c_str(), ios::in | ios::binary);
  if (!ifs)
  {
    return false;
  }
  ifs.close();
  
  return true;
}



// #################################################################
/**
 * @brief 2倍密格子の領域開始インデクス番号から、その領域が属する粗格子計算結果ファイル名と、その計算結果ファイルの開始インデクス番号を取得する
 * @param [in]  i                 密格子　開始インデクスi
 * @param [in]  j                        同j
 * @param [in]  k                        同k
 * @param [in]  coarse_dfi_fname  粗格子のdfiファイル名（どの変数のものでも良い）
 * @param [in]  coarse_prefix     粗格子計算結果ファイルプリフィクス e.g. "prs_16"
 * @param [in]  m_step            探索するステップ数
 * @param [out] coarse_sph_fname  ijk位置の結果を含む粗格子計算結果ファイル名
 * @param [out] c_size            粗格子の分割数
 * @param [out] coarse            粗格子　開始インデクス
 * @param [out] block             含まれるブロック数
 * return エラーコード
 */
bool FFV::getCoarseResult (int i, int j, int k,
                           std::string& coarse_dfi_fname,
                           std::string& coarse_prefix,
                           const int m_step,
                           std::string& coarse_sph_fname,
                           int* c_size,
                           int* coarse,
                           int* block
                           )
{
	// 密格子のi,j,kを粗格子のi0,j0,k0に変換
	int i0 = (i+1)/2; 
  int j0 = (j+1)/2; 
  int k0 = (k+1)/2;
  
  // ステップ数の文字列を生成
  char tmp[10+1]; // 10 digit + \0
  memset(tmp, 0, sizeof(char)*11);
  sprintf(tmp, "%010d", m_step);
  std::string step(tmp);
  
  
	// dfiファイルローダをインスタンス
  TPControl tp_dfi;
  
  tp_dfi.getTPinstance();
  
  // 入力ファイルをオープン
  int ierror = tp_dfi.readTPfile(coarse_dfi_fname);
  
  std::string str;
  std::string label;
  std::string label_base;
  std::string label_leaf;
  int ibuf;
  int iv[3];
  
  std::string fmt(C.file_fmt_ext);

	
  std::string buf;
	int rank = -1;
	int hi, hj, hk, ti, tj, tk;
  
  // サブドメインの分割数（粗格子）
  int Ci, Cj, Ck;
  

  // ノード情報を探索
  label_base = "/DistributedFileInfo/NodeInfo";
  
  if ( !tp_dfi.chkNode(label_base) )
  {
    Hostonly_ stamped_printf("\tParsing error : Missing '%s'\n", label_base.c_str());
    Exit(0);
  }
  
  // 粗格子 (i0,j0,k0) が含まれるランクを探す
  int n = 0;
  while ( n < numProc ) {
    
    if ( !tp_dfi.GetNodeStr(label_base, n+1, &str) ) // Node[@], @=0,...
    {
      Hostonly_ printf("\tParsing error : Missing 'Node'\n");
      Exit(0);
    }
    
    if ( !strcasecmp(str.substr(0,4).c_str(), "Node") )
    {
      label_leaf = label_base + "/" + str;
      
      // ランク
      label = label_leaf + "/RankID";
      
      if ( !tp_dfi.GetValue(label, &ibuf) )
      {
        Hostonly_ printf("\tParsing error : Invalid integer value for '%s'\n", label.c_str());
        Exit(0);
      }
      rank = ibuf;
      
      // VoxelSize
      label = label_leaf + "/VoxelSize";
      if ( !(tp_dfi.GetVector(label, iv, 3)) )
      {
        Hostonly_ printf("\tParsing error : Invalid integer value for '%s'\n", label.c_str());
        Exit(0);
      }
      Ci = iv[0];
      Cj = iv[1];
      Ck = iv[2];
      
      // HeadIndex
      label = label_leaf + "/HeadIndex";
      if ( !(tp_dfi.GetVector(label, iv, 3)) )
      {
        Hostonly_ printf("\tParsing error : Invalid integer value for '%s'\n", label.c_str());
        Exit(0);
      }
      hi = iv[0];
      hj = iv[1];
      hk = iv[2];
      
      // TailIndex
      label = label_leaf + "/TailIndex";
      if ( !(tp_dfi.GetVector(label, iv, 3)) )
      {
        Hostonly_ printf("\tParsing error : Invalid integer value for '%s'\n", label.c_str());
        Exit(0);
      }
      ti = iv[0];
      tj = iv[1];
      tk = iv[2];
      
      //printf("rank=%d size=(%d,%d,%d) head=(%d,%d,%d) tail=(%d,%d,%d) : (%d,%d,%d)\n",
      //       rank, Ci, Cj, Ck, hi, hj, hk, ti, tj, tk, i0, j0, k0);
      
      if ( i0>=hi && i0<=ti && j0>=hj && j0<=tj && k0>=hk && k0<=tk ) break; // found!
    }
    
    n++;
  }
  
  // 見つけられなかった
  if ( rank == -1 ) return false;
  
  
  // ランク番号のファイル名を生成する
  std::string target = DFI.GenerateFileName(coarse_prefix, fmt, m_step, rank, (bool)C.FIO.IOmode);
  if ( target.empty() ) return false;
  
  
  
  // 各方向に含まれるブロック数（dx_C/dx_F = 2）
  int bi = Ci * 2 / size[0];
  int bj = Cj * 2 / size[1];
  int bk = Ck * 2 / size[2];
  
  // 粗格子の読み込み開始のローカルインデクス
  int bh_i = i0 - hi + 1;
  int bh_j = j0 - hj + 1;
  int bh_k = k0 - hk + 1;
  
  
  // return value
  coarse_sph_fname = target;
  
  c_size[0] = Ci;
  c_size[1] = Cj;
  c_size[2] = Ck;
  
	coarse[0] = bh_i;
	coarse[1] = bh_j;
	coarse[2] = bh_k;
  
  block[0] = bi;
  block[1] = bj;
  block[2] = bk;
  
  //
  //printf("rk=%4d : fine(%4d %4d %4d) : coarse(%4d %4d %4d) : head(%4d %4d %4d) : cblk(%4d %4d %4d): block(%4d %4d %4d) : %s\n", 
  //       pn.myrank,  i,j,k,  i0,j0,k0,  bh_i,bh_j,bh_k,  Ci,Cj,Ck,  Fi,Fj,Fk,  coarse_sph_fname.c_str());
  //
  
	return true;
}



// #################################################################
/**
 * @brief 2倍密格子の領域開始インデクス番号から、その領域が属する粗格子計算結果ファイル名と、その計算結果ファイルの開始インデクス番号を取得する
 * @param [in]  i                 密格子　開始インデクスi
 * @param [in]  j                        同j
 * @param [in]  k                        同k
 * @param [in]  coarse_dfi_fname  粗格子のdfiファイル名（どの変数のものでも良い）
 * @param [in]  coarse_prefix     粗格子計算結果ファイルプリフィクス e.g. "prs_16"
 * @param [in]  m_step            探索するステップ数
 * @param [out] coarse_sph_fname  ijk位置の結果を含む粗格子計算結果ファイル名
 * @param [out] c_size            粗格子の分割数
 * @param [out] coarse            粗格子　開始インデクス
 * @param [out] block             含まれるブロック数
 * return エラーコード
 */
bool FFV::getCoarseResult2(int i, int j, int k,
                           std::string& coarse_dfi_fname,
                           std::string& coarse_prefix,
                           const int m_step,
                           std::string& coarse_sph_fname,
                           int* c_size,
                           int* coarse,
                           int* block
                           )
{
	// 密格子のi,j,kを粗格子のi0,j0,k0に変換
	int i0 = (i+1)/2;
  int j0 = (j+1)/2;
  int k0 = (k+1)/2;
  
  // ステップ数の文字列を生成
  char tmp[10+1]; // 10 digit + \0
  memset(tmp, 0, sizeof(char)*11);
  sprintf(tmp, "%010d", m_step);
  std::string step(tmp);
  
  
	// dfiファイルローダをインスタンス
  TPControl tp_dfi;
  
  tp_dfi.getTPinstance();
  
  // 入力ファイルをオープン
  int ierror = tp_dfi.readTPfile(coarse_dfi_fname);
  
  std::string str;
  std::string label;
  std::string label_base;
  std::string label_leaf;
  int ibuf;
  int iv[3];
  
  std::string fmt(C.file_fmt_ext);
  
  
  std::string buf;
	int rank = -1;
	int hi, hj, hk, ti, tj, tk;
  
  // サブドメインの分割数（粗格子）
  int Ci, Cj, Ck;
  
  
  // ノード情報を探索
  label_base = "/DistributedFileInfo/NodeInfo";
  
  if ( !tp_dfi.chkNode(label_base) )
  {
    Hostonly_ stamped_printf("\tParsing error : Missing '%s'\n", label_base.c_str());
    Exit(0);
  }
  
  const int* m_tail = paraMngr->GetVoxelTailIndex();
  int tail[3];
  tail[0] = m_tail[0] + 1;
  tail[1] = m_tail[1] + 1;
  tail[2] = m_tail[2] + 1;
  
	rank = myRank;
	Ci = size[0]/2;
	Cj = size[1]/2;
	Ck = size[2]/2;
	hi = (head[0]+1)/2;
	hj = (head[1]+1)/2;
	hk = (head[2]+1)/2;
	ti = (tail[0]+1)/2;
	tj = (tail[1]+1)/2;
	tk = (tail[2]+1)/2;
  
  
  // 見つけられなかった
  if ( rank == -1 ) return false;
  
  
  // ランク番号のファイル名を生成する
  std::string target = DFI.GenerateFileName(coarse_prefix, fmt, m_step, rank, (bool)C.FIO.IOmode);
  if ( target.empty() ) return false;
  
  
  
  // 各方向に含まれるブロック数（dx_C/dx_F = 2）
  int bi = Ci * 2 / size[0];
  int bj = Cj * 2 / size[1];
  int bk = Ck * 2 / size[2];
  
  // 粗格子の読み込み開始のローカルインデクス
  int bh_i = i0 - hi + 1;
  int bh_j = j0 - hj + 1;
  int bh_k = k0 - hk + 1;
  
  
  // return value
  coarse_sph_fname = target;
  
  c_size[0] = Ci;
  c_size[1] = Cj;
  c_size[2] = Ck;
  
	coarse[0] = bh_i;
	coarse[1] = bh_j;
	coarse[2] = bh_k;
  
  block[0] = bi;
  block[1] = bj;
  block[2] = bk;
  
  //
  //printf("rk=%4d : fine(%4d %4d %4d) : coarse(%4d %4d %4d) : head(%4d %4d %4d) : cblk(%4d %4d %4d): block(%4d %4d %4d) : %s\n",
  //       pn.myrank,  i,j,k,  i0,j0,k0,  bh_i,bh_j,bh_k,  Ci,Cj,Ck,  Fi,Fj,Fk,  coarse_sph_fname.c_str());
  //
  
	return true;
}



// 粗格子から密格子へ内挿
void FFV::Interpolation_from_coarse_initial(const int* m_st, const int* m_bk)
{
  int st[3], bk[3];
  st[0] = m_st[0];
  st[1] = m_st[1];
  st[2] = m_st[2];
  bk[0] = m_bk[0];
  bk[1] = m_bk[1];
  bk[2] = m_bk[2];
  
  fb_interp_coarse0_s_(d_p, size, &guide, d_r_p, st, bk);
  fb_interp_coarse0_v_(d_v, size, &guide, d_r_v, st, bk);
  
  if ( C.isHeatProblem() )
  {
    fb_interp_coarse0_s_(d_t, size, &guide, d_r_t, st, bk);
  }
  
}



// #################################################################
/**
 * @brief リスタートプロセス
 * @param [in]     fp     ファイルポインタ
 */
void FFV::Restart(FILE* fp)
{
  double flop_task;
  
  TIMING_start(tm_restart);
  
  if ( C.Start == initial_start) // 初期スタートのステップ，時間を設定する
  {
    Session_StartStep = CurrentStep = 0;
    Session_StartTime = CurrentTime = 0.0;
    
    // V00の値のセット．モードがONの場合はV00[0]=1.0に設定，そうでなければtmに応じた値
    if ( C.CheckParam == ON ) RF.setV00(CurrentTime, true);
    else                      RF.setV00(CurrentTime);
    
    double g[4];
    RF.copyV00(g);
    for (int i=0; i<4; i++) v00[i]=(REAL_TYPE)g[i];
    
  }
  
  else if ( C.Start == restart) // 同一解像度のリスタート
  {
    Hostonly_ fprintf(stdout, "\t>> Restart from Previous Calculated Results\n\n");
    Hostonly_ fprintf(fp, "\t>> Restart from Previous Calculated Results\n\n");
    
    flop_task = 0.0;
    Restart_std(fp, flop_task); 
  }
  
  else if ( C.Start == restart_refinement) // 粗い格子からのリスタート
  {
    Hostonly_ fprintf(stdout, "\t>> Restart from Previous Results on Coarse Mesh\n\n");
    Hostonly_ fprintf(fp, "\t>> Restart from Previous Results on Coarse Mesh\n\n");
    
    
    // 粗い格子のファイルをロードし、内挿処理を行う
    flop_task = 0.0;
    Restart_std(fp, flop_task);
    
    Hostonly_ fprintf(stdout,"\n");
    Hostonly_ fprintf(fp,"\n");
  }
  
  else if ( C.Start == restart_different_proc) // 異なる並列数でのリスタート
  {
    Hostonly_ fprintf(stdout, "\t>> Restart from Previous Calculated Results That Nproc Differ from\n\n");
    Hostonly_ fprintf(fp, "\t>> Restart from Previous Calculated Results That Nproc Differ from\n\n");
    
    flop_task = 0.0;
    Restart_std(fp, flop_task);
  }
  
  TIMING_stop(tm_restart);
}



// #################################################################
// リスタートの最大値と最小値の表示
void FFV::Restart_display_minmax(FILE* fp, double& flop)
{
  if ( (C.Start == restart) || (C.Start == restart_refinement) || (C.Start == restart_different_proc) )
  {
    Hostonly_ fprintf(stdout, "\tNon-dimensional value\n");
    Hostonly_ fprintf(fp, "\tNon-dimensional value\n");
    REAL_TYPE f_min, f_max, min_tmp, max_tmp;
    
    // Velocity
    fb_minmax_v_ (&f_min, &f_max, size, &guide, v00, d_v, &flop); // allreduceすること
    
    if ( numProc > 1 ) 
    {
      min_tmp = f_min;
      if( paraMngr->Allreduce(&min_tmp, &f_min, 1, MPI_MIN) != CPM_SUCCESS ) Exit(0);
      
      max_tmp = f_max;
      if( paraMngr->Allreduce(&max_tmp, &f_max, 1, MPI_MAX) != CPM_SUCCESS ) Exit(0);
    }
    
    Hostonly_ fprintf(stdout, "\t\tV : min=%13.6e / max=%13.6e\n", f_min, f_max);
    Hostonly_ fprintf(fp, "\t\tV : min=%13.6e / max=%13.6e\n", f_min, f_max);
    
    
    // Pressure
    fb_minmax_s_ (&f_min, &f_max, size, &guide, d_p, &flop);
    
    if ( numProc > 1 ) 
    {
      min_tmp = f_min;
      if( paraMngr->Allreduce(&min_tmp, &f_min, 1, MPI_MIN) != CPM_SUCCESS ) Exit(0);
      
      max_tmp = f_max;
      if( paraMngr->Allreduce(&max_tmp, &f_max, 1, MPI_MAX) != CPM_SUCCESS ) Exit(0);
    }
    
    Hostonly_ fprintf(stdout, "\t\tP : min=%13.6e / max=%13.6e\n", f_min, f_max);
    Hostonly_ fprintf(fp, "\t\tP : min=%13.6e / max=%13.6e\n", f_min, f_max);
    
    // temperature
    if ( C.isHeatProblem() ) 
    {
      fb_minmax_s_ (&f_min, &f_max, size, &guide, d_t, &flop);
      
      if ( numProc > 1 ) 
      {
        min_tmp = f_min;
        if( paraMngr->Allreduce(&min_tmp, &f_min, 1, MPI_MIN) != CPM_SUCCESS ) Exit(0);
        
        max_tmp = f_max;
        if( paraMngr->Allreduce(&max_tmp, &f_max, 1, MPI_MAX) != CPM_SUCCESS ) Exit(0);
      }
      
      Hostonly_ fprintf(stdout, "\t\tT : min=%13.6e / max=%13.6e\n", f_min, f_max);
      Hostonly_ fprintf(fp, "\t\tT : min=%13.6e / max=%13.6e\n", f_min, f_max);
    }
	}
}



// #################################################################
/**
 * @brief リスタート時の瞬時値ファイル読み込み
 * @param [in]  fp   ファイルポインタ
 * @param [out] flop 浮動小数点演算数
 */
void FFV::Restart_std(FILE* fp, double& flop)
{
  double time;
  unsigned step;
  std::string tmp, dtmp, fname;
  std::string fmt(C.file_fmt_ext);

  REAL_TYPE bp = ( C.Unit.Prs == Unit_Absolute ) ? C.BasePrs : 0.0;
  REAL_TYPE refD = C.RefDensity;
  REAL_TYPE refV = C.RefVelocity;
  int Dmode = C.Unit.File;
  
  // ガイド出力
  int gs = C.GuideOut;
  
  // dummy
  unsigned i_dummy=0;
  double f_dummy=0.0;
  
  // 入出力モード
  bool mio = (bool)C.FIO.IOmode;
  
  // 入力ディレクトリ
  dtmp = DFI.GenerateDirName(C.FIO.InDirPath, C.Restart_step, C.FIO.Slice);
  
  const int* m_div = paraMngr->GetDivNum();
  
  // 自身の領域終点インデックス
  int tail[3];
  for(int i=0;i<3;i++) tail[i]=head[i]+size[i]-1;
  
  REAL_TYPE r_time;
  DFI_IN_PRS->ReadData(C.Restart_step, guide, G_size, (int *)m_div, head, tail, d_p, r_time);
  if ( d_p == NULL ) Exit(0);
  time = (double)r_time;
  step = (unsigned)C.Restart_step;
  
  // ここでタイムスタンプを得る
  if (C.Unit.File == DIMENSIONAL) time /= C.Tscale;
  Session_StartStep = CurrentStep = step;
  Session_StartTime = CurrentTime = time;
  
  // v00[]に値をセット
  copyV00fromRF(Session_StartTime);
  
  DFI_IN_VEL->ReadData(C.Restart_step, guide, G_size, (int *)m_div, head, tail, d_wo, r_time);
  if( d_wo == NULL ) Exit(0);
  
  REAL_TYPE refv = (Dmode == DIMENSIONAL) ? refV : 1.0;
  REAL_TYPE scale = 1.0; // 瞬時値の時スケールは1.0
  REAL_TYPE u0[4];
  u0[0] = v00[0];
  u0[1] = v00[1];
  u0[2] = v00[2];
  u0[3] = v00[3];
  
  fb_shift_refv_in_(d_v, size, &guide, d_wo, u0, &scale, &refv, &flop);
  
  time = (double)r_time;
  step = (unsigned)C.Restart_step;
  
  if (C.Unit.File == DIMENSIONAL) time /= C.Tscale;
  
  if ( (step != Session_StartStep) || (time != Session_StartTime) )
  {
    Hostonly_ printf     ("\n\tTime stamp is different between files\n");
    Hostonly_ fprintf(fp, "\n\tTime stamp is different between files\n");
    Exit(0);
  }

  
  
  // Instantaneous Temperature fields
  if ( C.isHeatProblem() )
  {
    REAL_TYPE klv = ( C.Unit.Temp == Unit_KELVIN ) ? 0.0 : KELVIN;
    
    tmp = DFI.GenerateFileName(C.f_Temperature, fmt, C.Restart_step, myRank, mio);
    fname = dtmp + tmp;
    
    if ( !checkFile(fname) )
    {
      Hostonly_ printf("\n\tError : File open '%s'\n", fname.c_str());
      Exit(0);
    }
    F.readTemperature(fp, fname, size, guide, d_t, step, time, Dmode, C.BaseTemp, C.DiffTemp, klv, flop, gs, true, i_dummy, f_dummy);
    
    if (C.Unit.File == DIMENSIONAL) time /= C.Tscale;
    
    if ( (step != Session_StartStep) || (time != Session_StartTime) ) 
    {
      Hostonly_ printf     ("\n\tTime stamp is different between files\n");
      Hostonly_ fprintf(fp, "\n\tTime stamp is different between files\n");
      Exit(0);
    }
  }
  
}



// #################################################################
/**
 * @brief リスタート時の平均値ファイル読み込み
 * @param [in]  fp   ファイルポインタ
 * @param [out] flop 浮動小数点演算数
 */
void FFV::Restart_avrerage (FILE* fp, double& flop)
{
  std::string tmp, dtmp, fname;
  std::string fmt(C.file_fmt_ext);
  
  unsigned step = Session_StartStep;
  double   time = Session_StartTime;

  
  // ガイド出力
  int gs = C.GuideOut;
  
  if ( C.Interval[Interval_Manager::tg_average].isStep() )
  {
    if ( step >= C.Interval[Interval_Manager::tg_average].getStartStep() )
    {
      Hostonly_ printf     ("\tRestart from Previous Calculation Results of averaged field\n");
      Hostonly_ fprintf(fp, "\tRestart from Previous Calculation Results of averaged field\n");
      Hostonly_ printf     ("\tStep : base=%u current=%u\n", step, CurrentStep);
      Hostonly_ fprintf(fp, "\tStep : base=%u current=%u\n", step, CurrentStep);
    }
    else
    {
      return;
    }
  }
  else
  {
    if ( time >= C.Interval[Interval_Manager::tg_average].getStartTime() )
    {
      Hostonly_ printf     ("\tRestart from Previous Calculation Results of averaged field\n");
      Hostonly_ fprintf(fp, "\tRestart from Previous Calculation Results of averaged field\n");
      Hostonly_ printf     ("\tTime : base=%e[sec.]/%e[-] current=%e[-]\n", time*C.Tscale, time, CurrentTime);
      Hostonly_ fprintf(fp, "\tTime : base=%e[sec.]/%e[-] current=%e[-]\n", time*C.Tscale, time, CurrentTime);
    }
    else
    {
      return;
    }
  }
  
  unsigned step_avr = 0;
  double time_avr = 0.0;
  
  // 入出力モード
  bool mio = (bool)C.FIO.IOmode;
  
  // 入力ディレクトリ
  dtmp = DFI.GenerateDirName(C.FIO.InDirPath, C.Restart_stepAvr, C.FIO.Slice);
  
  
  // Pressure
  REAL_TYPE bp = ( C.Unit.Prs == Unit_Absolute ) ? C.BasePrs : 0.0;
  
  tmp = DFI.GenerateFileName(C.f_AvrPressure, fmt, C.Restart_stepAvr, myRank, mio);
  fname = dtmp + tmp;
  
  if ( !checkFile(fname) )
  {
    Hostonly_ printf("\n\tError : File open '%s'\n", fname.c_str());
    Exit(0);
  }
  F.readPressure(fp, fname, size, guide, d_ap, step, time, C.Unit.File, bp, C.RefDensity, C.RefVelocity, flop, gs, false, step_avr, time_avr);
  
  if ( (step != Session_StartStep) || (time != Session_StartTime) )
  {
    Hostonly_ printf     ("\n\tTime stamp is different between instantaneous and averaged files\n");
    Hostonly_ fprintf(fp, "\n\tTime stamp is different between instantaneous and averaged files\n");
    Exit(0);
  }
  
  CurrentStep_Avr = step_avr;
  CurrentTime_Avr = time_avr;
  
  
  // Velocity
  tmp = DFI.GenerateFileName(C.f_AvrVelocity, fmt, C.Restart_stepAvr, myRank, mio);
  fname = dtmp + tmp;
  
  if ( !checkFile(fname) )
  {
    Hostonly_ printf("\n\tError : File open '%s'\n", fname.c_str());
    Exit(0);
  }
  F.readVelocity(fp, fname, size, guide, d_av, d_wo, step, time, v00, C.Unit.File, C.RefVelocity, flop, gs, false, step_avr, time_avr);
  
  if ( (step_avr != CurrentStep_Avr) || (time_avr != CurrentTime_Avr) ) // 圧力とちがう場合
  {
    Hostonly_ printf     ("\n\tTime stamp is different between files\n");
    Hostonly_ fprintf(fp, "\n\tTime stamp is different between files\n");
    Exit(0);
  }
  
  
  if ( C.isHeatProblem() )
  {
    REAL_TYPE klv = ( C.Unit.Temp == Unit_KELVIN ) ? 0.0 : KELVIN;
    
    tmp = DFI.GenerateFileName(C.f_AvrTemperature, fmt, C.Restart_stepAvr, myRank, mio);
    fname = dtmp + tmp;
    
    if ( !checkFile(fname) )
    {
      Hostonly_ printf("\n\tError : File open '%s'\n", fname.c_str());
      Exit(0);
    }
    F.readTemperature(fp, fname, size, guide, d_at, step, time, C.Unit.File, C.BaseTemp, C.DiffTemp, klv, flop, gs, false, step_avr, time_avr);
    
    if ( (step_avr != CurrentStep_Avr) || (time_avr != CurrentTime_Avr) )
    {
      Hostonly_ printf     ("\n\tTime stamp is different between files\n");
      Hostonly_ fprintf(fp, "\n\tTime stamp is different between files\n");
      Exit(0);
    }
  }
}



// #################################################################
// 粗い格子を用いたリスタート
void FFV::Restart_coarse(FILE* fp, double& flop)
{
  std::string f_prs;
  std::string f_vel;
  std::string f_temp;
  std::string fmt(C.file_fmt_ext);
  
  // 出力モード
  bool mio = (bool)C.FIO.IOmode;
  
  // 粗格子の分割サイズ
  int r_size[3];
  
  // 粗格子サブドメインにおける読み込み開始インデクス
  int crs[3];
  
  // 粗格子サブドメインに含まれる密格子サブドメインの数
  int num_block[3];
  
  //並列時には各ランクに必要なファイル名と開始インデクスを取得
  if ( C.FIO.IOmode == IO_DISTRIBUTE )
  {
    int i, j, k;  // 密格子のグローバル開始インデクス
    i = head[0];
    j = head[1];
    k = head[2];
    
    // crs_i, _j, _kには同じ値が入る
    if ( !getCoarseResult2(i, j, k, C.f_dfi_prs, C.f_dfi_prfx_prs, C.Restart_step, f_prs, r_size, crs, num_block) )
    {
      Hostonly_ printf("\tError : Find invalid coarse sub-domain\n");
      Exit(0);
    }
    
    if ( !getCoarseResult2(i, j, k, C.f_dfi_vel, C.f_dfi_prfx_vel, C.Restart_step, f_vel, r_size, crs, num_block) )
    {
      Hostonly_ printf("\tError : Find invalid coarse sub-domain\n");
      Exit(0);
    }
    
    if ( C.isHeatProblem() )
    {
      if ( !getCoarseResult2(i, j, k, C.f_dfi_temp, C.f_dfi_prfx_temp, C.Restart_step, f_temp, r_size, crs, num_block) )
      {
        Hostonly_ printf("\tError : Find invalid coarse sub-domain\n");
        Exit(0);
      }
    }
  }
  else
  {
    crs[0] = 1;
    crs[1] = 1;
    crs[2] = 1;
    f_prs = DFI.GenerateFileName(C.f_dfi_prfx_prs, fmt, C.Restart_step, myRank, mio);
    f_vel = DFI.GenerateFileName(C.f_dfi_prfx_vel, fmt, C.Restart_step, myRank, mio);
    
    if ( C.isHeatProblem() )
    {
      f_temp= DFI.GenerateFileName(C.f_dfi_prfx_temp, fmt, C.Restart_step, myRank, mio);
    }
  }
  
  // テンポラリのファイルロード用メモリ領域
  double G_c_mem=0.0;
  double c_mem=0.0;
  
  allocArray_CoarseMesh(r_size, c_mem);
  
  display_memory_info(fp, G_c_mem, c_mem, "Coarse Mesh reading");
  
  
  // ガイド出力
  int gs = C.GuideOut;
  
  
  // dummy
  unsigned i_dummy=0;
  double f_dummy=0.0;
  
  unsigned step;
  double time;
  
  
  // 圧力の瞬時値　ここでタイムスタンプを得る
  REAL_TYPE bp = ( C.Unit.Prs == Unit_Absolute ) ? C.BasePrs : 0.0;
  
  if ( !checkFile(f_prs) )
  {
    Hostonly_ printf("\n\tError : File open '%s'\n", f_prs.c_str());
    Exit(0);
  }

  F.readPressure(fp, f_prs, r_size, guide, d_r_p, step, time, C.Unit.File, bp, C.RefDensity, C.RefVelocity, flop, gs, true, i_dummy, f_dummy);
  
  if (C.Unit.File == DIMENSIONAL) time /= C.Tscale;

  Session_StartStep = CurrentStep = step;
  Session_StartTime = CurrentTime = time;
  
  // v00[]に値をセット
  copyV00fromRF(Session_StartTime);

  
  // Instantaneous Velocity fields
  if ( !checkFile(f_vel) )
  {
    Hostonly_ printf("\n\tError : File open '%s'\n", f_vel.c_str());
    Exit(0);
  }
  // d_r_vはr_size, d_woはsize, rsize<sizeなのでバッファとして利用
  F.readVelocity(fp, f_vel, r_size, guide, d_r_v, d_wo, step, time, v00, C.Unit.File, C.RefVelocity, flop, gs, true, i_dummy, f_dummy);
  
  if (C.Unit.File == DIMENSIONAL) time /= C.Tscale;
  
  if ( (step != Session_StartStep) || (time != Session_StartTime) )
  {
    Hostonly_ printf     ("\n\tTime stamp is different between files\n");
    Hostonly_ fprintf(fp, "\n\tTime stamp is different between files\n");
    Exit(0);
  }
  
  // Instantaneous Temperature fields
  if ( C.isHeatProblem() )
  {
    REAL_TYPE klv = ( C.Unit.Temp == Unit_KELVIN ) ? 0.0 : KELVIN;
    
    if ( !checkFile(f_temp) )
    {
      Hostonly_ printf("\n\tError : File open '%s'\n", f_temp.c_str());
      Exit(0);
    }
    F.readTemperature(fp, f_temp, r_size, guide, d_r_t, step, time, C.Unit.File, C.BaseTemp, C.DiffTemp, klv, flop, gs, true, i_dummy, f_dummy);
    
    if (C.Unit.File == DIMENSIONAL) time /= C.Tscale;
    
    if ( (step != Session_StartStep) || (time != Session_StartTime) )
    {
      Hostonly_ printf     ("\n\tTime stamp is different between files\n");
      Hostonly_ fprintf(fp, "\n\tTime stamp is different between files\n");
      Exit(0);
    }
  }

  // 同期
  if ( paraMngr->BndCommS3D(d_r_p, r_size[0], r_size[1], r_size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
  if ( paraMngr->BndCommV3D(d_r_v, r_size[0], r_size[1], r_size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
  
  if ( C.isHeatProblem() )
  {
    if ( paraMngr->BndCommS3D(d_r_t, r_size[0], r_size[1], r_size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
  }
  
  // 内挿処理
  Interpolation_from_coarse_initial(crs, num_block);
  
  // Delete temporary array
  if ( d_r_p )
  {
    delete [] d_r_p;
    d_r_p=NULL;
  }
  
  if ( d_r_v )
  {
    delete [] d_r_v;
    d_r_v=NULL;
  }
  
  if ( C.isHeatProblem() )
  {
    if ( d_r_t )
    {
      delete [] d_r_t;
      d_r_t=NULL;
    }
  }
}



// #################################################################
// 並列分散時のファイル名の管理を行う
void FFV::setDFI()
{
  int* g_bbox_st = new int[3*numProc];
  int* g_bbox_ed = new int[3*numProc];
  
  // host nameの取得
  string host = paraMngr->GetHostName();
  if ( host.empty() ) Exit(0);
  
  
  // 並列時のみ
  if ( numProc > 1 )
  {
    const int* m_tail = paraMngr->GetVoxelTailIndex();
    int tail[3];
    
    tail[0] = m_tail[0] + 1;
    tail[1] = m_tail[1] + 1;
    tail[2] = m_tail[2] + 1;
    
    // 集約
    if ( paraMngr->Gather(head, 3, g_bbox_st, 3, 0) != CPM_SUCCESS ) Exit(0);
    if ( paraMngr->Gather(tail, 3, g_bbox_ed, 3, 0) != CPM_SUCCESS ) Exit(0);
    
  }
  else // serial
  {
    g_bbox_st[0] = 1;
    g_bbox_st[1] = 1;
    g_bbox_st[2] = 1;
    g_bbox_ed[0] = size[0];
    g_bbox_ed[1] = size[1];
    g_bbox_ed[2] = size[2];
  }
  
  
  std::string UnitL;
  
  switch (C.Unit.Length) {
    case LTH_ND:
      UnitL = "NonDimensional";
      break;
      
    case LTH_m:
      UnitL = "M";
      break;
      
    case LTH_cm:
      UnitL = "cm";
      break;
      
    case LTH_mm:
      UnitL = "mm";
      break;
      
    default:
      break;
  }
  
  
  std::string UnitV;
  
  if ( C.Unit.File == DIMENSIONAL )
  {
    UnitV = "m/s";
  }
  else
  {
    UnitV = "NonDimensional";
  }
  
  std::string UnitP;
  
  if ( C.Unit.File == DIMENSIONAL )
  {
    UnitP = "Pa";
  }
  else
  {
    UnitP = "NonDimensional";
  }
  
  // 無次元化の圧力差 != 動圧
  REAL_TYPE dp = C.RefDensity * C.RefVelocity * C.RefVelocity;
  
  // DFIクラスの初期化 >> 共通
  if ( !DFI.init(G_size,
                 paraMngr->GetDivNum(),
                 C.GuideOut,
                 C.Start,
                 C.RefLength,
                 C.RefVelocity,
                 C.BasePrs,
                 dp,
                 UnitL,
                 UnitV,
                 UnitP,
                 g_bbox_st,
                 g_bbox_ed,
                 host) ) Exit(0);
  
  
  // 後始末
  if ( g_bbox_st )
  {
    delete [] g_bbox_st;
    g_bbox_st = NULL;
  }
  
  if ( g_bbox_ed )
  {
    delete [] g_bbox_ed;
    g_bbox_ed = NULL;
  }
  
}



// #################################################################
// リスタート時の瞬時値ファイル読み込み（並列数が異なる場合）
void FFV::Restart_different(FILE* fp, double& flop)
{
  std::string fmt(C.file_fmt_ext);
  
  // 出力モード
  bool mio = (bool)C.FIO.IOmode;
  
  // 自身の領域終点インデックス
  int tail[3];
  for(int i=0;i<3;i++) tail[i]=head[i]+size[i]-1;
  
  // 読み込むdfiファイルの数
  int ndfi=2;
  if ( C.isHeatProblem() ) ndfi=3;
  
  // allocate dfi info class
  DfiInfo* DI;
  DI = new DfiInfo[ndfi];
  
  // set dfi info class
  vector<string> dfi_name;
  dfi_name.clear();
  if ( C.FIO.IOmode == IO_DISTRIBUTE ) // always IO_DISTRIBUTE
  {
    dfi_name.push_back(C.f_dfi_prs.c_str());
    dfi_name.push_back(C.f_dfi_vel.c_str());
    dfi_name.push_back(C.f_dfi_fvel.c_str());
    if ( C.isHeatProblem() ) dfi_name.push_back(C.f_dfi_temp.c_str());
  }
  
  // set dfi info class
  for(int ic=0;ic<ndfi;ic++){
    string fname = dfi_name[ic];
    DI[ic].ReadDfiFile(fname);
    ic++;
  }
  
  // 読み込むファイルを探しDRIに情報を格納する
  int nDRI=0;
  int max_nDRI=0;
  DifferentRestartInfo* DRI;
  if ( C.FIO.IOmode == IO_DISTRIBUTE ) // always IO_DISTRIBUTE
  {
    // 自身の領域が含まれる領域の数を数える
    for(int j=0; j< DI[0].NodeInfoSize; j++ ) {
      if( !(head[0] <= DI[0].Node[j].TailIndex[0] && //自身のすべての方向のheadが、読み込む側のすべての方向のtailより小さくなければ飛ばす
            head[1] <= DI[0].Node[j].TailIndex[1] &&
            head[2] <= DI[0].Node[j].TailIndex[2]) ) continue;
      if( !(tail[0] >= DI[0].Node[j].HeadIndex[0] && //自身のすべての方向のtailが、読み込む側のすべての方向のheadより大きくなければ飛ばす
            tail[1] >= DI[0].Node[j].HeadIndex[1] &&
            tail[2] >= DI[0].Node[j].HeadIndex[2]) ) continue;
      nDRI++;
    }
    
    // DRIのアロケート
    if(nDRI == 0 ){ //読み込む領域が見つからない場合はエラー  {
      Hostonly_ printf("\tError : cannot find area to read\n");
      Exit(0);
    }
    
    // 全ランクのnDRIの最大値を調べ、その最大値でDRIをアロケートする ---> STAGING対策
    max_nDRI=nDRI;
    if ( numProc > 1 )
    {
      MPI_Allreduce(&max_nDRI, &max_nDRI, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    }
    DRI = new DifferentRestartInfo[max_nDRI];
    
    // set DRI
    int ic=0;
    for(int j=0; j< DI[0].NodeInfoSize; j++ ) {
      if( !(head[0] <= DI[0].Node[j].TailIndex[0] && //自身のすべての方向のheadが、読み込む側のすべての方向のtailより小さくなければ飛ばす
            head[1] <= DI[0].Node[j].TailIndex[1] &&
            head[2] <= DI[0].Node[j].TailIndex[2]) ) continue;
      if( !(tail[0] >= DI[0].Node[j].HeadIndex[0] && //自身のすべての方向のtailが、読み込む側のすべての方向のheadより大きくなければ飛ばす
            tail[1] >= DI[0].Node[j].HeadIndex[1] &&
            tail[2] >= DI[0].Node[j].HeadIndex[2]) ) continue;
      
      int oh[3],ot[3];
      CalOverlap(oh,ot,head,tail,DI[0].Node[j].HeadIndex,DI[0].Node[j].TailIndex);
      
      DRI[ic].rank=DI[0].Node[j].RankID;// コンストラクタでrank=-1で初期化
      for(int i=0;i<3;i++) DRI[ic].overlap_head[i]=oh[i];
      for(int i=0;i<3;i++) DRI[ic].overlap_tail[i]=ot[i];
      for(int i=0;i<3;i++) DRI[ic].read_file_voxel_head[i]=DI[0].Node[j].HeadIndex[i];
      for(int i=0;i<3;i++) DRI[ic].read_file_voxel_tail[i]=DI[0].Node[j].TailIndex[i];
      for(int i=0;i<3;i++) DRI[ic].read_file_voxel_size[i]=DI[0].Node[j].VoxelSize[i];
      DRI[ic].f_prs  = DFI.GenerateFileName(C.f_different_nproc_pressure,  fmt, C.Restart_step, DI[0].Node[j].RankID, mio);
      DRI[ic].f_vel  = DFI.GenerateFileName(C.f_different_nproc_velocity,  fmt, C.Restart_step, DI[0].Node[j].RankID, mio);
      DRI[ic].f_fvel = DFI.GenerateFileName(C.f_different_nproc_fvelocity, fmt, C.Restart_step, DI[0].Node[j].RankID, mio);
      if ( C.isHeatProblem() ) DRI[ic].f_temp= DFI.GenerateFileName(C.f_different_nproc_temperature, fmt, C.Restart_step, DI[0].Node[j].RankID, mio);
      ic++;
    }
  }
  else
  {
    nDRI=1;
    max_nDRI=nDRI;
    DRI = new DifferentRestartInfo[max_nDRI];
    DRI[0].rank=0;
    for(int i=0;i<3;i++) DRI[0].overlap_head[i]=head[i];
    for(int i=0;i<3;i++) DRI[0].overlap_tail[i]=tail[i];
    for(int i=0;i<3;i++) DRI[0].read_file_voxel_head[i]=1;
    for(int i=0;i<3;i++) DRI[0].read_file_voxel_tail[i]=G_size[i];
    for(int i=0;i<3;i++) DRI[0].read_file_voxel_size[i]=G_size[i];
    DRI[0].f_prs  = DFI.GenerateFileName(C.f_different_nproc_pressure, fmt, C.Restart_step, myRank, mio);
    DRI[0].f_vel  = DFI.GenerateFileName(C.f_different_nproc_velocity, fmt, C.Restart_step, myRank, mio);
    DRI[0].f_fvel = DFI.GenerateFileName(C.f_different_nproc_fvelocity, fmt, C.Restart_step, myRank, mio);
    if ( C.isHeatProblem() ) DRI[0].f_temp= DFI.GenerateFileName(C.f_different_nproc_temperature, fmt, C.Restart_step, myRank, mio);
  }
  
  // ランクごとに前計算時のどのランクのファイルを必要とするかテーブルで持っておく
  int frank_size=numProc*max_nDRI;
  
  int* frank = new int[frank_size];
  
  for(int i=0;i<frank_size;i++) frank[i]=0;
  
  if ( C.FIO.IOmode == IO_DISTRIBUTE ) // Gather出力されていた場合は0ランクを見に行くのでMPI_SUM必要なし
  {
    for(int ic=0;ic<nDRI;ic++){
      frank[myRank*max_nDRI+ic]=DRI[ic].rank;
    }
    
    int* tmp = new int[frank_size];
    
    for(int i=0;i<frank_size;i++) tmp[i]=frank[i];
    
    if ( numProc > 1 )
    {
      MPI_Allreduce(tmp, frank, frank_size, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    }
    
    if ( tmp ) delete [] tmp;
  }
  
  // ランクごとに前計算時のどのファイルが自身のランクにステージングされているかリストしておく
  int numProc_last=DI[0].NumberOfRank;
  int nassign=numProc_last/numProc;
  if( (numProc_last%numProc) > 0 ) nassign++;
  
  int* assign = new int[nassign];
  
  for(int i=0;i<nassign;i++) assign[i]=-1;
  
  int irank=0;
  int iassign=0;
  for(int i=0;i<numProc_last;i++){
    if(irank==myRank) assign[iassign]=i;
    irank++;
    if(irank==numProc){
      iassign++;
      irank=0;
    }
  }
  
  // 読み込みようワークエリアの最大サイズを決定
  //REAL_TYPE *d_d_v;  ///< 速度
  //REAL_TYPE *d_d_p;  ///< 圧力
  //REAL_TYPE *d_d_t;  ///< 温度
  REAL_TYPE *d_wk;
  unsigned long long d_sz=0;
  for(int ic=0;ic<nDRI;ic++){// 自身が読むファイルのボクセルサイズの最大値
    if( d_sz < (DRI[ic].read_file_voxel_size[0]+2*guide)
       * (DRI[ic].read_file_voxel_size[1]+2*guide)
       * (DRI[ic].read_file_voxel_size[2]+2*guide) )
    {
      d_sz = (DRI[ic].read_file_voxel_size[0]+2*guide)
      * (DRI[ic].read_file_voxel_size[1]+2*guide)
      * (DRI[ic].read_file_voxel_size[2]+2*guide);
    }
  }
  for(int i=0;i<nassign;i++){// 自身のランクにステージングされているファイルのボクセルサイズの最大値
    if(assign[i]==-1) continue;
    if( d_sz < (DI[0].Node[assign[i]].VoxelSize[0]+2*guide)
       * (DI[0].Node[assign[i]].VoxelSize[1]+2*guide)
       * (DI[0].Node[assign[i]].VoxelSize[2]+2*guide) )
    {
      d_sz = (DI[0].Node[assign[i]].VoxelSize[0]+2*guide)
      * (DI[0].Node[assign[i]].VoxelSize[1]+2*guide)
      * (DI[0].Node[assign[i]].VoxelSize[2]+2*guide);
    }
  }
  d_sz=d_sz*3;
  
  
  // 読み込みwk用メモリ計算
  double G_d_mem=0.0;
  double mc=(double)d_sz;
  if(mc>(double)INT_MAX){// 整数値あふれ出しチェック --- 参考 894*894*894*3=2143550952 INT_MAX 2147483647
    printf("\tsize error : wkmaxsize>INT_MAX\n");
    Exit(0);
  }
  display_memory_info(fp, G_d_mem, mc, "Different Nproc Mesh reading");
  
  
  // ワークエリアの確保
  int wksize=(int)d_sz;
  if ( !(d_wk = new REAL_TYPE[ wksize ] ) ){
    Hostonly_  printf(    "\t>> cannot allocate work area : restart from different nproc\n\n");
    Exit(0);
  }
  
  
  
#if 0
  
  //並列デバッグ用ファイルポインタ
  FILE* ifdg;
  
  // 並列デバッグ用ランク別ファイルオープン
  int len = 6;
  char* tmp = new char[len];
  memset(tmp, 0, sizeof(char)*len);
  sprintf(tmp, "%06d", myRank);
  std::string buff(tmp);
  buff = buff+"_debug.txt";
  if( !(ifdg = fopen(buff.c_str(), "w")) )
  {
    fprintf(stderr, "Can't open file.(%s)\n", buff.c_str());
    Exit(0);
  }
  
  fprintf(ifdg,"\n");
  fprintf(ifdg,"*** debug file *** rank == %d\n", myRank);
  fprintf(ifdg,"\n");
  for(int ic=0;ic<nDRI;ic++){
    fprintf(ifdg,"\n");
    fprintf(ifdg,"\tDRI[%2d].rank = %d\n",ic,DRI[ic].rank);
    for(int i=0;i<3;i++){
      fprintf(ifdg,"\tDRI[%2d].overlap_head[%d] = %d\n",ic,i,DRI[ic].overlap_head[i]);
    }
    for(int i=0;i<3;i++){
      fprintf(ifdg,"\tDRI[%2d].overlap_tail[%d] = %d\n",ic,i,DRI[ic].overlap_tail[i]);
    }
    for(int i=0;i<3;i++){
      fprintf(ifdg,"\tDRI[%2d].read_file_voxel_head[%d] = %d\n",ic,i,DRI[ic].read_file_voxel_head[i]);
    }
    for(int i=0;i<3;i++){
      fprintf(ifdg,"\tDRI[%2d].read_file_voxel_tail[%d] = %d\n",ic,i,DRI[ic].read_file_voxel_tail[i]);
    }
    for(int i=0;i<3;i++){
      fprintf(ifdg,"\tDRI[%2d].read_file_voxel_size[%d] = %d\n",ic,i,DRI[ic].read_file_voxel_size[i]);
    }
    fprintf(ifdg,"\tDRI[%2d].f_prs = %s\n",ic,DRI[ic].f_prs.c_str());
    fprintf(ifdg,"\tDRI[%2d].f_vel = %s\n",ic,DRI[ic].f_vel.c_str());
    if ( C.isHeatProblem() ) fprintf(ifdg,"\tDRI[%2d].f_temp = %s\n",ic,DRI[ic].f_temp.c_str());
  }
  
  fprintf(ifdg,"\n");
  fprintf(ifdg,"*** frank table *** rank = %d\n", myRank);
  fprintf(ifdg,"\n");
  fprintf(ifdg,"\tnumProc    = %d\n", numProc);
  fprintf(ifdg,"\tnDRI       = %d\n", nDRI);
  fprintf(ifdg,"\tfrank_size = %d\n", frank_size);
  fprintf(ifdg,"\tmax_nDRI   = %d\n", max_nDRI);
  for(int i=0;i<frank_size;i++) fprintf(ifdg,"frank[%d] = %d\n",i,frank[i]);
  
  fprintf(ifdg,"\n");
  fprintf(ifdg,"*** assign list *** rank = %d\n", myRank);
  fprintf(ifdg,"\tnassign = %d\n", nassign);
  for(int i=0;i<nassign;i++) fprintf(ifdg,"\tassign[%d] = %d\n",i,assign[i]);
  
#endif
  
  if( C.Restart_staging )
  {
    for(int ic=0;ic<max_nDRI;ic++){
      
      // 各プロセスが今読みたいファイルの前計算時のリストを作る
      int* rank_list = new int[numProc];
      
      for(int i=0;i<numProc;i++) rank_list[i]=0;
      rank_list[myRank]=DRI[ic].rank;
      
      if ( C.FIO.IOmode == IO_DISTRIBUTE )
      {
        int *tmp = new int[numProc];
        for(int i=0;i<numProc;i++) tmp[i]=0;
        tmp[myRank]=DRI[ic].rank;
        if ( numProc > 1 )
        {
          MPI_Allreduce(tmp, rank_list, numProc, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        }
        if ( tmp) delete [] tmp;
      }
      else
      {
        //通信の必要なし ---> すべて0のはず
      }
      
      int recv_rank=DRI[ic].rank%numProc;// 自分が読みたいファイルがステージングされているランク
      if(DRI[ic].rank < 0) recv_rank=-1;
      
      // オーバーラップ分をファイルから読み込み
      ReadOverlap_Pressure(fp, flop, &DRI[ic], DI, d_wk, rank_list, recv_rank, assign, nassign);
      ReadOverlap_Velocity(fp, flop, &DRI[ic], DI, d_wk, rank_list, recv_rank, assign, nassign);
      ReadOverlap_FVelocity(fp, flop, &DRI[ic], DI, d_wk, rank_list, recv_rank, assign, nassign);
      if ( C.isHeatProblem() ) ReadOverlap_Temperature(fp, flop, &DRI[ic], DI, d_wk, rank_list, recv_rank, assign, nassign);
      
      if ( rank_list ) delete [] rank_list;
    }
  }
  else
  {
    // オーバーラップ分をファイルから読み込み
    for(int ic=0;ic<max_nDRI;ic++){
      ReadOverlap(fp, flop, &DRI[ic], d_wk);
    }
  }
  
  
#if 0
  if (ifdg) fclose(ifdg);
#endif
  
  if ( DI ) delete [] DI;
  if ( DRI ) delete [] DRI;
  if ( d_wk ) delete [] d_wk;
  if ( frank ) delete [] frank;
  if ( assign ) delete [] assign;
  
}


// #################################################################
void FFV::CalOverlap(int* overlap_h, int* overlap_t, int* h, int* t, int* head, int* tail)
{
  // pattern 1
  //        h|-----|t
  //   head|---|tail
  //
  // pattern 2
  //        h|-----|t
  //      head|---|tail
  //
  // pattern 3
  //        h|-----|t
  //         head|---|tail
  //
  // pattern 4
  //        h|-----|t
  //    head|-------|tail
  //
  
  int flag1,flag2;
  for(int i=0;i<3;i++){
    flag1=OFF;
    if(head[i] >= h[i]) flag1=ON;
    flag2=OFF;
    if(tail[i] <= t[i]) flag2=ON;
    
    if( flag1==OFF && flag2==ON)     // pattern 1
    {
      overlap_h[i]=h[i];
      overlap_t[i]=tail[i];
    }
    else if( flag1==ON && flag2==ON) // pattern 2
    {
      overlap_h[i]=head[i];
      overlap_t[i]=tail[i];
    }
    else if( flag1==ON && flag2==OFF) // pattern 3
    {
      overlap_h[i]=head[i];
      overlap_t[i]=t[i];
    }
    else if( flag1==OFF && flag2==OFF) // pattern 4
    {
      overlap_h[i]=h[i];
      overlap_t[i]=t[i];
    }
  }
  
}


// #################################################################
void FFV::SetOverlap(REAL_TYPE* write_wk, REAL_TYPE* read_wk, int dim, int gd,
                     int* h, int* s, int* overlap_h, int* overlap_t, int* head, int* size)
{
  int isx = overlap_h[0]-h[0];
  int isy = overlap_h[1]-h[1];
  int isz = overlap_h[2]-h[2];
  size_t ips = (size_t)dim * (size_t)( (s[0]+2*gd)*(s[1]+2*gd)*(isz+gd)+(s[0]+2*gd)*(isy+gd)+(isx+gd) );
  
  int isxr = overlap_h[0]-head[0];
  int isyr = overlap_h[1]-head[1];
  int iszr = overlap_h[2]-head[2];
  size_t ipsr = (size_t)dim * (size_t)( (size[0]+2*gd)*(size[1]+2*gd)*(iszr+gd) +(size[0]+2*gd)*(isyr+gd) +(isxr+gd) );
  
  int szx = overlap_t[0]-overlap_h[0]+1;
  int szy = overlap_t[1]-overlap_h[1]+1;
  int szz = overlap_t[2]-overlap_h[2]+1;
  
  for(int k=0;k<szz;k++){
    for(int j=0;j<szy;j++){
      for(int i=0;i<szx;i++){
        for(int id=0;id<dim;id++){
          size_t ip  = ips  + (size_t)dim*( _F_IDX_S3D(i+1, j+1, k+1, (s[0]+2*gd), (s[1]+2*gd), (s[2]+2*gd), 0) ) + (size_t)id;
          size_t ipr = ipsr + (size_t)dim*( _F_IDX_S3D(i+1, j+1, k+1, (size[0]+2*gd), (size[1]+2*gd), (size[2]+2*gd), 0) ) + (size_t)id;
          write_wk[ip] = read_wk[ipr];
        }
      }
    }
  }
}


// #################################################################
void FFV::ReadOverlap(FILE* fp, double& flop, DifferentRestartInfo* DRI, REAL_TYPE* d_wk)
{
  // ファイルの読み込み
  
  double time;
  unsigned step;
  //string tmp;
  
  // 読み込み領域の格子の分割サイズ
  int d_size[3];
  for(int i=0;i<3;i++){
    d_size[i] = DRI->read_file_voxel_size[i];
  }
  
  // 圧力の瞬時値
  REAL_TYPE bp = ( C.Unit.Prs == Unit_Absolute ) ? C.BasePrs : 0.0;
  
  // ガイド出力
  int gs = C.GuideOut;
  
  // dummy
  unsigned i_dummy=0;
  double f_dummy=0.0;
  
  
  // 圧力のオーバーラップ分を移す
  if ( !checkFile(DRI->f_prs) )
  {
    Hostonly_ printf("\n\tError : File open '%s'\n", DRI->f_prs.c_str());
    Exit(0);
  }
  F.readPressure(fp, DRI->f_prs, d_size, guide, d_wk, step, time, C.Unit.File, bp, C.RefDensity, C.RefVelocity, flop, gs, true, i_dummy, f_dummy);
  SetOverlap(d_p,d_wk,1,guide,head,size,DRI->overlap_head,DRI->overlap_tail,
             DRI->read_file_voxel_head,DRI->read_file_voxel_size);
  
  
  // ここでタイムスタンプを得る
  if (C.Unit.File == DIMENSIONAL) time /= C.Tscale;
  Session_StartStep = CurrentStep = step;
  Session_StartTime = CurrentTime = time;
  
  // v00[]に値をセット
  copyV00fromRF(Session_StartTime);
  
  
  // 速度のオーバーラップ分を移す
  if ( !checkFile(DRI->f_vel) )
  {
    Hostonly_ printf("\n\tError : File open '%s'\n", DRI->f_vel.c_str());
    Exit(0);
  }
  
  // ？？？　ベクトルなのになぜd_wk、多分バグ
  F.readVelocity(fp, DRI->f_vel, d_size, guide, d_wk, d_wo, step, time, v00, C.Unit.File, C.RefVelocity, flop, gs, true, i_dummy, f_dummy);
  SetOverlap(d_v,d_wk,3,guide,head,size,DRI->overlap_head,DRI->overlap_tail,
             DRI->read_file_voxel_head,DRI->read_file_voxel_size);
  
  if (C.Unit.File == DIMENSIONAL) time /= C.Tscale;
  
  if ( (step != Session_StartStep) || (time != Session_StartTime) )
  {
    Hostonly_ printf     ("\n\tTime stamp is different between files\n");
    Hostonly_ fprintf(fp, "\n\tTime stamp is different between files\n");
    Exit(0);
  }
  

  // Face Velocity のオーバーラップ分を移す
  if ( !checkFile(DRI->f_fvel) )
  {
    Hostonly_ printf("\n\tError : File open '%s'\n", DRI->f_fvel.c_str());
    Exit(0);
  }
  
  F.readVelocity(fp, DRI->f_fvel, d_size, guide, d_wk, d_wo, step, time, v00, C.Unit.File, C.RefVelocity, flop, gs, true, i_dummy, f_dummy);
  SetOverlap(d_vf,d_wk,3,guide,head,size,DRI->overlap_head,DRI->overlap_tail,
             DRI->read_file_voxel_head,DRI->read_file_voxel_size);

  if (C.Unit.File == DIMENSIONAL) time /= C.Tscale;
  
  if ( (step != Session_StartStep) || (time != Session_StartTime) )
  {
    Hostonly_ printf     ("\n\tTime stamp is different between files\n");
    Hostonly_ fprintf(fp, "\n\tTime stamp is different between files\n");
    Exit(0);
  }


  // Instantaneous Temperature fields
  if ( C.isHeatProblem() )
  {
    REAL_TYPE klv = ( C.Unit.Temp == Unit_KELVIN ) ? 0.0 : KELVIN;
    
    // 温度のオーバーラップ分を移す
    if ( !checkFile(DRI->f_temp) )
    {
      Hostonly_ printf("\n\tError : File open '%s'\n", DRI->f_temp.c_str());
      Exit(0);
    }
    F.readTemperature(fp, DRI->f_temp, d_size, guide, d_wk, step, time, C.Unit.File, C.BaseTemp, C.DiffTemp, klv, flop, gs, true, i_dummy, f_dummy);
    SetOverlap(d_t,d_wk,1,guide,head,size,DRI->overlap_head,DRI->overlap_tail,
               DRI->read_file_voxel_head,DRI->read_file_voxel_size);
    
    if (C.Unit.File == DIMENSIONAL) time /= C.Tscale;
    
    if ( (step != Session_StartStep) || (time != Session_StartTime) )
    {
      Hostonly_ printf     ("\n\tTime stamp is different between files\n");
      Hostonly_ fprintf(fp, "\n\tTime stamp is different between files\n");
      Exit(0);
    }
  }
  
}


// #################################################################
void FFV::ReadOverlap_Pressure(FILE* fp, double& flop, DifferentRestartInfo* DRI, DfiInfo* DI, REAL_TYPE* d_wk,
                               int* rank_list, int recv_rank, int* assign, int nassign)
{
  MPI_Request mpi_req;
  double time;
  unsigned step;
  std::string tmp;
  
  std::string fmt(C.file_fmt_ext);
  
  // 出力モード
  bool mio = (bool)C.FIO.IOmode;
  
  // set rank_dir
  int len = 6;
  char* buff = new char[len];
  memset(buff, 0, sizeof(char)*len);
  sprintf(buff, "%06d", myRank);
  std::string rank_dir(buff);
  if ( buff ) delete [] buff;
  
  // 圧力の瞬時値
  REAL_TYPE bp = ( C.Unit.Prs == Unit_Absolute ) ? C.BasePrs : 0.0;
  
  // ガイド出力
  int gs = C.GuideOut;
  
  // dummy
  unsigned i_dummy=0;
  double f_dummy=0.0;
  
  // 読み込み領域の格子の分割サイズ
  int d_size[3];
  
  // ファイルの送信
  for(int k=0;k<numProc;k++){
    if(rank_list[k]==-1) continue;
    for(int j=0;j<nassign;j++){
      if(assign[j]==-1) continue;
      
      if(assign[j]==rank_list[k]){ // 自身のランクにステージングされたファイルに一致するものがあれば
        int send_rank=k;
        int read_rank=rank_list[k];
        
        if(send_rank==myRank) continue; // 自身のランクであれば送る必要なし
        
        // 読み込み領域の格子の分割サイズ
        for(int i=0;i<3;i++){
          d_size[i] = DI[0].Node[read_rank].VoxelSize[i];
        }
        
        double check = (double)(d_size[0]+2*guide)*(double)(d_size[1]+2*guide)*(double)(d_size[2]+2*guide);
        if(check>(double)INT_MAX){
          printf("\tsize error : sendsize>INT_MAX\n");
          Exit(0);
        }
        
        unsigned sendsize=(d_size[0]+2*guide)*(d_size[1]+2*guide)*(d_size[2]+2*guide);
        for(unsigned i=0;i<sendsize;i++) d_wk[i] = 0.0;
        
        // ファイルの読み込み
        tmp = DFI.GenerateFileName(C.f_different_nproc_pressure, fmt, C.Restart_step, read_rank, mio);
        tmp = C.f_different_nproc_dir_prefix + rank_dir + "/" + tmp;
        if ( !checkFile(tmp) )
        {
          printf("\n\tError : File open '%s'\n", tmp.c_str());
          Exit(0);
        }
        F.readPressure(fp, tmp, d_size, guide, d_wk, step, time, C.Unit.File, bp, C.RefDensity, C.RefVelocity, flop, gs, true, i_dummy, f_dummy);
        
        // 送信
        if (paraMngr->Isend( &time, 1, send_rank, &mpi_req ) != CPM_SUCCESS) Exit(0);
        if (paraMngr->Wait( &mpi_req ) != MPI_SUCCESS) Exit(0);
        if (paraMngr->Isend( &step, 1, send_rank, &mpi_req ) != CPM_SUCCESS) Exit(0);
        if (paraMngr->Wait( &mpi_req ) != MPI_SUCCESS) Exit(0);
        if (paraMngr->Isend( d_wk, sendsize, send_rank, &mpi_req ) != CPM_SUCCESS) Exit(0);
        if (paraMngr->Wait( &mpi_req ) != MPI_SUCCESS) Exit(0);
      }
    }
  }
  
  // 送信だけしてリターン
  if(DRI->rank==-1) return;
  
  // ファイルの読み込み（受信）
  
  // 読み込み領域の格子の分割サイズ
  for(int i=0;i<3;i++){
    d_size[i] = DRI->read_file_voxel_size[i];
  }
  
  // 圧力のオーバーラップ分を移す
  if(recv_rank==myRank)
  {
    // ファイルから読み込み
    tmp = C.f_different_nproc_dir_prefix + rank_dir + "/" + DRI->f_prs;
    if ( !checkFile(tmp) )
    {
      printf("\n\tError : File open '%s'\n", tmp.c_str());
      Exit(0);
    }
    F.readPressure(fp, tmp, d_size, guide, d_wk, step, time, C.Unit.File, bp, C.RefDensity, C.RefVelocity, flop, gs, true, i_dummy, f_dummy);
  }
  else
  {
    // 受信
    unsigned recvsize=(d_size[0]+2*guide)*(d_size[1]+2*guide)*(d_size[2]+2*guide);
    for(unsigned i=0;i<recvsize;i++) d_wk[i] = 0.0;
    if (paraMngr->Irecv( &time, 1, recv_rank, &mpi_req ) != CPM_SUCCESS) Exit(0);
    if (paraMngr->Wait( &mpi_req ) != MPI_SUCCESS) Exit(0);
    if (paraMngr->Irecv( &step, 1, recv_rank, &mpi_req ) != CPM_SUCCESS) Exit(0);
    if (paraMngr->Wait( &mpi_req ) != MPI_SUCCESS) Exit(0);
    if (paraMngr->Irecv( d_wk, recvsize, recv_rank, &mpi_req ) != CPM_SUCCESS) Exit(0);
    if (paraMngr->Wait( &mpi_req ) != MPI_SUCCESS) Exit(0);
  }
  SetOverlap(d_p,d_wk,1,guide,head,size,DRI->overlap_head,DRI->overlap_tail,
             DRI->read_file_voxel_head,DRI->read_file_voxel_size);
  
  // ここでタイムスタンプを得る
  if (C.Unit.File == DIMENSIONAL) time /= C.Tscale;
  Session_StartStep = CurrentStep = step;
  Session_StartTime = CurrentTime = time;
  
  // v00[]に値をセット
  copyV00fromRF(Session_StartTime);
}


// #################################################################
void FFV::ReadOverlap_Velocity(FILE* fp, double& flop, DifferentRestartInfo* DRI, DfiInfo* DI, REAL_TYPE* d_wk,
                               int* rank_list, int recv_rank, int* assign, int nassign)
{
  MPI_Request mpi_req;
  double time;
  unsigned step;
  string tmp;
  std::string fmt(C.file_fmt_ext);
  
  // 出力モード
  bool mio = (bool)C.FIO.IOmode;
  
  // set rank_dir
  int len = 6;
  char* buff = new char[len];
  memset(buff, 0, sizeof(char)*len);
  sprintf(buff, "%06d", myRank);
  std::string rank_dir(buff);
  if ( buff ) delete [] buff;
  
  // ガイド出力
  int gs = C.GuideOut;
  
  // dummy
  unsigned i_dummy=0;
  double f_dummy=0.0;
  
  // 読み込み領域の格子の分割サイズ
  int d_size[3];
  
  // ファイルの送信
  for(int k=0;k<numProc;k++){
    if(rank_list[k]==-1) continue;
    for(int j=0;j<nassign;j++){
      if(assign[j]==-1) continue;
      
      if(assign[j]==rank_list[k]){ // 自身のランクにステージングされたファイルに一致するものがあれば
        int send_rank=k;
        int read_rank=rank_list[k];
        
        if(send_rank==myRank) continue; // 自身のランクであれば送る必要なし
        
        // 読み込み領域の格子の分割サイズ
        for(int i=0;i<3;i++){
          d_size[i] = DI[0].Node[read_rank].VoxelSize[i];
        }
        
        double check = (double)(d_size[0]+2*guide)*(double)(d_size[1]+2*guide)*(double)(d_size[2]+2*guide);
        if(check>(double)INT_MAX){
          printf("\tsize error : sendsize>INT_MAX\n");
          Exit(0);
        }
        
        unsigned sendsize=(d_size[0]+2*guide)*(d_size[1]+2*guide)*(d_size[2]+2*guide)*3;
        for(unsigned i=0;i<sendsize;i++) d_wk[i] = 0.0;
        
        // ファイルの読み込み
        tmp = DFI.GenerateFileName(C.f_different_nproc_velocity, fmt, C.Restart_step, read_rank, mio);
        tmp = C.f_different_nproc_dir_prefix + rank_dir + "/" + tmp;
        if ( !checkFile(tmp) )
        {
          printf("\n\tError : File open '%s'\n", tmp.c_str());
          Exit(0);
        }
        F.readVelocity(fp, tmp, d_size, guide, d_wk, d_wo, step, time, v00, C.Unit.File, C.RefVelocity, flop, gs, true, i_dummy, f_dummy);
        
        if (C.Unit.File == DIMENSIONAL) time /= C.Tscale;
        if ( (step != Session_StartStep) || (time != Session_StartTime) )
        {
          Hostonly_ printf     ("\n\tTime stamp is different between files\n");
          Hostonly_ fprintf(fp, "\n\tTime stamp is different between files\n");
          Exit(0);
        }
        
        // 送信
        if (paraMngr->Isend( &time, 1, send_rank, &mpi_req ) != CPM_SUCCESS) Exit(0);
        if (paraMngr->Wait( &mpi_req ) != MPI_SUCCESS) Exit(0);
        if (paraMngr->Isend( &step, 1, send_rank, &mpi_req ) != CPM_SUCCESS) Exit(0);
        if (paraMngr->Wait( &mpi_req ) != MPI_SUCCESS) Exit(0);
        if (paraMngr->Isend( d_wk, sendsize, send_rank, &mpi_req ) != CPM_SUCCESS) Exit(0);
        if (paraMngr->Wait( &mpi_req ) != MPI_SUCCESS) Exit(0);
      }
    }
  }
  
  // 送信だけしてリターン
  if(DRI->rank==-1) return;
  
  // ファイルの読み込み（受信）
  
  // 読み込み領域の格子の分割サイズ
  for(int i=0;i<3;i++){
    d_size[i] = DRI->read_file_voxel_size[i];
  }
  
  // 速度のオーバーラップ分を移す
  if(recv_rank==myRank)
  {
    // ファイルから読み込み
    tmp = C.f_different_nproc_dir_prefix + rank_dir + "/" + DRI->f_vel;
    if ( !checkFile(tmp) )
    {
      Hostonly_ printf("\n\tError : File open '%s'\n", tmp.c_str());
      Exit(0);
    }
    F.readVelocity(fp, tmp, d_size, guide, d_wk, d_wo, step, time, v00, C.Unit.File, C.RefVelocity, flop, gs, true, i_dummy, f_dummy);
    
    if (C.Unit.File == DIMENSIONAL) time /= C.Tscale;
    
    if ( (step != Session_StartStep) || (time != Session_StartTime) )
    {
      Hostonly_ printf     ("\n\tTime stamp is different between files\n");
      Hostonly_ fprintf(fp, "\n\tTime stamp is different between files\n");
      Exit(0);
    }
  }
  else
  {
    // 受信
    unsigned recvsize=(d_size[0]+2*guide)*(d_size[1]+2*guide)*(d_size[2]+2*guide)*3;
    for(unsigned i=0;i<recvsize;i++) d_wk[i] = 0.0;
    if (paraMngr->Irecv( &time, 1, recv_rank, &mpi_req ) != CPM_SUCCESS) Exit(0);
    if (paraMngr->Wait( &mpi_req ) != MPI_SUCCESS) Exit(0);
    if (paraMngr->Irecv( &step, 1, recv_rank, &mpi_req ) != CPM_SUCCESS) Exit(0);
    if (paraMngr->Wait( &mpi_req ) != MPI_SUCCESS) Exit(0);
    if (paraMngr->Irecv( d_wk, recvsize, recv_rank, &mpi_req ) != CPM_SUCCESS) Exit(0);
    if (paraMngr->Wait( &mpi_req ) != MPI_SUCCESS) Exit(0);
  }
  SetOverlap(d_v,d_wk,3,guide,head,size,DRI->overlap_head,DRI->overlap_tail,
             DRI->read_file_voxel_head,DRI->read_file_voxel_size);
  
}


// #################################################################
void FFV::ReadOverlap_FVelocity(FILE* fp, double& flop, DifferentRestartInfo* DRI, DfiInfo* DI, REAL_TYPE* d_wk,
                               int* rank_list, int recv_rank, int* assign, int nassign)
{
  MPI_Request mpi_req;
  double time;
  unsigned step;
  string tmp;
  
  std::string fmt(C.file_fmt_ext);
  
  // 出力モード
  bool mio = (bool)C.FIO.IOmode;
  
  // set rank_dir
  int len = 6;
  char* buff = new char[len];
  memset(buff, 0, sizeof(char)*len);
  sprintf(buff, "%06d", myRank);
  std::string rank_dir(buff);
  if ( buff ) delete [] buff;
  
  // ガイド出力
  int gs = C.GuideOut;
  
  // dummy
  unsigned i_dummy=0;
  double f_dummy=0.0;
  
  // 読み込み領域の格子の分割サイズ
  int d_size[3];
  
  // ファイルの送信
  for(int k=0;k<numProc;k++){
    if(rank_list[k]==-1) continue;
    for(int j=0;j<nassign;j++){
      if(assign[j]==-1) continue;
      
      if(assign[j]==rank_list[k]){ // 自身のランクにステージングされたファイルに一致するものがあれば
        int send_rank=k;
        int read_rank=rank_list[k];
        
        if(send_rank==myRank) continue; // 自身のランクであれば送る必要なし
        
        // 読み込み領域の格子の分割サイズ
        for(int i=0;i<3;i++){
          d_size[i] = DI[0].Node[read_rank].VoxelSize[i];
        }
        
        double check = (double)(d_size[0]+2*guide)*(double)(d_size[1]+2*guide)*(double)(d_size[2]+2*guide);
        if(check>(double)INT_MAX){
          printf("\tsize error : sendsize>INT_MAX\n");
          Exit(0);
        }
        
        unsigned sendsize=(d_size[0]+2*guide)*(d_size[1]+2*guide)*(d_size[2]+2*guide)*3;
        for(unsigned i=0;i<sendsize;i++) d_wk[i] = 0.0;
        
        // ファイルの読み込み
        tmp = DFI.GenerateFileName(C.f_different_nproc_fvelocity, fmt, C.Restart_step, read_rank, mio);
        tmp = C.f_different_nproc_dir_prefix + rank_dir + "/" + tmp;
        if ( !checkFile(tmp) )
        {
          printf("\n\tError : File open '%s'\n", tmp.c_str());
          Exit(0);
        }
        F.readVelocity(fp, tmp, d_size, guide, d_wk, d_wo, step, time, v00, C.Unit.File, C.RefVelocity, flop, gs, true, i_dummy, f_dummy);
        
        if (C.Unit.File == DIMENSIONAL) time /= C.Tscale;
        if ( (step != Session_StartStep) || (time != Session_StartTime) )
        {
          Hostonly_ printf     ("\n\tTime stamp is different between files\n");
          Hostonly_ fprintf(fp, "\n\tTime stamp is different between files\n");
          Exit(0);
        }
        
        // 送信
        if (paraMngr->Isend( &time, 1, send_rank, &mpi_req ) != CPM_SUCCESS) Exit(0);
        if (paraMngr->Wait( &mpi_req ) != MPI_SUCCESS) Exit(0);
        if (paraMngr->Isend( &step, 1, send_rank, &mpi_req ) != CPM_SUCCESS) Exit(0);
        if (paraMngr->Wait( &mpi_req ) != MPI_SUCCESS) Exit(0);
        if (paraMngr->Isend( d_wk, sendsize, send_rank, &mpi_req ) != CPM_SUCCESS) Exit(0);
        if (paraMngr->Wait( &mpi_req ) != MPI_SUCCESS) Exit(0);
      }
    }
  }
  
  // 送信だけしてリターン
  if(DRI->rank==-1) return;
  
  // ファイルの読み込み（受信）
  
  // 読み込み領域の格子の分割サイズ
  for(int i=0;i<3;i++){
    d_size[i] = DRI->read_file_voxel_size[i];
  }
  
  // 速度のオーバーラップ分を移す
  if(recv_rank==myRank)
  {
    // ファイルから読み込み
    tmp = C.f_different_nproc_dir_prefix + rank_dir + "/" + DRI->f_fvel;
    if ( !checkFile(tmp) )
    {
      Hostonly_ printf("\n\tError : File open '%s'\n", tmp.c_str());
      Exit(0);
    }
    F.readVelocity(fp, tmp, d_size, guide, d_wk, d_wo, step, time, v00, C.Unit.File, C.RefVelocity, flop, gs, true, i_dummy, f_dummy);
    
    if (C.Unit.File == DIMENSIONAL) time /= C.Tscale;
    
    if ( (step != Session_StartStep) || (time != Session_StartTime) )
    {
      Hostonly_ printf     ("\n\tTime stamp is different between files\n");
      Hostonly_ fprintf(fp, "\n\tTime stamp is different between files\n");
      Exit(0);
    }
  }
  else
  {
    // 受信
    unsigned recvsize=(d_size[0]+2*guide)*(d_size[1]+2*guide)*(d_size[2]+2*guide)*3;
    for(unsigned i=0;i<recvsize;i++) d_wk[i] = 0.0;
    if (paraMngr->Irecv( &time, 1, recv_rank, &mpi_req ) != CPM_SUCCESS) Exit(0);
    if (paraMngr->Wait( &mpi_req ) != MPI_SUCCESS) Exit(0);
    if (paraMngr->Irecv( &step, 1, recv_rank, &mpi_req ) != CPM_SUCCESS) Exit(0);
    if (paraMngr->Wait( &mpi_req ) != MPI_SUCCESS) Exit(0);
    if (paraMngr->Irecv( d_wk, recvsize, recv_rank, &mpi_req ) != CPM_SUCCESS) Exit(0);
    if (paraMngr->Wait( &mpi_req ) != MPI_SUCCESS) Exit(0);
  }
  SetOverlap(d_vf,d_wk,3,guide,head,size,DRI->overlap_head,DRI->overlap_tail,
             DRI->read_file_voxel_head,DRI->read_file_voxel_size);
  
}


// #################################################################
void FFV::ReadOverlap_Temperature(FILE* fp, double& flop, DifferentRestartInfo* DRI, DfiInfo* DI, REAL_TYPE* d_wk,
                                  int* rank_list, int recv_rank, int* assign, int nassign)
{
  MPI_Request mpi_req;
  double time;
  unsigned step;
  string tmp;
  
  std::string fmt(C.file_fmt_ext);
  
  // 出力モード
  bool mio = (bool)C.FIO.IOmode;
  
  // set rank_dir
  int len = 6;
  char* buff = new char[len];
  memset(buff, 0, sizeof(char)*len);
  sprintf(buff, "%06d", myRank);
  std::string rank_dir(buff);
  if ( buff ) delete [] buff;
  
  //
  REAL_TYPE klv = ( C.Unit.Temp == Unit_KELVIN ) ? 0.0 : KELVIN;
  
  // ガイド出力
  int gs = C.GuideOut;
  
  // dummy
  unsigned i_dummy=0;
  double f_dummy=0.0;
  
  // 読み込み領域の格子の分割サイズ
  int d_size[3];
  
  // ファイルの送信
  for(int k=0;k<numProc;k++){
    if(rank_list[k]==-1) continue;
    for(int j=0;j<nassign;j++){
      if(assign[j]==-1) continue;
      
      if(assign[j]==rank_list[k]){ // 自身のランクにステージングされたファイルに一致するものがあれば
        int send_rank=k;
        int read_rank=rank_list[k];
        
        if(send_rank==myRank) continue; // 自身のランクであれば送る必要なし
        
        // 読み込み領域の格子の分割サイズ
        for(int i=0;i<3;i++){
          d_size[i] = DI[0].Node[read_rank].VoxelSize[i];
        }
        
        double check = (double)(d_size[0]+2*guide)*(double)(d_size[1]+2*guide)*(double)(d_size[2]+2*guide);
        if(check>(double)INT_MAX){
          printf("\tsize error : sendsize>INT_MAX\n");
          Exit(0);
        }
        
        unsigned sendsize=(d_size[0]+2*guide)*(d_size[1]+2*guide)*(d_size[2]+2*guide);
        for(unsigned i=0;i<sendsize;i++) d_wk[i] = 0.0;
        
        // ファイルの読み込み
        tmp = DFI.GenerateFileName(C.f_different_nproc_temperature, fmt, C.Restart_step, read_rank, mio);
        tmp = C.f_different_nproc_dir_prefix + rank_dir + "/" + tmp;
        if ( !checkFile(tmp) )
        {
          printf("\n\tError : File open '%s'\n", tmp.c_str());
          Exit(0);
        }
        F.readTemperature(fp, tmp, d_size, guide, d_wk, step, time, C.Unit.File, C.BaseTemp, C.DiffTemp, klv, flop, gs, true, i_dummy, f_dummy);
        
        if (C.Unit.File == DIMENSIONAL) time /= C.Tscale;
        
        if ( (step != Session_StartStep) || (time != Session_StartTime) )
        {
          Hostonly_ printf     ("\n\tTime stamp is different between files\n");
          Hostonly_ fprintf(fp, "\n\tTime stamp is different between files\n");
          Exit(0);
        }
        
        // 送信
        if (paraMngr->Isend( &time, 1, send_rank, &mpi_req ) != CPM_SUCCESS) Exit(0);
        if (paraMngr->Wait( &mpi_req ) != MPI_SUCCESS) Exit(0);
        if (paraMngr->Isend( &step, 1, send_rank, &mpi_req ) != CPM_SUCCESS) Exit(0);
        if (paraMngr->Wait( &mpi_req ) != MPI_SUCCESS) Exit(0);
        if (paraMngr->Isend( d_wk, sendsize, send_rank, &mpi_req ) != CPM_SUCCESS) Exit(0);
        if (paraMngr->Wait( &mpi_req ) != MPI_SUCCESS) Exit(0);
      }
    }
  }
  
  // 送信だけしてリターン
  if(DRI->rank==-1) return;
  
  // ファイルの読み込み（受信）
  
  // 読み込み領域の格子の分割サイズ
  for(int i=0;i<3;i++){
    d_size[i] = DRI->read_file_voxel_size[i];
  }
  
  // 温度のオーバーラップ分を移す
  if(recv_rank==myRank)
  {
    // ファイルから読み込み
    tmp = C.f_different_nproc_dir_prefix + rank_dir + "/" + DRI->f_temp;
    if ( !checkFile(tmp) )
    {
      Hostonly_ printf("\n\tError : File open '%s'\n", tmp.c_str());
      Exit(0);
    }
    F.readTemperature(fp, tmp, d_size, guide, d_wk, step, time, C.Unit.File, C.BaseTemp, C.DiffTemp, klv, flop, gs, true, i_dummy, f_dummy);
    
    if (C.Unit.File == DIMENSIONAL) time /= C.Tscale;
    
    if ( (step != Session_StartStep) || (time != Session_StartTime) )
    {
      Hostonly_ printf     ("\n\tTime stamp is different between files\n");
      Hostonly_ fprintf(fp, "\n\tTime stamp is different between files\n");
      Exit(0);
    }
  }
  else
  {
    // 受信
    unsigned recvsize=(d_size[0]+2*guide)*(d_size[1]+2*guide)*(d_size[2]+2*guide);
    for(unsigned i=0;i<recvsize;i++) d_wk[i] = 0.0;
    if (paraMngr->Irecv( &time, 1, recv_rank, &mpi_req ) != CPM_SUCCESS) Exit(0);
    if (paraMngr->Wait( &mpi_req ) != MPI_SUCCESS) Exit(0);
    if (paraMngr->Irecv( &step, 1, recv_rank, &mpi_req ) != CPM_SUCCESS) Exit(0);
    if (paraMngr->Wait( &mpi_req ) != MPI_SUCCESS) Exit(0);
    if (paraMngr->Irecv( d_wk, recvsize, recv_rank, &mpi_req ) != CPM_SUCCESS) Exit(0);
    if (paraMngr->Wait( &mpi_req ) != MPI_SUCCESS) Exit(0);
  }
  SetOverlap(d_t,d_wk,1,guide,head,size,DRI->overlap_head,DRI->overlap_tail,
             DRI->read_file_voxel_head,DRI->read_file_voxel_size);
  
}
