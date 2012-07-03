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
 * @file   ffv_Restart.C
 * @brief  FFV Class
 * @author kero
 */

#include "ffv.h"


// ファイルのオープンチェック
bool FFV::checkFile(std::string fname)
{
  ifstream ifs(fname.c_str(), ios::in | ios::binary);
  if (!ifs) {
    return false;
  }
  ifs.close();
  
  return true;
}



// 2倍密格子の領域開始インデクス番号から、その領域が属する粗格子計算結果ファイル名と、その計算結果ファイルの開始インデクス番号を取得する
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
  char tmp[10]; // 10 digit
  memset(tmp, 0, sizeof(char)*10);
  sprintf(tmp, "%010d", m_step);
  std::string step(tmp);
  
	// dfiファイルを開いて
	ifstream ifs( coarse_dfi_fname.c_str() );
	if( !ifs ) return false;
	
	// 粗格子ijkが含まれるランクは？
  std::string buf;
	int rank = -1;
	int hi, hj, hk, ti, tj, tk;
  
  // サブドメインの分割数（粗格子）
  int Ci, Cj, Ck;
  
  
	while( getline(ifs, buf) ) {
    
    //if( buf.find("\"WorldNodeNum\"",0) != string::npos ) {
    //  np = get_intval( buf );	
    //}
    
		if( buf.find("\"GroupID\"",0) != string::npos ) {
			rank = get_intval( buf );
      
			while( getline(ifs, buf) ) {
        if( buf.find("\"VoxelSize\"",0) != string::npos ) {
					getline(ifs, buf);  Ci = get_intval( buf );	
					getline(ifs, buf);  Cj = get_intval( buf );	
					getline(ifs, buf);  Ck = get_intval( buf );	
				}
				if( buf.find("\"HeadIndex\"",0) != string::npos ) {
					getline(ifs, buf);  hi = get_intval( buf );	
					getline(ifs, buf);  hj = get_intval( buf );	
					getline(ifs, buf);  hk = get_intval( buf );	
				}
				if( buf.find("\"TailIndex\"",0) != string::npos ) {
					getline(ifs, buf);  ti = get_intval( buf );	
					getline(ifs, buf);  tj = get_intval( buf );	
					getline(ifs, buf);  tk = get_intval( buf );	
					break;
				}
			}
			if( i0>=hi && i0<=ti && j0>=hj && j0<=tj && k0>=hk && k0<=tk ) {
				// found!
				break;
			}
		}
	}
	if( rank == -1 ) return false;
  
  // このセッションの並列数
  //int mm = paraMngr->GetNumRank();
  //if ( np != mm ) {
  //  Hostonly_ printf("Error : The number of nodes in between previous[%d] and this[%d] session is different.\n", np, mm);
  //  Exit(0);
  //}
  
	// id=rankで、coarse_prefixをファイル名に含むsphファイルを探す
  std::string fname = "";
  std::string target = "";
	char id[32];
	sprintf(id, "id=\"%d\"", rank);
  
	while( getline(ifs, buf) ) {
		if( buf.find("\"FileName\"",0) != std::string::npos && buf.find(id,0) != string::npos ) {
			fname = get_strval( buf );
			if( (fname.find(coarse_prefix,0) != std::string::npos) && (fname.find(step,0) != std::string::npos) ) {
        target = fname;
        break;
			}
		}
	}
  
  if( target.empty() ) return false;
  
  
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



// ファイルから値をとりだす（整数）
int FFV::get_intval( string& buffer )
{
	int s = buffer.find( "value=\"", 0 ) + 7;
	int e = buffer.find( "\"", s );
	return atoi( buffer.substr( s, e-s ).c_str() );
}


// ファイルから値をとりだす（文字列）
string FFV::get_strval( string& buffer )
{
	int s = buffer.find( "value=\"", 0 ) + 7;
	int e = buffer.find( "\"", s );
  string result = buffer.substr( s, e-s );
	return result;
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
  
  fb_interp_coarse_s_(d_p, size, &guide, d_r_p, st, bk);
  fb_interp_coarse_v_(d_v, size, &guide, d_r_v, st, bk);
  
  if ( C.isHeatProblem() ) {
    fb_interp_coarse_s_(d_t, size, &guide, d_r_t, st, bk);
  }
  
}



// リスタートプロセス
void FFV::Restart()
{
  double flop_task;
  
  TIMING_start(tm_restart);
  
  if ( C.Start == initial_start) { // 初期スタートのステップ，時間を設定する
    
    Base_step = Current_step = Session_step = Total_step = 0;
    Base_time = Current_time = Session_time = Total_time = 0.0;
    
    // V00の値のセット．モードがONの場合はV00[0]=1.0に設定，そうでなければtmに応じた値
    if ( C.CheckParam == ON ) RF.setV00(Total_time, true);
    else                      RF.setV00(Total_time);
    
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
  else if ( C.Start == coarse_restart) // 粗い格子からのリスタート
  {
    Hostonly_ fprintf(stdout, "\t>> Restart from Previous Results on Coarse Mesh\n\n");
    Hostonly_ fprintf(fp, "\t>> Restart from Previous Results on Coarse Mesh\n\n");
    
    
    // 粗い格子のファイルをロードし、内挿処理を行う
    flop_task = 0.0;
    Restart_coarse(fp, flop_task);
    
    Hostonly_ fprintf(stdout,"\n");
    Hostonly_ fprintf(fp,"\n");
  }
  
  TIMING_stop(tm_restart);
}



// リスタートの最大値と最小値の表示
void FFV::Restart_display_minmax(FILE* fp, double& flop)
{
  if ( (C.Start == restart) || (C.Start == coarse_restart) ) {
    
    Hostonly_ fprintf(mp, "\tNon-dimensional value\n");
    Hostonly_ fprintf(fp, "\tNon-dimensional value\n");
    REAL_TYPE f_min, f_max, min_tmp, max_tmp, fpct;
    
    fpct = (REAL_TYPE)flop;
    
    // Velocity
    fb_minmax_v_ (&f_min, &f_max, size, guide, v00, d_v, &fpct); // allreduceすること
    
    if ( numProc > 1 ) {
      min_tmp = f_min;
      if( paraMngr->Allreduce(&min_tmp, &f_min, 1, MPI_MIN) != CPM_SUCCESS ) Exit(0);
      
      max_tmp = f_max;
      if( paraMngr->Allreduce(&max_tmp, &f_max, 1, SKL_MAX) != CPM_SUCCESS ) Exit(0);
    }
    
    Hostonly_ fprintf(stdout, "\t\tV : min=%13.6e / max=%13.6e\n", f_min, f_max);
    Hostonly_ fprintf(fp, "\t\tV : min=%13.6e / max=%13.6e\n", f_min, f_max);
    
    
    // Pressure
    fb_minmax_s_ (&f_min, &f_max, size, guide, d_p, &fpct);
    
    if ( numProc > 1 ) {
      min_tmp = f_min;
      if( !paraMngr->Allreduce(&min_tmp, &f_min, 1, MPI_MIN) != CPM_SUCCESS ) Exit(0);
      
      max_tmp = f_max;
      if( !paraMngr->Allreduce(&max_tmp, &f_max, 1, MPI_MAX) != CPM_SUCCESS ) Exit(0);
    }
    
    Hostonly_ fprintf(stdout, "\t\tP : min=%13.6e / max=%13.6e\n", f_min, f_max);
    Hostonly_ fprintf(fp, "\t\tP : min=%13.6e / max=%13.6e\n", f_min, f_max);
    
    // temperature
    if ( C.isHeatProblem() ) {
      fb_minmax_s_ (&f_min, &f_max, size, guide, d_t, &fpct);
      
      if ( numProc > 1 ) {
        min_tmp = f_min;
        if( !paraMngr->Allreduce(&min_tmp, &f_min, 1, MPI_MIN) != CPM_SUCCESS ) Exit(0);
        
        max_tmp = f_max;
        if( !paraMngr->Allreduce(&max_tmp, &f_max, 1, MPI_MAX) != CPM_SUCCESS ) Exit(0);
      }
      
      Hostonly_ fprintf(stdout, "\t\tT : min=%13.6e / max=%13.6e\n", f_min, f_max);
      Hostonly_ fprintf(fp, "\t\tT : min=%13.6e / max=%13.6e\n", f_min, f_max);
    }
    
    flop = (double)fpct;
	}
}


// リスタート時の瞬時値ファイル読み込み
void FFV::Restart_std(FILE* fp, double& flop)
{
  double time;
  unsigned step;
  string tmp;
  
  // 圧力の瞬時値　
  REAL_TYPE bp = ( C.Unit.Prs == Unit_Absolute ) ? C.BasePrs : 0.0;
  
  // ガイド出力
  int gs = C.GuideOut;
  
  // dummy
  unsigned i_dummy=0;
  double f_dummy=0.0;
  
  tmp = DFI.Generate_FileName(C.f_Pressure, C.Restart_step, myRank, (bool)C.FIO.IO_Input);
  
  if ( !checkFile(tmp) ) {
    Hostonly_ printf("\n\tError : File open '%s'\n", tmp.c_str());
    Exit(0);
  }
  
  F.readPressure(fp, tmp, size, guide, d_p, step, time, C.Unit.File, bp, C.RefDensity, C.RefVelocity, flop, gs, true, i_dummy, f_dummy);
  
  // ここでタイムスタンプを得る
  if (C.Unit.File == DIMENSIONAL) time /= C.Tscale;
  Base_step = step;
  Base_time = time;
  
  // v00[]に値をセット
  copyV00fromRF(Base_time);
  
  
  // Instantaneous Velocity fields
  tmp = DFI.Generate_FileName(C.f_Velocity, C.Restart_step, myRank, (bool)C.FIO.IO_Input);
  
  if ( !checkFile(tmp) ) {
    Hostonly_ printf("\n\tError : File open '%s'\n", tmp.c_str());
    Exit(0);
  }
  
  F.readVelocity(fp, tmp, size, guide, d_v, step, time, v00, C.Unit.File, C.RefVelocity, flop, gs, true, i_dummy, f_dummy);
  
  if (C.Unit.File == DIMENSIONAL) time /= (double)C.Tscale;
  
  if ( (step != Base_step) || (time != Base_time) ) 
  {
    Hostonly_ printf     ("\n\tTime stamp is different between files\n");
    Hostonly_ fprintf(fp, "\n\tTime stamp is different between files\n");
    Exit(0);
  }
  
  // Instantaneous Temperature fields
  if ( C.isHeatProblem() ) {
    
    REAL_TYPE klv = ( C.Unit.Temp == Unit_KELVIN ) ? 0.0 : KELVIN;
    
    tmp = DFI.Generate_FileName(C.f_Temperature, C.Restart_step, myRank, (bool)C.FIO.IO_Input);
    
    if ( !checkFile(tmp) ) {
      Hostonly_ printf("\n\tError : File open '%s'\n", tmp.c_str());
      Exit(0);
    }
    F.readTemperature(fp, tmp, size, guide, d_t, step, time, C.Unit.File, C.BaseTemp, C.DiffTemp, klv, flop, gs, true, i_dummy, f_dummy);
    
    if (C.Unit.File == DIMENSIONAL) time /= (double)C.Tscale;
    
    if ( (step != Base_step) || (time != Base_time) ) 
    {
      Hostonly_ printf     ("\n\tTime stamp is different between files\n");
      Hostonly_ fprintf(fp, "\n\tTime stamp is different between files\n");
      Exit(0);
    }
  }
  
}



// リスタート時の平均値ファイル読み込み
void FFV::Restart_avrerage (FILE* fp, double& flop)
{
  std::string tmp;
  
  unsigned step = Base_step;
  double time = Base_time;

  
  // ガイド出力
  int gs = C.GuideOut;
  
  if ( C.Interval[Interval_Manager::tg_avstart].isStep() ) {
    if ( step > C.Interval[Interval_Manager::tg_avstart].getIntervalStep() ) {
      Hostonly_ printf     ("\tRestart from Previous Calculation Results of averaged field\n");
      Hostonly_ fprintf(fp, "\tRestart from Previous Calculation Results of averaged field\n");
      Hostonly_ printf     ("\tStep : base=%d current=%d total=%d\n", step, Current_step, Total_step);
      Hostonly_ fprintf(fp, "\tStep : base=%d current=%d total=%d\n", step, Current_step, Total_step);
    }
    else {
      return;
    }
  }
  else {
    if ( time > C.Interval[Interval_Manager::tg_avstart].getIntervalTime() ) {
      Hostonly_ printf     ("\tRestart from Previous Calculation Results of averaged field\n");
      Hostonly_ fprintf(fp, "\tRestart from Previous Calculation Results of averaged field\n");
      Hostonly_ printf     ("\tTime : base=%e[sec.]/%e[-] current=%e[-] total=%e[-]\n", time*C.Tscale, time, Current_time, Total_time);
      Hostonly_ fprintf(fp, "\tTime : base=%e[sec.]/%e[-] current=%e[-] total=%e[-]\n", time*C.Tscale, time, Current_time, Total_time);
    }
    else {
      return;
    }
  }
  
  unsigned step_avr = 0;
  double time_avr = 0.0;
  
  // Pressure
  REAL_TYPE bp = ( C.Unit.Prs == Unit_Absolute ) ? C.BasePrs : 0.0;
  
  tmp = DFI.Generate_FileName(C.f_AvrPressure, C.Restart_step, myRank, (bool)C.FIO.IO_Input);
  if ( !checkFile(tmp) ) {
    Hostonly_ printf("\n\tError : File open '%s'\n", tmp.c_str());
    Exit(0);
  }
  F.readPressure(fp, tmp, size, guide, d_ap, step, time, C.Unit.File, bp, C.RefDensity, C.RefVelocity, flop, gs, false, step_avr, time_avr);
  
  if ( (step != Base_step) || (time != Base_time) ) {
    Hostonly_ printf     ("\n\tTime stamp is different between instantaneous and averaged files\n");
    Hostonly_ fprintf(fp, "\n\tTime stamp is different between instantaneous and averaged files\n");
    Exit(0);
  }
  
  if (C.Unit.File == DIMENSIONAL) time /= (double)C.Tscale;
  
  Total_step_avr = step_avr;
  Total_time_avr = time_avr;
  
  
  // Velocity
  tmp = DFI.Generate_FileName(C.f_AvrVelocity, C.Restart_step, myRank, (bool)C.FIO.IO_Input);
  if ( !checkFile(tmp) ) {
    Hostonly_ printf("\n\tError : File open '%s'\n", tmp.c_str());
    Exit(0);
  }
  F.readVelocity(fp, tmp, size, guide, d_av, step, time, v00, C.Unit.File, C.RefVelocity, flop, gs, false, step_avr, time_avr);
  
  if (C.Unit.File == DIMENSIONAL) time /= (double)C.Tscale;
  
  if ( (step_avr != Total_step_avr) || (time_avr != Total_time_avr) ) { // 圧力とちがう場合
    Hostonly_ printf     ("\n\tTime stamp is different between files\n");
    Hostonly_ fprintf(fp, "\n\tTime stamp is different between files\n");
    Exit(0);
  }
  
  
  if ( C.isHeatProblem() ) {
    
    REAL_TYPE klv = ( C.Unit.Temp == Unit_KELVIN ) ? 0.0 : KELVIN;
    
    tmp = DFI.Generate_FileName(C.f_AvrTemperature, C.Restart_step, myRank, (bool)C.FIO.IO_Input);
    if ( !checkFile(tmp) ) {
      Hostonly_ printf("\n\tError : File open '%s'\n", tmp.c_str());
      Exit(0);
    }
    F.readTemperature(fp, tmp, size, guide, d_at, step, time, C.Unit.File, C.BaseTemp, C.DiffTemp, klv, flop, gs, false, step_avr, time_avr);
    
    if (C.Unit.File == DIMENSIONAL) time /= (double)C.Tscale;
    
    if ( (step_avr != Total_step_avr) || (time_avr != Total_time_avr) ) {
      Hostonly_ printf     ("\n\tTime stamp is different between files\n");
      Hostonly_ fprintf(fp, "\n\tTime stamp is different between files\n");
      Exit(0);
    }
  }
}




// 粗い格子を用いたリスタート
void FFV::Restart_coarse(FILE* fp, double& flop)
{
  std::string f_prs;
  std::string f_vel;
  std::string f_temp;
  
  // 粗格子の分割サイズ
  int r_size[3];
  
  // 粗格子サブドメインにおける読み込み開始インデクス
  int crs[3];
  
  // 粗格子サブドメインに含まれる密格子サブドメインの数
  int num_block[3];
  
  //並列時には各ランクに必要なファイル名と開始インデクスを取得
  if ( C.FIO.IO_Input == IO_DISTRIBUTE ) {
    
    int i, j, k;  // 密格子のグローバル開始インデクス
    i = head[0];
    j = head[1];
    k = head[2];
    
    // crs_i, _j, _kには同じ値が入る 
    getCoarseResult(i, j, k, C.f_Coarse_dfi_prs, C.f_Coarse_pressure, C.Restart_step, f_prs, r_size, crs, num_block);
    
    getCoarseResult(i, j, k, C.f_Coarse_dfi_vel, C.f_Coarse_velocity, C.Restart_step, f_vel, r_size, crs, num_block);
    
    if ( C.isHeatProblem() ) {
      getCoarseResult(i, j, k, C.f_Coarse_dfi_temp, C.f_Coarse_temperature, C.Restart_step, f_temp, r_size, crs, num_block);
    }
  }
  else {
    crs[0] = 1;
    crs[1] = 1;
    crs[2] = 1;
    f_prs = DFI.Generate_FileName(C.f_Coarse_pressure, C.Restart_step, myRank);
    f_vel = DFI.Generate_FileName(C.f_Coarse_velocity, C.Restart_step, myRank);
    if ( C.isHeatProblem() )
      f_temp= DFI.Generate_FileName(C.f_Coarse_temperature, C.Restart_step, myRank);
  }
  
  // テンポラリのファイルロード
  allocArray_CoarseMesh(r_size, flop);
  
  
  // ガイド出力
  int gs = C.GuideOut;
  
  
  // dummy
  unsigned i_dummy=0;
  double f_dummy=0.0;
  
  unsigned step;
  double time;
  
  // 圧力の瞬時値　ここでタイムスタンプを得る
  REAL_TYPE bp = ( C.Unit.Prs == Unit_Absolute ) ? C.BasePrs : 0.0;
  
  if ( !checkFile(f_prs) ) {
    Hostonly_ printf("\n\tError : File open '%s'\n", f_prs.c_str());
    Exit(0);
  }
  
  F.readPressure(fp, f_prs, r_size, guide, d_r_p, step, time, C.Unit.File, bp, C.RefDensity, C.RefVelocity, flop, gs, true, i_dummy, f_dummy);
  
  if (C.Unit.File == DIMENSIONAL) time /= (double)C.Tscale;
  Base_step = step;
  Base_time = time;
  
  // v00[]に値をセット
  copyV00fromRF(Base_time);
  
  
  // Instantaneous Velocity fields
  if ( !checkFile(f_vel) ) {
    Hostonly_ printf("\n\tError : File open '%s'\n", f_vel.c_str());
    Exit(0);
  }
  F.readVelocity(fp, f_vel, r_size, guide, d_r_v, step, time, v00, C.Unit.File, C.RefVelocity, flop, gs, true, i_dummy, f_dummy);
  
  if (C.Unit.File == DIMENSIONAL) time /= (double)C.Tscale;
  
  if ( (step != Base_step) || (time != Base_time) )
  {
    Hostonly_ printf     ("\n\tTime stamp is different between files\n");
    Hostonly_ fprintf(fp, "\n\tTime stamp is different between files\n");
    Exit(0);
  }
  
  // Instantaneous Temperature fields
  if ( C.isHeatProblem() ) {
    
    REAL_TYPE klv = ( C.Unit.Temp == Unit_KELVIN ) ? 0.0 : KELVIN;
    
    if ( !checkFile(f_temp) ) {
      Hostonly_ printf("\n\tError : File open '%s'\n", f_temp.c_str());
      Exit(0);
    }
    F.readTemperature(fp, f_temp, r_size, guide, d_r_t, step, time, C.Unit.File, C.BaseTemp, C.DiffTemp, klv, flop, gs, true, i_dummy, f_dummy);
    
    if (C.Unit.File == DIMENSIONAL) time /= (double)C.Tscale;
    
    if ( (step != Base_step) || (time != Base_time) )
    {
      Hostonly_ printf     ("\n\tTime stamp is different between files\n");
      Hostonly_ fprintf(fp, "\n\tTime stamp is different between files\n");
      Exit(0);
    }
  }
  
  // 同期
  if ( paraMngr->BndCommV3DEx(d_r_v, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
  if ( paraMngr->BndCommS3D(d_r_p, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
  
  if ( C.isHeatProblem() ) {
    if ( paraMngr->BndCommS3D(d_r_t, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
  }
  
  // 内挿処理
  Interpolation_from_coarse_initial(crs, num_block);
}



// 並列分散時のファイル名の管理を行う
void FFV::setDFI()
{
  // 並列時のみ
  if ( numProc > 1 ) {
    int* g_bbox_st = new int[3*numProc];
    int* g_bbox_ed = new int[3*numProc];
    
    int st_i, st_j ,st_k, ed_i, ed_j, ed_k;
    
    // head and tail index
    for (int n=0; n<numProc; n++){
      
      // Fortran index
      st_i = head[0];
      st_j = head[1];
      st_k = head[2];
      
      const int* tail = paraMngr->GetVoxelTailIndex();
      ed_i = tail[0] + 1;
      ed_j = tail[1] + 1;
      ed_k = tail[2] + 1;
      
      if ( (g_bbox_st[3*n+0] = st_i) < 1 ) Exit(0);
      if ( (g_bbox_st[3*n+1] = st_j) < 1 ) Exit(0);
      if ( (g_bbox_st[3*n+2] = st_k) < 1 ) Exit(0);
      if ( (g_bbox_ed[3*n+0] = ed_i) < 1 ) Exit(0);
      if ( (g_bbox_ed[3*n+1] = ed_j) < 1 ) Exit(0);
      if ( (g_bbox_ed[3*n+2] = ed_k) < 1 ) Exit(0);
    }
    
    // DFIクラスの初期化
    if ( !DFI.init(G_size, paraMngr->GetDivNum(), C.GuideOut, C.Start, g_bbox_st, g_bbox_ed) ) Exit(0);
    
    // host name
    //for (int n=0; n<numProc; n++){
    //  const char* host = para_mng->GetHostName(n, procGrp);
    //  if ( !host ) Exit(0);
    
    //  DFI.copy_hostname(host, n);
    //}
    
    if ( g_bbox_st ) {
      delete [] g_bbox_st; 
      g_bbox_st = NULL;
    }
    if ( g_bbox_ed ) {
      delete [] g_bbox_ed; 
      g_bbox_ed = NULL;
    }
  }
  
}

