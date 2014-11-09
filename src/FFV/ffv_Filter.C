//##################################################################################
//
// FFV-C : Frontflow / violet Cartesian
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
 * @file   ffv_Filter.C
 * @brief  FFV filter
 * @author aics
 */

#include "ffv.h"

#define SafeDelClass(_C) if(_C){delete _C; _C=NULL;}
#define SafeDelArray(_A) if(_A){delete[] _A; _A=NULL;}
#define SafeDelArray_(_A, _S) if(_A){delete[] _A; _A=NULL; _S=0;}


int FFV::FilterLoop()
{
  int ret = 1;
  
  printf("\n\n\tData sampling....\n\n");
  
  C.Interval[Control::tg_compute].printInfo("tg_compute");
  C.Interval[Control::tg_console].printInfo("tg_console");
  C.Interval[Control::tg_history].printInfo("tg_history");
  C.Interval[Control::tg_basic].printInfo("tg_basic");
  C.Interval[Control::tg_statistic].printInfo("tg_statistic");
  C.Interval[Control::tg_derived].printInfo("tg_derived");
  C.Interval[Control::tg_accelra].printInfo("tg_accelra");
  C.Interval[Control::tg_sampled].printInfo("tg_sampled");
  C.Interval[Control::tg_END].printInfo("tg_END");
  
  //FFVのデータから取得できれば、ここでハードコードする必要がない。
  int interval = 32;
  int i = 0;
  
  for(;;)
  {
    printf("\n\n\tData sampling i=%d\n\n", i);
    
    if ( FFV_TerminateCtrl::getTerminateFlag() )
    {
      break;// forced terminate
    }
    
    /////////////////////////////////// FilterLoop(i)
    //
    int loop_ret = FilterLoop(i);
    //
    /////////////////////////////////// FilterLoop(i)
    
    if ( loop_ret == 0 ) break;
    
    i += interval;
  }
  
  printf("\n\n\tClose MO, End FilterLoop()\n\n");
  
  // サンプリングファイルのクローズ
  MO.closeFile();
  
  if ( fp_b ) fclose(fp_b);  ///< 基本情報
  if ( fp_w ) fclose(fp_w);  ///< 壁面情報
  if ( fp_c ) fclose(fp_c);  ///< コンポーネント情報
  if ( fp_d ) fclose(fp_d);  ///< 流量収支情報
  if ( fp_i ) fclose(fp_i);  ///< 反復履歴情報
  if ( fp_f ) fclose(fp_f);  ///< 力の履歴情報
  
  if ( !stepPost() ) ret = -1;
  
  return ret;
}

int FFV::FilterLoop(unsigned int step)
{
  
}

/*
int FFV::FilterLoop(unsigned int step)
{
  int ret=1;
  
  {
    char work_path[1024]="";
    getcwd( work_path, sizeof(work_path) );
    printf("\n\twork_path = %s\n", work_path);
    
    
    if ( !F.checkOutFile() )
    {
      printf("\n\t\tdfi is NULL, return.\n");
      return 0;
    }
    
    //bool is_gathered = (MO.getOutputType()==MonitorList::GATHER ? true : false);
    
    int istep = step;
    
    bool mio = (numProc > 1 ? true : false); //並列判定フラグ（逐次or並列の判定用）
    
    for( int rank_id=0; rank_id<numProc; rank_id++ )
    {
      cdm_Rank rank_v;   if(DFI_OUT_VEL!=NULL)  rank_v   = DFI_OUT_VEL->GetcdmProcess()->RankList[rank_id];
      cdm_Rank rank_p;   if(DFI_OUT_PRS!=NULL)  rank_p   = DFI_OUT_PRS->GetcdmProcess()->RankList[rank_id];
      cdm_Rank rank_t;   if(DFI_OUT_TEMP!=NULL) rank_t   = DFI_OUT_TEMP->GetcdmProcess()->RankList[rank_id];
      cdm_Rank rank_vrt; if(DFI_OUT_VRT!=NULL)  rank_vrt = DFI_OUT_VRT->GetcdmProcess()->RankList[rank_id];
      ;
      
      bool b_alloc = false; //false の場合、FFVに既にある配列を使用する。
      
      //----------------------------------------------------------------
      //デフォルトでは、REAL_TYPE=float ,コンパイル時オプション-D_REAL_IS_DOUBLE_
      //を付与することで, REAL_TYPE=doubleになる
      REAL_TYPE *p_arr_v=NULL, *p_arr_p=NULL, *p_arr_t=NULL, *p_arr_vrt=NULL;
      int        n_arr_v=0,     n_arr_p=0,     n_arr_t=0,     n_arr_vrt=0;
      
      if( b_alloc == false )
      {
        p_arr_v  = d_v;     //dv--セルセンター速度, d_wo--入出力のバッファワーク
        p_arr_p  = d_p;     //dp--圧力
        p_arr_t  = d_ws;    //d_ws--反復中に固定のソース, d_ie--内部エネルギー
        p_arr_vrt= d_vrt;   //d_vrt--渦度ベクトル
      }
      
      //フィールドデータの読込み
      
      int rc1 = FilterGetArrayFromSph(DFI_OUT_VEL,  &rank_v,  istep, &p_arr_v,  &n_arr_v );
      if( rc1 != CDM::E_CDM_SUCCESS ){ p_arr_v=NULL; n_arr_v=0;}
      
      int rc2 = FilterGetArrayFromSph(DFI_OUT_PRS,  &rank_p,  istep, &p_arr_p,  &n_arr_p );
      if( rc2 != CDM::E_CDM_SUCCESS ){ p_arr_p=NULL; n_arr_p=0;}
      
      int rc3 = FilterGetArrayFromSph(DFI_OUT_TEMP, &rank_t,  istep, &p_arr_t,  &n_arr_t );
      if( rc3 != CDM::E_CDM_SUCCESS ){ p_arr_t=NULL; n_arr_t=0;}
      
      int rc4 = FilterGetArrayFromSph(DFI_OUT_FVEL, &rank_vrt, istep, &p_arr_vrt, &n_arr_vrt);
      if( rc4 != CDM::E_CDM_SUCCESS ){ p_arr_vrt=NULL; n_arr_vrt=0;}
      
      if( rc1!=CDM::E_CDM_SUCCESS && rc2!=CDM::E_CDM_SUCCESS )
      {
        ret = CDM::E_CDM_ERROR;
        return ret;
      }
      
      if( rc1!=CDM::E_CDM_SUCCESS && rc2!=CDM::E_CDM_SUCCESS && rc3!=CDM::E_CDM_SUCCESS && rc4!=CDM::E_CDM_SUCCESS )
      {
        ret = CDM::E_CDM_ERROR;
        return ret;
      }
      
      //ここで、明示的にサンプリング元となるデータ配列(REAL_TYPE)の登録必要
      //   v  速度変数配列       p   圧力変数配列
      //   t  温度変数配列       vrt 渦度変数配列
      MO.setDataPtrs(p_arr_v, p_arr_p, p_arr_t, p_arr_vrt);
      
      //ランク毎に、サンプリングを行う
      MO.sampling();
      
      if( b_alloc == true )
      {
        SafeDelArray_(p_arr_v, n_arr_p);
        SafeDelArray_(p_arr_p, n_arr_p);
        SafeDelArray_(p_arr_t, n_arr_t);
        SafeDelArray_(p_arr_vrt, n_arr_vrt);
      }
    }
    
    //getherされた結果を出力する。
    MO.print( step, 0.0 );
  }
  
  return ret;
}*/



/* ZRM
int FFV::FilterLoop(unsigned int step)
{
  int ret=1;
  
  {
    char work_path[1024]="";
    getcwd( work_path, sizeof(work_path) );
    printf("\n\twork_path = %s\n", work_path);
    
    
    if ( !F.checkOutFile() )
    {
      printf("\n\t\tdfi is NULL, return.\n");
      return 0;
    }
    
    //bool is_gathered = (MO.getOutputType()==MonitorList::GATHER ? true : false);
    
    int istep = step;
    
    bool mio = (numProc > 1 ? true : false); //並列判定フラグ（逐次or並列の判定用）
    
    for( int rank_id=0; rank_id<numProc; rank_id++ )
    {
      cdm_Rank rank_v;   if(DFI_OUT_VEL!=NULL)  rank_v   = DFI_OUT_VEL->GetcdmProcess()->RankList[rank_id];
      cdm_Rank rank_p;   if(DFI_OUT_PRS!=NULL)  rank_p   = DFI_OUT_PRS->GetcdmProcess()->RankList[rank_id];
      cdm_Rank rank_t;   if(DFI_OUT_TEMP!=NULL) rank_t   = DFI_OUT_TEMP->GetcdmProcess()->RankList[rank_id];
      cdm_Rank rank_vrt; if(DFI_OUT_VRT!=NULL)  rank_vrt = DFI_OUT_VRT->GetcdmProcess()->RankList[rank_id];
      ;
      
      bool b_alloc = false; //false の場合、FFVに既にある配列を使用する。
      
      //----------------------------------------------------------------
      //デフォルトでは、REAL_TYPE=float ,コンパイル時オプション-D_REAL_IS_DOUBLE_
      //を付与することで, REAL_TYPE=doubleになる
      REAL_TYPE *p_arr_v=NULL, *p_arr_p=NULL, *p_arr_t=NULL, *p_arr_vrt=NULL;
      int        n_arr_v=0,     n_arr_p=0,     n_arr_t=0,     n_arr_vrt=0;
      
      if( b_alloc == false )
      {
        p_arr_v  = d_v;     //dv--セルセンター速度, d_wo--入出力のバッファワーク
        p_arr_p  = d_p;     //dp--圧力
        p_arr_t  = d_ws;    //d_ws--反復中に固定のソース, d_ie--内部エネルギー
        p_arr_vrt= d_vrt;   //d_vrt--渦度ベクトル
      }
      
      //フィールドデータの読込み
      
      int rc1 = FilterGetArrayFromSph(DFI_OUT_VEL,  &rank_v,  istep, &p_arr_v,  &n_arr_v );
      if( rc1 != CDM::E_CDM_SUCCESS ){ p_arr_v=NULL; n_arr_v=0;}
      
      int rc2 = FilterGetArrayFromSph(DFI_OUT_PRS,  &rank_p,  istep, &p_arr_p,  &n_arr_p );
      if( rc2 != CDM::E_CDM_SUCCESS ){ p_arr_p=NULL; n_arr_p=0;}
      
      int rc3 = FilterGetArrayFromSph(DFI_OUT_TEMP, &rank_t,  istep, &p_arr_t,  &n_arr_t );
      if( rc3 != CDM::E_CDM_SUCCESS ){ p_arr_t=NULL; n_arr_t=0;}
      
      int rc4 = FilterGetArrayFromSph(DFI_OUT_FVEL, &rank_vrt, istep, &p_arr_vrt, &n_arr_vrt);
      if( rc4 != CDM::E_CDM_SUCCESS ){ p_arr_vrt=NULL; n_arr_vrt=0;}
      
      if( rc1!=CDM::E_CDM_SUCCESS && rc2!=CDM::E_CDM_SUCCESS )
      {
        ret = CDM::E_CDM_ERROR;
        return ret;
      }
      
      if( rc1!=CDM::E_CDM_SUCCESS && rc2!=CDM::E_CDM_SUCCESS && rc3!=CDM::E_CDM_SUCCESS && rc4!=CDM::E_CDM_SUCCESS )
      {
        ret = CDM::E_CDM_ERROR;
        return ret;
      }
      
      //ここで、明示的にサンプリング元となるデータ配列(REAL_TYPE)の登録必要
      //   v  速度変数配列       p   圧力変数配列
      //   t  温度変数配列       vrt 渦度変数配列
      MO.setDataPtrs(p_arr_v, p_arr_p, p_arr_t, p_arr_vrt);
      
      //ランク毎に、サンプリングを行う
      MO.sampling();
      
      if( b_alloc == true )
      {
        SafeDelArray_(p_arr_v, n_arr_p);
        SafeDelArray_(p_arr_p, n_arr_p);
        SafeDelArray_(p_arr_t, n_arr_t);
        SafeDelArray_(p_arr_vrt, n_arr_vrt);
      }
    }
    
    //getherされた結果を出力する。
    MO.print( step, 0.0 );
  }
  
  return ret;
}
*/
int FFV::FilterGetArrayFromSph( cdm_DFI *dfi, cdm_Rank *rank, int step, REAL_TYPE **pArray, int *nArray)
{
  
}
/*
int FFV::FilterGetArrayFromSph( cdm_DFI *dfi, cdm_Rank *rank, int step, REAL_TYPE **pArray, int *nArray)
{
  CDM::E_CDM_ERRORCODE ret = CDM::E_CDM_SUCCESS;
  
  if( dfi == NULL || rank == NULL || step < 0 || pArray == NULL || nArray == NULL )
  {
    ret = CDM::E_CDM_ERROR;
    return (int)ret;
  }
  if( rank->RankID == 0 && rank->HostName == "" )
  {
    ret = CDM::E_CDM_ERROR;
    return (int)ret;
  }
  
  int istep = step;
  int outGc = C.GuideOut;
  
  int read_sta[3],read_end[3];
  for(int k=0; k<3; k++) {
    read_sta[k] = rank->HeadIndex[k] - outGc;
    read_end[k] = rank->TailIndex[k] + outGc;
  }
  
  float r_time=0.0, f_dummy=0.0;
  unsigned i_dummy=0;
  
  //計算空間の定義
  const int *Gdiv = paraMngr->GetDivNum();
  const int *Gvox = paraMngr->GetLocalVoxelSize();
  
  //読込み配列のサイズ
  int  ncomp=dfi->GetNumComponent();
  size_t  size=(Gvox[0]+2*outGc)*(Gvox[1]+2*outGc)*(Gvox[2]+2*outGc);
  if( size <= 0 ) return (int)ret;
  
  bool b_alloc=false;
  REAL_TYPE *p_array = *pArray;
  if( p_array == NULL ){ p_array = new REAL_TYPE[size*ncomp]; b_alloc=true;}
  
  //フィールドデータの読込み
  bool b_donot_avr = true;
  
  ret = dfi->ReadData(    p_array,    ///<読込み先配列のポインタ
                      istep,      ///<読込みフィールドデータのステップ番号
                      outGc,      ///<計算空間の仮想セル数
                      Gvox,       ///<計算空間全体のボクセルサイズ
                      (int*)Gdiv, ///<領域分割数
                      read_sta,   ///<計算領域の開始位置
                      read_end,   ///<計算領域の終了位置
                      r_time,     ///<dfi から読込んだ時間
                      b_donot_avr,///<平均を読込まない
                      i_dummy,
                      f_dummy  );
  
  if( ret != CDM::E_CDM_SUCCESS )
  {
    if( b_alloc == true ){ SafeDelArray_(p_array, size);}
  }
  else
  {
    ret = CDM::E_CDM_SUCCESS;
  }
  
  *pArray = p_array;
  *nArray = size;
  
  
  return (int)ret;
}*/



int FFV::FilterInitialize(int argc, char **argv)
{
  
}
// #################################################################
/* @brief フィルタ処理の初期化
 * @param [in] argc  main関数の引数の個数
 * @param [in] argv  main関数の引数リスト
 *
int FFV::FilterInitialize(int argc, char **argv)
{
  double TotalMemory   = 0.0;  ///< 計算に必要なメモリ量（ローカル）
  double PrepMemory    = 0.0;  ///< 初期化に必要なメモリ量（ローカル）
  double G_TotalMemory = 0.0;  ///< 計算に必要なメモリ量（グローバル）
  double G_PrepMemory  = 0.0;  ///< 初期化に必要なメモリ量（グローバル）
  double tmp_memory    = 0.0;  ///< 計算に必要なメモリ量（グローバル）
  double flop_task     = 0.0;  ///< flops計算用
  
  
  
  // cpm_ParaManagerのポインタをセット
  C.importCPM(paraMngr);
  V.importCPM(paraMngr);
  B.importCPM(paraMngr);
  BC.importCPM(paraMngr);
  MO.importCPM(paraMngr);
  
  // 入力ファイルの指定
  std::string input_file = argv[2];
  
  
  // ffvのパラメータローダのインスタンス生成
  TextParser tp_ffv;
  
  
  // パラメータのロードと保持
  {
    int ierror=0;
    
    if ( (ierror = tp_ffv.read(input_file)) != TP_NO_ERROR )
    {
      Hostonly_ stamped_printf("\tError at reading '%s' file : %d\n", input_file.c_str(), ierror);
      Exit(0);
    }
  }
  
  
  // TextParserクラスのポインタを各クラスに渡す
  C.importTP(&tp_ffv);
  B.importTP(&tp_ffv);
  M.importTP(&tp_ffv);
  MO.importTP(&tp_ffv);
  
  
  // 固定パラメータ
  fixedParameters();
  
  
  
  
  // ------------------------------------
  FILE* fp = NULL;
  
  // condition fileのオープン
  Hostonly_
  {
    if ( !(fp=fopen("filter-condition.txt", "w")) )
    {
      stamped_printf("\tSorry, can't open 'filter-condition.txt' file. Write failed.\n");
      return -1;
    }
  }
  
  // メッセージ表示
  Hostonly_
  {
    FBUtility::printVersion(fp,     "FFV Filter ", FFVC_VERSION_NO);
    FBUtility::printVersion(stdout, "FFV Filter ", FFVC_VERSION_NO);
  }
  
  
  // 反復制御クラスのインスタンス
  //C.getIteration();
  
  
  // 流体の解法アルゴリズムを取得
  C.getSolvingMethod4Flow();
  
  
  // 線形ソルバーの特定
  identifyLinearSolver(&tp_ffv);
  
  
  // 計算モデルの入力ソース情報を取得
  C.getGeometryModel();
  
  
  // Intrinsic classの同定
  identifyExample(fp);
  
  
  //CAUTION: setupDomain() --> DomainInitialize(tpf); --> getDomainInfo() TextParserから
  //         /DomainInfo/ActiveSubDomainFile の値を取得しています、値がstr.empty()==true
  //         EXEC_MODE = ffvc_solver;
  
  //EXEC_MODE 値を保存する
  int prev_EXEC_MODE = EXEC_MODE;
  
  //元の EXEC_MODE 値に戻る
  EXEC_MODE = prev_EXEC_MODE;
  
  
  // パラメータの取得と計算領域の初期化，並列モードを返す
  std::string str_para = setupDomain(&tp_ffv);
  
  
  // mat[], cmp[]の作成
  createTable(fp);
  
  
  // 媒質情報をパラメータファイルから読み込み，媒質リストを作成する
  setMediumList(fp);
  
  
  V.setControlVars(C.NoCompo, C.NoMedium, Ex);
  
  
  // CompoListの設定，外部境界条件の読み込み保持
  setBCinfo();
  
  
  
  // タイミング測定の初期化
  if ( C.Mode.Profiling != OFF )
  {
    ModeTiming = ON;
    TIMING__ PM.initialize( PM_NUM_MAX );
    TIMING__ PM.setRankInfo( paraMngr->GetMyRankID() );
    TIMING__ PM.setParallelMode(str_para, C.num_thread, C.num_process);
    set_timing_label();
  }
  
  
  // タイミング測定開始
  TIMING_start(tm_init_sct);
  
  
  
  // 前処理に用いるデータクラスのアロケート -----------------------------------------------------

  setArraySize();
  allocArray_Prep(PrepMemory, TotalMemory);
  
  
  
  TIMING_start(tm_voxel_prep_sct);
  
  
  
  // 各問題に応じてモデルを設定 >> Polylib + Cutlib
  // 外部境界面およびガイドセルのカットとIDの処理
  setModel(PrepMemory, TotalMemory, fp);
  
  
  
  // 定義点上に交点がある場合の処理 >> カットするポリゴンのエントリ番号でフィルする
  unsigned long fill_cut = V.modifyCutOnCellCenter(d_bid, d_cut, C.FillID);
  
  Hostonly_
  {
    printf(    "\tCut on cell center = %16ld\n\n", fill_cut);
    fprintf(fp,"\tCut on cell center = %16ld\n\n", fill_cut);
  }
  
  
  
  // 領域情報の表示
  Hostonly_
  {
    printf("\n----------\n");
    printf("\n\t>> Global Domain Information\n\n");
    C.printGlobalDomain(stdout, G_size, G_origin, G_region, pitch);
    
    fprintf(fp,"\n----------\n");
    fprintf(fp,"\n\t>> Global Domain Information\n\n");
    C.printGlobalDomain(fp, G_size, G_origin, G_region, pitch);
  }
  
  
  // メモリ消費量の情報を表示
  Hostonly_
  {
    printf(    "\n----------\n\n");
    fprintf(fp,"\n----------\n\n");
  }
  G_PrepMemory = PrepMemory;
  
  displayMemoryInfo(fp, G_PrepMemory, PrepMemory, "Preprocessor");
  
  
  
  
  // サンプリング準備
  //
  // In setMonitorList(), MO.getMonitor(&C, cmp); will be called.
  // All defined sampling objects were added into MonitorList.
  // Output file(s) were opened.
  //
  setMonitorList();
  
  
  
  // Fill
  Hostonly_
  {
    printf(    "\n----------\n\n");
    fprintf(fp,"\n----------\n\n");
    printf(    "\t>> Fill\n\n");
    fprintf(fp,"\t>> Fill\n\n");
  }
  
  fill(fp);
  
  
  // 全周カットのあるセルを固体セルIDで埋める > fill()で埋められているので不要．
  // V.replaceIsolatedFcell(d_bcd, C.FillID, d_bid);
  
  
  
  // ∆tの決め方とKindOfSolverの組み合わせで無効なものをはねる
  if ( !DT.chkDtSelect() )
  {
    Hostonly_
    {
      printf(    "\tCombination of specified 'TimeIncrement' and 'KindOfSolver' is not permitted.\n");
      fprintf(fp,"\tCombination of specified 'TimeIncrement' and 'KindOfSolver' is not permitted.\n");
    }
    return -1;
  }
  
  
  // パラメータファイルから得られた内部BCコンポーネント数を表示
  Hostonly_
  {
    printf("\n----------\n\n");
    printf("\t>> Components\n\n");
    C.printNoCompo(stdout);
    printf("\n"); fflush(stdout);
    
    fprintf(fp,"\n----------\n\n");
    fprintf(fp,"\t>> Components\n\n");
    C.printNoCompo(fp);
    fprintf(fp,"\n"); fflush(fp);
  }
  
  
  
  
  // HEX/FANコンポーネントの形状情報からBboxと体積率を計算
  if ( C.EnsCompo.fraction )
  {
    TIMING_start(tm_init_alloc);
    allocArray_CompoVF(PrepMemory, TotalMemory);
    TIMING_stop(tm_init_alloc);
    
    setComponentVF();
  }
  
  
  
  // 内部周期境界の場合のガイドセルのコピー処理
  V.adjMediumPrdcInner(d_bcd, cmp);
  
  
  // 媒質数とKindOfSolverの整合性をチェックする
  if ( !chkMediumConsistency() )
  {
    Hostonly_
    {
      stamped_printf(    "\tchkMediumConsistency()\n");
      stamped_fprintf(fp,"\tchkMediumConsistency()\n");
    }
    return -1;
  }
  
  
  
  
  // BCIndexにビット情報をエンコードとコンポーネントインデクスの再構築
  encodeBCindex(fp);
  
  
  
  // Polygonモニタの数をcmp[]にセット
  if ( C.SamplingMode )
  {
    MO.setMonitorNpoint(cmp, C.NoCompo);
  }
  
  
  
  // 体積力を使う場合のコンポーネント配列の確保
  TIMING_start(tm_init_alloc);
  allocArray_Forcing(PrepMemory, TotalMemory, fp, C.NoCompo);
  TIMING_stop(tm_init_alloc);
  
  
  // コンポーネントの体積率を8bitで量子化し，圧力損失コンポの場合にはFORCING_BITをON > bcdにエンコード
  V.setCmpFraction(cmp, d_bcd, d_cvf);
  
  
  
  
  // 周期境界条件が設定されている場合のBCIndexの周期条件の強制同期
  BC.setBCIperiodic(d_bcd, ensPeriodic);
  BC.setBCIperiodic(d_bcp, ensPeriodic);
  BC.setBCIperiodic(d_cdf, ensPeriodic);
  
  
  // bcd/bcp/cdfの同期
  if ( numProc > 1 )
  {
    if ( paraMngr->BndCommS3D(d_bcd, size[0], size[1], size[2], guide, 1) != CPM_SUCCESS ) Exit(0);
    if ( paraMngr->BndCommS3D(d_bcp, size[0], size[1], size[2], guide, 1) != CPM_SUCCESS ) Exit(0);
    if ( paraMngr->BndCommS3D(d_cdf, size[0], size[1], size[2], guide, 1) != CPM_SUCCESS ) Exit(0);
  }
  
  
  // サンプリング点の整合性をチェック
  if ( C.SamplingMode == ON ) MO.checkStatus();
  
  
  // 時間積分幅 deltaT や物理パラメータの設定
  setParameters();
  
  
  // 必要なパラメータをSetBC3Dクラスオブジェクトにコピーする >> setParameters()の後
  BC.setControlVars(&C, mat, &RF, Ex);
  
  
  
  
  // コンポーネントのグローバルインデクス情報を取得し，CompoListの内容とセル数の情報を表示する
  dispGlobalCompoInfo(fp);
  
  
  
  
  // 外部境界面の開口率を計算する
  V.countOpenAreaOfDomain(d_bcd, C.OpenDomain);
  
  
  Hostonly_
  {
    fprintf(fp,"\n----------\n\n\n");
    printf(    "\n----------\n\n\n");
  }
  
  if (C.FIO.IO_Voxel == ON)
  {
    Ex->writeSVX(d_bcd, &C);
    Hostonly_
    {
      fprintf(fp,"\tVoxel file 'example.svx' was written.\n");
      printf(    "\tVoxel file 'example.svx' was written.\n");
      fprintf(fp,"\n----------\n\n\n");
      printf(    "\n----------\n\n\n");
    }
  }
  
  // mid[]を解放する  ---------------------------
  if ( d_mid ) delete [] d_mid;
  
  
  
  // 各ノードの領域情報をファイル出力
  gatherDomainInfo();
  
  
  
  
  TIMING_stop(tm_voxel_prep_sct);
  // ここまでが準備の時間セクション
  
  
  
  // 計算に用いる配列のアロケート ----------------------------------------------------------------------------------
  allocate_Main(TotalMemory);
  
  
  
  // 初期値とリスタート処理 瞬時値と統計値に分けて処理　------------------
  Hostonly_
  {
    printf(    "\n----------\n\n");
    fprintf(fp,"\n----------\n\n");
  }
  
  // リスタートモードの選択
  if ( C.Start != initial_start)
  {
    selectRestartMode();
  }
  
  
  // 瞬時値のリスタート
  TIMING_start(tm_restart);
  Restart(fp);
  TIMING_stop(tm_restart);
  
  
  // 制御インターバルの初期化
  initInterval();
  
  
  // 統計値のリスタート
  if ( C.Mode.StatisticRestart == ON )
  {
    TIMING_start(tm_restart);
    RestartStatistic(fp, flop_task);
    TIMING_stop(tm_restart);
  }
  
  
  // リスタートの最大値と最小値の表示
  if ( C.Start != initial_start )
  {
    RestartDisplayMinmax(fp, flop_task);
  }
  
  
  
  // 利用ライブラリのバージョン番号取得
  C.ver_CPM = cpm_Base::getVersionInfo();
  C.ver_CDM = cdm_DFI::getVersionInfo();
  C.ver_Poly= PL->getVersionInfo();
  C.ver_PM  = PM.getVersionInfo();
  C.ver_CUT = cutlib_VersionInfo();
  C.ver_TP  = tp_ffv.getVersionInfo();
  
  
  // 制御パラメータ，物理パラメータの表示
  Hostonly_
  {
    displayParameters(fp);
  }
  
  
  
  // ドライバ条件のチェック
  BC.checkDriver(fp);
  
  
  
  // 初期条件の条件設定
  setInitialCondition();
  
  
  // サンプリング元となるデータ配列の登録
  if ( C.SamplingMode == ON )
  {
    MO.setDataPtrs(d_v, d_p, d_ie, d_vrt);
  }
  
  
  // 出力ファイルの初期化
  initFileOut();
  
  
  
  // メモリ使用量の表示
  Hostonly_
  {
    printf(    "\n----------\n\n");
    fprintf(fp,"\n----------\n\n");
  }
  
  G_TotalMemory = TotalMemory;
  
  displayMemoryInfo(fp, G_TotalMemory, TotalMemory, "Solver");
  
  
  
  // 履歴出力準備
  //prepHistoryOutput();
  
  
  Hostonly_ if ( fp ) fclose(fp);
  
  
  TIMING_stop(tm_init_sct);
  
  return 1;
}*/
