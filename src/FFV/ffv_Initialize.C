//##################################################################################
//
// FFV-C : Frontflow / violet Cartesian
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2013 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2013 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

/** 
 * @file   ffv_Initialize.C
 * @brief  FFV Class
 * @author kero
 */

#include "ffv.h"
#include "limits.h"


// #################################################################
/* @brief 初期化格子生成、ビットフラグ処理ほか
 * @param [in] argc  main関数の引数の個数
 * @param [in] argv  main関数の引数リスト
 */
int FFV::Initialize(int argc, char **argv)
{
  double TotalMemory   = 0.0;  ///< 計算に必要なメモリ量（ローカル）
  double PrepMemory    = 0.0;  ///< 初期化に必要なメモリ量（ローカル）
  double G_TotalMemory = 0.0;  ///< 計算に必要なメモリ量（グローバル）
  double G_PrepMemory  = 0.0;  ///< 初期化に必要なメモリ量（グローバル）
  double tmp_memory    = 0.0;  ///< 計算に必要なメモリ量（グローバル）？
  
  double flop_task     = 0.0;  ///< flops計算用

  
  
  // 固定パラメータ
  fixedParameters();
  

  
  // ------------------------------------
  FILE* fp = NULL;
  
  // condition fileのオープン
  Hostonly_
  {
    if ( !(fp=fopen("condition.txt", "w")) )
    {
      stamped_printf("\tSorry, can't open 'condition.txt' file. Write failed.\n");
      return -1;
    }
  }
  
  // メッセージ表示
  Hostonly_
  {
    FBUtility::printVersion(fp,     "Welcome to FFV  ", FFV_VERS);
    FBUtility::printVersion(stdout, "Welcome to FFV  ", FFV_VERS);
    
    FBUtility::printVersion(fp,     "FlowBase        ", FB_VERS);
    FBUtility::printVersion(stdout, "FlowBase        ", FB_VERS);
  }
  
  // 入力ファイルの指定
  std::string input_file = argv[1];

  
  // ffvのパラメータローダのインスタンス生成
  TPControl tp_ffv;
  
  tp_ffv.getTPinstance();
  
  
  // パラメータのロードと保持
  if ( tp_ffv.readTPfile(input_file) )
  {
    Hostonly_ stamped_printf("\tInput file '%s' can not find.\n", input_file.c_str());
    Exit(0);
  }
  
  
  
  // TPControlクラスのポインタを各クラスに渡す
  C.importTP(&tp_ffv);
  B.importTP(&tp_ffv);
  M.importTP(&tp_ffv);
  MO.importTP(&tp_ffv);
  
  
  // Intrinsic classの同定
  identifyExample(&tp_ffv, fp);
  
  
  // パラメータの取得と計算領域の初期化，並列モードを返す
  std::string str_para = setupDomain(&tp_ffv);

  
  
  // mat[], cmp[]の作成
  createTable(fp);
  
  
  
  
  // 媒質情報をパラメータファイルから読み込み，媒質リストを作成する
  setMediumList(fp);
  
  
  V.setControlVars(C.NoCompo, Ex);

  B.setControlVars(&C);

  B.countMedium(&C, mat);

  
  
  // CompoListの設定，外部境界条件の読み込み保持
  setBCinfo();
  
  
  
  // タイミング測定の初期化
  if ( C.Mode.Profiling != OFF )
  {
    ModeTiming = ON;
    TIMING__ PM.initialize( tm_END );
    TIMING__ PM.setRankInfo( paraMngr->GetMyRankID() );
    TIMING__ PM.setParallelMode(str_para, C.num_thread, C.num_process);
    set_timing_label();
  }
  
  
  // タイミング測定開始
  TIMING_start(tm_init_sct); 
  
  
  
  // 前処理に用いるデータクラスのアロケート -----------------------------------------------------
  TIMING_start(tm_init_alloc); 
  allocArray_Prep(PrepMemory, TotalMemory);
  TIMING_stop(tm_init_alloc);
  

  
  TIMING_start(tm_voxel_prep_sct);
  

  
  
  
  // 各問題に応じてモデルを設定
  setModel(PrepMemory, TotalMemory, fp);
  
  
  
  // 領域情報の表示
  Hostonly_
  {
    printf("\n---------------------------------------------------------------------------\n");
    printf("\n\t>> Global Domain Information\n\n");
    C.printGlobalDomain(stdout, G_size, G_origin, G_region, pitch);
    
    fprintf(fp,"\n---------------------------------------------------------------------------\n");
    fprintf(fp,"\n\t>> Global Domain Information\n\n");
    C.printGlobalDomain(fp, G_size, G_origin, G_region, pitch);
  }
  
  
  // メモリ消費量の情報を表示
  Hostonly_
  {
    printf(    "\n---------------------------------------------------------------------------\n\n");
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
  }
  G_PrepMemory = PrepMemory;
  
  displayMemoryInfo(fp, G_PrepMemory, PrepMemory, "Preprocessor");
  


  
  
  // Binaryの場合に，SOLIDセルを生成
  if ( !C.isCDS() && (C.Mode.Example == id_Polygon) )
  {
    generateSolid(fp);
  }
  

  
  // ガイドセル上にパラメータファイルで指定する媒質IDを代入する．周期境界の場合の処理も含む．
  for (int face=0; face<NOFACE; face++)
  {
    if ( (C.Mode.Example == id_Jet) && (face==0) )
    {
      // skip
    }
    else
    {
      V.adjMedium_on_GC(face, d_mid, BC.export_OBC(face)->get_Class(),
                      BC.export_OBC(face)->getGuideMedium(), BC.export_OBC(face)->get_PrdcMode());
    }
  }


  
  // Fill
  if ( (C.Mode.Example == id_Polygon) )
  {
    Hostonly_
    {
      printf(    "\n---------------------------------------------------------------------------\n\n");
      fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
      printf(    "\t>> Fill\n\n");
      fprintf(fp,"\t>> Fill\n\n");
    }
    
    fill(fp);
  }
  

  
  
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
    printf("\n---------------------------------------------------------------------------\n\n");
    printf("\t>> Components\n\n");
    C.printNoCompo(stdout);
    printf("\n"); fflush(stdout);
    
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
    fprintf(fp,"\t>> Components\n\n");
    C.printNoCompo(fp);
    fprintf(fp,"\n"); fflush(fp);
  }
  

  
  // セルモニターのモニタ位置をセット
  setupCellMonitor();
  
  
  
  // ボクセルのスキャン
  scanVoxel(fp);
  

  
  
  // CDSの場合，WALLとSYMMETRICのときに，カットを外部境界に接する内部セルに設定 fill以前には行わない
  if ( C.isCDS() ) 
  {
    V.setOBC_Cut(&BC, d_cut);
  }
  
  
  // セルIDのノード間同期
  if ( paraMngr->BndCommS3D(d_mid, size[0], size[1], size[2], guide, 1) != CPM_SUCCESS ) Exit(0);
  
  
  

  // HEX/FANコンポーネントの形状情報からBboxと体積率を計算
  if ( C.isVfraction() )
  {
    TIMING_start(tm_init_alloc);
    allocArray_CompoVF(PrepMemory, TotalMemory);
    TIMING_stop(tm_init_alloc); 
    
    setComponentVF();
  }

  
  
  // 内部周期境界の場合のガイドセルのコピー処理
  V.adjMediumPrdc_Inner(d_mid, cmp);
  
  
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

  
  
  // コンポーネント情報を保存
  //setLocalCmpIdx_Binary();

  
  
  // BCIndexにビット情報をエンコードとコンポーネントインデクスの再構築
  encodeBCindex();
  
  
  // カット情報からBCのコンポーネント情報を設定
  //BboxLocalBC();
  
  
  // 体積力を使う場合のコンポーネント配列の確保
  TIMING_start(tm_init_alloc);
  allocArray_Forcing(PrepMemory, TotalMemory, fp);
  TIMING_stop(tm_init_alloc); 
  

  // コンポーネントの体積率を8bitで量子化し，圧力損失コンポの場合にはFORCING_BITをON > bcdにエンコード
  V.setCmpFraction(cmp, d_bcd, d_cvf);


  
  // RefMediumがMediumList中にあるかどうかをチェックし、RefMatを設定
  if ( (C.RefMat = C.find_ID_from_Label(mat, C.NoMedium, C.RefMedium)) == 0 )
  {
    Hostonly_
    {
      printf(     "RefMat[%d] is not listed in MediumTable.\n", C.RefMat);
      fprintf(fp, "RefMat[%d] is not listed in MediumTable.\n", C.RefMat);
    }
    return -1;
  }


  // 周期境界条件が設定されている場合のBCIndexの周期条件の強制同期
  BC.setBCIperiodic(d_bcd);
  BC.setBCIperiodic(d_bcp);
  BC.setBCIperiodic(d_bcv);

  if ( C.isHeatProblem() )
  {
    BC.setBCIperiodic(d_bh1);
    BC.setBCIperiodic(d_bh2);
  }

  // bcd/bcp/bcv/bchの同期
  if ( numProc > 1 )
  {
    if ( paraMngr->BndCommS3D(d_bcd, size[0], size[1], size[2], guide, 1) != CPM_SUCCESS ) Exit(0);
    if ( paraMngr->BndCommS3D(d_bcp, size[0], size[1], size[2], guide, 1) != CPM_SUCCESS ) Exit(0);
    if ( paraMngr->BndCommS3D(d_bcv, size[0], size[1], size[2], guide, 1) != CPM_SUCCESS ) Exit(0);
    
    if ( C.isHeatProblem() )
    {
      if ( paraMngr->BndCommS3D(d_bh1, size[0], size[1], size[2], guide, 1) != CPM_SUCCESS ) Exit(0);
      if ( paraMngr->BndCommS3D(d_bh2, size[0], size[1], size[2], guide, 1) != CPM_SUCCESS ) Exit(0);
    }
  }

  
  // 法線計算
  getCompoArea();
  
  
  
  // 時間積分幅 deltaT や物理パラメータの設定
  setParameters();


  
  // 必要なパラメータをSetBC3Dクラスオブジェクトにコピーする >> setParameters()の後
  BC.setControlVars(&C, mat, &RF, Ex);
  
  
// ##########
#if 0
  // チェックのため，全計算セルのBCIndexの内容を表示する
  if ( !V.dbg_chkBCIndexP(bcd, bcp, "BCindex.txt") )
  {
    Hostonly_
    {
      stamped_printf("\tVoxInfo::DbgChkBCIndexP()\n");
    }
    return -1;
  }
#endif
// ##########
  
  
  
  // 温度計算の場合の初期値指定
  if ( C.isHeatProblem() )
  {
    B.get_Medium_InitTemp(cmp);
  }
  
  // set phase 
  if ( C.BasicEqs == INCMP_2PHASE )
  {
    B.get_Phase(cmp);
  }
  

  // コンポーネントのグローバルインデクス情報を取得
  setGlobalCmpIdx();
  
  
  
  // CompoListの内容とセル数の情報を表示する
  Hostonly_
  {
    printf(    "\n---------------------------------------------------------------------------\n\n");
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
    printf(    "\t>> Component Information\n\n");
    fprintf(fp,"\t>> Component Information\n\n");
  }
  
  displayCompoInfo(fp);
  

  
  // 外部境界面の開口率を計算する
  V.countOpenAreaOfDomain(d_bcd, C.OpenDomain);
  
  
  Hostonly_
  {
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n\n");
    printf(    "\n---------------------------------------------------------------------------\n\n\n");
  }
  
  
  
  // 性能測定モードがオフのときのみ，出力指定，あるいはMonitorListの場合に，svxファイルを出力する．
  if ( C.Hide.PM_Test == OFF )
  {
    if ( (C.Sampling.log == ON) || (C.FIO.IO_Voxel == Control::Sphere_SVX) )
    {
      Ex->writeSVX(d_mid, &C);
    }
  }
  
  
  // mid[]を解放する  ---------------------------
  if ( d_mid ) delete [] d_mid;
  

  // Monitor Listの処理 --------------------------------------------
  if (C.Sampling.log == ON) setMonitorList();
  
  
  
  // 各ノードの領域情報をファイル出力
  gatherDomainInfo();
  

  
  TIMING_stop(tm_voxel_prep_sct);
  // ここまでがボクセル準備の時間セクション
  
  
// ##########  
#if 0
  write_distance(cut);
#endif
// ##########  
  

  
  // 計算に用いる配列のアロケート ----------------------------------------------------------------------------------
  allocate_Main(TotalMemory);
  

  
  // 初期値とリスタート処理 瞬時値と平均値に分けて処理　------------------
  Hostonly_
  {
    printf(    "\n---------------------------------------------------------------------------\n\n");
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
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
  init_Interval();
  

  // 平均値のリスタート
  if ( (C.Start != initial_start) && (C.Mode.Average == ON) )
  {
    TIMING_start(tm_restart);
    Restart_avrerage(fp, flop_task);
    TIMING_stop(tm_restart);
  }
  

  // リスタートの最大値と最小値の表示
  if ( C.Start != initial_start)
  {
    Restart_display_minmax(fp, flop_task);
  }

  
  
  // 利用ライブラリのバージョン番号取得
  C.ver_CPM = cpm_Base::getVersionInfo();
  C.ver_CIO = cio_DFI::getVersionInfo();
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
  


  if ( C.Sampling.log == ON ) 
  {
    // サンプリング指定がある場合，モニタ結果出力ファイル群のオープン
    MO.openFile(C.HistoryMonitorName.c_str());
    
    // サンプリング元となるデータ配列の登録
    if ( C.isHeatProblem() ) 
    {
      MO.setDataPtrs(d_v, d_p, d_t);
    }
    else 
    {
      MO.setDataPtrs(d_v, d_p);
    }
  }
  
  
  /* 20130611 commentout PLOT3D形状データの書き出し
  PLT3D.Initialize(size, guide, deltaX, dfi_mng[var_Plot3D], &C, &FP3DW, &DFI, d_ws, d_p, d_wo, d_v, d_t, d_p0, d_wv, d_bcv, d_bcd);
  
  if (C.FIO.Format == plt3d_fmt)
  {
    PLT3D.setValuePlot3D();
    if (C.P3Op.IS_xyz == ON) PLT3D.OutputPlot3D_xyz(CurrentStep, origin, pitch);// ---> moving grid を考慮したときOutputPlot3D_postに組み込む
    if (C.P3Op.IS_DivideFunc == ON)
    {
      if (C.P3Op.IS_function_name == ON) Hostonly_ PLT3D.OutputPlot3D_function_name_divide();
    }
    else
    {
      if (C.P3Op.IS_function_name == ON) Hostonly_ PLT3D.OutputPlot3D_function_name();
    }
    if (C.P3Op.IS_fvbnd == ON) PLT3D.OutputPlot3D_fvbnd();
  }
  */

  
  // 出力ファイルの初期化
  init_FileOut();
  

  // セッションを開始したときに、初期値をファイル出力  リスタートと性能測定モードのときには出力しない
  if ( (C.Hide.PM_Test == OFF) && (0 == CurrentStep) )
  {
    flop_task = 0.0;
    FileOutput(flop_task);
    
    /* 20130611 commentout
    if (C.FIO.Format == plt3d_fmt) PLT3D.OutputPlot3D_post(CurrentStep, CurrentTime, v00, origin, pitch, dfi_mng[var_Plot3D], flop_task);
    */
    
    if ( (C.Mode.Average == ON) && (C.Start != initial_start) )
    {
      double flop_count=0.0;
      AverageOutput(flop_count);
    }
  }

  
  // 粗い格子を用いたリスタート時には出力
  if ( (C.Start == restart_sameDiv_refinement) || (C.Start == restart_diffDiv_refinement) )
  {
    flop_task = 0.0;
    FileOutput(flop_task, true);
    
    /* 20130611 commentout 
    if (C.FIO.Format == plt3d_fmt) PLT3D.OutputPlot3D_post(CurrentStep, CurrentTime, v00, origin, pitch, dfi_mng[var_Plot3D], flop_task);
     */
  }


  // SOR2SMA
  switch (IC[ItrCtl::ic_prs_pr].get_LS())
  {
    case SOR2SMA:
    case SOR2CMA:
    case GMRES:
    case RBGS:
		case PCG:
		case PBiCGSTAB:
      allocate_SOR2SMA_buffer(TotalMemory);
      break;
  }
  
  // Krylov subspace
  switch (IC[ItrCtl::ic_prs_pr].get_LS())
  {
    case GMRES:
      allocArray_Krylov(TotalMemory);
      break;
      
    case PCG:
			allocArray_PCG(TotalMemory);
			break;
      
		case PBiCGSTAB:
			allocArray_PBiCGSTAB(TotalMemory);
			break;
  }
  
  
  // メモリ使用量の表示
  Hostonly_
  {
    printf(    "\n---------------------------------------------------------------------------\n\n");
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
  }
  
  G_TotalMemory = TotalMemory;
  
  displayMemoryInfo(fp, G_TotalMemory, TotalMemory, "Solver");

  
  
  // 履歴出力準備
  prep_HistoryOutput();
  
  
  
  Hostonly_ if ( fp ) fclose(fp);
  
  TIMING_stop(tm_init_sct);
  
  // チェックモードの場合のコメント表示，前処理のみで中止---------------------------------------------------------
  if ( C.CheckParam == ON)
  {
		Hostonly_
    {
      printf(     "\n\tCheck mode --- Only pre-process\n\n");
      fprintf(fp, "\n\tCheck mode --- Only pre-process\n\n");
    }
    return 0;
	}

  return 1;
}


// #################################################################
/* @brief 主計算部分に用いる配列のアロケーション
 * @param [in,out] total ソルバーに使用するメモリ量
 */
void FFV::allocate_Main(double &total)
{
  TIMING_start(tm_init_alloc);
  allocArray_Main(total);
  
  //allocArray_Collocate (total);
  
  if ( C.LES.Calc == ON )
  {
    allocArray_LES (total);
  }
  
  if ( C.isHeatProblem() )
  {
    allocArray_Heat(total);
  }
  
  if ( (C.AlgorithmF == Control::Flow_FS_AB2) || (C.AlgorithmF == Control::Flow_FS_AB_CN) )
  {
    allocArray_AB2(total);
  }
  
  if ( C.BasicEqs == INCMP_2PHASE )
  {
    allocArray_Interface(total);
  }
  
  // 時間平均用の配列をアロケート
  if ( C.Mode.Average == ON )
  {
    allocArray_Average(total);
  }
  
  TIMING_stop(tm_init_alloc);
}



/* #################################################################
// ポリゴングループの座標値からboxを計算する
void FFV::calcBboxfromPolygonGroup()
{
  // float で定義すること
  float poly_org[3];
  float poly_dx[3];
  
  // 有次元に変換 Polylib: 並列計算領域情報　ポリゴンは実スケールで，ガイドセル領域部分も含めて指定する
  poly_dx[0]  = pitch[0] * C.RefLength;
  poly_dx[1]  = pitch[1] * C.RefLength;
  poly_dx[2]  = pitch[2] * C.RefLength;
  poly_org[0] = origin[0]* C.RefLength;
  poly_org[1] = origin[1]* C.RefLength;
  poly_org[2] = origin[2]* C.RefLength;
  
  
  // ポリゴン情報へのアクセス
  vector<PolygonGroup*>* pg_roots = PL->get_root_groups();
  vector<PolygonGroup*>::iterator it;
  
  PolylibNS::Vec3f m_min, m_max;
  PolylibNS::Vec3f t1(poly_org), t2(poly_dx), t3;
  
  t3.assign((float)size[0]*t2.t[0], (float)size[1]*t2.t[1], (float)size[2]*t2.t[2]);
  
  // サブドメインの1層外側までをサーチ対象とする
  m_min = t1 - t2;      
  m_max = t1 + t3 + t2;
  printf("Search area Bbox min : %f %f %f\n", m_min.t[0], m_min.t[1], m_min.t[2]);
  printf("Search area Bbox max : %f %f %f\n", m_max.t[0], m_max.t[1], m_max.t[2]);
  
  
  // ポリゴングループのループ
  int m=0;
  for (it = pg_roots->begin(); it != pg_roots->end(); it++)
  {
    string m_pg = PolyPP[m].get_Group();
    
    // false; ポリゴンが一部でもかかればピックアップ
    vector<Triangle*>* trias = PL->search_polygons(m_pg, m_min, m_max, false);
    
    PolylibNS::Vec3f *p;
    FB::Vec3f bbox_min( 1.0e6,  1.0e6,  1.0e6);
    FB::Vec3f bbox_max(-1.0e6, -1.0e6, -1.0e6);
    unsigned c=0;
    vector<Triangle*>::iterator it2;
    
    for (it2 = trias->begin(); it2 != trias->end(); it2++)
    {
      p = (*it2)->get_vertex();
      
      // PolulibNS::Vec3f >> FB::Vec3f
      FB::Vec3f p0(p[0].t[0], p[0].t[1], p[0].t[2]);
      FB::Vec3f p1(p[1].t[0], p[1].t[1], p[1].t[2]);
      FB::Vec3f p2(p[2].t[0], p[2].t[1], p[2].t[2]);
      
      CompoFraction::get_min(bbox_min, p0);
      CompoFraction::get_min(bbox_min, p1);
      CompoFraction::get_min(bbox_min, p2);
      
      CompoFraction::get_max(bbox_max, p0);
      CompoFraction::get_max(bbox_max, p1);
      CompoFraction::get_max(bbox_max, p2);
      
#if 0
      printf("%d : p0=(%6.3e %6.3e %6.3e)  p1=(%6.3e %6.3e %6.3e) p2=(%6.3e %6.3e %6.3e) n=(%6.3e %6.3e %6.3e)\n", c++,
             p[0].t[0], p[0].t[1], p[0].t[2],
             p[1].t[0], p[1].t[1], p[1].t[2],
             p[2].t[0], p[2].t[1], p[2].t[2],
             n.t[0], n.t[1], n.t[2]);
#endif
    }
    
    PolyPP[m].set_Min(bbox_min);
    PolyPP[m].set_Max(bbox_max);
    
    printf("%20s : (%6.3e %6.3e %6.3e) - (%6.3e %6.3e %6.3e)\n",
           m_pg.c_str(), bbox_min.x, bbox_min.y, bbox_min.z,
           bbox_max.x, bbox_max.y, bbox_max.z);
    
    delete trias; // 後始末
    m++;
  }
  
}
*/

// #################################################################
/* @brief 全Voxelモデルの媒質数とKOSの整合性をチェック
 * @retval エラーコード
 */
bool FFV::chkMediumConsistency()
{
  int nmSolid = C.NoMediumSolid;
  int nmFluid = C.NoMediumFluid;
  
  if ( numProc > 1 )
  {
    int nms = nmSolid;
    int nmf = nmFluid;
    paraMngr->Allreduce(&nms, &nmSolid, 1, MPI_SUM);
    paraMngr->Allreduce(&nmf, &nmFluid, 1, MPI_SUM);
  }
  
  if ( (nmFluid == 0) && (nmSolid == 0) )
  {
    Hostonly_ printf("\tError : No medium\n");
    return false;
  }
  
  switch (C.KindOfSolver)
  {
    case FLOW_ONLY:
    case THERMAL_FLOW:
    case THERMAL_FLOW_NATURAL:
      
      if ( nmFluid == 0 )
      {
        Hostonly_ printf("\tError : No FLUID medium\n");
        return false;
      }
      break;
      
    case CONJUGATE_HEAT_TRANSFER:
      if ( ( nmFluid == 0 ) || ( nmSolid == 0 ) )
      {
        Hostonly_ printf("\tError : Fluid/Solid should have at least one medium.\n");
        return false;
      }
      break;
      
    case SOLID_CONDUCTION:
      if ( nmSolid == 0 )
      {
        Hostonly_ printf("\tError : No Solid medium\n");
        return false;
      }
      break;
  };
  
  return true;
}



// #################################################################
/* @brief 時刻をRFクラスからv00[4]にコピーする
 * @param [in] time 設定する時刻
 */
void FFV::copyV00fromRF(double m_time)
{
  RF.setV00(m_time);
  
  double g[4];
  RF.copyV00(g);
  for (int i=0; i<4; i++) v00[i]=(REAL_TYPE)g[i];
}


// #################################################################
/* @brief mat[], cmp[]のテーブルを作成
 * @param [in] fp  file pointer
 */
void FFV::createTable(FILE* fp)
{
  // コンポーネント数，境界条件数，媒質数を取得し，配列をアロケートする
  C.getNoOfComponent();
  
  // 媒質リストをインスタンス [0]はダミーとして利用しないので，配列の大きさはプラス１する
  if ( !(mat = new MediumList[C.NoCompo+1]) ) Exit(0);
  
  // CompoListクラスをインスタンス [0]はダミーとして利用しないので，配列の大きさはプラス１する
  if ( !(cmp = new CompoList[C.NoCompo+1]) ) Exit(0);
  
  // CompoList, MediumListのポインタをセット
  BC.importCMP_MAT(cmp, mat);
  
  Hostonly_
  {
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
    fprintf(fp,"\n\t>> Tables\n\n");
    fprintf(fp,"\tNo. of Medium    = %3d\n", C.NoMedium);
    fprintf(fp,"\tNo. of LocalBC   = %3d\n", C.NoBC);
    fprintf(fp,"\t----------------------\n");
    fprintf(fp,"\tNo. of Component = %3d\n", C.NoCompo);
    fprintf(fp,"\n");
    
    printf(    "\n---------------------------------------------------------------------------\n\n");
    printf(    "\n\t>> Tables\n\n");
    printf(    "\tNo. of Medium    = %3d\n", C.NoMedium);
    printf(    "\tNo. of LocalBC   = %3d\n", C.NoBC);
    printf(    "\t----------------------\n");
    printf(    "\tNo. of Component = %3d\n", C.NoCompo);
    printf(    "\n");
  }
  
}



// #################################################################
/* @brief CompoListの情報を表示する
 * @param [in]  fp   ファイルポインタ
 */
void FFV::displayCompoInfo(FILE* fp)
{
  // サマリー
  Hostonly_
  {
    B.printCompoSummary(stdout, cmp, C.BasicEqs);
    B.printCompoSummary(fp, cmp, C.BasicEqs);
  }
  
  // セル数の情報を表示する
  Hostonly_
  {
    double cr = (double)G_Wcell/ ( (double)G_size[0] * (double)G_size[1] * (double)G_size[2]) *100.0;
    
    fprintf(stdout, "\tThis model includes %4d solid %s  [Solid cell ratio inside computational domain : %9.5f percent]\n\n", 
            C.NoMediumSolid, (C.NoMediumSolid>1) ? "IDs" : "ID", cr);
    fprintf(fp, "\tThis model includes %4d solid %s  [Solid cell ratio inside computational domain : %9.5f percent]\n\n", 
            C.NoMediumSolid, (C.NoMediumSolid>1) ? "IDs" : "ID", cr);
  }
  
  
  
  // 詳細情報
  if ( C.NoBC >0 )
  {
    Hostonly_
    {
      B.printCompo( stdout, compo_global_bbox, mat, cmp, BC.export_OBC() );
      B.printCompo( fp, compo_global_bbox, mat, cmp, BC.export_OBC() );
    }
  }
  
  
  // Check consistency of boundary condition
  for (int n=1; n<=C.NoCompo; n++)
  {
    if ( cmp[n].isKindCompo() )
    {
      if ( cmp[n].getType() == HT_SN )
      {
        if ( (C.KindOfSolver == FLOW_ONLY) || (C.KindOfSolver == THERMAL_FLOW) || (C.KindOfSolver == SOLID_CONDUCTION) )
        {
          Hostonly_ printf("\tInconsistent parameters of combination between Kind of Solver and Heat Transfer type SN. Check QBCF\n");
          fflush(stdout);
          Exit(0);
        }
      }
      if ( cmp[n].getType() == HT_SF )
      {
        if ( (C.KindOfSolver == FLOW_ONLY) || (C.KindOfSolver == THERMAL_FLOW_NATURAL) || (C.KindOfSolver == SOLID_CONDUCTION) )
        {
          Hostonly_ printf("\tInconsistent parameters of combination between Kind of Solver and Heat Transfer type SF. Check QBCF\n");
          fflush(stdout);
          Exit(0);
        }
      }
      
    }
  }
  
}

// #################################################################
/* @brief メモリ消費情報を表示
 * @param [in]     fp    ファイルポインタ
 * @param [in,out] G_mem グローバルメモリサイズ
 * @param [in]     L_mem ローカルメモリサイズ
 * @param [in]     str   表示用文字列
 */
void FFV::displayMemoryInfo(FILE* fp, double G_mem, double L_mem, const char* str)
{
  if ( numProc > 1 )
  {
    double tmp_memory = G_mem;
    if ( paraMngr->Allreduce(&tmp_memory, &G_mem, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
  }
  
  Hostonly_
  {
    FBUtility::MemoryRequirement(str, G_mem, L_mem, stdout);
    FBUtility::MemoryRequirement(str, G_mem, L_mem, fp);
  }
  
  Hostonly_
  {
    printf("\n\n");
    fprintf(fp, "\n\n");
  }
}


// #################################################################
/* @brief 制御パラメータ，物理パラメータの表示
 * @param [in]  fp   ファイルポインタ
 */
void FFV::displayParameters(FILE* fp)
{
  // 20130611 C.displayParams(stdout, fp, IC, &DT, &RF, mat, &FP3DW);
  C.displayParams(stdout, fp, IC, &DT, &RF, mat);
  Ex->printPara(stdout, &C);
  Ex->printPara(fp, &C);
  
  // 外部境界面の開口率を表示
  C.printOuterArea(stdout, G_Fcell, G_Acell, G_size);
  C.printOuterArea(fp, G_Fcell, G_Acell, G_size);
  
  // 境界条件のリストと外部境界面のBC設定を表示

  printf(    "\n---------------------------------------------------------------------------\n\n");
  fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
  printf(    "\t>> Outer Boundary Conditions\n\n");
  fprintf(fp,"\t>> Outer Boundary Conditions\n\n");
  
  B.printFaceOBC(stdout, G_region, BC.export_OBC(), mat);
  B.printFaceOBC(fp, G_region, BC.export_OBC(), mat);

  
  // モニタ情報の表示
  if ( C.Sampling.log == ON ) {
    
    MO.printMonitorInfo(stdout, C.HistoryMonitorName.c_str(), false); // ヘッダのみ
    
    FILE *fp_mon=NULL;
    Hostonly_
    {
      if ( !(fp_mon=fopen("sampling_info.txt", "w")) )
      {
        stamped_printf("\tSorry, can't open 'sampling_info.txt' file. Write failed.\n");
        //return -1;
        Exit(0);
      }
    }
    
    MO.printMonitorInfo(fp_mon, C.HistoryMonitorName.c_str(), true);  // 詳細モード
    Hostonly_ if ( fp_mon ) fclose(fp_mon);
  }
}


// #################################################################
/* @brief 計算領域情報を設定する
 * @param [in] tp_dom  TPControlクラス
 */
void FFV::DomainInitialize(TPControl* tp_dom)
{
  // メンバ変数にパラメータをロード : 分割指示 (1-with / 2-without)
  int div_type = getDomainInfo(tp_dom);

  
// ##########  
# if 0
  printDomainInfo();
  fflush(stdout);
#endif
// ##########
  

  // 袖通信の最大数
  size_t Nvc  = (size_t)C.guide;
  size_t Ncmp = 6; // カットを使うので6成分（面）
  
  int m_sz[3]  = {G_size[0], G_size[1], G_size[2]};
  int m_div[3] = {G_division[0], G_division[1], G_division[2]};
  

  // 有次元の場合に無次元化する　paraMngr->VoxelInit()で計算する基本量
  if (C.Unit.Param == DIMENSIONAL )
  {
    for (int i=0; i<3; i++) {
      pitch[i]    /= C.RefLength;
      G_origin[i] /= C.RefLength;
      G_region[i] /= C.RefLength;
    }
  }

  double m_org[3] = {(double)G_origin[0], (double)G_origin[1], (double)G_origin[2]};
  double m_reg[3] = {(double)G_region[0], (double)G_region[1], (double)G_region[2]};
  
  
  // 領域分割モードのパターン
  //     分割指定(G_div指定)         domain.txt
  // 1)  G_div指定なし |  G_orign + G_region + (G_pitch || G_voxel)
  // 2)  G_div指定あり |  G_orign + G_region + (G_pitch || G_voxel)
  // 3)  G_div指定なし |          + ActiveDomainInfo
  // 4)  G_div指定あり |          + ActiveDomainInfo
  
  /* Policy
   * G_regionは必須
   * G_voxelとG_pitchでは，G_voxelが優先．両方が指定されている場合はエラー
   * G_pitchが指定されており，割り切れない場合は，領域を拡大
   *        G_voxel = (int)ceil(G_region/G_pitch)
   *        G_region = G_pitch * G_voxel
   */
  
  
  // 分割数を元に分割する >> CPMlibの仕様
  switch (div_type) 
  {
    case 1: // 分割数が指示されている場合
      if ( paraMngr->VoxelInit(m_div, m_sz, m_org, m_reg, Nvc, Ncmp) != CPM_SUCCESS )
      {
        cout << "Domain decomposition error : " << endl;
        Exit(0);
      }
      break;
      
    case 2: // 分割数が指示されていない場
      if ( paraMngr->VoxelInit(m_sz, m_org, m_reg, Nvc, Ncmp) != CPM_SUCCESS )
      {
        cout << "Domain decomposition error : " << endl;
        Exit(0);
      }
      break;
      
    default:
      Exit(0);
      break;
  }


  // 分割後のパラメータをDomainInfoクラスメンバ変数に保持
  setNeighborInfo(C.guide);
  
  
  // チェック
  unsigned long tz = (unsigned long)size[0] * (unsigned long)size[1] * (unsigned long)size[2];
  if ( tz >= UINT_MAX)
  {
    Hostonly_ stamped_printf("\n\tError : Product of size[] exceeds UINT_MAX\n\n");
    Exit(0);
  }

}


// #################################################################
/* @brief BCIndexにビット情報をエンコードする
 */
void FFV::encodeBCindex()
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  // 基本ビット情報（Active, State, コンポ，媒質情報）を全領域についてエンコードする
  V.setBCIndexBase(d_bcd, d_mid, d_cvf, mat, cmp, L_Acell, G_Acell, C.KindOfSolver);
  
  
  // STATEとACTIVEビットのコピー
  V.copyBCIbase(d_bcp, d_bcd);
  V.copyBCIbase(d_bcv, d_bcd);
  
  if ( C.isHeatProblem() )
  {
    V.copyBCIbase(d_bh1, d_bcd);
    V.copyBCIbase(d_bh2, d_bcd);
  }
  
  
  
  // BCIndexP に圧力計算のビット情報をエンコードする -----
  if ( C.isCDS() )
  {
    C.NoWallSurface = V.setBCIndexP(d_bcd, d_bcp, d_mid, &BC, cmp, C.Mode.Example, d_cut, d_bid, true);
  }
  else // binary
  {
    C.NoWallSurface = V.setBCIndexP(d_bcd, d_bcp, d_mid, &BC, cmp, C.Mode.Example, d_cut, d_bid, false);
  }
  
#if 0
  V.dbg_chkBCIndexP(d_bcd, d_bcp, "BCindexP.txt", cmp);
#endif
  
  
  
  // BCIndexV に速度計算のビット情報をエンコードする -----
  //if ( C.isCDS() )
  //{
  V.setBCIndexV(d_bcv, d_mid, d_bcp, &BC, cmp, C.Mode.Example, true, d_cut, d_bid);
  //}
  //else // binary
  //{
  //  V.setBCIndexV(d_bcv, d_mid, d_bcp, &BC, cmp, C.Mode.Example);
  //}
  
  
  
  // ##########
#if 0
  V.dbg_chkBCIndexV(d_bcv, "BCindexV.txt");
#endif
  // ##########
  
  // BCIndexT に温度計算のビット情報をエンコードする -----
  if ( C.isHeatProblem() )
  {
    if ( C.isCDS() )
    {
      V.setBCIndexH(d_bcd, d_bh1, d_bh2, d_mid, &BC, C.KindOfSolver, cmp, true, d_cut, d_bid);
    }
    else // binary
    {
      V.setBCIndexH(d_bcd, d_bh1, d_bh2, d_mid, &BC, C.KindOfSolver, cmp);
    }
    
    
    // ##########
#if 0
    V.dbg_chkBCIndexH(d_bcv, "BCindexH.txt");
#endif
    // ##########
    
  }
  
  // 内部領域のFluid, Solidのセル数を数える C.Wcell(Local), G_Wcell(global)
  V.countCellState(L_Wcell, G_Wcell, d_bcd, SOLID);
  V.countCellState(L_Fcell, G_Fcell, d_bcd, FLUID);
  
  
  // getLocalCmpIdx()などで作成したコンポーネントのインデクスの再構築
  // ここまでのbboxはmid[]でサーチした結果，BCindex処理の過程で範囲が変わるので変更
  resizeCompoBbox();
}



// #################################################################
/* @brief 初期インデクスの情報を元に，一層拡大したインデクス値を返す
 * @param [in,out] m_st 拡大された開始点（Fortranインデクス）
 * @param [in,out] m_ed 拡大された終了点（Fortranインデクス）
 * @param [in]     st_i 開始点（Cインデクス）
 * @param [in]     len  コンポーネントの存在長さ
 * @param [in]     m_x  軸方向のサイズ
 * @param [in]     dir  方向
 * @param [in]     m_id キーID
 */
void FFV::EnlargeIndex(int& m_st, int& m_ed, const int st_i, const int len, const int m_x, const int dir, const int m_id)
{
  int ed_i = st_i + len - 1;
  int n_st = st_i - 1;
  int n_ed = ed_i + 1;
  int max_c1 = m_x + guide;
  
  int label_st, label_ed;
  
  switch (dir) 
  {
    case 0:
      label_st = X_MINUS;
      label_ed = X_PLUS;
      break;
      
    case 1:
      label_st = Y_MINUS;
      label_ed = Y_PLUS;
      break;
      
    case 2:
      label_st = Z_MINUS;
      label_ed = Z_PLUS;
      break;
      
    default:
      Hostonly_ stamped_printf("\tError : DIRECTION\n");
      Exit(0);
  }
  
  // BVが-方向のガイドセル内のみにある場合
  if ( ed_i < guide )
  {
    if( nID[label_st] < 0 ) // 計算領域の外部面に接する場合は，対象外
    {
      m_st = 0;
      m_ed = 0;
    }
    else // 計算領域内部にある場合（並列時）
    {
      if ( n_ed == guide ) // ガイドセル1層外側の場合
      {
        m_st = 1; // F index
        m_ed = 1; // F index
      }
      else
      {
        m_st = 0;
        m_ed = 0;
      }
    }
  }
  
  // BVが+方向のガイドセル内のみにある場合
  else if ( st_i >= max_c1 )
  {
    if( nID[label_ed] < 0 ) // 計算領域の外部面に接する場合は，対象外
    {
      m_st = 0;
      m_ed = 0;
    }
    else
    {
      if ( n_st == (max_c1 - 1) ) // ガイドセル1層外側の場合
      {
        m_st = m_x; // F index
        m_ed = m_x; // F index
      }
      else
      {
        m_st = 0;
        m_ed = 0;
      }
    }
    //debug; Hostonly_ printf("(2) dir=%d : st=%d ed=%d\n", dir, m_st, m_ed);
  }
  
  // BVが内部領域のみにある場合（逐次・並列で同じ処理）
  else if ( (st_i >= guide) && (ed_i < max_c1) )
  {
    if ( st_i == guide ) // 最外層セル
    {
      m_st = 1; // F index
    }
    else
    {
      m_st = n_st + 1 - guide; // F index
    }
    
    if ( ed_i == (max_c1 - 1) ) // 最外層セル
    {
      m_ed = m_x; // F index
    }
    else // 内部
    {
      m_ed = n_ed + 1 - guide; // F index
    }
    //debug; Hostonly_ printf("(3) dir=%d : st=%d ed=%d\n", dir, m_st, m_ed);
  }
  
  // BVが両方向のガイドセルにまたがる場合（逐次・並列で同じ処理）
  else if ( (st_i < guide) && (ed_i >= max_c1) )
  {
    m_st = 1; // F index
    m_ed = m_x; // F index
    //debug; Hostonly_ printf("(4) dir=%d : st=%d ed=%d\n", dir, m_st, m_ed);
  }
  
  // BVが-方向のガイドセルから内部領域にある場合
  else if ( (st_i < guide) && (ed_i < max_c1) )
  {
    m_st = 1; // F index
    
    if ( ed_i == (max_c1 - 1) ) // 最外層セル
    {
      m_ed = m_x; // F index
    }
    else // 内部
    {
      m_ed = n_ed + 1 - guide; // F index
    }
    //debug; Hostonly_ printf("(5) dir=%d : st=%d ed=%d\n", dir, m_st, m_ed);
  }
  
  // BVが+方向のガイドセルから内部領域にある場合
  else if ( (st_i < max_c1) && (ed_i >= max_c1) )
  {
    m_ed = m_x; // F index
    
    if ( st_i == guide ) // 端点
    {
      m_st = 1; // F index
    }
    else // 内部
    {
      m_st = n_st + 1 - guide; // F index
    }
    //debug; Hostonly_ printf("(6) dir=%d : st=%d ed=%d\n", dir, m_st, m_ed);
  }
  
  else
  {
    string m_dir;
    if      ( dir == 0 ) m_dir = "X";
    else if ( dir == 1 ) m_dir = "Y";
    else                 m_dir = "Z";
    
    Hostonly_
    {
      stamped_printf("\tError : Unexpected case for ID=%d, (%d - %d): %s\n", m_id, st_i, ed_i, m_dir.c_str());
    }
    Exit(0);
  }
  
}



// #################################################################
/* @brief ポリゴンの場合のフィル操作
 * @param [in] fp    ファイルポインタ
 */
void FFV::fill(FILE* fp)
{
  
  // 指定媒質の属性をチェック
  bool flag = false;
  
  for (int i=1; i<=C.NoCompo; i++)
  {
    if ( (i == C.Fill_Fluid) && (cmp[i].getState() == FLUID) )
    {
      flag = true;
    }
  }
  if ( !flag )
  {
    Hostonly_ printf("\tSpecified Medium of filling fluid is not FLUID\n");
    Exit(0);
  }
  
  
  flag = false;
  for (int i=1; i<=C.NoCompo; i++)
  {
    if ( (i == C.Fill_Solid) && (cmp[i].getState() == SOLID) )
    {
      flag = true;
    }
  }
  if ( !flag )
  {
    Hostonly_ printf("\tSpecified Medium of filling solid is not SOLID\n");
    Exit(0);
  }

  
  
  unsigned long target_count; ///< フィルの対象となるセル数
  unsigned long replaced;     ///< 置換された数
  unsigned long filled;       ///< FLUIDでフィルされた数
  
  
  // 最初にフィル対象のセル数を求める >> 全計算内部セル数
  unsigned long total_cell = (unsigned long)size[0] * (unsigned long)size[1] * (unsigned long)size[2];
  
  
  if ( numProc > 1 )
  {
    unsigned long tmp_fc = total_cell;
    if ( paraMngr->Allreduce(&tmp_fc, &total_cell, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
  }
  
  filled = V.countPainted(d_mid);
  
  target_count = total_cell - filled;
  
  Hostonly_
  {
    printf(    "\tFill initialize -----\n\n");
    fprintf(fp,"\tFill initialize -----\n\n");

    printf    ("\t\tWhole num. of cells    = %16ld\n", total_cell);
    fprintf(fp,"\t\tWhole num. of cells    = %16ld\n", total_cell);
    
    printf    ("\t\tGenerated solid        = %16ld\n", filled);
    fprintf(fp,"\t\tGenerated solid        = %16ld\n", filled);
    
    printf    ("\t\tInitial target count   = %16ld\n", target_count);
    fprintf(fp,"\t\tInitial target count   = %16ld\n", target_count);
  }
  
  
  
  
  
  // チェック用のリストのアロケート
  int* list = new int[C.NoCompo+1];
  
  // CellMonitorをリストアップ
  for (int n=1; n<=C.NoCompo; n++)
  {
    list[n] = 0;
    if ( cmp[n].isMONITOR() ) list[n] = n;
  }
  
  
  
  // 1st pass
  
  Hostonly_
  {
    printf(    "\t1st Fill -----\n\n");
    fprintf(fp,"\t1st Fill -----\n\n");
    printf    ("\t\tFilling medium         : %s\n", mat[C.Fill_Fluid].getAlias().c_str());
    fprintf(fp,"\t\tFilling medium         : %s\n", mat[C.Fill_Fluid].getAlias().c_str());
    printf(    "\t\tHint                   : %s\n", FBUtility::getDirection(C.Fill_Hint).c_str());
    fprintf(fp,"\t\thint                   : %s\n", FBUtility::getDirection(C.Fill_Hint).c_str());
  }
  
  // ヒントが与えられている場合
  // 確実に流体のセルのみをペイントする
  if ( C.Fill_Hint >= 0 )
  {

    filled = (unsigned long)V.fillSeed(d_mid, C.Fill_Hint, C.Fill_Fluid, d_cut);

    if ( numProc > 1 )
    {
      unsigned long tmp_fs = filled;
      if ( paraMngr->Allreduce(&tmp_fs, &filled, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
    }
    
    if ( filled == 0 )
    {
      Hostonly_
      {
        printf(    "No cells to paint\n");
        fprintf(fp,"No cells to paint\n");
      }
      Exit(0);
    }
    
    Hostonly_
    {
      printf(    "\t\tPainted cells          : %16ld\n", filled);
      fprintf(fp,"\t\tPainted cells          : %16ld\n", filled);
    }
  }

  // ペイントされたシードセル数をターゲットから差し引く
  target_count -= filled;
  
  
  Hostonly_
  {
    printf(    "\t\tRemaining target cells = %16ld\n\n", target_count);
    fprintf(fp,"\t\tRemaining target cells = %16ld\n\n", target_count);
  }

  
  
  // BIDによるフィル
  // 隣接する流体セルと接続しており，かつ固体セルに挟まれていないセルのみペイントする
  
  int c=-1; // iteration
  
  while (target_count > 0) {
    
    unsigned fs;
    filled = (unsigned long)V.fill_by_bid(d_bid, d_mid, d_cut, C.Fill_Fluid, fs, list);
    replaced = (unsigned long)fs;
    
    if ( numProc > 1 )
    {
      unsigned long t_fc = filled;
      if ( paraMngr->Allreduce(&t_fc, &filled, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
      
      t_fc = replaced;
      if ( paraMngr->Allreduce(&t_fc, &replaced, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
    }
    
    if ( filled == 0 ) break; // フィル対象がなくなったら終了
    
    target_count -= filled;
    target_count -= replaced;
    c++;
    
    // 同期
    if ( numProc > 1 )
    {
      if ( paraMngr->BndCommS4DEx(d_cut, 6, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
      if ( paraMngr->BndCommS3D(d_mid, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
      if ( paraMngr->BndCommS3D(d_bid, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
    }
  }
  
  Hostonly_
  {
    printf(    "\t\tBID Iteration          = %5d\n", c+1);
    fprintf(fp,"\t\tBID Iteration          = %5d\n", c+1);
    printf(    "\t\t\tFLUID filled   = %16ld\n", filled);
    fprintf(fp,"\t\t\tFLUID filled   = %16ld\n", filled);
    printf(    "\t\t\tSOLID replaced = %16ld\n", replaced);
    fprintf(fp,"\t\t\tSOLID replaced = %16ld\n", replaced);
    printf(    "\t\t\tRemaining cell = %16ld\n\n", target_count);
    fprintf(fp,"\t\t\tRemaining cell = %16ld\n\n", target_count);
  }
  
  
  
  // midによる穴埋め
  c = -1;
  replaced = 0;
  
  while (target_count > 0) {
    
    unsigned z1 = V.fill_by_mid(d_bid, d_mid, d_cut, C.Fill_Fluid, list);
    
    if ( numProc > 1 )
    {
      unsigned t_fc = z1;
      if ( paraMngr->Allreduce(&t_fc, &z1, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
    }
    
    replaced += (unsigned long)z1;
    
    if ( z1 == 0 ) break; // フィル対象がなくなったら終了

    c++;
    
    // 同期
    if ( numProc > 1 )
    {
      if ( paraMngr->BndCommS4DEx(d_cut, 6, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
      if ( paraMngr->BndCommS3D(d_mid, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
      if ( paraMngr->BndCommS3D(d_bid, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
    }
  }
  
  Hostonly_
  {
    printf(    "\t\tHoll filling Iteration = %5d\n", c+1);
    fprintf(fp,"\t\tHoll filling Iteration = %5d\n", c+1);
    printf(    "\t\t\tFilled by SOLID= %16ld\n", replaced);
    fprintf(fp,"\t\t\tFilled by SOLID= %16ld\n", replaced);
  }
  
  
  
  // 2nd pass
  
  Hostonly_
  {
    printf(    "\n\t2nd Fill -----\n\n");
    fprintf(fp,"\n\t2nd Fill -----\n\n");
  }
  
  // 既にペイントした流体セルをクリア
  unsigned long fz = (unsigned long)V.fillReplace(d_mid, C.Fill_Fluid, 0);
  
  if ( numProc > 1 )
  {
    unsigned long t_fc = fz;
    if ( paraMngr->Allreduce(&t_fc, &fz, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
  }
  
  target_count = fz;
  
  
  Hostonly_
  {
    printf(    "\t\tCleared cell           = %16ld\n", fz);
    fprintf(fp,"\t\tCleared cell           = %16ld\n", fz);
    printf(    "\t\tTarget cell            = %16ld\n", target_count);
    fprintf(fp,"\t\tTarget cell            = %16ld\n", target_count);
  }
  
  
  
  if ( C.Fill_Hint >= 0 ) // ヒントが与えられている場合
  {
    filled = V.fillSeed(d_mid, C.Fill_Hint, C.Fill_Fluid, d_cut);
    
    if ( numProc > 1 )
    {
      unsigned long tmp_fs = filled;
      if ( paraMngr->Allreduce(&tmp_fs, &filled, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
    }
    
    if ( filled == 0 )
    {
      Hostonly_
      {
        printf(    "Failed second painting (%s includes solid cell)\n", FBUtility::getDirection(C.Fill_Hint).c_str());
        fprintf(fp,"Failed second painting (%s includes solid cell)\n", FBUtility::getDirection(C.Fill_Hint).c_str());
      }
      Exit(0);
    }
    
    Hostonly_
    {
      printf(    "\t\tPainted cells          = %16ld\n", filled);
      fprintf(fp,"\t\tPainted cells          = %16ld\n", filled);
    }
  }
  
  // ペイントされたシードセル数をターゲットから差し引く
  target_count -= filled;
  
  Hostonly_
  {
    printf(    "\t\tRemaining target cells = %16ld\n\n", target_count);
    fprintf(fp,"\t\tRemaining target cells = %16ld\n\n", target_count);
  }
  
  
  // BIDによるフィル
  c = -1;
  
  while (target_count > 0) {
    
    unsigned fs;
    filled = (unsigned long)V.fill_by_bid(d_bid, d_mid, d_cut, C.Fill_Fluid, fs, list);
    replaced = (unsigned long)fs;
    
    if ( numProc > 1 )
    {
      unsigned long t_fc = filled;
      if ( paraMngr->Allreduce(&t_fc, &filled, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
      
      t_fc = replaced;
      if ( paraMngr->Allreduce(&t_fc, &replaced, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
    }
    
    if ( filled == 0 ) break; // フィル対象がなくなったら終了
    
    target_count -= filled;
    target_count -= replaced;
    c++;
    
    // 同期
    if ( numProc > 1 )
    {
      if ( paraMngr->BndCommS4DEx(d_cut, 6, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
      if ( paraMngr->BndCommS3D(d_mid, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
      if ( paraMngr->BndCommS3D(d_bid, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
    }
  }
  
  Hostonly_
  {
    printf(    "\t\tBID Iteration          = %5d\n", c+1);
    fprintf(fp,"\t\tBID Iteration          = %5d\n", c+1);
    printf(    "\t\t\tFLUID filled   = %16ld\n", filled);
    fprintf(fp,"\t\t\tFLUID filled   = %16ld\n", filled);
    printf(    "\t\t\tSOLID replaced = %16ld\n", replaced);
    fprintf(fp,"\t\t\tSOLID replaced = %16ld\n", replaced);
    printf(    "\t\t\tRemaining cell = %16ld\n\n", target_count);
    fprintf(fp,"\t\t\tRemaining cell = %16ld\n\n", target_count);
  }
  
  
  
  // 固体に変更
  // Allreduce時の桁あふれ対策のため、unsigned long で集約
  c = -1;
  while ( target_count > 0 ) {
    
    // 未ペイントのセルに対して、固体IDを与える
    replaced = (unsigned long)V.fillReplace(d_mid, 0, C.Fill_Solid);
    
    
    if ( numProc > 1 )
    {
      unsigned long t_fc = replaced;
      if ( paraMngr->Allreduce(&t_fc, &replaced, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
    }
    
    if ( replaced == 0 ) break;
    
    target_count -= replaced;
    c++;
  }
  
  Hostonly_
  {
    printf(    "\t\tSOLID filling Iteration= %5d\n", c+1);
    fprintf(fp,"\t\tSOLID filling Iteration= %5d\n", c+1);
    printf(    "\t\t\tFilled by SOLID= %16ld\n\n", target_count);
    fprintf(fp,"\t\t\tFilled by SOLID= %16ld\n\n", target_count);
  }
  

  
  // 確認 paintedは未ペイントセルがある場合に1
  unsigned long painted = (unsigned long)V.fill_check(d_mid);
  
  if ( numProc > 1 )
  {
    unsigned long t_painted = painted;
    if ( paraMngr->Allreduce(&t_painted, &painted, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
  }
  
  if ( painted  != 0 )
  {
    Hostonly_
    {
      printf(    "\tFill operation is done, but still remains unpainted cells.\n");
      fprintf(fp,"\tFill operation is done, but still remains unpainted cells.\n");
    }
    Exit(0);
  }

    if ( list ) delete [] list;
  
}



// #################################################################
/* @brief 固定パラメータの設定
 */
void FFV::fixedParameters()
{
  // 精度
  if ( sizeof(REAL_TYPE) == sizeof(double) )
  {
    C.Mode.Precision = sizeof(double);
  }
  else
  {
    C.Mode.Precision = sizeof(float);
  }
  
  // ログファイル名
  C.HistoryName        = "history_base.txt";
  C.HistoryCompoName   = "history_compo.txt";
  C.HistoryDomfxName   = "history_domainflux.txt";
  C.HistoryForceName   = "history_force.txt";
  C.HistoryWallName    = "history_log_wall.txt";
  C.HistoryItrName     = "history_iteration.txt";
  C.HistoryMonitorName = "sampling.txt";
  
  C.f_Pressure       = "prs";
  C.f_Velocity       = "vel";
  C.f_Fvelocity      = "fvel";
  C.f_Temperature    = "tmp";
  C.f_AvrPressure    = "prsa";
  C.f_AvrVelocity    = "vela";
  C.f_AvrTemperature = "tmpa";
  C.f_DivDebug       = "div";
  C.f_Helicity       = "hlt";
  C.f_TotalP         = "tp";
  C.f_I2VGT          = "i2vgt";
  C.f_Vorticity      = "vrt";
}



// #################################################################
/* @brief 並列処理時の各ノードの分割数を集めてファイルに保存する
 */
void FFV::gatherDomainInfo()
{
  // 統計処理の母数
  double d = 1.0 /(double)numProc;
  double r;
  
  if ( numProc > 1 )
  {
    r = 1.0 /(double)(numProc-1);
  }
  else
  {
    r = 1.0;
  }
  
  int* m_size = NULL;           ///< 領域分割数
  unsigned long* bf_fcl = NULL; ///< Fluid cell
  unsigned long* bf_wcl = NULL; ///< Solid cell
  unsigned long* bf_acl = NULL; ///< Active cell
  
  REAL_TYPE* m_org = NULL;      ///< 基点
  REAL_TYPE* m_reg = NULL;      ///< 領域サイズ
  double* bf_srf = NULL;        ///< 表面数
  
  int* st_buf = NULL; ///<
  int* ed_buf = NULL; ///<
  
  
  if( !(m_size = new int[numProc*3]) ) Exit(0);
  if( !(bf_fcl = new unsigned long[numProc]) )   Exit(0);
  if( !(bf_wcl = new unsigned long[numProc]) )   Exit(0);
  if( !(bf_acl = new unsigned long[numProc]) )   Exit(0);
  
  if( !(m_org  = new REAL_TYPE[numProc*3]) ) Exit(0);
  if( !(m_reg  = new REAL_TYPE[numProc*3]) ) Exit(0);
  if( !(bf_srf = new double[numProc]) )   Exit(0);
  
  if( !(st_buf = new int [numProc*3]) ) Exit(0);
  if( !(ed_buf = new int [numProc*3]) ) Exit(0);
  
  // 領域情報の収集
  if ( numProc > 1 )
  {
    if ( paraMngr->Gather(size, 3, m_size, 3, 0) != CPM_SUCCESS ) Exit(0);
    if ( paraMngr->Gather(origin, 3, m_org, 3, 0) != CPM_SUCCESS ) Exit(0);
    if ( paraMngr->Gather(region, 3, m_reg, 3, 0) != CPM_SUCCESS ) Exit(0);
    if ( paraMngr->Gather(&L_Fcell, 1, bf_fcl, 1, 0) != CPM_SUCCESS ) Exit(0);
    if ( paraMngr->Gather(&L_Wcell, 1, bf_wcl, 1, 0) != CPM_SUCCESS ) Exit(0);
    if ( paraMngr->Gather(&L_Acell, 1, bf_acl, 1, 0) != CPM_SUCCESS ) Exit(0);
  }
  else // serial
  {
    memcpy(m_size, size, 3*sizeof(int));
    bf_fcl[0] = L_Fcell;
    bf_wcl[0] = L_Wcell;
    bf_acl[0] = L_Acell;
    memcpy(m_org, origin, 3*sizeof(REAL_TYPE));
    memcpy(m_reg, region, 3*sizeof(REAL_TYPE));
  }
  
  
  // Info. of computational domain
  double vol = (double)( (double)G_size[0] * (double)G_size[1] * (double)G_size[2]);
  double srf = 2.0   * ( (double)G_size[0] * (double)G_size[1] 
                       + (double)G_size[1] * (double)G_size[2] 
                       + (double)G_size[2] * (double)G_size[0]);
  
  // ローカルノード
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  double m_srf = (double)((ix*jx + jx*kx + kx*ix)) * 2.0;
  
  if ( nID[X_MINUS] < 0 ) m_srf -= (double)(jx*kx);  // remove face which does not join communication
  if ( nID[Y_MINUS] < 0 ) m_srf -= (double)(ix*kx);
  if ( nID[Z_MINUS] < 0 ) m_srf -= (double)(ix*jx);
  if ( nID[X_PLUS]  < 0 ) m_srf -= (double)(jx*kx);
  if ( nID[Y_PLUS]  < 0 ) m_srf -= (double)(ix*kx);
  if ( nID[Z_PLUS]  < 0 ) m_srf -= (double)(ix*jx);
  
  if ( numProc > 1 )
  {
    if ( paraMngr->Gather(&m_srf, 1, bf_srf, 1, 0) != CPM_SUCCESS ) Exit(0);
  }
  else 
  {
    bf_srf[0] = m_srf;
  }
  
  // mean of domain
  double m_vol = 0;
  double m_efv = 0;
  m_srf = 0;
  
  for (int i=0; i<numProc; i++) {
    ix = m_size[3*i];
    jx = m_size[3*i+1];
    kx = m_size[3*i+2];
    m_vol += (double)(ix*jx*kx);
    m_srf += bf_srf[i];
    m_efv += (double)bf_acl[i];
  }
  
  double d_vol = m_vol * d;
  double d_srf = m_srf * d;
  double d_efv = m_efv * d;
  
  // std. deviation of domain
  double vol_dv = 0.0;
  double srf_dv = 0.0;
  double efv_dv = 0.0;
  
  double d1, d2, d3;
  
  for (int i=0; i<numProc; i++) {
    ix = m_size[3*i];
    jx = m_size[3*i+1];
    kx = m_size[3*i+2];
    d1 = (double)(ix*jx*kx) - d_vol;
    d2 = bf_srf[i] - d_srf;
    d3 = (double)bf_acl[i] - d_efv;
    vol_dv += d1 * d1;
    srf_dv += d2 * d2;
    efv_dv += d3 * d3;
  }
  vol_dv = sqrt(vol_dv*r);
  srf_dv = sqrt(srf_dv*r);
  efv_dv = sqrt(efv_dv*r);
  
  
  FILE *fp=NULL;
  if ( !(fp=fopen("DomainInfo.txt", "w")) )
  {
    stamped_printf("\tSorry, can't open 'DomainInfo.txt' file. Write failed.\n");
    Exit(0);
  }
  
  // 全体情報の表示
  C.printGlobalDomain(fp, G_size, G_origin, G_region, pitch);
  
  // ローカルノードの情報を表示
  for (int i=0; i<numProc; i++)
  {
    Hostonly_
    {
      fprintf(fp,"Domain %4d\n", i);
      fprintf(fp,"\t ix, jx,  kx        [-] =  %13ld %13ld %13ld\n",  m_size[i*3], m_size[i*3+1], m_size[i*3+2]);
      fprintf(fp,"\t(ox, oy, oz)  [m] / [-] = (%13.6e %13.6e %13.6e)  /  (%13.6e %13.6e %13.6e)\n", 
              m_org[i*3]*C.RefLength,  m_org[i*3+1]*C.RefLength,  m_org[i*3+2]*C.RefLength, m_org[i*3],  m_org[i*3+1],  m_org[i*3+2]);
      fprintf(fp,"\t(Lx, Ly, Lz)  [m] / [-] = (%13.6e %13.6e %13.6e)  /  (%13.6e %13.6e %13.6e)\n", 
              m_reg[i*3]*C.RefLength,  m_reg[i*3+1]*C.RefLength,  m_reg[i*3+2]*C.RefLength, m_reg[i*3],  m_reg[i*3+1],  m_reg[i*3+2]);
      
      if (C.NoCompo != 0) fprintf(fp, "\t no            Label   i_st    i_ed    j_st    j_ed    k_st    k_ed\n");
    }
    
    if ( numProc > 1 )
    {
      for (int n=1; n<=C.NoCompo; n++)
      {
        if ( !cmp[n].isKindMedium() )
        {
          if( paraMngr->Gather(cmp[n].getBbox_st(), 3, st_buf, 3, 0) != CPM_SUCCESS ) Exit(0);
          if( paraMngr->Gather(cmp[n].getBbox_ed(), 3, ed_buf, 3, 0) != CPM_SUCCESS ) Exit(0);
        }
        Hostonly_
        {
          fprintf(fp,"\t%3d %16s %7d %7d %7d %7d %7d %7d\n",
                  n, cmp[n].getAlias().c_str(), st_buf[i*3], ed_buf[i*3], st_buf[i*3+1], ed_buf[i*3+1], st_buf[i*3+2], ed_buf[i*3+2]);
        }
      }
    }
    else // serial
    {
      int *st, *ed;
      for (int n=1; n<=C.NoCompo; n++)
      {
        if ( !cmp[n].isKindMedium() )
        {
          st = cmp[n].getBbox_st();
          ed = cmp[n].getBbox_ed();
        
          Hostonly_
          {
            fprintf(fp,"\t%3d %16s %7d %7d %7d %7d %7d %7d\n",
                  n, cmp[n].getAlias().c_str(), st[0], ed[0], st[1], ed[1], st[2], ed[2]);
          }
        }
      }
    }
  }
  
  Hostonly_
  {
    fprintf(fp, "\n");
    fprintf(fp,"\n\t--------------------------------------------------\n");
    fprintf(fp,"\tReport of Whole Domain Statistics\n");
    fprintf(fp,"\tDomain size         = %7d %7d %7d\n", G_size[0], G_size[1], G_size[2]);
    fprintf(fp,"\tNumber of voxels    = %12.6e\n", vol);
    fprintf(fp,"\tNumber of surface   = %12.6e\n", srf);
    fprintf(fp,"\tEffective voxels    = %12.6e (%6.2f%%)\n", (REAL_TYPE)G_Acell, 100.0*(REAL_TYPE)G_Acell/vol);
    fprintf(fp,"\tFluid voxels        = %12.6e (%6.2f%%)\n", (REAL_TYPE)G_Fcell, 100.0*(REAL_TYPE)G_Fcell/vol);
    fprintf(fp,"\tWall  voxels        = %12.6e (%6.2f%%)\n", (REAL_TYPE)G_Wcell, 100.0*(REAL_TYPE)G_Wcell/vol);
    
    if ( numProc == 1 )
    {
      fprintf(fp,"\tDivision :          = %d : %s\n", numProc, "Serial");
    }
    else
    {
      fprintf(fp,"\tDivision :          = %d : %s\n", numProc, "Equal segregation");
    }
    
    fprintf(fp,"\n\t--------------------------------------------------\n");
    fprintf(fp,"\tDomain Statistics per MPI process\n");
    fprintf(fp,"\tMean volume in each domain           = %12.6e\n", d_vol);
    fprintf(fp,"\tStd. deviation of domain             = %12.6e\n", vol_dv);
    fprintf(fp,"\tMean comm. in each domain            = %12.6e\n", d_srf);
    fprintf(fp,"\tStd. deviation of surface            = %12.6e\n", srf_dv);
    fprintf(fp,"\tMean effective volume in each domain = %12.6e\n", d_efv);
    fprintf(fp,"\tStd. deviation of effective volume   = %12.6e\n", efv_dv);
    fprintf(fp,"\n");
    
    fprintf(fp,"\tDomain :     ix     jx     kx       Volume Vol_dv[%%]      Surface Srf_dv[%%] Fluid[%%] Solid[%%]      Eff_Vol Eff_Vol_dv[%%]      Eff_Srf Eff_srf_dv[%%]  Itr_scheme\n");
    fprintf(fp,"\t---------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
    
    double tmp_vol, tmp_acl, tmp_fcl, tmp_wcl;
    
    for (int i=0; i<numProc; i++)
    {
      ix = m_size[3*i];
      jx = m_size[3*i+1];
      kx = m_size[3*i+2];
      tmp_vol = (double)(ix*jx*kx);
      tmp_acl = (double)bf_acl[i];
      tmp_fcl = (double)bf_fcl[i];
      tmp_wcl = (double)bf_wcl[i];
      fprintf(fp,"\t%6d : %6d %6d %6d ", i, ix, jx, kx);
      fprintf(fp,"%12.4e  %8.3f ", tmp_vol, 100.0*(tmp_vol-d_vol)/d_vol);
      fprintf(fp,"%12.4e  %8.3f ", bf_srf[i], (d_srf == 0.0) ? 0.0 : 100.0*(bf_srf[i]-d_srf)/d_srf);
      fprintf(fp,"%8.3f %8.3f ", 100.0*tmp_fcl/tmp_vol, 100.0*tmp_wcl/tmp_vol);
      fprintf(fp,"%12.4e      %8.3f \n", tmp_acl, 100.0*(tmp_acl-d_efv)/d_efv);
    }
    fprintf(fp,"\t---------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
  }
  
  if (fp) fclose(fp);
  
  
  if( m_size ) { delete [] m_size; m_size=NULL; }
  if( m_org  ) { delete [] m_org;  m_org =NULL; }
  if( m_reg  ) { delete [] m_reg;  m_reg =NULL; }
  if( bf_srf ) { delete [] bf_srf; bf_srf=NULL; }
  if( bf_fcl ) { delete [] bf_fcl; bf_fcl=NULL; }
  if( bf_wcl ) { delete [] bf_wcl; bf_wcl=NULL; }
  if( bf_acl ) { delete [] bf_acl; bf_acl=NULL; }
  if( st_buf ) { delete [] st_buf; st_buf=NULL; }
  if( ed_buf ) { delete [] ed_buf; ed_buf=NULL; }
}



// #################################################################
/* @brief binaryの場合に，非BCポリゴンからSOLIDセルを生成
 * @param [in] fp stdout
 */
void FFV::generateSolid(FILE* fp)
{
  unsigned long zc=0;
  

  zc += V.SolidFromCut(d_mid, d_bid, d_cut, cmp);
      
  Hostonly_
  {
    printf(    "\n\tGenerated Solid cell from cut = %ld\n", zc);
    fprintf(fp,"\n\tGenerated Solid cell from cut = %ld\n", zc);
  }
  
  // midのガイドセル同期
  if ( paraMngr->BndCommS3D(d_mid, size[0], size[1], size[2], guide, 1) != CPM_SUCCESS ) Exit(0);
  
}



// #################################################################
// @brief コンポーネントの面積を計算
void FFV::getCompoArea()
{

  if ( C.NoBC == 0 ) return;
  if ( C.Mode.Example != id_Polygon ) return;
  
  
  float area=0.0;
  
  // コンポーネントで指定されるID面の法線を計算
  for (int n=1; n<=C.NoCompo; n++)
  {
    if ( cmp[n].isKindCompo() )
    {
      int type = cmp[n].getType();
      
      if ( type > OBSTACLE )
      {
        string label = cmp[n].getAlias();
        
        for (int i=0; i<C.num_of_polygrp; i++)
        {
          if ( FBUtility::compare(PolyPP[i].getGroup(), label) )
          {
            cmp[n].area = PolyPP[i].getGarea();
          }
        }
      }
    }

  }
}


// #################################################################
/* @brief グローバルな領域情報を取得
 * @param [in] tp_dom  TPControlクラス
 * @return 分割指示 (1-with / 2-without)
 */
int FFV::getDomainInfo(TPControl* tp_dom)
{
  // 領域分割モードのパターン
  //     分割指定(G_div指定)         domain.txt
  // 1)  G_div指定なし |  G_orign + G_region + (G_pitch || G_voxel)
  // 2)  G_div指定あり |  G_orign + G_region + (G_pitch || G_voxel)
  // 3)  G_div指定なし |          + ActiveDomainInfo
  // 4)  G_div指定あり |          + ActiveDomainInfo
  
  /* Policy
   * G_regionは必須
   * G_voxelとG_pitchでは，G_voxelが優先．両方が指定されている場合はエラー
   * G_pitchが指定されており，割り切れない場合は，領域を拡大
   *        G_voxel = (int)ceil(G_region/G_pitch)
   *        G_region = G_pitch * G_voxel
   */
  
  string label, str;
  int div_type = 1; // 指定分割 => 1
  
  
  // 長さの単位
  label = "/DomainInfo/UnitOfLength";
  
  if ( !tp_dom->GetValue(label, &str) )
  {
		Hostonly_ stamped_printf("\tParsing error : Invalid string for '%s'\n", label.c_str());
	  Exit(0);
  }
  
  if     ( !strcasecmp(str.c_str(), "NonDimensional") )  C.Unit.Length = LTH_ND;
  else if( !strcasecmp(str.c_str(), "M") )               C.Unit.Length = LTH_m;
  else if( !strcasecmp(str.c_str(), "cm") )              C.Unit.Length = LTH_cm;
  else if( !strcasecmp(str.c_str(), "mm") )              C.Unit.Length = LTH_mm;
  else
  {
    Hostonly_ stamped_printf("\tInvalid keyword is described at '%s'\n", label.c_str());
    Exit(0);
  }
  
  
  // G_origin　必須
  label = "/DomainInfo/GlobalOrigin";
  
  if ( !tp_dom->GetVector(label, G_origin, 3) )
  {
    cout << "ERROR : in parsing [" << label << "]" << endl;
    Exit(0);
  }
  
  // G_region 必須
  label = "/DomainInfo/GlobalRegion";
  
  if ( !tp_dom->GetVector(label, G_region, 3) )
  {
    Hostonly_ cout << "ERROR : in parsing [" << label << "]" << endl;
    Exit(0);
  }
  
  if ( (G_region[0]>0.0) && (G_region[1]>0.0) && (G_region[2]>0.0) )
  {
    ; // skip
  }
  else
  {
    Hostonly_ cout << "ERROR : in parsing [" << label << "]" << endl;
    Exit(0);
  }
  
  // 排他チェックフラグ
  bool g_flag = true;
  bool p_flag = true;
  
  
  // G_voxel
  label = "/DomainInfo/GlobalVoxel";
  
  if ( !tp_dom->GetVector(label, G_size, 3) ) g_flag = false;
  
  
  // 2D check
  if ( (Ex->mode == Intrinsic::dim_2d) && (G_size[2] != 3) )
  {
    Hostonly_ {
      printf("\tError : In case of 2 dimensional problem, kmax must be 3.\n");
      Exit(0);
    }
  }
  
  // 偶数のチェック
  if ( Ex->even == ON )
  {
    if ( G_size[0]/2*2 != G_size[0] )
    {
      printf("\tDimension size must be even for x direction (%d %d %d)\n", G_size[0], G_size[1], G_size[2]);
      Exit(0);
    }
    if ( G_size[1]/2*2 != G_size[1] )
    {
      printf("\tDimension size must be even for y direction (%d %d %d)\n", G_size[0], G_size[1], G_size[2]);
      Exit(0);
    }
    if ( (Ex->mode == Intrinsic::dim_3d) && (G_size[2]/2*2 != G_size[2]) )
    {
      printf("\tDimension size must be even for z direction (%d %d %d)\n", G_size[0], G_size[1], G_size[2]);
      Exit(0);
    }
  }
  
  
  // pitch
  label = "/DomainInfo/GlobalPitch";
  
  if ( !tp_dom->GetVector(label, pitch, 3) ) p_flag = false;
  
  
  // 排他チェック
  if ( !g_flag && !p_flag ) // 両方とも指定されていない
  {
    Hostonly_ cout << "\tNeither GlobalVoxel nor GlobalPitch is specified. You need to choose either parameter.\n" << endl;
    Exit(0);
  }
  else if ( g_flag ) // G_voxelの指定がある >> G_voxelが優先
  {
    if ( g_flag )
    {
      if ( (G_size[0]>0) && (G_size[1]>0) && (G_size[2]>0) )
      {
        
        pitch[0] = G_region[0] / (REAL_TYPE)G_size[0];
        pitch[1] = G_region[1] / (REAL_TYPE)G_size[1];
        pitch[2] = G_region[2] / (REAL_TYPE)G_size[2];
        
        // 等方性チェック
        if ( Ex->mode == Intrinsic::dim_3d )
        {
          if ( !( (pitch[0] == pitch[1]) && (pitch[1] == pitch[2]) ) )
          {
            Hostonly_ printf("\tGlobal Pitch must be same in all direction (%14.6e, %14.6e, %14.6e)\n", pitch[0], pitch[1], pitch[2]);
            Exit(0);
          }
        }
        else // dim_2d
        {
          if ( pitch[0] != pitch[1] )
          {
            Hostonly_ printf("\tGlobal Pitch must be same in X-Y direction (%14.6e, %14.6e)\n", pitch[0], pitch[1]);
            Exit(0);
          }
          pitch[2] = pitch[0];
          G_region[2] = (REAL_TYPE)G_size[2] * pitch[2];
          G_origin[2] = -0.5 * pitch[2];
        }
        
        
      }
      else
      {
        Hostonly_ printf("ERROR : in parsing [%s] >> (%d, %d, %d)\n", label.c_str(), G_size[0], G_size[1], G_size[2] );
        Exit(0);
      }
    }
  }
  else if ( !g_flag && p_flag ) // pitchのみ
  {
    if ( (pitch[0]>0.0) && (pitch[1]>0.0) && (pitch[2]>0.0) )
    {
      
      // 等方性チェック
      if ( !( (pitch[0] == pitch[1]) && (pitch[1] == pitch[2]) ) )
      {
        Hostonly_ printf("\tGlobal Pitch must be same in all direction (%14.6e, %14.6e, %14.6e)\n", pitch[0], pitch[1], pitch[2]);
        Exit(0);
      }
      
      // pitchを基準にして、全計算領域の分割数を計算する。切り上げ
      G_size[0] = (int)ceil(G_region[0]/pitch[0]);
      G_size[1] = (int)ceil(G_region[1]/pitch[1]);
      G_size[2] = (int)ceil(G_region[2]/pitch[2]);
      
      // 再計算された全計算領域の分割数から、全計算領域の大きさを再計算
      REAL_TYPE gr[3];
      gr[0] = (REAL_TYPE)G_size[0] * pitch[0];
      gr[1] = (REAL_TYPE)G_size[1] * pitch[1];
      gr[2] = (REAL_TYPE)G_size[2] * pitch[2];
      
      // 整合性チェック
      if ( Ex->mode == Intrinsic::dim_3d )
      {
        if ( (G_region[0] != gr[0]) || (G_region[1] != gr[1]) || (G_region[2] != gr[2]) )
        {
          Hostonly_ {
            printf("\tGlobal Region is modified due to maintain the consistency of domain parameters.\n\n");
            printf("\t[%12.6e  %12.6e  %12.6e] >> [%12.6e  %12.6e  %12.6e]\n\n",
                   G_region[0], G_region[1], G_region[2], gr[0], gr[1], gr[2]);
          }
          G_region[0] = gr[0];
          G_region[1] = gr[1];
          G_region[2] = gr[2];
        }
      }
      else // dim_2d
      {
        if ( (G_region[0] != gr[0]) || (G_region[1] != gr[1]) )
        {
          Hostonly_ {
            printf("\tGlobal Region is modified due to maintain the consistency of domain parameters.\n\n");
            printf("\t[%12.6e  %12.6e  %12.6e] >> [%12.6e  %12.6e  %12.6e]\n\n",
                   G_region[0], G_region[1], G_region[2], gr[0], gr[1], gr[2]);
          }
          G_region[0] = gr[0];
          G_region[1] = gr[1];
          G_region[2] = gr[2];
          G_origin[2] = -0.5 * pitch[2];
        }
      }
      
    }
    else // パラメータが無効の場合
    {
      Hostonly_ printf("ERROR : in parsing [%s] >> (%e, %e, %e)\n", label.c_str(), pitch[0], pitch[1], pitch[2] );
      Exit(0);
    }
  }
  else
  {
    Hostonly_ cout << "\tError : unknown\n" << endl;
    Exit(0);
  }
  

  
  
  // G_division オプション
  label = "/DomainInfo/GlobalDivision";
  
  if ( !tp_dom->GetVector(label, G_division, 3) )
  {
    Hostonly_ cout << "\tAutomatic domain division is selected." << endl;
    div_type = 2; // 自動分割
  }
  
  // プロセス分割数が指定されている場合のチェック
  if ( div_type == 1 )
  {
    if ( (G_division[0]>0) && (G_division[1]>0) && (G_division[2]>0) ) 
    {
      Hostonly_ printf("\tManual domain division is selected.\n");
    }
    else
    {
      Hostonly_ cout << "ERROR : in parsing [" << label << "]" << endl;
      Exit(0);
    }
  }

  
  // 有次元の場合には、単位をメートルに変換
  REAL_TYPE factor;
  switch ( C.Unit.Length )
  {
    case LTH_m:
    case LTH_ND:
      factor = 1.0; // 変更なし
      break;
      
    case LTH_cm:
      factor = 0.01;
      break;
      
    case LTH_mm:
      factor = 0.001;
      break;
      
    default:
      Exit(0);
  }
  
  for (int i=0; i<3; i++) {
    pitch[i]    *= factor;
    G_origin[i] *= factor;
    G_region[i] *= factor;
  }
  
  
  // ActiveSubdomainファイル名の取得
  label = "/DomainInfo/ActiveSubDomainFile";
  
  if ( !tp_dom->GetValue(label, &str ) )
  {
    Hostonly_ cout << "\tNo option : in parsing [" << label << "]" << endl;
  }
  //@todo  string hoge = str;
  
  
  return div_type;
}


// #################################################################
/* @brief Intrisic classの同定
 * @param [in] tpf ffvのパラメータを保持するTextParserのインスタンス
 * @param [in] fp  ファイル出力ポインタ
 */
void FFV::identifyExample(TPControl* tpf, FILE* fp)
{
  
  // 例題の種類を取得
  string keyword;
  string label;
  
  label = "/Steer/Example";
  
  if ( !(tpf->GetValue(label, &keyword )) )
  {
    Hostonly_ stamped_printf("\tError : '%s'\n", label.c_str());
    Exit(0);
  }
  
  if     ( FBUtility::compare(keyword, "ParallelPlate2D") )   C.Mode.Example = id_PPLT2D;
  else if( FBUtility::compare(keyword, "Duct") )              C.Mode.Example = id_Duct;
  else if( FBUtility::compare(keyword, "SHC1D") )             C.Mode.Example = id_SHC1D;
  else if( FBUtility::compare(keyword, "PerformanceTest") )   C.Mode.Example = id_PMT;
  else if( FBUtility::compare(keyword, "Rectangular") )       C.Mode.Example = id_Rect;
  else if( FBUtility::compare(keyword, "Cylinder") )          C.Mode.Example = id_Cylinder;
  else if( FBUtility::compare(keyword, "BackStep") )          C.Mode.Example = id_Step;
  else if( FBUtility::compare(keyword, "Polygon") )           C.Mode.Example = id_Polygon;
  else if( FBUtility::compare(keyword, "Sphere") )            C.Mode.Example = id_Sphere;
  else if( FBUtility::compare(keyword, "Jet") )               C.Mode.Example = id_Jet;
  else
  {
    Hostonly_
    {
      stamped_printf(     "\tInvalid keyword is described for '%s'\n", label.c_str());
      stamped_fprintf(fp, "\tInvalid keyword is described for '%s'\n", label.c_str());
    }
    Exit(0);
  }
  
  
  // 例題クラスの実体をインスタンスし，Exにポイントする
  if      ( C.Mode.Example == id_PPLT2D)   Ex = dynamic_cast<Intrinsic*>(new IP_PPLT2D);
  else if ( C.Mode.Example == id_Duct )    Ex = dynamic_cast<Intrinsic*>(new IP_Duct);
  else if ( C.Mode.Example == id_SHC1D)    Ex = dynamic_cast<Intrinsic*>(new IP_SHC1D);
  else if ( C.Mode.Example == id_PMT )     Ex = dynamic_cast<Intrinsic*>(new IP_PMT);
  else if ( C.Mode.Example == id_Rect )    Ex = dynamic_cast<Intrinsic*>(new IP_Rect);
  else if ( C.Mode.Example == id_Cylinder) Ex = dynamic_cast<Intrinsic*>(new IP_Cylinder);
  else if ( C.Mode.Example == id_Step )    Ex = dynamic_cast<Intrinsic*>(new IP_Step);
  else if ( C.Mode.Example == id_Polygon ) Ex = dynamic_cast<Intrinsic*>(new IP_Polygon);
  else if ( C.Mode.Example == id_Sphere )  Ex = dynamic_cast<Intrinsic*>(new IP_Sphere);
  else if ( C.Mode.Example == id_Jet )     Ex = dynamic_cast<Intrinsic*>(new IP_Jet);
  else
  {
    Hostonly_
    {
      stamped_printf(    "\tInvalid keyword is described for Exmple definition\n");
      stamped_fprintf(fp,"\tInvalid keyword is described for Exmple definition\n");
    }
    Exit(0);
  }
  
  
  // 組み込み例題クラス名を表示
  Hostonly_ Ex->printExample(fp, C.Mode.Example);
  
}


// #################################################################
/* @brief ファイル出力の初期化
 */
void FFV::init_FileOut()
{
  std::string hostname;
  hostname = paraMngr->GetHostName();
  
  
  // Format
  cio_DFI::E_CIO_FORMAT format;
  
  if ( C.FIO.Format == sph_fmt )
  {
    format = cio_DFI::E_CIO_FMT_SPH;
  }
  else if ( C.FIO.Format == bov_fmt )
  {
    format = cio_DFI::E_CIO_FMT_BOV;
  }
  
  // Datatype
  cio_DFI::E_CIO_DTYPE datatype;
  
  if ( sizeof(REAL_TYPE) == 4 )
  {
    datatype = cio_DFI::E_CIO_FLOAT32;
  }
  else if ( sizeof(REAL_TYPE) == 8 )
  {
    datatype = cio_DFI::E_CIO_FLOAT64;
  }
  else
  {
    Exit(0);
  }
  
  
  // 出力ファイルヘッダ
  int cio_tail[3], cio_div[3];
  for (int i=0; i<3; i++) cio_tail[i]=size[i];
  for (int i=0; i<3; i++) cio_div[i]=1;
  
  if ( numProc > 1)
  {
    const int* p_tail = paraMngr->GetVoxelTailIndex();
    for (int i=0; i<3; i++ ) cio_tail[i]=p_tail[i]+1;
    
    const int* p_div = paraMngr->GetDivNum();
    for (int i=0; i<3; i++ ) cio_div[i] = p_div[i];
  }
  

  int gc_out = C.GuideOut;
  REAL_TYPE cio_org[3], cio_pit[3];
  
  for (int i=0; i<3; i++)
  {
    cio_org[i] = origin[i] - pitch[i]*(REAL_TYPE)C.GuideOut; // ガイドセルによるオリジナルポイントの調整
    cio_pit[i] = pitch[i];
  }
  
  // セルセンター位置を基点とする
  for (int i=0; i<3; i++)
  {
    cio_org[i] += 0.5*cio_pit[i];
  }
  
  // 出力ファイルの指定が有次元の場合
  if ( C.Unit.File == DIMENSIONAL )
  {
    for (int i=0; i<3; i++)
    {
      cio_org[i] *= C.RefLength;
      cio_pit[i] *= C.RefLength;
    }
  }
  
  // make output directory
  std::string path = C.FIO.OutDirPath;

  
  // タイムスライス出力オプション
  cio_DFI::E_CIO_ONOFF TimeSliceDir;
  if ( C.FIO.Slice == ON )
  {
    TimeSliceDir = cio_DFI::E_CIO_ON;
  }
  else
  {
    TimeSliceDir = cio_DFI::E_CIO_OFF;
  }
  
  
  std::string process="./proc.dfi";
  std::string dfi_name;
  int comp = 1;
  
  
  // Divergence for Debug
  if ( C.FIO.Div_Debug == ON )
  {
    dfi_name = cio_DFI::Generate_DFI_Name(C.f_dfi_out_div);
    DFI_OUT_DIV = cio_DFI::WriteInit(MPI_COMM_WORLD,
                                     dfi_name,
                                     path,
                                     C.f_DivDebug,
                                     format,
                                     gc_out,
                                     datatype,
                                     cio_DFI::E_CIO_IJKN,
                                     comp,
                                     process,
                                     size,
                                     cio_pit,
                                     cio_org,
                                     cio_div,
                                     head,
                                     cio_tail,
                                     hostname,
                                     TimeSliceDir);
    
    if ( DFI_OUT_DIV == NULL )
    {
      Hostonly_ stamped_printf("\tFails to instance Divergence dfi.\n");
      Exit(0);
    }
  }
  

  
  // Pressure
  dfi_name = cio_DFI::Generate_DFI_Name(C.f_dfi_out_prs);
  DFI_OUT_PRS = cio_DFI::WriteInit(MPI_COMM_WORLD,
                                   dfi_name,
                                   path,
                                   C.f_Pressure,
                                   format,
                                   gc_out,
                                   datatype,
                                   cio_DFI::E_CIO_IJKN,
                                   comp,
                                   process,
                                   size,
                                   cio_pit,
                                   cio_org,
                                   cio_div,
                                   head,
                                   cio_tail,
                                   hostname,
                                   TimeSliceDir);
  
  if ( DFI_OUT_PRS == NULL )
  {
    Hostonly_ stamped_printf("\tFails to instance Pressure dfi.\n");
    Exit(0);
  }
  
  if( C.FIO.Slice == ON )
  {
    DFI_OUT_PRS->m_outSlice = true;
  }
  else
  {
    DFI_OUT_PRS->m_outSlice = false;
  }
  
  
  std::string UnitL;
  
  switch (C.Unit.Length)
  {
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
  
  DFI_OUT_PRS->SetUnitLength(true, UnitL, (double)C.RefLength);
  
  
  std::string UnitV;
  if ( C.Unit.File == DIMENSIONAL )
  {
    UnitV = "m/s";
  }
  else
  {
    UnitV = "NonDimensional";
  }
  
  DFI_OUT_PRS->SetUnitVelo(true, UnitV, (double)C.RefVelocity);
  
  
  std::string UnitP;
  if ( C.Unit.File == DIMENSIONAL )
  {
    UnitP = "Pa";
  }
  else
  {
    UnitP = "NonDimensional";
  }
  
  double DiffPrs;
  DiffPrs = (double)C.RefDensity * (double)C.RefVelocity * (double)C.RefVelocity;
  DFI_OUT_PRS->SetUnitPres(true, UnitP, (double)C.BasePrs, DiffPrs);
  
  cio_DFI::MakeDirectory(CIO::cioPath_DirName(dfi_name));
  std::string processdir = CIO::cioPath_DirName(dfi_name)+"/"+CIO::cioPath_FileName(process,".dfi");
  
  DFI_OUT_PRS->WriteProcDfiFile(MPI_COMM_WORLD,
                                processdir,
                                G_size,
                                cio_div,
                                head,
                                cio_tail,
                                G_origin,
                                cio_pit,
                                hostname,
                                true);
  
  
  // Velocity
  dfi_name = cio_DFI::Generate_DFI_Name(C.f_dfi_out_vel);
  comp = 3;
  DFI_OUT_VEL = cio_DFI::WriteInit(MPI_COMM_WORLD,
                                   dfi_name,
                                   path,
                                   C.f_Velocity,
                                   format,
                                   gc_out,
                                   datatype,
                                   cio_DFI::E_CIO_NIJK,
                                   comp,
                                   process,
                                   size,
                                   cio_pit,
                                   cio_org,
                                   cio_div,
                                   head,
                                   cio_tail,
                                   hostname,
                                   TimeSliceDir);
  
  if ( DFI_OUT_VEL == NULL )
  {
    Hostonly_ stamped_printf("\tFails to instance Velocity dfi.\n");
    Exit(0);
  }
  if ( C.FIO.Slice == ON )
  {
    DFI_OUT_VEL->m_outSlice = true;
  }
  else
  {
    DFI_OUT_VEL->m_outSlice = false;
  }
  
  DFI_OUT_VEL->SetUnitLength(true, UnitL, (double)C.RefLength);
  DFI_OUT_VEL->SetUnitVelo(  true, UnitV, (double)C.RefVelocity);
  DFI_OUT_VEL->SetUnitPres(  true, UnitP, (double)C.BasePrs, DiffPrs);
  
  
  // Fvelocity
  if (C.Mode.FaceV == ON )
  {
    dfi_name = cio_DFI::Generate_DFI_Name(C.f_dfi_out_fvel);
    comp = 3;
    DFI_OUT_FVEL = cio_DFI::WriteInit(MPI_COMM_WORLD,
                                      dfi_name,
                                      path,
                                      C.f_Fvelocity,
                                      format,
                                      gc_out,
                                      datatype,
                                      cio_DFI::E_CIO_NIJK,
                                      comp,
                                      process,
                                      size,
                                      cio_pit,
                                      cio_org,
                                      cio_div,
                                      head,
                                      cio_tail,
                                      hostname,
                                      TimeSliceDir);
    
    if( DFI_OUT_FVEL == NULL )
    {
      Hostonly_ stamped_printf("\tFails to instance Fvelocity dfi.\n");
      Exit(0);
    }
    if ( C.FIO.Slice == ON )
    {
      DFI_OUT_FVEL->m_outSlice = true;
    }
    else
    {
      DFI_OUT_FVEL->m_outSlice = false;
    }
    
    DFI_OUT_FVEL->SetUnitLength(true, UnitL, (double)C.RefLength);
    DFI_OUT_FVEL->SetUnitVelo(  true, UnitV, (double)C.RefVelocity);
    DFI_OUT_FVEL->SetUnitPres(  true, UnitP, (double)C.BasePrs, DiffPrs);
    
  }
  
  // Temperature
  if ( C.isHeatProblem() )
  {
    dfi_name = cio_DFI::Generate_DFI_Name(C.f_dfi_out_temp);
    comp = 1;
    DFI_OUT_TEMP = cio_DFI::WriteInit(MPI_COMM_WORLD,
                                      dfi_name,
                                      path,
                                      C.f_Temperature,
                                      format,
                                      gc_out,
                                      datatype,
                                      cio_DFI::E_CIO_IJKN,
                                      comp,
                                      process,
                                      size,
                                      cio_pit,
                                      cio_org,
                                      cio_div,
                                      head,
                                      cio_tail,
                                      hostname,
                                      TimeSliceDir);
    if( DFI_OUT_TEMP == NULL )
    {
      Hostonly_ stamped_printf("\tFails to instance Temperature dfi.\n");
      Exit(0);
    }
    if ( C.FIO.Slice == ON )
    {
      DFI_OUT_TEMP->m_outSlice = true;
    }
    else
    {
      DFI_OUT_TEMP->m_outSlice = false;
    }
    
    DFI_OUT_TEMP->SetUnitLength(true, UnitL, (double)C.RefLength);
    DFI_OUT_TEMP->SetUnitVelo(  true, UnitV, (double)C.RefVelocity);
    DFI_OUT_TEMP->SetUnitPres(  true, UnitP, (double)C.BasePrs, DiffPrs);
    
  }
  
  // 平均値
  if ( C.Mode.Average == ON ) {
    
    // Pressure
    dfi_name = cio_DFI::Generate_DFI_Name(C.f_dfi_out_prsa);
    comp = 1;
    DFI_OUT_PRSA = cio_DFI::WriteInit(MPI_COMM_WORLD,
                                      dfi_name,
                                      path,
                                      C.f_AvrPressure,
                                      format,
                                      gc_out,
                                      datatype,
                                      cio_DFI::E_CIO_IJKN,
                                      comp,
                                      process,
                                      size,
                                      cio_pit,
                                      cio_org,
                                      cio_div,
                                      head,
                                      cio_tail,
                                      hostname,
                                      TimeSliceDir);
    if( DFI_OUT_PRSA == NULL )
    {
      Hostonly_ stamped_printf("\tFails to instance AvrPressure dfi.\n");
      Exit(0);
    }
    if( C.FIO.Slice == ON )
    {
      DFI_OUT_PRSA->m_outSlice = true;
    }
    else
    {
      DFI_OUT_PRSA->m_outSlice = false;
    }
    
    DFI_OUT_PRSA->SetUnitLength(true, UnitL, (double)C.RefLength);
    DFI_OUT_PRSA->SetUnitVelo(  true, UnitV, (double)C.RefVelocity);
    DFI_OUT_PRSA->SetUnitPres(  true, UnitP, (double)C.BasePrs, DiffPrs);
    
    // Velocity
    dfi_name = cio_DFI::Generate_DFI_Name(C.f_dfi_out_vela);
    comp = 3;
    DFI_OUT_VELA = cio_DFI::WriteInit(MPI_COMM_WORLD,
                                      dfi_name,
                                      path,
                                      C.f_AvrVelocity,
                                      format,
                                      gc_out,
                                      datatype,
                                      cio_DFI::E_CIO_NIJK,
                                      comp,
                                      process,
                                      size,
                                      cio_pit,
                                      cio_org,
                                      cio_div,
                                      head,
                                      cio_tail,
                                      hostname,
                                      TimeSliceDir);
    if( DFI_OUT_VELA == NULL )
    {
      Hostonly_ stamped_printf("\tcan not instance AvrVelocity dfi\n");
      Exit(0);
    }
    if( C.FIO.Slice == ON )
    {
      DFI_OUT_VELA->m_outSlice = true;
    }
    else
    {
      DFI_OUT_VELA->m_outSlice = false;
    }
    
    DFI_OUT_VELA->SetUnitLength(true, UnitL, (double)C.RefLength);
    DFI_OUT_VELA->SetUnitVelo(  true, UnitV, (double)C.RefVelocity);
    DFI_OUT_VELA->SetUnitPres(  true, UnitP, (double)C.BasePrs, DiffPrs);
    
    // Temperature
    if ( C.isHeatProblem() )
    {
      dfi_name = cio_DFI::Generate_DFI_Name(C.f_dfi_out_tempa);
      comp = 1;
      DFI_OUT_TEMPA = cio_DFI::WriteInit(MPI_COMM_WORLD,
                                         dfi_name,
                                         path,
                                         C.f_AvrTemperature,
                                         format,
                                         gc_out,
                                         datatype,
                                         cio_DFI::E_CIO_IJKN,
                                         comp,
                                         process,
                                         size,
                                         cio_pit,
                                         cio_org,
                                         cio_div,
                                         head,
                                         cio_tail,
                                         hostname,
                                         TimeSliceDir);
      if( DFI_OUT_TEMPA == NULL )
      {
        Hostonly_ stamped_printf("\tFails to instance AvrTemperature dfi.\n");
        Exit(0);
      }
      if ( C.FIO.Slice == ON )
      {
        DFI_OUT_TEMPA->m_outSlice = true;
      }
      else
      {
        DFI_OUT_TEMPA->m_outSlice = false;
      }
      
      DFI_OUT_TEMPA->SetUnitLength(true, UnitL, (double)C.RefLength);
      DFI_OUT_TEMPA->SetUnitVelo(  true, UnitV, (double)C.RefVelocity);
      DFI_OUT_TEMPA->SetUnitPres(  true, UnitP, (double)C.BasePrs, DiffPrs);
    }
  }
  
  
  // Total Pressure
  if (C.Mode.TP == ON )
  {
    dfi_name = cio_DFI::Generate_DFI_Name(C.f_dfi_out_tp);
    comp = 1;
    DFI_OUT_TP = cio_DFI::WriteInit(MPI_COMM_WORLD,
                                    dfi_name,
                                    path,
                                    C.f_TotalP,
                                    format,
                                    gc_out,
                                    datatype,
                                    cio_DFI::E_CIO_IJKN,
                                    comp,
                                    process,
                                    size,
                                    cio_pit,
                                    cio_org,
                                    cio_div,
                                    head,
                                    cio_tail,
                                    hostname,
                                    TimeSliceDir);
    if( DFI_OUT_TP == NULL )
    {
      Hostonly_ stamped_printf("\tFails to instance TotalPressure dfi.\n");
      Exit(0);
    }
    if ( C.FIO.Slice == ON )
    {
      DFI_OUT_TP->m_outSlice = true;
    }
    else
    {
      DFI_OUT_TP->m_outSlice = false;
    }
    
    DFI_OUT_TP->SetUnitLength(true, UnitL, (double)C.RefLength);
    DFI_OUT_TP->SetUnitVelo(  true, UnitV, (double)C.RefVelocity);
    DFI_OUT_TP->SetUnitPres(  true, UnitP, (double)C.BasePrs, DiffPrs);
  }
  
  // Vorticity
  if (C.Mode.VRT == ON )
  {
    dfi_name = cio_DFI::Generate_DFI_Name(C.f_dfi_out_vrt);
    comp = 3;
    DFI_OUT_VRT = cio_DFI::WriteInit(MPI_COMM_WORLD,
                                     dfi_name,
                                     path,
                                     C.f_Vorticity,
                                     format,
                                     gc_out,
                                     datatype,
                                     cio_DFI::E_CIO_NIJK,
                                     comp,
                                     process,
                                     size,
                                     cio_pit,
                                     cio_org,
                                     cio_div,
                                     head,
                                     cio_tail,
                                     hostname,
                                     TimeSliceDir);
    if( DFI_OUT_VRT == NULL )
    {
      Hostonly_ stamped_printf("\tFails to instance Vorticity dfi.\n");
      Exit(0);
    }
    if( C.FIO.Slice == ON )
    {
      DFI_OUT_VRT->m_outSlice = true;
    }
    else
    {
      DFI_OUT_VRT->m_outSlice = false;
    }
    
    DFI_OUT_VRT->SetUnitLength(true, UnitL, (double)C.RefLength);
    DFI_OUT_VRT->SetUnitVelo(  true, UnitV, (double)C.RefVelocity);
    DFI_OUT_VRT->SetUnitPres(  true, UnitP, (double)C.BasePrs, DiffPrs);
  }
  
  // 2nd Invariant of VGT
  if (C.Mode.I2VGT == ON )
  {
    dfi_name = cio_DFI::Generate_DFI_Name(C.f_dfi_out_i2vgt);
    comp = 1;
    DFI_OUT_I2VGT = cio_DFI::WriteInit(MPI_COMM_WORLD,
                                       dfi_name,
                                       path,
                                       C.f_I2VGT,
                                       format,
                                       gc_out,
                                       datatype,
                                       cio_DFI::E_CIO_IJKN,
                                       comp,
                                       process,
                                       size,
                                       cio_pit,
                                       cio_org,
                                       cio_div,
                                       head,
                                       cio_tail,
                                       hostname,
                                       TimeSliceDir);
    if( DFI_OUT_I2VGT == NULL )
    {
      Hostonly_ stamped_printf("\tFails to instance 2nd Invariant of VGT dfi.\n");
      Exit(0);
    }
    if ( C.FIO.Slice == ON )
    {
      DFI_OUT_I2VGT->m_outSlice = true;
    }
    else
    {
      DFI_OUT_I2VGT->m_outSlice = false;
    }
    
    DFI_OUT_I2VGT->SetUnitLength(true, UnitL, (double)C.RefLength);
    DFI_OUT_I2VGT->SetUnitVelo(  true, UnitV, (double)C.RefVelocity);
    DFI_OUT_I2VGT->SetUnitPres(  true, UnitP, (double)C.BasePrs, DiffPrs);
  }
  
  // Helicity
  if (C.Mode.Helicity == ON )
  {
    dfi_name = cio_DFI::Generate_DFI_Name(C.f_dfi_out_hlt);
    comp = 1;
    DFI_OUT_HLT = cio_DFI::WriteInit(MPI_COMM_WORLD,
                                     dfi_name,
                                     path,
                                     C.f_Helicity,
                                     format,
                                     gc_out,
                                     datatype,
                                     cio_DFI::E_CIO_IJKN,
                                     comp,
                                     process,
                                     size,
                                     cio_pit,
                                     cio_org,
                                     cio_div,
                                     head,
                                     cio_tail,
                                     hostname,
                                     TimeSliceDir);
    if( DFI_OUT_HLT == NULL )
    {
      Hostonly_ stamped_printf("\tFails to instance Helicity dfi.\n");
      Exit(0);
    }
    if ( C.FIO.Slice == ON )
    {
      DFI_OUT_HLT->m_outSlice = true;
    }
    else
    {
      DFI_OUT_HLT->m_outSlice = false;
    }
    
    DFI_OUT_HLT->SetUnitLength(true, UnitL, (double)C.RefLength);
    DFI_OUT_HLT->SetUnitVelo(  true, UnitV, (double)C.RefVelocity);
    DFI_OUT_HLT->SetUnitPres(  true, UnitP, (double)C.BasePrs, DiffPrs);
  }
  
}



// #################################################################
/* @brief インターバルの初期化
 */
void FFV::init_Interval()
{
  
  // セッションの初期時刻をセット
  for (int i=0; i<Interval_Manager::tg_END; i++) {
    if ( i != Interval_Manager::tg_average)
    {
      C.Interval[i].setStart(Session_StartStep, Session_StartTime);
    }
    else
    {
      // tg_averageの場合は、get_Average_option()で指定済み
    }
  }
  
  // debug
  //printf("tg_instant start_step : %d\n",C.Interval[Interval_Manager::tg_instant].getStartStep());
  //printf("tg_average start_step : %d\n",C.Interval[Interval_Manager::tg_average].getStartStep());
  
  // 入力モードが有次元の場合に，無次元に変換
  if ( C.Unit.Param == DIMENSIONAL )
  {
    for (int i=0; i<Interval_Manager::tg_END; i++) {
      C.Interval[i].normalizeTime(C.Tscale);
    }
  }
  
  // Reference frame
  RF.setAccel( C.Interval[Interval_Manager::tg_accelra].getIntervalTime() );
  
  
  
  // インターバルの初期化
  double m_dt    = DT.get_DT();
  double m_tm    = CurrentTime;  // Restart()で設定
  unsigned m_stp = CurrentStep;

  
  if ( !C.Interval[Interval_Manager::tg_compute].initTrigger(m_stp, m_tm, m_dt, Interval_Manager::tg_compute) )
  {
    Hostonly_ printf("\t Error : Computation Period is assigned to zero.\n");
    Exit(0);
  }
  
  if ( !C.Interval[Interval_Manager::tg_console].initTrigger(m_stp, m_tm, m_dt, Interval_Manager::tg_console) )  // 基本履歴のコンソールへの出力
  {
    Hostonly_ printf("\t Error : Interval for Console output is assigned to zero.\n");
    Exit(0);
  }
  if ( !C.Interval[Interval_Manager::tg_history].initTrigger(m_stp, m_tm, m_dt, Interval_Manager::tg_history) )  // 履歴のファイルへの出力
  {
    Hostonly_ printf("\t Error : Interval for History output is assigned to zero.\n");
    Exit(0);
  }
  if ( !C.Interval[Interval_Manager::tg_instant].initTrigger(m_stp, m_tm, m_dt, Interval_Manager::tg_instant) )  // 瞬時値ファイル
  {
    Hostonly_ printf("\t Error : Interval for Instantaneous output is assigned to zero.\n");
    Exit(0);
  }
  
  if ( C.Mode.Average == ON )
  {
    if ( !C.Interval[Interval_Manager::tg_average].initTrigger(m_stp, m_tm, m_dt, Interval_Manager::tg_average) ) // 平均値ファイル
    {
      Hostonly_ printf("\t Error : Interval for Average output is assigned to zero.\n");
      Exit(0);
    }
  }
    
  if ( C.Sampling.log == ON )
  {
    if ( !C.Interval[Interval_Manager::tg_sampled].initTrigger(m_stp, m_tm, m_dt, Interval_Manager::tg_sampled) )  // サンプリング履歴
    {
      Hostonly_ printf("\t Error : Interval for Sampling output is assigned to zero.\n");
      Exit(0);
    }    
  }
  
  if (C.FIO.Format == plt3d_fmt)
  {
    if ( !C.Interval[Interval_Manager::tg_plot3d].initTrigger(m_stp, m_tm, m_dt, Interval_Manager::tg_plot3d) )  // 瞬時値ファイル
    {
      Hostonly_ printf("\t Error : Interval for plot3d output is assigned to zero.\n");
      Exit(0);
    }
  }
  
  Session_LastStep = C.Interval[Interval_Manager::tg_compute].getIntervalStep();
  
}



// #################################################################
/* @brief 距離の最小値を求める
 * @param [in,out] cut カット情報の配列
 * @param [in]     bid カット情報の配列
 * @param [in]     fp  file pointer
 */
void FFV::minDistance(float* cut, int* bid, FILE* fp)
{
  float global_min;
  float local_min = 1.0;
  float eps = 1.0/255.0; // 0.000392
  unsigned g=0;
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  TIMING_start(tm_cut_min);
  
#pragma omp parallel firstprivate(ix, jx, kx, eps, gd)
  {
    float th_min = 1.0;
    
#pragma omp for schedule(static) reduction(+:g)
    for (int k=1; k<=kx; k++) {
      for (int j=1; j<=jx; j++) {
        for (int i=1; i<=ix; i++) {
          
          size_t mp = _F_IDX_S4DEX(0, i, j, k, 6, ix, jx, kx, gd);
          size_t mb = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
          int bd = bid[mb];
          
          if ( TEST_BC(bd) ) // カットがあるか，IDによる判定
          {
            for (size_t n=0; n<6; n++) {
              float c = cut[mp+n];
              th_min = min(th_min, c);
              //th_min = min(th_min, cutPos->getPos(i+1,j+1,k+1,n)); slower than above inplementation
              
              if ( (c > 0.0) && (c <= eps) )
              {
                cut[mp+n] = eps;
                g++;
              }
            }
          }
          
          
// ##########
#if 0 // debug
          int b0 = (bd >> 0)  & MASK_5;
          int b1 = (bd >> 5)  & MASK_5;
          int b2 = (bd >> 10) & MASK_5;
          int b3 = (bd >> 15) & MASK_5;
          int b4 = (bd >> 20) & MASK_5;
          int b5 = (bd >> 25) & MASK_5;
          printf("%d %d %d %d %d %d ", b0, b1, b2, b3, b4, b5);
          /*
          printf("%3d %3d %3d : %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %d %d %d %d %d %d\n", i,j,k,
                 d_cut[mp+0],
                 d_cut[mp+1],
                 d_cut[mp+2],
                 d_cut[mp+3],
                 d_cut[mp+4],
                 d_cut[mp+5],
                 b0, b1, b2, b3, b4, b5);
           */
#endif
// ##########
        }
      }
    }
    
#pragma omp critical
    {
      local_min = min(local_min, th_min);
    }
  }
    
  global_min = local_min;
  
  // Allreduce時の桁あふれ対策のため、unsigned long で集約
  unsigned long gl = (unsigned long)g;
  
  if ( numProc > 1 )
  {
    float tmp = global_min;
    if ( paraMngr->Allreduce(&tmp, &global_min, 1, MPI_MIN) != CPM_SUCCESS ) Exit(0);
    
    unsigned long tmp_g = gl;
    if ( paraMngr->Allreduce(&tmp_g, &gl, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
  }
  
  TIMING_stop(tm_cut_min);

  if ( gl > 0 )
  {
    Hostonly_
    {
      fprintf(fp, "\tMinimum non-dimnensional distance is %e and replaced to %e : num = %ld\n\n", global_min, eps, gl);
      printf     ("\tMinimum non-dimnensional distance is %e and replaced to %e : num = %ld\n\n", global_min, eps, gl);
    }
  }

}


// #################################################################
/* @brief 履歴の出力準備
 */
void FFV::prep_HistoryOutput()
{
  // マスターノードでの履歴出力準備
  H = new History(&C);
  
  Hostonly_ {
    H->printHistoryTitle(stdout, IC, &C, true);
    
    // コンポーネント情報
    if ( C.Mode.Log_Base == ON ) 
    {
      // 基本情報　history.log, history_compo.log, history_domfx.log
      if ( !(fp_b=fopen(C.HistoryName.c_str(), "w")) ) 
      {
        stamped_printf("\tSorry, can't open '%s' file. Write failed.\n", C.HistoryName.c_str());
        Exit(0);
      }
      H->printHistoryTitle(fp_b, IC, &C, false);
      
      // コンポーネント履歴情報
      if ( !(fp_c=fopen(C.HistoryCompoName.c_str(), "w")) ) 
      {
        stamped_printf("\tSorry, can't open '%s' file. Write failed.\n", C.HistoryCompoName.c_str());
        Exit(0);
      }
      H->printHistoryCompoTitle(fp_c, cmp, &C);
      
      // 流量収支情報　
      if ( !(fp_d=fopen(C.HistoryDomfxName.c_str(), "w")) ) 
      {
        stamped_printf("\tSorry, can't open '%s' file. Write failed.\n", C.HistoryDomfxName.c_str());
        Exit(0);
      }
      H->printHistoryDomfxTitle(fp_d, &C);
      
      // 力の履歴情報　
      if ( !(fp_f=fopen(C.HistoryForceName.c_str(), "w")) ) 
      {
        stamped_printf("\tSorry, can't open '%s' file. Write failed.\n", C.HistoryForceName.c_str());
        Exit(0);
      }
      H->printHistoryForceTitle(fp_f);
    }
    
    // 反復履歴情報　history_itr.log
    if ( C.Mode.Log_Itr == ON ) 
    {
      if ( !(fp_i=fopen(C.HistoryItrName.c_str(), "w")) ) 
      {
				stamped_printf("\tSorry, can't open '%s' file.\n", C.HistoryItrName.c_str());
        Exit(0);
      }
    }
    
    // 壁面情報　history_wall.log
    if ( C.Mode.Log_Wall == ON ) 
    {
      if ( !(fp_w=fopen(C.HistoryWallName.c_str(), "w")) ) 
      {
				stamped_printf("\tSorry, can't open '%s' file.\n", C.HistoryWallName.c_str());
        Exit(0);
      }
      H->printHistoryWallTitle(fp_w);
    }
  }
  
  // サンプリング指定がある場合，モニタ結果出力ファイル群のオープン
  if ( C.Sampling.log == ON ) MO.openFile(C.HistoryMonitorName.c_str());
}


// #################################################################
/* @brief 読み込んだ領域情報のデバッグライト
 */
void FFV::printDomainInfo()
{
  cout << "\n####### read parameters ########" << endl;
  cout << " Gorg      = " << G_origin[0] << "," << G_origin[1] << "," << G_origin[2] << endl;
  cout << " Gvoxel    = " << G_size[0]   << "," << G_size[1]   << "," << G_size[2]   << endl;
  cout << " Gpitch    = " << pitch[0]    << "," << pitch[1]    << "," << pitch[2]    << endl;
  cout << " Gregion   = " << G_region[0] << "," << G_region[1] << "," << G_region[2] << endl;
  cout << " Gdiv      = " << G_division[0]  << "," << G_division[1]  << "," << G_division[2]  << endl;
}


// #################################################################
/* @brief セルBCのBbox情報をリサイズする
 * @param [in] order  CompoListのエントリ
 * @param [in] bx     BCindex
 */
void FFV::resizeBbox4Cell(const int order, const int* bx)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  int odr= order;
  
  // 初期値はローカルノードの大きさ
  int ist = ix;
  int jst = jx;
  int kst = kx;
  int ied = 0;
  int jed = 0;
  int ked = 0;
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd, odr) \
shared(ist, ied, jst, jed, kst, ked) schedule(static)
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        int s = bx[m];
        
        if ( DECODE_CMP(s) == odr )
        {
          if( i < ist ) ist = i;
          if( i > ied ) ied = i;
          if( j < jst ) jst = j;
          if( j > jed ) jed = j;
          if( k < kst ) kst = k;
          if( k > ked ) ked = k;
        }
      }
    }
  }
  
  // replace
  int nst[3] = {ist, jst, kst};
  int ned[3] = {ied, jed, ked};
  cmp[order].setBbox(nst, ned);
}


// #################################################################
/* @brief セルフェイスBCのBbox情報をリサイズする
 * @param [in] order  CompoListのエントリ
 * @param [in] bx     BCindex
 */
void FFV::resizeBbox4Face(const int order, const int* bx)
{
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  int odr= order;
  
  // 初期値はローカルノードの大きさ
  int ist = ix;
  int jst = jx;
  int kst = kx;
  int ied = 0;
  int jed = 0;
  int ked = 0;
  
  
#pragma omp parallel for firstprivate(ix, jx, kx, gd, odr) \
shared(ist, ied, jst, jed, kst, ked) schedule(static)
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        
        size_t m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        int s = bx[m];
        
        int flag = 0;
        
        if ( GET_FACE_BC(s, BC_FACE_W) == odr ) flag++;
        /*
        {
          if( i < nst[0] ) nst[0] = i;
          if( i > ned[0] ) ned[0] = i;
          if( j < nst[1] ) nst[1] = j;
          if( j > ned[1] ) ned[1] = j;
          if( k < nst[2] ) nst[2] = k;
          if( k > ned[2] ) ned[2] = k;
        }*/
        
        if ( GET_FACE_BC(s, BC_FACE_E) == odr ) flag++;
        /*
        {
          if( i < nst[0] ) nst[0] = i;
          if( i > ned[0] ) ned[0] = i;
          if( j < nst[1] ) nst[1] = j;
          if( j > ned[1] ) ned[1] = j;
          if( k < nst[2] ) nst[2] = k;
          if( k > ned[2] ) ned[2] = k;
        }*/
        
        if ( GET_FACE_BC(s, BC_FACE_S) == odr ) flag++;
        /*
        {
          if( i < nst[0] ) nst[0] = i;
          if( i > ned[0] ) ned[0] = i;
          if( j < nst[1] ) nst[1] = j;
          if( j > ned[1] ) ned[1] = j;
          if( k < nst[2] ) nst[2] = k;
          if( k > ned[2] ) ned[2] = k;
        }*/
        
        if ( GET_FACE_BC(s, BC_FACE_N) == odr ) flag++;
        /*
        {
          if( i < nst[0] ) nst[0] = i;
          if( i > ned[0] ) ned[0] = i;
          if( j < nst[1] ) nst[1] = j;
          if( j > ned[1] ) ned[1] = j;
          if( k < nst[2] ) nst[2] = k;
          if( k > ned[2] ) ned[2] = k;
        }*/
        
        if ( GET_FACE_BC(s, BC_FACE_B) == odr ) flag++;
        /*
        {
          if( i < nst[0] ) nst[0] = i;
          if( i > ned[0] ) ned[0] = i;
          if( j < nst[1] ) nst[1] = j;
          if( j > ned[1] ) ned[1] = j;
          if( k < nst[2] ) nst[2] = k;
          if( k > ned[2] ) ned[2] = k;
        }*/
        
        if ( GET_FACE_BC(s, BC_FACE_T) == odr ) flag++;
        /*
        {
          if( i < nst[0] ) nst[0] = i;
          if( i > ned[0] ) ned[0] = i;
          if( j < nst[1] ) nst[1] = j;
          if( j > ned[1] ) ned[1] = j;
          if( k < nst[2] ) nst[2] = k;
          if( k > ned[2] ) ned[2] = k;
        }*/
        
        if ( flag != 0 )
        {
          if( i < ist ) ist = i;
          if( i > ied ) ied = i;
          if( j < jst ) jst = j;
          if( j > jed ) jed = j;
          if( k < kst ) kst = k;
          if( k > ked ) ked = k;
        }
      }
    }
  }
  
  // replace
  int nst[3] = {ist, jst, kst};
  int ned[3] = {ied, jed, ked};
  cmp[order].setBbox(nst, ned);
  printf("%d %d %d %d %d %d\n", ist, ied, jst, jed, kst, ked);
}


// #################################################################
/* @brief コンポーネントリストに登録されたBbox情報をリサイズする
 */
void FFV::resizeCompoBbox()
{

  for (int n=1; n<=C.NoCompo; n++)
  {
    int typ = cmp[n].getType();
    
    // デフォルトで流れ計算パートのBC
    switch ( typ )
    {
      case SPEC_VEL:
      case SPEC_VEL_WH:
      case OUTFLOW:
        resizeBbox4Face(n, d_bcv); // 速度のBCindex　セルフェイス位置でのflux型BC
        break;
        
      case PERIODIC:
      case IBM_DF:
      case HEX:
      case FAN:
      case DARCY:
      case CELL_MONITOR:
        resizeBbox4Cell(n, d_bcd); // セルセンタ位置でのBC
        break;
    }
    
    // 熱計算のパート
    if ( C.isHeatProblem() )
    {
      switch ( typ )
      {
        case ADIABATIC:
        case HEATFLUX:
        case SPEC_VEL_WH:
        case OUTFLOW:
        case TRANSFER:
        case ISOTHERMAL:
          resizeBbox4Face(n, d_bh1);
          break;
          
        case RADIANT:
          break;
          
        case HEAT_SRC:
        case CNST_TEMP:
          resizeBbox4Cell(n, d_bh2);
          break;
      }
    }
    
  }
}


// #################################################################
/* @brief ボクセルをスキャンし情報を表示する
 * @param [in] fp file pointer
 */
void FFV::scanVoxel(FILE* fp)
{
  // 媒質IDリスト
  int cell_id[C.NoCompo+1];
  
  for (int i=1; i<=C.NoCompo; i++) cell_id[i]=0;
  
  
  // 外部境界面の媒質IDをセット
  for (int i=0; i<NOFACE; i++)
  {
    int m = BC.export_OBC(i)->getGuideMedium();
    
    if ( m<1 || m>C.NoMedium )
    {
      Hostonly_
      {
        stamped_printf (   "\tError : An ID of guide cell is out of range.\n");
        stamped_fprintf(fp,"\tError : An ID of guide cell is out of range.\n");
      }
      Exit(0);
    }
    
    cell_id[m] = 1;
    
  }
  
  
  // midに設定されたIDをスキャンし，IDの個数を返し，作業用のcolorList配列にIDを保持，midに含まれるIDの数をチェック
  V.scanCell(d_mid, cell_id, C.Hide.Change_ID, fp);
  
}


// #################################################################
/* @brief 外部境界条件を読み込み，Controlクラスに保持する
 */
void FFV::setBCinfo()
{

  // パラメータファイルをパースして，外部境界条件を保持する　>> VoxScan()以前に処理
  B.loadBC_Outer( BC.export_OBC(), mat, cmp );

  
  // パラメータファイルの情報を元にCompoListの情報を設定する
  B.loadBC_Local(&C, mat, cmp, PolyPP);
  
  
#if 0
  Hostonly_ B.checkList(mat, cmp);
#endif
  
  
  // 各コンポーネントが存在するかどうかを保持しておく
  setEnsComponent();
  
  // KOSと境界条件種類の整合性をチェック
  B.chkBCconsistency(C.KindOfSolver, cmp);
  
}



// #################################################################
/* @brief HEX,FANコンポーネントなどの体積率とbboxなどをセット
 * @note インデクスの登録と配列確保はVoxEncode()で、コンポーネント領域のリサイズ後に行う
 */
void FFV::setComponentVF()
{
  const int subsampling = 20; // 体積率のサブサンプリングの基数
  int f_st[3], f_ed[3];
  double flop;
  
  CompoFraction CF(size, guide, (float*)pitch, (float*)origin, subsampling);
  
  for (int n=1; n<=C.NoCompo; n++)
  {
    
    if ( cmp[n].isFORCING() )
    {
      // 形状パラメータのセット
      switch ( cmp[n].getType() ) 
      {
        case HEX:
          CF.setShapeParam((float*)cmp[n].nv, (float*)cmp[n].oc, (float*)cmp[n].dr, (float)cmp[n].depth, (float)cmp[n].shp_p1, (float)cmp[n].shp_p2);
          break;
          
        case FAN:
          CF.setShapeParam((float*)cmp[n].nv, (float*)cmp[n].oc, (float)cmp[n].depth, (float)cmp[n].shp_p1, (float)cmp[n].shp_p2);
          break;
          
        case DARCY:
          Exit(0);
          break;
          
        default:
          Exit(0);
          break;
      }
      
      // 回転角度の計算
      CF.get_angle(); 
      
      // bboxと投影面積の計算
      cmp[n].area = CF.get_BboxArea();
      
      // インデクスの計算 > あとで，VoxEncode()でresize
      CF.bbox_index(f_st, f_ed);
      
      // インデクスのサイズ登録と存在フラグ
      cmp[n].setBbox(f_st, f_ed);
      cmp[n].setEnsLocal(ON);
      
      // 体積率
      TIMING_start(tm_cmp_vertex8);
      flop = 0.0;
      CF.vertex8(f_st, f_ed, d_cvf, flop);
      TIMING_stop(tm_cmp_vertex8, flop);
      
      TIMING_start(tm_cmp_subdivision);
      flop = 0.0;
      CF.subdivision(f_st, f_ed, d_cvf, flop);
      TIMING_stop(tm_cmp_subdivision, flop);
    }
  }
  
// ########## 確認のための出力
#if 0
  REAL_TYPE org[3], pit[3];
  
  //  ガイドセルがある場合(GuideOut != 0)にオリジナルポイントを調整
  for (int i=0; i<3; i++) 
  {
    org[i] = C.org[i] - C.dx[i]*(REAL_TYPE)C.GuideOut;
    pit[i] = C.dx[i];
  }
  
  // 出力ファイルの指定が有次元の場合
  if ( C.Unit.File == DIMENSIONAL ) 
  {
    for (int i=0; i<3; i++) 
    {
      org[i] *= C.RefLength;
      pit[i] *= C.RefLength;
    }
  }
  F.writeRawSPH(cvf, size, guide, org, pit, sizeof(float));
#endif
// ##########
  
}



// #################################################################
/* @brief コンポーネントが存在するかを保持しておく
 */
void FFV::setEnsComponent()
{
  int c;
  
  // Vspec
  c = 0;
  for (int n=1; n<=C.NoCompo; n++)
  {
    if ( cmp[n].isVBC_IO() ) c++;
  }
  if ( c>0 ) C.EnsCompo.forcing = ON;
  
  // Forcing > HEX, FAN, DARCY
  c = 0;
  for (int n=1; n<=C.NoCompo; n++)
  {
    if ( cmp[n].isFORCING() ) c++;
  }
  if ( c>0 ) C.EnsCompo.forcing = ON;
  
  // Heat source > HEAT_SRC, CNST_TEMP
  c = 0;
  for (int n=1; n<=C.NoCompo; n++)
  {
    if ( cmp[n].isHsrc() ) c++;
  }
  if ( c>0 ) C.EnsCompo.hsrc = ON;
  
  // 周期境界 > PERIODIC
  c = 0;
  for (int n=1; n<=C.NoCompo; n++) 
  {
    if ( cmp[n].getType() == PERIODIC ) c++;
  }
  if ( c>0 ) C.EnsCompo.periodic = ON;
  
  // 流出境界 > OUTFLOW
  c = 0;
  for (int n=1; n<=C.NoCompo; n++) 
  {
    if ( cmp[n].getType() == OUTFLOW ) c++;
  }
  if ( c>0 ) C.EnsCompo.outflow = ON;
  
  // 体積率コンポーネント
  c = 0;
  for (int n=1; n<=C.NoCompo; n++) 
  {
    if ( cmp[n].isVFraction() ) c++;
  }
  if ( c>0 ) C.EnsCompo.fraction = ON;
  
  // トラクションフリー
  c = 0;
  for (int n=0; n<NOFACE; n++)
  {
    if ( BC.export_OBC(n)->get_Class()== OBC_TRC_FREE ) c++;
  }
  if ( c>0 ) C.EnsCompo.tfree = ON;
  
  // モニタ
  c = 0;
  for (int n=1; n<=C.NoCompo; n++) 
  {
    if ( cmp[n].isMONITOR() ) c++;
  }
  
  // チェック　MONITOR_LISTでCELL_MONITORが指定されている場合
  if ( (C.isMonitor() == ON) && (c < 1) ) 
  {
    Hostonly_ stamped_printf("\tError : CellMonitor in MonitorList is specified, however any Monitor can not be found.\n");
    Exit(0);
  }
  
  if ( (C.isMonitor() == OFF) && (c > 0) ) 
  {
    Hostonly_ stamped_printf("\tError : CellMonitor in MonitorList is NOT specified, however Monitor section is found in LocalBoundary.\n");
    Exit(0);
  }
  
  // @see C.EnsCompo.monitor==ONは，Control::getMonitorList()
  
}



// #################################################################
/* @brief コンポーネントのローカルなBbox情報からグローバルなBbox情報を求める
 */
void FFV::setGlobalCmpIdx()
{
  int st_i, st_j, st_k, ed_i, ed_j, ed_k;
  int node_st_i, node_st_j, node_st_k;
  int st[3], ed[3];
  
  // グローバルインデクスの配列インスタンス
  compo_global_bbox = new int[6*(C.NoCompo+1)];
  int* cgb = compo_global_bbox;
  
  // ローカルインデクスからグローバルインデクスに変換
  for (int m=1; m<=C.NoCompo; m++) 
  {
    
    if ( !cmp[m].isEnsLocal() ) // コンポーネントが存在しないノードはゼロを代入
    {
      cgb[6*m+0] = 0;
      cgb[6*m+1] = 0;
      cgb[6*m+2] = 0;
      cgb[6*m+3] = 0;
      cgb[6*m+4] = 0;
      cgb[6*m+5] = 0;
    }
    else // コンポーネントが存在する場合
    {
      cmp[m].getBbox(st, ed);

      st_i = st[0];
      st_j = st[1];
      st_k = st[2];
      ed_i = ed[0];
      ed_j = ed[1];
      ed_k = ed[2];
      
      if ( numProc > 1 ) 
      {
        node_st_i = head[0];
        node_st_j = head[1];
        node_st_k = head[2];
        
        cgb[6*m+0] = node_st_i + st_i;
        cgb[6*m+1] = node_st_j + st_j;
        cgb[6*m+2] = node_st_k + st_k;
        cgb[6*m+3] = node_st_i + ed_i;
        cgb[6*m+4] = node_st_j + ed_j;
        cgb[6*m+5] = node_st_k + ed_k;
      } 
      else 
      {
        cgb[6*m+0] = st_i;
        cgb[6*m+1] = st_j;
        cgb[6*m+2] = st_k;
        cgb[6*m+3] = ed_i;
        cgb[6*m+4] = ed_j;
        cgb[6*m+5] = ed_k;
      }

    }
  }
  
  // 領域全体のbboxを求める
  int* m_gArray = NULL;
  int* m_eArray = NULL;
  int array_size = 6*(C.NoCompo+1);
  int st_x, st_y, st_z, ed_x, ed_y, ed_z, es;
  
  if ( !(m_gArray = new int[numProc*6]) ) Exit(0);
  if ( !(m_eArray = new int[numProc]  ) ) Exit(0);
  
  for (int n=1; n<=C.NoCompo; n++)
  {
    if ( numProc > 1 ) 
    {
      es = ( cmp[n].isEnsLocal() ) ? 1 : 0;
      if ( paraMngr->Gather(&es, 1, m_eArray, 1, 0) != CPM_SUCCESS ) Exit(0);
      if ( paraMngr->Gather(&cgb[6*n], 6, m_gArray, 6, 0) != CPM_SUCCESS ) Exit(0);
      
      
      if (myRank == 0) { // マスターノードのみ
        
        // 初期値
        cgb[6*n+0] = 100000000;
        cgb[6*n+1] = 100000000;
        cgb[6*n+2] = 100000000;
        cgb[6*n+3] = 0;
        cgb[6*n+4] = 0;
        cgb[6*n+5] = 0;
        
        for (int m=0; m<numProc; m++) 
        {
          if ( m_eArray[m]==1 ) // コンポーネントの存在ランクのみを対象とする
          {
            st_x = m_gArray[6*m+0]; // 各ランクのコンポーネント存在範囲のインデクス
            st_y = m_gArray[6*m+1];
            st_z = m_gArray[6*m+2];
            ed_x = m_gArray[6*m+3];
            ed_y = m_gArray[6*m+4];
            ed_z = m_gArray[6*m+5];
            
            if( st_x < cgb[6*n+0] ) { cgb[6*n+0] = st_x; } // 最大値と最小値を求める
            if( st_y < cgb[6*n+1] ) { cgb[6*n+1] = st_y; }
            if( st_z < cgb[6*n+2] ) { cgb[6*n+2] = st_z; }
            if( ed_x > cgb[6*n+3] ) { cgb[6*n+3] = ed_x; }
            if( ed_y > cgb[6*n+4] ) { cgb[6*n+4] = ed_y; }
            if( ed_z > cgb[6*n+5] ) { cgb[6*n+5] = ed_z; }
          }
        }
      }
    }
    
  }
  
  // destroy
  if (m_gArray) 
  {
    delete[] m_gArray;
    m_gArray = NULL;
  }
  
  if (m_eArray) 
  {
    delete[] m_eArray;
    m_eArray = NULL;
  }
}



// #################################################################
// @brief 初期条件の設定
void FFV::setInitialCondition()
{
  double flop_task;
  Gemini_R* m_buf = new Gemini_R [C.NoCompo];
  
  double tm = CurrentTime * C.Tscale;
  
  if ( C.Start == initial_start )
  {
		REAL_TYPE U0[3];
    
		// 速度の初期条件の設定
    if (C.Unit.Param == DIMENSIONAL)
    {
      U0[0] = C.iv.VecU/C.RefVelocity;
      U0[1] = C.iv.VecV/C.RefVelocity;
      U0[2] = C.iv.VecW/C.RefVelocity;
    }
    else
    {
      U0[0] = C.iv.VecU;
      U0[1] = C.iv.VecV;
      U0[2] = C.iv.VecW;
    }
		fb_set_vector_(d_v, size, &guide, U0, d_bcd);
    fb_set_fvector_(d_vf, size, &guide, U0, d_bcd);
    
    // セルフェイスの設定　発散値は関係なし
    BC.mod_div(d_dv, d_bcv, CurrentTime, v00, m_buf, d_vf, d_v, &C, flop_task);
    
		// 外部境界面の流出流量と移流速度
    DomainMonitor( BC.export_OBC(), &C);
    
		// 外部境界面の移流速度を計算し，外部境界条件を設定
		BC.OuterVBC(d_v, d_vf, d_bcv, tm, &C, v00, flop_task);
    BC.InnerVBC(d_v, d_bcv, tm, v00, flop_task);
    BC.InnerVBC_Periodic(d_v, d_bcd);
    
		// 圧力
    REAL_TYPE ip;
    if (C.Unit.Param == DIMENSIONAL)
    {
      ip = FBUtility::convD2ND_P(C.iv.Pressure, C.BasePrs, C.RefDensity, C.RefVelocity, C.Unit.Prs);
    }
    else
    {
      ip = C.iv.Pressure;
    }
    
    U.xset(d_p, size, guide, ip, kind_scalar);
		BC.OuterPBC(d_p);
    
		// 温度
		if ( C.isHeatProblem() )
    {
      REAL_TYPE it;
      if (C.Unit.Param == DIMENSIONAL)
      {
        it = FBUtility::convK2ND(C.iv.Temperature, C.BaseTemp, C.DiffTemp);
      }
      else
      {
        it = C.iv.Temperature;
      }
      
      U.xset(d_t, size, guide, it, kind_scalar);
      
      // コンポーネントの初期値
      for (int m=1; m<=C.NoCompo; m++)
      {
        BC.setInitialTempCompo(m, d_bcd, d_t);
      }
      
			BC.OuterTBC(d_t);
		}
    
  }
  else // リスタート時
  {
    // 内部境界条件
    BC.InnerVBC(d_v, d_bcv, tm, v00, flop_task);
    BC.InnerVBC_Periodic(d_v, d_bcd);
    BC.InnerPBC_Periodic(d_p, d_bcd);
    
    // 外部境界条件
    BC.OuterVBC(d_v, d_vf, d_bcv, tm, &C, v00, flop_task);
    
    // 流出境界の流出速度の算出
    // dummy
    Gemini_R* m_buf = new Gemini_R [C.NoCompo];
    BC.mod_div(d_ws, d_bcv, tm, v00, m_buf, d_vf, d_v, &C, flop_task);
    if ( m_buf ) { delete [] m_buf; m_buf=NULL; }
    
    DomainMonitor(BC.export_OBC(), &C);
    
    //if ( C.isHeatProblem() ) BC.InnerTBC_Periodic()
    
  }

  
  // 初期解およびリスタート解の同期
  if ( numProc > 1 )
  {
    if ( paraMngr->BndCommV3D(d_v,  size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
    if ( paraMngr->BndCommV3D(d_vf, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
    if ( paraMngr->BndCommS3D(d_p,  size[0], size[1], size[2], guide, 1    ) != CPM_SUCCESS ) Exit(0);
    
    if ( C.isHeatProblem() ) 
    {
      if ( paraMngr->BndCommS3D(d_p, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
    }
  }

  // VOF
  if ( C.BasicEqs == INCMP_2PHASE )
  {
    setVOF();
    if ( numProc > 1 )
    {
      if ( paraMngr->BndCommS3D(d_vof, size[0], size[1], size[2], guide, guide) != CPM_SUCCESS ) Exit(0);
    }
  }
  
  // 後始末
  if ( m_buf ) delete [] m_buf;
}



// #################################################################
/* @brief midの情報から各BCコンポーネントのローカルなインデクスを取得する
 * @note  計算内部領域の境界と外部境界とでは，ガイドセル部分にあるコンポーネントIDの取り扱いが異なる
 *        外部境界に接する面では，幅はそのまま，始点はガイドセル部分を含む
 *        内部境界に接する面では，始点と幅はローカルノード内の計算内部領域に含まれるように調整
 */
void FFV::setLocalCmpIdx_Binary()
{
  int st_i[3], len[3];
  int m_st[3], m_ed[3];
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  for (int m=1; m<=C.NoCompo; m++)
  {
    switch ( cmp[m].getType() ) 
    {
      case HEX:
      case FAN:
        break; // 体積率でエンコード済み
        
      default:
        
        for (int d=0; d<3; d++) 
        {
          st_i[d] = 0;
          len[d] = 0;
        }
        
        // コンポーネント範囲
        //GetBndIndexExtGc()は自ノード内でのidのバウンディングボックスを取得．インデクスはローカルインデクスで，ガイドセルを含む配列の基点をゼロとするCのインデクス
        if ( !paraMngr->GetBndIndexExtGc(m, d_mid, ix, jx, kx, gd, st_i[0], st_i[1], st_i[2], len[0], len[1], len[2]) )
        {
          Hostonly_ stamped_printf("\tError : can not get component local index for ID[%d]\n", m);
          Exit(0);
        }
        
        // ノード内にコンポーネントがあるかどうかをチェック
        if ( (len[0]==0) || (len[1]==0) || (len[2]==0) ) { // コンポーネントなし
          cmp[m].setEnsLocal(OFF);
          // BV情報はCompoListのコンストラクタでゼロに初期化されているので，すべてゼロ
        }
        else 
        {
          
          for (int d=0; d<3; d++) 
          {
            int tmp_st=0;
            int tmp_ed=0;
            
            EnlargeIndex(tmp_st, tmp_ed, st_i[d], len[d], size[d], d, m);
            
            m_st[d] = tmp_st;
            m_ed[d] = tmp_ed;
          }
          cmp[m].setBbox(m_st, m_ed);
          cmp[m].setEnsLocal(ON); // コンポーネントあり
        }
        
        break;
    }
    
  }
}


// #################################################################
/* @brief ParseMatクラスをセットアップし，媒質情報を入力ファイルから読み込み，媒質リストを作成する
 * @param [in] fp  ファイルポインタ
 */
void FFV::setMediumList(FILE* fp)
{
  Hostonly_
  {
    printf(    "\n---------------------------------------------------------------------------\n\n");
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
    printf(    "\n\t>> Medium List\n\n");
    fprintf(fp,"\n\t>> Medium List\n\n");
  }
  
  
  // 媒質情報をロードし、 MediumTableタグ内のコンポーネント数，BC数を得る
  M.getMediumTable(C.NoMedium, mat);
  

  if ( !M.check(mat) )
  {
    Hostonly_ stamped_printf("Error : Duplicate label in Material Table\n");
    Exit(0);
  }

  // 媒質テーブルの表示
  Hostonly_
  {
    M.printMatList(stdout, mat);
    M.printMatList(fp, mat);
  }

}




// #################################################################
/* @brief 各種例題のモデルをセットアップ
 * @param [in] PrepMemory  前処理に必要なメモリ
 * @param [in] TotalMemory ソルバー実行に必要なメモリ
 * @param [in] fp          ファイルポインタ
 */
void FFV::setModel(double& PrepMemory, double& TotalMemory, FILE* fp)
{
  // d_midをゼロで初期化
  size_t mt = (size[0]+2*guide) * (size[1]+2*guide) *(size[2]+2*guide) * sizeof(int);
  memset(d_mid, 0, mt);
  
  
  switch (C.Mode.Example)
  {
    case id_Polygon: // ユーザ例題
      
      C.getGeometry(mat);
      
      
      // PolylibとCutlibのセットアップ
      setupPolygon2CutInfo(PrepMemory, TotalMemory, fp);
      break;
      
    case id_Sphere:
      if ( !C.isCDS() ) // binary
      {
        Ex->setup(d_mid, &C, G_origin, C.NoCompo, mat);
      }
      else // cut
      {
        // 初期値1.0をセット
        setupCutInfo4IP(PrepMemory, TotalMemory, fp);
        ((IP_Sphere*)Ex)->setup_cut(d_mid, &C, G_origin, C.NoCompo, mat, d_cut);
      }
      break;
      
    case id_SHC1D:
      setupCutInfo4IP(PrepMemory, TotalMemory, fp);
      Ex->setup(d_mid, &C, G_origin, C.NoCompo, mat);
      ((IP_SHC1D*)Ex)->setup_bc(d_bid);
      break;
      
    default: // ほかのIntrinsic problems
      if ( C.isCDS() ) // カットの場合
      {
        setupCutInfo4IP(PrepMemory, TotalMemory, fp);
      }
      
      Ex->setup(d_mid, &C, G_origin, C.NoCompo, mat);
      break;
  }
  
  // midのガイドセル同期
  if ( paraMngr->BndCommS3D(d_mid, size[0], size[1], size[2], guide, 1) != CPM_SUCCESS ) Exit(0);
  
}


// #################################################################
/**
 * @brief MonitorListの処理
 */
void FFV::setMonitorList()
{
  MO.setControlVars(d_bcd,
                    C.RefVelocity,
                    C.BaseTemp,
                    C.DiffTemp,
                    C.RefDensity,
                    C.RefLength,
                    C.BasePrs,
                    C.Unit.Temp,
                    C.Mode.Precision,
                    C.Unit.Prs,
                    C.num_process);
  
  
  //パラメータを取得し，セットの配列を確保する
  MO.getMonitor(&C);
  
  
  // 内部境界条件として指定されたモニタ設定を登録
  MO.setInnerBoundary(cmp, C.NoCompo);
}


// #################################################################
/* @brief 並列化と分割の方法を保持
 * @return 並列モード
 */
string FFV::setParallelism()
{
  string para_mode;
  
  C.num_thread  = omp_get_max_threads();
  
  // Serial or Parallel environment
  if( paraMngr->IsParallel() )
  {
    C.num_process = paraMngr->GetNumRank();
    
    if ( C.num_thread > 1 ) 
    {
      C.Parallelism = Control::Hybrid;
      para_mode = "Hybrid";
    }
    else 
    {
      C.Parallelism = Control::FlatMPI;
      para_mode = "FlatMPI";
    }
  }
  else 
  {
    C.num_process = 1;
    
    if ( C.num_thread > 1 ) 
    {
      C.Parallelism = Control::OpenMP;
      para_mode = "OpenMP";
    }
    else 
    {
      C.Parallelism = Control::Serial;
      para_mode = "Serial";
    }
  }
  return para_mode;
}


// #################################################################
/* @brief 時間積分幅や物理パラメータの設定
 */
void FFV::setParameters()
{
  // 無次元数などの計算パラメータを設定する．MediumListを決定した後，かつ，SetBC3Dクラスの初期化前に実施すること
  // 代表物性値をRefMatの示す媒質から取得
  // Δt=constとして，無次元の時間積分幅 deltaTを計算する．ただし，一定幅の場合に限られる．不定幅の場合には別途考慮の必要
  DT.set_Vars(C.KindOfSolver, C.Unit.Param, (double)deltaX, (double)C.Reynolds, (double)C.Peclet);
  
  
  // 無次元速度1.0を与えてdeltaTをセットし，エラーチェック
  switch ( DT.set_DT(1.0) ) 
  {
    case 0: // 成功
      break;
      
    case 1:
      Hostonly_ stamped_printf("\tdt selection error(1) : 'Kind of Solver' is solid conduction. Consider to specify other dt scheme or confirm 'Kind of Solver'.\n");
      Exit(0);
      break;
      
    case 2:
      Hostonly_ stamped_printf("\tdt selection error(2) : 'Kind of Solver' is solid conduction. Consider to specify other dt scheme or confirm 'Kind of Solver'.\n");
      Exit(0);
      break;
      
    case 3:
      Hostonly_ stamped_printf("\tdt selection error(3) : 'Kind of Solver' includes flow effect. Consider to specify other dt scheme or confirm 'Kind of Solver'.\n");
      Exit(0);
      break;
      
    case 4:
      Hostonly_ stamped_printf("\tdt selection error(4) : 'Kind of Solver' is solid conduction. Consider to specify other dt scheme or confirm 'Kind of Solver'.\n");
      Exit(0);
      break;
      
    case 5:
      Hostonly_ stamped_printf("\tdt selection error(5) : 'Kind of Solver' is solid conduction. Consider to specify other dt scheme or confirm 'Kind of Solver'.\n");
      Exit(0);
      break;
      
    case 6:
      Hostonly_ stamped_printf("\tdt selection error(6) : CFL max V and Diffusion ; Not implemented yet\n");
      Exit(0);
      break;
      
    case 7:
      Hostonly_ stamped_printf("\tdt selection error(7) : CFL max V for cp ;  Not implemented yet\n");
      Exit(0);
      break;
      
    default:
      Exit(0);
      break;
  }
  
  
  // 無次元時間積分幅
  deltaT = DT.get_DT();

  
  C.setParameters(mat, cmp, &RF, BC.export_OBC());
  
  
  // パラメータの無次元化（正規化）に必要な参照物理量の設定
  B.setRefMediumProperty(mat, C.RefMat);
}



// #################################################################
/* @brief BCで指定するCELL_MINITORのIDをmid[]にセットする
 */
void FFV::setupCellMonitor()
{
  
  // Cell_Monitorの指定がある場合，モニタ位置をセット
  if ( (C.Sampling.log == ON) && (C.isMonitor() == ON) )
  {
    
    // ShapeMonitorのインスタンス
    ShapeMonitor SM(size, guide, (float*)pitch, (float*)origin);
    
    for (int n=1; n<=C.NoCompo; n++)
    {
      
      if ( cmp[n].isMONITOR() )
      {
        switch ( cmp[n].get_Shape() )
        {
          case SHAPE_BOX:
          case SHAPE_CYLINDER:
            V.setMonitorShape(d_mid, n, &SM, cmp, C.RefLength);
            cmp[n].setEnsLocal(ON);
            break;
            
          case SHAPE_VOXEL:
            V.setMonitorCellID(d_mid, d_bid, d_cut, n, C.Fill_Fluid, cmp[n].getSamplingWidth());
            cmp[n].setEnsLocal(ON);
            break;
            
          default:
            Exit(0);
            break;
        }
      }
    }
  }
  
}



// #################################################################
/* @brief IP用にカット領域をアロケートする
 * @param [in,out] m_prep  前処理用のメモリサイズ
 * @param [in,out] m_total 本計算用のメモリリサイズ
 * @param [in]     fp      ファイルポインタ
 */
void FFV::setupCutInfo4IP(double& m_prep, double& m_total, FILE* fp)
{
  Hostonly_ 
  {
    fprintf(fp, "\n---------------------------------------------------------------------------\n\n");
    fprintf(fp, "\t>> Cut Info\n\n");
    printf("\n---------------------------------------------------------------------------\n\n");
    printf("\t>> Cut Info\n\n");
  }
  
  size_t n_cell[3];
  
  for (int i=0; i<3; i++) 
  {
    n_cell[i] = (size_t)(size[i] + 2*guide); // 分割数+ガイドセル
  }
  size_t size_n_cell = n_cell[0] * n_cell[1] * n_cell[2];
  
  
  // カット情報保持領域のアロケート
  TIMING_start(tm_init_alloc);
  allocArray_Cut(m_total);
  TIMING_stop(tm_init_alloc);
  
  // 使用メモリ量　
  double cut_mem, G_cut_mem;
  G_cut_mem = cut_mem = (double)size_n_cell * (double)(6*sizeof(float) + sizeof(int));
  m_prep += cut_mem;
  m_total+= cut_mem;
  
  
  displayMemoryInfo(fp, G_cut_mem, cut_mem, "Cut");
  

  
  // 初期値のセット
  for (size_t i=0; i<size_n_cell*6; i++) {
    d_cut[i] = 1.0f;
  }
  
  for (size_t i=0; i<size_n_cell; i++) {
    d_bid[i] = 0;
  }
  
}




// #################################################################
/* @brief パラメータのロードと計算領域を初期化し，並列モードを返す
 * @param [in] tpf ffvのパラメータを保持するTextParserインスタンス
 * @retval 並列モードの文字列
 */
string FFV::setupDomain(TPControl* tpf)
{
  
  // ランク情報をセット >> 各クラスでランク情報メンバ変数を利用する前にセットすること
  C.setRankInfo    (paraMngr, procGrp);
  B.setRankInfo    (paraMngr, procGrp);
  V.setRankInfo    (paraMngr, procGrp);
  F.setRankInfo    (paraMngr, procGrp);
  BC.setRankInfo   (paraMngr, procGrp);
  Ex->setRankInfo  (paraMngr, procGrp);
  MO.setRankInfo   (paraMngr, procGrp);
  /* 20130611
  FP3DR.setRankInfo(paraMngr, procGrp);
  FP3DW.setRankInfo(paraMngr, procGrp);
   */
  
  
  // 並列モードの取得
  string str = setParallelism();
  
  
  // 最初のパラメータの取得
  // 20130611 C.get_Steer_1(&DT, &FP3DR, &FP3DW);
  C.get1stParameter(&DT);
  
  // 代表パラメータをコピー
  Ex->setRefParameter(&C);
  
  
  // 例題クラス固有のパラメータを取得
  if ( !Ex->getTP(&C, tpf) ) Exit(0);
  
  
  // 領域設定 計算領域全体のサイズ，並列計算時のローカルのサイズ，コンポーネントのサイズなどを設定
  DomainInitialize(tpf);
  
  
  // 各例題の領域パラメータを設定する
  Ex->setDomainParameter(&C, size, origin, region, pitch);
  
  
  // 各クラスで領域情報を保持
  C.setNeighborInfo    (C.guide);
  B.setNeighborInfo    (C.guide);
  V.setNeighborInfo    (C.guide);
  F.setNeighborInfo    (C.guide);
  BC.setNeighborInfo   (C.guide);
  Ex->setNeighborInfo  (C.guide);
  MO.setNeighborInfo   (C.guide);
  /* 20130611
  FP3DR.setNeighborInfo(C.guide);
  FP3DW.setNeighborInfo(C.guide);
   */
  
  
  // 従属的なパラメータの取得
  C.get2ndParameter(IC, &RF);

  
  return str;
}


// #################################################################
/* @brief 幾何形状情報を準備し，交点計算を行う
 * @param [in,out] m_prep   前処理用のメモリサイズ
 * @param [in,out] m_total  本計算用のメモリリサイズ
 * @param [in]     fp       ファイルポインタ
 */
void FFV::setupPolygon2CutInfo(double& m_prep, double& m_total, FILE* fp)
{
  unsigned poly_gc[3] = {guide, guide, guide};
  
  // float で定義すること
  float poly_org[3];
  float poly_dx[3];

  
  // 有次元に変換 Polylib: 並列計算領域情報　ポリゴンは実スケールで，ガイドセル領域部分も含めて指定する
  poly_dx[0]  = pitch[0] * C.RefLength;
  poly_dx[1]  = pitch[1] * C.RefLength;
  poly_dx[2]  = pitch[2] * C.RefLength;
  poly_org[0] = origin[0]* C.RefLength;
  poly_org[1] = origin[1]* C.RefLength;
  poly_org[2] = origin[2]* C.RefLength;
  
  Hostonly_
  {
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
    fprintf(fp,"\t>> Polylib configuration\n\n");
    fprintf(fp,"\tfile name = %s\n\n", C.PolylibConfigName.c_str());
    
    printf("\n---------------------------------------------------------------------------\n\n");
    printf("\t>> Polylib configuration\n\n");
    printf("\tfile name = %s\n\n", C.PolylibConfigName.c_str());
  }
  
  // Polylib: インスタンス取得
  PL = MPIPolylib::get_instance();
  
  // Polylib: ポリゴングループのFactoryクラスを登録
  //PL->set_factory( new MyPolygonGroupFactory() );
  
  // Polylib: 並列計算領域情報を設定
  poly_stat = PL->init_parallel_info(MPI_COMM_WORLD,
                                     poly_org,         // 自ランクの基点座標
                                     (unsigned*)size,  // 自ランクの分割数
                                     poly_gc,          // ガイドセル数
                                     poly_dx           // 格子幅
                                     );
  if ( poly_stat != PLSTAT_OK )
  {
    Hostonly_
    {
      fprintf(fp,"\tRank [%6d]: p_polylib->init_parallel_info() failed.", myRank);
      printf    ("\tRank [%6d]: p_polylib->init_parallel_info() failed.", myRank);
    }
    Exit(0);
  }
  
  
  
  // Polylib: STLデータ読み込みとスケーリング
  TIMING_start(tm_polygon_load);
  
  // ロード
  poly_stat = PL->load_rank0( C.PolylibConfigName );
  
  if( poly_stat != PLSTAT_OK )
  {
    Hostonly_
    {
      printf    ("\tRank [%6d]: p_polylib->load_rank0() failed.", myRank);
      fprintf(fp,"\tRank [%6d]: p_polylib->load_rank0() failed.", myRank);
    }
    Exit(0);
  }
  
  // ポリゴンの修正
  //RepairPolygonData(PL);
  
  
  /* スケーリングする場合のみ表示
  if ( C.Scaling_Factor != 1.0 )
  {
    fprintf(fp,"\n\tScaling Factor           :   %f\n", C.Scaling_Factor);
  }
  
  // スケーリング
  poly_stat = PL->rescale_polygons(C.Scaling_Factor);
  
  if( poly_stat != PLSTAT_OK )
  {
    Hostonly_
    {
      printf    ("\tRank [%6d]: p_polylib->rescale_polygons() failed.", myRank);
      fprintf(fp,"\tRank [%6d]: p_polylib->rescale_polygons() failed.", myRank);
    }
    Exit(0);
  }*/
  
  
  TIMING_stop(tm_polygon_load);
  
  // 階層情報表示 debug brief hierarchy
// ##########
#if 0
  PL->show_group_hierarchy();
  PL->show_group_hierarchy(fp);
#endif
// ##########
  
  
  // ポリゴン情報へのアクセス
  vector<PolygonGroup*>* pg_roots = PL->get_root_groups();
  vector<PolygonGroup*>::iterator it;
  
  // Polygon Groupの数
  C.num_of_polygrp = pg_roots->size();
  
  if ( (C.num_of_polygrp < 1) || (C.num_of_polygrp > C.NoCompo) )
  {
    printf("\tError : Number of polygon group must be greater than 1 and less than NoCompo.\n");
    Exit(0);
  }
  
  // ポリゴングループの属性保持
  PolyPP = new PolygonProperty[C.num_of_polygrp];
  
  Hostonly_
  {
    printf(     "\tNumber of Polygon Group = %d\n\n", C.num_of_polygrp);
    fprintf(fp, "\tNumber of Polygon Group = %d\n\n", C.num_of_polygrp);
    
    printf(     "\t   Polygon Group Label     Medium Alias              Local BC      Element          Area\n");
    printf(     "\t   -------------------------------------------------------------------------------------\n");
    fprintf(fp, "\t   Polygon Group Label     Medium Alias              Local BC      Element          Area\n");
    fprintf(fp, "\t   -------------------------------------------------------------------------------------\n");
  }
  
  
  // ポリゴン情報の表示と管理クラスへのコピー
  int c = 0;
  for (it = pg_roots->begin(); it != pg_roots->end(); it++) {
    
    std::string m_pg = (*it)->get_name();     // グループラベル
    std::string m_mat = (*it)->get_label();   // 媒質ラベル
    std::string m_bc = (*it)->get_type();     // 境界条件ラベル
    int ntria= (*it)->get_group_num_tria();   // ローカルのポリゴン数
    REAL_TYPE area = (*it)->get_group_area(); // ローカルのポリゴン面積
    
    // mat[]の格納順を探す
    int mat_id;
    
    for (int i=1; i<=C.NoCompo; i++)
    {
      if ( FBUtility::compare(m_pg, mat[i].getAlias()) )
      {
        mat_id = i;
        break;
      }
    }
    
    // PolygonにIDを割り当てる
    poly_stat = (*it)->set_all_exid_of_trias(mat_id);
    
    if ( poly_stat != PLSTAT_OK )
    {
      Hostonly_
      {
        printf(     "\tError : Polylib::set_all_exid_of_trias()\n");
        fprintf(fp, "\tError : Polylib::set_all_exid_of_trias()\n");
        Exit(0);
      }
    }
    
    PolyPP[c].setGroup(m_pg);
    PolyPP[c].setMaterial(m_mat);
    PolyPP[c].setBClabel(m_bc);
    PolyPP[c].setLntria(ntria);
    PolyPP[c].setLarea(area);     
    
    if ( numProc > 1 )
    {
      int tmp = ntria;
      if ( paraMngr->Allreduce(&tmp, &ntria, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
      
      REAL_TYPE ta = area;
      if ( paraMngr->Allreduce(&ta, &area, 1, MPI_SUM) != CPM_SUCCESS ) Exit(0);
    }
    
    PolyPP[c].setGntria(ntria);   // グローバルなポリゴン数
    PolyPP[c].setGarea(area);     // グローバルなポリゴン面積
    
    Hostonly_
    {
      printf(    "\t  %20s %16s  %20s %12d  %e\n", PolyPP[c].getGroup().c_str(),
                                                   PolyPP[c].getMaterial().c_str(),
                                                   PolyPP[c].getBClabel().c_str(),
                                                   PolyPP[c].getGntria(),
                                                   PolyPP[c].getGarea());
      fprintf(fp,"\t  %20s %16s  %20s %12d  %e\n", PolyPP[c].getGroup().c_str(),
                                                   PolyPP[c].getMaterial().c_str(),
                                                   PolyPP[c].getBClabel().c_str(),
                                                   PolyPP[c].getGntria(),
                                                   PolyPP[c].getGarea());
    }
    
    c++;
    
// ########## show corrdinates and area
#if 0
    PL->show_group_info(m_pg); //debug
#endif
// ##########
    
  }
  
  delete pg_roots;
  
  Hostonly_
  {
    printf("\n");
    fprintf(fp, "\n");
  }
  
  
  
  
  // 使用メモリ量　基本クラスのみ
  double poly_mem, G_poly_mem;
  G_poly_mem = poly_mem = (double)PL->used_memory_size();
  m_prep += poly_mem;
  m_total+= poly_mem;
  
  displayMemoryInfo(fp, G_poly_mem, poly_mem, "Polygon");
  

  //calcBboxfromPolygonGroup();
  
  
  // Triangle display >> Debug
// ##########
#if 0
  PolylibNS::Vec3f m_min, m_max, t1(poly_org), t2(poly_dx), t3;
  t3.assign((float)size[0]*t2.t[0], (float)size[1]*t2.t[1], (float)size[2]*t2.t[2]);
  m_min = t1 - t2;      // 1層外側まで
  m_max = t1 + t3 + t2; //
  printf("min : %f %f %f\n", m_min.t[0], m_min.t[1], m_min.t[2]);
  printf("max : %f %f %f\n", m_max.t[0], m_max.t[1], m_max.t[2]);
  vector<Triangle*>* trias = PL->search_polygons("Tube", m_min, m_max, false); // false; ポリゴンが一部でもかかる場合
  
  PolylibNS::Vec3f *p, nrl, n;
  c=0;
  vector<Triangle*>::iterator it2;
  for (it2 = trias->begin(); it2 != trias->end(); it2++) {
    p = (*it2)->get_vertex();
    n = (*it2)->get_normal();
    printf("%d : p0=(%6.3e %6.3e %6.3e)  p1=(%6.3e %6.3e %6.3e) p2=(%6.3e %6.3e %6.3e) n=(%6.3e %6.3e %6.3e)\n", c++,
           p[0].t[0], p[0].t[1], p[0].t[2],
           p[1].t[0], p[1].t[1], p[1].t[2],
           p[2].t[0], p[2].t[1], p[2].t[2],
           n.t[0], n.t[1], n.t[2]);
  }
  
  delete trias;  //後始末
#endif
// ##########
  
  // Polylib: STLデータ書き出しテスト >> Debug
// ##########
#if 0
  unsigned poly_out_para = IO_DISTRIBUTE; // 逐次の場合と並列の場合で明示的に切り分けている．あとで，考慮
  string fname;
  
  if ( poly_out_para == IO_GATHER )
  {
    poly_stat = PL->save_rank0( &fname, "stl_b" );
    if ( poly_stat != PLSTAT_OK )
    {
      Hostonly_
      {
        printf("Rank [%d]: p_polylib->save_rank0() failed to write into '%s'.", myRank, fname.c_str());
        fprintf(fp,"Rank [%d]: p_polylib->save_rank0() failed to write into '%s'.", myRank, fname.c_str());
      }
      Exit(0);
    }
  }
  else
  {
    poly_stat = PL->save_parallel( &fname, "stl_b" );
    if ( poly_stat != PLSTAT_OK )
    {
      Hostonly_
      {
        printf("Rank [%d]: p_polylib->save_parallel() failed to write into '%s'.", myRank, fname.c_str());
        fprintf(fp,"Rank [%d]: p_polylib->save_parallel() failed to write into '%s'.", myRank, fname.c_str());
      }
      Exit(0);
    }
  }
#endif
// ##########  
  
  
  
  
  
  // Cutlib
  Hostonly_
  {
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
    fprintf(fp,"\t>> Cut Info\n\n");
    printf("\n---------------------------------------------------------------------------\n\n");
    printf("\t>> Cut Info\n\n");
  }
  
  size_t n_cell[3];
  
  unsigned cut_gc[3] = {guide, guide, guide};
  double cut_org[3];
  double cut_dx[3];
  
  // 有次元に変換 Polylib: 並列計算領域情報　ポリゴンは実スケール
  cut_dx[0]  = (double)(pitch[0] * C.RefLength);
  cut_dx[1]  = (double)(pitch[1] * C.RefLength);
  cut_dx[2]  = (double)(pitch[2] * C.RefLength);
  cut_org[0] = (double)(origin[0]* C.RefLength);
  cut_org[1] = (double)(origin[1]* C.RefLength);
  cut_org[2] = (double)(origin[2]* C.RefLength);

  /*
   
          Outer < | >     Inner      < | > Outer
        +----+----+----+----+-...-+----+----+----+
      0 |    |    |    |    |     |    |    |    |
        |    |    |    |    |     |    |    |    |
        +----+----+----+----+-...-+----+----+----+
        ^ -1   0    1    2    ...   ix  ix+1 ix+2
          sx                                  ex
   
   グリッド情報アクセッサクラスCell(const double org[3], const double pitch[3])
   のコンストラクタに渡す引数
     org[]はインデクス空間で(0,0,0)のセルの最小位置（上図の^の部分）
   
   */
  
  for ( int i=0; i<3; i++) {
    n_cell[i] = (size_t)(size[i] + 2*guide);  // 分割数+ガイドセル両側
    cut_org[i] -= cut_dx[i]*(double)guide;    // ガイドセルを含む領域のマイナス側頂点の座標
  }
  
  size_t size_n_cell = n_cell[0] * n_cell[1] * n_cell[2];
  
  if (size_n_cell*6  > UINT_MAX)
  {
    Hostonly_
    {
      stamped_printf("\n\tError : Product of size[]*6 exceeds UINT_MAX\n\n");
      stamped_fprintf(fp,"\n\tError : Product of size[]*6 exceeds UINT_MAX\n\n");
    }
    Exit(0);
  }
  
  
  
  // Cutlibの配列は各方向(引数)のサイズ
  TIMING_start(tm_init_alloc);
  cutPos = new CutPos32Array(n_cell); // 6*(float)
  cutBid = new CutBid5Array (n_cell); // (int32_t)
  TIMING_stop(tm_init_alloc);
  
  
  // グリッド情報アクセッサクラス
  GridAccessor* GRID = new Cell(cut_org, cut_dx);
  
  
  // 交点計算
  TIMING_start(tm_cutinfo);
  CalcCutInfo(GRID, PL, cutPos, cutBid);
  TIMING_stop(tm_cutinfo);
  
  
  // 使用メモリ量
  double cut_mem, G_cut_mem;
  G_cut_mem = cut_mem = (double)size_n_cell * (double)(6*sizeof(float) + sizeof(int)); // float
  m_prep += cut_mem;
  m_total+= cut_mem;
  
  
  displayMemoryInfo(fp, G_cut_mem, cut_mem, "Cut");
  
  
  // カットとID情報をポイント
  d_cut = (float*)cutPos->getDataPointer();
  d_bid = (int*)cutBid->getDataPointer();
  
  
  // カットの最小値
  minDistance(d_cut, d_bid, fp);
  
}


// #################################################################
// VOF値を気体(0.0)と液体(1.0)で初期化
void FFV::setVOF()
{
  size_t m;
  int s, odr;
  
  int ix = size[0];
  int jx = size[1];
  int kx = size[2];
  int gd = guide;
  
  for (int k=1; k<=kx; k++) {
    for (int j=1; j<=jx; j++) {
      for (int i=1; i<=ix; i++) {
        m = _F_IDX_S3D(i, j, k, ix, jx, kx, gd);
        s = d_bcd[m];
        odr = DECODE_CMP(s);
        
        if ( cmp[odr].getState() == FLUID ) 
        {
          d_vof[m] = ( cmp[odr].getPhase() == GAS ) ? 0.0 : 1.0;
        }
      }
    }
  }
}
