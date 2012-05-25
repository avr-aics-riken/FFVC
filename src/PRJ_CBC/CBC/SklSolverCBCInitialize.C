/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file SklSolverCBCInitialize.C
//@brief CBCソルバークラスのプリプロセッサ
//@author keno, FSI Team, VCAD, RIKEN

#include <stdlib.h>
#include <string>

#include "SklSolverCBC.h"
#include "CompoFraction.h"
#include "fileio/SklVoxDataSet.h"
#include <unistd.h> // for gethostname() and sleep()
#include "util_Path.h"


//@fn int SklSolverCBC::SklSolverInitialize()
//@brief 前処理
int
SklSolverCBC::SklSolverInitialize() {
/* USER WRITE
   Please write solver initialize code.
     Return value  ->   -1 : error
                         0 : forced termination
                         1 : normal return
*/

  SklParaComponent* para_cmp = SklGetParaComponent();
  const SklParaManager* para_mng = para_cmp->GetParaManager();
  
  unsigned long TotalMemory, PrepMemory, G_TotalMemory, G_PrepMemory, tmp_memory;
  REAL_TYPE  *ws, *v, *t, *p;
  int       *mid;
  unsigned  *bcv, *bh1, *bh2, *bcp, *bcd;
  unsigned  n=0;
  FILE* fp = NULL;
  REAL_TYPE flop_task=0.0;
  float *cvf=NULL; /// コンポーネントの体積率
  int  d_length;
  
  ws = v = t = p = NULL;
  mid = NULL;
  bcd = bcv = bh1 = bh2 = bcp = NULL;
  TotalMemory = PrepMemory = 0;
  G_TotalMemory = G_PrepMemory = 0;
  
  
  // 固定パラメータ
  fixed_parameters();
  
  
  // 前処理段階のみに使用するオブジェクトをインスタンス
  VoxInfo Vinfo;
  ParseBC     B;
  ParseMat    M;

  // 並列処理時のランク番号をセット．デフォルトで，procGrp = 0;
  pn.ID = para_mng->GetMyID(pn.procGrp);
  
  int m_np = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &m_np);
  
  // 並列処理モード
  setParallelism();
  std::string para_mode;
  
  switch (C.Parallelism) {
    case Control::Serial:
      para_mode = "Serial";
      break;
      
    case Control::OpenMP:
      para_mode = "OpenMP";
      break;
      
    case Control::FlatMPI:
      para_mode = "FlatMPI";
      break;
      
    case Control::Hybrid:
      para_mode = "Hybrid";
      break;
      
    default:
      Exit(0);
      break;
  }
  
  // condition fileのオープン
  Hostonly_ {
    if ( !(fp=fopen("condition.txt", "w")) ) {
      stamped_printf("\tSorry, can't open 'condition.txt' file. Write failed.\n");
      return -1;
    }
  }

  // メッセージ表示
  Hostonly_ {
    char buf[LABEL];
    memset(buf, 0, sizeof(char)*LABEL);
    strcpy(buf, "Welcome to V-SPHERE::CBC");
    FBUtility::printVersion(fp, buf, VERS_CBC);
    FBUtility::printVersion(mp, buf, VERS_CBC);
		memset(buf, 0, sizeof(char)*LABEL);
    strcpy(buf, "FlowBase                ");
    FBUtility::printVersion(fp, buf, FB_VERS);
    FBUtility::printVersion(mp, buf, FB_VERS);
  }
  
  // 例題の種類を取得し，C.Mode.Exampleにフラグをセットする
  getXMLExample(&C);
  
  // 組み込み例題クラスの実体をインスタンスし，*Exにポイントする
  connectExample(&C);
  
  // 組み込み例題クラス名を表示
  Hostonly_ {
    Ex->printExample(fp, Ex->getExampleName());
  }
  
  // コンフィギュレーション 各クラスで必要なオブジェクトのポインタを渡す -----------------------------------------------------
  if ( !C.receiveCfgPtr(m_solvCfg) ) {
    Hostonly_ stamped_printf("\tError during sending an object pointer of XML tree to Control class\n");
    return -1;
  }
  B.receiveCfgPtr(m_solvCfg);
  
  // XMLファイルのエントリタグをチェックする -----------------------------------------------------
  if ( !SklUtil::chkXMLTopTag(STEER) ) {
    Hostonly_ stamped_printf("\tSklUtil::chkXMLTopTag()\n");
    return -1;
  }
  if ( !SklUtil::chkXMLTopTag(PARAMETER) ) {
    Hostonly_ stamped_printf("\tSklUtil::chkXMLTopTag()\n");
    return -1;
  }
  
  // バージョン情報を取得し，ソルバークラスのバージョンと一致するかをチェックする
  C.getXML_Version();
  
  if ( C.version != (unsigned)VERS_CBC ) {
    Hostonly_ {
      fprintf(mp, "\t##### Version of XML description (%d) is NOT compliant with CBC ver. %d #####\n", C.version, VERS_CBC);
      fprintf(fp, "\t##### Version of XML description (%d) is NOT compliant with CBC ver. %d #####\n", C.version, VERS_CBC);
    }
    return -1;
  }
  if ( C.FB_version != (unsigned)FB_VERS ) {
    Hostonly_ {
      fprintf(mp, "\t##### Version of XML description (%d) is NOT compliant with FB ver. %d #####\n", C.FB_version, FB_VERS);
      fprintf(fp, "\t##### Version of XML description (%d) is NOT compliant with FB ver. %d #####\n", C.FB_version, FB_VERS);
    }
    return -1;
  }

  // 最初のXMLパラメータの取得
  C.getXML_Steer_1(&DT);
  
  // FiliIOのモードを修正
  if ( m_np == 1 ) {
    C.FIO.IO_Input  = IO_GATHER;
    C.FIO.IO_Output = IO_GATHER;
    Hostonly_ printf("\tMode of Parallel_Input/_Output was changed because of serial execution.\n");
  }

  // 計算領域全体のサイズ，並列計算時のローカルのサイズ，コンポーネントのサイズなどを設定する -----------------------------------------------------
  VoxelInitialize();

  // 組み込み例題クラスへ変数の値をコピー，固有のパラメータ取得，領域設定 -----------------------------------------------------
  Ex->setControlVars (&C);
  
  // 並列情報のコピー
  BC.setParallelInfo  (pn);
  F.setParallelInfo   (pn);
  MO.setParallelInfo  (pn);
  Ex->setParallelInfo (pn);  // 続く処理より先に初期化しておく

  // XMLパラメータの取得
  C.getXML_Steer_2(IC, &RF);
  
  // 組み込み例題の固有パラメータ
  if ( !Ex->getXML(m_solvCfg, &C) ) Exit(0);

  // ソルバークラスのノードローカルな変数の設定 -----------------------------------------------------
  ix      = (int*)&C.imax;
  jx      = (int*)&C.jmax;
  kx      = (int*)&C.kmax;
  x0      = &C.org[0];
  y0      = &C.org[1];
  z0      = &C.org[2];
  dh0     = &C.dh;
  dh      = &C.dh;
  guide   = C.guide;
  gc      = (int*)&C.guide;
  sz[0]   = ixc = (int)C.imax;
  sz[1]   = jxc = (int)C.jmax;
  sz[2]   = kxc = (int)C.kmax;
  size[0] = C.imax;
  size[1] = C.jmax;
  size[2] = C.kmax;
  
  // タイミング測定の初期化
  if ( C.Mode.Profiling != OFF ) {
    ModeTiming = ON;
    TIMING__ PM.initialize(tm_END);
    TIMING__ PM.setRankInfo(pn.ID);
    TIMING__ PM.setParallelMode(para_mode, C.num_thread, C.num_process);
    set_timing_label();
  }
  
  // タイミング測定開始
  TIMING_start(tm_init_sct); 

  // 前処理に用いるデータクラスのアロケート -----------------------------------------------------
  TIMING_start(tm_init_alloc); 
  allocArray_prep(PrepMemory, TotalMemory);
  TIMING_stop(tm_init_alloc); 
  
  if( !(ws  = dc_ws->GetData()) )   Exit(0);
  if( !(mid = dc_mid->GetData()) )  Exit(0);
  if( !(bcd = dc_bcd->GetData()) )  Exit(0);
  if( !(bcp = dc_bcp->GetData()) )  Exit(0);
  if( !(bcv = dc_bcv->GetData()) )  Exit(0);
  if ( C.isHeatProblem() ) {
    if( !(bh1 = dc_bh1->GetData()) )  Exit(0);
    if( !(bh2 = dc_bh2->GetData()) )  Exit(0);
  }
  
  
  // ファイルからIDを読み込む，または組み込み例題クラスでID情報を作成 --------------------------------------------------------
  Hostonly_ {
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
    fprintf(fp,"\t>> Voxel file information\n\n");
    fprintf(mp,"\n---------------------------------------------------------------------------\n\n");
    fprintf(mp,"\t>> Voxel file information\n\n");
  }
  
  TIMING_start(tm_voxel_prep_sct);
  
  REAL_TYPE* tmp_cut = NULL;

  switch (C.Mode.Example) {
      
    case id_Polygon:
      // パラメータの取得
      C.getXML_Polygon();
      
      // PolylibとCutlibのセットアップ
      setup_Polygon2CutInfo(PrepMemory, TotalMemory, fp);
      
      if( !(tmp_cut  = dc_cut->GetData()) ) Exit(0);
      
      if ( !C.isCDS() ) {
        unsigned zc = Solid_from_Cut(mid, tmp_cut, id_of_solid);
        Hostonly_ printf("\tGenerated Solid cell from cut = %d\n", zc);
      }
      break;
      
    case id_Users:
      // バイナリ
      if ( !C.isCDS() ) {
        TIMING_start(tm_voxel_load);
        if ( C.vxFormat == Control::Sphere_SVX ) {
          F.readSVX(this, fp, "SphereSVX", G_size, guide, dc_mid);
        }
        else {
          F.readSBX(this, fp, "SphereSBX", G_size, guide, dc_mid);
        }
        TIMING_stop(tm_voxel_load);
      }
      else { // カットなし
        Hostonly_ printf("\t>> Binary only for 'Users' mode\n\n");
        Exit(0);
      }
      break;
      
    case id_Sphere:
      if ( !C.isCDS() ) {
        Ex->setup(mid, &C, G_org);
      }
      else {
        // cutをアロケートし，初期値1.0をセット
        setup_CutInfo4IP(PrepMemory, TotalMemory, fp);
        Ex->setup_cut(mid, &C, G_org, dc_cut->GetData());
      }
      break;
      
    default: // ほかのIntrinsic problems
      if ( C.isCDS() ) {
        setup_CutInfo4IP(PrepMemory, TotalMemory, fp);
      }
      Ex->setup(mid, &C, G_org);
      break;
  }

  
  // midのガイドセル同期
  if( !dc_mid->CommBndCell(guide) ) return -1;
  

  // 領域情報の表示
  Hostonly_ {
    fprintf(mp,"\n---------------------------------------------------------------------------\n\n");
    fprintf(mp,"\n\t>> Global Domain Information\n\n");
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
    fprintf(fp,"\n\t>> Global Domain Information\n\n");
    C.printDomainInfo(mp, fp, G_size, G_org, G_Lbx);
  }
  
  // メモリ消費量の情報を表示
  Hostonly_ {
    fprintf(mp,"\n---------------------------------------------------------------------------\n\n");
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
  }
  G_PrepMemory = PrepMemory;
  if( para_cmp->IsParallel() ) {
    tmp_memory = G_PrepMemory;
    para_cmp->Allreduce(&tmp_memory, &G_PrepMemory, 1, SKL_ARRAY_DTYPE_ULONG, SKL_SUM, pn.procGrp);
  }
  Hostonly_  {
    FBUtility::displayMemory("prep", G_PrepMemory, PrepMemory, fp, mp);
  }
  
  // ∆tの決め方とKindOfSolverの組み合わせで無効なものをはねる
  if ( !DT.chkDtSelect() ) {
    Hostonly_ printf("\tCombination of specified 'Time_Increment' and 'Kind_of_Solver' is not permitted.\n");
    return -1;
  }
  
  // XMLから C.NoBC, C.NoID, C.NoCompoを取得
  // KOSオプションとパラメータの整合性をチェック
  // BCのテーブル取得と表示, NoMaterialをセットする
  Hostonly_ {
    fprintf(mp,"\n---------------------------------------------------------------------------\n\n");
    fprintf(mp, "\t>> Table Information\n\n");
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
    fprintf(fp, "\t>> Table Information\n\n");
  }
  setIDtables(&B, fp, mp);
  
  
  // VoxInfoクラスへ値をセット
  Vinfo.setNoCompo_BC(C.NoBC, C.NoCompo);
  
  
  // XMLから得られた内部BCコンポーネント数を表示
  Hostonly_ {
    fprintf(mp,"\n---------------------------------------------------------------------------\n\n");
    fprintf(mp,"\t>> Components\n\n");
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
    fprintf(fp,"\t>> Components\n\n");
    
    C.printNoCompo(fp);
    fprintf(fp,"\n"); fflush(fp);
    C.printNoCompo(mp);
    fprintf(mp,"\n"); fflush(mp);
  }

  // VoxInfoクラスで利用する変数をコピー
  Vinfo.setParallelInfo(pn);
  Vinfo.setControlVars(size, guide);
  
  
  // Fill
  if ( (C.Mode.Example == id_Polygon) ) { //&& C.isCDS() ) {
    
    Hostonly_ {
      fprintf(mp,"\n---------------------------------------------------------------------------\n\n");
      fprintf(mp,"\t>> Fill\n\n");
      fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
      fprintf(fp,"\t>> Fill\n\n");
    }
    
    int target_id = C.Mode.Base_Medium;
    unsigned long fill_count = size[0] * size[1] * size[2];
    
    //unsigned isolated_cell = Vinfo.test_opposite_cut(cut_id, mid, id_of_solid);
    //Hostonly_ printf("Filled cut = %d\n", isolated_cell);
    
    unsigned long tmp_fc = fill_count;
    if( para_cmp->IsParallel() ) {
      para_cmp->Allreduce(&tmp_fc, &fill_count, 1, SKL_ARRAY_DTYPE_ULONG, SKL_SUM, pn.procGrp);
    }
    Hostonly_ printf("\tInitial fill count     : %15ld\n\n", fill_count);
    
    int seed[3];
    seed[0] = 1;
    seed[1] = 1;
    seed[2] = 1;
    
    Hostonly_ {
      if ( !Vinfo.paint_first_seed(mid, seed, target_id) ) {
        printf("Failed first painting\n");
      }
    }
    
    // first fill
    fill_count--;
    
    int c=0;
    while (fill_count > 0) {
      
      unsigned fc = Vinfo.fill_cell_edge(dc_bid->GetData(), mid, dc_cut->GetData(), target_id, id_of_solid);
      
      unsigned t_fc = fc;

      if ( m_np > 1 ) {
        VoxInfo::uint_sum_Allreduce(&t_fc, &fc);
      }
      
      // 同期
      dc_mid->CommBndCell(guide);
      dc_bid->CommBndCell(guide);
      dc_cut->CommBndCell(guide);
      
      if ( fc == 0 ) break;
      
      fill_count -= (unsigned long)fc;
      
      Hostonly_ printf("\t%4d : Try edge fill with ID = %d : %ld\n", ++c, target_id, fill_count);
      
    }
    Hostonly_ printf("\tEdge  fill with ID = %d : %15ld\n\n", target_id, fill_count);
    
    // 固体に変更
    c = 0;
    while ( fill_count > 0 ) {
      
      unsigned fc = Vinfo.fill_inside(mid, id_of_solid);
      
      unsigned t_fc = fc;

      if ( m_np > 1 ) {
        VoxInfo::uint_sum_Allreduce(&t_fc, &fc);
      }
      
      if ( fc == 0 ) break;
      
      fill_count -= (unsigned long)fc;
      //Hostonly_ printf("\t%4d : Try solid fill with ID = %d : %ld\n", ++c, id_of_solid, fill_count);
      
    }
    Hostonly_ printf("\tSolid fill with ID = %d : %15ld\n\n", id_of_solid, fill_count);

    // チェック
    Ex->writeSVX(mid, &C);
    
    // 確認 paintedは未ペイントセルがある場合に1
    unsigned painted = Vinfo.check_fill(mid);
    unsigned t_painted = painted;

    if ( m_np > 1 ) {
      VoxInfo::uint_sum_Allreduce(&t_painted, &painted);
    }
    
    if ( painted  != 0 ) {
      Hostonly_ printf("\tFill operation is done, but still remains unpainted cells.\n");
      Exit(0);
    }
  }
  
  
  // バイナリボクセルのスキャン
  VoxScan(&Vinfo, &B, mid, fp);
  
  // スキャンしたセルIDの情報を表示する
  Hostonly_ {
    fprintf(mp, "\n---------------------------------------------------------------------------\n\n");
    fprintf(mp,"\t>> Information of Scaned Voxel\n\n");
    fprintf(fp, "\n---------------------------------------------------------------------------\n\n");
    fprintf(fp,"\t>> Information of Scaned Voxel\n\n");
    
		Vinfo.printScanedCell(fp);
		fflush(fp);
    Vinfo.printScanedCell(mp);
		fflush(mp);
  }

  // ボクセルモデルのIDがXMLに記述されたIDに含まれていること
	Hostonly_ {
		if ( !Vinfo.chkIDconsistency(B.get_IDtable_Ptr(), C.NoID) ) {
			stamped_printf("\tID in between XML and Voxel scaned is not consistent\n");
			return -1;
		}
	}

  // CompoList/MaterialListクラスをインスタンス．[0]はダミーとして利用しないので，配列の大きさはプラス１する
  cmp = new CompoList[C.NoCompo+1];
  mat = new MaterialList[C.NoMaterial+1];
  BC.setWorkList(cmp, mat);
  Vinfo.setWorkList(cmp, mat);

  
  // ParseBCクラスのセットアップ，CompoListの設定，外部境界条件の読み込み保持
  setBCinfo(&B);

  
  // ガイドセル上にXMLで指定するセルIDを代入する．周期境界の場合の処理も含む．
  for (int face=0; face<NOFACE; face++) {
    Vinfo.adjCellID_on_GC(face, dc_mid, BC.get_OBC_Ptr(face)->get_BCtype(), 
                         BC.get_OBC_Ptr(face)->get_GuideID(), BC.get_OBC_Ptr(face)->get_PrdcMode());
  }

  
  // Cell_Monitorの指定がある場合，モニタ位置をセット
  if ( (C.Sampling.log == ON) && (C.isMonitor() == ON) ) setShapeMonitor(mid);

  
  // CDSの場合，WALLとSYMMETRICのときに，カットを外部境界に接する内部セルに設定
  if ( C.isCDS() ) {
    Vinfo.setOBC_Cut(&BC, dc_cut->GetData());
  }

  // ParseMatクラスをセットアップし，媒質情報をXMLから読み込み，媒質リストを作成する
  Hostonly_  {
    fprintf(mp,"\n---------------------------------------------------------------------------\n\n");
    fprintf(mp,"\n\t>> Medium List\n\n");
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
    fprintf(fp,"\n\t>> Medium List\n\n");
  }
  setMaterialList(&B, &M, mp, fp);
  

#if 0
  // カット情報からボクセルモデルの固体をIDセット　デバッグ用
  if ( C.isCDS() ) {
    unsigned zc = Vinfo.Solid_from_Cut(mid, cut);
    Hostonly_ printf("\tCut cell = %d\n", zc);
  }
#endif
  
  // セルIDのノード間同期
  if( !dc_mid->CommBndCell(guide) ) return -1;
  

  // HEX/FANコンポーネントの形状情報からBboxと体積率を計算
  if ( C.isVfraction() ) {
    TIMING_start(tm_init_alloc);
    allocArray_compoVF(PrepMemory, TotalMemory);
    TIMING_stop(tm_init_alloc); 
    if( !(cvf = dc_cvf->GetData()) )  Exit(0);
    
    setComponentVF(cvf);
  }
  
  // コンポーネントのローカルインデクスを保存
  setLocalCmpIdx_Binary();
  
  
  // 内部周期境界の場合のガイドセルのコピー処理
  Vinfo.adjCellID_Prdc_Inner(dc_mid);

  // 媒質数とKindOfSolverの整合性をチェックする
  if ( !chkMediumConsistency() ) {
    Hostonly_ stamped_printf("\tchkMediumConsistency()\n");
    return -1;
  }


  // BCIndexへのエンコード処理
  Hostonly_  {
    fprintf(mp,"\n---------------------------------------------------------------------------\n\n");
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
  }
  
  // BCIndexにビット情報をエンコードとコンポーネントインデクスの再構築
  VoxEncode(&Vinfo, &M, mid, cvf);

  // CDSのとき，ポリゴンからVBCのコンポーネント情報を設定
  if ( C.isCDS() ) {
    setBbox_of_VIBC_Cut(&Vinfo, bcv);
  }
  
  // 体積力を使う場合のコンポーネント配列の確保
  allocComponentArray(PrepMemory, TotalMemory, fp);
  
  // コンポーネントの体積率を8bitで量子化し，圧力損失コンポの場合にはFORCING_BITをON > bcdにエンコード
  Vinfo.setCmpFraction(cmp, bcd, cvf);
  
#if 0
  // CompoListとMaterialListの関連を表示
  Hostonly_ M.dbg_printRelation(mp, fp, cmp);
#endif
  
  // CompoListとMaterialListのstate(Fluid / Solid)を比較しチェックする
  if ( !M.chkStateList(cmp)) {
    Hostonly_ stamped_printf("\tParseMat::chkStateList()\n");
    return -1;
  }

  // RefIDがCompoList中にあるかどうかをチェックする
  if ( B.isIDinCompo(C.RefID, C.NoCompo+1) ) { // falseのとき，重複がある，つまりIDがある
    Hostonly_ {
      fprintf(mp, "RefID[%d] is not listed in Medium XML file.\n", C.RefID);
      fprintf(fp, "RefID[%d] is not listed in Medium XML file.\n", C.RefID);
      stamped_printf("\tParseBC::isIDinCompo()\n");
    }
    return -1;
  }


  // 周期境界条件が設定されている場合のBCIndexの周期条件の強制同期
  BC.setBCIperiodic(dc_bcd);
  BC.setBCIperiodic(dc_bcp);
  BC.setBCIperiodic(dc_bcv);
  if ( C.isHeatProblem() ) {
    BC.setBCIperiodic(dc_bh1);
    BC.setBCIperiodic(dc_bh2);
  }
  
  // bcd/bcp/bcv/bchの同期
  dc_bcd->CommBndCell(guide);
  dc_bcp->CommBndCell(guide);
  dc_bcv->CommBndCell(guide);
  if ( C.isHeatProblem() ) {
    dc_bh1->CommBndCell(guide);
    dc_bh2->CommBndCell(guide);
  }
  
  // 法線計算のためのワーク配列 >> 不要にする
  if ( C.NoBC != 0 ) {
    Vinfo.alloc_voxel_nv((C.NoBC+1)*3);
    
    // コンポーネントで指定されるID面の法線を計算，向きはblowing/suctionにより決まる．　bcdをセットしたあとに処理
    for (unsigned n=1; n<=C.NoBC; n++) {
      if ( C.Mode.Example == id_Polygon ) {
        Vinfo.get_Compo_Area_Cut(n, PL);
      }
      else {
        //Vinfo.cal_Compo_Area_Normal(n, bcd, bcv, bh1, C.dh*C.RefLength, &compo_global_bbox[n*6]);
      }
    }
  }
  
  // 無次元数などの計算パラメータを設定する．MaterialListを決定した後，かつ，SetBC3Dクラスの初期化前に実施すること
  // 代表物性値をRefIDの示す媒質から取得
  // Δt=constとして，無次元の時間積分幅 deltaTを計算する．ただし，一定幅の場合に限られる．不定幅の場合には別途考慮の必要
  DT.set_Vars(C.KindOfSolver, C.Unit.Param, (double)C.dh, (double)C.Reynolds, (double)C.Peclet);

  // 無次元速度1.0を与えてdeltaTをセットし，エラーチェック
  switch ( DT.set_DT(1.0) ) {
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
  
  // Δtをフレームワークへセット
  SklSetDeltaT( DT.get_DT() );
  
  // SPHEREフレームワークにmaxStepを設定する >> V-Sphereがステップによる制御のみ．変化するΔtを用いる場合には注意
  if ( !C.Interval[Interval_Manager::tg_compute].initTrigger(SklGetTotalStep(), 
                                                             (double)SklGetTotalTime(), 
                                                             DT.get_DT(), Interval_Manager::tg_compute, 
                                                             (double)(C.RefLength/C.RefVelocity)) ) {
    Hostonly_ printf("\t Error : Computation Period is asigned to zero.\n");
    Exit(0);
  }
  C.LastStep = C.Interval[Interval_Manager::tg_compute].getIntervalStep();
  SklSetMaxStep( C.LastStep );
  
  // C.Interval[Interval_Manager::tg_compute].initTrigger()で初期化後
  C.setParameters(mat, cmp, B.get_NoBaseBC(), B.get_BaseBC_Ptr(), &RF);
  
  // 媒質による代表パラメータのコピー
  B.setRefValue(mat, cmp, &C);
  
  // 必要なパラメータをSetBC3Dクラスオブジェクトにコピーする >> C.setParameters()の後
  BC.setControlVars(&C, mat, cmp, &RF, Ex);
  
  if ( C.NoBC == 0 ) {
    Hostonly_ printf("\tNo Inner Boundary Conditions\n");
  }
  
  // パラメータの無次元表示（正規化）に必要な参照物理量の設定
  B.setRefMedium(mat, &C);
  
#if 0
  // チェックのため，全計算セルのBCIndexの内容を表示する
  if ( !Vinfo.dbg_chkBCIndexP(bcd, bcp, "BCindex.txt") ) {
    Hostonly_ stamped_printf("\tVoxInfo::dbg_chkBCIndexP()\n");
    return -1;
  }
#endif
  
  // 温度計算の場合の初期値指定
  if ( C.isHeatProblem() ) {
    B.getXML_Medium_InitTemp();
  }
  
  // set phase 
  if ( C.BasicEqs == INCMP_2PHASE ) {
    B.getXML_Phase();
  }
  
  // CompoListの内容を表示する
  Hostonly_  {
    fprintf(mp,"\n---------------------------------------------------------------------------\n\n");
    fprintf(mp,"\t>> Component List\n\n");
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
    fprintf(fp,"\t>> Component List\n\n");
    M.chkList(mp, fp, cmp, C.BasicEqs);
  }
  
  // セル数の情報を表示する
  Hostonly_ {
    REAL_TYPE cr = (REAL_TYPE)G_Wcell/ ( (REAL_TYPE)G_size[0] * (REAL_TYPE)G_size[1] * (REAL_TYPE)G_size[2]) *100.0;
    fprintf(mp, "\tThis voxel includes %4d solid %s  [Solid cell ratio inside computational domain : %8.4f percent]\n\n", 
            C.NoMediumSolid, (C.NoMediumSolid>1) ? "IDs" : "ID", cr);
    fprintf(fp, "\tThis voxel includes %4d solid %s  [Solid cell ratio inside computational domain : %8.4f percent]\n\n", 
            C.NoMediumSolid, (C.NoMediumSolid>1) ? "IDs" : "ID", cr);
  }
  
  // 外部境界面の開口率を計算する
  Vinfo.countOpenAreaOfDomain(bcd, C.OpenDomain);
  
  Hostonly_ {
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
    fprintf(mp,"\n---------------------------------------------------------------------------\n\n");
    fprintf(mp, "\n\n");
  }
  
  // Monitor Listの処理 --------------------------------------------
  MO.setControlVars(bcd, G_org, G_Lbx, C.org, C.dx, C.Lbx, size, guide,
                    C.RefVelocity, C.BaseTemp, C.DiffTemp, C.RefDensity, C.RefLength, C.BasePrs,
                    C.Unit.Temp, C.Mode.Precision, C.Unit.Prs);
  
  
  // モニタ機能がONの場合に，パラメータを取得し，セットの配列を確保する
  if ( C.Sampling.log == ON ) getXML_Monitor(m_solvCfg, &MO);
  
  
  // モニタリストが指定されている場合に，プローブ位置をID=255としてボクセルファイルに書き込む
  if (C.Sampling.log == ON ) {
    MO.write_ID(mid);
  }

  // 内部境界条件として指定されたモニタ設定を登録
   if ( (C.Sampling.log == ON) && (C.isMonitor() == ON) ) MO.setInnerBoundary(cmp, C.NoBC);

  // 組み込み例題 or MonitorListの場合に，svxファイルを出力する．
  if ( (C.Mode.Example != id_Users) || ( (C.Mode.Example == id_Users) && (C.Sampling.log == ON) ) ) {
    Hostonly_ printf("\n\twrite ID which includes Monitor List ID\n\n");
    
    // 性能測定モードがオフのときのみ出力
    if ( C.Hide.PM_Test == OFF ) Ex->writeSVX(mid, &C); // writeSVX(); ユーザ問題の場合には，単にtrueを返す
  }
  
  // mid[]を解放する  ---------------------------
  DataMngr.DeleteDataObj("mid");

  
  // コンポーネントのグローバルインデクス情報を取得
  setGlobalCmpIdx();
  
  // コンポーネントの内容リストを表示する
  Hostonly_ {
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
    fprintf(fp,"\t>> Component Information\n");
    fprintf(mp,"\n---------------------------------------------------------------------------\n\n");
    fprintf(mp,"\t>> Component Information\n");
    B.printCompoInfo(mp, fp, Vinfo.get_vox_nv_ptr(), compo_global_bbox, mat);
  }

  // コンポーネント数がゼロの場合のチェック
  for (unsigned n=1; n<=C.NoBC; n++) {
    if ( cmp[n].getElement() == 0 ) {
      Hostonly_ printf("\tError : No element was found in Component[%d]\n", n);
      fflush(stdout);
      Exit(0);
    }
  }

  // Check consistency of boundary condition
  for (unsigned n=1; n<=C.NoBC; n++) {
    if ( cmp[n].getType() == HT_SN ) {
      if ( (C.KindOfSolver == FLOW_ONLY) || (C.KindOfSolver == THERMAL_FLOW) || (C.KindOfSolver == SOLID_CONDUCTION) ) {
        Hostonly_ printf("\tInconsistent parameters of combination between Kind of Solver and Heat Transfer type SN. Check QBCF\n");
        fflush(stdout);
        Exit(0);
      }
    }
    if ( cmp[n].getType() == HT_SF ) {
      if ( (C.KindOfSolver == FLOW_ONLY) || (C.KindOfSolver == THERMAL_FLOW_NATURAL) || (C.KindOfSolver == SOLID_CONDUCTION) ) {
        Hostonly_ printf("\tInconsistent parameters of combination between Kind of Solver and Heat Transfer type SF. Check QBCF\n");
        fflush(stdout);
        Exit(0);
      }
    }
  }

  // 各ノードの領域情報をファイル出力
  gather_DomainInfo();

  
  TIMING_stop(tm_voxel_prep_sct);
  // ここまでがボクセル準備の時間セクション

#if 0
  write_distance(cut);
#endif
  
  
  // 計算に用いる配列のアロケート ----------------------------------------------------------------------------------
  TIMING_start(tm_init_alloc);
  allocArray_main (TotalMemory);
  
  //allocArray_Collocate (TotalMemory);
  
  if ( C.LES.Calc == ON ) {
    allocArray_LES (TotalMemory);
  }
  
  if ( C.isHeatProblem() ) {
    allocArray_heat (TotalMemory);
  }
  
  if ( C.AlgorithmF == Control::Flow_FS_RK_CN ) {
    allocArray_RK (TotalMemory);
  }
  
  if ( (C.AlgorithmF == Control::Flow_FS_AB2) || (C.AlgorithmF == Control::Flow_FS_AB_CN) ) {
    allocArray_AB2 (TotalMemory);
  }
  
  if ( C.BasicEqs == INCMP_2PHASE ) {
    allocArray_interface(TotalMemory);
  }

  // 時間平均用の配列をアロケート
  allocArray_average (TotalMemory, fp);


  
  TIMING_stop(tm_init_alloc);

  
  // サブドメインのbbox情報を集約する

  // 並列時のみ
  if ( m_np > 1 ) {
    int* g_bbox_st = new int[3*m_np];
    int* g_bbox_ed = new int[3*m_np];
    
    int st_i, st_j ,st_k, ed_i, ed_j, ed_k;

    // head and tail index
    for (int n=0; n<m_np; n++){
      
      // Fortran index
      st_i = para_mng->GetVoxelHeadIndex(n, 0, pn.procGrp) + 1;
      st_j = para_mng->GetVoxelHeadIndex(n, 1, pn.procGrp) + 1;
      st_k = para_mng->GetVoxelHeadIndex(n, 2, pn.procGrp) + 1;
      ed_i = para_mng->GetVoxelTailIndex(n, 0, pn.procGrp) + 1;
      ed_j = para_mng->GetVoxelTailIndex(n, 1, pn.procGrp) + 1;
      ed_k = para_mng->GetVoxelTailIndex(n, 2, pn.procGrp) + 1;
      
      if ( (g_bbox_st[3*n+0] = st_i) < 1 ) Exit(0);
      if ( (g_bbox_st[3*n+1] = st_j) < 1 ) Exit(0);
      if ( (g_bbox_st[3*n+2] = st_k) < 1 ) Exit(0);
      if ( (g_bbox_ed[3*n+0] = ed_i) < 1 ) Exit(0);
      if ( (g_bbox_ed[3*n+1] = ed_j) < 1 ) Exit(0);
      if ( (g_bbox_ed[3*n+2] = ed_k) < 1 ) Exit(0);
    }

    // DFIクラスの初期化
    if ( !DFI.init((int*)G_size, num_div_domain, (int)C.GuideOut, C.Start, g_bbox_st, g_bbox_ed) ) Exit(0);
    
    // host name
    for (int n=0; n<m_np; n++){
      const char* host = para_mng->GetHostName(n, pn.procGrp);
      if ( !host ) Exit(0);
      
      DFI.copy_hostname(host, n);
    }
      
    delete [] g_bbox_st; g_bbox_st = NULL;
    delete [] g_bbox_ed; g_bbox_ed = NULL;
  }
 

  
  
  
  
  // スタート処理 瞬時値と平均値に分けて処理　------------------
  Hostonly_ fprintf(mp,"\n---------------------------------------------------------------------------\n\n");
  Hostonly_ fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
  
  TIMING_start(tm_restart);
  
  if ( C.Start == initial_start) { // 初期スタートのステップ，時間を設定する
    
    SklSetBaseStep(0);
    SklSetBaseTime(0.0);
    
    // V00の値のセット．モードがtrueの場合はV00[0]=1.0に設定，そうでなければtmに応じた値
    if ( SklIsCheckMode() ) RF.setV00((double)SklGetTotalTime(), true);
    else                    RF.setV00((double)SklGetTotalTime());
    
    double g[4];
    RF.copyV00(g);
    for (int i=0; i<4; i++) v00[i]=(REAL_TYPE)g[i];
    
  }
  else if ( C.Start == restart) { // 同一解像度のリスタート
    
    Hostonly_ fprintf(mp, "\t>> Restart from Previous Calculated Results\n\n");
    Hostonly_ fprintf(fp, "\t>> Restart from Previous Calculated Results\n\n");
    
    Restart(fp, flop_task); 

  }
  else if ( C.Start == coarse_restart) { // 粗い格子からのリスタート
    
    Hostonly_ fprintf(mp, "\t>> Restart from Previous Results on Coarse Mesh\n\n");
    Hostonly_ fprintf(fp, "\t>> Restart from Previous Results on Coarse Mesh\n\n");
    

    // 粗い格子のファイルをロードし、内挿処理を行う
    Restart_coarse(fp, flop_task);
    
    Hostonly_ fprintf(mp,"\n");
    Hostonly_ fprintf(fp,"\n");
  }
  
  TIMING_stop(tm_restart);
  


  // セッションの初期時刻をセット
  for (int i=0; i<Interval_Manager::tg_END; i++) {
    C.Interval[i].setTime_init( (double)SklGetBaseTime() );
  }

  // インターバルの初期化
  double m_dt    = (double)SklGetDeltaT();
  double m_tm    = (double)SklGetTotalTime();
  unsigned m_stp = SklGetTotalStep();
  REAL_TYPE tm = (REAL_TYPE)m_tm;

  if ( !C.Interval[Interval_Manager::tg_console].initTrigger(m_stp, m_tm, m_dt, Interval_Manager::tg_console) ) {  // 基本履歴のコンソールへの出力
    Hostonly_ printf("\t Error : Interval for Console output is asigned to zero.\n");
    Exit(0);
  }
  if ( !C.Interval[Interval_Manager::tg_history].initTrigger(m_stp, m_tm, m_dt, Interval_Manager::tg_history) ) {  // 履歴のファイルへの出力
    Hostonly_ printf("\t Error : Interval for History output is asigned to zero.\n");
    Exit(0);
  }
  if ( !C.Interval[Interval_Manager::tg_instant].initTrigger(m_stp, m_tm, m_dt, Interval_Manager::tg_instant) ) {  // 瞬時値ファイル
    Hostonly_ printf("\t Error : Interval for Instantaneous output is asigned to zero.\n");
    Exit(0);
  }
  if ( C.Mode.Average == ON ) {
    //if ( !C.Interval[Interval_Manager::tg_average].initTrigger(m_stp, m_tm, m_dt, Interval_Manager::tg_average) ) {  // 平均値ファイル
    //  Hostonly_ printf("\t Error : Interval for Average output is asigned to zero.\n");
    //  Exit(0);
    //}
    // tg_averageの初期化はLoop中で行う．平均値開始時刻が未知のため．
    if ( !C.Interval[Interval_Manager::tg_avstart].initTrigger(m_stp, m_tm, m_dt, Interval_Manager::tg_avstart) ) {  // 平均値開始
      Hostonly_ printf("\t Error : Interval for Average start is asigned to zero.\n");
      Exit(0);
    }
  }
  if ( C.Sampling.log == ON ) {
    if ( !C.Interval[Interval_Manager::tg_sampled].initTrigger(m_stp, m_tm, m_dt, Interval_Manager::tg_sampled) ) {  // サンプリング履歴
      Hostonly_ printf("\t Error : Interval for Sampling output is asigned to zero.\n");
      Exit(0);
    }    
  }
  

  // 平均値のロード
  if ( C.Start == restart ) {
    TIMING_start(tm_restart);
    if ( C.Mode.Average == ON ) Restart_avrerage(fp, flop_task);
    TIMING_stop(tm_restart);
  }
  
  // データクラスのポイント
  if( !(v  = dc_v->GetData()) )     return -1;
  if( !(p  = dc_p->GetData()) )     return -1;
  
  if ( C.isHeatProblem() ) {
    if( !(t = dc_t->GetData()) )    return -1;
  }

  
	// リスタートの最大値と最小値の表示
  if ( (C.Start == restart) || (C.Start == coarse_restart) ) {
    
    Hostonly_ fprintf(mp, "\tNon-dimensional value\n");
    Hostonly_ fprintf(fp, "\tNon-dimensional value\n");
    REAL_TYPE f_min, f_max, min_tmp, max_tmp;
    
    // Velocity
    fb_minmax_v_ (&f_min, &f_max, sz, gc, v00, v, &flop_task); // allreduceすること
    
    if( para_mng->IsParallel() ){
      min_tmp = f_min;
      if( !para_mng->Allreduce(&min_tmp, &f_min, 1, SKL_ARRAY_DTYPE_REAL, SKL_MIN, pn.procGrp) ) Exit(0);
      
      max_tmp = f_max;
      if( !para_mng->Allreduce(&max_tmp, &f_max, 1, SKL_ARRAY_DTYPE_REAL, SKL_MAX, pn.procGrp) ) Exit(0);
    }
    
    Hostonly_ fprintf(mp, "\t\tV : min=%13.6e / max=%13.6e\n", f_min, f_max);
    Hostonly_ fprintf(fp, "\t\tV: min=%13.6e / max=%13.6e\n", f_min, f_max);
    
    
    // Pressure
    fb_minmax_s_ (&f_min, &f_max, sz, gc, p, &flop_task);
    
    if( para_mng->IsParallel() ){
      min_tmp = f_min;
      if( !para_mng->Allreduce(&min_tmp, &f_min, 1, SKL_ARRAY_DTYPE_REAL, SKL_MIN, pn.procGrp) ) Exit(0);
      
      max_tmp = f_max;
      if( !para_mng->Allreduce(&max_tmp, &f_max, 1, SKL_ARRAY_DTYPE_REAL, SKL_MAX, pn.procGrp) ) Exit(0);
    }
    
    Hostonly_ fprintf(mp, "\t\tP : min=%13.6e / max=%13.6e\n", f_min, f_max);
    Hostonly_ fprintf(fp, "\t\tP : min=%13.6e / max=%13.6e\n", f_min, f_max);
    
    // temperature
    if ( C.isHeatProblem() ) {
      fb_minmax_s_ (&f_min, &f_max, sz, gc, t, &flop_task);
      
      if( para_mng->IsParallel() ){
        min_tmp = f_min;
        if( !para_mng->Allreduce(&min_tmp, &f_min, 1, SKL_ARRAY_DTYPE_REAL, SKL_MIN, pn.procGrp) ) Exit(0);
        
        max_tmp = f_max;
        if( !para_mng->Allreduce(&max_tmp, &f_max, 1, SKL_ARRAY_DTYPE_REAL, SKL_MAX, pn.procGrp) ) Exit(0);
      }
      
      Hostonly_ fprintf(mp, "\t\tT : min=%13.6e / max=%13.6e\n", f_min, f_max);
      Hostonly_ fprintf(fp, "\t\tT : min=%13.6e / max=%13.6e\n", f_min, f_max);
    }
	}

  // 通信バッファの大きさの設定とアロケート
  if ( hasLinearSolver(SOR2SMA) ) {
    cf_sz[0] = (C.jmax+1) * (C.kmax+1) / 2;
    cf_sz[1] = (C.kmax+1) * (C.imax+1) / 2;
    cf_sz[2] = (C.imax+1) * (C.jmax+1) / 2;
    
    if( !( cf_x = SklAllocateSKL_REAL( cf_sz[0]*4, true ) ) ) return -1;
    if( !( cf_y = SklAllocateSKL_REAL( cf_sz[1]*4, true ) ) ) return -1;
    if( !( cf_z = SklAllocateSKL_REAL( cf_sz[2]*4, true ) ) ) return -1;
    TotalMemory += (unsigned long)( (cf_sz[0]*4+cf_sz[1]*4+cf_sz[2]*4)*sizeof(REAL_TYPE) );
  }

  
  // メモリ使用量の表示
  Hostonly_ {
    fprintf(mp,"\n---------------------------------------------------------------------------\n\n");
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
  }
  G_TotalMemory = TotalMemory;
  if( para_cmp->IsParallel() ) {
    tmp_memory = G_TotalMemory;
    para_cmp->Allreduce(&tmp_memory, &G_TotalMemory, 1, SKL_ARRAY_DTYPE_ULONG, SKL_SUM, pn.procGrp);
  }
  Hostonly_ {
    unsigned long mc = (unsigned long)( dc_wvex->GetArrayLength() * sizeof(REAL_TYPE) ); // temporaty array for vector output, see prepOutput();
    FBUtility::displayMemory("solver", G_TotalMemory, TotalMemory, fp, mp);
  }
  
  // 制御パラメータ，物理パラメータの表示
  Hostonly_ {
    C.displayParams(mp, fp, IC, &DT, &RF);
    Ex->printParaInfo(mp, fp, &C);
    
    // 外部境界面の開口率を表示
    C.printAreaInfo(mp, fp, G_Fcell, G_Acell, G_size);
    
    // 境界条件のリストと外部境界面のBC設定を表示
    Hostonly_ {
      fprintf(mp,"\n---------------------------------------------------------------------------\n\n");
      fprintf(mp,"\t>> Outer Boundary Conditions\n\n");
      fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
      fprintf(fp,"\t>> Outer Boundary Conditions\n\n");
      B.printOBCinfo(mp, fp, G_Lbx);
    }

    // モニタ情報の表示
    if ( C.Sampling.log == ON ) {
      
      MO.printMonitorInfo(mp, C.HistoryMonitorName, false); // ヘッダのみ
      
      FILE *fp_mon=NULL;
      Hostonly_ {
        if ( !(fp_mon=fopen("sampling_info.txt", "w")) ) {
          stamped_printf("\tSorry, can't open 'sampling_info.txt' file. Write failed.\n");
          return -1;
        }
      }
      
      MO.printMonitorInfo(fp_mon, C.HistoryMonitorName, true);  // 詳細モード
      Hostonly_ if ( fp_mon ) fclose(fp_mon);
    }
  }

  
  // ドライバ条件のチェック
  BC.checkDriver(fp);


  // 初期条件の条件設定
	if ( C.Start == initial_start ) {
		REAL_TYPE dt_init=1.0, tm_init=0.0;
		REAL_TYPE U0[3];
    
		// 速度の初期条件の設定
    if (C.Unit.Param == DIMENSIONAL) {
      U0[0] = C.iv.VecU/C.RefVelocity;
      U0[1] = C.iv.VecV/C.RefVelocity;
      U0[2] = C.iv.VecW/C.RefVelocity;
    }
    else {
      U0[0] = C.iv.VecU;
      U0[1] = C.iv.VecV;
      U0[2] = C.iv.VecW;
    }
		fb_set_vector_(v, sz, gc, U0, (int*)bcd);
    
		// 外部境界面の流出流量と移流速度
    DomainMonitor( BC.get_OBC_Ptr(), &C, flop_task);
    
		// 外部境界面の移流速度を計算し，外部境界条件を設定
    BC.OuterVBC_Periodic(dc_v);
		BC.OuterVBC(v, v, bcv, tm, (REAL_TYPE)m_dt, &C, v00, flop_task);
    BC.InnerVBC(v, bcv, tm, v00, flop_task);
    BC.InnerVBC_Periodic(dc_v, dc_bcd);
    
		// 圧力
    REAL_TYPE ip;
    int d_length;
    if (C.Unit.Param == DIMENSIONAL) {
      ip = FBUtility::convD2ND_P(C.iv.Pressure, C.BasePrs, C.RefDensity, C.RefVelocity, C.Unit.Prs);
    }
    else {
      ip = C.iv.Pressure;
    }

    d_length = (int)dc_p->GetArrayLength();
    fb_set_real_(dc_p->GetData(), &d_length, &ip);
		BC.OuterPBC(dc_p);
    
		// 温度
		if ( C.isHeatProblem() ) {
      REAL_TYPE it;
      if (C.Unit.Param == DIMENSIONAL) {
        it = FBUtility::convK2ND(C.iv.Temperature, C.BaseTemp, C.DiffTemp);
      }
      else {
        it = C.iv.Temperature;
      }

      d_length = (int)dc_t->GetArrayLength();
      fb_set_real_(dc_t->GetData(), &d_length, &it);
      
      // コンポーネントの初期値
      for (unsigned m=C.NoBC+1; m<=C.NoCompo; m++) {
        BC.setInitialTemp_Compo(m, bcd, dc_t);
      }
      
			BC.OuterTBC(dc_t);
		}
    
    // ユーザ例題のときに，速度の内部境界条件を設定する
    if ( C.Mode.Example == id_Users ) {
      ; // nothing
    }
    else {
      Ex->initCond(v, p);
    }
    
  }
  else { // リスタート時
    
    // 内部境界条件
    BC.InnerVBC(v, bcv, tm, v00, flop_task);
    BC.InnerVBC_Periodic(dc_v, dc_bcd);
    BC.InnerPBC_Periodic(dc_p, dc_bcd);
    
    // 外部境界条件
    BC.OuterVBC(v, v, bcv, tm, (REAL_TYPE)m_dt, &C, v00, flop_task);
    BC.OuterVBC_Periodic(dc_v);
    
    //流出境界の流出速度の算出
    REAL_TYPE coef = C.dh/(REAL_TYPE)m_dt;
    REAL_TYPE m_av[2];
    BC.mod_div(ws, bcv, coef, tm, v00, m_av, flop_task);
    DomainMonitor(BC.get_OBC_Ptr(), &C, flop_task);
    
    //if ( C.isHeatProblem() ) BC.InnerTBC_Periodic()
    
  }

  // VOF
  if ( C.BasicEqs == INCMP_2PHASE ) {
    REAL_TYPE* vof=NULL;
    if( !(vof = dc_vof->GetData()) )  return -1;
    setVOF(vof, bcd);
    if ( !dc_vof->CommBndCell(guide) ) return -1;
  }

  
  // 初期解およびリスタート解の同期
  if( !dc_v->CommBndCell(guide) ) return -1;
  if( !dc_p->CommBndCell(guide) ) return -1;
  if ( C.isHeatProblem() ) {
    if( !dc_t->CommBndCell(guide) ) return -1;
  }
  


  // マスターノードでの履歴出力準備
  H = new History(&C);
  
  Hostonly_ {
    H->printHistoryTitle(mp, IC, &C);
    
    // コンポーネント情報
    if ( C.Mode.Log_Base == ON ) {
      // 基本情報　history.log, history_compo.log, history_domfx.log
      if ( !(fp_b=fopen((char*)C.HistoryName, "w")) ) {
        stamped_printf("\tSorry, can't open '%s' file. Write failed.\n", (char*)C.HistoryName);
        return -1;
      }
      H->printHistoryTitle(fp_b, IC, &C);
      
      // コンポーネント履歴情報
      if ( !(fp_c=fopen((char*)C.HistoryCompoName, "w")) ) {
        stamped_printf("\tSorry, can't open '%s' file. Write failed.\n", (char*)C.HistoryCompoName);
        return -1;
      }
      H->printHistoryCompoTitle(fp_c, cmp, &C);

      // 流量収支情報　
      if ( !(fp_d=fopen((char*)C.HistoryDomfxName, "w")) ) {
        stamped_printf("\tSorry, can't open '%s' file. Write failed.\n", (char*)C.HistoryDomfxName);
        return -1;
      }
      H->printHistoryDomfxTitle(fp_d, &C);
      
      // 力の履歴情報　
      if ( !(fp_f=fopen((char*)C.HistoryForceName, "w")) ) {
        stamped_printf("\tSorry, can't open '%s' file. Write failed.\n", (char*)C.HistoryForceName);
        return -1;
      }
      H->printHistoryForceTitle(fp_f);
    }
    
    // 反復履歴情報　history_itr.log
    if ( C.Mode.Log_Itr == ON ) {
      if ( !(fp_i=fopen((char*)C.HistoryItrName, "w")) ) {
				stamped_printf("\tSorry, can't open '%s' file.\n", (char*)C.HistoryItrName);
        return -1;
      }
    }
    
    // 壁面情報　history_wall.log
    if ( C.Mode.Log_Wall == ON ) {
      if ( !(fp_w=fopen((char*)C.HistoryWallName, "w")) ) {
				stamped_printf("\tSorry, can't open '%s' file.\n", (char*)C.HistoryWallName);
        return -1;
      }
      H->printHistoryWallTitle(fp_w);
    }
  }

  // XMLによるサンプリング指定がある場合，モニタ結果出力ファイル群のオープン
  if ( C.Sampling.log == ON ) MO.openFile(C.HistoryMonitorName);
    
  // サンプリング元となるデータ配列の登録
  if ( C.Sampling.log == ON ) {
    if ( C.isHeatProblem() ) {
      MO.setDataPtrs(dc_v->GetData(), dc_p->GetData(), dc_t->GetData());
    }
    else {
      MO.setDataPtrs(dc_v->GetData(), dc_p->GetData());
    }
  }
  

  // 初期状態のファイル出力  リスタート時と性能測定モードのときには出力しない
	if ( (C.Hide.PM_Test == OFF) && (0 == SklGetTotalStep()) ) FileOutput(flop_task);
  
  
  // 粗い格子を用いたリスタート時には出力
  if ( C.Start == coarse_restart ) FileOutput(flop_task, true);
  
  
  // チェックモードの場合のコメント表示，前処理のみで中止---------------------------------------------------------
  if ( C.CheckParam == ON) {
		Hostonly_ fprintf(mp, "\n\tCheck mode --- Only pre-process\n\n");
		Hostonly_ fprintf(fp, "\n\tCheck mode --- Only pre-process\n\n");
    return 0;
	}

  // 組み込み例題の初期化
  Ex->PostInit(checkTime, &C);
  
  Hostonly_ if ( fp ) fclose(fp);
   
  TIMING_stop(tm_init_sct);

  return 1;
}

/**
 @fn void SklSolverCBC::allocArray_prep (unsigned long &prep, unsigned long &total)
 @brief 前処理に用いる配列のアロケーション
 @param total ソルバーに使用するメモリ量
 @param prep 前処理に使用するメモリ量
 */
void SklSolverCBC::allocArray_prep (unsigned long &prep, unsigned long &total)
{
  unsigned long mc=0;
  
  if ( !A.alloc_Real_S3D(this, dc_ws, "ws", size, guide, 0.0, mc) ) Exit(0);
  prep += mc;
  total+= mc;
  
  if ( !A.alloc_Int_S3D(this, dc_mid, "mid", size, guide, 0, mc) ) Exit(0);
  prep += mc;
  
  if ( !A.alloc_Uint_S3D(this, dc_bcd, "bcd", size, guide, 0, mc) ) Exit(0);
  prep += mc;
  total+= mc;
  
  if ( !A.alloc_Uint_S3D(this, dc_bcp, "bcp", size, guide, 0, mc) ) Exit(0);
  prep += mc;
  total+= mc;
  
  if ( !A.alloc_Uint_S3D(this, dc_bcv, "bcv", size, guide, 0, mc) ) Exit(0);
  prep += mc;
  total+= mc;
  
  if ( C.isHeatProblem() ) {
    if ( !A.alloc_Uint_S3D(this, dc_bh1, "bch1", size, guide, 0, mc) ) Exit(0);
    prep += mc;
    total+= mc;
    
    if ( !A.alloc_Uint_S3D(this, dc_bh2, "bch2", size, guide, 0, mc) ) Exit(0);
    prep += mc;
    total+= mc;
  }
}

/**
 @fn void SklSolverCBC::allocArray_compoVF (unsigned long &prep, unsigned long &total)
 @brief コンポーネント体積率の配列のアロケーション
 @param total ソルバーに使用するメモリ量
 @param prep 前処理に使用するメモリ量
 */
void SklSolverCBC::allocArray_compoVF (unsigned long &prep, unsigned long &total)
{
  unsigned long mc=0;
  
  if ( !A.alloc_Float_S3D(this, dc_cvf, "cvf", size, guide, 0.0, mc) ) Exit(0);
  prep += mc;
  total+= mc;
}

/**
 @fn void SklSolverCBC::allocArray_main (unsigned long &total)
 @brief 主計算部分に用いる配列のアロケーション
 @param total ソルバーに使用するメモリ量
 */
void SklSolverCBC::allocArray_main (unsigned long &total)
{
  unsigned long mc=0;
  
  if ( !A.alloc_Real_V3DEx(this, dc_v, "vel", size, guide, 0.0, mc) ) Exit(0);
  total+= mc;
  
  if ( !A.alloc_Real_V3DEx(this, dc_vc, "vc", size, guide, 0.0, mc) ) Exit(0);
  total+= mc;
  
  if ( !A.alloc_Real_V3DEx(this, dc_v0, "v0", size, guide, 0.0, mc) ) Exit(0);
  total+= mc;
  
  if ( !A.alloc_Real_V3DEx(this, dc_wv, "wv", size, guide, 0.0, mc) ) Exit(0);
  total+= mc;
  
  if ( !A.alloc_Real_V3DEx(this, dc_wvex, "wvex", size, guide, 0.0, mc) ) Exit(0);
  total+= mc;
  
  if ( !A.alloc_Real_S3D(this, dc_p, "prs", size, guide, 0.0, mc) ) Exit(0);
  total+= mc;
  
  if ( !A.alloc_Real_S3D(this, dc_p0, "prs0", size, guide, 0.0, mc) ) Exit(0);
  total+= mc;
  
  if ( !A.alloc_Real_S3D(this, dc_wk2, "wk2", size, guide, 0.0, mc) ) Exit(0);
  total+= mc;
}

/**
 @fn void SklSolverCBC::allocArray_interface (unsigned long &total)
 @brief allocation for interface equation
 @param total memory requirement in main solver
 */
void SklSolverCBC::allocArray_interface (unsigned long &total)
{
  unsigned long mc=0;
  
  if ( !A.alloc_Real_S3D(this, dc_vof, "vof", size, guide, 0.0, mc) ) Exit(0);
  total+= mc;
}

/**
 @fn void SklSolverCBC::allocArray_heat (unsigned long &total)
 @brief 熱の主計算部分に用いる配列のアロケーション
 @param total ソルバーに使用するメモリ量
 */
void SklSolverCBC::allocArray_heat (unsigned long &total)
{
  unsigned long mc=0;
  
  if ( !A.alloc_Real_S3D(this, dc_t, "tmp", size, guide, 0.0, mc) ) Exit(0);
  total+= mc;
  
  if ( !A.alloc_Real_S3D(this, dc_t0, "tmp0", size, guide, 0.0, mc) ) Exit(0);
  total+= mc;
  
  if ( !A.alloc_Real_V3DEx(this, dc_qbc, "qbc", size, guide, 0.0, mc) ) Exit(0);
  total+= mc;
}

/**
 @fn void SklSolverCBC::allocArray_AB2 (unsigned long &total)
 @brief Adams-Bashforth法に用いる配列のアロケーション
 @param total ソルバーに使用するメモリ量
 */
void SklSolverCBC::allocArray_AB2 (unsigned long &total)
{
  unsigned long mc=0;
  if ( !A.alloc_Real_V3DEx(this, dc_abf, "abf", size, guide, 0.0, mc) ) Exit(0);
  total+= mc;
}

/**
 @fn void SklSolverCBC::allocArray_RK (unsigned long &total)
 @brief Runge-Kutta法に用いる配列のアロケーション
 @param total ソルバーに使用するメモリ量
 */
void SklSolverCBC::allocArray_RK (unsigned long &total)
{
  unsigned long mc=0;
  
  if ( !A.alloc_Real_S3D(this, dc_dp, "dp", size, guide, 0.0, mc) ) Exit(0);
  total+= mc;
}


/**
 @fn void SklSolverCBC::allocArray_CoarseMesh(unsigned* r_size)
 @brief 前処理に用いる配列のアロケーション
 @param r_size 粗格子の領域サイズ

 */
void SklSolverCBC::allocArray_CoarseMesh(unsigned* r_size)
{
  unsigned long mc=0;

  if ( !A.alloc_Real_S3D(this, dc_r_p, "r_prs", r_size, guide, 0.0, mc) ) Exit(0);
  //prep += mc;
  
  if ( !A.alloc_Real_V3DEx(this, dc_r_v, "r_vel", r_size, guide, 0.0, mc) ) Exit(0);
  //prep += mc;
  
  
  if ( C.isHeatProblem() ) {
    if ( !A.alloc_Real_S3D(this, dc_r_t, "r_tmp", r_size, guide, 0.0, mc) ) Exit(0);
    //prep += mc;
  }
}


/**
 @fn void SklSolverCBC::allocArray_average (unsigned long &total, FILE* fp)
 @brief 平均値処理に用いる配列のアロケーション
 @param total ソルバーに使用するメモリ量
 @param fp ファイルポインタ
 @note
 - 配列長は計算内部領域のみなので，関連する処理に注意
 */
void SklSolverCBC::allocArray_average (unsigned long &total, FILE* fp)
{
  unsigned long mc=0;
  
  if ( C.Mode.Average == ON ) {
    if ( !SklCfgCheckOutFile("AvrPressure") || // 時間平均を指定しているが，出力ファイル記述がない場合
         !SklCfgCheckOutFile("AvrVelocity") || 
        (!SklCfgCheckOutFile("AvrTemperature") && C.isHeatProblem()) ) {
      Hostonly_ stamped_printf     ("\tAveraging, but there is no OutFile description for an average file. \n");
      Hostonly_ stamped_fprintf(fp, "\tAveraging, but there is no OutFile description for an average file. \n");
      Exit(0);
    }
    
    if ( !A.alloc_Real_S3D(this, dc_ap, "avtp", size, guide, 0.0, mc) ) Exit(0);
    total += mc;
    
    if ( !A.alloc_Real_V3DEx(this, dc_av, "avrv", size, guide, 0.0, mc) ) Exit(0);
    total += mc;
    
    if ( C.isHeatProblem() ) {
      if ( !A.alloc_Real_S3D(this, dc_at, "avrt", size, guide, 0.0, mc) ) Exit(0);
      total += mc;
    }
  }
}

/**
 @fn void SklSolverCBC::allocArray_LES (unsigned long &total)
 @brief LES計算に用いる配列のアロケーション
 @param total ソルバーに使用するメモリ量
 */
void SklSolverCBC::allocArray_LES (unsigned long &total)
{
  unsigned long mc=0;
  
  if ( !A.alloc_Real_S3D(this, dc_vt, "eddy", size, guide, 0.0, mc) ) Exit(0);
  total+= mc;
}

/**
 @fn void SklSolverCBC::allocArray_Collocate(unsigned long &total)
 @brief コロケート格子用のセルフェイスの配列
 @param total ソルバーに使用するメモリ量
 */
void SklSolverCBC::allocArray_Collocate(unsigned long &total)
{
  unsigned long mc=0;
  
  if ( !A.alloc_Real_V3DEx(this, dc_vf0, "vf0", size, guide, 0.0, mc) ) Exit(0);
  total+= mc;
}

/**
 @fn void SklSolverCBC::allocArray_Cut(unsigned long &total)
 @brief コロケート格子用のセルフェイスの配列
 @param total ソルバーに使用するメモリ量
 */
void SklSolverCBC::allocArray_Cut(unsigned long &total)
{
  unsigned long mc=0;

  if ( !A.alloc_Float_S4DEx(this, dc_cut, "cut", size, guide, 6, 0.0, mc) ) Exit(0);
  total+= mc;

  if ( !A.alloc_Int_S3D(this, dc_bid, "bid", size, guide, 0, mc) ) Exit(0);
  total += mc;
}

/**
 @fn void SklSolverCBC::allocComponentArray(unsigned long& m_prep, unsigned long& m_total, FILE* fp)
 @brief コンポーネントのワーク用配列のアロケート
 @param cutPos カット情報コンテナ
 @param m_prep  前処理用のメモリサイズ
 @param m_total 本計算用のメモリリサイズ
 @param fp
 @note
 - bcdに対して，共通の処理を行い，それをbcp, bcv, bchにコピーし，その後個別処理を行う．
 */
void SklSolverCBC::allocComponentArray(unsigned long& m_prep, unsigned long& m_total, FILE* fp)
{
  SklParaComponent* para_cmp = SklGetParaComponent();
  
  // 管理用のポインタ配列の確保
  component_array = new REAL_TYPE* [C.NoBC];
  
  for (int i=0; i<C.NoBC; i++) {
    component_array[i] = NULL;
  }
  
  
  // リサイズ後のインデクスサイズの登録と配列領域の確保
  int c_sz[3];
  int gd=2; // 両側それぞれ2セル
  unsigned long m_cmp_size=0;
  
  for (int n=1; n<=C.NoBC; n++) {

    if ( cmp[n].isFORCING() ) {

      // インデクスサイズをリサイズしたst[], ed[]から計算
      cmp[n].set_cmp_sz();
      cmp[n].get_cmp_sz(c_sz);
      
      // ワーク用の配列を確保
      size_t array_size = (c_sz[0]+2*gd) * (c_sz[1]+2*gd) * (c_sz[2]+2*gd) * 3;
      component_array[n] = new REAL_TYPE[array_size];
      m_cmp_size += array_size;
    }
  }
  
  // 使用メモリ量　
  unsigned long cmp_mem, G_cmp_mem;
  G_cmp_mem = cmp_mem = m_cmp_size * sizeof(REAL_TYPE);
  m_prep += cmp_mem;
  m_total+= cmp_mem;
  
  if( para_cmp->IsParallel() ) {
    para_cmp->Allreduce(&cmp_mem, &G_cmp_mem, 1, SKL_ARRAY_DTYPE_ULONG, SKL_SUM, pn.procGrp);
  }
  
  Hostonly_  {
    FBUtility::displayMemory("component", G_cmp_mem, cmp_mem, fp, mp);
  }
}

/**
 @fn bool SklSolverCBC::chkMediumConsistency(void)
 @brief 全Voxelモデルの媒質数とKOSの整合性をチェック
 @retval エラーコード
 @note NoMediumSolidなどは，ParseBC::setMedium()で取得
 */
bool SklSolverCBC::chkMediumConsistency(void)
{
  SklParaManager* para_mng = ParaCmpo->GetParaManager();
  unsigned long nmSolid = C.NoMediumSolid;
  unsigned long nmFluid = C.NoMediumFluid;
  
  if( para_mng->IsParallel() ){
    unsigned long nms = nmSolid;
    unsigned long nmf = nmFluid;
    if( !para_mng->Allreduce(&nms, &nmSolid, 1, SKL_ARRAY_DTYPE_ULONG, SKL_SUM, pn.procGrp) ) return false;
    if( !para_mng->Allreduce(&nmf, &nmFluid, 1, SKL_ARRAY_DTYPE_ULONG, SKL_SUM, pn.procGrp) ) return false;
  }
  
  if ( (nmFluid == 0) && (nmSolid == 0) ) {
    Hostonly_ printf("\tError : No medium\n");
    return false;
  }
  
  switch (C.KindOfSolver) {
    case FLOW_ONLY:
    case THERMAL_FLOW:
    case THERMAL_FLOW_NATURAL:
      
      if ( nmFluid == 0 ) {
        Hostonly_ printf("\tError : No FLUID medium\n");
        return false;
      }
      break;
      
    case CONJUGATE_HEAT_TRANSFER:
      if ( ( nmFluid == 0 ) || ( nmSolid == 0 ) ) {
        Hostonly_ printf("\tError : Fluid/Solid should have at least one medium.\n");
        return false;
      }
      break;
      
    case SOLID_CONDUCTION:
      if ( nmSolid == 0 ) {
        Hostonly_ printf("\tError : No Solid medium\n");
        return false;
      }
      break;
  };
  
  return true;
}

/**
 @fn void SklSolverCBC::connectExample(Control* Cref)
 @brief 組み込み例題のインスタンス
 @param Cref Controlクラスのポインタ
 */
void SklSolverCBC::connectExample(Control* Cref)
{
  if      ( Cref->Mode.Example == id_Users )   Ex = dynamic_cast<Intrinsic*>(new IP_Users);
  else if ( Cref->Mode.Example == id_PPLT2D)   Ex = dynamic_cast<Intrinsic*>(new IP_PPLT2D);
  else if ( Cref->Mode.Example == id_SHC1D)    Ex = dynamic_cast<Intrinsic*>(new IP_SHC1D);
  else if ( Cref->Mode.Example == id_Duct )    Ex = dynamic_cast<Intrinsic*>(new IP_Duct);
  else if ( Cref->Mode.Example == id_PMT )     Ex = dynamic_cast<Intrinsic*>(new IP_PMT);
  else if ( Cref->Mode.Example == id_Rect )    Ex = dynamic_cast<Intrinsic*>(new IP_Rect);
  else if ( Cref->Mode.Example == id_Cylinder) Ex = dynamic_cast<Intrinsic*>(new IP_Cylinder);
  else if ( Cref->Mode.Example == id_Step )    Ex = dynamic_cast<Intrinsic*>(new IP_Step);
  else if ( Cref->Mode.Example == id_Polygon ) Ex = dynamic_cast<Intrinsic*>(new IP_Polygon);
  else if ( Cref->Mode.Example == id_Sphere )  Ex = dynamic_cast<Intrinsic*>(new IP_Sphere);
  else {
    Hostonly_ stamped_printf("\tInvalid keyword is described for Exmple definition\n");
    Exit(0);
  }
}

/**
 @fn void SklSolverCBC::EnlargeIndex(int& m_st, int& m_ed, const unsigned st_i, const unsigned len, const unsigned mx, const unsigned dir, const int m_id)
 @brief 初期インデクスの情報を元に，一層拡大したインデクス値を返す
 @param m_st 拡大された開始点（Fortranインデクス）
 @param m_ed 拡大された終了点（Fortranインデクス）
 @param st_i 開始点（Cインデクス）
 @param len コンポーネントの存在長さ
 @param m_x 軸方向のサイズ
 @param dir 方向
 @param m_id キーID
 */
void SklSolverCBC::EnlargeIndex(int& m_st, int& m_ed, const unsigned st_i, const unsigned len, const unsigned m_x, const unsigned dir, const int m_id)
{
  unsigned ed_i = st_i + len - 1;
  unsigned n_st = st_i - 1;
  unsigned n_ed = ed_i + 1;
  unsigned max_c1 = m_x + guide;
  
  unsigned label_st, label_ed;
  
  switch (dir) {
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
  if ( ed_i < guide ) { 
    if( pn.nID[label_st] < 0 ){ // 計算領域の外部面に接する場合は，対象外
      m_st = 0;
      m_ed = 0;
    }
    else { // 計算領域内部にある場合（並列時）
      if ( n_ed == guide ) { // ガイドセル1層外側の場合
        m_st = 1; // F index
        m_ed = 1; // F index
      }
      else {
        m_st = 0;
        m_ed = 0;
      }
    }
  }
  
  // BVが+方向のガイドセル内のみにある場合
  else if ( st_i >= max_c1 ) {
    if( pn.nID[label_ed] < 0 ){ // 計算領域の外部面に接する場合は，対象外
      m_st = 0;
      m_ed = 0;
    }
    else {
      if ( n_st == (max_c1 - 1) ) { // ガイドセル1層外側の場合
        m_st = m_x; // F index
        m_ed = m_x; // F index
      }
      else {
        m_st = 0;
        m_ed = 0;
      }
    }
    //debug; Hostonly_ printf("(2) dir=%d : st=%d ed=%d\n", dir, m_st, m_ed);
  }
  
  // BVが内部領域のみにある場合（逐次・並列で同じ処理）
  else if ( (st_i >= guide) && (ed_i < max_c1) ) {
    if ( st_i == guide ) { // 最外層セル
      m_st = 1; // F index
    }
    else {
      m_st = n_st + 1 - guide; // F index
    }
    
    if ( ed_i == (max_c1 - 1) ) { // 最外層セル
      m_ed = m_x; // F index
    }
    else { // 内部
      m_ed = n_ed + 1 - guide; // F index
    }
    //debug; Hostonly_ printf("(3) dir=%d : st=%d ed=%d\n", dir, m_st, m_ed);
  }
  
  // BVが両方向のガイドセルにまたがる場合（逐次・並列で同じ処理）
  else if ( (st_i < guide) && (ed_i >= max_c1) ) {
    m_st = 1; // F index
    m_ed = m_x; // F index
    //debug; Hostonly_ printf("(4) dir=%d : st=%d ed=%d\n", dir, m_st, m_ed);
  }
  
  // BVが-方向のガイドセルから内部領域にある場合
  else if ( (st_i < guide) && (ed_i < max_c1) ) {
    m_st = 1; // F index
    
    if ( ed_i == (max_c1 - 1) ) { // 最外層セル
      m_ed = m_x; // F index
    }
    else { // 内部
      m_ed = n_ed + 1 - guide; // F index
    }
    //debug; Hostonly_ printf("(5) dir=%d : st=%d ed=%d\n", dir, m_st, m_ed);
  }
  
  // BVが+方向のガイドセルから内部領域にある場合
  else if ( (st_i < max_c1) && (ed_i >= max_c1) ) {
    m_ed = m_x; // F index
    
    if ( st_i == guide ) { // 端点
      m_st = 1; // F index
    }
    else { // 内部
      m_st = n_st + 1 - guide; // F index
    }
    //debug; Hostonly_ printf("(6) dir=%d : st=%d ed=%d\n", dir, m_st, m_ed);
  }
  
  else {
    string m_dir;
    if      ( dir == 0 ) m_dir = "X";
    else if ( dir == 1 ) m_dir = "Y";
    else                 m_dir = "Z";
    
    Hostonly_ stamped_printf("\tError : Unexpected case for ID=%d, (%d - %d): %s\n", m_id, st_i, ed_i, m_dir.c_str());
    Exit(0);
  }
  
}

/**
 @fn void SklSolverCBC::fixed_parameters(void)
 @brief 固定パラメータの設定
 */
void SklSolverCBC::fixed_parameters(void)
{
  // 次元
  C.NoDimension = 3;
  
  // 精度
  if ( sizeof(REAL_TYPE) == sizeof(double) ) {
    C.Mode.Precision = SPH_DOUBLE;
  }
  else {
    C.Mode.Precision = SPH_SINGLE;
  }
  
  // ログファイル名
  strcpy(C.HistoryName,        "history_base.txt");
  strcpy(C.HistoryCompoName,   "history_compo.txt");
  strcpy(C.HistoryDomfxName,   "history_domainflux.txt");
  strcpy(C.HistoryForceName,   "history_force.txt");
  strcpy(C.HistoryWallName,    "history_log_wall.txt");
  strcpy(C.HistoryItrName,     "history_iteration.txt");
  strcpy(C.HistoryMonitorName, "sample.log");
  
  C.f_Pressure       = "prs_";
  C.f_Velocity       = "vel_";
  C.f_Temperature    = "tmp_";
  C.f_AvrPressure    = "prsa_";
  C.f_AvrVelocity    = "vela_";
  C.f_AvrTemperature = "tmpa_";
  C.f_DivDebug       = "div_";
  C.f_Helicity       = "hlt_";
  C.f_TotalP         = "tp_";
  C.f_I2VGT          = "i2vgt_";
  C.f_Vorticity      = "vrt_";
  
}

/**
 @fn void SklSolverCBC::gather_DomainInfo(void)
 @brief 並列処理時の各ノードの分割数を集めてファイルに保存する
 */
void SklSolverCBC::gather_DomainInfo(void)
{
  SklParaComponent* para_cmp = SklGetParaComponent();
  SklParaManager* para_mng = para_cmp->GetParaManager();
  
  int nID[6], np=0;
  FILE *fp=NULL;
  REAL_TYPE vol, srf, m_vol, m_srf, vol_dv, srf_dv, m_efv, efv_dv;
  REAL_TYPE d1, d2, d3, d, r;
  int ix, jx, kx;
  
  np=para_mng->GetNodeNum(pn.procGrp);
  d = 1.0/(REAL_TYPE)np;
  if ( para_mng->IsParallel() ) {
    r = 1.0/(REAL_TYPE)(np-1);
  }
  else {
    r = 1.0;
  }
  
  unsigned* m_size=NULL; if( !(m_size = new unsigned[np*3]) ) Exit(0); // use new to assign variable array, and release at the end of this method
  REAL_TYPE* m_org=NULL;  if( !(m_org  = new REAL_TYPE[np*3]) ) Exit(0);
  REAL_TYPE* m_Lbx=NULL;  if( !(m_Lbx  = new REAL_TYPE[np*3]) ) Exit(0);
  unsigned* st_buf=NULL; if( !(st_buf = new unsigned[np*3]) ) Exit(0);
  unsigned* ed_buf=NULL; if( !(ed_buf = new unsigned[np*3]) ) Exit(0);
  REAL_TYPE *bf_srf=NULL; if( !(bf_srf = new REAL_TYPE[np]) )   Exit(0);
  unsigned* bf_fcl=NULL; if( !(bf_fcl = new unsigned[np]) )   Exit(0);
  unsigned* bf_wcl=NULL; if( !(bf_wcl = new unsigned[np]) )   Exit(0);
  unsigned* bf_acl=NULL; if( !(bf_acl = new unsigned[np]) )   Exit(0);
  
  // 領域情報の収集
  if ( para_mng->IsParallel() ) {
    if ( !para_mng->Gather(size,    3, SKL_ARRAY_DTYPE_UINT, 
                           m_size,  3, SKL_ARRAY_DTYPE_UINT, 0, pn.procGrp) ) Exit(0);
    if ( !para_mng->Gather(C.org,   3, SKL_ARRAY_DTYPE_REAL, 
                           m_org,   3, SKL_ARRAY_DTYPE_REAL, 0, pn.procGrp) ) Exit(0);
    if ( !para_mng->Gather(C.Lbx,   3, SKL_ARRAY_DTYPE_REAL, 
                           m_Lbx,   3, SKL_ARRAY_DTYPE_REAL, 0, pn.procGrp) ) Exit(0);
    if ( !para_mng->Gather(&C.Fcell,1, SKL_ARRAY_DTYPE_UINT, 
                           bf_fcl,  1, SKL_ARRAY_DTYPE_UINT, 0, pn.procGrp) ) Exit(0);
    if ( !para_mng->Gather(&C.Wcell,1, SKL_ARRAY_DTYPE_UINT, 
                           bf_wcl,  1, SKL_ARRAY_DTYPE_UINT, 0, pn.procGrp) ) Exit(0);
    if ( !para_mng->Gather(&C.Acell,1, SKL_ARRAY_DTYPE_UINT, 
                           bf_acl,  1, SKL_ARRAY_DTYPE_UINT, 0, pn.procGrp) ) Exit(0);
  }
  else { // serial
    memcpy(m_size, size, 3*sizeof(unsigned));
    bf_fcl[0] = C.Fcell;
    bf_wcl[0] = C.Wcell;
    bf_acl[0] = C.Acell;
    memcpy(m_org, C.org, 3*sizeof(REAL_TYPE));
    memcpy(m_Lbx, C.Lbx, 3*sizeof(REAL_TYPE));
  }
  
  // Info. of computational domain
  vol = (REAL_TYPE)(G_size[0]*G_size[1]*G_size[2]);
  srf = (REAL_TYPE)(2*(G_size[0]*G_size[1] + G_size[1]*G_size[2] + G_size[2]*G_size[0]));
  
  // amount of communication in each node
  ix = size[0];
  jx = size[1];
  kx = size[2];
  m_srf = (REAL_TYPE)(2*(ix*jx + jx*kx + kx*ix));
  if ( pn.nID[X_MINUS] < 0 ) m_srf -= (REAL_TYPE)(jx*kx);  // remove face which does not join communication
  if ( pn.nID[Y_MINUS] < 0 ) m_srf -= (REAL_TYPE)(ix*kx);
  if ( pn.nID[Z_MINUS] < 0 ) m_srf -= (REAL_TYPE)(ix*jx);
  if ( pn.nID[X_PLUS]  < 0 ) m_srf -= (REAL_TYPE)(jx*kx);
  if ( pn.nID[Y_PLUS]  < 0 ) m_srf -= (REAL_TYPE)(ix*kx);
  if ( pn.nID[Z_PLUS]  < 0 ) m_srf -= (REAL_TYPE)(ix*jx);
  
  if ( para_mng->IsParallel() ) {
    if ( !para_mng->Gather(&m_srf,  1,        SKL_ARRAY_DTYPE_REAL, 
                           bf_srf,  1,        SKL_ARRAY_DTYPE_REAL, 0, pn.procGrp) ) Exit(0);
  }
  else {
    bf_srf[0] = m_srf;
  }
  
  // mean of domain
  m_vol = m_srf = m_efv = 0.0;
  for (int i=0; i<np; i++) {
    ix = m_size[3*i];
    jx = m_size[3*i+1];
    kx = m_size[3*i+2];
    m_vol += (REAL_TYPE)(ix*jx*kx);
    m_srf += bf_srf[i];
    m_efv += bf_acl[i];
  }
  m_vol *= d;
  m_srf *= d;
  m_efv *= d;
  
  // std. deviation of domain
  vol_dv = srf_dv = efv_dv = 0.0;
  for (int i=0; i<np; i++) {
    ix = m_size[3*i];
    jx = m_size[3*i+1];
    kx = m_size[3*i+2];
    d1 = (REAL_TYPE)(ix*jx*kx) - m_vol;
    d2 = bf_srf[i] - m_srf;
    d3 = (REAL_TYPE)bf_acl[i] - m_efv;
    vol_dv += d1*d1;
    srf_dv += d2*d2;
    efv_dv += d3*d3;
  }
  vol_dv = sqrt(vol_dv*r);
  srf_dv = sqrt(srf_dv*r);
  efv_dv = sqrt(efv_dv*r);
  
  if ( !(fp=fopen("DomainInfo.txt", "w")) ) {
    stamped_printf("\tSorry, can't open 'DomainInfo.txt' file. Write failed.\n");
    Exit(0);
  }
  
  // 全体情報の表示
  C.printDomain(fp, G_size, G_org, G_Lbx);
  
  // ローカルノードの情報を表示
  for (int i=0; i<np; i++) {
    Hostonly_ {
      fprintf(fp,"Domain %4d\n", i);
      fprintf(fp,"\t ix, jx,  kx        [-] =  %13d %13d %13d\n",  m_size[i*3], m_size[i*3+1], m_size[i*3+2]);
      fprintf(fp,"\t(ox, oy, oz)  [m] / [-] = (%13.6e %13.6e %13.6e)  /  (%13.6e %13.6e %13.6e)\n", 
              m_org[i*3]*C.RefLength,  m_org[i*3+1]*C.RefLength,  m_org[i*3+2]*C.RefLength, m_org[i*3],  m_org[i*3+1],  m_org[i*3+2]);
      fprintf(fp,"\t(Lx, Ly, Lz)  [m] / [-] = (%13.6e %13.6e %13.6e)  /  (%13.6e %13.6e %13.6e)\n", 
              m_Lbx[i*3]*C.RefLength,  m_Lbx[i*3+1]*C.RefLength,  m_Lbx[i*3+2]*C.RefLength, m_Lbx[i*3],  m_Lbx[i*3+1],  m_Lbx[i*3+2]);
      
      if (C.NoBC != 0) fprintf(fp, "\t no            Label    ID    i_st    i_ed    j_st    j_ed    k_st    k_ed\n");
    }
    
    if( para_mng->IsParallel() ) {
      for (unsigned n=1; n<=C.NoBC; n++) {
        if( !para_mng->Gather(cmp[n].getBbox_st(), 3, SKL_ARRAY_DTYPE_UINT, st_buf, 3, SKL_ARRAY_DTYPE_UINT, 0, pn.procGrp) ) Exit(0);
        if( !para_mng->Gather(cmp[n].getBbox_ed(), 3, SKL_ARRAY_DTYPE_UINT, ed_buf, 3, SKL_ARRAY_DTYPE_UINT, 0, pn.procGrp) ) Exit(0);
        
        Hostonly_ {
          fprintf(fp,"\t%3d %16s %5d %7d %7d %7d %7d %7d %7d\n",
                  n, cmp[n].name, cmp[n].getID(), st_buf[i*3], ed_buf[i*3], st_buf[i*3+1], ed_buf[i*3+1], st_buf[i*3+2], ed_buf[i*3+2]);
        }
      }
    }
  }
  
  Hostonly_ {
    fprintf(fp, "\n");
    fprintf(fp,"\n\t--------------------------------------------------\n");
    fprintf(fp,"\tReport of Whole Domain Statistics\n");
    fprintf(fp,"\tDomain size         = %7d %7d %7d\n", G_size[0], G_size[1], G_size[2]);
    fprintf(fp,"\tNumber of voxels    = %12.6e\n", vol);
    fprintf(fp,"\tNumber of surface   = %12.6e\n", srf);
    fprintf(fp,"\tEffective voxels    = %12.6e (%6.2f%%)\n", (REAL_TYPE)G_Acell, 100.0*(REAL_TYPE)G_Acell/vol);
    fprintf(fp,"\tFluid voxels        = %12.6e (%6.2f%%)\n", (REAL_TYPE)G_Fcell, 100.0*(REAL_TYPE)G_Fcell/vol);
    fprintf(fp,"\tWall  voxels        = %12.6e (%6.2f%%)\n", (REAL_TYPE)G_Wcell, 100.0*(REAL_TYPE)G_Wcell/vol);
    if ( np == 1 ) {
      fprintf(fp,"\tDivision :          = %d : %s\n", np, "Serial");
    }
    else {
      fprintf(fp,"\tDivision :          = %d : %s\n", np, (para_mng->IsEv()) ? "Equal segregation" : "Multi-Box division");
    }
    fprintf(fp,"\n\t--------------------------------------------------\n");
    fprintf(fp,"\tDomain Statistics per MPI process\n");
    fprintf(fp,"\tMean volume in each domain           = %12.6e\n", m_vol);
    fprintf(fp,"\tStd. deviation of domain             = %12.6e\n", vol_dv);
    fprintf(fp,"\tMean comm. in each domain            = %12.6e\n", m_srf);
    fprintf(fp,"\tStd. deviation of surface            = %12.6e\n", srf_dv);
    fprintf(fp,"\tMean effective volume in each domain = %12.6e\n", m_efv);
    fprintf(fp,"\tStd. deviation of effective volume   = %12.6e\n", efv_dv);
    fprintf(fp,"\n");
    
    fprintf(fp,"\tDomain :     ix     jx     kx       Volume Vol_dv[%%]      Surface Srf_dv[%%] Fluid[%%] Solid[%%]      Eff_Vol Eff_Vol_dv[%%]      Eff_Srf Eff_srf_dv[%%]  Itr_scheme\n");
    fprintf(fp,"\t---------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
    
    REAL_TYPE tmp_vol, tmp_acl, tmp_fcl, tmp_wcl;
    for (int i=0; i<np; i++) {
      ix = m_size[3*i];
      jx = m_size[3*i+1];
      kx = m_size[3*i+2];
      tmp_vol = (REAL_TYPE)(ix*jx*kx);
      tmp_acl = (REAL_TYPE)bf_acl[i];
      tmp_fcl = (REAL_TYPE)bf_fcl[i];
      tmp_wcl = (REAL_TYPE)bf_wcl[i];
      fprintf(fp,"\t%6d : %6d %6d %6d ", i, ix, jx, kx);
      fprintf(fp,"%12.4e  %8.3f ", tmp_vol, 100.0*(tmp_vol-m_vol)/m_vol);
      fprintf(fp,"%12.4e  %8.3f ", bf_srf[i], (m_srf == 0.0) ? 0.0 : 100.0*(bf_srf[i]-m_srf)/m_srf);
      fprintf(fp,"%8.3f %8.3f ", 100.0*tmp_fcl/tmp_vol, 100.0*tmp_wcl/tmp_vol);
      fprintf(fp,"%12.4e      %8.3f \n", tmp_acl, 100.0*(tmp_acl-m_efv)/m_efv);
    }
    fprintf(fp,"\t---------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
  }
  
  if (fp) fclose(fp);
  
  if( m_size ) { delete [] m_size; m_size=NULL; }
  if( m_org  ) { delete [] m_org;  m_org =NULL; }
  if( m_Lbx  ) { delete [] m_Lbx;  m_Lbx =NULL; }
  if( st_buf ) { delete [] st_buf; st_buf=NULL; }
  if( ed_buf ) { delete [] ed_buf; ed_buf=NULL; }
  if( bf_srf ) { delete [] bf_srf; bf_srf=NULL; }
  if( bf_fcl ) { delete [] bf_fcl; bf_fcl=NULL; }
  if( bf_wcl ) { delete [] bf_wcl; bf_wcl=NULL; }
  if( bf_acl ) { delete [] bf_acl; bf_acl=NULL; }
}


/**
 @fn void SklSolverCBC::getXMLExample(Control* Cref)
 @brief 組み込み例題の設定
 @param Cref Controlクラスのポインタ
 */
void SklSolverCBC::getXMLExample(Control* Cref)
{
  const char *keyword=NULL;
  ParseSteer Tree(m_solvCfg);
  
  if ( !(keyword=Tree.getParam("Example")) ) Exit(0);
  
  if     ( !strcasecmp(keyword, "Users") )                    Cref->Mode.Example = id_Users;
  else if( !strcasecmp(keyword, "Parallel_Plate_2D") )        Cref->Mode.Example = id_PPLT2D;
  else if( !strcasecmp(keyword, "Duct") )                     Cref->Mode.Example = id_Duct;
  else if( !strcasecmp(keyword, "SHC1D"))                     Cref->Mode.Example = id_SHC1D;
  else if( !strcasecmp(keyword, "Performance_Test") )         Cref->Mode.Example = id_PMT;
  else if( !strcasecmp(keyword, "Rectangular") )              Cref->Mode.Example = id_Rect;
  else if( !strcasecmp(keyword, "Cylinder") )                 Cref->Mode.Example = id_Cylinder;
  else if( !strcasecmp(keyword, "Back_Step") )                Cref->Mode.Example = id_Step;
  else if( !strcasecmp(keyword, "Polygon") )                  Cref->Mode.Example = id_Polygon;
  else if( !strcasecmp(keyword, "Sphere") )                   Cref->Mode.Example = id_Sphere;
  else {
    Hostonly_ stamped_printf("\tInvalid keyword is described for Example definition\n");
    Exit(0);
  }
}

/**
 @fn void SklSolverCBC::getXML_Monitor(SklSolverConfig* CF, MonitorList* M)
 @brief XMLに記述されたモニタ座標情報を取得し，リストに保持する
 @param M MonitorList クラスオブジェクトのポインタ
 */
void SklSolverCBC::getXML_Monitor(SklSolverConfig* CF, MonitorList* M)
{
  const CfgElem *elemTop=NULL, *elmL1=NULL, *elmL2=NULL;
  const char* p=NULL;
  const char* str=NULL;
  const char* label=NULL;
  const char* method=NULL;
  const char* mode=NULL;
  const CfgParam * param=NULL;
  REAL_TYPE f_val=0.0;
  //int nvar=0;
  MonitorCompo::Type type;
  vector<string> variables;
  
  // Monitor_List section is already confirmed
  elemTop = CF->GetTop(STEER);
  elmL1 = elemTop->GetElemFirst("Monitor_List");
  
  // ログ出力のON/OFFはControl::getXML_Sampling()で取得済み
  
  // 集約モード
  if ( !elmL1->GetValue(CfgIdt("output_mode"), &str) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Output_Mode' in 'Monitor_List'\n");
    Exit(0);
  }
  if ( !strcasecmp("gather", str)) {
    C.Sampling.out_mode = MonitorList::GATHER;
    M->setOutputType(MonitorList::GATHER);
  }
  else if ( !strcasecmp("distribute", str)) {
    C.Sampling.out_mode = MonitorList::DISTRIBUTE;
    M->setOutputType(MonitorList::DISTRIBUTE);
  }
  else {
    Hostonly_ stamped_printf("\tParsing error : Invalid keyord for 'Output_Mode'\n");
    Exit(0);
  }
  
  // サンプリング間隔
  if ( !elmL1->GetValue("Sampling_Interval_Type", &str) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'Sampling_Interval_Type' in 'Monitor_List'\n");
    Exit(0);
  }
  else {
    if ( !strcasecmp(str, "step") ) {
      C.Interval[Interval_Manager::tg_sampled].setMode_Step();
    }
    else if ( !strcasecmp(str, "time") ) {
      C.Interval[Interval_Manager::tg_sampled].setMode_Time();
    }
    else {
      Hostonly_ stamped_printf("\tParsing error : Invalid keyword for 'Sampling_Interval_Type' in 'Monitor_List'\n");
      Exit(0);
    }
    
    if ( elmL1->GetValue("Sampling_Interval", &f_val) ) {
      C.Interval[Interval_Manager::tg_sampled].setInterval((double)f_val);
    }
    else {
      Hostonly_ stamped_printf("\tParsing error : fail to get 'Sampling_Interval' in 'Monitor_List'\n");
      Exit(0);
    }
  }
  
  // 単位指定
  if ( !elmL1->GetValue("Unit", &str) ) {
		Hostonly_ stamped_printf("\tParsing error : Invalid string for 'Unit' in 'Monitor_List'\n");
		Exit(0);
	}
  if     ( !strcasecmp(str, "Dimensional") ) {
    C.Sampling.unit = DIMENSIONAL;
    M->setSamplingUnit(DIMENSIONAL);
  }
  else if( !strcasecmp(str, "Non_Dimensional") ) {
    C.Sampling.unit = NONDIMENSIONAL;
    M->setSamplingUnit(NONDIMENSIONAL);
  }
  else {
    Hostonly_ stamped_printf("\tInvalid keyword is described at 'Unit' section\n");
    Exit(0);
  }
  
  // サンプリングの指定単位が有次元の場合に，無次元に変換
  if ( C.Sampling.unit == DIMENSIONAL ) {
    C.Interval[Interval_Manager::tg_sampled].normalizeInterval(C.Tscale);
  }
  
  // 指定モニタ個数のチェック
  if ( (elmL1->GetElemSize() < 1) && (C.isMonitor() == OFF) ) {
    Hostonly_ stamped_printf("\tError : No monitoring points. Please confirm 'Monitor_List' and 'InnerBoundary' in XML. \n");
    Exit(0);
  }
  
  
  // モニターリストの読み込み
  elmL2 = elmL1->GetElemFirst();
  
  while (elmL2) {
    
    // sampling type & param check
    p = elmL2->GetName();
    if ( !strcasecmp(p, "point_set") ) {
      type = MonitorCompo::POINT_SET;
      if ( 0 == elmL2->GetElemSize() ) {
        Hostonly_ stamped_printf("\tParsing error : At least, 1 elem of 'set' should be found in 'point_set'\n");
        Exit(0);
      }
    }
    else if ( !strcasecmp(p, "line") ) {
      type = MonitorCompo::LINE;
      if ( 2 != elmL2->GetElemSize() ) {
        Hostonly_ stamped_printf("\tParsing error : 2 elems (from/to) should be found in 'line'\n");
        Exit(0);
      }
    }
    else {
      Hostonly_ stamped_printf("\tParsing error : No valid keyword [point_set / line] in 'Monitor_List'\n");
      Exit(0);
    }
    
    // Labelの取得．ラベルなしでもエラーではない
    if ( !(label = elmL2->GetComment()) ) {
      Hostonly_ stamped_printf("\tParsing warning : No commnet in '%s'\n", p);
    }
    
    // variable
    variables.clear();
    param = elmL2->GetParamFirst("variable");
    while (param) {
      str=NULL;
      if ( !param->GetData(&str) ) {
        Hostonly_ stamped_printf("\tParsing error : fail to get 'variable' in 'Monitor_List'\n");
        Exit(0);
      }
      variables.push_back(str);
      param = elmL2->GetParamNext(param, "variable");
    }
    if (variables.size() == 0) {
      Hostonly_ stamped_printf("\tParsing error : No 'variable' in 'Monitor_List'\n");
      Exit(0);
    }
    
    // method
    if (!elmL2->GetValue("sampling_method", &method)) {
      Hostonly_ stamped_printf("\tParsing error : fail to get 'sampling_method' in 'Monitor_List'\n");
      Exit(0);
    }
    
    // mode
    if (!elmL2->GetValue("sampling_mode", &mode)) {
      Hostonly_ stamped_printf("\tParsing error : fail to get 'sampling_mode' in 'Monitor_List'\n");
      Exit(0);
    }
    
    // get coordinate
    if ( type == MonitorCompo::POINT_SET ) {
      vector<MonitorCompo::MonitorPoint> pointSet;
      getXML_Mon_Pointset(M, elmL2, pointSet);
      M->setPointSet(label, variables, method, mode, pointSet);
    }
    else {
      REAL_TYPE from[3], to[3];
      int nDivision;
      getXML_Mon_Line(M, elmL2, from, to, nDivision);
      M->setLine(label, variables, method, mode, from, to, nDivision);
    }
    
    elmL2 = elmL1->GetElemNext(elmL2); // ahead on the next pointer
  }
  
}

/**
 @fn void SklSolverCBC::getXML_Mon_Line(MonitorList* M, const CfgElem* elmL2, REAL_TYPE from[3], REAL_TYPE to[3], int& nDivision)
 @brief XMLに記述されたモニタ座標情報(Line)を取得
 @param M MonitorList クラスオブジェクトのポインタ
 @param elmL2 コンフィギュレーションツリーのポインタ
 @param from Line始点座標
 @param to   Line終点座標
 @param nDivision Line分割数
 @note データは無次元化して保持
 */
void SklSolverCBC::getXML_Mon_Line(MonitorList* M, const CfgElem* elmL2, REAL_TYPE from[3], REAL_TYPE to[3], int& nDivision)
{
  const CfgElem *elmL3=NULL;
  
  if ( !elmL2->GetValue("division", &nDivision) ) Exit(0);
  if ( nDivision == 0 ) Exit(0);
  
  // load parameter of 'from' and 'to'
  if ( !elmL2->GetValue("from", &elmL3) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'from' in 'line' >> %s\n", elmL3->GetName());
    Exit(0);
  }
  else {
    if ( !elmL3->GetVctValue("x", "y", "z", &from[0], &from[1], &from[2]) ) {
      Hostonly_ stamped_printf("\tParsing error : fail to get vec params in 'from'\n");
      Exit(0);
    }
    if (C.Sampling.unit == DIMENSIONAL) {
      normalizeCord(from);
    }
  }
  if ( !elmL2->GetValue("to", &elmL3) ) {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'to' in 'line' >> %s\n", elmL3->GetName());
    Exit(0);
  }
  else {
    if ( !elmL3->GetVctValue("x", "y", "z", &to[0], &to[1], &to[2]) ) {
      Hostonly_ stamped_printf("\tParsing error : fail to get vec params in 'to'\n");
      Exit(0);
    }
    if (C.Sampling.unit == DIMENSIONAL) {
      normalizeCord(to);
    }
  }
}

/**
 @fn void SklSolverCBC::getXML_Mon_Pointset(MonitorList* M, const CfgElem *elmL2, vector<MonitorCompo::MonitorPoint>& pointSet)
 @brief XMLに記述されたモニタ座標情報を取得(PointSet)
 @param M MonitorList クラスオブジェクトのポインタ
 @param m order
 @param elmL2 コンフィギュレーションツリーのポインタ
 @param pointSet PointSet配列
 @note データは無次元化して保持
 */
void SklSolverCBC::getXML_Mon_Pointset(MonitorList* M, const CfgElem *elmL2, vector<MonitorCompo::MonitorPoint>& pointSet)
{
  const CfgElem *elmL3=NULL;
  REAL_TYPE v[3];
  const char* str=NULL;
  char tmpstr[20];
  
  v[0] = v[1] = v[2] = 0.0;
  
  // load parameter for a set
  elmL3 = elmL2->GetElemFirst();
  
  for (unsigned j=0; j<elmL2->GetElemSize(); j++) {
    
    if ( strcasecmp("set", elmL3->GetName()) ) { // not agree
      Hostonly_ stamped_printf("\tParsing error : fail to get 'set' in 'point_set' >> %s\n", elmL3->GetName());
      Exit(0);
    }
    if ( !elmL3->GetVctValue("x", "y", "z", &v[0], &v[1], &v[2]) ) {
      Hostonly_ stamped_printf("\tParsing error : fail to get vec params in 'point_set'\n");
      Exit(0);
    }
    if (C.Sampling.unit == DIMENSIONAL) {
      normalizeCord(v);
    }
    
    // set Labelの取得．ラベルなしでもエラーではない
    if ( !(str = elmL3->GetComment()) ) {
      Hostonly_ stamped_printf("\tParsing warning : No commnet for 'set'\n");
    }
    if ( !str ) {
      sprintf(tmpstr, "point_%d",j);
      str = tmpstr;
    }
    
    pointSet.push_back(MonitorCompo::MonitorPoint(v, str));
    
    elmL3 = elmL2->GetElemNext(elmL3); // ahead on the next pointer
  }  
}

/**
 @fn bool SklSolverCBC::hasLinearSolver(unsigned L)
 @brief 種類Lの線形ソルバを利用する場合，trueを返す
 @param L 線形ソルバの種類
 */
bool SklSolverCBC::hasLinearSolver(unsigned L)
{
  for (int i=0; i<ItrCtl::ic_END; i++)
    if ( IC[i].get_LS() == L ) return true;
  
  return false;
}

/**
 @fn void SklSolverCBC::Restart (FILE* fp, REAL_TYPE& flop)
 @brief リスタート時の瞬時値ファイル読み込み
 @param fp ファイルポインタ
 @param flop
 */
void SklSolverCBC::Restart (FILE* fp, REAL_TYPE& flop)
{
  REAL_TYPE time;
  std::string tmp;
  
  REAL_TYPE *v = NULL;
  REAL_TYPE *p = NULL;
  
  if( !(v = dc_v->GetData()) )   Exit(0);
  if( !(p = dc_p->GetData()) )   Exit(0);
  
  // ステップ数 > ファイル名用
  int m_step = (int)C.Restart_step;
  int step;
  
  // 圧力の瞬時値　
  REAL_TYPE bp = ( C.Unit.Prs == Unit_Absolute ) ? C.BasePrs : 0.0;
  
  // ガイド出力
  int gs = (int)C.GuideOut;
  
  // dummy
  int i_dummy=0;
  REAL_TYPE f_dummy=0.0;
  
  tmp = DFI.Generate_FileName(C.f_Pressure, m_step, pn.ID, (bool)C.FIO.IO_Input);

  if ( !checkFile(tmp) ) {
    Hostonly_ printf("\n\tError : File open '%s'\n", tmp.c_str());
    Exit(0);
  }
  
  F.readPressure(fp, tmp, size, guide, p, step, time, C.Unit.File, bp, C.RefDensity, C.RefVelocity, flop, gs, true, i_dummy, f_dummy);

  // ここでタイムスタンプを得る
  if (C.Unit.File == DIMENSIONAL) time /= C.Tscale;
  SklSetBaseStep(step);
  SklSetBaseTime(time);
  
  // v00[]に値をセット
  copyV00fromRF((double)time);


  // Instantaneous Velocity fields
  tmp = DFI.Generate_FileName(C.f_Velocity, m_step, pn.ID, (bool)C.FIO.IO_Input);
  
  if ( !checkFile(tmp) ) {
    Hostonly_ printf("\n\tError : File open '%s'\n", tmp.c_str());
    Exit(0);
  }

  F.readVelocity(fp, tmp, size, guide, v, step, time, v00, C.Unit.File, C.RefVelocity, flop, gs, true, i_dummy, f_dummy);

  if (C.Unit.File == DIMENSIONAL) time /= C.Tscale;
  
  if ( (step != SklGetTotalStep()) || (time != (REAL_TYPE)SklGetTotalTime()) ) {
    Hostonly_ printf     ("\n\tTime stamp is different between files\n");
    Hostonly_ fprintf(fp, "\n\tTime stamp is different between files\n");
    Exit(0);
  }
  
  // Instantaneous Temperature fields
  if ( C.isHeatProblem() ) {
    
    REAL_TYPE *t = NULL;
    if( !(t = dc_t->GetData()) ) Exit(0);
    
    REAL_TYPE klv = ( C.Unit.Temp == Unit_KELVIN ) ? 0.0 : KELVIN;
    
    tmp = DFI.Generate_FileName(C.f_Temperature, m_step, pn.ID, (bool)C.FIO.IO_Input);
    
    if ( !checkFile(tmp) ) {
      Hostonly_ printf("\n\tError : File open '%s'\n", tmp.c_str());
      Exit(0);
    }
    F.readTemperature(fp, tmp, size, guide, t, step, time, C.Unit.File, C.BaseTemp, C.DiffTemp, klv, flop, gs, true, i_dummy, f_dummy);
    
    if (C.Unit.File == DIMENSIONAL) time /= C.Tscale;
    
    if ( (step != SklGetTotalStep()) || (time != (REAL_TYPE)SklGetTotalTime()) ) {
      Hostonly_ printf     ("\n\tTime stamp is different between files\n");
      Hostonly_ fprintf(fp, "\n\tTime stamp is different between files\n");
      Exit(0);
    }
  }
  
}


/**
 @fn void SklSolverCBC::Restart_coarse (FILE* fp, REAL_TYPE& flop)
 @brief 粗い格子を用いたリスタート
 @param fp ファイルポインタ
 @param flop
 */
void SklSolverCBC::Restart_coarse (FILE* fp, REAL_TYPE& flop)
{
  int step;
  REAL_TYPE time;
  
  
  // ステップ数 > ファイル名用
  int m_step = (int)C.Restart_step;
  
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
    i = pn.st_idx[0];
    j = pn.st_idx[1];
    k = pn.st_idx[2];
    
    // crs_i, _j, _kには同じ値が入る 
    getCoarseResult(i, j, k, C.f_Coarse_dfi_prs, C.f_Coarse_pressure, m_step, f_prs, r_size, crs, num_block);

    getCoarseResult(i, j, k, C.f_Coarse_dfi_vel, C.f_Coarse_velocity, m_step, f_vel, r_size, crs, num_block);
    
    if ( C.isHeatProblem() ) {
      getCoarseResult(i, j, k, C.f_Coarse_dfi_temp, C.f_Coarse_temperature, m_step, f_temp, r_size, crs, num_block);
    }
  }
  else {
    crs[0] = 1;
    crs[1] = 1;
    crs[2] = 1;
    f_prs = DFI.Generate_FileName(C.f_Coarse_pressure, m_step, pn.ID);
    f_vel = DFI.Generate_FileName(C.f_Coarse_velocity, m_step, pn.ID);
    if ( C.isHeatProblem() )
      f_temp= DFI.Generate_FileName(C.f_Coarse_temperature, m_step, pn.ID);
  }

  // テンポラリのファイルロード
  allocArray_CoarseMesh((unsigned*)r_size);

  
  REAL_TYPE *r_v = NULL;
  REAL_TYPE *r_p = NULL;
  
  if( !(r_v = dc_r_v->GetData()) )   Exit(0);
  if( !(r_p = dc_r_p->GetData()) )   Exit(0);
  
  // ガイド出力
  int gs = (int)C.GuideOut;
  
  
  // dummy
  int i_dummy=0;
  REAL_TYPE f_dummy=0.0;
  
  // 圧力の瞬時値　ここでタイムスタンプを得る
  REAL_TYPE bp = ( C.Unit.Prs == Unit_Absolute ) ? C.BasePrs : 0.0;

  if ( !checkFile(f_prs) ) {
    Hostonly_ printf("\n\tError : File open '%s'\n", f_prs.c_str());
    Exit(0);
  }

  F.readPressure(fp, f_prs, (unsigned*)r_size, guide, r_p, step, time, C.Unit.File, bp, C.RefDensity, C.RefVelocity, flop, gs, true, i_dummy, f_dummy);

  if (C.Unit.File == DIMENSIONAL) time /= C.Tscale;
  SklSetBaseStep(step);
  SklSetBaseTime(time);

  // v00[]に値をセット
  copyV00fromRF((double)time);
  
  
  // Instantaneous Velocity fields
  if ( !checkFile(f_vel) ) {
    Hostonly_ printf("\n\tError : File open '%s'\n", f_vel.c_str());
    Exit(0);
  }
  F.readVelocity(fp, f_vel, (unsigned*)r_size, guide, r_v, step, time, v00, C.Unit.File, C.RefVelocity, flop, gs, true, i_dummy, f_dummy);
  
  if (C.Unit.File == DIMENSIONAL) time /= C.Tscale;
  if ( (step != SklGetTotalStep()) || (time != (REAL_TYPE)SklGetTotalTime()) ) {
    Hostonly_ printf     ("\n\tTime stamp is different between files\n");
    Hostonly_ fprintf(fp, "\n\tTime stamp is different between files\n");
    Exit(0);
  }
  
  // Instantaneous Temperature fields
  if ( C.isHeatProblem() ) {
    
    REAL_TYPE *r_t = NULL;
    if( !(r_t = dc_r_t->GetData()) ) Exit(0);
    
    REAL_TYPE klv = ( C.Unit.Temp == Unit_KELVIN ) ? 0.0 : KELVIN;

    if ( !checkFile(f_temp) ) {
      Hostonly_ printf("\n\tError : File open '%s'\n", f_temp.c_str());
      Exit(0);
    }
    F.readTemperature(fp, f_temp, (unsigned*)r_size, guide, r_t, step, time, C.Unit.File, C.BaseTemp, C.DiffTemp, klv, flop, gs, true, i_dummy, f_dummy);
    
    if (C.Unit.File == DIMENSIONAL) time /= C.Tscale;
    if ( (step != SklGetTotalStep()) || (time != (REAL_TYPE)SklGetTotalTime()) ) {
      Hostonly_ printf     ("\n\tTime stamp is different between files\n");
      Hostonly_ fprintf(fp, "\n\tTime stamp is different between files\n");
      Exit(0);
    }
  }
  
  // 同期
  if( !dc_r_v->CommBndCell(guide) ) Exit(0);
  if( !dc_r_p->CommBndCell(guide) ) Exit(0);
  if ( C.isHeatProblem() ) {
    if( !dc_r_t->CommBndCell(guide) ) Exit(0);
  }
  
  // 内挿処理
  Interpolation_from_coarse_initial(crs, num_block);
}


/**
 @fn void SklSolverCBC::Restart_avrerage (FILE* fp, REAL_TYPE& flop)
 @brief リスタート時の平均値ファイル読み込み
 @param fp ファイルポインタ
 @param flop
 */
void SklSolverCBC::Restart_avrerage (FILE* fp, REAL_TYPE& flop)
{
  std::string tmp;
  
  int step = SklGetBaseStep();
  REAL_TYPE time = SklGetBaseTime();
  
  // ステップ数 > ファイル名用
  int m_step = (int)C.Restart_step;
  
  // ガイド出力
  int gs = (int)C.GuideOut;
  
  REAL_TYPE *av = NULL;
  REAL_TYPE *ap = NULL;
  
  if( !(av = dc_av->GetData()) )   Exit(0);
  if( !(ap = dc_ap->GetData()) )   Exit(0);

  if ( C.Interval[Interval_Manager::tg_avstart].isStep() ) {
    if ( step > C.Interval[Interval_Manager::tg_avstart].getIntervalStep() ) {
      Hostonly_ printf     ("\tRestart from Previous Calculation Results of averaged field\n");
      Hostonly_ fprintf(fp, "\tRestart from Previous Calculation Results of averaged field\n");
      Hostonly_ printf     ("\tStep : base=%d current=%d total=%d\n", step, SklGetCurrentStep(), SklGetTotalStep());
      Hostonly_ fprintf(fp, "\tStep : base=%d current=%d total=%d\n", step, SklGetCurrentStep(), SklGetTotalStep());
    }
    else {
      return;
    }
  }
  else {
    if ( time > C.Interval[Interval_Manager::tg_avstart].getIntervalTime() ) {
      Hostonly_ printf     ("\tRestart from Previous Calculation Results of averaged field\n");
      Hostonly_ fprintf(fp, "\tRestart from Previous Calculation Results of averaged field\n");
      Hostonly_ printf     ("\tTime : base=%e[sec.]/%e[-] current=%e[-] total=%e[-]\n", time*C.Tscale, time, SklGetCurrentTime(), SklGetTotalTime());
      Hostonly_ fprintf(fp, "\tTime : base=%e[sec.]/%e[-] current=%e[-] total=%e[-]\n", time*C.Tscale, time, SklGetCurrentTime(), SklGetTotalTime());
    }
    else {
      return;
    }
  }

  int step_avr = 0;
  REAL_TYPE time_avr = 0.0;
  
  // Pressure
  REAL_TYPE bp = ( C.Unit.Prs == Unit_Absolute ) ? C.BasePrs : 0.0;
  
  tmp = DFI.Generate_FileName(C.f_AvrPressure, m_step, pn.ID, (bool)C.FIO.IO_Input);
  if ( !checkFile(tmp) ) {
    Hostonly_ printf("\n\tError : File open '%s'\n", tmp.c_str());
    Exit(0);
  }
  F.readPressure(fp, tmp, size, guide, ap, step, time, C.Unit.File, bp, C.RefDensity, C.RefVelocity, flop, gs, false, step_avr, time_avr);
  
  if ( (step != SklGetBaseStep()) || ((REAL_TYPE)time != SklGetBaseTime()) ) {
    Hostonly_ printf     ("\n\tTime stamp is different between instantaneous and averaged files\n");
    Hostonly_ fprintf(fp, "\n\tTime stamp is different between instantaneous and averaged files\n");
    Exit(0);
  }
  
  if (C.Unit.File == DIMENSIONAL) time /= C.Tscale;
  
  SklSetAverageBaseStep(step_avr);
  SklSetAverageBaseTime(time_avr);
  
  
  // Velocity
  tmp = DFI.Generate_FileName(C.f_AvrVelocity, m_step, pn.ID, (bool)C.FIO.IO_Input);
  if ( !checkFile(tmp) ) {
    Hostonly_ printf("\n\tError : File open '%s'\n", tmp.c_str());
    Exit(0);
  }
  F.readVelocity(fp, tmp, size, guide, av, step, time, v00, C.Unit.File, C.RefVelocity, flop, gs, false, step_avr, time_avr);
  
  if (C.Unit.File == DIMENSIONAL) time /= C.Tscale;
  
  if ( (step_avr != SklGetAverageBaseStep()) || ((REAL_TYPE)time_avr != SklGetAverageBaseTime()) ) { // 圧力とちがう場合
    Hostonly_ printf     ("\n\tTime stamp is different between files\n");
    Hostonly_ fprintf(fp, "\n\tTime stamp is different between files\n");
    Exit(0);
  }
  

  if ( C.isHeatProblem() ) {
    
    REAL_TYPE *at = NULL;
    if( !(at = dc_at->GetData()) ) Exit(0);
    
    REAL_TYPE klv = ( C.Unit.Temp == Unit_KELVIN ) ? 0.0 : KELVIN;

    tmp = DFI.Generate_FileName(C.f_AvrTemperature, m_step, pn.ID, (bool)C.FIO.IO_Input);
    if ( !checkFile(tmp) ) {
      Hostonly_ printf("\n\tError : File open '%s'\n", tmp.c_str());
      Exit(0);
    }
    F.readTemperature(fp, tmp, size, guide, at, step, time, C.Unit.File, C.BaseTemp, C.DiffTemp, klv, flop, gs, false, step_avr, time_avr);
    
    if (C.Unit.File == DIMENSIONAL) time /= C.Tscale;
    
    if ( (step_avr != SklGetAverageBaseStep()) || ((REAL_TYPE)time_avr != SklGetAverageBaseTime()) ) {
      Hostonly_ printf     ("\n\tTime stamp is different between files\n");
      Hostonly_ fprintf(fp, "\n\tTime stamp is different between files\n");
      Exit(0);
    }
  }
}


/**
 @fn void SklSolverCBC::min_distance(float* cut, FILE* fp)
 @brief 距離の最小値を求める
 @param cut カット情報の配列
 @param fp file pointer
 @note cutlibのインスタンスを参照, windows does not support std::min() method?
 */
void SklSolverCBC::min_distance(float* cut, FILE* fp)
{
  //SklParaComponent* para_cmp = SklGetParaComponent();
  //SklParaManager* para_mng = para_cmp->GetParaManager();
  
  unsigned mm = (size[2]+2*guide)*(size[1]+2*guide)*(size[0]+2*guide)*6;
  
  if (mm >= UINT_MAX) {
    Hostonly_ stamped_printf("\n\tError : Product of size[]*6 exceeds UINT_MAX\n\n");
    Exit(0);
  }
  
  float global_min = 1.0;
  float c;
  float eps = 1.0/255.0; // 0.000392
  unsigned g=0;
  
#pragma omp parallel firstprivate(mm) private(c)
  {
    float local_min = 1.0;
  
#pragma omp for nowait schedule(static) reduction(+:g)
    for (unsigned i=0; i<mm; i++) {
      c = cut[i]; 
      if ( local_min > c ) local_min = c;
      if ( (c > 0.0f) && (c <= eps) ) {
        cut[i] = eps;
        g++;
      }
    }
#pragma omp critical
    {
      if ( local_min < global_min ) global_min = local_min;
    }
  } // end omp parallel

  int m_np;
  MPI_Comm_size(MPI_COMM_WORLD, &m_np);
  
  if ( m_np > 1 ) {
    
    float tmp = global_min;
    MPI_Allreduce(&tmp, &global_min, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
    
    unsigned tmp_g = g;
    VoxInfo::uint_sum_Allreduce(&tmp_g, &g);
  }

  Hostonly_ fprintf(fp, "\n\tMinimum non-dimnensional distance is %e and replaced to %e : num = %d\n\n", global_min, eps, g);
  Hostonly_ printf     ("\n\tMinimum non-dimnensional distance is %e and replaced to %e : num = %d\n\n", global_min, eps, g);
}



/**
 @fn void SklSolverCBC::resizeBVface(const int* st, const int* ed, const unsigned n, const unsigned* bx)
 @brief コンポーネントリストに登録されたセルフェイスBCのBV情報をリサイズする
 @param st 開始インデクス
 @param ed 終了インデクス
 @param n CompoListのエントリ
 @param bx BCindex
 */
void SklSolverCBC::resizeBVface(const int* st, const int* ed, const unsigned n, const unsigned* bx)
{
  int i,j,k;
  unsigned register s, m;
  int nst[3], ned[3];
  
  // 初期値はローカルノードの大きさ
  nst[0] = (int)size[0];
  nst[1] = (int)size[1];
  nst[2] = (int)size[2];
  ned[0] = 0;
  ned[1] = 0;
  ned[2] = 0;
  
  for (k=st[2]; k<=ed[2]; k++) {
    for (j=st[1]; j<=ed[1]; j++) {
      for (i=st[0]; i<=ed[0]; i++) {
        
        m = FBUtility::getFindexS3D(size, guide, i, j, k);
        s = bx[m];
        
        if ( GET_FACE_BC(s, BC_FACE_W) == n ) {
          if( i < nst[0] ) { nst[0] = i; }
          if( i > ned[0] ) { ned[0] = i; }
          if( j < nst[1] ) { nst[1] = j; }
          if( j > ned[1] ) { ned[1] = j; }
          if( k < nst[2] ) { nst[2] = k; }
          if( k > ned[2] ) { ned[2] = k; }
        }
        
        if ( GET_FACE_BC(s, BC_FACE_E) == n ) {
          if( i < nst[0] ) { nst[0] = i; }
          if( i > ned[0] ) { ned[0] = i; }
          if( j < nst[1] ) { nst[1] = j; }
          if( j > ned[1] ) { ned[1] = j; }
          if( k < nst[2] ) { nst[2] = k; }
          if( k > ned[2] ) { ned[2] = k; }
        }
        
        if ( GET_FACE_BC(s, BC_FACE_S) == n ) {
          if( i < nst[0] ) { nst[0] = i; }
          if( i > ned[0] ) { ned[0] = i; }
          if( j < nst[1] ) { nst[1] = j; }
          if( j > ned[1] ) { ned[1] = j; }
          if( k < nst[2] ) { nst[2] = k; }
          if( k > ned[2] ) { ned[2] = k; }
        }
        
        if ( GET_FACE_BC(s, BC_FACE_N) == n ) {
          if( i < nst[0] ) { nst[0] = i; }
          if( i > ned[0] ) { ned[0] = i; }
          if( j < nst[1] ) { nst[1] = j; }
          if( j > ned[1] ) { ned[1] = j; }
          if( k < nst[2] ) { nst[2] = k; }
          if( k > ned[2] ) { ned[2] = k; }
        }
        
        if ( GET_FACE_BC(s, BC_FACE_B) == n ) {
          if( i < nst[0] ) { nst[0] = i; }
          if( i > ned[0] ) { ned[0] = i; }
          if( j < nst[1] ) { nst[1] = j; }
          if( j > ned[1] ) { ned[1] = j; }
          if( k < nst[2] ) { nst[2] = k; }
          if( k > ned[2] ) { ned[2] = k; }
        }
        
        if ( GET_FACE_BC(s, BC_FACE_T) == n ) {
          if( i < nst[0] ) { nst[0] = i; }
          if( i > ned[0] ) { ned[0] = i; }
          if( j < nst[1] ) { nst[1] = j; }
          if( j > ned[1] ) { ned[1] = j; }
          if( k < nst[2] ) { nst[2] = k; }
          if( k > ned[2] ) { ned[2] = k; }
        }
      }
    }
  }
  
  // replace
  cmp[n].setBbox(nst, ned);
}

/**
 @fn void SklSolverCBC::resizeBVcell(const int* st, const int* ed, const unsigned n, const unsigned* bx)
 @brief コンポーネントリストに登録されたセル要素BCのBV情報をリサイズする
 @param st 開始インデクス
 @param ed 終了インデクス
 @param n CompoListのエントリ
 @param bx BCindex
 */
void SklSolverCBC::resizeBVcell(const int* st, const int* ed, const unsigned n, const unsigned* bx)
{
  int i,j,k;
  unsigned register s, m;
  int nst[3], ned[3];
  
  // 初期値はローカルノードの大きさ
  nst[0] = (int)size[0];
  nst[1] = (int)size[1];
  nst[2] = (int)size[2];
  ned[0] = 0;
  ned[1] = 0;
  ned[2] = 0;

  for (k=st[2]; k<=ed[2]; k++) {
    for (j=st[1]; j<=ed[1]; j++) {
      for (i=st[0]; i<=ed[0]; i++) {
        
        m = FBUtility::getFindexS3D(size, guide, i, j, k);
        s = bx[m];
        
        if ( ( s & MASK_6) == n ) {
          if( i < nst[0] ) { nst[0] = i; }
          if( i > ned[0] ) { ned[0] = i; }
          if( j < nst[1] ) { nst[1] = j; }
          if( j > ned[1] ) { ned[1] = j; }
          if( k < nst[2] ) { nst[2] = k; }
          if( k > ned[2] ) { ned[2] = k; }
        }
      }
    }
  }

  // replace
  cmp[n].setBbox(nst, ned);
}

/**
 @fn void SklSolverCBC::resizeCompoBV(const unsigned* bd, const unsigned* bv, const unsigned* bh1, const unsigned* bh2, const unsigned kos, const bool isHeat)
 @brief コンポーネントリストに登録されたBV情報をリサイズする
 @param bd  BCindex bcd
 @param bv  BCindex bcv
 @param bh1 BCindex bh1
 @param bh2 BCindex bh2
 @param kos KOS
 @param isHeat 熱問題のときtrue
 */
void SklSolverCBC::resizeCompoBV(const unsigned* bd, const unsigned* bv, const unsigned* bh1, const unsigned* bh2, const unsigned kos, const bool isHeat)
{
  int st[3], ed[3];
  unsigned typ, id;
  
  for (unsigned n=1; n<=C.NoBC; n++) {
    id = cmp[n].getID();
    cmp[n].getBbox(st, ed);
    typ = cmp[n].getType();
    
    // デフォルトで流れ計算パートのBC
    switch ( typ ) {
      case SPEC_VEL:
      case SPEC_VEL_WH:
      case OUTFLOW:
        resizeBVface(st, ed, n, bv); // 速度のBCindex　セルフェイス位置でのflux型BC
        break;
        
      case PERIODIC:
      case IBM_DF:
      case HEX:
      case FAN:
      case DARCY:
        resizeBVcell(st, ed, n, bd); // セルセンタ位置でのBC
        break;
    }
    
    // 熱計算のパート
    if ( isHeat ) {
      switch ( typ ) {
        case ADIABATIC:
        case HEATFLUX:
        case SPEC_VEL_WH:
        case OUTFLOW:
        case TRANSFER:
        case ISOTHERMAL:
          resizeBVface(st, ed, n, bh1);
          break;
          
        case RADIANT:
          break;
          
        case HEAT_SRC:
        case CNST_TEMP:
        case PERIODIC:
          resizeBVcell(st, ed, n, bh2);
          break;
      }            
    }
    
  }
}


/**
 @fn void SklSolverCBC::Interpolation_from_coarse_initial(const int* m_st)
 @param m_st 粗い格子の開始インデクス
 @param m_bk ブロック数
 @brief 粗い格子を用いたリスタート値の内挿
 */
void SklSolverCBC::Interpolation_from_coarse_initial(const int* m_st, const int* m_bk)
{
  REAL_TYPE *v = NULL;
  REAL_TYPE *p = NULL;
  REAL_TYPE *t = NULL;
  REAL_TYPE *r_v = NULL;
  REAL_TYPE *r_p = NULL;
  REAL_TYPE *r_t = NULL;
  
  int st[3], bk[3];
  st[0] = m_st[0];
  st[1] = m_st[1];
  st[2] = m_st[2];
  bk[0] = m_bk[0];
  bk[1] = m_bk[1];
  bk[2] = m_bk[2];
  
  if( !(v   = dc_v->GetData()) )     Exit(0);
  if( !(p   = dc_p->GetData()) )     Exit(0);
  if( !(r_v = dc_r_v->GetData()) )   Exit(0);
  if( !(r_p = dc_r_p->GetData()) )   Exit(0);
  
  fb_interp_coarse_s_(p, sz, gc, r_p, st, bk);
  fb_interp_coarse_v_(v, sz, gc, r_v, st, bk);
  
  if ( C.isHeatProblem() ) {
    if( !(t   = dc_t->GetData()) )   Exit(0);
    if( !(r_t = dc_r_t->GetData()) ) Exit(0);
    fb_interp_coarse_s_(t, sz, gc, r_t, st, bk);
  }

}


/**
 @fn void SklSolverCBC::setBCinfo(ParseBC* B)
 @brief ParseBCクラスをセットアップし，外部境界条件を読み込み，Controlクラスに保持する
 @param B
 */
void SklSolverCBC::setBCinfo(ParseBC* B)
{
  // CompoListクラスのオブジェクトポインタを渡す
  B->receiveCompoPtr(cmp);

  // XMLの情報を元にCompoListの情報を設定する
  B->setCompoList(&C);
  
  // 各コンポーネントが存在するかどうかを保持しておく
  setEnsComponent();

  // BoundaryOuterクラスのポインタを渡す
  B->setObcPtr(BC.get_OBC_Ptr());

  // XMLファイルをパースして，外部境界条件を保持する
  B->loadOuterBC();
  
  // KOSと境界条件種類の整合性をチェック
  B->chkBCconsistency(C.KindOfSolver);
}


/**
 @fn void SklSolverCBC::setComponentVF(float* cvf)
 @brief HEX,FANコンポーネントなどの体積率とbboxなどをセット
 @param cvf 体積率
 @note インデクスの登録と配列確保はVoxEncode()で、コンポーネント領域のリサイズ後に行う
 */
void SklSolverCBC::setComponentVF(float* cvf)
{
  int subsampling = 20; // 体積率のサブサンプリングの基数
  int f_st[3], f_ed[3];
  double flop;
  
  CompoFraction CF(size, guide, C.dx, C.org, subsampling);
  
  for (int n=1; n<=C.NoBC; n++) {
    
    if ( cmp[n].isFORCING() ) {
      // 形状パラメータのセット
      switch ( cmp[n].getType() ) {
        case HEX:
          CF.setShapeParam(cmp[n].nv, cmp[n].oc, cmp[n].dr, cmp[n].depth, cmp[n].shp_p1, cmp[n].shp_p2);
          break;
          
        case FAN:
          CF.setShapeParam(cmp[n].nv, cmp[n].oc, cmp[n].depth, cmp[n].shp_p1, cmp[n].shp_p2);
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
      cmp[n].setEns(ON);
      
      // 体積率
      TIMING_start(tm_cmp_vertex8);
      flop = 0.0;
      CF.vertex8    (f_st, f_ed, cvf, flop);
      TIMING_stop(tm_cmp_vertex8, (REAL_TYPE)flop);
      
      TIMING_start(tm_cmp_subdivision);
      flop = 0.0;
      CF.subdivision(f_st, f_ed, cvf, flop);
      TIMING_stop(tm_cmp_subdivision, (REAL_TYPE)flop);
    }
  }
  
  // 確認のための出力
#if 0
  REAL_TYPE org[3], pit[3];
  
  //  ガイドセルがある場合(GuideOut != 0)にオリジナルポイントを調整
  for (int i=0; i<3; i++) {
    org[i] = C.org[i] - C.dx[i]*(REAL_TYPE)C.GuideOut;
    pit[i] = C.dx[i];
  }
  
  // 出力ファイルの指定が有次元の場合
  if ( C.Unit.File == DIMENSIONAL ) {
    for (int i=0; i<3; i++) {
      org[i] *= C.RefLength;
      pit[i] *= C.RefLength;
    }
  }
  F.writeRawSPH(cvf, size, guide, org, pit, SPH_SINGLE);
#endif
  
}


/**
 @fn void SklSolverCBC::setMaterialList(ParseBC* B, ParseMat* M, FILE* mp, FILE* fp)
 @brief ParseMatクラスをセットアップし，媒質情報をXMLから読み込み，媒質リストを作成する
 @param B
 @param M 
 @param mp
 @param fp 
 */
void SklSolverCBC::setMaterialList(ParseBC* B, ParseMat* M, FILE* mp, FILE* fp)
{
  // ParseMatクラスの環境設定 
  M->setControlVars(&C, B->get_IDtable_Ptr(), mat, m_solvCfg);
  
  // Material情報の内容をXMLファイルをパースして，MaterialListクラスのオブジェクトBaseMatに保持する
  M->getXMLmaterial();
  
#if 0
  // Materialの基本リストを表示
  Hostonly_ M->dbg_printBaseMaterialList(mp, fp);
#endif
  
  // MaterialListを作成する
  M->makeMaterialList();
  
  // コンポーネントとMaterialリストの関連づけ（相互参照リスト）を作成する
  M->makeLinkCmpMat(cmp);
  
  // 媒質テーブルの表示
  Hostonly_ M->printMaterialList(mp, fp);
  
  // Model_Settingで指定した媒質とiTableのStateの不一致をチェック
  M->chkState_Mat_Cmp(cmp);
}


/**
 @fn void SklSolverCBC::setBbox_of_VIBC_Cut(VoxInfo* Vinfo, const unsigned* bv)
 @brief ポリゴンからVBCのboxをセット
 @param Vinfo
 @param bv BCindex V
 @note インデクスの登録はVoxEncode()で、コンポーネント領域のリサイズ後に行う
 */
void SklSolverCBC::setBbox_of_VIBC_Cut(VoxInfo* Vinfo, const unsigned* bv)
{
  int f_st[3], f_ed[3];
  
  for (int n=1; n<=C.NoBC; n++) {
    
    if ( cmp[n].isVBC_IO() ) { // SPEC_VEL || SPEC_VEL_WH || OUTFLOW

      // インデクスの計算 > あとで，VoxEncode()でresize
      Vinfo->findVIBCbbox(n, bv, f_st, f_ed);
      
      // インデクスのサイズ登録と存在フラグ
      cmp[n].setBbox(f_st, f_ed);
      cmp[n].setEns(ON);
    }
  }
  
}


/**
 @fn void SklSolverCBC::setEnsComponent(void)
 @brief コンポーネントが存在するかを保持しておく
 */
void SklSolverCBC::setEnsComponent(void)
{
  unsigned c;
  
  // Forcing > HEX, FAN, DARCY
  c = 0;
  for (int n=1; n<=C.NoBC; n++) {
    if ( cmp[n].isFORCING() ) c++;
  }
  if ( c>0 ) C.EnsCompo.forcing = ON;
  
  // Heat source > HEAT_SRC, CNST_TEMP
  c = 0;
  for (int n=1; n<=C.NoBC; n++) {
    if ( cmp[n].isHsrc() ) c++;
  }
  if ( c>0 ) C.EnsCompo.hsrc = ON;
  
  // 周期境界 > PERIODIC
  c = 0;
  for (int n=1; n<=C.NoBC; n++) {
    if ( cmp[n].getType() == PERIODIC ) c++;
  }
  if ( c>0 ) C.EnsCompo.periodic = ON;
  
  // 流出境界 > OUTFLOW
  c = 0;
  for (int n=1; n<=C.NoBC; n++) {
    if ( cmp[n].getType() == OUTFLOW ) c++;
  }
  if ( c>0 ) C.EnsCompo.outflow = ON;
  
  // 体積率コンポーネント
  c = 0;
  for (int n=1; n<=C.NoBC; n++) {
    if ( cmp[n].isVFraction() ) c++;
  }
  if ( c>0 ) C.EnsCompo.fraction = ON;
  
  // モニタ
  c = 0;
  for (int n=1; n<=C.NoBC; n++) {
    if ( cmp[n].isMONITOR() ) c++;
  }
  // MONITOR_LISTでCELL_MONITORが指定されている場合，C.EnsCompo.monitor==ON
  if ( (C.isMonitor() == ON) && (c < 1) ) {
    Hostonly_ stamped_printf("\tError : Cell_Monitor in MONITOR_LIST is specified, however any MONITOR can not be found.\n");
    Exit(0);
  }
  if ( (C.isMonitor() == OFF) && (c > 0) ) {
    Hostonly_ stamped_printf("\tError : Cell_Monitor in MONITOR_LIST is NOT specified, however MONITOR section is found in InnerBoundary.\n");
    Exit(0);
  }

}

/**
 @fn void SklSolverCBC::setIDtables(ParseBC* B, FILE* fp, FILE* mp)
 @brief XMLから境界条件数やID情報を取得し，表示する
 @param B
 @param fp
 @param mp 
 */
void SklSolverCBC::setIDtables(ParseBC* B, FILE* fp, FILE* mp)
{
  // C.NoBC, C.NoID, C.NoCompoを取得
  // NoID = scanXMLmodel();
  // NoCompo = NoBC + NoID;
	// ParseBCクラス内でiTable[NoID+1]を確保
  B->setControlVars(&C);
  
  // XMLファイルのModel_SettingからボクセルIDの情報を取得
  B->getXML_Model();
  
  // XMLから得られたIDテーブルを表示
  Hostonly_ {
    B->printTable(fp);
    fprintf(fp,"\n"); fflush(fp);
    
    B->printTable(mp);
    fprintf(mp,"\n"); fflush(mp);
  }
  
  // 流体と固体の媒質数，NoMaterialをセットする
  B->setMedium(&C);
}


/**
 @fn void SklSolverCBC::setGlobalCmpIdx(void)
 @brief コンポーネントのローカルなBbox情報からグローバルなBbox情報を求める
 */
void SklSolverCBC::setGlobalCmpIdx(void)
{
  SklParaComponent* para_cmp = SklGetParaComponent();
  SklParaManager* para_mng = para_cmp->GetParaManager();
  
  int st_i, st_j, st_k, ed_i, ed_j, ed_k;
  int node_st_i, node_st_j, node_st_k;
  int st[3], ed[3];
  
  // グローバルインデクスの配列インスタンス
  compo_global_bbox = new int[6*(C.NoCompo+1)];
  int* cgb = compo_global_bbox;
  
  // ローカルインデクスからグローバルインデクスに変換
  for (unsigned m=1; m<=C.NoCompo; m++) {
    
    if ( !cmp[m].isEns() ) { // コンポーネントが存在しないノードはゼロを代入
      cgb[6*m+0] = 0;
      cgb[6*m+1] = 0;
      cgb[6*m+2] = 0;
      cgb[6*m+3] = 0;
      cgb[6*m+4] = 0;
      cgb[6*m+5] = 0;
    }
    else { // コンポーネントが存在する場合
      cmp[m].getBbox(st, ed);
      st_i = st[0];
      st_j = st[1];
      st_k = st[2];
      ed_i = ed[0];
      ed_j = ed[1];
      ed_k = ed[2];
      
      if( para_mng->IsParallel() ){
        node_st_i = para_mng->GetVoxelHeadIndex(pn.ID, 0);
        node_st_j = para_mng->GetVoxelHeadIndex(pn.ID, 1);
        node_st_k = para_mng->GetVoxelHeadIndex(pn.ID, 2);
        
        cgb[6*m+0] = node_st_i + st_i;
        cgb[6*m+1] = node_st_j + st_j;
        cgb[6*m+2] = node_st_k + st_k;
        cgb[6*m+3] = node_st_i + ed_i;
        cgb[6*m+4] = node_st_j + ed_j;
        cgb[6*m+5] = node_st_k + ed_k;
      } 
      else {
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
  unsigned array_size = 6*(C.NoCompo+1);
  int np = para_mng->GetNodeNum(pn.procGrp);
  int st_x, st_y, st_z, ed_x, ed_y, ed_z, es;
  
  if ( !(m_gArray = new int[np*6]) ) Exit(0);
  if ( !(m_eArray = new int[np]  ) ) Exit(0);
  
  for (unsigned n=1; n<=C.NoBC; n++) {
    if (para_mng->IsParallel()) {
      es = ( cmp[n].isEns() ) ? 1 : 0;
      if (!para_mng->Gather(&es, 1, SKL_ARRAY_DTYPE_INT, m_eArray, 1, SKL_ARRAY_DTYPE_INT, 0, pn.procGrp)) Exit(0);
      if (!para_mng->Gather(&cgb[6*n], 6, SKL_ARRAY_DTYPE_INT, m_gArray, 6, SKL_ARRAY_DTYPE_INT, 0, pn.procGrp)) Exit(0);
      
      if (pn.ID == 0) { // マスターノードのみ
        
        // 初期値
        cgb[6*n+0] = 100000000;
        cgb[6*n+1] = 100000000;
        cgb[6*n+2] = 100000000;
        cgb[6*n+3] = 0;
        cgb[6*n+4] = 0;
        cgb[6*n+5] = 0;
        
        for (int m=0; m<np; m++) {
          if ( m_eArray[m]==1 ) { // コンポーネントの存在ランクのみを対象とする
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
  if (m_gArray) {
    delete[] m_gArray;
    m_gArray = NULL;
  }
  if (m_eArray) {
    delete[] m_eArray;
    m_eArray = NULL;
  }
}


/**
 @fn void SklSolverCBC::setLocalCmpIdx_Binary(void)
 @brief midの情報から各BCコンポーネントのローカルなインデクスを取得する
 @note
 - 計算内部領域の境界と外部境界とでは，ガイドセル部分にあるコンポーネントIDの取り扱いが異なる
 - 外部境界に接する面では，幅はそのまま，始点はガイドセル部分を含む
 - 内部境界に接する面では，始点と幅はローカルノード内の計算内部領域に含まれるように調整
 */
void SklSolverCBC::setLocalCmpIdx_Binary(void)
{
  unsigned st_i[3], len[3];
  int id;
  int m_st[3], m_ed[3];
  
  for (unsigned m=1; m<=C.NoBC; m++) {
    id = (int)cmp[m].getID();
    
    switch ( cmp[m].getType() ) {
      case HEX:
      case FAN:
        break; // 体積率でエンコード済み
        
      default:
        for (int d=0; d<3; d++) {
          st_i[d] = 0;
          len[d] = 0;
        }
        
        // コンポーネント範囲
        // GetBndIndexExtGc()は自ノード内でのidのバウンディングボックスを取得．インデクスはローカルインデクスで，ガイドセルを含む配列の基点をゼロとするCのインデクス
        if ( !dc_mid->GetBndIndexExtGc(id, st_i[0], st_i[1], st_i[2], len[0], len[1], len[2], 0) ) {
          Hostonly_ stamped_printf("\tError : can not get component local index for ID[%d]\n", id);
          Exit(0);
        }
        
        // ノード内にコンポーネントがあるかどうかをチェック
        if ( (len[0]==0) || (len[1]==0) || (len[2]==0) ) { // コンポーネントなし
          cmp[m].setEns(OFF);
          // BV情報はCompoListのコンストラクタでゼロに初期化されているので，すべてゼロ
        }
        else {
          
          for (unsigned d=0; d<3; d++) {
            int tmp_st=0;
            int tmp_ed=0;
            
            EnlargeIndex(tmp_st, tmp_ed, st_i[d], len[d], size[d], d, id);
            
            m_st[d] = tmp_st;
            m_ed[d] = tmp_ed;
          }
          cmp[m].setBbox(m_st, m_ed);
          cmp[m].setEns(ON); // コンポーネントあり
        }
        
        break;
    }
    
  }
}


/**
 @fn void SklSolverCBC::set_Parallel_Info(void)
 @brief 並列処理時のプロセッサの隣接ランクを計算する
 */
void SklSolverCBC::set_Parallel_Info(void)
{
  SklParaComponent* para_cmp = SklGetParaComponent();
  SklParaManager* para_mng = para_cmp->GetParaManager();
  int nID[6], sidx[3];
  
  if( para_mng->IsParallel() ){
    if (para_mng->IsEv()) {
      nID[X_MINUS] = para_mng->GetCommID(-1, 0, 0);
      nID[Y_MINUS] = para_mng->GetCommID( 0,-1, 0);
      nID[Z_MINUS] = para_mng->GetCommID( 0, 0,-1);
      nID[X_PLUS]  = para_mng->GetCommID( 1, 0, 0);
      nID[Y_PLUS]  = para_mng->GetCommID( 0, 1, 0);
      nID[Z_PLUS]  = para_mng->GetCommID( 0, 0, 1);
    }
    else if (para_mng->IsMb()) {
      nID[X_MINUS] = para_mng->IsCommID(-1, 0, 0);
      nID[Y_MINUS] = para_mng->IsCommID( 0,-1, 0);
      nID[Z_MINUS] = para_mng->IsCommID( 0, 0,-1);
      nID[X_PLUS]  = para_mng->IsCommID( 1, 0, 0);
      nID[Y_PLUS]  = para_mng->IsCommID( 0, 1, 0);
      nID[Z_PLUS]  = para_mng->IsCommID( 0, 0, 1);
    }
    else {
      stamped_printf("\tID %d : not parallel process.\n", pn.ID);
      Exit(0);
    }
    
    for(int i=0; i<6; i++) pn.nID[i] = nID[i];
  }
  
  if( para_mng->IsParallel() ){
    sidx[0] = para_mng->GetVoxelHeadIndex(pn.ID, 0, pn.procGrp) + 1;
    sidx[1] = para_mng->GetVoxelHeadIndex(pn.ID, 1, pn.procGrp) + 1;
    sidx[2] = para_mng->GetVoxelHeadIndex(pn.ID, 2, pn.procGrp) + 1;
  } 
  else {
    sidx[0] = sidx[1] = sidx[2] = 1;
  }
  
  for(int i=0; i<3; i++) pn.st_idx[i] = sidx[i];
  
}

/**
 @fn void SklSolverCBC::setParallelism(void)
 @brief 並列化と分割の方法を保持
 */
void SklSolverCBC::setParallelism(void)
{
  int m_np; ///< 並列ノード数
  MPI_Comm_size(MPI_COMM_WORLD, &m_np);
  
  C.num_thread  = omp_get_max_threads();
  
  // Serial or Parallel environment
  if( m_np > 1 ){
    C.num_process = m_np;

    if ( C.num_thread > 1 ) {
      C.Parallelism = Control::Hybrid;
    }
    else {
      C.Parallelism = Control::FlatMPI;
    }
  }
  else {
    C.num_process = 1;
    
    if ( C.num_thread > 1 ) {
      C.Parallelism = Control::OpenMP;
    }
    else {
      C.Parallelism = Control::Serial;
    }
  }

  
  // 空間分割
  C.Partition = Control::Equal;
  
}

/**
 @fn void SklSolverCBC::setShapeMonitor(int* mid)
 @brief Cell_Monitorで指定するIDでモニタ部分を指定するためのしかけ
 @param mid ID配列
 */
void SklSolverCBC::setShapeMonitor(int* mid)
{
  int f_st[3], f_ed[3];
  
  float ctr[3];
  float nv[3];
  float dr[3];
  float m_depth;
  float m_radius;
  float m_width;
  float m_height;
  
  // ShapeMonitorのインスタンス
  ShapeMonitor SM(size, guide, C.dx, C.org);
  
  for (int n=1; n<=C.NoBC; n++) {
    
    if (cmp[n].isMONITOR()) {

      // 形状パラメータのセット
      nv[0]   = cmp[n].nv[0];
      nv[1]   = cmp[n].nv[1];
      nv[2]   = cmp[n].nv[2];
      
      dr[0]   = cmp[n].dr[0];
      dr[1]   = cmp[n].dr[1];
      dr[2]   = cmp[n].dr[2];
      
      ctr[0]  = cmp[n].oc[0]/C.RefLength;
      ctr[1]  = cmp[n].oc[1]/C.RefLength;
      ctr[2]  = cmp[n].oc[2]/C.RefLength;
      
      m_depth = cmp[n].depth/C.RefLength;
      m_radius= cmp[n].shp_p1/C.RefLength;
      m_width = cmp[n].shp_p1/C.RefLength;
      m_height= cmp[n].shp_p2/C.RefLength;
      
      switch ( cmp[n].get_Shape() ) {
        case SHAPE_BOX:
          SM.setShapeParam(nv, ctr, dr, m_depth, m_width, m_height);
          break;
          
        case SHAPE_CYLINDER:
          SM.setShapeParam(nv, ctr, m_depth, m_radius);
          break;
          
        default:
          Exit(0);
          break;
      }
      
      // 回転角度の計算
      SM.get_angle(); 
      
      // bboxと投影面積の計算
      cmp[n].area = SM.get_BboxArea();
      
      // インデクスの計算 > あとで，VoxEncode()でresize
      SM.bbox_index(f_st, f_ed);
      
      // インデクスのサイズ登録と存在フラグ
      cmp[n].setBbox(f_st, f_ed);
      cmp[n].setEns(ON);
      
      int id = (int)cmp[n].getID();
      SM.setID(f_st, f_ed, mid, id);
      
    }
  }
  
}


/**
 @fn void SklSolverCBC::setup_Polygon2CutInfo(unsigned long& m_prep, unsigned long& m_total, FILE* fp)
 @brief Polylibを準備し，ポリゴンをロードする
 @param m_prep  前処理用のメモリサイズ
 @param m_total 本計算用のメモリリサイズ
 @param fp
 @note Polylib: 並列計算領域情報　ポリゴンは実スケールで，ガイドセル領域部分も含めて指定する
 */
void SklSolverCBC::setup_Polygon2CutInfo(unsigned long& m_prep, unsigned long& m_total, FILE* fp)
{
  //SklParaComponent* para_cmp = SklGetParaComponent();
  //SklParaManager* para_mng = para_cmp->GetParaManager();
  
  unsigned poly_gc[3];
  float    poly_org[3], poly_dx[3];
  poly_gc[0]  = poly_gc[1] = poly_gc[2] = guide;
  
  // 有次元スケールで処理する
  poly_dx[0]  = (float)C.dx[0]*C.RefLength;
  poly_dx[1]  = (float)C.dx[1]*C.RefLength;
  poly_dx[2]  = (float)C.dx[2]*C.RefLength;
  poly_org[0] = (float)C.org[0]*C.RefLength;
  poly_org[1] = (float)C.org[1]*C.RefLength;
  poly_org[2] = (float)C.org[2]*C.RefLength;
  
  Hostonly_ {
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
    fprintf(fp,"\t>> Polylib configuration\n\n");
    fprintf(fp,"\tfile name = %s\n\n", C.PolylibConfigName);
    fprintf(mp,"\n---------------------------------------------------------------------------\n\n");
    fprintf(mp,"\t>> Polylib configuration\n\n");
    fprintf(mp,"\tfile name = %s\n\n", C.PolylibConfigName);
  }
  
  // Polylib: インスタンス取得
  PL = MPIPolylib::get_instance();
  
  // Polylib: ポリゴングループのFactoryクラスを登録
  //PL->set_factory( new MyPolygonGroupFactory() );
  
  // Polylib: 並列計算領域情報を設定
  poly_stat = PL->init_parallel_info( MPI_COMM_WORLD,
                                     poly_org,   //myParaInfos[rank].bpos,
                                     size,       //myParaInfos[rank].bbsize,
                                     poly_gc,    //myParaInfos[rank].gcsize,
                                     poly_dx     //myParaInfos[rank].dx
                                     );
  if( poly_stat != PLSTAT_OK ) {
    fprintf(mp,"\tRank [%d]: p_polylib->init_parallel_info() failed.", pn.ID);
    Exit(0);
  }
  
  // Polylib: STLデータ読み込み
  TIMING_start(tm_polygon_load);
  poly_stat = PL->load_rank0( C.PolylibConfigName );
  if( poly_stat != PLSTAT_OK ) {
    fprintf(mp,"\tRank [%d]: p_polylib->load_rank0() failed.", pn.ID);
    Exit(0);
  }
  TIMING_stop(tm_polygon_load);
  
  // 階層情報表示 debug
#if 0
  PL->show_group_hierarchy();
  PL->show_group_hierarchy(fp);
#endif
  
  // IDの表示
  vector<PolygonGroup*>* pg_roots = PL->get_root_groups();
  vector<PolygonGroup*>::iterator it;
  
  Hostonly_ printf("\t  ID : Polygon Group\n");
  
  for (it = pg_roots->begin(); it != pg_roots->end(); it++) {
    std::string m_pg = (*it)->get_name();
    int m_id = (*it)->get_id();
    Hostonly_ printf("\t %3d : %s\n", m_id, m_pg.c_str());
#if 0
    PL->show_group_info(m_pg); debug
#endif
  }
  
  delete pg_roots;
  
  Hostonly_ {
    printf("\n");
    fprintf(fp, "\n");
  }
  
  // 使用メモリ量　基本クラスのみ
  unsigned long poly_mem, G_poly_mem;
  G_poly_mem = poly_mem = (unsigned long)PL->used_memory_size();
  m_prep += poly_mem;
  m_total+= poly_mem;
  
  int m_np;
  MPI_Comm_size(MPI_COMM_WORLD, &m_np);
  
  //if( para_cmp->IsParallel() ) {
  if ( m_np > 1 ) {
    // para_cmp->Allreduce(&poly_mem, &G_poly_mem, 1, SKL_ARRAY_DTYPE_ULONG, SKL_SUM, pn.procGrp);
    MPI_Allreduce(&poly_mem, &G_poly_mem, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  }
  
  Hostonly_  {
    FBUtility::displayMemory("Polygon", G_poly_mem, poly_mem, fp, mp);
  }
  
  // Triangle display
#if 0
   PolylibNS::Vec3f m_min, m_max, t1(poly_org), t2(poly_dx), t3;
   t3.assign((float)size[0]*t2.t[0], (float)size[1]*t2.t[1], (float)size[2]*t2.t[2]);
   m_min = t1 - t2;      // 1層外側まで
   m_max = t1 + t3 + t2; // 
   printf("min : %f %f %f\n", m_min.t[0], m_min.t[1], m_min.t[2]);
   printf("max : %f %f %f\n", m_max.t[0], m_max.t[1], m_max.t[2]);
   vector<Triangle*>* trias = PL->search_polygons("Tube", m_min, m_max, false); // false; ポリゴンが一部でもかかる場合
   
   PolylibNS::Vec3f *p, nrl, n;
   int c=0;
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
  
  // Polylib: STLデータ書き出しテスト
#if 0
   unsigned poly_out_para = IO_GATHER; // 逐次の場合と並列の場合で明示的に切り分けている．あとで，考慮
   string fname;
   
  if ( poly_out_para == IO_GATHER ) {
    poly_stat = PL->save_rank0( &fname, "stl_b" );
    if( poly_stat != PLSTAT_OK ) {
      Hostonly_ fprintf(mp, "Rank [%d]: p_polylib->save_rank0() failed to write into '%s'.", pn.ID, fname.c_str());
      Exit(0);
    }
  }
  else {
    poly_stat = PL->save_parallel( &fname, "stl_b" );
    if( poly_stat != PLSTAT_OK ) {
      fprintf(mp, "Rank [%d]: p_polylib->save_parallel() failed to write into '%s'.", pn.ID, fname.c_str());
      Exit(0);
    }
  }
#endif
  
  // Cutlib
  Hostonly_ {
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
    fprintf(fp,"\t>> Cut Info\n\n");
    fprintf(mp,"\n---------------------------------------------------------------------------\n\n");
    fprintf(mp,"\t>> Cut Info\n\n");
  }
  
  size_t n_cell[3];
  
  for ( unsigned i=0; i<3; i++) {
    n_cell[i] = size[i] + 2*guide;          // 分割数+ガイドセル両側
    poly_org[i] -= poly_dx[i]*(float)guide; // ガイドセル分だけシフト
  }
  
  unsigned size_n_cell = n_cell[0] * n_cell[1] * n_cell[2];
  
  
  // Cutlibの配列は各方向(引数)のサイズ
  TIMING_start(tm_init_alloc);
  cutPos = new CutPos32Array(n_cell); // 6*(float)
  cutBid = new CutBid5Array(n_cell);  // (int32_t)
  TIMING_stop(tm_init_alloc);
  
  TIMING_start(tm_cutinfo);
  CutInfoCell(poly_org, poly_dx, PL, cutPos, cutBid); // ガイドセルを含む全領域で交点を計算
  TIMING_stop(tm_cutinfo);
  
  // 使用メモリ量　
  unsigned long cut_mem, G_cut_mem;
  G_cut_mem = cut_mem = (unsigned long)size_n_cell *(6+1)*4;
  m_prep += cut_mem;
  m_total+= cut_mem;
  
  //if( para_cmp->IsParallel() ) {
  if ( m_np > 1 ) {
    //para_cmp->Allreduce(&cut_mem, &G_cut_mem, 1, SKL_ARRAY_DTYPE_ULONG, SKL_SUM, pn.procGrp);
    MPI_Allreduce(&cut_mem, &G_cut_mem, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  }
  
  Hostonly_  {
    FBUtility::displayMemory("Cut", G_cut_mem, cut_mem, fp, mp);
  }
  
  // カットとID情報をポイント
  float* cut = (float*)cutPos->getDataPointer();
  int* cut_id = (int*)cutBid->getDataPointer();
  
  // データクラスへのコピー
  allocArray_Cut(m_total);
  
  REAL_TYPE* tmp_cut = NULL;
  int* tmp_bid = NULL;
  if( !(tmp_cut  = dc_cut->GetData()) ) Exit(0);
  if( !(tmp_bid  = dc_bid->GetData()) ) Exit(0);

  int cut_len = dc_cut->GetArrayLength();
  int bid_len = dc_bid->GetArrayLength();
  fb_copy_real_ (tmp_cut, cut, &cut_len);
  fb_copy_int_  (tmp_bid, cut_id, &bid_len);
  
  // テスト
  int z=0;
  float *pos;
  float d_min=1.0;
  int bd, b0, b1, b2, b3, b4, b5;
  size_t mp, mb;
  
  for (int k=0; k<=size[2]+1; k++) {
    for (int j=0; j<=size[1]+1; j++) {
      for (int i=0; i<=size[0]+1; i++) {
        
        mp = FBUtility::getFindexS3Dcut(size, guide, 0, i, j, k);
        pos = &tmp_cut[mp];
        
        mb = FBUtility::getFindexS3D(size, guide, i, j, k);
        bd = tmp_bid[mb];
        
        if ( (pos[0]+pos[1]+pos[2]+pos[3]+pos[4]+pos[5]) < 6.0 ) { // 6方向のうちいずれかにカットがある
          
#if 0 // debug
          b0 = (bd >> 0)  & MASK_5;
          b1 = (bd >> 5)  & MASK_5;
          b2 = (bd >> 10) & MASK_5;
          b3 = (bd >> 15) & MASK_5;
          b4 = (bd >> 20) & MASK_5;
          b5 = (bd >> 25) & MASK_5;
#endif
          
          for (int n=0; n<6; n++) {
            d_min = min(d_min, pos[n]);
          }
          z++;
          
#if 0 // debug
          printf("%3d %3d %3d : %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %d %d %d %d %d %d\n", i,j,k,
                 pos[0], pos[1], pos[2], pos[3], pos[4], pos[5], b0, b1, b2, b3, b4, b5);
#endif
          
        }  
      }
    }
  }
  
  //if( para_mng->IsParallel() ) {
  if ( m_np > 1 ) {
    
    float tmp =d_min;
    MPI_Allreduce(&tmp, &d_min, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
    
    int tmp_g = z;
    VoxInfo::int_sum_Allreduce(&tmp_g, &z);
    
    unsigned tmp_u = size_n_cell;
    VoxInfo::uint_sum_Allreduce(&tmp_u, &size_n_cell);
    
  }
  
  Hostonly_ printf("\n\tMinimum dist. = %5.3e  : # of cut = %d : %f [percent]\n", d_min, z, (float)z/(float)size_n_cell*100.0);
  
  // カットの最小値
  min_distance(tmp_cut, fp);
  
}

/**
 @fn void SklSolverCBC::setup_CutInfo4IP(unsigned long& m_prep, unsigned long& m_total, FILE* fp)
 @brief IP用にカット領域をアロケートする
 @param m_prep  前処理用のメモリサイズ
 @param m_total 本計算用のメモリリサイズ
 @param fp
 */
void SklSolverCBC::setup_CutInfo4IP(unsigned long& m_prep, unsigned long& m_total, FILE* fp)
{
  SklParaComponent* para_cmp = SklGetParaComponent();
  SklParaManager* para_mng = para_cmp->GetParaManager();
  
  Hostonly_ {
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
    fprintf(fp,"\t>> Cut Info\n\n");
    fprintf(mp,"\n---------------------------------------------------------------------------\n\n");
    fprintf(mp,"\t>> Cut Info\n\n");
  }
  
  unsigned n_cell[3];
  size_t size_n_cell;
  
  for (int i=0; i<3; i++) {
    n_cell[i] = size[i] + 2*guide; // 分割数+ガイドセル
  }
  size_n_cell = n_cell[0] * n_cell[1] * n_cell[2];

  
  TIMING_start(tm_init_alloc);
  allocArray_Cut(m_total);
  TIMING_stop(tm_init_alloc);
  
  // 使用メモリ量　
  unsigned long cut_mem, G_cut_mem;
  G_cut_mem = cut_mem = (unsigned long)size_n_cell *6*4;
  m_prep += cut_mem;
  m_total+= cut_mem;
  
  if( para_cmp->IsParallel() ) {
    para_cmp->Allreduce(&cut_mem, &G_cut_mem, 1, SKL_ARRAY_DTYPE_ULONG, SKL_SUM, pn.procGrp);
  }
  
  Hostonly_  {
    FBUtility::displayMemory("Cut", G_cut_mem, cut_mem, fp, mp);
  }
  
  // 初期値のセット
  float* cut=NULL;
  cut = dc_cut->GetData();
  for (size_t i=0; i<size_n_cell*6; i++) {
    cut[i] = 1.0f;
  }
  
}


/**
 @fn void SklSolverCBC::setVOF(REAL_TYPE* vof, unsigned* bx)
 @brief VOF値を気体(0.0)と液体(1.0)で初期化
 @param vof VOF function
 @param bx BCindex ID
 */
void SklSolverCBC::setVOF(REAL_TYPE* vof, unsigned* bx)
{
  int i,j,k;
  unsigned s, odr, m0;
  
  for (k=1; k<=(int)size[2]; k++) {
    for (j=1; j<=(int)size[1]; j++) {
      for (i=1; i<=(int)size[0]; i++) {
        m0 = FBUtility::getFindexS3D(size, guide, i  , j  , k  );
        s = bx[m0];
        odr = DECODE_CMP(s);
        if ( cmp[odr].getState() == FLUID ) {
          vof[m0] = ( cmp[odr].getPhase() == GAS ) ? 0.0 : 1.0;
        }
      }
    }
  }
}


/**
 @fn unsigned SklSolverCBC::Solid_from_Cut(int* mid, float* cut, const int id)
 @brief ボクセルモデルにカット情報から得られた固体情報を転写する
 @param[in/out] mid セルID
 @param[in] cut 距離情報
 @param[in] id 固体ID 
 @retval 固体セル数
 */
unsigned SklSolverCBC::Solid_from_Cut(int* mid, float* cut, const int id)
{
  //SklParaManager* para_mng = ParaCmpo->GetParaManager();
  unsigned c=0;
  int q, m_id;
  unsigned mp, m;
  
  unsigned m_sz[3], gd;
  m_sz[0] = size[0];
  m_sz[1] = size[1];
  m_sz[2] = size[2];
  gd = guide;
  m_id= id;
  
#pragma omp parallel for firstprivate(m_sz, gd, m_id) private(q) schedule(static) reduction(+:c)
  for (int k=1; k<=(int)m_sz[2]; k++) {
    for (int j=1; j<=(int)m_sz[1]; j++) {
      for (int i=1; i<=(int)m_sz[0]; i++) {
        
        mp= FBUtility::getFindexS3D(m_sz, gd, i, j, k);
        m = FBUtility::getFindexS3Dcut(m_sz, gd, 0, i, j, k);
        q = 0;
        
        // セル内に交点があれば，壁
        if ( cut[m+0] < 1.0f ) q++;
        if ( cut[m+1] < 1.0f ) q++;
        if ( cut[m+2] < 1.0f ) q++;
        if ( cut[m+3] < 1.0f ) q++;
        if ( cut[m+4] < 1.0f ) q++;
        if ( cut[m+5] < 1.0f ) q++;
        
        if ( q > 0 ) {
          mid[mp] = m_id;
          c++;
        }
      }
    }
  }
  
  int m_np;
  MPI_Comm_size(MPI_COMM_WORLD, &m_np);
  
  if ( m_np > 1 ) {
    unsigned tmp = c;
    VoxInfo::uint_sum_Allreduce(&tmp, &c);
  }
  
  return c;
}


/**
 @fn void SklSolverCBC::VoxEncode(VoxInfo* Vinfo, ParseMat* M, int* mid, float* vf)
 @brief BCIndexにビット情報をエンコードする
 @param Vinfo
 @param M
 @param mid Voxel IDの配列
 @param vf コンポーネントの体積率 
 @note bcdに対して，共通の処理を行い，それをbcp, bcv, bchにコピーし，その後個別処理を行う．
 */
void SklSolverCBC::VoxEncode(VoxInfo* Vinfo, ParseMat* M, int* mid, float* vf)
{
  unsigned  *bcv, *bh1, *bh2, *bcp, *bcd;
  bcv = bh1 = bh2 = bcp = bcd = NULL;
  
  if( !(bcd = dc_bcd->GetData()) )  Exit(0);
  if( !(bcp = dc_bcp->GetData()) )  Exit(0);
  if( !(bcv = dc_bcv->GetData()) )  Exit(0);
  if ( C.isHeatProblem() ) {
    if( !(bh1 = dc_bh1->GetData()) )  Exit(0);
    if( !(bh2 = dc_bh2->GetData()) )  Exit(0);
  }

  // 基本ビット情報（Active, State, コンポ，媒質情報）を全領域についてエンコードする
  Vinfo->setBCIndex_base1(bcd, mid, vf);
  
  // bcdの同期
  if( !dc_bcd->CommBndCell(guide) ) Exit(0);

  // C.Acell > ノードローカルの有効セル数　
  // G_Acell > グローバルなアクティブセル数
  Vinfo->setBCIndex_base2(bcd, mid, &BC, C.Acell, G_Acell, C.KindOfSolver);
  
  // STATEとACTIVEビットのコピー
  Vinfo->copyBCIbase(bcp, bcd);
  Vinfo->copyBCIbase(bcv, bcd);
  if ( C.isHeatProblem() ) {
    Vinfo->copyBCIbase(bh1, bcd);
    Vinfo->copyBCIbase(bh2, bcd);
  }
  
  // BCIndexP に圧力計算のビット情報をエンコードする -----
  if ( C.isCDS() ) {
    C.NoWallSurface = Vinfo->setBCIndexP(bcd, bcp, mid, &BC, C.isCDS(), dc_cut->GetData());
  }
  else { // binary
    C.NoWallSurface = Vinfo->setBCIndexP(bcd, bcp, mid, &BC);
  }
  
#if 0
  Vinfo->dbg_chkBCIndexP(bcd, bcp, "BCindexP.txt");
#endif
  
  // BCIndexV に速度計算のビット情報をエンコードする -----
  if ( C.isCDS() ) {
    Vinfo->setBCIndexV(bcv, mid, &BC, bcp, C.isCDS(), dc_cut->GetData(), dc_bid->GetData());
  }
  else { // binary
    Vinfo->setBCIndexV(bcv, mid, &BC, bcp);
  }
  
#if 0
  Vinfo->dbg_chkBCIndexV(bcv, "BCindexV.txt");
#endif
  
  // BCIndexT に温度計算のビット情報をエンコードする -----
  if ( C.isHeatProblem() ) {
    Vinfo->setBCIndexH(bcd, bh1, bh2, mid, &BC, C.KindOfSolver);
#if 0
    Vinfo->dbg_chkBCIndexH(bcv, "BCindexH.txt");
#endif
  }

  // 内部領域のFluid, Solidのセル数を数える C.Wcell(Local), G_Wcell(global)
  Vinfo->countCellState(C.Wcell, G_Wcell, bcd, SOLID);
  Vinfo->countCellState(C.Fcell, G_Fcell, bcd, FLUID);
  
  // set local active cell ratio
  C.Eff_Cell_Ratio = (REAL_TYPE)C.Acell / C.getCellSize(size);
  
  // getLocalCmpIdx()などで作成したコンポーネントのインデクスの再構築
  resizeCompoBV(bcd, bcv, bh1, bh2, C.KindOfSolver, C.isHeatProblem());

}


//@fn void SklSolverCBC::VoxelInitialize(void)
//@brief 計算領域全体のサイズ，並列計算時のローカルのサイズ，コンポーネントのサイズなどを設定する
//@note パラメータScaling_Factorで強制的にスケールを変更可能
void SklSolverCBC::VoxelInitialize(void)
{
  SklParaComponent* para_cmp = SklGetParaComponent();
  SklParaManager* para_mng = para_cmp->GetParaManager();
  
  unsigned m_sz[3];
  REAL_TYPE m_org[3], m_pch[3], m_wth[3];
  REAL_TYPE scaling = C.Hide.Scaling_Factor;
  
  for (int i=0; i<3; i++) {
    m_org[i] = m_pch[i] = m_wth[i] =0.0;
    m_sz[i] = 0;
  }
  
  // ユーザ問題の場合
  if ( C.Mode.Example == id_Users ) {
    const char* fname=NULL;
    
    if ( !(fname = C.getVoxelFileName()) ) {
      Hostonly_ stamped_printf("\tRead error : A file name of voxel model is invalid.\n");
      Exit(0);
    }
    
    if ( C.vxFormat == Control::Sphere_SVX ) {
      unsigned type;
      float f_org[3], f_pch[3];
      
      if ( !getSVXHeaderInfo(fname, &type, m_sz, f_org, f_pch) ) {
        Hostonly_ stamped_printf("\tRead error : Header of the voxel file '%s' could not read.\n", fname);
        Exit(0);
      }
      
      for (int i=0; i<3; i++) {
        m_org[i] = (REAL_TYPE)f_org[i]*scaling;
        m_pch[i] = (REAL_TYPE)f_pch[i]*scaling;
      }
    }
    else { // SBX
      unsigned dtype, dims, vlen, gcell, rlen, crddef, aux;
      unsigned long long blksz, l_sz[3];
      double d_org[3], d_pch[3];
      
      if ( !getSBXHeaderInfo(fname, &dims, &vlen, &dtype, &gcell, &rlen, &crddef, &aux, &blksz, l_sz, d_org, d_pch) ) {
        Hostonly_ stamped_printf("\tRead error : Header of the voxel file '%s' could not read.\n", fname);
        Exit(0);
      }
      for (int i=0; i<3; i++) {
        m_sz[i]  = (unsigned)l_sz[i];
        m_org[i] = (REAL_TYPE)d_org[i]*scaling;
        m_pch[i] = (REAL_TYPE)d_pch[i]*scaling;
      }
    }
    
    // ファイルに記述されたヘッダは有次元であるため，無次元化する
    for (int i=0; i<3; i++) {
      m_org[i] /= C.RefLength;
      m_pch[i] /= C.RefLength;
    }
    
    // 全計算領域のbounding boxサイズ
    for (int i=0; i<3; i++) {
      m_wth[i] = (REAL_TYPE)m_sz[i] * m_pch[i];
    }
  }
  else { // その他の組み込み例題の場合
    
    // 分割数，基点，ピッチを取得する
    if ( !SklUtil::getCellInfo(C.NoDimension, m_sz, m_org, m_pch, m_wth) ) Exit(0);
    
    // 有次元の場合に無次元化
    if (C.Unit.Param == DIMENSIONAL ) {
      for (int i=0; i<3; i++) {
        m_org[i] /= C.RefLength;
        m_pch[i] /= C.RefLength;
        m_wth[i] /= C.RefLength;
      }
    }
    
    // sacaling もし指定があれば
    for (int i=0; i<3; i++) {
      m_org[i] *= scaling;
      m_pch[i] *= scaling;
      m_wth[i] *= scaling;
    }
    
    // 各例題のパラメータ設定
    Ex->setDomain(&C, m_sz, m_org, m_wth, m_pch);
  }
  
  // 以下は，共通処理
  int idiv, jdiv, kdiv;
  idiv = jdiv = kdiv = 0;
  
  // グローバルな値の保持 （無次元値）
  for (int i=0; i<3; i++) {
    G_size[i] = m_sz[i];
    G_org[i]  = m_org[i];
    G_Lbx[i]  = m_wth[i];
  }
  
  // 分割数をXMLから取得，指定なければ自動分割する
  if( para_mng->IsParallel() ){
    if( !GetCfgVoxelDivisionMethod(idiv, jdiv, kdiv) ) idiv = jdiv = kdiv = 0;
  }
  else {
    idiv = jdiv = kdiv = 1;
  }

  
  // 分割数を計算し，sz_paraで取得
  if( SKL_PARACMPO_SUCCESS != para_cmp->SklVoxelInit(m_sz[0], m_sz[1], m_sz[2], idiv, jdiv, kdiv, pn.procGrp)) {
    stamped_printf("\tID %d : Voxel Initialize error.\n", pn.ID);
    Exit(0);
  }
  
  const unsigned int* sz_para = para_mng->GetVoxelSize();
  if( !sz_para ){
    stamped_printf("\tID %d : Can't get voxel size.\n", pn.ID);
    Exit(0);
  }
  
  // 分割数を保持
  if( para_mng->IsParallel() ){
    const unsigned int* para_div = para_mng->GetVoxelDivInfo(pn.procGrp);
    if ( !para_div ) {
      stamped_printf("\tError : division size\n");
    }
    num_div_domain[0] = para_div[0];
    num_div_domain[1] = para_div[1];
    num_div_domain[2] = para_div[2];
  }
  
  
  // 並列処理時のランク情報と各ノードのスタートインデクスを計算
  set_Parallel_Info();
  
  // ノードローカルな値の設定　（無次元値）
  for(int i=0; i<3; i++) {
    m_sz[i] = sz_para[i];
    m_wth[i] = (REAL_TYPE)m_sz[i]*m_pch[i];
    m_org[i] += m_pch[i]*(REAL_TYPE)(pn.st_idx[i]-1); // セル（1，1，1）の小さい方のコーナー座標
  }
  
  //stamped_printf("Local : org(%e %e %e)\n", m_org[0], m_org[1], m_org[2]);
  
  // コントロールクラスのメンバ変数で値を保持
  C.setDomainInfo(m_sz, m_org, m_pch, m_wth);
}


/**
 @fn void SklSolverCBC::VoxScan(VoxInfo* Vinfo, ParseBC* B, int* mid, FILE* fp)
 @brief ボクセルをスキャンし情報を表示する
 @param Vinfo
 @param B
 @param mid 
 @param fp 
 @note
 - ボクセルデータに含まれるID数をカウント
 - XMLに記述されたパラメータと比較
 */
void SklSolverCBC::VoxScan(VoxInfo* Vinfo, ParseBC* B, int* mid, FILE* fp)
{
  // 外部境界面の媒質IDとその個数を取得
  int cell_id[NOFACE];
  for (int i=0; i<NOFACE; i++) cell_id[i] = 0;
  
  B->count_Outer_Cell_ID(cell_id);
  
  Hostonly_ {
    fprintf(fp, "\tCell IDs on Guide cell region\n");
    for ( int i=0; i<NOFACE; i++) {
      fprintf(fp, "\t\t%s = %d\n", FBUtility::getDirection(i).c_str(), cell_id[i]);
    }
    fprintf(mp, "\tCell IDs on Guide cell region\n");
    for ( int i=0; i<NOFACE; i++) {
      fprintf(mp, "\t\t%s = %d\n", FBUtility::getDirection(i).c_str(), cell_id[i]);
    }
  }
  
  // midにロードされたIDをスキャンし，IDの個数を返し，作業用のcolorList配列にIDを保持，midに含まれるIDの数をチェック
  unsigned sc=0;
  if ( (sc=Vinfo->scanCell(mid, cell_id, C.Hide.Change_ID)) > C.NoID ) {
    Hostonly_ stamped_printf("A number of IDs included in voxel model(%d) is grater than one described in 'Model_Setting'(%d)\n", 
                             sc, C.NoID);
    Exit(0);
  }
}

/**
 @fn void SklSolverCBC::write_distance(float* cut)
 @brief 距離の最小値を求め，閾値以上にする
 @retval 最小距離
 @param cut カット情報
 */
void SklSolverCBC::write_distance(float* cut)
{
  
  float* tmp =NULL;
  tmp = new float [(size[0]+2*guide)*(size[1]+2*guide)*(size[2]+2*guide)];
  
  size_t m, n;
  float s;
  
  for (int k=1; k<=size[2]; k++) { 
    for (int j=1; j<=size[1]; j++) {
      for (int i=1; i<=size[0]; i++) {
        n = FBUtility::getFindexS3D(size, guide, i, j, k);
        
        s = 1.0e6;
        for (int l=0; l<6; l++) {
          m = FBUtility::getFindexS3Dcut(size, guide, l, i, j, k);
          if ( cut[m] < s ) s = cut[m];
        }
        tmp[n] = s;
      }
    }
  }
  
  REAL_TYPE org[3], pit[3];
  
  
  //  ガイドセルがある場合(GuideOut != 0)にオリジナルポイントを調整
  for (int i=0; i<3; i++) {
    org[i] = C.org[i] - C.dx[i]*(REAL_TYPE)C.GuideOut;
    pit[i] = C.dx[i];
  }
  
  // 出力ファイルの指定が有次元の場合
  if ( C.Unit.File == DIMENSIONAL ) {
    for (int i=0; i<3; i++) {
      org[i] *= C.RefLength;
      pit[i] *= C.RefLength;
    }
  }
  F.writeRawSPH(tmp, size, guide, org, pit, SPH_SINGLE);
}


//@fn int SklSolverCBC::get_intval( std::string& buffer )
//@brief XMLファイルから値をとりだす
int SklSolverCBC::get_intval( std::string& buffer )
{
	int s = buffer.find( "value=\"", 0 ) + 7;
	int e = buffer.find( "\"", s );
	return atoi( buffer.substr( s, e-s ).c_str() );
}


//@fn std::string SklSolverCBC::get_strval( std::string& buffer )
//@brief XMLファイルから値をとりだす
std::string SklSolverCBC::get_strval( std::string& buffer )
{
	int s = buffer.find( "value=\"", 0 ) + 7;
	int e = buffer.find( "\"", s );
  std::string result = buffer.substr( s, e-s );
	return result;
}


//@fn bool SklSolverCBC::getCoarseResult()
//@note 2倍密格子の領域開始インデクス番号から、その領域が属する粗格子計算結果ファイル名と、その計算結果ファイルの開始インデクス番号を取得する
bool SklSolverCBC::getCoarseResult (
                                    int i,		                     // (in) 密格子　開始インデクスi
                                    int j,		                     // (in) 同j
                                    int k,		                     // (in) 同k
                                    std::string& coarse_dfi_fname, // (in) 粗格子のdfiファイル名（どの変数のものでも良い）
                                    std::string& coarse_prefix,	   // (in) 粗格子計算結果ファイルプリフィクス e.g. "prs_16"
                                    const int m_step,              // (in) 探索するステップ数
                                    std::string& coarse_sph_fname, // (out) ijk位置の結果を含む粗格子計算結果ファイル名
                                    int* c_size,                   // (out) 粗格子の分割数
                                    int* coarse,  	               // (out) 粗格子　開始インデクス
                                    int* block                     // (out) 含まれるブロック数
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
  
  //int mm;
  //MPI_Comm_size(MPI_COMM_WORLD, &mm); printf("nnum = %d\n", np);
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
  //       pn.ID,  i,j,k,  i0,j0,k0,  bh_i,bh_j,bh_k,  Ci,Cj,Ck,  Fi,Fj,Fk,  coarse_sph_fname.c_str());
  //
  
	return true;
}
