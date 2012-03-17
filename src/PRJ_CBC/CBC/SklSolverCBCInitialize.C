/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file SklSolverCBCInitialize.C
//@brief CBCソルバークラスのプリプロセッサ
//@author keno, FSI Team, VCAD, RIKEN

#include "SklSolverCBC.h"
#include "CompoFraction.h"

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
  
  ws = v = t = p = NULL;
  mid = NULL;
  bcd = bcv = bh1 = bh2 = bcp = NULL;
  TotalMemory = PrepMemory = 0;
  G_TotalMemory = G_PrepMemory = 0;
  
  // 前処理段階のみに使用するオブジェクトをインスタンス
  VoxInfo Vinfo;
  ParseBC     B;
  ParseMat    M;

  // 並列処理時のランク番号をセット．デフォルトで，procGrp = 0;
  pn.ID = para_mng->GetMyID(pn.procGrp);
  
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

  // 計算領域全体のサイズ，並列計算時のローカルのサイズ，コンポーネントのサイズなどを設定する -----------------------------------------------------
  VoxelInitialize();

  // 組み込み例題クラスへ変数の値をコピー，固有のパラメータ取得，領域設定 -----------------------------------------------------
  Ex->setControlVars (&C);
  
  // 並列情報のコピー
  B.setParallelInfo   (pn);
  M.setParallelInfo   (pn);
  BC.setParallelInfo  (pn);
  C.setParallelInfo   (pn);
  F.setParallelInfo   (pn);
  MO.setParallelInfo  (pn);
  TIMING__ PM.setParallelInfo  (pn);
  Ex->setParallelInfo (pn);  // 続く処理より先に初期化しておく

  // Serial or Parallel environment
  if( para_mng->IsParallel() ){
    int np=para_mng->GetNodeNum(pn.procGrp);
    Hostonly_ {
      fprintf(mp, "\tParallel Execution : %d\n", np);
      fprintf(fp, "\tParallel Execution : %d\n", np);
    }
  }
  else {
    Hostonly_ {
      fprintf(mp, "\tSerial Execution\n");
      fprintf(fp, "\tSerial Execution\n");
    }
  }

  // 領域情報の表示
  Hostonly_ {
    fprintf(mp,"\n---------------------------------------------------------------------------\n\n");
    fprintf(mp,"\n\t>> Global Domain Information\n\n");
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
    fprintf(fp,"\n\t>> Global Domain Information\n\n");
    C.printDomainInfo(mp, fp, G_size, G_org, G_Lbx);
  }

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
    set_timing_label();
  }
  
  // タイミング測定開始
  TIMING_start(tm_init_sct); 
  
  // 前処理に用いるデータクラスのアロケート -----------------------------------------------------
  TIMING_start(tm_init_alloc); 
  allocArray_prep(TotalMemory, PrepMemory);
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
  
  // BinaryとCut-Distanceの分岐
  if ( C.isCDS() ) { // Cut-Distance Scheme
    
#ifndef BINARY_VOXEL    
    // PolylibとCutlibのセットアップ
    if ( C.Mode.Example == id_Users ) {
      setup_Polygon2CutInfo(PrepMemory, TotalMemory, fp);
      
      // 媒質ファイルが指定された場合の処理
      TIMING_start(tm_voxel_load);
      if ( C.Mode.Medium_Spec == ON ) {
        if ( C.vxFormat == Control::Sphere_SVX ) {
          F.readSVX(this, fp, "SphereSVX", G_size, guide, dc_mid);
        }
        else {
          F.readSBX(this, fp, "SphereSBX", G_size, guide, dc_mid);
        }
      }
      else { // 媒質ファイルを使わない場合，指定された媒質番号で初期化
        fb_set_int_s_(dc_mid->GetData(), (int*)sz, (int*)&gc, (int*)mid);
      }
      TIMING_stop(tm_voxel_load);

    }
    else { // Intrinsic problem
      setup_CutInfo4IP(PrepMemory, TotalMemory, fp);
      Ex->setup(mid, &C, G_org);
    }
#endif
  }
  else { // Binary
    
    TIMING_start(tm_voxel_load);
    if ( C.Mode.Example == id_Users ) {
      if ( C.vxFormat == Control::Sphere_SVX ) {
        F.readSVX(this, fp, "SphereSVX", G_size, guide, dc_mid);
      }
      else {
        F.readSBX(this, fp, "SphereSBX", G_size, guide, dc_mid);
      }
    }
    else { // Intrinsic problem　ユーザ問題も組み込み例題クラスの一つとして実装されている
      Ex->setup(mid, &C, G_org);
    }
    TIMING_stop(tm_voxel_load);
    
  }
  
  // midのガイドセル同期
  if( !dc_mid->CommBndCell(guide) ) return -1;
  
  
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
  if ( !Vinfo.receiveCfgPtr(m_solvCfg) ) {
    Hostonly_ stamped_printf("\tError during sending an object pointer of XML tree to VoxInfo class\n");
    Exit(0);
  }
  Vinfo.setParallelInfo(pn);
  Vinfo.setControlVars(size, guide);
  
  // バイナリボクセルのスキャン
  VoxScan(&Vinfo, &B, mid, fp);
  
  // スキャンしたセルIDの情報を表示する
  Hostonly_ {
    fprintf(mp, "\n---------------------------------------------------------------------------\n\n");
    fprintf(mp,"\t>> Voxel Model Information\n\n");
    fprintf(fp, "\n---------------------------------------------------------------------------\n\n");
    fprintf(fp,"\t>> Voxel Model Information\n\n");
    
		Vinfo.printScanedCell(fp);
		fflush(fp);
    Vinfo.printScanedCell(mp);
		fflush(mp);
  }

  // ボクセルモデルのIDがXMLに記述されたIDに含まれていること
	Hostonly_ {
		if ( !Vinfo.chkIDconsistency(B.get_IDtable_Ptr(), C.NoID) ) {
			stamped_printf("\tID in between XML and scaned is not consistent\n");
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

#ifndef BINARY_VOXEL
  // CDSの場合，WALLとSYMMETRICのときに，カットを外部境界に接する内部セルに設定
  if ( C.isCDS() ) {
    Vinfo.setOBC_Cut(&BC, cutPos);
  }
#endif
  
#ifdef DEBUG
  // カット情報から壁面IDセット
  if ( C.isCDS() ) {
    unsigned zc = Vinfo.markSolid_from_Cut(mid, cutPos);
    Hostonly_ printf("\tCut cell = %d\n", zc);
  }
#endif
  
  // セルIDのノード間同期
  if( !dc_mid->CommBndCell(guide) ) return -1;
  

  // HEX/FANコンポーネントの形状情報からBboxと体積率を計算
  if ( C.isVfraction() ) {
    TIMING_start(tm_init_alloc);
    allocArray_compoVF(TotalMemory, PrepMemory);
    TIMING_stop(tm_init_alloc); 
    if( !(cvf = dc_cvf->GetData()) )  Exit(0);
    
    setComponentVF(cvf);
  }

  
  // コンポーネントのローカルインデクスをcmp.ciに保存
  getLocalCmpIdx();
  
  // コンポーネントのインデクス情報を取得
  getGlobalCmpIdx(&Vinfo);
  
  // 内部周期境界の場合のガイドセルのコピー処理
  Vinfo.adjCellID_Prdc_Inner(dc_mid);

  // 媒質数とKindOfSolverの整合性をチェックする
  if ( !C.chkMediumConsistency() ) {
    Hostonly_ stamped_printf("\tControl::chkMediumConsistency()\n");
    return -1;
  }
  
  // ParseMatクラスをセットアップし，媒質情報をXMLから読み込み，媒質リストを作成する
  Hostonly_  {
    fprintf(mp,"\n---------------------------------------------------------------------------\n\n");
    fprintf(mp,"\n\t>> Medium List\n\n");
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
    fprintf(fp,"\n\t>> Medium List\n\n");
  }
  setMaterialList(&B, &M, mp, fp);

  
  // BCIndexへのエンコード処理
  Hostonly_  {
    fprintf(mp,"\n---------------------------------------------------------------------------\n\n");
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
  }

  VoxEncode(&Vinfo, &M, mid, cvf, cutPos);

  // Select implementation of SOR iteration by Eff_Cell_Ratio
  C.select_Itr_Impl(IC);

  // コンポーネントの体積率を8bitで量子化し，圧力損失コンポの場合にはFORCING_BITをON > bcdにエンコード
  Vinfo.setCmpFraction(cmp, bcd, cvf);

#ifdef DEBUG
  // CompoListとMaterialListの関連を表示
  Hostonly_ M.printRelation(mp, fp, cmp);
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
  
  
#ifdef DEBUG
  // チェック comment out
  //FBUtility::chkGamma(C.imax, C.jmax, C.kmax, C.guide, bch);
#endif
  
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
      Vinfo.cal_Compo_Area_Normal(n, bcd, bcv, bh1, C.dh*C.RefLength, &GC_bv[n*6]);
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
  
  // パラメータの無次元表示（正規化）に必要な参照物理量の設定
  B.setRefMedium(mat, &C);
  
  // 流入出フラグを設定
  BC.set_InOut_flag();
  
#ifdef DEBUG
  // チェックのため，全計算セルのBCIndexの内容を表示する
  if ( !Vinfo.chkBCIndexP(bcd, bcp, "BCindex.txt") ) {
    Hostonly_ stamped_printf("\tVoxInfo::chkBCIndexP()\n");
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
    REAL_TYPE cr = (REAL_TYPE)G_Wcell/(REAL_TYPE)(G_size[0]*G_size[1]*G_size[2]) *100.0;
    fprintf(mp, "\tThis voxel includes %4d solid %s  [Solid cell ratio inside computational domain : %5.1f percent]\n\n", 
            C.NoMediumSolid, (C.NoMediumSolid>1) ? "IDs" : "ID", cr);
    fprintf(fp, "\tThis voxel includes %4d solid %s  [Solid cell ratio inside computational domain : %5.1f percent]\n\n", 
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
  // 機能がOFFの場合には直ちに戻る
  C.getXML_Monitor(&MO);
  
  // モニタリストが指定されている場合に，プローブ位置をID=255としてボクセルファイルに書き込んで出力する
  if (C.Sampling.log == ON  || MO.hasCellMonitor(cmp, C.NoBC)) {
    MO.write_ID(mid);
  }
  
  // History出力の指定
  C.getXML_History();
  
  // 内部境界条件として指定されたモニタ設定を登録
  MO.setInnerBoundary(cmp, C.NoBC);
    
  // 組み込み例題 or MonitorListの場合に，svxファイルを出力する．
  if ( (C.Mode.Example != id_Users) || ( (C.Mode.Example == id_Users) && (C.Sampling.log == ON  || MO.hasCellMonitor(cmp, C.NoBC)) ) ) {
    Hostonly_ printf("\n\twrite ID which includes Monitor List ID\n\n");
    
    // 性能測定モードがオフのときのみ出力
    if ( C.Hide.PM_Test == OFF ) Ex->writeSVX(mid, &C); // writeSVX(); ユーザ問題の場合には，単にtrueを返す
  }
  
  // mid[]を解放する  ---------------------------
  DataMngr.DeleteDataObj("mid");
  
  // コンポーネントの内容リストを表示する
  Hostonly_ {
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
    fprintf(fp,"\t>> Component Information\n");
    fprintf(mp,"\n---------------------------------------------------------------------------\n\n");
    fprintf(mp,"\t>> Component Information\n");
    B.printCompoInfo(mp, fp, Vinfo.get_vox_nv_ptr(), GC_bv, mat);
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
  
  // 各コンポーネントが存在するかどうかを保持しておく
  //setEnsComponent();
  
  // 各ノードの領域情報をファイル出力
  gather_DomainInfo();

  
  TIMING_stop(tm_voxel_prep_sct);
  // ここまでがボクセル準備の時間セクション

  
  // 計算に用いる配列のアロケート ----------------------------------------------------------------------------------
  TIMING_start(tm_init_alloc);
  allocArray_main (TotalMemory);
  
  allocArray_Collocate (TotalMemory);
  
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
  
  
  // リスタート 瞬時値と平均値に分けて処理　------------------
  if ( C.Start==Control::re_start) {
    Hostonly_ fprintf(mp,"\n---------------------------------------------------------------------------\n\n");
    Hostonly_ fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
    Hostonly_ fprintf(mp, "\t>> Restart from Previous Calculation Results\n\n");
    Hostonly_ fprintf(fp, "\t>> Restart from Previous Calculation Results\n\n");
    
    TIMING_start(tm_restart);
    load_Restart_file(fp); // 瞬時値のロード
    TIMING_stop(tm_restart);
    
    Hostonly_ fprintf(mp,"\n");
    Hostonly_ fprintf(fp,"\n");
  }
  else { // 初期スタートのステップ，時間を設定する
    SklSetBaseStep(0);
    SklSetBaseTime(0.0);
    
    // V00の値のセット．モードがtrueの場合はV00[0]=1.0に設定，そうでなければtmに応じた値
    if ( SklIsCheckMode() ) RF.setV00((double)SklGetTotalTime(), true);
    else                    RF.setV00((double)SklGetTotalTime());
    
    double g[4];
    RF.copyV00(g);
    for (int i=0; i<4; i++) v00[i]=(REAL_TYPE)g[i];
  }
  
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
    
  // V-Sphereに出力インターバルを通知
  C.tell_Interval_2_Sphere();
  
  // 平均値のロード
  if ( C.Start==Control::re_start) {
    TIMING_start(tm_restart);
    load_Restart_avr_file(fp);
    TIMING_stop(tm_restart);
  }
  
  // データクラスのポイント
  if( !(v  = dc_v->GetData()) )     return -1;
  if( !(p  = dc_p->GetData()) )     return -1;
  
  if ( C.isHeatProblem() ) {
    if( !(t = dc_t->GetData()) )    return -1;
  }
  
  // リスタート時，セルフェイスへの値の内挿 <v00をセットした後に処理>
  if ( C.Start==Control::re_start) {
    
    if (C.Mode.ShapeAprx == BINARY) {
      //BC.mod_Vec_CF(vf, bcv, bcd, tm, &C, v00, flop_task);
    }
    else if (C.Mode.ShapeAprx == CUT_INFO) {
      //BC.mod_Vec_CF(vf, bcv, bcd, tm, &C, v00, flop_task);
    }
    else {
      Hostonly_ stamped_printf("Invalid Shape Approximation\n");
      return -1;
    }
  }
  
	// リスタートの最大値と最小値の表示
  if ( C.Start == Control::re_start ) {
    Hostonly_ fprintf(mp, "\tNon-dimensional value\n");
    Hostonly_ fprintf(fp, "\tNon-dimensional value\n");
    REAL_TYPE f_min, f_max;
    fb_minmax_v_ (&f_min, &f_max, sz, gc, v00, v, &flop_task);
    Hostonly_ fprintf(mp, "\t\tV: min=%13.6e max=%13.6e\n", f_min, f_max);
    Hostonly_ fprintf(fp, "\t\tV: min=%13.6e max=%13.6e\n", f_min, f_max);
    
    fb_minmax_s_ (&f_min, &f_max, sz, gc, p, &flop_task);
    Hostonly_ fprintf(mp, "\t\tP: min=%13.6e max=%13.6e\n", f_min, f_max);
    Hostonly_ fprintf(fp, "\t\tP: min=%13.6e max=%13.6e\n", f_min, f_max);
    
    if ( C.isHeatProblem() ) {
      fb_minmax_s_ (&f_min, &f_max, sz, gc, t, &flop_task);
      Hostonly_ fprintf(mp, "\t\tT: min=%13.6e max=%13.6e\n", f_min, f_max);
      Hostonly_ fprintf(fp, "\t\tT: min=%13.6e max=%13.6e\n", f_min, f_max);
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
    if ( C.Sampling.log == ON || MO.hasCellMonitor(cmp, C.NoBC) ) {
      MO.printMonitorInfo(mp, C.HistoryMonitorName);
      MO.printMonitorInfo(fp, C.HistoryMonitorName);
    }
  }
  
  // ドライバ条件のチェック
  BC.checkDriver(fp);
  
  // 初期条件の条件設定
	if ( C.Start == Control::initial_start ) {
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
		fb_set_vector_(v, sz, gc, U0);
    
		// 外部境界面の流出流量と移流速度
    DomainMonitor( BC.get_OBC_Ptr(), &C, flop_task);
    
		// 外部境界面の移流速度を計算し，外部境界条件を設定
    BC.OuterVBC_Periodic(dc_v);
		BC.OuterVBC(v, v, bcv, tm, (REAL_TYPE)m_dt, &C, v00, flop_task);
    BC.InnerVBC(v, bcv, tm, v00, flop_task);
    BC.InnerVBC_Periodic(dc_v, dc_bcd);
    
		// 圧力
    REAL_TYPE ip;
    if (C.Unit.Param == DIMENSIONAL) {
      ip = FBUtility::convD2ND_P(C.iv.Pressure, C.BasePrs, C.RefDensity, C.RefVelocity, C.Unit.Prs);
    }
    else {
      ip = C.iv.Pressure;
    }

    fb_set_real_s_(dc_p->GetData(), sz, gc, &ip);
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

      fb_set_real_s_(dc_t->GetData(), sz, gc, &it);
      
      // コンポーネントの初期値
      for (unsigned m=C.NoBC+1; m<=C.NoCompo; m++) {
        BC.setInitialTemp_Compo(m, bcd, dc_t);
      }
      
			BC.OuterTBC(dc_t);
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

    
    // ユーザ例題のときに，速度の内部境界条件を設定する
    if ( C.Mode.Example == id_Users ) {
      ; // nothing
    }
    else {
      Ex->initCond(v, p);
    }
    
    // VOF
    if ( C.BasicEqs == INCMP_2PHASE ) {
      REAL_TYPE* vof=NULL;
      if( !(vof = dc_vof->GetData()) )  return -1;
      setVOF(vof, bcd);
      if ( !dc_vof->CommBndCell(guide) ) return -1;
    }
  }
  
  // 初期解およびリスタート解の同期
  if( !dc_v->CommBndCell(guide) ) return -1;
  if( !dc_p->CommBndCell(guide) ) return -1;
  if ( C.isHeatProblem() ) {
    if( !dc_t->CommBndCell(guide) ) return -1;
  }
  
  // 設定した初期値のファイル出力準備
  prepOutput();

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
  if ( (C.Sampling.log == ON) || MO.hasCellMonitor(cmp, C.NoBC) ) {
    if ( C.isHeatProblem() ) {
      MO.setDataPtrs(dc_v->GetData(), dc_p->GetData(), dc_t->GetData());
    }
    else {
      MO.setDataPtrs(dc_v->GetData(), dc_p->GetData());
    }
  }
  
  // 初期状態のファイル出力 性能測定モードのときには出力しない
	if ( (C.Hide.PM_Test == OFF) &&  (0 == SklGetTotalStep()) ) FileOutput(Control::IO_forced, flop_task);
  
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
 @fn void SklSolverCBC::allocArray_prep (unsigned long &total, unsigned long &prep)
 @brief 前処理に用いる配列のアロケーション
 @param total ソルバーに使用するメモリ量
 @param prep 前処理に使用するメモリ量
 */
void SklSolverCBC::allocArray_prep (unsigned long &total, unsigned long &prep)
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
 @fn void SklSolverCBC::allocArray_compoVF (unsigned long &total, unsigned long &prep)
 @brief コンポーネント体積率の配列のアロケーション
 @param total ソルバーに使用するメモリ量
 @param prep 前処理に使用するメモリ量
 */
void SklSolverCBC::allocArray_compoVF (unsigned long &total, unsigned long &prep)
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
 @fn void SklSolverCBC::allocArray_index3 (unsigned long &total)
 @brief Index計算に用いる配列のアロケーション
 @param total ソルバーに使用するメモリ量
 */
void SklSolverCBC::allocArray_index3 (unsigned long &total)
{
  unsigned long mc=0;
  unsigned idx_sz = C.Fcell * 3;
  
  if ( !A.alloc_Int_S1D(this, dc_index3, "index3", idx_sz, 0, 0.0, mc) ) Exit(0);
  total+= mc;
}

/**
 @fn void SklSolverCBC::allocArray_index (unsigned long &total)
 @brief Index計算に用いる配列のアロケーション
 @param total ソルバーに使用するメモリ量
 */
void SklSolverCBC::allocArray_index(unsigned long &total)
{
  unsigned long mc=0;
  unsigned idx_sz = C.Fcell;
  
  if ( !A.alloc_Uint_S1D(this, dc_index, "index", idx_sz, 0, 0, mc) ) Exit(0);
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
    if ( !SklCfgCheckInFile("AvrPressure") || // 時間平均を指定しているが，出力ファイル記述がない場合
        !SklCfgCheckInFile("AvrVelocity") ||
        (!SklCfgCheckInFile("AvrTemperature") && C.isHeatProblem()) ) {
      Hostonly_ stamped_printf     ("\tRestart mode and averaging, but there is no InFile description for an average file. \n");
      Hostonly_ stamped_fprintf(fp, "\tRestart mode and averaging, but there is no InFile description for an average file. \n");
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
  else {
    Hostonly_ stamped_printf("\tInvalid keyword is described for Exmple definition\n");
    Exit(0);
  }
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
      fprintf(fp,"%12.4e      %8.3f ", tmp_acl, 100.0*(tmp_acl-m_efv)/m_efv);
      fprintf(fp,"%12.4e      %8.3f  %s\n", 0.0, 0.0, (tmp_acl/tmp_vol>THRESHOLD_SOR_IMPLEMENTATION) ? "Mask-loop" : "Skip-loop");
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
 @fn void SklSolverCBC::getGlobalCmpIdx(VoxInfo* Vinfo)
 @brief コンポーネントのローカルなBV情報からグローバルなBV情報を求める
 */
void SklSolverCBC::getGlobalCmpIdx(VoxInfo* Vinfo)
{
  SklParaComponent* para_cmp = SklGetParaComponent();
  SklParaManager* para_mng = para_cmp->GetParaManager();
  
  int st_i, st_j, st_k, ed_i, ed_j, ed_k;
  int node_st_i, node_st_j, node_st_k;
  int st[3], ed[3];
  
  // グローバルインデクスの配列インスタンス
  GC_bv = new int[6*(C.NoCompo+1)];
  
  // ローカルインデクスからグローバルインデクスに変換
  for (unsigned m=1; m<=C.NoCompo; m++) {
    
    if ( !cmp[m].isEns() ) { // コンポーネントが存在しないノードはゼロを代入
      st_i = 0;
      st_j = 0;
      st_k = 0;
      ed_i = 0;
      ed_j = 0;
      ed_k = 0;
      
      GC_bv[6*m+0] = 0;
      GC_bv[6*m+1] = 0;
      GC_bv[6*m+2] = 0;
      GC_bv[6*m+3] = 0;
      GC_bv[6*m+4] = 0;
      GC_bv[6*m+5] = 0;
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
        
        GC_bv[6*m+0] = node_st_i + st_i;
        GC_bv[6*m+1] = node_st_j + st_j;
        GC_bv[6*m+2] = node_st_k + st_k;
        GC_bv[6*m+3] = node_st_i + ed_i;
        GC_bv[6*m+4] = node_st_j + ed_j;
        GC_bv[6*m+5] = node_st_k + ed_k;
      } 
      else {
        GC_bv[6*m+0] = st_i;
        GC_bv[6*m+1] = st_j;
        GC_bv[6*m+2] = st_k;
        GC_bv[6*m+3] = ed_i;
        GC_bv[6*m+4] = ed_j;
        GC_bv[6*m+5] = ed_k;
      }
    }
  }
  
  // コンポーネントBVの更新
  Vinfo->gatherGlobalIndex(GC_bv);
}


/**
 @fn void SklSolverCBC::getLocalCmpIdx(void)
 @brief midの情報から各BCコンポーネントのローカルなインデクスを取得する
 @note
 - 計算内部領域の境界と外部境界とでは，ガイドセル部分にあるコンポーネントIDの取り扱いが異なる
 - 外部境界に接する面では，幅はそのまま，始点はガイドセル部分を含む
 - 内部境界に接する面では，始点と幅はローカルノード内の計算内部領域に含まれるように調整
 */
void SklSolverCBC::getLocalCmpIdx(void)
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
          
          int m_st, m_ed;
          
          for (unsigned d=0; d<3; d++) {
            m_st = m_ed = 0;
            getEnlargedIndex(m_st, m_ed, st_i[d], len[d], size[d], d, id);
            cmp[m].setBbox_st(d, m_st);
            cmp[m].setBbox_ed(d, m_ed);
          }
          cmp[m].setEns(ON); // コンポーネントあり
        }
        
        break;
    }
    
  }
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
  else {
    Hostonly_ stamped_printf("\tInvalid keyword is described for Example definition\n");
    Exit(0);
  }
}


/**
 @fn void SklSolverCBC::getEnlargedIndex(int& m_st, int& m_ed, const unsigned st_i, const unsigned len, const unsigned mx, const unsigned dir, const int m_id)
 @brief 初期インデクスの情報を元に，一層拡大したインデクス値を返す
 @param m_st 拡大された開始点（Fortranインデクス）
 @param m_ed 拡大された終了点（Fortranインデクス）
 @param st_i 開始点（Cインデクス）
 @param len コンポーネントの存在長さ
 @param m_x 軸方向のサイズ
 @param dir 方向
 @param m_id キーID
 */
void SklSolverCBC::getEnlargedIndex(int& m_st, int& m_ed, const unsigned st_i, const unsigned len, const unsigned m_x, const unsigned dir, const int m_id)
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
 @fn void SklSolverCBC::load_Restart_file (FILE* fp)
 @brief リスタート時の瞬時値ファイル読み込み
 @param fp ファイルポインタ
 */
void SklSolverCBC::load_Restart_file (FILE* fp)
{
  int step;
  REAL_TYPE time;
  
  // 圧力の瞬時値　ここでタイムスタンプを得る
  F.loadSphScalar3D(this, fp, "Pressure", G_size, guide, dc_p, step, time, C.Unit.File);
  if (C.Unit.File == DIMENSIONAL) time /= C.Tscale;
  SklSetBaseStep(step);
  SklSetBaseTime(time);
  
  // v00[]に値をセット
  copyV00fromRF((double)time);
  
  // 有次元ファイルの場合，無次元に変換する
  if (C.Unit.File == DIMENSIONAL) {
    F.cnv_P_D2ND(dc_p, C.BasePrs, C.RefDensity, C.RefVelocity, C.Unit.Prs);
  }
  
  // Instantaneous Velocity fields
  F.loadSphVector3D(this, fp, "Velocity", G_size, guide, dc_v, step, time, v00, C.Unit.File);
  if (C.Unit.File == DIMENSIONAL) time /= C.Tscale;
  if ( (step != SklGetTotalStep()) || (time != (REAL_TYPE)SklGetTotalTime()) ) {
    Hostonly_ printf     ("\n\tTime stamp is different between files\n");
    Hostonly_ fprintf(fp, "\n\tTime stamp is different between files\n");
    Exit(0);
  }
  
  // convert from dimensional to non-dimensional, iff file is dimensional
  if (C.Unit.File == DIMENSIONAL) {
    F.cnv_V_D2ND(dc_v, C.RefVelocity);
  }
  
  // Instantaneous Temperature fields
  if ( C.isHeatProblem() ) {
    F.loadSphScalar3D(this, fp, "Temperature", G_size, guide, dc_t, step, time, C.Unit.File);
    if (C.Unit.File == DIMENSIONAL) time /= C.Tscale;
    if ( (step != SklGetTotalStep()) || (time != (REAL_TYPE)SklGetTotalTime()) ) {
      Hostonly_ printf     ("\n\tTime stamp is different between files\n");
      Hostonly_ fprintf(fp, "\n\tTime stamp is different between files\n");
      Exit(0);
    }
    
    // convert from dimensional to non-dimensional, iff file is dimensional
    if (C.Unit.File == DIMENSIONAL) {
      F.cnv_T_D2ND(dc_t, C.BaseTemp, C.DiffTemp, C.Unit.Temp);
    }
  }
  
  // VOF fields
  //F.loadSphScalar3D(this, fp, "VOF", G_size, guide, dc_vof, C.Unit.File);
}

/**
 @fn void SklSolverCBC::load_Restart_avr_file (FILE* fp)
 @brief リスタート時の平均値ファイル読み込み
 @param fp ファイルポインタ
 */
void SklSolverCBC::load_Restart_avr_file (FILE* fp)
{
  int step;
  REAL_TYPE time;
  
  step = SklGetBaseStep();
  time = SklGetBaseTime();
  
  if ( C.Mode.Average == OFF ) return;
  
  
  if ( C.Interval[Interval_Manager::tg_avstart].isStep() ) {
    if ( step > C.Interval[Interval_Manager::tg_avstart].getIntervalStep() ) {
      Hostonly_ printf     ("\tRestart from Previous Calculation Results of averaged field\n");
      Hostonly_ fprintf(fp, "\tRestart from Previous Calculation Results of averaged field\n");
      Hostonly_ printf     ("\tStep : base=%d current=%d total=%d\n", SklGetBaseStep(), SklGetCurrentStep(), SklGetTotalStep());
      Hostonly_ fprintf(fp, "\tStep : base=%d current=%d total=%d\n", SklGetBaseStep(), SklGetCurrentStep(), SklGetTotalStep());
    }
    else {
      return;
    }
  }
  else {
    if ( time > C.Interval[Interval_Manager::tg_avstart].getIntervalTime() ) {
      Hostonly_ printf     ("\tRestart from Previous Calculation Results of averaged field\n");
      Hostonly_ fprintf(fp, "\tRestart from Previous Calculation Results of averaged field\n");
      Hostonly_ printf     ("\tTime : base=%e[sec.]/%e[-] current=%e[-] total=%e[-]\n", SklGetBaseTime()*C.Tscale, SklGetBaseTime(), SklGetCurrentTime(), SklGetTotalTime());
      Hostonly_ fprintf(fp, "\tTime : base=%e[sec.]/%e[-] current=%e[-] total=%e[-]\n", SklGetBaseTime()*C.Tscale, SklGetBaseTime(), SklGetCurrentTime(), SklGetTotalTime());
    }
    else {
      return;
    }
  }
  
  // Pressure
  F.loadSphScalar3DAvr(this, fp, "AvrPressure", G_size, guide, dc_ap, step, time, C.Unit.File);
  if (C.Unit.File == DIMENSIONAL) time /= C.Tscale;
  SklSetAverageBaseStep(step);
  SklSetAverageBaseTime(time);
  
  // convert from dimensional to non-dimensional, iff file is dimensional
  if (C.Unit.File == DIMENSIONAL) {
    F.cnv_P_D2ND(dc_ap, C.BasePrs, C.RefDensity, C.RefVelocity, C.Unit.Prs);
  }
  
  // Velocity
  F.loadSphVector3DAvr(this, fp, "AvrVelocity", G_size, guide, dc_av, step, time, v00, C.Unit.File);
  if (C.Unit.File == DIMENSIONAL) time /= C.Tscale;
  if ( (step != SklGetAverageBaseStep()) || ((REAL_TYPE)time != SklGetAverageBaseTime()) ) {
    Hostonly_ printf     ("\n\tTime stamp is different between files\n");
    Hostonly_ fprintf(fp, "\n\tTime stamp is different between files\n");
    Exit(0);
  }
  
  // convert from dimensional to non-dimensional, iff file is dimensional
  if (C.Unit.File == DIMENSIONAL) {
    F.cnv_V_D2ND(dc_av, C.RefVelocity);
  }
  
  if ( C.isHeatProblem() ) {
    F.loadSphScalar3DAvr(this, fp, "AvrTemperature", G_size, guide, dc_at, step, time, C.Unit.File);
    if (C.Unit.File == DIMENSIONAL) time /= C.Tscale;
    if ( (step != SklGetAverageBaseStep()) || ((REAL_TYPE)time != SklGetAverageBaseTime()) ) {
      Hostonly_ printf     ("\n\tTime stamp is different between files\n");
      Hostonly_ fprintf(fp, "\n\tTime stamp is different between files\n");
      Exit(0);
    }
    
    // convert from dimensional to non-dimensional, iff file is dimensional
    if (C.Unit.File == DIMENSIONAL) {
      F.cnv_T_D2ND(dc_at, C.BaseTemp, C.DiffTemp, C.Unit.Temp);
    }
  }
}


/**
 @fn float SklSolverCBC::min_distance(float* cut, FILE* fp)
 @brief 距離の最小値を求め，閾値以上にする
 @retval 最小距離
 @param cut カット情報の配列
 @param fp file pointer
 @note cutlibのインスタンスを参照, windows does not support std::min() method?
 */
float SklSolverCBC::min_distance(float* cut, FILE* fp)
{
  float min_g = 1.0e6, min_l, c, eps;
  unsigned mm, g, i;
  mm = (size[2]+2*guide)*(size[1]+2*guide)*(size[0]+2*guide)*6;
  g = 0;
  eps = 1.0e-1;
  
  min_g = 1.0;
  for (i=0; i<mm; i++) {
    c = cut[i]; 
    if( min_g > c ) min_g = c;
    if ( c < eps ) {
      cut[i] = eps;
      g++;
    }
  }
  if ( g>0 ) {
    Hostonly_ fprintf(fp, "\t%d %s modified to %5.3e (non-dimnensional distance)\n", g, (g>1)?"cuts were":"cut was", eps);
    Hostonly_ printf("\t%d %s modified to %5.3e (non-dimnensional distance)\n", g, (g>1)?"cuts were":"cut was", eps);
  }
  return min_g;
}


/**
 @fn void SklSolverCBC::prepOutput(void)
 @brief ファイル出力の準備
 */
void SklSolverCBC::prepOutput (void)
{
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
  
  // for instantaneous fields
	if ( C.KindOfSolver != SOLID_CONDUCTION ) {
    
    // Divergence
    if ( C.FIO.Div_Debug == ON ) {
      if ( SklCfgCheckOutFile("Divergence") ) {
        if( !(m_outDiv = InitFile("Divergence", size, org, pit, dc_ws)) ) {
          Hostonly_ stamped_printf("\tInit failed.\n");
          Exit(0);
        }
      }
      else {
        Hostonly_ stamped_printf("\tMissing 'Divergence' in OutFile\n");
        Exit(0);
      }
    }
    
    // Pressure
		if ( SklCfgCheckOutFile("Pressure") ) {
			if( !(m_outPrs = InitFile("Pressure", size, org, pit, dc_ws)) ) {
				Hostonly_ stamped_printf("\tInit failed.\n");
				Exit(0);
			}
		}
		else {
			Hostonly_ stamped_printf("\tMissing 'Pressure' in OutFile\n");
			Exit(0);
		}
    
    // Velocity
		if ( SklCfgCheckOutFile("Velocity") ) {
			if( !(m_outUVW = InitFile("Velocity", size, org, pit, dc_wvex)) ) {
				Hostonly_ stamped_printf("\tInit failed.\n");
				Exit(0);
			}
		}
		else {
			Hostonly_ stamped_printf("\tMissing 'Velocity' in OutFile\n");
			Exit(0);
		}
		
		// Total pressure
		if (C.Mode.TP == ON) {
			if ( SklCfgCheckOutFile("TotalPressure") ) {
				if( !(m_outTP = InitFile("TotalPressure", size, org, pit, dc_ws)) ) {
					Hostonly_ stamped_printf("\tInit failed.\n");
					Exit(0);
				}
			}
			else {
				Hostonly_ stamped_printf("\tMissing 'TotalPressure' in OutFile\n");
				Exit(0);
			}
		}
		
		// Vorticity
		if (C.Mode.VRT == ON) {
			if ( SklCfgCheckOutFile("Vorticity") ) {
				if( !(m_outVrt = InitFile("Vorticity", size, org, pit, dc_wvex)) ) {
					Hostonly_ stamped_printf("\tInit failed.\n");
					Exit(0);
				}
			}
			else {
				Hostonly_ stamped_printf("\tMissing 'Vorticity' in OutFile\n");
				Exit(0);
			}
		}
    
    // 2nd Invariant of Velocity Gradient Tensor
		if (C.Mode.I2VGT == ON) {
			if ( SklCfgCheckOutFile("2ndInvrntVGT") ) {
				if( !(m_outI2VGT = InitFile("2ndInvrntVGT", size, org, pit, dc_ws)) ) {
					Hostonly_ stamped_printf("\tInit failed.\n");
					Exit(0);
				}
			}
			else {
				Hostonly_ stamped_printf("\tMissing '2ndInvrntVGT' in OutFile\n");
				Exit(0);
			}
		}
    
    // Helicity
		if (C.Mode.Helicity == ON) {
			if ( SklCfgCheckOutFile("Helicity") ) {
				if( !(m_outHlcty = InitFile("Helicity", size, org, pit, dc_ws)) ) {
					Hostonly_ stamped_printf("\tInit failed.\n");
					Exit(0);
				}
			}
			else {
				Hostonly_ stamped_printf("\tMissing 'Helicity' in OutFile\n");
				Exit(0);
			}
		}
    
    // Interface function
		if ( C.BasicEqs == INCMP_2PHASE ) {
			if ( SklCfgCheckOutFile("VOF") ) {
				if( !(m_outVrt = InitFile("VOF", size, org, pit, dc_vof)) ) {
					Hostonly_ stamped_printf("\tInit failed.\n");
					Exit(0);
				}
			}
			else {
				Hostonly_ stamped_printf("\tMissing 'VOF' in OutFile\n");
				Exit(0);
			}
		}
	}
  
  if( C.isHeatProblem() ){
    if ( SklCfgCheckOutFile("Temperature") ) {
      if( !(m_outTmp = InitFile("Temperature", size, org, pit, dc_ws)) ) {
        Hostonly_ stamped_printf("\tInit failed.\n");
        Exit(0);
      }
    }
    else {
      Hostonly_ stamped_printf("\tMissing 'Temperature' in OutFile\n");
      Exit(0);
    }
  }
  
  // for averaged fields
  if ( C.Mode.Average == OFF) return;
  
	if ( C.KindOfSolver != SOLID_CONDUCTION ) {
		if ( SklCfgCheckOutFile("AvrPressure") ) {
			if( !(m_outAvrPrs = InitFile("AvrPressure", size, org, pit, dc_ws)) ) {
				Hostonly_ stamped_printf("\tInit failed.\n");
				Exit(0);
			}
		}
		else {
			Hostonly_ stamped_printf("\tMissing 'AvrPressure' in OutFile\n");
			Exit(0);
		}
    
		if ( SklCfgCheckOutFile("AvrVelocity") ) {
			if( !(m_outAvrUVW = InitFile("AvrVelocity", size, org, pit, dc_wvex)) ) {
				Hostonly_ stamped_printf("\tInit failed.\n");
				Exit(0);
			}
		}
		else {
			Hostonly_ stamped_printf("\tMissing 'AvrVelocity' in OutFile\n");
			Exit(0);
		}
	}
  
  if( C.isHeatProblem() ){
    if ( SklCfgCheckOutFile("AvrTemperature") ) {
      if( !(m_outAvrTmp = InitFile("AvrTemperature", size, org, pit, dc_ws)) ) {
        Hostonly_ stamped_printf("\tInit failed.\n");
        Exit(0);
      }
    }
    else {
      Hostonly_ stamped_printf("\tMissing 'AvrTemperature' in OutFile\n");
      Exit(0);
    }
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
  
#ifdef DEBUG
  // Materialの基本リストを表示
  Hostonly_ M->printBaseMaterialList(mp, fp);
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
  
  //debug; printf("%d : [%d %d %d %d %d %d] [%d %d %d]\n", pn.ID, pn.nID[0],pn.nID[1],pn.nID[2],pn.nID[3],pn.nID[4],pn.nID[5],pn.st_idx[0], pn.st_idx[1],pn.st_idx[2]);
}



/**
 @fn void SklSolverCBC::setComponentVF(float* cvf)
 @brief HEX,FANコンポーネントなどの体積率とbboxなどをセット
 @param cvf 体積率
 */
void SklSolverCBC::setComponentVF(float* cvf)
{
  int subsampling = 20; // 体積率のサブサンプリングの基数
  int f_st[3], f_ed[3];
  double flop;
  
  CompoFraction CF(size, guide, C.dx, C.org, subsampling);
  CF.setParallelInfo(pn);
  
  for (int n=1; n<=C.NoBC; n++) {
    
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
        break;
    }
    
    // 回転角度の計算
    CF.get_angle(); 
    
    // bboxと投影面積の計算
    cmp[n].area = CF.get_BboxArea();
    
    // インデクスの計算
    CF.bbox_index(f_st, f_ed);
    
    // インデクスの登録
    cmp[n].setBbox_st(f_st);
    cmp[n].setBbox_ed(f_ed);
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
#ifdef DEBUG
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
}

#ifndef BINARY_VOXEL
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
  SklParaComponent* para_cmp = SklGetParaComponent();
  SklParaManager* para_mng = para_cmp->GetParaManager();
  
  unsigned poly_gc[3];
  float    poly_org[3], poly_dx[3];
  poly_gc[0]  = poly_gc[1] = poly_gc[2] = guide;
  poly_dx[0]  = (float)C.dx[0]*C.RefLength;
  poly_dx[1]  = (float)C.dx[1]*C.RefLength;
  poly_dx[2]  = (float)C.dx[2]*C.RefLength;
  poly_org[0] = (float)C.org[0]*C.RefLength;
  poly_org[1] = (float)C.org[1]*C.RefLength;
  poly_org[2] = (float)C.org[2]*C.RefLength;
  
  //stamped_printf("poly : org(%e %e %e) dimensional\n", poly_org[0], poly_org[1], poly_org[2]);
  
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
  
  // 階層情報表示
  PL->show_group_hierarchy();
  PL->show_group_hierarchy(fp);
  
  Hostonly_ {
    printf("\n");
    fprintf(fp, "\n");
  }
  
  // 使用メモリ量　基本クラスのみ
  unsigned long poly_mem, G_poly_mem;
  G_poly_mem = poly_mem = (unsigned long)PL->used_memory_size();
  m_prep += poly_mem;
  m_total+= poly_mem;
  
  if( para_cmp->IsParallel() ) {
    para_cmp->Allreduce(&poly_mem, &G_poly_mem, 1, SKL_ARRAY_DTYPE_ULONG, SKL_SUM, pn.procGrp);
  }
  
  Hostonly_  {
    FBUtility::displayMemory("Polygon", G_poly_mem, poly_mem, fp, mp);
  }
  
  /* Triangle display
   PolylibNS::Vec3f m_min, m_max, t1(poly_org), t2(poly_dx), t3;
   t3.assign((float)size[0]*t2.t[0], (float)size[1]*t2.t[1], (float)size[2]*t2.t[2]);
   m_min = t1 - t2;      // 1層外側まで
   m_max = t1 + t3 + t2; // 
   printf("min : %f %f %f\n", m_min.t[0], m_min.t[1], m_min.t[2]);
   printf("max : %f %f %f\n", m_max.t[0], m_max.t[1], m_max.t[2]);
   vector<Triangle*>* trias = PL->search_polygons("root", m_min, m_max, false); // false; ポリゴンが一部でもかかる場合
   
   PolylibNS::Vec3f *p, nrl, *n;
   int c=0;
   vector<Triangle*>::iterator it;
   for (it = trias->begin(); it != trias->end(); it++) {
   p = (*it)->get_vertex();
   n = &(*it)->get_normal();
   printf("%d : p0=(%6.3e %6.3e %6.3e)  p1=(%6.3e %6.3e %6.3e) p2=(%6.3e %6.3e %6.3e) n=(%6.3e %6.3e %6.3e)\n", c++, 
   p[0].t[0], p[0].t[1], p[0].t[2],
   p[1].t[0], p[1].t[1], p[1].t[2],
   p[2].t[0], p[2].t[1], p[2].t[2],
   (*n).t[0], (*n).t[1], (*n).t[2]);
   }
   
   delete trias; */ //後始末
  
  
  /* Polylib: STLデータ書き出しテスト
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
   */
  
  // Cutlib
  Hostonly_ {
    fprintf(fp,"\n---------------------------------------------------------------------------\n\n");
    fprintf(fp,"\t>> Cut Info\n\n");
    fprintf(mp,"\n---------------------------------------------------------------------------\n\n");
    fprintf(mp,"\t>> Cut Info\n\n");
  }
  
  size_t n_cell[3];
  size_t ista[3];
  size_t nlen[3];
  size_t size_n_cell;
  
  for ( unsigned i=0; i<3; i++) {
    n_cell[i] = size[i] + 2*guide; // 分割数+ガイドセル
    poly_org[i] -= poly_dx[i]*guide;   // ガイドセル分だけシフト
  }
  size_n_cell = n_cell[0] * n_cell[1] * n_cell[2];
  
  ista[0] = ista[1] = ista[2] = guide;
  nlen[0] = size[0];
  nlen[1] = size[1];
  nlen[2] = size[2];
  
  //stamped_printf("cut : org(%e %e %e) dimensional\n", poly_org[0], poly_org[1], poly_org[2]);
  //stamped_printf("cut : dx(%e %e %e) dimensional\n", poly_dx[0], poly_dx[1], poly_dx[2]);
  
  // Cutlibの配列は各方向(引数)のサイズ
  TIMING_start(tm_init_alloc);
  cutPos = new CutPos32Array(n_cell); // 6*(float)
  cutBid = new CutBid8Array(n_cell);  // 2*(int32_t)
  TIMING_stop(tm_init_alloc);
  
  TIMING_start(tm_cutinfo);
  CutInfoCell(poly_org, poly_dx, PL, cutPos, cutBid); // ガイドセルを含む全領域で交点を計算
  //CutInfoCell(ista, nlen, poly_org, poly_dx, PL, cutPos, cutBid);
  TIMING_stop(tm_cutinfo);
  
  // 使用メモリ量　
  unsigned long cut_mem, G_cut_mem;
  G_cut_mem = cut_mem = (unsigned long)size_n_cell *(6+2)*4;
  m_prep += cut_mem;
  m_total+= cut_mem;
  
  if( para_cmp->IsParallel() ) {
    para_cmp->Allreduce(&cut_mem, &G_cut_mem, 1, SKL_ARRAY_DTYPE_ULONG, SKL_SUM, pn.procGrp);
  }
  
  Hostonly_  {
    FBUtility::displayMemory("Cut", G_cut_mem, cut_mem, fp, mp);
  }
  
  // テスト
  int m, z=0;
  float pos6[6], d_min=1.0;
  BidType bid6[6];
  
  for (int k=ista[2]; k<=ista[2]+nlen[2]-1; k++) {
    for (int j=ista[1]; j<=ista[1]+nlen[1]-1; j++) {
      for (int i=ista[0]; i<=ista[0]+nlen[0]-1; i++) {
        
        cutPos->getPos(i+guide-1,j+guide-1,k+guide-1, pos6);
        cutBid->getBid(i+guide-1,j+guide-1,k+guide-1, bid6);
        
        if ( (pos6[0]+pos6[1]+pos6[2]+pos6[3]+pos6[4]+pos6[5]) < 6.0 ) { // 6方向のうちいずれかにカットがある
          //if ( (bid6[0]+bid6[1]+bid6[2]+bid6[3]+bid6[4]+bid6[5]) >0 ) { // これも有効
          for (int n=0; n<6; n++) {
            if (pos6[n] > 0.0) d_min = min(d_min, pos6[n]);
          }
          z++;
          //printf("%3d %3d %3d : %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %d %d %d %d %d %d\n", i,j,k,
          //     pos6[0], pos6[1], pos6[2], pos6[3], pos6[4], pos6[5],
          //     bid6[0], bid6[1], bid6[2], bid6[3], bid6[4], bid6[5]);
        }
      }
    }
  }
  Hostonly_ printf("\n\tMinimum dist. = %5.3e  : # of cut = %d : %f [percent]\n", d_min, z, (float)z/(float)size_n_cell*100.0);
  
  // カットの最小値を閾値以上にする
  cut = (float*)cutPos->getDataPointer();
  min_distance(cut, fp);
  
}

/**
 @fn void SklSolverCBC::setup_CutInfo4IP(unsigned long& m_prep, unsigned long& m_total, FILE* fp)
 @brief Polylibを準備し，ポリゴンをロードする
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
  
  size_t n_cell[3];
  size_t ista[3];
  size_t nlen[3];
  size_t size_n_cell;
  float  cut_org[3], cut_dx[3];
  
  for ( unsigned i=0; i<3; i++) {
    n_cell[i] = size[i] + 2*guide; // 分割数+ガイドセル
  }
  size_n_cell = n_cell[0] * n_cell[1] * n_cell[2];
  
  ista[0] = ista[1] = ista[2] = guide;
  nlen[0] = size[0];
  nlen[1] = size[1];
  nlen[2] = size[2];
  
  cut_dx[0]  = (float)C.dx[0]*C.RefLength;
  cut_dx[1]  = (float)C.dx[1]*C.RefLength;
  cut_dx[2]  = (float)C.dx[2]*C.RefLength;
  cut_org[0] = (float)( C.org[0] - C.dx[0]*guide )*C.RefLength;
  cut_org[1] = (float)( C.org[1] - C.dx[1]*guide )*C.RefLength;
  cut_org[2] = (float)( C.org[2] - C.dx[2]*guide )*C.RefLength;
  
  // Cutlibの配列は各方向(引数)のサイズ
  TIMING_start(tm_init_alloc);
  cutPos = new CutPos32Array(n_cell); // 6*(float)
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
  cut = (float*)cutPos->getDataPointer();
  
  for (size_t i=0; i<size_n_cell*6; i++) {
    cut[i] = 1.0f;
  }
  
}
#endif

/**
 @fn void SklSolverCBC::VoxEncode(VoxInfo* Vinfo, ParseMat* M, int* mid, float* vf, CutPos32Array* cutPos)
 @brief BCIndexにビット情報をエンコードする
 @param Vinfo
 @param M
 @param mid Voxel IDの配列
 @param vf コンポーネントの体積率 
 @param cutPos カット情報コンテナ
 @note
 - bcdに対して，共通の処理を行い，それをbcp, bcv, bchにコピーし，その後個別処理を行う．
 */
void SklSolverCBC::VoxEncode(VoxInfo* Vinfo, ParseMat* M, int* mid, float* vf, CutPos32Array* cutPos)
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
  C.NoWallSurface = Vinfo->setBCIndexP(bcd, bcp, mid, &BC, cutPos, C.isCDS()); 
  //Vinfo->chkBCIndexP(bcd, bcp, "BCindexP.txt");
  
  // BCIndexV に速度計算のビット情報をエンコードする -----
  Vinfo->setBCIndexV(bcv, mid, &BC, bcp, C.isCDS());
  // debug; Vinfo->chkBCIndexV(bcv, "BCindexV.txt");
  
  // BCIndexT に温度計算のビット情報をエンコードする -----
  if ( C.isHeatProblem() ) {
    Vinfo->setBCIndexH(bcd, bh1, bh2, mid, &BC, C.KindOfSolver);
    // debug; Vinfo->chkBCIndexH(bcv, "BCindexH.txt");
  }
  
  // 内部領域のFluid, Solidのセル数を数える C.Wcell(Local), G_Wcell(global)
  Vinfo->countCellState(C.Wcell, G_Wcell, bcd, SOLID);
  Vinfo->countCellState(C.Fcell, G_Fcell, bcd, FLUID);
  
  // set local active cell ratio
  C.Eff_Cell_Ratio = (REAL_TYPE)C.Acell / C.getCellSize(size);
  
  // コンポーネントのインデクスの再構築
  Vinfo->resizeCompoBV(bcd, bcv, bh1, bh2, C.KindOfSolver, C.isHeatProblem(), GC_bv);
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
  else { // 組み込み例題の場合
    
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
    
    // 各例題のパラメータ設定
    Ex->setDomain(&C, m_sz, m_org, m_wth, m_pch);
  }
  
  //stamped_printf("Model read : org(%e %e %e)\n", m_org[0], m_org[1], m_org[2]);
  
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
  if( para_cmp->IsParallel() ){
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
  
  // 並列処理時のランク情報と各ノードのスタートインデクスを計算
  set_Parallel_Info();
  
  // ノードローカルな値の設定　（無次元値）
  for(int i=0; i<3; i++) {
    m_sz[i] = sz_para[i];
    m_wth[i] = (REAL_TYPE)m_sz[i]*m_pch[i];
    m_org[i] += m_pch[i]*(REAL_TYPE)(pn.st_idx[i]-1);
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
  unsigned cell_id[NOFACE], count=0;
  for (int i=0; i<NOFACE; i++) cell_id[i] = 0;
  
  count = B->count_Outer_Cell_ID(cell_id);
  
  Hostonly_ {
    fprintf(fp, "\tCell IDs on Guide cell region\n");
    for ( int i=0; i<count; i++) {
      fprintf(fp, "\t\tID[%d] = %d\n", i+1, cell_id[i]);
    }
    fprintf(mp, "\tCell IDs on Guide cell region\n");
    for ( int i=0; i<count; i++) {
      fprintf(mp, "\t\tID[%d] = %d\n", i+1, cell_id[i]);
    }
  }
  
  // midにロードされたIDをスキャンし，IDの個数を返し，作業用のcolorList配列にIDを保持，midに含まれるIDの数をチェック
  unsigned sc=0;
  if ( (sc=Vinfo->scanCell(mid, count, cell_id, C.Hide.Change_ID)) > C.NoID ) {
    Hostonly_ stamped_printf("A number of IDs included in voxel model(%d) is grater than one described in 'Model_Setting'(%d)\n", 
                             sc, C.NoID);
    Exit(0);
  }
}

