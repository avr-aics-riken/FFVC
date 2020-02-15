//##################################################################################
//
// Copyright (c) 2016-2019 Research Institute for Information Technology, Kyushu university
// All rights researved.
//
//##################################################################################

/**
 * @file   Cloud.C
 * @brief  Cloud class
 * @author riit
 */

#include "Cloud.h"

//#############################################################################
// @brief ランタイム
bool Cloud::tracking(const unsigned step, const double time)
{
  int mflag = 0;
  
  nOutPart = 0;  // 領域外へ出た粒子数
  nPasPart = 0;  // 壁を通過した粒子数
  
  // 指定時刻が過ぎ、処理ステップになった場合
  if ( Interval.isStarted(step, time) )
  {
    TIMING_start("PT_Tracking");
    for(auto itr = chunkList.begin(); itr != chunkList.end(); ++itr)
    {
      // 保持している粒子を積分し、マイグレーション準備
      mflag += (*itr)->updatePosition(tr, scheme, EmitLife, dt, Rmap, nOutPart, nPasPart);
      
      // 放出タイミングで、登録した開始点から粒子を追加（マイグレーション先は除く）
      if ( Interval.isTriggered(step, time) )
      {
        (*itr)->addParticleFromOrigin();
      }
    }
    
    // 粒子の放出回数
    if ( Interval.isTriggered(step, time) ) nEmission++;
    
    // 全プロセスで共有  最大値をとる => どれか１プロセスでもマイグレーションが発生
    int tmp = mflag;
    if ( MPI_SUCCESS != MPI_Allreduce(&tmp,
                                      &mflag,
                                      1,
                                      MPI_INT,
                                      MPI_MAX,
                                      MPI_COMM_WORLD) ) return false;
    TIMING_stop("PT_Tracking");
    
    
    TIMING_start("PT_Migration");
    // マイグレーションが生じる場合のみ
    if ( mflag > 0 )
    {
      // 送信方向毎に全チャンクの送信粒子数の和をとる
      int* pInfo = PC.pInfo_ptr();
      memset(pInfo, 0, sizeof(int)*NDST*2);
      
      for(auto itr = chunkList.begin(); itr != chunkList.end(); ++itr)
      {
        int* p = (*itr)->pSend_ptr();
        
        for (int i=0; i<NDST; i++)
        {
          pInfo[2*i] += p[i];
        }
      }
      
      // 送信要素の最大値を保存  行き先毎には異なるが、最大数で考えておく
      unsigned max_particle = 0;
      
      for (int i=0; i<NDST; i++) {
        unsigned p = (unsigned)pInfo[2*i];
        if (max_particle < p) {
          max_particle = p;
        }
      }
      
      // バッファ長の更新が必要な場合、バッファ計算用の粒子数はBUF_UNIT単位で確保
      if (buf_max_particle < max_particle)
      {
        buf_max_particle = (max_particle / BUF_UNIT + 1) * BUF_UNIT;
        buf_updated = true;
      }
      
      // バッファ長を同期
      if ( !PC.getMax(buf_max_particle) ) return false;
      
      //Hostonly_ printf("buf_max_p= %ld\n", buf_max_particle);
      
      
      // マイグレーション処理
      if ( !migration() ) return false;
      
    } // mflag
    TIMING_stop("PT_Migration");
    
    
    // 出力
    if ( (step/OutInterval)*OutInterval == step )
    {
      // ログ出力
      TIMING_start("PT_Statistics");
      nParticle = getNparticle();
      
      if ( !PC.Statistics(nCommParticle,
                          nParticle,
                          gParticle,
                          nOutPart,
                          nPasPart) ) return false;
      
      sleepParticle +=  nOutPart + nPasPart;
      
      Hostonly_
      {
        logging(step);
      }
      TIMING_stop("PT_Statistics");
      
      
      
      // ファイル出力
      TIMING_start("PT_fileIO");
      if ( out_format == 0 )
      {
        if ( !write_ascii(step, time) ) return false;
      }
      else
      {
        if ( !write_binary(step, time) ) return false;
      }
      TIMING_stop("PT_fileIO");
    }
    
  } // Interval.isStarted
  
  return true;
}



//#############################################################################
// @brief 初期設定
bool Cloud::initCloud(FILE* fp)
{
  tr = new Tracking(size,
                    guide,
                    origin,
                    region,
                    pitch,
                    vSource,
                    bcd,
                    bid,
                    myRank);
  
  if ( !setPTinfo() ) return false;
  
  if ( !determineUniqueID() ) return false;
  
  displayParam(stdout);
  displayParam(fp);
  
  // リスタート時
  if ( restartFlag == 1 )
  {
    // 粒子をもつランク番号情報
    if ( !readRestartRank() ) return false;
    
    if ( restartForm == 0 ) // Ascii
    {
      if ( !readRestartParticleAscii() ) return  false;
    }
    else // Binary
    {
      if ( !readRestartParticleBinary() ) return  false;
    }
  }
  
  // ログ出力
  Hostonly_
  {
    if ( !(fpl=fopen("pt_log.txt", "w")) )
    {
      stamped_printf("\tSorry, can't open 'pt_log.txt' file. Write failed.\n");
      return false;
    }
  }
  
  // 通信クラスの設定
  Rmap = PC.setPtComm(G_division, myRank, numProc);
  
  
  return true;
}




//#############################################################################
// @brief マイグレーション処理
/* @note 方針
 - Chunk::updatePosition()で、マイグレーション候補にマーク、送信先毎の個数をカウント
 - 送信バッファは周囲のランク26方向に対してvectorで管理、毎回利用直前にクリアして使う
 - 粒子と情報をバッファに詰め込み、リストから削除、送信
 - 受信し、アンパック、情報を修正
 */
bool Cloud::migration()
{
  // マイグレーション用に確保したバッファ領域の準備
  if (buf_updated) {
    PC.resizeSendBuffer(buf_max_particle);
    buf_updated = false;
  }
  else
  {
    // バッファサイズの更新がない場合は再初期化のみ
    PC.initSendBuffer();
  }
  
  mark();
  
  TIMING_start("PT_ParticlePacking");
  // 各チャンクのマイグレーション候補を見つけ、行き先毎にパック
  int c[NDST]={0};
  for(auto itr = chunkList.begin(); itr != chunkList.end(); ++itr)
  {
    (*itr)->packParticle(PC.ps_ptr(),
                         PC.bs_ptr(),
                         buf_max_particle,
                         c);
  }
  TIMING_stop("PT_ParticlePacking");
  mark();
  
  // 周囲のランクと通信し、経路と送受信データ数を確定
  TIMING_start("PT_EstablishCommPath");
  if ( !PC.establishCommPath() ) return false;
  TIMING_stop("PT_EstablishCommPath");
  mark();
  
  // データ本体の送受信
  TIMING_start("PT_CommParticle");
  if ( !PC.commParticle(buf_max_particle) ) return false;
  TIMING_stop("PT_CommParticle");
  mark();
  
  // アンパック
  TIMING_start("PT_ParticleUnpacking");
  unpackParticle();
  TIMING_stop("PT_ParticleUnpacking");
  mark();
  
  return true;
}


//#############################################################################
// @brief 受信粒子のアンパック
void Cloud::unpackParticle()
{
  int* pInfo = PC.pInfo_ptr();
  REAL_TYPE* pr_buf = PC.pr_ptr();
  int* br_buf = PC.br_ptr();
  
  // 周囲26方向のバッファの区切り
  size_t bsz6 = buf_max_particle * 6; // 6要素
  size_t bsz2 = buf_max_particle * 2; // 2要素
  
  for (int i=0; i<NDST; i++)
  {
    const int cnt = pInfo[2*i+1] ;  // 受信粒子数
    
    int b, lid;
    Vec3r pos;
    Vec3r v;
    particle p;
    
    if (cnt > 0 && i != 13 )
    {
      for (int j=0; j<cnt; j++)
      {
        pos.x = pr_buf[bsz6*i + 6*j+0];
        pos.y = pr_buf[bsz6*i + 6*j+1];
        pos.z = pr_buf[bsz6*i + 6*j+2];
        v.x   = pr_buf[bsz6*i + 6*j+3];
        v.y   = pr_buf[bsz6*i + 6*j+4];
        v.z   = pr_buf[bsz6*i + 6*j+5];
        b     = br_buf[bsz2*i + 2*j+0];
        lid   = br_buf[bsz2*i + 2*j+1];
        
        p.pos = pos;
        p.bf  = Chunk::removeMigrate(b);
        p.vel = v;
        p.foo = lid;
        
        addParticle2ChunkList(p);
      }
    } // cnt
    
  } // NDST-loop
  
}



//#############################################################################
// @brief ログ出力
void Cloud::logging(const unsigned step)
{
  // マイグレーションで移動した粒子数と全粒子数
  fprintf(fpl, "step %ld migrated %d total %ld sleep %ld OutRegion %d PassWall %d\n",
          step, nCommParticle, gParticle, sleepParticle, nOutPart, nPasPart);
  
  // 各ランク毎の粒子数
  unsigned* p = PC.nPart_ptr();
  for (int i=0; i<numProc; i++) {
    fprintf(fpl, "%ld ", p[i]);
  }
  fprintf(fpl,"\n");
  
  // Number of total particles
  unsigned tp = nEmitParticle * (nEmission+1);
  if (gParticle != tp) fprintf(fpl, "Number of particles must be %ld %ld\n", tp, nEmission);
  
  fflush(fpl);
}


//#############################################################################
// @brief 開始点のユニークIDを割り振る
bool Cloud::determineUniqueID()
{
  unsigned* uid = new unsigned[numProc];
  memset(uid, 0, sizeof(unsigned)*numProc);
  unsigned sum=0;
  
  if (numProc > 1)
  {
    unsigned sbuf = nParticle;
    if ( MPI_SUCCESS != MPI_Allgather(&sbuf,
                                      1,
                                      MPI_UNSIGNED,
                                      uid,
                                      1,
                                      MPI_UNSIGNED,
                                      MPI_COMM_WORLD) ) return false;
    
    for (int i=0; i<numProc; i++) sum += uid[i];
  }
  else
  {
    sum = nParticle;
  }
  gParticle = sum;
  
  Hostonly_ printf("\tGlobal # of particles  = %ld\n", gParticle);
  
  // 各ランクのユニークIDの先頭番号
  unsigned* acc = new unsigned[numProc];
  
  acc[0] = 0;
  for (int i=1; i<numProc; i++)
  {
    acc[i] = acc[i-1] + uid[i-1];
  }
  
  
  // ユニークIDの振り直し
  unsigned nc = 0;
  for(auto itr = chunkList.begin(); itr != chunkList.end(); ++itr)
  {
    int id = acc[myRank]+nc;
    
    (*itr)->setUid(id);
    (*itr)->setPchunkUid(id);
    nc++;
    
  }
  
  if (nc != chunkList.size()) {
    printf("chunkList size is not proper\n");
    return false;
  }
  
#ifdef PT_DEBUG
  printf("*\trank=%d : nParticle=%d\n", myRank, nParticle);
#endif
  
  
  delete [] uid;
  delete [] acc;
  
  return true;
}


//#############################################################################
// @brief 粒子追跡情報を取得し，chunkに保持する
bool Cloud::setPTinfo()
{
  Monitor_Type mon_type;
  string str, label;
  string label_base, label_leaf;
  string name;
  double f_val=0.0;
  
  /* 粒子出力  >>  FFVの制御部で記述
   label = "/ParticleTracking/Tracking";
   
   if ( !(tpCntl->getInspectedValue(label, str )) )
   {
   Hostonly_ stamped_printf("\tError : '%s'\n", label.c_str());
   return false;
   }
   else
   {
   if ( !strcasecmp(str.c_str(), "on") ) tracking = ON;
   if ( !strcasecmp(str.c_str(), "off")) tracking = OFF;
   }
   */
  
  
  // 入力パラメータの次元
  label = "/unit/UnitOfInputParameter";
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
    Hostonly_ stamped_printf("\tParsing warning : No Label in '%s'\n", label.c_str());
    return false;
  }
  if ( !strcasecmp(str.c_str(), "NONDIMENSIONAL") ) {
    unit = NONDIMENSIONAL;
  }
  else if ( !strcasecmp(str.c_str(), "DIMENSIONAL") ) {
    unit = DIMENSIONAL;
  }
  else {
    Hostonly_ printf("Invalid unit\n");
    return false;
  }
  
  
  // 代表長
  label = "/reference/length";
  if ( !(tpCntl->getInspectedValue(label, f_val )) )
  {
    Hostonly_ stamped_printf("\tParsing warning : No Label in '%s'\n", label.c_str());
    return false;
  }
  refLen = (REAL_TYPE)f_val;
  
  
  // 代表速度
  label = "/reference/velocity";
  if ( !(tpCntl->getInspectedValue(label, f_val )) )
  {
    Hostonly_ stamped_printf("\tParsing warning : No Label in '%s'\n", label.c_str());
    return false;
  }
  refVel = (REAL_TYPE)f_val;
  
  
  
  label_base = "/ParticleTracking";
  label = label_base + "/method";
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
    Hostonly_ stamped_printf("\tParsing warning : No Label in '%s'\n", label.c_str());
    return false;
  }
  if ( !strcasecmp(str.c_str(), "euler") ) {
    scheme = pt_euler;
  }
  else if ( !strcasecmp(str.c_str(), "rk2") ) {
    scheme = pt_rk2;
  }
  else if ( !strcasecmp(str.c_str(), "rk4") ) {
    scheme = pt_rk4;
  }
  else {
    Hostonly_ printf("Invalid particle tracking method\n");
    return false;
  }
  
  
  // ファイル出力フォーマット
  label = label_base + "/file_format";
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
    Hostonly_ stamped_printf("\tParsing warning : No Label in '%s'\n", label.c_str());
    return false;
  }
  
  if ( !strcasecmp(str.c_str(), "ascii") ) {
    out_format = 0;
  }
  else if ( !strcasecmp(str.c_str(), "binary") ) {
    out_format = 1;
  }
  else {
    Hostonly_ printf("Invalid file format\n");
    return false;
  }
  
  
  // 粒子出力
  label = label_base + "/StartEmit";
  
  if ( !(tpCntl->getInspectedValue(label, f_val )) )
  {
    Hostonly_ stamped_printf("\tParsing error : Invalid string for '%s'\n", label.c_str());
    return false;
  }
  else
  {
    EmitStart = (int)f_val;
    Interval.setStart(EmitStart);
  }
  
  
  // 放出インターバル
  label = label_base + "/EmitInterval";
  
  if ( !(tpCntl->getInspectedValue(label, f_val )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  else
  {
    EmitInterval = (int)f_val;
    Interval.setInterval(EmitInterval);
  }
  
  
  // 出力インターバル
  label = label_base + "/OutInterval";
  
  if ( !(tpCntl->getInspectedValue(label, f_val )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  else
  {
    OutInterval = (int)f_val;
  }
  
  
  
  // 放出カウント
  label =  label_base + "/EmitLife";
  
  if ( !(tpCntl->getInspectedValue(label, f_val )) )
  {
    Hostonly_ stamped_printf("\tParsing error : No lifetime\n");
    return false;
  }
  else
  {
    int tmp = (int)f_val;
    if (tmp > MAX_LIFE) {
      Hostonly_ printf("Exceed MAX_LIFE\n");
      return false;
    }
    EmitLife = tmp;
  }
  
  
  // リスタート
  label = label_base + "/restart/step";
  
  if ( !(tpCntl->getInspectedValue(label, f_val )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  else
  {
    restartStep = (int)f_val;
  }
  
  if ( 0 != restartStep ) restartFlag = 1;
  
  
  // リスタート粒子ファイルのフォーマット
  label = label_base + "/restart/format";
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  else
  {
    if ( !strcasecmp(str.c_str(), "ascii") )
    {
      restartForm = 0;
    }
    else if ( !strcasecmp(str.c_str(), "binary") )
    {
      restartForm = 1;
    }
    else
    {
      Hostonly_ stamped_printf("\tParsing error : invalid keyword '%s'\n", str.c_str());
      return false;
    }
  }
  
  
  
  // 粒子放出開始グループ数のチェック
  int nnode=0;
  int nlist=0;
  
  nnode = tpCntl->countLabels(label_base);
  if ( nnode == 0 )
  {
    stamped_printf("\tcountLabels --- %s\n", label_base.c_str());
    return false;
  }
  
  for (int i=0; i<nnode; i++)
  {
    if ( !(tpCntl->getNodeStr(label_base, i+1, str)) )
    {
      printf("\tParsing error : No List[@]\n");
      return false;
    }
    
    if( strcasecmp(str.substr(0,4).c_str(), "List") ) continue;
    
    nlist++;
  }
  
  if ( nlist==0 )
  {
    Hostonly_ stamped_printf("\tError : No start points. Please confirm 'ParticleTracking' in Input parameter file. \n");
  }
  
  
  // ソルバー計算ステップ情報の取得
  
  // スタート
  label = "/TimeControl/Session/Start";
  if ( !(tpCntl->getInspectedValue(label, f_val )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  unsigned m_start = f_val;
  
  if ( !Interval.initTrigger(m_start,
                             (double)m_start * (double)dt,
                             (double)dt) ) return false;
  
  
  /* 終了
   label = "/TimeControl/Session/End";
   if ( !(tpCntl->getInspectedValue(label, f_val )) )
   {
   Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
   return false;
   }
   unsigned m_end = f_val;
   */
  
  
  // EmitGrp[]の作成
  nGrpEmit = nlist;
  Egrp = new EmitGroup[nGrpEmit];
  
#ifdef PT_DEBUG
  Hostonly_ printf("*\tnGrpEmit=%d\n", nGrpEmit);
#endif
  
  
  
  // リストの読み込み
  label_base = "/ParticleTracking";
  
  int odr =0;
  for (int i=0; i<nnode; i++)
  {
    if ( !(tpCntl->getNodeStr(label_base, i+1, str)) )
    {
      printf("\tParsing error : No List[@]\n");
      return false;
    }
    if( strcasecmp(str.substr(0,4).c_str(), "List") ) continue;
    
    
    label_leaf = label_base + "/" + str;
    
    // point distribution type
    label = label_leaf + "/Type";
    
    if ( !(tpCntl->getInspectedValue(label, str)) )
    {
      stamped_printf("\tParsing error : No entory '%s'\n", label.c_str());
      return false;
    }
    
    if ( !strcasecmp(str.c_str(), "PointSet") )
    {
      mon_type = mon_POINT_SET;
    }
    else if ( !strcasecmp(str.c_str(), "Line") )
    {
      mon_type = mon_LINE;
    }
    /*
     else if ( !strcasecmp(str.c_str(), "Cylinder") )
     {
     mon_type = mon_CYLINDER;
     }
     else if ( !strcasecmp(str.c_str(), "Box") )
     {
     mon_type = mon_BOX;
     }
     */
    else if ( !strcasecmp(str.c_str(), "Disc") )
    {
      mon_type = mon_DISC;
    }
    else
    {
      Hostonly_ stamped_printf("\tParsing error : Invalid keyword in label='%s', str=%s\n", label.c_str(), str.c_str());
      return false;
    }
    
    
    // Labelの取得
    label = label_leaf + "/Label";
    
    if ( !(tpCntl->getInspectedValue(label, str )) )
    {
      Hostonly_ stamped_printf("\tParsing warning : No Label in '%s'\n", label.c_str());
      return false;
    }
    else
    {
      Egrp[odr].setGrpName(str);
    }
    
    
    
    // 順番に粒子IDとチャンクIDを付与
    // 粒子IDは各サブドメイン毎の仮のもの
    switch(mon_type)
    {
      case mon_POINT_SET:
        Egrp[odr].setType(mon_POINT_SET);
        setPointset(label_leaf, odr);
        break;
        
      case mon_LINE:
        Egrp[odr].setType(mon_LINE);
        setLine(label_leaf, odr);
        break;
        
      case mon_DISC:
        Egrp[odr].setType(mon_DISC);
        setDisc(label_leaf, odr);
        break;
        
      default:
        return false;
        break;
    }
    
    odr++;
  } // nnode
  
  
  nParticle = getNparticle();
  
  return true;
}



//#############################################################################
/**
 * @brief 座標情報を取得し、chunkに保持(PointSet)
 * @param [in]  odr        グループ登録インデクス
 * @note データは無次元化して保持
 */
bool Cloud::setPointset(const string label_base, const int odr)
{
  char tmpstr[20];
  string str, label, label_leaf;
  
  // PointSet個数のチェック
  int nnode=0;
  int nlist=0;
  
  nnode = tpCntl->countLabels(label_base);
  
  if ( nnode == 0 )
  {
    stamped_printf("\tcountLabels --- %s\n", label_base.c_str());
    return false;
  }
  
  for (int i=0; i<nnode; i++)
  {
    if ( !(tpCntl->getNodeStr(label_base, i+1, str)) )
    {
      printf("\tParsing error : No node name\n");
      return false;
    }
    
    if ( strcasecmp(str.substr(0,3).c_str(), "set") ) continue;
    nlist++;
  }
  
  // PointSet取得
  int pc=0;
  for (int i=0; i<nnode; i++)
  {
    if ( !(tpCntl->getNodeStr(label_base, i+1, str)) )
    {
      printf("\tParsing error : No Elem name\n");
      return false;
    }
    
    if ( strcasecmp(str.substr(0,3).c_str(), "set") ) continue;
    pc++;
    
    label_leaf = label_base + "/" + str;
    
    // set coordinate
    label = label_leaf + "/Coordinate";
    REAL_TYPE v[3];
    for (int n=0; n<3; n++) v[n]=0.0;
    
    if ( !(tpCntl->getInspectedVector(label, v, 3 )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      return false;
    }
    
    // 入力パラメータの次元が有次元のとき，無次元化する
    if (unit == DIMENSIONAL) normalizeCord(v);
    
    
    // set tagの取得．ラベルなしでもエラーではない
    label = label_leaf + "/tag";
    
    if ( !(tpCntl->getInspectedValue(label, str)) )
    {
      Hostonly_ stamped_printf("\tParsing warning : No tag for '%s'\n", label.c_str());
    }
    if ( !strcasecmp(str.c_str(), "") )
    {
      sprintf(tmpstr, "point_%d", pc);
      str = tmpstr;
    }
    
    // registration
    Vec3r pos(v);
    registChunk(pos);
  }
  
  // 放出点の合計
  nEmitParticle += nnode;
  
  return true;
}


//#############################################################################
/**
 * @brief chunkを登録
 * @param [in]  pos    座標ベクトル
 */
void Cloud::registChunk(Vec3r pos)
{
  // 自領域内であれば、初期開始点として追加
  if ( tr->inOwnRegion(pos) )
  {
    if (restartFlag == 0) // 初期計算時は座標値を登録
    {
      Chunk* m = new Chunk(pos,
                           1,
                           EmitStart,
                           myRank,
                           EmitInterval);
      chunkList.push_back(m);
    }
    else // リスタート時は、座標点を登録しない、あとでリスタートファイルをロードする
    {
      Chunk* m = new Chunk(pos,
                           1,
                           EmitStart,
                           myRank,
                           EmitInterval,
                           false);
      chunkList.push_back(m);
    }
  }
  
}


//#############################################################################
/**
 * @brief 開始座標情報を取得し、chunkに保持(Line)
 * @param [in]  odr        グループ登録インデクス
 * @note データは無次元化して保持
 */
bool Cloud::setLine(const string label_base, const int odr)
{
  string str, label;
  REAL_TYPE v[3];
  int nDivision;
  REAL_TYPE from[3], to[3];
  
  label = label_base + "/Division";
  
  if ( !(tpCntl->getInspectedValue(label, nDivision )) )
  {
    Hostonly_ stamped_printf("\tParsing error : No Division\n");
    return false;
  }
  if ( nDivision == 0 ) return false;
  
  Egrp[odr].nDiv = nDivision;
  
  // 放出点の合計
  nEmitParticle += nDivision + 1;
  
  // load parameter of 'from' and 'to'
  label = label_base + "/From";
  
  for (int n=0; n<3; n++) v[n]=0.0;
  if ( !(tpCntl->getInspectedVector(label, v, 3 )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'From' in 'Line'\n");
    return false;
  }
  Egrp[odr].from[0] = from[0]=v[0];
  Egrp[odr].from[1] = from[1]=v[1];
  Egrp[odr].from[2] = from[2]=v[2];
  
  
  // 入力パラメータの次元が有次元のとき，無次元化する
  if (unit == DIMENSIONAL) normalizeCord(from);
  
  
  label=label_base+"/To";
  
  for (int n=0; n<3; n++) v[n]=0.0;
  if ( !(tpCntl->getInspectedVector(label, v, 3 )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get 'To' in 'Line'\n");
    return false;
  }
  Egrp[odr].to[0] = to[0]=v[0];
  Egrp[odr].to[1] = to[1]=v[1];
  Egrp[odr].to[2] = to[2]=v[2];
  
  // 入力パラメータの次元が有次元のとき，無次元化する
  if (unit == DIMENSIONAL) normalizeCord(to);
  
  
  // 点群の発生と登録
  
  int npnt = nDivision + 1;
  if (npnt < 2) return false;
  
  Vec3r st(from);
  Vec3r ed(to);
  Vec3r dd = ed - st;
  dd /= (REAL_TYPE)npnt - 1.0;
  
  
  for (int m = 0; m < npnt; m++)
  {
    Vec3r pos = st + dd * (REAL_TYPE)m;
    registChunk(pos);
  }
  return true;
}



//#############################################################################
/**
 * @brief 開始座標情報を取得し、chunkに保持(Disc)
 * @param [in]  odr        グループ登録インデクス
 * @note データは無次元化して保持
 */
bool Cloud::setDisc(const string label_base, const int odr)
{
  string str, label;
  int nSample;
  REAL_TYPE cnt[3], nv[3], radius;
  
  label=label_base+"/nSample";
  
  if ( !(tpCntl->getInspectedValue(label, nSample )) )
  {
    Hostonly_ stamped_printf("\tParsing error : No nSample\n");
    return false;
  }
  if ( nSample <= 0 ) return false;
  Egrp[odr].nSample = nSample;
  
  // 放出点の合計
  nEmitParticle += nSample;
  
  // load parameter
  label=label_base+"/center";
  
  if ( !(tpCntl->getInspectedVector(label, cnt, 3 )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  
  // 入力パラメータの次元が有次元のとき，無次元化する
  if (unit == DIMENSIONAL) normalizeCord(cnt);
  
  Egrp[odr].center[0] = cnt[0];
  Egrp[odr].center[1] = cnt[1];
  Egrp[odr].center[2] = cnt[2];
  
  
  label=label_base+"/normal";
  
  if ( !(tpCntl->getInspectedVector(label, nv, 3 )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  
  // 単位ベクトル化
  getUnitVector(nv);
  
  Egrp[odr].normal[0] = nv[0];
  Egrp[odr].normal[1] = nv[1];
  Egrp[odr].normal[2] = nv[2];
  
  
  label=label_base+"/radius";
  
  if ( !(tpCntl->getInspectedValue(label, radius )) )
  {
    Hostonly_ stamped_printf("\tParsing error : No Radius\n");
    return false;
  }
  if ( radius <= 0.0 ) return false;
  
  // 入力パラメータの次元が有次元のとき，無次元化する
  if (unit == DIMENSIONAL) radius /= refLen;
  
  Egrp[odr].radius = radius;
  
  
  // 点群の発生と登録
  samplingInCircle(cnt, nv, radius, nSample);
  
  return true;
}


//#############################################################################
/// @brief 半径r内のサンプリング
/// @param [in]  cnt        中心座標
/// @param [in]  nv         法線方向
/// @param [in]  radius     半径
/// @param [in]  nSample    サンプル数
/// @note 指定半径より少し内側（98%）にする
///       1) マスターが点群を発生
///       2) 全ランクにbroadcast
///       3) 各ランクで自領域のみ登録
/// @url https://teramonagi.hatenablog.com/entry/20141113/1415839321
bool Cloud::samplingInCircle(const REAL_TYPE* cnt,
                             const REAL_TYPE* nv,
                             const REAL_TYPE radius,
                             const int nSample)
{
  // 送受信バッファ
  REAL_TYPE* buf = new REAL_TYPE [3*nSample];
  
  // マスターランクで点群を発生
  Hostonly_ {
    REAL_TYPE theta, r, x, y;
    double pi = 2.0*asin(1.0);
    Vec3r q, t;
    Vec3r center(cnt);
    Vec3r normal(nv);
    vector<Vec3r> sampled;
    
    
    Vec3r angle = getAngle(normal);
    
    for (int i=0; i<nSample; i++)
    {
      theta = mts(0.0, 2.0*pi);
      r = sqrt( 2.0 * mts(0.0, (double)(0.5*radius*radius)) ) * 0.98;
      //r = mts(0.0, (double)radius) * 0.98;
      
      x = r * cos(theta);
      y = r * sin(theta);
      q = rotate_inv(angle, t.assign(x, y, 0.0)) + center;
      
      sampled.push_back(q);
    }
    
    // バッファに詰める
    int m = 0;
    for(auto itr = sampled.begin(); itr != sampled.end(); ++itr)
    {
      buf[3*m+0] = (*itr).x;
      buf[3*m+1] = (*itr).y;
      buf[3*m+2] = (*itr).z;
      m++;
    }
  } // Hostonly
  
  
  
  // 全ランクに送る
  if ( !PC.BcastParticles(3*nSample, buf) ) return false;
  
  
  // バッファから取り出し、登録
  for (int i=0; i<nSample; i++)
  {
    Vec3r p;
    p.x = buf[3*i+0];
    p.y = buf[3*i+1];
    p.z = buf[3*i+2];
    
    registChunk(p);
  }
  
  
  delete [] buf;
  buf = NULL;
  
  return false;
}



//#############################################################################
// 指定法線nvがz軸の方向ベクトルに向かう回転角を計算する
// 回転の符号はz軸に向かう回転が右ねじ方向の場合を正にとる
Vec3r Cloud::getAngle(Vec3r nv)
{
  REAL_TYPE alpha, beta, c, d, c_alp, c_bta;
  REAL_TYPE eps = 1.0e-5, f_yz, f_xz;
  Vec3r p, q, angle;
  Vec3r z(0.0, 0.0, 1.0);
  //Vec3r dir;         ///< 矩形の方向規定の参照ベクトル
  
  // 単位ベクトルnvがz軸の単位ベクトルと作る角度を返す
  // yz面への射影
  p.x = 0.0;
  p.y = nv.y;
  p.z = nv.z;
  c = p.length();
  
  if ( c != 0.0 )
  {
    c_alp = dot(z, p)/c;
    d = acos( c_alp );
    f_yz = c_alp+1.0;
    alpha = (nv.y >= 0.0) ? d : -d; // alpahはx軸周りの回転角
  }
  else
  {
    alpha = 0.0; // yz面への射影ベクトルがゼロの場合には回転しない
  }
  
  
  // 参照ベクトルをalphaだけ回転して評価ベクトルを生成 > xz面への射影
  q.assign(alpha, 0.0, 0.0);
  p = rotate(q, nv);
  c = p.length();
  
  if ( c != 0.0 )
  {
    c_bta = dot(z, p)/c;
    d = acos( c_bta );
    f_xz = c_bta+1.0;
    beta = (nv.x >= 0.0) ? -d : d; // betaはy軸まわりの回転角
  }
  else
  {
    beta = 0.0;
  }
  
  // pがz-方向の場合にだけ，y軸回りのみ回転させてz+方向にする
  if ( (f_yz<eps) && (f_xz<eps) )
  {
    alpha = 0.0;
    beta  = acos(-1.0);
  }
  
  angle.assign(alpha, beta, 0.0);
  
  
  /* 矩形の場合，単位ベクトルdirが回転した後，x軸の単位ベクトルへ回転する角度を計算
   if ( mon_type == mon_BOX )
   {
   REAL_TYPE c_gma, f_xy;
   Vec3r x(1.0, 0.0, 0.0);
   
   q = rotate(angle, dir); // 回転によりxy平面上に射影される > q.z=0.0
   c = q.length();
   
   if ( c != 0.0 )
   {
   c_gma = dot(x, q)/c;
   d = acos( c_gma );
   f_xy = c_gma+1.0;
   
   if ( f_xy<eps )
   {
   angle.z = 2.0*asin(1.0); // 反対方向なのでπ
   }
   else
   {
   angle.z = (q.y >= 0.0) ? -d : d;
   }
   }
   else
   {
   stamped_printf("\tInvalid Parameter of Heat exchanger : lateral vector is zero\n");
   Exit(0);
   }
   }
   */
  
  // ##########
#if 0
  stamped_printf("rotation angle = (%f %f %f)\n", angle.x, angle.y, angle.z);
#endif
  // ##########
  
  return angle;
}


//#############################################################################
/**
 * @brief 回転ベクトルp(alpha, beta, gamma)でベクトルuを回転する
 * @param [in] p 回転角度
 * @param [in] u 方向ベクトル
 * @return 角度
 * @note sin, cos = 5flop, dot=5flop, total 181flop
 */
Vec3r Cloud::rotate(const Vec3r p, const Vec3r u)
{
  Vec3r a, b, c;
  
  // line vector expression
  a.x =  cos(p.y)*cos(p.z);
  a.y =  sin(p.x)*sin(p.y)*cos(p.z) - cos(p.x)*sin(p.z);
  a.z =  cos(p.x)*sin(p.y)*cos(p.z) + sin(p.x)*sin(p.z);
  
  b.x =  cos(p.y)*sin(p.z);
  b.y =  sin(p.x)*sin(p.y)*sin(p.z) + cos(p.x)*cos(p.z);
  b.z =  cos(p.x)*sin(p.y)*sin(p.z) - sin(p.x)*cos(p.z);
  
  c.x = -sin(p.y);
  c.y =  sin(p.x)*cos(p.y);
  c.z =  cos(p.x)*cos(p.y);
  
  return Vec3r( dot(a, u), dot(b, u), dot(c, u) );
}


//#############################################################################
/**
 * @brief 回転ベクトルp(alpha, beta, gamma)に対して，-pでベクトルuを回転する
 * @param [in] p 回転角度
 * @param [in] u 方向ベクトル
 * @return 角度
 * @note sin, cos = 5flop, dot=5flop, total 181flop
 */
Vec3r Cloud::rotate_inv(const Vec3r p, const Vec3r u)
{
  Vec3r a, b, c;
  
  // line vector expression
  a.x =  cos(p.y)*cos(p.z);
  a.y =  cos(p.y)*sin(p.z);
  a.z = -sin(p.y);
  
  b.x =  sin(p.x)*sin(p.y)*cos(p.z) - cos(p.x)*sin(p.z);
  b.y =  sin(p.x)*sin(p.y)*sin(p.z) + cos(p.x)*cos(p.z);
  b.z =  sin(p.x)*cos(p.y);
  
  c.x =  cos(p.x)*sin(p.y)*cos(p.z) + sin(p.x)*sin(p.z);
  c.y =  cos(p.x)*sin(p.y)*sin(p.z) - sin(p.x)*cos(p.z);
  c.z =  cos(p.x)*cos(p.y);
  
  return Vec3r( dot(a, u), dot(b, u), dot(c, u) );
}


//#############################################################################
// @brief パラメータ表示
void Cloud::displayParam(FILE* fp)
{
  Hostonly_ {
    fprintf(fp,"\tDescribed unit      : %s\n", (unit==DIMENSIONAL) ? "DIMENSIONAL" : "NON-DIMENSIONAL");
    fprintf(fp,"\tReference Length    : %e\n", refLen);
    fprintf(fp,"\tReference Velocity  : %e\n", refVel);
    fprintf(fp,"\tIntegration Method  : ");
    switch (scheme) {
      case pt_euler:
        fprintf(fp,"Euler\n");
        break;
      case pt_rk2:
        fprintf(fp,"RK2\n");
        break;
      case pt_rk4:
        fprintf(fp,"RK4\n");
        break;
    }
    
    fprintf(fp,"\tStart step          : %d\n",EmitStart);
    fprintf(fp,"\tEmitInterval        : %d\n",EmitInterval);
    fprintf(fp,"\tLife time           : %d\n",EmitLife);
    fprintf(fp,"\tFile format         : %s\n", (out_format==0)? "ASCII":"BINARY");
    
    for (int i=0; i<nGrpEmit; i++)
    {
      fprintf(fp,"\n\tGroup               : %s\n",Egrp[i].getGrpName().c_str());
      fprintf(fp,"\t\tType                : ");
      switch (Egrp[i].getType()) {
        case mon_POINT_SET:
          fprintf(fp,"\t\tPointset\n");
          break;
          
        case mon_LINE:
          fprintf(fp,"\t\tLine\n");
          fprintf(fp,"\t\t# of division       : %d\n",Egrp[i].nDiv);
          fprintf(fp,"\t\tFrom                : (%14.6e %14.6e %14.6e) / (%14.6e %14.6e %14.6e) [m]\n",
                  Egrp[i].from[0],
                  Egrp[i].from[1],
                  Egrp[i].from[2],
                  Egrp[i].from[0] * refLen,
                  Egrp[i].from[1] * refLen,
                  Egrp[i].from[2] * refLen
                  );
          fprintf(fp,"\t\tTo                  : (%14.6e %14.6e %14.6e) / (%14.6e %14.6e %14.6e) [m]\n",
                  Egrp[i].to[0],
                  Egrp[i].to[1],
                  Egrp[i].to[2],
                  Egrp[i].to[0] * refLen,
                  Egrp[i].to[1] * refLen,
                  Egrp[i].to[2] * refLen
                  );
          break;
          
        case mon_DISC:
          fprintf(fp,"\t\tDisc\n");
          fprintf(fp,"\t\t# of Sample         : %d\n",Egrp[i].nSample);
          fprintf(fp,"\t\tRadius              : %14.6e / %14.6e [m]\n", Egrp[i].radius, Egrp[i].radius * refLen);
          fprintf(fp,"\t\tCenter              : (%14.6e %14.6e %14.6e) / (%14.6e %14.6e %14.6e) [m]\n",
                  Egrp[i].center[0],
                  Egrp[i].center[1],
                  Egrp[i].center[2],
                  Egrp[i].center[0] * refLen,
                  Egrp[i].center[1] * refLen,
                  Egrp[i].center[2] * refLen
                  );
          fprintf(fp,"\t\tNormal              : (%14.6e %14.6e %14.6e)\n",
                  Egrp[i].normal[0],
                  Egrp[i].normal[1],
                  Egrp[i].normal[2]);
          break;
      }
    } // nGrpEmit
    
    fprintf(fp, "\n\tTotal number of emission points : %ld\n", nEmitParticle);
    
    if ( restartFlag == 1 )
    {
      fprintf(fp,"\n\tRestart from %d step\n", restartStep);
    }
    
    fprintf(fp,"\n");
  } // Hostonly
}



//#############################################################################
// @brief タイミング測定区間にラベルを与える
// @note ffv側で登録しているので、これは未使用
void Cloud::set_timing_label()
{
  /*
   set_label("PT_Migration",         PerfMonitor::CALC, false);
   set_label("PT_ParticlePacking",   PerfMonitor::CALC);
   set_label("PT_EstablishCommPath", PerfMonitor::CALC);
   set_label("PT_CommParticle",      PerfMonitor::CALC);
   set_label("PT_ParticleUnpacking", PerfMonitor::CALC);
   set_label("PT_Tracking",          PerfMonitor::CALC);
   */
}


//#############################################################################
// @brief ascii output
bool Cloud::write_ascii(const unsigned step, const double time)
{
  if (nParticle == 0) return true;
  
  FILE* fp;
  char tmp_fname[50];
  
  sprintf( tmp_fname, "Particle");
  if ( !c_mkdir(tmp_fname) ) return false;
  
  sprintf( tmp_fname, "Particle/pt_%08ld_%06d.npt", step, myRank );
  
  if ( !(fp=fopen(tmp_fname, "w")) )
  {
    stamped_printf("\tSorry, can't open '%s' file. Write failed.\n", tmp_fname);
    return false;
  }
  
  fprintf(fp,"\n# step_time %ld %e\n", step, time);
  fprintf(fp,"# total_particles %ld\n", nParticle);
  fprintf(fp,"# no_of_chunks %d\n", (int)chunkList.size());
  
  int c=0;
  for(auto itr = chunkList.begin(); itr != chunkList.end(); ++itr)
  {
    fprintf(fp,"\n# Chunk %d\n", c++);
    (*itr)->write_ascii(fp, refLen, refVel);
  }
  
  fclose(fp);
  
  return true;
}


//#############################################################################
// @brief meta file output
bool Cloud::write_filelist(const unsigned step)
{
  if (nParticle == 0) return true;
  
  FILE* fp;
  char tmp_fname[50];
  sprintf( tmp_fname, "Particle/pt_files_%06d.lst", myRank );
  
  if ( !(fp=fopen(tmp_fname, "a")) )
  {
    stamped_printf("\tSorry, can't open 'pt_files_xxx.txt' file. Write failed.\n");
    return false;
  }
  
  sprintf( tmp_fname, "pt_%08ld_%06d.npt", step, myRank );
  fprintf(fp,"%s\n", tmp_fname);
  fclose(fp);
  
  return true;
}


// #################################################################
// @brief ディレクトリがなければ作成、既存なら何もしない（単一ディレクトリ）
// @param [in] path ディレクトリパス
bool Cloud::c_mkdir(const char* path)
{
  // 標準ライブラリ呼び出し
  // パーミッションはumaskで指定
  umask(022);
  
  int ret = mkdir(path, 0777); // rwx
  
  if ( 0 != ret )
  {
    // 既存以外のエラー
    if ( EEXIST != errno )
    {
      printf( "\tError(errno)=[%s]\n", strerror(errno) );
      return false;
    }
  }
  
  return true;
}


// #################################################################
// @brief binary output
bool Cloud::write_binary(const unsigned step, const double time)
{
  if (nParticle == 0) return true;
  
  char tmp_fname[50];
  
  sprintf( tmp_fname, "Particle");
  if ( !c_mkdir(tmp_fname) ) return false;
  
  sprintf( tmp_fname, "Particle/pt_%08ld_%06d.bpt", step, myRank );
  
  std::ofstream ofs(tmp_fname, std::ios::out | std::ios::binary);
  if (!ofs) {
    printf("\tCan't open '%s' file\n", tmp_fname);
    return false;
  }
  
  unsigned stp = step;
  double tm = time;
  unsigned csz = (unsigned)chunkList.size();
  
  ofs.write((char*)&stp, sizeof(unsigned));
  ofs.write((char*)&tm,  sizeof(double));
  ofs.write((char*)&nParticle, sizeof(unsigned));
  ofs.write((char*)&csz, sizeof(unsigned));
  
  
  for(auto itr = chunkList.begin(); itr != chunkList.end(); ++itr)
  {
    (*itr)->write_binary(ofs, refLen, refVel);
  }
  
  ofs.close();
  
  return true;
}


// #################################################################
// @brieaf リスタート時のランク番号情報
bool Cloud::readRestartRank()
{
  int flag = 0;
  char buf[20];
  int num_rank=0;
  int* rankno = NULL;
  FILE* fp;
  
  Hostonly_
  {
    if ( !(fp=fopen("./Particle/rpt_rank.txt", "r")) )
    {
      printf("\tSorry, can't open 'rpt_rank.txt' file.\n");
      flag = 1;
    }
  }
  
  int tmp = flag;
  if ( MPI_SUCCESS != MPI_Allreduce(&tmp,
                                    &flag,
                                    1,
                                    MPI_INT,
                                    MPI_MAX,
                                    MPI_COMM_WORLD) ) return false;
  
  if ( flag == 1 ) return false;
  
  
  Hostonly_ fscanf(fp, "%d %s", &num_rank, buf);
  
  if (MPI_Bcast(&num_rank, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS) {
    return false;
  }
  
  if ( num_rank < 1 ) return false;
  
  
  rankno = new int[num_rank];
  
  Hostonly_
  {
    for (int i=0; i<num_rank; i++)
    {
      fscanf(fp, "%d", &rankno[i]);
      //printf("restart rank = %d\n", rankno[i]);
    }
    
    fclose(fp);
  }
  
  if (MPI_Bcast(rankno, num_rank, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS) {
    return false;
  }
  
  for (int i=0; i< num_rank; i++)
  {
    if ( rankno[i] == myRank ) restartRankFlag = 1;
  }
  
  
  if ( rankno ) {
    delete [] rankno;
    rankno = NULL;
  }
  
  return true;
}


// #################################################################
// @brieaf リスタート時の粒子データ読み込み(Ascii)
bool Cloud::readRestartParticleAscii()
{
  int flag =  0;
  char buf[20];
  
  if ( restartRankFlag == 1 )
  {
    FILE* fp;
    char tmp_fname[80];
    
    sprintf( tmp_fname, "./Particle/pt_%08ld_%06d.npt", restartStep, myRank );
    
    if ( !(fp=fopen(tmp_fname, "r")) )
    {
      printf("\tSorry, can't open '%s' file.\n", tmp_fname);
      flag = 1;
    }
    
    int step;
    int npart;
    float time;
    int nchnk;
    
    if ( flag == 0 )
    {
      fscanf(fp, "%s %s %d %f", buf, buf, &step, &time);
      fscanf(fp, "%s %s %d", buf, buf, &npart);
      fscanf(fp, "%s %s %d", buf, buf, &nchnk);
      //printf("[%d] step=%d time=%f npart=%d nchunk=%d\n", myRank, step, time, npart, nchnk);
      
      for (int j=0; j<nchnk; j++)
      {
        particle p;
        Vec3r q, v, e;
        int lid, np, dd, sorg, semit, ch;
        int life;
        
        fscanf(fp, "%s %s %d", buf, buf, &ch);            // chunk no
        fscanf(fp, "%s %d", buf, &np);                    // particle
        fscanf(fp, "%s %d", buf, &dd);                    // emit_pnt_id
        fscanf(fp, "%s %e %e %e", buf, &e.x, &e.y, &e.z); // emit_origin
        fscanf(fp, "%s %d", buf, &sorg);                  // start_origin
        fscanf(fp, "%s %d", buf, &semit);                 // start_emit
        //printf("[%d] %d %d %e %e %e %d %d\n", myRank, ch, np, e.x, e.y, e.z, sorg, semit);
        
        for (int i=0; i<np; i++)
        {
          fscanf(fp, "%d %e %e %e %e %e %e %d %d",
                 &dd, &q.x, &q.y, &q.z, &v.x, &v.y, &v.z, &lid, &life);
          // lifeには前回の最後のライフカウントが入っている
          //printf("[%d] %e %e %e %d %e %e %e %d\n", myRank, q.x, q.y, q.z, life, v.x, v.y, v.z, lid);
          
          p.pos = q;
          p.bf  = Chunk::Activate(life);
          p.vel = v;
          p.foo = lid;
          
          addParticle2ChunkList(p);
        }
      }
    }
    
    fclose(fp);
    
  } // restartRankFlag
  
  
  int tmp = flag;
  if ( MPI_SUCCESS != MPI_Allreduce(&tmp,
                                    &flag,
                                    1,
                                    MPI_INT,
                                    MPI_MAX,
                                    MPI_COMM_WORLD) ) return false;
  
  return ( flag == 0 ) ? true : false;
}


// #################################################################
// @brieaf リスタート時の粒子データ読み込み(Binary)
bool Cloud::readRestartParticleBinary()
{
  int flag =  0;
  
  if ( restartRankFlag == 1 )
  {
    FILE* fp;
    char tmp_fname[80];
    
    sprintf( tmp_fname, "./Particle/restart_%08ld_%06d.brpt", restartStep, myRank );
    
    std::ifstream ifs(tmp_fname, std::ios::in | std::ios::binary);
    if (!ifs) {
      printf("\tCan't open %s file\n", tmp_fname);
      flag = 1;
    }
    
    unsigned step = 0;
    unsigned npart = 0;
    float time = 0.0;
    particle p;
    Vec3r q, v;
    int life, lid;
    
    if ( flag == 0 )
    {
      ifs.read((char*)&step, sizeof(unsigned));
      ifs.read((char*)&time,  sizeof(float));
      ifs.read((char*)&npart, sizeof(unsigned));
      fscanf(fp, "%d", &step);
      fscanf(fp, "%f", &time);
      fscanf(fp, "%d", &npart);
      
      for (int i=0; i<npart; i++)
      {
        ifs.read((char*)&q.x,  sizeof(float));
        ifs.read((char*)&q.y,  sizeof(float));
        ifs.read((char*)&q.z,  sizeof(float));
        ifs.read((char*)&v.x,  sizeof(float));
        ifs.read((char*)&v.y,  sizeof(float));
        ifs.read((char*)&v.z,  sizeof(float));
        ifs.read((char*)&lid,  sizeof(int));
        ifs.read((char*)&life, sizeof(int));
        
        p.pos = q;
        int b = life;
        p.bf  = Chunk::Activate(b);
        p.vel = v;
        p.foo = lid;
        
        addParticle2ChunkList(p);
      }
    }
    
    ifs.close();
    
  } // restartRankFlag
  
  
  int tmp = flag;
  if ( MPI_SUCCESS != MPI_Allreduce(&tmp,
                                    &flag,
                                    1,
                                    MPI_INT,
                                    MPI_MAX,
                                    MPI_COMM_WORLD) ) return false;
  
  return ( flag == 0 ) ? true : false;
}


// #################################################################
// @brieaf 粒子のchunkListへの追加
void Cloud::addParticle2ChunkList(particle p)
{
  // chunkList内で粒子IDをサーチし、存在しなければ、新たなチャンクを作る
  int check=0;
  for(auto itr = chunkList.begin(); itr != chunkList.end(); ++itr)
  {
    // IDがあれば、追加
    if ( p.foo == (*itr)->getUid() ) {
      (*itr)->addParticle(p);
      check=1;
      break;
    }
  }
  
  // uidのチャンクが存在しない場合
  if (check==0)
  {
    Chunk* m = new Chunk(p,
                         EmitStart,
                         myRank,
                         EmitInterval);
    chunkList.push_back(m);
  }
  
}
