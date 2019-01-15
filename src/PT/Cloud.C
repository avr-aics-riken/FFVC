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
  int max_part=0;  // 送受信バッファの計算に必要な送受信最大粒子数
  bool ret;
  
  TIMING_start("PT_Tracking");
  
  for(auto itr = chunkList.begin(); itr != chunkList.end(); ++itr)
  {
    int m = (*itr)->getGrp();
    if (m<0) {
      printf("Error : grp number\n");
      return false;
    }
    
    // 指定時刻が過ぎ、処理ステップになった場合
    if ( Interval[m].isStarted(step, time) && Interval[m].isTriggered(step, time) )
    {
      // 保持している粒子を積分し、マイグレーション準備
      (*itr)->updatePosition(tr, scheme, Egrp[m].getLife(), dt, max_part);
      
      // 登録した開始点から粒子を追加
      (*itr)->addParticleFromOrigin();
    }
  }
  
  // バッファ長の更新が必要な場合、バッファ計算用の粒子数はBUF_UNIT単位で確保
  if (buf_max_particle < max_part) {
    buf_max_particle = (max_part / BUF_UNIT + 1) * BUF_UNIT;
    buf_updated = true;
  }
  
  // バッファ長を同期
  ret = PC.migrateBuffer(buf_max_particle);

  TIMING_stop("PT_Tracking");
  if ( !ret ) return false;
  

  TIMING_start("PT_Migration");
  ret = migration();
  TIMING_stop("PT_Migration");
  if ( !ret ) return false;

  
  
  // 出力の場合
  if ( ((step/log_interval)*log_interval == step ) ||
       ((step/file_interval)*file_interval == step ) )
  {
    nParticle = 0;
    for(auto itr = chunkList.begin(); itr != chunkList.end(); ++itr) {
      nParticle += (unsigned)(*itr)->getNpoints();
    }
  }
  
  
  // ログ出力
  if ( (step/log_interval)*log_interval == step )
  {
    TIMING_start("PT_Statistics");
    ret = PC.Statistics(nCommParticle, nParticle, gParticle);
    Hostonly_ logging(step);
    TIMING_stop("PT_Statistics");
    if ( !ret ) return false;
  }
  
  
  // ファイル出力
  if ( (step/file_interval)*file_interval == step )
  {
    TIMING_start("PT_fileIO");
    ret = write_ascii(step, time);
    write_filelist(step);
    TIMING_stop("PT_fileIO");
    if ( !ret ) return false;
  }
  
  
  return true;
}



//#############################################################################
// @brief 初期設定
bool Cloud::initCloud(FILE* fp)
{
  if ( !setPTinfo() ) return false;
  
  if ( !determineUniqueID() ) return false;
  
  displayParam(stdout);
  displayParam(fp);
  
  tr = new Tracking(size,
                    guide,
                    origin,
                    region,
                    pitch,
                    vSource,
                    bcd,
                    myRank);
  
  // ログ出力
  if (log_interval > 0) {
    Hostonly_
    {
      if ( !(fpl=fopen("pt_log.txt", "w")) )
      {
        stamped_printf("\tSorry, can't open 'pt_log.txt' file. Write failed.\n");
        return false;
      }
    }
  }

  // 通信クラスの設定
  PC.setPtComm(G_division, myRank, numProc);

  
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
  bool ret;

  // マイグレーション用に確保したバッファ領域の準備
  if (buf_updated) {
    PC.resizeSendBuffer(buf_max_particle);
    buf_updated = false;
  }

  // バッファサイズの更新がない場合は再初期化のみ
  PC.initSendBuffer();


  // マイグレーション候補を見つけ、行き先毎にパック
  TIMING_start("PT_ParticlePacking");
  for(auto itr = chunkList.begin(); itr != chunkList.end(); ++itr)
  {
    (*itr)->packParticle(PC.ps_ptr(),
                         PC.bs_ptr(),
                         buf_max_particle,
                         PC.pInfo_ptr());
  }
  TIMING_stop("PT_ParticlePacking");

  // 周囲のランクと通信し、経路と送受信データ数を確定
  TIMING_start("PT_EstablishCommPath");
  ret = PC.establishCommPath();
  TIMING_stop("PT_EstablishCommPath");
  if ( !ret ) return false;

  // データ本体の送受信
  TIMING_start("PT_CommParticle");
  ret = PC.commParticle(buf_max_particle);
  TIMING_stop("PT_CommParticle");
  if ( !ret ) return false;


  // アンパック
  TIMING_start("PT_ParticleUnpacking");
  unpackParticle();
  TIMING_stop("PT_ParticleUnpacking");
  

  // マイグレーション終了
  flag_migration = false;

  return true;
}


//#############################################################################
// @brief 受信粒子のアンパック
void Cloud::unpackParticle()
{
  int* pInfo = PC.pInfo_ptr();
  REAL_TYPE* pr_buf = PC.pr_ptr();
  int* br_buf = PC.br_ptr();
  const int bsz = buf_max_particle;

  for (int i=0; i<NDST; i++)
  {
    const int gid = pInfo[4*i+1];   // グループID
    const int pid = pInfo[4*i+2];   // 粒子ID
    const int cnt = pInfo[4*i+3] ;  // 受信粒子数
    
    int b, foo;
    Vec3r pos;
    Vec3r v;
    particle p;

    for (int j=0; j<cnt; j++)
    {
      pos.x = pr_buf[bsz*i + 6*j+0];
      pos.y = pr_buf[bsz*i + 6*j+1];
      pos.z = pr_buf[bsz*i + 6*j+2];
      v.x   = pr_buf[bsz*i + 6*j+3];
      v.y   = pr_buf[bsz*i + 6*j+4];
      v.z   = pr_buf[bsz*i + 6*j+5];
      b     = br_buf[bsz*i + 2*j+0];
      foo   = br_buf[bsz*i + 2*j+1];
      b = Chunk::removeMigrate(b);

      p.pos = pos;
      p.bf  = b;
      p.vel = v;
      foo++;        // マイグレーションの回数
      p.foo = foo;
      //printf("unpack : %f %f %f %d %d %f %f %f %d\n",
      //       pos.x, pos.y, pos.z, br_buf[bsz*i + 2*j+0], b, v.x, v.y, v.z, foo);

      // chunkList内で粒子IDをサーチし、存在しなければ、新たなチャンクを作る
      int check=0;
      for(auto itr = chunkList.begin(); itr != chunkList.end(); ++itr)
      {
        // IDがあれば、追加
        if ( pid == (*itr)->getUid() ) {
          (*itr)->addParticle(p);
          check=1;
          break;
        }
      }

      // pidのチャンクが存在しない場合
      if (check==0)
      {
        REAL_TYPE v[3];
        Egrp[gid].getEmitOrg(v);
        Vec3r eo(v);
        
        Chunk* m = new Chunk(p,
                             gid,
                             pid,
                             Egrp[gid].getStart(),
                             myRank,
                             Egrp[gid].getInterval(),
                             eo);
        chunkList.push_back(m);
        Egrp[gid].incGroup();
      }

    }
  } // NDST-loop
  
}


//#############################################################################
// @brief 既に存在するグループIDの検索
// @param [in] c   検索対象ID
// @retval 存在すれば true
bool Cloud::searchGrp(const int c)
{
  for (int i=0; i<nGrpEmit; i++)
  {
    if ( Egrp[i].isEns(c) ) return true;
  }
  return false;
}



//#############################################################################
// @brief ログ出力
void Cloud::logging(const unsigned step)
{
  // マイグレーションで移動した粒子数と全粒子数
  fprintf(fpl, "step %ld migrated %d total %ld\n", step, nCommParticle, gParticle);
  
  /*
  for (int i=0; i<numProc; i++) {
    fprintf(fpl, "%d ", i);
  }
  fprintf(fpl,"\n");
  */
  
  // 各ランク毎の粒子数
  unsigned* p = PC.nPart_ptr();
  for (int i=0; i<numProc; i++) {
    fprintf(fpl, "%ld ", p[i]);
  }
  fprintf(fpl,"\n");
  
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
    (*itr)->setUid(acc[myRank]+nc);
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
  
  // ログ出力インターバル
  label = label_base + "/log_interval";
  
  if ( !(tpCntl->getInspectedValue(label, f_val )) )
  {
    Hostonly_ stamped_printf("\tParsing warning : No Label in '%s'\n", label.c_str());
    return false;
  }
  log_interval = (int)f_val;

  
  // ファイル出力インターバル
  label = label_base + "/file_interval";
  
  if ( !(tpCntl->getInspectedValue(label, f_val )) )
  {
    Hostonly_ stamped_printf("\tParsing warning : No Label in '%s'\n", label.c_str());
    return false;
  }
  file_interval = (int)f_val;
  
  
  // ファイル出力フォーマット
  label = label_base + "/file_format";
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
    Hostonly_ stamped_printf("\tParsing warning : No Label in '%s'\n", label.c_str());
    return false;
  }
  
  if ( !strcasecmp(str.c_str(), "ascii") ) {
    file_format = 0;
  }
  else if ( !strcasecmp(str.c_str(), "binary") ) {
    file_format = 1;
  }
  else {
    Hostonly_ printf("Invalid file format\n");
    return false;
  }
  


  // 粒子放出開始点の個数のチェック
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
      printf("\tParsing error : No No List[@]\n");
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

  // IntervalManagerのインスタンス
  Interval = new IntervalManager[nGrpEmit];

  for (int j=0; j<nGrpEmit; j++)
  {
    Interval[j].setMode(IntervalManager::By_step);
  }



  // リストの読み込み
  label_base = "/ParticleTracking";

  nlist=0;
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


    if ( !getTPparam(label_leaf, nlist) ) return false;


    for (int j=0; j<nGrpEmit; j++)
    {
      Interval[j].setStart(m_start);
    }


    // 順番に粒子IDとチャンクIDを付与
    // 粒子IDは各サブドメイン毎の仮のもの
    switch(mon_type)
    {
      case mon_POINT_SET:
        Egrp[nlist].setType(mon_POINT_SET);
        setPointset(label_leaf, nlist);
        break;

      case mon_LINE:
        Egrp[nlist].setType(mon_LINE);
        setLine(label_leaf, nlist);
        break;

      case mon_DISC:
        Egrp[nlist].setType(mon_DISC);
        setDisc(label_leaf, nlist);
        break;

      default:
        return false;
        break;
    }

    nlist++;
  } // nnode


  for (int j=0; j<nGrpEmit; j++)
  {
    if ( !Interval[j].initTrigger(m_start,
                                 (double)m_start * (double)dt,
                                 (double)dt) ) return false;
  }

  nParticle = chunkList.size();

  return true;
}


//#############################################################################
// @brief 粒子追跡パラメータ取得
bool Cloud::getTPparam(const string label_leaf, int odr)
{
  string str, label;
  double f_val=0.0;


  // 粒子出力
  label = label_leaf + "/StartEmit";

  if ( !(tpCntl->getInspectedValue(label, f_val )) )
  {
    Hostonly_ stamped_printf("\tParsing error : Invalid string for '%s'\n", label.c_str());
    return false;
  }
  else
  {
    Egrp[odr].setStart( (int)f_val );
  }

  label = label_leaf + "/Interval";

  if ( !(tpCntl->getInspectedValue(label, f_val )) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    return false;
  }
  else
  {
    Egrp[odr].setInterval( (int)f_val );
    Interval[odr].setInterval(f_val);
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


  label =  label_leaf + "/lifetime";

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
    Egrp[odr].setLife(tmp);
  }

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
    
    
    Egrp[odr].setEmitOrg(v);

    // 自領域内であれば、初期開始点として追加
    if ( !ModeTOOL ) {
      if ( inOwnRegion(v) )
      {
        Vec3r pos(v);
        Chunk* m = new Chunk(pos,
                             odr,
                             1,
                             Egrp[odr].getStart(),
                             myRank,
                             Egrp[odr].getInterval());
        chunkList.push_back(m);
        Egrp[odr].incGroup();
      }
    }
    
  }
  return true;
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
  
  REAL_TYPE ooo[3] = {dd.x, dd.y, dd.z};
  Egrp[odr].setEmitOrg(ooo);

  for (int m = 0; m < npnt; m++)
  {
    Vec3r pos = st + dd * (REAL_TYPE)m;

    printf("[%d] %d  tgt = (%14.6e %14.6e %14.6e)\n", myRank, m, pos.x, pos.y, pos.z);

    // 自領域内であれば、初期開始点として追加
    if ( !ModeTOOL ) {
      if ( inOwnRegion(pos) )
      {
        Chunk* m = new Chunk(pos,
                             odr,
                             1,
                             Egrp[odr].getStart(),
                             myRank,
                             Egrp[odr].getInterval());
        chunkList.push_back(m);
        Egrp[odr].incGroup();
      }
    }
    
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

  // load parameter of 'from' and 'to'
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

  // 入力パラメータの次元が有次元のとき，無次元化する
  if (unit == DIMENSIONAL) normalizeCord(nv);

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
  Egrp[odr].radius = radius;


  // 点群の発生と登録
  samplingInCircle(cnt, nv, radius, nSample, odr);

  return true;
}


//#############################################################################
/// @brief 半径r内のサンプリング
/// @param [in]  cnt        中心座標
/// @param [in]  nv         法線方向
/// @param [in]  radius     半径
/// @param [in]  nSample    サンプル数
/// @param [in]  odr        グループ登録インデクス
/// @note 指定半径より少し内側（98%）にする
/// @url https://teramonagi.hatenablog.com/entry/20141113/1415839321
void Cloud::samplingInCircle(const REAL_TYPE* cnt,
                             const REAL_TYPE* nv,
                             const REAL_TYPE radius,
                             const int nSample,
                             const int odr)
{
  REAL_TYPE theta, r, x, y;
  double pi = 2.0*asin(1.0);
  Vec3r q, t;
  Vec3r center(cnt);
  Vec3r normal(nv);

  Vec3r angle = getAngle(normal);

  for (int i=0; i<nSample; i++)
  {
    theta = mts(0.0, 2.0*pi);
    r = sqrt( 2.0*mts(0.0, (double)radius) );

    x = 0.98 * r * cos(theta);
    y = 0.98 * r * sin(theta);

    q = rotate_inv(angle, t.assign(x, y, 0.0)) + center;
    
    REAL_TYPE ooo[3] = {q.x, q.y, q.z};
    Egrp[odr].setEmitOrg(ooo);

    if ( !ModeTOOL ) {
      if ( inOwnRegion(q) )
      {
        Vec3r pos(q);
        Chunk* m = new Chunk(pos,
                             odr,
                             1,
                             Egrp[odr].getStart(),
                             myRank,
                             Egrp[odr].getInterval());
        chunkList.push_back(m);
        Egrp[odr].incGroup();
      }
    }
  }
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
  Vec3r dir;         ///< 矩形の方向規定の参照ベクトル

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

    for (int i=0; i<nGrpEmit; i++)
    {
      fprintf(fp,"\n\tGroup               : %s\n",Egrp[i].getGrpName().c_str());
      fprintf(fp,"\t\tStart step          : %d\n",Egrp[i].getStart());
      fprintf(fp,"\t\tInterval            : %d\n",Egrp[i].getInterval());
      fprintf(fp,"\t\tLife time           : %d\n",Egrp[i].getLife());
      fprintf(fp,"\t\tType                : ");
      switch (Egrp[i].getType()) {
        case mon_POINT_SET:
          fprintf(fp,"\t\tPointset\n");
          break;

        case mon_LINE:
          fprintf(fp,"\t\tLine\n");
          fprintf(fp,"\t\t# of division       : %d\n",Egrp[i].nDiv);
          fprintf(fp,"\t\tFrom                : (%14.6e %14.6e %14.6e)\n",
                                        Egrp[i].from[0],
                                        Egrp[i].from[1],
                                        Egrp[i].from[2]);
          fprintf(fp,"\t\tTo                  : (%14.6e %14.6e %14.6e)\n",
                                        Egrp[i].to[0],
                                        Egrp[i].to[1],
                                        Egrp[i].to[2]);
          break;

        case mon_DISC:
          fprintf(fp,"\t\tDisc\n");
          fprintf(fp,"\t\t# of Sample         : %d\n",Egrp[i].nSample);
          fprintf(fp,"\t\tRadius              : %14.6e\n",Egrp[i].radius);
          fprintf(fp,"\t\tCenter              : (%14.6e %14.6e %14.6e)\n",
                                        Egrp[i].center[0],
                                        Egrp[i].center[1],
                                        Egrp[i].center[2]);
          fprintf(fp,"\t\tNormal              : (%14.6e %14.6e %14.6e)\n",
                                        Egrp[i].normal[0],
                                        Egrp[i].normal[1],
                                        Egrp[i].normal[2]);
          break;
      }

      fprintf(fp,"\n");
    }
  }
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
  sprintf( tmp_fname, "Particle/pt_%08ld_%06d.npt", step, myRank );
  
  if ( !(fp=fopen(tmp_fname, "w")) )
  {
    stamped_printf("\tSorry, can't open 'pt_log.txt' file. Write failed.\n");
    return false;
  }
  
  
  fprintf(fp,"\n# step_time %ld %e\n", step, time);
  fprintf(fp,"# total_particles %ld\n", nParticle);
  fprintf(fp,"# no_of_chunks %d\n", (int)chunkList.size());
  
  int c=0;
  for(auto itr = chunkList.begin(); itr != chunkList.end(); ++itr)
  {
    fprintf(fp,"\n# Chunk %d\n", c++);
    (*itr)->write_ascii(fp);
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
    stamped_printf("\tSorry, can't open 'pt_file_xxx.txt' file. Write failed.\n");
    return false;
  }
  
  sprintf( tmp_fname, "pt_%08ld_%06d.npt", step, myRank );
  fprintf(fp,"%s\n", tmp_fname);
  fclose(fp);
  
  return true;
}
