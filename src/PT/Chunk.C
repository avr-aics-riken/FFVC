
//##################################################################################
//
// Copyright (c) 2016-2019 Research Institute for Information Technology, Kyushu university
// All rights researved.
//
//##################################################################################

/**
 * @file   Chunk.C
 * @brief  Chunk class
 * @author riit
 * @note 粒子単位の軌跡を管理
 */

#include "Chunk.h"


//#############################################################################
// @brief マイグレーションフラグの立っている粒子をパックし、リストから削除
// @param [in]     pbuf      粒子座標用送信バッファ
// @param [in]     fbuf      ビットフィールド用送信バッファ
// @param [in]     buf_size  最大粒子個数
// @param [in,out] c[NDST]   各方向の送信粒子数の全チャンクの総和
// @note bszは毎回異なる
void Chunk::packParticle(REAL_TYPE* ps_buf,
                         int* bs_buf,
                         const unsigned buf_size,
												 int* c)
{
  // マイグレーション粒子がない場合には処理をスキップ
  if ( !MigrationFlag ) return;

	// バッファには各チャンクのオフセット必要
    
  // 周囲26方向のバッファの区切り
  // バッファを長さbuf_sizeで26領域に分け、各領域の先頭アドレスから使う
  size_t bsz6 = buf_size * 6; // 6要素
  size_t bsz2 = buf_size * 2; // 2要素

  // 粒子データのパック
  for(auto itr = pchunk.begin(); itr != pchunk.end(); ++itr)
  {
    int b = (*itr).bf;
    int m = decMigrateDir(b); // m=[0,26]

    if ( IS_MIGRATE(b) ) {
      ps_buf[bsz6*m + 6*c[m]+0] = (*itr).pos.x;
      ps_buf[bsz6*m + 6*c[m]+1] = (*itr).pos.y;
      ps_buf[bsz6*m + 6*c[m]+2] = (*itr).pos.z;
      ps_buf[bsz6*m + 6*c[m]+3] = (*itr).vel.x;
			ps_buf[bsz6*m + 6*c[m]+4] = (*itr).vel.y;
      ps_buf[bsz6*m + 6*c[m]+5] = (*itr).vel.z;
      bs_buf[bsz2*m + 2*c[m]+0] = b;
      bs_buf[bsz2*m + 2*c[m]+1] = (*itr).foo;
			c[m]++;
      //printf("migrated %d to %d\n", c[m], m);
      //printf("pack : %f %f %f %d %d %f %f %f %d\n",
      //       (*itr).pos.x, (*itr).pos.y, (*itr).pos.z, b, getBit25(clrMigrateDir(b)), (*itr).vel.x, (*itr).vel.y, (*itr).vel.z, (*itr).foo);
    }
  }


  // 要素の削除
  for(auto itr = pchunk.begin(); itr != pchunk.end();)
  {
    if ( IS_MIGRATE( (*itr).bf) )
		{
      itr = pchunk.erase(itr);
			continue;
    }
		++itr;
  }
    
  // フラグをfalseに戻す
  MigrationFlag = false;

}


//#############################################################################
// @brief pchunkに保持している粒子を積分し、マイグレーションと寿命判定
// @param [in]  tr           Trackingクラスオブジェクトポインタ
// @param [in]  scheme       積分方法の指定
// @param [in]  life         粒子生存カウント
// @param [in]  dt           積分幅
// @param [in]  Rmap         3ｘ3の隣接ランクマップ
void Chunk::updatePosition(Tracking* tr,
                           const int scheme,
                           const int life,
                           const REAL_TYPE dt,
                           const int* Rmap)
{
  int pcnt[NDST] = {0};  // 行き先毎の送信粒子数
  int flag = 0;

  for(auto itr = pchunk.begin(); itr != pchunk.end(); ++itr)
  {
    Vec3r p = (*itr).pos;
    Vec3r v = (*itr).vel;
    Vec3i r;
    int c;
		bool wallFlag = false;
    
    if ( IS_ACTIVE( (*itr).bf) )
    {
      switch(scheme)
      {
        case pt_euler:
          r = tr->integrate_Euler(p, v, wallFlag, dt);
          break;
          
        case pt_rk2:
          r = tr->integrate_RK2(p, v, wallFlag, dt);
          break;
          
        case pt_rk4:
          break;
      }

      // d=[0,26], r = {-1, 0, +1}
      int d = getMigrateDir(r);
      
      if (d == 13) // 自領域 nothing
      {
        ;
      }
      else if (Rmap[d] == -1) // 計算領域外
      {
        (*itr).bf = Inactivate( (*itr).bf ); // 停止
        printf("Out of region %d %d %d\n", r.x, r.y, r.z);
      }
      else // Migration
      {
        //printf("b = %d\n", (*itr).bf );
        (*itr).bf = stampMigrate( (*itr).bf );  // マイグレーション候補
        //printf("stampMigrate(b) = %d\n", (*itr).bf );
        c = encMigrateDir((*itr).bf, r);        // 移動先をエンコード
        //printf("encMigrateDir(b) = %d\n", (*itr).bf );
        pcnt[c]++;
        flag++;
        //printf("migrate from %d > %d (%d %d %d)\n", myRank, c, r.x, r.y, r.z);
      }
			
			// 壁面を通過した場合
			if (wallFlag)
			{
				(*itr).bf = Inactivate( (*itr).bf ); // 停止
				printf("wall passing\n");
			}

      
      // 粒子位置情報をアップデート
      (*itr).pos = p;
      (*itr).vel = v;
      (*itr).bf++;
      //printf("inc(b) = %d\n", (*itr).bf );
      //printf("%f %f %f %d : %d %d %d\n", v.x, v.y, v.z, flag, r.x, r.y, r.z);
      
      // 寿命制御があるとき、指定値を過ぎたら停止 >> 削除するか？
      if (life > 0)
      {
        if ( life < getBit25( (*itr).bf) )
				{
					(*itr).bf = Inactivate( (*itr).bf );
					printf("life time\n");
				}
      }
      
    } // IS_ACTIVE
  } // itr=pchunk

	
	// 行き先毎の送信粒子数
	for (int i=0; i<NDST; i++)
	{
		pSend[i] = pcnt[i];
	}
  
  // packParticle()での番兵
	if ( flag > 0 ) {
		//printf("Migration [%d] %d\n", myRank, flag);
		MigrationFlag = true;
	}
  
}


//#############################################################################
// @brief pchunkに開始点から粒子を追加
// @note マイクレーション先は除く
void Chunk::addParticleFromOrigin()
{
	if ( startOrigin == 1 )
	{
		particle p;
		p.pos = EmitOrigin;
		p.bf  = 0;
		p.foo = uid;
		p.vel.assign(0.0, 0.0, 0.0);
		p.bf = Activate(p.bf);
		pchunk.push_front(p);
	}
}


//#############################################################################
// @brief pchunkのuidをセットする
// @param [in] m_uid   uid
void Chunk::setPchunkUid(const int m_uid)
{
	for(auto itr = pchunk.begin(); itr != pchunk.end(); ++itr)
	{
		(*itr).foo = m_uid;
	}
}



//#############################################################################
// @brief ascii format出力
void Chunk::write_ascii(FILE* fp, const REAL_TYPE refL, const REAL_TYPE refV)
{
  fprintf(fp,"particles %d\n",(int)pchunk.size());
  fprintf(fp,"emit_pnt_id %d\n", uid);
  fprintf(fp,"emit_origin %e %e %e\n", EmitOrigin.x, EmitOrigin.y, EmitOrigin.z);
  fprintf(fp,"start_origin %d\n", startOrigin);
  fprintf(fp,"start_emit %d\n", EmitStep);
  
  for(auto itr = pchunk.begin(); itr != pchunk.end(); ++itr)
  {
    Vec3r p = (*itr).pos;
    Vec3r v = (*itr).vel;
    int b   = (*itr).bf;
    int foo = (*itr).foo;
    fprintf(fp,"%d %e %e %e %d %e %e %e %d\n",
               BIT_SHIFT(b, ACTIVE_BIT),
               p.x * refL,
               p.y * refL,
               p.z * refL,
               getBit25(b),
               v.x * refV,
               v.y * refV,
               v.z * refV,
               foo
            );
  }
  fprintf(fp,"\n");
}


//#############################################################################
// @brief binary format出力
void Chunk::write_binary(std::ofstream &ofs, const REAL_TYPE refL, const REAL_TYPE refV)
{
	unsigned n = (unsigned)pchunk.size();
	unsigned est = (unsigned)EmitStep;
	
	ofs.write((char*)&n, sizeof(unsigned));
	ofs.write((char*)&uid, sizeof(unsigned));
	ofs.write((char*)&EmitOrigin.x, sizeof(REAL_TYPE));
	ofs.write((char*)&EmitOrigin.y, sizeof(REAL_TYPE));
	ofs.write((char*)&EmitOrigin.z, sizeof(REAL_TYPE));
	ofs.write((char*)&est, sizeof(unsigned));
	
	for(auto itr = pchunk.begin(); itr != pchunk.end(); ++itr)
	{
		Vec3r p = (*itr).pos * refL;
		Vec3r v = (*itr).vel * refV;
		int b   = (*itr).bf;
		int uid = (*itr).foo;
		int actv= BIT_SHIFT(b, ACTIVE_BIT);
		unsigned lc = getBit25(b);
		
		ofs.write((char*)&actv, sizeof(int));
		ofs.write((char*)&p.x, sizeof(REAL_TYPE));
		ofs.write((char*)&p.y, sizeof(REAL_TYPE));
		ofs.write((char*)&p.z, sizeof(REAL_TYPE));
		ofs.write((char*)&lc, sizeof(unsigned));
		ofs.write((char*)&v.x, sizeof(REAL_TYPE));
		ofs.write((char*)&v.y, sizeof(REAL_TYPE));
		ofs.write((char*)&v.z, sizeof(REAL_TYPE));
		ofs.write((char*)&uid, sizeof(int));
	}
}
