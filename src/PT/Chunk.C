
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
// @param [in]  pbuf   粒子座標用送信バッファ
// @param [in]  fbuf   ビットフィールド用送信バッファ
// @param [in]  max_p  最大粒子個数
// @param [in]  pInfo  送受信情報
// @note bszは毎回異なる
void Chunk::packParticle(REAL_TYPE* ps_buf,
                         int* bs_buf,
                         const int max_p,
                         int* pInfo)
{
  // 毎回クリア
  int c[NDST] = {}; // 0
  memset(pInfo, 0, sizeof(int)*NDST*4);
  
  // マイグレーション粒子がない場合には処理をスキップ
  if ( !mig ) return;

    
    
  int bsz6 = max_p * 6; // 6要素
  int bsz2 = max_p * 2; // 2要素

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
      //printf("pack : %f %f %f %d %f %f %f %d\n",
      //       (*itr).pos.x, (*itr).pos.y, (*itr).pos.z, b, (*itr).vel.x, (*itr).vel.y, (*itr).vel.z, (*itr).foo);

    }
  }

  // 粒子管理情報
  for (int i=0; i<NDST; i++)
  {
    pInfo[4*i + 0] = c[i];  // 送信粒子数
    pInfo[4*i + 1] = grp;   // グループID
    pInfo[4*i + 2] = uid;   // 開始点ID
  }


  // 要素の削除
  for(auto itr = pchunk.begin(); itr != pchunk.end(); ++itr)
  {
    if ( IS_MIGRATE( (*itr).bf) ) {
      pchunk.erase(itr);
    }
  }
    
  // フラグをfalseに戻す
  mig = false;

}


//#############################################################################
// @brief pchunkに保持している粒子を積分し、マイグレーションと寿命判定
// @param [in]  tr           Trackingクラスオブジェクトポインタ
// @param [in]  scheme       積分方法の指定
// @param [in]  life         粒子生存カウント
// @param [in]  dt           積分幅
// @param [out] max_particle 送受信バッファの計算に必要な送受信よ粒子数の最大値
void Chunk::updatePosition(Tracking* tr,
                           const int scheme,
                           const int life,
                           const REAL_TYPE dt,
                           int& max_particle)
{
  int pcnt[NDST] = {};  // 行き先毎の送信粒子数
  int flag = 0;

  for(auto itr = pchunk.begin(); itr != pchunk.end(); ++itr)
  {
    Vec3r p = (*itr).pos;
    Vec3r v = (*itr).vel;
    Vec3i r;
    int c;
    
    if ( IS_ACTIVE( (*itr).bf) )
    {
      switch(scheme)
      {
        case pt_euler:
          // 正 : 自領域外、マイグレーション
          // -1 : 計算領域外
          // -2 : 自領域内
          r = tr->integrate_Euler(p, v, dt);
          break;
          
        case pt_rk2:
          r = tr->integrate_RK2(p, v, dt);
          break;
          
        case pt_rk4:
          break;
      }
      //printf("r= %d %d %d\n", r.x, r.y, r.z);
      if (r.x == -2) {
        ; // nothing
      }
      else if (r.x == -1 || r.y == -1 || r.z == -1) {
        (*itr).bf = Inactivate( (*itr).bf ); // 停止
      }
      else {
        (*itr).bf = stampMigrate( (*itr).bf );  // マイグレーション候補
        c = encMigrateDir((*itr).bf, r);        // 移動先をエンコード
        pcnt[c]++;
        flag++;
        //printf("migrate from %d > %d (%d %d %d)\n", myRank, c, r.x, r.y, r.z);
      }
      
      // 粒子位置情報をアップデート
      (*itr).pos = p;
      (*itr).vel = v;
      (*itr).bf++;
      //printf("%f %f %f %d : %d %d %d\n", v.x, v.y, v.z, flag, r.x, r.y, r.z);
      
      // 寿命制御があるとき、指定値を過ぎたら停止 >> 削除するか？
      if (life > 0)
      {
        if ( life < getBit26( (*itr).bf) ) (*itr).bf = Inactivate( (*itr).bf );
      }
      
    } // IS_ACTIVE
  } // itr=pchunk


  // 送信要素の最大値を保存
  // 行き先毎には異なるが、最大数で考えておく
  for (int i=0; i<NDST; i++) {
    if (max_particle < pcnt[i]) {
      max_particle = pcnt[i];
    }
  }
  
  // packParticle()での番兵
  if ( flag > 0 ) mig = true;
  
}


//#############################################################################
// @brief pchunkに開始点から粒子を追加
// @note ライフカウントはゼロ
void Chunk::addParticleFromOrigin()
{
  particle p;
  p.pos = EmitOrigin;
  p.bf  = 0;
  p.foo = 0;
  p.vel.assign(0.0, 0.0, 0.0);
  p.bf = Activate(p.bf);
  pchunk.push_front(p);
}


//#############################################################################
// @brief ascii format出力
void Chunk::write_ascii(FILE* fp)
{
  fprintf(fp,"particles %d\n",(int)pchunk.size());
  fprintf(fp,"group_id %d\n", grp);
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
               p.x, p.y, p.z,
               getBit26(b),
               v.x, v.y, v.z,
               foo
            );
  }
  fprintf(fp,"\n");
}
