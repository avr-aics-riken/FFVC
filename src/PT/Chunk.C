
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
    }
  }

  // 粒子管理情報
  for (int i=0; i<NDST; i++)
  {
    pInfo[4*i + 0] = c[i];  // 送信粒子数
    pInfo[4*i + 1] = grp;   // グループID
    pInfo[4*i + 2] = uid;   // 粒子ID
  }


  // 要素の削除
  for(auto itr = pchunk.begin(); itr != pchunk.end(); ++itr)
  {
    if ( IS_MIGRATE( (*itr).bf) ) {
      pchunk.erase(itr);
    }
  }

}


//#############################################################################
// @brief pchunkに保持している粒子を積分し、マイグレーションと寿命判定
// @param [in]  tr           Trackingクラスオブジェクトポインタ
// @param [in]  scheme       積分方法の指定
// @param [in]  life         粒子寿命
// @param [in]  dt           積分幅
// @param [out] max_particle 送受信バッファの計算に必要な送受信よ粒子数の最大値
// @retval マイグレーションが発生する場合にはtrue
bool Chunk::updatePosition(Tracking* tr,
                           const int scheme,
                           const int life,
                           const REAL_TYPE dt,
                           int& max_particle)
{
  int flag=0;
  int pcnt[NDST] = {};  // 行き先毎の送信粒子数

  for(auto itr = pchunk.begin(); itr != pchunk.end(); ++itr)
  {
    Vec3r p = (*itr).pos;
    Vec3i r;
    int c;
    //printf("p:(%14.6f %14.6f %14.6f)\n", p.x, p.y, p.z);

    switch(scheme)
    {
      case pt_euler:
        // 正 : 自領域外、マイグレーション
        // -1 : 計算領域外
        // -2 : 自領域内
        r = tr->integrate_Euler(p, dt);
        break;

      case pt_rk2:
        r = tr->integrate_RK2(p, dt);
        break;

      case pt_rk4:
        break;
    }

    if (r.x == -2) {
      ; // nothing
    }
    else if (r.x == -1 || r.y == -1 || r.z == -1) {
      Inactivate( (*itr).bf ); // 停止
    }
    else {
      stampMigrate((*itr).bf);         // マイグレーション候補
      c = encMigrateDir((*itr).bf, r); // 移動先をエンコード
      flag++;
      pcnt[c]++;
    }

    // 粒子位置情報をアップデート
    (*itr).pos = p;

    // アクティブな粒子のみライフタイムをインクリメント
    if ( IS_ACTIVE( (*itr).bf) ) (*itr).bf++;

    // 指定値を過ぎたら停止 >> 削除するか？
    if ( life < getBit26((*itr).bf) ) Inactivate((*itr).bf);
  }


  // 送信要素の最大値を保存
  // 行き先毎には異なるが、最大数で考えておく
  for (int i=0; i<NDST; i++) {
    if (max_particle < pcnt[i]) {
      max_particle = pcnt[i];
    }
  }


  if (flag > 0) return true;

  return false;
}


//#############################################################################
// @brief pchunkに開始点から粒子を追加
// @note ライフタイムはゼロ
void Chunk::addParticleFromOrigin()
{
  particle p;
  p.pos = origin;
  p.bf  = 0;
  p.foo = 0;
  p.vel.assign(0.0, 0.0, 0.0);
  Activate(p.bf);
  pchunk.push_front(p);
}


//#############################################################################
// @brief pchunkの終端に粒子を追加
// @note ライフタイムはparticleが持っている値を継承
void Chunk::addParticle(particle p)
{
  Activate(p.bf);
  pchunk.push_back(p);
}


//#############################################################################
// @brief ascii format出力
void Chunk::write_ascii(FILE* fp)
{
  fprintf(fp,"# particles   %d\n",(int)pchunk.size());
  fprintf(fp,"# group_id    %d\n", grp);
  fprintf(fp,"# particle_id %d\n", uid);
  fprintf(fp,"# origin      %e %e %e\n", origin.x, origin.y, origin.z);
  fprintf(fp,"# start_point %d\n", start_point);
  fprintf(fp,"# start_step  %d\n\n", start_step);
  
  for(auto itr = pchunk.begin(); itr != pchunk.end(); ++itr)
  {
    Vec3r p = (*itr).pos;
    Vec3r v = (*itr).vel;
    int b   = (*itr).bf;
    fprintf(fp,"%d %e %e %e %d %e %e %e %d\n",
               (b >> ACTIVE_BIT) & 0x1,
               p.x, p.y, p.z,
               getBit26(b),
               v.x, v.y, v.z,
               (*itr).foo
            );
  }
}
