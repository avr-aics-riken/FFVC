
//##################################################################################
//
// Flow Base class
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2015 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2016 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
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

// @brief マイグレーションフラグの立っている粒子をパックし、リストから削除
// @param [in]  pbuf   粒子座標用送信バッファ
// @param [in]  fbuf   ビットフィールド用送信バッファ
// @param [in]  bsz    バッファサイズ
// @param [in]  pInfo  送受信情報
// @note bszは毎回異なる
void Chunk::packParticle(REAL_TYPE* ps_buf,
                         int* bs_buf,
                         const int bsz,
                         int* pInfo)
{
  // 毎回クリア
  int c[NDST] = {}; // 0
  memset(pInfo, 0, sizeof(int)*NDST*4);

  // 粒子データのパック
  for(auto itr = pchunk.begin(); itr != pchunk.end(); ++itr)
  {
    int   b = (*itr).bf;
    Vec3r p = (*itr).pos;
    int m = decMigrateDir(b); // m=[0,26]

    if ( IS_MIGRATE(b) ) {
      ps_buf[bsz*m + 3*c[m]+0] = p.x;
      ps_buf[bsz*m + 3*c[m]+1] = p.y;
      ps_buf[bsz*m + 3*c[m]+2] = p.z;
      bs_buf[bsz*m + c[m]] = b;
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


// @brief pchunkに保持している粒子を積分し、マイグレーションと寿命判定
// @param [in]  tr          Trackingクラスオブジェクトポインタ
// @param [in]  scheme      積分方法の指定
// @param [in]  life        粒子寿命
// @param [in]  dt          積分幅
// @param [out] buf_len     送受信バッファに必要な長さ
// @retval マイグレーションが発生する場合にはtrue
bool Chunk::updatePosition(Tracking* tr,
                           const int scheme,
                           const int life,
                           const REAL_TYPE dt,
                           int& buf_len)
{
  int flag=0;
  int pcnt[NDST] = {0};  // 行き先毎の送信粒子数

  for(auto itr = pchunk.begin(); itr != pchunk.end(); ++itr)
  {
    Vec3r p = (*itr).pos;
    Vec3i r;
    int c;

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
      stampMigrate((*itr).bf);        // マイグレーション候補
      c = encMigrateDir((*itr).bf, r); // 移動先をエンコード
      flag++;
      pcnt[c]++;
    }

    // 粒子位置情報をアップデート
    (*itr).pos = p;

    // アクティブな粒子のみライフタイムをインクリメント
    if ( IS_ACTIVE( (*itr).bf) ) (*itr).bf++;

    // 指定値を過ぎたら停止
    if ( life < getBit26((*itr).bf) ) Inactivate((*itr).bf);
  }


  // 送信バッファの最大値を保存
  // 行き先毎には異なるが、最大長で考えておく
  for (int i=0; i<NDST; i++) {
    if (buf_len < pcnt[i]) {
      buf_len = pcnt[i];
    }
  }


  if (flag > 0) return true;

  return false;
}


// @brief pchunkに開始点から粒子を追加
// @note ライフタイムはゼロ
void Chunk::addParticle()
{
  particle p;
  p.pos = origin;
  p.bf  = 0;
  Activate(p.bf);
  pchunk.push_front(p);
}


// @brief pchunkの終端に粒子を追加
// @note ライフタイムはparticleが持っている値を継承
void Chunk::addParticle(particle p)
{
  Activate(p.bf);
  pchunk.push_back(p);
}
