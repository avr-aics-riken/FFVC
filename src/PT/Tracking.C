//##################################################################################
//
// Copyright (c) 2019 Research Institute for Information Technology, Kyushu university
// All rights researved.
//
//##################################################################################


/**
 * @file   Tracking.C
 * @brief  Tracking class
 * @author riit
 */

#include "Tracking.h"


//#############################################################################
// @brief Euler陽解法による経路積分
// @param [in,out] p   粒子座標
// @param [out]    v   粒子並進速度
// @retval 正-計算領域内の他ランク, 負-領域外(-1)、自領域(-2)
Vec3i Tracking::integrate_Euler(Vec3r& p, Vec3r& v, const REAL_TYPE dt)
{
  //printf("[%d] p:(%14.6f %14.6f %14.6f)\n", myRank, p.x, p.y, p.z);
  // 予測値
  v = getV(p);
  //printf("[%d] p:(%14.6f %14.6f %14.6f)\n", myRank, v.x, v.y, v.z);
  Vec3r q = p + dt * v;

  // 格子幅より大きな移動の場合
  if (distance(q, p) > pch.length()) {
    // 積分幅を半分にして、2回やり直し
    Vec3r r = p + 0.5f*dt * v;
    q = r + 0.5f*dt * getV(r);

    Vec3i b = getBase(p);
    fprintf(stdout, "Retry at (%d, %d, %d)\n", b.x, b.y, b.z);
  }
  p = q;
  Vec3i dst(-2);
  if ( !inOwnRegion(q) ) dst = findRankDir(q);

  return dst;
}


//#############################################################################
// @brief RK2による経路積分
// @param [in,out] p   粒子座標
// @param [out]    v   粒子並進速度
// @retval 正-計算領域内の他ランク, 負-領域外(-1)、自領域(-2)
Vec3i Tracking::integrate_RK2(Vec3r& p, Vec3r& v, const REAL_TYPE dt)
{
  // 1st step
  v = getV(p);
  Vec3r q = p + 0.5f*dt * v;
  q = p + dt * getV(q);

  // 格子幅より大きな移動の場合
  if (distance(q, p) > pch.length()) {
    // 積分幅を半分にして、2回やり直し
    q = RK2(RK2(p, 0.5*dt), 0.5*dt);
    Vec3i b = getBase(p) ;
    fprintf(stdout, "Retry at (%d, %d, %d)\n", b.x, b.y, b.z);
  }
  p = q;
  Vec3i dst(-2);
  if ( !inOwnRegion(q) ) dst = findRankDir(q);
  return dst;
}


//#############################################################################
// @brief RK2
Vec3r Tracking::RK2(Vec3r p, const REAL_TYPE dt)
{
  Vec3r q = p + 0.5f*dt * getV(p);
  q = p + dt * getV(q);
  return q;
}


//#############################################################################
// @brief pの属するランクを探す 全領域外の判断も
// @param [in] p  空間座標（ローカル）
// @retval 正-計算領域内の他ランク, 負-領域外(-1)、自領域(-2)
// @note inOwnRegion()と同じ判断基準のこと
Vec3i Tracking::findRankDir(Vec3r p)
{
  Vec3r o = org;
  Vec3r e = reg;
  //printf("o %e %e %e / e %e %e %e / p %e %e %e\n",o.x, o.y, o.z, e.x, e.y, e.z, p.x, p.y, p.z);

  Vec3i d;
  if      (p.x < o.x)        d.x = -1;
  else if (o.x + e.x <= p.x) d.x =  1;
  else                       d.x =  0;

  if      (p.y < o.y)        d.y = -1;
  else if (o.y + e.y <= p.y) d.y =  1;
  else                       d.y =  0;

  if      (p.z < o.z)        d.z = -1;
  else if (o.z + e.z <= p.z) d.z =  1;
  else                       d.z =  0;

  return d;
}
