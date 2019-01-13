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
 * @file   Tracking.C
 * @brief  Tracking class
 * @author riit
 */

#include "Tracking.h"


//#############################################################################
// @brief Euler陽解法による経路積分
// @param [in,out] p   粒子座標
// @retval 正-計算領域内の他ランク, 負-領域外(-1)、自領域(-2)
Vec3i Tracking::integrate_Euler(Vec3r& p, const REAL_TYPE dt)
{
  //printf("[%d] p:(%14.6f %14.6f %14.6f)\n", myRank, p.x, p.y, p.z);
  // 予測値
  Vec3r q = Euler0(p, dt);

  // 格子幅より大きな移動の場合
  if (distance(q, p) > pch.length()) {
    // 積分幅を半分にして、2回やり直し

    q = Euler0(Euler0(p, 0.5*dt), 0.5*dt);

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
// @retval 正-計算領域内の他ランク, 負-領域外(-1)、自領域(-2)
Vec3i Tracking::integrate_RK2(Vec3r& p, const REAL_TYPE dt)
{
  // 1st step
  Vec3r q = RK2(p, dt);

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
// @brief Euler積分
Vec3r Tracking::Euler0(Vec3r p, const REAL_TYPE dt)
{
  Vec3r c = getRidx(p);
  //printf("[%d] c = (%14.6f %14.6f %14.6f) \n", myRank, c.x, c.y, c.z);
  Vec3i base = getBase(c);
  //printf("[%d] base = (%d %d %d) \n", myRank, base.x, base.y, base.z);
  Vec3r coef = getCoef(c, base);
  //printf("[%d] coef= (%14.6f %14.6f %14.6f) \n", myRank, coef.x, coef.y, coef.z);
  return p + dt * samplingV(coef, base);
}


//#############################################################################
// @brief Euler積分の一部分
Vec3r Tracking::Euler1(Vec3r p, const REAL_TYPE dt)
{
  Vec3r c = getRidx(p);
  Vec3i base = getBase(c) ;
  Vec3r coef = getCoef(c, base);
  return dt * samplingV(coef, base);
}


//#############################################################################
// @brief RK2
Vec3r Tracking::RK2(Vec3r p, const REAL_TYPE dt)
{
  Vec3r q = Euler0(p, 0.5*dt);
  q = p + Euler1(q, dt);
  return q;
}


//#############################################################################
// @brief pの属するランクを探す 全領域外の判断も
// @param [in] p  空間座標（ローカル）
// @retval 正-計算領域内の他ランク, 負-領域外(-1)、自領域(-2)
Vec3i Tracking::findRankDir(Vec3r p)
{
  Vec3r o = org;
  Vec3r e = reg;
  //printf("o %e %e %e / e %e %e %e / p %e %e %e\n",o.x, o.y, o.z, e.x, e.y, e.z, p.x, p.y, p.z);

  Vec3i d;
  if      (p.x < o.x)       d.x = -1;
  else if (o.x + e.x < p.x) d.x =  0;
  else                      d.x =  1;

  if      (p.y < o.y)       d.y = -1;
  else if (o.y + e.y < p.y) d.y =  0;
  else                      d.y =  1;

  if      (p.z < o.z)       d.z = -1;
  else if (o.z + e.z < p.z) d.z =  0;
  else                      d.z =  1;

  return d;
}
