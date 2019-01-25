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
// @param [in,out] p        粒子座標
// @param [out]    v        粒子並進速度
// @param [out]    wallFlag 壁を通り抜けた場合true
// @param [in]     dt       積分幅
// @retval 正-計算領域内の他ランク, 負-領域外(-1)、自領域(-2)
// @note 流体計算のΔtの設定から格子幅より大きなトラベルはない
Vec3i Tracking::integrate_Euler(Vec3r& p, Vec3r& v, bool& wallFlag, const REAL_TYPE dt)
{
  // 予測値
	bool oflag=false;
  v = getV(p, oflag);

  Vec3r q = p + dt * v;

  /* 格子幅より大きな移動の場合
  if (distance(q, p) > pch.length()) {
    // 積分幅を半分にして、2回やり直し
    Vec3r r = p + 0.5f*dt * v;
    q = r + 0.5f*dt * getV(r);

    Vec3i b = getBase(p);
    fprintf(stdout, "Retry at (%d, %d, %d)\n", b.x, b.y, b.z);
  }
	 */
	
	// 壁近傍の積分を行った場合
	if (oflag)
	{
		wallFlag = chkWallPassing(p, q);
	}
	
  p = q;
  Vec3i dst(0);
  if ( !inOwnRegion(q) ) dst = findRankDir(q);

  return dst;
}


//#############################################################################
// @brief RK2による経路積分
// @param [in,out] p        粒子座標
// @param [out]    v        粒子並進速度
// @param [out]    wallFlag 壁を通り抜けた場合true
// @param [in]     dt       積分幅
// @retval 正-計算領域内の他ランク, 負-領域外(-1)、自領域(-2)
Vec3i Tracking::integrate_RK2(Vec3r& p, Vec3r& v, bool& wallFlag, const REAL_TYPE dt)
{
  // 1st step
	bool oflag=false;
	v = getV(p, oflag);
  Vec3r q = p + 0.5f*dt * v;
  q = p + dt * getV(q, oflag);
	
	// 壁近傍の積分を行った場合
	if (oflag)
	{
		wallFlag = chkWallPassing(p, q);
	}
	
  p = q;
  Vec3i dst(0);
  if ( !inOwnRegion(q) ) dst = findRankDir(q);
  return dst;
}


//#############################################################################
// @brief RK2
Vec3r Tracking::RK2(Vec3r p, const REAL_TYPE dt)
{
	bool oflag=false;
  Vec3r q = p + 0.5f*dt * getV(p, oflag);
  q = p + dt * getV(q, oflag);
  return q;
}


//#############################################################################
// @brief pの属する隣接ランクの方向を返す
// @param [in] p  空間座標（ローカル）
// @retval 隣接方向{-1, 0, 1}
// @note inOwnRegion()と同じ判断基準のこと
Vec3i Tracking::findRankDir(Vec3r p)
{
  Vec3r o = org;
  Vec3r e = reg;

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


//#############################################################################
// @brief p点の速度をサンプリングする　境界セルでは0.5がけ
// @param [in]  p   粒子座標
// @param [out] of  判定フラグ
// @todo 速度は壁からの距離にするとよいが
Vec3r Tracking::getV(const Vec3r p, bool& of)
{
	Vec3r c = getRidx(p);
	Vec3i base = getBase(c) ;
	Vec3r coef = getCoef(c, base);
	
	REAL_TYPE xsi = coef.x;
	REAL_TYPE eta = coef.y;
	REAL_TYPE zta = coef.z;
	
	int i = base.x;
	int j = base.y;
	int k = base.z;
	
	int ix = size[0];
	int jx = size[1];
	int kx = size[2];
	
	// 点pを含むセルインデクス(fortran)
	Vec3i q = getInCellF(p);
	int bd = bcd[ _F_IDX_S3D(q.x, q.y , q.z , ix, jx, kx, gc) ];
	
	Vec3r vp = {0.0, 0.0, 0.0};
	
	
	if ( IS_FLUID(bd) )
	{
		int bx1 = bid[ _F_IDX_S3D(i  , j  , k  , ix, jx, kx, gc) ];
		int bx2 = bid[ _F_IDX_S3D(i+1, j  , k  , ix, jx, kx, gc) ];
		int bx3 = bid[ _F_IDX_S3D(i  , j+1, k  , ix, jx, kx, gc) ];
		int bx4 = bid[ _F_IDX_S3D(i+1, j+1, k  , ix, jx, kx, gc) ];
		int bx5 = bid[ _F_IDX_S3D(i  , j  , k+1, ix, jx, kx, gc) ];
		int bx6 = bid[ _F_IDX_S3D(i+1, j  , k+1, ix, jx, kx, gc) ];
		int bx7 = bid[ _F_IDX_S3D(i  , j+1, k+1, ix, jx, kx, gc) ];
		
		// いずれかの辺に交点がある (id != 0)
		if ( ( 0 != getBit5(bx1, X_plus) )
			|| ( 0 != getBit5(bx3, X_plus) )
			|| ( 0 != getBit5(bx5, X_plus) )
			|| ( 0 != getBit5(bx7, X_plus) )
			|| ( 0 != getBit5(bx1, Y_plus) )
			|| ( 0 != getBit5(bx2, Y_plus) )
			|| ( 0 != getBit5(bx5, Y_plus) )
			|| ( 0 != getBit5(bx6, Y_plus) )
			|| ( 0 != getBit5(bx1, Z_plus) )
			|| ( 0 != getBit5(bx2, Z_plus) )
			|| ( 0 != getBit5(bx3, Z_plus) )
			|| ( 0 != getBit5(bx4, Z_plus) ) )
		{
			// pを含むセルの速度に0.5がけ
			vp.x = 0.5 * vSrc[ _F_IDX_V3D(q.x, q.y , q.z , 0, ix, jx, kx, gc) ];
			vp.y = 0.5 * vSrc[ _F_IDX_V3D(q.x, q.y , q.z , 1, ix, jx, kx, gc) ];
			vp.z = 0.5 * vSrc[ _F_IDX_V3D(q.x, q.y , q.z , 2, ix, jx, kx, gc) ];
			of = true;
		}
		else // 交点のない全セル流体の場合、通常のtri-linear
		{
			REAL_TYPE u1 = vSrc[ _F_IDX_V3D(i  , j  , k  , 0, ix, jx, kx, gc) ];
			REAL_TYPE u2 = vSrc[ _F_IDX_V3D(i+1, j  , k  , 0, ix, jx, kx, gc) ];
			REAL_TYPE u3 = vSrc[ _F_IDX_V3D(i  , j+1, k  , 0, ix, jx, kx, gc) ];
			REAL_TYPE u4 = vSrc[ _F_IDX_V3D(i+1, j+1, k  , 0, ix, jx, kx, gc) ];
			REAL_TYPE u5 = vSrc[ _F_IDX_V3D(i  , j  , k+1, 0, ix, jx, kx, gc) ];
			REAL_TYPE u6 = vSrc[ _F_IDX_V3D(i+1, j  , k+1, 0, ix, jx, kx, gc) ];
			REAL_TYPE u7 = vSrc[ _F_IDX_V3D(i  , j+1, k+1, 0, ix, jx, kx, gc) ];
			REAL_TYPE u8 = vSrc[ _F_IDX_V3D(i+1, j+1, k+1, 0, ix, jx, kx, gc) ];
			
			REAL_TYPE v1 = vSrc[ _F_IDX_V3D(i  , j  , k  , 1, ix, jx, kx, gc) ];
			REAL_TYPE v2 = vSrc[ _F_IDX_V3D(i+1, j  , k  , 1, ix, jx, kx, gc) ];
			REAL_TYPE v3 = vSrc[ _F_IDX_V3D(i  , j+1, k  , 1, ix, jx, kx, gc) ];
			REAL_TYPE v4 = vSrc[ _F_IDX_V3D(i+1, j+1, k  , 1, ix, jx, kx, gc) ];
			REAL_TYPE v5 = vSrc[ _F_IDX_V3D(i  , j  , k+1, 1, ix, jx, kx, gc) ];
			REAL_TYPE v6 = vSrc[ _F_IDX_V3D(i+1, j  , k+1, 1, ix, jx, kx, gc) ];
			REAL_TYPE v7 = vSrc[ _F_IDX_V3D(i  , j+1, k+1, 1, ix, jx, kx, gc) ];
			REAL_TYPE v8 = vSrc[ _F_IDX_V3D(i+1, j+1, k+1, 1, ix, jx, kx, gc) ];
			
			REAL_TYPE w1 = vSrc[ _F_IDX_V3D(i  , j  , k  , 2, ix, jx, kx, gc) ];
			REAL_TYPE w2 = vSrc[ _F_IDX_V3D(i+1, j  , k  , 2, ix, jx, kx, gc) ];
			REAL_TYPE w3 = vSrc[ _F_IDX_V3D(i  , j+1, k  , 2, ix, jx, kx, gc) ];
			REAL_TYPE w4 = vSrc[ _F_IDX_V3D(i+1, j+1, k  , 2, ix, jx, kx, gc) ];
			REAL_TYPE w5 = vSrc[ _F_IDX_V3D(i  , j  , k+1, 2, ix, jx, kx, gc) ];
			REAL_TYPE w6 = vSrc[ _F_IDX_V3D(i+1, j  , k+1, 2, ix, jx, kx, gc) ];
			REAL_TYPE w7 = vSrc[ _F_IDX_V3D(i  , j+1, k+1, 2, ix, jx, kx, gc) ];
			REAL_TYPE w8 = vSrc[ _F_IDX_V3D(i+1, j+1, k+1, 2, ix, jx, kx, gc) ];
			
			vp.x = (1.0-xsi)*(1.0-eta)*(1.0-zta)* u1
					 +      xsi *(1.0-eta)*(1.0-zta)* u2
					 + (1.0-xsi)*     eta *(1.0-zta)* u3
					 +      xsi *     eta *(1.0-zta)* u4
					 + (1.0-xsi)*(1.0-eta)*     zta * u5
					 +      xsi *(1.0-eta)*     zta * u6
					 + (1.0-xsi)*     eta *     zta * u7
					 +      xsi *     eta *     zta * u8;
			
			vp.y = (1.0-xsi)*(1.0-eta)*(1.0-zta)* v1
					 +      xsi *(1.0-eta)*(1.0-zta)* v2
					 + (1.0-xsi)*     eta *(1.0-zta)* v3
					 +      xsi *     eta *(1.0-zta)* v4
					 + (1.0-xsi)*(1.0-eta)*     zta * v5
					 +      xsi *(1.0-eta)*     zta * v6
					 + (1.0-xsi)*     eta *     zta * v7
					 +      xsi *     eta *     zta * v8;
			
			vp.z = (1.0-xsi)*(1.0-eta)*(1.0-zta)* w1
					 +      xsi *(1.0-eta)*(1.0-zta)* w2
					 + (1.0-xsi)*     eta *(1.0-zta)* w3
					 +      xsi *     eta *(1.0-zta)* w4
					 + (1.0-xsi)*(1.0-eta)*     zta * w5
					 +      xsi *(1.0-eta)*     zta * w6
					 + (1.0-xsi)*     eta *     zta * w7
					 +      xsi *     eta *     zta * w8;
		}
	}
	// 固体の場合にはゼロを返す
	else
	{
		vp.assign(0.0, 0.0, 0.0);
	}
	
	return vp;
}


//#############################################################################
// @brief 壁を通過したかを簡易判定
// @param [in]  p   元の粒子座標
// @param [in]  q   移動後の粒子座標
// @retval true - 壁を通り抜けた
// @todo 本当は距離で判定する
bool Tracking::chkWallPassing(const Vec3r p, const Vec3r q)
{
	Vec3i ip = getInCellF(p);
	Vec3i iq = getInCellF(q);
	Vec3i d = iq - ip;
	int b = bid[ _F_IDX_S3D(ip.x, ip.y , ip.z , size[0], size[1], size[2], gc) ];
	
	// 各方向に移動し、かつ交点がある場合
	if ( d.x > 0 && ( 0 != getBit5(b, X_plus)) )  return true;
	if ( d.x < 0 && ( 0 != getBit5(b, X_minus)) ) return true;
	if ( d.y > 0 && ( 0 != getBit5(b, Y_plus)) )  return true;
	if ( d.y < 0 && ( 0 != getBit5(b, Y_minus)) ) return true;
	if ( d.z > 0 && ( 0 != getBit5(b, Z_plus)) )  return true;
	if ( d.z < 0 && ( 0 != getBit5(b, Z_minus)) ) return true;
	
	return false;
}
