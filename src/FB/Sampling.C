// #################################################################
//
// CAERU Library
//
// Copyright (c) All right reserved. 2012
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
//
// #################################################################

//@file   Sampling.C
//@brief  FlowBase Sampling class
//@author kero

#include "Sampling.h"

/// 渦度を計算
///
///   @param[in] v 速度配列
///   @param[in] index 計算位置のセルインデックス
///   @return 渦度ベクトル
///
Vec3r Sampling::calcVorticity(const REAL_TYPE* v, Vec3i index)
{
  if (isFluid(index)) {
  
    Vec3r v0 = getVector(v, index);
    Vec3r v1 = 2.0 * v0 - v00;

    Vec3r v_xm = isFluid(shift_xm(index)) ? getVector(v, shift_xm(index)) : v1;
    Vec3r v_xp = isFluid(shift_xp(index)) ? getVector(v, shift_xp(index)) : v1;
    Vec3r v_ym = isFluid(shift_ym(index)) ? getVector(v, shift_ym(index)) : v1;
    Vec3r v_yp = isFluid(shift_yp(index)) ? getVector(v, shift_yp(index)) : v1;
    Vec3r v_zm = isFluid(shift_zm(index)) ? getVector(v, shift_zm(index)) : v1;
    Vec3r v_zp = isFluid(shift_zp(index)) ? getVector(v, shift_zp(index)) : v1;

    REAL_TYPE hx = 0.5 / pch.x;
    REAL_TYPE hy = 0.5 / pch.y;
    REAL_TYPE hz = 0.5 / pch.z;

    Vec3r omg;
    omg.x = (v_yp.z - v_ym.z)*hy - (v_zp.y - v_zm.y)*hz;
    omg.y = (v_zp.x - v_zm.x)*hz - (v_xp.z - v_xm.z)*hx;
    omg.z = (v_xp.y - v_xm.y)*hx - (v_yp.x - v_ym.x)*hy;

    return omg;
  }
  else {
    return 0.0;
  }
}


/* -------- Nearest -------------------------------------------------------- */

/// コンストラクタ.
///
///   @param[in] mode サンプリングモード
///   @param[in] size,guide ローカルセル数，ガイドセル数
///   @param[in] crd  モニタ点座標
///   @param[in] org,pch  ローカル領域基点座標，セル幅
///   @param[in] v00  座標系移動速度
///   @param[in] bcd  BCindex ID
///
Nearest::Nearest(int mode, int size[], int guide,
                 Vec3r crd, Vec3r org, Vec3r pch, Vec3r v00, int* bcd) 
  : Sampling(mode, size, guide, crd, org, pch, v00, bcd)
{
}


/// スカラー変数をサンプリング.
///
///   @param[in] s サンプリング元スカラー変数配列
///
REAL_TYPE Nearest::samplingScalar(const REAL_TYPE* s)
{
  return getScalar(s, cIndex);
}


/// 速度をサンプリング.
///
///   @param[in] v サンプリング元速度配列
///
Vec3r Nearest::samplingVelocity(const REAL_TYPE* v)
{
  return getVector(v, cIndex);
}


/// 全圧をサンプリング.
///
///   @param[in] v サンプリング元速度配列
///   @param[in] p サンプリング元圧力配列
///
REAL_TYPE Nearest::samplingTotalPressure(const REAL_TYPE* v, const REAL_TYPE* p)
{
  return calcTotalPressure(getVector(v, cIndex), getScalar(p, cIndex));
}


/// 渦度をサンプリング.
///
///   @param[in] v サンプリング元速度配列
///
Vec3r Nearest::samplingVorticity(const REAL_TYPE* v)
{
  return calcVorticity(v, cIndex);
}


/* -------- Smoothing ------------------------------------------------------ */

/// コンストラクタ.
///
///   隣接セルに対して局所平均計算に使用するかどうかのフラグを立てる.
///
///   @param[in] mode サンプリングモード
///   @param[in] size,guide ローカルセル数，ガイドセル数
///   @param[in] crd  モニタ点座標
///   @param[in] org,pch  ローカル領域基点座標，セル幅
///   @param[in] v00  座標系移動速度
///   @param[in] bcd  BCindex ID
///
Smoothing::Smoothing(int mode, int size[], int guide,
                     Vec3r crd, Vec3r org, Vec3r pch, Vec3r v00, int* bcd) 
  : Sampling(mode, size, guide, crd, org, pch, v00, bcd)
{
  add_xm = permitToAdd(shift_xm(cIndex)) ? true : false;
  add_xp = permitToAdd(shift_xp(cIndex)) ? true : false;
  add_ym = permitToAdd(shift_ym(cIndex)) ? true : false;
  add_yp = permitToAdd(shift_yp(cIndex)) ? true : false;
  add_zm = permitToAdd(shift_zm(cIndex)) ? true : false;
  add_zp = permitToAdd(shift_zp(cIndex)) ? true : false;

  nAdd = 1;
  if (add_xm) nAdd++;
  if (add_xp) nAdd++;
  if (add_ym) nAdd++;
  if (add_yp) nAdd++;
  if (add_zm) nAdd++;
  if (add_zp) nAdd++;
}


/// スカラー変数をサンプリング.
///
///   @param[in] s サンプリング元スカラー変数配列
///
REAL_TYPE Smoothing::samplingScalar(const REAL_TYPE* s)
{
  REAL_TYPE sRet = getScalar(s, cIndex);
  if (add_xm) sRet += getScalar(s, shift_xm(cIndex));
  if (add_xp) sRet += getScalar(s, shift_xp(cIndex));
  if (add_ym) sRet += getScalar(s, shift_ym(cIndex));
  if (add_yp) sRet += getScalar(s, shift_yp(cIndex));
  if (add_zm) sRet += getScalar(s, shift_zm(cIndex));
  if (add_zp) sRet += getScalar(s, shift_zp(cIndex));
  sRet /= nAdd;

  return sRet;
}


/// 速度をサンプリング.
///
///   @param[in] v サンプリング元速度配列
///
Vec3r Smoothing::samplingVelocity(const REAL_TYPE* v)
{
  Vec3r vRet = getVector(v, cIndex);
  if (add_xm) vRet += getVector(v, shift_xm(cIndex));
  if (add_xm) vRet += getVector(v, shift_xp(cIndex));
  if (add_xm) vRet += getVector(v, shift_ym(cIndex));
  if (add_xm) vRet += getVector(v, shift_yp(cIndex));
  if (add_xm) vRet += getVector(v, shift_zm(cIndex));
  if (add_xm) vRet += getVector(v, shift_zp(cIndex));
  vRet /= nAdd;

  return vRet;
}


/// 全圧をサンプリング.
///
///   @param[in] v サンプリング元速度配列
///   @param[in] p サンプリング元圧力配列
///
REAL_TYPE Smoothing::samplingTotalPressure(const REAL_TYPE* v, const REAL_TYPE* p)
{
  REAL_TYPE tp = calcTotalPressure(getVector(v, cIndex), getScalar(p, cIndex));
  if (add_xm) {
    tp += calcTotalPressure(getVector(v, shift_xm(cIndex)), getScalar(p, shift_xm(cIndex)));
  }
  if (add_xp) {
    tp += calcTotalPressure(getVector(v, shift_xp(cIndex)), getScalar(p, shift_xp(cIndex)));
  }
  if (add_ym) {
    tp += calcTotalPressure(getVector(v, shift_ym(cIndex)), getScalar(p, shift_ym(cIndex)));
  }
  if (add_yp) {
    tp += calcTotalPressure(getVector(v, shift_yp(cIndex)), getScalar(p, shift_yp(cIndex)));
  }
  if (add_zm) {
    tp += calcTotalPressure(getVector(v, shift_zm(cIndex)), getScalar(p, shift_zm(cIndex)));
  }
  if (add_zp) {
    tp += calcTotalPressure(getVector(v, shift_zp(cIndex)), getScalar(p, shift_zp(cIndex)));
  }
  tp /= nAdd;

  return tp;
}


/// 渦度をサンプリング.
///
///   @param[in] v サンプリング元速度配列
///
Vec3r Smoothing::samplingVorticity(const REAL_TYPE* v)
{
  Vec3r omg = calcVorticity(v, cIndex);
  if (add_xm) omg += calcVorticity(v, shift_xm(cIndex));
  if (add_xp) omg += calcVorticity(v, shift_xp(cIndex));
  if (add_ym) omg += calcVorticity(v, shift_ym(cIndex));
  if (add_yp) omg += calcVorticity(v, shift_yp(cIndex));
  if (add_zm) omg += calcVorticity(v, shift_zm(cIndex));
  if (add_zp) omg += calcVorticity(v, shift_zp(cIndex));
  omg /= nAdd;

  return omg;
}

/* -------- Interpolation -------------------------------------------------- */

/// コンストラクタ.
///
///   補間計算の基準セルのインデックスおよび補間係数を計算.
///   流体-固体境界に接しているかをチェック.
///
///   @param[in] mode サンプリングモード
///   @param[in] size,guide ローカルセル数，ガイドセル数
///   @param[in] crd  モニタ点座標
///   @param[in] org,pch  ローカル領域基点座標，セル幅
///   @param[in] v00  座標系移動速度
///   @param[in] bcd  BCindex ID
///
Interpolation::Interpolation(int mode, int size[], int guide,
                             Vec3r crd, Vec3r org, Vec3r pch, Vec3r v00, int* bcd) 
  : Sampling(mode, size, guide, crd, org, pch, v00, bcd)
{
  Vec3r c = (crd - org) / pch;
  base.x = (int)(c.x + 0.5);
  base.y = (int)(c.y + 0.5);
  base.z = (int)(c.z + 0.5);
  coef[0] = c.x + 0.5 - (REAL_TYPE)base[0];
  coef[1] = c.y + 0.5 - (REAL_TYPE)base[1];
  coef[2] = c.z + 0.5 - (REAL_TYPE)base[2];

  onBoundary = false;
  if (checkBoundary(base))         onBoundary = true;
  if (checkBoundary(shift1(base))) onBoundary = true;
  if (checkBoundary(shift2(base))) onBoundary = true;
  if (checkBoundary(shift3(base))) onBoundary = true;
  if (checkBoundary(shift4(base))) onBoundary = true;
  if (checkBoundary(shift5(base))) onBoundary = true;
  if (checkBoundary(shift6(base))) onBoundary = true;
  if (checkBoundary(shift7(base))) onBoundary = true;
}


/// スカラー変数をサンプリング.
///
///   @param[in] s サンプリング元スカラー変数配列
///
REAL_TYPE Interpolation::samplingScalar(const REAL_TYPE* s)
{
  if (onBoundary) return getScalar(s, cIndex);  // nearest

  REAL_TYPE r[8];
  r[0] = getScalar(s, base);          // (0, 0, 0)
  r[1] = getScalar(s, shift1(base));  // (1, 0, 0)
  r[2] = getScalar(s, shift2(base));  // (0, 1, 0)
  r[3] = getScalar(s, shift3(base));  // (1, 1, 0)
  r[4] = getScalar(s, shift4(base));  // (0, 0, 1)
  r[5] = getScalar(s, shift5(base));  // (1, 0, 1)
  r[6] = getScalar(s, shift6(base));  // (0, 1, 1)
  r[7] = getScalar(s, shift7(base));  // (1, 1, 1)
  return Trilinear(coef, r);
}


/// 速度をサンプリング.
///
///   @param[in] v サンプリング元速度配列
///
Vec3r Interpolation::samplingVelocity(const REAL_TYPE* v)
{
  if (onBoundary) return getVector(v, cIndex);  // nearest

  Vec3r r[8];
  r[0] = getVector(v, base);
  r[1] = getVector(v, shift1(base));
  r[2] = getVector(v, shift2(base));
  r[3] = getVector(v, shift3(base));
  r[4] = getVector(v, shift4(base));
  r[5] = getVector(v, shift5(base));
  r[6] = getVector(v, shift6(base));
  r[7] = getVector(v, shift7(base));
  return Trilinear(coef, r);
}


/// 全圧をサンプリング.
///
///   @param[in] v サンプリング元速度配列
///   @param[in] p サンプリング元圧力配列
///
REAL_TYPE Interpolation::samplingTotalPressure(const REAL_TYPE* v, const REAL_TYPE* p)
{
  if (onBoundary) return calcTotalPressure(getVector(v, cIndex), getScalar(p, cIndex));  // nearest

  REAL_TYPE r[8];
  r[0] = calcTotalPressure(getVector(v, base), getScalar(p, base));
  r[1] = calcTotalPressure(getVector(v, shift1(base)), getScalar(p, shift1(base)));
  r[2] = calcTotalPressure(getVector(v, shift2(base)), getScalar(p, shift2(base)));
  r[3] = calcTotalPressure(getVector(v, shift3(base)), getScalar(p, shift3(base)));
  r[4] = calcTotalPressure(getVector(v, shift4(base)), getScalar(p, shift4(base)));
  r[5] = calcTotalPressure(getVector(v, shift5(base)), getScalar(p, shift5(base)));
  r[6] = calcTotalPressure(getVector(v, shift6(base)), getScalar(p, shift6(base)));
  r[7] = calcTotalPressure(getVector(v, shift7(base)), getScalar(p, shift7(base)));
  return Trilinear(coef, r);
}


/// 渦度をサンプリング.
///
///   @param[in] v サンプリング元速度配列
///
Vec3r Interpolation::samplingVorticity(const REAL_TYPE* v)
{
  if (onBoundary) return calcVorticity(v, cIndex);  // nearest

  Vec3r r[8];
  r[0] = calcVorticity(v, base);
  r[1] = calcVorticity(v, shift1(base));
  r[2] = calcVorticity(v, shift2(base));
  r[3] = calcVorticity(v, shift3(base));
  r[4] = calcVorticity(v, shift4(base));
  r[5] = calcVorticity(v, shift5(base));
  r[6] = calcVorticity(v, shift6(base));
  r[7] = calcVorticity(v, shift7(base));
  return Trilinear(coef, r);
}


/* -------- InterpolationStgV ---------------------------------------------- */

/// コンストラクタ.
///
///   @param[in] mode サンプリングモード
///   @param[in] size,guide ローカルセル数，ガイドセル数
///   @param[in] crd  モニタ点座標
///   @param[in] org,pch  ローカル領域基点座標，セル幅
///   @param[in] v00  座標系移動速度
///   @param[in] bcd  BCindex ID
///
InterpolationStgV::InterpolationStgV(int mode, int size[], int guide,
                                     Vec3r crd, Vec3r org, Vec3r pch, Vec3r v00, int* bcd) 
  : Interpolation(mode, size, guide, crd, org, pch, v00, bcd)
{
  Vec3r c = (crd - org) / pch;
  base_s[0] = (int)(c.x);
  base_s[1] = (int)(c.y);
  base_s[2] = (int)(c.z);
  coef_s[0] = c.x - (REAL_TYPE)base_s[0];
  coef_s[1] = c.y - (REAL_TYPE)base_s[1];
  coef_s[2] = c.z - (REAL_TYPE)base_s[2];
}


/// 速度をサンプリング.
///
///   @param[in] v サンプリング元速度配列
///
Vec3r InterpolationStgV::samplingVelocity(const REAL_TYPE* v)
{
  Vec3r vRet;
  REAL_TYPE t[3];
  REAL_TYPE r[8];
  int i, j, k;

  t[0] = coef_s[0];
  t[1] = coef[1];
  t[2] = coef[2];
  i = base_s.x;
  j = base.y;
  k = base.z;
  r[0] = v[ FBUtility::getFindexV3DEx(size, guide, 0, i  , j  , k  ) ];
  r[1] = v[ FBUtility::getFindexV3DEx(size, guide, 0, i+1, j  , k  ) ];
  r[2] = v[ FBUtility::getFindexV3DEx(size, guide, 0, i  , j+1, k  ) ];
  r[3] = v[ FBUtility::getFindexV3DEx(size, guide, 0, i+1, j+1, k  ) ];
  r[4] = v[ FBUtility::getFindexV3DEx(size, guide, 0, i  , j  , k+1) ];
  r[5] = v[ FBUtility::getFindexV3DEx(size, guide, 0, i+1, j  , k+1) ];
  r[6] = v[ FBUtility::getFindexV3DEx(size, guide, 0, i  , j+1, k+1) ];
  r[7] = v[ FBUtility::getFindexV3DEx(size, guide, 0, i+1, j+1, k+1) ];
  vRet.x = Trilinear(t, r);

  t[0] = coef[0];
  t[1] = coef_s[1];
  t[2] = coef[2];
  i = base.x;
  j = base_s.y;
  k = base.z;
  r[0] = v[ FBUtility::getFindexV3DEx(size, guide, 1, i  , j  , k  ) ];
  r[1] = v[ FBUtility::getFindexV3DEx(size, guide, 1, i+1, j  , k  ) ];
  r[2] = v[ FBUtility::getFindexV3DEx(size, guide, 1, i  , j+1, k  ) ];
  r[3] = v[ FBUtility::getFindexV3DEx(size, guide, 1, i+1, j+1, k  ) ];
  r[4] = v[ FBUtility::getFindexV3DEx(size, guide, 1, i  , j  , k+1) ];
  r[5] = v[ FBUtility::getFindexV3DEx(size, guide, 1, i+1, j  , k+1) ];
  r[6] = v[ FBUtility::getFindexV3DEx(size, guide, 1, i  , j+1, k+1) ];
  r[7] = v[ FBUtility::getFindexV3DEx(size, guide, 1, i+1, j+1, k+1) ];
  vRet.y = Trilinear(t, r);

  t[0] = coef[0];
  t[1] = coef[1];
  t[2] = coef_s[2];
  i = base.x;
  j = base.y;
  k = base_s.z;
  r[0] = v[ FBUtility::getFindexV3DEx(size, guide, 2, i  , j  , k  ) ];
  r[1] = v[ FBUtility::getFindexV3DEx(size, guide, 2, i+1, j  , k  ) ];
  r[2] = v[ FBUtility::getFindexV3DEx(size, guide, 2, i  , j+1, k  ) ];
  r[3] = v[ FBUtility::getFindexV3DEx(size, guide, 2, i+1, j+1, k  ) ];
  r[4] = v[ FBUtility::getFindexV3DEx(size, guide, 2, i  , j  , k+1) ];
  r[5] = v[ FBUtility::getFindexV3DEx(size, guide, 2, i+1, j  , k+1) ];
  r[6] = v[ FBUtility::getFindexV3DEx(size, guide, 2, i  , j+1, k+1) ];
  r[7] = v[ FBUtility::getFindexV3DEx(size, guide, 2, i+1, j+1, k+1) ];
  vRet.z = Trilinear(t, r);

  return vRet;
}
