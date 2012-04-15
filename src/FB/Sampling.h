#ifndef _SKL_FB_SAMPLING_H_
#define _SKL_FB_SAMPLING_H_

/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file MSampling.h
//@brief FlowBase Sampling class Header
//@author keno, FSI Team, VCAD, RIKEN

#include "FBDefine.h"
#include "FBUtility.h"
#include "vec3.h"
#include "basic_func.h"

using namespace FB;

/**
 * Samplingクラス
 *
 *   サンプリング基底クラス(アブストラクト)
 */
class Sampling {
public:

  /// モニタ点状態型
  enum PointStatus {
    POINT_STATUS_OK,   ///< 正常
    UNEXPECTED_SOLID,  ///< 流体セルを指定したが固体だった
    UNEXPECTED_FLUID,  ///< 固体セルを指定したが流体だった
  };

protected:
  int mode;          ///< サンプリングモード

  Vec3i cIndex;      ///< モニタ点を含むセルのインデックス

  unsigned size[3];  ///< セル(ローカル)サイズ
  unsigned guide;    ///< ガイドセル数
  Vec3r pch;         ///< セル幅

  Vec3r v00;         ///< 座標系移動速度
  unsigned* bcd;     ///< BCindex ID
  
public:
  /// ディフォルトコンストラクタ.
  Sampling() {}

  /// コンストラクタ.
  ///
  ///   @param[in] mode サンプリングモード
  ///   @param[in] size,guide ローカルセル数，ガイドセル数
  ///   @param[in] crd  モニタ点座標
  ///   @param[in] org,pch  ローカル領域基点座標，セル幅
  ///   @param[in] v00  座標系移動速度
  ///   @param[in] bcd  BCindex ID
  ///
  Sampling(int mode, unsigned size[], unsigned guide,
           Vec3r crd, Vec3r org, Vec3r pch, Vec3r v00, unsigned* bcd) {
    this->mode = mode;
    this->size[0] = size[0];
    this->size[1] = size[1];
    this->size[2] = size[2];
    this->guide = guide;
    this->pch = pch;
    this->v00 = v00;
    this->bcd = bcd;

    Vec3r c = (crd - org) / pch;
    cIndex.x = (int)(c.x) + 1;
    cIndex.y = (int)(c.y) + 1;
    cIndex.z = (int)(c.z) + 1;
  } 

  /// デストラクタ.
  virtual ~Sampling() {}

  /// モニタ点の状態を返す.
  PointStatus checkMonitorPoint() {
    if      (mode == SAMPLING_FLUID_ONLY && !isFluid(cIndex)) return UNEXPECTED_SOLID;
    else if (mode == SAMPLING_SOLID_ONLY &&  isFluid(cIndex)) return UNEXPECTED_FLUID;
    else return POINT_STATUS_OK;   // mode == ALL
  }

  /// 速度をサンプリング.
  ///
  ///   @param[in] v サンプリング元速度配列
  ///
  virtual Vec3r samplingVelocity(const REAL_TYPE* v) = 0;

  /// 圧力をサンプリング.
  ///
  ///   @param[in] p サンプリング元圧力配列
  ///
  virtual REAL_TYPE samplingPressure(const REAL_TYPE* p) = 0;

  /// 温度をサンプリング.
  ///
  ///   @param[in] t サンプリング元温度配列
  ///
  virtual REAL_TYPE samplingTemperature(const REAL_TYPE* t) = 0;

  /// 全圧をサンプリング.
  ///
  ///   @param[in] v サンプリング元速度配列
  ///   @param[in] p サンプリング元圧力配列
  ///
  virtual REAL_TYPE samplingTotalPressure(const REAL_TYPE* v, const REAL_TYPE* p) = 0;

  /// 渦度をサンプリング.
  ///
  ///   @param[in] v サンプリング元速度配列
  ///
  virtual Vec3r samplingVorticity(const REAL_TYPE* v) = 0;

protected:
  /// 全圧を計算.
  ///
  ///   @param[in] v  速度
  ///   @param[in] p  圧力
  ///   @return 全圧
  REAL_TYPE calcTotalPressure(const Vec3r v, REAL_TYPE p) {
    Vec3r v1 = v - v00;
    return 0.5 * v1.lengthSquared() + p;
  }

  /// 渦度を計算
  ///
  ///   @param[in] v 速度配列
  ///   @param[in] index セルインデックス
  ///   @return 渦度ベクトル
  ///
  Vec3r calcVorticity(const REAL_TYPE* v, Vec3i index);


  /// セルインデックスを(1,0,0)シフト.
  Vec3i shift1(Vec3i index) { return Vec3i(index.x+1, index.y  , index.z  ); }

  /// セルインデックスを(0,1,0)シフト.
  Vec3i shift2(Vec3i index) { return Vec3i(index.x  , index.y+1, index.z  ); }

  /// セルインデックスを(1,1,0)シフト.
  Vec3i shift3(Vec3i index) { return Vec3i(index.x+1, index.y+1, index.z  ); }

  /// セルインデックスを(0,0,1)シフト.
  Vec3i shift4(Vec3i index) { return Vec3i(index.x  , index.y  , index.z+1); }

  /// セルインデックスを(1,0,1)シフト.
  Vec3i shift5(Vec3i index) { return Vec3i(index.x+1, index.y  , index.z+1); }

  /// セルインデックスを(0,1,1)シフト.
  Vec3i shift6(Vec3i index) { return Vec3i(index.x  , index.y+1, index.z+1); }

  /// セルインデックスを(1,1,1)シフト.
  Vec3i shift7(Vec3i index) { return Vec3i(index.x+1, index.y+1, index.z+1); }

  /// セルインデックスを-xシフト.
  Vec3i shift_xm(Vec3i index) { return Vec3i(index.x-1, index.y  , index.z  ); }

  /// セルインデックスを+xシフト.
  Vec3i shift_xp(Vec3i index) { return Vec3i(index.x+1, index.y  , index.z  ); }

  /// セルインデックスを-yシフト.
  Vec3i shift_ym(Vec3i index) { return Vec3i(index.x  , index.y-1, index.z  ); }

  /// セルインデックスを+yシフト.
  Vec3i shift_yp(Vec3i index) { return Vec3i(index.x  , index.y+1, index.z  ); }

  /// セルインデックスを-zシフト.
  Vec3i shift_zm(Vec3i index) { return Vec3i(index.x  , index.y  , index.z-1); }

  /// セルインデックスを+zシフト.
  Vec3i shift_zp(Vec3i index) { return Vec3i(index.x  , index.y  , index.z+1); }

  /// セルが流体セルかどうか調べる.
  ///
  ///    @param[in] index セルインデックス
  ///
  bool isFluid(Vec3i index) {
    return IS_FLUID(bcd[FBUtility::getFindexS3D(size, guide, index.x, index.y, index.z)]);
  }

  /// セルでのスカラー値を取得.
  ///
  ///    @param[in] s スカラー配列
  ///    @param[in] index セルインデックス
  ///    @return セルでのスカラー値
  ///
  REAL_TYPE getScalar(const REAL_TYPE* s, Vec3i index) {
    return s[FBUtility::getFindexS3D(size, guide, index.x, index.y, index.z)];
  }

  /// セルでのベクトル値を取得.
  ///
  ///    @param[in] v ベクトル配列
  ///    @param[in] index セルインデックス
  ///    @return セルでのベクトル値
  ///
  Vec3r getVector(const REAL_TYPE* v, Vec3i index) {
    Vec3r vRet;
    vRet.x = v[FBUtility::getFindexV3DEx(size, guide, index.x, index.y, index.z, 0)];
    vRet.y = v[FBUtility::getFindexV3DEx(size, guide, index.x, index.y, index.z, 1)];
    vRet.z = v[FBUtility::getFindexV3DEx(size, guide, index.x, index.y, index.z, 2)];
    return vRet;
  }

};

/**
  * Nearestクラス.
  *
  *    - モニタ点を含むセルの値をサンプリングする．
  */
class Nearest : public Sampling {
protected:

  /// スカラー変数をサンプリング.
  ///
  ///   @param[in] s サンプリング元スカラー変数配列
  ///
  REAL_TYPE samplingScalar(const REAL_TYPE* s);

public:

  /// コンストラクタ.
  ///
  ///   @param[in] mode サンプリングモード
  ///   @param[in] size,guide ローカルセル数，ガイドセル数
  ///   @param[in] crd  モニタ点座標
  ///   @param[in] org,pch  ローカル領域基点座標，セル幅
  ///   @param[in] v00  座標系移動速度
  ///   @param[in] bcd  BCindex ID
  ///
  Nearest(int mode, unsigned size[], unsigned guide,
          Vec3r crd, Vec3r org, Vec3r pch, Vec3r v00, unsigned* bcd);

  /// デストラクタ.
  ~Nearest() {}

  /// 速度をサンプリング.
  ///
  ///   @param[in] v サンプリング元速度配列
  ///
  Vec3r samplingVelocity(const REAL_TYPE* v);

  /// 圧力をサンプリング.
  ///
  ///   @param[in] p サンプリング元圧力配列
  ///
  REAL_TYPE samplingPressure(const REAL_TYPE* p) { return samplingScalar(p); }

  /// 温度をサンプリング.
  ///
  ///   @param[in] t サンプリング元温度配列
  ///
  REAL_TYPE samplingTemperature(const REAL_TYPE* t) { return samplingScalar(t); }

  /// 全圧をサンプリング.
  ///
  ///   @param[in] v サンプリング元速度配列
  ///   @param[in] p サンプリング元圧力配列
  ///
  REAL_TYPE samplingTotalPressure(const REAL_TYPE*v, const REAL_TYPE* p);

  /// 渦度をサンプリング.
  ///
  ///   @param[in] v サンプリング元速度配列
  ///
  Vec3r samplingVorticity(const REAL_TYPE* v);
};


/**
 * Smoothingクラス
 *
 *   - モニタ点を含むセルとその隣接セルを用いて平均値を計算する.
 *   - モードがallの場合は、6個の全隣接セルを用いる.
 *   - モードがfluidの場合は、隣接セルとしては流体セルのみを用いる.
 *   - モードがsolidの場合は、隣接セルとしては固体セルのみを用いる.
 */
class Smoothing : public Sampling {

protected:

  bool add_xm;  ///< -xセルを平均値計算に参加させるかどうかのフラグ
  bool add_xp;  ///< +xセルを平均値計算に参加させるかどうかのフラグ
  bool add_ym;  ///< -yセルを平均値計算に参加させるかどうかのフラグ
  bool add_yp;  ///< +yセルを平均値計算に参加させるかどうかのフラグ
  bool add_zm;  ///< -zセルを平均値計算に参加させるかどうかのフラグ
  bool add_zp;  ///< +zセルを平均値計算に参加させるかどうかのフラグ

  int nAdd;  ///< 平均値計算参加セル数

  /// そのセルを平均値計算に使用するかどうか判定する.
  ///
  ///    @param[in] index セルインデックス
  ///
  bool permitToAdd(Vec3i index) {
    if      (mode == SAMPLING_FLUID_ONLY) return isFluid(index);
    else if (mode == SAMPLING_SOLID_ONLY) return !isFluid(index);
    else return true;  // モードがallの場合はセルのfluid/solidは考慮しない
  }

  /// スカラー変数をサンプリング.
  ///
  ///   @param[in] s サンプリング元スカラー変数配列
  ///
  REAL_TYPE samplingScalar(const REAL_TYPE* s);

public:

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
  Smoothing(int mode, unsigned size[], unsigned guide,
            Vec3r crd, Vec3r org, Vec3r pch, Vec3r v00, unsigned* bcd) ;

  /// デストラクタ.
  ~Smoothing() {}

  /// 速度をサンプリング.
  ///
  ///   @param[in] v サンプリング元速度配列
  ///
  Vec3r samplingVelocity(const REAL_TYPE* v);

  /// 圧力をサンプリング.
  ///
  ///   @param[in] p サンプリング元圧力配列
  ///
  REAL_TYPE samplingPressure(const REAL_TYPE* p) { return samplingScalar(p); }

  /// 温度をサンプリング.
  ///
  ///   @param[in] t サンプリング元温度配列
  ///
  REAL_TYPE samplingTemperature(const REAL_TYPE* t) { return samplingScalar(t); }

  /// 全圧をサンプリング.
  ///
  ///   @param[in] v サンプリング元速度配列
  ///   @param[in] p サンプリング元圧力配列
  ///
  REAL_TYPE samplingTotalPressure(const REAL_TYPE*v, const REAL_TYPE* s);

  /// 渦度をサンプリング.
  ///
  ///   @param[in] v サンプリング元速度配列
  ///
  Vec3r samplingVorticity(const REAL_TYPE* v);

};


/**
 * Interpolationクラス
 *
 *   - モニタ点をはさむ8セルから三重線形補間により計算する．
 *   - モードがfluidで、モニタ点を含むセル以外の7セルのどれかが固体セルの
 *     場合は，Nearest(モニタ点を含むセルの値)によりサンプリング.
 *   - モードがsolidで、モニタ点を含むセル以外の7セルのどれかが流体セルの
 *     場合は，Nearest(モニタ点を含むセルの値)によりサンプリング.
 */
class Interpolation : public Sampling {
protected:
  Vec3i base;       ///< 線形補間基準セルのインデックス
  REAL_TYPE coef[3]; ///< 線形補間係数

  bool onBoundary;  ///< 8セルにモードと違う態(流体/固体)があるかどうかのフラグ

  /// そのセルがモードと違うかどうか調べる.
  bool checkBoundary(Vec3i index) {
    if      (mode == SAMPLING_FLUID_ONLY) return !isFluid(index);
    else if (mode == SAMPLING_SOLID_ONLY) return isFluid(index);
    else return false;  // モードがallの場合はセルのfluid/solidは考慮しない
  }

  /// スカラー変数をサンプリング.
  ///
  ///   @param[in] s サンプリング元スカラー変数配列
  ///
  REAL_TYPE samplingScalar(const REAL_TYPE* s);

public:

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
  Interpolation(int mode, unsigned size[], unsigned guide,
                Vec3r crd, Vec3r org, Vec3r pch, Vec3r v00, unsigned* bcd);

  /// デストラクタ.
  ~Interpolation() {}

  /// 速度をサンプリング.
  ///
  ///   @param[in] v サンプリング元速度配列
  ///
  Vec3r samplingVelocity(const REAL_TYPE* v);

  /// 圧力をサンプリング.
  ///
  ///   @param[in] p サンプリング元圧力配列
  ///
  REAL_TYPE samplingPressure(const REAL_TYPE* p) { return samplingScalar(p); }

  /// 温度をサンプリング.
  ///
  ///   @param[in] t サンプリング元温度配列
  ///
  REAL_TYPE samplingTemperature(const REAL_TYPE* t) { return samplingScalar(t); }

  /// 全圧をサンプリング.
  ///
  ///   @param[in] v サンプリング元速度配列
  ///   @param[in] p サンプリング元圧力配列
  ///
  REAL_TYPE samplingTotalPressure(const REAL_TYPE* v, const REAL_TYPE* p);

  /// 渦度をサンプリング.
  ///
  ///   @param[in] v サンプリング元速度配列
  ///
  Vec3r samplingVorticity(const REAL_TYPE* v);
};


/**
 * InterpolationStgVクラス
 *
 *   - スタガード速度配置用の三重線形補間
 *   - Interpolationクラスを継承
 *   - 旧版との比較用(全圧と渦度には未対応)
 */
class InterpolationStgV : public Interpolation {
protected:
  Vec3i base_s;        ///< 線形補間基準セルのインデックス(スタガード配置)
  REAL_TYPE coef_s[3];  ///< 線形補間係数(スタガード配置)

public:

  /// コンストラクタ.
  ///
  ///   @param[in] mode サンプリングモード
  ///   @param[in] size,guide ローカルセル数，ガイドセル数
  ///   @param[in] crd  モニタ点座標
  ///   @param[in] org,pch  ローカル領域基点座標，セル幅
  ///   @param[in] v00  座標系移動速度
  ///   @param[in] bcd  BCindex ID
  ///
  InterpolationStgV(int mode, unsigned size[], unsigned guide,
                    Vec3r crd, Vec3r org, Vec3r pch, Vec3r v00, unsigned* bcd);

  /// デストラクタ.
  ~InterpolationStgV() {}

  /// 速度をサンプリング.
  ///
  ///   @param[in] v サンプリング元速度配列
  ///
  Vec3r samplingVelocity(const REAL_TYPE* v);

};

#endif // _SKL_FB_SAMPLING_H_
