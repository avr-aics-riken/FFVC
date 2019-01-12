#ifndef _PT_TRACKING_H_
#define _PT_TRACKING_H_

//##################################################################################
//
// Copyright (c) 2019 Research Institute for Information Technology, Kyushu university
// All rights researved.
//
//##################################################################################

// @file   Tracking.h
// @brief  Tracking class Header
// @author riit
// @note

#include "FB_Define.h"
#include <stdio.h>
#include "common/Vec3.h" // defined in Polylib
#include "PtDefine.h"

using namespace Vec3class;

/**
 *  粒子追跡クラス
 */
class Tracking {

protected:
  Vec3r pch;            ///< セル幅
  Vec3r org;            ///< ローカル起点座標
  Vec3r reg;            ///< サブドメイン領域サイズ
  int size[3];          ///< サブドメインサイズ
  int gc;               ///< ガイドセル数
  int myRank;

  REAL_TYPE* vSrc;      ///< 速度データ
  int* bcd;             ///< BCindex B


public:
  /// コンストラクタ
  Tracking()
  {
    size[0] = size[1] = size[2] = 0;
    gc = 0;
    vSrc = NULL;
    bcd = NULL;
    myRank = 0;
  }

  Tracking(const int m_sz[3],
           const int m_gc,
           const REAL_TYPE m_org[3],
           const REAL_TYPE m_reg[3],
           const REAL_TYPE m_pch[3],
           REAL_TYPE* m_Vsrc,
           int* m_bcd,
           const int m_rank)
  {
    size[0] = m_sz[0];
    size[1] = m_sz[1];
    size[2] = m_sz[2];
    gc = m_gc;

    org.x = m_org[0];
    org.y = m_org[1];
    org.z = m_org[2];

    reg.x = m_reg[0];
    reg.y = m_reg[1];
    reg.z = m_reg[2];

    pch.x = m_pch[0];
    pch.y = m_pch[1];
    pch.z = m_pch[2];

    myRank = m_rank;

    //printf("[%d] sz  (%14d %14d %14d) \n", myRank, size[0], size[1], size[2]);
    //printf("[%d] org (%14.6e %14.6e %14.6e) \n", myRank, org.x, org.y, org.z);
    //printf("[%d] reg (%14.6e %14.6e %14.6e) \n", myRank, reg.x, reg.y, reg.z);
    //printf("[%d] pch (%14.6e %14.6e %14.6e) \n", myRank, pch.x, pch.y, pch.z);

    vSrc = m_Vsrc;
    bcd = m_bcd;
  }


  /// デストラクタ
  ~Tracking() {}

public:


  // @brief Euler陽解法による経路積分、チェックあり
  // @param [in,out] p   粒子座標
  // @param [in]     dt  時間積分幅
  // @retval 移動先のランク番号 自領域の場合-1
  Vec3i integrate_Euler(Vec3r& p, const REAL_TYPE dt);


  // @brief RK2による経路積分
  // @param [in,out] p   粒子座標
  // @param [in]     dt  時間積分幅
  // @retval 移動先のランク番号 自領域の場合-1
  Vec3i integrate_RK2(Vec3r& p, const REAL_TYPE dt);


protected:

  // @brief ランクマップを作成 findRank()で利用
  //bool makeRankMap();


  // @brief pの属するランクを探す 全領域外の判断も
  // @param [in] p  空間座標
  // @retval 移動する方向
  Vec3i findRankDir(Vec3r p);


  // @brief pを含むセルセンター起点の小さい方のインデクス
  inline Vec3i getBase(const Vec3r p) {
    Vec3i b;
    b.x = (int)(p.x + 0.5);
    b.y = (int)(p.y + 0.5);
    b.z = (int)(p.z + 0.5);
    return b;
  }


  // @brief ゼロ起点の計算空間座標
  inline Vec3r getRidx(const Vec3r p) {
    return (p - org) / pch;
  }


  // @brief 内挿係数
  inline Vec3r getCoef(const Vec3r p, const Vec3i base) {
    Vec3r c;
    c.x = p.x + 0.5 - (REAL_TYPE)base.x;
    c.y = p.y + 0.5 - (REAL_TYPE)base.y;
    c.z = p.z + 0.5 - (REAL_TYPE)base.z;
    return c;
  }


  // @brief pを含むセルのインデクス（F)
  inline Vec3i getInCellF(const Vec3r p) {
    Vec3i c;
    c.x = (int)((p.x - org.x) / pch.x) + 1;
    c.y = (int)((p.y - org.y) / pch.y) + 1;
    c.z = (int)((p.z - org.z) / pch.z) + 1;
    return c;
  }


  Vec3i shift1(Vec3i index) { return Vec3i(index.x+1, index.y  , index.z  ); }
  Vec3i shift2(Vec3i index) { return Vec3i(index.x  , index.y+1, index.z  ); }
  Vec3i shift3(Vec3i index) { return Vec3i(index.x+1, index.y+1, index.z  ); }
  Vec3i shift4(Vec3i index) { return Vec3i(index.x  , index.y  , index.z+1); }
  Vec3i shift5(Vec3i index) { return Vec3i(index.x+1, index.y  , index.z+1); }
  Vec3i shift6(Vec3i index) { return Vec3i(index.x  , index.y+1, index.z+1); }
  Vec3i shift7(Vec3i index) { return Vec3i(index.x+1, index.y+1, index.z+1); }

  Vec3i shift_xm(Vec3i index) { return Vec3i(index.x-1, index.y  , index.z  ); }
  Vec3i shift_xp(Vec3i index) { return Vec3i(index.x+1, index.y  , index.z  ); }
  Vec3i shift_ym(Vec3i index) { return Vec3i(index.x  , index.y-1, index.z  ); }
  Vec3i shift_yp(Vec3i index) { return Vec3i(index.x  , index.y+1, index.z  ); }
  Vec3i shift_zm(Vec3i index) { return Vec3i(index.x  , index.y  , index.z-1); }
  Vec3i shift_zp(Vec3i index) { return Vec3i(index.x  , index.y  , index.z+1); }


  /**
   @brief Linear interpolation
   @param [in] t target dividing point [0,1] normalized
   @param [in] v values
   */
  template <typename T>
  inline T Linear(REAL_TYPE t, T v[2]) const
  {
    return ( v[0]*(1.0-t) + v[1]*t );
  }

  /**
   @brief Bilinear interpolation
   @param [in] t target dividing point [0,1][0,1] normalized
   @param [in] v values
   @note
   - v[0]; c[0,0]
   - v[1]; c[1,0]
   - v[2]; c[0,1]
   - v[3]; c[1,1]
   */
  template <typename T>
  inline T Bilinear(REAL_TYPE t[2], T v[4]) const
  {
    T r[2];
    r[0] = Linear( t[0], &v[0] );
    r[1] = Linear( t[0], &v[2] );
    return ( Linear( t[1], r ) );
  }

  /**
   @brief Trilinear interpolation
   @param [in] t target dividing point [0,1][0,1][0,1] normalized
   @param [in] v values
   @note
   - v[0]; c[0,0,0]
   - v[1]; c[1,0,0]
   - v[2]; c[0,1,0]
   - v[3]; c[1,1,0]
   - v[4]; c[0,0,1]
   - v[5]; c[1,0,1]
   - v[6]; c[0,1,1]
   - v[7]; c[1,1,1]
   */
  template <typename T>
  inline T Trilinear(Vec3r t, T v[8]) const
  {
    T r[2];
    r[0] = Bilinear( &t.x, &v[0] );
    r[1] = Bilinear( &t.x, &v[4] );
    return ( Linear( t.z, r ) );
  }


  /// セルが流体セルかどうか調べる
  ///
  ///    @param[in] index セルインデックス
  ///
  bool isFluid(Vec3i index) {
    size_t m =_F_IDX_S3D(index.x, index.y, index.z, size[0], size[1], size[2], gc);
    return IS_FLUID(bcd[m]);
  }


  /// セルでのスカラー値を取得
  ///
  ///    @param [in] s スカラー配列
  ///    @param [in] index セルインデックス
  ///    @return セルでのスカラー値
  ///
  REAL_TYPE getScalar(const REAL_TYPE* s, Vec3i index) {
    size_t m =_F_IDX_S3D(index.x, index.y, index.z, size[0], size[1], size[2], gc);
    return s[m];
  }


  /// セルでのベクトル値を取得
  ///
  ///    @param [in] idx セルインデックス
  ///    @return セルでのベクトル値
  ///
  Vec3r getVector(const Vec3i idx)
  {
    Vec3r vRet;
    size_t mx = _F_IDX_V3D(idx.x, idx.y, idx.z, 0, size[0], size[1], size[2], gc);
    size_t my = _F_IDX_V3D(idx.x, idx.y, idx.z, 1, size[0], size[1], size[2], gc);
    size_t mz = _F_IDX_V3D(idx.x, idx.y, idx.z, 2, size[0], size[1], size[2], gc);
    //printf("[%d] %d %d %d : %d %d %d : %d\n", myRank, idx.x, idx.y, idx.z,size[0], size[1], size[2], gc);
    vRet.x = vSrc[mx];
    vRet.y = vSrc[my];
    vRet.z = vSrc[mz];
    return vRet;
  }


  // @brief 固体セルを含む場合trueを返す
  bool isBCell(const Vec3i base)
  {
    int i = base.x;
    int j = base.y;
    int k = base.z;
    size_t m0 = _F_IDX_S3D(i  , j  , k  , size[0], size[1], size[2], gc);
    size_t m1 = _F_IDX_S3D(i+1, j  , k  , size[0], size[1], size[2], gc);
    size_t m2 = _F_IDX_S3D(i  , j+1, k  , size[0], size[1], size[2], gc);
    size_t m3 = _F_IDX_S3D(i+1, j+1, k  , size[0], size[1], size[2], gc);
    size_t m4 = _F_IDX_S3D(i  , j  , k+1, size[0], size[1], size[2], gc);
    size_t m5 = _F_IDX_S3D(i+1, j  , k+1, size[0], size[1], size[2], gc);
    size_t m6 = _F_IDX_S3D(i  , j+1, k+1, size[0], size[1], size[2], gc);
    size_t m7 = _F_IDX_S3D(i+1, j+1, k+1, size[0], size[1], size[2], gc);

    bool onBoundary=false;   ///< 8セルに固体があるかどうかのフラグ
    if ( !IS_FLUID(bcd[m0]) ||
         !IS_FLUID(bcd[m1]) ||
         !IS_FLUID(bcd[m2]) ||
         !IS_FLUID(bcd[m3]) ||
         !IS_FLUID(bcd[m4]) ||
         !IS_FLUID(bcd[m5]) ||
         !IS_FLUID(bcd[m6]) ||
         !IS_FLUID(bcd[m7]) ) onBoundary = true; // ある
    //printf("[%d] BC=%d : %d %d %d\n", myRank, onBoundary, i,j,k);
    return onBoundary;
  }


  Vec3r samplingV(Vec3r coef, const Vec3i base)
  {
    if ( isBCell(base) ) return getVector(base);  // nearest

    Vec3r r[8];
    r[0] = getVector(base);
    r[1] = getVector(shift1(base));
    r[2] = getVector(shift2(base));
    r[3] = getVector(shift3(base));
    r[4] = getVector(shift4(base));
    r[5] = getVector(shift5(base));
    r[6] = getVector(shift6(base));
    r[7] = getVector(shift7(base));

    return Trilinear(coef, r);
  }

  /// @brief 自領域内に存在するかどうかを判断
  /// @param [in] x  coordinate
  /// @retval true 自領域内
  bool inOwnRegion(const Vec3r p)
  {
    if ( p.x < org.x || p.x > org.x + reg.x
      || p.y < org.y || p.y > org.y + reg.y
      || p.z < org.z || p.z > org.z + reg.z ) return false;
    return true;
  }


  // @brief Euler陽解法による経路積分
  // @param [in,out] p   粒子座標
  // @param [in]     dt  時間積分幅
  Vec3r Euler0(Vec3r p, const REAL_TYPE dt);


  // @brief Euler陽解法の一部
  // @param [in,out] p   粒子座標
  // @param [in]     dt  時間積分幅
  Vec3r Euler1(Vec3r p, const REAL_TYPE dt);


  // @brief Runge-Kutta 2 stage
  // @param [in,out] p   粒子座標
  // @param [in]     dt  時間積分幅
  Vec3r RK2(Vec3r p, const REAL_TYPE dt);


};
#endif // _PT_TRACKING_H_
