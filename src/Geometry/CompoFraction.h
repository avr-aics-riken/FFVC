#ifndef _FB_FRACTION_H_
#define _FB_FRACTION_H_

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
//##################################################################################

/**
 * @file   CompoFraction.h
 * @brief  Component Fraction class Header
 * @author aics
 */

#ifndef DISABLE_MPI
  #include "mpi.h" // add header explicitly to avoid compile error for Intel MPI
#else
  #include "cpm_mpistub.h"
#endif

#include "FB_Define.h"
#include "common/Vec3.h" // defined in Polylib
#include <math.h>
#include "Geometry.h"

using namespace Vec3class;

class CompoFraction {

public:
  /** 形状 */
  enum ShapeList
  {
    shape_plate=1,
    shape_cylinder
  };

protected:
  // ローカル計算領域情報
  int size[3];     ///< セル(ローカル)サイズ
  int guide;       ///< ガイドセル数
  int myrank;      ///< rank
  int division;    ///< 細分化の分割数
  Vec3r pch;       ///< セル幅
  Vec3r org;       ///< 計算領域の基点
  Vec3r angle;     ///< 変換の回転角度

  // 形状パラメータ
  int smode;         ///< 形状モード
  int shape;         ///< 形状（1-plate, 2-cylinder）
  REAL_TYPE depth;   ///< 厚さ
  REAL_TYPE width;   ///< 矩形の幅（dir方向)
  REAL_TYPE height;  ///< 矩形の高さ
  REAL_TYPE R1;      ///< ファン半径
  REAL_TYPE R2;      ///< ボス半径
  Vec3r nv;          ///< 法線方向ベクトル（流出方向）
  Vec3r center;      ///< 形状の中心座標（前面の中心位置）
  Vec3r dir;         ///< 矩形の方向規定の参照ベクトル
  Vec3r box_min;     ///< Bounding boxの最小値
  Vec3r box_max;     ///< Bounding boxの最大値

public:
  /** デフォルトコンストラクタ */
  CompoFraction() {
    smode = 0;
    shape = 0;
    depth = 0.0;
    width = 0.0;
    height = 0.0;
    R1 = 0.0;
    R2 = 0.0;
  }

  /** コンストラクタ
   * @param [in] size   ローカルセル数
   * @param [in] guide  ガイドセル数
   * @param [in] m_rank rank
   * @param [in] pch    ローカル領域セル幅
   * @param [in] org    ローカル領域基点座標
   * @param [in] div    サブディビジョンの分割数
   */
  CompoFraction(const int* size, const int guide, const int m_rank, const REAL_TYPE* pch, const REAL_TYPE* org, const int div)
  {
    this->size[0]  = size[0];
    this->size[1]  = size[1];
    this->size[2]  = size[2];
    this->guide    = guide;
    this->myrank   = m_rank;
    this->pch.x    = pch[0];
    this->pch.y    = pch[1];
    this->pch.z    = pch[2];
    this->org.x    = org[0];
    this->org.y    = org[1];
    this->org.z    = org[2];
    this->division = div;
  }

  /** デストラクタ */
  virtual ~CompoFraction() {}

protected:
  /**
   * @brief 矩形領域のbboxを計算
   * @param [in] mn bboxの最小位置
   * @param [in] mx bboxの最大位置
   * @erturn 投影面積
   */
  REAL_TYPE calcBboxRect(Vec3r& mn, Vec3r& mx);


  /**
   * @brief 円筒領域のbboxを計算
   * @param [in] mn bboxの最小位置
   * @param [in] mx bboxの最大位置
   * @erturn 投影面積
   */
  REAL_TYPE calcBboxCircle(Vec3r& mn, Vec3r& mx);


  /**
   * @brief 円柱底面と線分の交点を求める
   * @param [in] st         bbox開始インデクス
   * @param [in] ed         bbox終了インデクス
   * @param [in,out] bid    境界ID（5ビット幅x6方向）
   * @param [in,out] cut    カット情報
   * @param [in]     pl     テストする平面方程式の係数
   * @param [in]     s_id   エンコードするID
   * @param [in]     Dsize  サイズ
   */
  int CylinderPlane(const int st[],
                    const int ed[],
                    int* bid,
                    long long* cut,
                    const REAL_TYPE pl[4],
                    const int s_id,
                    const int* Dsize);


  /**
   * @brief 円柱側面と線分の交点を求める
   * @param [in] st         bbox開始インデクス
   * @param [in] ed         bbox終了インデクス
   * @param [in,out] bid    境界ID（5ビット幅x6方向）
   * @param [in,out] cut    カット情報
   * @param [in]     s_id   エンコードするID
   * @param [in]     Dsize  サイズ
   */
  int CylinderSide(const int st[],
                   const int ed[],
                   int* bid,
                   long long* cut,
                   const int s_id,
                   const int* Dsize);


  /**
   * @brief 点pの属するセルインデクスを求める
   * @param [out] w インデクス
   * @param [in]  p テストする座標値
   * @note Fortran index
   */
  void findIndex(int* w, const Vec3r p)
  {
    Vec3r q = (p-org)/pch;

    w[0] = (int)ceil(q.x);
    w[1] = (int)ceil(q.y);
    w[2] = (int)ceil(q.z);

    if ( w[0] < 1 ) w[0] = 1;
    if ( w[1] < 1 ) w[1] = 1;
    if ( w[2] < 1 ) w[2] = 1;

    if ( w[0] > size[0] ) w[0] = size[0];
    if ( w[1] > size[1] ) w[1] = size[1];
    if ( w[2] > size[2] ) w[2] = size[2];
  }


  /**
   * @brief 円と線分の交点を求める
   * @param [out]  X      交点Xの座標
   * @param [in]   A      始点座標ベクトル
   * @param [in]   B      終点座標ベクトル
   * @retval  交点距離の小さい方を返す。負値の場合は交点が無い
   */
  REAL_TYPE intersectLineByCylinder(Vec3r& X, const Vec3r A, const Vec3r B);


  /**
   * @brief 平面と線分の交点を求める
   * @param [out]  X   平面P上の交点Xの座標
   * @param [in]   A   線分ABの端点
   * @param [in]   B   線分ABの端点
   * @param [in]   PL  平面PLの係数 ax+by+cz+d=0
   * @retval 平面P上の交点XのAからの距離, 負値の場合は交点が無い
   * @see http://www.sousakuba.com/Programming/gs_plane_line_intersect.html
   */
  REAL_TYPE intersectLineByPlane(Vec3r& X, const Vec3r A, const Vec3r B, const REAL_TYPE PL[4]);


  /**
   * @brief 円筒形状の内外判定
   * @param [in] p    テスト点座標
   * @param [in] mode 変換モード
   * @return 内部のときに1を返す
   * @note 186 flop
   */
  inline int judgeCylider(const Vec3r p, bool mode=false)
  {
    Vec3r q;

    if ( !mode )
    {
      q = rotate(angle, p-center); // 181 flop
    }
    else
    {
      q = p;
    }

    // palteの場合には、半径方向の判断のみでよい
    if ( shape == shape_cylinder )
    {
      if ( (q.z < 0.0) || (q.z > depth)  ) return 0;
    }

    REAL_TYPE r = sqrtf(q.x*q.x + q.y*q.y);

    return ( r<=R1 && r>=R2 ) ? 1 : 0;
  }


  /**
   * @brief 矩形形状の内外判定
   * @param [in] p    テスト点座標
   * @param [in] mode 変換モード
   * @return 内部のときに1を返す
   * @note 183 flop
   */
  inline int judgeRect(const Vec3r p, bool mode=false)
  {
    Vec3r q;

    if ( !mode )
    {
      q = rotate(angle, p-center); // 181 flop
    }
    else
    {
      q = p;
    }

    if ( (q.z < 0.0) || (q.z > depth)  ) return 0;
    if ( fabs(q.x) > 0.5*width )  return 0;
    if ( fabs(q.y) > 0.5*height ) return 0;

    return 1;
  }


  /**
   * @brief 回転ベクトルp(alpha, beta, gamma)でベクトルuを回転する
   * @param [in] p 回転角度
   * @param [in] u 方向ベクトル
   * @return 角度
   * @note sin, cos = 5flop, dot=5flop, total 181flop
   */
  inline Vec3r rotate(const Vec3r p, const Vec3r u)
  {
    Vec3r a, b, c;

    // line vector expression
    a.x =  cos(p.y)*cos(p.z);
    a.y =  sin(p.x)*sin(p.y)*cos(p.z) - cos(p.x)*sin(p.z);
    a.z =  cos(p.x)*sin(p.y)*cos(p.z) + sin(p.x)*sin(p.z);

    b.x =  cos(p.y)*sin(p.z);
    b.y =  sin(p.x)*sin(p.y)*sin(p.z) + cos(p.x)*cos(p.z);
    b.z =  cos(p.x)*sin(p.y)*sin(p.z) - sin(p.x)*cos(p.z);

    c.x = -sin(p.y);
    c.y =  sin(p.x)*cos(p.y);
    c.z =  cos(p.x)*cos(p.y);

    return Vec3r( dot(a, u), dot(b, u), dot(c, u) );
  }


  /**
   * @brief 回転ベクトルp(alpha, beta, gamma)に対して，-pでベクトルuを回転する
   * @param [in] p 回転角度
   * @param [in] u 方向ベクトル
   * @return 角度
   * @note sin, cos = 5flop, dot=5flop, total 181flop
   */
  inline Vec3r rotate_inv(const Vec3r p, const Vec3r u)
  {
    Vec3r a, b, c;

    // line vector expression
    a.x =  cos(p.y)*cos(p.z);
    a.y =  cos(p.y)*sin(p.z);
    a.z = -sin(p.y);

    b.x =  sin(p.x)*sin(p.y)*cos(p.z) - cos(p.x)*sin(p.z);
    b.y =  sin(p.x)*sin(p.y)*sin(p.z) + cos(p.x)*cos(p.z);
    b.z =  sin(p.x)*cos(p.y);

    c.x =  cos(p.x)*sin(p.y)*cos(p.z) + sin(p.x)*sin(p.z);
    c.y =  cos(p.x)*sin(p.y)*sin(p.z) - sin(p.x)*cos(p.z);
    c.z =  cos(p.x)*cos(p.y);

    return Vec3r( dot(a, u), dot(b, u), dot(c, u) );
  }


  /**
   * @brief インデックスを(1,0,0)シフト
   * @param [in] index 元のインデクス
   * @param [in] h     シフト幅
   */
  inline Vec3r shift_f1(const Vec3r index, const Vec3r h)
  {
    return Vec3r(index.x+h.x, index.y, index.z);
  }


  /**
   * @brief インデックスを(0,1,0)シフト
   * @param [in] index 元のインデクス
   * @param [in] h     シフト幅
   */
  inline Vec3r shift_f2 (const Vec3r index, const Vec3r h)
  {
    return Vec3r(index.x, index.y+h.y, index.z);
  }


  /**
   * @brief インデックスを(1,1,0)シフト
   * @param [in] index 元のインデクス
   * @param [in] h     シフト幅
   */
  inline Vec3r shift_f3 (const Vec3r index, const Vec3r h)
  {
    return Vec3r(index.x+h.x, index.y+h.y, index.z);
  }


  /**
   * @brief インデックスを(0,0,1)シフト
   * @param [in] index 元のインデクス
   * @param [in] h     シフト幅
   */
  inline Vec3r shift_f4 (const Vec3r index, const Vec3r h)
  {
    return Vec3r(index.x, index.y, index.z+h.z);
  }


  /**
   * @brief インデックスを(1,0,1)シフト
   * @param [in] index 元のインデクス
   * @param [in] h     シフト幅
   */
  inline Vec3r shift_f5 (const Vec3r index, const Vec3r h)
  {
    return Vec3r(index.x+h.x, index.y, index.z+h.z);
  }


  /**
   * @brief インデックスを(0,1,1)シフト
   * @param [in] index 元のインデクス
   * @param [in] h     シフト幅
   */
  inline Vec3r shift_f6 (const Vec3r index, const Vec3r h)
  {
    return Vec3r(index.x, index.y+h.y, index.z+h.z);
  }


  /**
   * @brief インデックスを(1,1,1)シフト
   * @param [in] index 元のインデクス
   * @param [in] h     シフト幅
   */
  inline Vec3r shift_f7 (const Vec3r index, const Vec3r h)
  {
    return Vec3r(index.x+h.x, index.y+h.y, index.z+h.z);
  }


  /**
   * @brief インデックスを(1,0,0)シフト
   * @param [in] index 元のインデクス
   * @param [in] h     シフト幅
   */
  inline Vec3r shift_E(const Vec3r index, const Vec3r h)
  {
    return shift_f1(index, h);
  }


  /**
   * @brief インデックスを(-1,0,0)シフト
   * @param [in] index 元のインデクス
   * @param [in] h     シフト幅
   */
  inline Vec3r shift_W(const Vec3r index, const Vec3r h)
  {
    return Vec3r(index.x-h.x, index.y, index.z  );
  }


  /**
   * @brief インデックスを(0,1,0)シフト
   * @param [in] index 元のインデクス
   * @param [in] h     シフト幅
   */
  inline Vec3r shift_N(const Vec3r index, const Vec3r h)
  {
    return shift_f2(index, h);
  }


  /**
   * @brief インデックスを(0,-1,0)シフト
   * @param [in] index 元のインデクス
   * @param [in] h     シフト幅
   */
  inline Vec3r shift_S(const Vec3r index, const Vec3r h)
  {
    return Vec3r(index.x, index.y-h.y, index.z);
  }


  /**
   * @brief インデックスを(0,0,1)シフト
   * @param [in] index 元のインデクス
   * @param [in] h     シフト幅
   */
  inline Vec3r shift_T(const Vec3r index, const Vec3r h)
  {
    return shift_f4(index, h);
  }


  /**
   * @brief インデックスを(0,0,-1)シフト
   * @param [in] index 元のインデクス
   * @param [in] h     シフト幅
   */
  inline Vec3r shift_B(const Vec3r index, const Vec3r h)
  {
    return Vec3r(index.x, index.y, index.z-h.z);
  }


  /**
   * @brief 交点情報をアップデート
   * @param [in]     A   線分ABの端点
   * @param [in]     B   線分ABの端点
   * @param [in,out] cut 量子化交点距離情報
   * @param [in,out] bid 交点ID情報
   * @param [in]     dir テストする方向
   * @param [in]     pid polygon id
   * @retval 新規交点の数
   * @note 短い距離を記録
   */
  inline unsigned updateCutCylinder(const Vec3r A,
                                    const Vec3r B,
                                    long long& cut,
                                    int& bid,
                                    const int dir,
                                    const int pid)
  {
    /* 内外判定：judgeCylinder()の戻り値が 1 のとき内側
     * テストする2点が内側の場合には交点がない、それ以外をテスト
     * 外側から内側に向かう線分に対して、交点距離とIDを設定
     */

    unsigned count = 0;

    // 交点座標
    Vec3r X;

    // 交点計算
    float t = intersectLineByCylinder(X, A, B);

    // 9bit幅の量子化
    int r = quantize9(t);

    bool record = false;

    if ( 0.0 <= t && judgeCylider(X, true) )
    {
      // 交点が記録されていない場合 >> 新規記録
      if ( ensCut(cut, dir) == 0 )
      {
        record = true;
        count = 1;
      }
      else // 交点が既に記録されている場合 >> 短い方を記録
      {
        if ( r < getBit9(cut, dir)) record = true;
      }
    }

    if ( record )
    {
      setBit5(bid, pid, dir);
      setCut9(cut, r, dir);
      //printf("%10.6f %6d dir=%d\n", t,r,dir);
    }

    return count;
  }


  /**
   * @brief 交点情報をアップデート
   * @param [in]     A   線分ABの端点
   * @param [in]     B   線分ABの端点
   * @param [in]     pl  平面PLの係数 ax+by+cz+d=0
   * @param [in,out] cut 量子化交点距離情報
   * @param [in,out] bid 交点ID情報
   * @param [in]     dir テストする方向
   * @param [in]     pid polygon id
   * @retval 新規交点の数
   * @note 短い距離を記録
   */
  unsigned updateCutPlane(const Vec3r A,
                          const Vec3r B,
                          const REAL_TYPE pl[4],
                          long long& cut,
                          int& bid,
                          const int dir,
                          const int pid)
  {
    /* 内外判定：judgeCylinder()の戻り値が 1 のとき内側
     * テストする2点が内側の場合には交点がない、それ以外をテスト
     * 外側から内側に向かう線分に対して、交点距離とIDを設定
     */

    unsigned count = 0;

    // 交点座標
    Vec3r X;

    // 交点計算
    float t = intersectLineByPlane(X, A, B, pl);

    if ( t < 0.0 || 1.0 < t ) return 0;


    // 9bit幅の量子化
    int r = quantize9(t);

    bool record = false;

    if ( 0.0 <= t && t <= 1.0 && judgeCylider(X, true) )
    {
      // 交点が記録されていない場合 >> 新規記録
      if ( ensCut(cut, dir) == 0 )
      {
        record = true;
        count = 1;
      }
      else // 交点が既に記録されている場合 >> 短い方を記録
      {
        if ( r < getBit9(cut, dir)) record = true;
      }
    }

    if ( record )
    {
      setBit5(bid, pid, dir);
      setCut9(cut, r, dir);
      //printf("%10.6f %6d dir=%d %d %d\n", t,r,dir, getBit5(bid, dir), getBit9(cut, dir));
    }

    return count;
  }



public:

  /**
   * @brief 形状のbboxと投影面積を求める
   * @param [out] st 開始インデクス
   * @param [out] ed 終了インデクス
   * @return 投影面積
   */
  REAL_TYPE getBboxArea (int* st, int* ed);


  /**
   * @brief 指定法線nvがz軸の方向ベクトルに向かう回転角を計算する
   */
  void getAngle();


  /** 円柱と線分の交点を求める
   * @param [in] st         bbox開始インデクス
   * @param [in] ed         bbox終了インデクス
   * @param [in,out] bid    境界ID（5ビット幅x6方向）
   * @param [in,out] cut    カット情報
   * @param [in]     tgt_id 固体媒質ID
   * @param [in]     Dsize  サイズ
   */
  bool intersectCylinder(const int st[],
                         const int ed[],
                         int* bid,
                         long long* cut,
                         const int tgt_id,
                         const int* Dsize=NULL);

  /**
   * @brief 矩形の形状パラメータをセットする
   * @param [in] m_nv     法線ベクトル
   * @param [in] m_ctr    中心座標
   * @param [in] m_dir    方向ベクトル
   * @param [in] m_depth  厚み
   * @param [in] m_width  幅
   * @param [in] m_height 高さ
   */
  void setShapeParam(const REAL_TYPE m_nv[3],
                     const REAL_TYPE m_ctr[3],
                     const REAL_TYPE m_dir[3],
                     const REAL_TYPE m_depth,
                     const REAL_TYPE m_width,
                     const REAL_TYPE m_height);

  /**
   * @brief 円筒の形状パラメータをセットする
   * @param [in] m_nv     法線ベクトル
   * @param [in] m_ctr    中心座標
   * @param [in] m_depth  厚み
   * @param [in] m_R1     外径
   * @param [in] m_R2     内径
   * @param [in] m_shape  形状 デフォルトでCYLINDER
   */
  void setShapeParam(const REAL_TYPE m_nv[3],
                     const REAL_TYPE m_ctr[3],
                     const REAL_TYPE m_depth,
                     const REAL_TYPE m_R1,
                     const REAL_TYPE m_R2=0.0f,
                     const int m_shape=shape_cylinder);


  /**
   * @brief 体積率が(0,1)の間のセルに対してサブディビジョンを実施
   * @param [in]     st    開始インデクス
   * @param [in]     ed    終了インデクス
   * @param [in,out] vf    フラクション
   * @param [in,out] flop  浮動小数点演算数
   */
  void subdivision(const int st[], const int ed[], REAL_TYPE* vf, double& flop);


  /**
   * @brief セルの8頂点の内外判定を行い，0, 1, otherに分類
   * @param [in]     st    開始インデクス
   * @param [in]     en    終了インデクス
   * @param [in,out] vf    フラクション
   * @param [in,out] flop  浮動小数点演算数
   */
  void vertex8(const int st[], const int ed[], REAL_TYPE* vf, double& flop);

};



class ShapeMonitor : public CompoFraction {

public:
  /** デフォルトコンストラクタ */
  ShapeMonitor() {}


  /** コンストラクタ
   * @param [in] size   ローカルセル数
   * @param [in] guide  ガイドセル数
   * @param [in] pch    ローカル領域セル幅
   * @param [in] org    ローカル領域基点座標
   */
  ShapeMonitor(const int* size, const int guide, const REAL_TYPE* pch, const REAL_TYPE* org)
  {
    this->size[0]  = size[0];
    this->size[1]  = size[1];
    this->size[2]  = size[2];
    this->guide    = guide;
    this->pch.x    = pch[0];
    this->pch.y    = pch[1];
    this->pch.z    = pch[2];
    this->org.x    = org[0];
    this->org.y    = org[1];
    this->org.z    = org[2];
  }


   /** デストラクタ */
  virtual ~ShapeMonitor() {}


public:
  /**
   * @brief セルの8頂点の内外判定より50%以上のセルにIDを設定する
   * @param [in]     st   開始インデクス
   * @param [in]     en   終了インデクス
   * @param [in,out] bcd  BCindex B
   * @param [in]     id   指定ID
   */
  void setID(const int st[], const int ed[], int* bcd, const int m_id);
};


#endif // _FB_FRACTION_H_
