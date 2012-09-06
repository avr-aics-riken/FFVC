/**@file
 * @brief Octreeセルデータへのアクセッサクラス
 */

#ifndef CUTINFOOCTREE_H
#define CUTINFOOCTREE_H

#include "CutInfo.h"

namespace cutlib {

/**
 * @defgroup CutInfoOctree 交点情報Octreeセルデータアクセッサクラス
 * @{
 */

/// Octree交点座標データアクセッサ仮想クラス
class CutPosOctree {
protected:
  /// SklCellデータ領域内での格納開始インデックス
  int index_;

public:
  /// コンストラクタ
  /**
   * @param[in] index SklCellデータ領域内での格納開始インデックス
   */
  CutPosOctree(int index) { index_ = index; }

  /// デストラクタ
  virtual ~CutPosOctree() {}

  /// SklCellデータ領域に結合
  /**
   * @param[in] data SklCellデータ領域ポインタ
   */
  virtual void assignData(float* data) = 0;

  /// 交点座標値を設定(d方向)
  /**
   * @param[in] d 交点探査方向(0〜5) 
   * @param[in] pos 交点座標値
   */
  virtual void setPos(int d, float pos) = 0;

  /// 交点座標値を設定(6方向まとめて)
  /**
   * @param[in] pos 交点座標値配列
   */
  virtual void setPos(const float pos[]) = 0;

  /// 交点座標値(d方向)を得る
  /**
   * @param[in] d 交点探査方向(0〜5) 
   * @return 交点座標値
   */
  virtual float getPos(int d) const = 0;

  /// 交点座標値(6方向まとめて)を得る
  /**
   * @param[out] pos 交点座標値配列
   */
  virtual void getPos(float pos[]) const = 0;

  /// 6方向の交点座標値を1.0でクリア
  virtual void clear() = 0;

  /// 必要な格納領域サイズをfloat単位で得る
  virtual unsigned getSizeInFloat() const = 0;
};

/// Octree境界IDデータアクセッサ仮想クラス
class CutBidOctree {
protected:
  /// SklCellデータ領域内での格納開始インデックス
  int index_;

public:
  /// コンストラクタ
  /**
   * @param[in] index SklCellデータ領域内での格納開始インデックス
   */
  CutBidOctree(int index) { index_ = index; }

  /// デストラクタ

  virtual ~CutBidOctree() {}

  /// SklCellデータ領域に結合
  /**
   * @param[in] data SklCellデータ領域ポインタ
   */
  virtual void assignData(float* data) = 0;

  /// 境界IDを設定(d方向)
  /**
   * @param[in] d 交点探査方向(0〜5) 
   * @param[in] bid 境界ID
   */
  virtual void setBid(int d, BidType bid) = 0;

  /// 境界IDを設定(6方向まとてて)
  /**
   * @param[in] bid 境界ID配列
   */
  virtual void setBid(const BidType bid[]) = 0;

  /// 境界ID(d方向)を得る
  /**
   * @param[in] d 交点探査方向(0〜5) 
   * @return 境界ID
   */
  virtual BidType getBid(int d) const = 0;

  /// 境界ID(6方向まとめて)を得る
  /**
   * @param[out] bid 境界ID配列
   */
  virtual void getBid(BidType bid[]) const = 0;

  /// 6方向の境界IDを0クリア
  virtual void clear() = 0;

  /// 必要な格納領域サイズをfloat単位で得る
  virtual unsigned getSizeInFloat() const = 0;
};

//-----------------------------------------------------------------------------

/// Octree交点座標データアクセサクラステンプレート
template<typename CUT_POS, unsigned SIZE_IN_FLOAT>
class CutPosOctreeTemplate : public CutPosOctree {
  /// SklCellデータ領域ポインタ
  CUT_POS* data_;

public:
  ///  必要な格納領域サイズ(float単位)
  /**
   * テンプレートパラメータとして指定
   */
  static const unsigned SizeInFloat = SIZE_IN_FLOAT;

  /// 必要な格納領域サイズをfloat単位で得る
  unsigned getSizeInFloat() const { return SizeInFloat; }

  /// コンストラクタ
  /**
   * @param[in] index SklCellデータ領域内での格納開始インデックス
   */
  CutPosOctreeTemplate(int index) : CutPosOctree(index) {}

  /// デストラクタ
  virtual ~CutPosOctreeTemplate() {}

  /// SklCellデータ領域に結合
  /**
   * @param[in] data SklCellデータ領域ポインタ
   */
  void assignData(float* data) { data_ = (CUT_POS*)&data[index_]; }

  /// 交点座標値を設定(d方向)
  /**
   * @param[in] d 交点探査方向(0〜5) 
   * @param[in] pos 交点座標値
   */
  void setPos(int d, float pos) { SetCutPos(*data_, d, pos); }

  /// 交点座標値を設定(6方向まとめて)
  /**
   * @param[in] pos 交点座標値配列
   */
  void setPos(const float pos[]) { SetCutPos(*data_, pos); }

  /// 交点座標値(d方向)を得る
  /**
   * @param[in] d 交点探査方向(0〜5) 
   * @return 交点座標値
   */
  float getPos(int d) const { return GetCutPos(*data_, d); }

  /// 交点座標値(6方向まとめて)を得る
  /**
   * @param[out] pos 交点座標値配列
   */
  void getPos(float pos[]) const { GetCutPos(*data_, pos); }

  /// 6方向の交点座標値を1.0でクリア
  void clear() { ClearCutPos(*data_); }
};

//-----------------------------------------------------------------------------

/// Octree境界IDデータアクセサクラステンプレート
template<typename CUT_BID, unsigned SIZE_IN_FLOAT>
class CutBidOctreeTemplate : public CutBidOctree {
  /// SklCellデータ領域ポインタ
  CUT_BID* data_;

public:
  ///  必要な格納領域サイズ(float単位)
  /**
   * テンプレートパラメータとして指定
   */
  static const unsigned SizeInFloat = SIZE_IN_FLOAT;

  /// 必要な格納領域サイズをfloat単位で得る
  unsigned getSizeInFloat() const { return SizeInFloat; }

  /// コンストラクタ
  /**
   * @param[in] index SklCellデータ領域内での格納開始インデックス
   */
  CutBidOctreeTemplate(int index) : CutBidOctree(index) {}

  /// デストラクタ
  virtual ~CutBidOctreeTemplate() {}

  /// SklCellデータ領域に結合
  /**
   * @param[in] data SklCellデータ領域ポインタ
   */
  void assignData(float* data) { data_ = (CUT_BID*)&data[index_]; }

  /// 境界IDを設定(d方向)
  /**
   * @param[in] d 交点探査方向(0〜5) 
   * @param[in] bid 境界ID
   */
  void setBid(int d, BidType bid) { SetCutBid(*data_, d, bid); }

  /// 境界IDを設定(6方向まとてて)
  /**
   * @param[in] bid 境界ID配列
   */
  void setBid(const BidType bid[]) { SetCutBid(*data_, bid); }

  /// 境界ID(d方向)を得る
  /**
   * @param[in] d 交点探査方向(0〜5) 
   * @return 境界ID
   */
  BidType getBid(int d) const { return GetCutBid(*data_, d); }

  /// 境界ID(6方向まとめて)を得る
  /**
   * @param[out] bid 境界ID配列
   */
  void getBid(BidType bid[]) const { GetCutBid(*data_, bid); }

  /// 6方向の境界IDを0クリア
  void clear() { ClearCutBid(*data_); }
};

//-----------------------------------------------------------------------------

/// CutPos32型交点座標データアクセッサクラス
typedef CutPosOctreeTemplate<CutPos32, 6> CutPos32Octree;

/// CutPos8型交点座標データアクセッサクラス
typedef CutPosOctreeTemplate<CutPos8, 2> CutPos8Octree;

/// CutBid8型境界IDデータアクセッサクラス
typedef CutBidOctreeTemplate<CutBid8, 2> CutBid8Octree;

/// CutBid5型境界IDデータアクセッサクラス
typedef CutBidOctreeTemplate<CutBid5, 1> CutBid5Octree;

/** @} */ // end gropu CutInfoOctree

} /* namespace cutlib */

#endif /* CUTINFOOCTREE_H */
