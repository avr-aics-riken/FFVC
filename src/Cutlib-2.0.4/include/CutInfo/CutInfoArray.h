/**@file
 * @brief 交点情報配列ラッパクラス
 */

#ifndef CUTINFOARRAY_H
#define CUTINFOARRAY_H

#include "CutInfo.h"

namespace cutlib {

/**
 * @defgroup CutInfoArray 交点情報一次元配列ラッパクラス
 * @{
 */

/// 一次元配列ラッパクラスの基底クラス
class CutInfoArray {
  /// サイズ(実際はnx_*ny_*nz_の一次元配列を使用する)
  size_t nx_, ny_, nz_;

public:
  /// コンストラクタ
  /**
   * @param[in] nx,ny,nz 配列サイズ(3次元で指定)
   */
  CutInfoArray(size_t nx, size_t ny, size_t nz) : nx_(nx), ny_(ny), nz_(nz) {}

  /// デストラクタ
  virtual ~CutInfoArray() {}

  /// x方向のサイズを得る
  size_t getSizeX() const { return nx_; }

  /// y方向のサイズを得る
  size_t getSizeY() const { return ny_; }

  /// z方向のサイズを得る
  size_t getSizeZ() const { return nz_; }

protected:
  /// 3次元インデックス(i,j,k)より1次元インデックスを計算
  /**
   * @param[in] i,j,k  3次元インデックス
   */
  size_t index(size_t i, size_t j, size_t k) const { return i+j*nx_+k*nx_*ny_; }
};


/// 交点座標一次元配列ラッパ仮想クラス
class CutPosArray : public CutInfoArray {
public:
  /// コンストラクタ
  /**
   * @param[in] nx,ny,nz 配列サイズ(3次元で指定)
   */
  CutPosArray(size_t nx, size_t ny, size_t nz) : CutInfoArray(nx, ny, nz) {}

  /// デストラクタ
  virtual ~CutPosArray() {}

  /// 交点座標値を設定(d方向)
  /**
   * @param[in] i,j,k 3次元インデックス
   * @param[in] d 交点探査方向(0〜5) 
   * @param[in] pos 交点座標値
   */
  virtual void setPos(size_t i, size_t j, size_t k, int d, float pos) = 0;

  /// 交点座標値を設定(6方向まとめて)
  /**
   * @param[in] i,j,k 3次元インデックス
   * @param[in] pos 交点座標配列
   */
  virtual void setPos(size_t i, size_t j, size_t k, const float pos[]) = 0;

  /// 交点座標値(d方向)を得る
  /**
   * @param[in] i,j,k 3次元インデックス
   * @param[in] d 交点探査方向(0〜5) 
   * @return 交点座標値
   */
  virtual float getPos(size_t i, size_t j, size_t k, int d) const = 0;

  /// 交点座標値(d方向)を得る(1次元インデックスで指定)
  /**
   * @param[in] ijk 1次元インデックス
   * @param[in] d 交点探査方向(0〜5) 
   * @return 交点座標値
   */
  virtual float getPos(size_t ijk, int d) const = 0;

  /// 交点座標値(6方向まとめて)を得る
  /**
   * @param[in] i,j,k 3次元インデックス
   * @param[out] pos 交点座標配列
   */
  virtual void getPos(size_t i, size_t j, size_t k, float pos[]) const = 0;

  /// 交点座標値(6方向まとめて)を得る(1次元インデックスで指定)
  /**
   * @param[in] ijk 1次元インデックス
   * @param[out] pos 交点座標配列
   */
  virtual void getPos(size_t ijk, float pos[]) const = 0;

  /// 全配列データを1.0でクリア
  virtual void clear() = 0;
};


/// 境界ID一次元配列ラッパ仮想クラス
class CutBidArray : public CutInfoArray {
public:
  /// コンストラクタ
  /**
   * @param[in] nx,ny,nz  配列サイズ(3次元で指定)
   */
  CutBidArray(size_t nx, size_t ny, size_t nz) : CutInfoArray(nx, ny, nz) {}

  /// デストラクタ
  virtual ~CutBidArray() {}

  /// 境界IDを設定(d方向)
  /**
   * @param[in] i,j,k 3次元インデックス
   * @param[in] d 交点探査方向(0〜5) 
   * @param[in] bid 境界ID
   */
  virtual void setBid(size_t i, size_t j, size_t k, int d, BidType bid) = 0;

  /// 境界IDを設定(6方向まとてて)
  /**
   * @param[in] i,j,k 3次元インデックス
   * @param[in] bid 境界ID配列
   */
  virtual void setBid(size_t i, size_t j, size_t k, const BidType bid[]) = 0;

  /// 境界ID(d方向)を得る
  /**
   * @param[in] i,j,k 3次元インデックス
   * @param[in] d 交点探査方向(0〜5) 
   * @return 境界ID
   */
  virtual BidType getBid(size_t i, size_t j, size_t k, int d) const = 0;

  /// 境界ID(d方向)を得る(1次元インデックスで指定)
  /**
   * @param[in] ijk 1次元インデックス
   * @param[in] d 交点探査方向(0〜5) 
   * @return 境界ID
   */
  virtual BidType getBid(size_t ijk, int d) const = 0;

  /// 境界ID(6方向まとめて)を得る
  /**
   * @param[in] i,j,k 3次元インデックス
   * @param[out] bid  境界ID配列
   */
  virtual void getBid(size_t i, size_t j, size_t k, BidType bid[]) const = 0;

  /// 境界ID(6方向まとめて)を得る(1次元インデックスで指定)
  /**
   * @param[in] ijk 1次元インデックス
   * @param[out] bid  境界ID配列
   */
  virtual void getBid(size_t ijk, BidType bid[]) const = 0;

  /// 全配列データを0クリア
  virtual void clear() = 0;
};

//-----------------------------------------------------------------------------

/// 交点座標一次元配列ラッパクラステンプレート
template<typename CUT_POS>
class CutPosArrayTemplate : public CutPosArray {
  size_t n_;        ///< 一次元データサイズ
  CUT_POS* data_;   ///< 一次元データポインタ
  bool allocated_;  ///< 一次元データ管理フラグ

public:
  /// コンストラクタ(自前で一次元データ領域を確保)
  /**
   * デストラクタで一次元データ領域を開放(allocated_=true)
   * @param[in] nx,ny,nz  配列サイズ(3次元で指定)
   */
  CutPosArrayTemplate(size_t nx, size_t ny, size_t nz) : CutPosArray(nx, ny, nz)
  {
    n_ = nx * ny * nz;
    data_ = new CUT_POS[n_];
    allocated_ = true;
  }

  /// コンストラクタ(自前で一次元データ領域を確保)
  /**
   * デストラクタで一次元データ領域を開放(allocated_=true)
   * @param[in] ndim  配列サイズ(3次元で指定)
   */
  CutPosArrayTemplate(const size_t ndim[]) : CutPosArray(ndim[0], ndim[1], ndim[2])
  {
    n_ = ndim[0] * ndim[1] * ndim[2];
    data_ = new CUT_POS[n_];
    allocated_ = true;
  }

  /// コンストラクタ(一次元データ領域をインポート)
  /**
   * デストラクタで一次元データ領域を開放しない(allocated_=false)
   * @param[in] data 交点座標基本型配列
   * @param[in] nx,ny,nz  配列サイズ(3次元で指定)
   */
  CutPosArrayTemplate(CUT_POS* data, size_t nx, size_t ny, size_t nz) : CutPosArray(nx, ny, nz)
  {
    n_ = nx * ny * nz;
    data_ = data;
    allocated_ = false;
  }

  /// ディストラクタ
  /**
   * allocated_=trueの場合のみ一次元データ領域を開放
   */
  ~CutPosArrayTemplate() { if (allocated_) delete [] data_; }

  /// 交点座標値を設定(d方向)
  /**
   * @param[in] i,j,k 3次元インデックス
   * @param[in] d 交点探査方向(0〜5) 
   * @param[in] pos 交点座標値
   */
  void setPos(size_t i, size_t j, size_t k, int d, float pos)
  {
    SetCutPos(data_[index(i,j,k)], d, pos);
  }

  /// 交点座標値を設定(6方向まとめて)
  /**
   * @param[in] i,j,k 3次元インデックス
   * @param[in] pos 交点座標配列
   */
  void setPos(size_t i, size_t j, size_t k, const float pos[])
  {
    SetCutPos(data_[index(i,j,k)], pos);
  }

  /// 交点座標値(d方向)を得る
  /**
   * @param[in] i,j,k 3次元インデックス
   * @param[in] d 交点探査方向(0〜5) 
   * @return 交点座標値
   */
  float getPos(size_t i, size_t j, size_t k, int d) const
  {
    return GetCutPos(data_[index(i,j,k)], d);
  }

  /// 交点座標値(d方向)を得る(1次元インデックスで指定)
  /**
   * @param[in] ijk 1次元インデックス
   * @param[in] d 交点探査方向(0〜5) 
   * @return 交点座標値
   */
  float getPos(size_t ijk, int d) const
  {
    return GetCutPos(data_[ijk], d);
  }

  /// 交点座標値(6方向まとめて)を得る
  /**
   * @param[in] i,j,k 3次元インデックス
   * @param[out] pos 交点座標配列
   */
  void getPos(size_t i, size_t j, size_t k, float pos[]) const
  {
    GetCutPos(data_[index(i,j,k)], pos);
  }

  /// 交点座標値(6方向まとめて)を得る(1次元インデックスで指定)
  /**
   * @param[in] ijk 1次元インデックス
   * @param[out] pos 交点座標配列
   */
  void getPos(size_t ijk, float pos[]) const
  {
    GetCutPos(data_[ijk], pos);
  }

  /// 一次元配列データへのポインタを得る
  CUT_POS* getDataPointer() const { return data_; }


  /// 一次元配列データのサイズを得る
  size_t getDataSize() const { return n_; }

  /// 全配列データを1.0でクリア
  void clear()
  {
    for (size_t i = 0; i < n_; i++) ClearCutPos(data_[i]);
  }

};

//-----------------------------------------------------------------------------

/// 境界ID一次元配列ラッパクラステンプレート
template<typename CUT_BID>
class CutBidArrayTemplate : public CutBidArray {
  size_t n_;        ///< 一次元データサイズ
  CUT_BID* data_;   ///< 一次元データポインタ
  bool allocated_;  ///< 一次元データ管理フラグ

public:
  /// コンストラクタ(自前で一次元データ領域を確保)
  /**
   * デストラクタで一次元データ領域を開放(allocated_=true)
   * @param[in] nx,ny,nz  配列サイズ(3次元で指定)
   */
  CutBidArrayTemplate(size_t nx, size_t ny, size_t nz) : CutBidArray(nx, ny, nz)
  {
    n_ = nx * ny * nz;
    data_ = new CUT_BID[n_];
    allocated_ = true;
  }

  /// コンストラクタ(自前で一次元データ領域を確保)
  /**
   * デストラクタで一次元データ領域を開放(allocated_=true)
   * @param[in] ndim  配列サイズ(3次元で指定)
   */
  CutBidArrayTemplate(const size_t ndim[]) : CutBidArray(ndim[0], ndim[1], ndim[2])
  {
    n_ = ndim[0] * ndim[1] * ndim[2];
    data_ = new CUT_BID[n_];
    allocated_ = true;
  }

  /// コンストラクタ(一次元データ領域をインポート)
  /**
   * デストラクタで一次元データ領域を開放しない(allocated_=false)
   * @param[in] data 境界ID基本型配列
   * @param[in] nx,ny,nz  配列サイズ(3次元で指定)
   */
  CutBidArrayTemplate(CUT_BID* data, size_t nx, size_t ny, size_t nz) : CutBidArray(nx, ny, nz)
  {
    n_ = nx * ny * nz;
    data_ = data;
    allocated_ = false;
  }

  /// デストラクタ
  /**
   * allocated_=trueの場合のみ一次元データ領域を開放
   */
  ~CutBidArrayTemplate() { if (allocated_) delete [] data_; }

  /// 境界IDを設定(d方向)
  /**
   * @param[in] i,j,k 3次元インデックス
   * @param[in] d 交点探査方向(0〜5) 
   * @param[in] bid 境界ID
   */
  void setBid(size_t i, size_t j, size_t k, int d, BidType bid)
  {
    SetCutBid(data_[index(i,j,k)], d, bid);
  }

  /// 境界IDを設定(6方向まとてて)
  /**
   * @param[in] i,j,k 3次元インデックス
   * @param[in] bid 境界ID配列
   */
  void setBid(size_t i, size_t j, size_t k, const BidType bid[])
  {
    SetCutBid(data_[index(i,j,k)], bid);
  }

  /// 境界ID(d方向)を得る
  /**
   * @param[in] i,j,k 3次元インデックス
   * @param[in] d 交点探査方向(0〜5) 
   * @return 境界ID
   */
  BidType getBid(size_t i, size_t j, size_t k, int d) const
  {
    return GetCutBid(data_[index(i,j,k)], d);
  }

  /// 境界ID(d方向)を得る(1次元インデックスで指定)
  /**
   * @param[in] ijk 1次元インデックス
   * @param[in] d 交点探査方向(0〜5) 
   * @return 境界ID
   */
  BidType getBid(size_t ijk, int d) const
  {
    return GetCutBid(data_[ijk], d);
  }

  /// 境界ID(6方向まとめて)を得る
  /**
   * @param[in] i,j,k 3次元インデックス
   * @param[out] bid  境界ID配列
   */
  void getBid(size_t i, size_t j, size_t k, BidType bid[]) const
  {
    GetCutBid(data_[index(i,j,k)], bid);
  }

  /// 境界ID(6方向まとめて)を得る(1次元インデックスで指定)
  /**
   * @param[in] ijk 1次元インデックス
   * @param[out] bid  境界ID配列
   */
  void getBid(size_t ijk, BidType bid[]) const
  {
    GetCutBid(data_[ijk], bid);
  }

  /// 一次元配列データへのポインタを得る
  CUT_BID* getDataPointer() const { return data_; }

  /// 一次元配列データのサイズを得る
  size_t getDataSize() const { return n_; }

  /// 全配列データを0クリア
  void clear()
  {
    for (size_t i = 0; i < n_; i++) ClearCutBid(data_[i]);
  }
};

//-----------------------------------------------------------------------------

/// CutPos32型交点座標配列ラッパクラス
typedef CutPosArrayTemplate<CutPos32> CutPos32Array;

/// CutPos8型交点座標配列ラッパクラス
typedef CutPosArrayTemplate<CutPos8> CutPos8Array;

/// CutBid8型境界ID配列ラッパクラス
typedef CutBidArrayTemplate<CutBid8> CutBid8Array;

/// CutBid5型境界ID配列ラッパクラス
typedef CutBidArrayTemplate<CutBid5> CutBid5Array;

/** @} */ // end gropu CutInfoArray

} /* namespace cutlib */

#endif /* CUTINFOARRAY_H */
