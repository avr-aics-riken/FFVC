/**@file
 * @brief 交点情報基本データ型
 */

#ifndef CUTINFO_H
#define CUTINFO_H

#include <cstddef>   /* for size_t */
#include <stdint.h>  /* for int32_t */

namespace cutlib {

/// 交点探査方向の順番
/**
 * X_M(x負方向), X_P(x正方向), ...
 * (2010/10/18 正負の順を変更)
 */
enum CutInfoOrder { X_M, X_P, Y_M, Y_P, Z_M, Z_P };


/// 境界ID型
typedef unsigned char BidType;

/**
 * @defgroup CutInfo 交点情報基本型
 * @{
 */
//-----------------------------------------------------------------------------

/// 交点座標基本型: 交点座標をfloat(32ビット)として格納
typedef float   CutPos32[6]; 

/// 交点座標基本型: 交点座標を8ビット量子化して，3つずつ，2つの32ビット整数に格納
typedef int32_t CutPos8[2];

/// 境界ID基本型: 0〜255の境界ID(8ビット)を3つずつ，2つの32ビット整数に格納
typedef int32_t CutBid8[2];

/// 境界ID基本型: 0〜31の境界ID(5ビット)を6つまとめて，1つの32ビット整数に格納
//typedef int32_t CutBid5[1];
typedef int32_t CutBid5;

//-----------------------------------------------------------------------------

/// 交点座標値を1.0でクリア
/**
 * @param[out] cp 交点座標基本型
 */
inline void ClearCutPos(CutPos32& cp)
{
  for (int d = 0; d < 6; d++) cp[d] = 1.0f;
}

/// 交点座標値を設定(d方向)
/**
 * @param[out] cp 交点座標基本型
 * @param[in] d 交点探査方向(0〜5)
 * @param[in] pos 交点座標値
 */
inline void SetCutPos(CutPos32& cp, int d, float pos)
{
  cp[d] = pos;
}

/// 交点座標値を設定(6方向まとめて)
/**
 * @param[out] cp 交点座標基本型
 * @param[in] pos 交点座標値配列
 */
inline void SetCutPos(CutPos32& cp, const float pos[])
{
  for (int d = 0; d < 6; d++) cp[d] = pos[d];
}

/// 交点座標値(d方向)を得る
/**
 * @param[in] cp 交点座標基本型
 * @param[in] d 交点探査方向(0〜5)
 * @return 交点座標値
 */
inline float GetCutPos(const CutPos32& cp, int d)
{
  return cp[d];
}

/// 交点座標値(6方向まとめて)を得る
/**
 * @param[in] cp 交点座標基本型
 * @param[out] pos 交点座標値配列
 */
inline void GetCutPos(const CutPos32& cp, float pos[])
{
  for (int d = 0; d < 6; d++) pos[d] = cp[d];
}

//-----------------------------------------------------------------------------

/// クリアマスク(方向別)
const int32_t CutPos8Clear0[3] = {
   ~(int32_t)255,
   ~((int32_t)255 << 8),
   ~((int32_t)255 << 16)
};

/// クリアマスク(6方向まとめて)
const int32_t CutPos8Clear = (int32_t)255 | ((int32_t)255<<8) | ((int32_t)255<<16);

/// 交点座標値を1.0でクリア
/**
 * @param[out] cp 交点座標基本型
 */
inline void ClearCutPos(CutPos8& cp)
{
  cp[0] = cp[1] = CutPos8Clear;
}

/// 交点座標値を設定(d方向)
/**
 * @param[out] cp 交点座標基本型
 * @param[in] d 交点探査方向(0〜5)
 * @param[in] pos 交点座標値
 */
inline void SetCutPos(CutPos8& cp, int d, float pos)
{
  int i =  (d < 3) ? 0 : 1;
  cp[i] &= CutPos8Clear0[d%3];
  cp[i] |= (int32_t)(pos*255) << (d%3)*8;
}

/// 交点座標値を設定(6方向まとめて)
/**
 * @param[out] cp 交点座標基本型
 * @param[in] pos 交点座標値配列
 */
inline void SetCutPos(CutPos8& cp, const float pos[])
{
  cp[0] = (int32_t)(pos[0]*255)
        | ((int32_t)(pos[1]*255) << 8)
        | ((int32_t)(pos[2]*255) << 16);
  cp[1] = (int32_t)(pos[3]*255)
        | ((int32_t)(pos[4]*255) << 8)
        | ((int32_t)(pos[5]*255) << 16);
}

/// 交点座標値(d方向)を得る
/**
 * @param[in] cp 交点座標基本型
 * @param[in] d 交点探査方向(0〜5)
 * @return 交点座標値
 */
inline float GetCutPos(const CutPos8& cp, int d)
{
  int i =  (d < 3) ? 0 : 1;
  return (float)((cp[i] >> (d%3)*8) & 255) / 255;
}

/// 交点座標値(6方向まとめて)を得る
/**
 * @param[in] cp 交点座標基本型
 * @param[out] pos 交点座標値配列
 */
inline void GetCutPos(const CutPos8& cp, float pos[])
{
  int32_t qpos3 = cp[0];
  pos[0] = (float)(qpos3 & 255) / 255; qpos3 >>= 8;
  pos[1] = (float)(qpos3 & 255) / 255; qpos3 >>= 8;
  pos[2] = (float)(qpos3 & 255) / 255;
  qpos3 = cp[1];
  pos[3] = (float)(qpos3 & 255) / 255; qpos3 >>= 8;
  pos[4] = (float)(qpos3 & 255) / 255; qpos3 >>= 8;
  pos[5] = (float)(qpos3 & 255) / 255;
}

//-----------------------------------------------------------------------------

/// クリアマスク(方向別)
const int32_t CutBid8Clear0[3] = {
   ~(int32_t)255,
   ~((int32_t)255 << 8),
   ~((int32_t)255 << 16)
};

/// 境界IDを0クリア
/**
 * @param[out] cb 境界ID基本型
 */
inline void ClearCutBid(CutBid8& cb)
{
  cb[0] = cb[1] = 0;
}

/// 境界IDを設定(d方向)
/**
 * @param[out] cb 境界ID基本型
 * @param[in] d 交点探査方向(0〜5)
 * @param[in] bid 境界ID
 */
inline void SetCutBid(CutBid8& cb, int d, BidType bid)
{
  int i =  (d < 3) ? 0 : 1;
  cb[i] &= CutBid8Clear0[d%3];
  cb[i] |= (int32_t)bid << (d%3)*8;
}

/// 境界IDを設定(6方向まとてて)
/**
 * @param[out] cb 境界ID基本型
 * @param[in] bid 境界ID配列
 */
inline void SetCutBid(CutBid8& cb, const BidType bid[])
{
  cb[0] = (int32_t)bid[0]
        | ((int32_t)bid[1] << 8)
        | ((int32_t)bid[2] << 16);
  cb[1] = (int32_t)bid[3]
        | ((int32_t)bid[4] << 8)
        | ((int32_t)bid[5] << 16);
}

/// 境界ID(d方向)を得る
/**
 * @param[in] cb 境界ID基本型
 * @param[in] d 交点探査方向(0〜5)
 * @return 境界ID
 */
inline BidType GetCutBid(const CutBid8& cb, int d)
{
  int i =  (d < 3) ? 0 : 1;
  return (cb[i] >> (d%3)*8) & 255;
}

/// 境界ID(6方向まとめて)を得る
/**
 * @param[in] cb 境界ID基本型
 * @param[out] bid  境界ID配列
 */
inline void GetCutBid(const CutBid8& cb, BidType bid[])
{
  int32_t bid3 = cb[0];
  bid[0] = bid3 & 255; bid3 >>= 8;
  bid[1] = bid3 & 255; bid3 >>= 8;
  bid[2] = bid3 & 255;
  bid3 = cb[1];
  bid[3] = bid3 & 255; bid3 >>= 8;
  bid[4] = bid3 & 255; bid3 >>= 8;
  bid[5] = bid3 & 255;
}

//-----------------------------------------------------------------------------

/// クリアマスク(方向別)
const int32_t CutBid5Clear0[6] = {
   ~(int32_t)31,
   ~((int32_t)31 << 5),
   ~((int32_t)31 << 10),
   ~((int32_t)31 << 15),
   ~((int32_t)31 << 20),
   ~((int32_t)31 << 25),
};

/// 境界IDを0クリア
/**
 * @param[out] cb 境界ID基本型
 */
inline void ClearCutBid(CutBid5& cb)
{
  cb = 0;
}

/// 境界IDを設定(d方向)
/**
 * @param[out] cb 境界ID基本型
 * @param[in] d 交点探査方向(0〜5)
 * @param[in] bid 境界ID
 */
inline void SetCutBid(CutBid5& cb, int d, BidType bid)
{
  cb &= CutBid5Clear0[d];
  cb |= (int32_t)bid << d*5;
}

/// 境界IDを設定(6方向まとてて)
/**
 * @param[out] cb 境界ID基本型
 * @param[in] bid 境界ID配列
 */
inline void SetCutBid(CutBid5& cb, const BidType bid[])
{
  cb = (int32_t)bid[0] | ((int32_t)bid[1] << 5)
     | ((int32_t)bid[2] << 10) | ((int32_t)bid[3] << 15)
     | ((int32_t)bid[4] << 20) | ((int32_t)bid[5] << 25);
}

/// 境界ID(d方向)を得る
/**
 * @param[in] cb 境界ID基本型
 * @param[in] d 交点探査方向(0〜5)
 * @return 境界ID
 */
inline BidType GetCutBid(const CutBid5& cb, int d)
{
  return (cb >> d*5) & 31;
}

/// 境界ID(6方向まとめて)を得る
/**
 * @param[in] cb 境界ID基本型
 * @param[out] bid  境界ID配列
 */
inline void GetCutBid(const CutBid5& cb, BidType bid[])
{
  int32_t bid6 = cb;
  for (int d = 0; d < 6; d++) {
     bid[d] = bid6 & 31;
     bid6 >>= 5;
  }
}

//-----------------------------------------------------------------------------

/// 一次元配列データから交点座標を得る(d方向)
/**
 * @param[in] cp 交点座標基本型配列
 * @param[in] ijk 一次元インデックス
 * @param[in] d 交点探査方向(0〜5)
 * @return 交点座標値
 */
template<typename CUT_POS>
inline float GetCutPos(const CUT_POS* cp, size_t ijk, int d)
{
  return GetCutPos(cp[ijk], d);
}

/// 一次元配列データから交点座標を得る(6方向まとめて)
/**
 * @param[in] cp 交点座標基本型配列
 * @param[in] ijk 一次元インデックス
 * @param[out] pos 交点座標値配列
 */
template<typename CUT_POS>
inline void GetCutPos(const CUT_POS* cp, size_t ijk, float pos[])
{
  GetCutPos(cp[ijk], pos);
}

/// 一次元配列データから境界IDを得る(d方向)
/**
 * @param[in] cb  境界ID基本型配列
 * @param[in] ijk 一次元インデックス
 * @param[in] d 交点探査方向(0〜5)
 * @return 境界ID
 */
template<typename CUT_BID>
inline BidType GetCutBid(const CUT_BID* cb, size_t ijk, int d)
{
  return GetCutBid(cb[ijk], d);
}

/// 一次元配列データから境界IDを得る(6方向まとめて)
/**
 * @param[in] cb  境界ID基本型配列
 * @param[in] ijk 一次元インデックス
 * @param[out] bid 境界ID配列
 */
template<typename CUT_BID>
inline void GetCutBid(const CUT_BID* cb, size_t ijk, BidType bid[])
{
  GetCutBid(cb[ijk], bid);
}

/** @} */ // end gropu CutInfo

} /* namespace cutlib */

#endif /* CUTINFOE_H */
