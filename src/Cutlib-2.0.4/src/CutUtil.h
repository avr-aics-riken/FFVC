/**@file
 * @brief 補助関数群 宣言
 */

#ifndef CUTUTIL_H
#define CUTUTIL_H

#include <cstdio>
#include <cstdarg>

#include "Cutlib.h"

namespace cutlib {
namespace util {

/// エラーメッセージ出力
void errMsg(const char* fmt, ...);

/// 三角形とz軸に平行な直線との交点を計算
/**
 * @param[in] t 三角形
 * @param[in] x,y 直線座標
 * @param[out] z  交点座標
 * @return true:交点あり/false:交点なし
 */
bool intersectXY(const Triangle* t, float x, float y, float& z);

/// 三角形とx軸に平行な直線との交点を計算
/**
 * @param[in] t 三角形
 * @param[in] y,z 直線座標
 * @param[out] x  交点座標
 * @return true:交点あり/false:交点なし
 */
bool intersectYZ(const Triangle* t, float y, float z, float& x);

/// 三角形とy軸に平行な直線との交点を計算
/**
 * @param[in] t 三角形
 * @param[in] z,x 直線座標
 * @param[out] y  交点座標
 * @return true:交点あり/false:交点なし
 */
bool intersectZX(const Triangle* t, float z, float x, float& y);

} /* namespace util */
} /* namespace cutlib */

#endif /* CUTUTIL_H */
