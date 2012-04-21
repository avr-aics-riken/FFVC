/**@file
 * @brief 補助関数群 実装
 */

#include "CutUtil.h"

namespace cutlib {
namespace util {

namespace {

enum { X, Y, Z };

/// 2次元外積計算
inline float cross2d(const Vec2f& v1, const Vec2f& v2)
{
  return v1.x * v2.y - v1.y * v2.x;
}

/// 点pが三角形の内部にあるかの判定
/**
 * @param[in] n 正:v0,v1,v2が反時計周り/負:時計周り
 * @param[in] v0,v1,v2 三角形の頂点座標
 * @param[in] p 判定対象点
 * @return true:内部/false:外部
 */
bool existIntersection(float n, const Vec2f& v0, const Vec2f& v1, const Vec2f& v2, const Vec2f& p)
{
//#ifdef IGNORE_NORMAL_DIRECTION
#if (defined IGNORE_NORMAL_DIRECTION && !defined IGNORE_STL_NORMAL)
  if (cross2d(v1-v0,p-v0) > 0.0) {
    if (cross2d(v2-v1,p-v1) < 0.0 ||
        cross2d(v0-v2,p-v2) < 0.0) return false;
  }
  else {
    if (cross2d(v2-v1,p-v1) > 0.0 ||
        cross2d(v0-v2,p-v2) > 0.0) return false;
  }
#else
  if (n > 0.0) {
     if (cross2d(v1-v0,p-v0) < 0.0 ||
         cross2d(v2-v1,p-v1) < 0.0 ||
         cross2d(v0-v2,p-v2) < 0.0) return false;
  }
  else {
     if (cross2d(v1-v0,p-v0) > 0.0 ||
         cross2d(v2-v1,p-v1) > 0.0 ||
         cross2d(v0-v2,p-v2) > 0.0) return false;
  }
#endif
  return true;
}

#ifdef IGNORE_STL_NORMAL

Vec3f n;

Vec3f calcNormal(const Triangle* t)
{
  Vec3f* v = t->get_vertex();
  Vec3f a = v[1] - v[0];
  Vec3f b = v[2] - v[0];
  return cross(a,b).normalize();
}
#endif

} /* namespace ANONYMOUS */


/// 三角形とz軸に平行な直線との交点を計算
/**
 * @param[in] t 三角形
 * @param[in] x,y 直線座標
 * @param[out] z  交点座標
 * @return true:交点あり/false:交点なし
 */
bool intersectXY(const Triangle* t, float x, float y, float& z)
{
#ifndef IGNORE_STL_NORMAL
  Vec3f n = t->get_normal();
#endif
  Vec3f* v = t->get_vertex();
  z = 1.0;

  if (n[Z] == 0.0) return false;

  Vec2f v0(v[0][X], v[0][Y]);
  Vec2f v1(v[1][X], v[1][Y]);
  Vec2f v2(v[2][X], v[2][Y]);
  Vec2f p(x, y);
  if (!existIntersection(n[Z], v0, v1, v2, p)) return false;

  z = (dot(n,v[0]) - n[X]*x - n[Y]*y) / n[Z];

  return true;
}


/// 三角形とx軸に平行な直線との交点を計算
/**
 * @param[in] t 三角形
 * @param[in] y,z 直線座標
 * @param[out] x  交点座標
 * @return true:交点あり/false:交点なし
 */
bool intersectYZ(const Triangle* t, float y, float z, float& x)
{
#ifndef IGNORE_STL_NORMAL
  Vec3f n = t->get_normal();
#else
  n = calcNormal(t);
#endif
  Vec3f* v = t->get_vertex();
  x = 1.0;

  if (n[X] == 0.0) return false;

  Vec2f v0(v[0][Y], v[0][Z]);
  Vec2f v1(v[1][Y], v[1][Z]);
  Vec2f v2(v[2][Y], v[2][Z]);
  Vec2f p(y, z);
  if (!existIntersection(n[X], v0, v1, v2, p)) return false;

  x = (dot(n,v[0]) - n[Y]*y - n[Z]*z) / n[X];

  return true;
}


/// 三角形とy軸に平行な直線との交点を計算
/**
 * @param[in] t 三角形
 * @param[in] z,x 直線座標
 * @param[out] y  交点座標
 * @return true:交点あり/false:交点なし
 */
bool intersectZX(const Triangle* t, float z, float x, float& y)
{
#ifndef IGNORE_STL_NORMAL
  Vec3f n = t->get_normal();
#endif
  Vec3f* v = t->get_vertex();
  y = 1.0;

  if (n[Y] == 0.0) return false;

  Vec2f v0(v[0][Z], v[0][X]);
  Vec2f v1(v[1][Z], v[1][X]);
  Vec2f v2(v[2][Z], v[2][X]);
  Vec2f p(z, x);
  if (!existIntersection(n[Y], v0, v1, v2, p)) return false;

  y = (dot(n,v[0]) - n[Z]*z - n[X]*x) / n[Y];

  return true;
}


void errMsg(const char* fmt, ...)
{
  va_list ap;
  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);
  va_end(ap);
}


} /* namespace util */
} /* namespace cutlib */
