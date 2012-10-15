/**@file
 * @brief カスタムポリゴンリスト用クラス 実装
 */

#include "CutTriangle.h"

#ifdef CUTLIB_TIMING
#include "CutTiming.h"
#endif

#include <algorithm>

namespace cutlib {

namespace {

enum { X, Y, Z };

} /* namespace ANONYMOUS */

/// コンストラクタ
/**
 * @param[in] _t Polylib三角形ポリゴンクラス
 * @param[in] _bid 境界ID
 */
CutTriangle::CutTriangle(Triangle* _t, BidType _bid) : t(_t), bid(_bid)
{
  Vec3f* v = t->get_vertex();
  bboxMin[X] = std::min(std::min(v[0][X], v[1][X]), v[2][X]);
  bboxMin[Y] = std::min(std::min(v[0][Y], v[1][Y]), v[2][Y]);
  bboxMin[Z] = std::min(std::min(v[0][Z], v[1][Z]), v[2][Z]);
  bboxMax[X] = std::max(std::max(v[0][X], v[1][X]), v[2][X]);
  bboxMax[Y] = std::max(std::max(v[0][Y], v[1][Y]), v[2][Y]);
  bboxMax[Z] = std::max(std::max(v[0][Z], v[1][Z]), v[2][Z]);
}

/// 三角形が直方体領域と交わるかの判定
/**
 * @param[in] min,max 直方体頂点座標
 * @return true:交わる/false:交わらない
 */
bool CutTriangle::intersectBox(const Vec3f& min, const Vec3f& max)
{
  if (bboxMin[X] > max[X] || bboxMax[X] < min[X]) return false;
  if (bboxMin[Y] > max[Y] || bboxMax[Y] < min[Y]) return false;
  if (bboxMin[Z] > max[Z] || bboxMax[Z] < min[Z]) return false;
  return true;
}


/// Polylib検索メソッドの結果をカスタムリストに追加
void CutTriangle::AppendCutTriangles(CutTriangles& ctList, 
                                     const Polylib* pl, const CutBoundaries* bList,
                                     const Vec3f& min, const Vec3f& max)
{
  CutBoundaries::const_iterator b;
  for (b = bList->begin(); b != bList->end(); b++) {

#ifdef CUTLIB_TIMING
    Timing::Start("Polylib::search_polygons");
#endif
    Triangles* tList = pl->search_polygons(b->name, min, max, false);
#ifdef CUTLIB_TIMING
    Timing::Stop("Polylib::search_polygons");
#endif

    Triangles::const_iterator t;
    for (t = tList->begin(); t != tList->end(); t++) {
      CutTriangle* ct = new CutTriangle(*t, b->id);
      ctList.push_back(ct);
    }
  delete tList;
  }
}


/// 直方体領域と交わる三角形のリストをコピー
/**
 * @param[in] ctListFrom 三角形リスト コピー元
 * @param[out] ctListTo 三角形リスト コピー先
 * @param[in] center 直方体領域中心
 * @param[in] d 直方体領域幅の1/2
 */
void CutTriangle::CopyCutTriangles(const CutTriangles& ctListFrom, 
                                   CutTriangles& ctListTo,
                                   const Vec3f& center, const Vec3f& d)
{
  Vec3f min = center - d;
  Vec3f max = center + d;
  ctListTo.clear();
  CutTriangles::const_iterator ct;
  for (ct = ctListFrom.begin(); ct != ctListFrom.end(); ct++) {
    if ((*ct)->intersectBox(min, max)) ctListTo.push_back(*ct);
  }
}
                                          

/// リスト内の三角形オブジェクトを消去
void CutTriangle::DeleteCutTriangles(CutTriangles& ctList)
{
  CutTriangles::iterator ct;
  for (ct = ctList.begin(); ct != ctList.end(); ct++) {
    delete *ct;
  }
}


} /* namespace cutlib */
