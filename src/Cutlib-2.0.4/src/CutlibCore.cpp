/**@file
 * @brief 境界情報計算関数(コア部分) 実装
 */

#include "CutlibCore.h"
#include "CutUtil.h"
#include "CutTiming.h"

namespace cutlib {
namespace core {

namespace {

enum { X, Y, Z };

/// 指定されたポリゴンに対して交点情報を計算
void calcCutInofoTriangle(const Vec3f& c, const Vec3f& d,
                          const Triangle* t, BidType bid,
                          float pos6[], BidType bid6[])
{
  float p, pos;

  if (util::intersectYZ(t, c[Y], c[Z], p)) {
    if (p >= c[X]) {
      pos = (p - c[X])/ d[X];
      if (pos < 1.0 && pos < pos6[X_P]) { pos6[X_P] = pos; bid6[X_P] = bid; }
    }
    if (p <= c[X]) {
      pos = (c[X] - p) / d[X];
      if (pos < 1.0 && pos < pos6[X_M]) { pos6[X_M] = pos; bid6[X_M] = bid; }
    }
  }

  if (util::intersectZX(t, c[Z], c[X], p)) {
    if (p >= c[Y]) {
      pos = (p - c[Y])/ d[Y];
      if (pos < 1.0 && pos < pos6[Y_P]) { pos6[Y_P] = pos; bid6[Y_P] = bid; }
    }
    if (p <= c[Y]) {
      pos = (c[Y]- p) / d[Y];
      if (pos < 1.0 && pos < pos6[Y_M]) { pos6[Y_M] = pos; bid6[Y_M] = bid; }
    }
  }

  if (util::intersectXY(t, c[X], c[Y], p)) {
    if (p >= c[Z]) {
      pos = (p - c[Z])/ d[Z];
      if (pos < 1.0 && pos < pos6[Z_P]) { pos6[Z_P] = pos; bid6[Z_P] = bid; }
    }
    if (p <= c[Z]) {
      pos = (c[Z] - p) / d[Z];
      if (pos < 1.0 && pos < pos6[Z_M]) { pos6[Z_M] = pos; bid6[Z_M] = bid; }
    }
  }
}

} /* namespace ANONYMOUS */


/// 指定された計算点(セル中心orノード)で交点情報を計算
/**
 * @param[in] c 計算点座標
 * @param[in] d セルピッチ
 * @param[in] pl Polylibクラスオブジェクト
 * @param[in] bList (境界ID,ポリゴングループ名)対応リスト
 * @param[out] pos6 交点座標値配列
 * @param[out] bid6 境界ID配列
 */
void calcCutInfo(const Vec3f& c, const Vec3f& d,
                 const Polylib* pl, const CutBoundaries* bList,
                 float pos6[], BidType bid6[])
{
  Vec3f min = c - d;
  Vec3f max = c + d;
  for (int i = 0; i < 6; i++) {
    pos6[i] = 1.0;
    bid6[i] = 0;
  }

  CutBoundaries::const_iterator b;
  for (b = bList->begin(); b != bList->end(); b++) {

    BidType bid = b->id;
    std::string name = b->name;

#ifdef CUTLIB_TIMING
    Timing::Start("Polylib::search_polygons");
#endif
    Triangles* tList = pl->search_polygons(name, min, max, false);
#ifdef CUTLIB_TIMING
    Timing::Stop("Polylib::search_polygons");
#endif

    Triangles::const_iterator t;
    for (t = tList->begin(); t != tList->end(); t++) {
       calcCutInofoTriangle(c, d, *t, bid, pos6, bid6);
    } 
    delete tList;

  } /* b */
}


#ifdef CUTLIB_OCTREE

/// Octree上のセルでの交点情報を計算
/**
 * 再帰的に呼び出される
 * @param[in,out]  cell  SklCellセル
 * @param[in] center セル中心座標
 * @param[in] d セルピッチ
 * @param cutPos 交点座標データアクセッサ
 * @param cutBid 境界IDデータアクセッサ
 * @param[in] ctList ポリゴンリスト
 */
void calcCutInfoTreeCell(SklCell* cell, const Vec3f& center, const Vec3f& d,
                         CutPosOctree* cutPos, CutBidOctree* cutBid,
                         const CutTriangles& ctList)
{
#ifdef CUTLIB_DEBUG
  std::cout << ctList.size() << "@" << cell->GetMyLevel() << std::endl;
#endif

  cutPos->assignData(cell->GetData());
  cutBid->assignData(cell->GetData());

  float pos6[6];
  BidType bid6[6];
  for (int i = 0; i < 6; i++) {
    pos6[i] = 1.0;
    bid6[i] = 0;
  }
  CutTriangles::const_iterator ct;
  for (ct = ctList.begin(); ct != ctList.end(); ct++) {
    calcCutInofoTriangle(center, d, (*ct)->t, (*ct)->bid, pos6, bid6);
  }

  cutPos->setPos(pos6);
  cutBid->setBid(bid6);

  if (cell->hasChild()) {
    for (TdPos p = 0; p < 8; p++) {
      SklCell* cellChild = cell->GetChildCell(p);
      Vec3f orgChild, dChild;
      cellChild->GetOrigin(orgChild[0], orgChild[1], orgChild[2]);
      cellChild->GetPitch(dChild[0], dChild[1], dChild[2]);
      Vec3f centerChild = orgChild + dChild * 0.5;

      CutTriangles ctListChild;
      CutTriangle::CopyCutTriangles(ctList, ctListChild, centerChild, dChild);
      calcCutInfoTreeCell(cellChild, centerChild, dChild,
                          cutPos, cutBid, ctListChild);
    }
  }
}


/// Octree上のセルでの交点情報を計算(デバッグ用)
/**
 * Polylibの検索メソッドを使用，再帰的に呼び出される
 * @param[in,out]  cell  SklCellセル
 * @param cutPos 交点座標データアクセッサ
 * @param cutBid 境界IDデータアクセッサ
 * @param[in] pl Polylibクラスオブジェクト
 * @param[in] bList (境界ID,ポリゴングループ名)対応リスト
 */
void calcCutInfoTreeCell0(SklCell* cell,
                          CutPosOctree* cutPos, CutBidOctree* cutBid,
                          const Polylib* pl, const CutBoundaries* bList)
{

  cutPos->assignData(cell->GetData());
  cutBid->assignData(cell->GetData());

  Vec3f orig, d;
  cell->GetOrigin(orig[0], orig[1], orig[2]);
  cell->GetPitch(d[0], d[1], d[2]);
  Vec3f center = orig + d * 0.5;

  float pos6[6];
  BidType bid6[6];
  calcCutInfo(center, d, pl, bList, pos6, bid6);

  cutPos->setPos(pos6);
  cutBid->setBid(bid6);

//cout << "L(" << cell->GetMyLevel() << ")P(" << cell->GetMyPos() << ") ";

  if (cell->hasChild()) {
    for (TdPos pos = 0; pos < 8; pos++) {
      SklCell* childCell = cell->GetChildCell(pos);
      calcCutInfoTreeCell0(childCell, cutPos, cutBid, pl, bList);
    }
  }

}

#endif //CUTLIB_OCTREE


} /* namespace core */
} /* namespace cutlib */
