/**@file
 * @brief 境界情報計算関数(コア部分) 宣言
 */

#ifndef CUTLIBCORE_H
#define CUTLIBCORE_H

#include "Cutlib.h"
#include "CutBoundary.h"
#include "CutTriangle.h"

#include "Polylib.h"
using namespace PolylibNS;

namespace cutlib {
namespace core {

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
                 float pos6[], BidType bid6[]);

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
                         const CutTriangles& ctList);


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
                          const Polylib* pl, const CutBoundaries* bList);

#endif //CUTLIB_OCTREE

} /* namespace core */
} /* namespace cutlib */

#endif /* CUTLIBCORE_H */
