/**@file
 * @brief カスタムポリゴンリスト用クラス 宣言
 */

#ifndef CUTTRIANGLE_H
#define CUTTRIANGLE_H

#include "Cutlib.h"
#include "CutBoundary.h"

namespace cutlib {

class CutTriangle;

/// 三角形リスト
typedef std::vector<CutTriangle*> CutTriangles;

/// BBox(binding box)情報と境界IDを持つカスタムポリゴンクラス
class CutTriangle {
public:
  Triangle* t;    ///< Polylib三角形ポリゴンクラス
  BidType bid;    ///< 境界ID
  Vec3f bboxMin; ///< BBox最小値
  Vec3f bboxMax; ///< BBox最大値

  /// コンストラクタ
  /**
   * @param[in] _t Polylib三角形ポリゴンクラス
   * @param[in] _bid 境界ID
   */
  CutTriangle(Triangle* _t, BidType _bid);

  /// 三角形が直方体領域と交わるかの判定
  /**
   * @param[in] min,max 直方体頂点座標
   * @return true:交わる/false:交わらない
   */
  bool intersectBox(const Vec3f& min, const Vec3f& max);

  /// Polylib検索メソッドの結果をカスタムリストに追加
  /**
   * @param[in,out] ctList 三角形リスト
   * @param[in] pl Polylibクラスオブジェクト
   * @param[in] bList (境界ID,ポリゴングループ名)対応リスト
   * @param[in] min,max 検索領域
   */
  static void AppendCutTriangles(CutTriangles& ctList, 
                                 const Polylib* pl, const CutBoundaries* bList,
                                 const Vec3f& min, const Vec3f& max);

  /// 直方体領域と交わる三角形のリストをコピー
  /**
   * @param[in] ctListFrom 三角形リスト コピー元
   * @param[out] ctListTo 三角形リスト コピー先
   * @param[in] center 直方体領域中心
   * @param[in] d 直方体領域幅の1/2
   */
  static void CopyCutTriangles(const CutTriangles& ctListFrom, 
                               CutTriangles& ctListTo,
                               const Vec3f& center, const Vec3f& d);

  /// リスト内の三角形オブジェクトを消去
  /**
   * @param[in,out] ctList 三角形リスト
   */
  static void DeleteCutTriangles(CutTriangles& ctList);
};


} /* namespace cutlib */

#endif /* CUTTRIANGLE_H */
