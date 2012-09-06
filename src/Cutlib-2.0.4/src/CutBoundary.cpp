/**@file
 * @brief Polylibポリゴングループ--境界ID 対応付けクラス 実装
 */

#ifdef CUTLIB_DEBUG
#include <iostream>
#endif
#include "CutBoundary.h"

namespace cutlib {

/// Polylibポリゴングループ--境界ID 対応リスト生成
/**
 * @param pl Polylibクラスオブジェクト
 * @return Polylibポリゴングループ--境界ID 対応リスト
 */
CutBoundaries* CutBoundary::createCutBoundaryList(const Polylib* pl)
{
  CutBoundaries* bList = new CutBoundaries;
  std::vector<PolygonGroup *>* rootGroups = pl->get_root_groups();
  std::vector<PolygonGroup*>::iterator it = rootGroups->begin();;

  for (it = rootGroups->begin(); it != rootGroups->end(); it++) {
    checkPolygonGroup(*it, bList);
  }

  delete rootGroups;

  return bList;
}


/// ポリゴングループのチェック
/**
 * ポリゴングループのツリーを再帰的にたどり,
 * リーフノードかつIDが1以上のグループをリストに登録
 *
 * @param[in] pg Polylibポリゴングループ
 * @param[in,out] bList Polylibポリゴングループ--境界ID 対応リスト
 */
void CutBoundary::checkPolygonGroup(PolygonGroup *pg, CutBoundaries* bList)
{
  std::vector<PolygonGroup *>& childGroups = pg->get_children();
  if (childGroups.size() == 0) {
    if (pg->get_id() > 0) {
      bList->push_back(CutBoundary(pg->acq_fullpath(), pg->get_id()));
#ifdef CUTLIB_DEBUG
      std::cout << "add boundary: name:" << pg->acq_fullpath()
           << " id:" << pg->get_id() << std::endl;
#endif
    }
    return;
  }
  
  // 子グループに対して再帰呼び出し
  std::vector<PolygonGroup*>::iterator it;
  for (it = childGroups.begin(); it != childGroups.end(); it++) {
    checkPolygonGroup(*it, bList);
  }
}


} /* namespace cutlib */
