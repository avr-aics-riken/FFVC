/**@file
 * @brief 境界情報計算関数 実装
 */

#include "Cutlib.h"
#include "CutlibCore.h"
#include "CutUtil.h"
#include "CutTriangle.h"
#include "CutBoundary.h"

#ifdef CUTLIB_TIMING
#include "CutTiming.h"
#endif

#include <iostream>
using namespace std;

namespace cutlib {

namespace {

/// 配列サイズ, 計算対象領域のチェック
CutlibReturn checkSize(const char* func_name, const size_t ista[], const size_t nlen[],
              const CutPosArray* cutPosArray, const CutBidArray* cutBidArray)
{
  if (ista[0]+nlen[0] > cutPosArray->getSizeX() ||
      ista[1]+nlen[1] > cutPosArray->getSizeY() ||
      ista[2]+nlen[2] > cutPosArray->getSizeZ()) {
     util::errMsg("*** %s: ista[]+nlen[] exceeds array size.\n", func_name);
    return CL_SIZE_EXCEED;
  }
  return CL_SUCCESS;
}


#ifdef CUTLIB_OCTREE

/// SklTreeのチェック
CutlibReturn checkTree(const char* func_name, SklTree* tree)
{
  if (tree == 0) {
    util::errMsg("*** %s: SklTree not initialized.", func_name);
    return CL_BAD_SKLTREE;
  }
  return CL_SUCCESS;
}

#endif //CUTLIB_OCTREE

} /* namespace ANONYMOUS */


/// 交点情報計算: セル中心間, 配列ラッパクラス, 計算領域指定
/**
 * @param[in] ista 計算対象領域開始セル位置
 * @param[in] nlen 計算対象領域セル数
 * @param[in] org 領域原点座標
 * @param[in] d セル間隔
 * @param[in] pl Polylibクラスオブジェクト
 * @param[in,out] cutPos 交点座標配列ラッパ
 * @param[in,out] cutBid 境界ID配列ラッパ
 */
CutlibReturn CutInfoCell(const size_t ista[], const size_t nlen[],
                         const float org[], const float d[],
                         const Polylib* pl,
                         CutPosArray* cutPos, CutBidArray* cutBid)
{
  { // check input parameters
    CutlibReturn ret;
    ret = checkSize("CutInfoCell", ista, nlen, cutPos, cutBid);
    if (ret != CL_SUCCESS) return ret;
  }

  CutBoundaries* bList = CutBoundary::createCutBoundaryList(pl);

#ifdef CUTLIB_TIMING
  Timing::Start("Total");
#endif

  cutPos->clear();
  cutBid->clear();

  for (size_t k = ista[2]; k < ista[2]+nlen[2]; k++) {
    for (size_t j = ista[1]; j < ista[1]+nlen[1]; j++) {
      for (size_t i = ista[0]; i < ista[0]+nlen[0]; i++) {
        float pos6[6];
        BidType bid6[6];
        Vec3f center(org[0]+(i+0.5)*d[0], org[1]+(j+0.5)*d[1], org[2]+(k+0.5)*d[2]);

        core::calcCutInfo(center, d, pl, bList, pos6, bid6);

        cutPos->setPos(i, j, k, pos6);
        cutBid->setBid(i, j, k, bid6);
      }
    }
  }

#ifdef CUTLIB_TIMING
  Timing::Stop("Total");
  cout << endl;
  Timing::Print("Polylib::search_polygons");
  Timing::Print("Total");
#endif

  delete bList;

  return CL_SUCCESS;
}


/// 交点情報計算: ノード間, 配列ラッパクラス, 計算領域指定
/**
 * @param[in] ista 計算対象領域開始ノード位置
 * @param[in] nlen 計算対象領域ノード数
 * @param[in] org 領域原点座標
 * @param[in] d セル間隔
 * @param[in] pl Polylibクラスオブジェクト
 * @param[in,out] cutPos 交点座標配列ラッパ
 * @param[in,out] cutBid 境界ID配列ラッパ
 */
CutlibReturn CutInfoNode(const size_t ista[], const size_t nlen[],
                         const float org[], const float d[],
                         const Polylib* pl,
                         CutPosArray* cutPos, CutBidArray* cutBid)
{
  { // check input parameters
    CutlibReturn ret;
    ret = checkSize("CutInfoNode", ista, nlen, cutPos, cutBid);
    if (ret != CL_SUCCESS) return ret;
  }

  CutBoundaries* bList = CutBoundary::createCutBoundaryList(pl);

#ifdef CUTLIB_TIMING
  Timing::Start("Total");
#endif

  cutPos->clear();
  cutBid->clear();

  for (size_t k = ista[2]; k < ista[2]+nlen[2]; k++) {
    for (size_t j = ista[1]; j < ista[1]+nlen[1]; j++) {
      for (size_t i = ista[0]; i < ista[0]+nlen[0]; i++) {
        float pos6[6];
        BidType bid6[6];
        Vec3f center(org[0]+i*d[0], org[1]+j*d[1], org[2]+k*d[2]);

        core::calcCutInfo(center, d, pl, bList, pos6, bid6);

        cutPos->setPos(i, j, k, pos6);
        cutBid->setBid(i, j, k, bid6);
      }
    }
  }

#ifdef CUTLIB_TIMING
  Timing::Stop("Total");
  cout << endl;
  Timing::Print("Polylib::search_polygons");
  Timing::Print("Total");
#endif

  delete bList;

  return CL_SUCCESS;
}


/// 交点情報計算: セル中心間, 配列ラッパクラス, 全領域
/**
 * @param[in] org 領域原点座標
 * @param[in] d セル間隔
 * @param[in] pl Polylibクラスオブジェクト
 * @param[in,out] cutPos 交点座標配列ラッパ
 * @param[in,out] cutBid 境界ID配列ラッパ
 */
CutlibReturn CutInfoCell(const float org[], const float d[],
                         const Polylib* pl,
                         CutPosArray* cutPos, CutBidArray* cutBid)
{
  size_t ista[3], nlen[3];
  ista[0] = ista[1] = ista[2] = 0;
  nlen[0] = cutPos->getSizeX();
  nlen[1] = cutPos->getSizeY();
  nlen[2] = cutPos->getSizeZ();

  return CutInfoCell(ista, nlen, org, d, pl, cutPos, cutBid);
}


/// 交点情報計算: ノード間, 配列ラッパクラス, 全領域
/**
 * @param[in] org 領域原点座標
 * @param[in] d セル間隔
 * @param[in] pl Polylibクラスオブジェクト
 * @param[in,out] cutPos 交点座標配列ラッパ
 * @param[in,out] cutBid 境界ID配列ラッパ
 */
CutlibReturn CutInfoNode(const float org[], const float d[],
                         const Polylib* pl,
                         CutPosArray* cutPos, CutBidArray* cutBid)
{
  size_t ista[3], nlen[3];
  ista[0] = ista[1] = ista[2] = 0;
  nlen[0] = cutPos->getSizeX();
  nlen[1] = cutPos->getSizeY();
  nlen[2] = cutPos->getSizeZ();

  return CutInfoNode(ista, nlen, org, d, pl, cutPos, cutBid);
}


#ifdef CUTLIB_OCTREE

/// 交点情報計算: Octree, 計算対象セルタイプ指定
/**
 * @param[in,out] tree SklTreeクラスオブジェクト
 * @param[in] pl Polylibクラスオブジェクト
 * @param cutPos 交点座標データアクセッサ
 * @param[in] leafCellOnly 計算対象セルタイプフラグ
 */
CutlibReturn CutInfoOctree(SklTree* tree,
                           const Polylib* pl,
                           CutPosOctree* cutPos, CutBidOctree* cutBid,
                           bool leafCellOnly)
{
  if (leafCellOnly) {
    return CutInfoOctreeLeafCell(tree, pl, cutPos, cutBid);
  }
  else {
    return CutInfoOctreeAllCell(tree, pl, cutPos, cutBid);
  }
}


/// 交点情報計算: Octree, リーフセルのみ
/**
 * @param[in,out] tree SklTreeクラスオブジェクト
 * @param[in] pl Polylibクラスオブジェクト
 * @param cutPos 交点座標データアクセッサ
 * @param cutBid 境界IDデータアクセッサ
 */
CutlibReturn CutInfoOctreeLeafCell(SklTree* tree,
                                   const Polylib* pl,
                                   CutPosOctree* cutPos, CutBidOctree* cutBid)
{
  { // check input parameters
    CutlibReturn ret;
    ret = checkTree("CutInfoOctreeLeafCell", tree);
    if (ret != CL_SUCCESS) return ret;
  }

  CutBoundaries* bList = CutBoundary::createCutBoundaryList(pl);

#ifdef CUTLIB_TIMING
  Timing::Start("Total");
#endif

  for (SklCell* cell = tree->GetLeafCellFirst(); cell != 0;
       cell = tree->GetLeafCellNext(cell)) {

    cutPos->assignData(cell->GetData());
    cutBid->assignData(cell->GetData());

    Vec3f orig, d;
    cell->GetOrigin(orig[0], orig[1], orig[2]);
    cell->GetPitch(d[0], d[1], d[2]);
    Vec3f center = orig + d * 0.5;

    float pos6[6];
    BidType bid6[6];
    core::calcCutInfo(center, d, pl, bList, pos6, bid6);

    cutPos->setPos(pos6);
    cutBid->setBid(bid6);
  }

#ifdef CUTLIB_TIMING
  Timing::Stop("Total");
  cout << endl;
  Timing::Print("Polylib::search_polygons");
  Timing::Print("Total");
#endif

  delete bList;

  return CL_SUCCESS;
}


/// 交点情報計算: Octree, 全セル計算
/**
 * @param[in,out] tree SklTreeクラスオブジェクト
 * @param[in] pl Polylibクラスオブジェクト
 * @param cutPos 交点座標データアクセッサ
 * @param cutBid 境界IDデータアクセッサ
 */
CutlibReturn CutInfoOctreeAllCell(SklTree* tree,
                                  const Polylib* pl,
                                  CutPosOctree* cutPos, CutBidOctree* cutBid)
{
  { // check input parameters
    CutlibReturn ret;
    ret = checkTree("CutInfoOctreeAllCell", tree);
    if (ret != CL_SUCCESS) return ret;
  }

  CutBoundaries* bList = CutBoundary::createCutBoundaryList(pl);

#ifdef CUTLIB_TIMING
  Timing::Start("Total");
#endif

  size_t nx, ny, nz;
  tree->GetSize(nx, ny, nz);

  for (size_t k = 0; k < nz; k++) {
    for (size_t j = 0; j < ny; j++) {
      for (size_t i = 0; i < nx; i++) {
        SklCell* cell = tree->GetRootCell(i, j, k);
        Vec3f org, d;
        cell->GetOrigin(org[0], org[1], org[2]);
        cell->GetPitch(d[0], d[1], d[2]);
        Vec3f center = org + d * 0.5;
        Vec3f min = center - d;
        Vec3f max = center + d;

        CutTriangles ctList;
        CutTriangle::AppendCutTriangles(ctList, pl, bList, min, max);

        core::calcCutInfoTreeCell(cell, center, d, cutPos, cutBid, ctList);

        CutTriangle::DeleteCutTriangles(ctList);
      }
    }
  }

#ifdef CUTLIB_TIMING
  Timing::Stop("Total");
  cout << endl;
  Timing::Print("Polylib::search_polygons");
  Timing::Print("Total");
#endif

  delete bList;

  return CL_SUCCESS;
}


/// 交点情報計算: Octree, 全セル計算, デバッグ用
/**
 * 全セルでPolylibの検索メソッドを使用
 *
 * @param[in,out] tree SklTreeクラスオブジェクト
 * @param[in] pl Polylibクラスオブジェクト
 * @param cutPos 交点座標データアクセッサ
 * @param cutBid 境界IDデータアクセッサ
 */
CutlibReturn CutInfoOctreeAllCell0(SklTree* tree,
                                   const Polylib* pl,
                                   CutPosOctree* cutPos, CutBidOctree* cutBid)
{
  { // check input parameters
    CutlibReturn ret;
    ret = checkTree("CutInfoOctreeAllCell", tree);
    if (ret != CL_SUCCESS) return ret;
  }

  CutBoundaries* bList = CutBoundary::createCutBoundaryList(pl);

#ifdef CUTLIB_TIMING
  Timing::Start("Total");
#endif

  size_t nx, ny, nz;
  tree->GetSize(nx, ny, nz);

  for (size_t k = 0; k < nz; k++) {
    for (size_t j = 0; j < ny; j++) {
      for (size_t i = 0; i < nx; i++) {
        SklCell* cell = tree->GetRootCell(i, j, k);
        core::calcCutInfoTreeCell0(cell, cutPos, cutBid, pl, bList);
      }
    }
  }

#ifdef CUTLIB_TIMING
  Timing::Stop("Total");
  cout << endl;
  Timing::Print("Polylib::search_polygons");
  Timing::Print("Total");
#endif

  delete bList;

  return CL_SUCCESS;
}

#endif //CUTLIB_OCTREE

} /* namespace cutlib */
