/*
 * Cutlib Version 2.0.2,  Nov.10.2010
 */

/**@file
 * @brief 境界情報計算関数 宣言
 */

#ifndef CUTLIB_H
#define CUTLIB_H

#include "Polylib.h"
using namespace PolylibNS;

#include "CutInfo/CutInfoArray.h"

#ifdef CUTLIB_OCTREE
#include "SklCompatibility.h"
#include "CutInfo/CutInfoOctree.h"
#endif

#include <string>
#include <vector>
#include <list>

namespace cutlib {

/// Polylib検索結果格納リスト
//typedef std::list<Triangle*> Triangles;
typedef std::vector<Triangle*> Triangles;

/**
 * @defgroup CutInfoFunc 交点情報計算関数
 * @{
 */
/// 交点情報計算関数リターンコード
enum CutlibReturn {
  CL_SUCCESS = 0,         ///< 成功
  CL_BAD_GROUP_LIST = 1,  ///< (境界ID,ポリゴングループ名)対応リストが不正
  CL_BAD_POLYLIB = 2,     ///< Polylibオブジェクトが不正(未初期化等)
  CL_BAD_SKLTREE = 3,     ///< SklTreeオブジェクトが不正(未初期化等)
  CL_SIZE_EXCEED = 4,     ///< ista[]+nlen[]が配列サイズを越えている
  CL_OTHER_ERROR = 10,    ///< その他のエラー
};


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
                         CutPosArray* cutPos, CutBidArray* cutBid);


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
                         CutPosArray* cutPos, CutBidArray* cutBid);


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
                         CutPosArray* cutPos, CutBidArray* cutBid);


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
                         CutPosArray* cutPos, CutBidArray* cutBid);

#ifdef CUTLIB_OCTREE

/// 交点情報計算: Octree, 計算対象セルタイプ指定
/**
 * @param[in,out] tree SklTreeクラスオブジェクト
 * @param[in] pl Polylibクラスオブジェクト
 * @param cutPos 交点座標データアクセッサ
 * @param cutBid 境界IDデータアクセッサ
 * @param[in] leafCellOnly 計算対象セルタイプフラグ
 */
CutlibReturn CutInfoOctree(SklTree* tree,
                           const Polylib* pl,
                           CutPosOctree* cutPos, CutBidOctree* cutBid,
                           bool leafCellOnly = true);


/// 交点情報計算: Octree, リーフセルのみ
/**
 * @param[in,out] tree SklTreeクラスオブジェクト
 * @param[in] pl Polylibクラスオブジェクト
 * @param cutPos 交点座標データアクセッサ
 * @param cutBid 境界IDデータアクセッサ
 */
CutlibReturn CutInfoOctreeLeafCell(SklTree* tree,
                                   const Polylib* pl,
                                   CutPosOctree* cutPos, CutBidOctree* cutBid);


/// 交点情報計算: Octree, 全セル計算
/**
 * @param[in,out] tree SklTreeクラスオブジェクト
 * @param[in] pl Polylibクラスオブジェクト
 * @param cutPos 交点座標データアクセッサ
 * @param cutBid 境界IDデータアクセッサ
 */
CutlibReturn CutInfoOctreeAllCell(SklTree* tree,
                          const Polylib* pl,
                          CutPosOctree* cutPos, CutBidOctree* cutBid);


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
                                   CutPosOctree* cutPos, CutBidOctree* cutBid);

#endif //CUTLIB_OCTREE


/// 交点情報計算: セル中心間, 一次元配列, 計算領域指定
/**
 * @param[in] ndim 全領域セル数
 * @param[in] ista 計算対象領域開始セル位置
 * @param[in] nlen 計算対象領域セル数
 * @param[in] org 領域原点座標
 * @param[in] d セル間隔
 * @param[in] pl Polylibクラスオブジェクト
 * @param[out] cutPos 交点座標配列
 * @param[out] cutBid 境界ID配列
 */
template<typename CUT_POS, typename CUT_BID>
CutlibReturn CutInfoCell(const size_t ndim[], const size_t ista[], const size_t nlen[],
                         const float org[], const float d[],
                         const Polylib* pl,
                         CUT_POS cutPos[], CUT_BID cutBid[])
{
  CutPosArrayTemplate<CUT_POS>* cutPosArray
      =  new CutPosArrayTemplate<CUT_POS>(cutPos, ndim[0], ndim[1], ndim[2]);
  CutBidArrayTemplate<CUT_BID>* cutBidArray
      =  new CutBidArrayTemplate<CUT_BID>(cutBid, ndim[0], ndim[1], ndim[2]);

  CutlibReturn ret = CutInfoCell(ista, nlen, org, d, pl, cutPosArray, cutBidArray);

  delete cutPosArray;
  delete cutBidArray;

  return ret;
}


/// 交点情報計算: セル中心間, 一次元配列, 全領域
/**
 * @param[in] ndim 全領域セル数
 * @param[in] org 領域原点座標
 * @param[in] d セル間隔
 * @param[in] pl Polylibクラスオブジェクト
 * @param[out] cutPos 交点座標配列
 * @param[out] cutBid 境界ID配列
 */
template<typename CUT_POS, typename CUT_BID>
CutlibReturn CutInfoCell(const size_t ndim[], const float org[], const float d[],
                         const Polylib* pl,
                         CUT_POS cutPos[], CUT_BID cutBid[])
{
  CutPosArrayTemplate<CUT_POS>* cutPosArray
      =  new CutPosArrayTemplate<CUT_POS>(cutPos, ndim[0], ndim[1], ndim[2]);
  CutBidArrayTemplate<CUT_BID>* cutBidArray
      =  new CutBidArrayTemplate<CUT_BID>(cutBid, ndim[0], ndim[1], ndim[2]);

  CutlibReturn ret = CutInfoCell(org, d, pl, cutPosArray, cutBidArray);

  delete cutPosArray;
  delete cutBidArray;

  return ret;
}


/// 交点情報計算: ノード間, 一次元配列, 計算領域指定
/**
 * @param[in] ndim 全領域セル数
 * @param[in] ista 計算対象領域開始ノード位置
 * @param[in] nlen 計算対象領域ノード数
 * @param[in] org 領域原点座標
 * @param[in] d セル間隔
 * @param[in] pl Polylibクラスオブジェクト
 * @param[out] cutPos 交点座標配列
 * @param[out] cutBid 境界ID配列
 */
template<typename CUT_POS, typename CUT_BID>
CutlibReturn CutInfoNode(const size_t ndim[], const size_t ista[], const size_t nlen[],
                         const float org[], const float d[],
                         const Polylib* pl,
                         CUT_POS cutPos[], CUT_BID cutBid[])
{
  CutPosArrayTemplate<CUT_POS>* cutPosArray
      =  new CutPosArrayTemplate<CUT_POS>(cutPos, ndim[0]+1, ndim[1]+1, ndim[2]+1);
  CutBidArrayTemplate<CUT_BID>* cutBidArray
      =  new CutBidArrayTemplate<CUT_BID>(cutBid, ndim[0]+1, ndim[1]+1, ndim[2]+1);

  CutlibReturn ret = CutInfoNode(ista, nlen, org, d, pl, cutPosArray, cutBidArray);

  delete cutPosArray;
  delete cutBidArray;

  return ret;
}


/// 交点情報計算: ノード間, 一次元配列, 全領域
/**
 * @param[in] ndim 全領域セル数
 * @param[in] org 領域原点座標
 * @param[in] d セル間隔
 * @param[in] pl Polylibクラスオブジェクト
 * @param[out] cutPos 交点座標配列
 * @param[out] cutBid 境界ID配列
 */
template<typename CUT_POS, typename CUT_BID>
CutlibReturn CutInfoNode(const size_t ndim[], const float org[], const float d[],
                         const Polylib* pl,
                         CUT_POS cutPos[], CUT_BID cutBid[])
{
  CutPosArrayTemplate<CUT_POS>* cutPosArray
      =  new CutPosArrayTemplate<CUT_POS>(cutPos, ndim[0]+1, ndim[1]+1, ndim[2]+1);
  CutBidArrayTemplate<CUT_BID>* cutBidArray
      =  new CutBidArrayTemplate<CUT_BID>(cutBid, ndim[0]+1, ndim[1]+1, ndim[2]+1);

  CutlibReturn ret = CutInfoNode(org, d, pl, cutPosArray, cutBidArray);

  delete cutPosArray;
  delete cutBidArray;

  return ret;
}

/** @} */ // end group CutInfoFunc

} /* namespace cutlib */

#endif /* CUTLIB_H */
