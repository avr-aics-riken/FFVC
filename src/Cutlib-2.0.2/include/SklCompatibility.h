/**@file
 * @brief pFTTクラスライブラリとの互換性のためのヘッダ
 *
 * Compatibility.hの内容を「Sph〜」を「Skl〜」に修正して作成
 */

#ifndef SKLCOMPATIBILITY_H_
#define SKLCOMPATIBILITY_H_

#include "Tree.h"
#include "Node.h"
#include "Cell.h"

#include "DimPolicy.h"
#include "DataPolicy.h"
#include "Timing.h"

typedef TD_NAMESPACE::Tree< 
    TD_NAMESPACE::ThreeDimPolicy, 
    TD_NAMESPACE::DefaultDataPolicy< float > 
> SklTree;
typedef SklTree::cell_type SklCell;

typedef SklCell::NeighborSet SklNghbrSet;
typedef SklTree::SplitList SklSplitList;

inline
double
SklGetTime()
{
    return TD_NAMESPACE::GetTime();
}

#endif /*SKLCOMPATIBILITY_H_*/
