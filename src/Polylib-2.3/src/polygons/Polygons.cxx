/*
 * Polylib - Polygon Management Library.
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2010-
 *
 */
#include <iostream>
#include "common/PolylibCommon.h"
#include "polygons/Polygons.h"

namespace PolylibNS {

using namespace std;

/************************************************************************
 *
 * Polygonsクラス
 *  @attension 三角形ポリゴン集合を管理する純粋仮想クラス
 *
 ***********************************************************************/
///
/// デストラクタ
///  @attention 継承しているクラスから呼び出されるために必要。
///
Polygons::~Polygons() {}

} //namespace PolylibNS
