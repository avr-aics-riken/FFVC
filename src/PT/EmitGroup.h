#ifndef _PT_GROUP_H_
#define _PT_GROUP_H_

//##################################################################################
//
// Copyright (c) 2019 Research Institute for Information Technology, Kyushu university
// All rights researved.
//
//##################################################################################

/**
 * @file   EmitGroup.h
 * @brief  Cloud class Header
 * @author riit
 * @note
 */

#include "PtDefine.h"
#include <string.h>


/**
 * @brief 粒子管理グループクラス
 * @note
 */
class EmitGroup {
public:
  // Line
  int nDiv;            ///< 分割数
  REAL_TYPE from[3];   ///< 線分
  REAL_TYPE to[3];     ///< 線分

  // Disc
  int nSample;         ///< サンプル数
  REAL_TYPE center[3]; ///< 円の中心座標
  REAL_TYPE normal[3]; ///< 法線ベクトル
  REAL_TYPE radius;    ///< 円の半径

protected:
  int type;            ///< 開始点タイプ (Pointset, Line, Disc,..)
  std::string name;    ///< グループ名


public:
  EmitGroup() {
    type = 0;
    nDiv = 0;
    nSample = 0;
    radius = 0.0;
    for (int i=0; i<3; i++) {
      from[i] = 0.0;
      to[i] = 0.0;
      center[i]=0.0;
      normal[i]=0.0;
    }
  }

  ~EmitGroup() {}

  
  void setGrpName(const std::string mm) {
    name = mm;
  }


  std::string getGrpName() {
    return name;
  }
  
  
  void setType(const int m_type) {
    type = m_type;
  }


  int getType() const {
    return type;
  }

};

#endif // _PT_GROUP_H_
