#ifndef _FB_INTRNSC_H_
#define _FB_INTRNSC_H_

//##################################################################################
//
// Flow Base class
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2013 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2013 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

/** 
 * @file   Intrinsic.h
 * @brief  FlowBase Intrinsic class Header
 * @author kero
 */

#include <math.h>
#include <fstream>

#include "cpm_ParaManager.h"

#include "DomainInfo.h"
#include "FB_Define.h"
#include "BndOuter.h"
#include "Control.h"
#include "Component.h"
#include "vec3.h"
#include "Medium.h"
#include "TextParser.h"


class Intrinsic : public DomainInfo {
  
public:
  
  REAL_TYPE RefL;      ///< 代表長さ [m]
  REAL_TYPE RefV;      ///< 代表速度 [m/s]
  
  int even;            ///< 偶数分割のチェック
  int mode;            ///< 次元数
  
  std::string m_fluid; ///< 流体のラベル
  std::string m_solid; ///< 固体のラベル
  
  /** 次元のモード */
  enum dim_mode 
  {
    dim_2d = 1,
    dim_3d
  };
  
  /** コンストラクタ */
  Intrinsic() 
  { 
    RefL = 0.0;
    RefV = 0.0;
    even = 0;
    mode = dim_3d;
  }
  
  /**　デストラクタ */
  virtual ~Intrinsic() {}
    
  
public:
  
  // 例題クラス固有のパラメータをロードする
  virtual bool getTP(Control* R, TextParser* tpCntl) 
  { 
    return true; 
  }
  
  
  // 例題名称の表示
  virtual void printExample(FILE* fp, const int m_id);
  
  
  // 例題クラス固有のパラメータの表示
  virtual void printPara(FILE* fp, const Control* R);
  
  
  // 領域パラメータを設定する
  virtual void setDomainParameter(Control* R, const int* sz, REAL_TYPE* org, REAL_TYPE* reg, REAL_TYPE* pch) {};
  
  
  // 計算領域の媒質情報を設定する
  virtual void setup(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, const MediumList* mat) {};
  
  
  
  
  // ユーザー例題の名称を返す
  const char* getExampleName(int m_id);
  
  
  // 代表パラメータのセット
  void setRefParameter(Control* Cref);
  
  
  // モデルIDをsphフォーマット(float)で出力する
  void writeSPH(const int *mid, const Control* R);
  
  
  // 例題のモデルをsvxフォーマットで出力する(体積率とID)
  void writeSVX(REAL_TYPE *vf, int *id, Control* R);
  
  
  // 例題のモデルをsvxフォーマットで出力する(ID)
  void writeSVX(int *id, Control* R);
  
};

#endif // _FB_INTRNSC_H_
