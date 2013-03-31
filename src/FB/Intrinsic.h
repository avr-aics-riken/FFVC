#ifndef _FB_INTRNSC_H_
#define _FB_INTRNSC_H_

// #################################################################
//
// CAERU Library
//
// Copyright (c) 2012-2013  All right reserved.
//
// Institute of Industrial Science, University of Tokyo, Japan. 
//
// #################################################################

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
#include "TPControl.h"
#include "Medium.h"


/** 組み込み例題のID */
enum Intrinsic_class 
{
  id_Duct,
  id_PPLT2D,
  id_SHC1D,
  id_PMT,
  id_Rect,
  id_Cylinder,
  id_Step,
  id_Polygon,
  id_Sphere,
  id_Jet
};


class Intrinsic : public DomainInfo {
  
public:
  REAL_TYPE RefL; ///< 代表長さ
  
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
  }
  
  /**　デストラクタ */
  virtual ~Intrinsic() {}
    
  
public:
  
  // ユーザー例題の名称を返す
  virtual const char* getExampleName() 
  { 
    return NULL; 
  }
  
  // パラメータをロードする
  virtual bool getTP(Control* R, TPControl* tpCntl) 
  { 
    return true; 
  }
  
  
  // 
  virtual void initCond(REAL_TYPE* v, REAL_TYPE* p) {};
  
  
  // 例題名称の表示
  virtual void printExample(FILE* fp, const char* str);
  
  
  // パラメータの表示
  virtual void printPara(FILE* fp, const Control* R);
  
  
  // 領域パラメータを設定する
  virtual void setDomainParameter(Control* R, const int* sz, REAL_TYPE* org, REAL_TYPE* reg, REAL_TYPE* pch) {};
  
  
  // 計算領域の媒質情報を設定する
  virtual void setup(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, MediumList* mat) {};
  
  
  // モデルIDをsphフォーマット(float)で出力する
  void writeSPH(const int *mid, const Control* R);
  
  
  // 例題のモデルをsvxフォーマットで出力する(体積率とID)
  void writeSVX(REAL_TYPE *vf, int *id, Control* R);
  
  
  // 例題のモデルをsvxフォーマットで出力する(ID)
  void writeSVX(int *id, Control* R);
  
};

#endif // _FB_INTRNSC_H_
