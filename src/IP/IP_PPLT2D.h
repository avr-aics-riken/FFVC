#ifndef _IP_PPLT2D_H_
#define _IP_PPLT2D_H_

// #################################################################
//
// FFV : Frontflow / violet Cartesian
//
// Copyright (c) 2012-2013  All right reserved.
//
// Institute of Industrial Science, University of Tokyo, Japan. 
//
// #################################################################

/** 
 * @file   IP_PPLT2D.h
 * @brief  IP_PPLT2D class Header
 * @author kero
 */

#include "Intrinsic.h"
#include "IP_Define.h"

class IP_PPLT2D : public Intrinsic {
public:
  std::string m_fluid; ///< 流体のラベル
  std::string m_solid; ///< 固体のラベル
  
public:
  /** コンストラクタ */
  IP_PPLT2D(){}
  
  /**　デストラクタ */
  ~IP_PPLT2D() {}
  
protected:

public:
  
  // パラメータを取得する
  virtual bool getTP(Control* R, TPControl* tpCntl);
  
  
  // 領域パラメータを設定する
  virtual void setDomainParameter(Control* R, const int* sz, REAL_TYPE* org, REAL_TYPE* reg, REAL_TYPE* pch);
  
  
  // PPLT2Dの計算領域のセルIDを設定する
  virtual void setup(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, MediumList* mat);
  
  
  /** 
   * @brief 例題の名称を返す
   */
  virtual const char* getExampleName() 
  {
    return ("Parallel Plate 2D");
  }
};
#endif // _IP_PPLT2D_H_
