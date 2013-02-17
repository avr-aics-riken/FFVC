#ifndef _IP_PMT_H_
#define _IP_PMT_H_

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
 * @file   IP_PMT.h
 * @brief  IP_PMT class Header
 * @author kero
 */

#include "Intrinsic.h"
#include "IP_Define.h"

class IP_PMT : public Intrinsic {
public:
  /** コンストラクタ */
  IP_PMT(){}
  
  /**　デストラクタ */
  ~IP_PMT() {}

public:
  std::string m_fluid; ///< 流体のラベル
  std::string m_solid; ///< 固体のラベル
  
protected:

public:
  // パラメータを取得する
  virtual bool getTP(Control* R, TPControl* tpCntl);
  
  
  // 領域情報を設定する
  virtual void setDomain(Control* R, const int* sz, REAL_TYPE* org, REAL_TYPE* reg, REAL_TYPE* pch);
  
  
  // 計算領域のセルIDを設定する
  virtual void setup(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, MediumList* mat);
  
  
  /** 
   * @brief 例題の名称を返す
   */
  virtual const char* getExampleName() 
  {
    return ("Performance Test");
  }
  
};
#endif // _IP_PMT_H_
