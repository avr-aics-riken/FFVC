#ifndef _IP_STEP_H_
#define _IP_STEP_H_

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
 * @file   IP_Step.h
 * @brief  IP_Step class Header
 * @author kero
 */

#include "Intrinsic.h"
#include "IP_Define.h"

class IP_Step : public Intrinsic {
  
protected:
  REAL_TYPE width;           ///< 流路の幅
  REAL_TYPE height;          ///< ドライバ部の高さ
  REAL_TYPE drv_length;      ///< ドライバの長さ
  std::string m_driver;      ///< ドライバ部分のラベル
  std::string m_driver_face; ///< ドライバ指定面のラベル
  
public:
  /** コンストラクタ */
  IP_Step(){
    width  = 0.0;
    height = 0.0;
    drv_length = 0.0;
  }
  
  /**　デストラクタ */
  ~IP_Step() {}

public:

  virtual bool getTP(Control* R, TPControl* tpCntl);
  
  virtual void printPara(FILE* fp, const Control* R);
  
  virtual void setup(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, MediumList* mat);

};
#endif // _IP_STEP_H_
