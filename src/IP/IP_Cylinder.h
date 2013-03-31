#ifndef _IP_CYL_H_
#define _IP_CYL_H_

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
 * @file   IP_Cylinder.h
 * @brief  IP_Cylinder class Header
 * @author kero
 */

#include "Intrinsic.h"
#include "IP_Define.h"

class IP_Cylinder : public Intrinsic {
protected:
  int mode;                  ///< モード（2D or 3D）
  REAL_TYPE width;           ///< 流路の幅
  REAL_TYPE height;          ///< ドライバ部の高さ
  REAL_TYPE drv_length;      ///< ドライバの長さ
  int drv_mode;              ///< ドライバのON/OFF
  std::string m_fluid;       ///< 流体のラベル
  std::string m_solid;       ///< 固体のラベル
  std::string m_driver;      ///< ドライバ部分のラベル
  std::string m_driver_face; ///< ドライバ指定面のラベル
  
public:
  /** コンストラクタ */
  IP_Cylinder(){
    mode   = 0;
    width  = 0.0;
    height = 0.0;
    drv_length = 0.0;
  }
  
  /**　デストラクタ */
  ~IP_Cylinder() {}

public:
  // パラメータを取得する
  virtual bool getTP(Control* R, TPControl* tpCntl);
  
  
  // 領域パラメータを設定する
  virtual void setDomainParameter(Control* R, const int* sz, REAL_TYPE* org, REAL_TYPE* reg, REAL_TYPE* pch);
  
  
  // Cylinderの計算領域のセルIDを設定する
  virtual void setup(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, MediumList* mat);
  
  
  // パラメータの表示
  virtual void printPara(FILE* fp, const Control* R);
  
  
  /** 
   * @brief 例題の名称を返す
   */
  virtual const char* getExampleName()
  {
    return ("Back Step");
  }
  
};
#endif // _IP_CYL_H_
