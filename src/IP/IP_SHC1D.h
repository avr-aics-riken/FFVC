#ifndef _IP_SHC1D_H_
#define _IP_SHC1D_H_

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
 @file   IP_SHC1D.h
 @brief  IP_SHC1D class Header
 @author kero
 */

#include "Intrinsic.h"
#include "IP_Define.h"

class IP_SHC1D : public Intrinsic {
  
public:   
  std::string m_fluid;      ///< 流体のラベル
  std::string m_solid;      ///< 固体のラベル
  std::string m_inactive;   ///< 固体で不活性セルのラベル
  
public:
  /** コンストラクタ */
  IP_SHC1D(){}
  
  /**　デストラクタ */
  ~IP_SHC1D() {}

public:
  // パラメータを取得する
  virtual bool getTP(Control* R, TPControl* tpCntl);
  
  
  // 領域情報を設定する
  virtual void setDomain(Control* R, const int* sz, REAL_TYPE* org, REAL_TYPE* reg, REAL_TYPE* pch);
  
  
  // モデルIDのセットアップ
  virtual void setup(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, MediumList* mat);
  
  
  /** 
   @brief 例題の名称を返す
   */
  virtual const char* getExampleName() 
  {
    return ("Steady 1D Heat Conduction");
  }
  
  
  // 計算領域のセルIDとカット情報を設定する
  void setup_bc(int* bid);
  
};

#endif // _IP_SHC1D_H_
