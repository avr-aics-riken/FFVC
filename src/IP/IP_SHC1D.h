#ifndef _IP_SHC1D_H_
#define _IP_SHC1D_H_

// #################################################################
//
// CAERU Library
//
// Copyright (c) All right reserved. 2012
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
//
// #################################################################

/**
 @file IP_SHC1D.h
 @brief IP_SHC1D class Header
 @author kero
 */

#include "Intrinsic.h"
#include "IP_Define.h"

class IP_SHC1D : public Intrinsic {
  
public:
  std::string m_inactive;
  std::string m_fin;
  std::string m_isothermal;
  std::string m_adiabatic;
  std::string m_fluid;      ///< 流体のラベル
  std::string m_solid;      ///< 固体のラベル
  
public:
  /** コンストラクタ */
  IP_SHC1D(){}
  
  /**　デストラクタ */
  ~IP_SHC1D() {}
  
protected:

public:
  /** パラメータをロード
   * @apram[in] R      Controlクラス
   * @param[in] tpCntl テキストパーサクラス
   * @return true-成功, false-エラー
   */
  virtual bool getTP(Control* R, TPControl* tpCntl);
  
  
  /** 領域を設定する
   * @param[in] R   Controlクラスのポインタ
   * @param[in] sz  分割数
   * @param[in] org 計算領域の基点
   * @param[in] wth 計算領域のbounding boxサイズ
   * @param[in] pch セル幅
   */
  virtual void setDomain(Control* R, unsigned sz[3], REAL_TYPE org[3], REAL_TYPE wth[3], REAL_TYPE pch[3]);
  
  
  /** 矩形の計算領域のセルIDを設定する
   * @param[in/out] mid   媒質情報の配列
   * @param[in]     R     Controlクラスのポインタ
   * @param[in]     G_org グローバルな原点（無次元）
   * @param[in]     Nmax  Controlクラスのポインタ
   * @param[in]     mat   MediumListクラスのポインタ
   */
  virtual void setup(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, MediumList* mat);
  
  
  /** 
   @brief 例題の名称を返す
   */
  virtual const char* getExampleName(void) {
    return ("Steady 1D Heat Conduction");
  }
};

#endif // _IP_SHC1D_H_
