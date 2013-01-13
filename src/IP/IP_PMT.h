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
  /** パラメータをロード
   * @param [in] R      Controlクラス
   * @param [in] tpCntl テキストパーサクラス
   * @return true-成功, false-エラー
   */
  virtual bool getTP(Control* R, TPControl* tpCntl);
  
  
  /** 領域を設定する
   * @param [in]     R   Controlクラスのポインタ
   * @param [in]     sz  分割数
   * @param [in,out] org 計算領域の基点
   * @param [in,out] reg 計算領域のbounding boxサイズ
   * @param [in,out] pch セル幅
   */
  virtual void setDomain(Control* R, const int* sz, REAL_TYPE* org, REAL_TYPE* reg, REAL_TYPE* pch);
  
  
  /** Cavityの計算領域のセルIDを設定する
   * @param [in,out] mid   媒質情報の配列
   * @param [in]     R     Controlクラスのポインタ
   * @param [in]     G_org グローバルな原点（無次元）
   * @param [in]     Nmax  Controlクラスのポインタ
   * @param [in]     mat   MediumListクラスのポインタ
   */
  virtual void setup(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, MediumList* mat);
  
  
  /** 
   * @brief 例題の名称を返す
   */
  virtual const char* getExampleName(void) 
  {
    return ("Performance Test");
  }
  
};
#endif // _IP_PMT_H_
