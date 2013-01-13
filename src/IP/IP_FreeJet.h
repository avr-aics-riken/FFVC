#ifndef _IP_FREEJET_H_
#define _IP_FREEJET_H_

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
 @file IP_FreeJet.h
 @brief IP_FreeJet class Header
 @author keno, FSI Team, RIKEN
 */

#include "Intrinsic.h"
#include "IP_Define.h"

class FreeJet : public Intrinsic {
public:
  /** コンストラクタ */
  FreeJet(){}
  
  /**　デストラクタ */
  ~FreeJet() {}
  
protected:
  virtual bool printPara(FILE* fp, const Control* R);
  
  
  /** 矩形の計算領域のセルIDを設定する
   * @param [in,out] mid   媒質情報の配列
   * @param [in]     R     Controlクラスのポインタ
   * @param [in]     G_org グローバルな原点（無次元）
   * @param [in]     Nmax  Controlクラスのポインタ
   * @param [in]     mat   MediumListクラスのポインタ
   */
  virtual void setup(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, MediumList* mat);
  
public:
  virtual bool getParaXML(SklSolverConfig* CF, Control* R);
  virtual bool initVars(Control* R);
  
  /** 領域を設定する
   * @param [in]     R   Controlクラスのポインタ
   * @param [in]     sz  分割数
   * @param [in,out] org 計算領域の基点
   * @param [in,out] reg 計算領域のbounding boxサイズ
   * @param [in,out] pch セル幅
   */
  virtual bool setDomain(Control* R, const int* sz, REAL_TYPE* org, REAL_TYPE* reg, REAL_TYPE* pch);
  
  
  /** 
   @brief 例題の名称を返す
   */
  virtual const char* getExampleName(void) 
  {
    return ("FreeJet");
  }
  
};

#endif // _IP_FREEJET_H_
