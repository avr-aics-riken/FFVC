#ifndef _IP_FREEJET_H_
#define _IP_FREEJET_H_

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
  virtual bool printPara(FILE* fp, Control* R);
  
  
  /** 矩形の計算領域のセルIDを設定する
   * @param[in/out] mid   媒質情報の配列
   * @param[in]     R     Controlクラスのポインタ
   * @param[in]     G_org グローバルな原点（無次元）
   * @param[in]     Nmax  Controlクラスのポインタ
   * @param[in]     mat   MediumListクラスのポインタ
   */
  virtual void setup(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, MediumList* mat);
  
public:
  virtual bool getParaXML(SklSolverConfig* CF, Control* R);
  virtual bool initVars(Control* R);
  
  /** 領域を設定する
   * @param[in] R   Controlクラスのポインタ
   * @param[in] sz  分割数
   * @param[in] org 計算領域の基点
   * @param[in] reg 計算領域のbounding boxサイズ
   * @param[in] pch セル幅
   */
  virtual bool setDomain(Control* R, const int* sz, const REAL_TYPE* org, const REAL_TYPE* reg, const REAL_TYPE* pch);
  
  
  /** 
   @brief 例題の名称を返す
   */
  virtual const char* getExampleName(void) {
    return ("FreeJet");
  }
  
  virtual void PostInit(REAL_TYPE &checkTime, Control* R);
};

#endif // _IP_FREEJET_H_
