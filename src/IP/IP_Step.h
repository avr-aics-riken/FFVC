#ifndef _IP_STEP_H_
#define _IP_STEP_H_

//##################################################################################
//
// FFV-C : Frontflow / violet Cartesian
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2015 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2016 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

/** 
 * @file   IP_Step.h
 * @brief  IP_Step class Header
 * @author aics
 */

#include "Intrinsic.h"
#include "IP_Define.h"
#include "FBUtility.h"


class IP_Step : public Intrinsic {
  
protected:
  REAL_TYPE width;           ///< 流路の幅
  REAL_TYPE height;          ///< ステップ部の高さ
  
public:
  /** コンストラクタ */
  IP_Step(){
    width  = 0.0;
    height = 0.0;
  }
  
  /**　デストラクタ */
  ~IP_Step() {}

public:

  /**
   * @brief パラメータをロード
   * @param [in] R      Controlクラス
   * @param [in] tpCntl テキストパーサクラス
   * @return true-成功, false-エラー
   */
  virtual bool getTP(Control* R, TextParser* tpCntl);
  
  /**
   * @brief パラメータの表示
   * @param [in] fp ファイルポインタ
   * @param [in] R  コントロールクラスのポインタ
   */
  virtual void printPara(FILE* fp, const Control* R);
  
  
  /**
   * @brief 計算領域のセルIDを設定する
   * @param [in,out] bcd      BCindex B
   * @param [in]     R        Controlクラスのポインタ
   * @param [in]     NoMedium 媒質数
   * @param [in]     mat      MediumListクラスのポインタ
   * @param [in]     NoCompo  コンポーネント数
   * @param [in]     cmp      CompoListクラスのポインタ
   * @param [out]    cut      カット情報
   * @param [out]    bid      境界ID
   */
  virtual void setup(int* bcd,
                     Control* R,
                     const int NoMedium,
                     const MediumList* mat,
                     const int NoCompo,
                     const CompoList* cmp,
                     int* cutL,
                     int* cutU,
                     int* bid);

};
#endif // _IP_STEP_H_
