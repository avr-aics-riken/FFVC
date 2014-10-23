#ifndef _FB_INTRNSC_H_
#define _FB_INTRNSC_H_

//##################################################################################
//
// Flow Base class
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2014 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

/** 
 * @file   Intrinsic.h
 * @brief  FlowBase Intrinsic class Header
 * @author aics
 */

#include "cpm_ParaManager.h" // 最初に必要
#include <math.h>
#include <fstream>
#include "DomainInfo.h"
#include "FB_Define.h"
#include "BndOuter.h"
#include "Control.h"
#include "Component.h"
#include "Medium.h"
#include "TextParser.h"


class Intrinsic : public DomainInfo {
  
public:
  
  REAL_TYPE RefL;      ///< 代表長さ [m]
  REAL_TYPE RefV;      ///< 代表速度 [m/s]
  
  int even;            ///< 偶数分割のチェック
  int mode;            ///< 次元数
  
  std::string m_fluid; ///< 流体のラベル
  std::string m_solid; ///< 固体のラベル
  
  /** 次元のモード */
  enum dim_mode 
  {
    dim_2d = 1,
    dim_3d
  };
  
  /** 形状 */
  enum shape_type {
    id_circular = 1,
    id_rectangular
  };
  
  /** コンストラクタ */
  Intrinsic() 
  { 
    RefL = 0.0;
    RefV = 0.0;
    even = 0;
    mode = dim_3d;
  }
  
  /**　デストラクタ */
  virtual ~Intrinsic() {}
    
  
public:
  
  /**
   * @brief 例題クラス固有のパラメータをロードする
   * @param [in] R      Controlクラス
   * @param [in] tpCntl テキストパーサクラス
   * @return true-成功, false-エラー
   */
  virtual bool getTP(Control* R, TextParser* tpCntl) 
  { 
    return true; 
  }
  
  
  /**
   * @brief 例題名称の表示
   * @param [in] fp   出力ファイルのファイルポインタ
   * @param [in] m_id 例題ID
   */
  virtual void printExample(FILE* fp, const int m_id);
  
  
  
  /**
   * @brief パラメータの表示
   * @param [in] fp ファイルポインタ
   * @param [in] R  コントロールクラスのポインタ
   */
  virtual void printPara(FILE* fp, const Control* R);
  
  
  /**
   * @brief 領域パラメータを設定する
   * @param [in]     R     Controlクラスのポインタ
   * @param [in]     sz    分割数
   * @param [in,out] m_org 計算領域の基点
   * @param [in,out] m_reg 計算領域のbounding boxサイズ
   * @param [in,out] m_pch セル幅
   */
  virtual void setDomainParameter(Control* R,
                                  const int* sz,
                                  REAL_TYPE* org,
                                  REAL_TYPE* reg,
                                  REAL_TYPE* pch) {};
  
  
  /**
   * @brief 外部境界の設定
   * @param [in]     face     面番号
   * @param [in,out] bcd      BCindexx B
   * @param [in]     R        Controlクラスのポインタ
   * @param [in]     G_org    グローバルな原点（無次元）
   * @param [in]     NoMedium 媒質数
   * @param [in]     mat      MediumListクラスのポインタ
   * @param [out]    cut      カット情報
   * @param [out]    bid      境界ID
   */
  virtual void setOBC(const int face,
                      int* bcd,
                      Control* R,
                      REAL_TYPE* G_org,
                      const int NoMedium,
                      const MediumList* mat,
                      float* cut,
                      int* bid) {};
  
  
  
  /**
   * @brief 代表パラメータの設定
   * @param [in] Cref  コントロールクラスのポインタ
   */
  void setRefParameter(Control* Cref);
  
  
  /**
   * @brief 計算領域のセルIDとカット情報を設定する
   * @param [in,out] bcd      BCindex B
   * @param [in]     R        Controlクラスのポインタ
   * @param [in]     NoMedium 媒質数
   * @param [in]     mat      MediumListクラスのポインタ
   * @param [in]     cut      交点情報
   * @param [in]     bid      境界ID
   */
  virtual void setup(int* bcd,
                     Control* R,
                     const int NoMedium,
                     const MediumList* mat,
                     float* cut,
                     int* bid) {};
  
  
  
  /**
   * @brief モデルIDをsphフォーマット(float)で出力する
   * @param [in] bcd BCindex B
   * @param [in] R   コントロールクラスのポインタ
   */
  void writeSPH(const int* bcd, const Control* R);
  
  
  
  /**
   * @brief 例題のモデルをsvxフォーマットで出力する(体積率とID)
   * @param [in] vf 体積占有率
   * @param [in] id ID情報
   * @param [in] R  コントロールクラスのポインタ
   */
  void writeSVX(REAL_TYPE *vf, int *id, Control* R);
  
  
  /**
   * @brief 例題のモデルをsvxフォーマットで出力する(ID)
   * @param [in] bcd BCindex B
   * @param [in] R   コントロールクラスのポインタ
   */
  void writeSVX(const int* bcd, Control* R);
  
};

#endif // _FB_INTRNSC_H_
