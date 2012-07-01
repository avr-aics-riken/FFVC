#ifndef _FB_INTRNSC_H_
#define _FB_INTRNSC_H_

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
 * @file   Intrinsic.h
 * @brief  FlowBase Intrinsic class Header
 * @author kero
 */

#include <math.h>
#include <fstream>

#include "cpm_Define.h"
#include "cpm_ParaManager.h"

#include "DomainInfo.h"
#include "FB_Define.h"
#include "BndOuter.h"
#include "Control.h"
#include "Component.h"
#include "vec3.h"
#include "TPControl.h"
#include "Medium.h"


/** 組み込み例題のID */
enum Intrinsic_class 
{
  id_Duct,
  id_PPLT2D,
  id_SHC1D,
  id_PMT,
  id_Rect,
  id_Cylinder,
  id_Step,
  id_Polygon,
  id_Sphere
};


class Intrinsic : public DomainInfo {
  
public:
  REAL_TYPE RefL; ///< 代表長さ
  
  /** 次元のモード */
  enum dim_mode 
  {
    dim_2d = 1,
    dim_3d
  };
  
  /** コンストラクタ */
  Intrinsic() 
  { 
    RefL = 0.0;
  }
  
  /**　デストラクタ */
  virtual ~Intrinsic() {}
    
  
public:
  
  /** 
   @brief ユーザー例題の名称を返す
   */
  virtual const char* getExampleName() 
  { 
    return NULL; 
  }
  
  /**
   * @param [in] R       Controlクラスのポインタ
   * @param [in] tpCntl  TPControlクラスのポインタ
   * @return true-成功、false-エラー
   */
  virtual bool getTP(Control* R, TPControl* tpCntl) 
  { 
    return true; 
  }
  
  virtual void initCond(REAL_TYPE* v, REAL_TYPE* p) {};
  
  virtual void PostInit(REAL_TYPE &checkTime, Control* R) {};
  
  
  /**
   @brief 例題名称の表示
   @param [in] fp   出力ファイルのファイルポインタ
   @param [in] str  表示文字列
   */
  virtual void printExample(FILE* fp, const char* str);
  
  
  /**
   @brief パラメータの表示
   @param [in] fp ファイルポインタ
   @param [in] R  コントロールクラスのポインタ
   */
  virtual void printPara(FILE* fp, const Control* R);
  
  
  /** 領域を設定する
   * @param [in]     R   Controlクラスのポインタ
   * @param [in]     sz  分割数
   * @param [in/out] org 計算領域の基点
   * @param [in/out] reg 計算領域の大きさ
   * @param [in/out] pch セル幅
   */
  virtual void setDomain(Control* R, const int* sz, REAL_TYPE* org, REAL_TYPE* reg, REAL_TYPE* pch) {};
  
  
  /** 計算領域の媒質情報を設定する
   * @param [in/out] mid   媒質情報の配列
   * @param [in]     R     Controlクラスのポインタ
   * @param [in]     G_org Controlクラスのポインタ
   * @param [in]     Nmax  Controlクラスのポインタ
   * @param [in]     mat   MediumListクラスのポインタ
   */
  virtual void setup(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, MediumList* mat) {};
  
  
  
  virtual void setup_cut(int* mid, Control* R, REAL_TYPE* G_org, const int Nmax, MediumList* mat, float* cut) {};
  
  
  
  /**
   * @brief 例題のモデルをsvxフォーマットで出力する(体積率とID)
   * @param [in] vf 体積占有率
   * @param [in] id ID情報
   * @param [in] R  コントロールクラスのポインタ
   */
  void writeSVX(REAL_TYPE *vf, int *id, Control* R);
  
  
  /**
   * @brief 例題のモデルをsvxフォーマットで出力する(ID)
   * @param [in] id ID情報
   * @param [in] R  コントロールクラスのポインタ
   */
  void writeSVX(int *id, Control* R);
  
  
  /** MediumList中に登録されているkeyに対するIDを返す
   * @param [in] mat  MediumListクラス
   * @param [in] Namx リストの最大数
   * @param [in] key  探査するラベル
   * @return keyに対するIDを返す。発見できない場合はzero
   */
  int find_ID_from_Label(MediumList* mat, const int Nmax, const std::string key);
  
};

#endif // _FB_INTRNSC_H_
