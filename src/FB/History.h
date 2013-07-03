#ifndef _FB_HISTORY_H_
#define _FB_HISTORY_H_

//##################################################################################
//
// Flow Base class
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2013 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2013 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

/** 
 * @file   History.h
 * @brief  FlowBase History class Header
 * @author kero
 */


#include "Control.h"
#include "Component.h"
#include "FBUtility.h"
#include "vec3.h"

class History {
protected:
  REAL_TYPE Tscale;          ///< 時間スケール
  REAL_TYPE RefVelocity;     ///< 代表速度
  REAL_TYPE BaseTemp;        ///< 基準温度
  REAL_TYPE DiffTemp;        ///< 代表温度差
  REAL_TYPE RefDensity;      ///< 代表密度
  REAL_TYPE RefLength;       ///< 代表長さ
  REAL_TYPE BasePrs;         ///< 基準圧力
  REAL_TYPE RefSpecificHeat; ///< 代表比熱
  REAL_TYPE time;            ///< 時刻
  REAL_TYPE v_max;           ///< 速度成分の最大値
  REAL_TYPE dhd;             ///< 有次元の格子幅
  REAL_TYPE dh;              ///< 無次元の格子幅
  REAL_TYPE rhocp;           ///< 熱流束計算の無次元化量
  REAL_TYPE dynamic_p;       ///< 動圧
  REAL_TYPE base_mf;         ///< 流量の基準値
  int step;                  ///< ステップ数
  int Unit_Temp;             ///< 温度単位
  int Unit_Prs;              ///< 圧力基準モード
  int Unit_Log;              ///< ログ出力の単位
  
public:
  
  /** デフォルトコンストラクタ */
  History() {}
  
  /** コンストラクタ */
  History(const Control* Cref) 
  {
    RefVelocity     = Cref->RefVelocity;
    BaseTemp        = Cref->BaseTemp;
    DiffTemp        = Cref->DiffTemp;
    RefDensity      = Cref->RefDensity;
    RefLength       = Cref->RefLength;
    BasePrs         = Cref->BasePrs;
    RefSpecificHeat = Cref->RefSpecificHeat;
    Unit_Temp       = Cref->Unit.Temp;
    Unit_Prs        = Cref->Unit.Prs;
    Unit_Log        = Cref->Unit.Log;
    dh              = Cref->deltaX;
    Tscale          = RefLength / RefVelocity;
    dhd             = dh*RefLength;
    rhocp           = RefVelocity * DiffTemp * RefDensity * RefSpecificHeat;
    dynamic_p       = RefVelocity * RefVelocity * RefDensity;
    base_mf         = RefVelocity * RefLength * RefLength;

    time  = 0.0;
    v_max = 0.0;
    step  = 0;
  }
  
  /**　デストラクタ */
  ~History() {}
  
  
protected:

  /** @brief モードに対応する時刻を返す */
  REAL_TYPE printTime() const
  {
    return ( (Unit_Log == DIMENSIONAL) ? time*Tscale : time );
  }
  

  /** @brief モードに対応する速度最大値を返す */
  REAL_TYPE printVmax() const
  {
    return ( (Unit_Log == DIMENSIONAL) ? v_max*RefVelocity : v_max );
  }
  

  /**
   * @brief モードに対応する長さを返す
   * @param [in] var  長さ
   */
  REAL_TYPE printLen(const REAL_TYPE var) const
  {
    return ( (Unit_Log == DIMENSIONAL) ? var*RefLength : var );
  }
  

  /** 
   * @brief モードに対応する速度値を返す
   * @param [in] var  速度
   */
  REAL_TYPE printVel(const REAL_TYPE var) const
  {
    return ( (Unit_Log == DIMENSIONAL) ? var*RefVelocity : var );
  }
  

  /**
   * @brief 圧力の次元変換
   * @param [in] var  圧力
   */
  REAL_TYPE printPrs(const REAL_TYPE var) const
  {
    if (Unit_Log != DIMENSIONAL) return var;
    const REAL_TYPE a = var * dynamic_p;
    return ( (Unit_Prs==Unit_Absolute) ? BasePrs+a : a );
  }
  

  /** 
   * @brief 全圧の次元変換
   * @param [in] var  圧力
   */
  REAL_TYPE printTP(const REAL_TYPE var) const
  {
    return ( (Unit_Log == DIMENSIONAL) ? var*dynamic_p : var );
  }
  

  /**
   * @brief 全圧の次元変換
   * @param [in] var  全圧
   */
  REAL_TYPE printForce(const REAL_TYPE var) const
  {
    return ( (Unit_Log == DIMENSIONAL) ? var*dynamic_p*RefLength*RefLength : var );
  }
  

  /** 
   * @brief モードに対応する流量を返す
   * @param [in] var  流量
   */
  REAL_TYPE printMF(const REAL_TYPE var) const
  {
    return ( (Unit_Log == DIMENSIONAL) ? var*base_mf : var );
  }
  
  
  /** 
   * @brief モードに対応する温度を返す
   * @param [in] var  温度
   */
  REAL_TYPE printTmp(const REAL_TYPE var) const
  {
    if (Unit_Log != DIMENSIONAL) return var;
    const REAL_TYPE a = BaseTemp + DiffTemp*var; // Kelvin
    return ( (Unit_Temp==Unit_KELVIN) ? a : a-KELVIN  );
  }
  
  

  /**
   * @brief モードに対応する熱量を返す(面要素)
   * @param [in] var  熱量
   */
  REAL_TYPE printQF(const REAL_TYPE var) const
  {
    return ( (Unit_Log == DIMENSIONAL) ? var*dhd*dhd*rhocp : var*dh*dh );
  }
  

  /**
   * @brief モードに対応する熱量を返す(体積要素)
   * @param [in] var  熱量
   */
  REAL_TYPE printQV(const REAL_TYPE var) const
  {
    return ( (Unit_Log == DIMENSIONAL) ? var*dhd*dhd*dhd*rhocp : var*dh*dh*dh );
  }
  

  /** 
   * @brief モードに対応する熱流束を返す
   * @param [in] var  熱流束
   */
  REAL_TYPE printHflux(const REAL_TYPE var) const
  {
    return ( (Unit_Log == DIMENSIONAL) ? var*RefVelocity : var );
  }
  
public:
  
  // 標準履歴の出力
  void printHistory(FILE* fp, const double* avr, const double* rms, const ItrCtl* IC, const Control* C, const double stptm, const bool disp);
  
  
  // 標準履歴モニタのヘッダー出力
  void printHistoryTitle(FILE* fp, const ItrCtl* IC, const Control* C, const bool disp);
  
  
  // コンポーネントモニタの履歴出力(dimensional value)
  void printHistoryCompo(FILE* fp, const CompoList* cmp, const Control* C);
  
  
  // コンポーネントモニタのヘッダー出力
  void printHistoryCompoTitle(FILE* fp, const CompoList* cmp, const Control* C);
  
  
  // 計算領域の流束履歴の出力
  void printHistoryDomfx(FILE* fp, const Control* C);
  
  
  // 計算領域の流束履歴のヘッダー出力
  void printHistoryDomfxTitle(FILE* fp, const Control* C);
  
  
  // 物体に働く力の履歴の出力
  void printHistoryForce(FILE* fp, const REAL_TYPE* force);
  
  
  // 物体に働く力の履歴のヘッダー出力
  void printHistoryForceTitle(FILE* fp);
  
  
  // 反復履歴出力
  void printHistoryItr(FILE* fp, const ItrCtl* IC, const FB::Vec3i idx);
  
  
  // 反復過程の状況モニタのヘッダー出力
  void printHistoryItrTitle(FILE* fp);
  
  
  // 壁面履歴の出力
  void printHistoryWall(FILE* fp, const REAL_TYPE* range_Yp, const REAL_TYPE* range_Ut);
  
  
  // 壁面モニタのヘッダー出力
  void printHistoryWallTitle(FILE* fp);
  
  
  // タイムスタンプの更新
  void updateTimeStamp(const int m_stp, const REAL_TYPE m_tm, const REAL_TYPE vMax);
};

#endif // _FB_HISTORY_H_
