#ifndef _SKL_FB_HISTORY_H_
#define _SKL_FB_HISTORY_H_

/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2010
 *
 */

//@file History.h
//@brief FlowBase History class Header
//@author keno, FSI Team, VCAD, RIKEN

#include "IterationCtl.h"
#include "Control.h"
#include "Component.h"
#include "FBUtility.h"

class History {
protected:
  SKL_REAL Tscale;          /// 時間スケール
  SKL_REAL RefVelocity;     /// 代表速度
  SKL_REAL BaseTemp;        /// 基準温度
  SKL_REAL DiffTemp;        /// 代表温度差
  SKL_REAL RefDensity;      /// 代表密度
  SKL_REAL RefLength;       /// 代表長さ
  SKL_REAL BasePrs;         /// 基準圧力
  SKL_REAL RefSpecificHeat; /// 代表比熱
  SKL_REAL time;            /// 時刻
  SKL_REAL v_max;           /// 速度成分の最大値
  SKL_REAL dhd;             /// 有次元の格子幅
  SKL_REAL dh;              /// 無次元の格子幅
  SKL_REAL rhocp;           /// 熱流束計算の無次元化量
  unsigned step;            /// ステップ数
  unsigned Unit_Temp;       /// 温度単位
  unsigned Unit_Prs;        /// 圧力基準モード
  unsigned Unit_Log;        /// ログ出力の単位
  
public:
  History() {}
  History(const Control* Cref) {
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
    dh              = Cref->dh;
    Tscale          = RefLength / RefVelocity;
    dhd             = dh*RefLength;
    rhocp           = RefVelocity * DiffTemp * RefDensity * RefSpecificHeat;
    
    time  = 0.0;
    v_max = 0.0;
    step  = 0;
  }
  ~History() {}
  
protected:
  //@fn SKL_REAL printTime(void) const
  //@brief モードに対応する時刻を返す
  SKL_REAL printTime(void) const {
    return ( (Unit_Log == DIMENSIONAL) ? time*Tscale : time );
  }
  
  //@fn SKL_REAL printVmax(void) const
  //@brief モードに対応する速度最大値を返す
  SKL_REAL printVmax(void) const {
    return ( (Unit_Log == DIMENSIONAL) ? v_max*RefVelocity : v_max );
  }
  
  //@fn SKL_REAL printLen(const SKL_REAL var) const
  //@brief モードに対応する長さを返す
  SKL_REAL printLen(const SKL_REAL var) const {
    return ( (Unit_Log == DIMENSIONAL) ? var*RefLength : var );
  }
  
  //@fn SKL_REAL printVel(const SKL_REAL var) const
  //@brief モードに対応する速度値を返す
  SKL_REAL printVel(const SKL_REAL var) const {
    return ( (Unit_Log == DIMENSIONAL) ? var*RefVelocity : var );
  }
  
  //@fn SKL_REAL printMF(const SKL_REAL var) const
  //@brief モードに対応する流量を返す
  SKL_REAL printMF(const SKL_REAL var) const {
    const SKL_REAL cf = RefVelocity * RefLength * RefLength;
    return ( (Unit_Log == DIMENSIONAL) ? var*cf : var );
  }
  
  //@fn SKL_REAL printQF(const SKL_REAL var) const
  //@brief モードに対応する熱量を返す(面表素)
  SKL_REAL printQF(const SKL_REAL var) const {
    return ( (Unit_Log == DIMENSIONAL) ? var*dhd*dhd*rhocp : var*dh*dh );
  }
  
  //@fn SKL_REAL printQV(const SKL_REAL var) const
  //@brief モードに対応する熱量を返す(体積要素)
  SKL_REAL printQV(const SKL_REAL var) const {
    return ( (Unit_Log == DIMENSIONAL) ? var*dhd*dhd*dhd*rhocp : var*dh*dh*dh );
  }
  
  //@fn SKL_REAL printHflux(const SKL_REAL var) const
  //@brief モードに対応する熱流束を返す
  SKL_REAL printHflux(const SKL_REAL var) const {
    return ( (Unit_Log == DIMENSIONAL) ? var*RefVelocity : var );
  }
  
public:
  void printHistory          (FILE* fp, const SKL_REAL* delta, const ItrCtl* IC, const Control* C) const;
  void printHistoryTitle     (FILE* fp, const ItrCtl* IC, const Control* C) const;
  void printHistoryCompo     (FILE* fp, const CompoList* cmp, const Control* C) const;
  void printHistoryCompoTitle(FILE* fp, const CompoList* cmp, const Control* C) const;
  void printHistoryDomfx     (FILE* fp, const Control* C) const;
  void printHistoryDomfxTitle(FILE* fp, const Control* C) const;
  void printHistoryItr       (FILE* fp, const unsigned itr, const SKL_REAL nrm, const int* idx) const;
  void printHistoryItrTitle  (FILE* fp) const;
  void printHistoryWall      (FILE* fp, SKL_REAL* range_Yp, SKL_REAL* range_Ut) const;
  void printHistoryWallTitle (FILE* fp) const;
  void updateTimeStamp       (const unsigned m_stp, const SKL_REAL m_tm, const SKL_REAL vMax);
};

#endif // _SKL_FB_HISTORY_H_