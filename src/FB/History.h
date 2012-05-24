#ifndef _FB_HISTORY_H_
#define _FB_HISTORY_H_

/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file History.h
//@brief FlowBase History class Header
//@author keno, FSI Team, VCAD, RIKEN

#include "Control.h"
#include "Component.h"
#include "FBUtility.h"

class History {
protected:
  REAL_TYPE Tscale;          /// 時間スケール
  REAL_TYPE RefVelocity;     /// 代表速度
  REAL_TYPE BaseTemp;        /// 基準温度
  REAL_TYPE DiffTemp;        /// 代表温度差
  REAL_TYPE RefDensity;      /// 代表密度
  REAL_TYPE RefLength;       /// 代表長さ
  REAL_TYPE BasePrs;         /// 基準圧力
  REAL_TYPE RefSpecificHeat; /// 代表比熱
  REAL_TYPE time;            /// 時刻
  REAL_TYPE v_max;           /// 速度成分の最大値
  REAL_TYPE dhd;             /// 有次元の格子幅
  REAL_TYPE dh;              /// 無次元の格子幅
  REAL_TYPE rhocp;           /// 熱流束計算の無次元化量
  REAL_TYPE dynamic_p;       /// 動圧
  REAL_TYPE base_mf;         /// 流量の基準値
  unsigned step;             /// ステップ数
  unsigned Unit_Temp;        /// 温度単位
  unsigned Unit_Prs;         /// 圧力基準モード
  unsigned Unit_Log;         /// ログ出力の単位
  
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
    dynamic_p       = RefVelocity * RefVelocity * RefDensity;
    base_mf         = RefVelocity * RefLength * RefLength;

    time  = 0.0;
    v_max = 0.0;
    step  = 0;
  }
  ~History() {}
  
protected:
  //@fn REAL_TYPE printTime(void)
  //@brief モードに対応する時刻を返す
  REAL_TYPE printTime(void) {
    return ( (Unit_Log == DIMENSIONAL) ? time*Tscale : time );
  }
  
  //@fn REAL_TYPE printVmax(void)
  //@brief モードに対応する速度最大値を返す
  REAL_TYPE printVmax(void) {
    return ( (Unit_Log == DIMENSIONAL) ? v_max*RefVelocity : v_max );
  }
  
  //@fn REAL_TYPE printLen(const REAL_TYPE var)t
  //@brief モードに対応する長さを返す
  REAL_TYPE printLen(const REAL_TYPE var) {
    return ( (Unit_Log == DIMENSIONAL) ? var*RefLength : var );
  }
  
  //@fn REAL_TYPE printVel(const REAL_TYPE var)
  //@brief モードに対応する速度値を返す
  REAL_TYPE printVel(const REAL_TYPE var) {
    return ( (Unit_Log == DIMENSIONAL) ? var*RefVelocity : var );
  }
  
  //@fn REAL_TYPE printPrs(const REAL_TYPE var)
  //@brief 圧力の次元変換
  REAL_TYPE printPrs(const REAL_TYPE var) {
    if (Unit_Log != DIMENSIONAL) return var;
    const REAL_TYPE a = var * dynamic_p;
    return ( (Unit_Prs==Unit_Absolute) ? BasePrs+a : a );
  }
  
  //@fn REAL_TYPE printTP(const REAL_TYPE var)
  //@brief 全圧の次元変換
  REAL_TYPE printTP(const REAL_TYPE var) {
    return ( (Unit_Log == DIMENSIONAL) ? var*dynamic_p : var );
  }
  
  //@fn REAL_TYPE printForce(const REAL_TYPE var)
  //@brief 全圧の次元変換
  REAL_TYPE printForce(const REAL_TYPE var) {
    return ( (Unit_Log == DIMENSIONAL) ? var*dynamic_p*RefLength*RefLength : var );
  }
  
  //@fn REAL_TYPE printMF(const REAL_TYPE var)
  //@brief モードに対応する流量を返す
  REAL_TYPE printMF(const REAL_TYPE var) {
    return ( (Unit_Log == DIMENSIONAL) ? var*base_mf : var );
  }
  
  
  REAL_TYPE printTmp(const REAL_TYPE var) {
    if (Unit_Log != DIMENSIONAL) return var;
    const REAL_TYPE a = BaseTemp + DiffTemp*var; // Kelvin
    return ( (Unit_Temp==Unit_KELVIN) ? a : a-KELVIN  );
  }
  
  
  //@fn REAL_TYPE printQF(const REAL_TYPE var)
  //@brief モードに対応する熱量を返す(面表素)
  REAL_TYPE printQF(const REAL_TYPE var) {
    return ( (Unit_Log == DIMENSIONAL) ? var*dhd*dhd*rhocp : var*dh*dh );
  }
  
  //@fn REAL_TYPE printQV(const REAL_TYPE var)
  //@brief モードに対応する熱量を返す(体積要素)
  REAL_TYPE printQV(const REAL_TYPE var) {
    return ( (Unit_Log == DIMENSIONAL) ? var*dhd*dhd*dhd*rhocp : var*dh*dh*dh );
  }
  
  //@fn REAL_TYPE printHflux(const REAL_TYPE var)
  //@brief モードに対応する熱流束を返す
  REAL_TYPE printHflux(const REAL_TYPE var) {
    return ( (Unit_Log == DIMENSIONAL) ? var*RefVelocity : var );
  }
  
public:
  void printHistory          (FILE* fp, const REAL_TYPE* avr, const REAL_TYPE* rms, const ItrCtl* IC, const Control* C);
  void printHistoryTitle     (FILE* fp, const ItrCtl* IC, const Control* C);
  void printHistoryCompo     (FILE* fp, const CompoList* cmp, const Control* C);
  void printHistoryCompoTitle(FILE* fp, const CompoList* cmp, const Control* C);
  void printHistoryDomfx     (FILE* fp, const Control* C);
  void printHistoryForceTitle(FILE* fp, const Control* C);
  void printHistoryForce     (FILE* fp, const Control* C, REAL_TYPE* force);
  void printHistoryDomfxTitle(FILE* fp, const Control* C);
  void printHistoryItr       (FILE* fp, const unsigned itr, const REAL_TYPE nrm, const int* idx);
  void printHistoryItrTitle  (FILE* fp);
  void printHistoryWall      (FILE* fp, REAL_TYPE* range_Yp, REAL_TYPE* range_Ut);
  void printHistoryWallTitle (FILE* fp);
  void updateTimeStamp       (const unsigned m_stp, const REAL_TYPE m_tm, const REAL_TYPE vMax);
};

#endif // _FB_HISTORY_H_
