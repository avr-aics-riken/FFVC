#ifndef _SKL_FB_BND_OUTER_H_
#define _SKL_FB_BND_OUTER_H_

/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file BndOuter.h
//@brief FlowBase BoundaryOuter class Header
//@author keno, FSI Team, VCAD, RIKEN
//@todo 平均速度のみ，1/2(min+max)は未処理
//@note メンバ変数に追加したら，dataCopy()処理にも加えること

#include "Skl.h"
#include "FBDefine.h"

class BoundaryOuter {
protected:
  int HTref;            /// 熱伝達境界の参照モード(Bulk, Local)
  int BC_ID;            /// 基本境界条件リストのID
  int BCtype;           /// 境界条件の種類
  int drv_dir;          /// ドライバーの方向
  int drv_lid;          /// ドライバフェイスIDの位置
  int gc_id;            /// ガイドセルのID
  int gc_medium;        /// ガイドセルの媒質ID
  unsigned Face_mode;   /// 周期境界のときの面の状況指定（upstream, downstream）
  unsigned hType;       /// 熱境界条件の種別
  unsigned HTmode;      /// 熱伝達境界の種別(HT_N, HT_B, HT_S, HT_SN, HT_SF)
  unsigned oflowType;   /// 流出，流入出境界条件のときの流出対流速度の評価モード（average, minmax）
  unsigned Prdc_mode;   /// 周期境界のモード（simple, directional, driver）
  unsigned pType;       /// 外部境界の圧力指定(ディリクレ，勾配ゼロ)
  unsigned vType;       /// 速度プロファイル（constant, harmonic, zero）
  unsigned valid_cell;  /// 境界面で流量計算に有効なセル数（Fluid cell）
  REAL_TYPE var1;       /// 多目的用の変数(熱流束，熱伝達係数を共用するので排他的に使用)
  REAL_TYPE var2;       /// 多目的用の変数(温度)
  REAL_TYPE dm[3];      /// ローカルな計算領域境界面のモニタ値 (0-sum, 1-min, 2-max) コピー不要
  
public: 
  unsigned mon_ref;     /// IN_OUT境界条件のときのBC格納番号
  unsigned Face_inout;  /// 
  REAL_TYPE nv[3];      /// 
  REAL_TYPE ca[5];      /// 
  REAL_TYPE cb[5];      /// 
  REAL_TYPE p;          ///  
  
  enum periodic_dir {
    prdc_upstream,
    prdc_downstream
  };
  
  enum periodic_kind {
    prdc_Simple,
    prdc_Directional,
    prdc_Driver
  };
  
  BoundaryOuter() {
    BCtype = BC_ID = drv_dir = HTref = 0;
    drv_lid = 0;
    mon_ref = pType = vType = hType = oflowType = 0;
    HTmode = gc_medium = gc_id = Prdc_mode = Face_mode = Face_inout = 0;
    p = var1 = var2 = 0.0;
    valid_cell = 0;
		for (int i=0; i<5; i++) ca[i] = cb[i] = 0.0;
    for (int i=0; i<3; i++) nv[i] = 0.0;
    for (int i=0; i<3; i++) dm[i]=0.0;
  }
  ~BoundaryOuter() {}
  
public:
  int get_BC_ID(void)         const { return BC_ID; };
  int get_BCtype(void)        const { return BCtype; };
  int get_DriverDir(void)     const { return drv_dir; };
  int get_DriverIndex(void)   const { return drv_lid; };
  int get_GuideID(void)       const { return gc_id; };
  int get_GuideMedium(void)   const { return gc_medium; };
  int get_HTmodeRef(void)     const { return HTref; };
  
  REAL_TYPE get_CoefHT(void)   const { return var1; };
  REAL_TYPE get_Heatflux(void) const { return var1; };
  REAL_TYPE get_Temp(void)     const { return var1; };
  
  unsigned get_FaceMode(void) const { return Face_mode; };
  unsigned get_HTmode(void)   const { return HTmode; };
  unsigned get_hType(void)    const { return hType; };
  unsigned get_MonRef(void)   const { return mon_ref; };
  unsigned get_ofv(void)      const { return oflowType; };
  unsigned get_PrdcMode(void) const { return Prdc_mode; };
  unsigned get_pType(void)    const { return pType; };
  unsigned get_vType(void)    const { return vType; };
  unsigned get_ValidCell(void)const { return valid_cell; };
  
  //@fn void getDomainV(REAL_TYPE* vv)
  //@brief ローカルのモニタ積算値
  //@ret vv[3] (0-sum, 1-min, 2-max)
  //@param face 外部境界面番号
  void getDomainV(REAL_TYPE* vv) {
    vv[0] = dm[0];
    vv[1] = dm[1];
    vv[2] = dm[2];
  }
  
  void addVec         (REAL_TYPE* vec);
  void dataCopy       (BoundaryOuter* src);
  void set_BC_ID      (int key);
  void set_BCtype     (int key);
  void set_CoefHT     (REAL_TYPE val);
  void set_DomainV    (REAL_TYPE* vv, int face, bool mode=false);
  void set_DriverDir  (int key);
  void set_DriverIndex(int key);
  void set_FaceMode   (unsigned key);
  void set_GuideID    (int key);
  void set_GuideMedium(int key);
  void set_Heatflux   (REAL_TYPE val);
  void set_HTmode     (unsigned key);
  void set_HTmodeRef  (int key);
  void set_hType      (unsigned key);
  void set_MonRef     (unsigned key);
  void set_ofv        (unsigned key);
  void set_PrdcMode   (unsigned key);
  void set_pType      (unsigned key);
  void set_Temp       (REAL_TYPE val);
  void set_vType      (unsigned key);
  void set_ValidCell  (unsigned val);
  
};

#endif // _SKL_FB_BND_OUTER_H_
