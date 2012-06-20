#ifndef _FB_BND_OUTER_H_
#define _FB_BND_OUTER_H_

// #################################################################
//
// CAERU Library
//
// Copyright (c) All right reserved. 2012
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
//
// #################################################################

//@file BndOuter.h
//@brief FlowBase BoundaryOuter Class Header
//@author kero
//@note メンバ変数に追加したら，dataCopy()処理にも加えること

#include "FB_Define.h"
#include <string>

class BoundaryOuter {
protected:
  int BCclass;       /// 境界条件の種類
  int subType;       /// サブタイプ（汎用）
                     /// outflow >> 流出対流速度の評価モード（average, minmax）
                     /// wall >> (fixed, slide)
  int drv_dir;       /// ドライバーの方向
  int drv_lid;       /// ドライバフェイスIDの位置
  int gc_medium;     /// ガイドセルの媒質インデクス
  int v_profile;     /// 速度プロファイル（constant, harmonic, zero）
  int Face_mode;     /// 周期境界のときの面の状況指定（upstream, downstream）
  int hType;         /// 熱境界条件の種別
  int HTref;         /// 熱伝達境界の参照モード(Bulk, Local)
  int HTmode;        /// 熱伝達境界の種別(HT_N, HT_B, HT_S, HT_SN, HT_SF)
  int Prdc_mode;     /// 周期境界のモード（simple, directional, driver）
  int pType;         /// 外部境界の圧力指定(ディリクレ，勾配ゼロ)
  int valid_cell;    /// 境界面で流量計算に有効なセル数（Fluid cell）
  REAL_TYPE var1;    /// 多目的用の変数(熱流束，熱伝達係数を共用するので排他的に使用)
  REAL_TYPE var2;    /// 多目的用の変数(温度)
  REAL_TYPE dm[3];   /// ローカルな計算領域境界面のモニタ値 (0-sum, 1-min, 2-max) コピー不要
  std::string label; /// ラベル
  std::string alias; /// 別名
  
public: 
  int mon_ref;       /// IN_OUT境界条件のときのBC格納番号
  REAL_TYPE nv[3];   /// 
  REAL_TYPE ca[5];   /// 
  REAL_TYPE cb[5];   /// 
  REAL_TYPE p;       ///  
  
  enum periodic_dir {
    prdc_upstream,
    prdc_downstream
  };
  
  enum periodic_kind {
    prdc_Simple,
    prdc_Directional,
    prdc_Driver
  };
  
  enum wall_kind {
    fixed,
    slide
  };
  
  BoundaryOuter() {
    BCclass = drv_dir = HTref = subType = 0;
    drv_lid = 0;
    mon_ref = pType = v_profile = hType = 0;
    HTmode = gc_medium = Prdc_mode = Face_mode = 0;
    p = var1 = var2 = 0.0;
    valid_cell = 0;
		for (int i=0; i<5; i++) ca[i] = cb[i] = 0.0;
    for (int i=0; i<3; i++) nv[i] = 0.0;
    for (int i=0; i<3; i++) dm[i]=0.0;
  }
  ~BoundaryOuter() {}
  
public:
  int get_Class(void)       const { return BCclass; };
  int get_DriverDir(void)   const { return drv_dir; };
  int get_DriverIndex(void) const { return drv_lid; };
  int get_GuideMedium(void) const { return gc_medium; };
  int get_HTmodeRef(void)   const { return HTref; };
  int get_FaceMode(void)    const { return Face_mode; };
  int get_HTmode(void)      const { return HTmode; };
  int get_hType(void)       const { return hType; };
  int get_MonRef(void)      const { return mon_ref; };
  int get_ofv(void)         const { return subType; };
  int get_PrdcMode(void)    const { return Prdc_mode; };
  int get_pType(void)       const { return pType; };
  int get_Type(void)        const { return subType; };
  int get_V_Profile(void)   const { return v_profile; };
  int get_ValidCell(void)   const { return valid_cell; };
  
  REAL_TYPE get_CoefHT(void)   const { return var1; };
  REAL_TYPE get_Heatflux(void) const { return var1; };
  REAL_TYPE get_Temp(void)     const { return var1; };
  
  std::string get_Label(void) const { return label; };
  std::string get_Alias(void) const { return alias; };
  
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
  void set_Alias      (std::string key);
  void set_Class      (const int key);
  void set_CoefHT     (REAL_TYPE val);
  void set_DomainV    (REAL_TYPE* vv, int face, bool mode=false);
  void set_DriverDir  (int key);
  void set_DriverIndex(int key);
  void set_FaceMode   (int key);
  void set_GuideMedium(int key);
  void set_Heatflux   (REAL_TYPE val);
  void set_HTmode     (int key);
  void set_HTmodeRef  (int key);
  void set_hType      (int key);
  void set_Label      (std::string key);
  void set_MonRef     (int key);
  void set_ofv        (int key);
  void set_PrdcMode   (int key);
  void set_pType      (int key);
  void set_Temp       (REAL_TYPE val);
  void set_Type       (const int key);
  void set_V_Profile  (const int key);
  void set_ValidCell  (const int val);
  
};

#endif // _FB_BND_OUTER_H_
