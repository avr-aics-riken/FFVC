#ifndef _FB_BND_OUTER_H_
#define _FB_BND_OUTER_H_

// #################################################################
//
// CAERU Library
//
// Copyright (c) 2012 All right reserved.
//
// Institute of Industrial Science, University of Tokyo, Japan. 
//
// #################################################################

/** 
 * @file   BndOuter.h
 * @brief  FlowBase BoundaryOuter Class Header
 * @author kero
 * @note メンバ変数に追加したら，dataCopy()処理にも加えること
 */

#include <string>
#include "FB_Define.h"
#include "cpm_Define.h"

class BoundaryOuter {
private:
  int BCclass;       ///< 境界条件の種類
  int wallType;      ///< wall >> (fixed, slide)
  int outType;       /// outflow >> 流出対流速度の評価モード（average, minmax）
  int drv_dir;       ///< ドライバーの方向
  int drv_lid;       ///< ドライバフェイスIDの位置
  int gc_medium;     ///< ガイドセルの媒質インデクス
  int v_profile;     ///< 速度プロファイル（constant, harmonic, zero）
  int Face_mode;     ///< 周期境界のときの面の状況指定（upstream, downstream）
  int hType;         ///< 熱境界条件の種別
  int HTref;         ///< 熱伝達境界の参照モード(Bulk, Local)
  int HTmode;        ///< 熱伝達境界の種別(HT_N, HT_B, HT_S, HT_SN, HT_SF)
  int Prdc_mode;     ///< 周期境界のモード（simple, directional, driver）
  int pType;         ///< 外部境界の圧力指定(ディリクレ，勾配ゼロ)
  int valid_cell;    ///< 境界面で流量計算に有効なセル数（Fluid cell）
  REAL_TYPE var1;    ///< 多目的用の変数(熱流束，熱伝達係数を共用するので排他的に使用)
  REAL_TYPE var2;    ///< 多目的用の変数(温度)
  REAL_TYPE dm[3];   ///< ローカルな計算領域境界面のモニタ値 (0-sum, 1-min, 2-max) コピー不要
  std::string label; ///< ラベル
  std::string alias; ///< 別名
  
public: 
  int mon_ref;       ///< IN_OUT境界条件のときのBC格納番号
  REAL_TYPE nv[3];   ///< 法線
  REAL_TYPE ca[5];   ///< 係数
  REAL_TYPE cb[5];   ///< 係数
  REAL_TYPE p;       ///< ワーク
  
  /** 周期境界の方向 */
  enum periodic_dir 
  {
    prdc_upstream,
    prdc_downstream
  };
  
  /** 周期境界の種類 */
  enum periodic_kind 
  {
    prdc_Simple,
    prdc_Directional,
    prdc_Driver
  };
  
  /** 壁面の種類 */
  enum wall_kind 
  {
    fixed,
    slide
  };
  
  /** コンストラクタ */
  BoundaryOuter() 
  {
    BCclass = drv_dir = HTref = wallType = 0;
    drv_lid = 0;
    outType = 0;
    mon_ref = pType = v_profile = hType = 0;
    HTmode = gc_medium = Prdc_mode = Face_mode = 0;
    p = var1 = var2 = 0.0;
    valid_cell = 0;
		for (int i=0; i<5; i++) ca[i] = cb[i] = 0.0;
    for (int i=0; i<3; i++) nv[i] = 0.0;
    for (int i=0; i<3; i++) dm[i]=0.0;
  }
  
  /**　デストラクタ */
  ~BoundaryOuter() {}
  
  
public:
  int get_Class() const
  { 
    return BCclass; 
  }
  
  int get_DriverDir() const
  { 
    return drv_dir;
  }
  
  int get_DriverIndex() const
  { 
    return drv_lid; 
  }
  
  int get_GuideMedium() const
  { 
    return gc_medium; 
  }
  
  int get_HTmodeRef() const
  { 
    return HTref; 
  }
  
  int get_FaceMode() const
  { 
    return Face_mode; 
  }
  
  int get_HTmode() const
  { 
    return HTmode; 
  }
  
  int get_hType() const
  { 
    return hType;
  }
  
  int get_MonRef() const
  {
    return mon_ref;
  }
  
  int get_ofv() const
  { 
    return outType; 
  }
  
  int get_PrdcMode() const
  { 
    return Prdc_mode;
  }
  
  int get_pType() const 
  { 
    return pType; 
  }
  
  int get_wallType() const 
  { 
    return wallType; 
  }
  
  int get_V_Profile() const 
  { 
    return v_profile; 
  }
  
  int get_ValidCell() const
  { 
    return valid_cell;
  }
  
  REAL_TYPE get_CoefHT() const 
  { 
    return var1; 
  }
  
  REAL_TYPE get_Heatflux() const
  {
    return var1;
  }
  
  REAL_TYPE get_Temp() const 
  { 
    return var1; 
  }
  
  std::string get_Label() const 
  { 
    return label; 
  }
  
  std::string get_Alias() const 
  { 
    return alias;
  }
  
  /**
   * @brief ローカルのモニタ積算値
   * @return dm[] (0-sum, 1-min, 2-max)
   */
  REAL_TYPE* getDomainV() 
  {
    return dm;
  }
  
  void addVec         (REAL_TYPE* vec);
  void dataCopy       (BoundaryOuter* src);
  void set_Alias      (std::string key);
  void set_Class      (const int key);
  void set_CoefHT     (REAL_TYPE val);
  
  
  // モニタ値を保持する
  void set_DomainV(const REAL_TYPE* vv, const int face, const int mode);
  
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
  void set_wallType   (const int key);
  void set_V_Profile  (const int key);
  void set_ValidCell  (const int val);
  
};

#endif // _FB_BND_OUTER_H_
