#ifndef _FB_COMPO_H_
#define _FB_COMPO_H_

// #################################################################
//
// CAERU Library
//
// Copyright (c) All right reserved. 2012
//
// Institute of Industrial Science, The University of Tokyo, Japan. 
//
// #################################################################

//@file Component.h
//@brief FlowBase CompoList class Header
//@author kero

#include <string>
#include <string.h>
#include "FB_Define.h"


class CompoList {
public:
  /// 速度指定方法
  enum velocity_spec_policy {
    BC_type_velocity,
    BC_type_massflow
  };
  
  /// 発熱項の指定の種別
  enum heat_src_type {
    hsrc_watt=1,
    hsrc_density
  };
  
  /// 熱伝達流束計算時の参照値モード
  enum HeatTransfer_Ref_mode {
    HT_mode_local=1,
    HT_mode_bulk
  };
  
  /// 変動速度指定時の変数指定子
  enum coef_ID {
    amplitude = 0,
    frequency,
    initphase,
    bias
  };
  
  /// 速度変化のモード
  enum velocity_variation {
    vel_constant,
    vel_harmonic,
    vel_zero
  };
  
  /// 熱伝達境界条件 自然対流　Vertical
  enum Natural_Vertical_HT {
    vert_laminar_alpha,
    vert_laminar_beta,
    vert_turbulent_alpha,
    vert_turbulent_beta,
    vert_Ra_critial
  };
  
  /// 熱伝達境界条件 自然対流　Horizontal
  enum Natural_Horizontal_HT {
    lower_laminar_alpha,
    lower_laminar_beta,
    lower_turbulent_alpha,
    lower_turbulent_beta,
    lower_Ra_critial
  };
  
  /// 熱伝達境界条件 強制対流
  enum Forced_HT {
    alpha,
    beta,
    gamma
  };
  
  /// 境界条件の入力単位指定
  enum BC_unit_policy {
    unit_mmAq=1,
    unit_mmHg,
    unit_Pa,
    unit_NonDimensional
  };
  
  /// BC位置
  enum BC_Location {
    same_direction=1,
    opposite_direction
  };

  
protected:
  unsigned long element;    /// 要素数
  unsigned type;            /// 
  unsigned attrb;           /// 
  unsigned h_type;          /// 
  unsigned variable;        ///
  unsigned ens;             /// 
  unsigned phase;           /// 
  unsigned var_u1;          /// 内部周期境界の方向，圧力単位指定，流出速度のタイプ，セルモニタの状態
  unsigned bc_dir;          /// VBCの指定方向
  
  int state;           /// Fluid(1) or Solid(0)
  int st[3];           /// コンポーネントインデクスBV範囲の始点
  int ed[3];           /// コンポーネントインデクスBV範囲の終点
  int c_size[3];       /// コンポーネントワーク配列の大きさ
  int def;             /// BC指定時の面を挟む相手先のセルID
  int mat_odr;         /// 媒質リストのインデクス
  int shape;           /// 形状パラメータ
  int sampling_method; /// サンプリングの方法（NEAREST, INTERPOLATION, SMOOTHING）
  int sampling_mode;   /// サンプリングモード（ALL, FLOW, SOLID）
  int usw;             /// 汎用変数
  
  REAL_TYPE var1;           /// パラメータ保持 (Velocity, Pressure, Massflow, Epsiolon of Radiation)
  REAL_TYPE var2;           /// パラメータ保持 (Heat Value, Heat flux, Heat Transfer, Pressure loss, Projection of Radiation)
  REAL_TYPE var3;           /// パラメータ保持 (Heat Density, Temperature)
  REAL_TYPE var_m;          /// モニタの値を保持
  REAL_TYPE temp_init;      /// 温度の初期値
  
  std::string name;         ///< ラベル
  
  
public:
  REAL_TYPE area;           ///< 断面積
  REAL_TYPE nv[3];          ///< 法線方向ベクトル（流出方向）
  REAL_TYPE oc[3];          ///< 形状の中心座標（前面の中心位置）
  REAL_TYPE dr[3];          ///< 補助方向ベクトル（nvと直交）
  REAL_TYPE depth;          ///< 厚さ（nv方向）
  REAL_TYPE shp_p1;         ///< 矩形の幅（dr方向), ファン半径
  REAL_TYPE shp_p2;         ///< 矩形の高さ, ボス半径
  REAL_TYPE val[var_END];   ///< データ保持用ワーク
  REAL_TYPE ca[6];          ///< 係数セット a
  REAL_TYPE cb[6];          ///< 係数セット b
  
  CompoList() {
    element = 0;
    type = variable = attrb = bc_dir = 0;
    mat_odr = 0;
    h_type = 0;
    state = shape = -1;
    def = 0;
    ens = OFF;
    area = 0.0;
    usw = 0;
    var_u1 = 0;
    phase = 0;
    for (int i=0; i<3; i++) {
      nv[i] = 0.0;
      st[i] = 0;
      ed[i] = 0;
      oc[i] = 0.0;
      dr[i] = 0.0;
      c_size[i] = 0;
    }
    for (int i=0; i<var_END; i++) val[i]=0.0;
    for (int i=0; i<6; i++) ca[i] = cb[i] = 0.0;
    var1 = var2 = var3 = var_m = temp_init = 0.0;
    depth = shp_p1 = shp_p2 = 0.0;
  }
  ~CompoList() {}
  
public:
  bool isFORCING   (void); 
  bool isHBC       (void);
  bool isHsrc      (void);
  bool isMONITOR   (void);
  bool isVBC       (void);
  bool isVBC_IO    (void);
  bool isVecForcing(void);
  bool isVFraction (void);
  
  REAL_TYPE get_CoefRadEps(void)  const { return var1; };
  REAL_TYPE get_CoefRadPrj(void)  const { return var2; };
  REAL_TYPE get_CoefMassflow(void)const { return var1; };
  REAL_TYPE get_CoefPrsLoss(void) const { return var2; };
  REAL_TYPE get_CoefHT(void)      const { return var2; };
  REAL_TYPE get_Heatflux(void)    const { return var2; };
  REAL_TYPE get_HeatDensity(void) const { return var3; };
  REAL_TYPE get_HeatValue(void)   const { return var2; };
  REAL_TYPE get_Massflow(void)    const { return var1; };
  REAL_TYPE get_Mon_Calorie(void) const { return var_m; };
  REAL_TYPE get_Mon_Heatflux(void)const { return var_m; };
  REAL_TYPE get_Mon_Temp(void)    const { return var_m; };
  REAL_TYPE get_Pressure(void)    const { return var1; };
  REAL_TYPE get_Temp(void)        const { return var3; };
  REAL_TYPE get_Velocity(void)    const { return var1; };
  
  std::string getVarStr (void);
  std::string getBCstr  (void);
  std::string getLabel  (void) { return name; }
  
  void setAttrb            (const unsigned key);
  void setBClocation       (const unsigned key);
  void setBbox             (const int m_st[], const int m_ed[]);
  void setDef              (const int key);
  void setElement          (const unsigned long key);
  void setEns              (const unsigned key);
  void setHtype            (const unsigned key);
  void setInitTemp         (const REAL_TYPE var);
  void setMatOdr           (const int key);
  void setLabel            (const std::string pnt);
  void setOutflowType      (const unsigned key);
  void setPeriodicDir      (const unsigned key);
  void setPrsUnit          (const unsigned key);
  void setPhase            (const unsigned m_phase);
  void setState            (const int key);
  void setStateCellMonitor (const unsigned key);
  void setType             (const unsigned key);
  void set_CoefMassflow    (const REAL_TYPE var);
  void set_CoefHT          (const REAL_TYPE var);
  void set_CoefPrsLoss     (const REAL_TYPE var);
  void set_CoefRadEps      (const REAL_TYPE var);
  void set_CoefRadPrj      (const REAL_TYPE var);
  void set_cmp_sz          (void);
  void set_Heatflux        (const REAL_TYPE var);
  void set_HeatDensity     (const REAL_TYPE var);
  void set_HeatValue       (const REAL_TYPE var);
  void set_HSRC_policy     (const bool kind);
  void set_Massflow        (const REAL_TYPE var);
  void set_Mon_Calorie     (const REAL_TYPE var);
  void set_Mon_Heatflux    (const REAL_TYPE var);
  void set_Mon_Temp        (const REAL_TYPE var);
  void set_Pressure        (const REAL_TYPE var);
  void set_SamplingMethod  (const int key);
  void set_SamplingMode    (const int key);
  void set_Shape           (const int key);
  void set_sw_Heatgen      (const int key);
  void set_sw_HexDir       (const int key);
  void set_sw_HTmodeRef    (const int key);
  void set_P_BCtype        (const int key);
  void set_Temp            (const REAL_TYPE var);
  void set_V_profile       (const int key);
  void set_VBC_policy      (const bool kind);
  void set_Velocity        (const REAL_TYPE var);
  
  inline REAL_TYPE getInitTemp(void) const        { return temp_init; }
  
  int get_sw_Heatgen(void)  const { return usw; };
  int get_sw_HexDir(void)   const { return usw; };
  int get_sw_HTmodeRef(void)const { return usw; };
  int get_P_BCtype(void)    const { return usw; };
  int get_V_Profile(void)   const { return usw; };
  inline int getDef(void) const                   { return def; };
  inline int getState(void) const                 { return state; }
  inline int get_Shape(void) const                { return shape; }
  inline int get_SamplingMethod(void) const       { return sampling_method; }
  inline int get_SamplingMode(void) const         { return sampling_mode; }
  inline int getMatOdr(void) const                { return mat_odr; }
  
  inline unsigned getBClocation(void) const       { return bc_dir; }
  inline unsigned getPeriodicDir(void) const      { return var_u1; }
  inline unsigned getPrsUnit(void) const          { return var_u1; }
  inline unsigned getOutflowType(void) const      { return var_u1; }
  inline unsigned getStateCellMonitor(void) const { return var_u1; }
  inline unsigned getType (void) const            { return type; }
  inline unsigned getHtype (void) const           { return h_type; }
  inline unsigned getAttrb(void) const            { return attrb; }
  
  inline unsigned long getElement(void) const     { return element; }
  
  inline bool isPolicy_Massflow(void) const       { return (attrb==BC_type_massflow) ? true : false; }
  inline bool isPolicy_HeatDensity(void) const    { return (usw==hsrc_density) ? true : false; }
  
  inline bool isEns(void) const                   { return ( (ens == ON) ? true : false ); }
  
  //@fn unsigned getPhase(void) const
  //@brief return pahse ID (SOLID=0, FLUID=1, GAS=2, LIQUID=3)
  inline unsigned getPhase(void) const            { return phase; }
  
  //@fn void encodeVarType (unsigned var)
  //@brief メンバ変数variableに変数の種類をエンコードする
  //@param var 変数タイプを表すビット数
  inline void encodeVarType (const unsigned var)  { variable |= (0x1 << var); }
  
  //@fn bool isVarEncoded (unsigned var) const
  //@brief メンバ変数variableにエンコードされた変数タイプとの照合を行う
  //@param var 変数タイプのビット数
  inline bool isVarEncoded (const unsigned var) const   { return ( ( 0x1 & (variable >> var)) ? true : false ); }

  //@fn void getBbox(int* m_st, int* m_ed)
  //@brief コンポーネントのBV情報を返す
  inline void getBbox(int* m_st, int* m_ed) {
    m_st[0] = st[0];
    m_st[1] = st[1];
    m_st[2] = st[2];
    m_ed[0] = ed[0];
    m_ed[1] = ed[1];
    m_ed[2] = ed[2];
  }
  
  //@fn int* getBbox_st(void)
  //@brief コンポーネントのBV情報stのアドレスを返す
  int* getBbox_st(void) { return (st); }
  
  //@fn int* getBbox_ed(void)
  //@brief コンポーネントのBV情報edのアドレスを返す
  int* getBbox_ed(void) { return (ed); }
  
  //@fn void get_cmp_sz(int* m_sz)
  //@brief コンポーネントのサイズを返す
  void get_cmp_sz(int* m_sz) { 
    m_sz[0] = c_size[0];
    m_sz[1] = c_size[1];
    m_sz[2] = c_size[2];
  }
};

#endif // _FB_COMPO_H_
