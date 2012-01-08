#ifndef _SKL_FB_COMPO_H_
#define _SKL_FB_COMPO_H_

/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2010
 *
 */

//@file Component.h
//@brief FlowBase CompoList class Header
//@author keno, FSI Team, VCAD, RIKEN

#include <string>

#include "Skl.h"
#include "SklSolverBase.h"
#include "FBDefine.h"
#include "mydebug.h"

using namespace std;

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
  
  /// 圧力単位
  enum Unit_Pressure {
    Gauge=1,
    Absolute
  };
  
  /// 温度単位
  enum Temp_Unit {
    Unit_KELVIN=1,
    Unit_CELSIUS
  };
  
  /// 境界条件の入力単位指定
  enum BC_unit_policy {
    unit_mmAq=1,
    unit_mmHg,
    unit_Pa,
    unit_NonDimensional
  };
  
protected:
  unsigned ID;       /// セルID
  unsigned type;     /// 
  unsigned element;  /// 要素数
  unsigned attrb;    /// 
  unsigned h_type;   /// 
  unsigned variable; ///
  unsigned mat_odr;  /// 
  unsigned ens;      /// 
  unsigned phase;    /// 
  unsigned var_u1;   /// 内部周期境界の方向，圧力単位指定，流出速度のタイプ，セルモニタの状態
  unsigned usw;      ///
  int      state;    ///
  int      st[3];    /// コンポーネントインデクスBV範囲の始点
  int      ed[3];    /// コンポーネントインデクスBV範囲の終点
  int      def;      /// BC指定時の面を挟む相手先のセルID
  SKL_REAL var1;     /// パラメータ保持 (Velocity, Pressure, Massflow, Epsiolon of Radiation)
  SKL_REAL var2;     /// パラメータ保持 (Heat Value, Heat flux, Heat Transfer, Pressure loss, Projection of Radiation)
  SKL_REAL var3;     /// パラメータ保持 (Heat Density, Temperature)
  SKL_REAL var_m;    /// モニタの値を保持
  SKL_REAL temp_init;/// 温度の初期値
  
public:
  SKL_REAL dir;
  SKL_REAL area;
  SKL_REAL nv[3];
  SKL_REAL val[var_END];
  SKL_REAL ca[6];
  SKL_REAL cb[6];
  unsigned sub_domain;
  char     name[LABEL];
  
  CompoList() {
    ID = type = element = variable = mat_odr = attrb = 0;
    h_type = 0;
    state = -1;
    def = 0;
    ens = OFF;
    dir = area = 0.0;
    usw = OFF;
    var_u1 = 0;
    sub_domain = 0;
    phase = 0;
    for (int i=0; i<3; i++) {
      nv[i]=0.0;
      st[i]=0.0;
      ed[i]=0.0;
    }
    for (int i=0; i<var_END; i++) val[i]=0.0;
    for (int i=0; i<6; i++) ca[i] = cb[i] = 0.0;
    for (int n=0; n<LABEL; n++) name[n]='\0';
    var1 = var2 = var3 = var_m = temp_init = 0.0;
  }
  ~CompoList() {}
  
public:
  bool isFORCING   (void) const; 
  bool isHBC       (void) const;
  bool isHsrc      (void) const;
  bool isMONITOR   (void) const;
  bool isVBC       (void) const;
  bool isVecForcing(void) const;
  bool isVFraction (void) const;
  
  SKL_REAL get_CoefRadEps(void)  const { return var1; };
  SKL_REAL get_CoefRadPrj(void)  const { return var2; };
  SKL_REAL get_CoefMassflow(void)const { return var1; };
  SKL_REAL get_CoefPrsLoss(void) const { return var2; };
  SKL_REAL get_CoefHT(void)      const { return var2; };
  SKL_REAL get_Heatflux(void)    const { return var2; };
  SKL_REAL get_HeatDensity(void) const { return var3; };
  SKL_REAL get_HeatValue(void)   const { return var2; };
  SKL_REAL get_Massflow(void)    const { return var1; };
  SKL_REAL get_Mon_Calorie(void) const { return var_m; };
  SKL_REAL get_Mon_Heatflux(void)const { return var_m; };
  SKL_REAL get_Mon_Temp(void)    const { return var_m; };
  SKL_REAL get_Pressure(void)    const { return var1; };
  SKL_REAL get_Temp(void)        const { return var3; };
  SKL_REAL get_Velocity(void)    const { return var1; };
  
  unsigned get_sw_Heatgen(void)  const { return usw; };
  unsigned get_sw_HexDir(void)   const { return usw; };
  unsigned get_sw_HTmodeRef(void)const { return usw; };
  unsigned get_sw_V_profile(void)const { return usw; };
  unsigned get_sw_P_BCtype(void) const { return usw; };
  
  string getVarStr (void) const;
  string getBCstr  (void) const;
  
  void addVec              (SKL_REAL* vec);
  void setAttrb            (unsigned key);
  void setCompoBV_ed       (unsigned odr, int val);
  void setCompoBV_st       (unsigned odr, int val);
  void setDef              (int key);
  void setElement          (unsigned key);
  void setEns              (unsigned key);
  void setHtype            (unsigned key);
  void setID               (unsigned key);
  void setInitTemp         (SKL_REAL var);
  void setMatOdr           (unsigned key);
  void setName             (const char* pnt);
  void setOutflowType      (unsigned key);
  void setPeriodicDir      (unsigned key);
  void setPrsUnit          (unsigned key);
  void setPhase            (unsigned m_phase);
  void setState            (int key);
  void setStateCellMonitor (unsigned key);
  void setType             (unsigned key);
  void set_CoefMassflow    (SKL_REAL var);
  void set_CoefHT          (SKL_REAL var);
  void set_CoefPrsLoss     (SKL_REAL var);
  void set_CoefRadEps      (SKL_REAL var);
  void set_CoefRadPrj      (SKL_REAL var);
  void set_Heatflux        (SKL_REAL var);
  void set_HeatDensity     (SKL_REAL var);
  void set_HeatValue       (SKL_REAL var);
  void set_HSRC_policy     (bool kind);
  void set_Massflow        (SKL_REAL var);
  void set_Mon_Calorie     (SKL_REAL var);
  void set_Mon_Heatflux    (SKL_REAL var);
  void set_Mon_Temp        (SKL_REAL var);
  void set_Pressure        (SKL_REAL var);
  void set_sw_Heatgen      (unsigned key);
  void set_sw_HexDir       (unsigned key);
  void set_sw_HTmodeRef    (unsigned key);
  void set_sw_P_BCtype     (unsigned key);
  void set_sw_V_profile    (unsigned key);
  void set_Temp            (SKL_REAL var);
  void set_VBC_policy      (bool kind);
  void set_Velocity        (SKL_REAL var);
  
  //@fn SKL_REAL getInitTemp(void) const
  inline SKL_REAL getInitTemp(void) const {
    return temp_init;
  }
  
  //@fn int getDef(void) const
  inline int getDef(void) const { return def; };
  
  //@fn unsigned getPeriodicDir(void) const
  inline unsigned getPeriodicDir(void) const {
    return var_u1;
  }
  
  //@fn unsigned getPrsUnit(void) const
  inline unsigned getPrsUnit(void) const {
    return var_u1;
  }
  
  //@fn unsigned getOutflowType(void) const
  inline unsigned getOutflowType(void) const {
    return var_u1;
  }
  
  //@fn unsigned getStateCellMonitor(void) const
  inline unsigned getStateCellMonitor(void) const {
    return var_u1;
  }
  
  //@fn unsigned getPhase(void) const
  //@brief return pahse ID (SOLID=0, FLUID=1, GAS=2, LIQUID=3)
  inline unsigned getPhase(void) const {
    return phase;
  }
  
  /**
   @fn void encodeVarType (unsigned var)
   @brief メンバ変数variableに変数の種類をエンコードする
   @param var 変数タイプを表すビット数
   */
  inline void encodeVarType (unsigned var) { 
    variable |= (0x1 << var); 
  }
  
  /**
   @fn bool isVarEncoded (unsigned var) const
   @brief メンバ変数variableにエンコードされた変数タイプとの照合を行う
   @param var 変数タイプのビット数
   @retval 一致すればtrue
   */
  inline bool isVarEncoded (unsigned var) const { 
    return ( ( 0x1 & (variable >> var)) ? true : false );
  }
  
  //@fn int getState(void) const
  //@brief セルの状態 SOLID/FLUIDを返す
  inline int getState(void) const { return state; }
  
  //@fn unsigned getElement(void) const
  //@brief 要素数を返す
  inline unsigned getElement(void) const { return element; }
  
  //@fn unsigned getID(void) const
  //@brief IDを返す
  inline unsigned getID(void) const { return ID; }
  
  //@fn unsigned getMatOdr(void) const
  //@brief MaterialListのエントリを返す
  inline unsigned getMatOdr(void) const { return mat_odr; }
  
  //@fn unsigned getType(void) const
  //@brief 境界条件の種類typeを返す
  inline unsigned getType (void) const { return type; }
  
  //@fn unsigned getHtype(void) const
  //@brief 熱伝達境界条件の種類h_typeを返す
  inline unsigned getHtype (void) const { return h_type; }
  
  //@fn unsigned getAttrb(void) const
  inline unsigned getAttrb(void) const { return attrb; }
  
  //@fn bool isPolicy_Massflow(void) const
  //@brief ポリシーが流量の場合trueを返す
  inline bool isPolicy_Massflow(void) const { return (attrb==BC_type_massflow) ? true : false; }
  
  //@fn bool isPolicy_HeatDensity(void) const
  //@brief ポリシーが発熱密度の場合trueを返す
  inline bool isPolicy_HeatDensity(void) const { return (usw==hsrc_density) ? true : false; }
  
  //@fn bool isEns(void) const
  //@brief コンポーネントが自ノードに存在しているかどうかを返す
  inline bool isEns(void) const {
    return ( (ens == ON) ? true : false );
  }
  
  //@fn void getCompoBV(int* m_st, int* m_ed)
  //@brief コンポーネントのBV情報を返す
  inline void getCompoBV(int* m_st, int* m_ed) {
    m_st[0] = st[0];
    m_st[1] = st[1];
    m_st[2] = st[2];
    m_ed[0] = ed[0];
    m_ed[1] = ed[1];
    m_ed[2] = ed[2];
  }
  
  //@fn int getCompoBV_st_x(void)
  //@brief コンポーネントのBV情報st_xを返す
  inline int getCompoBV_st_x(void) {
    return (st[0]);
  }
  
  //@fn int getCompoBV_st_y(void)
  //@brief コンポーネントのBV情報st_yを返す
  inline int getCompoBV_st_y(void) {
    return (st[1]);
  }
  
  //@fn int getCompoBV_st_z(void)
  //@brief コンポーネントのBV情報st_zを返す
  inline int getCompoBV_st_z(void) {
    return (st[2]);
  }
  
  //@fn int getCompoBV_ed_x(void)
  //@brief コンポーネントのBV情報ed_xを返す
  inline int getCompoBV_ed_x(void) {
    return (ed[0]);
  }
  
  //@fn int getCompoBV_ed_y(void)
  //@brief コンポーネントのBV情報ed_yを返す
  inline int getCompoBV_ed_y(void) {
    return (ed[1]);
  }
  
  //@fn int getCompoBV_ed_z(void)
  //@brief コンポーネントのBV情報ed_zを返す
  inline int getCompoBV_ed_z(void) {
    return (ed[2]);
  }
  
  //@fn int* getCompoBV_adrs_st(void)
  //@brief コンポーネントのBV情報stのアドレスを返す
  int* getCompoBV_adrs_st(void) {
    return (st);
  }
  
  //@fn int* getCompoBV_adrs_ed(void)
  //@brief コンポーネントのBV情報edのアドレスを返す
  int* getCompoBV_adrs_ed(void) {
    return (ed);
  }
};

#endif // _SKL_FB_COMPO_H_