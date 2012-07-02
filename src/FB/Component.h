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

/** 
 * @file Component.h
 * @brief FlowBase CompoList class Header
 * @author kero
 */

#include <string>
#include <string.h>
#include "FB_Define.h"
#include "cpm_Define.h"


class CompoList {
public:
  /** 速度指定方法 */
  enum velocity_spec_policy 
  {
    BC_type_velocity,
    BC_type_massflow
  };
  
  /** 発熱項の指定の種別 */
  enum heat_src_type 
  {
    hsrc_watt=1,
    hsrc_density
  };
  
  /** 熱伝達流束計算時の参照値モード */
  enum HeatTransfer_Ref_mode 
  {
    HT_mode_local=1,
    HT_mode_bulk
  };
  
  /** 変動速度指定時の変数指定子 */
  enum coef_ID 
  {
    amplitude = 0,
    frequency,
    initphase,
    bias
  };
  
  /** 速度変化のモード */
  enum velocity_variation 
  {
    vel_constant,
    vel_harmonic,
    vel_zero
  };
  
  /** 熱伝達境界条件 自然対流　Vertical */
  enum Natural_Vertical_HT 
  {
    vert_laminar_alpha,
    vert_laminar_beta,
    vert_turbulent_alpha,
    vert_turbulent_beta,
    vert_Ra_critial
  };
  
  /** 熱伝達境界条件 自然対流　Horizontal */
  enum Natural_Horizontal_HT 
  {
    lower_laminar_alpha,
    lower_laminar_beta,
    lower_turbulent_alpha,
    lower_turbulent_beta,
    lower_Ra_critial
  };
  
  /** 熱伝達境界条件 強制対流 */
  enum Forced_HT 
  {
    alpha,
    beta,
    gamma
  };
  
  /** 境界条件の入力単位指定 */
  enum BC_unit_policy 
  {
    unit_mmAq=1,
    unit_mmHg,
    unit_Pa,
    unit_NonDimensional
  };
  
  /** BC位置 */
  enum BC_Location 
  {
    same_direction=1,
    opposite_direction
  };

  
private:
  unsigned long element; ///< 要素数
  
  int type;              /// 
  int attrb;             /// 
  int h_type;            /// 
  int variable;          ///
  int ens;               /// 
  int phase;             /// 
  int var_u1;            /// 内部周期境界の方向，圧力単位指定，流出速度のタイプ，セルモニタの状態
  int bc_dir;            /// VBCの指定方向
  int state;             /// Fluid(1) or Solid(0)
  int st[3];             /// コンポーネントインデクスBV範囲の始点
  int ed[3];             /// コンポーネントインデクスBV範囲の終点
  int c_size[3];         /// コンポーネントワーク配列の大きさ
  int def;               /// BC指定時の面を挟む相手先のセルID
  int mat_odr;           /// 媒質リストのインデクス
  int shape;             /// 形状パラメータ
  int sampling_method;   /// サンプリングの方法（NEAREST, INTERPOLATION, SMOOTHING）
  int sampling_mode;     /// サンプリングモード（ALL, FLOW, SOLID）
  int usw;               /// 汎用変数
  
  REAL_TYPE var1;        /// パラメータ保持 (Velocity, Pressure, Massflow, Epsiolon of Radiation)
  REAL_TYPE var2;        /// パラメータ保持 (Heat Value, Heat flux, Heat Transfer, Pressure loss, Projection of Radiation)
  REAL_TYPE var3;        /// パラメータ保持 (Heat Density, Temperature)
  REAL_TYPE var_m;       /// モニタの値を保持
  REAL_TYPE temp_init;   /// 温度の初期値
  
  std::string name;      ///< ラベル
  
  
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
  
  /** コンストラクタ */
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
  
  /**　デストラクタ */
  ~CompoList() {}
  
public:
  bool isFORCING   (); 
  bool isHBC       ();
  bool isHsrc      ();
  bool isMONITOR   ();
  bool isVBC       ();
  bool isVBC_IO    ();
  bool isVecForcing();
  bool isVFraction ();
  
  REAL_TYPE get_CoefRadEps() const
  { 
    return var1; 
  }
  
  REAL_TYPE get_CoefRadPrj() const
  { 
    return var2; 
  }
  
  REAL_TYPE get_CoefMassflow() const
  { 
    return var1; 
  }
  
  REAL_TYPE get_CoefPrsLoss() const
  { 
    return var2; 
  }
  
  REAL_TYPE get_CoefHT() const
  { 
    return var2; 
  }
  
  REAL_TYPE get_Heatflux() const
  { 
    return var2; 
  }
  
  REAL_TYPE get_HeatDensity() const
  { 
    return var3; 
  }
  
  REAL_TYPE get_HeatValue() const
  { 
    return var2; 
  }
  
  REAL_TYPE get_Massflow() const
  { 
    return var1; 
  }
  
  REAL_TYPE get_Mon_Calorie() const
  { 
    return var_m; 
  }
  
  REAL_TYPE get_Mon_Heatflux() const
  { 
    return var_m; 
  }
  
  REAL_TYPE get_Mon_Temp() const 
  { 
    return var_m; 
  }
  
  REAL_TYPE get_Pressure() const
  { 
    return var1; 
  }
  
  REAL_TYPE get_Temp() const
  { 
    return var3; 
  }
  
  REAL_TYPE get_Velocity() const
  { 
    return var1; 
  }
  
  
  /**
   * @brief 変数名を返す
   */
  const char* getVarStr();
  
  
  /**
   * @brief BCのラベル名を返す
   */
  std::string getBCstr();
  
  
  /**
   * @brief ラベル名を返す
   */
  std::string getLabel() const
  { 
    return name; 
  }
  
  void setAttrb            (const int key);
  void setBClocation       (const int key);
  void setBbox             (const int m_st[], const int m_ed[]);
  void setDef              (const int key);
  void setElement          (const unsigned long key);
  void setEns              (const int key);
  void setHtype            (const int key);
  void setInitTemp         (const REAL_TYPE var);
  void setMatOdr           (const int key);
  void setLabel            (const std::string pnt);
  void setOutflowType      (const int key);
  void setPeriodicDir      (const int key);
  void setPrsUnit          (const int key);
  void setPhase            (const int m_phase);
  void setState            (const int key);
  void setStateCellMonitor (const int key);
  void setType             (const int key);
  void set_CoefMassflow    (const REAL_TYPE var);
  void set_CoefHT          (const REAL_TYPE var);
  void set_CoefPrsLoss     (const REAL_TYPE var);
  void set_CoefRadEps      (const REAL_TYPE var);
  void set_CoefRadPrj      (const REAL_TYPE var);
  void set_cmp_sz          ();
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
  
  REAL_TYPE getInitTemp() const
  { 
    return temp_init; 
  }
  
  int get_sw_Heatgen() const
  { 
    return usw; 
  }
  
  int get_sw_HexDir() const
  { 
    return usw; 
  }
  
  int get_sw_HTmodeRef() const
  { 
    return usw; 
  }
  
  int get_P_BCtype() const
  {
    return usw; 
  }
  
  int get_V_Profile() const
  { 
    return usw; 
  }
  
  int getDef() const
  { 
    return def; 
  }
  
  int getState() const
  { 
    return state; 
  }
  
  int get_Shape() const
  { 
    return shape; 
  }
  
  int get_SamplingMethod() const
  { 
    return sampling_method; 
  }
  
  int get_SamplingMode() const
  { 
    return sampling_mode; 
  }
  
  int getMatOdr() const
  { 
    return mat_odr;
  }
  
  int getBClocation() const
  { 
    return bc_dir; 
  }
  
  int getPeriodicDir() const
  { 
    return var_u1; 
  }
  
  int getPrsUnit() const
  { 
    return var_u1; 
  }
  
  int getOutflowType() const
  { 
    return var_u1; 
  }
  
  int getStateCellMonitor() const 
  { 
    return var_u1; 
  }
  
  int getType () const
  { 
    return type; 
  }
  
  int getHtype () const
  { 
    return h_type; 
  }
  
  int getAttrb() const
  { 
    return attrb; 
  }
  
  unsigned long getElement() const
  { 
    return element; 
  }
  
  bool isPolicy_Massflow() const
  { 
    return (attrb==BC_type_massflow) ? true : false; 
  }
  
  bool isPolicy_HeatDensity() const
  { 
    return (usw==hsrc_density) ? true : false; 
  }
  
  bool isEns() const
  { 
    return ( (ens == ON) ? true : false ); 
  }
  

  //@brief return pahse ID (SOLID=0, FLUID=1, GAS=2, LIQUID=3)
  int getPhase() const
  { 
    return phase; 
  }
  

  //@brief メンバ変数variableに変数の種類をエンコードする
  //@param var 変数タイプを表すビット数
  void encodeVarType (const int var)
  { 
    variable |= (0x1 << var); 
  }
  

  //@brief メンバ変数variableにエンコードされた変数タイプとの照合を行う
  //@param var 変数タイプのビット数
  bool isVarEncoded (const int var) const
  { 
    return ( ( 0x1 & (variable >> var)) ? true : false ); 
  }


  //@brief コンポーネントのBV情報を返す
  inline void getBbox(int* m_st, int* m_ed) 
  {
    m_st[0] = st[0];
    m_st[1] = st[1];
    m_st[2] = st[2];
    m_ed[0] = ed[0];
    m_ed[1] = ed[1];
    m_ed[2] = ed[2];
  }
  

  //@brief コンポーネントのBV情報stのアドレスを返す
  int* getBbox_st()
  { 
    return st; 
  }
  

  //@brief コンポーネントのBV情報edのアドレスを返す
  int* getBbox_ed() 
  { 
    return ed; 
  }
  

  //@brief コンポーネントのサイズを返す
  void get_cmp_sz(int* m_sz) 
  { 
    m_sz[0] = c_size[0];
    m_sz[1] = c_size[1];
    m_sz[2] = c_size[2];
  }
};

#endif // _FB_COMPO_H_
