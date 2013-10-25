#ifndef _FB_COMPO_H_
#define _FB_COMPO_H_

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
 * @file   Component.h
 * @brief  FlowBase CompoList class Header
 * @author kero
 */

#include "cpm_Define.h"
#include <string>
#include <string.h>
#include "FB_Define.h"


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
  
  int type;              ///< 境界条件の種類 媒質の場合には0
  int attrb;             /// 
  int h_type;            ///
  int ens;               ///
  int phase;             ///
  int var_u1;            /// 内部周期境界の方向，圧力単位指定，セルモニタの状態
  int bc_dir;            /// VBCの指定方向
  int state;             ///< Fluid(1) or Solid(0)
  int st[3];             /// コンポーネントインデクスBbox範囲の始点
  int ed[3];             /// コンポーネントインデクスBbox範囲の終点
  int c_size[3];         /// コンポーネントワーク配列の大きさ
  int usw;               /// 汎用変数
  int heatmode;          /// 熱輸送のときON
  
  REAL_TYPE var1;        /// パラメータ保持 (Velocity, Pressure, Massflow, Epsiolon of Radiation)
  REAL_TYPE var2;        /// パラメータ保持 (Heat Value, Heat flux, Heat Transfer, Pressure loss, Projection of Radiation)
  REAL_TYPE var3;        /// パラメータ保持 (Heat Density, Temperature)
  REAL_TYPE var_m;       /// モニタの値を保持
  REAL_TYPE temp_init;   /// 温度の初期値
  
  std::string alias;     ///< 局所境界条件の別名
  std::string medium;    ///< Medium名
  
public:
  REAL_TYPE area;           ///< 断面積（有次元）
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
    type = 0;
    attrb = bc_dir = 0;
    h_type = 0;
    state = -1;
    ens = OFF;
    area = 0.0;
    usw = 0;
    var_u1 = 0;
    phase = 0;
    heatmode = OFF;
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
  
  /**
   * @brief ラベル名を返す
   */
  std::string getAlias() const
  {
    return alias;
  }

  
  int getAttrb() const
  {
    return attrb;
  }
  
  
  //@brief コンポーネントのBbox情報を返す
  inline void getBbox(int* m_st, int* m_ed)
  {
    m_st[0] = st[0];
    m_st[1] = st[1];
    m_st[2] = st[2];
    m_ed[0] = ed[0];
    m_ed[1] = ed[1];
    m_ed[2] = ed[2];
  }
  
  
  //@brief コンポーネントのBbox情報edのアドレスを返す
  int* getBbox_ed()
  {
    return ed;
  }
  
  
  //@brief コンポーネントのBbox情報stのアドレスを返す
  int* getBbox_st()
  {
    return st;
  }
  
  
  /**
   * @brief BCのラベル名を返す
   */
  std::string getBCstr();
  
  
  //@brief コンポーネントのサイズを返す
  void get_cmp_sz(int* m_sz)
  {
    m_sz[0] = c_size[0];
    m_sz[1] = c_size[1];
    m_sz[2] = c_size[2];
  }
  
  
  // サブドメインにコンポーネントが存在するかどうかを返す
  bool existLocal() const
  {
    return ( (ens == ON) ? true : false );
  }
  
  
  
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
  
  REAL_TYPE getCoefHT() const
  { 
    return var2; 
  }
  
  REAL_TYPE getHeatflux() const
  { 
    return var2; 
  }
  
  REAL_TYPE getHeatDensity() const
  { 
    return var3; 
  }
  
  REAL_TYPE get_HeatValue() const
  { 
    return var2; 
  }
  
  REAL_TYPE getInitTemp() const
  {
    return temp_init;
  }
  
  
  REAL_TYPE get_Massflow() const
  { 
    return var1; 
  }
  
  
  /**
   * @brief 媒質名を返す
   */
  std::string getMedium() const
  {
    return medium;
  }
  
  
  REAL_TYPE getMonCalorie() const
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
  
  //@brief return pahse ID (SOLID=0, FLUID=1, GAS=2, LIQUID=3)
  int getPhase() const
  {
    return phase;
  }
  
  
  REAL_TYPE get_Pressure() const
  { 
    return var1; 
  }
  
  int get_sw_Heatgen() const
  {
    return usw;
  }
  
  int get_sw_HexDir() const
  {
    return usw;
  }
  
  
  REAL_TYPE getTemp() const
  { 
    return var3; 
  }
  
  REAL_TYPE get_Velocity() const
  { 
    return var1; 
  }
  
  
  
  
  int get_P_BCtype() const
  {
    return usw;
  }
  
  int get_V_Profile() const
  {
    return usw;
  }
  
  int getState() const
  {
    return state;
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
  
  int getType () const
  {
    return type;
  }
  
  int getHtype () const
  {
    return h_type;
  }
  
  
  unsigned long getElement() const
  {
    return element;
  }
  
  
  /**
   * @brief 流体要素であればtrue
   */
  bool isFluid() const
  {
    return (state == FLUID) ? true : false;
  }
  
  
  /**
   @brief 境界条件タイプがFORCINGかどうかを調べる
   */
  bool isFORCING() const
  {
    if ((type == HEX) ||
        (type == FAN) ||
        (type == DARCY) ) return true;
    return false;
  }
  
  
  // @brief 熱問題のときture
  bool isHeatMode() const
  {
    return (heatmode==ON) ? true : false;
  }
  
  
  // @brief 境界条件タイプが熱源かどうかを調べる
  bool isHsrc() const
  {
    if ((type == HEAT_SRC) ||
        (type == CNST_TEMP) ) return true;
    return false;
  }
  
  
  /**
   * @brief typeがコンポーネントであればtrue
   * @note OBSTACLE=1, BCは2以上，非設定=-1（媒質）
   */
  bool isKindCompo() const
  {
    return (type > OBSTACLE) ? true : false;
  }
  
  
  /**
   * @brief typeが媒質であればtrue
   */
  bool isKindMedium() const
  {
    return (type == 0) ? true : false;
  }
  
  
  /**
   * @brief typeがOBSTACLEであればtrue
   */
  bool isKindObstacle() const
  {
    return (type == OBSTACLE) ? true : false;
  }
  
  
  bool isPolicy_HeatDensity() const
  {
    return (usw==hsrc_density) ? true : false;
  }
  
  
  bool isPolicy_Massflow() const
  {
    return (attrb==BC_type_massflow) ? true : false;
  }
  
  
  // @brief 内部境界条件タイプがコンポーネントモニタかどうかを調べる
  bool isCompoMonitor() const
  {
    if ((type == SPEC_VEL) ||
        (type == OUTFLOW) ||
        (type == HEX) ||
        (type == DARCY) ||
        (type == FAN) ||
        (type == DARCY) ||
        (type == HEATFLUX) ||
        (type == TRANSFER) ||
        (type == ISOTHERMAL) ||
        (type == RADIANT)
        ) return true;
    return false;
  }
  
  
  // @brief 内部境界条件タイプが速度指定かどうかを調べる
  bool isVBC() const
  {
    if ((type == SPEC_VEL) ||
        (type == OUTFLOW) ||
        (type == IBM_DF) ||
        (type == HEX) ||
        (type == FAN) ||
        (type == DARCY) ) return true;
    return false;
  }
  
  // @brief コンポーネントが速度規定
  bool isVBC_IO() const
  {
    if ((type == SPEC_VEL) ||
        (type == OUTFLOW) ) return true;
    return false;
  }
  
  
  // @brief 体積率の必要なコンポーネントかどうか
  bool isVFraction() const
  {
    if ((type == HEAT_SRC) ||
        (type == CNST_TEMP) ||
        (type == IBM_DF) ||
        (type == HEX) ||
        (type == FAN) ||
        (type == DARCY) )  return true;
    return false;
  }
  
  
  // @brief ベクトル強制をするかどうかを調べる
  bool isVecForcing() const
  {
    if ( isFORCING() && usw==ON ) return true;
    
    return false;
  }
  

  void setAlias            (const std::string pnt);
  
  void setAttrb            (const int key);
  void setBClocation       (const int key);
  void setBbox             (const int m_st[], const int m_ed[]);
  void set_CoefMassflow    (const REAL_TYPE var);
  void setCoefHT           (const REAL_TYPE var);
  void set_CoefPrsLoss     (const REAL_TYPE var);
  void set_CoefRadEps      (const REAL_TYPE var);
  void set_CoefRadPrj      (const REAL_TYPE var);
  void set_cmp_sz          ();
  void setElement          (const unsigned long key);
  
  // サブドメインにコンポーネントが存在するかどうかを設定
  void setEnsLocal         (const int key);
  
  void setHeatflux         (const REAL_TYPE var);
  
  
  // @brief 吸発熱密度の保持
  void setHeatDensity (const REAL_TYPE var);
  
  
  // @brief 熱問題の指定
  void setHeatmode (const int mode)
  {
    heatmode = mode;
  }
  
  
  // @brief 吸発熱量の保持
  void setHeatValue (const REAL_TYPE var);
  
  
  // @brief 発熱項の指定ポリシーを指定する
  // @param [in] kind ポリシー種別　true-発熱量, false-発熱密度
  void setHsrcPolicy (const bool kind);
  
  void setHtype            (const int key);
  void setInitTemp         (const REAL_TYPE var);
  
  void set_Massflow        (const REAL_TYPE var);
  void setMedium           (const std::string pnt);
  
  // @brief モニター値を保持する
  void setMonitorValue (const REAL_TYPE var);
  
  
  void set_P_BCtype        (const int key);
  void setPeriodicDir      (const int key);
  void setPhase            (const int m_phase);
  void set_Pressure        (const REAL_TYPE var);
  void setPrsUnit          (const int key);
  void setSamplingWidth    (const int key);
  void setState            (const int key);
  void setType             (const int key);

  void set_sw_Heatgen      (const int key);
  void set_sw_HexDir       (const int key);
  
  
  // @brief 温度の保持
  void setTemp (const REAL_TYPE var);
  
  
  void set_V_profile       (const int key);
  void set_VBC_policy      (const bool kind);
  void set_Velocity        (const REAL_TYPE var);
  

};

#endif // _FB_COMPO_H_
