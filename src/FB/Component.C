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
 * @file   Component.C
 * @brief  FlowBase CompoList class
 * @author kero
 */

#include "Component.h"


// #################################################################
// 変数名を返す
const char* CompoList::getVarStr()
{
  std::string var;
  
  if ( isVarEncoded(var_Velocity) )         var += "Velocity ";              // 0
  if ( isVarEncoded(var_Pressure) )         var += "Pressure ";              // 1
  if ( isVarEncoded(var_Temperature) )      var += "Temperature ";           // 2
  if ( isVarEncoded(var_Density) )          var += "MassDensity ";           // 3
  if ( isVarEncoded(var_TotalP) )           var += "TotalPressure ";         // 4
  if ( isVarEncoded(var_Velocity_Avr) )     var += "Averaged Velocity ";     // 5
  if ( isVarEncoded(var_Pressure_Avr) )     var += "Averaged Pressure ";     // 6
  if ( isVarEncoded(var_Temperature_Avr) )  var += "Averaged Temperature ";  // 7
  if ( isVarEncoded(var_Density_Avr) )      var += "Averaged MassDensity ";  // 8
  if ( isVarEncoded(var_TotalP_Avr) )       var += "Averaged TotalPressure ";// 9
  
  return var.c_str();
}


// #################################################################
//BCのラベル名を返す
std::string CompoList::getBCstr() const
{
  std::string bc;
  
  if      ( type == ADIABATIC )     bc = "Adiabatic";
  else if ( type == OBSTACLE )      bc = "Obstacle";
  else if ( type == HEATFLUX )      bc = "Direct Heat Flux";
  else if ( type == ISOTHERMAL )    bc = "Isothermal";
  else if ( type == RADIANT )       bc = "Radiant";
  else if ( type == SPEC_VEL)       bc = "Specified Velocity";
  else if ( type == SPEC_VEL_WH)    bc = "Specified Velocity with Temperature";
  else if ( type == OUTFLOW)        bc = "Outflow";
  else if ( type == IBM_DF )        bc = "Immersed Boundary Method(Direct Forcing)";
  else if ( type == HEAT_SRC )      bc = "Heat Source";
  else if ( type == CNST_TEMP )     bc = "Constant Temperature";
  else if ( type == HEX )           bc = "Pressure Loss";
  else if ( type == FAN )           bc = "Fan";
  else if ( type == DARCY )         bc = "Darcy";
  else if ( type == CELL_MONITOR )  bc = "Cell Monitor";
  else if ( type == PERIODIC )      bc = "Periodic";
  else if ( type == INACTIVE )      bc = "Inactive";
  else if ( type == TRANSFER )
  {
    if      ( h_type == HT_N )      bc = "Heat Transfer type N";
    else if ( h_type == HT_S )      bc = "Heat Transfer type S";
    else if ( h_type == HT_SN)      bc = "Heat Transfer type SN (Natural convection)";
    else if ( h_type == HT_SF)      bc = "Heat Transfer type SF (Forced convection)";
    else if ( h_type == HT_B )      bc = "Heat Transfer type B";
  }
  else                              bc = "Medium";
  
  return bc;
}


// #################################################################
/**
 @brief 境界条件タイプが熱境界条件かどうかを調べる
 @retval HBCであればtrue
 */
bool CompoList::isHBC() const
{
  if ((type == ADIABATIC)  || 
      (type == HEATFLUX)   ||
      (type == TRANSFER)   ||
      (type == ISOTHERMAL) ||
      (type == RADIANT)    ||
      (type == SPEC_VEL_WH)||
      (type == HEAT_SRC)   ||
      (type == CNST_TEMP) ) return true;
  return false;
}

// #################################################################
/**
 @brief 境界条件タイプが熱源かどうかを調べる
 @retval FORCINGであればtrue
 */
bool CompoList::isHsrc() const
{
  if ((type == HEAT_SRC) || 
      (type == CNST_TEMP) ) return true;
  return false;
}



// #################################################################
/**

 @brief 内部境界条件タイプが速度指定かどうかを調べる
 @retval VBCであればtrue
 */
bool CompoList::isVBC() const
{
  if ((type == SPEC_VEL) ||
      (type == SPEC_VEL_WH) ||
      (type == OUTFLOW) ||
      (type == IBM_DF) ||
      (type == HEX) ||
      (type == FAN) ||
      (type == DARCY) ) return true;
  return false;
}



// #################################################################
//@brief 体積率の必要なコンポーネントかどうか
bool CompoList::isVFraction() const 
{
  if ((type == HEAT_SRC) ||
      (type == CNST_TEMP) ||
      (type == IBM_DF) ||
      (type == HEX) ||
      (type == FAN) ||
      (type == DARCY) )  return true;
  return false;
}


// #################################################################
//@brief エイリアス名をセットする
void CompoList::setAlias(const std::string pnt)
{
  alias = pnt;
}



// #################################################################
//@brief attrbをセットする
//@param key アトリビュート
void CompoList::setAttrb(const int key) 
{ 
  attrb = key; 
}


// #################################################################
//@brief BCの方向をセットする
//@param key CDSの場合のBCの方向
void CompoList::setBClocation(const int key) 
{ 
  bc_dir = key; 
}


// #################################################################
//@brief コンポーネントのBV情報を設定する
void CompoList::setBbox(const int m_st[], const int m_ed[]) 
{
  st[0] = m_st[0];
  st[1] = m_st[1];
  st[2] = m_st[2];
  ed[0] = m_ed[0];
  ed[1] = m_ed[1];
  ed[2] = m_ed[2];
}


// #################################################################
//@brief コンポーネントのサイズを計算
void CompoList::set_cmp_sz(void)
{
  c_size[0] = ed[0] - st[0] + 1;
  c_size[1] = ed[1] - st[1] + 1;
  c_size[2] = ed[2] - st[2] + 1;
}


// #################################################################
//@brief 熱伝達係数の保持
void CompoList::set_CoefHT(const REAL_TYPE var) 
{
  var2 = var;
}


// #################################################################
//@brief 輻射のイプシロンの保持
void CompoList::set_CoefRadEps(const REAL_TYPE var) 
{
  var1 = var;
}


// #################################################################
//@brief 輻射の射出係数の保持
void CompoList::set_CoefRadPrj(const REAL_TYPE var) 
{
  var2 = var;
}


// #################################################################
//@brief 流量の有次元化係数の保持
void CompoList::set_CoefMassflow(const REAL_TYPE var) 
{
  var1 = var;
}


// #################################################################
//@brief 圧力損失の有次元化係数の保持
void CompoList::set_CoefPrsLoss(const REAL_TYPE var) 
{
  var2 = var;
}


// #################################################################
//@brief 指定セルを保持する
void CompoList::setDef(const int key)
{
  def = key;
}


// #################################################################
//@brief elementをセットする
//@param key 要素数 element
void CompoList::setElement(const unsigned long key)
{
  element = key;
}


// #################################################################
//@brief コンポーネントが自ノードに存在しているかどうかをセットする
void CompoList::setEnsLocal(const int key)
{
  ens = key;
}



// #################################################################
//@brief 吸発熱密度の保持
void CompoList::set_HeatDensity(const REAL_TYPE var) 
{
  var3 = var;
}


// #################################################################
//@brief 熱流束の保持
void CompoList::set_Heatflux(const REAL_TYPE var) 
{
  var2 = var;
}


// #################################################################
//@brief 吸発熱量の保持
void CompoList::set_HeatValue(const REAL_TYPE var) 
{
  var2 = var;
}


// #################################################################
//@brief 発熱項の指定ポリシーをセットする
//@param kind ポリシー種別　true-発熱量指定, false-発熱密度指定
void CompoList::set_HSRC_policy(const bool kind) 
{
  usw = ( kind ) ? hsrc_watt : hsrc_density;
}


// #################################################################
//@brief h_typeをセットする
//@param key 境界条件の種類
void CompoList::setHtype(const int key)
{
  h_type = key;
}


// #################################################################
//@brief 初期温度の指定
void CompoList::setInitTemp(const REAL_TYPE key)
{
  temp_init = key;
}



// #################################################################
//@brief 流量の保持
void CompoList::set_Massflow(const REAL_TYPE var) 
{
  var1 = var;
}



// #################################################################
//@brief 媒質名をセットする
void CompoList::setMedium(const std::string pnt)
{
  medium = pnt;
}


// #################################################################
//@brief モニタ温度の保持
void CompoList::set_Mon_Temp(const REAL_TYPE var) 
{
  var_m = var;
}


// #################################################################
//@brief モニタ熱流束の保持
void CompoList::set_Mon_Heatflux(const REAL_TYPE var) 
{
  var_m = var;
}


// #################################################################
//@brief モニタ熱量の保持
void CompoList::set_Mon_Calorie(const REAL_TYPE var) 
{
  var_m = var;
}


// #################################################################
//@brief 圧力境界条件タイプ指定モードの保持
void CompoList::set_P_BCtype(const int key)
{
  usw = key;
}


// #################################################################
//@brief 周期境界の上流方向を保持する
void CompoList::setPeriodicDir(const int key)
{
  var_u1 = key;
}


// #################################################################
//@brief set pahse ID (SOLID=0, FLUID=1, GAS=2, LIQUID=3)
void CompoList::setPhase(const int m_phase)
{
  phase = m_phase;
}


// #################################################################
//@brief 圧力値の保持
void CompoList::set_Pressure(const REAL_TYPE var) 
{
  var1 = var;
}


// #################################################################
//@brief 圧力の単位を指定する
void CompoList::setPrsUnit(const int key)
{
  var_u1 = key;
}


// #################################################################
void CompoList::setSamplingWidth(const int key)
{
  sampling_width = key;
}


// #################################################################
void CompoList::set_Shape(const int key)
{
  shape = key;
}


// #################################################################
//@brief stateをセットする
//@param key セルの状態 SOLID/FLUID
void CompoList::setState(const int key)
{
  state = key;
}


// #################################################################
//@brief セルモニタスイッチ ON/OFF
void CompoList::setStateCellMonitor(const int key)
{
  var_u1 = key;
}


// #################################################################
//@brief 熱伝達の参照指定モードの保持
void CompoList::set_sw_HTmodeRef(const int key) 
{
  usw = key;
}


// #################################################################
//@brief 熱交換機の方向指定モードの保持
void CompoList::set_sw_HexDir (const int key) 
{
  usw = key;
}


// #################################################################
//@brief 発熱量指定モードの保持
void CompoList::set_sw_Heatgen(const int key) 
{
  usw = key;
}


// #################################################################
//@brief 温度の保持
void CompoList::set_Temp(const REAL_TYPE var) 
{
  var3 = var;
}


// #################################################################
//@brief typeをセットする
//@param key 境界条件の種類
void CompoList::setType(const int key)
{
  type = key;
}


// #################################################################
//@brief 速度指定ポリシーをセットする
//@param kind ポリシー種別　true-速度指定, false-流量指定
void CompoList::set_VBC_policy(const bool kind)
{
  attrb = ( kind ) ? BC_type_velocity : BC_type_massflow;
}


// #################################################################
//@brief 速度の保持
void CompoList::set_Velocity(const REAL_TYPE var) 
{
  var1 = var;
}


// #################################################################
//@brief 速度プロファイル指定モードの保持
void CompoList::set_V_profile(const int key)
{
  usw = key;
}

