//##################################################################################
//
// Flow Base class
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2015 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2016 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
// Copyright (c) 2016-2020 Research Institute for Information Research, Kyushu University.
// All rights reserved.
//
//##################################################################################

/** 
 * @file   Component.C
 * @brief  FlowBase CompoList class
 * @author aics
 */

#include "Component.h"


// #################################################################
// BCのラベル名を返す
// @see ParseBC::setKeywordLBC()
std::string CompoList::getBCstr()
{
  std::string bc;
  
  if      ( type == ADIABATIC )     bc = "Adiabatic";
  else if ( type == OBSTACLE )      bc = "Obstacle";
  else if ( type == HEATFLUX )      bc = "Direct Heat Flux";
  else if ( type == ISOTHERMAL )    bc = "Isothermal";
  else if ( type == RADIANT )       bc = "Radiant";
  else if ( type == SPEC_VEL)
  {
    if ( heatmode == OFF ) bc = "Specified Velocity";
    else                   bc = "Specified Velocity with Temperature";
  }  
  else if ( type == OUTFLOW)        bc = "Outflow";
  else if ( type == IBM_DF )        bc = "Immersed Boundary Method(Direct Forcing)";
  else if ( type == HEAT_SRC )      bc = "Heat Source";
  else if ( type == CNST_TEMP )     bc = "Constant Temperature";
  else if ( type == HEX )           bc = "Pressure Loss";
  else if ( type == FAN )           bc = "Fan";
  else if ( type == SOLIDREV )      bc = "Solid Revolution";
  else if ( type == DARCY )         bc = "Darcy";
  else if ( type == PERIODIC )      bc = "Periodic";
  else if ( type == TRANSFER )
  {
    if      ( h_type == HT_S )      bc = "Heat Transfer type S";
    else if ( h_type == HT_SN)      bc = "Heat Transfer type SN (Natural convection)";
    else if ( h_type == HT_SF)      bc = "Heat Transfer type SF (Forced convection)";
  }
  else if ( type == OUTER_BC )      bc = "Outer BC";
  else                              bc = "Medium";
  
  return bc;
}


// #################################################################
// Polylibファイル用のBCのラベル名を返す
// @see ParseBC::setKeywordLBC()
std::string CompoList::getBCstr2Polylib()
{
  std::string bc;
  
  if      ( type == ADIABATIC )     bc = "Adiabatic";
  else if ( type == HEATFLUX )      bc = "DirectHeatFlux";
  else if ( type == TRANSFER )
  {
    if      ( h_type == HT_S )      bc = "HeatTransferS";
    else if ( h_type == HT_SN)      bc = "HeatTransferSN";
    else if ( h_type == HT_SF)      bc = "HeatTransferSF";
  }
  else if ( type == ISOTHERMAL )    bc = "Isothermal";
  else if ( type == RADIANT )       bc = "Radiation";
  else if ( type == SPEC_VEL)       bc = "SpecifiedVelocity";
  else if ( type == OUTFLOW)        bc = "Outflow";
  else if ( type == IBM_DF )        bc = "Forcing";
  else if ( type == HEAT_SRC )      bc = "HeatSource";
  else if ( type == CNST_TEMP )     bc = "SpecifiedTemperature";
  else if ( type == HEX )           bc = "PressureLoss";
  else if ( type == FAN )           bc = "Fan";
  else if ( type == DARCY )         bc = "Darcy";
  else if ( type == PERIODIC )      bc = "Periodic";
  else if ( type == OBSTACLE )      bc = "Obstacle";
  else if ( type == SOLIDREV )      bc = "SolidRevolution";

  return bc;
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
//@brief コンポーネントのBbox情報を設定する
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
void CompoList::setCoefHT(const REAL_TYPE var) 
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
//@brief コンポーネントが自ノードに存在しているかどうかをセットする
void CompoList::setEnsLocal(const int key)
{
  ens = key;
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
//@brief モニタ値の保持
void CompoList::setMonitorValue( const REAL_TYPE var) 
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
//@brief stateをセットする
//@param key セルの状態 SOLID/FLUID
void CompoList::setState(const int key)
{
  state = key;
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
// 温度の保持
void CompoList::setTemp(const REAL_TYPE var) 
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

