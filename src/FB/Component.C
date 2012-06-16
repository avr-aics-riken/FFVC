/*
 * SPHERE - Skeleton for PHysical and Engineering REsearch
 *
 * Copyright (c) RIKEN, Japan. All right reserved. 2004-2012
 *
 */

//@file Component.C
//@brief FlowBase CompoList class
//@author keno, FSI Team, VCAD, RIKEN

#include "Component.h"

/**
 @fn std::string CompoList::getVarStr(void)
 @brief 変数名を返す
 @retval 変数の文字列
 */
std::string CompoList::getVarStr(void)
{
  std::string var;
  
  if ( isVarEncoded(var_Velocity) )         var += "Velocity ";              // 0
  if ( isVarEncoded(var_Pressure) )         var += "Pressure ";              // 1
  if ( isVarEncoded(var_Temperature) )      var += "Temperature ";           // 2
  if ( isVarEncoded(var_Density) )          var += "Density ";               // 3
  if ( isVarEncoded(var_TotalP) )           var += "TotalPressure ";         // 4
  if ( isVarEncoded(var_Velocity_Avr) )     var += "Averaged Velocity ";     // 5
  if ( isVarEncoded(var_Pressure_Avr) )     var += "Averaged Pressure ";     // 6
  if ( isVarEncoded(var_Temperature_Avr) )  var += "Averaged Temperature ";  // 7
  if ( isVarEncoded(var_Density_Avr) )      var += "Averaged Density ";      // 8
  if ( isVarEncoded(var_TotalP_Avr) )       var += "Averaged TotalPressure ";// 9
  
  return var;
}

/**
 @fn std::string CompoList::getBCstr(void)
 @brief BCのラベル名を返す
 */
std::string CompoList::getBCstr(void)
{
  std::string bc;
  if      ( type == ADIABATIC )     bc = "Adiabatic";
  else if ( type == HEATFLUX )      bc = "Direct Heat Flux";
  else if ( type == ISOTHERMAL )    bc = "Isothermal";
  else if ( type == RADIANT )       bc = "Radiant";
  else if ( type == SPEC_VEL)       bc = "Specified_Velocity";
  else if ( type == SPEC_VEL_WH)    bc = "Specified_Velocity with Temperature";
  else if ( type == OUTFLOW)        bc = "Outflow";
  else if ( type == IBM_DF )        bc = "Immersed Boundary Method(Direct Forcing)";
  else if ( type == HEAT_SRC )      bc = "Heat Source";
  else if ( type == CNST_TEMP )     bc = "Constant Temperature";
  else if ( type == HEX )           bc = "Pressure Loss";
  else if ( type == FAN )           bc = "Fan";
  else if ( type == DARCY )         bc = "Darcy";
  else if ( type == CELL_MONITOR )  bc = "Cell_Monitor";
  else if ( type == PERIODIC )      bc = "Periodic";
  else if ( type == INACTIVE )      bc = "Inactive";
  else if ( type == TRANSFER ) {
    if      ( h_type == HT_N )      bc = "Heat Transfer type N";
    else if ( h_type == HT_S )      bc = "Heat Transfer type S";
    else if ( h_type == HT_SN)      bc = "Heat Transfer type SN (Natural convection)";
    else if ( h_type == HT_SF)      bc = "Heat Transfer type SF (Forced convection)";
    else if ( h_type == HT_B )      bc = "Heat Transfer type B";
  }
  else                              bc = "";
  return bc;
}

/**
 @fn bool CompoList::isFORCING(void)
 @brief 境界条件タイプがFORCINGかどうかを調べる
 @retval FORCINGであればtrue
 */
bool CompoList::isFORCING(void)
{
  if ((type == HEX) || 
      (type == FAN) || 
      (type == DARCY) ) return true;
  return false;
}

/**
 @fn bool CompoList::isHBC(void)
 @brief 境界条件タイプが熱境界条件かどうかを調べる
 @retval HBCであればtrue
 */
bool CompoList::isHBC(void)
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

/**
 @fn bool CompoList::isHsrc(void)
 @brief 境界条件タイプが熱源かどうかを調べる
 @retval FORCINGであればtrue
 */
bool CompoList::isHsrc(void)
{
  if ((type == HEAT_SRC) || 
      (type == CNST_TEMP) ) return true;
  return false;
}

//@fn bool CompoList::isVBC_IO(void)
bool CompoList::isVBC_IO(void)
{
  if ((type == SPEC_VEL) ||
      (type == SPEC_VEL_WH) ||
      (type == OUTFLOW) ) return true;
  return false;
}

/**
 @fn bool CompoList::isMONITOR(void)
 @brief コンポーネントタイプがモニタかどうかを調べる
 @retval モニタであればtrue
 */
bool CompoList::isMONITOR(void) {
  return ( (type == CELL_MONITOR) ? true : false );
}

/**
 @fn bool CompoList::isVBC(void)
 @brief 内部境界条件タイプが速度指定かどうかを調べる
 @retval VBCであればtrue
 */
bool CompoList::isVBC(void)
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

/**
 @fn bool CompoList::isVecForcing(void)
 @brief ベクトル強制をするかどうかを調べる
 @retval ベクトルを強制する場合true
 */
bool CompoList::isVecForcing(void)
{
  if ( isFORCING() ) {
    if ( usw==ON) return true;
  }
  return false;
}

//@fn bool CompoList::isVFraction(void)
//@brief 体積率の必要なコンポーネントかどうか
bool CompoList::isVFraction(void) {
  if ((type == HEAT_SRC) ||
      (type == CNST_TEMP) ||
      (type == IBM_DF) ||
      (type == HEX) ||
      (type == FAN) ||
      (type == DARCY) )  return true;
  return false;
}

//@fn void CompoList::setAttrb(const unsigned key)
//@brief attrbをセットする
//@param key アトリビュート
void CompoList::setAttrb(const unsigned key) { attrb = key; }

//@fn void CompoList::setBClocation(const unsigned key)
//@brief BCの方向をセットする
//@param key CDSの場合のBCの方向
void CompoList::setBClocation(const unsigned key) { bc_dir = key; }

//@fn void CompoList::setBbox(int m_st[], int m_ed[])
//@brief コンポーネントのBV情報を設定する
void CompoList::setBbox(const int m_st[], const int m_ed[]) {
  st[0] = m_st[0];
  st[1] = m_st[1];
  st[2] = m_st[2];
  ed[0] = m_ed[0];
  ed[1] = m_ed[1];
  ed[2] = m_ed[2];
}

//@fn void CompoList::setDef(const int key)
//@brief 指定セルを保持する
void CompoList::setDef(const int key) {
  def = key;
}

//@fn void CompoList::setElement(const unsigned long key)
//@brief elementをセットする
//@param key 要素数 element
void CompoList::setElement(const unsigned long key) { 
  element = key; 
}

//@fn void CompoList::void setEns(const bool key)
//@brief コンポーネントが自ノードに存在しているかどうかをセットする
void CompoList::setEns(const unsigned key) {
  ens = key;
}

//@fn void CompoList::setHtype(const int key)
//@brief h_typeをセットする
//@param key 境界条件の種類
void CompoList::setHtype(const unsigned key) { 
  h_type = key; 
}

//@fn void CompoList::setID(const unsigned key)
//@brief IDをセットする
//@param key 媒質ID ID
void CompoList::setID(const unsigned key) { 
  ID = key; 
}

//@fn void CompoList::setInitTemp(const REAL_TYPE key)
//@brief 初期温度の指定
void CompoList::setInitTemp(const REAL_TYPE key) {
  temp_init = key;
}

//@fn void CompoList::setMatOdr (const unsigned key)
//@brief mat_odrをセットする
//@param key MediumListのエントリ番号
void CompoList::setMatOdr(const unsigned key) { 
  mat_odr = key; 
}

/**
 @fn void CompoList::setName(const char* pnt)
 @brief ラベル名をセットする
 @param pnt ラベル名のアドレス
 @attention NULL check
 */
void CompoList::setName(const char* pnt) { 
  strcpy(name, pnt); 
}

//@fn void CompoList::setOutflowType(const unsigned key)
//@brief 流出速度のタイプを指定する
//@note V_AVERAGE | V_MINMAX
void CompoList::setOutflowType(const unsigned key) {
  var_u1 = key;
}

//@fn void CompoList::setPeriodicDir(const unsigned key)
//@brief 周期境界の上流方向を保持する
void CompoList::setPeriodicDir(const unsigned key) {
  var_u1 = key;
}

//@fn void CompoList::setPhase(const unsigned m_phase)
//@brief set pahse ID (SOLID=0, FLUID=1, GAS=2, LIQUID=3)
void CompoList::setPhase(const unsigned m_phase) {
  phase = m_phase;
}

//@fn void CompoList::setsetPrsUnit(const unsigned key)
//@brief 圧力の単位を指定する
void CompoList::setPrsUnit(const unsigned key) {
  var_u1 = key;
}

//@fn void CompoList::set_SamplingMethod(const int key)
void CompoList::set_SamplingMethod(const int key) { 
  sampling_method = key; 
}

//@fn void CompoList::set_SamplingMode(const int key)
void CompoList::set_SamplingMode(const int key) { 
  sampling_mode = key; 
}

//@fn void CompoList::set_Shape(const int key)
void CompoList::set_Shape(const int key) { 
  shape = key; 
}

//@fn void CompoList::setState  (const int key)
//@brief stateをセットする
//@param key セルの状態 SOLID/FLUID
void CompoList::setState(const int key) { 
  state = key; 
}

//@fn void CompoList::setStateCellMonitor(const unsigned key)
//@brief セルモニタスイッチ ON/OFF
void CompoList::setStateCellMonitor(const unsigned key) {
  var_u1 = key;
}

//@fn void CompoList::setType(const int key)
//@brief typeをセットする
//@param key 境界条件の種類
void CompoList::setType(const unsigned key) { 
  type = key; 
}


//@fn void CompoList::set_cmp_sz(void)
//@brief コンポーネントのサイズを計算
//@note ガイドセルは片側2層を仮定(c_sizeには含めない)
void CompoList::set_cmp_sz(void) {
  c_size[0] = ed[0] - st[0] + 1;
  c_size[1] = ed[1] - st[1] + 1;
  c_size[2] = ed[2] - st[2] + 1;
}


//@fn void CompoList::set_CoefHT(const REAL_TYPE var)
//@brief 熱伝達係数の保持
void CompoList::set_CoefHT(const REAL_TYPE var) {
  var2 = var;
}

//@fn void CompoList::set_CoefRadEps(const REAL_TYPE var)
//@brief 輻射のイプシロンの保持
void CompoList::set_CoefRadEps(const REAL_TYPE var) {
  var1 = var;
}

//@fn void CompoList::CoefRadPrj(const REAL_TYPE var)
//@brief 輻射の射出係数の保持
void CompoList::set_CoefRadPrj(const REAL_TYPE var) {
  var2 = var;
}

//@fn void CompoList::set_CoefMassflow(const REAL_TYPE var)
//@brief 流量の有次元化係数の保持
void CompoList::set_CoefMassflow(const REAL_TYPE var) {
  var1 = var;
}

//@fn void CompoList::CoefPrsLoss(const REAL_TYPE var)
//@brief 圧力損失の有次元化係数の保持
void CompoList::set_CoefPrsLoss(const REAL_TYPE var) {
  var2 = var;
}

//@fn void CompoList::set_HeatDensity(const REAL_TYPE var)
//@brief 吸発熱密度の保持
void CompoList::set_HeatDensity(const REAL_TYPE var) {
  var3 = var;
}

//@fn void CompoList::set_Heatflux(const REAL_TYPE var)
//@brief 熱流束の保持
void CompoList::set_Heatflux(const REAL_TYPE var) {
  var2 = var;
}

//@fn void CompoList::set_HeatValue(const REAL_TYPE var)
//@brief 吸発熱量の保持
void CompoList::set_HeatValue(const REAL_TYPE var) {
  var2 = var;
}

//@fn void CompoList::set_HSRC_policy(const bool kind)
//@brief 発熱項の指定ポリシーをセットする
//@param kind ポリシー種別　true-発熱量指定, false-発熱密度指定
void CompoList::set_HSRC_policy(const bool kind) {
  usw = ( kind ) ? hsrc_watt : hsrc_density;
}

//@fn void CompoList::set_Massflow(const REAL_TYPE var)
//@brief 流量の保持
void CompoList::set_Massflow(const REAL_TYPE var) {
  var1 = var;
}

//@fn void CompoList::set_Mon_Temp(const REAL_TYPE var)
//@brief モニタ温度の保持
void CompoList::set_Mon_Temp(const REAL_TYPE var) {
  var_m = var;
}

//@fn void CompoList::set_Mon_Heatflux(const REAL_TYPE var)
//@brief モニタ熱流束の保持
void CompoList::set_Mon_Heatflux(const REAL_TYPE var) {
  var_m = var;
}

//@fn void CompoList::set_Mon_Calorie(const REAL_TYPE var)
//@brief モニタ熱量の保持
void CompoList::set_Mon_Calorie(const REAL_TYPE var) {
  var_m = var;
}

//@fn void CompoList::set_Pressure(const REAL_TYPE var)
//@brief 圧力値の保持
void CompoList::set_Pressure(const REAL_TYPE var) {
  var1 = var;
}

//@fn void CompoList::set_sw_V_profile(const unsigned key)
//@brief 速度プロファイル指定モードの保持
void CompoList::set_sw_V_profile(const unsigned key) {
  usw = key;
}

//@fn void CompoList::set_sw_P_BCtype(const unsigned key)
//@brief 圧力境界条件タイプ指定モードの保持
void CompoList::set_sw_P_BCtype(const unsigned key) {
  usw = key;
}

//@fn void CompoList::set_sw_HTmodeRef(const unsigned key)
//@brief 熱伝達の参照指定モードの保持
void CompoList::set_sw_HTmodeRef(const unsigned key) {
  usw = key;
}

//@fn void CompoList::set_sw_HexDir (const unsigned key)
//@brief 熱交換機の方向指定モードの保持
void CompoList::set_sw_HexDir (const unsigned key) {
  usw = key;
}

//@fn void CompoList::set_sw_Heatgen(const unsigned key)
//@brief 発熱量指定モードの保持
void CompoList::set_sw_Heatgen(const unsigned key) {
  usw = key;
}

//@fn void CompoList::set_Temp(const REAL_TYPE var)
//@brief 温度の保持
void CompoList::set_Temp(const REAL_TYPE var) {
  var3 = var;
}

//@fn void CompoList::set_Velocity(const REAL_TYPE var)
//@brief 速度の保持
void CompoList::set_Velocity(const REAL_TYPE var) {
  var1 = var;
}

//@fn void CompoList::set_VBC_policy(const bool kind)
//@brief 速度指定ポリシーをセットする
//@param kind ポリシー種別　true-速度指定, false-流量指定
void CompoList::set_VBC_policy(const bool kind)
{
  attrb = ( kind ) ? BC_type_velocity : BC_type_massflow;
}
