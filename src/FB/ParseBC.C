//##################################################################################
//
// Flow Base class
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2014 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

//@file   ParseBC.C
//@brief  FlowBase ParseBC class
//@author aics

#include "ParseBC.h"
#include <math.h>


// #################################################################
// KOSと境界条件数の整合性をチェックする
void ParseBC::checkList(const MediumList* mat, const CompoList* cmp)
{
  printf("\n#############################\n");
  printf("Medium             Alias          State\n");
  for (int i=1; i<=NoCompo; i++)
  {
    printf("%6d\t%16s \t%7s\n", i, mat[i].getAlias().c_str(), (mat[i].getState()==FLUID)?"FLUID":"SOLID");
  }
  
  printf("\n\nCompo.             Alias          State       Medium\n");
  for (int i=1; i<=NoCompo; i++)
  {
    printf("%6d\t%16s \t%7s %12s\n", i, cmp[i].getAlias().c_str(), (cmp[i].getState()==FLUID)?"FLUID":"SOLID", cmp[i].getMedium().c_str());
  }
  printf("#############################\n\n");
}


// #################################################################
// KOSと境界条件数の整合性をチェックする 
void ParseBC::chkBCconsistency(const int kos, const CompoList* cmp)
{
  if (kos == FLOW_ONLY) 
  {
    for (int n=1; n<=NoCompo; n++)
    {
      if ( cmp[n].isHeatMode() )
      {
        Hostonly_ stamped_printf("\t'KindOfSolver' is FLOW_ONLY, but 'LocalBoundary' has heat boundary condition.\n");
        Exit(0);
      }
    }
  }
  else if (kos == SOLID_CONDUCTION)
  {
    for (int n=1; n<=NoCompo; n++)
    {
      if ( cmp[n].isVBC() ) 
      {
        Hostonly_ stamped_printf("\t'KindOfSolver' is SOLID_CONDUCTION, but 'LocalBoundary' has velocity boundary condition.\n");
        Exit(0);
      }
    }
  }
}


// #################################################################
/**
 * @brief ラベルの重複を調べる
 * @param [in] n       テストするBaseBcの格納番号の最大値
 * @param [in] m_label テストラベル
 */
bool ParseBC::chkDuplicate(const int n, const string m_label)
{
	for (int i=0; i<n; i++)
  {
    if ( BaseBc[i].getAlias() == m_label ) return false;
	}
	return true;
}

// #################################################################
/**
 * @brief 外部境界条件の候補を探し，内容をコピーする
 * @param [in]     str  string
 * @param [in]     face face number
 * @param [in,out] bc   array of bc
 */
bool ParseBC::findOBClist(const string str, const int face, BoundaryOuter* bc)
{
  for (int i=1; i<=NoBaseBC; i++)
  {
    if ( !strcasecmp( str.c_str(), BaseBc[i].getAlias().c_str() ) )
    {
      bc[face].dataCopy( &BaseBc[i] );
      break;
    }
    else
    {
      if ( i == NoBaseBC ) // 最後までみつからない
      {
        return false;
      }
    }
  }
  return true;
}

// #################################################################
// KOSと媒質の状態の整合性をチェックし，媒質数をカウント，C.NoMediumFluid, C.NoMediumSolidをセット
void ParseBC::countMedium(Control* Cref, const MediumList* mat)
{
  // check at least one fluid
  if ( KindOfSolver != SOLID_CONDUCTION )
  {
    bool check=false;
    for (int i=1; i<=NoMedium; i++)
    {
      if ( mat[i].getState() == FLUID ) check = true;
    }
    if ( !check )
    {
      Hostonly_ stamped_printf("\tAnalysis model should have at least one FLUID medium in MediumTable.\n");
      Exit(0);
    }
  }
  
  // 流体と固体の媒質数をセット
  int m_fluid=0, m_solid=0;
  
  for (int i=1; i<=NoMedium; i++)
  {
    if      ( mat[i].getState() == SOLID ) m_solid++;
    else if ( mat[i].getState() == FLUID ) m_fluid++;
  }
  
  Cref->NoMediumFluid = m_fluid;
  Cref->NoMediumSolid = m_solid;
}




// #################################################################
/**
 * @brief ConstTemperatureのパラメータを取得する
 * @param [in]  label_base ラベルディレクトリ
 * @param [in]  n          コンポーネントリストの格納番号
 * @param [out] cmp        CompoList
 */
void ParseBC::getIbcCnstTemp(const string label_base, const int n, CompoList* cmp)
{
  string label = label_base + "/Temperature";
  
  REAL_TYPE tmp = getValueReal(label);
  cmp[n].setTemp( tmp );
}



// #################################################################
/**
 * @brief Direct_Fluxのパラメータを取得する
 * @param [in]  label_base ラベルディレクトリ
 * @param [in]  n          コンポーネントリストの格納番号
 * @param [out] cmp        CompoList
 */
void ParseBC::getIbcHeatFlux(const string label_base, const int n, CompoList* cmp)
{
  string label = label_base + "/HeatFlux";
  
  cmp[n].setHeatflux( getValueReal(label) ); /// @note [W/m^2]
}


// #################################################################
/**
 * @brief HeatGenerationのパラメータを取得する
 * @param [in]  label_base ラベルディレクトリ
 * @param [in]  n          コンポーネントリストの格納番号
 * @param [out] cmp        CompoList
 */
void ParseBC::getIbcHeatSrc(const string label_base, const int n, CompoList* cmp)
{
  REAL_TYPE ct=0.0f;
  string str;
  string label;
  
  
  // 発熱量は発熱密度の形式で保持
  label = label_base + "/HeatReleaseValue";
  
  if ( tpCntl->chkLabel(label) ) // 発熱量 [W]
  {
    if ( tpCntl->getInspectedValue(label, ct) )
    {
      cmp[n].setHsrcPolicy(CompoList::hsrc_watt);
      cmp[n].setHeatReleaseValue(ct);
    }
    else
    {
      stamped_printf("\tParsing error : Invalid value for '%s\n", label.c_str());
      Exit(0);
    }
  }
  else // 発熱密度 [W/m^3]
  {
    label = label_base + "/HeatGenerationDensity";
    
    if ( tpCntl->chkLabel(label) )
    {
      if ( tpCntl->getInspectedValue(label, ct) )
      {
        cmp[n].setHsrcPolicy(CompoList::hsrc_density);
        cmp[n].setHeatGenerationDensity(ct);
      }
      else
      {
        stamped_printf("\tParsing error : Invalid value for '%s\n", label.c_str());
        Exit(0);
      }
    }
    else
    {
      stamped_printf("\tParsing error : Invalid STRING value for '%s\n", label.c_str());
      Exit(0);
    }
  }
  
}


// #################################################################
/**
 * @brief HeatTransferSのパラメータを取得
 * @param [in]  label_base ラベルディレクトリ
 * @param [in]  n          コンポーネントリストの格納番号
 * @param [out] cmp        CompoList
 */
void ParseBC::getIbcHT_S(const string label_base, const int n, CompoList* cmp)
{
  string label;
  
  // 熱伝達係数
  label = label_base + "/CoefOfHeatTransfer";
  
  cmp[n].setCoefHT( getValueReal(label) );
  
  // 表面温度
  label = label_base + "/BulkTemperature";
  
  REAL_TYPE st = getValueReal(label);
  cmp[n].setTemp( st );

}


// #################################################################
/**
 * @brief HeatTransferSFのパラメータを取得する
 * @param [in]  label_base ラベルディレクトリ
 * @param [in]  n          コンポーネントリストの格納番号
 * @param [out] cmp        CompoList
 */
void ParseBC::getIbcHT_SF(const string label_base, const int n, CompoList* cmp)
{
  string str;
  string label;
  
  // 表面温度
  label = label_base+"/BulkTemperature";
  REAL_TYPE st = getValueReal(label);
  cmp[n].setTemp( st );
  
  // coefficients
  label = label_base+"/Alpha";
  cmp[n].ca[CompoList::alpha] = getValueReal(label);
  
  label = label_base+"/Beta";
  cmp[n].ca[CompoList::beta]  = getValueReal(label);
  
  label = label_base+"/Gamma";
  cmp[n].ca[CompoList::gamma] = getValueReal(label);
}



// #################################################################
/**
 * @brief HeatTransferSNのパラメータを取得する
 * @param [in]  label_base ラベルディレクトリ
 * @param [in]  n          コンポーネントリストの格納番号
 * @param [out] cmp        CompoList
 */
void ParseBC::getIbcHT_SN(const string label_base, const int n, CompoList* cmp)
{
  string str;
  string label;
  
  // 表面温度
  label = label_base + "/BulkTemperature";
  
  REAL_TYPE st = getValueReal(label);
  cmp[n].setTemp( st );
  
  // Vertical and upper face values
  label = label_base + "/VerticalLaminarAlpha";
  cmp[n].ca[CompoList::vert_laminar_alpha]    = getValueReal(label);
  
  label = label_base + "/VerticalLaminarBeta";
  cmp[n].ca[CompoList::vert_laminar_beta]     = getValueReal(label);
  
  label = label_base + "/VerticalTurbulentAlpha";
  cmp[n].ca[CompoList::vert_turbulent_alpha]  = getValueReal(label);
  
  label = label_base + "/VerticalTurbulentBeta";
  cmp[n].ca[CompoList::vert_turbulent_beta]   = getValueReal(label);
  
  label = label_base + "/VerticalRaCritial";
  cmp[n].ca[CompoList::vert_Ra_critial]       = getValueReal(label);
  
  
  // Lower face values
  label = label_base + "/LowerLaminarAlpha";
  cmp[n].cb[CompoList::lower_laminar_alpha]   = getValueReal(label);
  
  label = label_base + "/LowerLaminarBeta";
  cmp[n].cb[CompoList::lower_laminar_beta]    = getValueReal(label);
  
  label = label_base + "/LowerTurbulentAlpha";
  cmp[n].cb[CompoList::lower_turbulent_alpha] = getValueReal(label);
  
  label = label_base + "/LowerTurbulentBeta";
  cmp[n].cb[CompoList::lower_turbulent_beta]  = getValueReal(label);
  
  label = label_base + "/LowerRaCritial";
  cmp[n].cb[CompoList::lower_Ra_critial]      = getValueReal(label);
}



// #################################################################
/**
 * @brief 境界条件IsoThermalのパラメータを取得し保持する
 * @param [in]  label_base ラベルディレクトリ
 * @param [in]  n          コンポーネントリストの格納番号
 * @param [out] cmp        CompoList
 */
void ParseBC::getIbcIsoTherm(const string label_base, const int n, CompoList* cmp)
{
  string label;
  
  // 表面温度
  label = label_base+"/Temperature";
  REAL_TYPE tmp = getValueReal(label);
  cmp[n].setTemp( tmp );
}


// #################################################################
/**
 * @brief 内部の流出境界のパラメータを取得する
 * @param [in]  label_base ラベルディレクトリ
 * @param [in]  n          コンポーネントリストの格納番号
 * @param [out] cmp        CompoList
 */
void ParseBC::getIbcOutflow(const string label_base, const int n, CompoList* cmp)
{
  string str;
  string label;
  REAL_TYPE dv[3];
  
  // 圧力境界のタイプ default
  cmp[n].set_P_BCtype( P_GRAD_ZERO );
  
  // Hidden parameter
  label = label_base + "/PressureType";
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
    ; // なくてもOK
  }
  else {
    if ( !strcasecmp("Dirichlet", str.c_str()) )
    {
      cmp[n].set_P_BCtype( P_DIRICHLET );
    }
    else if ( !strcasecmp("GradZero", str.c_str()) )
    {
      cmp[n].set_P_BCtype( P_GRAD_ZERO );
    }
    else
    {
      printf("\tParsing error : Invalid string value for 'PressureType' : %s\n", str.c_str());
      Exit(0);
    }
    
    if ( cmp[n].get_P_BCtype() == P_DIRICHLET )
    {
      label = label_base + "/PressureValue";
      REAL_TYPE tmp=0.0;
      if ( !(tpCntl->getInspectedValue(label, tmp)) )
      {
        printf("\tParsing error : Missing for 'PressureValue'\n");
        Exit(0);
      }
      cmp[n].set_Pressure(tmp);
    }
  }
  
  // 法線ベクトルの取得
  label = label_base + "/OrientationVector";
  if ( !Control::getVec(label, dv, tpCntl, true) ) Exit(0);
  cmp[n].nv[0] = dv[0];
  cmp[n].nv[1] = dv[1];
  cmp[n].nv[2] = dv[2];
  
  // 出口方向を指定
  cmp[n].setBClocation(CompoList::opposite_direction);
  
}



// #################################################################
/**
 * @brief 内部の周期境界のパラメータを取得する
 * @param [in]  label_base ラベルディレクトリ
 * @param [in]  n          コンポーネントリストの格納番号
 * @param [out] cmp        CompoList
 */
void ParseBC::getIbcPeriodic(const string label_base, const int n, CompoList* cmp)
{
  int dir=0;
  REAL_TYPE ct=0.0;
  string str;
  string label;
  
  // 上流側の方向
  label = label_base + "/UpstreamDirection";
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
	  printf("\tParsing error : fail to get 'upstream_direction' in 'LocalBoundary > Periodic'\n");
	  Exit(0);
  }
  if ( !strcasecmp("xminus", str.c_str()) )
  {
		dir = X_minus;
	}
  else if ( !strcasecmp("xplus", str.c_str()) )
  {
    dir = X_plus;
  }
  else if ( !strcasecmp("yminus", str.c_str()) )
  {
    dir = Y_minus;
  }
  else if ( !strcasecmp("yplus", str.c_str()) )
  {
    dir = Y_plus;
  }
  else if ( !strcasecmp("zminus", str.c_str()) )
  {
    dir = Z_minus;
  }
  else if ( !strcasecmp("zplus", str.c_str()) )
  {
    dir = Z_plus;
  }
  else
  {
    stamped_printf("\tParsing error : Invalid direction in 'LocalBoundary > Periodic'\n");
    Exit(0);
  }
	cmp[n].setPeriodicDir((int)dir);
  
  
  // 圧力差
  label = label_base + "/PressureDifference";
  
  if ( !(tpCntl->getInspectedValue(label, ct )) )
  {
    printf("\tParsing error : Invalid value of 'Pressure difference' in 'LocalBoundary > Periodic'\n");
    Exit(0);
  }
  else
  {
    cmp[n].ca[0] = ct;
  }
}


// #################################################################
/**
 * @brief 回転体のパラメータを取得する
 * @param [in]  label_base ラベルディレクトリ
 * @param [in]  n          コンポーネントリストの格納番号
 * @param [out] cmp        CompoList
 */
void ParseBC::getIbcSolidRev(const string label_base, const int n, CompoList* cmp)
{
  string str;
  string label;
  REAL_TYPE dv[3], ct;
  
  
  // 法線ベクトルの取得 trueで単位ベクトル化
  label = label_base + "/OrientationVector";
  if ( !Control::getVec(label, dv, tpCntl, true) ) Exit(0);
  cmp[n].nv[0] = dv[0];
  cmp[n].nv[1] = dv[1];
  cmp[n].nv[2] = dv[2];
  
  
  // 中心座標の取得
  label = label_base + "/Center";
  if ( !Control::getVec(label, dv, tpCntl, false) ) Exit(0);
  if ( Unit_Param == DIMENSIONAL )
  {
    cmp[n].oc[0] = dv[0];
    cmp[n].oc[1] = dv[1];
    cmp[n].oc[2] = dv[2];
  }
  else
  {
    cmp[n].oc[0] = dv[0] * RefLength;
    cmp[n].oc[1] = dv[1] * RefLength;
    cmp[n].oc[2] = dv[2] * RefLength;
  }
  
  
  // 形状パラメータ
  label = label_base + "/Depth";
  ct = getValueReal(label);
  if ( Unit_Param == DIMENSIONAL )
  {
    cmp[n].depth = ct;
  }
  else
  {
    cmp[n].depth = ct * RefLength;
  }
  
  
  
  label = label_base + "/Radius";
  ct = getValueReal(label);
  if ( Unit_Param == DIMENSIONAL )
  {
    cmp[n].shp_p1 = ct;
  }
  else
  {
    cmp[n].shp_p1 = ct * RefLength;
  }
  
  
  // 回転数
  label = label_base + "/RotationFrequency";
  if ( !(tpCntl->getInspectedValue(label, ct)) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  REAL_TYPE pi = 2.0*asin(1.0);
  REAL_TYPE omg = 2.0*pi*ct/60.0; // [rad/s]
  
  cmp[n].ca[0] = ( Unit_Param == DIMENSIONAL ) ? omg : ct;
  
}




// #################################################################
/**
 * @brief 内部の流入境界のパラメータを取得する
 * @param [in]  label_base ラベルディレクトリ
 * @param [in]  n          コンポーネントリストの格納番号
 * @param [out] cmp        CompoList
 * @note Control::setparameters()でcmp[].ca[]に値をセットする
 */
void ParseBC::getIbcSpecVel(const string label_base, const int n, CompoList* cmp)
{
  string str;
  string label;
  REAL_TYPE ct=0;
  REAL_TYPE dv[3];
  
  // 指定タイプの特定
  label = label_base + "/Type";
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : Invalid SpecifiedType in '%s'\n", label.c_str());
    Exit(0);
  }
  if ( !strcasecmp("Velocity", str.c_str()) )
  {
    cmp[n].set_VBC_policy(true);
  }
  else if ( !strcasecmp("Massflow", str.c_str()) )
  {
    cmp[n].set_VBC_policy(false);
  }
  else
  {
    Hostonly_ printf("\tParsing error : Invalid string value '%s' for 'Type'\n", str.c_str());
    Exit(0);
  }
  
  // 速度指定タイプ
  cmp[n].set_V_profile( getVprofile(label_base) );
  
  
  // 速度パラメータの読み込み
  getVelocity(label_base, cmp[n].get_V_Profile(), cmp[n].ca, "LocalBoundary", cmp[n].isPolicy_Massflow());
  
  
  // 法線ベクトル
  label = label_base + "/OrientationVector";
  if ( !Control::getVec(label, dv, tpCntl, true) ) Exit(0);
  cmp[n].nv[0] = dv[0];
  cmp[n].nv[1] = dv[1];
  cmp[n].nv[2] = dv[2];
  
  
  // 境界条件の方向
  label = label_base + "/InOut";
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
    stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  if ( !strcasecmp("in", str.c_str()) )
  {
    cmp[n].setBClocation(CompoList::same_direction);
  }
  else if ( !strcasecmp("out", str.c_str()) )
  {
    cmp[n].setBClocation(CompoList::opposite_direction);
  }
  else
  {
    printf("\tParsing error : Invalid string value '%s' for 'InOut'\n", str.c_str());
    Exit(0);
  }
  
  
  // heat problem
  if ( HeatProblem )
  {
    label = label_base + "/Temperature";
    
    if ( !(tpCntl->getInspectedValue(label, ct )) )
    {
      Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    
    cmp[n].setTemp( ct );
    
    if ( Unit_Param != DIMENSIONAL )
    {
      Hostonly_ stamped_printf("\tWarning: Heat condition must be a dimensional value\n");
      Exit(0);
    }
    
  }
  
}



// #################################################################
/**
 * @brief
 * @param [in]  label_base ラベルディレクトリ
 * @param [in]  n          コンポーネントリストの格納番号
 * @param [out] cmp        CompoList
 */
void ParseBC::get_Darcy(const string label_base, const int n, CompoList* cmp)
{
  REAL_TYPE v[3];
  string label;
  
  for (int n=0; n<3; n++) v[n]=0.0;
  
  
  // 透過率の取得
  label = label_base + "/Permeability";
  
  if( !(tpCntl->getInspectedVector(label, v, 3)) ) {
    stamped_printf("\tParsing error : fail to get permeability params in 'Darcy'\n");
    Exit(0);
  }
  cmp[n].ca[0] = v[0]; // 透過率[m^2]は境界条件設定時に無次元化する
  cmp[n].ca[1] = v[1];
  cmp[n].ca[2] = v[2];
}



// #################################################################
/**
 * @brief Fanのパラメータを取得する
 * @param [in]  label_base ラベルディレクトリ
 * @param [in]  n          コンポーネントリストの格納番号
 * @param [out] cmp        CompoList
 */
void ParseBC::get_IBC_Fan(const string label_base, const int n, CompoList* cmp)
{
  string str,str_u;
  string label;
  REAL_TYPE dv[3];
  
  // 入力単位の指定
  label=label_base+"/Unit";//
  //cout <<  "label : " << label << endl;
  if ( !(tpCntl->getInspectedValue(label, str_u )) ) {
    stamped_printf("\tParsing error : Invalid float value for 'unit' in 'Pressure_Loss'\n");
    Exit(0);
  }
  if ( !strcasecmp("mmaq", str_u.c_str()) ) {
    cmp[n].setPrsUnit(CompoList::unit_mmAq);
  }
  else if ( !strcasecmp("mmhg", str_u.c_str()) ) {
    cmp[n].setPrsUnit(CompoList::unit_mmHg);
  }
  else if ( !strcasecmp("pa", str_u.c_str()) ) {
    cmp[n].setPrsUnit(CompoList::unit_Pa);
  }
  else if ( !strcasecmp("nondimension", str_u.c_str()) ) {
    cmp[n].setPrsUnit(CompoList::unit_NonDimensional);
  }
  else {
    stamped_printf("\tDescribed unit is out of scope.\n");
    Exit(0);
  }
  
  // 法線ベクトルの取得
  label = label_base + "/OrientationVector";
  if ( !Control::getVec(label, dv, tpCntl, true) ) Exit(0);
  cmp[n].nv[0] = dv[0];
  cmp[n].nv[1] = dv[1];
  cmp[n].nv[2] = dv[2];
  
  // 中心座標の取得
  label = label_base + "/Center";
  if ( !Control::getVec(label, dv, tpCntl, false) ) Exit(0);
  cmp[n].oc[0] = dv[0];
  cmp[n].oc[1] = dv[1];
  cmp[n].oc[2] = dv[2];
  
  // 形状パラメータ
  label=label_base+"/Depth";
  cmp[n].depth  = getValueReal(label);
  
  label=label_base+"/FanRadius";
  cmp[n].shp_p1 = getValueReal(label);
  
  label=label_base+"/BossRadius";
  cmp[n].shp_p2 = getValueReal(label);
  
  if ( cmp[n].shp_p1 <= cmp[n].shp_p2 )
  {
    stamped_printf("\tError : Radius of boss is greater than fan.\n");
    Exit(0);
  }
  
}



// #################################################################
/**
 * @brief Direct Forcingのパラメータを取得する
 * @param [in]  label_base ラベルディレクトリ
 * @param [in]  n          コンポーネントリストの格納番号
 * @param [out] cmp        CompoList
 */
void ParseBC::get_IBC_IBM_DF(const string label_base, const int n, CompoList* cmp)
{
  string str;
  string label;
  REAL_TYPE dv[3];
  
  // 法線ベクトル
  label = label_base + "/OrientationVector";
  if ( !Control::getVec(label, dv, tpCntl, true) ) Exit(0);
  cmp[n].nv[0] = dv[0];
  cmp[n].nv[1] = dv[1];
  cmp[n].nv[2] = dv[2];
  
  // Velocity
  label=label_base+"/Velocity";
  
  REAL_TYPE ct = getValueReal(label);
  
  if ( Unit_Param == DIMENSIONAL )
  {
    cmp[n].set_Velocity( ct );
  }
  else {
    cmp[n].set_Velocity( ct * RefVelocity );
  }
}




// #################################################################
/**
 * @brief HeatExchangerのパラメータを取得する
 * @param [in]  label_base ラベルディレクトリ
 * @param [in]  n          コンポーネントリストの格納番号
 * @param [out] cmp        CompoList
 * @note この時点ではRefDensityの値が未定なので，あとでパラメータ処理
 * @see Control::setParameters()
 */
void ParseBC::get_IBC_PrsLoss(const string label_base, const int n, CompoList* cmp)
{
  string str,str_u;
  string label;
  REAL_TYPE v[4];
  REAL_TYPE dv[3];
  
  // 入力単位の指定
  label = label_base + "/Unit";
  
  if ( !(tpCntl->getInspectedValue(label, str_u )) ) {
		stamped_printf("\tParsing error : Invalid float value for 'Unit' in 'PressureLoss'\n");
		Exit(0);
  }
  if ( !strcasecmp("mmaq", str_u.c_str()) ) {
    cmp[n].setPrsUnit(CompoList::unit_mmAq);
  }
  else if ( !strcasecmp("mmhg", str_u.c_str()) ) {
    cmp[n].setPrsUnit(CompoList::unit_mmHg);
  }
  else if ( !strcasecmp("pa", str_u.c_str()) ) {
    cmp[n].setPrsUnit(CompoList::unit_Pa);
  }
  else if ( !strcasecmp("NonDimension", str_u.c_str()) ) {
    cmp[n].setPrsUnit(CompoList::unit_NonDimensional);
  }
  else {
    stamped_printf("\tDescribed unit is out of scope.\n");
    Exit(0);
  }
  
  // 法線ベクトルの取得
  label = label_base + "/OrientationVector";
  if ( !Control::getVec(label, dv, tpCntl, true) ) Exit(0);
  cmp[n].nv[0] = dv[0];
  cmp[n].nv[1] = dv[1];
  cmp[n].nv[2] = dv[2];
  
  // 方向ベクトルの取得
  label = label_base + "/Dir";
  if ( !Control::getVec(label, dv, tpCntl, true) ) Exit(0);
  cmp[n].dr[0] = dv[0];
  cmp[n].dr[1] = dv[1];
  cmp[n].dr[2] = dv[2];
  
  // 中心座標の取得
  label = label_base + "/Center";
  if ( !Control::getVec(label, dv, tpCntl, false) ) Exit(0);
  cmp[n].oc[0] = dv[0];
  cmp[n].oc[1] = dv[1];
  cmp[n].oc[2] = dv[2];
  
  // 形状パラメータ
  label = label_base + "/Depth";
  cmp[n].depth  = getValueReal(label);
  
  label = label_base + "/Width";
  cmp[n].shp_p1 = getValueReal(label);
  
  label = label_base + "/Height";
  cmp[n].shp_p2 = getValueReal(label);
  
  // 圧力損失パラメータ
  label = label_base + "/C";
  if( !(tpCntl->getInspectedVector(label, v, 4)) )
  {
    stamped_printf("\tParsing error : fail to get vec params in '%s\n", label.c_str());
    Exit(0);
  }
  cmp[n].ca[0]=v[0];
  cmp[n].ca[1]=v[1];
  cmp[n].ca[2]=v[2];
  cmp[n].ca[3]=v[3];
  
  label = label_base + "/Uthreshold";
  cmp[n].ca[4] = getValueReal(label);
  
  label = label_base + "/Thickness";
  cmp[n].ca[5] = getValueReal(label);
  
  // 熱交換器の方向強制オプション
  
  label = label_base + "/Vector";
  if ( !(tpCntl->getInspectedValue(label, str )) ) {
    stamped_printf("\tParsing error : Invalid string for 'vector' in 'PressureLoss'\n");
    Exit(0);
  }
  if ( !strcasecmp("directional", str.c_str()) ) {
    cmp[n].set_sw_HexDir( ON );
  }
  else {
    cmp[n].set_sw_HexDir( OFF );
  }
  
}



// #################################################################
/** 境界条件Radiantのパラメータを取得し保持する
 * @brief
 * @param [in]  label_base ラベルディレクトリ
 * @param [in]  n          コンポーネントリストの格納番号
 * @param [out] cmp        CompoList
 * @note 境界条件自体は未実装
 */
void ParseBC::get_IBC_Radiant(const string label_base, const int n, CompoList* cmp)
{
  string label;
  
  // 係数
  label=label_base+"/Epsilon";
  cmp[n].set_CoefRadEps( getValueReal(label) );
  
  // 射出率
  label=label_base+"/Projection";
  cmp[n].set_CoefRadPrj( getValueReal(label) );
}




// #################################################################
// 温度計算の場合の各媒質の初期値を取得する
void ParseBC::getInitTempOfMedium(CompoList* cmp, Control* C)
{  
  string label, label_base;
  string str;
  int no_list=0;
  
  label_base = "/StartCondition/InitialState/MediumTemperature";
  
  if ( !(tpCntl->chkNode(label_base)) )
  {
    Hostonly_ stamped_printf("\tParsing error : fail to get '%s'\n", label_base.c_str());
    Exit(0);
  }
  else // label_base直下のノード数を得る
  {
	  no_list = tpCntl->countLabels(label_base);
  }

  if ( no_list != NoMedium )
  {
    Hostonly_ stamped_printf("\tNo of medium [%d] listed in '%s' is not agree with one [%d] in MediumTable.\n", no_list, label_base.c_str(), NoMedium);
    Exit(0);
  }
  
  
  for (int i=1; i<=no_list; i++)
  {
    if ( !(tpCntl->getNodeStr(label_base, i, str)) )
    {
      Hostonly_ stamped_printf("\tGetNodeStr error\n");
      Exit(0);
    }
    
    // medium name
    label = label_base + "/" + str;
    REAL_TYPE ct=0.0;
    
    if ( !(tpCntl->getInspectedValue(label, ct)) )
    {
      Hostonly_ stamped_printf("\tParsing error : Invalid keyword '%s'\n", label.c_str());
      Exit(0);
    }
    
    int flag = 0;
    
    for (int m=1; m<=NoCompo; m++)
    {
      if ( FBUtility::compare(cmp[m].getAlias(), str) )
      {
        if ( cmp[m].isKindMedium())
        {
          cmp[m].setInitTemp( ct );
          flag++;
        }
        else
        {
          Hostonly_ stamped_printf("\tError : Label '%s' is not a medium\n", cmp[m].getAlias().c_str());
          Exit(0);
        }
      }
    }
    
    if ( flag == 0 )
    {
      Hostonly_ stamped_printf("\tParsing error : Invalid keyword in '%s'\n", label_base.c_str());
      Exit(0);
    }
    
  }
}



// #################################################################
/**
 * @brief 外部境界条件のキーワードを照合し， BCの文字列を返す
 * @param [in] id
 */
string ParseBC::getOBCstr(const int id)
{
  string bc;
  if     ( id == OBC_WALL )      bc = "Wall";
  else if( id == OBC_OUTFLOW )   bc = "Outflow";
  else if( id == OBC_SPEC_VEL )  bc = "SpecifiedVelocity";
  else if( id == OBC_SYMMETRIC ) bc = "Symmetric";
  else if( id == OBC_PERIODIC )  bc = "Periodic";
  else if( id == OBC_TRC_FREE )  bc = "TractionFree";
  else if( id == OBC_FAR_FIELD ) bc = "FarField";
  else if( id == OBC_INTRINSIC ) bc = "Intrinsic";
  else                           bc = "";
  return bc;
}



// #################################################################
/**
 * @brief 外部境界の遠方境界のパラメータを取得する
 * @param [in] label_base ラベルディレクトリ
 * @param [in] n          面番号
 */
void ParseBC::getObcFarField(const string label_base, const int n)
{
  string str;
  string label;
  
  BaseBc[n].set_pType(P_GRAD_ZERO);
  BaseBc[n].p = 0.0; // ダミー値
  
  // 外部雰囲気温
  if ( HeatProblem )
  {
    if ( Unit_Param != DIMENSIONAL )
    {
      stamped_printf("\tError: Heat condition must be given by dimensional value\n");
      Exit(0);
    }
    
    label = label_base + "/AmbientTemperature";
    REAL_TYPE ct=0.0;
    
    if ( !(tpCntl->getInspectedValue(label, ct )) )
    {
      stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }

    BaseBc[n].setTemp( ct );
  }
  
  // 圧力境界のタイプ  default
  BaseBc[n].set_pType(P_GRAD_ZERO);
  BaseBc[n].p = 0.0; // ダミー値
  
  // Hidden option
  label = label_base + "/PressureType";
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
    ; // エラーではない
  }
  else
  {
    if ( !strcasecmp("dirichlet", str.c_str()) )
    {
      BaseBc[n].set_pType(P_DIRICHLET);
    }
    else if ( !strcasecmp("gradzero", str.c_str()) )
    {
      BaseBc[n].set_pType(P_GRAD_ZERO);
    }
    else
    {
      stamped_printf("\tParsing error : Invalid string value for 'PressureType' : %s\n", str.c_str());
      Exit(0);
    }
  }
  
  
  // 圧力の値
  if ( BaseBc[n].get_pType() == P_DIRICHLET )
  {
    if ( !strcasecmp("dirichlet", str.c_str()) )
    {
      label = label_base + "/PrsValue";
      REAL_TYPE ct=0.0;
      
      if ( !(tpCntl->getInspectedValue(label, ct )) )
      {
        stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
        Exit(0);
      }
      else
      {
        BaseBc[n].p = ct;
      }
    }
  }
  
}


// #################################################################
/**
 * @brief 外部の壁面熱伝達境界のパラメータを取得する
 * @param [in] label_base ラベルディレクトリ
 * @param [in] n          面番号
 * @param [in] kind       熱伝達境界の種類
 */
void ParseBC::getObcHeatTransfer(const string label_base, const int n, const string kind)
{
  string label;
  string str;
  REAL_TYPE ct;
  

  if ( !strcasecmp(kind.c_str(), "HeatTransferS") )
  {
    BaseBc[n].set_HTmode(HT_S);
    label = label_base + "/CoefOfHeatTransfer";
    BaseBc[n].setCoefHT( getValueReal(label) );
    
    label = label_base + "/SurfaceTemperature";
    ct = getValueReal(label);
    BaseBc[n].setTemp( ct );
  }
  else if ( !strcasecmp(kind.c_str(), "HeatTransferSF") )
  {
    BaseBc[n].set_HTmode(HT_SF);
    label = label_base + "/SurfaceTemperature";
    BaseBc[n].setTemp( getValueReal(label) );
    
    
    label=label_base+"/RefTempMode";
    BaseBc[n].set_HTmodeRef( CompoList::HT_mode_bulk );
    if ( tpCntl->getInspectedValue(label, str ) )
    {
      if ( !strcasecmp("Bulk", str.c_str()) )
      {
        BaseBc[n].set_HTmodeRef( CompoList::HT_mode_bulk );
      }
      else if ( !strcasecmp("Local", str.c_str()) )
      {
        BaseBc[n].set_HTmodeRef( CompoList::HT_mode_local );
      }
    }

    // coefficients
    label = label_base + "/Alpha";
    BaseBc[n].ca[0] = getValueReal(label);
    
    label = label_base + "/Beta";
    BaseBc[n].ca[1] = getValueReal(label);
    
    label = label_base + "/Gamma";
    BaseBc[n].ca[2] = getValueReal(label);
  }
  else if ( !strcasecmp(kind.c_str(), "HeatTransferSN") )
  {
    BaseBc[n].set_HTmode(HT_SN);
    label = label_base + "/SurfaceTemperature";
    BaseBc[n].setTemp( getValueReal(label) );
    
    // reference mode
    label = label_base + "/RefTempMode";
    BaseBc[n].set_HTmodeRef( CompoList::HT_mode_bulk );
    if ( tpCntl->getInspectedValue(label, str ) )
    {
      if ( !strcasecmp("Bulk", str.c_str()) )
      {
        BaseBc[n].set_HTmodeRef( CompoList::HT_mode_bulk );
      }
      else if ( !strcasecmp("Local", str.c_str()) )
      {
        BaseBc[n].set_HTmodeRef( CompoList::HT_mode_local );
      }
    }
    
    // Vertical and upper face values
    label = label_base + "/VertivalLaminarAlpha";
    BaseBc[n].ca[0] = getValueReal(label);
    
    label = label_base + "/VertivalLaminarBeta";
    BaseBc[n].ca[1] = getValueReal(label);
    
    label = label_base + "/vertivalTurbulentAlpha";
    BaseBc[n].ca[2] = getValueReal(label);
    
    label = label_base + "/VertivalTurbulentBeta";
    BaseBc[n].ca[3] = getValueReal(label);
    
    label = label_base + "/VertivalRaCritial";
    BaseBc[n].ca[4] = getValueReal(label);
    
    // Lower face values
    label = label_base + "/LowerLaminarAlpha";
    BaseBc[n].cb[0] = getValueReal(label);
    
    label = label_base + "/LowerLaminarBeta";
    BaseBc[n].cb[1] = getValueReal(label);
    
    label = label_base + "/LowerTurbulentAlpha";
    BaseBc[n].cb[2] = getValueReal(label);
    
    label = label_base + "/LowerTurbulentBeta";
    BaseBc[n].cb[3] = getValueReal(label);
    
    label = label_base + "/LowerRaCritial";
    BaseBc[n].cb[4] = getValueReal(label);
  }
  else {
    stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
}



// #################################################################
/**
 * @brief 外部境界の流出条件のパラメータを取得する
 * @param [in] label_base ラベルディレクトリ
 * @param [in] n          面番号
 * @note 圧力の値は，Control::setParameters()で無次元化する
 */
void ParseBC::getObcOutflow(const string label_base, const int n)
{
  string str;
  string label;
  REAL_TYPE ct;
  
  
  // 圧力境界のタイプ  default
  BaseBc[n].set_pType(P_GRAD_ZERO);
  BaseBc[n].p = 0.0; // ダミー値
  
  label = label_base + "/PressureType";
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
    printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  else
  {
    if ( !strcasecmp("dirichlet", str.c_str()) )
    {
      BaseBc[n].set_pType(P_DIRICHLET);
    }
    else if ( !strcasecmp("neumann", str.c_str()) )
    {
      BaseBc[n].set_pType(P_GRAD_ZERO);
    }
    else
    {
      stamped_printf("\tParsing error : Invalid string value for 'PressureType' : %s\n", str.c_str());
      Exit(0);
    }
  }
  
  
  // 値
  label = label_base + "/PrsValue";
  
  if ( !(tpCntl->getInspectedValue(label, ct )) )
  {
    stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  else
  {
    BaseBc[n].p = ct;
  }
}


// #################################################################
/**
 * @brief 外部境界の周期条件のパラメータを取得する
 * @param [in] label_base ラベル
 * @param [in] n          面番号
 * @note 圧力の値は，Control::setParameters()で無次元化する
 */
void ParseBC::getObcPeriodic(const string label_base, const int n)
{
  REAL_TYPE ct;
  int def;
  string str;
  string label;
  
  // モード
  label = label_base + "/Mode";
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
    printf("\tParsing error : fail to get '%s'\n", label.c_str());
    Exit(0);
  }
  else
  {
    if ( !strcasecmp(str.c_str(), "SimpleCopy") )
    {
      BaseBc[n].set_PrdcMode(BoundaryOuter::prdc_Simple);
    }
    else if ( !strcasecmp(str.c_str(), "directional") )
    {
      BaseBc[n].set_PrdcMode(BoundaryOuter::prdc_Directional);
    }
    else if ( !strcasecmp(str.c_str(), "driver") )
    {
      BaseBc[n].set_PrdcMode(BoundaryOuter::prdc_Driver);
    }
  }
  
  // Directional
  if ( BaseBc[n].getPrdcMode() == BoundaryOuter::prdc_Directional )
  {
    label = label_base + "/PressureDifference";
    if ( !(tpCntl->getInspectedValue(label, ct )) )
    {
      printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    else
    {
      BaseBc[n].p = ct;
    }
    
    label = label_base + "/FlowDirection";
    if ( !(tpCntl->getInspectedValue(label, str )) )
    {
      printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    else
    {
      if ( !strcasecmp(str.c_str(), "Upstream") )
      {
        BaseBc[n].set_FaceMode(BoundaryOuter::prdc_upstream);
      }
      else if ( !strcasecmp(str.c_str(), "Downstream") )
      {
        BaseBc[n].set_FaceMode(BoundaryOuter::prdc_downstream);
      }
      else
      {
        printf("\tParsing error : Invalid keyword in '%s'\n", label.c_str());
        Exit(0);
      }
    }
  }
  
  // Driver
  if ( BaseBc[n].getPrdcMode() == BoundaryOuter::prdc_Driver )
  {
    
    label = label_base + "/DriverDirection";
    if ( !(tpCntl->getInspectedValue(label, str )) )
    {
      printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    else
    {
      if ( !strcasecmp("xminus", str.c_str()) )
      {
        def = X_minus;
      }
      else if ( !strcasecmp("xplus", str.c_str()) )
      {
        def = X_plus;
      }
      else if ( !strcasecmp("yminus", str.c_str()) )
      {
        def = Y_minus;
      }
      else if ( !strcasecmp("yplus", str.c_str()) )
      {
        def = Y_plus;
      }
      else if ( !strcasecmp("zminus", str.c_str()) )
      {
        def = Z_minus;
      }
      else if ( !strcasecmp("zplus", str.c_str()) )
      {
        def = Z_plus;
      }
      else
      {
        printf("\tParsing error : Invalid keyword in '%s'\n", label.c_str());
        Exit(0);
      }
      BaseBc[n].set_DriverDir(def);
    }
    
    
    label = label_base + "/DriverLidIndex";
    if ( !(tpCntl->getInspectedValue(label, def )) )
    {
      printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    else
    {
      BaseBc[n].setDriverIndex(def);
    }
  }
}



// #################################################################
/**
 * @brief 外部境界の流入条件のパラメータを取得する
 * @param [in] label_base ラベルディレクトリ
 * @param [in] n          面番号
 */
void ParseBC::getObcSpecVH(const string label_base, const int n)
{
  
  string label;
  
  // 速度境界条件のプロファイル
  BaseBc[n].set_V_Profile( getVprofile(label_base) );
  
  // 法線ベクトル
  if ( BaseBc[n].get_V_Profile() != CompoList::vel_zero )
  {
    REAL_TYPE v[3];
    label = label_base + "/OrientationVector";
    if ( !Control::getVec(label, v, tpCntl, true) ) Exit(0);
    BaseBc[n].addVec(v);
  }
  
  // 速度のパラメータ読み込み
  getVelocity(label_base, BaseBc[n].get_V_Profile(), BaseBc[n].ca, "OuterBoundary");

  
  // heat problem
  if ( HeatProblem )
  {
    if ( Unit_Param != DIMENSIONAL )
    {
      stamped_printf("\tError: Heat condition must be given by dimensional value\n");
      Exit(0);
    }
    
    label = label_base + "/Temperature";
    REAL_TYPE ct;
    
    if ( !(tpCntl->getInspectedValue(label, ct )) )
    {
      stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    else
    {
      BaseBc[n].setTemp( ct );
      BaseBc[n].set_hType(CNST_TEMP);
    }
  }
}


// #################################################################
/**
 * @brief トラクションフリーの外部境界のパラメータを取得する
 * @param [in] label_base ラベルディレクトリ
 * @param [in] n          面番号
 */
void ParseBC::getObcTrcfree(const string label_base, const int n)
{
  BaseBc[n].set_pType(P_DIRICHLET);
  BaseBc[n].p = 0.0; // ゲージ圧zero 固定
  
  /* 外部雰囲気温
  if ( HeatProblem )
  {
    if ( Unit_Param != DIMENSIONAL )
    {
      stamped_printf("\tError: Heat condition must be given by dimensional value\n");
      Exit(0);
    }
      
    string label = label_base + "/AmbientTemperature";
    REAL_TYPE ct;
    
    if ( !(tpCntl->getInspectedValue(label, ct )) )
    {
      stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }

    BaseBc[n].setTemp( ct );
  }
  */
}


// #################################################################
/**
 * @brief 外部境界の壁面条件のパラメータを取得する
 * @param [in] label_base ラベルディレクトリ
 * @param [in] n          面番号
 */
void ParseBC::getObcWall(const string label_base, const int n)
{
  REAL_TYPE ct=0.0;
  REAL_TYPE dv[3];
  string str;
  string label, label2;
  
  // 速度のタイプの特定
  label = label_base + "/Type";
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
    stamped_printf("\tParsing error : Invalid Specified_Type in '%s'\n", label.c_str());
    Exit(0);
  }
  
  if ( !strcasecmp("fixed", str.c_str()) )
  {
	  BaseBc[n].set_wallType(BoundaryOuter::fixed);
    BaseBc[n].set_V_Profile( CompoList::vel_zero );
  }
  else if ( !strcasecmp("slide", str.c_str()) )
  {
	  BaseBc[n].set_wallType(BoundaryOuter::slide);
    BaseBc[n].set_V_Profile( getVprofile(label_base) );
  }
  else
  {
	  printf("\tParsing error : Invalid string value '%s' for 'Type'\n", str.c_str());
	  Exit(0);
  }
  
  // 法線ベクトル
  if ( BaseBc[n].get_V_Profile() != CompoList::vel_zero )
  {
    label = label_base + "/OrientationVector";
    if ( !Control::getVec(label, dv, tpCntl, true) ) Exit(0);
    BaseBc[n].addVec(dv);
  }
  
  // 速度のパラメータ読み込み
  getVelocity(label_base, BaseBc[n].get_V_Profile(), BaseBc[n].ca, "OuterBoundary");
  
  
  
  // heat problem
  if ( HeatProblem )
  {
    label = label_base + "/ThermalOption";
    
    if ( !(tpCntl->getInspectedValue(label, str )) )
    {
      stamped_printf("\tParsing error : fail to get '%s'\n", label.c_str());
      Exit(0);
    }
    if ( Unit_Param != DIMENSIONAL )
    {
      stamped_printf("\tError: Heat condition must be given by dimensional value\n");
      Exit(0);
    }
    
    if ( !strcasecmp(str.c_str(), "Adiabatic") )
    {
      BaseBc[n].set_hType(ADIABATIC);
      BaseBc[n].p = 0.0;
    }
    else if( !strcasecmp(str.substr(0,12).c_str(), "HeatTransfer") )
    {
      BaseBc[n].set_hType(TRANSFER);
      getObcHeatTransfer(label, n, str);
    }
    else if( !strcasecmp(str.c_str(), "HeatFlux") )
    {
      BaseBc[n].set_hType(HEATFLUX);
      label2 = label_base + "/Flux";
      BaseBc[n].setHeatflux( getValueReal(label2) ); // 正符号は流入
    }
    else if( !strcasecmp(str.c_str(), "Isothermal") )
    {
      BaseBc[n].set_hType(ISOTHERMAL);
      label2 = label_base + "/Temperature";
      ct = getValueReal(label2); // 表面温度
      BaseBc[n].setTemp( ct );
    }
    else if( !strcasecmp(str.c_str(), "ConstantTemperature") )
    {
      BaseBc[n].set_hType(CNST_TEMP);
      label = label_base + "/Temperature";
      ct = getValueReal(label); // 指定温度
      BaseBc[n].setTemp( ct );
    }
    else
    {
      stamped_printf("\tParsing error : Invalid string value for 'ThermalOption' : %s\n", str.c_str());
      Exit(0);
    }
  }
}



// #################################################################
// 2相流問題で気相か液相かを取得する
void ParseBC::get_Phase(CompoList* cmp)
{
  string str,p;
  string label,label_base;
  int NoParam;
  
  ///////////////////////////////////////////////////////////////////////////////
  stamped_printf("\tWARNING not yet\n");
  Exit(0);
  ///////////////////////////////////////////////////////////////////////////////
  
  
  label_base="/PhaseIdetification";
  //cout <<  "label : " << label << endl;
  if ( !(tpCntl->chkNode(label_base)) ) {
    stamped_printf("\tParsing error : Missing the section of 'PhaseIdetification'\n");
    Exit(0);
  }
  
  // load statement list
  NoParam=tpCntl->countLabels(label_base);
  if ( NoParam < 0) {
    stamped_printf("\tcountLabels --- %s\n",label_base.c_str());
    Exit(0);
  }
  
  for (int i=1; i<=NoParam; i++) {
    if ( !(tpCntl->getNodeStr(label_base,i,str)) )
    {
      stamped_printf("\tParsing error : No Leaf Node \n");
      Exit(0);
    }
    label=label_base+"/"+str;//
    
    /*
     // state
     p = param->GetName();
     if ( !strcasecmp(p, "gas") ) {
     m_phase = GAS;
     }
     else if ( !strcasecmp(p, "liquid") ) {
     m_phase = LIQUID;
     }
     else {
     stamped_printf("\tParsing error : No valid keyword '%s' in 'Phase_Idetification'\n", p);
     Exit(0);
     }
     
     // ID
     if ( !param->isSetID() ) {
     stamped_printf("\tParsing error : No ID for statement in 'Phase_Idetification'\n");
     Exit(0);
     }
     if ( -1 == (id=param->GetID()) ) {
     stamped_printf("\tParsing error : No valid ID for statement in 'Phase_Idetification'\n");
     Exit(0);
     }
     
     // IDがiTableの中にリストアップされているかを調べる
     if ( (0 >= def) || (def > NoMedium) ) {
     stamped_printf("\tParsing error : ID[%d] described in 'id' is not included in 'Model_Setting'\n", id);
     Exit(0);
     }
     
     // set phase of FLUID
     for (int n=1; n<=NoMedium; n++) {
     if ( cmp[n].getID() == id ) {
     if ( cmp[n].getState() == FLUID ) cmp[n].setPhase(m_phase);
     }
     }
     */
    
  }
  
  // check Phase of Fluid
  int tmp;
  bool sw=true;
  for (int n=1; n<=NoCompo; n++)
  {
    if ( (cmp[n].getState() == FLUID) && cmp[n].isKindMedium() )
    {
      tmp = cmp[n].getPhase();
      if ( (tmp!=GAS) && (tmp!=LIQUID) )
      {
        stamped_printf("\tcomponent [%d] is fluid, but not identified by gas or liquid.\n", n);
        sw = false;
      }
    }
  }
  if ( sw == false ) Exit(0);
}



// #################################################################
//@brief 外部境界の速度境界条件のタイプを取得し，返す
int ParseBC::getVprofile(const string label_base)
{
  string label, str;
  
  label = label_base + "/Profile";
  
  if ( !(tpCntl->getInspectedValue(label, str )) )
  {
    printf("\tParsing error : fail to get 'Profile' in '%s'\n", str.c_str());
    Exit(0);
  }
  if ( !strcasecmp("constant", str.c_str()) )
  {
		return CompoList::vel_constant;
  }
  else if ( !strcasecmp("harmonic", str.c_str()) )
  {
		return CompoList::vel_harmonic;
  }
  else
  {
	  printf("\tParsing error : Invalid string value for '%s' : %s\n", label.c_str(), str.c_str());
	  Exit(0);
  }
  return 0;
}


// #################################################################
/**
 @brief 速度のパラメータを取得する
 @param label_base 
 @param ca 係数パラメータの配列
 @param policy 流量指定(true) or 速度指定(false) > 外部境界は速度のみ
 @note 
 - 値は，Control::setParameters()で無次元化する
 - 速度プロファイルは単振動と一定値の場合で係数の保持パターンが異なる
 - 内部境界の場合には，流量指定と速度指定があるので分岐処理（未実装）
 */
void ParseBC::getVelocity(const string label_base, const int prof, REAL_TYPE* ca, const char* str, const bool policy)
{
  string label;
  REAL_TYPE ct=0.0, vel;
  
  if ( prof == CompoList::vel_zero )
  {
    ca[CompoList::amplitude] = 0.0;
    ca[CompoList::frequency] = 0.0;
    ca[CompoList::initphase] = 0.0;
    ca[CompoList::bias]      = 0.0;
  }
  else if ( prof == CompoList::vel_constant )
  {
    label = label_base + "/Velocity";
    
    if ( !(tpCntl->getInspectedValue(label, ct )) )
    {
      stamped_printf("\tParsing error : Invalid value in '%s > Velocity'\n", str);
      Exit(0);
    }
    if ( policy ) // 流量指定の場合 LocalBoundaryの場合のみ
    {
      vel = ct; // 有次元でも無次元でも，モデルの断面積を計算して，後ほどパラメータ計算 >> Control::setParameters()
    }
    else
    {
      vel = ( Unit_Param == DIMENSIONAL ) ? ct : ct * RefVelocity;
    }
    
    ca[CompoList::amplitude] = 0.0;
    ca[CompoList::frequency] = 0.0;
    ca[CompoList::initphase] = 0.0;
    ca[CompoList::bias]      = vel;
  }
  else if ( prof == CompoList::vel_harmonic)
  {
    label = label_base + "/Amplitude";
    
    if ( !(tpCntl->getInspectedValue(label, ct )) )
    {
      stamped_printf("\tParsing error : fail to get 'Amplitude' in '%s'\n", label_base.c_str());
      Exit(0);
    }
    ca[CompoList::amplitude] = ct;
    
    label = label_base + "/Frequency";
    
    if ( !(tpCntl->getInspectedValue(label, ct )) )
    {
      stamped_printf("\tParsing error : fail to get 'Frequency' in '%s'\n", label_base.c_str());
      Exit(0);
    }
    ca[CompoList::frequency] = ct;
    
    label = label_base + "/InitialPhase";
    
    if ( !(tpCntl->getInspectedValue(label, ct )) )
    {
      stamped_printf("\tParsing error : fail to get 'InitialPhase' in '%s'\n", label_base.c_str());
      Exit(0);
    }
    ca[CompoList::initphase] = ct;
    
    
    label = label_base + "/ConstantBias";
    
    if ( !(tpCntl->getInspectedValue(label, ct )) )
    {
      stamped_printf("\tParsing error : fail to get 'ConstantBias' in '%s'\n", label_base.c_str());
      Exit(0);
    }
    ca[CompoList::bias] = ct;
  }
  else
  {
    Exit(0);
  }
  
}


// #################################################################
// TPのポインタを受け取る
void ParseBC::importTP(TextParser* tp)
{ 
  if ( !tp ) Exit(0);
  tpCntl = tp;
}



// #################################################################
// CompoListに内部境界条件の情報を設定する
/**
 * @note 最初にBCの情報を登録，その後IDの情報を登録
 * @note パラメータファイルから各内部BCのidをパースし，cmpに保持する
 * @note 格納番号は1からスタート
 */
void ParseBC::loadLocalBC(Control* C, MediumList* mat, CompoList* cmp)
{ 
  string str, label;
  string label_base, label_ename, label_leaf;
  

  label_base = "/BcTable/LocalBoundary";
  

  //
  // 内部境界の条件設定 --- NoBC = 内部境界の数
  //
  
  for (int k=1; k<=NoBC; k++)
  {
    int m = NoMedium + k;
    
    if ( !(tpCntl->getNodeStr(label_base, k, str)) )
    {
      stamped_printf("\tParsing error : No Leaf Node \n");
      Exit(0);
    }
    
    cmp[m].setAlias(str);
    
    
    label_leaf = label_base + "/" + str;
    
    
    // medium of alias
    label = label_leaf + "/Medium";
    
    if ( !(tpCntl->getInspectedValue(label, str )) )
    {
      Hostonly_ stamped_printf("\tParsing error : No '%s'\n", label.c_str());
      Exit(0);
    }
    cmp[m].setMedium(str);
    
    
    // class
    label = label_leaf + "/Class";
    
    if ( !(tpCntl->getInspectedValue(label, str )) )
    {
      Hostonly_ stamped_printf("\tParsing error : No '%s'\n", label.c_str());
      Exit(0);
    }
    
    // check consitency between Polylib file and input parameter file
    if ( C->Mode.Example == id_Polygon )
    {
      string m_class, m_medium;
      loadLocalBCfromPolylibFile(C, cmp[m].getAlias(), m_class, m_medium);
      
      if ( !strcasecmp(m_class.c_str(), str.c_str()) )
      {
        ;
      }
      else
      {
        Hostonly_ stamped_printf("\tError : Check consistency of LocalBC class name in between Polylib file and input parameter file. [%s , %s]\n", m_class.c_str(), str.c_str());
        Exit(0);
      }
      
      if ( !strcasecmp(m_medium.c_str(), cmp[m].getMedium().c_str()) )
      {
        ;
      }
      else
      {
        Hostonly_ stamped_printf("\tError : Check consistency of Medium name in between Polylib file and input parameter file. [%s , %s]\n", m_medium.c_str(), cmp[m].getMedium().c_str());
        Exit(0);
      }
    }
    
    
    // cmp[].type, h_typeのセット ---> setType
    setKeywordLBC(str, m, cmp);
    
    
    // 各BCの処理
    int tp = cmp[m].getType();
    
    if      ( tp == OBSTACLE )      { ; }
    else if ( tp == SPEC_VEL )      getIbcSpecVel(label_leaf, m, cmp);
    else if ( tp == OUTFLOW )       getIbcOutflow(label_leaf, m, cmp);
    else if ( tp == IBM_DF )        get_IBC_IBM_DF(label_leaf, m, cmp);
    else if ( tp == HEX )           get_IBC_PrsLoss(label_leaf, m, cmp);
    else if ( tp == FAN )           get_IBC_Fan(label_leaf, m, cmp);
    else if ( tp == DARCY )         get_Darcy(label_leaf, m, cmp);
    else if ( tp == PERIODIC )      getIbcPeriodic(label_leaf, m, cmp);
    else if ( tp == MONITOR )       { ; }
    else if ( tp == SOLIDREV )      getIbcSolidRev(label_leaf, m, cmp);
    
    
    if ( HeatProblem )
    {
      cmp[m].setHeatmode(ON);
      
      if ( C->KindOfSolver == FLOW_ONLY )
      {
        Hostonly_ stamped_printf("Parse Error : Heat BC is not allowed on FLOWONLY mode.\n");
        Exit(0);
      }
      
      if ( Unit_Param != DIMENSIONAL )
      {
        Hostonly_ stamped_printf("\tError : Heat condition must be given by dimensional value\n");
        Exit(0);
      }
      
      if      ( tp == OBSTACLE )   cmp[m].setHeatflux( 0.0 );
      else if ( tp == ADIABATIC )  cmp[m].setHeatflux( 0.0 );
      else if ( tp == HEATFLUX )   getIbcHeatFlux(label_leaf, m, cmp);
      else if ( tp == TRANSFER ) 
      {
        switch ( cmp[m].getHtype() )
        {
          case HT_S:
            getIbcHT_S(label_leaf, m, cmp);
            break;
            
          case HT_SN:
            getIbcHT_SN(label_leaf, m, cmp);
            break;
            
          case HT_SF:
            getIbcHT_SF(label_leaf, m, cmp);
            break;
        }        
      }
      else if ( tp == ISOTHERMAL )   getIbcIsoTherm(label_leaf, m, cmp);
      else if ( tp == RADIANT )      get_IBC_Radiant(label_leaf, m, cmp);
      else if ( tp == HEAT_SRC )     getIbcHeatSrc(label_leaf, m, cmp);
      else if ( tp == CNST_TEMP )    getIbcCnstTemp(label_leaf, m, cmp);
      else if ( tp == SPEC_VEL ) { ; } // 既に読み込み済み
      else if ( tp == SOLIDREV ) { ; } // 既に読み込み済み
      else 
      {
        Hostonly_ printf("\tError : Invalid Local BC keyword [%d]\n", tp);
        Exit(0);
      }
    }
  }
  
  
  // この時点まで，mat[]の媒質情報は媒質を保持しているオーダーにしかないので，媒質情報をLBCのオーダーにコピーする
  for (int k=1; k<=NoBC; k++)
  {
    int m = NoMedium + k; // cmp[]のLBCのインデクス
    
    // 媒質ラベルからmat[]の格納番号のサーチ
    int odr = -1;
    
    for (int i=1; i<=NoMedium; i++)
    {
      if ( !strcasecmp( cmp[m].getMedium().c_str(), mat[i].getAlias().c_str()) )
      {
        odr = i;
        break;
      }
    }
    
    if ( (odr < 1) || (odr > NoMedium) )
    {
      Hostonly_ stamped_printf("\tSomthing wrong %d\n", odr);
      Exit(0);
    }
    
    // LocalBC分の媒質情報のコピー
    if ( mat[odr].getState() == FLUID )
    {
      mat[m].P[p_density]              = mat[odr].P[p_density];
      mat[m].P[p_kinematic_viscosity]  = mat[odr].P[p_kinematic_viscosity];
      mat[m].P[p_viscosity]            = mat[odr].P[p_viscosity];
      mat[m].P[p_thermal_conductivity] = mat[odr].P[p_thermal_conductivity];
      mat[m].P[p_specific_heat]        = mat[odr].P[p_specific_heat];
      mat[m].P[p_speed_of_sound]       = mat[odr].P[p_speed_of_sound];
      mat[m].P[p_vol_expansion]        = mat[odr].P[p_vol_expansion];
    }
    else  // solid
    {
      mat[m].P[p_density]              = mat[odr].P[p_density];
      mat[m].P[p_specific_heat]        = mat[odr].P[p_specific_heat];
      mat[m].P[p_thermal_conductivity] = mat[odr].P[p_thermal_conductivity];
    }
    
    // cmp[]のaliasをmat[]へコピーしておく
    mat[m].setAlias( cmp[m].getAlias() );
    
    // cmp[]の状態を設定
    mat[m].setState(mat[odr].getState());
    cmp[m].setState(mat[odr].getState());
  }
  
  
  for (int k=1; k<=NoMedium; k++)
  {
    cmp[k].setState(mat[k].getState());
    cmp[k].setAlias(mat[k].getAlias());
    cmp[k].setMedium(mat[k].getAlias());
  }
  
#if 0
  for (int k=0; k<=NoCompo; k++)
  {
    printf("%3d : %3d %16s %16s %e %e %e %e %e %e %e\n",
           k,
           cmp[k].getState(),
           cmp[k].getAlias().c_str(),
           cmp[k].getMedium().c_str(),
           mat[k].P[p_density],
           mat[k].P[p_kinematic_viscosity],
           mat[k].P[p_viscosity],
           mat[k].P[p_thermal_conductivity],
           mat[k].P[p_specific_heat],
           mat[k].P[p_speed_of_sound],
           mat[k].P[p_vol_expansion]
           );
  }
#endif
}


// #################################################################
/** 
 * @brief 内部境界条件のalias名に対するclassとmediumをpolylib.tpから取得する
 * @param [in]  C        pointer to Control class
 * @param [in]  obj_name object name
 * @param [out] m_lass   class name
 * @param [out] m_medium medium name
 */
void ParseBC::loadLocalBCfromPolylibFile(Control* C, const string obj_name, string& m_class, string& m_medium)
{
  // ffvのパラメータローダのインスタンス生成
  TextParser tp;
  int ierror=0;
  
  if ( (ierror = tp.read(C->PolylibConfigName)) != TP_NO_ERROR )
  {
    Hostonly_ stamped_printf("\tError at reading '%s' file : %d\n", C->PolylibConfigName.c_str(), ierror);
    Exit(0);
  }
  
  // Number of polygon objects
  vector<string> nodes;
  tp.getLabelVector("/Polylib", nodes);
  
  int m_obj = nodes.size();
  
  if ( m_obj != C->NoBC )
  {
    stamped_printf("\tError : Number of LocalBC does not agree between '%s' and input parameter file. \n",
                   C->PolylibConfigName.c_str());
    Exit(0);
  }
  
  
  string label, str;
  
  // Does obj_name exist?
  bool flag=false;
  for (vector<string>::iterator it = nodes.begin(); it != nodes.end(); it++)
  {
    if ( !strcasecmp((*it).c_str(), obj_name.c_str()) ) flag = true;
  }
  if ( !flag )
  {
    Hostonly_ stamped_printf("\tParsing error : No '%s'\n", label.c_str());
    Exit(0);
  }
  
  label = "/Polylib/" + obj_name;
  
  // class name
  string label2 = label + "/type";
  
  if ( !(tp.getInspectedValue(label2, str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : No '%s'\n", label2.c_str());
    Exit(0);
  }
  m_class = str;

  // medium name
  string label3 = label + "/label";
  
  if ( !(tp.getInspectedValue(label3, str )) )
  {
    Hostonly_ stamped_printf("\tParsing error : No '%s'\n", label3.c_str());
    Exit(0);
  }
  m_medium = str;
  
}


// #################################################################
// 外部境界条件を取得，保持する
void ParseBC::loadOuterBC(BoundaryOuter* bc, const MediumList* mat, CompoList* cmp, int* ensPrdc)
{
  string label_base, label_leaf, label;
  string str;
  
  // Basic Outer BCリストの読み込み
  label_base = "/BcTable/OuterBoundary";
  
  if ( !(tpCntl->chkNode(label_base)) )
  {
    Hostonly_ stamped_printf("\tParsing error : Missing '%s'\n", label.c_str());
    Exit(0);
  }
  
  vector<string> nodes;
  
  // 直下のラベルを取得
  tpCntl->getLabelVector(label_base, nodes);
  
  
  // ラベルの重複チェックとセット
  // ラベル名はパラメータファイルの出現順になっている
  int m = 1;
  for (vector<string>::iterator it = nodes.begin(); it != nodes.end(); it++)
  {
    if ( !chkDuplicate(m, *it) ) Exit(0);
    
    if ( strcasecmp( (*it).c_str(), "FaceBC") )
    {
      BaseBc[m].setAlias(*it);
      m++;
    }
  }
  
  
  //BasicBCs
  for (int i=1; i<=NoBaseBC; i++)
  {
    if ( !(tpCntl->getNodeStr(label_base, i, str)) )
    {
      Hostonly_ printf("\tParsing error : Missing candidates of outer boundary\n");
      Exit(0);
    }
    
    label_leaf = label_base + "/" + str;
    
    // Classに境界条件の種別をセットする
    label = label_leaf + "/Class";
    
    if ( !(tpCntl->getInspectedValue(label, str )) )
    {
      Hostonly_ printf("\tParsing error : No Class\n");
      Exit(0);
    }
    setKeywordOBC(str, i);
  
    
    // ガイドセルの媒質ラベルを取得   周期境界は不要
    if ( BaseBc[i].getClass() != OBC_PERIODIC )
    {
      label = label_leaf + "/Medium";
      
      if ( !(tpCntl->getInspectedValue(label, str )) )
      {
        Hostonly_ printf("\tParsing error : No entory '/BCTable/OUterBoundary/*/Medium'\n");
        Exit(0);
      }
      
      
      // ラベル名が一致するエントリ番号をセットする 
      for (int m=1; m<=NoCompo; m++)
      {
        if ( !strcasecmp( str.c_str(), mat[m].getAlias().c_str() ) )
        {
          BaseBc[i].setGuideMedium(m);
          break;
        }
      }
    }

    
    // 各条件に応じたパラメータをロード
    switch ( BaseBc[i].getClass() )
    {
      case OBC_WALL:
        getObcWall(label_leaf, i);
        break;
        
      case OBC_OUTFLOW:
        getObcOutflow(label_leaf, i);
        break;
        
      case OBC_SPEC_VEL:
        getObcSpecVH(label_leaf, i);
        break;
        
      case OBC_TRC_FREE:
        getObcTrcfree(label_leaf, i);
        break;
        
      case OBC_FAR_FIELD:
        getObcFarField(label_leaf, i);
        break;
        
      case OBC_PERIODIC:
        getObcPeriodic(label_leaf, i);
        break;
        
      // nothing to do
      //case OBC_SYMMETRIC:
      //case OBC_INTRINSIC:
      //break;
    }
  }
  
  
  
  // 各フェイスに境界条件を設定する
  label_base = "/BcTable/OuterBoundary/FaceBC";
  
  if ( !(tpCntl->chkNode(label_base)) )
  {
    Hostonly_ printf("\tParsing error : Missing OuterBoundary '%s'\n", label_base.c_str());
    Exit(0);
  }
  
  // check
  int nnode = tpCntl->countLabels(label_base);
  if ( nnode != NOFACE ) 
  {
    Hostonly_ printf("\tParsing error : OuterBoundary FaceBC count != 6\n");
    Exit(0);
  }
  
  
  // 各面に与える境界条件番号を取得し，BaseBcから境界情報リストに内容をコピー
  label = label_base + "/Xminus";
  if ( !(tpCntl->getInspectedValue(label, str )) ) Exit(0);
  if ( !findOBClist(str, X_minus, bc) )
  {
    Hostonly_ printf("\tParsing error : Alias cannot be found : FaceBC/Xminus\n");
    Exit(0);
  }
  
  label = label_base + "/Xplus";
  if ( !(tpCntl->getInspectedValue(label, str )) ) Exit(0);
  if ( !findOBClist(str, X_plus, bc) )
  {
    Hostonly_ printf("\tParsing error : Alias cannot be found : FaceBC/Xplus\n");
    Exit(0);
  }
  
  label = label_base + "/Yminus";
  if ( !(tpCntl->getInspectedValue(label, str )) ) Exit(0);
  if ( !findOBClist(str, Y_minus, bc) )
  {
    Hostonly_ printf("\tParsing error : Alias cannot be found : FaceBC/Yminus\n");
    Exit(0);
  }
  
  label = label_base + "/Yplus";
  if ( !(tpCntl->getInspectedValue(label, str )) ) Exit(0);
  if ( !findOBClist(str, Y_plus, bc) )
  {
    Hostonly_ printf("\tParsing error : Alias cannot be found : FaceBC/Yplus\n");
    Exit(0);
  }
  
  label = label_base + "/Zminus";
  if ( !(tpCntl->getInspectedValue(label, str )) ) Exit(0);
  if ( !findOBClist(str, Z_minus, bc) )
  {
    Hostonly_ printf("\tParsing error : Alias cannot be found : FaceBC/Zminus\n");
    Exit(0);
  }
  
  label = label_base + "/Zplus";
  if ( !(tpCntl->getInspectedValue(label, str )) ) Exit(0);
  if ( !findOBClist(str, Z_plus, bc) )
  {
    Hostonly_ printf("\tParsing error : Alias cannot be found : FaceBC/Zplus\n");
    Exit(0);
  }
  

  
  // 周期境界条件の整合性のチェック
  
  // 部分周期境界の数
  int p_flag=0;
  for (int n=0; n<NOFACE; n++)
  {
    if (bc[n].getPrdcMode() == BoundaryOuter::prdc_Driver) p_flag++;
  }
  
  // 部分周期条件を使わない場合，対になる外部境界のチェック
  if ( p_flag == 0 )
  {
    int n_pair=0;
    
    // 周期境界条件の指定をチェック
    for (int n=0; n<NOFACE; n++)
    {
      if ( bc[n].getClass() == OBC_PERIODIC )
      {
        n_pair = oppositeDir(n);
        if ( bc[n_pair].getClass() != OBC_PERIODIC )
        {
          Hostonly_ printf("\tFace BC[rank=%d] : No consistent Periodic Bnoudary in %s direction\n", myRank, FBUtility::getDirection(n_pair).c_str());
          Exit(0);
        }
      }
    }
    
    // 対になるモードのチェック
    for (int n=0; n<NOFACE; n++)
    {
      if ( bc[n].getClass() == OBC_PERIODIC )
      {
        n_pair = oppositeDir(n);
        
        switch (bc[n].getPrdcMode())
        {
          case BoundaryOuter::prdc_Simple:
            if ( bc[n_pair].getPrdcMode() != BoundaryOuter::prdc_Simple )
            {
              Hostonly_ printf("\tFace BC : No consistent SIMPLE Periodic Bnoudary in %s direction\n", FBUtility::getDirection(n_pair).c_str());
              Exit(0);
            }
            break;
            
          case BoundaryOuter::prdc_Directional:
            if ( bc[n_pair].getPrdcMode() != BoundaryOuter::prdc_Directional )
            {
              Hostonly_ printf("\tFace BC : No consistent DIRECTIONAL Periodic Bnoudary in %s direction\n", FBUtility::getDirection(n_pair).c_str());
              Exit(0);
            }
            if ( bc[n].p != bc[n_pair].p ) // 同じ値が入っていること
            {
              Hostonly_ printf("\tFace BC : Pressure difference value is not same in %s direction\n", FBUtility::getDirection(n_pair).c_str());
              Exit(0);
            }
            if ( (bc[n].get_FaceMode() == BoundaryOuter::prdc_upstream) && (bc[n_pair].get_FaceMode() != BoundaryOuter::prdc_downstream) )
            {
              Hostonly_ printf("\tFace BC : No consistent Upstream/Downstream relation in %s direction\n", FBUtility::getDirection(n).c_str());
              Exit(0);
            }
            if ( (bc[n].get_FaceMode() == BoundaryOuter::prdc_downstream) && (bc[n_pair].get_FaceMode() != BoundaryOuter::prdc_upstream) )
            {
              Hostonly_ printf("\tFace BC : No consistent Upstream/Downstream relation in %s direction\n", FBUtility::getDirection(n).c_str());
              Exit(0);
            }
            
            break;
        }        
      }
    }
    
    // 周期境界の方向を保存
    for (int n=0; n<NOFACE; n++)
    {
      if ( bc[n].getClass() == OBC_PERIODIC )
      {
        if ( n == X_minus )      ensPrdc[0]=ON;
        else if ( n == X_plus )  ensPrdc[0]=ON;
        else if ( n == Y_minus ) ensPrdc[1]=ON;
        else if ( n == Y_plus )  ensPrdc[1]=ON;
        else if ( n == Z_minus ) ensPrdc[2]=ON;
        else if ( n == Z_plus )  ensPrdc[2]=ON;
      }
    }
    
  }
  else { // Driverが指定された場合の内部境界との整合性をチェック
    int n_pair=0;
    for (int n=0; n<NOFACE; n++)
    {
      n_pair = oppositeDir(n);
      if ( bc[n].getClass() == OBC_PERIODIC )
      {
        if (bc[n].getPrdcMode() == BoundaryOuter::prdc_Driver)
        {
          // 他方は周期境界以外であること
          if ( bc[n_pair].getClass() == OBC_PERIODIC )
          {
            Hostonly_ printf("\tFace BC : %s direction should be non periodic BC\n", FBUtility::getDirection(n_pair).c_str());
            Exit(0);
          }
          
          int cflag=0;
          for (int c=1; c<=NoCompo; c++)
          {
            if ( (cmp[c].getType() == PERIODIC) )
            {
              if ( (int)cmp[c].getPeriodicDir() != bc[n].get_DriverDir() )
              {
                Hostonly_ printf("\tPeriodic Driver BC : No consistent Periodic Bnoudary in %s direction\n", FBUtility::getDirection(n_pair).c_str());
                Exit(0);
              }
              else
              {
                cflag++;
              }
            }
          }
          if (cflag != 1)
          {
            Hostonly_ printf("\tPeriodic Driver BC can not detemine uniquely\n");
            Exit(0);
          }
        }        
      }
    }
  }
  
}


// #################################################################
/**
 * @brief 外部境界面の反対方向を返す
 * @param [in] dir 評価する方向
 * @return dirと反対方向
 */
int ParseBC::oppositeDir(const int dir)
{
  int n_pair=0;
  
  if      ( dir == X_minus ) n_pair=X_plus;
  else if ( dir == X_plus )  n_pair=X_minus;
  else if ( dir == Y_minus ) n_pair=Y_plus;
  else if ( dir == Y_plus )  n_pair=Y_minus;
  else if ( dir == Z_minus ) n_pair=Z_plus;
  else if ( dir == Z_plus )  n_pair=Z_minus;
  
  return n_pair;
}


// #################################################################
// コンポーネントの情報を表示する
void ParseBC::printCompo(FILE* fp, const int* gci, const MediumList* mat, CompoList* cmp, const BoundaryOuter* bc)
{
  int m;
  bool flag;
  Vec3i st, ed;
  
  // VBC ---------------------------------------------------
  if ( existComponent(SPEC_VEL, cmp) )
  {
    fprintf(fp, "\n\t[Specified_Velocity]\n");
    fprintf(fp, "\t no                    Label    i_st    i_ed    j_st    j_ed    k_st    k_ed   Area[m*m]    Elements\n");
    
    for (int n=1; n<=NoCompo; n++)
    {
      if ( cmp[n].getType() == SPEC_VEL )
      {
        st = getCmpGbbox_st(n, gci);
        ed = getCmpGbbox_ed(n, gci);
        
        fprintf(fp,"\t%3d %24s %7d %7d %7d %7d %7d %7d %11.4e %11ld\n\n",
                n, cmp[n].getAlias().c_str(),
                st.x, ed.x, st.y, ed.y, st.z, ed.z,
                cmp[n].area, cmp[n].getElement());
        fprintf(fp, "\t                     vector_x    vector_y    vector_z   Direction");
        
        if ( cmp[n].get_V_Profile() == CompoList::vel_zero )
        {
          fprintf(fp, "\n");
        }
        else if ( cmp[n].get_V_Profile() == CompoList::vel_constant )
        {
          fprintf(fp, "    Q[m^3/s]   Vel.[m/s]\n");
        }
        else
        {
          fprintf(fp, "    Q[m^3/s]   Amp.[m/s]   Freq.[Hz]  Phase[rad] Intcpt[m/s] Strauhal[-]\n");
        }
      }
    }
    
    
    for (int n=1; n<=NoCompo; n++)
    {
      if ( cmp[n].getType() == SPEC_VEL )
      {
        fprintf(fp,"\t\t\t   %10.3e  %10.3e  %10.3e  %10s ",
                cmp[n].nv[0], cmp[n].nv[1], cmp[n].nv[2],
                (cmp[n].getBClocation()==CompoList::same_direction) ? "In" : "Out");
        
        if ( cmp[n].get_V_Profile() == CompoList::vel_zero )
        {
          fprintf(fp,"\n");
        }
        else if ( cmp[n].get_V_Profile() == CompoList::vel_constant )
        {
          fprintf(fp,"%11.3e %11.3e\n", cmp[n].ca[CompoList::bias]*cmp[n].area, cmp[n].ca[CompoList::bias] );
        }
        else
        {
          fprintf(fp,"%11.3e %11.3e %11.3e %11.3e %11.3e %11.4e\n",
                  cmp[n].ca[CompoList::amplitude]*cmp[n].area,
                  cmp[n].ca[CompoList::amplitude],
                  cmp[n].ca[CompoList::frequency],
                  cmp[n].ca[CompoList::initphase],
                  cmp[n].ca[CompoList::bias],
                  cmp[n].ca[CompoList::frequency]*RefLength/RefVelocity);
        }
      }
    }
    
    // with constant temperature
    if ( existComponent(SPEC_VEL, cmp) && HeatProblem )
    {
      fprintf(fp, "\n\t[Specified_Velocity with Constant Temperature]\n");
      fprintf(fp, "\t no                    Label   # of faces    i_st    i_ed    j_st    j_ed    k_st    k_ed      Temp(C)      Temp[-]\n");
      
      for (int n=1; n<=NoCompo; n++)
      {
        if ( cmp[n].getType() == SPEC_VEL )
        {
          st = getCmpGbbox_st(n, gci);
          ed = getCmpGbbox_ed(n, gci);
          
          fprintf(fp, "\t%3d %24s %12ld %7d %7d %7d %7d %7d %7d %12.4e %12.4e\n", 
                  n, cmp[n].getAlias().c_str(), cmp[n].getElement(),
                  st.x, ed.x, st.y, ed.y, st.z, ed.z, 
									cmp[n].getTemp(), 
									FBUtility::convTempD2ND(cmp[n].getTemp(), BaseTemp, DiffTemp)); // 保持されている温度はCelsius
        }
      }
      fprintf(fp, "\n");
    }
  }
  
  if ( existComponent(OUTFLOW, cmp) )
  {
    fprintf(fp, "\n\t[Outflow]\n");
    fprintf(fp, "\t no                    Label   # of faces    i_st    i_ed    j_st    j_ed    k_st    k_ed  outflow_vel  pressure\n");
    
    for (int n=1; n<=NoCompo; n++)
    {
      if ( cmp[n].getType() == OUTFLOW )
      {
        st = getCmpGbbox_st(n, gci);
        ed = getCmpGbbox_ed(n, gci);
        
        fprintf(fp,"\t%3d %24s %12ld %7d %7d %7d %7d %7d %7d ",
                n, cmp[n].getAlias().c_str(), cmp[n].getElement(),
                st.x, ed.x, st.y, ed.y, st.z, ed.z);
        
        if (cmp[n].get_P_BCtype() == P_DIRICHLET)
        {
          fprintf(fp, "%12.4e;", FBUtility::convPrsND2D(cmp[n].get_Pressure(), BasePrs, RefDensity, RefVelocity, Unit_Prs) );
        }
        else
        {
          fprintf(fp," Grad_p = 0   ---");
        }
        fprintf(fp, "\n");
      }
    }
    fprintf(fp, "\n");
  }
  
  // Forcing
  if ( existComponent(IBM_DF, cmp) )
  {
    fprintf(fp, "\n\t[IBM_DF]\n");
    fprintf(fp, "\t no                    Label   i_st    i_ed    j_st    j_ed    k_st    k_ed\n");
    
    for (int n=1; n<=NoCompo; n++)
    {
      if ( cmp[n].getType() == IBM_DF )
      {
        st = getCmpGbbox_st(n, gci);
        ed = getCmpGbbox_ed(n, gci);
        
        fprintf(fp,"\t%3d %24s %7d %7d %7d %7d %7d %7d\n",
                n, cmp[n].getAlias().c_str(),
                st.x, ed.x, st.y, ed.y, st.z, ed.z);
      }
    }
    fprintf(fp, "\n");
    
    fprintf(fp, "\t no                    Label   vector_x   vector_y   vector_z    Vel.[m/s]      Vel.[-]\n");
    
    for (int n=1; n<=NoCompo; n++)
    {
      if ( cmp[n].getType() == IBM_DF )
      {
        fprintf(fp,"\t%3d %24s %10.3e %10.3e %10.3e %12.4e %12.4e \n",
                n, cmp[n].getAlias().c_str(), cmp[n].nv[0], cmp[n].nv[1], cmp[n].nv[2],
                cmp[n].get_Velocity(), FBUtility::convVelD2ND(cmp[n].get_Velocity(), RefVelocity));
      }
    }
  }
  
  // Forcing ---------------------------------------------------
	// Heat Exchanger
  if ( existComponent(HEX, cmp) )
  {
    fprintf(fp, "\n\t[Heat Exchanger]\n");
    
    fprintf(fp, "\t no                    Label   vector_x   vector_y   vector_z     O_x[m]     O_y[m]     O_z[m]      dir_x      dir_y      dir_z\n");
    for (int n=1; n<=NoCompo; n++)
    {
      if ( cmp[n].getType() == HEX )
      {
        fprintf(fp, "\t%3d %24s %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n", 
                n, cmp[n].getAlias().c_str(),
								cmp[n].nv[0], cmp[n].nv[1], cmp[n].nv[2],
                cmp[n].oc[0], cmp[n].oc[1], cmp[n].oc[2],
                cmp[n].dr[0], cmp[n].dr[1], cmp[n].dr[2]);
      }
    }
    fprintf(fp, "\n");
    
    fprintf(fp, "\t                                     Depth[m]   Width[m]  Height[m]  Area[m*m]\n");
    for (int n=1; n<=NoCompo; n++)
    {
      if ( cmp[n].getType() == HEX )
      {
        fprintf(fp, "\t%3d %24s %10.3e %10.3e %10.3e %10.3e\n", 
                n, cmp[n].getAlias().c_str(),
								cmp[n].depth, cmp[n].shp_p1, cmp[n].shp_p2, cmp[n].area);
      }
    }
    fprintf(fp, "\n");
    
    fprintf(fp, "\t                                      i_st    i_ed    j_st    j_ed    k_st    k_ed\n");
    for (int n=1; n<=NoCompo; n++)
    {
      if ( cmp[n].getType() == HEX )
      {
        st = getCmpGbbox_st(n, gci);
        ed = getCmpGbbox_ed(n, gci);
        
        fprintf(fp, "\t%3d %24s %7d %7d %7d %7d %7d %7d\n", 
                n, cmp[n].getAlias().c_str(),
                st.x, ed.x, st.y, ed.y, st.z, ed.z);
      }
    }
    fprintf(fp, "\n");
    
    fprintf(fp, "\t no                    Label   c1         c2         c3         c4  u_th[m/s]  thick[mm]     vec_forcing\n");
    for (int n=1; n<=NoCompo; n++)
    {
      if ( cmp[n].getType() == HEX )
      {
        fprintf(fp, "\t%3d %24s %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e     %s\n", 
                n, cmp[n].getAlias().c_str(),
                cmp[n].ca[0], cmp[n].ca[1], cmp[n].ca[2], cmp[n].ca[3], cmp[n].ca[4]*RefVelocity, cmp[n].ca[5]*RefLength*1000.0,
                (cmp[n].get_sw_HexDir()==ON) ? "Directional":"Non-directional");
      }
    }
  }
  
  // Fan
  if ( existComponent(FAN, cmp) )
  {
    fprintf(fp, "\n\t[Fan]\n");
    
    fprintf(fp, "\t no                    Label   vector_x   vector_y   vector_z      O_x[m]     O_y[m]     O_z[m]\n");
    for (int n=1; n<=NoCompo; n++)
    {
      if ( cmp[n].getType() == FAN )
      {
        fprintf(fp, "\t%3d %24s %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n", 
                n, cmp[n].getAlias().c_str(),
								cmp[n].nv[0], cmp[n].nv[1], cmp[n].nv[2],
                cmp[n].oc[0], cmp[n].oc[1], cmp[n].oc[2]);
      }
    }
    fprintf(fp, "\n");
    
    fprintf(fp, "\t                                     Depth[m]     Fan[m]    Boss[m]  Area[m*m]\n");
    for (int n=1; n<=NoCompo; n++)
    {
      if ( cmp[n].getType() == FAN )
      {
        fprintf(fp, "\t%3d %24s %10.3e %10.3e %10.3e %10.3e\n", 
                n, cmp[n].getAlias().c_str(),
								cmp[n].depth, cmp[n].shp_p1, cmp[n].shp_p2, cmp[n].area);
      }
    }
    fprintf(fp, "\n");
    
    fprintf(fp, "\t                                      i_st    i_ed    j_st    j_ed    k_st    k_ed\n");
    for (int n=1; n<=NoCompo; n++)
    {
      if ( cmp[n].getType() == FAN )
      {
        st = getCmpGbbox_st(n, gci);
        ed = getCmpGbbox_ed(n, gci);
        
        fprintf(fp, "\t%3d %24s %7d %7d %7d %7d %7d %7d\n", 
                n, cmp[n].getAlias().c_str(),
                st.x, ed.x, st.y, ed.y, st.z, ed.z);
      }
    }
    fprintf(fp, "\n");
    /*
     fprintf(fp, "\t no                    Label    ID         c1         c2         c3         c4  u_th[m/s]  thick[mm]     vec_forcing\n");
     for(n=1; n<=NoBC; n++) {
     if ( cmp[n].getType() == FAN ) {
     fprintf(fp, "\t%3d %24s %5d %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e     %s\n", 
     n, cmp[n].getAlias().c_str(), cmp[n].getOdrOdr(),
     cmp[n].ca[0], cmp[n].ca[1], cmp[n].ca[2], cmp[n].ca[3], cmp[n].ca[4]*RefVelocity, cmp[n].ca[5]*RefLength*1000.0,
     (cmp[n].get_sw_HexDir()==ON) ? "Directional":"Non-directional");
     }
     }*/
  }
  
  // Solid Revolution
  if ( existComponent(SOLIDREV, cmp) )
  {
    fprintf(fp, "\n\t[Solid Revolution]\n");
    
    fprintf(fp, "\t no                    Label        n_x        n_y        n_z     O_x[m]     O_y[m]     O_z[m]\n");
    for (int n=1; n<=NoCompo; n++)
    {
      if ( cmp[n].getType() == SOLIDREV )
      {
        fprintf(fp, "\t%3d %24s %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",
                n, cmp[n].getAlias().c_str(),
                cmp[n].nv[0], cmp[n].nv[1], cmp[n].nv[2],
                cmp[n].oc[0], cmp[n].oc[1], cmp[n].oc[2]);
      }
    }
    fprintf(fp, "\n");
    
    fprintf(fp, "\t                               Depth[m]  Radius[m]  Rot. Frequency[Hz]\n");
    for (int n=1; n<=NoCompo; n++)
    {
      if ( cmp[n].getType() == SOLIDREV )
      {
        fprintf(fp, "\t%3d %24s %10.3e %10.3e          %10.3e\n",
                n, cmp[n].getAlias().c_str(), cmp[n].depth, cmp[n].shp_p1, cmp[n].ca[0]);
      }
    }
    fprintf(fp, "\n");
    
    fprintf(fp, "\t                                i_st    i_ed    j_st    j_ed    k_st    k_ed\n");
    for (int n=1; n<=NoCompo; n++)
    {
      if ( cmp[n].getType() == SOLIDREV )
      {
        st = getCmpGbbox_st(n, gci);
        ed = getCmpGbbox_ed(n, gci);
        
        fprintf(fp, "\t%3d %24s %7d %7d %7d %7d %7d %7d\n",
                n, cmp[n].getAlias().c_str(),
                st.x, ed.x, st.y, ed.y, st.z, ed.z);
      }
    }
    fprintf(fp, "\n");
  }
  
    
	// Darcy Law
  if ( existComponent(DARCY, cmp) )
  {
    fprintf(fp, "\n\t[Darcy medium]\n");
    
    fprintf(fp, "\t no                    Label   Area[m*m]   Prmblty_x   Prmblty_y   Prmblty_z[m^2]    Prmblty_x   Prmblty_y   Prmblty_z[-]\n");
    for (int n=1; n<=NoCompo; n++)
    {
      if ( cmp[n].getType() == DARCY )
      {
        fprintf(fp, "\t%3d %24s %10.3e %11.4e %11.4e %11.4e       %11.4e %11.4e %11.4e\n", 
                n, cmp[n].getAlias().c_str(),
                cmp[n].area, 
								cmp[n].ca[0], cmp[n].ca[1], cmp[n].ca[2], cmp[n].ca[3], cmp[n].ca[4], cmp[n].ca[5]);
      }
    }
    
    fprintf(fp, "\t                                      i_st    i_ed    j_st    j_ed    k_st    k_ed\n");
    for (int n=1; n<=NoCompo; n++)
    {
      if ( cmp[n].getType() == DARCY )
      {
        st = getCmpGbbox_st(n, gci);
        ed = getCmpGbbox_ed(n, gci);
        
        fprintf(fp, "\t%3d %24s %7d %7d %7d %7d %7d %7d\n", 
                n, cmp[n].getAlias().c_str(),
                st.x, ed.x, st.y, ed.y, st.z, ed.z);
      }
    }
  }
  
  // Heat Face ---------------------------------------------------
  // Adiabatic
  if ( HeatProblem && existComponent(ADIABATIC, cmp) )
  {
    fprintf(fp, "\n\t[Adiabatic]\n");
    fprintf(fp, "\t no                    Label   # of faces    i_st    i_ed    j_st    j_ed    k_st    k_ed\n");
    
    for (int n=1; n<=NoCompo; n++)
    {
      if ( cmp[n].getType() == ADIABATIC )
      {
        st = getCmpGbbox_st(n, gci);
        ed = getCmpGbbox_ed(n, gci);
        
        fprintf(fp, "\t%3d %24s %12ld %7d %7d %7d %7d %7d %7d\n", n, cmp[n].getAlias().c_str(), cmp[n].getElement(),
                st.x, ed.x, st.y, ed.y, st.z, ed.z);
      }
    }
  }
  
  // Direct Heat Flux
  if ( HeatProblem && existComponent(HEATFLUX, cmp) )
  {
    fprintf(fp, "\n\t[Direct Heat Flux]\n");
    fprintf(fp, "\t no                    Label   # of faces    i_st    i_ed    j_st    j_ed    k_st    k_ed   Area[m*m]   flux(W/m^2)        q[-]\n");
    
    for (int n=1; n<=NoCompo; n++)
    {
      if ( cmp[n].getType() == HEATFLUX )
      {
        st = getCmpGbbox_st(n, gci);
        ed = getCmpGbbox_ed(n, gci);
        
        fprintf(fp, "\t%3d %24s %12ld %7d %7d %7d %7d %7d %7d %11.4e %12.4e %12.4e\n",
                n, cmp[n].getAlias().c_str(), cmp[n].getElement(),
                st.x, ed.x, st.y, ed.y, st.z, ed.z,
                cmp[n].area, cmp[n].getHeatflux(), cmp[n].getHeatflux()/(RefVelocity*DiffTemp*RefDensity*RefSpecificHeat));
      }
    }
  }
  
  // Heat Transfer S
  if ( HeatProblem && existCompoTransfer(HT_S, cmp) )
  {
    fprintf(fp, "\n\t[Heat Transfer : Type S]\n");
    fprintf(fp, "\t no                    Label   # of faces    i_st    i_ed    j_st    j_ed    k_st    k_ed   Area[m*m]    H[W/m^2K]      Temp[C]      Temp[-]\n");
    
    for (int n=1; n<=NoCompo; n++)
    {
      if ( cmp[n].getHtype() == HT_S )
      {
        st = getCmpGbbox_st(n, gci);
        ed = getCmpGbbox_ed(n, gci);
        
        fprintf(fp, "\t%3d %24s %12ld %7d %7d %7d %7d %7d %7d %11.4e %12.4e %12.4e %12.4e\n", 
                n, cmp[n].getAlias().c_str(), cmp[n].getElement(),
                st.x, ed.x, st.y, ed.y, st.z, ed.z, 
                cmp[n].area, cmp[n].getCoefHT(), cmp[n].getTemp(),
                FBUtility::convTempD2ND(cmp[n].getTemp(), BaseTemp, DiffTemp) );
      }
    }
  }
  
  // Heat Transfer SN
  if ( HeatProblem && existCompoTransfer(HT_SN, cmp) )
  {
    fprintf(fp, "\n\t[Heat Transfer : Type SN]\n"); 
    fprintf(fp, "\t no                    Label   # of faces    i_st    i_ed    j_st    j_ed    k_st    k_ed   Area[m*m]      Temp[C]      Temp[-]   Type\n");
    
    for (int n=1; n<=NoCompo; n++)
    {
      if ( cmp[n].getHtype() == HT_SN )
      {
        st = getCmpGbbox_st(n, gci);
        ed = getCmpGbbox_ed(n, gci);
        
        fprintf(fp, "\t%3d %24s %12ld %7d %7d %7d %7d %7d %7d %11.4e %12.4e %12.4e\n",
                n, cmp[n].getAlias().c_str(), cmp[n].getElement(),
                st.x, ed.x, st.y, ed.y, st.z, ed.z,
                cmp[n].area, cmp[n].getTemp(),
                FBUtility::convTempD2ND(cmp[n].getTemp(), BaseTemp, DiffTemp));
      }
    }
    fprintf(fp, "\n");
    
    fprintf(fp, "\t no                    Label   vert_lam_a  vert_lam_b vert_turb_a vert_turb_b  vert_Ra_cr   lwr_lam_a   lwr_lam_b  lwr_turb_a  lwr_turb_b   lwr_Ra_cr\n");
    for (int n=1; n<=NoCompo; n++)
    {
      if ( cmp[n].getHtype() == HT_SN )
      {
        fprintf(fp, "\t%3d %24s %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e\n", 
                n, 
                cmp[n].getAlias().c_str(),
                cmp[n].ca[CompoList::vert_laminar_alpha], 
                cmp[n].ca[CompoList::vert_laminar_beta], 
                cmp[n].ca[CompoList::vert_turbulent_alpha], 
                cmp[n].ca[CompoList::vert_turbulent_beta], 
                cmp[n].ca[CompoList::vert_Ra_critial], 
                cmp[n].cb[CompoList::lower_laminar_alpha], 
                cmp[n].cb[CompoList::lower_laminar_beta], 
                cmp[n].cb[CompoList::lower_turbulent_alpha], 
                cmp[n].cb[CompoList::lower_turbulent_beta],
                cmp[n].cb[CompoList::lower_Ra_critial] );
      }
    }
  }
  
  // Heat Transfer SF
  if ( HeatProblem && existCompoTransfer(HT_SF, cmp) )
  {
    fprintf(fp, "\n\t[Heat Transfer : Type SF]\n");
    fprintf(fp, "\t no                    Label   # of faces    i_st    i_ed    j_st    j_ed    k_st    k_ed   Area[m*m]      Temp[C]      Temp[-]   Type\n");
    
    for (int n=1; n<=NoCompo; n++)
    {
      if ( cmp[n].getHtype() == HT_SF )
      {
        st = getCmpGbbox_st(n, gci);
        ed = getCmpGbbox_ed(n, gci);
        
        fprintf(fp, "\t%3d %24s %12ld %7d %7d %7d %7d %7d %7d %11.4e %12.4e %12.4e\n",
                n, 
                cmp[n].getAlias().c_str(), 
                cmp[n].getElement(),
                st.x, ed.x, st.y, ed.y, st.z, ed.z,
                cmp[n].area, cmp[n].getTemp(),
                FBUtility::convTempD2ND(cmp[n].getTemp(), BaseTemp, DiffTemp));
      }
    }
    fprintf(fp, "\t no                    Label   alpha        beta       gamma\n");
    for (int n=1; n<=NoCompo; n++)
    {
      if ( cmp[n].getHtype() == HT_SN )
      {
        fprintf(fp, "\t%3d %24s %11.4e %11.4e %11.4e\n", 
                n, 
                cmp[n].getAlias().c_str(), 
                cmp[n].ca[CompoList::alpha], 
                cmp[n].ca[CompoList::beta], 
                cmp[n].ca[CompoList::gamma]);
      }
    }
  }
  
  // Iso Thermal
  if ( HeatProblem && existComponent(ISOTHERMAL, cmp))
  {
    fprintf(fp, "\n\t[Iso-Thermal]\n");
    fprintf(fp, "\t no                    Label   # of faces    i_st    i_ed    j_st    j_ed    k_st    k_ed   Area[m*m]   Sf.Temp[C]   Sf.Temp[-]\n");
    
    for (int n=1; n<=NoCompo; n++)
    {
      if ( cmp[n].getType() == ISOTHERMAL )
      {
        st = getCmpGbbox_st(n, gci);
        ed = getCmpGbbox_ed(n, gci);
        
        fprintf(fp, "\t%3d %24s %12ld %7d %7d %7d %7d %7d %7d %11.4e %12.4e %12.4e \n", 
                n, cmp[n].getAlias().c_str(), cmp[n].getElement(),
                st.x, ed.x, st.y, ed.y, st.z, ed.z,
								cmp[n].area, cmp[n].getTemp(),
                FBUtility::convTempD2ND(cmp[n].getTemp(), BaseTemp, DiffTemp) );
      }
    }
  }
  
  // Radiant
  if ( HeatProblem && existComponent(RADIANT, cmp))
  {
    fprintf(fp, "\n\t[Radiant]\n");
    fprintf(fp, "\t no                    Label   # of faces    i_st    i_ed    j_st    j_ed    k_st    k_ed   Area[m*m]   ep[-]   pj[-]\n");
    
    for (int n=1; n<=NoCompo; n++)
    {
      if ( cmp[n].getType() == RADIANT )
      {
        st = getCmpGbbox_st(n, gci);
        ed = getCmpGbbox_ed(n, gci);
        
        fprintf(fp, "\t%3d %24s %12ld %7d %7d %7d %7d %7d %7d %11.4e %12.4e %12.4e\n", 
                n, cmp[n].getAlias().c_str(), cmp[n].getElement(),
                st.x, ed.x, st.y, ed.y, st.z, ed.z,
                cmp[n].area, cmp[n].get_CoefRadEps(), cmp[n].get_CoefRadPrj());
      }
    }
  }
  
  // Heat source ---------------------------------------------------
  // Heat generation
  if ( HeatProblem && existComponent(HEAT_SRC, cmp))
  {
    fprintf(fp, "\n\t[Heat Generation]\n");
    fprintf(fp, "\t no                    Label   # of faces    i_st    i_ed    j_st    j_ed    k_st    k_ed     Q-Vol[W]     Q[W/m^3]  Q-nrmlzd[-]\n");
    
    for (int n=1; n<=NoCompo; n++)
    {
      if ( cmp[n].getType() == HEAT_SRC )
      {
        st = getCmpGbbox_st(n, gci);
        ed = getCmpGbbox_ed(n, gci);
        
        fprintf(fp, "\t%3d %24s %12ld %7d %7d %7d %7d %7d %7d %12.4e %12.4e %12.4e\n", 
                n, cmp[n].getAlias().c_str(), cmp[n].getElement(),
                st.x, ed.x, st.y, ed.y, st.z, ed.z,
                cmp[n].getHeatValue(),
                cmp[n].getHeatDensity(),
                FBUtility::convHsrcD2ND(cmp[n].getHeatDensity(), RefVelocity, RefLength, DiffTemp, mat[n].P[p_density], mat[n].P[p_specific_heat]));
      }
    }
  }
  
  // Constant Temperature
  if ( HeatProblem && existComponent(CNST_TEMP, cmp))
  {
    fprintf(fp, "\n\t[Constant Temperature]\n");
    fprintf(fp, "\t no                    Label   # of faces    i_st    i_ed    j_st    j_ed    k_st    k_ed      Temp[C]      Temp[-]\n");
    
    for (int n=1; n<=NoCompo; n++)
    {
      if ( cmp[n].getType() == CNST_TEMP )
      {
        st = getCmpGbbox_st(n, gci);
        ed = getCmpGbbox_ed(n, gci);
        
        fprintf(fp, "\t%3d %24s %12ld %7d %7d %7d %7d %7d %7d %12.4e %12.4e\n", 
                n, cmp[n].getAlias().c_str(), cmp[n].getElement(),
                st.x, ed.x, st.y, ed.y, st.z, ed.z,
                cmp[n].getTemp(), 
                FBUtility::convTempD2ND(cmp[n].getTemp(), BaseTemp, DiffTemp));
      }
    }
  }

  
  // Obstacle ---------------------------------------------------
  if ( existComponent(OBSTACLE, cmp) )
  {
    fprintf(fp, "\n\t[Obstacle]\n");
    fprintf(fp, "\t no                    Label    i_st    i_ed    j_st    j_ed    k_st    k_ed   Area[m*m]    Elements\n");
    
    for (int n=1; n<=NoCompo; n++)
    {
      if ( cmp[n].getType() == OBSTACLE )
      {
        st = getCmpGbbox_st(n, gci);
        ed = getCmpGbbox_ed(n, gci);
        
        fprintf(fp, "\t%3d %24s %7d %7d %7d %7d %7d %7d %11.4e %11ld\n",
                n, cmp[n].getAlias().c_str(),
                st.x, ed.x, st.y, ed.y, st.z, ed.z,
								cmp[n].area, cmp[n].getElement());
      }
    }
    fprintf(fp, "\n");
  }
  
  
  // Periodic ---------------------------------------------------
  if ( existComponent(PERIODIC, cmp) )
  {
    fprintf(fp, "\n\t[Periodic]\n");
    fprintf(fp, "\t no                    Label   i_st    i_ed    j_st    j_ed    k_st    k_ed    Pressure Difference [Pa]/[-]  Driver\n");
    
    int dir_in=0, dir_out=0, pp_in=0, pp_out=0;
    
    for (int n=1; n<=NoCompo; n++)
    {
      if ( cmp[n].getType() == PERIODIC )
      {
        dir_in = cmp[n].getPeriodicDir();
        
        st = getCmpGbbox_st(n, gci);
        ed = getCmpGbbox_ed(n, gci);
        
        fprintf(fp, "\t%3d %24s %7d %7d %7d %7d %7d %7d     ", 
                n, cmp[n].getAlias().c_str(),
                st.x, ed.x, st.y, ed.y, st.z, ed.z);
        fprintf(fp,"%12.6e / %12.6e ", cmp[n].ca[0], FBUtility::convPrsD2ND(cmp[n].ca[0], BasePrs, RefDensity, RefVelocity, Unit_Prs));
        fprintf(fp, "%7s\n", FBUtility::getDirection(dir_in).c_str());
        
        // ドライバの方向が周期境界であるかをチェック
        if ( bc[dir_in].getClass() != OBC_PERIODIC )
        {
          fprintf(fp, "\tError : Specified driver direction[%s] by component is not PERIODIC.", FBUtility::getDirection(dir_in).c_str());
          Exit(0);
        }
        
        // ドライバの方向が一致しているかどうかをチェック
        dir_out = bc[dir_in].get_DriverDir();
        if ( dir_in != dir_out )
        {
          fprintf(fp, "\tError : The specification of driver direction is different between Outer(%d) and Inner(%d) in XML.", dir_out, dir_in);
          Exit(0);
        }
        
        // 入力のインデクス値とボクセルの位置が一致しているかどうかをチェック
        switch ( dir_in )
        {
          case X_minus:
          case X_plus:
            pp_in = getCmpGbbox_st(n, gci).x;
            break;
            
          case Y_minus:
          case Y_plus:
            pp_in = getCmpGbbox_st(n, gci).y;
            break;
            
          case Z_minus:
          case Z_plus:
            pp_in = getCmpGbbox_st(n, gci).z;
            break;
        }
        
        pp_out = bc[dir_in].get_DriverIndex();
        if ( pp_in != pp_out )
        {
          fprintf(fp, "\tError : The index of inner periodic cell is different in between XML(%d) and voxel model(%d).", pp_out, pp_in);
          Exit(0);
        }
      }
    }
    fprintf(fp, "\n");
  }
  
  fflush(fp);
}


// #################################################################
// 取得したcmpList[]の内容を表示する
void ParseBC::printCompoSummary(FILE* fp, CompoList* cmp, const int basicEq)
{
  if( !fp ) Exit(0);
  
  if ( basicEq == INCMP_2PHASE )
  {
    fprintf(fp,"\t  No :  # of Faces/Cells       State   Phase                     Label : Class                                 Medium\n");
    fprintf(fp,"\t------------------------------------------------------------------------------------------------------------------------------\n");
    
    for (int i=1; i<=NoCompo; i++)
    {
      fprintf(fp,"\t%4d : %18ld ", i, cmp[i].getElement());
      ( cmp[i].getState() == FLUID ) ? fprintf(fp, "      Fluid ") : fprintf(fp, "      Solid ") ;
      ( cmp[i].getPhase() == GAS )   ? fprintf(fp, "        Gas ") : fprintf(fp, "     Liquid ") ;
      fprintf(fp, " %24s : %-36s  %-20s", (cmp[i].getAlias().empty()) ? "" : cmp[i].getAlias().c_str(),
              cmp[i].getBCstr().c_str(),
              cmp[i].getMedium().c_str() );
    }
  }
  else
  {
    if ( KindOfSolver == FLOW_ONLY )
    {
      fprintf(fp,"\t  No :   # of Faces/Cells       State                    Label : Class                                 Medium\n");
      fprintf(fp,"\t------------------------------------------------------------------------------------------------------------------------------\n");
      
      for (int i=1; i<=NoCompo; i++)
      {
        fprintf(fp,"\t%4d : %18ld ", i, cmp[i].getElement());
        ( cmp[i].getState() == FLUID ) ? fprintf(fp, "      Fluid ") : fprintf(fp, "      Solid ") ;
        fprintf(fp, "%24s : %-36s  %-20s", (cmp[i].getAlias().empty()) ? "" : cmp[i].getAlias().c_str(),
                cmp[i].getBCstr().c_str(),
                cmp[i].getMedium().c_str() );
        fprintf(fp,"\n");
      }
    }
    else
    {
      fprintf(fp,"\t  No :   # of Faces/Cells       State   Init.Temp[C]                    Label : Class                                 Medium\n");
      fprintf(fp,"\t------------------------------------------------------------------------------------------------------------------------------\n");
      
      for (int i=1; i<=NoCompo; i++)
      {
        fprintf(fp,"\t%4d : %18ld ", i, cmp[i].getElement());
        ( cmp[i].getState() == FLUID ) ? fprintf(fp, "      Fluid ") : fprintf(fp, "      Solid ") ;
        fprintf(fp, "%14.4e %24s : %-36s  %-20s",
                cmp[i].getInitTemp(),
                (cmp[i].getAlias().empty()) ? "" : cmp[i].getAlias().c_str(),
                cmp[i].getBCstr().c_str(),
                cmp[i].getMedium().c_str() );
        fprintf(fp,"\n");
      }
    }
  }
  fprintf(fp,"\n");
}



// #################################################################
// 外部境界条件の各面の情報を表示する
void ParseBC::printFaceOBC(FILE* fp, const REAL_TYPE* G_reg, const BoundaryOuter* bc, const MediumList* mat)
{
  for (int i=0; i<NOFACE; i++)
  {
    fprintf(fp,"\t      Set %s up as %s : < %s >\n", 
            FBUtility::getDirection(i).c_str(), 
            getOBCstr(bc[i].getClass()).c_str(),
            bc[i].getAlias().c_str());
    
    printOBC(fp, &bc[i], mat, G_reg, i);
    
    fprintf(fp,"\n");
  }
  fflush(fp);
}



// #################################################################
/**
 * @brief 速度の外部境界条件処理の表示
 * @param [in] fp    ファイルポインタ
 * @param [in] ref   BoundaryOuter
 * @param [in] mat   MediumList
 * @param [in] G_reg グローバルの領域の大きさ
 * @param [in] face  面番号
 */
void ParseBC::printOBC(FILE* fp, const BoundaryOuter* ref, const MediumList* mat, const REAL_TYPE* G_reg, const int face)
{
  REAL_TYPE a, b, c;
  
  if ( ref->getClass() == OBC_INTRINSIC)
  {
    fprintf(fp,"\t\t\tGuide Cell Medium = Intrinsic\n");
  }
  else if ( ref->getClass() == OBC_PERIODIC)
  {
    ; // 表示なし
  }
  else
  {
    fprintf(fp,"\t\t\tGuide Cell Medium = %s\n", mat[ref->getGuideMedium()].getAlias().c_str());
  }
  
  
  switch ( ref->getClass() )
  {
    case OBC_WALL:
      if ( ref->get_V_Profile() == CompoList::vel_harmonic )
      {
        fprintf(fp,"\t\t\tVelocity  Harmonic Oscillation\n");
        fprintf(fp,"\t\t\t          nrml vec. V(%10.3e, %10.3e, %10.3e) [-]\n", ref->nv[0], ref->nv[1], ref->nv[2]);
        fprintf(fp,"\t\t\t          amp.    %12.6e [m/s] / %12.6e [-]\n", ref->ca[0], ref->ca[0]/RefVelocity);
        fprintf(fp,"\t\t\t          freq.   %12.6e [Hz]  / %12.6e [-]\n", ref->ca[1], ref->ca[1]*RefLength/RefVelocity);
        fprintf(fp,"\t\t\t          phase   %12.6e [rad] / %12.6e [-]\n", ref->ca[2], ref->ca[2]);
        fprintf(fp,"\t\t\t          intcpt. %12.6e [m/s] / %12.6e [-]\n", ref->ca[3], ref->ca[3]/RefVelocity);
      }
      else // vel_constant, vel_zero
      {
        c = ref->ca[CompoList::bias]; // 有次元値
        fprintf(fp,"\t\t\tVelocity V(%10.3e, %10.3e, %10.3e) [m/s] / V(%10.3e, %10.3e, %10.3e) [-]\n", 
                ref->nv[0]*c, ref->nv[1]*c, ref->nv[2]*c, ref->nv[0]*c/RefVelocity, ref->nv[1]*c/RefVelocity, ref->nv[2]*c/RefVelocity); 
      }
      
      fprintf(fp,"\t\t\tPressure Gradient is %s\n", ( Mode_Gradp == P_GRAD_ZERO) ? "zero" : "calculated from Navier-Stokes eqs.");
      
      if ( HeatProblem )
      {
        int htp = ref->getHtype();
        
        if ( htp == ADIABATIC )
        {
          fprintf(fp,"\t\t\tAdiabatic\n");
        }
        else if ( htp == TRANSFER)
        {
          int ht_mode = ref->getHTmode();
          
          if ( ht_mode == HT_S )
          {
            fprintf(fp, "\tHeat Transfer Type S  : H. T. Coef. = %e \n", ref->getCoefHT());
            fprintf(fp, "\t                        Surf. temp. = %e \n", ref->getTemp());
          }
          else if ( ht_mode == HT_SN )
          {
            fprintf(fp, "\tHeat Transfer Type SN : Surf. temp. = %e \n", ref->getTemp());
            fprintf(fp, "\t                        Ref. Temp.  = %s \n", (ref->getHTmodeRef()==CompoList::HT_mode_bulk) ? "Bulk" : "Local");
            fprintf(fp, "\t                        vertival_laminar_alpha   = %12.6e \n", ref->ca[0]);
            fprintf(fp, "\t                        vertival_laminar_beta    = %12.6e \n", ref->ca[1]);
            fprintf(fp, "\t                        vertival_turbulent_alpha = %12.6e \n", ref->ca[2]);
            fprintf(fp, "\t                        vertival_turbulent_beta  = %12.6e \n", ref->ca[3]);
            fprintf(fp, "\t                        vertival_Ra_critial      = %12.6e \n", ref->ca[4]);
            fprintf(fp, "\t                        lower_laminar_alpha      = %12.6e \n", ref->cb[0]);
            fprintf(fp, "\t                        lower_laminar_beta       = %12.6e \n", ref->cb[1]);
            fprintf(fp, "\t                        lower_turbulent_alpha    = %12.6e \n", ref->cb[2]);
            fprintf(fp, "\t                        lower_turbulent_beta     = %12.6e \n", ref->cb[3]);
            fprintf(fp, "\t                        lower_Ra_critial         = %12.6e \n", ref->cb[4]);
          }
          else if ( ht_mode == HT_SF )
          {
            fprintf(fp, "\tHeat Transfer Type SF : Surf. temp. = %e \n", ref->getTemp());
            fprintf(fp, "\t                        Ref. Temp.  = %s \n", (ref->getHTmodeRef()==CompoList::HT_mode_bulk) ? "Bulk" : "Local");
            fprintf(fp, "\t                        alpha       = %12.6e \n", ref->ca[0]);
            fprintf(fp, "\t                        beta        = %12.6e \n", ref->ca[1]);
            fprintf(fp, "\t                        gamma       = %12.6e \n", ref->ca[2]);
          }
        }
        else if ( htp == HEATFLUX )
        {
          fprintf(fp,"\t\t\tHeat Flux  = %12.6e [W/m^2] / %12.6e [-]\n", ref->getHeatflux(), ref->getHeatflux());
        }
        else if ( htp == ISOTHERMAL )
        {
          fprintf(fp,"\t\t\tIsothermal   = %12.6e [C] / %12.6e [-]\n",
                  ref->getTemp(),
                  FBUtility::convTempD2ND(ref->getTemp(), BaseTemp, DiffTemp));
        }
        //else if ( htp == CNST_TEMP ) {
        //  fprintf(fp,"Dirichlet %e [C]\n", ref->T1.c);
        //}        
      }
      break;
      
    case OBC_SPEC_VEL:
      if ( ref->get_V_Profile() == CompoList::vel_harmonic )
      {
        fprintf(fp,"\t\t\tVelocity  Harmonic Oscillation\n");
        fprintf(fp,"\t\t\t          nrml vec. V(%10.3e, %10.3e, %10.3e) [-]\n", ref->nv[0], ref->nv[1], ref->nv[2]);
        fprintf(fp,"\t\t\t          amp.    %12.6e [m/s] / %12.6e [-]\n", ref->ca[0], ref->ca[0]/RefVelocity);
        fprintf(fp,"\t\t\t          freq.   %12.6e [Hz]  / %12.6e [-]\n", ref->ca[1], ref->ca[1]*RefLength/RefVelocity);
        fprintf(fp,"\t\t\t          phase   %12.6e [rad] / %12.6e [-]\n", ref->ca[2], ref->ca[2]);
        fprintf(fp,"\t\t\t          intcpt. %12.6e [m/s] / %12.6e [-]\n", ref->ca[3], ref->ca[3]/RefVelocity);
      }
      else // vel_zero, vel_constant
      {
        c = ref->ca[CompoList::bias];
        fprintf(fp,"\t\t\tVelocity V(%10.3e, %10.3e, %10.3e) [m/s] / V(%10.3e, %10.3e, %10.3e) [-]\n", 
                ref->nv[0]*c, ref->nv[1]*c, ref->nv[2]*c, ref->nv[0]*c/RefVelocity, ref->nv[1]*c/RefVelocity, ref->nv[2]*c/RefVelocity);
      }
      
      fprintf(fp,"\t\t\tPressure Gradient is zero.\n");
      
      if ( HeatProblem )
      {
        fprintf(fp, "\t\t\tSpecified Temperature  = %12.6e [C] / %12.6e [-] \n", 
                ref->getTemp(),
                FBUtility::convTempD2ND(ref->getTemp(), BaseTemp, DiffTemp));
      }
      break;
      
      
    case OBC_OUTFLOW:
      if (ref->get_pType() == P_DIRICHLET)
      {
        fprintf(fp,"\t\t\t%12.6e [Pa]  /  %12.6e\n", ref->p, FBUtility::convPrsD2ND(ref->p, BasePrs, RefDensity, RefVelocity, Unit_Prs));        
      }
      else
      {
        fprintf(fp,"\t\t\tPressure Gradient is zero\n");
      }
      if ( HeatProblem )
      {
        // ?
      }
      fprintf(fp,"\t\t\tNumber of Element = %d\n",ref->getValidCell());
      break;
      
      
    case OBC_TRC_FREE:
      fprintf(fp,"\t\t\t%12.6e [Pa]  /  %12.6e [-]\n", ref->p, FBUtility::convPrsD2ND(ref->p, BasePrs, RefDensity, RefVelocity, Unit_Prs));
      
      if ( HeatProblem )
      {
        fprintf(fp, "\t\t\t    Ambient Temperature  = %12.6e \n", ref->getTemp());
      }
      break;
      
      
    case OBC_PERIODIC:
      switch (ref->getPrdcMode())
      {
        case BoundaryOuter::prdc_Simple:
          fprintf(fp, "\t\t\tSimple periodic copy\n");
          break;
          
        case BoundaryOuter::prdc_Directional:
          fprintf(fp,"\t\t\tPressure Difference = %12.6e [Pa]    /  %12.6e [-]\n", 
                  ref->p, FBUtility::convPrsD2ND(ref->p, BasePrs, RefDensity, RefVelocity, Unit_Prs));
          
          switch (face)
          {
            case X_minus:
            case X_plus:
              fprintf(fp,"\t\t\tPressure Gradient   = %12.6e [Pa/m]  /  %12.6e [-]\n", 
                      ref->p/(G_reg[0]*RefLength), FBUtility::convPrsD2ND(ref->p, BasePrs, RefDensity, RefVelocity, Unit_Prs)/G_reg[0]);
              break;
              
            case Y_minus:
            case Y_plus:
              fprintf(fp,"\t\t\tPressure Gradient   = %12.6e [Pa/m]  /  %12.6e [-]\n", 
                      ref->p/(G_reg[1]*RefLength), FBUtility::convPrsD2ND(ref->p, BasePrs, RefDensity, RefVelocity, Unit_Prs)/G_reg[1]);
              break;
              
            case Z_minus:
            case Z_plus:
              fprintf(fp,"\t\t\tPressure Gradient   = %12.6e [Pa/m]  /  %12.6e [-]\n", 
                      ref->p/(G_reg[2]*RefLength), FBUtility::convPrsD2ND(ref->p, BasePrs, RefDensity, RefVelocity, Unit_Prs)/G_reg[2]);
              break;
          }          
          break;
          
        case BoundaryOuter::prdc_Driver:
          fprintf(fp,"\t\t\tDriver direction    : %s\n", FBUtility::getDirection(ref->get_DriverDir()).c_str());
          break;
          
        default:
          Exit(0);
      }
      break;
      
      
    case OBC_FAR_FIELD:
      fprintf(fp, "\t\t\tVelocity : Neumann, Pressure : %s\n", (ref->get_pType() == P_DIRICHLET) ? "Dirichlet" : "Neumann");
      if (ref->get_pType() == P_DIRICHLET)
      {
        fprintf(fp,"\t\t\t%12.6e [Pa]  /  %12.6e [-]\n", ref->p, FBUtility::convPrsD2ND(ref->p, BasePrs, RefDensity, RefVelocity, Unit_Prs));
      }
      else
      {
        fprintf(fp,"\t\t\tPressure Gradient is zero\n");
      }
      break;
      
    // nothing to do
    case OBC_SYMMETRIC:
    case OBC_INTRINSIC:
      break;
      
      
    default:
      printf("\n\tError : OuterBC\n");
      Exit(0);
      break;
  }
  
  fflush(fp);
}



// #################################################################
// 必要な変数をセットする
void ParseBC::setControlVars(Control* Cref)
{
  RefVelocity = Cref->RefVelocity;
  BaseTemp    = Cref->BaseTemp;
  DiffTemp    = Cref->DiffTemp;
  RefDensity  = Cref->RefDensity;
  HeatProblem = Cref->isHeatProblem();
  Unit_Param  = Cref->Unit.Param;
  RefLength   = Cref->RefLength;
  KindOfSolver= Cref->KindOfSolver;
  Unit_Prs    = Cref->Unit.Prs;
	BasePrs     = Cref->BasePrs;
  Mode_Gradp  = Cref->Mode.PrsNeuamnnType;
  NoCompo     = Cref->NoCompo;
  NoBC        = Cref->NoBC;
  NoMedium    = Cref->NoMedium;
  
  int m;
  double s, two=2.0;
  
  s = (double)MASK_5; // bit幅マスクは2^(bit幅)-1を表し，ちょうど0を除いた個数となっている
  m = (int)(log10(s+1.0)/log10(two) );
  
  if ( NoCompo > s )
  {
    printf("Error : No. of Component (= NoBC + NoMedium) must be less or equal %d(%dbit-width)\n", (int)s, m);
    Exit(0);
  }
  
  s = (double)(MASK_5-1); // 0と31を除く
  m = (int)( log10(s+2.0)/log10(two) );
  
  if ( NoBC > s )
  {
    printf("Error : No. of BC must be less or equal %d(%dbit-width)\n", (int)s, m);
    Exit(0);
  }
  
  // 外部境界条件の種類数を取得し，内部保持配列をインスタンス
  string label;
  string str;
  
  label = "/BcTable/OuterBoundary";
  
  if ( !(tpCntl->chkNode(label)) )
  {
    Hostonly_ stamped_printf("\tParsing error : Missing '%s\n", label.c_str());
    Exit(0);
  }
  
  vector<string> nodes;
  
  // OuterBoundary直下のラベルを取得
  tpCntl->getLabelVector(label, nodes);
  
  // ラベル数から"FaceBC"を除いた数
  NoBaseBC = nodes.size() -1;
  BaseBc = new BoundaryOuter[NoBaseBC+1];
}



// #################################################################
/**
 * @brief 内部境界条件の照合を行う
 * @param [in]  keyword テストキーワード
 * @param [in]  m       BaseBcの格納番号
 * @param [out] cmp     CompoList
 */
void ParseBC::setKeywordLBC(const string keyword, const int m, CompoList* cmp)
{
  if     ( FBUtility::compare(keyword, "Adiabatic") )            cmp[m].setType(ADIABATIC);
  else if( FBUtility::compare(keyword, "DirectHeatFlux") )       cmp[m].setType(HEATFLUX);
  else if( FBUtility::compare(keyword, "HeatTransferS") ) {      cmp[m].setType(TRANSFER); cmp[m].setHtype(HT_S); }
  else if( FBUtility::compare(keyword, "HeatTransferSF") ){      cmp[m].setType(TRANSFER); cmp[m].setHtype(HT_SF); }
  else if( FBUtility::compare(keyword, "HeatTransferSN") ){      cmp[m].setType(TRANSFER); cmp[m].setHtype(HT_SN); }
  else if( FBUtility::compare(keyword, "IsoThermal") )           cmp[m].setType(ISOTHERMAL);
  else if( FBUtility::compare(keyword, "Radiation") )            cmp[m].setType(RADIANT);
  else if( FBUtility::compare(keyword, "SpecifiedVelocity") )    cmp[m].setType(SPEC_VEL);
  else if( FBUtility::compare(keyword, "Outflow") )              cmp[m].setType(OUTFLOW);
  else if( FBUtility::compare(keyword, "Forcing") )              cmp[m].setType(IBM_DF);
  else if( FBUtility::compare(keyword, "HeatSource") )           cmp[m].setType(HEAT_SRC);
  else if( FBUtility::compare(keyword, "SpecifiedTemperature") ) cmp[m].setType(CNST_TEMP);
  else if( FBUtility::compare(keyword, "PressureLoss") )         cmp[m].setType(HEX);
  else if( FBUtility::compare(keyword, "Fan") )                  cmp[m].setType(FAN);
  else if( FBUtility::compare(keyword, "Darcy") )                cmp[m].setType(DARCY);
  else if( FBUtility::compare(keyword, "Periodic") )             cmp[m].setType(PERIODIC);
  else if( FBUtility::compare(keyword, "Obstacle") )             cmp[m].setType(OBSTACLE);
  else if( FBUtility::compare(keyword, "Monitor") )              cmp[m].setType(MONITOR);
  else if( FBUtility::compare(keyword, "SolidRevolution") )      cmp[m].setType(SOLIDREV);
  else {
    Hostonly_ stamped_printf("\tInvalid keyword is described '%s'\n", keyword.c_str());
    Exit(0);
  }
}


// #################################################################
/**
 * @brief 外部境界条件の照合を行う
 * @param [in] keyword テストキーワード
 * @param [in] m       コンポーネントの格納番号
 */
void ParseBC::setKeywordOBC(const string keyword, const int m)
{
  if     ( FBUtility::compare(keyword, "Wall") )              BaseBc[m].setClass(OBC_WALL);
  else if( FBUtility::compare(keyword, "Outflow") )           BaseBc[m].setClass(OBC_OUTFLOW);
  else if( FBUtility::compare(keyword, "SpecifiedVelocity") ) BaseBc[m].setClass(OBC_SPEC_VEL);
  else if( FBUtility::compare(keyword, "Symmetric") )         BaseBc[m].setClass(OBC_SYMMETRIC);
  else if( FBUtility::compare(keyword, "Periodic") )          BaseBc[m].setClass(OBC_PERIODIC);
  else if( FBUtility::compare(keyword, "TractionFree") )      BaseBc[m].setClass(OBC_TRC_FREE);
  else if( FBUtility::compare(keyword, "FarField") )          BaseBc[m].setClass(OBC_FAR_FIELD);
  else if( FBUtility::compare(keyword, "intrinsic") )         BaseBc[m].setClass(OBC_INTRINSIC);
  else {
    Hostonly_ stamped_printf("\tParsing error : Invalid keyword is described '%s'\n", keyword.c_str());
    Exit(0);
  }
}


// #################################################################
// 参照物理量を設定する
void ParseBC::setRefMediumProperty(const REAL_TYPE m_rho, const REAL_TYPE m_cp)
{  
  RefDensity      = m_rho;
  RefSpecificHeat = m_cp;
}
